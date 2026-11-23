#ifndef REGION_READ_CACHE_H
#define REGION_READ_CACHE_H

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "sam_utils.h"

using region_read_interval_t = std::pair<hts_pos_t, hts_pos_t>;

struct cached_region_read_t {
    bam1_t read;
    size_t data_offset; // An offset stays valid when the packed read data grows.
    hts_pos_t end, prefix_max_end;
};

inline void merge_region_read_windows(std::vector<region_read_interval_t>& windows) {
    constexpr hts_pos_t MAX_MERGE_GAP = 1024;
    std::sort(windows.begin(), windows.end());
    size_t n_merged = 0;
    for (region_read_interval_t window : windows) {
        if (n_merged > 0 && window.first <= windows[n_merged-1].second + MAX_MERGE_GAP) {
            windows[n_merged-1].second = std::max(windows[n_merged-1].second, window.second);
        } else {
            windows[n_merged++] = window;
        }
    }
    windows.resize(n_merged);
}

inline bool region_read_windows_cover(const std::vector<region_read_interval_t>& windows, hts_pos_t begin, hts_pos_t end) {
    auto next = std::upper_bound(windows.begin(), windows.end(), region_read_interval_t(begin, HTS_POS_MAX));
    return next != windows.begin() && std::prev(next)->second >= end;
}

struct region_read_cache_t {
    open_samFile_t* bam_file;
    int tid;
    std::vector<region_read_interval_t> coverage;
    std::vector<cached_region_read_t> records;
    std::vector<uint8_t> read_data;
};

struct region_read_itr_t {
    open_samFile_t* bam_file;
    const region_read_cache_t* cache = nullptr;
    std::vector<region_read_interval_t> regions;
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> hts_iter{nullptr, &hts_itr_destroy};
    size_t region_idx = 0, record_idx = 0;
};

inline region_read_cache_t*& active_region_read_cache() {
    static thread_local region_read_cache_t* cache = nullptr;
    return cache;
}

inline void release_region_read_cache(region_read_cache_t* cache) {
    if (active_region_read_cache() == cache) active_region_read_cache() = nullptr;
    delete cache;
}

using region_read_cache_ptr = std::unique_ptr<region_read_cache_t, decltype(&release_region_read_cache)>;

// Decode the breakpoint windows once, before a pass restricts the CRAM fields.
// The caller owns the cache for the duration of the task.
template<typename Variant>
region_read_cache_ptr load_region_read_cache(open_samFile_t* bam_file, const std::string& contig_name,
                                       const std::vector<Variant*>& calls, hts_pos_t padding) {
    int tid = sam_hdr_name2tid(bam_file->header, contig_name.c_str());
    region_read_cache_ptr cache(nullptr, &release_region_read_cache);
    if (tid < 0 || calls.empty()) return cache;

    cache.reset(new region_read_cache_t);
    cache->bam_file = bam_file;
    cache->tid = tid;
    hts_pos_t contig_len = sam_hdr_tid2len(bam_file->header, tid);
    for (Variant* call : calls) {
        for (auto pos : {call->start, call->end}) {
            hts_pos_t begin = std::max<hts_pos_t>(0, pos-padding);
            hts_pos_t end = std::min(contig_len, pos+padding+1);
            if (begin < end) cache->coverage.emplace_back(begin, end);
        }
    }
    if (cache->coverage.empty()) return {nullptr, &release_region_read_cache};
    merge_region_read_windows(cache->coverage);

    std::vector<std::string> region_names;
    for (const region_read_interval_t& interval : cache->coverage) {
        region_names.push_back(contig_name + ":" + std::to_string(interval.first+1) + "-" + std::to_string(interval.second));
    }
    std::vector<char*> regions;
    for (std::string& name : region_names) regions.push_back(&name[0]);
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iter(
        sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size()), &hts_itr_destroy);
    if (!iter) return {nullptr, &release_region_read_cache};

    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), &bam_destroy1);
    if (!read) throw std::bad_alloc();
    constexpr size_t MAX_CACHE_BYTES = 48u << 20;
    int status;
    while ((status = sam_itr_next(bam_file->file, iter.get(), read.get())) >= 0) {
        size_t data_offset = cache->read_data.size();
        size_t cache_bytes = data_offset + read->l_data + (cache->records.size()+1)*sizeof(cached_region_read_t);
        // Discard partial coverage. The caller will query HTSlib directly at extreme depth.
        if (cache_bytes > MAX_CACHE_BYTES) return {nullptr, &release_region_read_cache};

        hts_pos_t end = bam_endpos(read.get());
        hts_pos_t prefix_max_end = cache->records.empty() ? end : std::max(end, cache->records.back().prefix_max_end);
        cached_region_read_t record = {*read, data_offset, end, prefix_max_end};
        record.read.data = nullptr;
        cache->records.push_back(record);
        cache->read_data.insert(cache->read_data.end(), read->data, read->data+read->l_data);
    }
    if (status < -1) throw std::runtime_error("Failed to read genotype block on " + contig_name);
    active_region_read_cache() = cache.get();
    return cache;
}

inline void seek_read_region(region_read_itr_t& iter) {
    // A long read may start before the region. Skip only reads known to end at or before its start.
    const auto& records = iter.cache->records;
    auto first = std::upper_bound(records.begin()+iter.record_idx, records.end(), iter.regions[iter.region_idx].first,
        [](hts_pos_t pos, const cached_region_read_t& record) { return pos < record.prefix_max_end; });
    iter.record_idx = first - records.begin();
}

inline bool query_cached_read_regions(region_read_itr_t& iter, char** region_names, unsigned n_regions) {
    const region_read_cache_t* cache = active_region_read_cache();
    if (!cache || cache->bam_file != iter.bam_file || n_regions == 0) return false;

    for (unsigned i = 0; i < n_regions; i++) {
        int tid;
        hts_pos_t begin, end;
        if (!sam_parse_region(iter.bam_file->header, region_names[i], &tid, &begin, &end, HTS_PARSE_THOUSANDS_SEP)) return false;
        if (tid != cache->tid || !region_read_windows_cover(cache->coverage, begin, end) || begin >= end) return false;
        iter.regions.emplace_back(begin, end);
    }
    std::sort(iter.regions.begin(), iter.regions.end());
    iter.cache = cache;
    seek_read_region(iter);
    return true;
}

inline region_read_itr_t query_read_regions(open_samFile_t* bam_file, char** region_names, unsigned n_regions) {
    region_read_itr_t iter;
    iter.bam_file = bam_file;
    if (!query_cached_read_regions(iter, region_names, n_regions)) {
        iter.hts_iter.reset(sam_itr_regarray(bam_file->idx, bam_file->header, region_names, n_regions));
    }
    return iter;
}

inline region_read_itr_t query_read_region(open_samFile_t* bam_file, std::string region_name) {
    region_read_itr_t iter;
    iter.bam_file = bam_file;
    char* region = &region_name[0];
    if (!query_cached_read_regions(iter, &region, 1)) {
        iter.hts_iter.reset(sam_itr_querys(bam_file->idx, bam_file->header, region));
    }
    return iter;
}

inline int next_region_read(region_read_itr_t& iter, bam1_t* out) {
    if (!iter.cache) return sam_itr_next(iter.bam_file->file, iter.hts_iter.get(), out);

    const auto& records = iter.cache->records;
    while (iter.record_idx < records.size() && iter.region_idx < iter.regions.size()) {
        const auto& record = records[iter.record_idx];
        if (record.read.core.pos >= iter.regions[iter.region_idx].second) {
            iter.region_idx++;
            if (iter.region_idx < iter.regions.size()) seek_read_region(iter);
            continue;
        }

        iter.record_idx++;
        if (record.end > iter.regions[iter.region_idx].first) {
            bam1_t source = record.read;
            source.data = const_cast<uint8_t*>(iter.cache->read_data.data()+record.data_offset);
            if (!bam_copy1(out, &source)) throw std::bad_alloc();
            return 0;
        }
    }
    return -1;
}

#endif
