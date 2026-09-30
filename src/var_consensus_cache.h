#ifndef VAR_CONSENSUS_CACHE_H
#define VAR_CONSENSUS_CACHE_H

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_map>
#include "types.h"
#include "sam_utils.h"

namespace sv_consensus_cache {

inline uint64_t hash_bytes(const void* data, size_t size, uint64_t hash = UINT64_C(14695981039346656037)) {
    const unsigned char* bytes = static_cast<const unsigned char*>(data);
    for (size_t i = 0; i < size; i++) hash = (hash ^ bytes[i]) * UINT64_C(1099511628211);
    return hash;
}

inline void write_metrics(std::ostream& out, const consensus_alignment_metrics_t& m) {
    out << '\t' << m.length << '\t' << m.alt_score << '\t' << m.ref_score << '\t' << m.aux_ref_score;
    out << '\t' << m.inferred_missing_aux << '\t' << m.covered_edit_distance << '\t' << m.local_alt_ref_edit_distance << '\t' << m.expected_alt_gain;
    out << '\t' << m.main_edit_covered << '\t' << m.alt_ref_begin << '\t' << m.alt_ref_end;
    for (int i = 0; i < 2; i++) out << '\t' << m.split_ref_lengths[i] << '\t' << m.split_sizes[i] << '\t' << m.split_scores[i] << '\t' << m.independent_ref_scores[i];
}

inline void read_metrics(std::istream& in, consensus_alignment_metrics_t& m) {
    in >> m.length >> m.alt_score >> m.ref_score >> m.aux_ref_score;
    in >> m.inferred_missing_aux >> m.covered_edit_distance >> m.local_alt_ref_edit_distance >> m.expected_alt_gain;
    in >> m.main_edit_covered >> m.alt_ref_begin >> m.alt_ref_end;
    for (int i = 0; i < 2; i++) in >> m.split_ref_lengths[i] >> m.split_sizes[i] >> m.split_scores[i] >> m.independent_ref_scores[i];
}

enum anchor_ownership_t { NO_ANCHORS = 0, LEFT_ANCHOR = 1, RIGHT_ANCHOR = 2, BOTH_ANCHORS = LEFT_ANCHOR | RIGHT_ANCHOR };

struct build_t {
    uint64_t reads_hash = 0, allele_hash = 0;
    int anchors = NO_ANCHORS;
    double avg = 0, stddev = 0;
    // One digit per read: bit 0 = consistent, bit 1 = exact. '-' means no reads.
    std::string flags = "-";
};

struct anchor_t {
    hts_pos_t start = 0, end = 0;
    int seq_len = 0;
    void capture(const std::shared_ptr<sv_t::anchor_aln_t>& anchor) {
        start = anchor->start; end = anchor->end; seq_len = anchor->seq_len;
    }
    void restore(const std::shared_ptr<sv_t::anchor_aln_t>& anchor) const {
        anchor->start = start; anchor->end = end; anchor->seq_len = seq_len;
    }
};

struct record_t {
    hts_pos_t start = 0, end = 0;
    uint64_t variant_hash = 0;
    std::vector<build_t> builds;
    consensus_alignment_metrics_t alt1, alt2, extended1, extended2;
    anchor_t left, right;

    void write(std::ostream& out) const {
        out << '\t' << start << '\t' << end << '\t' << variant_hash << '\t' << builds.size();
        for (const auto& b : builds) out << '\t' << b.reads_hash << '\t' << b.allele_hash << '\t' << b.avg << '\t' << b.stddev << '\t' << b.flags << '\t' << b.anchors;
        write_metrics(out, alt1); write_metrics(out, alt2); write_metrics(out, extended1); write_metrics(out, extended2);
        for (const auto& a : {left, right}) out << '\t' << a.start << '\t' << a.end << '\t' << a.seq_len;
        out << '\n';
    }
    void read(std::istream& in) {
        size_t n;
        in >> start >> end >> variant_hash >> n;
        builds.resize(n);
        for (auto& b : builds) in >> b.reads_hash >> b.allele_hash >> b.avg >> b.stddev >> b.flags >> b.anchors;
        read_metrics(in, alt1); read_metrics(in, alt2); read_metrics(in, extended1); read_metrics(in, extended2);
        for (auto* a : {&left, &right}) in >> a->start >> a->end >> a->seq_len;
    }
    void capture(sv_t* sv) {
        alt1 = sv->sample_info.alt_consensus1_metrics; alt2 = sv->sample_info.alt_consensus2_metrics;
        extended1 = sv->sample_info.ext_alt_consensus1_metrics; extended2 = sv->sample_info.ext_alt_consensus2_metrics;
        left.capture(sv->left_anchor_aln); right.capture(sv->right_anchor_aln);
    }
    void restore_alt(sv_t* sv, size_t index) const {
        if (index == 0) {
            sv->sample_info.alt_consensus1_metrics = alt1;
            sv->sample_info.ext_alt_consensus1_metrics = extended1;
        } else {
            sv->sample_info.alt_consensus2_metrics = alt2;
            sv->sample_info.ext_alt_consensus2_metrics = extended2;
        }
        if (builds[index].anchors & LEFT_ANCHOR) left.restore(sv->left_anchor_aln);
        if (builds[index].anchors & RIGHT_ANCHOR) right.restore(sv->right_anchor_aln);
    }
};

class cache_t {
    struct contig_cache_t {
        int id;
        std::unordered_map<std::string, record_t> by_id;
    };
    std::string directory;
    bool reading = false;
    std::unordered_map<std::string, contig_cache_t> records;
public:
    bool active() const { return !directory.empty(); }
    void initialize(const std::string& path, bool read) { directory = path; reading = read; records.clear(); }

    // Called before any workers start. Map keys stay fixed while workers run;
    // each SV is processed by one worker, so separate records need no locking.
    void prepare(const std::string& contig, int contig_id, const std::vector<sv_t*>& svs) {
        if (!active()) return;
        auto& cache = records[contig];
        cache.id = contig_id;
        auto& by_id = cache.by_id;
        if (reading) {
            std::ifstream in(directory+"/"+std::to_string(contig_id)+".tsv");
            std::string line;
            while (std::getline(in, line)) {
                std::istringstream row(line);
                std::string id; row >> id;
                by_id[id].read(row);
            }
        } else {
            for (auto* sv : svs) by_id.emplace(sv->id, record_t());
        }
    }
    const record_t* find(sv_t* sv) const {
        const auto& by_id = records.at(sv->chr).by_id;
        auto it = by_id.find(sv->id);
        return it == by_id.end() ? nullptr : &it->second;
    }
    void save(sv_t* sv, record_t& record) {
        if (!reading) records.at(sv->chr).by_id.at(sv->id) = std::move(record);
    }
    // Only the main thread writes, after all workers have finished.
    void write() const {
        if (!active() || reading) return;
        for (const auto& contig : records) {
            std::ofstream out;
            out << std::setprecision(std::numeric_limits<double>::max_digits10);
            for (const auto& entry : contig.second.by_id) {
                if (entry.second.builds.empty()) continue;
                if (!out.is_open()) out.open(directory+"/"+std::to_string(contig.second.id)+".tsv");
                out << entry.first;
                entry.second.write(out);
            }
        }
    }
};

inline cache_t& cache() { static cache_t instance; return instance; }

struct build_input_t {
    const char* allele;
    const std::vector<std::shared_ptr<bam1_t>>& reads;
    int anchors; // Alt builds declare their anchors; reference builds own none.
};

class session_t {
    record_t pending;
    const record_t* cached = nullptr;
    bool enabled;
public:
    session_t(sv_t* sv, const std::vector<build_input_t>& inputs) : enabled(cache().active()) {
        if (!enabled) return;
        pending.start = sv->start; pending.end = sv->end;
        std::string variant = sv->unique_key();
        pending.variant_hash = hash_bytes(variant.data(), variant.size());
        for (const auto& input : inputs) {
            build_t build;
            build.anchors = input.anchors;
            build.reads_hash = hash_bytes(nullptr, 0);
            for (const auto& read : input.reads) {
                std::string name = std::string(bam_get_qname(read.get())) + (is_first_read(read.get()) ? "/1" : "/2");
                // Include the terminator to separate names unambiguously.
                build.reads_hash = hash_bytes(name.c_str(), name.size()+1, build.reads_hash);
            }
            build.allele_hash = hash_bytes(input.allele, strlen(input.allele));
            pending.builds.push_back(build);
        }
        cached = cache().find(sv);
        if (cached && (cached->start != pending.start || cached->end != pending.end || cached->variant_hash != pending.variant_hash || cached->builds.size() != pending.builds.size())) cached = nullptr;
        for (size_t i = 0; i < inputs.size(); i++) {
            if (inputs[i].anchors != NO_ANCHORS && hit(i)) cached->restore_alt(sv, i);
        }
    }
    bool hit(size_t index) const {
        if (!cached) return false;
        const auto& a = pending.builds[index]; const auto& b = cached->builds[index];
        return a.reads_hash == b.reads_hash && a.allele_hash == b.allele_hash && a.anchors == b.anchors;
    }
    template<class Compute> std::vector<bool> build(size_t index, double& avg, double& stddev, std::vector<bool>& exact, Compute compute) {
        if (hit(index)) {
            const build_t& b = cached->builds[index];
            avg = b.avg; stddev = b.stddev;
            size_t n = b.flags == "-" ? 0 : b.flags.size();
            std::vector<bool> consistent(n); exact.resize(n);
            for (size_t i = 0; i < n; i++) {
                consistent[i] = (b.flags[i]-'0') & 1; exact[i] = (b.flags[i]-'0') & 2;
            }
            return consistent;
        }
        auto consistent = compute();
        if (enabled) {
            build_t& b = pending.builds[index]; b.avg = avg; b.stddev = stddev;
            if (!consistent.empty()) b.flags.resize(consistent.size());
            for (size_t i = 0; i < consistent.size(); i++) b.flags[i] = '0' + unsigned(consistent[i]) + 2*unsigned(exact[i]);
        }
        return consistent;
    }
    void finish(sv_t* sv) { if (enabled && !cached) { pending.capture(sv); cache().save(sv, pending); } }
};

} // namespace sv_consensus_cache
#endif
