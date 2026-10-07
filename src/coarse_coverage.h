#ifndef COARSE_COVERAGE_H_
#define COARSE_COVERAGE_H_

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class coarse_coverage_builder_t {
    int64_t contig_len_, window_size_, position_ = 0;
    uint64_t depth_ = 0;
    std::priority_queue<int64_t, std::vector<int64_t>, std::greater<int64_t>> ends_;

    void record_segment(int64_t begin, int64_t end) {
        begin = std::max(int64_t(0), begin);
        end = std::min(contig_len_, end);
        if (depth_ == 0 || begin >= end) return;
        uint32_t depth = depth_ > std::numeric_limits<uint32_t>::max() ? std::numeric_limits<uint32_t>::max() : depth_;
        for (size_t window = begin/window_size_; window <= size_t((end-1)/window_size_); window++) maxima[window] = std::max(maxima[window], depth);
    }

    void advance_to(int64_t position) {
        if (position < position_) throw std::runtime_error("Coverage alignments are not coordinate sorted");
        while (!ends_.empty() && ends_.top() <= position) {
            int64_t end = ends_.top();
            record_segment(position_, end);
            uint64_t ended = 0;
            while (!ends_.empty() && ends_.top() == end) {
                ends_.pop();
                ended++;
            }
            depth_ -= ended;
            position_ = end;
        }
        record_segment(position_, position);
        position_ = position;
    }

public:
    std::vector<uint32_t> maxima;

    coarse_coverage_builder_t(int64_t contig_len, int64_t window_size) : contig_len_(contig_len), window_size_(window_size) {
        if (contig_len < 0 || window_size <= 0) throw std::invalid_argument("Invalid coarse coverage dimensions");
        maxima.resize((contig_len+window_size-1)/window_size);
    }

    void add_alignment(int64_t begin, int64_t end) {
        begin = std::max(int64_t(0), begin);
        end = std::min(contig_len_, end);
        if (begin >= end) return;
        advance_to(begin);
        ends_.push(end);
        depth_++;
    }

    void finish() {
        advance_to(contig_len_);
    }
};

class coarse_coverage_track_t {
    int64_t window_size_ = 1;

public:
    std::vector<uint32_t> maxima;

    coarse_coverage_track_t() {}
    coarse_coverage_track_t(int64_t window_size, std::vector<uint32_t> maxima) : window_size_(window_size), maxima(std::move(maxima)) {
        if (window_size <= 0) throw std::invalid_argument("Invalid coarse coverage window size");
    }

    bool exceeds(int64_t begin, int64_t end, uint64_t threshold) const {
        begin = std::max(int64_t(0), begin);
        if (begin >= end || maxima.empty()) return false;
        size_t first = begin/window_size_, last = (end-1)/window_size_;
        if (first >= maxima.size()) return false;
        last = std::min(last, maxima.size()-1);
        for (size_t window = first; window <= last; window++) if (maxima[window] > threshold) return true;
        return false;
    }
};

inline void write_coarse_coverage(const std::string& path, const std::vector<uint32_t>& maxima) {
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw std::runtime_error("Failed to open coarse coverage track " + path);
    if (!maxima.empty()) output.write(reinterpret_cast<const char*>(maxima.data()), maxima.size()*sizeof(uint32_t));
    if (!output) throw std::runtime_error("Failed to write coarse coverage track " + path);
}

inline coarse_coverage_track_t read_coarse_coverage(const std::string& path, int64_t contig_len, int64_t window_size) {
    if (contig_len < 0 || window_size <= 0) throw std::invalid_argument("Invalid coarse coverage dimensions");
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) return {};
    size_t expected_windows = (contig_len+window_size-1)/window_size;
    std::streamsize expected_bytes = expected_windows*sizeof(uint32_t);
    if (input.tellg() != expected_bytes) throw std::runtime_error("Invalid coarse coverage track size in " + path);
    input.seekg(0);
    std::vector<uint32_t> maxima(expected_windows);
    if (!maxima.empty()) input.read(reinterpret_cast<char*>(maxima.data()), expected_bytes);
    if (!input) throw std::runtime_error("Failed to read coarse coverage track " + path);
    return coarse_coverage_track_t(window_size, std::move(maxima));
}

#endif
