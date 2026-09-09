#ifndef HP_MISMATCH_RATE_THRESHOLDS_H
#define HP_MISMATCH_RATE_THRESHOLDS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "hp_read_info.h"
#include "consensus.h"

static const char HP_MISMATCH_RATE_THRESHOLDS_FILENAME[] = "hp_3p_mismatch_rate_thresholds.tsv";

struct hp_mismatch_rate_thresholds_t {
    std::array<std::vector<double>, 4> thresholds_by_base;

    explicit hp_mismatch_rate_thresholds_t(const std::string& fname) { load(fname); }

    void load(const std::string& fname) {
        std::ifstream fin(fname);
        if (!fin) throw std::runtime_error("Unable to open " + fname + ".");

        std::string hp_len_header, a_header, c_header, g_header, t_header;
        if (!(fin >> hp_len_header >> a_header >> c_header >> g_header >> t_header) || hp_len_header != "HP_LEN" || a_header != "A" || c_header != "C" || g_header != "G" || t_header != "T") {
            throw std::runtime_error("Invalid header in " + fname + ".");
        }

        for (std::vector<double>& thresholds : thresholds_by_base) thresholds.clear();
        int hp_len;
        double thresholds[4];
        int expected_hp_len = 0;
        while (fin >> hp_len) {
            if (!(fin >> thresholds[0] >> thresholds[1] >> thresholds[2] >> thresholds[3])) throw std::runtime_error("Invalid data in " + fname + ".");
            if (hp_len != expected_hp_len) throw std::runtime_error("Invalid HP length in " + fname + ".");
            for (int hp_base_idx = 0; hp_base_idx < 4; hp_base_idx++) {
                if (!std::isfinite(thresholds[hp_base_idx]) || thresholds[hp_base_idx] < 0 || thresholds[hp_base_idx] > 1) throw std::runtime_error("Invalid threshold in " + fname + ".");
                thresholds_by_base[hp_base_idx].push_back(thresholds[hp_base_idx]);
            }
            expected_hp_len++;
        }
        if (!fin.eof()) throw std::runtime_error("Invalid data in " + fname + ".");
        if (expected_hp_len == 0) throw std::runtime_error("No HP mismatch-rate thresholds in " + fname + ".");
    }

    double get_threshold(int hp_len, char hp_base) const {
        if (hp_len < 0) throw std::invalid_argument("HP length cannot be negative.");
        int hp_base_idx = get_base_idx(hp_base);
        if (hp_base_idx < 0) throw std::invalid_argument("HP base must be A, C, G, or T.");
        const std::vector<double>& thresholds = thresholds_by_base[hp_base_idx];
        return thresholds[std::min<size_t>(hp_len, thresholds.size() - 1)];
    }

private:
    static int get_base_idx(char base) {
        if (base == 'A' || base == 'a') return 0;
        if (base == 'C' || base == 'c') return 1;
        if (base == 'G' || base == 'g') return 2;
        if (base == 'T' || base == 't') return 3;
        return -1;
    }
};

struct hp_read_observation_t {
    std::string seq;
    std::vector<uint8_t> quals;
    int left_tail_len;
    bool is_reverse;
    double ref_3p_mismatch_rate;
    int right_tail_len; // Includes any HP bases reassigned to this tail by gap correction.
    std::vector<uint8_t> original_quals; // Preserves input qualities when normalization or tail recalibration changes the working qualities.
};

// Key: sequenced HP base, corrected HP length, 1-based distance from HP into the 3' tail, original Phred quality.
using hp_tail_error_table_t = std::map<std::array<int, 4>, std::pair<uint64_t, uint64_t>>; // Compared bases, mismatches.
using hp_tail_quality_table_t = std::map<std::array<int, 4>, int>; // Cached Phred qualities; -1 means no calibration data.
using hp_tail_error_bin_t = std::pair<std::array<int, 4>, std::pair<uint64_t, uint64_t>>;

struct hp_tail_quality_model_t {
    std::map<std::array<int, 2>, std::vector<hp_tail_error_bin_t>> bins_by_base_and_pos;
    uint64_t min_observations = 100;
};

inline hp_tail_quality_model_t make_hp_tail_quality_model(const hp_tail_error_table_t& error_table) {
    hp_tail_quality_model_t model;
    for (const auto& entry : error_table) {
        if (entry.second.first > 0) model.bins_by_base_and_pos[{entry.first[0], entry.first[2]}].push_back(entry);
    }
    return model;
}

inline hp_tail_quality_model_t read_hp_tail_quality_model(const std::string& fname) {
    std::ifstream fin(fname);
    if (!fin) throw std::runtime_error("Unable to open HP tail calibration table " + fname + ".");
    hp_tail_error_table_t table;
    std::string line;
    bool header_read = false;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (!header_read) {
            if (line != "HP_BASE\tHP_LEN\tTAIL_POS\tBASE_QUAL\tOBSERVATIONS\tERRORS\tERROR_PROBABILITY") throw std::runtime_error("Invalid HP tail calibration header in " + fname + ".");
            header_read = true;
            continue;
        }
        std::istringstream row(line);
        std::string base, extra;
        int hp_len, tail_pos, qual;
        int64_t observations, errors;
        double probability;
        if (!(row >> base >> hp_len >> tail_pos >> qual >> observations >> errors >> probability) || row >> extra || base.size() != 1 || base_to_index(base[0]) < 0 || hp_len < 0 || tail_pos < 1 || qual < 0 || qual >= 255 || observations <= 0 || errors < 0 || errors > observations || !std::isfinite(probability) || probability < 0 || probability > 1) throw std::runtime_error("Invalid HP tail calibration row in " + fname + ".");
        // Counts retain full precision and supply the weights when neighboring bins are pooled.
        if (!table.emplace(std::array<int, 4>{base_to_index(base[0]), hp_len, tail_pos, qual}, std::make_pair(uint64_t(observations), uint64_t(errors))).second) throw std::runtime_error("Duplicate HP tail calibration bin in " + fname + ".");
    }
    if (!header_read || fin.bad()) throw std::runtime_error("Unable to read HP tail calibration table " + fname + ".");
    return make_hp_tail_quality_model(table);
}

inline int estimate_hp_tail_quality(const std::array<int, 4>& key, const hp_tail_quality_model_t& model) {
    auto group = model.bins_by_base_and_pos.find({key[0], key[2]});
    if (group == model.bins_by_base_and_pos.end()) return -1;
    std::vector<std::pair<int, const hp_tail_error_bin_t*>> neighbors;
    for (const auto& bin : group->second) {
        if (bin.first[1] < (key[1] + 1) / 2 || bin.first[1] > 2 * key[1]) continue;
        int distance = std::abs(bin.first[1] - key[1]) + std::abs(bin.first[3] - key[3]);
        neighbors.push_back({distance, &bin});
    }
    std::sort(neighbors.begin(), neighbors.end(), [](const std::pair<int, const hp_tail_error_bin_t*>& a, const std::pair<int, const hp_tail_error_bin_t*>& b) { return a.first < b.first; });
    uint64_t observations = 0, errors = 0;
    // Pool equally near HP-length/quality bins together, stopping once the sample target is reached.
    for (size_t i = 0; i < neighbors.size() && observations < model.min_observations;) {
        int distance = neighbors[i].first;
        do {
            observations += neighbors[i].second->second.first;
            errors += neighbors[i].second->second.second;
            i++;
        } while (i < neighbors.size() && neighbors[i].first == distance);
    }
    if (observations < model.min_observations) return -1;
    // Phred Q = -10 log10(P(error)); zero observed errors receive the existing Q40 cap.
    double qual = errors == 0 ? 40.0 : -10.0 * std::log10(double(errors) / observations);
    return std::max(0, std::min(40, (int) std::lround(qual)));
}

inline void recalibrate_hp_3p_tail_qualities(std::vector<hp_read_observation_t>& observations, int hp_len, char hp_base, const hp_tail_quality_model_t& model, hp_tail_quality_table_t& quality_cache) {
    int hp_base_idx = base_to_index(hp_base);
    if (hp_base_idx < 0 || hp_len < 0 || model.bins_by_base_and_pos.empty()) return;
    for (hp_read_observation_t& observation : observations) {
        int sequenced_hp_base_idx = observation.is_reverse ? 3 - hp_base_idx : hp_base_idx;
        int tail_beg = observation.is_reverse ? 0 : (int) observation.seq.length() - observation.right_tail_len;
        int tail_end = observation.is_reverse ? observation.left_tail_len : observation.seq.length();
        for (int qpos = tail_beg; qpos < tail_end; qpos++) {
            int original_qual = observation.original_quals.empty() ? observation.quals[qpos] : observation.original_quals[qpos];
            if (original_qual == 255 || base_to_index(observation.seq[qpos]) < 0) continue;
            int tail_pos = observation.is_reverse ? observation.left_tail_len - qpos : qpos - tail_beg + 1;
            std::array<int, 4> key = {sequenced_hp_base_idx, hp_len, tail_pos, original_qual};
            auto it = quality_cache.find(key);
            if (it == quality_cache.end()) it = quality_cache.emplace(key, estimate_hp_tail_quality(key, model)).first;
            if (it->second < 0) continue;
            if (observation.original_quals.empty()) observation.original_quals = observation.quals;
            observation.quals[qpos] = it->second;
        }
    }
}

const int MIN_REF_HP_LEN = 5;

struct hp_run_context_t {
    hts_pos_t beg, end;
    char base;
    int left_side_5p_reads = 0, right_side_5p_reads = 0;
    int left_side_5p_indel_reads = 0, right_side_5p_indel_reads = 0;

    hp_run_context_t(hts_pos_t beg, hts_pos_t end, char base) : beg(beg), end(end), base(base) {}
};

// Correct direct HP expansions using a clipped 3' tail of at least min_tail_len bases.
// Reassign the subtracted HP bases to the 3' tail; keep the 5' boundary and read sequence unchanged.
inline void subtract_hp_3p_tail_gap(hp_read_info_t& hp_read_info, const hp_run_context_t& hp_run, char* contig_seq, hts_pos_t contig_len,
    bool is_reverse, bool clipped_3p, int min_tail_len) {

    int tail_len = hp_read_info.tail_3p_len;
    if (hp_read_info.hp_len_iteratively_resolved || hp_read_info.hp_len <= hp_run.end - hp_run.beg || !clipped_3p || tail_len < min_tail_len || tail_len <= 0) return;
    // Block if >=10% of qualifying 5' reads on the 3' side carry an indel: left for reverse, right for forward.
    // No qualifying reads on that side means no supported blocker.
    int reads = is_reverse ? hp_run.left_side_5p_reads : hp_run.right_side_5p_reads;
    int indel_reads = is_reverse ? hp_run.left_side_5p_indel_reads : hp_run.right_side_5p_indel_reads;
    if (reads > 0 && 10 * indel_reads >= reads) return;

    const std::string& read_seq = hp_read_info.read.seq;
    int tail_beg = is_reverse ? 0 : read_seq.length() - tail_len;
    int best_gap = 0, best_mismatches = tail_len + 1;
    // Remap the entire 3' tail ungapped, leaving 0-10 reference bases between it and the HP.
    // Minimize mismatches; ascending gaps and strict improvement favor the smallest gap on ties.
    for (int gap = 0; gap <= 10; gap++) {
        hts_pos_t ref_beg = is_reverse ? hp_run.beg - gap - tail_len : hp_run.end + gap;
        if (ref_beg < 0 || ref_beg + tail_len > contig_len) continue;
        int mismatches = number_of_mismatches_fast(read_seq.c_str() + tail_beg, contig_seq + ref_beg, tail_len, best_mismatches - 1);
        if (mismatches < best_mismatches) {
            best_mismatches = mismatches;
            best_gap = gap;
            if (best_mismatches == 0) break;
        }
    }
    if (best_gap == 0) return;
    hp_read_info.hp_len -= best_gap;
    // Negative estimates are undefined; zero remains a valid corrected length.
    if (hp_read_info.hp_len < 0) {
        hp_read_info.hp_len = UNDEFINED_HP_LEN;
        return;
    }
    hp_read_info.tail_3p_len += best_gap;
    // Fill the selected reference gap with the newly included query bases; do not search again.
    int enlarged_tail_beg = is_reverse ? 0 : tail_beg - best_gap;
    hts_pos_t enlarged_ref_beg = is_reverse ? hp_run.beg - hp_read_info.tail_3p_len : hp_run.end;
    hp_read_info.tail_3p_mismatches = number_of_mismatches_fast(read_seq.c_str() + enlarged_tail_beg, contig_seq + enlarged_ref_beg, hp_read_info.tail_3p_len, hp_read_info.tail_3p_len);
}

inline bool is_usable_hp_read(const hp_read_info_t& hp_read_info, int min_clip_len, double max_seq_error) {

    if (hp_read_info.hp_len == UNDEFINED_HP_LEN) return false;
    if (hp_read_info.hp_deletion_extends_outside_hp || hp_read_info.hp_insertion_has_non_hp_bases) return false;
    if (hp_read_info.tail_5p_len < min_clip_len || hp_read_info.tail_3p_len < min_clip_len) return false;
    if (double(hp_read_info.tail_5p_mismatches) / hp_read_info.tail_5p_len > max_seq_error) return false;
    return true;
}

template<typename HpRun = hp_run_context_t>
std::vector<HpRun> find_hp_runs(char* contig_seq, hts_pos_t contig_len, hts_pos_t chunk_beg, hts_pos_t chunk_end) {

    std::vector<HpRun> hp_runs;
    hts_pos_t beg = chunk_beg;
    if (beg > 0 && beg < contig_len && contig_seq[beg] == contig_seq[beg - 1]) {
        char base = contig_seq[beg];
        while (beg < contig_len && contig_seq[beg] == base) beg++;
    }

    for (; beg < chunk_end;) {
        char base = contig_seq[beg];
        hts_pos_t end = beg + 1;
        while (end < contig_len && contig_seq[end] == base) end++;

        if ((base == 'A' || base == 'C' || base == 'G' || base == 'T') && end - beg >= MIN_REF_HP_LEN) {
            hp_runs.emplace_back(beg, end, base);
        }
        beg = end;
    }
    return hp_runs;
}

inline void add_hp_5p_blocker_evidence(bam1_t* read, hp_run_context_t& hp_run, const config_t& config, const stats_t& stats) {
    if (is_unmapped(read) || !is_primary(read) || read->core.l_qseq <= 0 || !is_proper_pair(read, stats.min_is, stats.max_is)) return;
    hp_adjacent_indel_info_t indel_info = get_adjacent_indel_info(read, {hp_run.beg, hp_run.end}, stats.read_len/2);
    bool is_reverse = bam_is_rev(read);
    const hp_side_indel_info_t& five_p_info = is_reverse ? indel_info.right : indel_info.left;
    int& side_5p_reads = is_reverse ? hp_run.right_side_5p_reads : hp_run.left_side_5p_reads;
    int& side_5p_indel_reads = is_reverse ? hp_run.right_side_5p_indel_reads : hp_run.left_side_5p_indel_reads;
    if (five_p_info.aligned_len >= config.min_clip_len) {
        side_5p_reads++;
        if (five_p_info.indel_len != 0) side_5p_indel_reads++;
    }
}

inline hp_read_info_t estimate_hp_read_for_calibration(bam1_t* read, const hp_run_context_t& hp_run, char* contig_seq, hts_pos_t contig_len, const config_t& config) {
    hts_pair_pos_t hp_range = {hp_run.beg, hp_run.end};
    hts_pos_t extend = read->core.l_qseq - 1;
    hts_pos_t ref_allele_beg = std::max((hts_pos_t) 0, hp_run.beg - extend);
    hts_pos_t ref_allele_end = std::min(contig_len, hp_run.end + extend);
    hts_pair_pos_t ref_allele_hp_range = {hp_run.beg - ref_allele_beg, hp_run.end - ref_allele_beg};
    bool has_no_left_indel = five_p_evidence_permits_iterative_hp_len_estimation(hp_run.left_side_5p_reads, hp_run.left_side_5p_indel_reads);
    bool has_no_right_indel = five_p_evidence_permits_iterative_hp_len_estimation(hp_run.right_side_5p_reads, hp_run.right_side_5p_indel_reads);
    hp_read_info_t info = calculate_hp_read_info(read, hp_range, hp_run.base, contig_seq, contig_len, contig_seq + ref_allele_beg, ref_allele_end - ref_allele_beg, ref_allele_hp_range, has_no_left_indel, has_no_right_indel, 0, config.max_seq_error);
    // Keep the original BAM clipping flag even when HP estimation uses a secondary alignment.
    bool is_reverse = bam_is_rev(read);
    subtract_hp_3p_tail_gap(info, hp_run, contig_seq, contig_len, is_reverse, is_reverse ? get_left_clip_size(read) > 0 : get_right_clip_size(read) > 0, config.min_clip_len);
    return info;
}

inline void collect_hp_5p_blocker_evidence(open_samFile_t* alignment_file, const std::string& contig_name, std::vector<hp_run_context_t>& hp_runs, const config_t& config, const stats_t& stats) {
    if (hp_runs.empty()) return;
    std::string region = contig_name + ":" + std::to_string(hp_runs.front().beg + 1) + "-" + std::to_string(hp_runs.back().end);
    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iter(sam_itr_querys(alignment_file->idx, alignment_file->header, region.c_str()), &hts_itr_destroy);
    if (!iter) throw std::runtime_error("Unable to query alignments for " + region + ".");
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), &bam_destroy1);
    int status;
    while ((status = sam_itr_next(alignment_file->file, iter.get(), read.get())) >= 0) {
        auto first = std::upper_bound(hp_runs.begin(), hp_runs.end(), read->core.pos, [](hts_pos_t pos, const hp_run_context_t& run) { return pos < run.end; });
        for (auto run = first; run != hp_runs.end() && run->beg < bam_endpos(read.get()); run++) add_hp_5p_blocker_evidence(read.get(), *run, config, stats);
    }
    if (status < -1) throw std::runtime_error("Error while reading alignments for " + region + ".");
}

// HP interpretation selects a calibration row only; it never changes the clip read's sequence or placement.
inline std::vector<uint8_t> recalibrate_clip_read_qualities(bam1_t* read, const std::vector<hp_run_context_t>& hp_runs, char* contig_seq, hts_pos_t contig_len, const config_t& config, const hp_tail_quality_model_t& model, hp_tail_quality_table_t& quality_cache) {
    const uint8_t* bam_quals = bam_get_qual(read);
    std::vector<uint8_t> quals(bam_quals, bam_quals + read->core.l_qseq);
    std::replace(quals.begin(), quals.end(), uint8_t(255), uint8_t(0));
    if (model.bins_by_base_and_pos.empty() || is_unmapped(read) || !is_primary(read) || read->core.l_qseq <= 0) return quals;
    hp_read_info_t selected;
    const hp_run_context_t* selected_run = nullptr;
    auto first = std::upper_bound(hp_runs.begin(), hp_runs.end(), read->core.pos, [](hts_pos_t pos, const hp_run_context_t& run) { return pos < run.end; });
    for (auto run = first; run != hp_runs.end() && run->beg < bam_endpos(read); run++) {
        hp_read_info_t info = estimate_hp_read_for_calibration(read, *run, contig_seq, contig_len, config);
        if (!is_usable_hp_read(info, config.min_clip_len, config.max_seq_error)) continue;
        if (selected_run != nullptr) return quals; // Conflicting HP contexts do not define a unique calibration row.
        selected = std::move(info);
        selected_run = &*run;
    }
    if (selected_run == nullptr) return quals;
    bool is_reverse = bam_is_rev(read);
    int left_tail_len = is_reverse ? selected.tail_3p_len : selected.tail_5p_len;
    int right_tail_len = is_reverse ? selected.tail_5p_len : selected.tail_3p_len;
    std::vector<hp_read_observation_t> observations{{selected.read.seq, quals, left_tail_len, is_reverse, 0.0, right_tail_len, std::vector<uint8_t>(bam_quals, bam_quals + read->core.l_qseq)}};
    recalibrate_hp_3p_tail_qualities(observations, selected.hp_len, selected_run->base, model, quality_cache);
    return std::move(observations[0].quals);
}

#endif // HP_MISMATCH_RATE_THRESHOLDS_H
