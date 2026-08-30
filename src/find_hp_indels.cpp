#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "htslib/sam.h"
#include "htslib/tbx.h"
#include "htslib/vcf.h"

#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "utils.h"
#include "sam_utils.h"
#include "hp_read_info.h"
#include "hp_mismatch_rate_thresholds.h"
#include "vcf_utils.h"
#include "consensus.h"

const int MIN_REF_HP_LEN = 5;
const int MIN_SUPPORTING_READS = 3;
const int MIN_READS_FOR_HP_MM_THRESHOLD = 100;
const double HP_MM_OUTLIER_MAD_MULTIPLIER = 4.0;
const double MAX_HP_MM_THRESHOLD = 0.75;
const hts_pos_t CHUNK_SIZE = 1000000;

struct mate_info_t {
    std::string seq, qual;
    int mapq;
};

using mate_map_t = std::unordered_map<std::string, mate_info_t>;
struct hp_read_observation_t {
    std::string seq;
    std::vector<uint8_t> quals;
    int left_tail_len;
    bool is_reverse;
    double ref_3p_mismatch_rate;
};

struct hp_run_t {
    hts_pos_t beg, end;
    char base;
    int usable_reads = 0;
    std::unordered_map<int, int> hp_len_counts;
    std::unordered_map<int, std::vector<hp_read_observation_t>> observations_by_hp_len;
    int left_side_5p_reads = 0, right_side_5p_reads = 0;
    int left_side_5p_indel_reads = 0, right_side_5p_indel_reads = 0;

    hp_run_t(hts_pos_t beg, hts_pos_t end, char base) : beg(beg), end(end), base(base) {}
};

bool is_usable_hp_read(const hp_read_info_t& hp_read_info, int min_clip_len, double max_seq_error) {

    if (hp_read_info.hp_len == UNDEFINED_HP_LEN) return false;
    if (hp_read_info.hp_deletion_extends_outside_hp || hp_read_info.hp_insertion_has_non_hp_bases) return false;
    if (hp_read_info.tail_5p_len < min_clip_len || hp_read_info.tail_3p_len < min_clip_len) return false;
    if (double(hp_read_info.tail_5p_mismatches) / hp_read_info.tail_5p_len > max_seq_error) return false;
    return true;
}

std::vector<hp_run_t> find_hp_runs(char* contig_seq, hts_pos_t contig_len, hts_pos_t chunk_beg, hts_pos_t chunk_end) {

    std::vector<hp_run_t> hp_runs;
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

std::string workdir;
std::vector<mate_map_t> mateseqs_w_mapq;
std::vector<int> active_threads_per_chr;
std::vector<std::mutex> mutex_per_chr;

void read_mates(int contig_id) {
    std::lock_guard<std::mutex> lock(mutex_per_chr[contig_id]);
    if (active_threads_per_chr[contig_id] == 0) {
        std::string fname = workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt";
        std::ifstream fin(fname);
        std::string qname, read_seq, qual;
        int mapq;
        while (fin >> qname >> read_seq >> qual >> mapq) {
            mateseqs_w_mapq[contig_id][qname] = {read_seq, qual, mapq};
        }
    }
    active_threads_per_chr[contig_id]++;
}

void release_mates(int contig_id) {
    std::lock_guard<std::mutex> lock(mutex_per_chr[contig_id]);
    active_threads_per_chr[contig_id]--;
    if (active_threads_per_chr[contig_id] == 0) {
        mateseqs_w_mapq[contig_id].clear();
    }
}

struct hp_positional_consensus_t {
    std::string seq;
    std::vector<int> coverage;
    std::vector<int> read_offsets;
    int trim_beg = 0;
    int hp_beg = 0, hp_end = 0;
    int callable_beg = 0, callable_end = 0;
};

using hp_mismatch_rates_by_len_t = std::vector<std::vector<double>>;
using hp_mismatch_rates_by_len_and_base_t = std::array<hp_mismatch_rates_by_len_t, 4>;

struct hp_read_mismatch_rate_t {
    double rate;
    int sequenced_hp_base_idx;
};

struct hp_chunk_result_t {
    std::vector<std::shared_ptr<sv_t>> hp_indels;
    std::vector<std::vector<hp_read_mismatch_rate_t>> hp_indel_3p_mismatch_rates;
    hp_mismatch_rates_by_len_and_base_t threshold_estimation_rates_by_hp_len_and_base;
};

hp_positional_consensus_t build_hp_positional_consensus(const std::vector<hp_read_observation_t>& observations, int hp_len) {

    hp_positional_consensus_t consensus;
    if (observations.empty()) return consensus;

    int hp_anchor = 0;
    for (const auto& observation : observations) {
        hp_anchor = std::max(hp_anchor, observation.left_tail_len);
    }
    consensus.hp_beg = hp_anchor;
    consensus.hp_end = hp_anchor + hp_len;
    consensus.callable_beg = consensus.hp_beg;
    consensus.callable_end = consensus.hp_end;

    int consensus_len = 0;
    std::vector<std::string> seqs;
    std::vector<const uint8_t*> quals;
    std::vector<hts_pos_t> read_offsets;
    for (const auto& observation : observations) {
        int offset = hp_anchor - observation.left_tail_len;
        consensus.read_offsets.push_back(offset);
        read_offsets.push_back(offset);
        seqs.push_back(observation.seq);
        quals.push_back(observation.quals.data());
        consensus_len = std::max(consensus_len, offset + (int) observation.seq.length());
        if (observation.is_reverse) {
            consensus.callable_end = std::max(consensus.callable_end, offset + (int) observation.seq.length());
        } else {
            consensus.callable_beg = std::min(consensus.callable_beg, offset);
        }
    }

    int left_5p_reads = 0, right_5p_reads = 0;
    for (const hp_read_observation_t& observation : observations) {
        if (observation.is_reverse) right_5p_reads++;
        else left_5p_reads++;
    }
    for (int read_idx = 0; read_idx < observations.size(); read_idx++) {
        int offset = consensus.read_offsets[read_idx];
        if (left_5p_reads >= MIN_SUPPORTING_READS && observations[read_idx].is_reverse) {
            int left_flank_end = std::min((int) seqs[read_idx].length(), consensus.hp_beg - offset);
            std::fill(seqs[read_idx].begin(), seqs[read_idx].begin() + std::max(0, left_flank_end), 'N');
        }
        if (right_5p_reads >= MIN_SUPPORTING_READS && !observations[read_idx].is_reverse) {
            int right_flank_beg = std::max(0, consensus.hp_end - offset);
            std::fill(seqs[read_idx].begin() + std::min((int) seqs[read_idx].length(), right_flank_beg), seqs[read_idx].end(), 'N');
        }
    }

    positional_consensus_t positional_consensus = build_positional_consensus(seqs, quals, read_offsets);
    consensus.seq = positional_consensus.seq;
    consensus.coverage = positional_consensus.coverage;

    int trim_end = consensus_len;
    while (consensus.trim_beg < trim_end && consensus.coverage[consensus.trim_beg] < MIN_SUPPORTING_READS) {
        consensus.trim_beg++;
    }
    while (trim_end > consensus.trim_beg && consensus.coverage[trim_end - 1] < MIN_SUPPORTING_READS) {
        trim_end--;
    }
    consensus.seq = consensus.seq.substr(consensus.trim_beg, trim_end - consensus.trim_beg);
    consensus.coverage = std::vector<int>(consensus.coverage.begin() + consensus.trim_beg, consensus.coverage.begin() + trim_end);
    consensus.hp_beg -= consensus.trim_beg;
    consensus.hp_end -= consensus.trim_beg;
    consensus.callable_beg = std::max(0, consensus.callable_beg - consensus.trim_beg);
    consensus.callable_end = std::min((int) consensus.seq.length(), consensus.callable_end - consensus.trim_beg);
    return consensus;
}

struct hp_read_mismatch_rates_t {
    std::vector<hp_read_mismatch_rate_t> candidate_rates;
    std::vector<hp_read_mismatch_rate_t> threshold_estimation_rates;
};

hp_read_mismatch_rates_t calculate_hp_consensus_mismatch_rates(const hp_positional_consensus_t& consensus, const std::vector<hp_read_observation_t>& observations,
    int hp_len, char hp_base) {

    hp_read_mismatch_rates_t mismatch_rates;
    if (hp_len < 0) return mismatch_rates;
    int ref_hp_base_idx = base_to_index(hp_base);
    if (ref_hp_base_idx < 0) return mismatch_rates;

    int left_5p_reads = 0, right_5p_reads = 0;
    for (const hp_read_observation_t& observation : observations) {
        if (observation.is_reverse) right_5p_reads++;
        else left_5p_reads++;
    }

    for (int read_idx = 0; read_idx < observations.size(); read_idx++) {
        const hp_read_observation_t& observation = observations[read_idx];
        bool has_independent_consensus = observation.is_reverse ? left_5p_reads >= MIN_SUPPORTING_READS : right_5p_reads >= MIN_SUPPORTING_READS;
        int sequenced_hp_base_idx = observation.is_reverse ? 3 - ref_hp_base_idx : ref_hp_base_idx;
        if (!has_independent_consensus) {
            mismatch_rates.candidate_rates.push_back({observation.ref_3p_mismatch_rate, sequenced_hp_base_idx});
            continue;
        }
        if (consensus.seq.empty()) continue;

        int tail_beg = observation.is_reverse ? 0 : observation.left_tail_len + hp_len;
        int tail_end = observation.is_reverse ? observation.left_tail_len : observation.seq.length();
        int consensus_offset = consensus.read_offsets[read_idx] - consensus.trim_beg;
        tail_beg = std::max(tail_beg, -consensus_offset);
        tail_end = std::min(tail_end, (int) consensus.seq.length() - consensus_offset);
        int compared_bases = tail_end - tail_beg;
        if (compared_bases <= 0) continue;
        int consensus_beg = consensus_offset + tail_beg;
        int mismatches = number_of_mismatches_fast(observation.seq.c_str() + tail_beg, consensus.seq.c_str() + consensus_beg, compared_bases, compared_bases);
        double mismatch_rate = double(mismatches) / compared_bases;
        mismatch_rates.candidate_rates.push_back({mismatch_rate, sequenced_hp_base_idx});
        mismatch_rates.threshold_estimation_rates.push_back({mismatch_rate, sequenced_hp_base_idx});
    }
    return mismatch_rates;
}

double median_double(std::vector<double> values) {
    if (values.empty()) return -1.0;
    std::sort(values.begin(), values.end());
    size_t middle = values.size() / 2;
    return values.size() % 2 == 1 ? values[middle] : (values[middle - 1] + values[middle]) / 2;
}

std::vector<double> estimate_hp_mismatch_rate_thresholds(const hp_mismatch_rates_by_len_t& mismatch_rates_by_hp_len, double default_threshold) {
    std::vector<double> thresholds(mismatch_rates_by_hp_len.size(), -1.0);
    std::vector<bool> threshold_computed(mismatch_rates_by_hp_len.size(), false);
    for (int hp_len = 0; hp_len < mismatch_rates_by_hp_len.size(); hp_len++) {
        if (hp_len < MIN_REF_HP_LEN) continue;
        std::vector<double> rates = mismatch_rates_by_hp_len[hp_len];
        int min_neighbor_hp_len = (hp_len + 1) / 2;
        int max_neighbor_hp_len = std::min<int>(mismatch_rates_by_hp_len.size() - 1, 2 * hp_len);
        for (int distance = 1; rates.size() < MIN_READS_FOR_HP_MM_THRESHOLD && (hp_len - distance >= min_neighbor_hp_len || hp_len + distance <= max_neighbor_hp_len); distance++) {
            int lower_hp_len = hp_len - distance;
            int upper_hp_len = hp_len + distance;
            if (lower_hp_len >= min_neighbor_hp_len) rates.insert(rates.end(), mismatch_rates_by_hp_len[lower_hp_len].begin(), mismatch_rates_by_hp_len[lower_hp_len].end());
            if (upper_hp_len <= max_neighbor_hp_len) rates.insert(rates.end(), mismatch_rates_by_hp_len[upper_hp_len].begin(), mismatch_rates_by_hp_len[upper_hp_len].end());
        }
        if (rates.size() < MIN_READS_FOR_HP_MM_THRESHOLD) continue;

        double center = median_double(rates);
        std::vector<double> absolute_deviations;
        absolute_deviations.reserve(rates.size());
        for (double rate : rates) absolute_deviations.push_back(std::abs(rate - center));
        double robust_sigma = 1.4826 * median_double(absolute_deviations);
        thresholds[hp_len] = std::min(MAX_HP_MM_THRESHOLD, std::max(default_threshold, center + HP_MM_OUTLIER_MAD_MULTIPLIER * robust_sigma));
        threshold_computed[hp_len] = true;
    }

    int nearest_smaller_computed_hp_len = -1;
    for (int hp_len = 0; hp_len < thresholds.size(); hp_len++) {
        if (threshold_computed[hp_len]) {
            nearest_smaller_computed_hp_len = hp_len;
        } else if (nearest_smaller_computed_hp_len >= 0) {
            thresholds[hp_len] = thresholds[nearest_smaller_computed_hp_len];
        } else {
            thresholds[hp_len] = default_threshold;
        }
    }
    return thresholds;
}

void write_hp_mismatch_rate_thresholds(const std::array<std::vector<double>, 4>& thresholds_by_base, const std::string& fname) {
    std::ofstream fout(fname);
    if (!fout) throw std::runtime_error("Unable to open " + fname + " for writing.");
    fout << "HP_LEN\tA\tC\tG\tT\n" << std::setprecision(17);
    size_t hp_len_count = 0;
    for (const std::vector<double>& thresholds : thresholds_by_base) hp_len_count = std::max(hp_len_count, thresholds.size());
    for (size_t hp_len = 0; hp_len < hp_len_count; hp_len++) {
        fout << hp_len;
        for (int hp_base_idx = 0; hp_base_idx < 4; hp_base_idx++) fout << "\t" << thresholds_by_base[hp_base_idx][hp_len];
        fout << "\n";
    }
    if (!fout) throw std::runtime_error("Unable to write " + fname + ".");
}

void call_aux_from_hp_consensus(std::shared_ptr<sv_t>& hp_indel, const hp_run_t& hp_run, int alt_hp_len, const hp_positional_consensus_t& consensus,
    const std::string& contig_name, char* contig_seq, hts_pos_t contig_len, config_t& config, stats_t& stats, StripedSmithWaterman::Aligner& aligner,
    const StripedSmithWaterman::Filter& filter) {

    if (consensus.seq.empty() || consensus.hp_beg < 0 || consensus.hp_end > (int) consensus.seq.size()) return;

    hts_pos_t extend = std::max<hts_pos_t>(stats.read_len, consensus.seq.size());
    hts_pos_t ref_beg = std::max<hts_pos_t>(0, hp_run.beg - extend);
    hts_pos_t ref_end = std::min<hts_pos_t>(contig_len, hp_run.end + extend);
    int left_flank_len = hp_run.beg - ref_beg;
    int right_flank_len = ref_end - hp_run.end;
    std::string alt_reference(contig_seq + ref_beg, left_flank_len);
    alt_reference += std::string(alt_hp_len, hp_run.base);
    alt_reference.append(contig_seq + hp_run.end, right_flank_len);
    int alt_hp_beg = left_flank_len;
    int alt_hp_end = alt_hp_beg + alt_hp_len;

    StripedSmithWaterman::Alignment aln;
    if (!aligner.Align(consensus.seq.c_str(), alt_reference.c_str(), alt_reference.size(), filter, &aln, 0) ||
        aln.cigar.empty() || aln.ref_begin < 0) {
        return;
    }
    if (is_clipped(aln, config.min_clip_len)) return;
    if (aln.ref_begin > alt_hp_beg || aln.ref_end < alt_hp_end - 1) return;

    auto alt_base_to_genomic = [&](int alt_pos) -> hts_pos_t {
        if (alt_pos < alt_hp_beg) return ref_beg + alt_pos;
        return hp_run.end + (alt_pos - alt_hp_end);
    };

    std::shared_ptr<sv_t> alt_hp = std::make_shared<deletion_t>(contig_name, ref_beg + alt_hp_beg - 1, ref_beg + alt_hp_end - 1,
        "", nullptr, nullptr, nullptr, nullptr);
    detect_svs_from_aln(aln, contig_name, ref_beg, consensus.seq, alt_hp, 0, 0, stats, config);

    hp_alignment_summary_t aln_summary = summarize_hp_alignment(normalized_cigar(aln), aln.ref_begin, consensus.seq,
        {alt_hp_beg, alt_hp_end}, hp_run.base);
    std::pair<hts_pos_t, hts_pos_t> callable_range = get_highq_ref_range(aln.cigar, aln.ref_begin,
        consensus.seq.length(), consensus.callable_beg, consensus.seq.length() - consensus.callable_end);
    callable_range.first = std::min<hts_pos_t>(callable_range.first, alt_hp_beg);
    callable_range.second = std::max<hts_pos_t>(callable_range.second, alt_hp_end);
    hp_indel->junction_remap_ref_beg = ref_beg + callable_range.first;
    hp_indel->junction_remap_ref_end = hp_run.end + callable_range.second - alt_hp_end;

    for (snp_t snp : alt_hp->aux_snps) {
        int alt_pos = snp.pos - ref_beg;
        if (alt_hp_beg <= alt_pos && alt_pos < alt_hp_end) continue;
        auto qpos_it = std::find(aln_summary.qpos_to_rpos.begin(), aln_summary.qpos_to_rpos.end(), alt_pos);
        if (qpos_it == aln_summary.qpos_to_rpos.end()) continue;
        int qpos = qpos_it - aln_summary.qpos_to_rpos.begin();
        snp.pos = alt_base_to_genomic(alt_pos);
        if (consensus.callable_beg <= qpos && qpos < consensus.callable_end) {
            hp_indel->aux_snps.push_back(snp);
        }
    }
    for (std::shared_ptr<sv_t> aux_indel : alt_hp->aux_indels) {
        int alt_beg = aux_indel->start + 1 - ref_beg;
        int alt_end = aux_indel->end + 1 - ref_beg;
        int qpos = aux_indel->left_anchor_aln->seq_len;
        if (aux_indel->svtype() == "INS") {
            if (alt_hp_beg <= alt_beg && alt_beg <= alt_hp_end) continue;
            hts_pos_t genomic_boundary = alt_base_to_genomic(alt_beg);
            if (genomic_boundary > 0 && qpos >= consensus.callable_beg && qpos + aux_indel->ins_seq.length() <= consensus.callable_end) {
                aux_indel->start = aux_indel->end = genomic_boundary - 1;
                hp_indel->aux_indels.push_back(aux_indel);
            }
        } else {
            if (alt_beg < alt_hp_end && alt_end > alt_hp_beg) continue;
            hts_pos_t genomic_beg = alt_base_to_genomic(alt_beg);
            hts_pos_t genomic_end = alt_base_to_genomic(alt_end - 1) + 1;
            if (genomic_beg > 0 && consensus.callable_beg < qpos && qpos < consensus.callable_end) {
                aux_indel->start = genomic_beg - 1;
                aux_indel->end = genomic_end - 1;
                hp_indel->aux_indels.push_back(aux_indel);
            }
        }
    }
}

void add_rescued_hp_read(hp_run_t& hp_run, const mate_info_t& mate,
    bool anchor_is_reverse, char* contig_seq, hts_pos_t contig_len,
    int read_len, int min_clip_len, double max_seq_error, StripedSmithWaterman::Aligner& aligner,
    const StripedSmithWaterman::Filter& filter) {

    hts_pos_t extend = std::max(0, read_len - 1);
    hts_pos_t allele_beg = std::max((hts_pos_t) 0, hp_run.beg - extend);
    hts_pos_t allele_end = std::min(contig_len, hp_run.end + extend);
    std::string ref_allele(contig_seq + allele_beg, allele_end - allele_beg);
    hts_pair_pos_t allele_hp_range = {hp_run.beg - allele_beg, hp_run.end - allele_beg};

    std::string mate_seq = mate.seq;
    std::vector<uint8_t> mate_quals = decode_qualities(mate.qual);
    if (!anchor_is_reverse) {
        rc(mate_seq);
        std::reverse(mate_quals.begin(), mate_quals.end());
    }
    StripedSmithWaterman::Alignment aln;
    if (!aligner.Align(mate_seq.c_str(), ref_allele.c_str(), ref_allele.length(), filter, &aln, 0)) {
        return;
    }
    bool aln_as_rev = !anchor_is_reverse;
    int left_clip = aln.query_begin;
    int right_clip = mate_seq.length() - aln.query_end - 1;
    if (aln.ref_begin < 0 || aln.cigar.empty() || (!aln_as_rev && left_clip > 0) || (aln_as_rev && right_clip > 0)) {
        return;
    }

    bool has_no_left_indel = five_p_evidence_permits_iterative_hp_len_estimation(hp_run.left_side_5p_reads, hp_run.left_side_5p_indel_reads);
    bool has_no_right_indel = five_p_evidence_permits_iterative_hp_len_estimation(hp_run.right_side_5p_reads, hp_run.right_side_5p_indel_reads);
    bp_support_read_t rescued_read;
    rescued_read.mapq = mate.mapq;
    rescued_read.seq = mate_seq;
    rescued_read.mate_is_reverse = anchor_is_reverse;
    hp_read_info_t hp_read_info = calculate_hp_read_info(aln, mate_seq,
        allele_hp_range, hp_run.base, &ref_allele[0], ref_allele.length(), aln_as_rev, rescued_read, 0,
        has_no_left_indel, has_no_right_indel, max_seq_error);
    if (!is_usable_hp_read(hp_read_info, min_clip_len, max_seq_error)) return;
    hp_run.usable_reads++;
    hp_run.hp_len_counts[hp_read_info.hp_len]++;
    bool is_reverse = !hp_read_info.read.mate_is_reverse;
    int left_tail_len = is_reverse ? hp_read_info.tail_3p_len : hp_read_info.tail_5p_len;
    hp_run.observations_by_hp_len[hp_read_info.hp_len].push_back({mate_seq, mate_quals, left_tail_len, is_reverse, double(hp_read_info.tail_3p_mismatches) / hp_read_info.tail_3p_len});
}

hp_chunk_result_t find_hp_indels_for_chunk(int id, size_t contig_id, std::string contig_name,
    char* contig_seq, hts_pos_t contig_len, hts_pos_t chunk_beg, hts_pos_t chunk_end,
    config_t* config, stats_t* stats, bam_pool_t* bam_pool,
    StripedSmithWaterman::Aligner& aligner, const StripedSmithWaterman::Filter& filter) {

    hp_chunk_result_t result;
    if (contig_len == 0) return result;

    std::vector<hp_run_t> hp_runs = find_hp_runs(contig_seq, contig_len, chunk_beg, chunk_end);
    if (hp_runs.empty()) return result;

    std::vector<hts_pos_t> hp_run_ends;
    std::vector<hts_pos_t> hp_run_begs;
    hp_run_ends.reserve(hp_runs.size());
    hp_run_begs.reserve(hp_runs.size());
    for (const hp_run_t& hp_run : hp_runs) {
        hp_run_begs.push_back(hp_run.beg);
        hp_run_ends.push_back(hp_run.end);
    }

    open_samFile_t* alignment_file = bam_pool->get_bam_reader(id);
    read_mates(contig_id);
    try {
        const mate_map_t& mateseqs_w_mapq_chr = mateseqs_w_mapq[contig_id];
        int min_anchor_mapq = config->high_confidence_mapq;
        hts_pos_t query_beg = std::max((hts_pos_t) 0, hp_runs.front().beg - stats->max_is);
        hts_pos_t query_end = std::min(contig_len, hp_runs.back().end + stats->max_is);
        std::stringstream region_ss;
        region_ss << contig_name << ":" << query_beg + 1 << "-" << query_end;
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iter(sam_itr_querys(alignment_file->idx, alignment_file->header, region_ss.str().c_str()), &hts_itr_destroy);
        if (!iter) throw std::runtime_error("Unable to query alignments for " + contig_name + ".");

        std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), &bam_destroy1);
        if (!read) throw std::runtime_error("Unable to allocate an alignment record.");

        int read_status;
        while ((read_status = sam_itr_next(alignment_file->file, iter.get(), read.get())) >= 0) {
            if (is_unmapped(read.get()) || !is_primary(read.get()) || read->core.l_qseq <= 0) continue;
            if (!is_proper_pair(read.get(), stats->min_is, stats->max_is)) continue;

            hts_pos_t read_beg = read->core.pos;
            hts_pos_t read_end = bam_endpos(read.get());
            size_t hp_idx = std::upper_bound(hp_run_ends.begin(), hp_run_ends.end(), read_beg) - hp_run_ends.begin();
            for (; hp_idx < hp_runs.size() && hp_runs[hp_idx].beg < read_end; hp_idx++) {
                hp_run_t& hp_run = hp_runs[hp_idx];
                hp_adjacent_indel_info_t indel_info = get_adjacent_indel_info(read.get(), {hp_run.beg, hp_run.end}, stats->read_len/2);
                bool is_reverse = bam_is_rev(read.get());
                const hp_side_indel_info_t& five_p_info = is_reverse ? indel_info.right : indel_info.left;
                int& side_5p_reads = is_reverse ? hp_run.right_side_5p_reads : hp_run.left_side_5p_reads;
                int& side_5p_indel_reads = is_reverse ? hp_run.right_side_5p_indel_reads : hp_run.left_side_5p_indel_reads;
                if (five_p_info.aligned_len >= config->min_clip_len) {
                    side_5p_reads++;
                    if (five_p_info.indel_len != 0) side_5p_indel_reads++;
                }
            }
        }
        if (read_status < -1) throw std::runtime_error("Error while reading alignments for " + contig_name + ".");

        iter.reset(sam_itr_querys(alignment_file->idx, alignment_file->header, region_ss.str().c_str()));
        if (!iter) throw std::runtime_error("Unable to query alignments for " + contig_name + ".");
        while ((read_status = sam_itr_next(alignment_file->file, iter.get(), read.get())) >= 0) {
            if (is_unmapped(read.get()) || !is_primary(read.get()) || read->core.l_qseq <= 0) continue;

            hts_pos_t read_beg = read->core.pos;
            hts_pos_t read_end = bam_endpos(read.get());
            size_t hp_idx = std::upper_bound(hp_run_ends.begin(), hp_run_ends.end(), read_beg) - hp_run_ends.begin();
            for (; hp_idx < hp_runs.size() && hp_runs[hp_idx].beg < read_end; hp_idx++) {
                hp_run_t& hp_run = hp_runs[hp_idx];

                hts_pair_pos_t hp_range = {hp_run.beg, hp_run.end};
                hts_pos_t extend = read->core.l_qseq - 1;
                hts_pos_t ref_allele_beg = std::max((hts_pos_t) 0, hp_run.beg - extend);
                hts_pos_t ref_allele_end = std::min(contig_len, hp_run.end + extend);
                hts_pair_pos_t ref_allele_hp_range = {hp_run.beg - ref_allele_beg, hp_run.end - ref_allele_beg};
                bool has_no_left_indel = five_p_evidence_permits_iterative_hp_len_estimation(hp_run.left_side_5p_reads, hp_run.left_side_5p_indel_reads);
                bool has_no_right_indel = five_p_evidence_permits_iterative_hp_len_estimation(hp_run.right_side_5p_reads, hp_run.right_side_5p_indel_reads);
                hp_read_info_t hp_read_info = calculate_hp_read_info(read.get(), hp_range, hp_run.base,
                    contig_seq, contig_len, contig_seq + ref_allele_beg, ref_allele_end - ref_allele_beg,
                    ref_allele_hp_range, has_no_left_indel, has_no_right_indel, 0, config->max_seq_error);
                if (!is_usable_hp_read(hp_read_info, config->min_clip_len, config->max_seq_error)) continue;
                hp_run.usable_reads++;
                hp_run.hp_len_counts[hp_read_info.hp_len]++;
                bool is_reverse = !hp_read_info.read.mate_is_reverse;
                int left_tail_len = is_reverse ? hp_read_info.tail_3p_len : hp_read_info.tail_5p_len;
                const uint8_t* bam_quals = bam_get_qual(read.get());
                std::vector<uint8_t> quals(bam_quals, bam_quals + read->core.l_qseq);
                std::replace(quals.begin(), quals.end(), uint8_t(255), uint8_t(0));
                hp_run.observations_by_hp_len[hp_read_info.hp_len].push_back({hp_read_info.read.seq, quals, left_tail_len, is_reverse, double(hp_read_info.tail_3p_mismatches) / hp_read_info.tail_3p_len});
            }

            if (read->core.qual < min_anchor_mapq || !is_dc_pair(read.get()) || mateseqs_w_mapq_chr.empty()) continue;
            std::string mate_qname = get_mate_lookup_qname(read.get());
            auto mate_it = mateseqs_w_mapq_chr.find(mate_qname);
            if (mate_it == mateseqs_w_mapq_chr.end()) continue;

            if (!bam_is_rev(read.get())) {
                hts_pos_t min_hp_beg = read->core.pos + stats->read_len / 2;
                hts_pos_t max_hp_beg = read->core.pos + stats->max_is;
                size_t rescue_hp_idx = std::lower_bound(hp_run_begs.begin(), hp_run_begs.end(), min_hp_beg) - hp_run_begs.begin();
                for (; rescue_hp_idx < hp_runs.size() && hp_runs[rescue_hp_idx].beg <= max_hp_beg; rescue_hp_idx++) {
                    add_rescued_hp_read(hp_runs[rescue_hp_idx], mate_it->second, false,
                        contig_seq, contig_len, stats->read_len, config->min_clip_len, config->max_seq_error,
                        aligner, filter);
                }
            } else {
                hts_pos_t min_hp_end = read_end - stats->max_is;
                hts_pos_t max_hp_end = read_end - stats->read_len / 2;
                size_t rescue_hp_idx = std::lower_bound(hp_run_ends.begin(), hp_run_ends.end(), min_hp_end) - hp_run_ends.begin();
                for (; rescue_hp_idx < hp_runs.size() && hp_runs[rescue_hp_idx].end <= max_hp_end; rescue_hp_idx++) {
                    add_rescued_hp_read(hp_runs[rescue_hp_idx], mate_it->second, true,
                        contig_seq, contig_len, stats->read_len, config->min_clip_len, config->max_seq_error,
                        aligner, filter);
                }
            }
        }
        if (read_status < -1) throw std::runtime_error("Error while reading alignments for " + contig_name + ".");

        int max_hp_reads_per_locus = 2 * stats->get_max_depth(contig_name);
        for (const hp_run_t& hp_run : hp_runs) {
            if (hp_run.usable_reads > max_hp_reads_per_locus) continue;
            int ref_hp_len = hp_run.end - hp_run.beg;
            for (const auto& hp_len_count : hp_run.hp_len_counts) {
                int alt_hp_len = hp_len_count.first;
                int supporting_reads = hp_len_count.second;
                if (supporting_reads < MIN_SUPPORTING_READS) continue;

                auto observations_it = hp_run.observations_by_hp_len.find(alt_hp_len);
                if (observations_it == hp_run.observations_by_hp_len.end()) continue;
                hp_positional_consensus_t consensus = build_hp_positional_consensus(observations_it->second, alt_hp_len);
                hp_read_mismatch_rates_t read_mismatch_rates = calculate_hp_consensus_mismatch_rates(consensus, observations_it->second, alt_hp_len, hp_run.base);
                for (const hp_read_mismatch_rate_t& read_rate : read_mismatch_rates.candidate_rates) {
                    hp_mismatch_rates_by_len_t& rates_by_hp_len = result.threshold_estimation_rates_by_hp_len_and_base[read_rate.sequenced_hp_base_idx];
                    if (rates_by_hp_len.size() <= alt_hp_len) rates_by_hp_len.resize(alt_hp_len + 1);
                }
                for (const hp_read_mismatch_rate_t& read_rate : read_mismatch_rates.threshold_estimation_rates) result.threshold_estimation_rates_by_hp_len_and_base[read_rate.sequenced_hp_base_idx][alt_hp_len].push_back(read_rate.rate);

                if (hp_run.beg == 0 || alt_hp_len == ref_hp_len) continue;

                int hp_len_diff = alt_hp_len - ref_hp_len;
                if (std::abs(hp_len_diff) < config->min_sv_size) continue;
                hts_pos_t anchor_pos = hp_run.beg - 1;
                std::shared_ptr<sv_t> hp_indel;
                if (hp_len_diff > 0) {
                    hp_indel = std::make_shared<insertion_t>(contig_name, anchor_pos, anchor_pos,
                        std::string(hp_len_diff, hp_run.base), nullptr, nullptr, nullptr, nullptr);
                } else {
                    hp_indel = std::make_shared<deletion_t>(contig_name, anchor_pos,
                        anchor_pos - hp_len_diff, "", nullptr, nullptr, nullptr, nullptr);
                }
                hp_indel->source = "HP";
                hp_indel->junction_remap_ref_beg = hp_indel->start;
                hp_indel->junction_remap_ref_end = hp_indel->end + 1;
                hp_indel->hp_ref_beg = hp_run.beg;
                hp_indel->hp_ref_end = hp_run.end;
                hp_indel->sample_info.alt1_hp_len_mode = alt_hp_len;
                call_aux_from_hp_consensus(hp_indel, hp_run, alt_hp_len, consensus, contig_name, contig_seq, contig_len, *config, *stats, aligner, filter);
                result.hp_indel_3p_mismatch_rates.push_back(std::move(read_mismatch_rates.candidate_rates));
                result.hp_indels.push_back(hp_indel);
            }
        }
    } catch (...) {
        release_mates(contig_id);
        throw;
    }
    release_mates(contig_id);
    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: find_hp_indels <workdir> <reference.fa> <alignments.bam|cram> <output.vcf.gz>\n";
        return 1;
    }

    workdir = argv[1];
    std::string reference_fname = argv[2];
    std::string alignment_fname = argv[3];
    std::string out_vcf_fname = argv[4];

    config_t config;
    config.parse(workdir + "/config.txt");

    stats_t stats;
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    contig_map_t contig_map(workdir);

    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(reference_fname);

    bam_pool_t bam_pool(config.threads, alignment_fname, reference_fname);
    mateseqs_w_mapq.resize(contig_map.size());
    active_threads_per_chr = std::vector<int>(contig_map.size());
    mutex_per_chr = std::vector<std::mutex>(contig_map.size());
    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Filter filter;

    std::vector<std::future<hp_chunk_result_t>> futures;
    ctpl::thread_pool thread_pool(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        hts_pos_t contig_len = chr_seqs.get_len(contig_name);
        for (hts_pos_t chunk_beg = 0; chunk_beg < contig_len; chunk_beg += CHUNK_SIZE) {
            futures.push_back(thread_pool.push(find_hp_indels_for_chunk, contig_id, contig_name, chr_seqs.get_seq(contig_name), contig_len,
                chunk_beg, std::min(contig_len, chunk_beg + CHUNK_SIZE), &config, &stats, &bam_pool, std::ref(aligner), std::cref(filter)));
        }
    }
    thread_pool.stop(true);

    std::vector<std::shared_ptr<sv_t>> hp_indels;
    std::vector<std::vector<hp_read_mismatch_rate_t>> hp_indel_3p_mismatch_rates;
    hp_mismatch_rates_by_len_and_base_t threshold_estimation_rates_by_hp_len_and_base;
    for (std::future<hp_chunk_result_t>& future : futures) {
        hp_chunk_result_t chunk_result = future.get();
        hp_indels.insert(hp_indels.end(), chunk_result.hp_indels.begin(), chunk_result.hp_indels.end());
        hp_indel_3p_mismatch_rates.insert(hp_indel_3p_mismatch_rates.end(), std::make_move_iterator(chunk_result.hp_indel_3p_mismatch_rates.begin()), std::make_move_iterator(chunk_result.hp_indel_3p_mismatch_rates.end()));
        for (int hp_base_idx = 0; hp_base_idx < 4; hp_base_idx++) {
            hp_mismatch_rates_by_len_t& threshold_estimation_rates_by_hp_len = threshold_estimation_rates_by_hp_len_and_base[hp_base_idx];
            const hp_mismatch_rates_by_len_t& chunk_threshold_estimation_rates_by_hp_len = chunk_result.threshold_estimation_rates_by_hp_len_and_base[hp_base_idx];
            if (threshold_estimation_rates_by_hp_len.size() < chunk_threshold_estimation_rates_by_hp_len.size()) threshold_estimation_rates_by_hp_len.resize(chunk_threshold_estimation_rates_by_hp_len.size());
            for (int hp_len = 0; hp_len < chunk_threshold_estimation_rates_by_hp_len.size(); hp_len++) {
                threshold_estimation_rates_by_hp_len[hp_len].insert(threshold_estimation_rates_by_hp_len[hp_len].end(), chunk_threshold_estimation_rates_by_hp_len[hp_len].begin(), chunk_threshold_estimation_rates_by_hp_len[hp_len].end());
            }
        }
    }
    size_t hp_len_count = 0;
    for (const hp_mismatch_rates_by_len_t& rates_by_hp_len : threshold_estimation_rates_by_hp_len_and_base) hp_len_count = std::max(hp_len_count, rates_by_hp_len.size());
    for (hp_mismatch_rates_by_len_t& rates_by_hp_len : threshold_estimation_rates_by_hp_len_and_base) rates_by_hp_len.resize(hp_len_count);
    std::array<std::vector<double>, 4> hp_mismatch_rate_thresholds_by_base;
    for (int hp_base_idx = 0; hp_base_idx < 4; hp_base_idx++) hp_mismatch_rate_thresholds_by_base[hp_base_idx] = estimate_hp_mismatch_rate_thresholds(threshold_estimation_rates_by_hp_len_and_base[hp_base_idx], config.max_seq_error);
    write_hp_mismatch_rate_thresholds(hp_mismatch_rate_thresholds_by_base, workdir + "/" + HP_MISMATCH_RATE_THRESHOLDS_FILENAME);
    std::vector<std::shared_ptr<sv_t>> filtered_hp_indels;
    for (int hp_indel_idx = 0; hp_indel_idx < hp_indels.size(); hp_indel_idx++) {
        std::shared_ptr<sv_t>& hp_indel = hp_indels[hp_indel_idx];
        int hp_len = hp_indel->sample_info.alt1_hp_len_mode;
        if (hp_len < 0) continue;

        int passing_reads = 0;
        for (const hp_read_mismatch_rate_t& read_rate : hp_indel_3p_mismatch_rates[hp_indel_idx]) {
            int hp_base_idx = read_rate.sequenced_hp_base_idx;
            if (hp_base_idx < 0 || hp_len >= hp_mismatch_rate_thresholds_by_base[hp_base_idx].size()) continue;
            if (read_rate.rate <= hp_mismatch_rate_thresholds_by_base[hp_base_idx][hp_len]) {
                passing_reads++;
            }
        }
        if (passing_reads < MIN_SUPPORTING_READS) continue;
        filtered_hp_indels.push_back(hp_indel);
    }
    hp_indels = std::move(filtered_hp_indels);
    std::string command;
    for (int i = 0; i < argc; i++) {
        if (!command.empty()) command += " ";
        command += argv[i];
    }
    std::unique_ptr<bcf_hdr_t, decltype(&bcf_hdr_destroy)> hdr( generate_vcf_header(chr_seqs, "", config, command), &bcf_hdr_destroy);
    std::sort(hp_indels.begin(), hp_indels.end(), [&hdr](const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) {
        int a_rid = bcf_hdr_name2id(hdr.get(), a->chr.c_str());
        int b_rid = bcf_hdr_name2id(hdr.get(), b->chr.c_str());
        if (a_rid != b_rid) return a_rid < b_rid;
        return sv_output_order(a, b);
    });
    std::unique_ptr<htsFile, decltype(&hts_close)> out(bcf_open(out_vcf_fname.c_str(), "wz"), &hts_close);
    if (!out) throw std::runtime_error("Unable to open " + out_vcf_fname + " for writing.");
    if (hts_set_threads(out.get(), std::max(1, std::min(config.threads, 4))) != 0) {
        throw std::runtime_error("Unable to enable threaded VCF compression.");
    }
    if (bcf_hdr_write(out.get(), hdr.get()) != 0) {
        throw std::runtime_error("Unable to write the VCF header to " + out_vcf_fname + ".");
    }
    std::unique_ptr<bcf1_t, decltype(&bcf_destroy)> record(bcf_init(), &bcf_destroy);
    if (!record) throw std::runtime_error("Unable to allocate a VCF record.");
    for (size_t i = 0; i < hp_indels.size(); i++) {
        std::shared_ptr<sv_t>& hp_indel = hp_indels[i];
        hp_indel->id = "HP_" + hp_indel->svtype() + "_" + std::to_string(i + 1);
        hp_indel->sample_info.gt.clear();
        sv2bcf(hdr.get(), record.get(), hp_indel.get(), chr_seqs.get_seq(hp_indel->chr));
        int hp_ref_range[] = {(int) hp_indel->hp_ref_beg + 1, (int) hp_indel->hp_ref_end + 1};
        bcf_update_info_int32(hdr.get(), record.get(), "HP_REF_RANGE", hp_ref_range, 2);
        if (bcf_write(out.get(), hdr.get(), record.get()) != 0) {
            throw std::runtime_error("Unable to write VCF record " + hp_indel->id + ".");
        }
    }
    out.reset();

    if (tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf) != 0) {
        throw std::runtime_error("Unable to index " + out_vcf_fname + ".");
    }
}
