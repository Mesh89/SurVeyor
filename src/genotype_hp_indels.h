#ifndef GENOTYPE_HP_INDELS_H
#define GENOTYPE_HP_INDELS_H

#include <memory>

#include "htslib/sam.h"
#include "types.h"
#include "sam_utils.h"
#include "utils.h"
#include "stat_tests.h"
#include "../libs/ssw_cpp.h"

#include "genotype.h"
#include "hp_read_info.h"
#include "var_utils.h"

const double MAX_TAIL_MISMATCH_RATE = 0.2; 

// The goal of this function is to identify reads that are unlikely to have a correct hp len
// For now, let's use a simple heuristic: if a read has a hp run len <= 5 bp, it must be a good read
std::vector<bool> get_valid_reads_mask(const std::vector<hp_read_info_t>& hp_read_infos, int min_tail_len, double max_mismatch_rate) {

    std::vector<bool> valid_reads_mask;
    valid_reads_mask.reserve(hp_read_infos.size());
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (hp_read_info.hp_len == UNDEFINED_HP_LEN) {
            valid_reads_mask.push_back(false);
        } else if (hp_read_info.aligned_5p_tail_len < min_tail_len) {
            valid_reads_mask.push_back(false);
        } else if (hp_read_info.hp_len <= 5) {
            valid_reads_mask.push_back(hp_read_info.is_good_read(min_tail_len, max_mismatch_rate));
        } else {
            valid_reads_mask.push_back(true);
        }
    }
    return valid_reads_mask;
}


// Find the mode of HP lengths, optionally restricted to good reads
// Good reads = <20% mismatches on both tails
int find_hp_len_mode(const std::vector<hp_read_info_t>& hp_read_infos, int min_tail_len, double max_mismatch_rate, bool good_reads_only) {
    std::unordered_map<int, int> hp_len_counts;
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (!good_reads_only || hp_read_info.is_good_read(min_tail_len, max_mismatch_rate)) {
            hp_len_counts[hp_read_info.hp_len]++;
        }
    }

    int mode_hp_len = -1, max_count = 0;
    for (const auto& kv : hp_len_counts) {
        if (kv.second > max_count || (kv.second == max_count && kv.first < mode_hp_len)) {
            mode_hp_len = kv.first;
            max_count = kv.second;
        }
    }
    if (max_count > 0) {
        return mode_hp_len;
    }
    return -1; // no good reads
}

double quantile_linear_interp(std::vector<int>& sorted_values, double q) {
    if (sorted_values.empty()) {
        return -1.0;
    }
    if (sorted_values.size() == 1) {
        return sorted_values[0];
    }

    double pos = q * (sorted_values.size() - 1);
    int lo = std::floor(pos);
    int hi = std::ceil(pos);
    if (lo == hi) {
        return sorted_values[lo];
    }

    double frac = pos - lo;
    return sorted_values[lo] + frac * (sorted_values[hi] - sorted_values[lo]);
}

double find_hp_len_iqr(const std::vector<hp_read_info_t>& hp_read_infos, int min_tail_len, double max_mismatch_rate, bool good_reads_only) {
    std::vector<int> hp_lens;
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (!good_reads_only || hp_read_info.is_good_read(min_tail_len, max_mismatch_rate)) {
            hp_lens.push_back(hp_read_info.hp_len);
        }
    }

    if (hp_lens.empty()) {
        return -1.0;
    }

    std::sort(hp_lens.begin(), hp_lens.end());
    return quantile_linear_interp(hp_lens, 0.75) - quantile_linear_interp(hp_lens, 0.25);
}

void set_hp_tail_mismatch_rates(const std::vector<hp_read_info_t>& hp_read_infos, int min_tail_len, double max_mismatch_rate,
    bool good_reads_only, double& avg_5p_mismatch_rate, double& avg_3p_mismatch_rate) {
    
    double sum_5p_mismatch_rate = 0.0;
    double sum_3p_mismatch_rate = 0.0;
    int n_reads = 0;

    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (good_reads_only && !hp_read_info.is_good_read(min_tail_len, max_mismatch_rate)) continue;

        sum_5p_mismatch_rate += double(hp_read_info.tail_5p_mismatches) / hp_read_info.tail_5p_len;
        sum_3p_mismatch_rate += double(hp_read_info.tail_3p_mismatches) / hp_read_info.tail_3p_len;
        n_reads++;
    }

    if (n_reads == 0) {
        return;
    }

    avg_5p_mismatch_rate = sum_5p_mismatch_rate / n_reads;
    avg_3p_mismatch_rate = sum_3p_mismatch_rate / n_reads;
}

void set_hp_read_mapq_stats(const std::vector<bp_support_read_t>& reads, int& min_mapq, int& max_mapq, double& avg_mapq, double& stddev_mapq) {
    if (reads.empty()) {
        return;
    }

    std::vector<int> mapqs;
    mapqs.reserve(reads.size());
    for (const bp_support_read_t& read : reads) {
        mapqs.push_back(read.mapq);
    }

    min_mapq = *std::min_element(mapqs.begin(), mapqs.end());
    max_mapq = *std::max_element(mapqs.begin(), mapqs.end());
    avg_mapq = mean(mapqs);
    stddev_mapq = stddev(mapqs);
}

void set_hp_read_mapq_stats(const std::vector<std::shared_ptr<bam1_t>>& reads, int& min_mapq, int& max_mapq, double& avg_mapq, double& stddev_mapq) {
    std::vector<bp_support_read_t> support_reads;
    support_reads.reserve(reads.size());
    for (const std::shared_ptr<bam1_t>& read : reads) {
        support_reads.emplace_back(read.get());
    }
    set_hp_read_mapq_stats(support_reads, min_mapq, max_mapq, avg_mapq, stddev_mapq);
}

std::vector<std::vector<hp_read_info_t>> cluster_reads_by_two_modes(const std::vector<hp_read_info_t>& hp_read_infos,
    int min_tail_len, double max_mismatch_rate) {
    int n = hp_read_infos.size();
    if (n == 0) {
        return {};
    }

    std::unordered_map<int, int> good_counts;
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (hp_read_info.is_good_read(min_tail_len, max_mismatch_rate)) {
            good_counts[hp_read_info.hp_len]++;
        }
    }

    if (good_counts.empty()) {
        return {};
    }

    // Pick seed modes from good reads only, breaking count ties toward the smaller length.
    std::vector<std::pair<int, int>> count_items(good_counts.begin(), good_counts.end());
    std::sort(count_items.begin(), count_items.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;
    });

    if (count_items.size() == 1) {
        return {hp_read_infos};
    }

    int a_mode = count_items[0].first, b_mode = count_items[1].first;
    if (a_mode > b_mode) std::swap(a_mode, b_mode);

    // Absorb every read into the nearest good-read mode; distance ties go to the larger good-read seed.
    std::vector<int> a_idx, b_idx;
    for (int i = 0; i < n; i++) {
        int x = hp_read_infos[i].hp_len;
        int a_dist = abs(x - a_mode);
        int b_dist = abs(x - b_mode);
        if (a_dist < b_dist) {
            a_idx.push_back(i);
        } else if (b_dist < a_dist) {
            b_idx.push_back(i);
        } else if (good_counts[a_mode] >= good_counts[b_mode]) {
            a_idx.push_back(i);
        } else {
            b_idx.push_back(i);
        }
    }

    std::vector<hp_read_info_t> a_reads, b_reads;
    for (int i : a_idx) a_reads.push_back(hp_read_infos[i]);
    for (int i : b_idx) b_reads.push_back(hp_read_infos[i]);
    return {a_reads, b_reads};
}

struct hp_allele_cluster_t {
    int allele_idx;
    int allele_len;
    std::vector<hp_read_info_t> reads;

    hp_allele_cluster_t(int allele_idx, int allele_len) : allele_idx(allele_idx), allele_len(allele_len) {}
};

int nearest_hp_allele_idx(int hp_len, const std::string& read_seq, 
    std::vector<std::unique_ptr<char[]>>& alt_alleles, std::vector<int>& alt_allele_lens,
    std::vector<int>& hp_run_lens, 
    StripedSmithWaterman::Aligner& aligner, const std::vector<int>* allowed_allele_idxs = NULL) {
    
    std::vector<int> all_allele_idxs;
    if (allowed_allele_idxs == NULL) {
        for (int i = 0; i < hp_run_lens.size(); i++) all_allele_idxs.push_back(i);
        allowed_allele_idxs = &all_allele_idxs;
    }

    int best_dist = INT32_MAX;
    std::vector<int> best_idxs;
    for (int idx : *allowed_allele_idxs) {
        int dist = abs(hp_len - hp_run_lens[idx]);
        if (dist < best_dist) {
            best_dist = dist;
            best_idxs.clear();
            best_idxs.push_back(idx);
        } else if (dist == best_dist) {
            best_idxs.push_back(idx);
        }
    }

    if (best_idxs.empty()) return -1;
    if (best_idxs.size() == 1) return best_idxs[0];

    // If there are multiple alleles with the same distance, choose shortest HP run length
    // This is because empirically, I have observed that sequencing tend to overestimate the HP run length rather than underestimate it
    bool same_hp_run_len = true;
    for (int idx : best_idxs) {
        if (hp_run_lens[idx] != hp_run_lens[best_idxs[0]]) {
            same_hp_run_len = false;
            break;
        }
    }
    if (!same_hp_run_len) {
        int best_idx = best_idxs[0];
        for (int idx : best_idxs) {
            if (hp_run_lens[idx] < hp_run_lens[best_idx]) {
                best_idx = idx;
            }
        }
        return best_idx;
    }

    // If there optimal HP run len has multiple candidate hp indels, pick the one with the best alignment score
    if (read_seq.empty()) {
        throw std::runtime_error("Missing read sequence for HP allele tie-break.");
    }

    StripedSmithWaterman::Alignment aln;
    StripedSmithWaterman::Filter filter_score_only(false, false, 0, 32767);

    int best_idx = -1;
    int best_score = INT32_MIN;
    for (int idx : best_idxs) {
        aligner.Align(read_seq.c_str(), alt_alleles[idx].get(), alt_allele_lens[idx], filter_score_only, &aln, 0);
        if (aln.sw_score > best_score ||
            (aln.sw_score == best_score && (best_idx == -1 || hp_run_lens[idx] < hp_run_lens[best_idx]))) {
            best_score = aln.sw_score;
            best_idx = idx;
        }
    }
    return best_idx;
}

bool candidate_haplotypes_support_exact_ref(const std::string& read_seq,
    const char* ref_haplotype, int ref_haplotype_len, hts_pair_pos_t ref_hp_range,
    const char* alt_haplotype, int alt_haplotype_len, hts_pair_pos_t alt_hp_range) {

    ungapped_aln_t ref_aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(), ref_haplotype, ref_haplotype_len, 0);
    ungapped_aln_t alt_aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(), alt_haplotype, alt_haplotype_len, 0);

    bool ref_spans_hp = ref_aln.ref_begin <= ref_hp_range.beg && ref_aln.ref_end >= ref_hp_range.end;
    bool alt_spans_hp = alt_aln.ref_begin <= alt_hp_range.beg && alt_aln.ref_end >= alt_hp_range.end;
    return ref_spans_hp && (!alt_spans_hp || ref_aln.score > alt_aln.score);
}

std::vector<hp_allele_cluster_t> cluster_reads_by_nearest_allele_len(const std::vector<hp_read_info_t>& hp_read_infos,
    std::vector<std::unique_ptr<char[]>>& alt_alleles, std::vector<int>& alt_allele_lens, std::vector<int>& hp_run_lens, StripedSmithWaterman::Aligner& aligner) {

    if (hp_read_infos.empty() || hp_run_lens.empty()) return {};

    // First, assign each read to its nearest allele length
    std::vector<int> allele_counts(hp_run_lens.size(), 0);
    std::vector<long long> allele_dist_sums(hp_run_lens.size(), 0);
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        int allele_idx = nearest_hp_allele_idx(hp_read_info.hp_len, hp_read_info.read.seq, alt_alleles, alt_allele_lens, hp_run_lens, aligner);
        for (int i = 0; i < hp_run_lens.size(); i++) {
            if (hp_run_lens[i] == hp_run_lens[allele_idx]) {
                allele_counts[i]++;
                allele_dist_sums[i] += abs(hp_read_info.hp_len - hp_run_lens[i]);
            }
        }
    }

    // Then pick the (max 2) hp run lens with the most assigned reads, breaking count ties toward smaller total distance, then smaller length
    std::vector<int> allele_idxs;
    for (int i = 0; i < hp_run_lens.size(); i++) {
        if (allele_counts[i] > 0) allele_idxs.push_back(i);
    }
    std::sort(allele_idxs.begin(), allele_idxs.end(), [&](int a, int b) {
        if (allele_counts[a] != allele_counts[b]) return allele_counts[a] > allele_counts[b];
        if (allele_dist_sums[a] != allele_dist_sums[b]) return allele_dist_sums[a] < allele_dist_sums[b];
        if (hp_run_lens[a] != hp_run_lens[b]) return hp_run_lens[a] < hp_run_lens[b];
        return a < b;
    });

    std::vector<int> chosen_allele_idxs;
    chosen_allele_idxs.push_back(allele_idxs[0]);
    for (int allele_idx : allele_idxs) {
        if (hp_run_lens[allele_idx] != hp_run_lens[chosen_allele_idxs[0]]) {
            chosen_allele_idxs.push_back(allele_idx);
            break;
        }
    }

    std::sort(chosen_allele_idxs.begin(), chosen_allele_idxs.end(), [&](int a, int b) {
        if (hp_run_lens[a] != hp_run_lens[b]) return hp_run_lens[a] < hp_run_lens[b];
        return a < b;
    });

    // Finally, form max 2 clusters of reads by assigning reads to the max 2 retained alleles
    std::vector<hp_allele_cluster_t> clusters;
    for (int allele_idx : chosen_allele_idxs) {
        clusters.push_back(hp_allele_cluster_t(allele_idx, hp_run_lens[allele_idx]));
    }
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        int allele_idx = nearest_hp_allele_idx(hp_read_info.hp_len, hp_read_info.read.seq, alt_alleles, alt_allele_lens, hp_run_lens, aligner, &chosen_allele_idxs);
        for (hp_allele_cluster_t& cluster : clusters) {
            if (cluster.allele_idx == allele_idx) {
                cluster.reads.push_back(hp_read_info);
                break;
            }
        }
    }

    return clusters;
}

std::vector<hp_allele_cluster_t> cluster_reads_by_nearest_allele_len_reassign(const std::vector<hp_read_info_t>& hp_read_infos,
    std::vector<int>& hp_run_lens, std::vector<sv_t*>& hp_indels, evidence_map_t* evidence_map, int seed) {

    std::vector<hp_allele_cluster_t> clusters;
    for (int allele_idx = 0; allele_idx < hp_indels.size(); allele_idx++) {
        clusters.push_back(hp_allele_cluster_t(allele_idx, hp_run_lens[allele_idx]));
    }
    int ref_allele_idx = hp_indels.size();
    clusters.push_back(hp_allele_cluster_t(ref_allele_idx, hp_run_lens[ref_allele_idx]));

    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        bool assigned = false;
        for (int allele_idx = 0; allele_idx < hp_indels.size(); allele_idx++) {
            if (evidence_map->is_read_assigned_to_this_sv(hp_read_info.read, hp_indels[allele_idx])) {
                clusters[allele_idx].reads.push_back(hp_read_info);
                assigned = true;
                break;
            }
        }
        if (!assigned) {
            clusters[ref_allele_idx].reads.push_back(hp_read_info);
        }
    }

    std::stable_sort(clusters.begin(), clusters.end(), [](const hp_allele_cluster_t& a, const hp_allele_cluster_t& b) {
        return a.reads.size() > b.reads.size();
    });
    clusters.erase(std::remove_if(clusters.begin(), clusters.end(), [](const hp_allele_cluster_t& cluster) {
        return cluster.reads.empty();
    }), clusters.end());
    if (clusters.size() > 2) {
        std::vector<hp_read_info_t> discarded_reads;
        for (int i = 2; i < clusters.size(); i++) {
            discarded_reads.insert(discarded_reads.end(), clusters[i].reads.begin(), clusters[i].reads.end());
        }
        clusters.erase(clusters.begin()+2, clusters.end());

        if (clusters[0].allele_len != clusters[1].allele_len) {
            for (const hp_read_info_t& hp_read_info : discarded_reads) {
                int dist0 = abs(hp_read_info.hp_len - clusters[0].allele_len);
                int dist1 = abs(hp_read_info.hp_len - clusters[1].allele_len);
                int cluster_idx = 0;
                if (dist1 < dist0 || (dist1 == dist0 && clusters[1].allele_len < clusters[0].allele_len)) {
                    cluster_idx = 1;
                }
                clusters[cluster_idx].reads.push_back(hp_read_info);
            }
        } else {
            // This mimicks the behavior of evidence_map.load
            // Ambiguous reads are randomly assigned to one of the two clusters, weighted by the number of unique reads in each cluster
            // Clusters with less than 3 unique reads are treated as having 0 weight, and all ambiguous reads are assigned to the other cluster
            // However, if both clusters have less than 3 unique reads, all ambiguous reads are assigned to the cluster (i.e., variant) with the higher EPR
            int n0 = clusters[0].reads.size(), n1 = clusters[1].reads.size();
            if (n0 < 3) n0 = 0;
            if (n1 < 3) n1 = 0;
            int total_n = n0+n1;
            if (total_n == 0) {
                int cluster_idx = hp_indels[clusters[0].allele_idx]->sample_info.epr >= hp_indels[clusters[1].allele_idx]->sample_info.epr ? 0 : 1;
                for (const hp_read_info_t& hp_read_info : discarded_reads) {
                    clusters[cluster_idx].reads.push_back(hp_read_info);
                }
            } else {
                std::mt19937 gen(seed);
                std::uniform_int_distribution<> dis(1, total_n);
                for (const hp_read_info_t& hp_read_info : discarded_reads) {
                    int cluster_idx = dis(gen) <= n0 ? 0 : 1;
                    clusters[cluster_idx].reads.push_back(hp_read_info);
                }
            }
        }
    }

    return clusters;
}



std::vector<int> calculate_aln_scores(const std::vector<hp_read_info_t>& hp_read_infos, char* ref_seq, int ref_len, 
    StripedSmithWaterman::Aligner& aligner) {
    std::vector<int> scores;
    scores.reserve(hp_read_infos.size());

    StripedSmithWaterman::Alignment aln;
    StripedSmithWaterman::Filter filter_score_only(false, false, 0, 32767);
    for (const hp_read_info_t& hp_read_info : hp_read_infos) {
        if (hp_read_info.read.seq.empty()) {
            scores.push_back(0);
            continue;
        }

        aligner.Align(hp_read_info.read.seq.c_str(), ref_seq, ref_len, filter_score_only, &aln, 0);
        scores.push_back(aln.sw_score);
    }

    return scores;
}

int choose_best_allele_idx_for_hp_len(int allele_len, std::vector<hp_read_info_t>& good_hp_read_infos,
    const std::vector<int>& hp_run_lens, const std::vector<std::unique_ptr<char[]>>& alt_alleles,
    const std::vector<int>& alt_allele_lens, StripedSmithWaterman::Aligner& aligner) {

    std::vector<int> candidate_allele_idxs;
    for (int allele_idx = 0; allele_idx < alt_alleles.size(); allele_idx++) {
        if (hp_run_lens[allele_idx] == allele_len) {
            candidate_allele_idxs.push_back(allele_idx);
        }
    }
    if (candidate_allele_idxs.size() == 1) return candidate_allele_idxs[0];
    if (candidate_allele_idxs.empty()) return hp_run_lens.size() - 1;

    int best_allele_idx = -1;
    int best_score_sum = INT32_MIN;

    for (int allele_idx : candidate_allele_idxs) {
        std::vector<int> scores = calculate_aln_scores(good_hp_read_infos,
            alt_alleles[allele_idx].get(), alt_allele_lens[allele_idx], aligner);
        int score_sum = 0;
        for (int score : scores) score_sum += score;

        if (score_sum > best_score_sum) {
            best_score_sum = score_sum;
            best_allele_idx = allele_idx;
        }
    }

    return best_allele_idx;
}

void genotype_hp_indels_group(std::vector<sv_t*>& hp_indels, hts_pair_pos_t ref_hp_range, open_samFile_t* bam_file, char* contig_seq, hts_pos_t contig_len,
    stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner,
    std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, evidence_logger_t* evidence_logger,
    bool reassign_evidence, evidence_map_t* evidence_map) {

    if (hp_indels.empty()) return;
    for (sv_t* hp_indel : hp_indels) {
        hp_indel->hp_genotyped = true;
        hp_indel->hp_ref_beg = ref_hp_range.beg;
        hp_indel->hp_ref_end = ref_hp_range.end;
    }

    char hp_base = get_homopolymer_base(hp_indels[0], contig_seq);

    std::vector<hp_read_info_t> hp_read_infos, hp_read_infos_assigned_outside_group;

    std::stringstream region_ss;
    region_ss << hp_indels[0]->chr << ":" << ref_hp_range.beg-1 << "-" << ref_hp_range.end+1;
    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region_ss.str().c_str());

    // Build alt alleles
    hts_pos_t extend = stats.read_len - 1;
    hts_pos_t ref_hp_len = ref_hp_range.end - ref_hp_range.beg;
    hts_pos_t alt_start = std::max(hts_pos_t(0), ref_hp_range.beg - extend);
    hts_pos_t alt_end = std::min(contig_len, ref_hp_range.end + extend);
    hts_pos_t left_flank_len = ref_hp_range.beg - alt_start;
    hts_pos_t right_flank_len = alt_end - ref_hp_range.end;
    hts_pos_t ref_allele_len = left_flank_len + ref_hp_len + right_flank_len;

    std::unique_ptr<char[]> ref_allele(new char[ref_allele_len + 1]);
    strncpy(ref_allele.get(), contig_seq + alt_start, left_flank_len);
    memset(ref_allele.get() + left_flank_len, hp_base, ref_hp_len);
    strncpy(ref_allele.get() + left_flank_len + ref_hp_len, contig_seq + ref_hp_range.end, right_flank_len);
    ref_allele[ref_allele_len] = '\0';
    hts_pair_pos_t ref_allele_hp_range = {left_flank_len, left_flank_len + ref_hp_len};

    std::vector<std::unique_ptr<char[]>> candidate_ref_alleles;
    std::vector<int> candidate_ref_allele_lens;
    std::vector<hts_pair_pos_t> candidate_ref_allele_hp_ranges;
    std::vector<std::unique_ptr<char[]>> alt_alleles;
    std::vector<int> alt_allele_lens;
    std::vector<hts_pair_pos_t> alt_allele_hp_ranges;
    for (sv_t* hp_indel : hp_indels) {
        char* lh_seq = generate_haplotype_left(contig_seq, ref_hp_range.beg - 1, extend, hp_indel->aux_indels, hp_indel->aux_snps);
        char* rh_seq = generate_haplotype_right(contig_seq, contig_len, ref_hp_range.end, extend, hp_indel->aux_indels, hp_indel->aux_snps);
        hts_pos_t alt_lh_len = strlen(lh_seq);
        hts_pos_t alt_rh_len = strlen(rh_seq);

        // Exact-REF validation compares haplotypes with identical candidate AUX
        // context, changing only the HP length between REF and ALT.
        hts_pos_t candidate_ref_len = alt_lh_len + ref_hp_len + alt_rh_len;
        std::unique_ptr<char[]> candidate_ref_seq(new char[candidate_ref_len + 1]);
        strncpy(candidate_ref_seq.get(), lh_seq, alt_lh_len);
        memset(candidate_ref_seq.get() + alt_lh_len, hp_base, ref_hp_len);
        strncpy(candidate_ref_seq.get() + alt_lh_len + ref_hp_len, rh_seq, alt_rh_len);
        candidate_ref_seq[candidate_ref_len] = '\0';
        candidate_ref_alleles.push_back(std::move(candidate_ref_seq));
        candidate_ref_allele_lens.push_back(candidate_ref_len);
        candidate_ref_allele_hp_ranges.push_back({alt_lh_len, alt_lh_len + ref_hp_len});

        hts_pos_t alt_hp_len = ref_hp_len + hp_indel->svlen();
        hts_pos_t alt_len = alt_lh_len + alt_hp_len + alt_rh_len;
        std::unique_ptr<char[]> alt_seq(new char[alt_len + 1]);
        strncpy(alt_seq.get(), lh_seq, alt_lh_len);
        memset(alt_seq.get() + alt_lh_len, hp_base, alt_hp_len);
        strncpy(alt_seq.get() + alt_lh_len + alt_hp_len, rh_seq, alt_rh_len);
        alt_seq[alt_len] = '\0';
        alt_alleles.push_back(std::move(alt_seq));
        alt_allele_lens.push_back(alt_len);
        alt_allele_hp_ranges.push_back({alt_lh_len, alt_lh_len + alt_hp_len});
        delete[] lh_seq;
        delete[] rh_seq;
    }

    std::vector<char*> ref_seqs = {ref_allele.get()};
    std::vector<hts_pos_t> ref_lens = {ref_allele_len};
    for (int i = 0; i < hp_indels.size(); i++) {
        std::vector<hts_pos_t> alt_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(
            ref_seqs, ref_lens, alt_alleles[i].get(), alt_allele_lens[i], stats.read_len);
        hp_indels[i]->expected_alt1_reads_frac = (double) alt_ref_diff_reads_expected_positions.size() /
            std::max(hts_pos_t(1), hts_pos_t(alt_allele_lens[i]) - stats.read_len + 1);
    }

    // For each read overlapping the reference HP, calculate its observed HP length,
    // 5'/3' tail lengths and 5'/3' tail mismatch counts.
    StripedSmithWaterman::Alignment alt_aln;
    StripedSmithWaterman::Filter filter_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Filter filter_default;

    bam1_t* read = bam_init1();
    std::vector<std::shared_ptr<bam1_t>> collected_reads;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (!is_proper_pair(read, stats.min_is, stats.max_is)) continue;

        collected_reads.emplace_back(bam_dup1(read), bam_destroy1);
    }
    hts_itr_destroy(iter);

    int left_side_5p_reads = 0, right_side_5p_reads = 0;
    int left_side_5p_indel_reads = 0, right_side_5p_indel_reads = 0;
    for (const std::shared_ptr<bam1_t>& collected_read : collected_reads) {
        hp_adjacent_indel_info_t indel_info = get_adjacent_indel_info(collected_read.get(), ref_hp_range, stats.read_len/2);
        bool is_reverse = bam_is_rev(collected_read.get());
        const hp_side_indel_info_t& five_p_info = is_reverse ? indel_info.right : indel_info.left;
        int& side_5p_reads = is_reverse ? right_side_5p_reads : left_side_5p_reads;
        int& side_5p_indel_reads = is_reverse ? right_side_5p_indel_reads : left_side_5p_indel_reads;
        if (five_p_info.aligned_len >= config.min_clip_len) {
            side_5p_reads++;
            if (five_p_info.indel_len != 0) side_5p_indel_reads++;
        }
    }
    bool has_no_left_indel = five_p_evidence_permits_iterative_hp_len_estimation(left_side_5p_reads, left_side_5p_indel_reads);
    bool has_no_right_indel = five_p_evidence_permits_iterative_hp_len_estimation(right_side_5p_reads, right_side_5p_indel_reads);

    for (const std::shared_ptr<bam1_t>& collected_read : collected_reads) {
        bam1_t* read = collected_read.get();
        bool assigned_outside_group = reassign_evidence;
        for (int i = 0; i < hp_indels.size(); i++) {
             if (!reassign_evidence || !evidence_map->is_read_assigned_to_different_sv(read, hp_indels[i])) {
                assigned_outside_group = false;
                break;
            }
        }

        hp_read_info_t hp_read_info = calculate_hp_read_info(
            read, ref_hp_range, hp_base, contig_seq, contig_len, ref_allele.get(), ref_allele_len,
            ref_allele_hp_range, has_no_left_indel, has_no_right_indel, 10, config.max_seq_error);

        if (hp_read_info.tail_3p_len < config.min_clip_len || hp_read_info.tail_5p_len < config.min_clip_len) {
            // Even if we are discarding this reads because the tails are too short, we still want to prevent it from being used as evidence for non-HP indels
            // when they are aligned better to one of the HP indels
            std::string seq = get_sequence(read);
            for (int i = 0; i < hp_indels.size(); i++) {
                // if the read is assigned to a different SV, no need to align it, just count and continue

                aligner.Align(seq.c_str(), alt_alleles[i].get(), alt_allele_lens[i], filter_score_only, &alt_aln, 0);
                std::vector<std::shared_ptr<bam1_t>> read_v = { std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1) };
                std::vector<int> alt_aln_scores = { alt_aln.sw_score };
                if (evidence_logger) evidence_logger->log_reads_associations(hp_indels[i]->id, 1, read_v, alt_aln_scores);
            }
            continue;
        }
        if (assigned_outside_group) {
            hp_read_infos_assigned_outside_group.push_back(hp_read_info);
        } else {
            hp_read_infos.push_back(hp_read_info);
        }
    }
    collected_reads.clear();

    // Realign reads that "escaped" to a different locus, but their mate betrays them
    StripedSmithWaterman::Alignment ref_aln;
    std::stringstream possible_mates_ss;
    possible_mates_ss << hp_indels[0]->chr << ":" << std::max(hts_pos_t(0), ref_hp_range.beg - stats.max_is)
        << "-" << std::min(contig_len, ref_hp_range.end + stats.max_is);
    iter = sam_itr_querys(bam_file->idx, bam_file->header, possible_mates_ss.str().c_str());
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (!is_dc_pair(read)) continue;

        std::string mate_seq;
        int mate_mapq;
        hts_pos_t endpos = bam_endpos(read);
        if (!bam_is_rev(read) && read->core.pos >= ref_hp_range.beg-stats.max_is && read->core.pos <= ref_hp_range.beg-stats.read_len/2) {
            std::string qname = get_mate_lookup_qname(read);
            if (!mateseqs_w_mapq_chr.count(qname)) continue;
            mate_seq = mateseqs_w_mapq_chr[qname].first;
            mate_mapq = mateseqs_w_mapq_chr[qname].second;
            rc(mate_seq);
        } else if (bam_is_rev(read) && endpos >= ref_hp_range.end+stats.read_len/2 && endpos <= ref_hp_range.end+stats.max_is) {
            std::string qname = get_mate_lookup_qname(read);
            if (!mateseqs_w_mapq_chr.count(qname)) continue;
            mate_seq = mateseqs_w_mapq_chr[qname].first;
            mate_mapq = mateseqs_w_mapq_chr[qname].second;
        } else {
            continue;
        }

        aligner.Align(mate_seq.c_str(), ref_allele.get(), ref_allele_len, filter_default, &ref_aln, 0);

        // reject 5p-clipped reads
        bool aln_as_rev = !bam_is_rev(read); // may seem confusing, but read is the *mate*, not the read we are realigning
        if ((!aln_as_rev && get_left_clip_size(ref_aln) > 0) ||
             (aln_as_rev && get_right_clip_size(ref_aln) > 0)) {
            continue;
        }

        // The rescued read information mostly comes from its mate
        bp_support_read_t rescued_read;
        rescued_read.read_name = bam_get_qname(read);
        rescued_read.mapq = mate_mapq;
        rescued_read.mate_mapq = read->core.qual;
        rescued_read.seq = mate_seq;
        rescued_read.mate_is_reverse = bam_is_rev(read);
        rescued_read.mate_pos = read->core.pos;
        rescued_read.mate_endpos = endpos;
        rescued_read.is_first_in_pair = !is_first_read(read);
        hp_read_info_t hp_read_info = calculate_hp_read_info(ref_aln, mate_seq, ref_allele_hp_range,
            hp_base, ref_allele.get(), ref_allele_len, aln_as_rev, rescued_read, 10,
            has_no_left_indel, has_no_right_indel, config.max_seq_error);

        if (hp_read_info.tail_3p_len < config.min_clip_len || hp_read_info.tail_5p_len < config.min_clip_len) {
            continue;
        }

        // We need this to feed into set_bp_consensus_info
        // However, we don't want to feed it to set_hp_read_mapq_stats
        hp_read_info.rescued = true;
        bool assigned_outside_group = reassign_evidence;
        for (sv_t* hp_indel : hp_indels) {
            if (!reassign_evidence || !evidence_map->is_read_assigned_to_different_sv(hp_read_info.read, hp_indel)) {
                assigned_outside_group = false;
                break;
            }
        }
        if (assigned_outside_group) {
            hp_read_infos_assigned_outside_group.push_back(hp_read_info);
        } else {
            hp_read_infos.push_back(hp_read_info);
        }
    }

    bam_destroy1(read);
    hts_itr_destroy(iter);

    auto remove_invalid_hp_reads = [&](std::vector<hp_read_info_t>& reads) {
        std::vector<bool> valid_reads_mask = get_valid_reads_mask(reads, config.min_clip_len, MAX_TAIL_MISMATCH_RATE);
        std::vector<hp_read_info_t> valid_reads;
        valid_reads.reserve(reads.size());
        for (int i = 0; i < reads.size(); i++) {
            if (valid_reads_mask[i]) valid_reads.push_back(reads[i]);
        }
        reads.swap(valid_reads);
    };
    remove_invalid_hp_reads(hp_read_infos);
    remove_invalid_hp_reads(hp_read_infos_assigned_outside_group);

    std::vector<int> hp_run_lens;
    for (sv_t* hp_indel : hp_indels) {
        hp_run_lens.push_back(ref_hp_range.end - ref_hp_range.beg + hp_indel->svlen());
    }
    hp_run_lens.push_back(ref_hp_range.end - ref_hp_range.beg); // hp_run_lens[hp_indels.size()] corresponds to the reference allele

    // Cluster reads into up to two clusters around the candidate allele HP lengths
    std::vector<hp_allele_cluster_t> clusters = reassign_evidence ?
        cluster_reads_by_nearest_allele_len_reassign(hp_read_infos, hp_run_lens, hp_indels, evidence_map, config.seed) :
        cluster_reads_by_nearest_allele_len(hp_read_infos, alt_alleles, alt_allele_lens, hp_run_lens, aligner);

    std::vector<std::vector<hp_read_info_t>> ref_assigned_hp_read_infos(hp_indels.size());
    std::vector<std::vector<bp_support_read_t>> ref_all_reads(hp_indels.size());
    std::vector<std::vector<bp_support_read_t>> ref_good_reads(hp_indels.size()), ref_good_reads_non_rescued(hp_indels.size());
    std::vector<std::vector<bool>> ref_is_consistent(hp_indels.size());
    std::vector<std::vector<bool>> ref_is_exact_match(hp_indels.size());

    auto candidate_supports_exact_ref = [&](int allele_idx, const hp_read_info_t& hp_read_info) {
        return candidate_haplotypes_support_exact_ref(hp_read_info.read.seq,
            candidate_ref_alleles[allele_idx].get(), candidate_ref_allele_lens[allele_idx],
            candidate_ref_allele_hp_ranges[allele_idx], alt_alleles[allele_idx].get(),
            alt_allele_lens[allele_idx], alt_allele_hp_ranges[allele_idx]);
    };

    auto add_ref_support = [&](int allele_idx, const hp_read_info_t& hp_read_info) {
        ref_assigned_hp_read_infos[allele_idx].push_back(hp_read_info);
        ref_all_reads[allele_idx].push_back(hp_read_info.read);
        bool is_consistent = hp_read_info.is_good_read(config.min_clip_len, MAX_TAIL_MISMATCH_RATE);
        ref_is_consistent[allele_idx].push_back(is_consistent);
        bool is_exact_match = hp_read_info.hp_len == hp_run_lens.back() &&
            !hp_read_info.original_alignment_has_indel_outside_hp && !hp_read_info.hp_deletion_extends_outside_hp &&
            !hp_read_info.hp_insertion_has_non_hp_bases && !hp_read_info.hp_run_extends_into_3p_tail;
        if (is_exact_match && hp_read_info.ref_hp_has_non_hp_read_base) {
            is_exact_match = candidate_supports_exact_ref(allele_idx, hp_read_info);
        }
        ref_is_exact_match[allele_idx].push_back(is_exact_match);
        if (!is_consistent) return;
        ref_good_reads[allele_idx].push_back(hp_read_info.read);
        if (!hp_read_info.rescued) ref_good_reads_non_rescued[allele_idx].push_back(hp_read_info.read);
    };

    std::vector<int> cluster_allele_idxs;
    for (const hp_allele_cluster_t& cluster : clusters) {
        cluster_allele_idxs.push_back(cluster.allele_idx);
    }

    // For each read assigned outside this HP group, check which cluster it would belong to. 
    // If it belongs to the reference allele, incremenet RR and ORR counts for all indels in this group.
    // If it belongs to an indel allele, incrememt OAR counts for all indels in this group.
    for (const hp_read_info_t& hp_read_info : hp_read_infos_assigned_outside_group) {
        int allele_idx = nearest_hp_allele_idx(hp_read_info.hp_len, hp_read_info.read.seq, alt_alleles, alt_allele_lens, hp_run_lens, aligner, &cluster_allele_idxs);
        if (allele_idx == -1) continue;

        bool is_ref_allele = allele_idx == hp_run_lens.size() - 1;
        if (is_ref_allele) {
            for (int i = 0; i < hp_indels.size(); i++) {
                hp_indels[i]->sample_info.orr_bp1_reads++;
                evidence_map->register_orr_support(hp_indels[i]->sample_info, 1, hp_read_info.read);
                add_ref_support(i, hp_read_info);
            }
        } else {
            for (int i = 0; i < hp_indels.size(); i++) {
                hp_indels[i]->sample_info.oar_bp1_reads++;
                evidence_map->register_oar_support(hp_indels[i]->sample_info, 1, hp_read_info.read);
            }
        }
    }
 
    std::vector<std::vector<hp_read_info_t>> alt_assigned_hp_read_infos(hp_indels.size());
    std::vector<std::vector<bp_support_read_t>> alt_all_reads(hp_indels.size());
    std::vector<std::vector<bp_support_read_t>> alt_good_reads(hp_indels.size()),alt_good_reads_non_rescued(hp_indels.size());
    std::vector<std::vector<bool>> alt_is_consistent(hp_indels.size());
    std::vector<std::vector<bool>> alt_is_exact_match(hp_indels.size());

    // Associate each allele-aware cluster to its selected allele
    for (const hp_allele_cluster_t& allele_cluster : clusters) {
        const std::vector<hp_read_info_t>& cluster = allele_cluster.reads;

        // Mark good reads
        std::vector<bp_support_read_t> all_reads, good_reads, good_reads_non_rescued;
        std::vector<bool> is_consistent_read;
        for (const hp_read_info_t& hp_read_info : cluster) {
            all_reads.push_back(hp_read_info.read);
            bool is_consistent = hp_read_info.is_good_read(config.min_clip_len, MAX_TAIL_MISMATCH_RATE);
            is_consistent_read.push_back(is_consistent);
            if (is_consistent) {
                good_reads.push_back(hp_read_info.read);
                if (!hp_read_info.rescued) {
                    good_reads_non_rescued.push_back(hp_read_info.read);
                }
            }
        }

        std::vector<int> best_allele_idxs = {allele_cluster.allele_idx};
        if (!reassign_evidence) {
            best_allele_idxs.clear();
            for (int allele_idx = 0; allele_idx < hp_run_lens.size(); allele_idx++) {
                if (hp_run_lens[allele_idx] == allele_cluster.allele_len) {
                    best_allele_idxs.push_back(allele_idx);
                }
            }
        }

        int ref_allele_idx = hp_run_lens.size() - 1;
        bool is_ref_allele = allele_cluster.allele_idx == ref_allele_idx;

        for (int best_allele_idx : best_allele_idxs) {
            const char* allele_seq = is_ref_allele ? ref_allele.get() : alt_alleles[best_allele_idx].get();
            hts_pos_t allele_len = is_ref_allele ? ref_allele_len : alt_allele_lens[best_allele_idx];
            hts_pair_pos_t allele_hp_range = is_ref_allele ? ref_allele_hp_range : alt_allele_hp_ranges[best_allele_idx];
            bool allele_has_aux_indels = is_ref_allele ? 0 : !hp_indels[best_allele_idx]->aux_indels.empty();

            std::vector<bool> is_exact_match;
            for (const hp_read_info_t& hp_read_info : cluster) {
                bool hp_len_match = hp_read_info.hp_len == hp_run_lens[best_allele_idx];
                bool has_unexplained_outside_hp = has_unexplained_indel_outside_hp(
                    hp_read_info, allele_has_aux_indels, allele_seq, allele_len, allele_hp_range, aligner);
                bool exact_match = hp_len_match && !has_unexplained_outside_hp &&
                    (!is_ref_allele || !hp_read_info.hp_run_extends_into_3p_tail);
                is_exact_match.push_back(exact_match);
            }

            if (is_ref_allele) {
                // Cluster best matches the reference allele
                for (int i = 0; i < hp_indels.size(); i++) {
                    for (const hp_read_info_t& hp_read_info : cluster) {
                        add_ref_support(i, hp_read_info);
                    }
                }
            } else {
                // Cluster best matches an indel allele, so assign its reads to that allele
                alt_assigned_hp_read_infos[best_allele_idx].insert(alt_assigned_hp_read_infos[best_allele_idx].end(), cluster.begin(), cluster.end());
                alt_all_reads[best_allele_idx].insert(alt_all_reads[best_allele_idx].end(), all_reads.begin(), all_reads.end());
                alt_good_reads[best_allele_idx].insert(alt_good_reads[best_allele_idx].end(), good_reads.begin(), good_reads.end());
                alt_good_reads_non_rescued[best_allele_idx].insert(alt_good_reads_non_rescued[best_allele_idx].end(), good_reads_non_rescued.begin(), good_reads_non_rescued.end());
                alt_is_consistent[best_allele_idx].insert(alt_is_consistent[best_allele_idx].end(), is_consistent_read.begin(), is_consistent_read.end());
                alt_is_exact_match[best_allele_idx].insert(alt_is_exact_match[best_allele_idx].end(), is_exact_match.begin(), is_exact_match.end());

                std::vector<int> alt_scores = calculate_aln_scores(cluster, alt_alleles[best_allele_idx].get(), alt_allele_lens[best_allele_idx], aligner);
                if (evidence_logger) evidence_logger->log_reads_associations(hp_indels[best_allele_idx]->id, 1, all_reads, alt_scores);
            }
        }

        // We incremenet OAR counts for indels that are not the best match for this cluster
        if (!is_ref_allele) {
            int supporting_hpid = hp_indels[allele_cluster.allele_idx]->hpid;
            const std::string& supporting_sv_id = hp_indels[allele_cluster.allele_idx]->id;
            for (int target_allele_idx = 0; target_allele_idx < hp_indels.size(); target_allele_idx++) {
                if (std::find(best_allele_idxs.begin(), best_allele_idxs.end(), target_allele_idx) != best_allele_idxs.end()) {
                    continue;
                }
                for (const hp_read_info_t& hp_read_info : cluster) {
                    hp_indels[target_allele_idx]->sample_info.oar_bp1_reads++;
                    evidence_map->register_oar_support(hp_indels[target_allele_idx]->sample_info, 1,
                        hp_read_info.read, supporting_hpid, supporting_sv_id);
                }
            }
        }
    }

    for (int i = 0; i < hp_indels.size(); i++) {
        for (int j = 0; j < alt_all_reads[i].size(); j++) {
            if (!alt_is_consistent[i][j]) continue;
            const bp_support_read_t& read = alt_all_reads[i][j];
            evidence_map->record_assigned_read_consistency(read, read.mate_mapq >= config.high_confidence_mapq, alt_is_exact_match[i][j]);
        }

        set_bp_consensus_info(hp_indels[i]->sample_info.alt_bp1.reads_info, alt_all_reads[i], alt_is_consistent[i], alt_is_exact_match[i], 0.0, 0.0);
        hp_indels[i]->sample_info.alt1_hp_len_mode = find_hp_len_mode(alt_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, false);
        hp_indels[i]->sample_info.alt1_consistent_hp_len_mode = find_hp_len_mode(alt_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, true);
        hp_indels[i]->sample_info.alt1_consistent_hp_len_iqr = find_hp_len_iqr(alt_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, true);
        set_hp_tail_mismatch_rates(alt_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, false, 
            hp_indels[i]->sample_info.alt1_hp_5p_mismatch_rate, hp_indels[i]->sample_info.alt1_hp_3p_mismatch_rate);

        // We shouldn't use the mapping qualities of rescued reads since they reflect the confidence
        // of the original mapping rather than the remapping to the HP alleles
        set_hp_read_mapq_stats(alt_good_reads_non_rescued[i], hp_indels[i]->sample_info.alt1_hp_min_mapq, hp_indels[i]->sample_info.alt1_hp_max_mapq, 
            hp_indels[i]->sample_info.alt1_hp_avg_mapq, hp_indels[i]->sample_info.alt1_hp_stddev_mapq);

        set_bp_consensus_info(hp_indels[i]->sample_info.ref_bp1.reads_info, ref_all_reads[i], ref_is_consistent[i], ref_is_exact_match[i], 0.0, 0.0);
        hp_indels[i]->sample_info.ref1_hp_len_mode = find_hp_len_mode(ref_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, false);
        hp_indels[i]->sample_info.ref1_consistent_hp_len_mode = find_hp_len_mode(ref_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, true);
        hp_indels[i]->sample_info.ref1_consistent_hp_len_iqr = find_hp_len_iqr(ref_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, true);
        set_hp_tail_mismatch_rates(ref_assigned_hp_read_infos[i], config.min_clip_len, MAX_TAIL_MISMATCH_RATE, false,
            hp_indels[i]->sample_info.ref1_hp_5p_mismatch_rate, hp_indels[i]->sample_info.ref1_hp_3p_mismatch_rate);
        set_hp_read_mapq_stats(ref_good_reads_non_rescued[i], hp_indels[i]->sample_info.ref1_hp_min_mapq, hp_indels[i]->sample_info.ref1_hp_max_mapq,
            hp_indels[i]->sample_info.ref1_hp_avg_mapq, hp_indels[i]->sample_info.ref1_hp_stddev_mapq);
    }
}

// hp_indels are guaranteed to be on the same chromosome
void genotype_hp_indels(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<sv_t*> hp_indels,
    stats_t& stats, config_t& config, contig_map_t& contig_map, bam_pool_t* bam_pool,
    std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq_chr, 
    std::vector<double>* global_crossing_isize_dist,
    evidence_logger_t* evidence_logger, bool reassign_evidence, evidence_map_t* evidence_map) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    std::unordered_map<std::string, std::vector<sv_t*>> hp_indels_by_ref_hp_range;
    for (sv_t* hp_indel : hp_indels) {
        hts_pair_pos_t ref_hp_range = find_ref_hp_range_for_indel(hp_indel, contig_seq, contig_len);
        std::string ref_hp_range_key = std::to_string(ref_hp_range.beg) + ":" + std::to_string(ref_hp_range.end);
        hp_indels_by_ref_hp_range[ref_hp_range_key].push_back(hp_indel);
    }

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    for (auto& kv : hp_indels_by_ref_hp_range) {
        std::vector<sv_t*>& hp_indels_in_range = kv.second;
        hts_pair_pos_t ref_hp_range = find_ref_hp_range_for_indel(hp_indels_in_range[0], contig_seq, contig_len);
        genotype_hp_indels_group(hp_indels_in_range, ref_hp_range, bam_file, contig_seq, contig_len, stats, config, aligner,
            *mateseqs_w_mapq_chr, evidence_logger, reassign_evidence, evidence_map);
    }

    depth_filter_indel(contig_name, hp_indels, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, hp_indels, bam_file, config, stats);

    release_mates(contig_id);
}

#endif // GENOTYPE_HP_INDELS_H
