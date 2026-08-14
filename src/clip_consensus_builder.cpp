#include <cstdint>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "htslib/sam.h"
#include "sam_utils.h"
#include "sw_utils.h"
#include "hsr_utils.h"
#include "vcf_utils.h"
#include "utils.h"
#include "assemble.h"
#include "../libs/cptl_stl.h"

std::mutex mtx;
config_t config;
stats_t stats;
std::string workdir, workspace;
chr_seqs_map_t contigs;

std::unordered_map<std::string, int> detected_svs_count;
std::unordered_set<std::string> detected_svs_count_is_hq;

struct sync_hts_reader_t {
    std::vector<open_samFile_t*> files;
    std::vector<hts_itr_t*> iters;
    int read_len;
    
    struct cmp_reads {
        bool operator()(bam1_t* r1, bam1_t* r2) {
            return get_unclipped_start(r1) > get_unclipped_start(r2);
        }
    };
    std::priority_queue<bam1_t*, std::vector<bam1_t*>, cmp_reads> read_queue;

    sync_hts_reader_t(std::vector<std::string> fnames, std::string region, int read_len) : read_len(read_len) {
        bam1_t* read = bam_init1();
        for (std::string fname : fnames) {
            if (!file_exists(fname)) continue;
            open_samFile_t* file = new open_samFile_t(fname, true);
            hts_itr_t* iter = sam_itr_querys(file->idx, file->header, region.c_str());
            if (sam_itr_next(file->file, iter, read) >= 0) {
                read_queue.push(bam_dup1(read));
            }
            files.push_back(file);
            iters.push_back(iter);
        }
        bam_destroy1(read);
        fill_reads();
    }

    void fill_reads() {
        bam1_t* read = bam_init1();
        for (size_t i = 0; i < files.size(); i++) {
            open_samFile_t* file = files[i];
            hts_itr_t* iter = iters[i];
            while (sam_itr_next(file->file, iter, read) >= 0) {
                read_queue.push(bam_dup1(read));
                if (read->core.pos-read_len > get_unclipped_start(read_queue.top())) break;
            }
        }
        bam_destroy1(read);
    }

    bool next_read(bam1_t*& next_read) {
        if (read_queue.empty()) return false;
        next_read = read_queue.top();
        read_queue.pop();
        fill_reads();
        return true;
    }

    ~sync_hts_reader_t() {
        while (!read_queue.empty()) {
            bam_destroy1(read_queue.top());
            read_queue.pop();
        }
        for (hts_itr_t* iter : iters) {
            sam_itr_destroy(iter);
        }
        for (open_samFile_t* file : files) {
            delete file;
        }
    }
};

std::pair<int, int> get_dels_ins_in_first_n_chars(std::vector<uint32_t>& cigar, int n) {

    if (n < 0) return {0, 0};

    int dels = 0, inss = 0;
    int offset = 0;
    for (uint32_t c : cigar) {
        int len = bam_cigar_oplen(c);
        char op = bam_cigar_opchr(c);

        // since we are "unrolling" soft-clipped bases, they must be accounted for
        bool consumes_ref = bam_cigar_type(c) & 2 || op == 'S';
        if (consumes_ref && offset + len > n) {
            len = n-offset;
        }

        if (op == 'D') {
            dels += len;
        } else if (op == 'I') {
            inss += len;
        }

        if (consumes_ref) {
            offset += len;
            if (offset == n) break;
        }
    }
    return {dels, inss};
}

hts_pos_t get_start_offset(bam1_t* r1, bam1_t* r2) {
    hts_pos_t offset = get_unclipped_start(r2) - get_unclipped_start(r1);
    /* suppose the difference in starting position between R1 and R2 is N, we may be tempted to align the start
     * of R2 to position N+1 of R1
     * However, if R1 has insertions or deletions in the first N bps, the difference in starting positions
     * is no longer accurate to decide where the start of R2 aligns on R1.
     * If R1 has I bps inserted in the first N bps, then R2 will align to position N+1+I.
     * Conversely, if D bps are deleted, R2 will align to position N+1-D
     */
    std::vector<uint32_t> cigar(bam_get_cigar(r1), bam_get_cigar(r1)+r1->core.n_cigar);
    std::pair<int, int> del_ins = get_dels_ins_in_first_n_chars(cigar, offset);
    return offset + del_ins.second - del_ins.first;
}

hts_pos_t get_end_offset(bam1_t* r1, bam1_t* r2) {
    hts_pos_t offset = get_unclipped_end(r2) - get_unclipped_end(r1);
    uint32_t* cigar_array = bam_get_cigar(r2);
    std::vector<uint32_t> rev_cigar;
    for (int i = r2->core.n_cigar-1; i >= 0; i--) {
        rev_cigar.push_back(cigar_array[i]);
    }
    std::pair<int, int> del_ins = get_dels_ins_in_first_n_chars(rev_cigar, offset);
    return offset + del_ins.second - del_ins.first;
}

void make_offsets_nonnegative(std::vector<hts_pos_t>& read_start_offsets) {
    hts_pos_t min_offset = 0;
    for (hts_pos_t offset : read_start_offsets) {
        if (offset < min_offset) min_offset = offset;
    }
    if (min_offset < 0) {
        for (hts_pos_t& offset : read_start_offsets) {
            offset -= min_offset;
        }
    }
}

std::vector<hts_pos_t> get_read_start_offsets(std::deque<bam1_t*>& reads) {
    std::vector<hts_pos_t> read_start_offsets;
    if (reads.empty()) return read_start_offsets;

    std::string r0_seq = get_sequence(reads[0]);
    for (bam1_t* r : reads) {
        std::string r_seq = get_sequence(r);
        int offset1 = get_start_offset(reads[0], r);
        int overlap1_len = std::min(r0_seq.length()-offset1, r_seq.length());
        int offset1_mismatches = number_of_mismatches_fast(r0_seq.c_str()+offset1, r_seq.c_str(), overlap1_len, INT32_MAX);
        int offset1_matches = overlap1_len - offset1_mismatches;
        int offset1_score = offset1_matches - 4*offset1_mismatches; // each mismatch costs 4 points
        
        int end_offset;
        if (get_unclipped_end(r) >= get_unclipped_end(reads[0])) {
            end_offset = get_end_offset(reads[0], r);
        } else {
            // otherwise indels between the ends of the two reads would be ignored
            end_offset = -get_end_offset(r, reads[0]);
        }
        int offset2 = end_offset - (r->core.l_qseq - reads[0]->core.l_qseq);
        int offset2_score = INT32_MIN;
        if (offset2 >= 0 && offset2 < r0_seq.length()) {
            int overlap2_len = std::min(r0_seq.length()-offset2, r_seq.length());
            int offset2_mismatches = number_of_mismatches_fast(r0_seq.c_str()+offset2, r_seq.c_str(), overlap2_len, INT32_MAX);
            int offset2_matches = overlap2_len - offset2_mismatches;
            offset2_score = offset2_matches - 4*offset2_mismatches;
        }
    
        read_start_offsets.push_back(offset1_score >= offset2_score ? offset1 : offset2);
    }
    make_offsets_nonnegative(read_start_offsets);
    return read_start_offsets;
}

// mismatch_score, gap_open_score and gap_extend_score must be negative
int compute_read_score(bam1_t* r, int match_score, int mismatch_score, int gap_open_score, int gap_extend_score) {
	int m_bases = 0;
	int eq_bases = 0;
	int x_bases = 0;
	int indel_bases = 0;
	int score = 0;
	uint32_t* cigar = bam_get_cigar(r);
	for (int i = 0; i < r->core.n_cigar; i++) {
		char op = bam_cigar_opchr(cigar[i]);
		int len = bam_cigar_oplen(cigar[i]);
		if (op == 'M') m_bases += len;
		else if (op == '=') eq_bases += len;
		else if (op == 'X') x_bases += len;
		else if (op == 'I' || op == 'D') {
			score += gap_open_score + (len-1)*gap_extend_score;
			indel_bases += len;
		}
	}

	int nm = get_nm(r);
	int mismatches_in_m = std::max(0, nm - indel_bases - x_bases);
	int matches = eq_bases + std::max(0, m_bases - mismatches_in_m);
	int mismatches = x_bases + mismatches_in_m;
	score += match_score*matches + mismatch_score*mismatches;
	return score;
}

struct base_score_t {
    int freq = 0, qual = 0;
    char base;

    base_score_t(char base) : base(base) {}
};
bool operator < (const base_score_t& bs1, const base_score_t& bs2) {
    if (bs1.freq != bs2.freq) return bs1.freq < bs2.freq;
    return bs1.qual < bs2.qual;
}

std::string build_full_consensus_seq(std::vector<std::string>& seqs, std::vector<uint8_t*>& quals, 
    std::vector<hts_pos_t> read_start_offsets, int& lowq_prefix, int& lowq_suffix) {

    // if not already sorted, sort by start offset
    if (!std::is_sorted(read_start_offsets.begin(), read_start_offsets.end())) {
        std::vector<int> order(seqs.size());
        for (int i = 0; i < order.size(); i++) order[i] = i;
        std::sort(order.begin(), order.end(), [&read_start_offsets](int i1, int i2) {
            return read_start_offsets[i1] < read_start_offsets[i2];
        });
        std::vector<std::string> sorted_seqs;
        std::vector<uint8_t*> sorted_quals;
        std::vector<hts_pos_t> sorted_read_start_offsets;
        for (int i : order) {
            sorted_seqs.push_back(seqs[i]);
            sorted_quals.push_back(quals[i]);
            sorted_read_start_offsets.push_back(read_start_offsets[i]);
        }
        seqs = sorted_seqs;
        quals = sorted_quals;
        read_start_offsets = sorted_read_start_offsets;
    }

    hts_pos_t consensus_len = 0;
    for (int i = 0; i < read_start_offsets.size(); i++) {
        if (consensus_len < read_start_offsets[i] + seqs[i].length()) {
            consensus_len = read_start_offsets[i] + seqs[i].length();
        }
    }
    std::string consensus(consensus_len, 'N');

    std::vector<hts_pos_t> read_end_offsets;
    for (int i = 0; i < seqs.size(); i++) {
        hts_pos_t start_offset = read_start_offsets[i];
        hts_pos_t end_offset = start_offset + seqs[i].length() - 1;
        read_end_offsets.push_back(end_offset);
    }

    int s = 0;
    lowq_prefix = 0, lowq_suffix = 0;
    bool low_prefix_done = false;
    for (int i = 0; i < consensus_len; i++) {
        while (s < seqs.size() && read_end_offsets[s] < i) s++;

        base_score_t base_scores[5] = { base_score_t('A'), base_score_t('C'), base_score_t('G'), base_score_t('T'), base_score_t('N') };
        for (int j = s; j < seqs.size() && read_start_offsets[j] <= i; j++) {
            if (read_end_offsets[j] < i) continue;

            char nucl = seqs[j][i - read_start_offsets[j]];
            uint8_t qual = quals[j][i - read_start_offsets[j]];
            base_scores[nt_map[(uint8_t)nucl]].freq++;
            base_scores[nt_map[(uint8_t)nucl]].qual += qual;
        }

        base_score_t best_base_score = max(base_scores[0], base_scores[1], base_scores[2], base_scores[3]);
        consensus[i] = best_base_score.base;
        
        // determine length of low-quality prefix and suffix
        // low-quality prefix is the last position from the start where coverage is < 3 AND max base freq is < 2
        if (base_scores[0].freq + base_scores[1].freq + base_scores[2].freq + base_scores[3].freq < 3) {
            // not enough coverage
            if (!low_prefix_done && best_base_score.freq < 2) {
                lowq_prefix = i + 1;
            } else if (lowq_suffix == 0 && best_base_score.freq < 2) {
                lowq_suffix = consensus_len - i;
            }
        } else {
            low_prefix_done = true;
        }
    }
    return consensus;
}

// Use kmers to select reads that are likely to be part of the same haplotype
std::vector<int> select_reads_by_kmer(std::vector<std::string>& seqs, std::vector<hts_pos_t>& read_start_offsets) {

    const int K = sizeof(uint32_t)*8/2; // 16-mers

    uint64_t nucl_bm[256] = { 0 };
	nucl_bm['A'] = nucl_bm['a'] = 0;
	nucl_bm['C'] = nucl_bm['c'] = 1;
	nucl_bm['G'] = nucl_bm['g'] = 2;
	nucl_bm['T'] = nucl_bm['t'] = 3;
	nucl_bm['N'] = 0;

    int consensus_len = 0;
    for (int i = 0; i < read_start_offsets.size(); i++) {
        if (consensus_len < read_start_offsets[i] + seqs[i].length()) {
            consensus_len = read_start_offsets[i] + seqs[i].length();
        }
    }
    std::vector<std::vector<std::pair<uint32_t, int>>> kmer_counts_by_pos(consensus_len);

    for (int i = 0; i < seqs.size(); i++) {
        std::string& seq = seqs[i];
        if (seq.length() < K) continue;

        uint32_t kmer = 0;
        for (int j = 0; j < seq.length(); j++) {
            kmer = ((kmer << 2) | nucl_bm[seq[j]]);

            if (j >= K-1) {
                std::vector<std::pair<uint32_t, int>>& kmer_counts = kmer_counts_by_pos[read_start_offsets[i]+j];
                bool found = false;
                for (int i = 0; i < kmer_counts.size(); i++) {
                    if (kmer_counts[i].first == kmer) {
                        kmer_counts[i].second++;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    kmer_counts.push_back({kmer, 1});
                }
            }
        }
    }

    // select pos that maximizes the product of the frequencies of the 1st and 2nd most frequent kmers
    // select most frequent kmer at that pos as mandatory kmer
    int chosen_pos = 0;
    uint32_t chosen_kmer = 0;
    uint32_t chosen_freq1 = 0, chosen_freq2 = 0;
    for (int i = 0; i < kmer_counts_by_pos.size(); i++) {
        // find 1st and 2nd most frequent kmers
        std::sort(kmer_counts_by_pos[i].begin(), kmer_counts_by_pos[i].end(), [](const std::pair<uint32_t, int>& p1, const std::pair<uint32_t, int>& p2) {
            return p1.second > p2.second;
        });

        if (kmer_counts_by_pos[i].empty()) continue;
        int kmer1_freq = kmer_counts_by_pos[i][0].second;
        int kmer2_freq = kmer_counts_by_pos[i].size() <= 1 ? 1 : kmer_counts_by_pos[i][1].second;
        if (kmer1_freq*kmer2_freq > chosen_freq1*chosen_freq2 || 
            kmer1_freq*kmer2_freq == chosen_freq1*chosen_freq2 && kmer2_freq > chosen_freq2) {
            chosen_pos = i;
            chosen_freq1 = kmer1_freq;
            chosen_freq2 = kmer2_freq;
            chosen_kmer = kmer_counts_by_pos[i][0].first;
        }
    }

    std::vector<int> selected_idxs;
    if (chosen_freq1 < 3) { // if not enough reads to form a cluster, select all the reads
        for (int i = 0; i < seqs.size(); i++) {
            selected_idxs.push_back(i);
        }
    } else {
        for (int i = 0; i < seqs.size(); i++) {
            std::string& seq = seqs[i];
            if (seq.length() < K) continue;

            uint32_t kmer = 0;
            int kmer_start = chosen_pos - K + 1 - read_start_offsets[i], kmer_end = chosen_pos - read_start_offsets[i];
            if (kmer_start < 0 || kmer_end >= seq.length()) continue;
            for (int j = kmer_start; j <= kmer_end; j++) {
                kmer = ((kmer << 2) | nucl_bm[seq[j]]);
            }

            if (kmer == chosen_kmer) {
                selected_idxs.push_back(i);
            }
        }
    }
    return selected_idxs;
}

// vector<int> contains the number of mismatches for each read, negative if the read is rejected, non-negative if accepted
std::vector<int> find_accepted_reads(std::string& consensus_seq, std::deque<bam1_t*>& reads, std::vector<hts_pos_t>& read_start_offsets) {

    std::vector<int> accepted(reads.size(), 0);
    for (int i = 0; i < reads.size(); i++) {
        bam1_t* r = reads[i];
        hts_pos_t offset = read_start_offsets[i];
        int mm = 0;

        // filter reads with too many differences from the consensus_seq
        uint8_t* seq_array = bam_get_seq(r);
        for (int j = 0; j < r->core.l_qseq; j++) {
            if (j + offset >= consensus_seq.length()) {
                throw std::runtime_error("Consensus/read offset invariant violated in find_accepted_reads. Please report this to the developers.");
            }
            if (consensus_seq[j + offset] != get_base(seq_array, j)) {
                mm++;
            }
        }

        // Keep rejected reads strictly negative, including reads with zero
        // mismatches to the consensus. Accepted reads store their mismatch count.
        accepted[i] = -(mm + 1);
        if (mm <= std::ceil(config.max_seq_error * r->core.l_qseq)) {
            // The read should either map much better to the consensus than to the reference,
            // or contain an SV-sized net indel with respect to the reference.
            int orig_score = compute_read_score(r, 1, -4, -6, -1);
            int new_score = (r->core.l_qseq-mm)*1 - mm*4;
            indel_summary_t indels = get_indel_summary(r);
            bool improved_alignment = new_score - orig_score >= config.min_diff_hsr*5; // each mismatch costs 5 points
            bool has_sv_sized_indel = std::abs(indels.dels-indels.inss) >= config.min_sv_size;
            bool indel_explained = has_sv_sized_indel && new_score > orig_score;
            if (improved_alignment || indel_explained) {
                accepted[i] = mm;
            }
        }
    }
    return accepted;
}

std::string build_full_consensus_seq(std::deque<bam1_t*>& clipped, bool use_kmer_selection,
                                     std::vector<bool>& accepted, int& lowq_prefix, int& lowq_suffix) {

    std::vector<std::string> seqs;
    std::vector<uint8_t*> quals;
    std::vector<hts_pos_t> read_start_offsets = get_read_start_offsets(clipped);

    for (bam1_t* r : clipped) {
        seqs.push_back(get_sequence(r));
        quals.push_back(bam_get_qual(r));
    }

    std::deque<bam1_t*> selected_clipped;
    std::vector<int> selected_idxs;
    if (use_kmer_selection) { // let's try partitioning the sequences according to kmer
        selected_idxs = select_reads_by_kmer(seqs, read_start_offsets);
        if (selected_idxs.size() >= 3) {
            std::vector<std::string> selected_seqs;
            std::vector<uint8_t*> selected_quals;
            std::vector<hts_pos_t> selected_read_start_offsets;
            int min_offset = INT32_MAX;
            for (int i : selected_idxs) {
                selected_seqs.push_back(seqs[i]);
                selected_quals.push_back(quals[i]);
                selected_read_start_offsets.push_back(read_start_offsets[i]);
                selected_clipped.push_back(clipped[i]);
                if (min_offset > read_start_offsets[i]) min_offset = read_start_offsets[i];
            }
            seqs = selected_seqs;
            quals = selected_quals;
            read_start_offsets = selected_read_start_offsets;
            for (int i = 0; i < read_start_offsets.size(); i++) {
                read_start_offsets[i] -= min_offset;
            }
        }
    } else {
        selected_clipped = clipped;
        selected_idxs.resize(clipped.size());
        for (int i = 0; i < clipped.size(); i++) selected_idxs[i] = i;
    }

    std::string consensus_seq = build_full_consensus_seq(seqs, quals, read_start_offsets, lowq_prefix, lowq_suffix);

    std::vector<int> selected_accepted = find_accepted_reads(consensus_seq, selected_clipped, read_start_offsets);

    int n_accepted = 0;
    accepted = std::vector<bool>(clipped.size(), false);
    for (int i = 0; i < selected_accepted.size(); i++) {
        if (selected_accepted[i] >= 0) {
            accepted[selected_idxs[i]] = true;
            n_accepted++;
        }
    }

    return consensus_seq;
}

void dedup_cluster(std::deque<bam1_t*>& cluster) {
    std::unordered_map<std::string, bam1_t*> seen;
    for (bam1_t* r : cluster) {
        std::string key = std::to_string(r->core.pos) + " " +  get_sequence(r) + " " + 
        std::to_string(r->core.mpos) + " " + std::string(get_mc(r));
        if (!seen.count(key) || seen[key]->core.qual < r->core.qual) {
            seen[key] = r;
        }
    }

    std::deque<bam1_t*> unique_cluster;
    for (bam1_t* r : cluster) {
        std::string key = std::to_string(r->core.pos) + " " +  get_sequence(r) + " " + 
        std::to_string(r->core.mpos) + " " + std::string(get_mc(r));
        if (seen[key] == r) {
            unique_cluster.push_back(r);
        }
    }
    cluster.swap(unique_cluster);
}

std::vector<consensus_t*> build_full_consensus(std::string contig_name, std::deque<bam1_t*> clipped, std::deque<bool>& used) {

    if (clipped.size() <= 2 || clipped.size() > 20*stats.get_max_depth(contig_name)) {
        return {};
    }

    std::deque<bam1_t*> orig_clipped = clipped;

    std::sort(clipped.begin(), clipped.end(), [](bam1_t* r1, bam1_t* r2) {
        return get_unclipped_start(r1) < get_unclipped_start(r2);
    });
    if (get_unclipped_start(clipped[0]) < 0) return {}; // exclude clusters made of reads that are clipped due to hitting the beginning of the chromosome

    dedup_cluster(clipped);

    std::unordered_set<bam1_t*> used_reads; // reads used to build a consensus

    std::vector<consensus_t*> consensuses;
    while (clipped.size() >= 3) {
        std::vector<bool> accepted;
        int lowq_prefix, lowq_suffix;
        std::string consensus_seq = build_full_consensus_seq(clipped, true, accepted, lowq_prefix, lowq_suffix);

        int accepted_reads_n = std::count(accepted.begin(), accepted.end(), true);
        if (accepted_reads_n < 3) {
            consensus_seq = build_full_consensus_seq(clipped, false, accepted, lowq_prefix, lowq_suffix);
        }
        accepted_reads_n = std::count(accepted.begin(), accepted.end(), true);

        std::deque<bam1_t*> accepted_reads, rejected_reads;
        for (size_t i = 0; i < clipped.size(); i++) {
            if (accepted[i]) accepted_reads.push_back(clipped[i]);
            else rejected_reads.push_back(clipped[i]);
        }

        if (accepted_reads.size() < 3) {
            // there are noisy regions where reads are clipped due to noisy tails, and they "drown" correct HSRs
            // from creating a meaningful consensus. If possible, try removing them
            const auto old_size = clipped.size();
            clipped.erase(std::remove_if(clipped.begin(), clipped.end(), [](bam1_t* r) {
                return is_left_clipped(r, config.min_clip_len) || is_right_clipped(r, config.min_clip_len);
            }), clipped.end());
            if (clipped.size() == old_size) {
                break; // if no clipped read can be removed, give up on this cluster
            } else {
                continue; // skip clipped.swap(rejected_reads) in this case
            }
        }

        if (accepted_reads.size() >= 3) {
            for (bam1_t* r : accepted_reads) used_reads.insert(r);

            // rebuild consensus sequence using only accepted reads
            consensus_seq = build_full_consensus_seq(accepted_reads, false, accepted, lowq_prefix, lowq_suffix);

            hts_pos_t start = get_unclipped_start(accepted_reads[0]), end = 0;
            for (bam1_t* r : accepted_reads) end = std::max(end, get_unclipped_end(r));
            std::string highq_consensus_seq = consensus_seq.substr(lowq_prefix, consensus_seq.length()-lowq_prefix-lowq_suffix);
            consensus_ref_classification_t ref_classification = classify_consensus_against_ref(
                contigs.get_seq(contig_name), contigs.get_len(contig_name), start, end, highq_consensus_seq, config);
            if (ref_classification == consensus_ref_classification_t::ALIGNS_WELL) {
                clipped.swap(rejected_reads);
                continue;
            }
            bool left_clipped = ref_classification == consensus_ref_classification_t::LEFT;

            hts_pos_t breakpoint = left_clipped ? INT32_MAX : 0; // the current HTS_POS_MAX does not compile on some compilers
            hts_pos_t other_bp_lower_boundary = consensus_t::LOWER_BOUNDARY_NON_CALCULATED;
            hts_pos_t other_bp_upper_boundary = consensus_t::UPPER_BOUNDARY_NON_CALCULATED;
            int fwd_clipped = 0, rev_clipped = 0;
            uint8_t max_mapq = 0;
            for (bam1_t* r : accepted_reads) {
                if (bam_is_rev(r)) rev_clipped++;
                else fwd_clipped++;

                max_mapq = std::max(max_mapq, r->core.qual);

                if (!is_samechr(r) || is_samestr(r)) continue;
                if (left_clipped && bam_is_rev(r) && !is_mate_left_clipped(r)) {
                    other_bp_lower_boundary = std::max(other_bp_lower_boundary, r->core.mpos);
                    other_bp_upper_boundary = std::min(other_bp_upper_boundary, r->core.mpos+stats.max_is);
                } else if (!left_clipped && !bam_is_rev(r) && !is_mate_right_clipped(r)) {
                    hts_pos_t mate_endpos = get_mate_endpos(r);
                    other_bp_lower_boundary = std::max(other_bp_lower_boundary, mate_endpos-stats.max_is);
                    other_bp_upper_boundary = std::min(other_bp_upper_boundary, get_mate_endpos(r));
                }
            }

            if (other_bp_lower_boundary >= other_bp_upper_boundary) {
                clipped.swap(rejected_reads);
                continue;
            }

            // Calculate breakpoint as the most supported clipped position among accepted reads
            // if ties, choose smallest for left-clipped, largest for right-clipped
            // if no support, set breakpoint to start (left-clipped) or end (right-clipped)
            std::unordered_map<hts_pos_t, int> start_counts, end_counts;
            for (bam1_t* r : accepted_reads) {
                if (get_left_clip_size(r) >= config.min_clip_len) start_counts[r->core.pos]++;
                if (get_right_clip_size(r) >= config.min_clip_len) end_counts[bam_endpos(r)]++;
            }
            bool is_hsr = true;
            if (left_clipped) {
                // breakpoint is the most common r->start among accepted reads. If ties, choose the smallest
                int max_count = 0;
                for (auto& p : start_counts) {
                    if (p.second > max_count || (p.second == max_count && p.first < breakpoint)) {
                        max_count = p.second;
                        breakpoint = p.first;
                    }
                }
                if (!max_count) {
                    breakpoint = start;
                } else if (max_count >= 3) {
                    is_hsr = false;
                }
            } else {
                // breakpoint is the most common r->end among accepted reads. If ties, choose the largest
                int max_count = 0;
                for (auto& p : end_counts) {
                    if (p.second > max_count || (p.second == max_count && p.first > breakpoint)) {
                        max_count = p.second;
                        breakpoint = p.first;
                    }
                }
                if (!max_count) {
                    breakpoint = end;
                } else if (max_count >= 3) {
                    is_hsr = false;
                }
            }
            
            int clip_len = left_clipped ? breakpoint-start : end-breakpoint;
            if (is_hsr) clip_len = 0;

            consensus_t* consensus = new consensus_t(left_clipped, start, breakpoint, end, consensus_seq,
                fwd_clipped, rev_clipped, clip_len, max_mapq, lowq_prefix, lowq_suffix);
            consensus->other_bp_lower_boundary = other_bp_lower_boundary;
            consensus->other_bp_upper_boundary = other_bp_upper_boundary;
            consensus->is_hsr = is_hsr;
            consensuses.push_back(consensus);
        }
        clipped.swap(rejected_reads);
    }

    for (size_t i = 0; i < used.size(); i++) {
        if (used_reads.count(orig_clipped[i])) {
            used[i] = true;
        }
    }

    return consensuses;
}

void route_consensuses(const std::vector<consensus_t*>& consensuses,
                       std::vector<consensus_t*>& lc_consensuses, std::vector<consensus_t*>& rc_consensuses) {
    for (consensus_t* consensus : consensuses) {
        if (consensus->left_clipped) lc_consensuses.push_back(consensus);
        else rc_consensuses.push_back(consensus);
    }
}

void process_unused_read(bam1_t* read, const std::string& contig_name) {
    std::vector<std::shared_ptr<sv_t>> svs = detect_svs_from_aln(read, contig_name, get_sequence(read), nullptr, 0, 0, stats, config);
    std::lock_guard<std::mutex> lock(mtx);
    for (auto& sv : svs) {
        detected_svs_count[sv->unique_key(false)]++;
        if (read->core.qual >= config.high_confidence_mapq) {
            detected_svs_count_is_hq.insert(sv->unique_key(false));
        }
    }
}

bool reads_belong_to_same_cluster(bam1_t* r1, bam1_t* r2) {
    return overlap(get_unclipped_start(r1), get_unclipped_end(r1), get_unclipped_start(r2), get_unclipped_end(r2)) >= std::min(r1->core.l_qseq, r2->core.l_qseq)/2;
}

void drop_invalid_other_bp_intervals(std::vector<consensus_t*>& consensuses) {
    for (consensus_t*& consensus : consensuses) {
        if (consensus->other_bp_lower_boundary > consensus->other_bp_upper_boundary) {
            delete consensus;
            consensus = nullptr;
        }
    }
    consensuses.erase(std::remove(consensuses.begin(), consensuses.end(), nullptr), consensuses.end());
}

void build_consensuses(int id, std::string contig_name, std::vector<std::string> bam_fnames, std::string clip_fname) {
    
    std::ofstream clip_fout(clip_fname);
    std::deque<bam1_t*> cluster;
    std::deque<bool> used_for_consensus;

    sync_hts_reader_t sync_reader(bam_fnames, contig_name, stats.read_len);
    std::vector<consensus_t*> lc_consensuses, rc_consensuses;

    bam1_t* read = nullptr;
    while (sync_reader.next_read(read)) {

        if (is_left_clipped(read, config.min_clip_len) && is_right_clipped(read, config.min_clip_len)) {
            bam_destroy1(read);
            continue;
        }
        if (!is_clipped(read, config.min_clip_len) && !is_hidden_split_read(read, config)) {
            bam_destroy1(read);
            continue;
        }

        if (cluster.size() >= 3 && !reads_belong_to_same_cluster(cluster.front(), read)) { // candidate cluster complete
            std::vector<consensus_t*> consensuses = build_full_consensus(contig_name, cluster, used_for_consensus);
            route_consensuses(consensuses, lc_consensuses, rc_consensuses);
        }
        while (!cluster.empty() && !reads_belong_to_same_cluster(cluster.front(), read)) {
            if (!used_for_consensus.front()) process_unused_read(cluster.front(), contig_name);
            bam_destroy1(cluster.front());
            cluster.pop_front();
            used_for_consensus.pop_front();
        }
        cluster.push_back(read);
        used_for_consensus.push_back(false);
    }

    if (cluster.size() >= 3) {
        std::vector<consensus_t*> consensuses = build_full_consensus(contig_name, cluster, used_for_consensus);
        route_consensuses(consensuses, lc_consensuses, rc_consensuses);
    }
    for (int i = 0; i < used_for_consensus.size(); i++) {
        if (!used_for_consensus[i]) process_unused_read(cluster[i], contig_name);
    }
    for (bam1_t* r : cluster) bam_destroy1(r);

    filter_fully_contained(rc_consensuses);
    filter_fully_contained(lc_consensuses);

    merge_overlapping_clusters(rc_consensuses, stats.read_len/2);
    merge_overlapping_clusters(lc_consensuses, stats.read_len/2);

    drop_invalid_other_bp_intervals(rc_consensuses);
    drop_invalid_other_bp_intervals(lc_consensuses);

    enforce_max_ploidy(rc_consensuses, 4);
    enforce_max_ploidy(lc_consensuses, 4);

    if (lc_consensuses.empty() && rc_consensuses.empty()) return;

    std::vector<consensus_t*> all_consensuses;
    all_consensuses.insert(all_consensuses.end(), lc_consensuses.begin(), lc_consensuses.end());
    all_consensuses.insert(all_consensuses.end(), rc_consensuses.begin(), rc_consensuses.end());
    
    std::vector<std::string> consensus_strs;
    for (consensus_t* consensus : all_consensuses) {
        consensus_strs.push_back(consensus->to_string());
        delete consensus;
    }
    std::sort(consensus_strs.begin(), consensus_strs.end());
    for (const std::string& consensus_str : consensus_strs) {
        clip_fout << consensus_str << std::endl;
    }
    clip_fout.close();
}

int main(int argc, char* argv[]) {
    workdir = argv[1];
    workspace = workdir + "/workspace";

    std::string reference_fname = argv[2];
    std::string sample_name = argv[3];

    contig_map_t contig_map(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    contigs.read_fasta_into_map(reference_fname);

    std::vector<std::future<void> > futures;
    ctpl::thread_pool thread_pool(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future;

        std::string sr_bam_fname = workspace + "/sr/" + std::to_string(contig_id) + ".bam";
        std::string hsr_bam_fname = workspace + "/hsr/" + std::to_string(contig_id) + ".bam";
        future = thread_pool.push(build_consensuses, contig_name, std::vector<std::string>{sr_bam_fname, hsr_bam_fname}, 
            workspace + "/consensuses/" + std::to_string(contig_id) + ".txt");
        futures.push_back(std::move(future));
    }
    thread_pool.stop(true);
    for (size_t i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

    // Write detected SVs to VCF
    std::unordered_map<std::string, std::vector<std::shared_ptr<sv_t>>> svs_by_chr;
    for (const std::string& sv_str : detected_svs_count_is_hq) {
        if (detected_svs_count[sv_str] >= 2) { // SV detected by at least 2 reads
            std::string chr, insseq, svtype;
            hts_pos_t start, end;
            size_t pos1 = sv_str.find(':');
            size_t pos2 = sv_str.find(':', pos1+1);
            size_t pos3 = sv_str.find(':', pos2+1);
            size_t pos4 = sv_str.find(':', pos3+1);
            chr = sv_str.substr(0, pos1);
            start = std::stol(sv_str.substr(pos1+1, pos2-pos1-1));
            end = std::stol(sv_str.substr(pos2+1, pos3-pos2-1));
            svtype = sv_str.substr(pos3+1, pos4-pos3-1);
            insseq = sv_str.substr(pos4+1);

            if (svtype == "DEL") {
                std::shared_ptr<deletion_t> del = std::make_shared<deletion_t>(chr, start, end, "", nullptr, nullptr, nullptr, nullptr);
                del->source = "READ";
                del->junction_remap_ref_beg = start;
                del->junction_remap_ref_end = end + 1;
                svs_by_chr[chr].push_back(del);
            } else if (svtype == "INS") {
                std::shared_ptr<insertion_t> ins = std::make_shared<insertion_t>(chr, start, end, insseq, nullptr, nullptr, nullptr, nullptr);
                ins->source = "READ";
                ins->junction_remap_ref_beg = start;
                ins->junction_remap_ref_end = end + 1;
                svs_by_chr[chr].push_back(ins);
            }
        }
    }

    std::string full_cmd_fname = workdir + "/cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
	std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);
    
    bcf_hdr_t* sv_vcf_header = generate_vcf_header(contigs, sample_name, config, full_cmd_str);
    std::string sv_vcf_fname = workdir + "/intermediate_results/read_svs.vcf.gz";
    htsFile* sv_vcf_fout = hts_open(sv_vcf_fname.c_str(), "wz");
    if (bcf_hdr_write(sv_vcf_fout, sv_vcf_header) != 0) {
        throw std::runtime_error("Failed to write to " + sv_vcf_fname + ".");
    }

    bcf1_t* sv_vcf_record = bcf_init();
    for (const std::string& contig_name : contigs.ordered_contigs) {
        std::vector<std::shared_ptr<sv_t>>& svs = svs_by_chr[contig_name];
        std::sort(svs.begin(), svs.end(), sv_output_order);
        for (const auto& sv : svs) {
            sv2bcf(sv_vcf_header, sv_vcf_record, sv.get(), contigs.get_seq(sv->chr));
            if (bcf_write(sv_vcf_fout, sv_vcf_header, sv_vcf_record) != 0) {
                throw std::runtime_error("Failed to write to " + sv_vcf_fname + ".");
            }
        }
    }

    bcf_close(sv_vcf_fout);
    bcf_destroy(sv_vcf_record);
    bcf_hdr_destroy(sv_vcf_header);
}
