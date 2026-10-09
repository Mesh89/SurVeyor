#include <cstdint>
#include <iostream>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <future>
#include <functional>
#include <memory>
#include <string>
#include <thread>
#include <tuple>
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
#include "consensus.h"
#include "hp_mismatch_rate_thresholds.h"
#include "coarse_coverage.h"
#include "../libs/cptl_stl.h"

std::mutex mtx;
config_t config;
stats_t stats;
std::string workdir, workspace;
chr_seqs_map_t contigs;
hp_tail_quality_model_t hp_tail_quality_model;

std::unordered_map<std::string, int> detected_svs_count;
std::unordered_set<std::string> detected_svs_count_is_hq;
std::unordered_map<std::string, coarse_coverage_track_t> coverage_tracks;

// BAM records stay immutable and alive throughout a window's work. Cache only
// within that task, even when adjacent windows share the same BAM pointers.
struct clip_read_cache_t {
    struct read_data_t {
        std::string seq, dedup_key;
        std::vector<uint8_t> quals;
        bool has_dedup_key = false, has_quals = false;

        explicit read_data_t(bam1_t* read) : seq(get_sequence(read)) {}
    };

    std::unordered_map<const bam1_t*, read_data_t> entries;
    hp_tail_quality_table_t& quality_cache;
    char* contig_seq;
    hts_pos_t contig_len;

    explicit clip_read_cache_t(hp_tail_quality_table_t& quality_cache, char* contig_seq = nullptr, hts_pos_t contig_len = 0) : quality_cache(quality_cache), contig_seq(contig_seq), contig_len(contig_len) {}

    read_data_t& get(bam1_t* read) {
        auto entry = entries.find(read);
        if (entry == entries.end()) entry = entries.emplace(read, read_data_t(read)).first;
        return entry->second;
    }

    const std::string& sequence(bam1_t* read) {
        return get(read).seq;
    }

    const std::string& dedup_key(bam1_t* read) {
        read_data_t& data = get(read);
        if (!data.has_dedup_key) {
            data.dedup_key = std::to_string(read->core.pos) + " " + data.seq + " " + std::to_string(read->core.mpos) + " " + std::string(get_mc(read));
            data.has_dedup_key = true;
        }
        return data.dedup_key;
    }

    std::vector<uint8_t>& qualities(bam1_t* read) {
        read_data_t& data = get(read);
        if (!data.has_quals) {
            data.quals = recalibrate_clip_read_qualities(read, config, hp_tail_quality_model, quality_cache, &data.seq, contig_seq, contig_len);
            data.has_quals = true;
        }
        return data.quals;
    }
};

bool cluster_touches_excessive_coverage(const std::string& contig_name, const std::deque<bam1_t*>& clipped) {
    const auto& tracks = coverage_tracks;
    auto track = tracks.find(contig_name);
    if (track == tracks.end() || clipped.empty()) return false;
    uint64_t threshold = uint64_t(20)*stats.get_max_depth(contig_name);
    for (bam1_t* read : clipped) if (track->second.exceeds(read->core.pos, bam_endpos(read), threshold)) return true;
    return false;
}

struct sync_hts_reader_t {
    std::vector<open_samFile_t*> files;
    std::vector<hts_itr_t*> iters;
    std::vector<hts_pos_t> last_positions;
    std::vector<bool> finished;
    std::vector<uint64_t> read_orders;
    int read_len;
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> scratch_read{bam_init1(), &bam_destroy1};

    struct queued_read_t {
        bam1_t* read;
        size_t stream;
        uint64_t order;
    };
    struct cmp_reads {
        bool operator()(const queued_read_t& r1, const queued_read_t& r2) {
            hts_pos_t start1 = get_unclipped_start(r1.read), start2 = get_unclipped_start(r2.read);
            if (start1 != start2) return start1 > start2;
            // Equal starts must not change order when the amount of lookahead changes.
            if (r1.stream != r2.stream) return r1.stream > r2.stream;
            return r1.order > r2.order;
        }
    };
    std::priority_queue<queued_read_t, std::vector<queued_read_t>, cmp_reads> read_queue;

    sync_hts_reader_t(std::vector<std::string> fnames, std::string region, int read_len) : read_len(read_len) {
        for (std::string fname : fnames) {
            if (!file_exists(fname)) continue;
            open_samFile_t* file = new open_samFile_t(fname);
            files.push_back(file);
            hts_itr_t* iter = sam_itr_querys(file->idx, file->header, region.c_str());
            if (!iter) throw std::runtime_error("Unable to query " + region + " in " + fname);
            iters.push_back(iter);
            last_positions.push_back(-1);
            finished.push_back(false);
            read_orders.push_back(0);
        }
        fill_reads();
    }

    void fill_reads() {
        bam1_t* read = scratch_read.get();
        for (size_t i = 0; i < files.size(); i++) {
            open_samFile_t* file = files[i];
            hts_itr_t* iter = iters[i];
            // Unread BAM positions cannot precede last_positions[i], but their
            // unclipped starts can be up to read_len bases earlier.
            while (!finished[i] && (read_queue.empty() || last_positions[i]-read_len <= get_unclipped_start(read_queue.top().read))) {
                int status = sam_itr_next(file->file, iter, read);
                if (status < -1) {
                    throw std::runtime_error("Failed to read " + std::string(file->file->fn));
                }
                if (status < 0) {
                    finished[i] = true;
                    break;
                }
                last_positions[i] = read->core.pos;
                read_queue.push({bam_dup1(read), i, read_orders[i]++});
            }
        }
    }

    bool next_read(bam1_t*& next_read) {
        if (read_queue.empty()) return false;
        next_read = read_queue.top().read;
        read_queue.pop();
        fill_reads();
        return true;
    }

    ~sync_hts_reader_t() {
        while (!read_queue.empty()) {
            bam_destroy1(read_queue.top().read);
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
        bool consumes_ref = bam_cigar_type(bam_cigar_op(c)) & 2 || op == 'S';
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

std::vector<hts_pos_t> get_read_start_offsets(std::deque<bam1_t*>& reads, clip_read_cache_t& read_cache) {
    std::vector<hts_pos_t> read_start_offsets;
    if (reads.empty()) return read_start_offsets;

    const std::string& r0_seq = read_cache.sequence(reads[0]);
    for (bam1_t* r : reads) {
        const std::string& r_seq = read_cache.sequence(r);
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

// Use kmers to select reads that are likely to be part of the same haplotype
std::vector<int> select_reads_by_kmer(std::vector<std::string>& seqs, const std::vector<uint8_t*>& quals, std::vector<hts_pos_t>& read_start_offsets, size_t selection_rank = 0, std::vector<int>* runner_up_idxs = nullptr) {

    const int K = sizeof(uint32_t)*8/2; // 16-mers
    struct kmer_support_t {
        uint32_t kmer;
        int count = 0;
        std::array<int, K> qual_sums{};
        bool confident = false;

        kmer_support_t(uint32_t kmer) : kmer(kmer) {}
    };

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
    std::vector<std::vector<kmer_support_t>> kmer_counts_by_pos(consensus_len);

    for (int i = 0; i < seqs.size(); i++) {
        std::string& seq = seqs[i];
        if (seq.length() < K) continue;

        uint32_t kmer = 0;
        for (int j = 0; j < seq.length(); j++) {
            kmer = ((kmer << 2) | nucl_bm[seq[j]]);

            if (j >= K-1) {
                std::vector<kmer_support_t>& kmer_counts = kmer_counts_by_pos[read_start_offsets[i]+j];
                size_t group = 0;
                while (group < kmer_counts.size() && kmer_counts[group].kmer != kmer) group++;
                if (group == kmer_counts.size()) kmer_counts.push_back(kmer_support_t(kmer));
                kmer_support_t& support = kmer_counts[group];
                support.count++;
                // Confidence is monotonic; keep counting reads after all bases reach the threshold.
                if (!support.confident) {
                    support.confident = true;
                    for (int base = 0; base < K; base++) {
                        int qual = quals[i][j-K+1+base];
                        support.qual_sums[base] += qual == 255 ? 0 : qual;
                        if (support.qual_sums[base] < 40) support.confident = false;
                    }
                }
            }
        }
    }

    // Prefer two supported groups, then confident bases in both groups, then the existing frequency ranking.
    // select the requested kmer rank at that pos as mandatory kmer
    int chosen_pos = 0;
    uint32_t chosen_kmer = 0;
    uint32_t chosen_freq1 = 0, chosen_freq2 = 0;
    bool chosen_supported = false, chosen_confident = false;
    for (int i = 0; i < kmer_counts_by_pos.size(); i++) {
        // find 1st and 2nd most frequent kmers
        std::sort(kmer_counts_by_pos[i].begin(), kmer_counts_by_pos[i].end(), [](const kmer_support_t& p1, const kmer_support_t& p2) {
            return p1.count > p2.count;
        });

        if (kmer_counts_by_pos[i].empty()) continue;
        const auto& kmers = kmer_counts_by_pos[i];
        int kmer1_freq = kmers[0].count;
        int kmer2_freq = kmers.size() <= 1 ? 1 : kmers[1].count;
        bool supported = kmers.size() > 1 && kmer1_freq >= 3 && kmer2_freq >= 3;
        bool confident = kmers.size() > 1 && kmers[0].confident && kmers[1].confident;
        if (std::make_tuple(supported, confident, kmer1_freq*kmer2_freq, kmer2_freq) >
            std::make_tuple(chosen_supported, chosen_confident, chosen_freq1*chosen_freq2, chosen_freq2)) {
            chosen_pos = i;
            chosen_freq1 = kmer1_freq;
            chosen_freq2 = kmer2_freq;
            chosen_kmer = kmers[0].kmer;
            chosen_supported = supported;
            chosen_confident = confident;
        }
    }

    if (runner_up_idxs) runner_up_idxs->clear();
    if (selection_rank > 0) {
        if (chosen_freq1 < 3) return {};
        const auto& kmers = kmer_counts_by_pos[chosen_pos];
        if (selection_rank >= kmers.size() || kmers[selection_rank].count < 3) return {};
        chosen_kmer = kmers[selection_rank].kmer;
    }

    std::vector<int> selected_idxs;
    if (chosen_freq1 < 3) { // if not enough reads to form a cluster, select all the reads
        for (int i = 0; i < seqs.size(); i++) {
            selected_idxs.push_back(i);
        }
    } else {
        const auto& kmers = kmer_counts_by_pos[chosen_pos];
        bool collect_runner_up = runner_up_idxs && kmers.size() > 1 && kmers[1].count >= 3;
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
            if (collect_runner_up && kmer == kmers[1].kmer) {
                runner_up_idxs->push_back(i);
            }
        }
    }
    return selected_idxs;
}

// Returns an acceptance mask in the same order as reads.
std::vector<bool> find_accepted_reads(std::string& consensus_seq, std::deque<bam1_t*>& reads, std::vector<hts_pos_t>& read_start_offsets, const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds, clip_read_cache_t& read_cache) {

    std::vector<bool> accepted(reads.size(), false);
    std::vector<consensus_hp_region_t> hp_regions = find_consensus_hp_regions(consensus_seq);
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

        const std::string& read_seq = read_cache.sequence(r);
        ungapped_aln_t aln(0, r->core.l_qseq, offset, offset + r->core.l_qseq, mm, r->core.l_qseq - mm);
        if (passes_consensus_mismatch_filter(read_seq, bam_is_rev(r), consensus_seq, aln, hp_regions, hp_mismatch_rate_thresholds, config, true)) {
            // The read should either map much better to the consensus than to the reference,
            // or contain an SV-sized net indel with respect to the reference.
            int orig_score = compute_read_score(r, 1, -4, -6, -1);
            int new_score = fixed_ungapped_aln(read_seq.c_str(), read_seq.length(), consensus_seq.c_str(), consensus_seq.length(), offset, 1, -4, get_left_clip_size(r), get_right_clip_size(r)).score;
            indel_summary_t indels = get_indel_summary(r);
            bool improved_alignment = new_score - orig_score >= config.min_diff_hsr*5; // each mismatch costs 5 points
            bool has_sv_sized_indel = std::abs(indels.dels-indels.inss) >= config.min_sv_size;
            bool indel_explained = has_sv_sized_indel && new_score > orig_score;
            if (improved_alignment || indel_explained) {
                accepted[i] = true;
            }
        }
    }
    return accepted;
}

std::string build_full_consensus_seq(std::deque<bam1_t*>& clipped, bool use_kmer_selection, std::vector<bool>& accepted, int& lowq_prefix, int& lowq_suffix, 
    std::string& consensus_qual, clip_read_cache_t& read_cache,
    const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds = nullptr, size_t selection_rank = 0) {

    std::vector<std::string> seqs;
    std::vector<uint8_t*> quals;
    std::vector<hts_pos_t> read_start_offsets = get_read_start_offsets(clipped, read_cache);

    for (bam1_t* r : clipped) {
        seqs.push_back(read_cache.sequence(r));
        quals.push_back(read_cache.qualities(r).data());
    }

    std::vector<int> selected_idxs, runner_up_idxs;
    if (use_kmer_selection) { // let's try partitioning the sequences according to kmer
        selected_idxs = select_reads_by_kmer(seqs, quals, read_start_offsets, selection_rank, &runner_up_idxs);
    } else {
        selected_idxs.resize(clipped.size());
        for (int i = 0; i < clipped.size(); i++) selected_idxs[i] = i;
    }

    int n_accepted = 0;
    auto build_selected_consensus = [&](const std::vector<int>& idxs, bool require_min_support) {
        n_accepted = 0;
        if (require_min_support && idxs.size() < 3) {
            accepted.assign(clipped.size(), false);
            consensus_qual.clear();
            lowq_prefix = lowq_suffix = 0;
            return std::string();
        }

        std::deque<bam1_t*> selected_clipped;
        std::vector<std::string> selected_seqs;
        std::vector<uint8_t*> selected_quals;
        std::vector<hts_pos_t> selected_read_start_offsets;
        bool select_subset = use_kmer_selection && idxs.size() >= 3;
        if (select_subset) {
            int min_offset = INT32_MAX;
            for (int i : idxs) {
                selected_seqs.push_back(seqs[i]);
                selected_quals.push_back(quals[i]);
                selected_read_start_offsets.push_back(read_start_offsets[i]);
                selected_clipped.push_back(clipped[i]);
                if (min_offset > read_start_offsets[i]) min_offset = read_start_offsets[i];
            }
            for (int i = 0; i < selected_read_start_offsets.size(); i++) {
                selected_read_start_offsets[i] -= min_offset;
            }
        } else if (!use_kmer_selection) {
            selected_clipped = clipped;
        }

        auto& consensus_seqs = select_subset ? selected_seqs : seqs;
        auto& consensus_quals = select_subset ? selected_quals : quals;
        auto& consensus_offsets = select_subset ? selected_read_start_offsets : read_start_offsets;
        std::string consensus_seq = build_full_consensus_seq(consensus_seqs, consensus_quals, consensus_offsets, lowq_prefix, lowq_suffix, consensus_qual, true);

        std::vector<bool> selected_accepted = find_accepted_reads(consensus_seq, selected_clipped, consensus_offsets, hp_mismatch_rate_thresholds, read_cache);

        accepted = std::vector<bool>(clipped.size(), false);
        for (int i = 0; i < selected_accepted.size(); i++) {
            if (selected_accepted[i]) {
                accepted[idxs[i]] = true;
                n_accepted++;
            }
        }
        return consensus_seq;
    };

    std::string consensus_seq = build_selected_consensus(selected_idxs, use_kmer_selection && selection_rank > 0);
    if (use_kmer_selection && selection_rank == 0 && n_accepted < 3) {
        // Retry only the runner-up at the same chosen position before the caller falls back to all reads.
        return build_selected_consensus(runner_up_idxs, true);
    }
    return consensus_seq;
}

void dedup_cluster(std::deque<bam1_t*>& cluster, clip_read_cache_t& read_cache) {
    std::unordered_map<std::string, bam1_t*> seen;
    for (bam1_t* r : cluster) {
        const std::string& key = read_cache.dedup_key(r);
        if (!seen.count(key) || seen[key]->core.qual < r->core.qual) {
            seen[key] = r;
        }
    }

    std::deque<bam1_t*> unique_cluster;
    for (bam1_t* r : cluster) {
        const std::string& key = read_cache.dedup_key(r);
        if (seen[key] == r) {
            unique_cluster.push_back(r);
        }
    }
    cluster.swap(unique_cluster);
}

std::set<int> construction_read_indel_lengths(const std::deque<bam1_t*>& accepted_reads) {
    std::set<int> indel_lengths;
    for (const bam1_t* r : accepted_reads) {
        const uint32_t* cigar = bam_get_cigar(r);
        for (uint32_t i = 0; i < r->core.n_cigar; i++) {
            int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CINS) indel_lengths.insert(len);
            else if (op == BAM_CDEL) indel_lengths.insert(-len);
        }
    }
    return indel_lengths;
}

std::vector<consensus_t*> build_full_consensus(std::string contig_name, std::deque<bam1_t*> clipped, std::deque<bool>& used, clip_read_cache_t& read_cache, const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {

    if (clipped.size() <= 2 || clipped.size() > 20*stats.get_max_depth(contig_name)) {
        return {};
    }
    if (cluster_touches_excessive_coverage(contig_name, clipped)) return {};

    std::deque<bam1_t*> orig_clipped = clipped;

    std::sort(clipped.begin(), clipped.end(), [](bam1_t* r1, bam1_t* r2) {
        return get_unclipped_start(r1) < get_unclipped_start(r2);
    });
    if (get_unclipped_start(clipped[0]) < 0) return {}; // exclude clusters made of reads that are clipped due to hitting the beginning of the chromosome

    dedup_cluster(clipped, read_cache);

    std::unordered_set<bam1_t*> used_reads; // reads used to build a consensus

    std::vector<consensus_t*> consensuses;
    while (clipped.size() >= 3) {
        std::vector<bool> accepted;
        int lowq_prefix, lowq_suffix;
        std::string consensus_qual;
        std::string consensus_seq = build_full_consensus_seq(clipped, true, accepted, lowq_prefix, lowq_suffix, consensus_qual, read_cache, hp_mismatch_rate_thresholds);

        int accepted_reads_n = std::count(accepted.begin(), accepted.end(), true);
        if (accepted_reads_n < 3) {
            consensus_seq = build_full_consensus_seq(clipped, false, accepted, lowq_prefix, lowq_suffix, consensus_qual, read_cache, hp_mismatch_rate_thresholds);
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
            consensus_seq = build_full_consensus_seq(accepted_reads, false, accepted, lowq_prefix, lowq_suffix, consensus_qual, read_cache, hp_mismatch_rate_thresholds);

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

                if (!is_proper_pair(r, stats.min_is, stats.max_is)) continue;
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
            
            int seq_bp_idx = -1, clip_len = 0;
            if (!is_hsr) {
                // Use the same sequence placements as the final accepted-read consensus.
                std::vector<hts_pos_t> read_start_offsets = get_read_start_offsets(accepted_reads, read_cache);
                std::unordered_map<int, int> seq_bp_counts;
                for (size_t i = 0; i < accepted_reads.size(); i++) {
                    bam1_t* r = accepted_reads[i];
                    if (left_clipped && is_left_clipped(r, config.min_clip_len) && r->core.pos == breakpoint) {
                        seq_bp_counts[read_start_offsets[i] + get_left_clip_size(r)]++;
                    } else if (!left_clipped && is_right_clipped(r, config.min_clip_len) && bam_endpos(r) == breakpoint) {
                        seq_bp_counts[read_start_offsets[i] + r->core.l_qseq - get_right_clip_size(r)]++;
                    }
                }
                int max_count = 0;
                for (auto& p : seq_bp_counts) {
                    if (p.second > max_count || (p.second == max_count && (left_clipped ? p.first < seq_bp_idx : p.first > seq_bp_idx))) {
                        max_count = p.second;
                        seq_bp_idx = p.first;
                    }
                }
                clip_len = left_clipped ? seq_bp_idx : (int) consensus_seq.length() - seq_bp_idx;
            }

            consensus_t* consensus = new consensus_t(left_clipped, start, breakpoint, end, consensus_seq,
                consensus_qual, fwd_clipped, rev_clipped, clip_len, max_mapq, lowq_prefix, lowq_suffix);
            consensus->other_bp_lower_boundary = other_bp_lower_boundary;
            consensus->other_bp_upper_boundary = other_bp_upper_boundary;
            consensus->is_hsr = is_hsr;
            consensus->seq_bp_idx = seq_bp_idx;
            consensus->cigar_indel_lengths = construction_read_indel_lengths(accepted_reads);
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

void process_unused_read(bam1_t* read, const std::string& contig_name, clip_read_cache_t* read_cache = nullptr) {
    std::string decoded_seq;
    if (!read_cache) decoded_seq = get_sequence(read);
    const std::string& read_seq = read_cache ? read_cache->sequence(read) : decoded_seq;
    std::vector<std::shared_ptr<sv_t>> svs = detect_svs_from_aln(read, contig_name, read_seq, get_qual_ascii(read), nullptr, 0, 0, stats, config);
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

struct clip_read_t {
    bam1_t* read;
    // Adjacent windows share read-only BAM records and only ever set this flag.
    std::atomic<bool> used_for_consensus{false};
    // The reader holds one use until eviction/EOF, preventing finalization before
    // all windows containing this read have been submitted. Each window adds one.
    std::atomic<size_t> remaining_uses{1};

    explicit clip_read_t(bam1_t* read) : read(read) {}
    ~clip_read_t() { bam_destroy1(read); }

    void finish_use(const std::string& contig_name, clip_read_cache_t* read_cache = nullptr) {
        // The last user sees all consensus-usage flags, regardless of completion order.
        if (remaining_uses.fetch_sub(1) == 1 && !used_for_consensus) process_unused_read(read, contig_name, read_cache);
    }
};

struct consensus_window_t {
    int contig_id;
    std::string contig_name;
    hts_pos_t end;
    std::vector<std::shared_ptr<clip_read_t>> reads;
    size_t initial_cluster_size = 0, retired_reads = 0;
    size_t window_id = 0;
    size_t cluster_read_visits = 0;
    bool last_window = false;
    std::vector<std::unique_ptr<consensus_t>> consensuses;

    consensus_window_t(int contig_id, std::string contig_name, hts_pos_t end) : contig_id(contig_id), contig_name(contig_name), end(end) {}
};

std::shared_ptr<consensus_window_t> build_consensuses(int id, std::shared_ptr<consensus_window_t> window, const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {
    // The sample model is immutable; only calibration-bin lookups outlive this window.
    thread_local hp_tail_quality_table_t quality_cache;
    clip_read_cache_t read_cache(quality_cache, contigs.get_seq(window->contig_name), contigs.get_len(window->contig_name));
    std::deque<bam1_t*> cluster;
    std::deque<bool> used_for_consensus;

    for (size_t i = 0; i < window->reads.size(); i++) {
        bam1_t* read = window->reads[i]->read;
        // The initial cluster is the exact active window left by the preceding task.
        if (i >= window->initial_cluster_size) {
            if (cluster.size() >= 3 && !reads_belong_to_same_cluster(cluster.front(), read)) { // candidate cluster complete
                std::vector<consensus_t*> consensuses = build_full_consensus(window->contig_name, cluster, used_for_consensus, read_cache, hp_mismatch_rate_thresholds);
                for (consensus_t* consensus : consensuses) window->consensuses.emplace_back(consensus);
            }
            while (!cluster.empty() && !reads_belong_to_same_cluster(cluster.front(), read)) {
                if (used_for_consensus.front()) window->reads[window->retired_reads]->used_for_consensus = true;
                window->retired_reads++;
                cluster.pop_front();
                used_for_consensus.pop_front();
            }
        }
        cluster.push_back(read);
        used_for_consensus.push_back(false);
    }

    // Other tasks stop after evicting their last owned cluster; the surviving cluster
    // belongs to the next window and must not be flushed at this artificial boundary.
    if (window->last_window && cluster.size() >= 3) {
        std::vector<consensus_t*> consensuses = build_full_consensus(window->contig_name, cluster, used_for_consensus, read_cache, hp_mismatch_rate_thresholds);
        for (consensus_t* consensus : consensuses) window->consensuses.emplace_back(consensus);
    }
    for (size_t i = 0; i < used_for_consensus.size(); i++) {
        if (used_for_consensus[i]) window->reads[window->retired_reads+i]->used_for_consensus = true;
    }
    for (auto& read : window->reads) read->finish_use(window->contig_name, &read_cache);
    // Completed windows retain only consensus results, not BAM records or read slots.
    std::vector<std::shared_ptr<clip_read_t>>().swap(window->reads);
    return window;
}

void write_consensuses(std::string contig_name, std::string clip_fname, std::vector<consensus_t*>& lc_consensuses, std::vector<consensus_t*>& rc_consensuses) {
    std::ofstream clip_fout(clip_fname);

    filter_fully_contained(rc_consensuses);
    filter_fully_contained(lc_consensuses);

    merge_overlapping_clusters(rc_consensuses, stats.read_len/2);
    merge_overlapping_clusters(lc_consensuses, stats.read_len/2);

    filter_poly_g_tail_consensuses(rc_consensuses, contigs.get_seq(contig_name), contigs.get_len(contig_name), config);
    filter_poly_g_tail_consensuses(lc_consensuses, contigs.get_seq(contig_name), contigs.get_len(contig_name), config);

    enforce_max_ploidy(rc_consensuses, 8);
    enforce_max_ploidy(lc_consensuses, 8);

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
    if (argc != 4) {
        std::cerr << "Usage: clip_consensus_builder <workdir> <reference.fa> <sample>\n";
        return 1;
    }
    workdir = argv[1];
    workspace = workdir + "/workspace";

    std::string reference_fname = argv[2];
    std::string sample_name = argv[3];

    contig_map_t contig_map(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    contigs.read_fasta_into_map(reference_fname, true, config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        coverage_tracks.emplace(contig_name, read_coarse_coverage(workspace + "/coverage/" + std::to_string(contig_id) + ".bin", contigs.get_len(contig_name), config.coverage_window_size));
    }
    hp_tail_quality_model = read_hp_tail_quality_model(workdir + "/hp_3p_tail_error_probabilities.txt");

    hp_mismatch_rate_thresholds_t hp_mismatch_rate_thresholds(workdir + "/" + HP_MISMATCH_RATE_THRESHOLDS_FILENAME);

    const hts_pos_t window_size = 100000;
    const size_t max_cluster_read_visits = 20000;
    const int consensus_workers = config.threads > 1 ? config.threads - 1 : 1;
    const size_t max_pending = 2*size_t(consensus_workers);
    struct contig_windows_t {
        std::vector<std::shared_ptr<consensus_window_t>> windows;
        size_t completed = 0;
        bool all_submitted = false;
    };
    std::vector<contig_windows_t> contig_windows(contig_map.size());
    std::mutex work_mutex;
    std::condition_variable work_cv;
    std::deque<std::shared_ptr<consensus_window_t>> ready_windows, completed_windows;
    std::deque<std::future<void>> postprocessing_futures;
    std::atomic<bool> cancelled{false};
    bool scan_done = false;
    std::exception_ptr failure;
    auto report_failure = [&](std::exception_ptr error) {
        {
            std::lock_guard<std::mutex> lock(work_mutex);
            if (!failure) failure = error;
            cancelled = true;
        }
        work_cv.notify_all();
    };
    size_t pending = 0;
    // Destroy the pool before the queues and failure callback, including on exceptions.
    ctpl::thread_pool thread_pool(consensus_workers);
    auto collect_next = [&]() {
        std::unique_lock<std::mutex> lock(work_mutex);
        work_cv.wait(lock, [&]() { return failure || !completed_windows.empty(); });
        if (failure) std::rethrow_exception(failure);
        std::shared_ptr<consensus_window_t> window = std::move(completed_windows.front());
        completed_windows.pop_front();
        lock.unlock();
        pending--;

        contig_windows_t& results = contig_windows[window->contig_id];
        results.windows[window->window_id] = window;
        results.completed++;
        if (!results.all_submitted || results.completed != results.windows.size()) return;

        if (postprocessing_futures.size() == 2) {
            postprocessing_futures.front().get();
            postprocessing_futures.pop_front();
        }
        // Every window on this contig has finished. Restore the original order for
        // consensus postprocessing; unused reads were handled by their last user.
        std::vector<consensus_t*> lc_consensuses, rc_consensuses;
        for (auto& result : results.windows) {
            for (auto& consensus : result->consensuses) {
                if (consensus->left_clipped) lc_consensuses.push_back(consensus.release());
                else rc_consensuses.push_back(consensus.release());
            }
            result.reset();
        }
        postprocessing_futures.push_back(thread_pool.push([&](int id, std::string contig_name, std::string clip_fname, std::vector<consensus_t*>& lc_consensuses, std::vector<consensus_t*>& rc_consensuses) {
            try {
                write_consensuses(contig_name, clip_fname, lc_consensuses, rc_consensuses);
            } catch (...) {
                report_failure(std::current_exception());
                throw;
            }
        }, window->contig_name, workspace + "/consensuses/" + std::to_string(window->contig_id) + ".txt", std::move(lc_consensuses), std::move(rc_consensuses)));
        results.windows.clear();
    };
    auto submit_window = [&](std::shared_ptr<consensus_window_t> window) {
        contig_windows_t& results = contig_windows[window->contig_id];
        window->window_id = results.windows.size();
        results.windows.push_back(nullptr);
        if (window->last_window) results.all_submitted = true;
        // Completed results have separate storage and do not block new work behind
        // a slow earlier window. Bound only work awaiting completion collection.
        if (pending >= max_pending) collect_next();
        pending++;
        thread_pool.push([&, window](int id) {
            try {
                build_consensuses(id, window, &hp_mismatch_rate_thresholds);
                {
                    std::lock_guard<std::mutex> lock(work_mutex);
                    completed_windows.push_back(window);
                }
                work_cv.notify_all();
            } catch (...) {
                report_failure(std::current_exception());
            }
        });
    };
    auto scan_windows = [&](const std::function<bool(std::shared_ptr<consensus_window_t>)>& emit_window) {
        for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
            if (cancelled) return;
            std::string contig_name = contig_map.get_name(contig_id);
            std::string sr_bam_fname = workspace + "/sr/" + std::to_string(contig_id) + ".bam";
            std::string hsr_bam_fname = workspace + "/hsr/" + std::to_string(contig_id) + ".bam";
            sync_hts_reader_t sync_reader({sr_bam_fname, hsr_bam_fname}, contig_name, stats.read_len);
            std::shared_ptr<consensus_window_t> window = std::make_shared<consensus_window_t>(contig_id, contig_name, window_size);
            std::deque<std::shared_ptr<clip_read_t>> active_reads;
            bam1_t* read = nullptr;
            while (sync_reader.next_read(read)) {
                std::shared_ptr<clip_read_t> owned_read = std::make_shared<clip_read_t>(read);
                if (cancelled) return;
                if (is_left_clipped(read, config.min_clip_len) && is_right_clipped(read, config.min_clip_len)) continue;
                if (!is_clipped(read, config.min_clip_len) && !is_hidden_split_read(read, config)) continue;

                window->reads.push_back(owned_read);
                owned_read->remaining_uses++;
                if (active_reads.size() >= 3 && !reads_belong_to_same_cluster(active_reads.front()->read, read)) {
                    window->cluster_read_visits += active_reads.size();
                }
                while (!active_reads.empty() && !reads_belong_to_same_cluster(active_reads.front()->read, read)) {
                    active_reads.front()->finish_use(contig_name);
                    active_reads.pop_front();
                }
                active_reads.push_back(owned_read);
                hts_pos_t cluster_start = get_unclipped_start(active_reads.front()->read);
                if (cluster_start >= window->end || window->cluster_read_visits >= max_cluster_read_visits) {
                    // Include the triggering read so the worker completes its final cluster,
                    // even beyond either scheduling limit. Carry the surviving cluster intact.
                    if (window->reads.size() > active_reads.size()) {
                        if (!emit_window(window)) return;
                    } else {
                        for (auto& window_read : window->reads) window_read->finish_use(contig_name);
                    }
                    // A new window resets the work counter while retaining the 100 kb grid.
                    window = std::make_shared<consensus_window_t>(contig_id, contig_name, (cluster_start/window_size+1)*window_size);
                    window->reads.assign(active_reads.begin(), active_reads.end());
                    for (auto& active_read : active_reads) active_read->remaining_uses++;
                    window->initial_cluster_size = active_reads.size();
                }
            }
            window->last_window = true;
            if (!emit_window(window)) return;
            for (auto& active_read : active_reads) active_read->finish_use(contig_name);
        }
    };
    std::thread scanner;
    try {
        if (config.threads > 1) {
            scanner = std::thread([&]() {
                try {
                    scan_windows([&](std::shared_ptr<consensus_window_t> window) {
                        std::unique_lock<std::mutex> lock(work_mutex);
                        work_cv.wait(lock, [&]() { return cancelled || ready_windows.size() < max_pending; });
                        if (cancelled) return false;
                        ready_windows.push_back(std::move(window));
                        lock.unlock();
                        work_cv.notify_all();
                        return true;
                    });
                } catch (...) {
                    report_failure(std::current_exception());
                }
                {
                    std::lock_guard<std::mutex> lock(work_mutex);
                    scan_done = true;
                }
                work_cv.notify_all();
            });
            while (true) {
                if (pending >= max_pending) {
                    collect_next();
                    continue;
                }
                std::unique_lock<std::mutex> lock(work_mutex);
                work_cv.wait(lock, [&]() { return failure || scan_done || !ready_windows.empty() || !completed_windows.empty(); });
                if (failure) std::rethrow_exception(failure);
                if (!completed_windows.empty()) {
                    lock.unlock();
                    collect_next();
                    continue;
                }
                if (ready_windows.empty()) break; // Scanning has finished and every window was submitted.
                auto window = std::move(ready_windows.front());
                ready_windows.pop_front();
                lock.unlock();
                work_cv.notify_all();
                submit_window(window);
            }
        } else {
            // Preserve inline scanning and the single consensus worker for low thread counts.
            scan_windows([&](std::shared_ptr<consensus_window_t> window) {
                if (cancelled) return false;
                submit_window(window);
                return true;
            });
        }
        while (pending) collect_next();
    } catch (...) {
        report_failure(std::current_exception());
    }
    if (scanner.joinable()) scanner.join();
    while (!postprocessing_futures.empty()) {
        try {
            postprocessing_futures.front().get();
        } catch (...) {
            report_failure(std::current_exception());
        }
        postprocessing_futures.pop_front();
    }
    thread_pool.stop(true);
    if (failure) std::rethrow_exception(failure);

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
