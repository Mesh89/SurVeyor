#ifndef GENOTYPE_DELS_H
#define GENOTYPE_DELS_H

#include "htslib/sam.h"
#include "types.h"
#include "sam_utils.h"
#include "utils.h"
#include "stat_tests.h"
#include "hp_mismatch_rate_thresholds.h"
#include "../libs/ssw_cpp.h"

#include "genotype.h"


struct del_read_alignment_targets_t {
    char* alt_seq;
    hts_pos_t alt_len, alt_lh_len;
    char* ref_bp1_seq;
    hts_pos_t ref_bp1_len, ref_bp1_pos;
    char* ref_bp2_seq;
    hts_pos_t ref_bp2_len, ref_bp2_pos;
    hts_pos_t alt_start, alt_end, ref_bp1_end, ref_bp2_start;
};

struct del_read_alignment_sequences_t {
    std::unique_ptr<char[]> lh_seq, rh_seq, alt_seq, ref_bp1_seq, ref_bp1_w_aux_seq, ref_bp2_seq, ref_bp2_w_aux_seq;
    hts_pos_t alt_len, alt_lh_len, alt_rh_len;
    hts_pos_t ref_bp1_len, ref_bp1_w_aux_len, ref_bp2_len, ref_bp2_w_aux_len;
    del_read_alignment_targets_t targets;
};

del_read_alignment_sequences_t build_del_read_alignment_sequences(deletion_t* del, char* contig_seq, hts_pos_t contig_len, int del_start, int del_end, stats_t& stats) {
    del_read_alignment_sequences_t result;
    hts_pos_t extend = stats.read_len-1;

    result.lh_seq.reset(generate_haplotype_left(contig_seq, del_start-1, extend, del->aux_indels, del->aux_snps));
    result.rh_seq.reset(generate_haplotype_right(contig_seq, contig_len, del_end, extend, del->aux_indels, del->aux_snps));
    hts_pos_t alt_start = std::max(hts_pos_t(0), del_start-extend);
    hts_pos_t alt_end = std::min(del_end+extend, contig_len);
    result.alt_lh_len = strlen(result.lh_seq.get());
    result.alt_rh_len = strlen(result.rh_seq.get());
    result.alt_len = result.alt_lh_len + del->ins_seq.length() + result.alt_rh_len;
    result.alt_seq.reset(new char[result.alt_len + 1]);
    strncpy(result.alt_seq.get(), result.lh_seq.get(), result.alt_lh_len);
    strncpy(result.alt_seq.get()+result.alt_lh_len, del->ins_seq.c_str(), del->ins_seq.length());
    strncpy(result.alt_seq.get()+result.alt_lh_len+del->ins_seq.length(), result.rh_seq.get(), result.alt_rh_len);
    result.alt_seq[result.alt_len] = 0;
    to_uppercase(result.alt_seq.get());

    hts_pos_t ref_bp1_start = alt_start;
    hts_pos_t ref_bp1_end = std::min(del_start+extend, contig_len);
    hts_pos_t ref_bp1_pos = del_start - ref_bp1_start;
    result.ref_bp1_len = ref_bp1_end - ref_bp1_start;
    result.ref_bp1_seq.reset(new char[result.ref_bp1_len + 1]);
    strncpy(result.ref_bp1_seq.get(), contig_seq+ref_bp1_start, result.ref_bp1_len);
    result.ref_bp1_seq[result.ref_bp1_len] = 0;

    hts_pos_t del_len = del_end - del_start;
    hts_pos_t ref_bp1_w_aux_del_len = std::min(extend, del_len);
    hts_pos_t ref_bp1_w_aux_rh_len = std::min(result.alt_rh_len, extend - ref_bp1_w_aux_del_len);
    result.ref_bp1_w_aux_len = result.alt_lh_len + ref_bp1_w_aux_del_len + ref_bp1_w_aux_rh_len;
    result.ref_bp1_w_aux_seq.reset(new char[result.ref_bp1_w_aux_len + 1]);
    strncpy(result.ref_bp1_w_aux_seq.get(), result.lh_seq.get(), result.alt_lh_len);
    strncpy(result.ref_bp1_w_aux_seq.get() + result.alt_lh_len, contig_seq + del_start, ref_bp1_w_aux_del_len);
    strncpy(result.ref_bp1_w_aux_seq.get() + result.alt_lh_len + ref_bp1_w_aux_del_len, result.rh_seq.get(), ref_bp1_w_aux_rh_len);
    result.ref_bp1_w_aux_seq[result.ref_bp1_w_aux_len] = 0;

    hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), del_end-extend);
    hts_pos_t ref_bp2_end = alt_end;
    hts_pos_t ref_bp2_pos = del_end - ref_bp2_start;
    result.ref_bp2_len = ref_bp2_end - ref_bp2_start;
    result.ref_bp2_seq.reset(new char[result.ref_bp2_len + 1]);
    strncpy(result.ref_bp2_seq.get(), contig_seq+ref_bp2_start, result.ref_bp2_len);
    result.ref_bp2_seq[result.ref_bp2_len] = 0;

    hts_pos_t ref_bp2_w_aux_del_len = std::min(extend, del_len);
    hts_pos_t ref_bp2_w_aux_lh_len = std::min(result.alt_lh_len, extend - ref_bp2_w_aux_del_len);
    result.ref_bp2_w_aux_len = ref_bp2_w_aux_lh_len + ref_bp2_w_aux_del_len + result.alt_rh_len;
    result.ref_bp2_w_aux_seq.reset(new char[result.ref_bp2_w_aux_len + 1]);
    strncpy(result.ref_bp2_w_aux_seq.get(), result.lh_seq.get() + result.alt_lh_len - ref_bp2_w_aux_lh_len, ref_bp2_w_aux_lh_len);
    strncpy(result.ref_bp2_w_aux_seq.get() + ref_bp2_w_aux_lh_len, contig_seq + del_end - ref_bp2_w_aux_del_len, ref_bp2_w_aux_del_len);
    strncpy(result.ref_bp2_w_aux_seq.get() + ref_bp2_w_aux_lh_len + ref_bp2_w_aux_del_len, result.rh_seq.get(), result.alt_rh_len);
    result.ref_bp2_w_aux_seq[result.ref_bp2_w_aux_len] = 0;

    result.targets.alt_seq = result.alt_seq.get();
    result.targets.alt_len = result.alt_len;
    result.targets.alt_lh_len = result.alt_lh_len;
    result.targets.ref_bp1_seq = result.ref_bp1_w_aux_seq.get();
    result.targets.ref_bp1_len = result.ref_bp1_w_aux_len;
    result.targets.ref_bp1_pos = ref_bp1_pos;
    result.targets.ref_bp2_seq = result.ref_bp2_w_aux_seq.get();
    result.targets.ref_bp2_len = result.ref_bp2_w_aux_len;
    result.targets.ref_bp2_pos = ref_bp2_pos;
    result.targets.alt_start = alt_start;
    result.targets.alt_end = alt_end;
    result.targets.ref_bp1_end = ref_bp1_end;
    result.targets.ref_bp2_start = ref_bp2_start;
    return result;
}

struct del_read_evidence_t {
    std::vector<std::shared_ptr<bam1_t>> alt_reads, ref_bp1_reads, ref_bp2_reads;
    std::vector<int> alt_positions, alt_scores;
    std::vector<bool> alt_spans_bp1, alt_spans_bp2;
    std::vector<std::string> er_read_names;
    int er = 0, er_hq = 0;
    bool too_deep = false;
};

del_read_evidence_t collect_cached_del_read_evidence(deletion_t* del, open_samFile_t* bam_file, int del_start, int del_end, const del_read_alignment_targets_t& targets, stats_t& stats, config_t& config, evidence_map_t* evidence_map) {
    del_read_evidence_t result;
    if (del->sample_info.too_deep) {
        result.too_deep = true;
        return result;
    }

    std::stringstream l_region, r_region;
    l_region << del->chr << ":" << targets.alt_start << "-" << targets.ref_bp1_end;
    r_region << del->chr << ":" << targets.ref_bp2_start << "-" << targets.alt_end;
    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
    bam1_t* read = bam_init1();

    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < del_start || del_end < get_unclipped_start(read)) continue;
        if (del_start < get_unclipped_start(read) && get_unclipped_end(read) < del_end) continue;
        if (!is_samechr(read) || is_samestr(read)) continue;
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < del_start-stats.max_is) continue;
            if (read->core.mpos > del_start && (abs(read->core.pos-read->core.mpos) > 5 || !is_left_clipped(read, config.min_clip_len))) continue;
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > del_end+stats.max_is) continue;
            if (mate_endpos < del_end && (abs(mate_endpos-bam_endpos(read)) > 5 || !is_right_clipped(read, config.min_clip_len))) continue;
        }

        int cached_alt_bp1_pos = evidence_map->get_read_alt_pos(read, del, 1);
        int cached_alt_bp2_pos = evidence_map->get_read_alt_pos(read, del, 2);
        if (cached_alt_bp1_pos >= 0 && cached_alt_bp2_pos >= 0 && cached_alt_bp1_pos != cached_alt_bp2_pos) throw std::runtime_error("Different cached ALT alignment positions for deletion " + del->id + " and read " + read_name_with_suffix(read) + ".");
        int cached_alt_pos = cached_alt_bp1_pos >= 0 ? cached_alt_bp1_pos : cached_alt_bp2_pos;
        bool cached_alt_read = cached_alt_pos >= 0;
        bool cached_ref_bp1_read = evidence_map->is_read_ref(read, del, 1);
        bool cached_ref_bp2_read = evidence_map->is_read_ref(read, del, 2);
        bool cached_er_read = evidence_map->is_read_er(read, del);

        if (cached_er_read) {
            result.er++;
            if (read->core.qual >= config.high_confidence_mapq) result.er_hq++;
        } else if (cached_alt_read) {
            result.alt_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            result.alt_scores.push_back(0);
            result.alt_positions.push_back(cached_alt_pos);
            result.alt_spans_bp1.push_back(cached_alt_bp1_pos >= 0);
            result.alt_spans_bp2.push_back(cached_alt_bp2_pos >= 0);
        } else {
            if (cached_ref_bp1_read) result.ref_bp1_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            if (cached_ref_bp2_read) result.ref_bp2_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        }
    }

    free(regions[0]);
    free(regions[1]);
    bam_destroy1(read);
    hts_itr_destroy(iter);
    return result;
}

void write_aligned_del_read_evidence(deletion_t* del, open_samFile_t* bam_file, char* contig_seq, hts_pos_t contig_len, stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t& evidence_logger) {
    int del_start = del->start + 1, del_end = del->end + 1;
    del_read_alignment_sequences_t sequences = build_del_read_alignment_sequences(del, contig_seq, contig_len, del_start, del_end, stats);
    const del_read_alignment_targets_t& targets = sequences.targets;
    std::vector<std::shared_ptr<bam1_t>> alt_bp1_reads, alt_bp2_reads;
    std::vector<std::shared_ptr<bam1_t>> ref_bp1_reads, ref_bp2_reads;
    std::vector<int> alt_bp1_scores, alt_bp2_scores, alt_bp1_positions, alt_bp2_positions;
    std::vector<std::string> er_read_names;
    int evidence_count = 0;

    std::stringstream l_region, r_region;
    l_region << del->chr << ":" << targets.alt_start << "-" << targets.ref_bp1_end;
    r_region << del->chr << ":" << targets.ref_bp2_start << "-" << targets.alt_end;
    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
    bam1_t* read = bam_init1();

    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Filter filter_with_pos_and_cigar(true, true, 0, 32767);
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < del_start || del_end < get_unclipped_start(read)) continue;
        if (del_start < get_unclipped_start(read) && get_unclipped_end(read) < del_end) continue;

        std::string seq;
        if (!is_samechr(read) || is_samestr(read)) continue;
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < del_start-stats.max_is) continue;
            if (read->core.mpos > del_start && (abs(read->core.pos-read->core.mpos) > 5 || !is_left_clipped(read, config.min_clip_len))) continue;
            seq = get_sequence(read, true);
            rc(seq);
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > del_end+stats.max_is) continue;
            if (mate_endpos < del_end && (abs(mate_endpos-bam_endpos(read)) > 5 || !is_right_clipped(read, config.min_clip_len))) continue;
            seq = get_sequence(read, true);
        }

        uint16_t ref_aln_score = 0;
        bool ref_bp1_better = false, ref_bp2_better = false;
        bool ref_is_exact_match = is_perfectly_aligned(read);
        if (ref_is_exact_match) {
            ref_aln_score = read->core.l_qseq;
            if (read->core.pos < del_start && bam_endpos(read) > del_start) ref_bp1_better = true;
            if (read->core.pos < del_end && bam_endpos(read) > del_end) ref_bp2_better = true;
        } else {
            aligner.Align(seq.c_str(), targets.ref_bp1_seq, targets.ref_bp1_len, filter_with_pos, &ref1_aln, 0);
            aligner.Align(seq.c_str(), targets.ref_bp2_seq, targets.ref_bp2_len, filter_with_pos, &ref2_aln, 0);
            ref_aln_score = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln.sw_score : ref2_aln.sw_score;
            ref_bp1_better = ref1_aln.sw_score >= ref2_aln.sw_score && ref1_aln.ref_begin < targets.ref_bp1_pos && ref1_aln.ref_end >= targets.ref_bp1_pos;
            ref_bp2_better = ref2_aln.sw_score >= ref1_aln.sw_score && ref2_aln.ref_begin < targets.ref_bp2_pos && ref2_aln.ref_end >= targets.ref_bp2_pos;
        }

        alt_aln = align_fast(aligner, seq.c_str(), targets.alt_seq, targets.alt_len, filter_with_pos_and_cigar, ref_is_exact_match);
        hts_pos_t alt_right_flank_pos = targets.alt_lh_len + del->ins_seq.length();
        bool alt_spans_bp1 = alt_aln.ref_begin < targets.alt_lh_len && alt_aln.ref_end >= targets.alt_lh_len;
        bool alt_spans_bp2 = alt_aln.ref_begin < alt_right_flank_pos && alt_aln.ref_end >= alt_right_flank_pos;
        bool alt_better = alt_aln.sw_score > ref_aln_score && (alt_spans_bp1 || alt_spans_bp2);
        if (alt_better) {
            std::shared_ptr<bam1_t> alt_read(bam_dup1(read), bam_destroy1);
            if (alt_spans_bp1) {
                alt_bp1_reads.push_back(alt_read);
                alt_bp1_scores.push_back(alt_aln.sw_score);
                alt_bp1_positions.push_back(alt_aln.ref_begin);
            }
            if (alt_spans_bp2) {
                alt_bp2_reads.push_back(alt_read);
                alt_bp2_scores.push_back(alt_aln.sw_score);
                alt_bp2_positions.push_back(alt_aln.ref_begin);
            }
            evidence_count++;
        } else if (ref_aln_score > alt_aln.sw_score) {
            if (ref_bp1_better) {
                ref_bp1_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                evidence_count++;
            }
            if (ref_bp2_better) {
                ref_bp2_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                evidence_count++;
            }
        } else {
            er_read_names.push_back(read_name_with_suffix(read));
            evidence_count++;
        }

        if (evidence_count > 4 * stats.get_max_depth(del->chr)) {
            del->sample_info.too_deep = true;
            break;
        }
    }

    free(regions[0]);
    free(regions[1]);
    bam_destroy1(read);
    hts_itr_destroy(iter);
    if (del->sample_info.too_deep) return;

    evidence_logger.log_reads_associations(del->chr, del->id, 1, alt_bp1_reads, alt_bp1_scores, alt_bp1_positions);
    evidence_logger.log_reads_associations(del->chr, del->id, 2, alt_bp2_reads, alt_bp2_scores, alt_bp2_positions);
    evidence_logger.log_ref_reads_associations(del->chr, del->id, 1, ref_bp1_reads);
    evidence_logger.log_ref_reads_associations(del->chr, del->id, 2, ref_bp2_reads);
    evidence_logger.log_er_reads_associations(del->chr, del->id, er_read_names);
}

inline void genotype_del(deletion_t* del, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_map_t* evidence_map,
                const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {
    int del_start = del->start, del_end = del->end;

    // build alt allele
    /* POS in VCF is the base BEFORE the deletion - i.e., the first deleted base is POS+1.
     * Therefore, we want the ALT allele to *include* base POS
     * (note that POS is 1-based in the VCF file, but htslib kindly returns the 0-based coordinate here).
     * As for the END coordinate, my current understanding (which may change) is that it represents the last base deleted.
     * Therefore, the ALT allele should NOT include base END, i.e. it should start at END+1.
     * Here we shift both coordinates by 1, to make them the base immediately AFTER the breakpoints, which is a bit more intuitive for me. */
    del_start++; del_end++;

    // all ranges will be start-inclusive and end-exclusive, i.e. [a,b)
    del_read_alignment_sequences_t sequences = build_del_read_alignment_sequences(del, contig_seq, contig_len, del_start, del_end, stats);
    char* lh_seq = sequences.lh_seq.get();
    char* rh_seq = sequences.rh_seq.get();
    char* alt_seq = sequences.alt_seq.get();
    char* ref_bp1_seq = sequences.ref_bp1_seq.get();
    char* ref_bp2_seq = sequences.ref_bp2_seq.get();
    hts_pos_t alt_lh_len = sequences.alt_lh_len;
    hts_pos_t alt_rh_len = sequences.alt_rh_len;
    hts_pos_t alt_len = sequences.alt_len;
    std::vector<char*> ref_seqs = {sequences.ref_bp1_w_aux_seq.get(), sequences.ref_bp2_w_aux_seq.get()};
    std::vector<hts_pos_t> ref_lens = {sequences.ref_bp1_w_aux_len, sequences.ref_bp2_w_aux_len};
    std::vector<hts_pos_t> alt_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_seq, alt_len, stats.read_len);
    del->sample_info.expected_alt1_reads_frac = (double) alt_ref_diff_reads_expected_positions.size() / std::max(hts_pos_t(1), alt_len - stats.read_len + 1);
    del->sample_info.max_feasible_alt1_len = get_max_feasible_alt_len(alt_ref_diff_reads_expected_positions, stats.read_len);
    del_read_alignment_targets_t& read_targets = sequences.targets;

    del_read_evidence_t read_evidence = collect_cached_del_read_evidence(del, bam_file, del_start, del_end, read_targets, stats, config, evidence_map);
    std::vector<std::shared_ptr<bam1_t>>& alt_better_reads = read_evidence.alt_reads;
    std::vector<std::shared_ptr<bam1_t>>& ref_bp1_better_seqs = read_evidence.ref_bp1_reads;
    std::vector<std::shared_ptr<bam1_t>>& ref_bp2_better_seqs = read_evidence.ref_bp2_reads;
    std::vector<int>& alt_better_read_positions = read_evidence.alt_positions;
    std::vector<int>& alt_better_read_scores = read_evidence.alt_scores;
    del->sample_info.alt_ref_equal_reads = read_evidence.er;
    del->sample_info.alt_ref_equal_reads_highmq = read_evidence.er_hq;
    if (read_evidence.too_deep) evidence_map->clear_other_read_support_for_too_deep(del->sample_info);


    std::string alt_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bool> alt_is_exact_read, ref_bp1_is_exact_read, ref_bp2_is_exact_read;
    auto alt_is_consistent_read = gen_consensus_and_classify_seqs(alt_seq, alt_better_reads, std::vector<bool>(), alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read, hp_mismatch_rate_thresholds);
    auto ref_bp1_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp1_seq, ref_bp1_better_seqs, std::vector<bool>(), ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score, ref_bp1_is_exact_read, hp_mismatch_rate_thresholds);
    auto ref_bp2_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp2_seq, ref_bp2_better_seqs, std::vector<bool>(), ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score, ref_bp2_is_exact_read, hp_mismatch_rate_thresholds);

    std::vector<int> alt_better_read_positions_consistent = get_consistent_reads_start_positions(alt_is_consistent_read, alt_better_read_positions);
    del->sample_info.alt1_occ_ratio = occ_ratio(alt_better_read_positions_consistent, alt_ref_diff_reads_expected_positions.size());

    for (int i = 0; i < alt_better_reads.size(); i++) {
        std::shared_ptr<bam1_t>& r = alt_better_reads[i];
        evidence_map->record_assigned_alt_read(del, r.get(), alt_is_consistent_read[i], get_mq(r.get()) >= config.high_confidence_mapq, alt_is_exact_read[i]);
    }

    auto score_del_consensus = [&](const std::string& consensus_seq) {
        char* lh_seq = generate_haplotype_left(contig_seq, del_start-1, consensus_seq.length(), del->aux_indels, del->aux_snps);
        hts_pos_t lh_len = strlen(lh_seq);
        char* rh_seq = generate_haplotype_right(contig_seq, contig_len, del_end, consensus_seq.length(), del->aux_indels, del->aux_snps);
        hts_pos_t rh_len = strlen(rh_seq);
        alignment_targets_t targets;
        targets.alt_len = lh_len+del->ins_seq.length()+rh_len;
        targets.alt_seq = new char[targets.alt_len+1];
        strncpy(targets.alt_seq, lh_seq, lh_len);
        strncpy(targets.alt_seq+lh_len, del->ins_seq.c_str(), del->ins_seq.length());
        strncpy(targets.alt_seq+lh_len+del->ins_seq.length(), rh_seq, rh_len);
        targets.alt_seq[targets.alt_len] = 0;
        targets.left_flank_end = lh_len;
        targets.right_flank_start = lh_len+del->ins_seq.length();

        hts_pos_t lh_start = std::max(hts_pos_t(0), del_start-lh_len);
        hts_pos_t rh_end = std::min(del_end+rh_len, contig_len);
        hts_pos_t lbp_end = std::min(del_start+hts_pos_t(consensus_seq.length()), contig_len);
        hts_pos_t rbp_start = std::max(hts_pos_t(0), del_end-hts_pos_t(consensus_seq.length()));
        targets.ref_seqs.push_back(contig_seq+lh_start);
        targets.ref_lens.push_back(lbp_end-lh_start);
        targets.ref_seqs.push_back(contig_seq+rbp_start);
        targets.ref_lens.push_back(rh_end-rbp_start);

        targets.left_independent_ref_seq = concat2(lh_seq, contig_seq+del_start, lh_len, lbp_end-del_start);
        targets.left_independent_ref_len = lh_len+lbp_end-del_start;
        targets.right_independent_ref_seq = concat2(contig_seq+rbp_start, rh_seq, del_end-rbp_start, rh_len);
        targets.right_independent_ref_len = del_end-rbp_start+rh_len;

        consensus_alignment_metrics_t metrics = score_consensus_alignment(consensus_seq, targets, aligner);
        delete[] targets.alt_seq;
        delete[] targets.left_independent_ref_seq;
        delete[] targets.right_independent_ref_seq;
        delete[] lh_seq;
        delete[] rh_seq;
        return metrics;
    };

    if (alt_consensus_seq.length() >= 2*config.min_clip_len) {
        consensus_alignment_metrics_t unextended_metrics = score_del_consensus(alt_consensus_seq);
        del->sample_info.alt_consensus1_metrics = unextended_metrics;
        std::shared_ptr<consensus_t> alt_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, std::max<hts_pos_t>(0, del_start-GENOTYPE_CONSENSUS_EXTENSION), del_start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, del_end, std::min<hts_pos_t>(contig_len, del_end+GENOTYPE_CONSENSUS_EXTENSION), contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        
        del->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        del->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        del->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        del->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;
        
        consensus_alignment_metrics_t extended_metrics = score_del_consensus(alt_consensus_seq);
        del->sample_info.ext_alt_consensus1_metrics = extended_metrics;

        int lf_aln_rlen = extended_metrics.split_ref_lengths[0];
        int rf_aln_rlen = extended_metrics.split_ref_lengths[1];
        del->left_anchor_aln->start = del_start-lf_aln_rlen;
        del->left_anchor_aln->end = del_start;
        del->left_anchor_aln->seq_len = lf_aln_rlen;
        del->right_anchor_aln->start = del_end;
        del->right_anchor_aln->end = del_end+rf_aln_rlen;
        del->right_anchor_aln->seq_len = rf_aln_rlen;
    }

    set_bp_consensus_info(del->sample_info.alt_bp1.reads_info, alt_better_reads, alt_is_consistent_read, alt_is_exact_read, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(del->sample_info.ref_bp1.reads_info, ref_bp1_better_seqs, ref_bp1_is_consistent_read, ref_bp1_is_exact_read, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(del->sample_info.ref_bp2.reads_info, ref_bp2_better_seqs, ref_bp2_is_consistent_read, ref_bp2_is_exact_read, ref_bp2_avg_score, ref_bp2_stddev_score);

}

inline void genotype_dels(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<deletion_t*> dels,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t& stats, config_t& config, contig_map_t& contig_map,
    bam_pool_t* bam_pool, std::string workdir, std::vector<double>* global_crossing_isize_dist, const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

    int contig_id = contig_map.get_id(contig_name);
    auto chromosome_data = acquire_chromosome_data(contig_id);
    auto mateseqs_w_mapq_chr = chromosome_data.first;
    evidence_map_t* evidence_map = chromosome_data.second;

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (deletion_t* del : dels) {
        target_ivals.push_back({std::max<hts_pos_t>(0, del->start-GENOTYPE_CONSENSUS_EXTENSION), std::min<hts_pos_t>(contig_len, del->start+GENOTYPE_CONSENSUS_EXTENSION)});
        target_ivals.push_back({std::max<hts_pos_t>(0, del->end-GENOTYPE_CONSENSUS_EXTENSION), std::min<hts_pos_t>(contig_len, del->end+GENOTYPE_CONSENSUS_EXTENSION)});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_file, candidate_reads_for_extension);

    std::vector<deletion_t*> small_deletions, large_deletions;       
    std::vector<sv_t*> small_svs;  
    for (deletion_t* del : dels) {
        genotype_del(del, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_map, hp_mismatch_rate_thresholds);
        if (-del->svlen() >= stats.max_is) {
            large_deletions.push_back(del);
        } else {
            small_deletions.push_back(del);
            small_svs.push_back(del);
        }
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    depth_filter_del(contig_name, dels, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, small_svs, bam_file, config, stats);
    std::string mates_nms_file = workdir + "/workspace/long-pairs/" + std::to_string(contig_id) + ".txt";
    calculate_ptn_ratio(contig_name, dels, bam_file, config, stats, false, evidence_map, mates_nms_file);
    count_stray_pairs(contig_name, dels, bam_file, config, stats);
    release_chromosome_data(contig_id);
}
#endif // GENOTYPE_DELS_H
