#ifndef GENOTYPE_INSS_H
#define GENOTYPE_INSS_H

#include "stat_tests.h"

#include "genotype.h"

struct ins_alignment_targets_t {
    hts_pos_t ins_start, ins_end, alt_start, alt_end;
    hts_pos_t alt_lf_len, alt_rf_len;
    hts_pos_t alt_bp1_len, alt_bp1_pos, alt_bp2_len, alt_bp2_pos;
    hts_pos_t ref_bp1_start, ref_bp1_end, ref_bp1_len;
    hts_pos_t ref_bp2_start, ref_bp2_end, ref_bp2_len;
    hts_pos_t ref_bp1_len_with_aux, ref_bp1_pos_with_aux, ref_bp2_len_with_aux, ref_bp2_pos_with_aux;
    std::unique_ptr<char[]> alt_bp1_seq, alt_bp2_seq, ref_bp1_seq_with_aux, ref_bp2_seq_with_aux;
};

ins_alignment_targets_t build_ins_alignment_targets(insertion_t* ins, char* contig_seq, hts_pos_t contig_len, stats_t& stats) {
    ins_alignment_targets_t targets;
    hts_pos_t extend = stats.read_len + 20;
    targets.ins_start = ins->start + 1;
    targets.ins_end = ins->end + 1;
    targets.alt_start = std::max(hts_pos_t(0), targets.ins_start-extend);
    targets.alt_end = std::min(targets.ins_end+extend, contig_len);
    hts_pos_t ins_lh_len = std::min(extend, hts_pos_t(ins->ins_seq.length()));
    hts_pos_t ins_rh_len = ins_lh_len;

    std::unique_ptr<char[]> lf_seq(generate_haplotype_left(contig_seq, targets.ins_start-1, extend, ins->aux_indels, ins->aux_snps));
    std::unique_ptr<char[]> rf_seq(generate_haplotype_right(contig_seq, contig_len, targets.ins_end, extend, ins->aux_indels, ins->aux_snps));
    targets.alt_lf_len = strlen(lf_seq.get());
    targets.alt_rf_len = strlen(rf_seq.get());

    if (ins_lh_len < extend) {
        hts_pos_t extra = std::min(extend-ins_lh_len, targets.alt_rf_len);
        targets.alt_bp1_len = targets.alt_lf_len + ins_lh_len + extra;
        targets.alt_bp1_seq.reset(new char[targets.alt_bp1_len+1]);
        strncpy(targets.alt_bp1_seq.get(), lf_seq.get(), targets.alt_lf_len);
        strncpy(targets.alt_bp1_seq.get()+targets.alt_lf_len, ins->ins_seq.c_str(), ins_lh_len);
        strncpy(targets.alt_bp1_seq.get()+targets.alt_lf_len+ins_lh_len, rf_seq.get(), extra);
    } else {
        targets.alt_bp1_len = targets.alt_lf_len + ins_lh_len;
        targets.alt_bp1_seq.reset(new char[targets.alt_bp1_len+1]);
        strncpy(targets.alt_bp1_seq.get(), lf_seq.get(), targets.alt_lf_len);
        strncpy(targets.alt_bp1_seq.get()+targets.alt_lf_len, ins->ins_seq.c_str(), ins_lh_len);
    }
    targets.alt_bp1_seq[targets.alt_bp1_len] = 0;

    if (ins_rh_len < extend) {
        hts_pos_t extra = std::min(extend-ins_rh_len, targets.alt_lf_len);
        targets.alt_bp2_len = extra + ins_rh_len + targets.alt_rf_len;
        targets.alt_bp2_seq.reset(new char[targets.alt_bp2_len+1]);
        strncpy(targets.alt_bp2_seq.get(), lf_seq.get()+(targets.alt_lf_len-extra), extra);
        strncpy(targets.alt_bp2_seq.get()+extra, ins->ins_seq.c_str(), ins_rh_len);
        strncpy(targets.alt_bp2_seq.get()+extra+ins_rh_len, rf_seq.get(), targets.alt_rf_len);
    } else {
        targets.alt_bp2_len = ins_rh_len + targets.alt_rf_len;
        targets.alt_bp2_seq.reset(new char[targets.alt_bp2_len+1]);
        strncpy(targets.alt_bp2_seq.get(), ins->ins_seq.c_str()+(ins->ins_seq.length()-ins_rh_len), ins_rh_len);
        strncpy(targets.alt_bp2_seq.get()+ins_rh_len, rf_seq.get(), targets.alt_rf_len);
    }
    targets.alt_bp2_seq[targets.alt_bp2_len] = 0;
    to_uppercase(targets.alt_bp1_seq.get());
    to_uppercase(targets.alt_bp2_seq.get());

    targets.ref_bp1_start = targets.alt_start;
    targets.ref_bp1_end = std::min(targets.ins_start+extend, contig_len);
    targets.ref_bp1_len = targets.ref_bp1_end - targets.ref_bp1_start;
    targets.ref_bp2_start = std::max(hts_pos_t(0), targets.ins_end-extend);
    targets.ref_bp2_end = targets.alt_end;
    targets.ref_bp2_len = targets.ref_bp2_end - targets.ref_bp2_start;

    hts_pos_t ref_allele_len = targets.ins_end - targets.ins_start;
    hts_pos_t ref_bp1_ref_len = std::min(extend, ref_allele_len);
    hts_pos_t ref_bp1_rh_len = std::min(targets.alt_rf_len, extend-ref_bp1_ref_len);
    targets.ref_bp1_len_with_aux = targets.alt_lf_len + ref_bp1_ref_len + ref_bp1_rh_len;
    targets.ref_bp1_seq_with_aux.reset(new char[targets.ref_bp1_len_with_aux+1]);
    strncpy(targets.ref_bp1_seq_with_aux.get(), lf_seq.get(), targets.alt_lf_len);
    strncpy(targets.ref_bp1_seq_with_aux.get()+targets.alt_lf_len, contig_seq+targets.ins_start, ref_bp1_ref_len);
    strncpy(targets.ref_bp1_seq_with_aux.get()+targets.alt_lf_len+ref_bp1_ref_len, rf_seq.get(), ref_bp1_rh_len);
    targets.ref_bp1_seq_with_aux[targets.ref_bp1_len_with_aux] = 0;
    targets.ref_bp1_pos_with_aux = targets.alt_lf_len;

    hts_pos_t ref_bp2_ref_len = std::min(extend, ref_allele_len);
    hts_pos_t ref_bp2_lh_len = std::min(targets.alt_lf_len, extend-ref_bp2_ref_len);
    targets.ref_bp2_len_with_aux = ref_bp2_lh_len + ref_bp2_ref_len + targets.alt_rf_len;
    targets.ref_bp2_seq_with_aux.reset(new char[targets.ref_bp2_len_with_aux+1]);
    strncpy(targets.ref_bp2_seq_with_aux.get(), lf_seq.get()+targets.alt_lf_len-ref_bp2_lh_len, ref_bp2_lh_len);
    strncpy(targets.ref_bp2_seq_with_aux.get()+ref_bp2_lh_len, contig_seq+targets.ins_end-ref_bp2_ref_len, ref_bp2_ref_len);
    strncpy(targets.ref_bp2_seq_with_aux.get()+ref_bp2_lh_len+ref_bp2_ref_len, rf_seq.get(), targets.alt_rf_len);
    targets.ref_bp2_seq_with_aux[targets.ref_bp2_len_with_aux] = 0;
    targets.ref_bp2_pos_with_aux = ref_bp2_lh_len + ref_bp2_ref_len;
    targets.alt_bp1_pos = targets.alt_lf_len;
    targets.alt_bp2_pos = targets.alt_bp2_len - targets.alt_rf_len;
    return targets;
}

struct ins_read_evidence_t {
    bool alt_bp1 = false, alt_bp2 = false, ref_bp1 = false, ref_bp2 = false, er = false;
    int alt_bp1_score = 0, alt_bp2_score = 0, alt_bp1_pos = -1, alt_bp2_pos = -1;
};

ins_read_evidence_t get_cached_ins_read_evidence(insertion_t* ins, bam1_t* read, evidence_map_t* evidence_map) {
    ins_read_evidence_t result;
    result.alt_bp1_pos = evidence_map->get_read_alt_pos(read, ins, 1);
    result.alt_bp2_pos = evidence_map->get_read_alt_pos(read, ins, 2);
    result.alt_bp1 = result.alt_bp1_pos >= 0;
    result.alt_bp2 = result.alt_bp2_pos >= 0;
    if (!result.alt_bp1 && !result.alt_bp2) {
        result.ref_bp1 = evidence_map->is_read_ref(read, ins, 1);
        result.ref_bp2 = evidence_map->is_read_ref(read, ins, 2);
        if (!result.ref_bp1 && !result.ref_bp2) result.er = evidence_map->is_read_er(read, ins);
    }
    return result;
}

void write_aligned_ins_read_evidence(insertion_t* ins, open_samFile_t* bam_file, char* contig_seq, hts_pos_t contig_len, stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t& evidence_logger) {
    ins_alignment_targets_t targets = build_ins_alignment_targets(ins, contig_seq, contig_len, stats);
    std::stringstream l_region, r_region;
    l_region << ins->chr << ":" << targets.ref_bp1_start << "-" << targets.ref_bp1_end;
    r_region << ins->chr << ":" << targets.ref_bp2_start << "-" << targets.ref_bp2_end;
    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> alt_bp1_reads, alt_bp2_reads, ref_bp1_reads, ref_bp2_reads;
    std::vector<std::string> er_read_names;
    std::vector<int> alt_bp1_positions, alt_bp2_positions, alt_bp1_scores, alt_bp2_scores;
    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Alignment alt1_aln, alt2_aln, ref1_aln, ref2_aln;

    int aln_reads = 0;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (targets.ins_start < get_unclipped_start(read) && get_unclipped_end(read) < targets.ins_end) continue;

        bool only_allow_alt = get_unclipped_end(read) < targets.ins_start || targets.ins_end < get_unclipped_start(read);
        bool can_align_to_alt = can_align_to(read, stats.max_is, read->core.tid, targets.ins_start) || can_align_to(read, stats.max_is, read->core.tid, targets.ins_end);
        if (only_allow_alt && (!can_align_to_alt || !is_hidden_split_read(read, config))) continue;

        aln_reads++;
        if (aln_reads > 8 * stats.get_max_depth(ins->chr)) {
            ins->sample_info.too_deep = true;
            break;
        }

        std::string seq = get_sequence(read);
        bool ref_is_exact_match = is_perfectly_aligned(read);
        alt1_aln = align_fast(aligner, seq.c_str(), targets.alt_bp1_seq.get(), targets.alt_bp1_len, filter_with_pos, ref_is_exact_match);
        alt2_aln = align_fast(aligner, seq.c_str(), targets.alt_bp2_seq.get(), targets.alt_bp2_len, filter_with_pos, ref_is_exact_match);
        bool alt1_covers_bp1 = alt1_aln.ref_begin <= targets.alt_bp1_pos && alt1_aln.ref_end >= targets.alt_bp1_pos;
        bool alt2_covers_bp2 = alt2_aln.ref_begin <= targets.alt_bp2_pos && alt2_aln.ref_end >= targets.alt_bp2_pos;
        bool ref1_covers_bp1, ref2_covers_bp2;
        if (ref_is_exact_match) {
            ref1_aln.Clear();
            ref2_aln.Clear();
            ref1_aln.sw_score = read->core.l_qseq;
            ref2_aln.sw_score = read->core.l_qseq;
            ref1_covers_bp1 = read->core.pos <= targets.ins_start && bam_endpos(read) > targets.ins_start;
            ref2_covers_bp2 = read->core.pos <= targets.ins_end && bam_endpos(read) > targets.ins_end;
        } else {
            aligner.Align(seq.c_str(), targets.ref_bp1_seq_with_aux.get(), targets.ref_bp1_len_with_aux, filter_with_pos, &ref1_aln, 0);
            aligner.Align(seq.c_str(), targets.ref_bp2_seq_with_aux.get(), targets.ref_bp2_len_with_aux, filter_with_pos, &ref2_aln, 0);
            ref1_covers_bp1 = ref1_aln.ref_begin <= targets.ref_bp1_pos_with_aux && ref1_aln.ref_end >= targets.ref_bp1_pos_with_aux;
            ref2_covers_bp2 = ref2_aln.ref_begin <= targets.ref_bp2_pos_with_aux && ref2_aln.ref_end >= targets.ref_bp2_pos_with_aux;
        }
        StripedSmithWaterman::Alignment& alt_aln = alt1_aln.sw_score >= alt2_aln.sw_score ? alt1_aln : alt2_aln;
        StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
        if (alt_aln.sw_score > ref_aln.sw_score) {
            if (alt1_aln.sw_score >= alt2_aln.sw_score && alt1_covers_bp1) {
                alt_bp1_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                alt_bp1_scores.push_back(alt1_aln.sw_score);
                alt_bp1_positions.push_back(alt1_aln.ref_begin);
            }
            if (alt1_aln.sw_score <= alt2_aln.sw_score && alt2_covers_bp2) {
                alt_bp2_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                alt_bp2_scores.push_back(alt2_aln.sw_score);
                alt_bp2_positions.push_back(alt2_aln.ref_begin);
            }
        } else if (!only_allow_alt && alt_aln.sw_score < ref_aln.sw_score && !is_clipped(ref_aln, config.min_clip_len)) {
            if (ref1_aln.sw_score >= ref2_aln.sw_score && ref1_covers_bp1) ref_bp1_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            if (ref1_aln.sw_score <= ref2_aln.sw_score && ref2_covers_bp2) ref_bp2_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        } else if (!only_allow_alt && alt_aln.sw_score == ref_aln.sw_score) {
            er_read_names.push_back(read_name_with_suffix(read));
        }
    }

    free(regions[0]);
    free(regions[1]);
    bam_destroy1(read);
    hts_itr_destroy(iter);
    if (ins->sample_info.too_deep) return;

    evidence_logger.log_reads_associations(ins->chr, ins->id, 1, alt_bp1_reads, alt_bp1_scores, alt_bp1_positions);
    evidence_logger.log_reads_associations(ins->chr, ins->id, 2, alt_bp2_reads, alt_bp2_scores, alt_bp2_positions);
    evidence_logger.log_ref_reads_associations(ins->chr, ins->id, 1, ref_bp1_reads);
    evidence_logger.log_ref_reads_associations(ins->chr, ins->id, 2, ref_bp2_reads);
    evidence_logger.log_er_reads_associations(ins->chr, ins->id, er_read_names);
}

inline void genotype_ins(insertion_t* ins, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_map_t* evidence_map,
                const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {

    ins_alignment_targets_t alignment_targets = build_ins_alignment_targets(ins, contig_seq, contig_len, stats);
    hts_pos_t ins_start = alignment_targets.ins_start, ins_end = alignment_targets.ins_end;
    hts_pos_t alt_lf_len = alignment_targets.alt_lf_len, alt_rf_len = alignment_targets.alt_rf_len;
    hts_pos_t alt_bp1_len = alignment_targets.alt_bp1_len, alt_bp2_len = alignment_targets.alt_bp2_len;
    char* alt_bp1_seq = alignment_targets.alt_bp1_seq.get();
    char* alt_bp2_seq = alignment_targets.alt_bp2_seq.get();
    ins->alt1_hp_len = longest_homopolymer_len(alt_bp1_seq, alt_bp1_len);
    ins->alt2_hp_len = longest_homopolymer_len(alt_bp2_seq, alt_bp2_len);
    hts_pos_t ref_bp1_start = alignment_targets.ref_bp1_start, ref_bp1_end = alignment_targets.ref_bp1_end, ref_bp1_len = alignment_targets.ref_bp1_len;
    hts_pos_t ref_bp2_start = alignment_targets.ref_bp2_start, ref_bp2_end = alignment_targets.ref_bp2_end, ref_bp2_len = alignment_targets.ref_bp2_len;
    char* ref_bp1_w_aux_seq = alignment_targets.ref_bp1_seq_with_aux.get();
    char* ref_bp2_w_aux_seq = alignment_targets.ref_bp2_seq_with_aux.get();
    hts_pos_t ref_bp1_w_aux_len = alignment_targets.ref_bp1_len_with_aux, ref_bp1_w_aux_pos = alignment_targets.ref_bp1_pos_with_aux;
    hts_pos_t ref_bp2_w_aux_len = alignment_targets.ref_bp2_len_with_aux, ref_bp2_w_aux_pos = alignment_targets.ref_bp2_pos_with_aux;
    ins->ref1_hp_len = longest_homopolymer_len(contig_seq+ref_bp1_start, ref_bp1_len);
    ins->ref2_hp_len = longest_homopolymer_len(contig_seq+ref_bp2_start, ref_bp2_len);
    if (ins->sample_info.too_deep) return;

    std::vector<char*> ref_seqs = {ref_bp1_w_aux_seq, ref_bp2_w_aux_seq};
    std::vector<hts_pos_t> ref_lens = {ref_bp1_w_aux_len, ref_bp2_w_aux_len};
    std::vector<hts_pos_t> alt1_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_bp1_seq, alt_bp1_len, stats.read_len);
    std::vector<hts_pos_t> alt2_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_bp2_seq, alt_bp2_len, stats.read_len);
    ins->sample_info.expected_alt1_reads_frac = (double) alt1_ref_diff_reads_expected_positions.size() / std::max<hts_pos_t>(1, alt_bp1_len - stats.read_len + 1);
    ins->sample_info.expected_alt2_reads_frac = (double) alt2_ref_diff_reads_expected_positions.size() / std::max<hts_pos_t>(1, alt_bp2_len - stats.read_len + 1);
    ins->sample_info.max_feasible_alt1_len = get_max_feasible_alt_len(alt1_ref_diff_reads_expected_positions, stats.read_len);
    ins->sample_info.max_feasible_alt2_len = get_max_feasible_alt_len(alt2_ref_diff_reads_expected_positions, stats.read_len);

    std::stringstream l_region, r_region;
    l_region << ins->chr << ":" << ref_bp1_start << "-" << ref_bp1_end;
    r_region << ins->chr << ":" << ref_bp2_start << "-" << ref_bp2_end;

    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> alt_bp1_better_reads, alt_bp2_better_reads, ref_bp1_better_reads, ref_bp2_better_reads;
    std::vector<int> alt_bp1_better_read_positions, alt_bp2_better_read_positions;
    std::vector<int> alt_bp1_better_scores, alt_bp2_better_scores;
    hts_pos_t alt_bp1_pos = alignment_targets.alt_bp1_pos;
    hts_pos_t alt_bp2_pos = alignment_targets.alt_bp2_pos;

    int aln_reads = 0;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (ins_start < get_unclipped_start(read) && get_unclipped_end(read) < ins_end) continue;

        bool only_allow_alt = get_unclipped_end(read) < ins_start || ins_end < get_unclipped_start(read);
        bool can_align_to_alt = can_align_to(read, stats.max_is, read->core.tid, ins_start) || can_align_to(read, stats.max_is, read->core.tid, ins_end);
        if (only_allow_alt && (!can_align_to_alt || !is_hidden_split_read(read, config))) continue;

        aln_reads++;
        if (aln_reads > 8 * stats.get_max_depth(ins->chr)) {
            alt_bp1_better_reads.clear();
            alt_bp1_better_read_positions.clear();
            alt_bp1_better_scores.clear();
            alt_bp2_better_reads.clear();
            alt_bp2_better_read_positions.clear();
            alt_bp2_better_scores.clear();
            ref_bp1_better_reads.clear();
            ref_bp2_better_reads.clear();
            ins->sample_info.alt_ref_equal_reads = 0;
            ins->sample_info.alt_ref_equal_reads_highmq = 0;
            evidence_map->clear_other_read_support_for_too_deep(ins->sample_info);
            break;
        }

        ins_read_evidence_t read_evidence = get_cached_ins_read_evidence(ins, read, evidence_map);

        if (read_evidence.alt_bp1) {
            alt_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_bp1_better_scores.push_back(read_evidence.alt_bp1_score);
            alt_bp1_better_read_positions.push_back(read_evidence.alt_bp1_pos);
        }
        if (read_evidence.alt_bp2) {
            alt_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_bp2_better_scores.push_back(read_evidence.alt_bp2_score);
            alt_bp2_better_read_positions.push_back(read_evidence.alt_bp2_pos);
        }
        if (read_evidence.ref_bp1) ref_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        if (read_evidence.ref_bp2) ref_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        if (read_evidence.er) {
            ins->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) ins->sample_info.alt_ref_equal_reads_highmq++;
        }
    }
    std::string alt_bp1_consensus_seq, alt_bp2_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_bp1_avg_score, alt_bp2_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_bp1_stddev_score, alt_bp2_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bool> alt_bp1_is_exact_read, alt_bp2_is_exact_read, ref_bp1_is_exact_read, ref_bp2_is_exact_read;
    auto alt_bp1_is_consistent_read = gen_consensus_and_classify_seqs(alt_bp1_seq, alt_bp1_better_reads, std::vector<bool>(), alt_bp1_consensus_seq, alt_bp1_avg_score, alt_bp1_stddev_score, alt_bp1_is_exact_read, hp_mismatch_rate_thresholds);
    auto alt_bp2_is_consistent_read = gen_consensus_and_classify_seqs(alt_bp2_seq, alt_bp2_better_reads, std::vector<bool>(), alt_bp2_consensus_seq, alt_bp2_avg_score, alt_bp2_stddev_score, alt_bp2_is_exact_read, hp_mismatch_rate_thresholds);

    std::vector<int> alt_bp1_better_read_positions_consistent = get_consistent_reads_start_positions(alt_bp1_is_consistent_read, alt_bp1_better_read_positions);
    ins->sample_info.alt1_occ_ratio = occ_ratio(alt_bp1_better_read_positions_consistent, alt1_ref_diff_reads_expected_positions.size());

    std::vector<int> alt_bp2_better_read_positions_consistent = get_consistent_reads_start_positions(alt_bp2_is_consistent_read, alt_bp2_better_read_positions);
    ins->sample_info.alt2_occ_ratio = occ_ratio(alt_bp2_better_read_positions_consistent, alt2_ref_diff_reads_expected_positions.size());

    for (int i = 0; i < alt_bp1_better_reads.size(); i++) {
        const auto& r = alt_bp1_better_reads[i];
        evidence_map->record_assigned_alt_read(ins, r.get(), alt_bp1_is_consistent_read[i], get_mq(r.get()) >= config.high_confidence_mapq, alt_bp1_is_exact_read[i]);
    }
    for (int i = 0; i < alt_bp2_better_reads.size(); i++) {
        const auto& r = alt_bp2_better_reads[i];
        evidence_map->record_assigned_alt_read(ins, r.get(), alt_bp2_is_consistent_read[i], get_mq(r.get()) >= config.high_confidence_mapq, alt_bp2_is_exact_read[i]);
    }

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    auto ref_bp1_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp1_seq, ref_bp1_better_reads, std::vector<bool>(), ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score, ref_bp1_is_exact_read, hp_mismatch_rate_thresholds);

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    auto ref_bp2_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp2_seq, ref_bp2_better_reads, std::vector<bool>(), ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score, ref_bp2_is_exact_read, hp_mismatch_rate_thresholds);

    auto score_ins_consensus = [&](const std::string& consensus_seq, bool bp1) {
        char* lf_seq = generate_haplotype_left(contig_seq, ins_start-1, consensus_seq.length(), ins->aux_indels, ins->aux_snps);
        int alt_lf_len = strlen(lf_seq);
        char* rf_seq = generate_haplotype_right(contig_seq, contig_len, ins_end, consensus_seq.length(), ins->aux_indels, ins->aux_snps);
        int alt_rf_len = strlen(rf_seq);
        int ins_seq_portion_len = std::min(ins->ins_seq.length(), consensus_seq.length());
        int extra_len = std::max(0, int(consensus_seq.length())-int(ins->ins_seq.length()));

        alignment_targets_t targets;
        if (bp1) {
            extra_len = std::min(extra_len, alt_rf_len);
            targets.alt_len = alt_lf_len+ins_seq_portion_len+extra_len;
            targets.alt_seq = new char[targets.alt_len+1];
            strncpy(targets.alt_seq, lf_seq, alt_lf_len);
            strncpy(targets.alt_seq+alt_lf_len, ins->ins_seq.c_str(), ins_seq_portion_len);
            strncpy(targets.alt_seq+alt_lf_len+ins_seq_portion_len, rf_seq, extra_len);
            targets.left_flank_end = alt_lf_len;
            targets.right_flank_start = alt_lf_len+ins_seq_portion_len;
        } else {
            extra_len = std::min(extra_len, alt_lf_len);
            targets.alt_len = extra_len+ins_seq_portion_len+alt_rf_len;
            targets.alt_seq = new char[targets.alt_len+1];
            strncpy(targets.alt_seq, lf_seq+(alt_lf_len-extra_len), extra_len);
            strncpy(targets.alt_seq+extra_len, ins->ins_seq.c_str()+(ins->ins_seq.length()-ins_seq_portion_len), ins_seq_portion_len);
            strncpy(targets.alt_seq+extra_len+ins_seq_portion_len, rf_seq, alt_rf_len);
            targets.left_flank_end = extra_len;
            targets.right_flank_start = extra_len+ins_seq_portion_len;
        }
        targets.alt_seq[targets.alt_len] = 0;

        hts_pos_t breakpoint = bp1 ? ins_start : ins_end;
        hts_pos_t ref_start = std::max(hts_pos_t(0), breakpoint-hts_pos_t(consensus_seq.length()));
        hts_pos_t ref_end = std::min(breakpoint+hts_pos_t(consensus_seq.length()), contig_len);
        targets.ref_seqs.push_back(contig_seq+ref_start);
        targets.ref_lens.push_back(ref_end-ref_start);

        char* independent_ref_seq = concat3(lf_seq, contig_seq+ins_start, rf_seq, alt_lf_len, ins_end-ins_start, alt_rf_len);
        targets.left_independent_ref_seq = targets.right_independent_ref_seq = independent_ref_seq;
        targets.left_independent_ref_len = targets.right_independent_ref_len = alt_lf_len+ins_end-ins_start+alt_rf_len;

        consensus_alignment_metrics_t metrics = score_consensus_alignment(consensus_seq, targets, aligner);
        delete[] targets.alt_seq;
        delete[] independent_ref_seq;
        delete[] lf_seq;
        delete[] rf_seq;
        return metrics;
    };

    ins->sample_info.ext_alt_consensus1_metrics.length = alt_bp1_consensus_seq.length();
    if (alt_bp1_consensus_seq.length() >= 2*config.min_clip_len) {
        consensus_alignment_metrics_t unextended_metrics = score_ins_consensus(alt_bp1_consensus_seq, true);
        ins->sample_info.alt_consensus1_metrics = unextended_metrics;

        // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_bp1_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_bp1_consensus_seq, std::string(alt_bp1_consensus_seq.length(), '!'), 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp1_consensus, candidate_reads_for_extension_itree, std::max<hts_pos_t>(0, ins_start-GENOTYPE_CONSENSUS_EXTENSION), ins_start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        extend_consensus_to_right(alt_bp1_consensus, candidate_reads_for_extension_itree, ins_start, std::min<hts_pos_t>(contig_len, ins_start+GENOTYPE_CONSENSUS_EXTENSION), contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        ins->sample_info.alt_lext_reads = alt_bp1_consensus->left_ext_reads;
        ins->sample_info.alt_rext_reads = alt_bp1_consensus->right_ext_reads;
        ins->sample_info.hq_alt_lext_reads = alt_bp1_consensus->hq_left_ext_reads;
        ins->sample_info.hq_alt_rext_reads = alt_bp1_consensus->hq_right_ext_reads;
        alt_bp1_consensus_seq = alt_bp1_consensus->sequence;

        consensus_alignment_metrics_t extended_metrics = score_ins_consensus(alt_bp1_consensus_seq, true);
        ins->sample_info.ext_alt_consensus1_metrics = extended_metrics;

        int lf_aln_rlen = extended_metrics.split_ref_lengths[0];
        ins->left_anchor_aln->start = ins_start-lf_aln_rlen;
        ins->left_anchor_aln->end = ins_start;
        ins->left_anchor_aln->seq_len = lf_aln_rlen;
    }

    ins->sample_info.ext_alt_consensus2_metrics.length = alt_bp2_consensus_seq.length();
    if (alt_bp2_consensus_seq.length() >= 2*config.min_clip_len) {
        consensus_alignment_metrics_t unextended_metrics = score_ins_consensus(alt_bp2_consensus_seq, false);
        ins->sample_info.alt_consensus2_metrics = unextended_metrics;

        std::shared_ptr<consensus_t> alt_bp2_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_bp2_consensus_seq, std::string(alt_bp2_consensus_seq.length(), '!'), 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp2_consensus, candidate_reads_for_extension_itree, std::max<hts_pos_t>(0, ins_end-GENOTYPE_CONSENSUS_EXTENSION), ins_end, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        extend_consensus_to_right(alt_bp2_consensus, candidate_reads_for_extension_itree, ins_end, std::min<hts_pos_t>(contig_len, ins_end+GENOTYPE_CONSENSUS_EXTENSION), contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        ins->sample_info.alt_lext_reads += alt_bp2_consensus->left_ext_reads;
        ins->sample_info.alt_rext_reads += alt_bp2_consensus->right_ext_reads;
        ins->sample_info.hq_alt_lext_reads += alt_bp2_consensus->hq_left_ext_reads;
        ins->sample_info.hq_alt_rext_reads += alt_bp2_consensus->hq_right_ext_reads;
        alt_bp2_consensus_seq = alt_bp2_consensus->sequence;

        consensus_alignment_metrics_t extended_metrics = score_ins_consensus(alt_bp2_consensus_seq, false);
        ins->sample_info.ext_alt_consensus2_metrics = extended_metrics;

        int rf_aln_rlen = extended_metrics.split_ref_lengths[1];
        ins->right_anchor_aln->start = ins_end;
        ins->right_anchor_aln->end = ins_end+rf_aln_rlen;
        ins->right_anchor_aln->seq_len = rf_aln_rlen;
    }

    delete[] ref_bp1_seq;
    delete[] ref_bp2_seq;

    set_bp_consensus_info(ins->sample_info.alt_bp1.reads_info, alt_bp1_better_reads, alt_bp1_is_consistent_read, alt_bp1_is_exact_read, alt_bp1_avg_score, alt_bp1_stddev_score);
    set_bp_consensus_info(ins->sample_info.alt_bp2.reads_info, alt_bp2_better_reads, alt_bp2_is_consistent_read, alt_bp2_is_exact_read, alt_bp2_avg_score, alt_bp2_stddev_score);
    set_bp_consensus_info(ins->sample_info.ref_bp1.reads_info, ref_bp1_better_reads, ref_bp1_is_consistent_read, ref_bp1_is_exact_read, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(ins->sample_info.ref_bp2.reads_info, ref_bp2_better_reads, ref_bp2_is_consistent_read, ref_bp2_is_exact_read, ref_bp2_avg_score, ref_bp2_stddev_score);

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

inline void genotype_inss(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<insertion_t*> inss,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t& stats, config_t& config, contig_map_t& contig_map,
    bam_pool_t* bam_pool, std::vector<double>* global_crossing_isize_dist, const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

    int contig_id = contig_map.get_id(contig_name);
    auto chromosome_data = acquire_chromosome_data(contig_id);
    auto mateseqs_w_mapq_chr = chromosome_data.first;
    evidence_map_t* evidence_map = chromosome_data.second;

    std::vector<hts_pair_pos_t> target_ivals;
    for (insertion_t* ins : inss) {
        target_ivals.push_back({std::max<hts_pos_t>(0, ins->start-GENOTYPE_CONSENSUS_EXTENSION), std::min<hts_pos_t>(contig_len, ins->start+GENOTYPE_CONSENSUS_EXTENSION)});
        target_ivals.push_back({std::max<hts_pos_t>(0, ins->end-GENOTYPE_CONSENSUS_EXTENSION), std::min<hts_pos_t>(contig_len, ins->end+GENOTYPE_CONSENSUS_EXTENSION)});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_pool->get_bam_reader(id), candidate_reads_for_extension);

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);
    for (insertion_t* ins : inss) {
        genotype_ins(ins, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_map, hp_mismatch_rate_thresholds);
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    depth_filter_ins(contig_name, inss, bam_file, config, stats);
    calculate_ptn_ratio(contig_name, inss, bam_file, config, stats, false, evidence_map, *mateseqs_w_mapq_chr);

    std::vector<sv_t*> inss_sv(inss.begin(), inss.end());
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, inss_sv, bam_file, config, stats);

    release_chromosome_data(contig_id);

}
#endif // GENOTYPE_INSS_H
