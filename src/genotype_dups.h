#ifndef GENOTYPE_DUPS_H
#define GENOTYPE_DUPS_H

#include "stat_tests.h"

#include "genotype.h"

struct small_dup_alignment_targets_t {
    hts_pos_t dup_start, dup_end, ref_start, ref_end, ref_len, svlen;
    std::unique_ptr<char[]> ref_seq;
    std::vector<std::unique_ptr<char[]>> alt_alleles;
    std::vector<char*> alt_seqs;
};

small_dup_alignment_targets_t build_small_dup_alignment_targets(duplication_t* dup, char* contig_seq, hts_pos_t contig_len, stats_t& stats) {
    small_dup_alignment_targets_t targets;
    targets.dup_start = dup->start + 1;
    targets.dup_end = std::min(dup->end + 1, contig_len);
    targets.svlen = dup->svlen();
    hts_pos_t extend = stats.read_len + 20;
    targets.ref_start = std::max(hts_pos_t(0), targets.dup_start-extend);
    targets.ref_end = std::min(targets.dup_end+extend, contig_len);
    targets.ref_len = targets.ref_end - targets.ref_start;
    targets.ref_seq.reset(new char[targets.ref_len + 1]);
    strncpy(targets.ref_seq.get(), contig_seq+targets.ref_start, targets.ref_len);
    targets.ref_seq[targets.ref_len] = 0;

    for (int copies = 1; copies*targets.svlen < stats.read_len; copies++) {
        int alt_len = targets.ref_len + copies*targets.svlen;
        std::unique_ptr<char[]> alt_seq(new char[alt_len+1]);
        int pos = 0;
        strncpy(alt_seq.get(), contig_seq+targets.ref_start, targets.dup_end-targets.ref_start);
        pos += targets.dup_end-targets.ref_start;
        for (int i = 0; i < copies; i++) {
            strncpy(alt_seq.get()+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
            strncpy(alt_seq.get()+pos, contig_seq+targets.dup_start, targets.dup_end-targets.dup_start);
            pos += targets.dup_end-targets.dup_start;
        }
        strncpy(alt_seq.get()+pos, contig_seq+targets.dup_end, targets.ref_end-targets.dup_end);
        pos += targets.ref_end-targets.dup_end;
        alt_seq[pos] = 0;
        targets.alt_seqs.push_back(alt_seq.get());
        targets.alt_alleles.push_back(std::move(alt_seq));
    }
    return targets;
}

struct large_dup_alignment_targets_t {
    hts_pos_t dup_start, dup_end;
    hts_pos_t ref_bp1_start, ref_bp1_end, ref_bp1_len, ref_bp1_pos;
    hts_pos_t ref_bp2_start, ref_bp2_end, ref_bp2_len, ref_bp2_pos;
    hts_pos_t alt_len, alt_lh_len, alt_rh_len;
    std::unique_ptr<char[]> alt_seq;
};

large_dup_alignment_targets_t build_large_dup_alignment_targets(duplication_t* dup, char* contig_seq, hts_pos_t contig_len, stats_t& stats) {
    large_dup_alignment_targets_t targets;
    targets.dup_start = dup->start + 1;
    targets.dup_end = std::min(dup->end + 1, contig_len);
    hts_pos_t extend = stats.read_len + 20;
    targets.ref_bp1_start = std::max(hts_pos_t(0), targets.dup_start-extend);
    targets.ref_bp1_end = std::min(targets.dup_start+extend, contig_len);
    targets.ref_bp1_pos = targets.dup_start - targets.ref_bp1_start;
    targets.ref_bp1_len = targets.ref_bp1_end - targets.ref_bp1_start;
    targets.ref_bp2_start = std::max(hts_pos_t(0), targets.dup_end-extend);
    targets.ref_bp2_end = std::min(targets.dup_end+extend, contig_len);
    targets.ref_bp2_pos = targets.dup_end - targets.ref_bp2_start;
    targets.ref_bp2_len = targets.ref_bp2_end - targets.ref_bp2_start;
    targets.alt_lh_len = targets.dup_end - targets.ref_bp2_start;
    targets.alt_rh_len = targets.ref_bp1_end - targets.dup_start;
    targets.alt_len = targets.alt_lh_len + dup->ins_seq.length() + targets.alt_rh_len;
    targets.alt_seq.reset(new char[targets.alt_len + 1]);
    strncpy(targets.alt_seq.get(), contig_seq+targets.ref_bp2_start, targets.alt_lh_len);
    strncpy(targets.alt_seq.get()+targets.alt_lh_len, dup->ins_seq.c_str(), dup->ins_seq.length());
    strncpy(targets.alt_seq.get()+targets.alt_lh_len+dup->ins_seq.length(), contig_seq+targets.dup_start, targets.alt_rh_len);
    targets.alt_seq[targets.alt_len] = 0;
    to_uppercase(targets.alt_seq.get());
    return targets;
}

struct small_dup_read_evidence_t {
    std::vector<int> alt_idxs, alt_scores, alt_positions;
    bool ref = false, er = false;
};

small_dup_read_evidence_t get_cached_small_dup_read_evidence(duplication_t* dup, bam1_t* read, size_t n_alt_seqs, evidence_map_t* evidence_map) {
    small_dup_read_evidence_t result;
    const std::vector<evidence_map_t::read_alt_association_t>& associations = evidence_map->get_read_alt_associations(read);
    std::string read_name = read_name_with_suffix(read);
    for (const auto& association : associations) {
        if (association.sv != dup || association.bp != 1 || association.alt_idx < 0 || !evidence_map->is_read_alt_association(read_name, association)) continue;
        if (association.alt_idx >= n_alt_seqs) throw std::runtime_error("Invalid cached ALT index for small duplication " + dup->id + ".");
        result.alt_idxs.push_back(association.alt_idx);
        result.alt_scores.push_back(0);
        result.alt_positions.push_back(association.pos);
    }
    if (result.alt_idxs.empty()) result.ref = evidence_map->is_read_ref(read, dup, 1);
    if (result.alt_idxs.empty() && !result.ref) result.er = evidence_map->is_read_er(read, dup);
    return result;
}

struct large_dup_read_evidence_t {
    bool alt = false, ref_bp1 = false, ref_bp2 = false, er = false;
    int alt_score = 0, alt_position = -1;
    bool alt_spans_bp1 = false, alt_spans_bp2 = false;
};

large_dup_read_evidence_t get_cached_large_dup_read_evidence(duplication_t* dup, bam1_t* read, evidence_map_t* evidence_map) {
    large_dup_read_evidence_t result;
    int bp0_pos = evidence_map->get_read_alt_pos(read, dup, 0), bp1_pos = evidence_map->get_read_alt_pos(read, dup, 1), bp2_pos = evidence_map->get_read_alt_pos(read, dup, 2);
    if (bp1_pos >= 0 && bp2_pos >= 0 && bp1_pos != bp2_pos) throw std::runtime_error("Different cached ALT alignment positions for duplication " + dup->id + " and read " + read_name_with_suffix(read) + ".");
    result.alt = bp0_pos >= 0 || bp1_pos >= 0 || bp2_pos >= 0;
    if (result.alt) {
        result.alt_position = bp1_pos >= 0 ? bp1_pos : (bp2_pos >= 0 ? bp2_pos : bp0_pos);
        result.alt_spans_bp1 = bp1_pos >= 0;
        result.alt_spans_bp2 = bp2_pos >= 0;
    } else {
        result.ref_bp1 = evidence_map->is_read_ref(read, dup, 1);
        result.ref_bp2 = evidence_map->is_read_ref(read, dup, 2);
        if (!result.ref_bp1 && !result.ref_bp2) result.er = evidence_map->is_read_er(read, dup);
    }
    return result;
}

void write_aligned_small_dup_read_evidence(duplication_t* dup, open_samFile_t* bam_file, char* contig_seq, hts_pos_t contig_len, stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t& evidence_logger) {
    small_dup_alignment_targets_t targets = build_small_dup_alignment_targets(dup, contig_seq, contig_len, stats);
    std::stringstream region;
    region << dup->chr << ":" << targets.ref_start << "-" << targets.ref_end;
    hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region.str().c_str());
    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> ref_reads;
    std::vector<std::string> er_read_names;
    std::vector<std::vector<std::shared_ptr<bam1_t>>> alt_reads(targets.alt_seqs.size());
    std::vector<std::vector<int>> alt_positions(targets.alt_seqs.size()), alt_scores(targets.alt_seqs.size());
    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Filter filter_with_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Alignment alt_aln, ref_aln;

    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < targets.dup_start || targets.dup_end < get_unclipped_start(read)) continue;
        if (targets.dup_start < get_unclipped_start(read) && get_unclipped_end(read) < targets.dup_end) continue;
        if (!is_samechr(read) || is_samestr(read)) continue;

        bool reverse_seq = !bam_is_mrev(read);
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < targets.dup_start-stats.max_is || read->core.mpos > targets.dup_end) continue;
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > targets.dup_end+stats.max_is || mate_endpos < targets.dup_start) continue;
        }

        std::string seq = get_sequence(read, true);
        if (reverse_seq) rc(seq);
        aligner.Align(seq.c_str(), targets.ref_seq.get(), targets.ref_len, filter_with_score_only, &ref_aln, 0);
        uint16_t best_aln_score = 0;
        std::vector<uint16_t> alt_aln_scores(targets.alt_seqs.size());
        std::vector<int> alt_aln_positions(targets.alt_seqs.size());
        for (int i = 0; i < targets.alt_seqs.size(); i++) {
            aligner.Align(seq.c_str(), targets.alt_seqs[i], strlen(targets.alt_seqs[i]), filter_with_pos, &alt_aln, 0);
            alt_aln_scores[i] = alt_aln.sw_score;
            alt_aln_positions[i] = alt_aln.ref_begin;
            if (alt_aln.sw_score > best_aln_score) best_aln_score = alt_aln.sw_score;
        }
        if (best_aln_score > ref_aln.sw_score) {
            for (int i = 0; i < targets.alt_seqs.size(); i++) {
                if (alt_aln_scores[i] != best_aln_score) continue;
                alt_reads[i].push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                alt_scores[i].push_back(alt_aln_scores[i]);
                alt_positions[i].push_back(alt_aln_positions[i]);
            }
        } else if (best_aln_score < ref_aln.sw_score) {
            ref_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        } else {
            er_read_names.push_back(read_name_with_suffix(read));
        }
    }

    int alt_with_most_reads = 0;
    for (int i = 1; i < alt_reads.size(); i++) {
        if (alt_reads[i].size() > alt_reads[alt_with_most_reads].size()) alt_with_most_reads = i;
    }
    if (alt_reads[alt_with_most_reads].size() > 20*stats.get_max_depth(dup->chr) || ref_reads.size() > 20*stats.get_max_depth(dup->chr)) dup->sample_info.too_deep = true;

    bam_destroy1(read);
    hts_itr_destroy(iter);
    if (dup->sample_info.too_deep) return;

    for (int i = 0; i < alt_reads.size(); i++) evidence_logger.log_reads_associations(dup->chr, dup->id, 1, alt_reads[i], alt_scores[i], alt_positions[i], i);
    evidence_logger.log_ref_reads_associations(dup->chr, dup->id, 1, ref_reads);
    evidence_logger.log_er_reads_associations(dup->chr, dup->id, er_read_names);
}

void write_aligned_large_dup_read_evidence(duplication_t* dup, open_samFile_t* bam_file, char* contig_seq, hts_pos_t contig_len, stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t& evidence_logger) {
    large_dup_alignment_targets_t targets = build_large_dup_alignment_targets(dup, contig_seq, contig_len, stats);
    std::stringstream l_region, r_region;
    l_region << dup->chr << ":" << targets.ref_bp1_start << "-" << targets.ref_bp1_end;
    r_region << dup->chr << ":" << targets.ref_bp2_start << "-" << targets.ref_bp2_end;
    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> alt_reads, ref_bp1_reads, ref_bp2_reads;
    std::vector<std::string> er_read_names;
    std::vector<int> alt_positions, alt_scores;
    std::vector<bool> alt_spans_bp1, alt_spans_bp2;
    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;

    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < targets.dup_start || targets.dup_end < get_unclipped_start(read)) continue;
        if (targets.dup_start < get_unclipped_start(read) && get_unclipped_end(read) < targets.dup_end) continue;
        if (!is_samechr(read) || is_samestr(read)) continue;

        bool reverse_seq = !bam_is_mrev(read);
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < targets.dup_start-stats.max_is || read->core.mpos > targets.dup_end) continue;
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > targets.dup_end+stats.max_is || mate_endpos < targets.dup_start) continue;
        }

        std::string seq = get_sequence(read, true);
        if (reverse_seq) rc(seq);
        uint16_t ref_aln_score = 0;
        bool ref_bp1_better = false, ref_bp2_better = false;
        bool ref_is_exact_match = is_perfectly_aligned(read);
        if (ref_is_exact_match) {
            ref_aln_score = read->core.l_qseq;
            if (read->core.pos < targets.dup_start && bam_endpos(read) > targets.dup_start) ref_bp1_better = true;
            if (read->core.pos < targets.dup_end && bam_endpos(read) > targets.dup_end) ref_bp2_better = true;
        } else {
            aligner.Align(seq.c_str(), contig_seq+targets.ref_bp1_start, targets.ref_bp1_len, filter_with_pos, &ref1_aln, 0);
            aligner.Align(seq.c_str(), contig_seq+targets.ref_bp2_start, targets.ref_bp2_len, filter_with_pos, &ref2_aln, 0);
            ref_aln_score = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln.sw_score : ref2_aln.sw_score;
            ref_bp1_better = ref1_aln.sw_score >= ref2_aln.sw_score && ref1_aln.ref_begin <= targets.ref_bp1_pos && ref1_aln.ref_end >= targets.ref_bp1_pos;
            ref_bp2_better = ref2_aln.sw_score >= ref1_aln.sw_score && ref2_aln.ref_begin <= targets.ref_bp2_pos && ref2_aln.ref_end >= targets.ref_bp2_pos;
        }
        alt_aln = align_fast(aligner, seq.c_str(), targets.alt_seq.get(), targets.alt_len, filter_with_pos, ref_is_exact_match);
        hts_pos_t alt_right_flank_pos = targets.alt_lh_len + dup->ins_seq.length();
        bool spans_bp1 = alt_aln.ref_begin < targets.alt_lh_len && alt_aln.ref_end >= targets.alt_lh_len;
        bool spans_bp2 = alt_aln.ref_begin < alt_right_flank_pos && alt_aln.ref_end >= alt_right_flank_pos;
        if (alt_aln.sw_score > ref_aln_score) {
            alt_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_positions.push_back(alt_aln.ref_begin);
            alt_scores.push_back(alt_aln.sw_score);
            alt_spans_bp1.push_back(spans_bp1);
            alt_spans_bp2.push_back(spans_bp2);
        } else {
            if (alt_aln.sw_score < ref_aln_score && ref_bp1_better) ref_bp1_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            if (alt_aln.sw_score < ref_aln_score && ref_bp2_better) ref_bp2_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            if (alt_aln.sw_score == ref_aln_score) er_read_names.push_back(read_name_with_suffix(read));
        }

        if (ref_bp1_reads.size() + ref_bp2_reads.size() + er_read_names.size() > 4*stats.get_max_depth(dup->chr)) {
            dup->sample_info.too_deep = true;
            break;
        }
    }

    free(regions[0]);
    free(regions[1]);
    bam_destroy1(read);
    hts_itr_destroy(iter);
    if (dup->sample_info.too_deep) return;

    std::vector<std::shared_ptr<bam1_t>> alt_bp0_reads, alt_bp1_reads, alt_bp2_reads;
    std::vector<int> alt_bp0_scores, alt_bp1_scores, alt_bp2_scores, alt_bp0_positions, alt_bp1_positions, alt_bp2_positions;
    for (size_t i = 0; i < alt_reads.size(); i++) {
        if (!alt_spans_bp1[i] && !alt_spans_bp2[i]) {
            alt_bp0_reads.push_back(alt_reads[i]);
            alt_bp0_scores.push_back(alt_scores[i]);
            alt_bp0_positions.push_back(alt_positions[i]);
        }
        if (alt_spans_bp1[i]) {
            alt_bp1_reads.push_back(alt_reads[i]);
            alt_bp1_scores.push_back(alt_scores[i]);
            alt_bp1_positions.push_back(alt_positions[i]);
        }
        if (alt_spans_bp2[i]) {
            alt_bp2_reads.push_back(alt_reads[i]);
            alt_bp2_scores.push_back(alt_scores[i]);
            alt_bp2_positions.push_back(alt_positions[i]);
        }
    }
    evidence_logger.log_reads_associations(dup->chr, dup->id, 0, alt_bp0_reads, alt_bp0_scores, alt_bp0_positions);
    evidence_logger.log_reads_associations(dup->chr, dup->id, 1, alt_bp1_reads, alt_bp1_scores, alt_bp1_positions);
    evidence_logger.log_reads_associations(dup->chr, dup->id, 2, alt_bp2_reads, alt_bp2_scores, alt_bp2_positions);
    evidence_logger.log_ref_reads_associations(dup->chr, dup->id, 1, ref_bp1_reads);
    evidence_logger.log_ref_reads_associations(dup->chr, dup->id, 2, ref_bp2_reads);
    evidence_logger.log_er_reads_associations(dup->chr, dup->id, er_read_names);
}

inline void genotype_small_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_map_t* evidence_map,
                const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {
    small_dup_alignment_targets_t alignment_targets = build_small_dup_alignment_targets(dup, contig_seq, contig_len, stats);
    hts_pos_t dup_start = alignment_targets.dup_start, dup_end = alignment_targets.dup_end;
    hts_pos_t ref_start = alignment_targets.ref_start, ref_end = alignment_targets.ref_end, ref_len = alignment_targets.ref_len;
    hts_pos_t svlen = alignment_targets.svlen;
    char* ref_seq = alignment_targets.ref_seq.get();
    std::vector<char*>& alt_seqs = alignment_targets.alt_seqs;
    dup->ref1_hp_len = longest_homopolymer_len(ref_seq, ref_len);
    if (dup->sample_info.too_deep) {
        dup->alt1_hp_len = longest_homopolymer_len(alt_seqs[0], strlen(alt_seqs[0]));
        return;
    }

    std::stringstream region;
    region << dup->chr << ":" << ref_start << "-" << ref_end;

	hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region.str().c_str());
	
    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> ref_better_reads;
    std::vector<std::vector<std::shared_ptr<bam1_t>>> alt_better_reads(alt_seqs.size());
    std::vector<std::vector<int>> alt_better_read_positions(alt_seqs.size());
    std::vector<std::vector<int>> alt_better_reads_scores(alt_seqs.size());

    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;
        if (!is_samechr(read) || is_samestr(read)) continue;

        bool reverse_seq = !bam_is_mrev(read);
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < dup_start-stats.max_is || read->core.mpos > dup_end) continue;
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > dup_end+stats.max_is || mate_endpos < dup_start) continue;
        }

        small_dup_read_evidence_t read_evidence = get_cached_small_dup_read_evidence(dup, read, alt_seqs.size(), evidence_map);

        for (size_t i = 0; i < read_evidence.alt_idxs.size(); i++) {
            int alt_idx = read_evidence.alt_idxs[i];
            alt_better_reads[alt_idx].push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_better_reads_scores[alt_idx].push_back(read_evidence.alt_scores[i]);
            alt_better_read_positions[alt_idx].push_back(read_evidence.alt_positions[i]);
        }
        if (read_evidence.ref) ref_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        if (read_evidence.er) {
            dup->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) dup->sample_info.alt_ref_equal_reads_highmq++;
        }
    }

    int alt_with_most_reads = 0;
    for (int i = 1; i < alt_better_reads.size(); i++) {
        if (alt_better_reads[i].size() > alt_better_reads[alt_with_most_reads].size()) alt_with_most_reads = i;
    }

    if (alt_better_reads[alt_with_most_reads].size() > 20*stats.get_max_depth(dup->chr) || ref_better_reads.size() > 20*stats.get_max_depth(dup->chr)) {
        alt_better_reads[alt_with_most_reads].clear();
        alt_better_read_positions[alt_with_most_reads].clear();
        alt_better_reads_scores[alt_with_most_reads].clear();
        ref_better_reads.clear();
        dup->sample_info.alt_ref_equal_reads = 0;
        dup->sample_info.alt_ref_equal_reads_highmq = 0;
        evidence_map->clear_other_read_support_for_too_deep(dup->sample_info);
    }

    std::vector<char*> ref_seqs = {ref_seq};
    std::vector<hts_pos_t> ref_lens = {ref_len};
    int alt_len = strlen(alt_seqs[alt_with_most_reads]);
    dup->alt1_hp_len = longest_homopolymer_len(alt_seqs[alt_with_most_reads], alt_len);
    std::vector<hts_pos_t> alt_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_seqs[alt_with_most_reads], alt_len, stats.read_len);
    dup->sample_info.expected_alt1_reads_frac = (double) alt_ref_diff_reads_expected_positions.size() / std::max(1, alt_len - stats.read_len + 1);
    dup->sample_info.max_feasible_alt1_len = get_max_feasible_alt_len(alt_ref_diff_reads_expected_positions, stats.read_len);

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    double alt_stddev_score, ref_stddev_score;
    std::vector<bool> alt_is_exact_read, ref_is_exact_read;
    auto alt_is_consistent_read = gen_consensus_and_classify_seqs(alt_seqs[alt_with_most_reads], alt_better_reads[alt_with_most_reads], std::vector<bool>(), alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read, hp_mismatch_rate_thresholds);
    auto ref_is_consistent_read = gen_consensus_and_classify_seqs(ref_seq, ref_better_reads, std::vector<bool>(), ref_consensus_seq, ref_avg_score, ref_stddev_score, ref_is_exact_read, hp_mismatch_rate_thresholds);

    std::vector<int> alt_better_read_positions_consistent = get_consistent_reads_start_positions(alt_is_consistent_read, alt_better_read_positions[alt_with_most_reads]);
    dup->sample_info.alt1_occ_ratio = occ_ratio(alt_better_read_positions_consistent, alt_ref_diff_reads_expected_positions.size());

    auto score_dup_consensus = [&](const std::string& consensus_seq) {
        alignment_targets_t targets;
        hts_pos_t ref_start = std::max(hts_pos_t(0), dup->start-hts_pos_t(consensus_seq.length()));
        hts_pos_t ref_end = std::min(dup->end+hts_pos_t(consensus_seq.length()), contig_len);
        int n_extra_copies = alt_with_most_reads+1;
        targets.alt_len = ref_end-ref_start+n_extra_copies*svlen;
        targets.alt_seq = new char[targets.alt_len+1];
        int pos = 0;
        strncpy(targets.alt_seq, contig_seq+ref_start, dup_end-ref_start);
        pos += dup_end-ref_start;
        for (int i = 0; i < n_extra_copies; i++) {
            strncpy(targets.alt_seq+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
            strncpy(targets.alt_seq+pos, contig_seq+dup_start, dup_end-dup_start);
            pos += dup_end-dup_start;
        }
        strncpy(targets.alt_seq+pos, contig_seq+dup_end, ref_end-dup_end);
        pos += ref_end-dup_end;
        targets.alt_seq[pos] = 0;
        targets.ref_seqs.push_back(contig_seq+ref_start);
        targets.ref_lens.push_back(ref_end-ref_start);
        targets.left_flank_end = dup_start-ref_start;
        targets.right_flank_start = pos-(ref_end-dup_end);
        targets.right_flank_end_offset = 0;
        targets.left_independent_ref_seq = targets.right_independent_ref_seq = contig_seq+ref_start;
        targets.left_independent_ref_len = targets.right_independent_ref_len = ref_end-ref_start;
        consensus_alignment_metrics_t metrics = score_consensus_alignment(consensus_seq, targets, aligner);
        delete[] targets.alt_seq;
        return metrics;
    };

    if (alt_consensus_seq.length() >= 2*config.min_clip_len) {
        consensus_alignment_metrics_t unextended_metrics = score_dup_consensus(alt_consensus_seq);
        dup->sample_info.alt_consensus1_metrics = unextended_metrics;

        std::shared_ptr<consensus_t> alt_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, std::max<hts_pos_t>(0, dup->start-GENOTYPE_CONSENSUS_EXTENSION), dup->start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, dup->end, std::min<hts_pos_t>(contig_len, dup->end+GENOTYPE_CONSENSUS_EXTENSION), contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        dup->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        dup->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        dup->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        dup->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;

        consensus_alignment_metrics_t extended_metrics = score_dup_consensus(alt_consensus_seq);
        dup->sample_info.ext_alt_consensus1_metrics = extended_metrics;

        int lf_aln_rlen = extended_metrics.split_ref_lengths[0];
        int rf_aln_rlen = extended_metrics.split_ref_lengths[1];
        dup->left_anchor_aln->start = dup_end-lf_aln_rlen;
        dup->left_anchor_aln->end = dup_end;
        dup->left_anchor_aln->seq_len = lf_aln_rlen;
        dup->right_anchor_aln->start = dup_start;
        dup->right_anchor_aln->end = dup_start+rf_aln_rlen;
        dup->right_anchor_aln->seq_len = rf_aln_rlen;
    }

    alt_is_consistent_read = classify_seqs_with_ref_seq(alt_consensus_seq, alt_better_reads[alt_with_most_reads],
        alt_is_consistent_read, alt_avg_score, alt_stddev_score, alt_is_exact_read);

    for (int i = 0; i < alt_better_reads[alt_with_most_reads].size(); i++) {
        std::shared_ptr<bam1_t>& r = alt_better_reads[alt_with_most_reads][i];
        evidence_map->record_assigned_alt_read(dup, r.get(), alt_is_consistent_read[i], get_mq(r.get()) >= config.high_confidence_mapq, alt_is_exact_read[i]);
    }

    set_bp_consensus_info(dup->sample_info.alt_bp1.reads_info, alt_better_reads[alt_with_most_reads], alt_is_consistent_read, alt_is_exact_read, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp1.reads_info, ref_better_reads, ref_is_consistent_read, ref_is_exact_read, ref_avg_score, ref_stddev_score);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

inline void genotype_large_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_map_t* evidence_map,
                const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {

    large_dup_alignment_targets_t alignment_targets = build_large_dup_alignment_targets(dup, contig_seq, contig_len, stats);
    hts_pos_t dup_start = alignment_targets.dup_start, dup_end = alignment_targets.dup_end;
    hts_pos_t ref_bp1_start = alignment_targets.ref_bp1_start, ref_bp1_end = alignment_targets.ref_bp1_end, ref_bp1_len = alignment_targets.ref_bp1_len, ref_bp1_pos = alignment_targets.ref_bp1_pos;
    hts_pos_t ref_bp2_start = alignment_targets.ref_bp2_start, ref_bp2_end = alignment_targets.ref_bp2_end, ref_bp2_len = alignment_targets.ref_bp2_len, ref_bp2_pos = alignment_targets.ref_bp2_pos;
    hts_pos_t alt_lh_len = alignment_targets.alt_lh_len, alt_rh_len = alignment_targets.alt_rh_len, alt_len = alignment_targets.alt_len;
    char* alt_seq = alignment_targets.alt_seq.get();
    dup->ref1_hp_len = longest_homopolymer_len(contig_seq+ref_bp1_start, ref_bp1_len);
    dup->ref2_hp_len = longest_homopolymer_len(contig_seq+ref_bp2_start, ref_bp2_len);
    dup->alt1_hp_len = dup->alt2_hp_len = longest_homopolymer_len(alt_seq, alt_len);
    if (dup->sample_info.too_deep) return;

    std::vector<char*> ref_seqs = {contig_seq+ref_bp1_start, contig_seq+ref_bp2_start};
    std::vector<hts_pos_t> ref_lens = {ref_bp1_len, ref_bp2_len};
    std::vector<hts_pos_t> alt_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_seq, alt_len, stats.read_len);
    dup->sample_info.expected_alt1_reads_frac = (double) alt_ref_diff_reads_expected_positions.size() / std::max(hts_pos_t(1), alt_len - stats.read_len + 1);
    dup->sample_info.max_feasible_alt1_len = get_max_feasible_alt_len(alt_ref_diff_reads_expected_positions, stats.read_len);

    std::stringstream l_region, r_region;
    l_region << dup->chr << ":" << ref_bp1_start << "-" << ref_bp1_end;
    r_region << dup->chr << ":" << ref_bp2_start << "-" << ref_bp2_end;
    
    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);

    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> alt_better_reads, ref_bp1_better_reads, ref_bp2_better_reads;
    std::vector<int> alt_better_read_positions;
    std::vector<int> alt_better_reads_scores;
    std::vector<bool> alt_better_read_spans_bp1, alt_better_read_spans_bp2;

    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;
        if (!is_samechr(read) || is_samestr(read)) continue;

        bool reverse_seq = !bam_is_mrev(read);
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < dup_start-stats.max_is || read->core.mpos > dup_end) continue;
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > dup_end+stats.max_is || mate_endpos < dup_start) continue;
        }

        large_dup_read_evidence_t read_evidence = get_cached_large_dup_read_evidence(dup, read, evidence_map);

        if (read_evidence.alt) {
            alt_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_better_read_positions.push_back(read_evidence.alt_position);
            alt_better_reads_scores.push_back(read_evidence.alt_score);
            alt_better_read_spans_bp1.push_back(read_evidence.alt_spans_bp1);
            alt_better_read_spans_bp2.push_back(read_evidence.alt_spans_bp2);
        } else {
            if (read_evidence.ref_bp1) ref_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            if (read_evidence.ref_bp2) ref_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            if (read_evidence.er) {
                dup->sample_info.alt_ref_equal_reads++;
                if (read->core.qual >= config.high_confidence_mapq) dup->sample_info.alt_ref_equal_reads_highmq++;
            }
        }

        if (ref_bp1_better_reads.size() + ref_bp2_better_reads.size() + dup->sample_info.alt_ref_equal_reads > 4*stats.get_max_depth(dup->chr)) {
            alt_better_reads.clear();
            alt_better_read_positions.clear();
            alt_better_reads_scores.clear();
            alt_better_read_spans_bp1.clear();
            alt_better_read_spans_bp2.clear();
            ref_bp1_better_reads.clear();
            ref_bp2_better_reads.clear();
            dup->sample_info.alt_ref_equal_reads = 0;
            dup->sample_info.alt_ref_equal_reads_highmq = 0;
            evidence_map->clear_other_read_support_for_too_deep(dup->sample_info);
            break;
        }
    }

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;

    std::string alt_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bool> alt_is_exact_read, ref_bp1_is_exact_read, ref_bp2_is_exact_read;
    auto alt_is_consistent_read = gen_consensus_and_classify_seqs(alt_seq, alt_better_reads, std::vector<bool>(), alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read, hp_mismatch_rate_thresholds);
    auto ref_bp1_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp1_seq, ref_bp1_better_reads, std::vector<bool>(), ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score, ref_bp1_is_exact_read, hp_mismatch_rate_thresholds);
    auto ref_bp2_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp2_seq, ref_bp2_better_reads, std::vector<bool>(), ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score, ref_bp2_is_exact_read, hp_mismatch_rate_thresholds);

    std::vector<int> alt_better_read_positions_consistent = get_consistent_reads_start_positions(alt_is_consistent_read, alt_better_read_positions);
    dup->sample_info.alt1_occ_ratio = occ_ratio(alt_better_read_positions_consistent, alt_ref_diff_reads_expected_positions.size());

    for (int i = 0; i < alt_better_reads.size(); i++) {
        std::shared_ptr<bam1_t>& r = alt_better_reads[i];
        evidence_map->record_assigned_alt_read(dup, r.get(), alt_is_consistent_read[i], get_mq(r.get()) >= config.high_confidence_mapq, alt_is_exact_read[i]);
    }

    auto score_dup_consensus = [&](const std::string& consensus_seq) {
        alignment_targets_t targets;
        hts_pos_t lh_start = std::max(hts_pos_t(0), dup->end-hts_pos_t(consensus_seq.length()));
        hts_pos_t lh_len = dup->end-lh_start;
        hts_pos_t rh_start = dup->start;
        hts_pos_t rh_end = std::min(dup->start+hts_pos_t(consensus_seq.length()), contig_len);
        hts_pos_t rh_len = rh_end-dup->start;
        targets.alt_len = lh_len+dup->ins_seq.length()+rh_len;
        targets.alt_seq = new char[targets.alt_len+1];
        strncpy(targets.alt_seq, contig_seq+lh_start, lh_len);
        strncpy(targets.alt_seq+lh_len, dup->ins_seq.c_str(), dup->ins_seq.length());
        strncpy(targets.alt_seq+lh_len+dup->ins_seq.length(), contig_seq+rh_start, rh_len);
        targets.alt_seq[targets.alt_len] = 0;
        hts_pos_t ref_bp1_start = std::max(hts_pos_t(0), dup->start-hts_pos_t(consensus_seq.length()));
        hts_pos_t ref_bp1_end = std::min(dup->start+hts_pos_t(consensus_seq.length()), contig_len);
        hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), dup->end-hts_pos_t(consensus_seq.length()));
        hts_pos_t ref_bp2_end = std::min(dup->end+hts_pos_t(consensus_seq.length()), contig_len);
        targets.ref_seqs.push_back(contig_seq+ref_bp1_start);
        targets.ref_lens.push_back(ref_bp1_end-ref_bp1_start);
        targets.ref_seqs.push_back(contig_seq+ref_bp2_start);
        targets.ref_lens.push_back(ref_bp2_end-ref_bp2_start);
        targets.left_flank_end = lh_len;
        targets.right_flank_start = lh_len+dup->ins_seq.length();
        targets.right_flank_end_offset = 0;
        targets.left_independent_ref_seq = contig_seq+ref_bp2_start;
        targets.left_independent_ref_len = ref_bp2_end-ref_bp2_start;
        targets.right_independent_ref_seq = contig_seq+ref_bp1_start;
        targets.right_independent_ref_len = ref_bp1_end-ref_bp1_start;
        consensus_alignment_metrics_t metrics = score_consensus_alignment(consensus_seq, targets, aligner);
        delete[] targets.alt_seq;
        return metrics;
    };

    if (alt_consensus_seq.length() >= 2*config.min_clip_len) {
        consensus_alignment_metrics_t unextended_metrics = score_dup_consensus(alt_consensus_seq);
        dup->sample_info.alt_consensus1_metrics = unextended_metrics;

       // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, std::max<hts_pos_t>(0, dup->end-GENOTYPE_CONSENSUS_EXTENSION), dup->end, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, dup->start, std::min<hts_pos_t>(contig_len, dup->start+GENOTYPE_CONSENSUS_EXTENSION), contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        dup->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        dup->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        dup->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        dup->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;

        consensus_alignment_metrics_t extended_metrics = score_dup_consensus(alt_consensus_seq);
        dup->sample_info.ext_alt_consensus1_metrics = extended_metrics;

        int lf_aln_rlen = extended_metrics.split_ref_lengths[0];
        int rf_aln_rlen = extended_metrics.split_ref_lengths[1];
        dup->left_anchor_aln->start = dup->end-lf_aln_rlen;
        dup->left_anchor_aln->end = dup->end;
        dup->left_anchor_aln->seq_len = lf_aln_rlen;
        dup->right_anchor_aln->start = dup->start;
        dup->right_anchor_aln->end = dup->start+rf_aln_rlen;
        dup->right_anchor_aln->seq_len = rf_aln_rlen;
    }

    delete[] ref_bp1_seq;
    delete[] ref_bp2_seq;

    set_bp_consensus_info(dup->sample_info.alt_bp1.reads_info, alt_better_reads, alt_is_consistent_read, alt_is_exact_read, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp1.reads_info, ref_bp1_better_reads, ref_bp1_is_consistent_read, ref_bp1_is_exact_read, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp2.reads_info, ref_bp2_better_reads, ref_bp2_is_consistent_read, ref_bp2_is_exact_read, ref_bp2_avg_score, ref_bp2_stddev_score);

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

inline void genotype_dups(int id, std::string contig_name, char* contig_seq, hts_pos_t contig_len, std::vector<duplication_t*> dups,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t& stats, config_t& config, contig_map_t& contig_map,
    bam_pool_t* bam_pool, std::string workdir, std::vector<double>* global_crossing_isize_dist, const hp_mismatch_rate_thresholds_t* hp_mismatch_rate_thresholds) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

    int contig_id = contig_map.get_id(contig_name);
    auto chromosome_data = acquire_chromosome_data(contig_id);
    auto mateseqs_w_mapq_chr = chromosome_data.first;
    evidence_map_t* evidence_map = chromosome_data.second;

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (duplication_t* dup : dups) {
        target_ivals.push_back({std::max<hts_pos_t>(0, dup->start-GENOTYPE_CONSENSUS_EXTENSION), std::min<hts_pos_t>(contig_len, dup->start+GENOTYPE_CONSENSUS_EXTENSION)});
        target_ivals.push_back({std::max<hts_pos_t>(0, dup->end-GENOTYPE_CONSENSUS_EXTENSION), std::min<hts_pos_t>(contig_len, dup->end+GENOTYPE_CONSENSUS_EXTENSION)});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_file, candidate_reads_for_extension);

    std::vector<sv_t*> small_dups;
    for (duplication_t* dup : dups) {
        if (dup->svlen() <= stats.read_len-2*config.min_clip_len) {
			genotype_small_dup(dup, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_map, hp_mismatch_rate_thresholds);
            small_dups.push_back(dup);
		} else {
			genotype_large_dup(dup, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_map, hp_mismatch_rate_thresholds);
		}
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    depth_filter_dup(contig_name, dups, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, small_dups, bam_file, config, stats);
    std::string mates_nms_file = workdir + "/workspace/outward-pairs/" + std::to_string(contig_id) + ".txt";
    calculate_ptn_ratio(contig_name, dups, bam_file, config, stats, contig_len, false, evidence_map, mates_nms_file);
    count_stray_pairs(contig_name, dups, bam_file, config, stats);
    release_chromosome_data(contig_id);
}
#endif // GENOTYPE_DUPS_H
