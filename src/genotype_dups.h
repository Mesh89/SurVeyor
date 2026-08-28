#ifndef GENOTYPE_DUPS_H
#define GENOTYPE_DUPS_H

#include "stat_tests.h"

#include "genotype.h"

void genotype_small_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t* evidence_logger,
                bool reassign_evidence, evidence_map_t* evidence_map) {

	hts_pos_t dup_start = dup->start, dup_end = dup->end;
    
	hts_pos_t extend = stats.read_len + 20;
	hts_pos_t svlen = dup->svlen();
    
	// See comments for relative code in genotype_del
	dup_start++; dup_end++;
    if (dup_end > contig_len) {
        dup_end = contig_len;
    }

    hts_pos_t ref_start = std::max(hts_pos_t(0), dup_start-extend), ref_end = std::min(dup_end+extend, contig_len);
	hts_pos_t ref_len = ref_end - ref_start;
	char* ref_seq = new char[ref_len + 1];
	strncpy(ref_seq, contig_seq+ref_start, ref_len);
	ref_seq[ref_len] = 0;

    std::vector<char*> alt_seqs;
	for (int copies = 1; copies*svlen < stats.read_len; copies++) {
		int alt_len = ref_len + copies*svlen;

		char* alt_seq = new char[alt_len+1];
		int pos = 0;
		strncpy(alt_seq, contig_seq+ref_start, dup_end-ref_start);
		pos += dup_end - ref_start;
		for (int i = 0; i < copies; i++) {
            strncpy(alt_seq+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
			strncpy(alt_seq+pos, contig_seq+dup_start, dup_end-dup_start);
			pos += dup_end-dup_start;
		}
		strncpy(alt_seq+pos, contig_seq+dup_end, ref_end-dup_end);
		pos += ref_end - dup_end;
		alt_seq[pos] = 0;
		alt_seqs.push_back(alt_seq);
	}

    std::stringstream region;
    region << dup->chr << ":" << ref_start << "-" << ref_end;

	hts_itr_t* iter = sam_itr_querys(bam_file->idx, bam_file->header, region.str().c_str());
	
    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> ref_better_reads;
    std::vector<std::string> er_reads;
    std::vector<std::vector<std::shared_ptr<bam1_t>>> alt_better_reads(alt_seqs.size());
    std::vector<std::vector<int>> alt_better_read_positions(alt_seqs.size());
    std::vector<std::vector<int>> alt_better_reads_scores(alt_seqs.size());

    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Filter filter_with_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Alignment alt_aln, ref_aln;
    static const std::vector<evidence_map_t::read_alt_association_t> no_cached_alt_associations;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;
        if (!is_samechr(read) || is_samestr(read)) continue;

        std::string seq;
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < dup_start-stats.max_is || read->core.mpos > dup_end) continue;
            seq = get_sequence(read, true);
            rc(seq);
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > dup_end+stats.max_is || mate_endpos < dup_start) continue;
            seq = get_sequence(read, true);
        }

        const std::vector<evidence_map_t::read_alt_association_t>& cached_alt_associations = reassign_evidence ? evidence_map->get_read_alt_associations(read) : no_cached_alt_associations;
        bool cached_alt_read = false;
        for (const auto& association : cached_alt_associations) {
            if (association.sv == dup && association.bp == 1 && association.alt_idx >= 0) cached_alt_read = true;
        }
        bool cached_ref_read = reassign_evidence && evidence_map->is_read_ref(read, dup, 1);
        bool cached_er_read = reassign_evidence && evidence_map->is_read_er(read, dup);
        if (cached_alt_read) {
            if (evidence_map->is_read_assigned_to_different_sv(read, dup)) {
                continue;
            }
            for (const auto& association : cached_alt_associations) {
                if (association.sv != dup || association.bp != 1 || association.alt_idx < 0) continue;
                int cached_alt_idx = association.alt_idx;
                if (cached_alt_idx < 0 || cached_alt_idx >= alt_seqs.size()) throw std::runtime_error("Invalid cached ALT index for small duplication " + dup->id + ".");
                alt_better_reads[cached_alt_idx].push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                alt_better_reads_scores[cached_alt_idx].push_back(0);
                alt_better_read_positions[cached_alt_idx].push_back(association.pos);
            }
            continue;
        }
        if (cached_ref_read) {
            ref_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            continue;
        }
        if (cached_er_read) {
            dup->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) dup->sample_info.alt_ref_equal_reads_highmq++;
            continue;
        }

        aligner.Align(seq.c_str(), ref_seq, ref_len, filter_with_score_only, &ref_aln, 0);

        uint16_t best_aln_score = 0;
        std::vector<uint16_t> alt_aln_scores(alt_seqs.size());
        std::vector<int> alt_aln_start_positions(alt_seqs.size());
        for (int i = 0; i < alt_seqs.size(); i++) {
            aligner.Align(seq.c_str(), alt_seqs[i], strlen(alt_seqs[i]), filter_with_pos, &alt_aln, 0);
            alt_aln_scores[i] = alt_aln.sw_score;
            alt_aln_start_positions[i] = alt_aln.ref_begin;
            if (alt_aln.sw_score > best_aln_score) {
                best_aln_score = alt_aln.sw_score;
            }
        }
        
        bool add_alt_better_read = best_aln_score > ref_aln.sw_score;
        bool add_ref_better_read = best_aln_score < ref_aln.sw_score;

        // OAR/ORR are propagated by the read owner; assigned-away ALT reads remain excluded from AR and REF reads continue into RR.
        if ((add_alt_better_read || add_ref_better_read) && reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, dup)) {
            if (add_alt_better_read) continue;
        }

        if (add_alt_better_read) {
            for (int i = 0; i < alt_seqs.size(); i++) {
                if (alt_aln_scores[i] == best_aln_score) {
                    alt_better_reads[i].push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                    alt_better_reads_scores[i].push_back(alt_aln_scores[i]);
                    alt_better_read_positions[i].push_back(alt_aln_start_positions[i]);
                }
            }
        } else if (add_ref_better_read) {
            ref_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        } else {
            dup->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) {
                dup->sample_info.alt_ref_equal_reads_highmq++;
            }
            if (evidence_logger) er_reads.push_back(read_name_with_suffix(read));
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
        er_reads.clear();
        dup->sample_info.alt_ref_equal_reads = 0;
        dup->sample_info.alt_ref_equal_reads_highmq = 0;
        evidence_map->clear_other_read_support_for_too_deep(dup->sample_info);
    }

    if (evidence_logger && !dup->sample_info.too_deep) {
        for (int i = 0; i < alt_better_reads.size(); i++) evidence_logger->log_reads_associations(dup->chr, dup->id, 1, alt_better_reads[i], alt_better_reads_scores[i], alt_better_read_positions[i], i);
        evidence_logger->log_ref_reads_associations(dup->chr, dup->id, 1, ref_better_reads);
        evidence_logger->log_er_reads_associations(dup->chr, dup->id, er_reads);
    }

    std::vector<char*> ref_seqs = {ref_seq};
    std::vector<hts_pos_t> ref_lens = {ref_len};
    int alt_len = strlen(alt_seqs[alt_with_most_reads]);
    std::vector<hts_pos_t> alt_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_seqs[alt_with_most_reads], alt_len, stats.read_len);
    dup->sample_info.expected_alt1_reads_frac = (double) alt_ref_diff_reads_expected_positions.size() / std::max(1, alt_len - stats.read_len + 1);
    dup->sample_info.max_feasible_alt1_len = get_max_feasible_alt_len(alt_ref_diff_reads_expected_positions, stats.read_len);

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    double alt_stddev_score, ref_stddev_score;
    std::vector<bool> alt_is_exact_read, ref_is_exact_read;
    auto alt_is_consistent_read = gen_consensus_and_classify_seqs(alt_seqs[alt_with_most_reads], alt_better_reads[alt_with_most_reads], std::vector<bool>(), alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read);
    auto ref_is_consistent_read = gen_consensus_and_classify_seqs(ref_seq, ref_better_reads, std::vector<bool>(), ref_consensus_seq, ref_avg_score, ref_stddev_score, ref_is_exact_read);

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

    if (reassign_evidence) {
        for (int i = 0; i < alt_better_reads[alt_with_most_reads].size(); i++) {
            std::shared_ptr<bam1_t>& r = alt_better_reads[alt_with_most_reads][i];
            evidence_map->record_assigned_alt_read(dup, r.get(), alt_is_consistent_read[i], get_mq(r.get()) >= config.high_confidence_mapq, alt_is_exact_read[i]);
        }
    }

    set_bp_consensus_info(dup->sample_info.alt_bp1.reads_info, alt_better_reads[alt_with_most_reads], alt_is_consistent_read, alt_is_exact_read, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp1.reads_info, ref_better_reads, ref_is_consistent_read, ref_is_exact_read, ref_avg_score, ref_stddev_score);

    delete[] ref_seq;
    for (char* alt_seq : alt_seqs) {
        delete[] alt_seq;
    }

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_large_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t* evidence_logger,
                bool reassign_evidence, evidence_map_t* evidence_map) {

    hts_pos_t dup_start = dup->start, dup_end = dup->end;

	hts_pos_t extend = stats.read_len + 20;

	// See comments for relative code in genotype_del
	dup_start++; dup_end++;
    if (dup_end > contig_len) {
        dup_end = contig_len;
    }

	// all ranges will be start-inclusive and end-exclusive, i.e. [a,b)

	hts_pos_t ref_bp1_start = std::max(hts_pos_t(0), dup_start-extend), ref_bp1_end = std::min(dup_start+extend, contig_len);
    hts_pos_t ref_bp1_pos = dup_start - ref_bp1_start;
	hts_pos_t ref_bp1_len = ref_bp1_end - ref_bp1_start;
	hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), dup_end-extend), ref_bp2_end = std::min(dup_end+extend, contig_len);
	hts_pos_t ref_bp2_pos = dup_end - ref_bp2_start;
	hts_pos_t ref_bp2_len = ref_bp2_end - ref_bp2_start;

	// build alt allele
	hts_pos_t alt_lh_len = dup_end - ref_bp2_start;
	hts_pos_t alt_rh_len = ref_bp1_end - dup_start;
	hts_pos_t alt_len = alt_lh_len + dup->ins_seq.length() + alt_rh_len;
	char* alt_seq = new char[alt_len + 1];
	strncpy(alt_seq, contig_seq+ref_bp2_start, alt_lh_len);
    strncpy(alt_seq+alt_lh_len, dup->ins_seq.c_str(), dup->ins_seq.length());
    strncpy(alt_seq+alt_lh_len+dup->ins_seq.length(), contig_seq+dup_start, alt_rh_len);
	alt_seq[alt_len] = 0;
    to_uppercase(alt_seq);

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
    std::vector<std::string> er_reads;
    std::vector<int> alt_better_read_positions;
    std::vector<int> alt_better_reads_scores;
    std::vector<bool> alt_better_read_spans_bp1, alt_better_read_spans_bp2;

    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;

        std::string seq;
        
        if (!is_samechr(read) || is_samestr(read)) continue;
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < dup_start-stats.max_is || read->core.mpos > dup_end) continue;
            seq = get_sequence(read, true);
            rc(seq);
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > dup_end+stats.max_is || mate_endpos < dup_start) continue;
            seq = get_sequence(read, true);
        }

        int cached_alt_bp1_pos = reassign_evidence ? evidence_map->get_read_alt_pos(read, dup, 1) : -1;
        int cached_alt_bp2_pos = reassign_evidence ? evidence_map->get_read_alt_pos(read, dup, 2) : -1;
        if (cached_alt_bp1_pos >= 0 && cached_alt_bp2_pos >= 0 && cached_alt_bp1_pos != cached_alt_bp2_pos) throw std::runtime_error("Different cached ALT alignment positions for duplication " + dup->id + " and read " + read_name_with_suffix(read) + ".");
        int cached_alt_pos = cached_alt_bp1_pos >= 0 ? cached_alt_bp1_pos : cached_alt_bp2_pos;
        bool cached_alt_read = cached_alt_pos >= 0;
        bool cached_ref_bp1_read = reassign_evidence && evidence_map->is_read_ref(read, dup, 1);
        bool cached_ref_bp2_read = reassign_evidence && evidence_map->is_read_ref(read, dup, 2);
        bool cached_ref_read = cached_ref_bp1_read || cached_ref_bp2_read;
        bool cached_er_read = reassign_evidence && evidence_map->is_read_er(read, dup);

        uint16_t ref_aln_score = 0;
        bool increase_ref_bp1_better = false, increase_ref_bp2_better = false;
        hts_pos_t alt_right_flank_pos = alt_lh_len + dup->ins_seq.length();
        bool alt_spans_bp1 = cached_alt_bp1_pos >= 0;
        bool alt_spans_bp2 = cached_alt_bp2_pos >= 0;
        bool add_alt_better_read = cached_alt_read;
        bool add_ref_bp1_better_read = cached_ref_bp1_read;
        bool add_ref_bp2_better_read = cached_ref_bp2_read;
        bool ref_better_read = cached_ref_read;
        if (cached_er_read) {
            dup->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) dup->sample_info.alt_ref_equal_reads_highmq++;
        } else if (!cached_alt_read && !cached_ref_read) {
            // align to REF (two breakpoints)
            bool ref_is_exact_match = is_perfectly_aligned(read);
            if (ref_is_exact_match) {
                ref_aln_score = read->core.l_qseq;
                if (read->core.pos < dup_start && bam_endpos(read) > dup_start) increase_ref_bp1_better = true;
                if (read->core.pos < dup_end && bam_endpos(read) > dup_end) increase_ref_bp2_better = true;
            } else {
                aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter_with_pos, &ref1_aln, 0);
                aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter_with_pos, &ref2_aln, 0);
                ref_aln_score = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln.sw_score : ref2_aln.sw_score;
                if (ref1_aln.sw_score >= ref2_aln.sw_score && ref1_aln.ref_begin <= ref_bp1_pos && ref1_aln.ref_end >= ref_bp1_pos) increase_ref_bp1_better = true;
                if (ref2_aln.sw_score >= ref1_aln.sw_score && ref2_aln.ref_begin <= ref_bp2_pos && ref2_aln.ref_end >= ref_bp2_pos) increase_ref_bp2_better = true;
            }

            // align to ALT
            alt_aln = align_fast(aligner, seq.c_str(), alt_seq, alt_len, filter_with_pos, ref_is_exact_match);
            alt_spans_bp1 = alt_aln.ref_begin < alt_lh_len && alt_aln.ref_end >= alt_lh_len;
            alt_spans_bp2 = alt_aln.ref_begin < alt_right_flank_pos && alt_aln.ref_end >= alt_right_flank_pos;
            add_alt_better_read = alt_aln.sw_score > ref_aln_score;
            ref_better_read = alt_aln.sw_score < ref_aln_score;
            add_ref_bp1_better_read = ref_better_read && increase_ref_bp1_better;
            add_ref_bp2_better_read = ref_better_read && increase_ref_bp2_better;
        }

        // OAR/ORR are propagated by the read owner; assigned-away ALT reads remain excluded from AR and REF reads continue into RR.
        if ((add_alt_better_read || add_ref_bp1_better_read || add_ref_bp2_better_read) && reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, dup)) {
            if (add_alt_better_read) continue;
        }

        if (add_alt_better_read) {
            alt_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_better_read_positions.push_back(cached_alt_read ? cached_alt_pos : alt_aln.ref_begin);
            alt_better_reads_scores.push_back(cached_alt_read ? 0 : alt_aln.sw_score);
            alt_better_read_spans_bp1.push_back(alt_spans_bp1);
            alt_better_read_spans_bp2.push_back(alt_spans_bp2);
        } else if (ref_better_read) {
            if (add_ref_bp1_better_read) {
                ref_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            }
            if (add_ref_bp2_better_read) {
                ref_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            }
        } else if (!cached_er_read) {
            dup->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) {
                dup->sample_info.alt_ref_equal_reads_highmq++;
            }
            if (evidence_logger) er_reads.push_back(read_name_with_suffix(read));
        }

        if (ref_bp1_better_reads.size() + ref_bp2_better_reads.size() + dup->sample_info.alt_ref_equal_reads > 4*stats.get_max_depth(dup->chr)) {
            alt_better_reads.clear();
            alt_better_read_positions.clear();
            alt_better_reads_scores.clear();
            alt_better_read_spans_bp1.clear();
            alt_better_read_spans_bp2.clear();
            ref_bp1_better_reads.clear();
            ref_bp2_better_reads.clear();
            er_reads.clear();
            dup->sample_info.alt_ref_equal_reads = 0;
            dup->sample_info.alt_ref_equal_reads_highmq = 0;
            evidence_map->clear_other_read_support_for_too_deep(dup->sample_info);
            break;
        }
    }

    if (evidence_logger) {
        std::vector<std::shared_ptr<bam1_t>> alt_bp0_reads, alt_bp1_reads, alt_bp2_reads;
        std::vector<int> alt_bp0_scores, alt_bp1_scores, alt_bp2_scores, alt_bp0_positions, alt_bp1_positions, alt_bp2_positions;
        for (size_t i = 0; i < alt_better_reads.size(); i++) {
            if (!alt_better_read_spans_bp1[i] && !alt_better_read_spans_bp2[i]) {
                alt_bp0_reads.push_back(alt_better_reads[i]);
                alt_bp0_scores.push_back(alt_better_reads_scores[i]);
                alt_bp0_positions.push_back(alt_better_read_positions[i]);
            }
            if (alt_better_read_spans_bp1[i]) {
                alt_bp1_reads.push_back(alt_better_reads[i]);
                alt_bp1_scores.push_back(alt_better_reads_scores[i]);
                alt_bp1_positions.push_back(alt_better_read_positions[i]);
            }
            if (alt_better_read_spans_bp2[i]) {
                alt_bp2_reads.push_back(alt_better_reads[i]);
                alt_bp2_scores.push_back(alt_better_reads_scores[i]);
                alt_bp2_positions.push_back(alt_better_read_positions[i]);
            }
        }
        evidence_logger->log_reads_associations(dup->chr, dup->id, 0, alt_bp0_reads, alt_bp0_scores, alt_bp0_positions);
        evidence_logger->log_reads_associations(dup->chr, dup->id, 1, alt_bp1_reads, alt_bp1_scores, alt_bp1_positions);
        evidence_logger->log_reads_associations(dup->chr, dup->id, 2, alt_bp2_reads, alt_bp2_scores, alt_bp2_positions);
        evidence_logger->log_ref_reads_associations(dup->chr, dup->id, 1, ref_bp1_better_reads);
        evidence_logger->log_ref_reads_associations(dup->chr, dup->id, 2, ref_bp2_better_reads);
        evidence_logger->log_er_reads_associations(dup->chr, dup->id, er_reads);
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
    auto alt_is_consistent_read = gen_consensus_and_classify_seqs(alt_seq, alt_better_reads, std::vector<bool>(), alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read);
    auto ref_bp1_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp1_seq, ref_bp1_better_reads, std::vector<bool>(), ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score, ref_bp1_is_exact_read);
    auto ref_bp2_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp2_seq, ref_bp2_better_reads, std::vector<bool>(), ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score, ref_bp2_is_exact_read);

    std::vector<int> alt_better_read_positions_consistent = get_consistent_reads_start_positions(alt_is_consistent_read, alt_better_read_positions);
    dup->sample_info.alt1_occ_ratio = occ_ratio(alt_better_read_positions_consistent, alt_ref_diff_reads_expected_positions.size());

    if (reassign_evidence) {
        for (int i = 0; i < alt_better_reads.size(); i++) {
            std::shared_ptr<bam1_t>& r = alt_better_reads[i];
            evidence_map->record_assigned_alt_read(dup, r.get(), alt_is_consistent_read[i], get_mq(r.get()) >= config.high_confidence_mapq, alt_is_exact_read[i]);
        }
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

    delete[] alt_seq;

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_dups(int id, std::string contig_name, char* contig_seq, hts_pos_t contig_len, std::vector<duplication_t*> dups,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t& stats, config_t& config, contig_map_t& contig_map,
    bam_pool_t* bam_pool, std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq_chr,
    std::string workdir, std::vector<double>* global_crossing_isize_dist, evidence_logger_t* evidence_logger,
    bool reassign_evidence, evidence_map_t* evidence_map) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

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
			genotype_small_dup(dup, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_logger, reassign_evidence, evidence_map);
            small_dups.push_back(dup);
		} else {
			genotype_large_dup(dup, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_logger, reassign_evidence, evidence_map);
		}
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    release_mates(contig_id);

    depth_filter_dup(contig_name, dups, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, small_dups, bam_file, config, stats);
    std::string mates_nms_file = workdir + "/workspace/outward-pairs/" + std::to_string(contig_id) + ".txt";
    calculate_ptn_ratio(contig_name, dups, bam_file, config, stats, contig_len, evidence_logger, false, evidence_map, mates_nms_file);
    count_stray_pairs(contig_name, dups, bam_file, config, stats);
}

#endif // GENOTYPE_DUPS_H
