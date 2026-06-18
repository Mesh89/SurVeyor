#ifndef GENOTYPE_DUPS_H
#define GENOTYPE_DUPS_H

#include "stat_tests.h"

#include "genotype.h"

void genotype_small_dup(duplication_t* dup, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t* evidence_logger,
                bool reassign_evidence, evidence_map_t* evidence_map, std::unordered_map<std::string, std::shared_ptr<sv_t>>& sv_map) {

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
    std::vector<std::vector<std::shared_ptr<bam1_t>>> alt_better_reads(alt_seqs.size());
    std::vector<std::vector<int>> alt_better_read_positions(alt_seqs.size());
    std::vector<std::vector<int>> alt_better_reads_scores(alt_seqs.size());

    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Filter filter_with_pos_and_cigar(true, true, 0, 32767);
    StripedSmithWaterman::Filter filter_with_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Alignment alt_aln, ref_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;
        if (!is_samechr(read) || is_samestr(read)) continue;

        // if the read is assigned to a different SV, no need to align it, just count and continue
        if (reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, dup)) {
            dup->sample_info.assigned_to_other_sv_bp1_reads++;
            continue;
        }

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
        
        if (best_aln_score > ref_aln.sw_score) {
            for (int i = 0; i < alt_seqs.size(); i++) {
                if (alt_aln_scores[i] == best_aln_score) {
                    alt_better_reads[i].push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                    alt_better_reads_scores[i].push_back(alt_aln_scores[i]);
                    alt_better_read_positions[i].push_back(alt_aln_start_positions[i]);
                }
            }
        } else if (best_aln_score < ref_aln.sw_score) {
            ref_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        } else {
            dup->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) {
                dup->sample_info.alt_ref_equal_reads_highmq++;
            }
        }
    }

    int alt_with_most_reads = 0;
    for (int i = 1; i < alt_better_reads.size(); i++) {
        if (alt_better_reads[i].size() > alt_better_reads[alt_with_most_reads].size()) {
            alt_with_most_reads = i;
        }
    }

    if (alt_better_reads[alt_with_most_reads].size() > 20*stats.get_max_depth(dup->chr) || ref_better_reads.size() > 20*stats.get_max_depth(dup->chr)) {
        alt_better_reads[alt_with_most_reads].clear();
        ref_better_reads.clear();
        dup->sample_info.alt_ref_equal_reads = 0;
        dup->sample_info.alt_ref_equal_reads_highmq = 0;
        dup->sample_info.too_deep = true;
    }

    if (evidence_logger) evidence_logger->log_reads_associations(dup->id, 1, alt_better_reads[alt_with_most_reads], alt_better_reads_scores[alt_with_most_reads]);

    std::vector<char*> ref_seqs = {ref_seq};
    std::vector<hts_pos_t> ref_lens = {ref_len};
    int alt_len = strlen(alt_seqs[alt_with_most_reads]);
    std::vector<hts_pos_t> alt_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_seqs[alt_with_most_reads], alt_len, stats.read_len);
    dup->expected_alt1_reads_frac = (double) alt_ref_diff_reads_expected_positions.size() / std::max(1, alt_len - stats.read_len + 1);

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    double alt_stddev_score, ref_stddev_score;
    std::vector<bool> alt_is_exact_read, ref_is_exact_read;
    auto alt_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(alt_seqs[alt_with_most_reads], alt_better_reads[alt_with_most_reads], std::vector<bool>(), alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read);
    auto ref_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_seq, ref_better_reads, std::vector<bool>(), ref_consensus_seq, ref_avg_score, ref_stddev_score, ref_is_exact_read);

    std::vector<int> alt_better_read_positions_consistent = get_consistent_reads_start_positions(alt_better_reads_consistent, alt_better_reads[alt_with_most_reads], alt_better_read_positions[alt_with_most_reads]);
    dup->sample_info.alt1_occ_ratio = occ_ratio(alt_better_read_positions_consistent, alt_ref_diff_reads_expected_positions.size());

    if (alt_consensus_seq.length() >= 2*config.min_clip_len) {
       // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, dup->start-stats.max_is, dup->start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, dup->end, dup->end+stats.max_is, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        dup->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        dup->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        dup->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        dup->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;

        hts_pos_t ref_start = dup->start-alt_consensus_seq.length();
        if (ref_start < 0) ref_start = 0;
        hts_pos_t ref_end = dup->end+alt_consensus_seq.length();
        if (ref_end > contig_len) ref_end = contig_len;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+ref_start, ref_end-ref_start, filter_with_pos_and_cigar, &ref_aln, 0);
        dup->sample_info.ext_alt_consensus1_to_ref_score = ref_aln.query_end - ref_aln.query_begin - ref_aln.mismatches;

        int n_extra_copies = alt_with_most_reads+1;
        int alt_len = ref_end - ref_start + n_extra_copies*svlen;
		char* alt_seq = new char[alt_len+1];
		int pos = 0;
		strncpy(alt_seq, contig_seq+ref_start, dup_end-ref_start);
		pos += dup_end - ref_start;
		for (int i = 0; i < n_extra_copies; i++) {
            strncpy(alt_seq+pos, dup->ins_seq.c_str(), dup->ins_seq.length());
            pos += dup->ins_seq.length();
			strncpy(alt_seq+pos, contig_seq+dup_start, dup_end-dup_start);
			pos += dup_end-dup_start;
		}
		strncpy(alt_seq+pos, contig_seq+dup_end, ref_end-dup_end);
		pos += ref_end - dup_end;
		alt_seq[pos] = 0;
        aligner.Align(alt_consensus_seq.c_str(), alt_seq, alt_len, filter_with_pos_and_cigar, &alt_aln, 0);
        dup->sample_info.ext_alt_consensus1_to_alt_score = alt_aln.query_end - alt_aln.query_begin - alt_aln.mismatches;
        delete[] alt_seq;

        int lf_seq_end = dup_start - ref_start;
        int lf_aln_rlen = std::max(0, lf_seq_end-alt_aln.ref_begin);
        int rf_seq_len = ref_end - dup_end;
        int rf_seq_start = pos - rf_seq_len;
        int rf_aln_rlen = std::max(0, alt_aln.ref_end-rf_seq_start);

        int temp;
        auto query_lh_aln_score = find_aln_prefix_score(alt_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        dup->sample_info.alt_consensus1_split_size1 = query_lh_aln_score.second;
        dup->sample_info.alt_consensus1_split_size2 = query_rh_aln_score.second;
        dup->sample_info.alt_consensus1_split_score1 = query_lh_aln_score.first;
        dup->sample_info.alt_consensus1_split_score2 = query_rh_aln_score.first;

        dup->left_anchor_aln->start = dup_end - lf_aln_rlen;
        dup->left_anchor_aln->end = dup_end;
        dup->left_anchor_aln->seq_len = lf_aln_rlen;
        dup->right_anchor_aln->start = dup_start;
        dup->right_anchor_aln->end = dup_start + rf_aln_rlen;
        dup->right_anchor_aln->seq_len = rf_aln_rlen;

        dup->sample_info.ext_alt_consensus1_length = alt_consensus_seq.length();

        ref_aln.Clear();
        std::string lh_query = alt_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), contig_seq+ref_start, ref_end-ref_start, filter_with_score_only, &ref_aln, 0);
        dup->sample_info.alt_consensus1_split_score1_ind_aln = ref_aln.sw_score;

        ref_aln.Clear();
        std::string rh_query = alt_consensus_seq.substr(alt_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), contig_seq+ref_start, ref_end-ref_start, filter_with_score_only, &ref_aln, 0);
        dup->sample_info.alt_consensus1_split_score2_ind_aln = ref_aln.sw_score;
    }

    alt_better_reads_consistent = find_seqs_consistent_with_ref_seq(alt_consensus_seq, alt_better_reads_consistent, alt_avg_score, alt_stddev_score, alt_is_exact_read);

    if (reassign_evidence) { // increment ORC counters for other SVs that lost support from these reads
        for (int i = 0; i < alt_better_reads_consistent.size(); i++) {
            std::shared_ptr<bam1_t>& r = alt_better_reads_consistent[i];
            for (std::pair<std::string, int>& ov : evidence_map->get_non_chosen_svs_for_read(r.get())) {
                increase_orc(sv_map, ov.first, ov.second, r.get(), get_mq(r.get()) >= config.high_confidence_mapq, alt_is_exact_read[i]);
            }
        }
    }

    set_bp_consensus_info(dup->sample_info.alt_bp1.reads_info, alt_better_reads[alt_with_most_reads].size(), alt_better_reads_consistent, alt_is_exact_read, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp1.reads_info, ref_better_reads.size(), ref_better_reads_consistent, ref_is_exact_read, ref_avg_score, ref_stddev_score);

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
                bool reassign_evidence, evidence_map_t* evidence_map, std::unordered_map<std::string, std::shared_ptr<sv_t>>& sv_map) {

    hts_pos_t dup_start = dup->start, dup_end = dup->end;

	hts_pos_t extend = stats.read_len + 20;

	// See comments for relative code in genotype_del
	dup_start++; dup_end++;

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

    std::vector<char*> ref_seqs = {contig_seq+ref_bp1_start, contig_seq+ref_bp2_start};
    std::vector<hts_pos_t> ref_lens = {ref_bp1_len, ref_bp2_len};
    std::vector<hts_pos_t> alt_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_seq, alt_len, stats.read_len);
    dup->expected_alt1_reads_frac = (double) alt_ref_diff_reads_expected_positions.size() / std::max(hts_pos_t(1), alt_len - stats.read_len + 1);

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

    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Filter filter_with_pos_and_cigar(true, true, 0, 32767);
    StripedSmithWaterman::Filter filter_with_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < dup_start || dup_end < get_unclipped_start(read)) continue;
        if (dup_start < get_unclipped_start(read) && get_unclipped_end(read) < dup_end) continue;

        // if the read is assigned to a different SV, no need to align it, just count and continue
        if (reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, dup)) {
            dup->sample_info.assigned_to_other_sv_bp1_reads++;
            dup->sample_info.assigned_to_other_sv_bp2_reads++;
            continue;
        }

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

        // align to REF (two breakpoints)
        uint16_t ref_aln_score = 0;
        bool increase_ref_bp1_better = false, increase_ref_bp2_better = false;
        if (is_perfectly_aligned(read)) {
            ref_aln_score = read->core.l_qseq;
            if (read->core.pos < dup_start && bam_endpos(read) > dup_start) {
                increase_ref_bp1_better = true;
            }
            if (read->core.pos < dup_end && bam_endpos(read) > dup_end) {
                increase_ref_bp2_better = true;
            }
        } else {
            aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter_with_pos, &ref1_aln, 0);
            aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter_with_pos, &ref2_aln, 0);
            ref_aln_score = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln.sw_score : ref2_aln.sw_score;
            if (ref1_aln.sw_score >= ref2_aln.sw_score && ref1_aln.ref_begin <= ref_bp1_pos && ref1_aln.ref_end >= ref_bp1_pos) {
                increase_ref_bp1_better = true;
            }
            if (ref2_aln.sw_score >= ref1_aln.sw_score && ref2_aln.ref_begin <= ref_bp2_pos && ref2_aln.ref_end >= ref_bp2_pos) {
                increase_ref_bp2_better = true;
            }
        }

        // align to ALT
        aligner.Align(seq.c_str(), alt_seq, alt_len, filter_with_pos, &alt_aln, 0);

        if (alt_aln.sw_score > ref_aln_score) {
            alt_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_better_read_positions.push_back(alt_aln.ref_begin);
            alt_better_reads_scores.push_back(alt_aln.sw_score);
        } else if (alt_aln.sw_score < ref_aln_score) {
            if (increase_ref_bp1_better) {
                ref_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            }
            if (increase_ref_bp2_better) {
                ref_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            }
        } else {
            dup->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) {
                dup->sample_info.alt_ref_equal_reads_highmq++;
            }
        }

        if (ref_bp1_better_reads.size() + ref_bp2_better_reads.size() + dup->sample_info.alt_ref_equal_reads > 4*stats.get_max_depth(dup->chr)) {
            alt_better_reads.clear();
            ref_bp1_better_reads.clear();
            ref_bp2_better_reads.clear();
            dup->sample_info.alt_ref_equal_reads = 0;
            dup->sample_info.alt_ref_equal_reads_highmq = 0;
            dup->sample_info.too_deep = true;
            break;
        }
    }

    if (evidence_logger) evidence_logger->log_reads_associations(dup->id, 1, alt_better_reads, alt_better_reads_scores);

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
    auto alt_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(alt_seq, alt_better_reads, std::vector<bool>(), alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read);
    auto ref_bp1_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_reads, std::vector<bool>(), ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score, ref_bp1_is_exact_read);
    auto ref_bp2_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_reads, std::vector<bool>(), ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score, ref_bp2_is_exact_read);
    delete[] ref_bp1_seq;
    delete[] ref_bp2_seq;

    std::vector<int> alt_better_read_positions_consistent = get_consistent_reads_start_positions(alt_better_reads_consistent, alt_better_reads, alt_better_read_positions);
    dup->sample_info.alt1_occ_ratio = occ_ratio(alt_better_read_positions_consistent, alt_ref_diff_reads_expected_positions.size());

    if (reassign_evidence) { // increment ORC counters for other SVs that lost support from these reads
        for (int i = 0; i < alt_better_reads_consistent.size(); i++) {
            std::shared_ptr<bam1_t>& r = alt_better_reads_consistent[i];
            for (std::pair<std::string, int>& ov : evidence_map->get_non_chosen_svs_for_read(r.get())) {
                increase_orc(sv_map, ov.first, ov.second, r.get(), get_mq(r.get()) >= config.high_confidence_mapq, alt_is_exact_read[i]);
            }
        }
    }

    if (alt_consensus_seq.length() >= 2*config.min_clip_len) {
       // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, dup->end-stats.max_is, dup->end, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, dup->start, dup->start+stats.max_is, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        dup->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        dup->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        dup->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        dup->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;

        hts_pos_t lh_start = dup->end-alt_consensus_seq.length();
        if (lh_start < 0) lh_start = 0;
        hts_pos_t lh_len = dup->end-lh_start;
        hts_pos_t rh_start = dup->start;
        hts_pos_t rh_end = dup->start+alt_consensus_seq.length();
        if (rh_end > contig_len) rh_end = contig_len;
        hts_pos_t rh_len = rh_end-dup->start;
    
        delete[] alt_seq;
        alt_seq = new char[lh_len + dup->ins_seq.length() + rh_len + 1];
        strncpy(alt_seq, contig_seq+lh_start, lh_len);
        strncpy(alt_seq+lh_len, dup->ins_seq.c_str(), dup->ins_seq.length());
        strncpy(alt_seq+lh_len+dup->ins_seq.length(), contig_seq+rh_start, rh_len);
        alt_seq[lh_len+dup->ins_seq.length()+rh_len] = 0;

        // align to ref+SV
        aligner.Align(alt_consensus_seq.c_str(), alt_seq, lh_len+dup->ins_seq.length()+rh_len, filter_with_pos_and_cigar, &alt_aln, 0);

        // length of the left and right flanking regions of the deletion covered by alt_consensus_seq
        int lf_aln_rlen = std::max(hts_pos_t(0), lh_len - alt_aln.ref_begin);
        int rf_aln_rlen = std::max(hts_pos_t(0), alt_aln.ref_end - lh_len - (int) dup->ins_seq.length());

        // length of the alt_consensus_seq covering left and right flanking regions of the deletion
        // note that this may be different from lf_aln_rlen and rf_aln_rlen, since the aln can include indels
        int temp;
        auto query_lh_aln_score = find_aln_prefix_score(alt_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        dup->sample_info.alt_consensus1_split_size1 = query_lh_aln_score.second;
        dup->sample_info.alt_consensus1_split_size2 = query_rh_aln_score.second;
        dup->sample_info.alt_consensus1_split_score1 = query_lh_aln_score.first;
        dup->sample_info.alt_consensus1_split_score2 = query_rh_aln_score.first;

        dup->left_anchor_aln->start = dup->end - lf_aln_rlen;
        dup->left_anchor_aln->end = dup->end;
        dup->left_anchor_aln->seq_len = lf_aln_rlen;
        dup->right_anchor_aln->start = dup->start;
        dup->right_anchor_aln->end = dup->start + rf_aln_rlen;
        dup->right_anchor_aln->seq_len = rf_aln_rlen;

        // align to ref
        hts_pos_t ref_bp1_start = dup->start - alt_consensus_seq.length(), ref_bp1_end = dup->start + alt_consensus_seq.length();
        if (ref_bp1_start < 0) ref_bp1_start = 0;
        if (ref_bp1_end > contig_len) ref_bp1_end = contig_len;
        hts_pos_t ref_bp2_start = dup->end - alt_consensus_seq.length(), ref_bp2_end = dup->end + alt_consensus_seq.length();
        if (ref_bp2_start < 0) ref_bp2_start = 0;
        if (ref_bp2_end > contig_len) ref_bp2_end = contig_len;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter_with_pos_and_cigar, &ref1_aln, 0);
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter_with_pos_and_cigar, &ref2_aln, 0);

        dup->sample_info.ext_alt_consensus1_length = alt_consensus_seq.length();
        dup->sample_info.ext_alt_consensus1_to_alt_score = alt_aln.query_end - alt_aln.query_begin - alt_aln.mismatches;
        dup->sample_info.ext_alt_consensus1_to_ref_score = std::max(
            ref1_aln.query_end - ref1_aln.query_begin - ref1_aln.mismatches,
            ref2_aln.query_end - ref2_aln.query_begin - ref2_aln.mismatches
        );

        ref1_aln.Clear();
        std::string lh_query = alt_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter_with_score_only, &ref1_aln, 0);
        dup->sample_info.alt_consensus1_split_score1_ind_aln = ref1_aln.sw_score;

        ref2_aln.Clear();
        std::string rh_query = alt_consensus_seq.substr(alt_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter_with_score_only, &ref2_aln, 0);
        dup->sample_info.alt_consensus1_split_score2_ind_aln = ref2_aln.sw_score;
    }

    set_bp_consensus_info(dup->sample_info.alt_bp1.reads_info, alt_better_reads.size(), alt_better_reads_consistent, alt_is_exact_read, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp1.reads_info, ref_bp1_better_reads.size(), ref_bp1_better_reads_consistent, ref_bp1_is_exact_read, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(dup->sample_info.ref_bp2.reads_info, ref_bp2_better_reads.size(), ref_bp2_better_reads_consistent, ref_bp2_is_exact_read, ref_bp2_avg_score, ref_bp2_stddev_score);

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
    bool reassign_evidence, evidence_map_t* evidence_map, std::unordered_map<std::string, std::shared_ptr<sv_t>>* sv_map) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (duplication_t* dup : dups) {
        target_ivals.push_back({dup->start-stats.max_is, dup->start+stats.max_is});
        target_ivals.push_back({dup->end-stats.max_is, dup->end+stats.max_is});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_file, candidate_reads_for_extension);

    std::vector<sv_t*> small_dups;
    for (duplication_t* dup : dups) {
        if (dup->svlen() <= stats.read_len-2*config.min_clip_len) {
			genotype_small_dup(dup, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_logger, reassign_evidence, evidence_map, *sv_map);
            small_dups.push_back(dup);
		} else {
			genotype_large_dup(dup, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_logger, reassign_evidence, evidence_map, *sv_map);
		}
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    release_mates(contig_id);

    depth_filter_dup(contig_name, dups, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, small_dups, bam_file, config, stats);
    std::string mates_nms_file = workdir + "/workspace/outward-pairs/" + std::to_string(contig_id) + ".txt";
    calculate_ptn_ratio(contig_name, dups, bam_file, config, stats, evidence_logger, false, evidence_map, mates_nms_file);
    count_stray_pairs(contig_name, dups, bam_file, config, stats);
}

#endif // GENOTYPE_DUPS_H
