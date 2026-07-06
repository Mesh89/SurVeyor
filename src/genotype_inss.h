#ifndef GENOTYPE_INSS_H
#define GENOTYPE_INSS_H

#include "stat_tests.h"

#include "genotype.h"

void genotype_ins(insertion_t* ins, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner, evidence_logger_t* evidence_logger,
                bool reassign_evidence, evidence_map_t* evidence_map, 
                std::unordered_map<std::string, std::shared_ptr<sv_t>>& sv_map) {

    hts_pos_t ins_start = ins->start, ins_end = ins->end;

	hts_pos_t extend = stats.read_len + 20;

	// build alt allele
	/*
	 * POS in VCF is the base BEFORE the insertion
	 * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" insertion (no deletion),POS == END, otherwise the last base deleted
	 * As usual, in order to make intervals [ ), we increase the coordinates by 1
	 */
	ins_start++; ins_end++;

	hts_pos_t alt_start = std::max(hts_pos_t(0), ins_start-extend);
	hts_pos_t alt_end = std::min(ins_end+extend, contig_len);
	int alt_lf_len = ins_start-alt_start, alt_rf_len = alt_end-ins_end;
    hts_pos_t ins_lh_len = std::min(extend, (hts_pos_t) ins->ins_seq.length());
	hts_pos_t ins_rh_len = ins_lh_len;

    char* lf_seq = generate_haplotype_left(contig_seq, ins_start-1, extend, ins->aux_indels, ins->aux_snps);
    alt_lf_len = strlen(lf_seq);

    char* rf_seq = generate_haplotype_right(contig_seq, contig_len, ins_end, extend,  ins->aux_indels, ins->aux_snps);
    alt_rf_len = strlen(rf_seq);

    char* alt_bp1_seq;
    int alt_bp1_len;
    if (ins_lh_len < extend) {
        // in this case, alt_bp1 = left flank + insertion + enough right flank to make sure that
        // there are 'extend' bases on each side of the breakpoint (if possible)
        int extra = std::min(extend-ins_lh_len, (hts_pos_t) alt_rf_len); // there may not be enough bp if near the end of contig
        alt_bp1_len = alt_lf_len + ins_lh_len + extra;
        alt_bp1_seq = new char[alt_bp1_len+1];
        strncpy(alt_bp1_seq, lf_seq, alt_lf_len);
        strncpy(alt_bp1_seq+alt_lf_len, ins->ins_seq.c_str(), ins_lh_len);
        strncpy(alt_bp1_seq+alt_lf_len+ins_lh_len, rf_seq, extra);
        alt_bp1_seq[alt_bp1_len] = 0;
    } else {
        alt_bp1_len = alt_lf_len + ins_lh_len;
        alt_bp1_seq = new char[alt_bp1_len+1];
        strncpy(alt_bp1_seq, lf_seq, alt_lf_len);
        strncpy(alt_bp1_seq+alt_lf_len, ins->ins_seq.c_str(), ins_lh_len);
        alt_bp1_seq[alt_bp1_len] = 0;
    }

	char* alt_bp2_seq;
    int alt_bp2_len;
    if (ins_rh_len < extend) {
        // in this case, alt_bp2 = enough left flank + insertion + right flank
        int extra = std::min(extend-ins_rh_len, (hts_pos_t) alt_lf_len); // there may not be enough bp if near the start of contig
        alt_bp2_len = extra + ins_rh_len + alt_rf_len;
        alt_bp2_seq = new char[alt_bp2_len+1];
        strncpy(alt_bp2_seq, lf_seq+(alt_lf_len - extra), extra);
        strncpy(alt_bp2_seq+extra, ins->ins_seq.c_str(), ins_rh_len);
        strncpy(alt_bp2_seq+extra+ins_rh_len, rf_seq, alt_rf_len);
        alt_bp2_seq[alt_bp2_len] = 0;
    } else {
        alt_bp2_len = ins_rh_len + alt_rf_len;
        alt_bp2_seq = new char[alt_bp2_len+1];
        strncpy(alt_bp2_seq, ins->ins_seq.c_str()+(ins->ins_seq.length()-ins_rh_len), ins_rh_len);
        strncpy(alt_bp2_seq+ins_rh_len, rf_seq, alt_rf_len);
        alt_bp2_seq[alt_bp2_len] = 0;
    }

    hts_pos_t ref_bp1_start = alt_start, ref_bp1_end = std::min(ins_start+extend, contig_len);
    hts_pos_t ref_bp1_len = ref_bp1_end - ref_bp1_start;
    hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), ins_end-extend), ref_bp2_end = alt_end;
    hts_pos_t ref_bp2_len = ref_bp2_end - ref_bp2_start;

    hts_pos_t ref_allele_len = ins_end - ins_start;
    hts_pos_t ref_bp1_w_aux_ref_len = std::min(extend, ref_allele_len);
    hts_pos_t ref_bp1_w_aux_rh_len = std::min((hts_pos_t) alt_rf_len, extend - ref_bp1_w_aux_ref_len);
    hts_pos_t ref_bp1_w_aux_len = alt_lf_len + ref_bp1_w_aux_ref_len + ref_bp1_w_aux_rh_len;
    char* ref_bp1_w_aux_seq = new char[ref_bp1_w_aux_len + 1];
    strncpy(ref_bp1_w_aux_seq, lf_seq, alt_lf_len);
    strncpy(ref_bp1_w_aux_seq + alt_lf_len, contig_seq + ins_start, ref_bp1_w_aux_ref_len);
    strncpy(ref_bp1_w_aux_seq + alt_lf_len + ref_bp1_w_aux_ref_len, rf_seq, ref_bp1_w_aux_rh_len);
    ref_bp1_w_aux_seq[ref_bp1_w_aux_len] = 0;

    hts_pos_t ref_bp2_w_aux_ref_len = std::min(extend, ref_allele_len);
    hts_pos_t ref_bp2_w_aux_lh_len = std::min((hts_pos_t) alt_lf_len, extend - ref_bp2_w_aux_ref_len);
    hts_pos_t ref_bp2_w_aux_len = ref_bp2_w_aux_lh_len + ref_bp2_w_aux_ref_len + alt_rf_len;
    char* ref_bp2_w_aux_seq = new char[ref_bp2_w_aux_len + 1];
    strncpy(ref_bp2_w_aux_seq, lf_seq + alt_lf_len - ref_bp2_w_aux_lh_len, ref_bp2_w_aux_lh_len);
    strncpy(ref_bp2_w_aux_seq + ref_bp2_w_aux_lh_len, contig_seq + ins_end - ref_bp2_w_aux_ref_len, ref_bp2_w_aux_ref_len);
    strncpy(ref_bp2_w_aux_seq + ref_bp2_w_aux_lh_len + ref_bp2_w_aux_ref_len, rf_seq, alt_rf_len);
    ref_bp2_w_aux_seq[ref_bp2_w_aux_len] = 0;

    delete[] lf_seq;
    delete[] rf_seq;

    hts_pos_t ref_bp1_w_aux_pos = alt_lf_len;
    hts_pos_t ref_bp2_w_aux_pos = ref_bp2_w_aux_lh_len + ref_bp2_w_aux_ref_len;

    std::vector<char*> ref_seqs = {ref_bp1_w_aux_seq, ref_bp2_w_aux_seq};
    std::vector<hts_pos_t> ref_lens = {ref_bp1_w_aux_len, ref_bp2_w_aux_len};
    std::vector<hts_pos_t> alt1_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_bp1_seq, alt_bp1_len, stats.read_len);
    std::vector<hts_pos_t> alt2_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_bp2_seq, alt_bp2_len, stats.read_len);
    ins->expected_alt1_reads_frac = (double) alt1_ref_diff_reads_expected_positions.size() / std::max(1, alt_bp1_len - stats.read_len + 1);
    ins->expected_alt2_reads_frac = (double) alt2_ref_diff_reads_expected_positions.size() / std::max(1, alt_bp2_len - stats.read_len + 1);

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
    hts_pos_t alt_bp1_pos = alt_lf_len;
    hts_pos_t alt_bp2_pos = alt_bp2_len - alt_rf_len;

    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Filter filter_with_pos_and_cigar(true, true, 0, 32767);
    StripedSmithWaterman::Filter filter_with_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Alignment alt1_aln, alt2_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (ins_start < get_unclipped_start(read) && get_unclipped_end(read) < ins_end) continue;

        // sometimes, reads that support an indel do not overlap the indel
        // in this case, we consider HSR reads that, according to their mate, can align to the an insertion bp
        // we only allow them to support the alt allele
        bool only_allow_alt = false;
        if (get_unclipped_end(read) < ins_start || ins_end < get_unclipped_start(read)) only_allow_alt = true;

        bool can_align_to_alt = can_align_to(read, stats.max_is, read->core.tid, ins_start) || can_align_to(read, stats.max_is, read->core.tid, ins_end);
        if (only_allow_alt && (!can_align_to_alt || !is_hidden_split_read(read, config))) continue;

        std::string seq = get_sequence(read);
        
        // align to ALT
        aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter_with_pos, &alt1_aln, 0);
        aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter_with_pos, &alt2_aln, 0);
        bool alt1_covers_bp1 = alt1_aln.ref_begin <= alt_bp1_pos && alt1_aln.ref_end >= alt_bp1_pos;
        bool alt2_covers_bp2 = alt2_aln.ref_begin <= alt_bp2_pos && alt2_aln.ref_end >= alt_bp2_pos;
        
        // align to REF (two breakpoints)
        aligner.Align(seq.c_str(), ref_bp1_w_aux_seq, ref_bp1_w_aux_len, filter_with_pos, &ref1_aln, 0);
        aligner.Align(seq.c_str(), ref_bp2_w_aux_seq, ref_bp2_w_aux_len, filter_with_pos, &ref2_aln, 0);

        StripedSmithWaterman::Alignment& alt_aln = alt1_aln.sw_score >= alt2_aln.sw_score ? alt1_aln : alt2_aln;
        StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
        bool alt_or_ref_read = false;
        bool add_alt_bp1_better_read = false;
        bool add_alt_bp2_better_read = false;
        bool add_ref_bp1_better_read = false;
        bool add_ref_bp2_better_read = false;
        if (alt_aln.sw_score > ref_aln.sw_score) {
            if (alt1_aln.sw_score >= alt2_aln.sw_score && alt1_covers_bp1) {
                add_alt_bp1_better_read = true;
                alt_or_ref_read = true;
            }
            if (alt1_aln.sw_score <= alt2_aln.sw_score && alt2_covers_bp2) {
                add_alt_bp2_better_read = true;
                alt_or_ref_read = true;
            }
        } else if (!only_allow_alt && alt_aln.sw_score < ref_aln.sw_score && !is_clipped(ref_aln, config.min_clip_len)) {
            if (ref1_aln.sw_score >= ref2_aln.sw_score && ref1_aln.ref_begin <= ref_bp1_w_aux_pos && ref1_aln.ref_end >= ref_bp1_w_aux_pos) {
                add_ref_bp1_better_read = true;
                alt_or_ref_read = true;
            }
            if (ref1_aln.sw_score <= ref2_aln.sw_score && ref2_aln.ref_begin <= ref_bp2_w_aux_pos && ref2_aln.ref_end >= ref_bp2_w_aux_pos) {
                add_ref_bp2_better_read = true;
                alt_or_ref_read = true;
            }
        } else if (!only_allow_alt && alt_aln.sw_score == ref_aln.sw_score) {
            ins->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) {
                ins->sample_info.alt_ref_equal_reads_highmq++;
            }
        }

        if (alt_or_ref_read && reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, ins)) {
            if (add_alt_bp1_better_read || add_ref_bp1_better_read) {
                ins->sample_info.assigned_to_other_sv_bp1_reads++;
            }
            if (add_alt_bp2_better_read || add_ref_bp2_better_read) {
                ins->sample_info.assigned_to_other_sv_bp2_reads++;
            }
            continue;
        }

        if (add_alt_bp1_better_read) {
            alt_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_bp1_better_scores.push_back(alt1_aln.sw_score);
            alt_bp1_better_read_positions.push_back(alt1_aln.ref_begin);
        }
        if (add_alt_bp2_better_read) {
            alt_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_bp2_better_scores.push_back(alt2_aln.sw_score);
            alt_bp2_better_read_positions.push_back(alt2_aln.ref_begin);
        }
        if (add_ref_bp1_better_read) {
            ref_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        }
        if (add_ref_bp2_better_read) {
            ref_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        }

        if (alt_bp1_better_reads.size() + alt_bp2_better_reads.size() + ref_bp1_better_reads.size() + ref_bp2_better_reads.size() + ins->sample_info.alt_ref_equal_reads > 4 * stats.get_max_depth(ins->chr)) {
            alt_bp1_better_reads.clear();
            alt_bp2_better_reads.clear();
            ref_bp1_better_reads.clear();
            ref_bp2_better_reads.clear();
            ins->sample_info.alt_ref_equal_reads = 0;
            ins->sample_info.alt_ref_equal_reads_highmq = 0;
            ins->sample_info.too_deep = true;
            break;
        }
    }

    if (evidence_logger) evidence_logger->log_reads_associations(ins->id, 1, alt_bp1_better_reads, alt_bp1_better_scores);
    if (evidence_logger) evidence_logger->log_reads_associations(ins->id, 2, alt_bp2_better_reads, alt_bp2_better_scores);

    std::string alt_bp1_consensus_seq, alt_bp2_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_bp1_avg_score, alt_bp2_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_bp1_stddev_score, alt_bp2_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bool> alt_bp1_is_exact_read, alt_bp2_is_exact_read, ref_bp1_is_exact_read, ref_bp2_is_exact_read;
    auto alt_bp1_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(alt_bp1_seq, alt_bp1_better_reads, std::vector<bool>(), alt_bp1_consensus_seq, alt_bp1_avg_score, alt_bp1_stddev_score, alt_bp1_is_exact_read);
    delete[] alt_bp1_seq;
    auto alt_bp2_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(alt_bp2_seq, alt_bp2_better_reads, std::vector<bool>(), alt_bp2_consensus_seq, alt_bp2_avg_score, alt_bp2_stddev_score, alt_bp2_is_exact_read);
    delete[] alt_bp2_seq;

    std::vector<int> alt_bp1_better_read_positions_consistent = get_consistent_reads_start_positions(alt_bp1_better_reads_consistent, alt_bp1_better_reads, alt_bp1_better_read_positions);
    ins->sample_info.alt1_occ_ratio = occ_ratio(alt_bp1_better_read_positions_consistent, alt1_ref_diff_reads_expected_positions.size());

    std::vector<int> alt_bp2_better_read_positions_consistent = get_consistent_reads_start_positions(alt_bp2_better_reads_consistent, alt_bp2_better_reads, alt_bp2_better_read_positions);
    ins->sample_info.alt2_occ_ratio = occ_ratio(alt_bp2_better_read_positions_consistent, alt2_ref_diff_reads_expected_positions.size());

    if (reassign_evidence) {  // increment ORC counters for other SVs that lost support from these reads
        std::vector<bam1_t*> alt_better_reads_consistent;
        std::unordered_map<std::string, bool> alt_better_read_is_exact;
        std::unordered_set<std::string> seen;
        for (int i = 0; i < alt_bp1_better_reads_consistent.size(); i++) {
            const auto& r = alt_bp1_better_reads_consistent[i];
            std::string read_name = read_name_with_suffix(r.get());
            alt_better_read_is_exact[read_name] |= alt_bp1_is_exact_read[i];
            if (!seen.count(read_name)) {
                seen.insert(read_name);
                alt_better_reads_consistent.push_back(r.get());
            }
        }
        for (int i = 0; i < alt_bp2_better_reads_consistent.size(); i++) {
            const auto& r = alt_bp2_better_reads_consistent[i];
            std::string read_name = read_name_with_suffix(r.get());
            alt_better_read_is_exact[read_name] |= alt_bp2_is_exact_read[i];
            if (!seen.count(read_name)) {
                seen.insert(read_name);
                alt_better_reads_consistent.push_back(r.get());
            }
        }
        for (bam1_t* r : alt_better_reads_consistent) {
            for (std::pair<std::string, int>& ov : evidence_map->get_non_chosen_svs_for_read(r)) {
                increase_orc(sv_map, ov.first, ov.second, r, get_mq(r) >= config.high_confidence_mapq, alt_better_read_is_exact[read_name_with_suffix(r)]);
            }
        }
    }

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    auto ref_bp1_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_reads, std::vector<bool>(), ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score, ref_bp1_is_exact_read);
    delete[] ref_bp1_seq;

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    auto ref_bp2_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_reads, std::vector<bool>(), ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score, ref_bp2_is_exact_read);
    delete[] ref_bp2_seq;

    if (alt_bp1_consensus_seq.length() >= 2*config.min_clip_len) {
        // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_bp1_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_bp1_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp1_consensus, candidate_reads_for_extension_itree, ins_start-stats.max_is, ins_start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_bp1_consensus, candidate_reads_for_extension_itree, ins_start, ins_start+stats.max_is, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        ins->sample_info.alt_lext_reads = alt_bp1_consensus->left_ext_reads;
        ins->sample_info.alt_rext_reads = alt_bp1_consensus->right_ext_reads;
        ins->sample_info.hq_alt_lext_reads = alt_bp1_consensus->hq_left_ext_reads;
        ins->sample_info.hq_alt_rext_reads = alt_bp1_consensus->hq_right_ext_reads;
        alt_bp1_consensus_seq = alt_bp1_consensus->sequence;
    }
    if (alt_bp2_consensus_seq.length() >= 2*config.min_clip_len) {
        std::shared_ptr<consensus_t> alt_bp2_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_bp2_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp2_consensus, candidate_reads_for_extension_itree, ins_end-stats.max_is, ins_end, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        extend_consensus_to_right(alt_bp2_consensus, candidate_reads_for_extension_itree, ins_end, ins_end+stats.max_is, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        ins->sample_info.alt_lext_reads += alt_bp2_consensus->left_ext_reads;
        ins->sample_info.alt_rext_reads += alt_bp2_consensus->right_ext_reads;
        ins->sample_info.hq_alt_lext_reads += alt_bp2_consensus->hq_left_ext_reads;
        ins->sample_info.hq_alt_rext_reads += alt_bp2_consensus->hq_right_ext_reads;
        alt_bp2_consensus_seq = alt_bp2_consensus->sequence;
    }

    ins->sample_info.ext_alt_consensus1_length = alt_bp1_consensus_seq.length();
    ins->sample_info.ext_alt_consensus2_length = alt_bp2_consensus_seq.length();

    ref1_aln.Clear();
    alt1_aln.Clear();
    if (alt_bp1_consensus_seq.length() >= 2*config.min_clip_len) {
        hts_pos_t ref_bp1_start = ins_start-alt_bp1_consensus_seq.length();
        if (ref_bp1_start < 0) ref_bp1_start = 0;
        hts_pos_t ref_bp1_end = ins_start+alt_bp1_consensus_seq.length();
        if (ref_bp1_end > contig_len) ref_bp1_end = contig_len;
        aligner.Align(alt_bp1_consensus_seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_end-ref_bp1_start, filter_with_pos_and_cigar, &ref1_aln, 0);

        char* lf_seq = generate_haplotype_left(contig_seq, ins_start-1, alt_bp1_consensus_seq.length(), ins->aux_indels, ins->aux_snps);
        int alt_lf_len = strlen(lf_seq);

        char* rf_seq = generate_haplotype_right(contig_seq, contig_len, ins_end, alt_bp1_consensus_seq.length(),  ins->aux_indels, ins->aux_snps);
        int alt_rf_len = strlen(rf_seq);

        char* alt_bp1_seq = new char[2*alt_bp1_consensus_seq.length()+1];
        strncpy(alt_bp1_seq, lf_seq, alt_lf_len);
        int ins_seq_portion_len = std::min(ins->ins_seq.length(), alt_bp1_consensus_seq.length());
        strncpy(alt_bp1_seq+alt_lf_len, ins->ins_seq.c_str(), ins_seq_portion_len);
        int extra_len = alt_bp1_consensus_seq.length() - ins->ins_seq.length();
        if (extra_len > 0) {
            if (extra_len > alt_rf_len) extra_len = alt_rf_len; // there may not be enough bp if near the end of contig
            strncpy(alt_bp1_seq+alt_lf_len+ins_seq_portion_len, rf_seq, extra_len);
        } else {
            extra_len = 0;
        }
        int alt_bp1_seq_len = alt_lf_len + ins_seq_portion_len + extra_len;
        alt_bp1_seq[alt_bp1_seq_len] = 0;

        aligner.Align(alt_bp1_consensus_seq.c_str(), alt_bp1_seq, alt_bp1_seq_len, filter_with_pos_and_cigar, &alt1_aln, 0);
        delete[] alt_bp1_seq;

        // length of the left and right flanking regions of the insertion covered by alt_bp1_consensus_seq
        int lf_aln_rlen = std::max(0, alt_lf_len - alt1_aln.ref_begin);
        int rf_aln_rlen = std::max(0, alt1_aln.ref_end+1 - alt_lf_len - ins_seq_portion_len);

        // length of the alt_consensus_seq covering left and right flanking regions of the deletion
        // note that this may be different from lf_aln_rlen and rf_aln_rlen, since the aln can include indels
        auto query_lh_aln_score = find_aln_prefix_score(alt1_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt1_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        ins->sample_info.alt_consensus1_split_size1 = query_lh_aln_score.second;
        ins->sample_info.alt_consensus1_split_size2 = query_rh_aln_score.second;
        ins->sample_info.alt_consensus1_split_score1 = query_lh_aln_score.first;
        ins->sample_info.alt_consensus1_split_score2 = query_rh_aln_score.first;

        ins->left_anchor_aln->start = ins_start - lf_aln_rlen;
        ins->left_anchor_aln->end = ins_start;
        ins->left_anchor_aln->seq_len = lf_aln_rlen;

        StripedSmithWaterman::Alignment aln;

        char* ref_seq = concat3(lf_seq, contig_seq+ins_start, rf_seq, alt_lf_len, ins_end-ins_start, alt_rf_len);

        aln.Clear();
        std::string lh_query = alt_bp1_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), ref_seq, alt_lf_len+ins_end-ins_start+alt_rf_len, filter_with_score_only, &aln, 0);
        ins->sample_info.alt_consensus1_split_score1_ind_aln = aln.sw_score;

        aln.Clear();
        std::string rh_query = alt_bp1_consensus_seq.substr(alt_bp1_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), ref_seq, alt_lf_len+ins_end-ins_start+alt_rf_len, filter_with_score_only, &aln, 0);
        ins->sample_info.alt_consensus1_split_score2_ind_aln = aln.sw_score;

        delete[] ref_seq;
        delete[] lf_seq;
        delete[] rf_seq;
    }

    ref2_aln.Clear();
    alt2_aln.Clear();
    if (alt_bp2_consensus_seq.length() >= 2*config.min_clip_len) {
        hts_pos_t ref_bp2_start = ins_end-alt_bp2_consensus_seq.length();
        if (ref_bp2_start < 0) ref_bp2_start = 0;
        hts_pos_t ref_bp2_end = ins_end+alt_bp2_consensus_seq.length();
        if (ref_bp2_end > contig_len) ref_bp2_end = contig_len; 
        aligner.Align(alt_bp2_consensus_seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_end-ref_bp2_start, filter_with_pos_and_cigar, &ref2_aln, 0);

        char* lf_seq = generate_haplotype_left(contig_seq, ins_start-1, alt_bp2_consensus_seq.length(), ins->aux_indels, ins->aux_snps);
        int alt_lf_len = strlen(lf_seq);

        char* rf_seq = generate_haplotype_right(contig_seq, contig_len, ins_end, alt_bp2_consensus_seq.length(),  ins->aux_indels, ins->aux_snps);
        int alt_rf_len = strlen(rf_seq);

        alt_bp2_seq = new char[2*alt_bp2_consensus_seq.length()+1];
        int extra_len = alt_bp2_consensus_seq.length() - ins->ins_seq.length();
        if (extra_len > 0) {
            if (extra_len > alt_lf_len) extra_len = alt_lf_len; // there may not be enough bp if near the start of contig
            strncpy(alt_bp2_seq, lf_seq+(alt_lf_len - extra_len), extra_len);
        } else {
            extra_len = 0;
        }
        int ins_seq_portion_len = std::min(ins->ins_seq.length(), alt_bp2_consensus_seq.length());
        strncpy(alt_bp2_seq+extra_len, ins->ins_seq.c_str()+(ins->ins_seq.length()-ins_seq_portion_len), ins_seq_portion_len);
        strncpy(alt_bp2_seq+extra_len+ins_seq_portion_len, rf_seq, alt_rf_len);
        int alt_bp2_seq_len = extra_len + ins_seq_portion_len + alt_rf_len;
        alt_bp2_seq[alt_bp2_seq_len] = 0;

        aligner.Align(alt_bp2_consensus_seq.c_str(), alt_bp2_seq, alt_bp2_seq_len, filter_with_pos_and_cigar, &alt2_aln, 0);
        delete[] alt_bp2_seq;

        // length of the left and right flanking regions of the insertion covered by alt_consensus_seq
        int lf_aln_rlen = std::max(0, extra_len - alt2_aln.ref_begin);
        int rf_aln_rlen = std::max(0, alt2_aln.ref_end+1 - (extra_len + ins_seq_portion_len));

        // length of the alt_consensus_seq covering left and right flanking regions of the insertion
        // note that this may be different from lf_aln_rlen and rf_aln_rlen, since the aln can include indels
        auto query_lh_aln_score = find_aln_prefix_score(alt2_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt2_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        ins->sample_info.alt_consensus2_split_size1 = query_lh_aln_score.second;
        ins->sample_info.alt_consensus2_split_size2 = query_rh_aln_score.second;
        ins->sample_info.alt_consensus2_split_score1 = query_lh_aln_score.first;
        ins->sample_info.alt_consensus2_split_score2 = query_rh_aln_score.first;

        ins->right_anchor_aln->start = ins_end;
        ins->right_anchor_aln->end = ins_end + rf_aln_rlen;
        ins->right_anchor_aln->seq_len = rf_aln_rlen;

        StripedSmithWaterman::Alignment aln;

        char* ref_seq = concat3(lf_seq, contig_seq+ins_start, rf_seq, alt_lf_len, ins_end-ins_start, alt_rf_len);

        aln.Clear();
        std::string lh_query = alt_bp2_consensus_seq.substr(0, query_lh_aln_score.second);
        aligner.Align(lh_query.c_str(), ref_seq, alt_lf_len+ins_end-ins_start+alt_rf_len, filter_with_score_only, &aln, 0);
        ins->sample_info.alt_consensus2_split_score1_ind_aln = aln.sw_score;

        aln.Clear();
        std::string rh_query = alt_bp2_consensus_seq.substr(alt_bp2_consensus_seq.length()-query_rh_aln_score.second);
        aligner.Align(rh_query.c_str(), ref_seq, alt_lf_len+ins_end-ins_start+alt_rf_len, filter_with_score_only, &aln, 0);
        ins->sample_info.alt_consensus2_split_score2_ind_aln = aln.sw_score;

        delete[] ref_seq;
        delete[] lf_seq;
        delete[] rf_seq;
    }

    ins->sample_info.ext_alt_consensus1_to_ref_score = ref1_aln.query_end - ref1_aln.query_begin - ref1_aln.mismatches;
    ins->sample_info.ext_alt_consensus2_to_ref_score = ref2_aln.query_end - ref2_aln.query_begin - ref2_aln.mismatches;
    ins->sample_info.ext_alt_consensus1_to_alt_score = alt1_aln.query_end - alt1_aln.query_begin - alt1_aln.mismatches;
    ins->sample_info.ext_alt_consensus2_to_alt_score = alt2_aln.query_end - alt2_aln.query_begin - alt2_aln.mismatches;

    set_bp_consensus_info(ins->sample_info.alt_bp1.reads_info, alt_bp1_better_reads.size(), alt_bp1_better_reads_consistent, alt_bp1_is_exact_read, alt_bp1_avg_score, alt_bp1_stddev_score);
    set_bp_consensus_info(ins->sample_info.alt_bp2.reads_info, alt_bp2_better_reads.size(), alt_bp2_better_reads_consistent, alt_bp2_is_exact_read, alt_bp2_avg_score, alt_bp2_stddev_score);
    set_bp_consensus_info(ins->sample_info.ref_bp1.reads_info, ref_bp1_better_reads.size(), ref_bp1_better_reads_consistent, ref_bp1_is_exact_read, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(ins->sample_info.ref_bp2.reads_info, ref_bp2_better_reads.size(), ref_bp2_better_reads_consistent, ref_bp2_is_exact_read, ref_bp2_avg_score, ref_bp2_stddev_score);

    delete[] ref_bp1_w_aux_seq;
    delete[] ref_bp2_w_aux_seq;

    free(regions[0]);
    free(regions[1]);

    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_inss(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<insertion_t*> inss,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t& stats, config_t& config, contig_map_t& contig_map,
    bam_pool_t* bam_pool, std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq_chr,
    std::vector<double>* global_crossing_isize_dist, evidence_logger_t* evidence_logger,
    bool reassign_evidence, evidence_map_t* evidence_map, std::unordered_map<std::string, std::shared_ptr<sv_t>>* sv_map) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (insertion_t* ins : inss) {
        target_ivals.push_back({ins->start-stats.max_is, ins->start+stats.max_is});
        target_ivals.push_back({ins->end-stats.max_is, ins->end+stats.max_is});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_pool->get_bam_reader(id), candidate_reads_for_extension);
             
    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);
    for (insertion_t* ins : inss) { 
        genotype_ins(ins, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner, evidence_logger, reassign_evidence, evidence_map, *sv_map);
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    depth_filter_ins(contig_name, inss, bam_file, config, stats);
    calculate_ptn_ratio(contig_name, inss, bam_file, config, stats, evidence_logger, false, evidence_map, *mateseqs_w_mapq_chr);
    std::vector<sv_t*> inss_sv(inss.begin(), inss.end());
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, inss_sv, bam_file, config, stats);

    release_mates(contig_id);
}

#endif // GENOTYPE_INSS_H
