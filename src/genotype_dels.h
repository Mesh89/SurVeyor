#ifndef GENOTYPE_DELS_H
#define GENOTYPE_DELS_H

#include "htslib/sam.h"
#include "types.h"
#include "sam_utils.h"
#include "utils.h"
#include "stat_tests.h"
#include "../libs/ssw_cpp.h"

#include "genotype.h"


void genotype_del(deletion_t* del, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                std::unordered_map<std::string, std::pair<std::string, int> >& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner,
                evidence_logger_t* evidence_logger, bool reassign_evidence, evidence_map_t* evidence_map) {
    int del_start = del->start, del_end = del->end;

    hts_pos_t extend = stats.read_len-1;

    // build alt allele
    /* POS in VCF is the base BEFORE the deletion - i.e., the first deleted base is POS+1.
     * Therefore, we want the ALT allele to *include* base POS
     * (note that POS is 1-based in the VCF file, but htslib kindly returns the 0-based coordinate here).
     * As for the END coordinate, my current understanding (which may change) is that it represents the last base deleted.
     * Therefore, the ALT allele should NOT include base END, i.e. it should start at END+1.
     * Here we shift both coordinates by 1, to make them the base immediately AFTER the breakpoints, which is a bit more intuitive for me. */
    del_start++; del_end++;

    // all ranges will be start-inclusive and end-exclusive, i.e. [a,b)
    char* lh_seq = generate_haplotype_left(contig_seq, del_start-1, extend, del->aux_indels, del->aux_snps);
    char* rh_seq = generate_haplotype_right(contig_seq, contig_len, del_end, extend, del->aux_indels, del->aux_snps);
    hts_pos_t alt_start = std::max(hts_pos_t(0), del_start-extend);
    hts_pos_t alt_end = std::min(del_end+extend, contig_len);
    hts_pos_t alt_lh_len = strlen(lh_seq), alt_rh_len = strlen(rh_seq);
    hts_pos_t alt_len = alt_lh_len + del->ins_seq.length() + alt_rh_len;
    char* alt_seq = new char[alt_len + 1];
    strncpy(alt_seq, lh_seq, alt_lh_len);
    strncpy(alt_seq+alt_lh_len, del->ins_seq.c_str(), del->ins_seq.length());
    strncpy(alt_seq+alt_lh_len+del->ins_seq.length(), rh_seq, alt_rh_len);
    alt_seq[alt_len] = 0;
    to_uppercase(alt_seq);

    // extract ref alleles - will be useful for consensus generation
    hts_pos_t ref_bp1_start = alt_start, ref_bp1_end = std::min(del_start+extend, contig_len);
    hts_pos_t ref_bp1_pos = del_start - ref_bp1_start;
    hts_pos_t ref_bp1_len = ref_bp1_end - ref_bp1_start;
    char* ref_bp1_seq = new char[ref_bp1_len + 1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;

    hts_pos_t del_len = del_end - del_start;
    hts_pos_t ref_bp1_w_aux_del_len = std::min(extend, del_len);
    hts_pos_t ref_bp1_w_aux_rh_len = std::min(alt_rh_len, extend - ref_bp1_w_aux_del_len);
    hts_pos_t ref_bp1_w_aux_len = alt_lh_len + ref_bp1_w_aux_del_len + ref_bp1_w_aux_rh_len;
    char* ref_bp1_w_aux_seq = new char[ref_bp1_w_aux_len + 1];
    strncpy(ref_bp1_w_aux_seq, lh_seq, alt_lh_len);
    strncpy(ref_bp1_w_aux_seq + alt_lh_len, contig_seq + del_start, ref_bp1_w_aux_del_len);
    strncpy(ref_bp1_w_aux_seq + alt_lh_len + ref_bp1_w_aux_del_len, rh_seq, ref_bp1_w_aux_rh_len);
    ref_bp1_w_aux_seq[ref_bp1_w_aux_len] = 0;

    hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), del_end-extend), ref_bp2_end = alt_end;
    hts_pos_t ref_bp2_pos = del_end - ref_bp2_start;
    hts_pos_t ref_bp2_len = ref_bp2_end - ref_bp2_start;
    char* ref_bp2_seq = new char[ref_bp2_len + 1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;

    hts_pos_t ref_bp2_w_aux_del_len = std::min(extend, del_len);
    hts_pos_t ref_bp2_w_aux_lh_len = std::min(alt_lh_len, extend - ref_bp2_w_aux_del_len);
    hts_pos_t ref_bp2_w_aux_len = ref_bp2_w_aux_lh_len + ref_bp2_w_aux_del_len + alt_rh_len;
    char* ref_bp2_w_aux_seq = new char[ref_bp2_w_aux_len + 1];
    strncpy(ref_bp2_w_aux_seq, lh_seq + alt_lh_len - ref_bp2_w_aux_lh_len, ref_bp2_w_aux_lh_len);
    strncpy(ref_bp2_w_aux_seq + ref_bp2_w_aux_lh_len, contig_seq + del_end - ref_bp2_w_aux_del_len, ref_bp2_w_aux_del_len);
    strncpy(ref_bp2_w_aux_seq + ref_bp2_w_aux_lh_len + ref_bp2_w_aux_del_len, rh_seq, alt_rh_len);
    ref_bp2_w_aux_seq[ref_bp2_w_aux_len] = 0;

    std::vector<char*> ref_seqs = {ref_bp1_w_aux_seq, ref_bp2_w_aux_seq};
    std::vector<hts_pos_t> ref_lens = {ref_bp1_w_aux_len, ref_bp2_w_aux_len};
    std::vector<hts_pos_t> alt_ref_diff_reads_expected_positions = get_diff_reads_expected_positions(ref_seqs, ref_lens, alt_seq, alt_len, stats.read_len);
    del->expected_alt1_reads_frac = (double) alt_ref_diff_reads_expected_positions.size() / std::max(hts_pos_t(1), alt_len - stats.read_len + 1);

    std::vector<std::shared_ptr<bam1_t>> alt_better_reads, ref_bp1_better_seqs, ref_bp2_better_seqs;
    std::vector<int> alt_better_read_positions;
    std::vector<int> alt_better_read_scores;
    
    bam1_t* read = bam_init1();

    std::stringstream l_region, r_region;
    l_region << del->chr << ":" << alt_start << "-" << ref_bp1_end;
    r_region << del->chr << ":" << ref_bp2_start << "-" << alt_end;

    char* regions[2];
    regions[0] = strdup(l_region.str().c_str());
    regions[1] = strdup(r_region.str().c_str());
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
    
    StripedSmithWaterman::Filter filter_with_pos(true, false, 0, 32767);
    StripedSmithWaterman::Filter filter_with_pos_and_cigar(true, true, 0, 32767);
    StripedSmithWaterman::Filter filter_with_score_only(false, false, 0, 32767);

    StripedSmithWaterman::Alignment alt_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < del_start || del_end < get_unclipped_start(read)) continue;
        if (del_start < get_unclipped_start(read) && get_unclipped_end(read) < del_end) continue;

        std::string seq;
        if (!is_samechr(read) || is_samestr(read)) continue;
        if (!bam_is_mrev(read)) {
            if (read->core.mpos < del_start-stats.max_is) continue; // positive mate and after deletion, potentially discordant...
            // however, there is an exception: both reads in a pair can be left-clipped on the end-side of the deletion (we allow 5bp tolerance)
            if (read->core.mpos > del_start && (abs(read->core.pos-read->core.mpos) > 5 || !is_left_clipped(read, config.min_clip_len))) continue;
            seq = get_sequence(read, true);
            rc(seq);
        } else {
            hts_pos_t mate_endpos = get_mate_endpos(read);
            if (mate_endpos > del_end+stats.max_is) continue; // negative mate and before deletion, potentially discordant...
            // however, there is an exception: both reads in a pair can be right-clipped on the start-side of the deletion (we allow 5bp tolerance)
            if (mate_endpos < del_end && (abs(mate_endpos-bam_endpos(read)) > 5 || !is_right_clipped(read, config.min_clip_len))) continue;
            seq = get_sequence(read, true);
        }

        // align to REF (two breakpoints)
        uint16_t ref_aln_score = 0;
        bool increase_ref_bp1_better = false, increase_ref_bp2_better = false;
        bool ref_is_exact_match = is_perfectly_aligned(read);
        if (ref_is_exact_match) {
            ref_aln_score = read->core.l_qseq;
            if (read->core.pos < del_start && bam_endpos(read) > del_start) {
                increase_ref_bp1_better = true;
            }
            if (read->core.pos < del_end && bam_endpos(read) > del_end) {
                increase_ref_bp2_better = true;
            }
        } else {
            aligner.Align(seq.c_str(), ref_bp1_w_aux_seq, ref_bp1_w_aux_len, filter_with_pos, &ref1_aln, 0);
            aligner.Align(seq.c_str(), ref_bp2_w_aux_seq, ref_bp2_w_aux_len, filter_with_pos, &ref2_aln, 0);
            ref_aln_score = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln.sw_score : ref2_aln.sw_score;
            if (ref1_aln.sw_score >= ref2_aln.sw_score && ref1_aln.ref_begin < ref_bp1_pos && ref1_aln.ref_end >= ref_bp1_pos) {
                increase_ref_bp1_better = true;
            }
            if (ref2_aln.sw_score >= ref1_aln.sw_score && ref2_aln.ref_begin < ref_bp2_pos && ref2_aln.ref_end >= ref_bp2_pos) {
                increase_ref_bp2_better = true;
            }
        }

        // align to ALT
        alt_aln = align_fast(aligner, seq.c_str(), alt_seq, alt_len, filter_with_pos_and_cigar, ref_is_exact_match);
        hts_pos_t alt_right_flank_pos = alt_lh_len + del->ins_seq.length();
        bool alt_spans_left_bp = alt_aln.ref_begin < alt_lh_len && alt_aln.ref_end >= alt_lh_len;
        bool alt_spans_right_bp = alt_aln.ref_begin < alt_right_flank_pos && alt_aln.ref_end >= alt_right_flank_pos;
        bool alt_spans_sv = alt_spans_left_bp || alt_spans_right_bp;
        bool alt_or_ref_read = false;
        bool add_alt_better_read = false;
        bool add_ref_bp1_better_seq = false;
        bool add_ref_bp2_better_seq = false;
        if (alt_aln.sw_score > ref_aln_score && alt_spans_sv) {
            if (is_clipped(alt_aln, config.min_clip_len)) continue;
            add_alt_better_read = true;
            alt_or_ref_read = true;
        } else if (ref_aln_score > alt_aln.sw_score) {
            if (increase_ref_bp1_better) {
                add_ref_bp1_better_seq = true;
                alt_or_ref_read = true;
            }
            if (increase_ref_bp2_better) {
                add_ref_bp2_better_seq = true;
                alt_or_ref_read = true;
            }
        } else {
            del->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) {
                del->sample_info.alt_ref_equal_reads_highmq++;
            }
        }

        // OAR reads remain excluded from AR; ORR reads continue into the regular RR statistics below.
        if (alt_or_ref_read && reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, del)) {
            if (add_alt_better_read && alt_spans_left_bp) {
                del->sample_info.oar_bp1_reads++;
                evidence_map->register_oar_support(del->sample_info, 1, read);
            }
            if (add_alt_better_read && alt_spans_right_bp) {
                del->sample_info.oar_bp2_reads++;
                evidence_map->register_oar_support(del->sample_info, 2, read);
            }
            if (add_ref_bp1_better_seq) {
                del->sample_info.orr_bp1_reads++;
                evidence_map->register_orr_support(del->sample_info, 1, read);
            }
            if (add_ref_bp2_better_seq) {
                del->sample_info.orr_bp2_reads++;
                evidence_map->register_orr_support(del->sample_info, 2, read);
            }
            if (add_alt_better_read) continue;
        }

        if (add_alt_better_read) {
            alt_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_better_read_scores.push_back(alt_aln.sw_score);
            alt_better_read_positions.push_back(alt_aln.ref_begin);
        }
        if (add_ref_bp1_better_seq) {
            ref_bp1_better_seqs.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        }
        if (add_ref_bp2_better_seq) {
            ref_bp2_better_seqs.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
        }

        if (alt_better_reads.size() + ref_bp1_better_seqs.size() + ref_bp2_better_seqs.size() + del->sample_info.alt_ref_equal_reads > 4 * stats.get_max_depth(del->chr)) {
            alt_better_reads.clear();
            ref_bp1_better_seqs.clear();
            ref_bp2_better_seqs.clear();
            del->sample_info.alt_ref_equal_reads = 0;
            del->sample_info.alt_ref_equal_reads_highmq = 0;
            evidence_map->clear_other_read_support_for_too_deep(del->sample_info);
            break;
        }
    }

    if (evidence_logger) evidence_logger->log_reads_associations(del->id, 1, alt_better_reads, alt_better_read_scores);

    std::string alt_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bool> alt_is_exact_read, ref_bp1_is_exact_read, ref_bp2_is_exact_read;
    auto alt_better_reads_consistent = gen_consensus_and_find_consistent_seqs_subset(alt_seq, alt_better_reads, std::vector<bool>(), alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read);
    auto ref_bp1_better_seqs_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_bp1_seq, ref_bp1_better_seqs, std::vector<bool>(), ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score, ref_bp1_is_exact_read);
    auto ref_bp2_better_seqs_consistent = gen_consensus_and_find_consistent_seqs_subset(ref_bp2_seq, ref_bp2_better_seqs, std::vector<bool>(), ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score, ref_bp2_is_exact_read);

    std::vector<int> alt_better_read_positions_consistent = get_consistent_reads_start_positions(alt_better_reads_consistent, alt_better_reads, alt_better_read_positions);
    del->sample_info.alt1_occ_ratio = occ_ratio(alt_better_read_positions_consistent, alt_ref_diff_reads_expected_positions.size());

    if (reassign_evidence) {
        for (int i = 0; i < alt_better_reads_consistent.size(); i++) {
            std::shared_ptr<bam1_t>& r = alt_better_reads_consistent[i];
            evidence_map->record_assigned_read_consistency(r.get(), get_mq(r.get()) >= config.high_confidence_mapq, alt_is_exact_read[i]);
        }
    }

    if (alt_consensus_seq.length() >= 2*config.min_clip_len) {
        // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_consensus_seq, 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, del_start-stats.max_is, del_start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr); 
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, del_end, del_end+stats.max_is, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr);
        
        del->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        del->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        del->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        del->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;
        
        char* lh_seq = generate_haplotype_left(contig_seq, del_start-1, alt_consensus_seq.length(), del->aux_indels, del->aux_snps);
        hts_pos_t lh_len = strlen(lh_seq);
        char* rh_seq = generate_haplotype_right(contig_seq, contig_len, del_end, alt_consensus_seq.length(), del->aux_indels, del->aux_snps);
        hts_pos_t rh_len = strlen(rh_seq);
        
        hts_pos_t lh_start = del_start - lh_len;
        if (lh_start < 0) lh_start = 0;
        hts_pos_t rh_end = del_end + rh_len;
        if (rh_end > contig_len) rh_end = contig_len;

        delete[] alt_seq;
        alt_seq = new char[lh_len + rh_len + del->ins_seq.length() + 1];
        strncpy(alt_seq, lh_seq, lh_len);
        strncpy(alt_seq+lh_len, del->ins_seq.c_str(), del->ins_seq.length());
        strncpy(alt_seq+lh_len+del->ins_seq.length(), rh_seq, rh_len);
        alt_seq[lh_len+rh_len+del->ins_seq.length()] = 0;

        // align to ref+SV
        aligner.Align(alt_consensus_seq.c_str(), alt_seq, lh_len+del->ins_seq.length()+rh_len, filter_with_pos_and_cigar, &alt_aln, 0);

        // length of the left and right flanking regions of the deletion covered by alt_consensus_seq
        int lf_aln_rlen = std::max(hts_pos_t(0), lh_len - alt_aln.ref_begin);
        int rf_aln_rlen = std::max(hts_pos_t(0), alt_aln.ref_end + 1 - lh_len - hts_pos_t(del->ins_seq.length()));

        // length of the alt_consensus_seq covering left and right flanking regions of the deletion
        // note that this may be different from lf_aln_rlen and rf_aln_rlen, since the aln can include indels
        auto query_lh_aln_score = find_aln_prefix_score(alt_aln.cigar, lf_aln_rlen, 1, -4, -6, -1);
        auto query_rh_aln_score = find_aln_suffix_score(alt_aln.cigar, rf_aln_rlen, 1, -4, -6, -1);
        int split_query_covered = query_lh_aln_score.second + query_rh_aln_score.second;

        del->sample_info.alt_consensus1_split_size1 = query_lh_aln_score.second;
        del->sample_info.alt_consensus1_split_size2 = query_rh_aln_score.second;
        del->sample_info.alt_consensus1_split_score1 = query_lh_aln_score.first;
        del->sample_info.alt_consensus1_split_score2 = query_rh_aln_score.first;

        del->left_anchor_aln->start = del_start - lf_aln_rlen;
        del->left_anchor_aln->end = del_start;
        del->left_anchor_aln->seq_len = lf_aln_rlen;
        del->right_anchor_aln->start = del_end;
        del->right_anchor_aln->end = del_end + rf_aln_rlen;
        del->right_anchor_aln->seq_len = rf_aln_rlen;

        // align to ref
        hts_pos_t lbp_start = lh_start, lbp_end = del_start + alt_consensus_seq.length();
        hts_pos_t rbp_start = del_end - alt_consensus_seq.length(), rbp_end = rh_end;
        if (lbp_end > contig_len) lbp_end = contig_len;
        if (rbp_start < 0) rbp_start = 0;
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+lbp_start, lbp_end-lbp_start, filter_with_pos_and_cigar, &ref1_aln, 0);
        aligner.Align(alt_consensus_seq.c_str(), contig_seq+rbp_start, rbp_end-rbp_start, filter_with_pos_and_cigar, &ref2_aln, 0);

        del->sample_info.ext_alt_consensus1_length = alt_consensus_seq.length();
        del->sample_info.ext_alt_consensus1_to_alt_score = alt_aln.query_end - alt_aln.query_begin - alt_aln.mismatches;
        del->sample_info.ext_alt_consensus1_to_ref_score = std::max(
            ref1_aln.query_end - ref1_aln.query_begin - ref1_aln.mismatches,
            ref2_aln.query_end - ref2_aln.query_begin - ref2_aln.mismatches
        );

        ref1_aln.Clear();
        std::string lh_query = alt_consensus_seq.substr(0, query_lh_aln_score.second);
        char* ref_seq_for_lh_aln = concat2(lh_seq, contig_seq+del_start, lh_len, lbp_end-del_start);
        aligner.Align(lh_query.c_str(), ref_seq_for_lh_aln, lh_len+lbp_end-del_start, filter_with_score_only, &ref1_aln, 0);
        del->sample_info.alt_consensus1_split_score1_ind_aln = ref1_aln.sw_score;
        delete[] ref_seq_for_lh_aln;
        delete[] lh_seq;

        ref2_aln.Clear();
        std::string rh_query = alt_consensus_seq.substr(alt_consensus_seq.length()-query_rh_aln_score.second);
        char* ref_seq_for_rh_aln = concat2(contig_seq+rbp_start, rh_seq, del_end-rbp_start, rh_len);
        aligner.Align(rh_query.c_str(), ref_seq_for_rh_aln, del_end-rbp_start+rh_len, filter_with_score_only, &ref2_aln, 0);
        del->sample_info.alt_consensus1_split_score2_ind_aln = ref2_aln.sw_score;
        delete[] ref_seq_for_rh_aln;
        delete[] rh_seq;
    }

    set_bp_consensus_info(del->sample_info.alt_bp1.reads_info, alt_better_reads.size(), alt_better_reads_consistent, alt_is_exact_read, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(del->sample_info.ref_bp1.reads_info, ref_bp1_better_seqs.size(), ref_bp1_better_seqs_consistent, ref_bp1_is_exact_read, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(del->sample_info.ref_bp2.reads_info, ref_bp2_better_seqs.size(), ref_bp2_better_seqs_consistent, ref_bp2_is_exact_read, ref_bp2_avg_score, ref_bp2_stddev_score);

    delete[] alt_seq;
    delete[] ref_bp1_seq;
    delete[] ref_bp1_w_aux_seq;
    delete[] ref_bp2_seq;
    delete[] ref_bp2_w_aux_seq;
    delete[] lh_seq;
    delete[] rh_seq;

    free(regions[0]);
    free(regions[1]);
    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_dels(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<deletion_t*> dels,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t& stats, config_t& config, contig_map_t& contig_map,
    bam_pool_t* bam_pool, std::unordered_map<std::string, std::pair<std::string, int> >* mateseqs_w_mapq_chr,
    std::string workdir, std::vector<double>* global_crossing_isize_dist, evidence_logger_t* evidence_logger,
    bool reassign_evidence, evidence_map_t* evidence_map) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

    int contig_id = contig_map.get_id(contig_name);
    read_mates(contig_id);

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);

    std::vector<hts_pair_pos_t> target_ivals;
    for (deletion_t* del : dels) {
        target_ivals.push_back({del->start-stats.max_is, del->start+stats.max_is});
        target_ivals.push_back({del->end-stats.max_is, del->end+stats.max_is});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_file, candidate_reads_for_extension);

    std::vector<deletion_t*> small_deletions, large_deletions;       
    std::vector<sv_t*> small_svs;  
    for (deletion_t* del : dels) {
        genotype_del(del, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, 
            stats, config, aligner, evidence_logger, reassign_evidence, evidence_map);
        if (-del->svlen() >= stats.max_is) {
            large_deletions.push_back(del);
        } else {
            small_deletions.push_back(del);
            small_svs.push_back(del);
        }
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    release_mates(contig_id);

    depth_filter_del(contig_name, dels, bam_file, config, stats);
    calculate_confidence_interval_size(contig_name, *global_crossing_isize_dist, small_svs, bam_file, config, stats);
    std::string mates_nms_file = workdir + "/workspace/long-pairs/" + std::to_string(contig_id) + ".txt";
    calculate_ptn_ratio(contig_name, dels, bam_file, config, stats, evidence_logger, false, evidence_map, mates_nms_file);
    count_stray_pairs(contig_name, dels, bam_file, config, stats);
}

#endif // GENOTYPE_DELS_H
