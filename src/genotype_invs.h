#ifndef GENOTYPE_INVS_H
#define GENOTYPE_INVS_H

#include "utils.h"
#include "stat_tests.h"

#include "genotype.h"

void genotype_small_inv(inversion_t* inv, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                ext_mate_map_t& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner) {
    hts_pos_t sv_start = inv->start, sv_end = inv->end;

    hts_pos_t extend = stats.read_len + 20;

    // build alt allele
    /*
     * POS in VCF is the base BEFORE the inversion
     * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" inversion (no deletion),POS == END, otherwise the last base deleted
     * As usual, in order to make intervals [ ), we increase the coordinates by 1
     */
    sv_start++; sv_end++;

    hts_pos_t alt_start = std::max(hts_pos_t(0), sv_start-extend);
    hts_pos_t alt_end = std::min(sv_end+extend, contig_len);
    int alt_lf_len = sv_start-alt_start, alt_rf_len = alt_end-sv_end;
    char* inv_seq;
    if (inv->ins_seq.empty()) {
        inv_seq = new char[inv->inv_end-inv->inv_start+1];
        strncpy(inv_seq, contig_seq+inv->inv_start, inv->inv_end-inv->inv_start);
        inv_seq[inv->inv_end-inv->inv_start] = 0;
        rc(inv_seq);
    } else {
        inv_seq = new char[inv->ins_seq.length()+1];
        strncpy(inv_seq, inv->ins_seq.c_str(), inv->ins_seq.length());
        inv_seq[inv->ins_seq.length()] = '\0';
    }
    int inv_len = strlen(inv_seq);
    int alt_len = alt_lf_len + inv_len + alt_rf_len;
    char* alt_seq = new char[alt_len+1];
    strncpy(alt_seq, contig_seq+alt_start, alt_lf_len);
    strncpy(alt_seq+alt_lf_len, inv_seq, inv_len);
    strncpy(alt_seq+alt_lf_len+inv_len, contig_seq+sv_end, alt_rf_len);
    alt_seq[alt_len] = 0;

    hts_pos_t ref_start = std::max(hts_pos_t(0), sv_start-extend);
    hts_pos_t ref_end = std::min(sv_end+extend, contig_len);
    int ref_len = ref_end - ref_start;
    char* ref_seq = new char[ref_len+1];
    strncpy(ref_seq, contig_seq+ref_start, ref_len);
    ref_seq[ref_len] = 0;
    inv->ref1_hp_len = longest_homopolymer_len(ref_seq, ref_len);
    inv->alt1_hp_len = longest_homopolymer_len(alt_seq, alt_len);

    std::stringstream region;
    region << inv->chr << ":" << ref_start << "-" << ref_end;
    char* regions[1];
    regions[0] = strdup(region.str().c_str());

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 1);

    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> alt_better_seqs, ref_better_seqs;
    std::vector<bool> alt_better_seqs_isrc, ref_better_seqs_isrc;

    StripedSmithWaterman::Filter filter_with_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Alignment alt_aln, ref_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;
        if (get_unclipped_end(read) < sv_start || sv_end < get_unclipped_start(read)) continue;
        if (sv_start < get_unclipped_start(read) && get_unclipped_end(read) < sv_end) continue;

        std::string seq = get_sequence(read);
        bool is_rc = false;
        if (bam_is_rev(read) && bam_is_mrev(read) && sv_end+stats.read_len/2 <= get_mate_endpos(read) || 
            !bam_is_rev(read) && !bam_is_mrev(read) && read->core.mpos <= sv_start-stats.read_len/2) {
            rc(seq);
            is_rc = true;
        }

        // align to ALT
        aligner.Align(seq.c_str(), alt_seq, alt_len, filter_with_score_only, &alt_aln, 0);

        // // align to REF
        aligner.Align(seq.c_str(), ref_seq, ref_len, filter_with_score_only, &ref_aln, 0);

        if (alt_aln.sw_score > ref_aln.sw_score) {
            alt_better_seqs.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            alt_better_seqs_isrc.push_back(is_rc);
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            ref_better_seqs.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
            ref_better_seqs_isrc.push_back(is_rc);
        } else {
            inv->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) {
                inv->sample_info.alt_ref_equal_reads_highmq++;
            }
        }

        if (alt_better_seqs.size() + ref_better_seqs.size() + inv->sample_info.alt_ref_equal_reads > 4 * stats.get_max_depth(inv->chr)) {
            alt_better_seqs.clear();
            ref_better_seqs.clear();
            inv->sample_info.alt_ref_equal_reads = 0;
            inv->sample_info.alt_ref_equal_reads_highmq = 0;
            inv->sample_info.too_deep = true;
            break;
        }
    }

    std::string alt_consensus_seq, ref_consensus_seq;
    double alt_avg_score, ref_avg_score;
    double alt_stddev_score, ref_stddev_score;
    std::vector<bool> alt_is_exact_read, ref_is_exact_read;
    auto alt_is_consistent_read = gen_consensus_and_classify_seqs(alt_seq, alt_better_seqs, alt_better_seqs_isrc, alt_consensus_seq, alt_avg_score, alt_stddev_score, alt_is_exact_read, nullptr);
    auto ref_is_consistent_read = gen_consensus_and_classify_seqs(ref_seq, ref_better_seqs, ref_better_seqs_isrc, ref_consensus_seq, ref_avg_score, ref_stddev_score, ref_is_exact_read, nullptr);
    auto score_inv_consensus = [&](const std::string& consensus_seq) {
        alignment_targets_t targets;
        hts_pos_t target_start = std::max<hts_pos_t>(0, sv_start-GENOTYPE_CONSENSUS_EXTENSION);
        int target_lf_len = sv_start-target_start;
        hts_pos_t target_end = std::min<hts_pos_t>(sv_end+GENOTYPE_CONSENSUS_EXTENSION, contig_len);
        int target_rf_len = target_end-sv_end;
        targets.alt_len = target_lf_len+inv_len+target_rf_len;
        targets.alt_seq = new char[targets.alt_len+1];
        strncpy(targets.alt_seq, contig_seq+target_start, target_lf_len);
        strncpy(targets.alt_seq+target_lf_len, inv_seq, inv_len);
        strncpy(targets.alt_seq+target_lf_len+inv_len, contig_seq+sv_end, target_rf_len);
        targets.alt_seq[targets.alt_len] = 0;
        targets.ref_seqs.push_back(contig_seq+target_start);
        targets.ref_lens.push_back(target_end-target_start);
        targets.ref_starts.push_back(target_start);
        targets.edits.push_back({allele_edit_kind_t::INDEL, sv_start, sv_end, target_lf_len, target_lf_len+inv_len, int(sv_end-sv_start+inv_len), true});
        targets.left_flank_end = target_lf_len;
        targets.right_flank_start = target_lf_len+inv_len;
        targets.left_independent_ref_seq = targets.right_independent_ref_seq = contig_seq+target_start;
        targets.left_independent_ref_len = targets.right_independent_ref_len = target_end-target_start;
        consensus_alignment_metrics_t metrics = score_consensus_alignment(consensus_seq, targets, aligner);
        delete[] targets.alt_seq;
        return metrics;
    };

    if (alt_consensus_seq.length() >= 2*config.min_clip_len) {
        consensus_alignment_metrics_t unextended_metrics = score_inv_consensus(alt_consensus_seq);
        inv->sample_info.alt_consensus1_metrics = unextended_metrics;

        // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_consensus_seq, std::string(alt_consensus_seq.length(), '!'), 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_consensus, candidate_reads_for_extension_itree, std::max<hts_pos_t>(0, sv_start-GENOTYPE_CONSENSUS_EXTENSION), sv_start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        extend_consensus_to_right(alt_consensus, candidate_reads_for_extension_itree, sv_end, std::min<hts_pos_t>(contig_len, sv_end+GENOTYPE_CONSENSUS_EXTENSION), contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        inv->sample_info.alt_lext_reads = alt_consensus->left_ext_reads;
        inv->sample_info.alt_rext_reads = alt_consensus->right_ext_reads;
        inv->sample_info.hq_alt_lext_reads = alt_consensus->hq_left_ext_reads;
        inv->sample_info.hq_alt_rext_reads = alt_consensus->hq_right_ext_reads;
        alt_consensus_seq = alt_consensus->sequence;

        consensus_alignment_metrics_t extended_metrics = score_inv_consensus(alt_consensus_seq);
        inv->sample_info.ext_alt_consensus1_metrics.length = extended_metrics.length;
        inv->sample_info.ext_alt_consensus1_metrics.alt_score = extended_metrics.alt_score;
        inv->sample_info.ext_alt_consensus1_metrics.ref_score = extended_metrics.ref_score;
        inv->sample_info.ext_alt_consensus1_metrics.covered_edit_distance = extended_metrics.covered_edit_distance;
        inv->sample_info.ext_alt_consensus1_metrics.main_edit_covered = extended_metrics.main_edit_covered;

        hts_pos_t target_start = std::max<hts_pos_t>(0, sv_start-GENOTYPE_CONSENSUS_EXTENSION);
        int target_lf_len = sv_start-target_start;
        int lf_aln_rlen = std::max(0, target_lf_len-extended_metrics.alt_ref_begin);
        int rf_aln_rlen = std::max(0, extended_metrics.alt_ref_end-target_lf_len);
        inv->left_anchor_aln->start = sv_start-lf_aln_rlen;
        inv->left_anchor_aln->end = sv_start;
        inv->left_anchor_aln->seq_len = lf_aln_rlen;
        inv->right_anchor_aln->start = sv_end;
        inv->right_anchor_aln->end = sv_end+rf_aln_rlen;
        inv->right_anchor_aln->seq_len = rf_aln_rlen;
    }

    set_bp_consensus_info(inv->sample_info.alt_bp1.reads_info, alt_better_seqs, alt_is_consistent_read, alt_is_exact_read, alt_avg_score, alt_stddev_score);
    set_bp_consensus_info(inv->sample_info.ref_bp1.reads_info, ref_better_seqs, ref_is_consistent_read, ref_is_exact_read, ref_avg_score, ref_stddev_score);

    delete[] alt_seq;
    delete[] ref_seq;
    delete[] inv_seq;
    
    free(regions[0]);
    bam_destroy1(read);
    hts_itr_destroy(iter);
}

void genotype_large_inv(inversion_t* inv, open_samFile_t* bam_file, IntervalTree<ext_read_t*>& candidate_reads_for_extension_itree, 
                ext_mate_map_t& mateseqs_w_mapq_chr, char* contig_seq, hts_pos_t contig_len,
                stats_t& stats, config_t& config, StripedSmithWaterman::Aligner& aligner) {

    hts_pos_t sv_start = inv->start, sv_end = inv->end;

    hts_pos_t extend = stats.read_len + 20;

    // build alt allele
    /*
     * POS in VCF is the base BEFORE the inversion
     * END seems to be the base BEFORE the reference resumes - i.e., for a "clean" inversion (no deletion),POS == END, otherwise the last base deleted
     * As usual, in order to make intervals [ ), we increase the coordinates by 1
     */
    sv_start++; sv_end++;

    hts_pos_t alt_start = std::max(hts_pos_t(0), sv_start-extend);
    hts_pos_t alt_end = std::min(sv_end+extend, contig_len);
    int alt_lf_len = sv_start-alt_start, alt_rf_len = alt_end-sv_end;
    hts_pos_t inv_border_len = std::min(extend, inv->inv_end-inv->inv_start);

    char* inv_prefix = new char[inv_border_len+1];
    strncpy(inv_prefix, contig_seq+inv->inv_start, inv_border_len);
    inv_prefix[inv_border_len] = 0;

    char* inv_prefix_rc = strdup(inv_prefix);
    rc(inv_prefix_rc);
    delete[] inv_prefix;

    char* inv_suffix = new char[inv_border_len+1];
    strncpy(inv_suffix, contig_seq+inv->inv_end-inv_border_len, inv_border_len);
    inv_suffix[inv_border_len] = 0;

    char* inv_suffix_rc = strdup(inv_suffix);
    rc(inv_suffix_rc);
    delete[] inv_suffix;

    int alt_bp1_len = alt_lf_len + inv_border_len;
    char* alt_bp1_seq = new char[alt_bp1_len+1];
    strncpy(alt_bp1_seq, contig_seq+alt_start, alt_lf_len);
    strncpy(alt_bp1_seq+alt_lf_len, inv_suffix_rc, inv_border_len);
    alt_bp1_seq[alt_bp1_len] = 0;

    int alt_bp2_len = inv_border_len + alt_rf_len;
    char* alt_bp2_seq = new char[alt_bp2_len+1];
    strncpy(alt_bp2_seq, inv_prefix_rc, inv_border_len);
    strncpy(alt_bp2_seq+inv_border_len, contig_seq+sv_end, alt_rf_len);
    alt_bp2_seq[alt_bp2_len] = 0;
    free(inv_prefix_rc);
    free(inv_suffix_rc);

    hts_pos_t ref_bp1_start = std::max(hts_pos_t(0), sv_start-extend);
    hts_pos_t ref_bp1_end = std::min(sv_start+extend, contig_len);
    hts_pos_t ref_bp1_len = ref_bp1_end-ref_bp1_start;
    hts_pos_t ref_bp2_start = std::max(hts_pos_t(0), sv_end-extend);
    hts_pos_t ref_bp2_end = std::min(sv_end+extend, contig_len);
    hts_pos_t ref_bp2_len = ref_bp2_end-ref_bp2_start;

    char* regions[4];
    std::stringstream ss;
    ss << inv->chr << ":" << ref_bp1_start << "-" << ref_bp1_end;
    regions[0] = strdup(ss.str().c_str());
    ss.str("");
    
    ss << inv->chr << ":" << ref_bp2_start << "-" << ref_bp2_end;
    regions[1] = strdup(ss.str().c_str());
    ss.str("");

    ss << inv->chr << ":" << std::max(hts_pos_t(0), inv->inv_start-extend) << "-" << inv->inv_start+extend;
    regions[2] = strdup(ss.str().c_str());
    ss.str("");

    ss << inv->chr << ":" << std::max(hts_pos_t(0), inv->inv_end-extend) << "-" << inv->inv_end+extend;
    regions[3] = strdup(ss.str().c_str());
    ss.str("");

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 4);

    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> alt_bp1_better_reads, alt_bp2_better_reads, ref_bp1_better_reads, ref_bp2_better_reads;
    std::vector<bool> alt_bp1_better_reads_isrc, alt_bp2_better_reads_isrc, ref_bp1_better_reads_isrc, ref_bp2_better_reads_isrc;

    StripedSmithWaterman::Filter filter_with_score_only(false, false, 0, 32767);
    StripedSmithWaterman::Alignment alt1_aln, alt2_aln, ref1_aln, ref2_aln;
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        hts_pos_t read_start = get_unclipped_start(read), read_end = get_unclipped_end(read);
        if (overlap(read_start, read_end, sv_start, sv_start+1) == 0 && 
            overlap(read_start, read_end, sv_end, sv_end+1) == 0 &&
            overlap(read_start, read_end, inv->inv_start, inv->inv_start+1) == 0 &&
            overlap(read_start, read_end, inv->inv_end, inv->inv_end+1) == 0) continue;

        alt1_aln.Clear();
        alt2_aln.Clear();
        ref1_aln.Clear();
        ref2_aln.Clear();

        std::string seq = get_sequence(read);
        bool alt_is_rc = false, ref_is_rc = false;
        // if mate is outside the inversion and pointing towards it
        if (bam_is_mrev(read) && sv_end+stats.read_len/2 <= get_mate_endpos(read) ||
            !bam_is_mrev(read) && read->core.mpos <= sv_start-stats.read_len/2) {
            if (is_samestr(read)) { // and mate is in the same orientation, we need to rc
                rc(seq);
                alt_is_rc = true;
                ref_is_rc = true;
            }

            // align to ALT
            aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter_with_score_only, &alt1_aln, 0);
            aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter_with_score_only, &alt2_aln, 0);

            // align to REF
            aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter_with_score_only, &ref1_aln, 0);
            aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter_with_score_only, &ref2_aln, 0);
        } else if (sv_start-10 <= read->core.mpos && get_mate_endpos(read) <= sv_end+10) { // if mate is inside the inversion

            // if read and mate point towards each other, it means that if the inversion is true, they were RC together
            // therefore, we need to align it as it is to REF but reverse-complemented to ALT
            if (is_proper_pair(read, stats.min_is, stats.max_is)) {
                if (bam_is_mrev(read)) {
                    aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter_with_score_only, &ref1_aln, 0);
                    rc(seq);
                    aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter_with_score_only, &alt2_aln, 0);
                    alt_is_rc = true;
                } else {
                    aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter_with_score_only, &ref2_aln, 0);
                    rc(seq);
                    aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter_with_score_only, &alt1_aln, 0);
                    alt_is_rc = true;
                }
            } else if (is_samestr(read)) { // if both point in the same direction, the mate was RC by itself, and we can align the read as it is
                if (bam_is_rev(read)) {
                    aligner.Align(seq.c_str(), contig_seq+ref_bp2_start, ref_bp2_len, filter_with_score_only, &ref2_aln, 0);
                    aligner.Align(seq.c_str(), alt_bp2_seq, alt_bp2_len, filter_with_score_only, &alt2_aln, 0);
                } else {
                    aligner.Align(seq.c_str(), contig_seq+ref_bp1_start, ref_bp1_len, filter_with_score_only, &ref1_aln, 0);
                    aligner.Align(seq.c_str(), alt_bp1_seq, alt_bp1_len, filter_with_score_only, &alt1_aln, 0);
                }
            }
        }

        StripedSmithWaterman::Alignment& alt_aln = alt1_aln.sw_score >= alt2_aln.sw_score ? alt1_aln : alt2_aln;
        StripedSmithWaterman::Alignment& ref_aln = ref1_aln.sw_score >= ref2_aln.sw_score ? ref1_aln : ref2_aln;
        if (alt_aln.sw_score > ref_aln.sw_score) {
            if (alt1_aln.sw_score >= alt2_aln.sw_score) {
                alt_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                alt_bp1_better_reads_isrc.push_back(alt_is_rc);
            }
            if (alt1_aln.sw_score <= alt2_aln.sw_score) {
                alt_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                alt_bp2_better_reads_isrc.push_back(alt_is_rc);
            }
        } else if (alt_aln.sw_score < ref_aln.sw_score) {
            if (ref1_aln.sw_score >= ref2_aln.sw_score) {
                ref_bp1_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                ref_bp1_better_reads_isrc.push_back(ref_is_rc);
            } 
            if (ref1_aln.sw_score <= ref2_aln.sw_score) {
                ref_bp2_better_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
                ref_bp2_better_reads_isrc.push_back(ref_is_rc);
            }
        } else if (alt_aln.sw_score == ref_aln.sw_score && alt_aln.sw_score != 0) {
            inv->sample_info.alt_ref_equal_reads++;
            if (read->core.qual >= config.high_confidence_mapq) {
                inv->sample_info.alt_ref_equal_reads_highmq++;
            }
        }

        if (alt_bp1_better_reads.size() + alt_bp2_better_reads.size() + ref_bp1_better_reads.size() + ref_bp2_better_reads.size() + inv->sample_info.alt_ref_equal_reads > 4 * stats.get_max_depth(inv->chr)) {
            alt_bp1_better_reads.clear();
            alt_bp2_better_reads.clear();
            ref_bp1_better_reads.clear();
            ref_bp2_better_reads.clear();
            inv->sample_info.alt_ref_equal_reads = 0;
            inv->sample_info.alt_ref_equal_reads_highmq = 0;
            inv->sample_info.too_deep = true;
            break;
        }
    }

    free(regions[0]);
    free(regions[1]);
    free(regions[2]);
    free(regions[3]);
    bam_destroy1(read);
    hts_itr_destroy(iter);

    std::string alt_bp1_consensus_seq, alt_bp2_consensus_seq, ref_bp1_consensus_seq, ref_bp2_consensus_seq;
    double alt_bp1_avg_score, alt_bp2_avg_score, ref_bp1_avg_score, ref_bp2_avg_score;
    double alt_bp1_stddev_score, alt_bp2_stddev_score, ref_bp1_stddev_score, ref_bp2_stddev_score;
    std::vector<bool> alt_bp1_is_exact_read, alt_bp2_is_exact_read, ref_bp1_is_exact_read, ref_bp2_is_exact_read;
    auto alt_bp1_is_consistent_read = gen_consensus_and_classify_seqs(alt_bp1_seq, alt_bp1_better_reads, alt_bp1_better_reads_isrc, alt_bp1_consensus_seq, alt_bp1_avg_score, alt_bp1_stddev_score, alt_bp1_is_exact_read, nullptr);
    auto alt_bp2_is_consistent_read = gen_consensus_and_classify_seqs(alt_bp2_seq, alt_bp2_better_reads, alt_bp2_better_reads_isrc, alt_bp2_consensus_seq, alt_bp2_avg_score, alt_bp2_stddev_score, alt_bp2_is_exact_read, nullptr);

    char* ref_bp1_seq = new char[ref_bp1_len+1];
    strncpy(ref_bp1_seq, contig_seq+ref_bp1_start, ref_bp1_len);
    ref_bp1_seq[ref_bp1_len] = 0;
    auto ref_bp1_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp1_seq, ref_bp1_better_reads, ref_bp1_better_reads_isrc, ref_bp1_consensus_seq, ref_bp1_avg_score, ref_bp1_stddev_score, ref_bp1_is_exact_read, nullptr);

    char* ref_bp2_seq = new char[ref_bp2_len+1];
    strncpy(ref_bp2_seq, contig_seq+ref_bp2_start, ref_bp2_len);
    ref_bp2_seq[ref_bp2_len] = 0;
    inv->ref1_hp_len = longest_homopolymer_len(ref_bp1_seq, ref_bp1_len);
    inv->ref2_hp_len = longest_homopolymer_len(ref_bp2_seq, ref_bp2_len);
    inv->alt1_hp_len = longest_homopolymer_len(alt_bp1_seq, alt_bp1_len);
    inv->alt2_hp_len = longest_homopolymer_len(alt_bp2_seq, alt_bp2_len);
    auto ref_bp2_is_consistent_read = gen_consensus_and_classify_seqs(ref_bp2_seq, ref_bp2_better_reads, ref_bp2_better_reads_isrc, ref_bp2_consensus_seq, ref_bp2_avg_score, ref_bp2_stddev_score, ref_bp2_is_exact_read, nullptr);

    auto score_inv_consensus = [&](const std::string& consensus_seq, bool bp1) {
        hts_pos_t breakpoint = bp1 ? sv_start : sv_end;
        hts_pos_t target_start = std::max(hts_pos_t(0), breakpoint-hts_pos_t(consensus_seq.length()));
        hts_pos_t target_end = std::min(breakpoint+hts_pos_t(consensus_seq.length()), contig_len);
        hts_pos_t inv_border_len = std::min(hts_pos_t(consensus_seq.length()), inv->inv_end-inv->inv_start);
        int extra_len = consensus_seq.length()-inv_border_len;

        alignment_targets_t targets;
        targets.alt_len = target_end-target_start;
        targets.alt_seq = new char[targets.alt_len+1];
        if (bp1) {
            int alt_lf_len = sv_start-target_start;
            strncpy(targets.alt_seq, contig_seq+target_start, alt_lf_len);
            char* inv_suffix_rc = new char[inv_border_len+1];
            strncpy(inv_suffix_rc, contig_seq+inv->inv_end-inv_border_len, inv_border_len);
            inv_suffix_rc[inv_border_len] = 0;
            rc(inv_suffix_rc);
            strncpy(targets.alt_seq+alt_lf_len, inv_suffix_rc, inv_border_len);
            extra_len = std::min(extra_len, int(contig_len-sv_end));
            strncpy(targets.alt_seq+alt_lf_len+inv_border_len, contig_seq+sv_end, extra_len);
            delete[] inv_suffix_rc;
            targets.left_flank_end = alt_lf_len;
            targets.has_right_split = false;
        } else {
            extra_len = std::min(extra_len, int(sv_start));
            strncpy(targets.alt_seq, contig_seq+sv_start-extra_len, extra_len);
            char* inv_prefix_rc = new char[inv_border_len+1];
            strncpy(inv_prefix_rc, contig_seq+inv->inv_start, inv_border_len);
            inv_prefix_rc[inv_border_len] = 0;
            rc(inv_prefix_rc);
            strncpy(targets.alt_seq+extra_len, inv_prefix_rc, inv_border_len);
            delete[] inv_prefix_rc;
            strncpy(targets.alt_seq+extra_len+inv_border_len, contig_seq+sv_end, target_end-sv_end);
            targets.right_flank_start = extra_len+inv_border_len;
            targets.has_left_split = false;
        }
        targets.alt_seq[targets.alt_len] = 0;

        targets.ref_seqs.push_back(contig_seq+target_start);
        targets.ref_lens.push_back(target_end-target_start);
        targets.ref_starts.push_back(target_start);
        if (inv_border_len == inv->inv_end-inv->inv_start) {
            int alt_edit_begin = bp1 ? int(sv_start-target_start) : extra_len;
            targets.edits.push_back({allele_edit_kind_t::INDEL, sv_start, sv_end, alt_edit_begin, alt_edit_begin+int(inv_border_len), int(sv_end-sv_start+inv_border_len), true});
        }
        if (bp1) {
            targets.left_independent_ref_seq = contig_seq+target_start;
            targets.left_independent_ref_len = target_end-target_start;
        } else {
            targets.right_independent_ref_seq = contig_seq+target_start;
            targets.right_independent_ref_len = target_end-target_start;
        }

        consensus_alignment_metrics_t metrics = score_consensus_alignment(consensus_seq, targets, aligner);
        delete[] targets.alt_seq;
        return metrics;
    };

    inv->sample_info.ext_alt_consensus1_metrics.length = alt_bp1_consensus_seq.length();
    if (alt_bp1_consensus_seq.length() >= 2*config.min_clip_len) {
        consensus_alignment_metrics_t unextended_metrics = score_inv_consensus(alt_bp1_consensus_seq, true);
        inv->sample_info.alt_consensus1_metrics = unextended_metrics;

        // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_bp1_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_bp1_consensus_seq, std::string(alt_bp1_consensus_seq.length(), '!'), 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp1_consensus, candidate_reads_for_extension_itree, std::max<hts_pos_t>(0, sv_start-GENOTYPE_CONSENSUS_EXTENSION), sv_start, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        extend_consensus_to_right(alt_bp1_consensus, candidate_reads_for_extension_itree, sv_start, std::min<hts_pos_t>(contig_len, sv_start+GENOTYPE_CONSENSUS_EXTENSION), contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        inv->sample_info.alt_lext_reads = alt_bp1_consensus->left_ext_reads;
        inv->sample_info.alt_rext_reads = alt_bp1_consensus->right_ext_reads;
        inv->sample_info.hq_alt_lext_reads = alt_bp1_consensus->hq_left_ext_reads;
        inv->sample_info.hq_alt_rext_reads = alt_bp1_consensus->hq_right_ext_reads;
        alt_bp1_consensus_seq = alt_bp1_consensus->sequence;

        consensus_alignment_metrics_t extended_metrics = score_inv_consensus(alt_bp1_consensus_seq, true);
        inv->sample_info.ext_alt_consensus1_metrics.length = extended_metrics.length;
        inv->sample_info.ext_alt_consensus1_metrics.alt_score = extended_metrics.alt_score;
        inv->sample_info.ext_alt_consensus1_metrics.ref_score = extended_metrics.ref_score;
        inv->sample_info.ext_alt_consensus1_metrics.covered_edit_distance = extended_metrics.covered_edit_distance;
        inv->sample_info.ext_alt_consensus1_metrics.main_edit_covered = extended_metrics.main_edit_covered;

        hts_pos_t target_start = std::max(hts_pos_t(0), sv_start-hts_pos_t(alt_bp1_consensus_seq.length()));
        int lf_aln_rlen = std::max(0, int(sv_start-target_start)-extended_metrics.alt_ref_begin);
        inv->left_anchor_aln->start = sv_start-lf_aln_rlen;
        inv->left_anchor_aln->end = sv_start;
        inv->left_anchor_aln->seq_len = lf_aln_rlen;
    }

    inv->sample_info.ext_alt_consensus2_metrics.length = alt_bp2_consensus_seq.length();
    if (alt_bp2_consensus_seq.length() >= 2*config.min_clip_len) {
        consensus_alignment_metrics_t unextended_metrics = score_inv_consensus(alt_bp2_consensus_seq, false);
        inv->sample_info.alt_consensus2_metrics = unextended_metrics;

        // all we care about is the consensus sequence
        std::shared_ptr<consensus_t> alt_bp2_consensus = std::make_shared<consensus_t>(false, 0, 0, 0, alt_bp2_consensus_seq, std::string(alt_bp2_consensus_seq.length(), '!'), 0, 0, 0, 0, 0, 0);
        extend_consensus_to_left(alt_bp2_consensus, candidate_reads_for_extension_itree, std::max<hts_pos_t>(0, sv_end-GENOTYPE_CONSENSUS_EXTENSION), sv_end, contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        extend_consensus_to_right(alt_bp2_consensus, candidate_reads_for_extension_itree, sv_end, std::min<hts_pos_t>(contig_len, sv_end+GENOTYPE_CONSENSUS_EXTENSION), contig_len, config.high_confidence_mapq, stats, mateseqs_w_mapq_chr, GENOTYPE_CONSENSUS_EXTENSION);
        inv->sample_info.alt_lext_reads += alt_bp2_consensus->left_ext_reads;
        inv->sample_info.alt_rext_reads += alt_bp2_consensus->right_ext_reads;
        inv->sample_info.hq_alt_lext_reads += alt_bp2_consensus->hq_left_ext_reads;
        inv->sample_info.hq_alt_rext_reads += alt_bp2_consensus->hq_right_ext_reads;
        alt_bp2_consensus_seq = alt_bp2_consensus->sequence;

        consensus_alignment_metrics_t extended_metrics = score_inv_consensus(alt_bp2_consensus_seq, false);
        inv->sample_info.ext_alt_consensus2_metrics.length = extended_metrics.length;
        inv->sample_info.ext_alt_consensus2_metrics.alt_score = extended_metrics.alt_score;
        inv->sample_info.ext_alt_consensus2_metrics.ref_score = extended_metrics.ref_score;
        inv->sample_info.ext_alt_consensus2_metrics.covered_edit_distance = extended_metrics.covered_edit_distance;
        inv->sample_info.ext_alt_consensus2_metrics.main_edit_covered = extended_metrics.main_edit_covered;

        hts_pos_t target_start = std::max(hts_pos_t(0), sv_end-hts_pos_t(alt_bp2_consensus_seq.length()));
        hts_pos_t target_end = std::min(sv_end+hts_pos_t(alt_bp2_consensus_seq.length()), contig_len);
        int target_len = target_end-target_start;
        int rf_aln_rlen = std::max(0, int(target_end-sv_end)-(target_len-extended_metrics.alt_ref_end));
        inv->right_anchor_aln->start = sv_end;
        inv->right_anchor_aln->end = sv_end+rf_aln_rlen;
        inv->right_anchor_aln->seq_len = rf_aln_rlen;
    }

    delete[] alt_bp1_seq;
    delete[] alt_bp2_seq;
    delete[] ref_bp1_seq;
    delete[] ref_bp2_seq;

    set_bp_consensus_info(inv->sample_info.alt_bp1.reads_info, alt_bp1_better_reads, alt_bp1_is_consistent_read, alt_bp1_is_exact_read, alt_bp1_avg_score, alt_bp1_stddev_score);
    set_bp_consensus_info(inv->sample_info.alt_bp2.reads_info, alt_bp2_better_reads, alt_bp2_is_consistent_read, alt_bp2_is_exact_read, alt_bp2_avg_score, alt_bp2_stddev_score);
    set_bp_consensus_info(inv->sample_info.ref_bp1.reads_info, ref_bp1_better_reads, ref_bp1_is_consistent_read, ref_bp1_is_exact_read, ref_bp1_avg_score, ref_bp1_stddev_score);
    set_bp_consensus_info(inv->sample_info.ref_bp2.reads_info, ref_bp2_better_reads, ref_bp2_is_consistent_read, ref_bp2_is_exact_read, ref_bp2_avg_score, ref_bp2_stddev_score);
}

bool is_small_inv(inversion_t* inv, stats_t& stats, config_t& config) {
    return inv->end-inv->start+inv->svlen() < stats.read_len-2*config.min_clip_len;
}

void genotype_invs(int id, std::string contig_name, char* contig_seq, int contig_len, std::vector<inversion_t*> invs,
    bcf_hdr_t* in_vcf_header, bcf_hdr_t* out_vcf_header, stats_t& stats, config_t& config, contig_map_t& contig_map,
    bam_pool_t* bam_pool) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);

    int contig_id = contig_map.get_id(contig_name);
    auto chromosome_data = acquire_chromosome_data(contig_id);
    auto mateseqs_w_mapq_chr = chromosome_data.first;

    std::vector<hts_pair_pos_t> target_ivals;
    for (inversion_t* inv : invs) {
        target_ivals.push_back({std::max<hts_pos_t>(0, inv->start-GENOTYPE_CONSENSUS_EXTENSION), std::min<hts_pos_t>(contig_len, inv->start+GENOTYPE_CONSENSUS_EXTENSION)});
        target_ivals.push_back({std::max<hts_pos_t>(0, inv->end-GENOTYPE_CONSENSUS_EXTENSION), std::min<hts_pos_t>(contig_len, inv->end+GENOTYPE_CONSENSUS_EXTENSION)});
    }
    std::vector<ext_read_t*> candidate_reads_for_extension;
    IntervalTree<ext_read_t*> candidate_reads_for_extension_itree = get_candidate_reads_for_extension_itree(contig_name, contig_len, target_ivals, bam_pool->get_bam_reader(id), candidate_reads_for_extension);

    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);
    for (inversion_t* inv : invs) {
        if (is_small_inv(inv, stats, config)) {
            genotype_small_inv(inv, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner);
        } else {
            genotype_large_inv(inv, bam_file, candidate_reads_for_extension_itree, *mateseqs_w_mapq_chr, contig_seq, contig_len, stats, config, aligner);
        }
    }

    for (ext_read_t* ext_read : candidate_reads_for_extension) delete ext_read;

    calculate_ptn_ratio(contig_name, invs, bam_file, config, stats);
    depth_filter_inv(contig_name, invs, bam_file, config, stats);

    // apply hard filters to inversions and calculate GT
    for (inversion_t* inv : invs) {
        bool fail_dp = false, fail_sr = false;
        inv->sample_info.filters.clear();
        if (inv->sample_info.left_flanking_cov > stats.get_max_depth(inv->chr) || inv->sample_info.right_flanking_cov > stats.get_max_depth(inv->chr) ||
            inv->sample_info.left_flanking_cov < stats.get_min_depth(inv->chr) || inv->sample_info.right_flanking_cov < stats.get_min_depth(inv->chr) ||
            inv->sample_info.left_anchor_cov > stats.get_max_depth(inv->chr) || inv->sample_info.right_anchor_cov > stats.get_max_depth(inv->chr)) {
            inv->sample_info.filters.push_back("ANOMALOUS_FLANKING_DEPTH");
            fail_dp = fail_sr = true;
        }
        if (inv->sample_info.alt_bp1.reads_info.consistent_reads() > stats.get_max_depth(inv->chr) 
            || inv->sample_info.alt_bp2.reads_info.consistent_reads() > stats.get_max_depth(inv->chr)) {
            inv->sample_info.filters.push_back("ANOMALOUS_SC_NUMBER");
            fail_dp = fail_sr = true;
        }

        if (inv->sample_info.alt_bp1.pairs_info.pos_max_mq < config.high_confidence_mapq || inv->sample_info.alt_bp2.pairs_info.pos_max_mq < config.high_confidence_mapq ||
            inv->sample_info.alt_bp1.pairs_info.neg_max_mq < config.high_confidence_mapq || inv->sample_info.alt_bp2.pairs_info.neg_max_mq < config.high_confidence_mapq) {
            inv->sample_info.filters.push_back("LOW_MAPQ_DISC_PAIRS");
            fail_dp = true;
        }
        if (inv->sample_info.alt_bp1.pairs_info.pairs < stats.get_min_disc_pairs_by_insertion_size(inv->svlen())/2 ||
            inv->sample_info.alt_bp2.pairs_info.pairs < stats.get_min_disc_pairs_by_insertion_size(inv->svlen())/2) {
            inv->sample_info.filters.push_back("NOT_ENOUGH_DISC_PAIRS");
            fail_dp = true;
        }
        double ptn_ratio_bp1 = double(inv->sample_info.alt_bp1.pairs_info.pairs)/std::max(1, inv->sample_info.alt_bp1.pairs_info.pairs+inv->sample_info.ref_bp1.pairs_info.pairs);
        double ptn_ratio_bp2 = double(inv->sample_info.alt_bp2.pairs_info.pairs)/std::max(1, inv->sample_info.alt_bp2.pairs_info.pairs+inv->sample_info.ref_bp2.pairs_info.pairs);
        if (ptn_ratio_bp1 < 0.25 || ptn_ratio_bp2 < 0.25) {
            inv->sample_info.filters.push_back("LOW_PTN_RATIO");
            fail_dp = true;
        }

        if (inv->left_anchor_aln->end-inv->left_anchor_aln->start < stats.max_is/2 || inv->right_anchor_aln->end-inv->right_anchor_aln->start < stats.max_is/2) {
            inv->sample_info.filters.push_back("SHORT_ANCHOR");
            fail_sr = true;
        }
        if (inv->sample_info.ext_alt_consensus1_metrics.alt_score <= inv->sample_info.ext_alt_consensus1_metrics.ref_score || 
            !is_small_inv(inv, stats, config) && inv->sample_info.ext_alt_consensus2_metrics.alt_score <= inv->sample_info.ext_alt_consensus2_metrics.ref_score) {
            inv->sample_info.filters.push_back("LOW_ALT_CONSENSUS_SCORE");
            fail_sr = true;
        }

        if (!fail_sr || !fail_dp) {
            inv->sample_info.filters.clear();
            inv->sample_info.filters.push_back("PASS");
        }

        double pairs_ptn = fail_dp ? 0 : std::min(ptn_ratio_bp1, ptn_ratio_bp2);

        double reads_ptn1 = double(inv->sample_info.alt_bp1.reads_info.consistent_reads())/std::max(1, inv->sample_info.alt_bp1.reads_info.consistent_reads()+inv->sample_info.ref_bp1.reads_info.consistent_reads());
        double reads_ptn2 = reads_ptn1;
        if (inv->sample_info.alt_bp2.reads_info.computed) {
            reads_ptn2 = double(inv->sample_info.alt_bp2.reads_info.consistent_reads())/std::max(1, inv->sample_info.alt_bp2.reads_info.consistent_reads()+inv->sample_info.ref_bp2.reads_info.consistent_reads());
        }
        double reads_ptn = fail_sr ? 0 : std::min(reads_ptn1, reads_ptn2);

        double ptn = std::max(pairs_ptn, reads_ptn);
        if (ptn >= 0.75) {
            inv->sample_info.gt = {bcf_gt_unphased(1), bcf_gt_unphased(1)};
        } else if (ptn <= 0.25) {
            inv->sample_info.gt = {bcf_gt_unphased(0), bcf_gt_unphased(0)};
        } else {
            inv->sample_info.gt = {bcf_gt_unphased(0), bcf_gt_unphased(1)};
        }
    }

    release_chromosome_data(contig_id);
}

#endif // GENOTYPE_INVS_H
