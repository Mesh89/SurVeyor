#include "genotype.h"
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <new>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "assemble.h"
#include "htslib/hts_endian.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/tbx.h"

#include "sw_utils.h"
#include "types.h"
#include "utils.h"
#include "sam_utils.h"
#include "extend_1sr_consensus.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "vcf_utils.h"
#include "stat_tests.h"
#include "reference_guided_assembly.h"
#include "SegTree.h"
#include "../libs/IntervalTree.h"

#include "genotype_dels.h"
#include "genotype_dups.h"
#include "genotype_inss.h"
#include "genotype_invs.h"
#include "genotype_hp_indels.h"

chr_seqs_map_t chr_seqs;
config_t config;
contig_map_t contig_map;
stats_t stats;

std::string bam_fname, reference_fname, workdir;
bam_pool_t* bam_pool;

std::vector<hts_pos_t> global_isize_dist;

std::vector<double> global_crossing_isize_dist;

StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
StripedSmithWaterman::Aligner harsh_aligner(1, 4, 100, 1, false);

std::vector<std::unordered_map<std::string, std::pair<std::string, int> > > mateseqs_w_mapq;
std::vector<int> active_threads_per_chr;
std::vector<std::mutex> mutex_per_chr;
std::vector<std::vector<sv_t*>> evidence_svs_by_contig;
std::vector<std::unique_ptr<evidence_map_t>> evidence_maps_by_contig;
std::vector<bool> evidence_map_load_attempted;
evidence_map_t empty_evidence_map;
bool load_evidence_maps = false;
evidence_mode_t evidence_mode = evidence_mode_t::CACHED;
std::string reads_association_dir;

void update_record_bp_reads_info(bcf_hdr_t* out_hdr, bcf1_t* b, sv_t::bp_reads_info_t bp_reads_info, std::string prefix, int bp_number) {
    if (!bp_reads_info.computed) return;

    std::string bp_number_str = std::to_string(bp_number);
    std::string read_fmt_prefix = prefix + "R" + bp_number_str;

    bcf_update_format_int32(out_hdr, b, read_fmt_prefix.c_str(), &(bp_reads_info.reads), 1);

    int consistent_reads = bp_reads_info.consistent_reads();
    int consistent_exact_reads = bp_reads_info.consistent_exact_reads();
    int exact_reads = bp_reads_info.exact_reads();
    if (bp_reads_info.reads) {
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "C").c_str(), &(consistent_reads), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "E").c_str(), &exact_reads, 1);
    }
    if (consistent_reads) {
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CE").c_str(), &consistent_exact_reads, 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CF").c_str(), &(bp_reads_info.consistent_fwd), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CR").c_str(), &(bp_reads_info.consistent_rev), 1);
        float rcas = bp_reads_info.consistent_avg_score;
        bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CAS").c_str(), &rcas, 1);
        float rcss = bp_reads_info.consistent_stddev_score;
        bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CSS").c_str(), &rcss, 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CHQ").c_str(), &(bp_reads_info.consistent_high_mq), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CmQ").c_str(), &(bp_reads_info.consistent_min_mq), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMQ").c_str(), &(bp_reads_info.consistent_max_mq), 1);
        float rcavmq = bp_reads_info.consistent_avg_mq;
        bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CAQ").c_str(), &rcavmq, 1);
        float rcstdmq = bp_reads_info.consistent_stddev_mq;
        bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CSQ").c_str(), &rcstdmq, 1);
        int mate_cov_bps[] = {bp_reads_info.fwd_mate_cov_bps, bp_reads_info.rev_mate_cov_bps};
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMSPAN").c_str(), mate_cov_bps, 2);
        int hq_mate_cov_bps[] = {bp_reads_info.fwd_hq_mate_cov_bps, bp_reads_info.rev_hq_mate_cov_bps};
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMHQSPAN").c_str(), hq_mate_cov_bps, 2);

        if (consistent_exact_reads) {
            bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CEHQ").c_str(), &(bp_reads_info.consistent_exact_high_mq), 1);
            bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CEmQ").c_str(), &(bp_reads_info.consistent_exact_min_mq), 1);
            bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CEMQ").c_str(), &(bp_reads_info.consistent_exact_max_mq), 1);
            bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CEF").c_str(), &(bp_reads_info.consistent_exact_fwd), 1);
            bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CER").c_str(), &(bp_reads_info.consistent_exact_rev), 1);
        }
    }
    if (exact_reads) {
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "EHQ").c_str(), &(bp_reads_info.exact_high_mq), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "EmQ").c_str(), &(bp_reads_info.exact_min_mq), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "EMQ").c_str(), &(bp_reads_info.exact_max_mq), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "EF").c_str(), &(bp_reads_info.exact_fwd), 1);
        bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "ER").c_str(), &(bp_reads_info.exact_rev), 1);
    }
}

void reset_record_bp_reads_info(bcf_hdr_t* out_hdr, bcf1_t* b, std::string prefix, int bp_number) {
    std::string bp_number_str = std::to_string(bp_number);
    std::string read_fmt_prefix = prefix + "R" + bp_number_str;

    bcf_update_format_int32(out_hdr, b, read_fmt_prefix.c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "C").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CF").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CR").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CE").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CEHQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CEmQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CEMQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CEF").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CER").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "E").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "EHQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "EmQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "EMQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "EF").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "ER").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CAS").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CSS").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CHQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CmQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMQ").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CAQ").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (read_fmt_prefix + "CSQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMSPAN").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (read_fmt_prefix + "CMHQSPAN").c_str(), NULL, 0);
}

void update_record_bp_pairs_info(bcf_hdr_t* out_hdr, bcf1_t* b, sv_t::bp_pairs_info_t bp_pairs_info, std::string prefix, int bp_number) {
    if (!bp_pairs_info.computed) return;

    std::string bp_number_str = std::to_string(bp_number);
    std::string pairs_fmt_prefix = prefix + "SP" + bp_number_str;

    bcf_update_format_int32(out_hdr, b, pairs_fmt_prefix.c_str(), &(bp_pairs_info.pairs), 1);

    if (bp_pairs_info.pairs) {
        int pairs_hq[] = {bp_pairs_info.pos_high_mapq, bp_pairs_info.neg_high_mapq};
        bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "HQ").c_str(), pairs_hq, 2);
        int pairs_min_mq[] = {bp_pairs_info.pos_min_mq, bp_pairs_info.neg_min_mq};
        bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "mQ").c_str(), pairs_min_mq, 2);
        int pairs_max_mq[] = {bp_pairs_info.pos_max_mq, bp_pairs_info.neg_max_mq};
        bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "MQ").c_str(), pairs_max_mq, 2);
        float pairs_avg_mq[] = {(float) bp_pairs_info.pos_avg_mq, (float) bp_pairs_info.neg_avg_mq};
        bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "AQ").c_str(), pairs_avg_mq, 2);
        float pairs_stddev_mq[] = {(float) bp_pairs_info.pos_stddev_mq, (float) bp_pairs_info.neg_stddev_mq};
        bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "SQ").c_str(), pairs_stddev_mq, 2);
        int pairs_span[] = {bp_pairs_info.lf_span, bp_pairs_info.rf_span};
        bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "SPAN").c_str(), pairs_span, 2);
        float pairs_avg_nm[] = {(float) bp_pairs_info.pos_avg_nm, (float) bp_pairs_info.neg_avg_nm};
        bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "NMA").c_str(), pairs_avg_nm, 2);
        float pairs_stddev_nm[] = {(float) bp_pairs_info.pos_stddev_nm, (float) bp_pairs_info.neg_stddev_nm};
        bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "NMS").c_str(), pairs_stddev_nm, 2);
    }
}

void reset_record_bp_pairs_info(bcf_hdr_t* out_hdr, bcf1_t* b, std::string prefix, int bp_number) {
    std::string bp_number_str = std::to_string(bp_number);
    std::string pairs_fmt_prefix = prefix + "SP" + bp_number_str;

    bcf_update_format_int32(out_hdr, b, pairs_fmt_prefix.c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "HQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "mQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "MQ").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "AQ").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "SQ").c_str(), NULL, 0);
    bcf_update_format_int32(out_hdr, b, (pairs_fmt_prefix + "SPAN").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "NMA").c_str(), NULL, 0);
    bcf_update_format_float(out_hdr, b, (pairs_fmt_prefix + "NMS").c_str(), NULL, 0);
}

void update_record_bp_consensus_info(bcf_hdr_t* out_hdr, bcf1_t* b, sv_t::bp_consensus_info_t bp_consensus_info, std::string prefix, int bp_number) {
    update_record_bp_reads_info(out_hdr, b, bp_consensus_info.reads_info, prefix, bp_number);
    update_record_bp_pairs_info(out_hdr, b, bp_consensus_info.pairs_info, prefix, bp_number);
}

void reset_record_bp_consensus_info(bcf_hdr_t* out_hdr, bcf1_t* b, std::string prefix, int bp_number) {
    reset_record_bp_reads_info(out_hdr, b, prefix, bp_number);
    reset_record_bp_pairs_info(out_hdr, b, prefix, bp_number);
}

int calculate_mh_len(sv_t* sv) {
    hts_pos_t chr_len = chr_seqs.get_len(sv->chr);
    char* chr_seq = chr_seqs.get_seq(sv->chr);

    // if complex indel, we do not consider microhomology
    if ((sv->svtype() == "DEL" || sv->svtype() == "DUP") && !sv->ins_seq.empty()) {
        return 0;
    } else if (sv->svtype() == "INS" && sv->start != sv->end) {
        return 0;
    }

    // cap microhomology length at 1000 bp to avoid excessive runtime
    int rf_len_cap = std::min(hts_pos_t(1000), sv->svsize());
    char* right_flanking = generate_haplotype_right(chr_seq, chr_len, sv->end+1, rf_len_cap, sv->aux_indels, sv->aux_snps);
    int rf_len = strlen(right_flanking);

    int mh_len = 0;
    if (sv->svtype() == "DEL" || sv->svtype() == "DUP") {
        while (mh_len < rf_len && toupper(chr_seq[sv->start+mh_len+1]) == toupper(right_flanking[mh_len])) {
            mh_len++;
        }
    } else if (sv->svtype() == "INS") {
        while (mh_len < sv->ins_seq.length() && sv->end+mh_len+1 < chr_len && toupper(sv->ins_seq[mh_len]) == toupper(chr_seq[sv->end+mh_len+1])) {
            mh_len++;
        }
    }

    delete[] right_flanking;

    return mh_len;
}

bool genotype_when_reassigning_evidence(const sv_t* sv) {
    return config.training_mode || sv->sample_info.epr >= MIN_EPR;
}

void update_other_read_support_fields(bcf_hdr_t* out_hdr, bcf1_t* record, const std::string& prefix, int bp_number, bool computed, int reads, int consistent, int consistent_hq, int exact) {
    const std::string tag = prefix + std::to_string(bp_number);
    const std::string consistent_tag = tag + "C";
    const std::string consistent_hq_tag = tag + "CHQ";
    const std::string exact_tag = tag + "E";

    if (!computed) {
        bcf_update_format_int32(out_hdr, record, tag.c_str(), NULL, 0);
        if (prefix == "OAR") bcf_update_format_int32(out_hdr, record, (tag + "ALL").c_str(), NULL, 0);
        bcf_update_format_int32(out_hdr, record, consistent_tag.c_str(), NULL, 0);
        bcf_update_format_int32(out_hdr, record, consistent_hq_tag.c_str(), NULL, 0);
        bcf_update_format_int32(out_hdr, record, exact_tag.c_str(), NULL, 0);
        return;
    }

    bcf_update_format_int32(out_hdr, record, tag.c_str(), &reads, 1);
    if (prefix == "OAR") bcf_update_format_int32(out_hdr, record, (tag + "ALL").c_str(), &reads, 1);
    bcf_update_format_int32(out_hdr, record, consistent_tag.c_str(), &consistent, 1);
    bcf_update_format_int32(out_hdr, record, consistent_hq_tag.c_str(), &consistent_hq, 1);
    bcf_update_format_int32(out_hdr, record, exact_tag.c_str(), &exact, 1);
}

void update_oar_max_fields(bcf_hdr_t* out_hdr, bcf1_t* record, int bp_number, bool computed,
    const std::unordered_map<int, int>& reads_by_hpid,
    const std::unordered_map<std::string, int>& reads_by_vid) {
    const std::string tag_prefix = "OAR" + std::to_string(bp_number);

    if (!computed) {
        bcf_update_format_int32(out_hdr, record, (tag_prefix + "HPID").c_str(), NULL, 0);
        bcf_update_format_string(out_hdr, record, (tag_prefix + "VID").c_str(), NULL, 0);
        bcf_update_format_int32(out_hdr, record, (tag_prefix + "MAX").c_str(), NULL, 0);
        return;
    }

    int32_t best_hpid = bcf_int32_missing;
    int32_t best_hpid_count = 0;

    for (const auto& entry : reads_by_hpid) {
        int hpid = entry.first;
        int count = entry.second;
        if (count > best_hpid_count || (count == best_hpid_count && hpid < best_hpid)) {
            best_hpid = hpid;
            best_hpid_count = count;
        }
    }

    if (reads_by_hpid.empty()) {
        bcf_update_format_int32(out_hdr, record, (tag_prefix + "HPID").c_str(), NULL, 0);
    } else {
        bcf_update_format_int32(out_hdr, record, (tag_prefix + "HPID").c_str(), &best_hpid, 1);
    }

    int32_t best_vid_count = 0;
    std::vector<std::string> best_vids;
    for (const auto& entry : reads_by_vid) {
        const std::string& vid = entry.first;
        int count = entry.second;
        if (count > best_vid_count) {
            best_vid_count = count;
            best_vids.clear();
            best_vids.push_back(vid);
        } else if (count == best_vid_count) {
            best_vids.push_back(vid);
        }
    }

    if (best_vids.empty()) {
        bcf_update_format_string(out_hdr, record, (tag_prefix + "VID").c_str(), NULL, 0);
    } else {
        std::sort(best_vids.begin(), best_vids.end());
        std::string joined_vids = best_vids[0];
        for (size_t i = 1; i < best_vids.size(); i++) {
            joined_vids += "," + best_vids[i];
        }
        const char* joined_vids_cstr = joined_vids.c_str();
        bcf_update_format_string(out_hdr, record, (tag_prefix + "VID").c_str(), &joined_vids_cstr, 1);
    }
    bcf_update_format_int32(out_hdr, record, (tag_prefix + "MAX").c_str(), &best_vid_count, 1);
}

void update_record(bcf_hdr_t* in_hdr, bcf_hdr_t* out_hdr, sv_t* sv, char* chr_seq, hts_pos_t chr_len, int sample_idx) {
    
    bcf_translate(out_hdr, in_hdr, sv->vcf_entry);

    bcf_subset(out_hdr, sv->vcf_entry, 1, &sample_idx);

    // update INFO fields
    int sv_end = sv->end+1;
    bcf_update_info_int32(out_hdr, sv->vcf_entry, "END", &sv_end, 1);

    if (sv->ins_seq.find("-") == std::string::npos) {
        if (!sv->ins_seq.empty()) {
            int svinslen = sv->ins_seq.length();
            bcf_update_info_int32(out_hdr, sv->vcf_entry, "SVINSLEN", &svinslen, 1);
        }
        if (sv->svtype() != "DUP") {
            int svlen = sv->svlen();
            bcf_update_info_int32(out_hdr, sv->vcf_entry, "SVLEN", &svlen, 1);
        } else {
            bcf_update_info_int32(out_hdr, sv->vcf_entry, "SVLEN", NULL, 0);
        }
    }
    hts_pos_t la_aln_start = bring_within_range(sv->left_anchor_aln->start, hts_pos_t(0), chr_len);
    hts_pos_t la_aln_end = bring_within_range(sv->left_anchor_aln->end, hts_pos_t(0), chr_len);
    hts_pos_t ra_aln_start = bring_within_range(sv->right_anchor_aln->start, hts_pos_t(0), chr_len);
    hts_pos_t ra_aln_end = bring_within_range(sv->right_anchor_aln->end, hts_pos_t(0), chr_len);
    
    base_frequencies_t left_anchor_base_freqs = get_base_frequencies(chr_seq+la_aln_start, la_aln_end-la_aln_start);
	int labc[] = {left_anchor_base_freqs.a, left_anchor_base_freqs.c, left_anchor_base_freqs.g, left_anchor_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "LEFT_ANCHOR_BASE_COUNT", labc, 4);
    
	std::tuple<base_frequencies_t, base_frequencies_t, base_frequencies_t> left_flanking_bfs = get_left_flanking_base_frequencies_50_100_500(chr_seq, sv->start);
	int labc50[] = {std::get<0>(left_flanking_bfs).a, std::get<0>(left_flanking_bfs).c, std::get<0>(left_flanking_bfs).g, std::get<0>(left_flanking_bfs).t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "LEFT_FLANKING_BASE_COUNT_50", labc50, 4);
	int labc100[] = {std::get<1>(left_flanking_bfs).a, std::get<1>(left_flanking_bfs).c, std::get<1>(left_flanking_bfs).g, std::get<1>(left_flanking_bfs).t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "LEFT_FLANKING_BASE_COUNT_100", labc100, 4);
	int labc500[] = {std::get<2>(left_flanking_bfs).a, std::get<2>(left_flanking_bfs).c, std::get<2>(left_flanking_bfs).g, std::get<2>(left_flanking_bfs).t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "LEFT_FLANKING_BASE_COUNT_500", labc500, 4);

	base_frequencies_t right_anchor_base_freqs = get_base_frequencies(chr_seq+ra_aln_start, ra_aln_end-ra_aln_start);
	int rabc[] = {right_anchor_base_freqs.a, right_anchor_base_freqs.c, right_anchor_base_freqs.g, right_anchor_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "RIGHT_ANCHOR_BASE_COUNT", rabc, 4);

	std::tuple<base_frequencies_t, base_frequencies_t, base_frequencies_t> right_flanking_bfs = get_right_flanking_base_frequencies_50_100_500(chr_seq, sv->end, chr_len);
	int rabc50[] = {std::get<0>(right_flanking_bfs).a, std::get<0>(right_flanking_bfs).c, std::get<0>(right_flanking_bfs).g, std::get<0>(right_flanking_bfs).t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "RIGHT_FLANKING_BASE_COUNT_50", rabc50, 4);
	int rabc100[] = {std::get<1>(right_flanking_bfs).a, std::get<1>(right_flanking_bfs).c, std::get<1>(right_flanking_bfs).g, std::get<1>(right_flanking_bfs).t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "RIGHT_FLANKING_BASE_COUNT_100", rabc100, 4);
	int rabc500[] = {std::get<2>(right_flanking_bfs).a, std::get<2>(right_flanking_bfs).c, std::get<2>(right_flanking_bfs).g, std::get<2>(right_flanking_bfs).t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "RIGHT_FLANKING_BASE_COUNT_500", rabc500, 4);

	base_frequencies_t prefix_ref_base_freqs = get_base_frequencies(chr_seq+sv->start, std::min(sv->end-sv->start, hts_pos_t(5000)));
	int svrefpbc[] = {prefix_ref_base_freqs.a, prefix_ref_base_freqs.c, prefix_ref_base_freqs.g, prefix_ref_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "SV_REF_PREFIX_BASE_COUNT", svrefpbc, 4);

	base_frequencies_t suffix_ref_base_freqs = get_base_frequencies(chr_seq+sv->end-std::min(sv->end-sv->start, hts_pos_t(5000)), std::min(sv->end-sv->start, hts_pos_t(5000)));
	int svrefsbc[] = {suffix_ref_base_freqs.a, suffix_ref_base_freqs.c, suffix_ref_base_freqs.g, suffix_ref_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "SV_REF_SUFFIX_BASE_COUNT", svrefsbc, 4);

    int d = sv->ins_seq.find("-");
	std::string ins_seq_fh = sv->ins_seq.substr(0, d);
	std::string ins_seq_sh = sv->ins_seq.substr(d+1);
	base_frequencies_t prefix_base_freqs = get_base_frequencies(ins_seq_fh.c_str(), ins_seq_fh.length());
	base_frequencies_t suffix_base_freqs = ins_seq_sh.empty() ? prefix_base_freqs : get_base_frequencies(ins_seq_sh.c_str(), ins_seq_sh.length());
	int pbc[] = {prefix_base_freqs.a, prefix_base_freqs.c, prefix_base_freqs.g, prefix_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "INS_PREFIX_BASE_COUNT", pbc, 4);
	int sbc[] = {suffix_base_freqs.a, suffix_base_freqs.c, suffix_base_freqs.g, suffix_base_freqs.t};
	bcf_update_info_int32(out_hdr, sv->vcf_entry, "INS_SUFFIX_BASE_COUNT", sbc, 4);

    int mh_len = calculate_mh_len(sv);
    if (mh_len > 0) {
        bcf_update_info_int32(out_hdr, sv->vcf_entry, "MH_LEN", &mh_len, 1);
    }
    bcf_update_info_flag(out_hdr, sv->vcf_entry, "HP_GENOTYPED", "", sv->hp_genotyped);
    if (sv->hp_ref_beg != HTS_POS_MIN && sv->hp_ref_end != HTS_POS_MIN) {
        int hp_ref_range[] = {(int) sv->hp_ref_beg + 1, (int) sv->hp_ref_end + 1}; // convert to 1-based coordinates for VCF
        bcf_update_info_int32(out_hdr, sv->vcf_entry, "HP_REF_RANGE", hp_ref_range, 2);
    } else {
        bcf_update_info_int32(out_hdr, sv->vcf_entry, "HP_REF_RANGE", NULL, 0);
    }
    // update FORMAT fields
    bcf_update_genotypes(out_hdr, sv->vcf_entry, sv->sample_info.gt.data(), sv->sample_info.gt.size());

    if (sv->sample_info.expected_alt1_reads_frac != sv_t::EXPECTED_ALT_READS_FREQ_NOT_COMPUTED) {
        float earf[] = {(float) sv->sample_info.expected_alt1_reads_frac, 0.0f};
        if (sv->sample_info.expected_alt2_reads_frac == sv_t::EXPECTED_ALT_READS_FREQ_NOT_COMPUTED) {
            bcf_float_set_missing(earf[1]);
        } else {
            earf[1] = (float) sv->sample_info.expected_alt2_reads_frac;
        }
        bcf_update_format_float(out_hdr, sv->vcf_entry, "EARF", earf, 2);
    }
    const bool mfal1_computed = sv->sample_info.max_feasible_alt1_len != sv_t::sample_info_t::NOT_COMPUTED;
    const bool mfal2_computed = sv->sample_info.max_feasible_alt2_len != sv_t::sample_info_t::NOT_COMPUTED;
    if (mfal1_computed || mfal2_computed) {
        int32_t mfal[] = {bcf_int32_missing, bcf_int32_missing};
        if (mfal1_computed) mfal[0] = sv->sample_info.max_feasible_alt1_len;
        if (mfal2_computed) mfal[1] = sv->sample_info.max_feasible_alt2_len;
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MFAL", mfal, 2);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MFAL", NULL, 0);
    }

    reset_record_bp_consensus_info(out_hdr, sv->vcf_entry, "A", 1);
    reset_record_bp_consensus_info(out_hdr, sv->vcf_entry, "A", 2);
    reset_record_bp_consensus_info(out_hdr, sv->vcf_entry, "R", 1);
    reset_record_bp_consensus_info(out_hdr, sv->vcf_entry, "R", 2);
    reset_record_bp_pairs_info(out_hdr, sv->vcf_entry, "N", 1);
    reset_record_bp_pairs_info(out_hdr, sv->vcf_entry, "N", 2);
    reset_record_bp_pairs_info(out_hdr, sv->vcf_entry, "S", 1);
    reset_record_bp_pairs_info(out_hdr, sv->vcf_entry, "S", 2);

    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.alt_bp1, "A", 1);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.alt_bp2, "A", 2);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.ref_bp1, "R", 1);
    update_record_bp_consensus_info(out_hdr, sv->vcf_entry, sv->sample_info.ref_bp2, "R", 2);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.neutral_bp1_pairs, "N", 1);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.neutral_bp2_pairs, "N", 2);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.bp1_stray_pairs, "S", 1);
    update_record_bp_pairs_info(out_hdr, sv->vcf_entry, sv->sample_info.bp2_stray_pairs, "S", 2);

    int er = sv->sample_info.alt_ref_equal_reads;
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ER", &er, 1);

    int erhq = sv->sample_info.alt_ref_equal_reads_highmq;
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "ERHQ", &erhq, 1);

    update_other_read_support_fields(out_hdr, sv->vcf_entry, "OAR", 1, true, sv->sample_info.oar_bp1_reads, sv->sample_info.oar_bp1_consistent_reads, sv->sample_info.oar_bp1_consistent_hq_reads, sv->sample_info.oar_bp1_exact_reads);
    update_oar_max_fields(out_hdr, sv->vcf_entry, 1, true,
        sv->sample_info.oar_bp1_reads_by_hpid, sv->sample_info.oar_bp1_reads_by_vid);
    update_other_read_support_fields(out_hdr, sv->vcf_entry, "ORR", 1, true, sv->sample_info.orr_bp1_reads, sv->sample_info.orr_bp1_consistent_reads, sv->sample_info.orr_bp1_consistent_hq_reads, sv->sample_info.orr_bp1_exact_reads);

    const bool bp2_computed = sv->sample_info.ref_bp2.reads_info.computed;
    update_other_read_support_fields(out_hdr, sv->vcf_entry, "OAR", 2, bp2_computed, sv->sample_info.oar_bp2_reads, sv->sample_info.oar_bp2_consistent_reads, sv->sample_info.oar_bp2_consistent_hq_reads, sv->sample_info.oar_bp2_exact_reads);
    update_oar_max_fields(out_hdr, sv->vcf_entry, 2, bp2_computed,
        sv->sample_info.oar_bp2_reads_by_hpid, sv->sample_info.oar_bp2_reads_by_vid);
    update_other_read_support_fields(out_hdr, sv->vcf_entry, "ORR", 2, bp2_computed, sv->sample_info.orr_bp2_reads, sv->sample_info.orr_bp2_consistent_reads, sv->sample_info.orr_bp2_consistent_hq_reads, sv->sample_info.orr_bp2_exact_reads);

    int td = sv->sample_info.too_deep;
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "TD", &td, 1);

    int median_depths[] = {sv->sample_info.left_flanking_cov, sv->sample_info.indel_left_cov, sv->sample_info.indel_right_cov, sv->sample_info.right_flanking_cov};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MD", median_depths, 4);

    int median_depths_highmq[] = {sv->sample_info.left_flanking_cov_highmq, sv->sample_info.indel_left_cov_highmq, sv->sample_info.indel_right_cov_highmq, sv->sample_info.right_flanking_cov_highmq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "MDHQ", median_depths_highmq, 4);

    int cluster_depths[] = {sv->sample_info.left_anchor_cov, sv->sample_info.right_anchor_cov};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "CLMD", cluster_depths, 2);

    int cluster_depths_highmq[] = {sv->sample_info.left_anchor_cov_highmq, sv->sample_info.right_anchor_cov_highmq};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "CLMDHQ", cluster_depths_highmq, 2);

    if (sv->min_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MINSIZE", &(sv->min_conf_size), 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MINSIZE", NULL, 0);
    }
    if (sv->max_conf_size != deletion_t::SIZE_NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MAXSIZE", &(sv->max_conf_size), 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "MAXSIZE", NULL, 0);
    }
    if (sv->ks_pval != deletion_t::KS_PVAL_NOT_COMPUTED) {
        float ks_pval = sv->ks_pval;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "KSPVAL", &ks_pval, 1);
    } else {
        bcf_update_format_float(out_hdr, sv->vcf_entry, "KSPVAL", NULL, 0);
    }

    int ext_reads[] = {sv->sample_info.alt_lext_reads, sv->sample_info.alt_rext_reads};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AXR", ext_reads, 2);
    int hq_ext_reads[] = {sv->sample_info.hq_alt_lext_reads, sv->sample_info.hq_alt_rext_reads};
    bcf_update_format_int32(out_hdr, sv->vcf_entry, "AXRHQ", hq_ext_reads, 2);

    if (sv->sample_info.alt_consensus1_metrics.length > 0) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AL", &(sv->sample_info.alt_consensus1_metrics.length), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AAS", &(sv->sample_info.alt_consensus1_metrics.alt_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARS", &(sv->sample_info.alt_consensus1_metrics.ref_score), 1);
        int ass[] = {sv->sample_info.alt_consensus1_metrics.split_sizes[0], sv->sample_info.alt_consensus1_metrics.split_sizes[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASS", ass, 2);
        int assc[] = {sv->sample_info.alt_consensus1_metrics.split_scores[0], sv->sample_info.alt_consensus1_metrics.split_scores[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASSC", assc, 2);
        int asscia[] = {sv->sample_info.alt_consensus1_metrics.independent_ref_scores[0], sv->sample_info.alt_consensus1_metrics.independent_ref_scores[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASSCIA", asscia, 2);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AL", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AAS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASSC", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASSCIA", NULL, 0);
    }
    if (sv->sample_info.alt_consensus2_metrics.length > 0) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AL2", &(sv->sample_info.alt_consensus2_metrics.length), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AAS2", &(sv->sample_info.alt_consensus2_metrics.alt_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARS2", &(sv->sample_info.alt_consensus2_metrics.ref_score), 1);
        int ass2[] = {sv->sample_info.alt_consensus2_metrics.split_sizes[0], sv->sample_info.alt_consensus2_metrics.split_sizes[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASS2", ass2, 2);
        int assc2[] = {sv->sample_info.alt_consensus2_metrics.split_scores[0], sv->sample_info.alt_consensus2_metrics.split_scores[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASSC2", assc2, 2);
        int assc2ia[] = {sv->sample_info.alt_consensus2_metrics.independent_ref_scores[0], sv->sample_info.alt_consensus2_metrics.independent_ref_scores[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASSC2IA", assc2ia, 2);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AL2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AAS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ARS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASSC2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "ASSC2IA", NULL, 0);
    }

    if (sv->sample_info.ext_alt_consensus1_metrics.length > 0) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XAAS", &(sv->sample_info.ext_alt_consensus1_metrics.alt_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XARS", &(sv->sample_info.ext_alt_consensus1_metrics.ref_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XAL", &(sv->sample_info.ext_alt_consensus1_metrics.length), 1);
        int xass[] = {sv->sample_info.ext_alt_consensus1_metrics.split_sizes[0], sv->sample_info.ext_alt_consensus1_metrics.split_sizes[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASS", xass, 2);
        int xassc[] = {sv->sample_info.ext_alt_consensus1_metrics.split_scores[0], sv->sample_info.ext_alt_consensus1_metrics.split_scores[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASSC", xassc, 2);
        int xasscia[] = {sv->sample_info.ext_alt_consensus1_metrics.independent_ref_scores[0], sv->sample_info.ext_alt_consensus1_metrics.independent_ref_scores[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASSCIA", xasscia, 2);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XAAS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XARS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XAL", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASS", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASSC", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASSCIA", NULL, 0);
    }

    if (sv->sample_info.ext_alt_consensus2_metrics.length > 0) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XAAS2", &(sv->sample_info.ext_alt_consensus2_metrics.alt_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XARS2", &(sv->sample_info.ext_alt_consensus2_metrics.ref_score), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XAL2", &(sv->sample_info.ext_alt_consensus2_metrics.length), 1);
        int xass2[] = {sv->sample_info.ext_alt_consensus2_metrics.split_sizes[0], sv->sample_info.ext_alt_consensus2_metrics.split_sizes[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASS2", xass2, 2);
        int xassc2[] = {sv->sample_info.ext_alt_consensus2_metrics.split_scores[0], sv->sample_info.ext_alt_consensus2_metrics.split_scores[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASSC2", xassc2, 2);
        int xassc2ia[] = {sv->sample_info.ext_alt_consensus2_metrics.independent_ref_scores[0], sv->sample_info.ext_alt_consensus2_metrics.independent_ref_scores[1]};
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASSC2IA", xassc2ia, 2);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XAAS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XARS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XAL2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASS2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASSC2", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "XASSC2IA", NULL, 0);
    }

    if (sv->sample_info.alt1_occ_ratio != sv_t::sample_info_t::NOT_COMPUTED) {
        float occ_ratio1 = sv->sample_info.alt1_occ_ratio;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1C_OCCR", &occ_ratio1, 1);
    } else {
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1C_OCCR", NULL, 0);
    }
    if (sv->sample_info.alt2_occ_ratio != sv_t::sample_info_t::NOT_COMPUTED) {
        float occ_ratio2 = sv->sample_info.alt2_occ_ratio;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR2C_OCCR", &occ_ratio2, 1);
    } else {
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR2C_OCCR", NULL, 0);
    }

    if (sv->sample_info.alt1_hp_len_mode != sv_t::sample_info_t::NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1HPMODE", &(sv->sample_info.alt1_hp_len_mode), 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1HPMODE", NULL, 0);
    }
    if (sv->sample_info.alt1_consistent_hp_len_mode != sv_t::sample_info_t::NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1CHPMODE", &(sv->sample_info.alt1_consistent_hp_len_mode), 1);
        float ar1_hp_iqr = sv->sample_info.alt1_consistent_hp_len_iqr;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1CHPIQR", &ar1_hp_iqr, 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1CHPMODE", NULL, 0);
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1CHPIQR", NULL, 0);
    }
    if (sv->sample_info.alt1_hp_min_mapq != sv_t::sample_info_t::NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1CHPmQ", &(sv->sample_info.alt1_hp_min_mapq), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1CHPMQ", &(sv->sample_info.alt1_hp_max_mapq), 1);
        float ar1_hp_avg_mapq = sv->sample_info.alt1_hp_avg_mapq;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1CHPAQ", &ar1_hp_avg_mapq, 1);
        float ar1_hp_stddev_mapq = sv->sample_info.alt1_hp_stddev_mapq;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1CHPSQ", &ar1_hp_stddev_mapq, 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1CHPmQ", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "AR1CHPMQ", NULL, 0);
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1CHPAQ", NULL, 0);
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1CHPSQ", NULL, 0);
    }
    if (sv->sample_info.alt1_hp_5p_mismatch_rate != sv_t::sample_info_t::NOT_COMPUTED) {
        float ar1_hp_5p_mismatch_rate = sv->sample_info.alt1_hp_5p_mismatch_rate;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1HP5PMR", &ar1_hp_5p_mismatch_rate, 1);
        float ar1_hp_3p_mismatch_rate = sv->sample_info.alt1_hp_3p_mismatch_rate;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1HP3PMR", &ar1_hp_3p_mismatch_rate, 1);
    } else {
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1HP5PMR", NULL, 0);
        bcf_update_format_float(out_hdr, sv->vcf_entry, "AR1HP3PMR", NULL, 0);
    }
    if (sv->sample_info.ref1_hp_len_mode != sv_t::sample_info_t::NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1HPMODE", &(sv->sample_info.ref1_hp_len_mode), 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1HPMODE", NULL, 0);
    }
    if (sv->sample_info.ref1_consistent_hp_len_mode != sv_t::sample_info_t::NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1CHPMODE", &(sv->sample_info.ref1_consistent_hp_len_mode), 1);
        float rr1_hp_iqr = sv->sample_info.ref1_consistent_hp_len_iqr;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1CHPIQR", &rr1_hp_iqr, 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1CHPMODE", NULL, 0);
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1CHPIQR", NULL, 0);
    }
    if (sv->sample_info.ref1_hp_min_mapq != sv_t::sample_info_t::NOT_COMPUTED) {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1CHPmQ", &(sv->sample_info.ref1_hp_min_mapq), 1);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1CHPMQ", &(sv->sample_info.ref1_hp_max_mapq), 1);
        float rr1_hp_avg_mapq = sv->sample_info.ref1_hp_avg_mapq;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1CHPAQ", &rr1_hp_avg_mapq, 1);
        float rr1_hp_stddev_mapq = sv->sample_info.ref1_hp_stddev_mapq;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1CHPSQ", &rr1_hp_stddev_mapq, 1);
    } else {
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1CHPmQ", NULL, 0);
        bcf_update_format_int32(out_hdr, sv->vcf_entry, "RR1CHPMQ", NULL, 0);
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1CHPAQ", NULL, 0);
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1CHPSQ", NULL, 0);
    }
    if (sv->sample_info.ref1_hp_5p_mismatch_rate != sv_t::sample_info_t::NOT_COMPUTED) {
        float rr1_hp_5p_mismatch_rate = sv->sample_info.ref1_hp_5p_mismatch_rate;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1HP5PMR", &rr1_hp_5p_mismatch_rate, 1);
        float rr1_hp_3p_mismatch_rate = sv->sample_info.ref1_hp_3p_mismatch_rate;
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1HP3PMR", &rr1_hp_3p_mismatch_rate, 1);
    } else {
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1HP5PMR", NULL, 0);
        bcf_update_format_float(out_hdr, sv->vcf_entry, "RR1HP3PMR", NULL, 0);
    }

    std::string filters;
    for (size_t i = 0; i < sv->sample_info.filters.size(); ++i) {
        filters += sv->sample_info.filters[i];
        if (i + 1 < sv->sample_info.filters.size()) {
            filters += ";";
        }
    }
    if (filters.empty()) {
        filters = "PASS";
    }

    const char* filters_cstr = filters.c_str();
    bcf_update_format_char(out_hdr, sv->vcf_entry, "FT", filters_cstr, strlen(filters_cstr));
}

hts_pos_t get_covered_bps(std::vector<hts_pair_pos_t>& pos_pairs) {
    if (pos_pairs.empty()) {
        return 0;
    }
    
    std::sort(pos_pairs.begin(), pos_pairs.end(), [](const hts_pair_pos_t &a, const hts_pair_pos_t &b) {
        return a.beg < b.beg;
    });
    
    hts_pos_t total = 0;
    hts_pos_t current_start = pos_pairs[0].beg;
    hts_pos_t current_end = pos_pairs[0].end;
    
    for (size_t i = 1; i < pos_pairs.size(); ++i) {
        if (pos_pairs[i].beg <= current_end) {
            current_end = std::max(current_end, pos_pairs[i].end);
        } else {
            total += current_end - current_start;
            current_start = pos_pairs[i].beg;
            current_end = pos_pairs[i].end;
        }
    }
    total += current_end - current_start;
    
    return total;
}

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, std::vector<bp_support_read_t>& reads,
    std::vector<bool>& is_consistent_read, std::vector<bool>& is_exact_read,
    double consistent_avg_score, double consistent_stddev_score) {

    if (reads.size() != is_consistent_read.size() || reads.size() != is_exact_read.size()) {
        throw std::runtime_error("Read classification vector size mismatch.");
    }

    bp_reads_info.computed = true;
    bp_reads_info.reads = reads.size();

    std::vector<int> mqs;
    std::vector<hts_pair_pos_t> fwd_mate_positions, rev_mate_positions;
    std::vector<hts_pair_pos_t> fwd_hq_mate_positions, rev_hq_mate_positions;

    double sum_mq = 0;
    int n_consistent_reads = 0;
    for (size_t i = 0; i < reads.size(); ++i) {
        const bp_support_read_t& read = reads[i];
        int mq = read.mate_mapq;
        bool is_consistent = is_consistent_read[i];
        bool is_exact = is_exact_read[i];

        if (read.mate_is_reverse) {
            if (is_consistent) bp_reads_info.consistent_fwd++;
            if (is_exact) bp_reads_info.exact_fwd++;
            if (is_consistent && is_exact) bp_reads_info.consistent_exact_fwd++;
            if (is_consistent) {
                rev_mate_positions.push_back({read.mate_pos, read.mate_endpos});
                if (mq >= config.high_confidence_mapq) {
                    rev_hq_mate_positions.push_back({read.mate_pos, read.mate_endpos});
                }
            }
        } else {
            if (is_consistent) bp_reads_info.consistent_rev++;
            if (is_exact) bp_reads_info.exact_rev++;
            if (is_consistent && is_exact) bp_reads_info.consistent_exact_rev++;
            if (is_consistent) {
                fwd_mate_positions.push_back({read.mate_pos, read.mate_endpos});
                if (mq >= config.high_confidence_mapq) {
                    fwd_hq_mate_positions.push_back({read.mate_pos, read.mate_endpos});
                }
            }
        }

        if (is_consistent) {
            bp_reads_info.consistent_min_mq = std::min(bp_reads_info.consistent_min_mq, mq);
            bp_reads_info.consistent_max_mq = std::max(bp_reads_info.consistent_max_mq, mq);
            if (mq >= config.high_confidence_mapq) {
                bp_reads_info.consistent_high_mq++;
            }
            sum_mq += mq;
            mqs.push_back(mq);
            n_consistent_reads++;
        }
        if (is_consistent && is_exact) {
            bp_reads_info.consistent_exact_min_mq = std::min(bp_reads_info.consistent_exact_min_mq, mq);
            bp_reads_info.consistent_exact_max_mq = std::max(bp_reads_info.consistent_exact_max_mq, mq);
            if (mq >= config.high_confidence_mapq) {
                bp_reads_info.consistent_exact_high_mq++;
            }
        }
        if (is_exact) {
            bp_reads_info.exact_min_mq = std::min(bp_reads_info.exact_min_mq, mq);
            bp_reads_info.exact_max_mq = std::max(bp_reads_info.exact_max_mq, mq);
            if (mq >= config.high_confidence_mapq) {
                bp_reads_info.exact_high_mq++;
            }
        }
    }
    bp_reads_info.consistent_avg_score = consistent_avg_score;
    bp_reads_info.consistent_stddev_score = consistent_stddev_score;

    bp_reads_info.consistent_avg_mq = sum_mq/std::max(1, n_consistent_reads);
    bp_reads_info.consistent_stddev_mq = stddev(mqs);

    bp_reads_info.fwd_mate_cov_bps = get_covered_bps(fwd_mate_positions);
    bp_reads_info.rev_mate_cov_bps = get_covered_bps(rev_mate_positions);
    bp_reads_info.fwd_hq_mate_cov_bps = get_covered_bps(fwd_hq_mate_positions);
    bp_reads_info.rev_hq_mate_cov_bps = get_covered_bps(rev_hq_mate_positions);
}

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, std::vector<std::shared_ptr<bam1_t>>& reads,
    std::vector<bool>& is_consistent_read, std::vector<bool>& is_exact_read,
    double consistent_avg_score, double consistent_stddev_score) {

    std::vector<bp_support_read_t> bp_reads;
    bp_reads.reserve(reads.size());
    for (const std::shared_ptr<bam1_t>& read : reads) {
        bp_reads.emplace_back(read.get());
    }
    set_bp_consensus_info(bp_reads_info, bp_reads, is_consistent_read, is_exact_read,
        consistent_avg_score, consistent_stddev_score);
}

std::vector<std::string> gen_consensus_seqs(std::string ref_seq, std::vector<std::string>& seqs, const std::vector<const uint8_t*>& quals) {
    std::vector<std::string> temp1, temp2;
    std::vector<StripedSmithWaterman::Alignment> consensus_contigs_alns;

    std::vector<std::string> consensus_seqs; 

    consensus_seqs = generate_reference_guided_consensus(ref_seq, temp1, seqs, temp2, aligner, harsh_aligner, consensus_contigs_alns, config, stats, false);
    consensus_seqs.push_back("");

    std::vector<seq_w_pp_t> seqs_w_pp, temp3, temp4;
    for (std::string& seq : seqs) {
        seqs_w_pp.push_back({seq, true, true});
    }
    std::vector<std::string> consensus_seqs2 = assemble_reads(temp3, seqs_w_pp, temp4, config, stats);
    consensus_seqs.insert(consensus_seqs.end(), consensus_seqs2.begin(), consensus_seqs2.end());

    for (std::string& consensus_seq : consensus_seqs) {
        if (!consensus_seq.empty() && consensus_seq != "HAS_CYCLE") correct_contig(consensus_seq, seqs, config.max_seq_error, config.min_clip_len, quals);
    }

    return consensus_seqs;
}

// Returns a consistency mask over reads; is_exact_match uses the same index space.
std::vector<bool> gen_consensus_and_classify_seqs(std::string ref_seq,
    std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<bool> revcomp_read, std::string& consensus_seq, double& avg_score, double& stddev_score, 
    std::vector<bool>& is_exact_match) {

    if (reads.empty()) {
        avg_score = 0;
        stddev_score = 0;
        consensus_seq = "";
        is_exact_match.assign(reads.size(), false);
        return {};
    }

    if (revcomp_read.empty()) {
        revcomp_read.resize(reads.size(), false);
    }

    std::vector<std::string> seqs;
    std::vector<std::vector<uint8_t>> quals_storage;
    for (int i = 0; i < reads.size(); i++) {
        std::shared_ptr<bam1_t> read = reads[i];
        std::string seq = get_sequence(read.get());
        const uint8_t* read_quals = bam_get_qual(read.get());
        quals_storage.emplace_back(read_quals, read_quals + read->core.l_qseq);
        if (revcomp_read[i]) {
            rc(seq);
            std::reverse(quals_storage.back().begin(), quals_storage.back().end());
        }
        seqs.push_back(seq);
    }

    std::vector<const uint8_t*> quals;
    for (const std::vector<uint8_t>& read_quals : quals_storage) {
        quals.push_back(read_quals.data());
    }

    avg_score = 0;
    std::vector<std::string> consensus_seqs = gen_consensus_seqs(ref_seq, seqs, quals);

    std::vector<std::shared_ptr<bam1_t>> consistent_reads;
    std::vector<int> start_positions, end_positions;
    std::vector<int> chosen_seqs_idxs;
    double cum_score = 0;
    std::vector<double> aln_scores;
    int separation = 0;
    int chosen_cseq_idx = -1;
    for (int i = 0; i < consensus_seqs.size(); i++) {
        std::string cseq = consensus_seqs[i];
        if (cseq.empty()) {
            separation = i;
            continue;
        }
        if (cseq == "HAS_CYCLE") continue;

        std::vector<std::shared_ptr<bam1_t>> curr_consistent_reads;
        std::vector<int> curr_start_positions, curr_end_positions;
        std::vector<int> curr_seqs_idxs;
        double curr_cum_score = 0;
        std::vector<double> curr_aln_scores;
        for (int j = 0; j < reads.size(); j++) {
            std::shared_ptr<bam1_t> read = reads[j];
            const std::string& read_seq = seqs[j];

            ungapped_aln_t ungapped_aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(), cseq.c_str(), cseq.length(), std::max(0, config.min_clip_len - 1));

            if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
                curr_seqs_idxs.push_back(j);
                curr_consistent_reads.push_back(read);
                curr_cum_score += double(ungapped_aln.score)/read_seq.length();
                curr_aln_scores.push_back(double(ungapped_aln.score)/read_seq.length());
                curr_start_positions.push_back(ungapped_aln.ref_begin);
                curr_end_positions.push_back(ungapped_aln.ref_end);
            }
        }

        curr_cum_score /= log(cseq.length());
        if (curr_cum_score > cum_score) {
            consistent_reads = curr_consistent_reads;
            chosen_seqs_idxs = curr_seqs_idxs;
            cum_score = curr_cum_score;
            aln_scores = curr_aln_scores;
            start_positions = curr_start_positions;
            end_positions = curr_end_positions;
            consensus_seq = cseq;
            chosen_cseq_idx = i;
        }
    }

    std::sort(start_positions.begin(), start_positions.end());
    std::sort(end_positions.begin(), end_positions.end(), std::greater<int>());
    if (!consistent_reads.empty()) avg_score = cum_score/consistent_reads.size();
    else avg_score = 0;
    stddev_score = stddev(aln_scores);

    std::vector<std::string> chosen_seqs;
    for (int idx : chosen_seqs_idxs) {
        chosen_seqs.push_back(seqs[idx]);
    }
    std::string evidence_consensus_seq;
    if (start_positions.size() >= 2) {
        correct_contig(consensus_seq, chosen_seqs, config.max_seq_error, config.min_clip_len);
        evidence_consensus_seq = consensus_seq;

        consensus_seq = consensus_seq.substr(start_positions[1], end_positions[1]-start_positions[1]);

        int p;
        int window_len = 2*stats.read_len/3;
        while ((p = find_char_in_str(consensus_seq, 'N', 0, window_len)) != -1) {
            consensus_seq = consensus_seq.substr(p+1);
        }

        while ((p = find_char_in_str(consensus_seq, 'N', consensus_seq.length()-window_len, consensus_seq.length())) != -1) {
            consensus_seq = consensus_seq.substr(0, p);
        }
    } else {
        consensus_seq = "";
    }

    // Correction and trimming can change which reads are consistent with the consensus, as well as which of those reads match it exactly.  
    // Recompute the returned evidence against the corrected consensus sequence.
    std::vector<bool> is_consistent_read(reads.size(), false);
    is_exact_match.assign(reads.size(), false);
    aln_scores.clear();
    cum_score = 0;
    int n_consistent_reads = 0;
    if (!evidence_consensus_seq.empty()) {
        for (int i = 0; i < reads.size(); i++) {
            const std::string& read_seq = seqs[i];
            ungapped_aln_t ungapped_aln = best_ungapped_aln(read_seq.c_str(), read_seq.length(),
                evidence_consensus_seq.c_str(), evidence_consensus_seq.length(), std::max(0, config.min_clip_len - 1));

            if (ungapped_aln.mismatch_rate() <= config.max_seq_error) {
                is_consistent_read[i] = true;
                is_exact_match[i] = ungapped_aln.mismatches == 0;
                n_consistent_reads++;
                double aln_score = double(ungapped_aln.score) / read_seq.length();
                cum_score += aln_score;
                aln_scores.push_back(aln_score);
            }
        }
    }

    if (n_consistent_reads > 0) {
        avg_score = cum_score / log(evidence_consensus_seq.length()) / n_consistent_reads;
    } else {
        avg_score = 0;
    }
    stddev_score = stddev(aln_scores);

    return is_consistent_read;
}

std::vector<bool> classify_seqs_with_ref_seq(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads,
    const std::vector<bool>& is_eligible_read, double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read) {

    if (reads.size() != is_eligible_read.size()) {
        throw std::runtime_error("Read eligibility vector size mismatch.");
    }

    if (reads.empty()) {
        avg_score = 0;
        stddev_score = 0;
        is_exact_read.assign(reads.size(), false);
        return {};
    }

    std::vector<bool> is_consistent_read(reads.size(), false);
    is_exact_read.assign(reads.size(), false);

    std::vector<double> aln_scores;
    for (int i = 0; i < reads.size(); i++) {
        if (!is_eligible_read[i]) continue;
        std::shared_ptr<bam1_t> read = reads[i];
        std::string seq = get_sequence(read.get(), true);
        if (!bam_is_mrev(read)) rc(seq);

        ungapped_aln_t ungapped_aln = best_ungapped_aln(seq.c_str(), seq.length(), ref_seq.c_str(), ref_seq.length(), config.min_clip_len - 1);
        if (ungapped_aln.query_end - ungapped_aln.query_begin <= 0) continue;

        double mismatch_rate = double(ungapped_aln.mismatches)/(ungapped_aln.query_end-ungapped_aln.query_begin);
        if (mismatch_rate <= config.max_seq_error) {
            is_consistent_read[i] = true;
            is_exact_read[i] = mismatch_rate == 0;
            aln_scores.push_back(double(ungapped_aln.score)/seq.length());
        }
    }

    avg_score = mean(aln_scores);
    stddev_score = stddev(aln_scores);

    return is_consistent_read;
}

std::pair<std::unordered_map<std::string, std::pair<std::string, int>>*, evidence_map_t*> acquire_chromosome_data(int contig_id) {
    mutex_per_chr[contig_id].lock();
    if (active_threads_per_chr[contig_id] == 0) {
		std::string fname = workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt";
		std::ifstream fin(fname);
		std::string qname, read_seq, qual; int mapq;
		while (fin >> qname >> read_seq >> qual >> mapq) {
			mateseqs_w_mapq[contig_id][qname] = {read_seq, mapq};
		}
    }
	if (load_evidence_maps && !evidence_map_load_attempted[contig_id]) {
        evidence_map_load_attempted[contig_id] = true;
        std::string contig_name = contig_map.get_name(contig_id);
        std::string alt_fname = reads_association_dir + "/" + contig_name + ".alt.txt";
        std::string ref_fname = reads_association_dir + "/" + contig_name + ".ref.txt";
        std::string er_fname = reads_association_dir + "/" + contig_name + ".er.txt";
        std::string td_fname = reads_association_dir + "/" + contig_name + ".td.txt";
        std::ifstream alt_fin(alt_fname), ref_fin(ref_fname), er_fin(er_fname), td_fin(td_fname);
        bool has_alt = bool(alt_fin), has_ref = bool(ref_fin), has_er = bool(er_fin);
        if (td_fin) {
            std::unordered_map<std::string, sv_t*> sv_by_id;
            for (sv_t* sv : evidence_svs_by_contig[contig_id]) sv_by_id[sv->id] = sv;
            std::string sv_id;
            while (td_fin >> sv_id) sv_by_id.at(sv_id)->sample_info.too_deep = true;
        }
        if (has_alt || has_ref || has_er) {
            evidence_maps_by_contig[contig_id].reset(new evidence_map_t());
            evidence_maps_by_contig[contig_id]->load(has_alt ? alt_fname : "", has_ref ? ref_fname : "", has_er ? er_fname : "", evidence_svs_by_contig[contig_id], config, reassigns_evidence(evidence_mode));
        }
	}
	active_threads_per_chr[contig_id]++;
	evidence_map_t* evidence_map = evidence_maps_by_contig[contig_id] ? evidence_maps_by_contig[contig_id].get() : &empty_evidence_map;
	mutex_per_chr[contig_id].unlock();
	return {&mateseqs_w_mapq[contig_id], evidence_map};
}

void release_chromosome_data(int contig_id) {
    mutex_per_chr[contig_id].lock();
	active_threads_per_chr[contig_id]--;
	if (active_threads_per_chr[contig_id] == 0) {
		mateseqs_w_mapq[contig_id].clear();
		evidence_maps_by_contig[contig_id].reset();
		evidence_map_load_attempted[contig_id] = false;
	}
	mutex_per_chr[contig_id].unlock();
}

IntervalTree<ext_read_t*> get_candidate_reads_for_extension_itree(std::string contig_name, hts_pos_t contig_len, std::vector<hts_pair_pos_t> target_ivals, open_samFile_t* bam_file,
                                                                  std::vector<ext_read_t*>& candidate_reads_for_extension) {
    int contig_id = contig_map.get_id(contig_name);
    candidate_reads_for_extension = get_extension_reads(contig_name, target_ivals, contig_len, config, stats, bam_file);
    std::vector<Interval<ext_read_t*>> it_ivals;
    for (ext_read_t* ext_read : candidate_reads_for_extension) {
        Interval<ext_read_t*> it_ival(ext_read->start, ext_read->end, ext_read);
        it_ivals.push_back(it_ival);
    }
    return IntervalTree<ext_read_t*>(it_ivals);
}

void clear_invalid_stat_tests(bcf_hdr_t* hdr, std::vector<std::shared_ptr<deletion_t>>& dels) {
    std::vector<std::pair<std::shared_ptr<deletion_t>, float>> dels_w_epr;
    hts_pos_t chr_len = 0;
    for (const auto& del : dels) {
        dels_w_epr.push_back({del, get_sv_epr(hdr, del->vcf_entry)});
        chr_len = std::max(chr_len, del->end);
    }

    std::sort(dels_w_epr.begin(), dels_w_epr.end(), 
    [](const std::pair<std::shared_ptr<deletion_t>, float>& a, const std::pair<std::shared_ptr<deletion_t>, float>& b) {
        return a.second > b.second; // Sort by EPR in descending order
    });

    SegTree seg_tree(chr_len+1);
    for (const auto& del_epr : dels_w_epr) {
        const auto& del = del_epr.first;
        hts_pos_t midpoint = (del->start + del->end) / 2;
        if (seg_tree.any_ge(midpoint-stats.max_is, midpoint, 1)) {
            del->ks_pval = deletion_t::KS_PVAL_NOT_COMPUTED;
            del->min_conf_size = deletion_t::SIZE_NOT_COMPUTED;
            del->max_conf_size = deletion_t::SIZE_NOT_COMPUTED;
        } else {
            seg_tree.add(midpoint-stats.max_is, midpoint, 1);
        }
    }
}

void rebalance_depth_cov(hts_pos_t strong_start, hts_pos_t strong_end, hts_pos_t weak_start, hts_pos_t weak_end, 
    int per_base_cov_delta, int& cov_to_update) {
    int ov = overlap(strong_start, strong_end, weak_start, weak_end);
    if (ov <= 0 || weak_end-weak_start <= 0) {
        return;
    }

    uint64_t weak_total_cov = (weak_end - weak_start) * cov_to_update;
    weak_total_cov += per_base_cov_delta * ov;
    cov_to_update = weak_total_cov / (weak_end - weak_start);
}

void rebalance_covs(bcf_hdr_t* hdr, std::vector<std::shared_ptr<deletion_t>>& dels) {
    std::vector<std::pair<std::shared_ptr<deletion_t>, float>> dels_w_epr;
    for (const auto& del : dels) {
        dels_w_epr.push_back({del, get_sv_epr(hdr, del->vcf_entry)});
    }

    std::vector<Interval<std::pair<std::shared_ptr<deletion_t>, float>>> it_ivals;
    for (const auto& del_epr : dels_w_epr) {
        const auto& del = del_epr.first;
        float epr = del_epr.second;
        Interval<std::pair<std::shared_ptr<deletion_t>, float>> it_ival(del->start, del->end, {del, epr});
        it_ivals.push_back(it_ival);
    }
    IntervalTree<std::pair<std::shared_ptr<deletion_t>, float>> it_tree(it_ivals);

    for (auto& curr_del : dels_w_epr) {
        auto ov_dels = it_tree.findOverlapping(curr_del.first->start, curr_del.first->end);
        for (const auto& it_del : ov_dels) {
            auto& ov_del = it_del.value;
            int ov_del_alt_ac = count_alt_alleles(hdr, ov_del.first->vcf_entry);
            if (ov_del.second <= curr_del.second || ov_del.first == curr_del.first || ov_del_alt_ac == 0) continue;

            deletion_t* strong_del = ov_del.first.get();
            deletion_t* weak_del = curr_del.first.get();
            int strong_del_alt_ac = ov_del_alt_ac;

            int strong_del_avg_flanking_cov = (strong_del->sample_info.left_flanking_cov + strong_del->sample_info.right_flanking_cov)/2;
            int strong_del_avg_indel_cov = (strong_del->sample_info.indel_left_cov + strong_del->sample_info.indel_right_cov)/2;
            int avg_depth_delta = std::max(0, (strong_del_avg_flanking_cov-strong_del_avg_indel_cov)/2 * strong_del_alt_ac);

            // int strong_del_avg_flanking_cov_hq = (strong_del->sample_info.left_flanking_cov_highmq + strong_del->sample_info.right_flanking_cov_highmq)/2;
            // int strong_del_avg_indel_cov_hq = (strong_del->sample_info.indel_left_cov_highmq + strong_del->sample_info.indel_right_cov_highmq)/2;
            // int avg_depth_delta_hq = std::max(0, (strong_del_avg_flanking_cov_hq-strong_del_avg_indel_cov_hq)/2 * strong_del_alt_ac);

            if (weak_del->end - weak_del->start <= config.indel_tested_region_size) {
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start, weak_del->end, avg_depth_delta, weak_del->sample_info.indel_left_cov);
                weak_del->sample_info.indel_right_cov = weak_del->sample_info.indel_left_cov;
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start, weak_del->end, avg_depth_delta, weak_del->sample_info.indel_left_cov_highmq);
                weak_del->sample_info.indel_right_cov_highmq = weak_del->sample_info.indel_left_cov_highmq;
            } else {
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start, weak_del->start+config.indel_tested_region_size, avg_depth_delta, weak_del->sample_info.indel_left_cov);
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->end-config.indel_tested_region_size, weak_del->end, avg_depth_delta, weak_del->sample_info.indel_right_cov);
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start, weak_del->start+config.indel_tested_region_size, avg_depth_delta, weak_del->sample_info.indel_left_cov_highmq);
                rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->end-config.indel_tested_region_size, weak_del->end, avg_depth_delta, weak_del->sample_info.indel_right_cov_highmq);
                
            }
            rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->start-config.flanking_size, weak_del->start, avg_depth_delta, weak_del->sample_info.left_flanking_cov);
            rebalance_depth_cov(strong_del->start, strong_del->end, weak_del->end, weak_del->end+config.flanking_size, avg_depth_delta, weak_del->sample_info.right_flanking_cov);
        }
    }
}

int main(int argc, char* argv[]) {

    if (argc < 7) throw std::runtime_error("Usage: genotype <input.vcf.gz> <output.vcf.gz> <alignment.bam> <reference.fa> <workdir> <sample> [options]");
    std::string in_vcf_fname = argv[1];
    std::string out_vcf_fname = argv[2];
    bam_fname = argv[3];
    reference_fname = argv[4];
    workdir = argv[5];
    std::string sample_name = argv[6];

    std::string alt_reads_association_dir = workdir + "/reads_to_sv_associations";
    for (int i = 7; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--reassign-evidence") {
            evidence_mode = evidence_mode_t::REASSIGN;
        } else if (arg == "--cached-evidence") {
            evidence_mode = evidence_mode_t::CACHED;
        } else if (arg == "--alt-read-associations") {
            if (++i >= argc) {
                throw std::runtime_error("Missing path after --alt-read-associations.");
            }
            alt_reads_association_dir = argv[i];
        } else {
            throw std::runtime_error("Unknown genotype option: " + arg);
        }
    }
    bool reassign_evidence = reassigns_evidence(evidence_mode);
    load_evidence_maps = true;
    reads_association_dir = alt_reads_association_dir;

    contig_map.load(workdir);
    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    chr_seqs.read_fasta_into_map(reference_fname);
    bam_pool = new bam_pool_t(config.threads, bam_fname, reference_fname);

    // read crossing isize distribution
    std::ifstream crossing_isizes_dist_fin(workdir + "/crossing_isizes.txt");
	int isize, count;
	while (crossing_isizes_dist_fin >> isize >> count) {
		for (int i = 0; i < count; i++) global_crossing_isize_dist.push_back(isize);
	}
	std::random_shuffle(global_crossing_isize_dist.begin(), global_crossing_isize_dist.end());
	if (global_crossing_isize_dist.size() > 100000) global_crossing_isize_dist.resize(100000);
    std::sort(global_crossing_isize_dist.begin(), global_crossing_isize_dist.end());
    global_crossing_isize_dist.shrink_to_fit();
	crossing_isizes_dist_fin.close();

    std::string full_cmd_fname = workdir + "/cmd.txt";
	std::ifstream full_cmd_fin(full_cmd_fname);
    std::string full_cmd_str;
	std::getline(full_cmd_fin, full_cmd_str);

    mateseqs_w_mapq.resize(contig_map.size());
    active_threads_per_chr = std::vector<int>(contig_map.size());
	mutex_per_chr = std::vector<std::mutex>(contig_map.size());
    evidence_svs_by_contig.resize(contig_map.size());
    evidence_maps_by_contig.resize(contig_map.size());
    evidence_map_load_attempted = std::vector<bool>(contig_map.size());

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    if (in_vcf_file == NULL) {
        throw std::runtime_error("Unable to open file " + in_vcf_fname + ".");
    }

    bcf_hdr_t* in_vcf_header = bcf_hdr_read(in_vcf_file);
    if (in_vcf_header == NULL) {
        throw std::runtime_error("Failed to read the VCF header.");
    }

    bcf1_t* vcf_record = bcf_init();
    std::unordered_map<std::string, std::vector<std::shared_ptr<sv_t>>> hp_by_chr;
    std::unordered_map<std::string, std::vector<std::shared_ptr<deletion_t>>> dels_by_chr;
    std::unordered_map<std::string, std::vector<std::shared_ptr<duplication_t>>> dups_by_chr;
    std::unordered_map<std::string, std::vector<std::shared_ptr<insertion_t>>> inss_by_chr;
    std::unordered_map<std::string, std::vector<std::shared_ptr<inversion_t>>> invs_by_chr;
    while (bcf_read(in_vcf_file, in_vcf_header, vcf_record) == 0) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(in_vcf_header, vcf_record);
        if (sv == nullptr) {
            std::cout << "Ignoring SV of unsupported type: " << vcf_record->d.id << std::endl; 
            continue;
        }
        if (sv->start > sv->end) {
            std::cout << "Discarding SV with invalid coordinates: " << sv->id << std::endl;
            continue;
        }

        sv->vcf_entry = bcf_dup(vcf_record);
        bool retained = true;
        if (should_genotype_as_hp_indel(sv.get(), chr_seqs.get_seq(sv->chr), chr_seqs.get_len(sv->chr))) {
            hp_by_chr[sv->chr].push_back(sv);
        } else if (sv->svtype() == "DEL") {
            dels_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<deletion_t>(sv));
        } else if (sv->svtype() == "DUP") {
            dups_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<duplication_t>(sv));
        } else if (sv->svtype() == "INS") {
        	inss_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<insertion_t>(sv));
        } else if (sv->svtype() == "INV") {
            invs_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<inversion_t>(sv));
        } else {
            retained = false;
        }
        if (retained) evidence_svs_by_contig[contig_map.get_id(sv->chr)].push_back(sv.get());
    }

    int* imap = NULL;
    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    bcf_hdr_t* out_vcf_header = bcf_subset_header(in_vcf_header, sample_name, imap);
    add_fmt_tags(out_vcf_header);
    if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
	        throw std::runtime_error("Failed to read the VCF header.");
    }

    // genotype chrs in descending order of svs
    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    const int BLOCK_SIZE = 20;

    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
	    	std::string contig_name = contig_map.get_name(contig_id);

        std::vector<std::shared_ptr<sv_t>>& hps = hp_by_chr[contig_name];
        std::vector<hts_pair_pos_t> ref_hp_ranges;
        for (std::shared_ptr<sv_t>& hp : hps) {
            ref_hp_ranges.push_back(find_ref_hp_range_for_indel(hp.get(), chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name)));
        }

        std::vector<sv_t*> block_hps;
        for (int i = 0; i < hps.size(); i++) {
            block_hps.push_back(hps[i].get());
            if ((i == hps.size()-1 && !block_hps.empty()) 
            || (block_hps.size() >= BLOCK_SIZE && ref_hp_ranges[i].beg != ref_hp_ranges[i+1].beg)) {
                std::future<void> future = thread_pool.push(genotype_hp_indels, contig_name, chr_seqs.get_seq(contig_name),
                        chr_seqs.get_len(contig_name), block_hps, std::ref(stats), std::ref(config), std::ref(contig_map), bam_pool,
                        &global_crossing_isize_dist, evidence_mode);
                block_hps.clear();
            }
        }

        std::vector<std::shared_ptr<deletion_t>>& dels = dels_by_chr[contig_name];
        std::vector<deletion_t*> block_dels;
        for (int i = 0; i < dels.size(); i++) {
            if (!reassign_evidence || genotype_when_reassigning_evidence(dels[i].get())) {
                block_dels.push_back(dels[i].get());
            }
            if (block_dels.size() == BLOCK_SIZE || (i == dels.size()-1 && !block_dels.empty())) {
                std::future<void> future = thread_pool.push(genotype_dels, contig_name, chr_seqs.get_seq(contig_name),
                        chr_seqs.get_len(contig_name), block_dels, in_vcf_header, out_vcf_header, std::ref(stats), std::ref(config),
                        std::ref(contig_map), bam_pool, workdir, &global_crossing_isize_dist);
                futures.push_back(std::move(future));
                block_dels.clear();
            }
        }

        std::vector<std::shared_ptr<duplication_t>>& dups = dups_by_chr[contig_name];
        std::vector<duplication_t*> block_dups;
        for (int i = 0; i < dups.size(); i++) {
            if (!reassign_evidence || genotype_when_reassigning_evidence(dups[i].get())) {
                block_dups.push_back(dups[i].get());
            }
            if (block_dups.size() == BLOCK_SIZE || (i == dups.size()-1 && !block_dups.empty())) {
                std::future<void> future = thread_pool.push(genotype_dups, contig_name, chr_seqs.get_seq(contig_name),
                        chr_seqs.get_len(contig_name), block_dups, in_vcf_header, out_vcf_header, std::ref(stats), std::ref(config),
                        std::ref(contig_map), bam_pool, workdir, &global_crossing_isize_dist);
                futures.push_back(std::move(future));
                block_dups.clear();
            }
        }

        std::vector<std::shared_ptr<insertion_t>>& inss = inss_by_chr[contig_name];
        std::vector<insertion_t*> block_inss;
        for (int i = 0; i < inss.size(); i++) {
            if (!reassign_evidence || genotype_when_reassigning_evidence(inss[i].get())) {
                block_inss.push_back(inss[i].get());
            }
            if (block_inss.size() == BLOCK_SIZE || (i == inss.size()-1 && !block_inss.empty())) {
                std::future<void> future = thread_pool.push(genotype_inss, contig_name, chr_seqs.get_seq(contig_name),
                        chr_seqs.get_len(contig_name), block_inss, in_vcf_header, out_vcf_header, std::ref(stats), std::ref(config),
                        std::ref(contig_map), bam_pool, &global_crossing_isize_dist);
                futures.push_back(std::move(future));
                block_inss.clear();
            }
        }

        std::vector<std::shared_ptr<inversion_t>>& invs = invs_by_chr[contig_name];
        std::vector<inversion_t*> block_invs;
        for (int i = 0; i < invs.size(); i++) {
            if (!reassign_evidence || genotype_when_reassigning_evidence(invs[i].get())) {
                block_invs.push_back(invs[i].get());
            }
            if (block_invs.size() == BLOCK_SIZE || (i == invs.size()-1 && !block_invs.empty())) {
                std::future<void> future = thread_pool.push(genotype_invs, contig_name, chr_seqs.get_seq(contig_name),
                        chr_seqs.get_len(contig_name), block_invs, in_vcf_header, out_vcf_header, std::ref(stats), std::ref(config),
                        std::ref(contig_map), bam_pool);
                futures.push_back(std::move(future));
                block_invs.clear();
            }
        }
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();
    // if (reassign_evidence) {
    //     for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
    //         std::string contig_name = contig_map.get_name(contig_id);
    //         std::vector<std::shared_ptr<deletion_t>>& dels = dels_by_chr[contig_name];
    //         clear_invalid_stat_tests(in_vcf_header, dels);
    //         rebalance_covs(in_vcf_header, dels);
    //     }
    // }

    // print contigs in vcf order
    int n_seqs;
    const char** seqnames = bcf_hdr_seqnames(in_vcf_header, &n_seqs);
    for (int i = 0; i < n_seqs; i++) {
        	std::string contig_name = seqnames[i];
        	std::vector<std::shared_ptr<sv_t>> contig_svs;
            if (hp_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), hp_by_chr[contig_name].begin(), hp_by_chr[contig_name].end());
        	if (dels_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), dels_by_chr[contig_name].begin(), dels_by_chr[contig_name].end());
        	if (dups_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), dups_by_chr[contig_name].begin(), dups_by_chr[contig_name].end());
        	if (inss_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), inss_by_chr[contig_name].begin(), inss_by_chr[contig_name].end());
            if (invs_by_chr.count(contig_name) > 0) contig_svs.insert(contig_svs.end(), invs_by_chr[contig_name].begin(), invs_by_chr[contig_name].end());
        	std::sort(contig_svs.begin(), contig_svs.end(), [](const std::shared_ptr<sv_t>& sv1, const std::shared_ptr<sv_t>& sv2) {
                return std::tie(sv1->start, sv1->end, sv1->id) < std::tie(sv2->start, sv2->end, sv2->id);
            });

			for (auto& sv : contig_svs) {
				bcf_update_info_int32(out_vcf_header, vcf_record, "AC", NULL, 0);
				bcf_update_info_int32(out_vcf_header, vcf_record, "AN", NULL, 0);
                if (!reassign_evidence || genotype_when_reassigning_evidence(sv.get())) update_record(in_vcf_header, out_vcf_header, sv.get(), chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), imap[0]);
				if (bcf_write(out_vcf_file, out_vcf_header, sv->vcf_entry) != 0) throw std::runtime_error("Failed to write VCF record to " + out_vcf_fname);
			}
	}
    free(seqnames);
    delete[] imap;

    bcf_destroy(vcf_record);
    bcf_hdr_destroy(out_vcf_header);
    bcf_hdr_destroy(in_vcf_header);
    bcf_close(in_vcf_file);
    bcf_close(out_vcf_file);
    delete bam_pool;

    tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
