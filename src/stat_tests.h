#ifndef STAT_TESTS_H
#define STAT_TESTS_H

#include <cstdint>
#include <cstring>
#include <cmath>
#include <htslib/sam.h>
#include <unordered_map>

#include "../libs/ssw_cpp.h"
#include "../libs/ks-test.h"
#include "genotype.h"
#include "htslib/hts.h"
#include "types.h"
#include "utils.h"
#include "sam_utils.h"
#include "sw_utils.h"

struct read_w_cached_info_t {
    bam1_t* read;
    hts_pos_t start, end, isize;
    int64_t as;
    int aln_len;
    bool left_clipped, right_clipped, is_rev, is_mrev;
    int references = 0;

    read_w_cached_info_t(bam1_t* read) : read(bam_dup1(read)), as(get_AS_tag(read)), start(read->core.pos), end(bam_endpos(read)),
                                         aln_len(get_aligned_portion_len(read)), left_clipped(is_left_clipped(read, 0)),
                                         right_clipped(is_right_clipped(read, 0)), is_rev(bam_is_rev(read)), is_mrev(bam_is_mrev(read)),
                                         isize(read->core.isize) {}

    ~read_w_cached_info_t() { bam_destroy1(read); }
};

struct region_depth_t {
	hts_pos_t start, end;
	std::vector<uint32_t> depths, highmq_depths;

	region_depth_t(hts_pos_t start, hts_pos_t end) : start(start), end(end), depths(std::vector<uint32_t>(end-start+1)),
			highmq_depths(std::vector<uint32_t>(end-start+1)) {}
};
struct region_w_median_t {
	hts_pos_t start, end;
	int* median_target;
	int* median_target_highmq_only;

	region_w_median_t(hts_pos_t start, hts_pos_t end, int* median_target, int* median_target_highmq_only) :
		start(start), end(end), median_target(median_target), median_target_highmq_only(median_target_highmq_only) {}
};
void depth_filter_indel(std::string contig_name, std::vector<sv_t*>& svs, open_samFile_t* bam_file, config_t& config, stats_t& stats) {

	if (svs.empty()) return;
	std::sort(svs.begin(), svs.end(), [](const sv_t* s1, const sv_t* s2) {
		return s1->start < s2->start;
	});

    std::vector<char*> regions;
    for (sv_t* sv : svs) {
        std::stringstream ss;
        ss << contig_name << ":" << std::max(hts_pos_t(1), sv->start - config.flanking_size) << "-" << std::min(sv->start + config.indel_tested_region_size, sv->end);
        char* region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);

        ss.str(std::string());
        ss << contig_name << ":" << std::max(sv->end - config.indel_tested_region_size, sv->start) << "-" << sv->end + config.flanking_size;
        region = new char[1000];
        strcpy(region, ss.str().c_str());
        regions.push_back(region);
    }

    std::vector<region_w_median_t> regions_of_interest;
    for (sv_t* sv : svs) {
    	regions_of_interest.emplace_back(std::max(hts_pos_t(0), sv->start - config.flanking_size), sv->start, &(sv->sample_info.left_flanking_cov), &(sv->sample_info.left_flanking_cov_highmq));
    	regions_of_interest.emplace_back(sv->start, std::min(sv->start + config.indel_tested_region_size, sv->end), &(sv->sample_info.indel_left_cov), &(sv->sample_info.indel_left_cov_highmq));
    	regions_of_interest.emplace_back(std::max(sv->end - config.indel_tested_region_size, sv->start), sv->end, &(sv->sample_info.indel_right_cov), &(sv->sample_info.indel_right_cov_highmq));
    	regions_of_interest.emplace_back(sv->end, sv->end + config.flanking_size, &(sv->sample_info.right_flanking_cov), &(sv->sample_info.right_flanking_cov_highmq));
		regions_of_interest.emplace_back(sv->left_anchor_aln->start, sv->left_anchor_aln->end, &(sv->sample_info.left_anchor_cov), &(sv->sample_info.left_anchor_cov_highmq));
		regions_of_interest.emplace_back(sv->right_anchor_aln->start, sv->right_anchor_aln->end, &(sv->sample_info.right_anchor_cov), &(sv->sample_info.right_anchor_cov_highmq));
	}
    std::sort(regions_of_interest.begin(), regions_of_interest.end(), [](region_w_median_t& r1, region_w_median_t& r2) {return r1.start < r2.start;});

    std::vector<region_depth_t> merged_regions_of_interest;
    auto curr_region = regions_of_interest[0];
    for (const auto& region : regions_of_interest) {
    	if (region.start - curr_region.end <= stats.read_len) {
    		curr_region.end = std::max(region.end, curr_region.end);
    	} else {
    		merged_regions_of_interest.emplace_back(curr_region.start, curr_region.end);
    		curr_region = region;
        }
    }
    merged_regions_of_interest.emplace_back(curr_region.start, curr_region.end);

    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	int curr_pos = 0;
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        hts_pos_t rs = read->core.pos, re = bam_endpos(read);
        while (curr_pos < merged_regions_of_interest.size() && rs > merged_regions_of_interest[curr_pos].end) curr_pos++;
        if (curr_pos == merged_regions_of_interest.size()) break;

        auto& curr_region = merged_regions_of_interest[curr_pos];

        if (overlap(curr_region.start, curr_region.end, rs, re)) {
			int ov_start = std::max(0, int(rs-curr_region.start));
			int ov_end = std::min(int(curr_region.depths.size()), int(re-curr_region.start));
			for (int i = ov_start; i < ov_end; i++) {
				curr_region.depths[i]++;
			}
		}

        if (read->core.qual < config.high_confidence_mapq) continue;

        if (overlap(curr_region.start, curr_region.end, rs, re)) {
        	int ov_start = std::max(0, int(rs-curr_region.start));
        	int ov_end = std::min(int(curr_region.depths.size()), int(re-curr_region.start));
        	for (int i = ov_start; i < ov_end; i++) {
        		curr_region.highmq_depths[i]++;
        	}
        }
	}

	curr_pos = 0;
	for (auto& region : regions_of_interest) {
		while (curr_pos < merged_regions_of_interest.size() && region.start > merged_regions_of_interest[curr_pos].end) curr_pos++;

		auto& merged_region = merged_regions_of_interest[curr_pos];
		std::vector<uint32_t> depths(merged_region.depths.begin() + (region.start - merged_region.start),
		                             merged_region.depths.begin() + (region.end - merged_region.start));
		if (depths.empty()) continue;
		std::sort(depths.begin(), depths.end());
		*region.median_target = depths[depths.size()/2];

		std::vector<uint32_t> highmq_depths(merged_region.highmq_depths.begin() + (region.start - merged_region.start),
		                                    merged_region.highmq_depths.begin() + (region.end - merged_region.start));
		std::sort(highmq_depths.begin(), highmq_depths.end());
		*region.median_target_highmq_only = highmq_depths[highmq_depths.size()/2];
	}

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

void depth_filter_del(std::string contig_name, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	std::vector<sv_t*> testable_dels(deletions.begin(), deletions.end());
    depth_filter_indel(contig_name, testable_dels, bam_file, config, stats);
}
void depth_filter_dup(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	std::vector<sv_t*> testable_dups(duplications.begin(), duplications.end());
    depth_filter_indel(contig_name, testable_dups, bam_file, config, stats);
}
void depth_filter_ins(std::string contig_name, std::vector<insertion_t*>& insertions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	std::vector<sv_t*> testable_inss(insertions.begin(), insertions.end());
	depth_filter_indel(contig_name, testable_inss, bam_file, config, stats);
}
void depth_filter_inv(std::string contig_name, std::vector<inversion_t*>& invertions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	std::vector<sv_t*> testable_invs(invertions.begin(), invertions.end());
	depth_filter_indel(contig_name, testable_invs, bam_file, config, stats);
}

struct pairs_data_t {
	std::vector<int> pairs_pos_mqs, pairs_neg_mqs;
	std::vector<int> pairs_pos_nms, pairs_neg_nms;
	hts_pos_t bp_lf_start = INT32_MAX, bp_lf_end = 0;
	hts_pos_t bp_rf_start = INT32_MAX, bp_rf_end = 0;
};

void set_bp_pairs_info(sv_t::bp_pairs_info_t& bp_pairs_info, std::vector<int>& supp_pairs_pos_mqs, std::vector<int>& supp_pairs_neg_mqs,
					   std::vector<int>& supp_pairs_pos_nms, std::vector<int>& supp_pairs_neg_nms, config_t& config) {

	bp_pairs_info.computed = true;
	bp_pairs_info.pairs = supp_pairs_pos_mqs.size();
	if (supp_pairs_pos_mqs.empty()) return;

	bp_pairs_info.pos_high_mapq = std::count_if(supp_pairs_pos_mqs.begin(), supp_pairs_pos_mqs.end(), [&](int mq) { return mq >= config.high_confidence_mapq; });
	bp_pairs_info.neg_high_mapq = std::count_if(supp_pairs_neg_mqs.begin(), supp_pairs_neg_mqs.end(), [&](int mq) { return mq >= config.high_confidence_mapq; });
	bp_pairs_info.pos_avg_mq = mean(supp_pairs_pos_mqs);
	bp_pairs_info.pos_stddev_mq = stddev(supp_pairs_pos_mqs);
	bp_pairs_info.neg_avg_mq = mean(supp_pairs_neg_mqs);
	bp_pairs_info.neg_stddev_mq = stddev(supp_pairs_neg_mqs);
	bp_pairs_info.pos_avg_nm = mean(supp_pairs_pos_nms);
	bp_pairs_info.pos_stddev_nm = stddev(supp_pairs_pos_nms);
	bp_pairs_info.neg_avg_nm = mean(supp_pairs_neg_nms);
	bp_pairs_info.neg_stddev_nm = stddev(supp_pairs_neg_nms);
	bp_pairs_info.pos_min_mq = *std::min_element(supp_pairs_pos_mqs.begin(), supp_pairs_pos_mqs.end());
	bp_pairs_info.neg_min_mq = *std::min_element(supp_pairs_neg_mqs.begin(), supp_pairs_neg_mqs.end());
	bp_pairs_info.pos_max_mq = *std::max_element(supp_pairs_pos_mqs.begin(), supp_pairs_pos_mqs.end());
	bp_pairs_info.neg_max_mq = *std::max_element(supp_pairs_neg_mqs.begin(), supp_pairs_neg_mqs.end());
}
void set_bp_pairs_info(sv_t::bp_pairs_info_t& bp_pairs_info, pairs_data_t& pairs_data, config_t& config) {
	set_bp_pairs_info(bp_pairs_info, pairs_data.pairs_pos_mqs, pairs_data.pairs_neg_mqs, pairs_data.pairs_pos_nms, pairs_data.pairs_neg_nms, config);
	bp_pairs_info.lf_span = pairs_data.bp_lf_end - pairs_data.bp_lf_start;
	bp_pairs_info.rf_span = pairs_data.bp_rf_end - pairs_data.bp_rf_start;
}

void count_stray_pairs(std::string contig_name, std::vector<deletion_t*> deletions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {

	if (deletions.empty()) return;

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->start < d2->start;
	});

	std::vector<char*> l_cluster_regions, r_cluster_regions;
	for (deletion_t* deletion : deletions) {
		std::stringstream ss;
		ss << contig_name << ":" << std::max(hts_pos_t(0), deletion->start-stats.max_is) << "-" << deletion->start;
		char* region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		l_cluster_regions.push_back(region);

		ss.str("");
		ss << contig_name << ":" << deletion->end << "-" << deletion->end+stats.max_is;
		region = new char[ss.str().length()+1];
		strcpy(region, ss.str().c_str());
		r_cluster_regions.push_back(region);
	}

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->start < d2->start;
	});

	std::vector<std::vector<int> > pos_mqs(deletions.size()), neg_mqs(deletions.size());
	std::vector<std::vector<int> > pos_nms(deletions.size()), neg_nms(deletions.size());

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, l_cluster_regions.data(), l_cluster_regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read) || bam_is_rev(read)) continue;

		while (curr_pos < deletions.size() && deletions[curr_pos]->start < read->core.pos) curr_pos++;

		// if the pair is discordant and it overlaps left_anchor_aln, increase l_cluster_region_disc_pair
		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_outward(read)) {
			for (int i = curr_pos; i < deletions.size() && deletions[i]->start-stats.max_is <= bam_endpos(read); i++) {
				if (read->core.pos <= deletions[i]->start) {
					pos_mqs[i].push_back(read->core.qual);
					neg_mqs[i].push_back(get_mq(read));

					pos_nms[i].push_back(get_nm(read));
					neg_nms[i].push_back(0);
				}
				// TODO: try to remap the mate instead?
			}
		}
	}

	for (int i = 0; i < deletions.size(); i++) {
		set_bp_pairs_info(deletions[i]->sample_info.bp1_stray_pairs, pos_mqs[i], neg_mqs[i], pos_nms[i], neg_nms[i], config);
	}

	std::sort(deletions.begin(), deletions.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->end < d2->end;
	});

	pos_mqs = std::vector<std::vector<int> >(deletions.size());
	neg_mqs = std::vector<std::vector<int> >(deletions.size());
	pos_nms = std::vector<std::vector<int> >(deletions.size());
	neg_nms = std::vector<std::vector<int> >(deletions.size());

	curr_pos = 0;
	iter = sam_itr_regarray(bam_file->idx, bam_file->header, r_cluster_regions.data(), r_cluster_regions.size());
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read) || !bam_is_rev(read)) continue;

		while (curr_pos < deletions.size() && deletions[curr_pos]->end+stats.max_is < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_outward(read)) {
			for (int i = curr_pos; i < deletions.size() && deletions[i]->end <= bam_endpos(read); i++) {
				if (read->core.pos <= deletions[i]->end+stats.max_is) {
					pos_mqs[i].push_back(get_mq(read));
					neg_mqs[i].push_back(read->core.qual);

					pos_nms[i].push_back(0);
					neg_nms[i].push_back(get_nm(read));
				}
			}
		}
	}

	for (int i = 0; i < deletions.size(); i++) {
		set_bp_pairs_info(deletions[i]->sample_info.bp2_stray_pairs, pos_mqs[i], neg_mqs[i], pos_nms[i], neg_nms[i], config);
	}

	for (char* region : l_cluster_regions) {
		delete[] region;
	}
	for (char* region : r_cluster_regions) {
		delete[] region;
	}
}

void count_stray_pairs(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, config_t& config, stats_t& stats) {

	if (duplications.empty()) return;

	std::vector<char*> lc_cluster_regions, rc_cluster_regions;
	for (duplication_t* dup : duplications) {
		std::stringstream ss;
		ss << contig_name << ":" << dup->start << "-" << dup->start+stats.max_is;
		rc_cluster_regions.push_back(strdup(ss.str().c_str()));

		ss.str("");
		ss << contig_name << ":" << std::max(hts_pos_t(0), dup->end-stats.max_is) << "-" << dup->end;
		lc_cluster_regions.push_back(strdup(ss.str().c_str()));
	}

	std::sort(duplications.begin(), duplications.end(), [](const duplication_t* d1, const duplication_t* d2) {
		return d1->start < d2->start;
	});

	std::vector<std::vector<int> > pos_mqs(duplications.size()), neg_mqs(duplications.size());
	std::vector<std::vector<int> > pos_nms(duplications.size()), neg_nms(duplications.size());

	bam1_t* read = bam_init1();
	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, rc_cluster_regions.data(), rc_cluster_regions.size());
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read) || !bam_is_rev(read)) continue;

		while (curr_pos < duplications.size() && duplications[curr_pos]->start+stats.max_is < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_long(read, stats.max_is)) {
			for (int i = curr_pos; i < duplications.size() && duplications[i]->start <= bam_endpos(read); i++) {
				if (read->core.pos <= duplications[i]->start+stats.max_is) {
					pos_mqs[i].push_back(get_mq(read));
					neg_mqs[i].push_back(read->core.qual);

					pos_nms[i].push_back(0);
					neg_nms[i].push_back(get_nm(read));
				}
			}
		}
	}
	hts_itr_destroy(iter);

	for (int i = 0; i < duplications.size(); i++) {
		set_bp_pairs_info(duplications[i]->sample_info.bp1_stray_pairs, pos_mqs[i], neg_mqs[i], pos_nms[i], neg_nms[i], config);
	}
	
	for (char* region : rc_cluster_regions) {
		free(region);
	}

	std::sort(duplications.begin(), duplications.end(), [](const duplication_t* d1, const duplication_t* d2) {
		return d1->end < d2->end;
	});

	pos_mqs = std::vector<std::vector<int> >(duplications.size());
	neg_mqs = std::vector<std::vector<int> >(duplications.size());
	pos_nms = std::vector<std::vector<int> >(duplications.size());
	neg_nms = std::vector<std::vector<int> >(duplications.size());

	curr_pos = 0;
	iter = sam_itr_regarray(bam_file->idx, bam_file->header, lc_cluster_regions.data(), lc_cluster_regions.size());
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read) || bam_is_rev(read)) continue;

		while (curr_pos < duplications.size() && duplications[curr_pos]->end < read->core.pos) curr_pos++;

		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || is_long(read, stats.max_is)) {
			for (int i = curr_pos; i < duplications.size() && duplications[i]->end-stats.max_is <= bam_endpos(read); i++) {
				if (read->core.pos <= duplications[i]->end) {
					pos_mqs[i].push_back(read->core.qual);
					neg_mqs[i].push_back(get_mq(read));

					pos_nms[i].push_back(get_nm(read));
					neg_nms[i].push_back(0);
				}
			}
		}
	}
	hts_itr_destroy(iter);

	for (int i = 0; i < duplications.size(); i++) {
		set_bp_pairs_info(duplications[i]->sample_info.bp2_stray_pairs, pos_mqs[i], neg_mqs[i], pos_nms[i], neg_nms[i], config);
	}

	for (char* region : lc_cluster_regions) {
		free(region);
	}
}

void calculate_confidence_interval_size(std::string contig_name, const std::vector<double>& global_crossing_isize_dist,
										std::vector<sv_t*>& svs, open_samFile_t* bam_file, 
										config_t& config, stats_t& stats) {

	if (svs.empty()) return;

	std::sort(svs.begin(), svs.end(), [](const sv_t* d1, const sv_t* d2) {
		return (d1->start+d1->end)/2 < (d2->start+d2->end)/2;
	});

    std::vector<hts_pos_t> midpoints, sizes;
    std::vector<uint64_t> sums(svs.size());
	std::vector<uint64_t> sq_sums(svs.size());
    std::vector<uint32_t> ns(svs.size());
    std::vector<char*> regions;
    std::vector<std::pair<hts_pos_t, hts_pos_t> > regions_coos;
    for (sv_t* sv : svs) {
        hts_pos_t midpoint = (sv->start+sv->end)/2, size = sv->end-sv->start;
        midpoints.push_back(midpoint);
        sizes.push_back(size);

        regions_coos.push_back({std::max(hts_pos_t(1), sv->start-stats.max_is), sv->start});
        regions_coos.push_back({std::max(hts_pos_t(1), midpoint-stats.max_is), midpoint});
    }
	std::sort(regions_coos.begin(), regions_coos.end());
    for (auto coos : regions_coos) {
    	std::stringstream ss;
    	ss << contig_name << ":" << coos.first << "-" << coos.second;
    	char* region = new char[1000];
    	strcpy(region, ss.str().c_str());
		regions.push_back(region);
    }

	int curr_pos = 0;
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    std::vector<std::vector<double> > local_dists(svs.size());
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {

    	if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
        if (!is_samechr(read) || is_samestr(read) || bam_is_rev(read) || read->core.isize <= 0) continue;

        while (curr_pos < svs.size() && midpoints[curr_pos] < read->core.pos) curr_pos++;

        hts_pos_t start = read->core.pos + read->core.l_qseq/2;
        hts_pos_t end = read->core.pos + read->core.isize - read->core.l_qseq/2;
        for (int i = curr_pos; i < midpoints.size() && midpoints[i] < read->core.pos+read->core.isize; i++) {
            if (start <= midpoints[i] && midpoints[i] <= end && read->core.isize <= stats.max_is+sizes[i]) {
                sums[i] += read->core.isize;
                sq_sums[i] += read->core.isize*read->core.isize;
                ns[i]++;
                local_dists[i].push_back(read->core.isize);
            }
        }
    }

    for (int i = 0; i < svs.size(); i++) {
		sv_t* sv = svs[i];
        int n = ns[i];
		uint64_t sum = sums[i], sq_sum = sq_sums[i];
        if (n >= 4) {
            int avg_is = sum/n;
            int var_is = (sq_sum - sum*sum/n)/(n-1);
            int confidence_ival = 2.576 * sqrt(var_is/n);
			sv->min_conf_size = abs(avg_is - stats.pop_avg_crossing_is) - confidence_ival;
            sv->max_conf_size = abs(avg_is - stats.pop_avg_crossing_is) + confidence_ival;

            if (!global_crossing_isize_dist.empty()) {
				double ks_pval = ks_test(global_crossing_isize_dist, local_dists[i]);
				if (std::isfinite(ks_pval)) sv->ks_pval = ks_pval;
				else sv->ks_pval = sv_t::KS_PVAL_NOT_COMPUTED;
            }
        }
    }

    for (char* region : regions) {
        delete[] region;
    }
    hts_itr_destroy(iter);
    bam_destroy1(read);
}

struct pairs_count_t {
	hts_pos_t pos;
	sv_t::bp_pairs_info_t* pairs_info;
	std::vector<int> supp_pairs_pos_mqs, supp_pairs_neg_mqs;
	std::vector<int> supp_pairs_pos_nms, supp_pairs_neg_nms;

	pairs_count_t(hts_pos_t pos, sv_t::bp_pairs_info_t* pairs_info) : pos(pos), pairs_info(pairs_info) {}
};

void calculate_ptn_ratio(std::string contig_name, std::vector<sv_t*>& svs, open_samFile_t* bam_file, config_t& config, stats_t& stats) {

	if (svs.empty()) return;

	std::vector<pairs_count_t> bkp_with_conc_pairs_count;
	for (sv_t* sv : svs) {
		bkp_with_conc_pairs_count.push_back({sv->start, &(sv->sample_info.ref_bp1.pairs_info)});
		bkp_with_conc_pairs_count.push_back({sv->end, &(sv->sample_info.ref_bp2.pairs_info)});
	}

	std::sort(bkp_with_conc_pairs_count.begin(), bkp_with_conc_pairs_count.end(), [](const pairs_count_t& p1, const pairs_count_t& p2) {
		return p1.pos < p2.pos;
	});

	std::vector<char*> regions;
	for (auto& b : bkp_with_conc_pairs_count) {
		std::stringstream ss;
		ss << contig_name << ":" << std::max(hts_pos_t(1), b.pos-stats.max_is) << "-" << b.pos+stats.max_is;
		char* region = strdup(ss.str().c_str());
		regions.push_back(region);
	}

	int curr_pos = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		while (curr_pos < bkp_with_conc_pairs_count.size() && bkp_with_conc_pairs_count[curr_pos].pos < read->core.pos-stats.max_is) curr_pos++;

		if (is_unmapped(read) || is_mate_unmapped(read) || !is_primary(read)) continue;
		if (!is_samechr(read) || is_samestr(read) || is_outward(read)) continue;
		if (std::abs(read->core.isize) < stats.min_is || std::abs(read->core.isize) > stats.max_is) continue;

		hts_pos_t start, end;
		if (bam_is_rev(read)) {
			start = read->core.mpos + stats.read_len/2;
			end = bam_endpos(read) - stats.read_len/2;
		} else {
			start = read->core.pos + stats.read_len/2;
			end = read->core.pos + read->core.isize - stats.read_len/2;
		}
		for (int i = curr_pos; i < bkp_with_conc_pairs_count.size() && bkp_with_conc_pairs_count[i].pos <= end; i++) {
			if (start <= bkp_with_conc_pairs_count[i].pos && bkp_with_conc_pairs_count[i].pos <= end) {
				if (!bam_is_rev(read)) {
					bkp_with_conc_pairs_count[i].supp_pairs_pos_mqs.push_back(read->core.qual);
					bkp_with_conc_pairs_count[i].supp_pairs_pos_nms.push_back(get_nm(read));
				} else {
					bkp_with_conc_pairs_count[i].supp_pairs_neg_mqs.push_back(read->core.qual);
					bkp_with_conc_pairs_count[i].supp_pairs_neg_nms.push_back(get_nm(read));
				}
			}
		}
	}

	for (char* region : regions) {
		free(region);
	}

	for (auto& b : bkp_with_conc_pairs_count) {
		set_bp_pairs_info(*(b.pairs_info), b.supp_pairs_pos_mqs, b.supp_pairs_neg_mqs, b.supp_pairs_pos_nms, b.supp_pairs_neg_nms, config);
	}
}

void add_read(pairs_data_t& pairs_data, bam1_t* read, int64_t mate_nm) {
	if (bam_is_rev(read)) {
		if (read->core.mpos < pairs_data.bp_lf_start) pairs_data.bp_lf_start = read->core.mpos;
		if (get_mate_endpos(read) > pairs_data.bp_lf_end) pairs_data.bp_lf_end = get_mate_endpos(read);
		if (read->core.pos < pairs_data.bp_rf_start) pairs_data.bp_rf_start = read->core.pos;
		if (bam_endpos(read) > pairs_data.bp_rf_end) pairs_data.bp_rf_end = bam_endpos(read);
		
		pairs_data.pairs_pos_mqs.push_back(get_mq(read));
		pairs_data.pairs_neg_mqs.push_back(read->core.qual);

		pairs_data.pairs_pos_nms.push_back(mate_nm); // qname_to_mate_nm[std::string(bam_get_qname(read))]);
		pairs_data.pairs_neg_nms.push_back(get_nm(read));
	} else {
		if (read->core.pos < pairs_data.bp_lf_start) pairs_data.bp_lf_start = read->core.pos;
		if (bam_endpos(read) > pairs_data.bp_lf_end) pairs_data.bp_lf_end = bam_endpos(read);
		if (read->core.mpos < pairs_data.bp_rf_start) pairs_data.bp_rf_start = read->core.mpos;
		if (get_mate_endpos(read) > pairs_data.bp_rf_end) pairs_data.bp_rf_end = get_mate_endpos(read);
		
		pairs_data.pairs_pos_mqs.push_back(read->core.qual);
		pairs_data.pairs_neg_mqs.push_back(get_mq(read));

		pairs_data.pairs_pos_nms.push_back(get_nm(read));
		pairs_data.pairs_neg_nms.push_back(mate_nm); // qname_to_mate_nm[std::string(bam_get_qname(read))]);
	}
}

void calculate_ptn_ratio(std::string contig_name, std::vector<deletion_t*>& deletions, open_samFile_t* bam_file, 
	config_t& config, stats_t& stats, bool reassign_evidence, evidence_map_t* evidence_map, std::string nm_file = "") {
	
	if (deletions.empty()) return;

	std::vector<char*> regions;
	for (deletion_t* del : deletions) {
		std::stringstream ss;
		ss << contig_name << ":" << std::max(hts_pos_t(1), del->start-stats.max_is) << "-" << del->start;
		regions.push_back(strdup(ss.str().c_str()));
		ss.str("");
		ss << contig_name << ":" << std::max(hts_pos_t(1), del->end-stats.max_is) << "-" << del->end;
		regions.push_back(strdup(ss.str().c_str()));
	}

	std::vector<deletion_t*> deletions_by_start(deletions.begin(), deletions.end());
	std::sort(deletions_by_start.begin(), deletions_by_start.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->start < d2->start;
	});

	std::vector<deletion_t*> deletions_by_end(deletions.begin(), deletions.end());
	std::sort(deletions_by_end.begin(), deletions_by_end.end(), [](const deletion_t* d1, const deletion_t* d2) {
		return d1->end < d2->end;
	});

	std::unordered_map<std::string, int64_t> qname_to_mate_nm; 
	std::ifstream mateseqs_fin(nm_file);
	std::string qname, seq;
	int64_t nm;
	while (mateseqs_fin >> qname >> seq >> nm) {
		qname_to_mate_nm[qname] = nm;
	}

	std::vector<pairs_data_t> alt_pairs_data(deletions.size());
	std::vector<pairs_data_t> ref_bp1_pairs_data(deletions.size()), ref_bp2_pairs_data(deletions.size());
	std::vector<pairs_data_t> neutral_bp1_pairs_data(deletions.size()), neutral_bp2_pairs_data(deletions.size());

	int curr_del_bystart_idx = 0, curr_del_byend_idx = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read) || bam_is_rev(read)) continue;
		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read) || read->core.isize <= 0) continue;

		while (curr_del_bystart_idx < deletions_by_start.size() && deletions_by_start[curr_del_bystart_idx]->start < read->core.pos) curr_del_bystart_idx++;
		while (curr_del_byend_idx < deletions_by_end.size() && deletions_by_end[curr_del_byend_idx]->end < read->core.pos) curr_del_byend_idx++;

		hts_pos_t pair_start = std::min(read->core.pos + read->core.l_qseq/2, bam_endpos(read));
		hts_pos_t pair_end = std::max(read->core.pos + read->core.isize - read->core.l_qseq/2, read->core.mpos);
		if (read->core.isize > stats.max_is) { // supports alt allele
			for (int i = curr_del_bystart_idx; i < deletions.size() && deletions_by_start[i]->start <= pair_end; i++) {
				if (pair_start <= deletions_by_start[i]->start+5 && deletions_by_start[i]->end-5 <= pair_end && 
					deletions_by_start[i]->start-pair_start + pair_end-deletions_by_start[i]->end <= stats.max_is) {
					if (reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, deletions_by_start[i])) continue;
					add_read(alt_pairs_data[i], read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
				}
			}
		} else if (read->core.isize >= stats.min_is) {
			for (int i = curr_del_bystart_idx; i < deletions.size() && deletions_by_start[i]->start <= pair_end; i++) {
				if (pair_start <= deletions_by_start[i]->start && deletions_by_start[i]->start <= pair_end) {
					pairs_data_t& pairs_data = (read->core.isize < stats.min_is+(-deletions_by_start[i]->svlen()) ? ref_bp1_pairs_data[i] : neutral_bp1_pairs_data[i]);
					add_read(pairs_data, read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
				}
			}
			for (int i = curr_del_byend_idx; i < deletions.size() && deletions_by_end[i]->end <= pair_end; i++) {
				if (pair_start <= deletions_by_end[i]->end && deletions_by_end[i]->end <= pair_end) {
					pairs_data_t& pairs_data = (read->core.isize < stats.min_is+(-deletions_by_end[i]->svlen()) ? ref_bp2_pairs_data[i] : neutral_bp2_pairs_data[i]);
					add_read(pairs_data, read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
				}
			}
		}
	}
	for (int i = 0; i < deletions.size(); i++) {
		set_bp_pairs_info(deletions_by_start[i]->sample_info.alt_bp1.pairs_info, alt_pairs_data[i], config);
		set_bp_pairs_info(deletions_by_start[i]->sample_info.ref_bp1.pairs_info, ref_bp1_pairs_data[i], config);
		set_bp_pairs_info(deletions_by_start[i]->sample_info.neutral_bp1_pairs, neutral_bp1_pairs_data[i], config);
		set_bp_pairs_info(deletions_by_end[i]->sample_info.ref_bp2.pairs_info, ref_bp2_pairs_data[i], config);
		set_bp_pairs_info(deletions_by_end[i]->sample_info.neutral_bp2_pairs, neutral_bp2_pairs_data[i], config);
	}
}
void calculate_ptn_ratio(std::string contig_name, std::vector<duplication_t*>& duplications, open_samFile_t* bam_file, 
	config_t& config, stats_t& stats, hts_pos_t contig_len, bool reassign_evidence, evidence_map_t* evidence_map, std::string nm_file = "") {
	
	if (duplications.empty()) return;

	std::vector<char*> regions;
	for (duplication_t* dup : duplications) {
		hts_pos_t region_start = std::max(hts_pos_t(1), dup->start-stats.max_is);
		hts_pos_t region_end = std::min(contig_len, dup->start+stats.max_is);

		std::stringstream ss;
		ss << contig_name << ":" << region_start << "-" << region_end;
		regions.push_back(strdup(ss.str().c_str()));
	}

	std::vector<duplication_t*> duplications_by_start(duplications.begin(), duplications.end());
	std::sort(duplications_by_start.begin(), duplications_by_start.end(), [](const duplication_t* d1, const duplication_t* d2) {
		return d1->start < d2->start;
	});

	std::vector<duplication_t*> duplications_by_end(duplications.begin(), duplications.end());
	std::sort(duplications_by_end.begin(), duplications_by_end.end(), [](const duplication_t* d1, const duplication_t* d2) {
		return d1->end < d2->end;
	});

	std::unordered_map<std::string, int64_t> qname_to_mate_nm; 
	std::ifstream mateseqs_fin(nm_file);
	std::string qname, seq;
	int64_t nm;
	while (mateseqs_fin >> qname >> seq >> nm) {
		qname_to_mate_nm[qname] = nm;
	}

	std::unordered_map<duplication_t*, size_t> duplication_to_idx;
	for (size_t i = 0; i < duplications.size(); i++) duplication_to_idx[duplications[i]] = i;

	std::vector<pairs_data_t> alt_pairs_data(duplications.size()), ref_pairs_data(duplications.size());
	std::vector<pairs_data_t> neutral_bp1_pairs_data(duplications.size()), neutral_bp2_pairs_data(duplications.size());

	int curr_dup_dp_idx = 0;
	int curr_dup_bystart_idx = 0, curr_dup_byend_idx = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read)) continue;
		if (is_mate_unmapped(read) || !is_samechr(read) || is_samestr(read)) continue;

		while (curr_dup_dp_idx < duplications.size() && duplications_by_start[curr_dup_dp_idx]->start < read->core.pos-stats.max_is) curr_dup_dp_idx++;
		while (curr_dup_bystart_idx < duplications.size() && duplications_by_start[curr_dup_bystart_idx]->start < read->core.pos) curr_dup_bystart_idx++;
		while (curr_dup_byend_idx < duplications.size() && duplications_by_end[curr_dup_byend_idx]->end < read->core.pos) curr_dup_byend_idx++;

		hts_pos_t pair_start = read->core.pos + read->core.l_qseq/2, pair_end = read->core.mpos + read->core.l_qseq/2;
		if (bam_is_rev(read) && is_outward(read)) {
			for (int i = curr_dup_dp_idx; i < duplications.size() && bam_endpos(read) > duplications_by_start[i]->start; i++) {
				duplication_t* dup = duplications_by_start[i];
				int idx = duplication_to_idx[dup];
				if (dup->start < pair_start && pair_start < dup->start+stats.max_is && dup->end-stats.max_is < pair_end && pair_end < dup->end) {
					if (reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, dup)) continue;
					add_read(alt_pairs_data[idx], read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
				}
			}
		} else if (!bam_is_rev(read) && read->core.isize >= stats.min_is && read->core.isize <= stats.max_is) {
			for (int i = curr_dup_bystart_idx; i < duplications.size() && duplications_by_start[i]->start < pair_end; i++) {
				duplication_t* dup = duplications_by_start[i];
				int idx = duplication_to_idx[dup];
				if (pair_start <= dup->start && dup->end <= pair_end) {
					if (read->core.isize > stats.max_is-dup->svlen()) {
						add_read(ref_pairs_data[idx], read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
					} else {
						add_read(neutral_bp1_pairs_data[idx], read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
						add_read(neutral_bp2_pairs_data[idx], read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
					}
				} else if (pair_start <= dup->start && dup->start <= pair_end) {
					add_read(neutral_bp1_pairs_data[idx], read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
				}
			}
			for (int i = curr_dup_byend_idx; i < duplications.size() && duplications_by_end[i]->end < pair_end; i++) {
				duplication_t* dup = duplications_by_end[i];
				int idx = duplication_to_idx[dup];
				if (pair_start <= dup->start && dup->end <= pair_end) {
					continue;
				} else if (pair_start <= dup->end && dup->end <= pair_end) {
					add_read(neutral_bp2_pairs_data[idx], read, qname_to_mate_nm[std::string(bam_get_qname(read))]);
				}
			}
		}
	}
	for (int i = 0; i < duplications.size(); i++) {
		duplication_t* dup = duplications[i];
		int idx = duplication_to_idx[dup];
		set_bp_pairs_info(dup->sample_info.alt_bp1.pairs_info, alt_pairs_data[idx], config);
		set_bp_pairs_info(dup->sample_info.ref_bp1.pairs_info, ref_pairs_data[idx], config);
		set_bp_pairs_info(dup->sample_info.neutral_bp1_pairs, neutral_bp1_pairs_data[idx], config);
		set_bp_pairs_info(dup->sample_info.neutral_bp2_pairs, neutral_bp2_pairs_data[idx], config);
	}
}

void calculate_ptn_ratio(std::string contig_name, std::vector<insertion_t*>& insertions, open_samFile_t* bam_file, config_t& config, stats_t& stats,
	bool reassign_evidence, evidence_map_t* evidence_map, ext_mate_map_t& mateseqs_w_mapq_chr) {
	
	if (insertions.empty()) return;

	std::vector<char*> regions;
    for (insertion_t* ins : insertions) {
        std::stringstream ss;
        ss << ins->chr << ":" << std::max(hts_pos_t(1), ins->start-stats.max_is) << "-" << ins->end+stats.max_is;
        regions.push_back(strdup(ss.str().c_str()));
    }

	std::vector<insertion_t*> insertions_by_start(insertions.begin(), insertions.end());
	std::sort(insertions_by_start.begin(), insertions_by_start.end(), [](const insertion_t* i1, const insertion_t* i2) {
		return i1->start < i2->start;
	});

	std::vector<insertion_t*> insertions_by_end(insertions.begin(), insertions.end());
	std::sort(insertions_by_end.begin(), insertions_by_end.end(), [](const insertion_t* i1, const insertion_t* i2) {
		return i1->end < i2->end;
	});

	std::unordered_map<insertion_t*, size_t> insertion_to_idx;
	for (size_t i = 0; i < insertions.size(); i++) insertion_to_idx[insertions[i]] = i;

	std::vector<pairs_data_t> alt_bp1_pairs_data(insertions.size()), alt_bp2_pairs_data(insertions.size());
	std::vector<pairs_data_t> ref_bp1_pairs_data(insertions.size());
	std::vector<pairs_data_t> neutral_bp1_pairs_data(insertions.size());
	std::vector<pairs_data_t> stray_bp1_pairs_data(insertions.size()), stray_bp2_pairs_data(insertions.size());	

	std::unordered_map<std::string, std::vector<std::pair<pairs_data_t*, bam1_t*> > > pairs_waiting_for_neg_nm;

	int curr_del_bystart_idx = 0, curr_del_byend_idx = 0;
    hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
    bam1_t* read = bam_init1();
    while (sam_itr_next(bam_file->file, iter, read) >= 0) {
        if (is_unmapped(read) || !is_primary(read)) continue;

        while (curr_del_bystart_idx < insertions_by_start.size() && insertions_by_start[curr_del_bystart_idx]->start < read->core.pos) curr_del_bystart_idx++;
		while (curr_del_byend_idx < insertions_by_end.size() && insertions_by_end[curr_del_byend_idx]->end < read->core.pos-stats.max_is) curr_del_byend_idx++;

        std::string qname = get_mate_lookup_qname(read);
		
		if (mateseqs_w_mapq_chr.count(qname) > 0) {
			std::string mate_seq = mateseqs_w_mapq_chr[qname].seq;
			int mate_seq_len = mate_seq.length();
			int min_aligned_len = (mate_seq_len+1)/2;
			
			if (!bam_is_rev(read)) {
				rc(mate_seq);
				for (int i = curr_del_bystart_idx; i < insertions_by_start.size() && insertions_by_start[i]->start < read->core.pos+stats.max_is; i++) {
					insertion_t* ins = insertions_by_start[i];
					int idx = insertion_to_idx[ins];

					ungapped_aln_t aln = best_ungapped_aln(mate_seq.c_str(), mate_seq_len, ins->ins_seq.c_str(), ins->ins_seq.length(), mate_seq_len/2, 1, -4);
					int lc_size = aln.query_begin;
					int rc_size = mate_seq_len - aln.query_end;
					int aligned_len = aln.query_end - aln.query_begin;

					if (aligned_len >= min_aligned_len && aln.mismatch_rate() <= config.max_seq_error &&
						(lc_size < config.min_clip_len || aln.ref_begin == 0) &&
						(rc_size < config.min_clip_len || aln.ref_end >= ins->ins_seq.length())) {
						if (reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, ins)) continue;
						add_read(alt_bp1_pairs_data[idx], read, aln.mismatches);
					} else {
						add_read(stray_bp1_pairs_data[idx], read, aln.mismatches);
					}
				}
			} else {
				for (int i = curr_del_byend_idx; i < insertions.size() && insertions_by_end[i]->end < read->core.pos; i++) {
					insertion_t* ins = insertions_by_end[i];
					int idx = insertion_to_idx[ins];

					ungapped_aln_t aln = best_ungapped_aln(mate_seq.c_str(), mate_seq_len, ins->ins_seq.c_str(), ins->ins_seq.length(), mate_seq_len/2, 1, -4);
					int lc_size = aln.query_begin;
					int rc_size = mate_seq_len - aln.query_end;
					int aligned_len = aln.query_end - aln.query_begin;

					if (aligned_len >= min_aligned_len && aln.mismatch_rate() <= config.max_seq_error &&
						(lc_size < config.min_clip_len || aln.ref_begin == 0) &&
						(rc_size < config.min_clip_len || aln.ref_end >= ins->ins_seq.length())) {
						if (reassign_evidence && evidence_map->is_read_assigned_to_different_sv(read, ins)) continue;
						add_read(alt_bp2_pairs_data[idx], read, aln.mismatches);
					} else {
						add_read(stray_bp2_pairs_data[idx], read, aln.mismatches);
					}
				}	
			}
		} else if (is_proper_pair(read, stats.min_is, stats.max_is)) {
			if (!bam_is_rev(read)) {
				hts_pos_t pair_start = read->core.pos + read->core.l_qseq/2, pair_end = read->core.mpos + read->core.l_qseq/2;
				for (int i = curr_del_bystart_idx; i < insertions.size() && insertions_by_start[i]->start <= pair_end; i++) {
					if (pair_start <= insertions_by_start[i]->start && insertions_by_start[i]->start <= pair_end) {
						insertion_t* ins = insertions_by_start[i];
						int idx = insertion_to_idx[ins];
						pairs_data_t* pairs_data = (read->core.isize > stats.max_is-ins->svlen() ? &ref_bp1_pairs_data[idx] : &neutral_bp1_pairs_data[idx]);
						pairs_waiting_for_neg_nm[std::string(bam_get_qname(read))].push_back({pairs_data, bam_dup1(read)});
					}
				}
			} else {
				std::vector<std::pair<pairs_data_t*, bam1_t*> >& pairs = pairs_waiting_for_neg_nm[std::string(bam_get_qname(read))];
				for (int i = 0; i < pairs.size(); i++) {
					add_read(*pairs[i].first, pairs[i].second, get_nm(read));
					bam_destroy1(pairs[i].second);
				}
				pairs_waiting_for_neg_nm.erase(std::string(bam_get_qname(read)));
			}
		}
    }
    for (int i = 0; i < insertions.size(); i++) {
        insertion_t* ins = insertions[i];
		int idx = insertion_to_idx[ins];
		set_bp_pairs_info(ins->sample_info.alt_bp1.pairs_info, alt_bp1_pairs_data[idx], config);
		set_bp_pairs_info(ins->sample_info.alt_bp2.pairs_info, alt_bp2_pairs_data[idx], config);
		set_bp_pairs_info(ins->sample_info.bp1_stray_pairs, stray_bp1_pairs_data[idx], config);
		set_bp_pairs_info(ins->sample_info.bp2_stray_pairs, stray_bp2_pairs_data[idx], config);
		set_bp_pairs_info(ins->sample_info.ref_bp1.pairs_info, ref_bp1_pairs_data[idx], config);
		set_bp_pairs_info(ins->sample_info.neutral_bp1_pairs, neutral_bp1_pairs_data[idx], config);
    }

    for (char* region : regions) free(region);
    hts_itr_destroy(iter);

}
void calculate_ptn_ratio(std::string contig_name, std::vector<inversion_t*>& inversions, open_samFile_t* bam_file, config_t& config, stats_t& stats) {
	if (inversions.empty()) return;

	std::vector<sv_t*> svs(inversions.begin(), inversions.end());
	calculate_ptn_ratio(contig_name, svs, bam_file, config, stats);

	std::vector<char*> regions;
	for (inversion_t* inv : inversions) {
		std::stringstream ss;
		ss << contig_name << ":" << inv->inv_start << "-" << inv->inv_start+stats.max_is;
		regions.push_back(strdup(ss.str().c_str()));
		ss.str("");

		ss << contig_name << ":" << std::max(hts_pos_t(1), inv->inv_end-stats.max_is) << "-" << inv->inv_end;
		regions.push_back(strdup(ss.str().c_str()));
		ss.str("");
	}

	std::vector<inversion_t*> inversions_by_start(inversions.begin(), inversions.end());
	std::sort(inversions_by_start.begin(), inversions_by_start.end(), [](const inversion_t* i1, const inversion_t* i2) {
		return i1->inv_start < i2->inv_start;
	});
	std::vector<inversion_t*> inversions_by_end(inversions.begin(), inversions.end());
	std::sort(inversions_by_end.begin(), inversions_by_end.end(), [](const inversion_t* i1, const inversion_t* i2) {
		return i1->inv_end < i2->inv_end;
	});

	std::vector<std::vector<int> > supp_pairs_bp1_pos_mqs(inversions.size()), supp_pairs_bp1_neg_mqs(inversions.size());
	std::vector<std::vector<int> > supp_pairs_bp2_pos_mqs(inversions.size()), supp_pairs_bp2_neg_mqs(inversions.size());
	
	int by_start_idx = 0, by_end_idx = 0;
	hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions.data(), regions.size());
	bam1_t* read = bam_init1();
	while (sam_itr_next(bam_file->file, iter, read) >= 0) {
		if (is_unmapped(read) || !is_primary(read) || !is_samechr(read) || !is_samestr(read)) continue;

		if (bam_is_rev(read)) {
			while (by_start_idx < inversions.size() && inversions_by_start[by_start_idx]->inv_start+stats.max_is < bam_endpos(read)) by_start_idx++;
		} else {
			while (by_end_idx < inversions.size() && inversions_by_end[by_end_idx]->inv_end < read->core.pos) by_end_idx++;
		}

		int curr_idx = bam_is_rev(read) ? by_start_idx : by_end_idx;
		for (int i = curr_idx; i < inversions.size(); i++) {
			inversion_t* inv = bam_is_rev(read) ? inversions_by_start[i] : inversions_by_end[i];

			if (bam_is_rev(read) && inv->inv_start+stats.max_is < read->core.pos) break;
			if (!bam_is_rev(read) && inv->inv_end < read->core.pos) break;

			hts_pos_t mate_startpos = read->core.mpos, mate_endpos = get_mate_endpos(read);
			if (!bam_is_mrev(read) && 
				overlap(read->core.pos, bam_endpos(read), inv->inv_end-stats.max_is, inv->inv_end) >= stats.read_len/2 &&
				overlap(mate_startpos, mate_endpos, inv->inv_start-stats.max_is, inv->inv_start) >= stats.read_len/2 &&
				read->core.mpos < read->core.pos) {
				supp_pairs_bp1_pos_mqs[i].push_back(get_mq(read));
				supp_pairs_bp1_neg_mqs[i].push_back(read->core.qual);
			}
			if (bam_is_mrev(read) && 
				overlap(read->core.pos, bam_endpos(read), inv->inv_start, inv->inv_start+stats.max_is) >= stats.read_len/2 &&
				overlap(mate_startpos, mate_endpos, inv->inv_end, inv->inv_end+stats.max_is) >= stats.read_len/2 &&
				read->core.mpos > read->core.pos) {
				supp_pairs_bp2_pos_mqs[i].push_back(read->core.qual);
				supp_pairs_bp2_neg_mqs[i].push_back(get_mq(read));
			}
		}
	}

	std::vector<int> empty_nms;
	for (int i = 0; i < inversions.size(); i++) {
		set_bp_pairs_info(inversions_by_end[i]->sample_info.alt_bp1.pairs_info, supp_pairs_bp1_pos_mqs[i], supp_pairs_bp1_neg_mqs[i], empty_nms, empty_nms, config);
		set_bp_pairs_info(inversions_by_start[i]->sample_info.alt_bp2.pairs_info, supp_pairs_bp2_pos_mqs[i], supp_pairs_bp2_neg_mqs[i], empty_nms, empty_nms, config);
	}
}


#endif //STAT_TESTS_H
