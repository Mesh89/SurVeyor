#ifndef REMAPPING_H
#define REMAPPING_H

#include <memory>
#include <vector>
#include <string>
#include <random>
#include <htslib/sam.h>

#include "../libs/ssw_selected.h"
#include "sam_utils.h"
#include "dc_remapper.h"

struct region_score_t {
    int total_score = 0;
    int remap_start = 0, remap_end = 0;
    int ro_accepted = 0, lo_accepted = 0;
};
bool operator > (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score > s2.total_score;
}
bool operator >= (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score >= s2.total_score;
}
bool operator < (const region_score_t& s1, const region_score_t& s2) {
    return s1.total_score < s2.total_score;
}

struct region_t {
    int contig_id; // id in our own mapping
    int original_bam_id; // id in the bam file
    int start, end;
    region_score_t score;

    region_t(int contig_id, int original_bam_id, int start, int end)
            : contig_id(contig_id), original_bam_id(original_bam_id), start(start), end(end) {}
};
bool operator > (const region_t& r1, const region_t& r2) {
	if (r1.score.total_score != r2.score.total_score) {
		return r1.score.total_score > r2.score.total_score;
	}
	return r1.end-r1.start < r2.end-r2.start;
}

struct remap_info_t {
    int start = 0, end = 0;
    int score = 0;
    std::string cigar;
    bool accepted = false, left_clipped = false, right_clipped = false;

    remap_info_t() {}

    remap_info_t(StripedSmithWaterman::Alignment& aln, bool accepted, int min_clip_len) : start(aln.ref_begin), end(aln.ref_end),
    score(aln.sw_score), cigar(aln.cigar_string), accepted(accepted) {
        if (aln.cigar.empty()) return;
        uint32_t f = aln.cigar[0], l = aln.cigar[aln.cigar.size()-1];
        left_clipped = cigar_int_to_op(f) == 'S' && cigar_int_to_len(f) >= min_clip_len;
        right_clipped = cigar_int_to_op(l) == 'S' && cigar_int_to_len(l) >= min_clip_len;
    }
};

struct remap_side_t {
    std::vector<remap_info_t> read_remap_infos;
    std::unique_ptr<remap_info_t> clip_remap_info;
};

region_t get_candidate_region(std::vector<std::shared_ptr<bam1_t>> subcluster, std::string& m_contig_name, int m_contig_id, hts_pos_t m_contig_len, int max_is, int max_trans_size) {
    std::shared_ptr<bam1_t> leftmost_reverse_mate = NULL, rightmost_reverse_mate = NULL;
    std::shared_ptr<bam1_t> leftmost_forward_mate = NULL, rightmost_forward_mate = NULL;
    /* gets two regions
     * 1. from reverse mates, max_is to the left of right-most read and max_insertion_size to the right of the left-most
     * 2. from right-most forward read, max_is to the right and max_insertion_size to the left
     */
    for (std::shared_ptr<bam1_t> r : subcluster) { // get leftmost reverse read
        if (bam_is_mrev(r)) {
            if (leftmost_reverse_mate == NULL || leftmost_reverse_mate->core.mpos > r->core.mpos) {
                leftmost_reverse_mate = r;
            }
            if (rightmost_reverse_mate == NULL || rightmost_reverse_mate->core.mpos < r->core.mpos) {
                rightmost_reverse_mate = r;
            }
        }
    }
    for (std::shared_ptr<bam1_t> r : subcluster) { // get rightmost forward read
        if (!bam_is_mrev(r)) {
            if (leftmost_forward_mate == NULL || leftmost_forward_mate->core.mpos > r->core.mpos) {
                leftmost_forward_mate = r;
            }
            if (rightmost_forward_mate == NULL || rightmost_forward_mate->core.mpos < r->core.mpos) {
                rightmost_forward_mate = r;
            }
        }
    }

    hts_pos_t start = INT_MAX;
    hts_pos_t end = 0;
    if (leftmost_reverse_mate != NULL) {
        start = std::min(start, rightmost_reverse_mate->core.mpos-max_is);
        end = std::max(end, leftmost_reverse_mate->core.mpos+max_trans_size);
    }
    if (rightmost_forward_mate != NULL) {
        start = std::min(start, rightmost_forward_mate->core.mpos-max_trans_size);
        end = std::max(end, leftmost_forward_mate->core.mpos+max_is);
    }

    // hts_pos_t contig_len = contigs.get_len(m_contig_name);
    return region_t(m_contig_id, subcluster[0]->core.mtid, std::max(hts_pos_t(0), start), std::min(end, m_contig_len));
}

remap_info_t remap_seq(char* region, int region_len, std::string& seq, StripedSmithWaterman::Aligner& aligner,
                StripedSmithWaterman::Filter& filter, int min_clip_len) {
    if (region_len <= 0) {
        return remap_info_t();
    }

    StripedSmithWaterman::Alignment alignment;
    aligner.Align(seq.c_str(), region, region_len, filter, &alignment, 0);
    bool accepted = 2*alignment.sw_score >= seq.length() && alignment.query_end-alignment.query_begin >= min_clip_len; 

    return remap_info_t(alignment, accepted, min_clip_len);
}

// if remapped reads are too far away, find the max_is bp window where the sum of read scores is maximised, and remap all the reads there
void remap_to_maxis_window(std::vector<remap_info_t>& remap_infos, remap_info_t& clip_remap_info, char* region_ptr,
		int remapped_region_start, int remapped_region_end, std::vector<std::string>& read_seqs, std::string& clip_seq, char clip_dir, int max_is, int min_clip_len,
		StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Filter& filter) {
	std::sort(remap_infos.begin(), remap_infos.end(),
			[](const remap_info_t& ri1, const remap_info_t& ri2) { return ri1.start < ri2.start; });
	int s = 0, score = 0;
	int best_start = 0, best_end = 0, best_score = 0;
	for (int e = 0; e < remap_infos.size(); e++) {
		while (remap_infos[e].start-remap_infos[s].start > max_is) {
			score -= remap_infos[s].score;
			s++;
		}
		score += remap_infos[e].score;
		if (score > best_score) {
			best_score = score;
			best_start = s, best_end = e;
		}
	}

	int strictest_region_start = remap_infos[best_start].start, strictest_region_end = remap_infos[best_end].end;
	int strictest_region_len = strictest_region_end - strictest_region_start;
	remap_infos.clear();


	// add some additional space to remap the clip
	int padding = max_is - strictest_region_len;
	if (!clip_seq.empty() && padding > 0) {
		if (clip_dir == 'L') {
			strictest_region_start -= padding;
			// we do not allow the new strictest region to cross the original remapped region - this is the simplest way to make sure
			// it does not go out of the contig boundary and crash
			strictest_region_start = std::max(strictest_region_start, remapped_region_start);
		}
		else if (clip_dir == 'R') {
			strictest_region_end += padding;
			strictest_region_end = std::min(strictest_region_end, remapped_region_end);
		}
		strictest_region_len = strictest_region_end - strictest_region_start;
	}

	for (std::string& s : read_seqs) {
		remap_info_t remap_info = remap_seq(region_ptr+strictest_region_start, strictest_region_len, s, aligner, filter, min_clip_len);
		remap_info.start += strictest_region_start;
		remap_info.end += strictest_region_start;
		remap_infos.push_back(remap_info);
	}
	if (!clip_seq.empty()) {
		clip_remap_info = remap_seq(region_ptr+strictest_region_start, strictest_region_len, clip_seq, aligner, filter, min_clip_len);
		clip_remap_info.start += strictest_region_start;
		clip_remap_info.end += strictest_region_start;
	}
}

int get_strict_region_start(std::vector<remap_info_t>& ro_remap_infos, std::vector<remap_info_t>& lo_remap_infos) {
    std::vector<remap_info_t> remap_infos;
    for (remap_info_t& ri : ro_remap_infos) {
        if (ri.left_clipped) remap_infos.push_back(ri);
    }
    for (remap_info_t& ri : lo_remap_infos) {
        if (ri.left_clipped) remap_infos.push_back(ri);
    }
    if (remap_infos.empty()) return 0;

    remap_info_t guard; guard.start = INT32_MAX;
    remap_infos.push_back(guard);

    std::sort(remap_infos.begin(), remap_infos.end(),
              [](remap_info_t& ri1, remap_info_t& ri2) { return ri1.start < ri2.start; });

    int strict_region_start = remap_infos[0].start, max_freq = 1, curr_freq = 1;
    for (int i = 1; i < remap_infos.size(); i++) {
        if (remap_infos[i-1].start != remap_infos[i].start) {
            if (curr_freq > max_freq) {
                strict_region_start = remap_infos[i-1].start;
                max_freq = curr_freq;
            }
            curr_freq = 0;
        }
        curr_freq++;
    }

    return strict_region_start;
}
int get_strict_region_end(std::vector<remap_info_t>& rc_remap_infos, std::vector<remap_info_t>& lc_remap_infos) {
    std::vector<remap_info_t> remap_infos;
    for (remap_info_t& ri : rc_remap_infos) {
        if (ri.right_clipped) remap_infos.push_back(ri);
    }
    for (remap_info_t& ri : lc_remap_infos) {
        if (ri.right_clipped) remap_infos.push_back(ri);
    }
    if (remap_infos.empty()) return 0;

    remap_info_t guard; guard.end = 0;
    remap_infos.push_back(guard);

    std::sort(remap_infos.begin(), remap_infos.end(),
              [](remap_info_t& ri1, remap_info_t& ri2) { return ri1.end > ri2.end; });

    int strict_region_end = remap_infos[0].end, max_freq = 1, curr_freq = 1;
    for (int i = 1; i < remap_infos.size(); i++) {
        if (remap_infos[i-1].end != remap_infos[i].end) {
            if (curr_freq > max_freq) {
                strict_region_end = remap_infos[i-1].end;
                max_freq = curr_freq;
            }
            curr_freq = 0;
        }
        curr_freq++;
    }

    return strict_region_end;
}


region_score_t compute_score_supp(region_t& region, char* contig_seq, std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
                       std::unordered_map<std::string, std::string>& mateseqs,
                       remap_side_t& ro_remap_side, remap_side_t& lo_remap_side,
                       StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner,
                       bool do_rc, config_t& config, stats_t& stats) {

    StripedSmithWaterman::Filter filter;

    // mateseqs contains seqs are stored as in fasta/q file
    std::vector<std::string> r_mates, l_mates;
    for (std::shared_ptr<bam1_t> read : r_cluster->cluster->reads) {
        std::string s = get_mate_seq(read.get(), mateseqs);
        if (!do_rc) rc(s);   // by default, we map mates of fwd reads on the rev strand
        r_mates.push_back(s);
    }
    for (std::shared_ptr<bam1_t> read : l_cluster->cluster->reads) {
        std::string s = get_mate_seq(read.get(), mateseqs);
        if (do_rc) rc(s);   // by default, we map mates of rev reads on the fwd strand
        l_mates.push_back(s);
    }

    char* region_ptr = contig_seq + region.start;
    int strict_region_start = 0, strict_region_end = region.end - region.start;

    region_score_t score;
    remap_info_t rc_remap_info, lc_remap_info;

    // map the clips (if available)
    if (r_cluster->clip_consensus && r_cluster->clip_consensus->clip_len >= config.min_clip_len) {
        std::string s = r_cluster->clip_consensus->clip_sequence();
        if (do_rc) rc(s);
        remap_info_t remap_info = remap_seq(region_ptr, region.end-region.start, s, aligner, filter, config.min_clip_len);
        if (remap_info.accepted) {
            if (!do_rc) {
                strict_region_start = remap_info.start;
            } else if (do_rc) {
                strict_region_end = remap_info.end;
            }
        }
        rc_remap_info = remap_info;
    }
    if (l_cluster->clip_consensus && l_cluster->clip_consensus->clip_len >= config.min_clip_len) {
        std::string s = l_cluster->clip_consensus->clip_sequence();
        if (do_rc) rc(s);
        remap_info_t remap_info = remap_seq(region_ptr, region.end-region.start, s, aligner, filter, config.min_clip_len);
        if (remap_info.accepted) {
            if (do_rc) {
                strict_region_start = remap_info.start;
            } else if (!do_rc) {
                strict_region_end = remap_info.end;
            }
        }
        lc_remap_info = remap_info;
    }

    // map the reads
    int strict_region_len = strict_region_end - strict_region_start;

    int min_start = INT32_MAX, max_end = 0;
    for (std::string& s : r_mates) {
        remap_info_t remap_info = remap_seq(region_ptr+strict_region_start, strict_region_len, s, aligner, filter, config.min_clip_len);
        remap_info.start += strict_region_start;
        min_start = std::min(min_start, remap_info.start);
        remap_info.end += strict_region_start;
        max_end = std::max(max_end, remap_info.end);
        ro_remap_side.read_remap_infos.push_back(remap_info);
    }
    if (max_end-min_start > stats.max_is) { // all the mappings need to be within a config.max_is bp window
    	std::string clip_seq = r_cluster->clip_consensus ? r_cluster->clip_consensus->clip_sequence() : "";
    	char clip_dir = do_rc ? 'R' : 'L'; // R means that the clip determines the right-most boundary of the inserted region/sequence, L the left-most
    	remap_to_maxis_window(ro_remap_side.read_remap_infos, rc_remap_info, region_ptr, strict_region_start, strict_region_end,
    			r_mates, clip_seq, clip_dir, stats.max_is, config.min_clip_len, aligner, filter);
    }

    min_start = INT32_MAX, max_end = 0;
    for (std::string& s : l_mates) {
        remap_info_t remap_info = remap_seq(region_ptr+strict_region_start, strict_region_len, s, aligner, filter, config.min_clip_len);
        remap_info.start += strict_region_start;
        min_start = std::min(min_start, remap_info.start);
        remap_info.end += strict_region_start;
        max_end = std::max(max_end, remap_info.end);
        lo_remap_side.read_remap_infos.push_back(remap_info);
    }
    if (max_end-min_start > stats.max_is) { // all the mappings need to be within a config.max_is bp window
    	std::string clip_seq = l_cluster->clip_consensus ? l_cluster->clip_consensus->clip_sequence() : "";
    	char clip_dir = do_rc ? 'L' : 'R';
    	remap_to_maxis_window(lo_remap_side.read_remap_infos, lc_remap_info, region_ptr, strict_region_start, strict_region_end,
    			l_mates, clip_seq, clip_dir, stats.max_is, config.min_clip_len, aligner, filter);
    }

    // if strict limits not identified through remapping of the clips, choose most common start and most common ends
    if (strict_region_start == 0) {
        strict_region_start = get_strict_region_start(ro_remap_side.read_remap_infos, lo_remap_side.read_remap_infos);
    }
    if (strict_region_end == region.end-region.start) {
        int temp = get_strict_region_end(ro_remap_side.read_remap_infos, lo_remap_side.read_remap_infos);
        if (temp > 0) strict_region_end = temp;
    }

    // un-accept all remappings that fall outside the strict limits or they are not clipped correctly
    for (remap_info_t& remap_info : ro_remap_side.read_remap_infos) {
        if (remap_info.left_clipped && strict_region_start < remap_info.start-5) remap_info.accepted = false;
        if (remap_info.start < strict_region_start-5) remap_info.accepted = false;
        if (remap_info.right_clipped && remap_info.end+5 < strict_region_end) remap_info.accepted = false;
        if (strict_region_end < remap_info.start) remap_info.accepted = false;
    }
    for (remap_info_t& remap_info : lo_remap_side.read_remap_infos) {
        if (remap_info.left_clipped && strict_region_start < remap_info.start-5) remap_info.accepted = false;
        if (remap_info.end < strict_region_start) remap_info.accepted = false;
        if (remap_info.right_clipped && remap_info.end+5 < strict_region_end) remap_info.accepted = false;
        if (strict_region_end < remap_info.start) remap_info.accepted = false;
    }

    for (remap_info_t& remap_info : ro_remap_side.read_remap_infos) {
        if (remap_info.accepted) {
            score.total_score += remap_info.score;
            score.ro_accepted++;
        }
    }
    for (remap_info_t& remap_info : lo_remap_side.read_remap_infos) {
        if (remap_info.accepted) {
            score.total_score += remap_info.score;
            score.lo_accepted++;
        }
    }

    // find inserted sequence coordinates
    int rc_start = INT32_MAX, rc_end = 0;
		if (rc_remap_info.accepted) {
			rc_start = rc_remap_info.start;
			rc_end = rc_remap_info.end;
		} else {
			for (remap_info_t& ri : ro_remap_side.read_remap_infos) {
				if (ri.accepted) {
					rc_start = std::min(rc_start, ri.start);
					rc_end = std::max(rc_end, ri.end);
			}
		}
	}
	int lc_start = INT32_MAX, lc_end = 0;
		if (lc_remap_info.accepted) {
			lc_start = lc_remap_info.start;
			lc_end = lc_remap_info.end;
		} else {
			for (remap_info_t& ri : lo_remap_side.read_remap_infos) {
				if (ri.accepted) {
					lc_start = std::min(lc_start, ri.start);
					lc_end = std::max(lc_end, ri.end);
			}
		}
	}

	if (!do_rc) {
		score.remap_start = rc_start;
		score.remap_end = lc_end;
	} else {
		score.remap_start = lc_start;
		score.remap_end = rc_end;
	}

    if (r_cluster->clip_consensus) {
        ro_remap_side.clip_remap_info.reset(new remap_info_t(rc_remap_info));
        if (rc_remap_info.accepted) {
            score.total_score += rc_remap_info.score;
            score.ro_accepted++;
        }
    }
    if (l_cluster->clip_consensus) {
        lo_remap_side.clip_remap_info.reset(new remap_info_t(lc_remap_info));
        if (lc_remap_info.accepted) {
            score.total_score += lc_remap_info.score;
            score.lo_accepted++;
        }
    }
    return score;
}

void compute_score(region_t& region, char* contig_seq, std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
                   std::unordered_map<std::string, std::string>& mateseqs,
                   remap_side_t* ro_remap_side, remap_side_t* lo_remap_side,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner,
                   bool& is_rc, config_t& config, stats_t& stats) {

    remap_side_t fwd_ro_remap_side, fwd_lo_remap_side, rc_ro_remap_side, rc_lo_remap_side;
    region_score_t score = compute_score_supp(region, contig_seq, r_cluster, l_cluster, mateseqs, fwd_ro_remap_side, fwd_lo_remap_side,
                                              aligner,permissive_aligner, false, config, stats);
    region_score_t rc_score = compute_score_supp(region, contig_seq, r_cluster, l_cluster, mateseqs, rc_ro_remap_side, rc_lo_remap_side,
                                                 aligner, permissive_aligner, true, config, stats);
    if (score >= rc_score) {
        is_rc = false;
        region.score = score;
        if (ro_remap_side != NULL) {
            ro_remap_side->read_remap_infos.swap(fwd_ro_remap_side.read_remap_infos);
            ro_remap_side->clip_remap_info.swap(fwd_ro_remap_side.clip_remap_info);
        }
        if (lo_remap_side != NULL) {
            lo_remap_side->read_remap_infos.swap(fwd_lo_remap_side.read_remap_infos);
            lo_remap_side->clip_remap_info.swap(fwd_lo_remap_side.clip_remap_info);
        }
    } else {
        is_rc = true;
        region.score = rc_score;
        if (ro_remap_side != NULL) {
            ro_remap_side->read_remap_infos.swap(rc_ro_remap_side.read_remap_infos);
            ro_remap_side->clip_remap_info.swap(rc_ro_remap_side.clip_remap_info);
        }
        if (lo_remap_side != NULL) {
            lo_remap_side->read_remap_infos.swap(rc_lo_remap_side.read_remap_infos);
            lo_remap_side->clip_remap_info.swap(rc_lo_remap_side.clip_remap_info);
        }
    }
}


std::shared_ptr<insertion_cluster_t> subsample_cluster(std::shared_ptr<insertion_cluster_t> reads_cluster, int size, std::default_random_engine& rng) {
    std::vector<std::shared_ptr<bam1_t>> subset(reads_cluster->cluster->reads);
    if (subset.size() > size) {
        std::shuffle(subset.begin(), subset.end(), rng);
        subset.erase(subset.begin() + size, subset.end());
    }
    auto subsampled_cluster = std::make_shared<insertion_cluster_t>(std::make_shared<cluster_t>());
    for (std::shared_ptr<bam1_t> read : subset) subsampled_cluster->add_stable_read(read);
    subsampled_cluster->add_clip_cluster(reads_cluster->clip_consensus);
    return subsampled_cluster;
}

std::pair<region_t, region_t> compute_best_and_base_regions(int contig_id, std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster, std::vector<region_t>& regions,
                                                            std::unordered_map<std::string, std::string>& mateseqs, chr_seqs_map_t& contigs, contig_map_t& contig_map, 
                                                            StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner, StripedSmithWaterman::Aligner& aligner_to_base,
                                                            config_t& config, stats_t& stats) {

    const int SMALL_SAMPLE_SIZE = 15;
    const int CLUSTER_CANDIDATES = 3;

    bool is_rc;

    // if too many regions and too many reads, subsample
    if (r_cluster->cluster->count + l_cluster->cluster->count > 2*SMALL_SAMPLE_SIZE && regions.size() > CLUSTER_CANDIDATES) {
    	std::default_random_engine rng(config.seed);
        std::shared_ptr<insertion_cluster_t> subsampled_lc = subsample_cluster(l_cluster, SMALL_SAMPLE_SIZE, rng);
        std::shared_ptr<insertion_cluster_t> subsampled_rc = subsample_cluster(r_cluster, SMALL_SAMPLE_SIZE, rng);

        // compute best score
        for (int i = 0; i < regions.size(); i++) {
            char* contig_seq = contigs.get_seq(contig_map.get_name(regions[i].contig_id));
            compute_score(regions[i], contig_seq, subsampled_rc, subsampled_lc, mateseqs, NULL, NULL,
                          aligner, permissive_aligner, is_rc, config, stats);
        }
        std::sort(regions.begin(), regions.end(), std::greater<region_t>());

        regions.erase(regions.begin()+CLUSTER_CANDIDATES, regions.end());
    }

    // compute best score
    for (region_t& region : regions) {
        char* contig_seq = contigs.get_seq(contig_map.get_name(region.contig_id));
        compute_score(region, contig_seq, r_cluster, l_cluster, mateseqs, NULL, NULL,
            aligner, permissive_aligner, is_rc, config, stats);
    }
    std::sort(regions.begin(), regions.end(), std::greater<region_t>());

    int start = r_cluster->start - stats.max_is;
    int end = l_cluster->end + stats.max_is;
    int contig_len = contigs.get_len(contig_map.get_name(contig_id));
    region_t base_region(contig_id, 0, std::max(0, start), std::min(end, contig_len));

    char* contig_seq = contigs.get_seq(contig_map.get_name(base_region.contig_id));
    compute_score(base_region, contig_seq, r_cluster, l_cluster, mateseqs, NULL, NULL,
                  aligner_to_base, permissive_aligner, is_rc, config, stats);

    return {regions[0], base_region};
}

#endif // REMAPPING_H
