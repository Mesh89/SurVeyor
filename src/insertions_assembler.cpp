#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <algorithm>
#include <memory>
#include <unistd.h>

#include <htslib/sam.h>
#include <htslib/kseq.h>
#include <vector>

#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "../libs/IntervalTree.h"
#include "dc_remapper.h"
#include "htslib/vcf.h"
#include "sam_utils.h"
#include "clustering_utils.h"
#include "remapping.h"
#include "assemble.h"
#include "reference_guided_assembly.h"
#include "vcf_utils.h"


config_t config;
stats_t stats;
std::string workdir;
std::mutex mtx;

contig_map_t contig_map;
chr_seqs_map_t contigs;

std::ofstream assembly_failed_no_seq, assembly_failed_cycle_writer, assembly_failed_too_many_reads_writer;

std::vector<std::shared_ptr<sv_t>> insertions;

const double BASE_ACCEPTANCE_THRESHOLD = 0.95;
const int TOO_MANY_READS = 1000;

struct cc_v_distance_t {
    std::shared_ptr<insertion_cluster_t> c1, c2;
    int distance;

    cc_v_distance_t(std::shared_ptr<insertion_cluster_t> c1, std::shared_ptr<insertion_cluster_t> c2, int distance) : c1(c1), c2(c2), distance(distance) {}
};
bool operator < (const cc_v_distance_t& ccd1, const cc_v_distance_t& ccd2) { // reverse op for priority queue
    return ccd1.distance < ccd2.distance;
}

bool find_insertion_from_cluster_pair(std::shared_ptr<insertion_cluster_t> r_cluster, std::shared_ptr<insertion_cluster_t> l_cluster,
                   int contig_id, bam_hdr_t* header, std::unordered_map<std::string, std::string>& mateseqs,
				   std::unordered_map<std::string, std::string>& matequals,
                   StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Aligner& permissive_aligner,
                   StripedSmithWaterman::Aligner& aligner_to_base, StripedSmithWaterman::Aligner& harsh_aligner,
                   stats_t& stats) {

	if (l_cluster->cluster->count + r_cluster->cluster->count > TOO_MANY_READS) return false;

    std::vector<region_t> regions;

    /* == Find candidate regions for remapping == */
    std::vector<std::shared_ptr<bam1_t>> full_cluster;
    full_cluster.insert(full_cluster.end(), l_cluster->cluster->reads.begin(), l_cluster->cluster->reads.end());
    full_cluster.insert(full_cluster.end(), r_cluster->cluster->reads.begin(), r_cluster->cluster->reads.end());
    std::sort(full_cluster.begin(), full_cluster.end(), [] (std::shared_ptr<bam1_t> r1, std::shared_ptr<bam1_t> r2) {
        if (r1->core.mtid != r2->core.mtid) return r1->core.mtid < r2->core.mtid;
        else return r1->core.mpos < r2->core.mpos;
    });

    std::vector<std::shared_ptr<bam1_t>> subcluster;
    for (std::shared_ptr<bam1_t> r : full_cluster) {
        if (!subcluster.empty() && (subcluster[0]->core.mtid != r->core.mtid ||
                r->core.mpos-subcluster[0]->core.mpos > stats.max_is)) {
            std::string m_contig_name = header->target_name[subcluster[0]->core.mtid];
            regions.push_back(get_candidate_region(subcluster, m_contig_name, contig_map.get_id(m_contig_name), contigs.get_len(m_contig_name), stats.max_is, config.max_trans_size));
            subcluster.clear();
        }
        if (!is_mate_unmapped(r.get())) subcluster.push_back(r);
    }
    if (!subcluster.empty()) {
        std::string m_contig_name = header->target_name[subcluster[0]->core.mtid];
        regions.push_back(get_candidate_region(subcluster, m_contig_name, contig_map.get_id(m_contig_name), contigs.get_len(m_contig_name), stats.max_is, config.max_trans_size));
    }

    std::string contig_name = contig_map.get_name(contig_id);
    if (regions.empty()) {
		std::shared_ptr<sv_t> ins = detect_de_novo_insertion(contig_name, contigs, r_cluster, l_cluster, mateseqs, matequals, assembly_failed_no_seq, assembly_failed_cycle_writer, assembly_failed_too_many_reads_writer,
				aligner_to_base, harsh_aligner, config, stats);
                
        if (ins != NULL) {
			mtx.lock();
            insertions.push_back(ins);
			mtx.unlock();
            return true;
		}
		return false;
    }

    /* == Find best candidate region == */
    auto best_and_base_region = compute_best_and_base_regions(contig_id, r_cluster, l_cluster, regions, mateseqs, contigs, 
        contig_map, aligner, permissive_aligner, aligner_to_base, config, stats);
    region_t best_region = best_and_base_region.first;
    region_t base_region = best_and_base_region.second;

    remap_side_t ro_remap_side, lo_remap_side;
    char* contig_seq = contigs.get_seq(contig_map.get_name(best_region.contig_id));
    bool is_rc;
    compute_score(best_region, contig_seq, r_cluster, l_cluster, mateseqs, &ro_remap_side, &lo_remap_side,
                  aligner, permissive_aligner, is_rc, config, stats);

    // build corrected consensus sequence and realign reads to it
    bool left_bp_precise = false, right_bp_precise = false;
    std::string corrected_consensus_sequence = generate_consensus_sequences(contig_name, contigs, contig_map, r_cluster, l_cluster, ro_remap_side, lo_remap_side,
					best_region, is_rc, left_bp_precise, right_bp_precise, mateseqs, matequals, aligner, harsh_aligner, config, stats);

    /* == Try assembly == */
	int rc_accepted_reads = 0, lc_accepted_reads = 0;
	for (int i = 0; i < r_cluster->cluster->reads.size(); i++) {
        rc_accepted_reads += ro_remap_side.read_remap_infos[i].accepted;
    }
    for (int i = 0; i < l_cluster->cluster->reads.size(); i++) {
        lc_accepted_reads += lo_remap_side.read_remap_infos[i].accepted;
    }
	int tot_reads = r_cluster->cluster->reads.size() + l_cluster->cluster->reads.size();
	if (rc_accepted_reads == 0 || lc_accepted_reads == 0 || double(rc_accepted_reads+lc_accepted_reads)/tot_reads < 0.5) {
		std::shared_ptr<sv_t> ins = detect_de_novo_insertion(contig_name, contigs, r_cluster, l_cluster, mateseqs, matequals, assembly_failed_no_seq, assembly_failed_cycle_writer, assembly_failed_too_many_reads_writer,
				aligner_to_base, harsh_aligner, config, stats);

		if (ins != NULL) {
            mtx.lock();
            insertions.push_back(ins);
			mtx.unlock();
			return true;
		}
	}

	if (corrected_consensus_sequence.empty()) return false;

    bool success = false;
    std::shared_ptr<insertion_t> insertion = detect_reference_guided_assembly_insertion(contig_name, contigs.get_seq(contig_name), contigs.get_len(contig_name), 
        corrected_consensus_sequence, r_cluster, l_cluster, ro_remap_side, lo_remap_side, best_region, is_rc, left_bp_precise, right_bp_precise, aligner, stats, config);
    if (insertion != NULL) {
        mtx.lock();
        insertions.push_back(insertion);
        mtx.unlock();
        success = true;
    }
    return success;
}

std::vector<std::shared_ptr<insertion_cluster_t> > cluster_reads(open_samFile_t* dc_file, int contig_id, std::vector<std::shared_ptr<consensus_t>>& clip_consensuses) {

	std::string contig_name = contig_map.get_name(contig_id);
    hts_itr_t* iter = sam_itr_querys(dc_file->idx, dc_file->header, contig_name.c_str());
    bam1_t* read = bam_init1();

    std::vector<std::shared_ptr<bam1_t>> reads;
    std::vector<std::shared_ptr<cluster_t>> clusters;
    while (sam_itr_next(dc_file->file, iter, read) >= 0) {
        std::shared_ptr<cluster_t> cluster = std::make_shared<cluster_t>(read, config.high_confidence_mapq);
        cluster->id = clusters.size();
        cluster->ra_start = cluster->la_start;
        cluster->ra_end = cluster->la_end;
        cluster->ra_rev = cluster->la_rev;

        clusters.push_back(cluster);
        reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
    }
    sam_itr_destroy(iter);
    bam_destroy1(read);

    if (clusters.empty()) return std::vector<std::shared_ptr<insertion_cluster_t> >();

    // union-find data structure
    int n_reads = clusters.size();    
    int max_cluster_size = (stats.get_max_depth(contig_name) * stats.max_is)/stats.read_len;
    cluster_clusters(clusters, reads, stats.max_is, max_cluster_size, true);

    std::vector<std::shared_ptr<insertion_cluster_t> > read_clusters;
    for (int i = 0; i < n_reads; i++) {
        if (clusters[i] != NULL) read_clusters.push_back(std::make_shared<insertion_cluster_t>(clusters[i]));
    }

    if (!clip_consensuses.empty()) { // TODO: temporary, enhance logic
        std::sort(clip_consensuses.begin(), clip_consensuses.end(), [](std::shared_ptr<consensus_t> c1, std::shared_ptr<consensus_t> c2) {
            return c1->breakpoint < c2->breakpoint;
        });

        int curr_j = 0;
        std::vector<std::shared_ptr<consensus_t>> clip_consensus_per_cluster(read_clusters.size(), NULL);

        bool left_facing = clip_consensuses[0]->left_clipped;
        if (left_facing) {
            std::sort(read_clusters.begin(), read_clusters.end(), [](std::shared_ptr<insertion_cluster_t> rc1, std::shared_ptr<insertion_cluster_t> rc2) {
                return rc1->end < rc2->end;
            });

            for (std::shared_ptr<consensus_t> clip_consensus : clip_consensuses) {
                while (curr_j < read_clusters.size() && read_clusters[curr_j]->end < clip_consensus->breakpoint) curr_j++;

                for (int j = curr_j; j < read_clusters.size(); j++) {
                    if (read_clusters[j]->end - clip_consensus->breakpoint <= stats.max_is) {
                        if (clip_consensus_per_cluster[j] == NULL || clip_consensus_per_cluster[j]->reads() < clip_consensus->reads()) {
                            clip_consensus_per_cluster[j] = clip_consensus;
                        }
                    } else break;
                }
            }
        } else {
            std::sort(read_clusters.begin(), read_clusters.end(), [](std::shared_ptr<insertion_cluster_t> rc1, std::shared_ptr<insertion_cluster_t> rc2) {
                return rc1->start < rc2->start;
            });

            for (std::shared_ptr<consensus_t> clip_consensus : clip_consensuses) {
                while (curr_j < read_clusters.size() && clip_consensus->breakpoint-read_clusters[curr_j]->start > stats.max_is) curr_j++;

                for (int j = curr_j; j < read_clusters.size(); j++) {
                    if (clip_consensus->breakpoint - read_clusters[j]->start >= 0) {
                        if (clip_consensus_per_cluster[j] == NULL || clip_consensus_per_cluster[j]->reads() < clip_consensus->reads()) {
                            clip_consensus_per_cluster[j] = clip_consensus;
                        }
                    } else break;
                }
            }
        }

        for (int j = 0; j < read_clusters.size(); j++) {
            read_clusters[j]->add_clip_cluster(clip_consensus_per_cluster[j]);
        }
    }

    std::sort(read_clusters.begin(), read_clusters.end(), [](const std::shared_ptr<insertion_cluster_t> rc1, const std::shared_ptr<insertion_cluster_t> rc2) {
    	return rc1->cluster->count > rc2->cluster->count;
    });

    return read_clusters;
}

bool is_semi_mapped(bam1_t* read) {
	return is_proper_pair(read, stats.min_is, stats.max_is) && !is_mate_clipped(read);
}

void add_semi_mapped_pairs(std::string clipped_fname, int contig_id, std::vector<std::shared_ptr<insertion_cluster_t> >& r_clusters, std::vector<std::shared_ptr<insertion_cluster_t> >& l_clusters) {

	if (!file_exists(clipped_fname)) return;

	std::vector<bam1_t*> l_semi_mapped_pairs, r_semi_mapped_pairs;

	open_samFile_t* clipped_file = new open_samFile_t(clipped_fname.c_str(), true);
	std::string contig_name = contig_map.get_name(contig_id);
	hts_itr_t* iter = sam_itr_querys(clipped_file->idx, clipped_file->header, contig_name.c_str());
	bam1_t* read = bam_init1();
	while (sam_itr_next(clipped_file->file, iter, read) >= 0) {
		if (is_semi_mapped(read) && get_sequence(read).find("N") == std::string::npos) {
            if (read->core.flag & BAM_FMREVERSE) {
				l_semi_mapped_pairs.push_back(bam_dup1(read));
			} else {
				r_semi_mapped_pairs.push_back(bam_dup1(read));
			}
		}
	}
	delete clipped_file;
    hts_itr_destroy(iter);
    bam_destroy1(read);

	// add semimapped pairs
	std::vector<Interval<bam1_t*>> r_iv, l_iv;
	for (bam1_t* r : r_semi_mapped_pairs) {
		r_iv.push_back(Interval<bam1_t*>(r->core.mpos, get_mate_endpos(r), r));
	}
	for (bam1_t* r : l_semi_mapped_pairs) {
		l_iv.push_back(Interval<bam1_t*>(r->core.mpos, get_mate_endpos(r), r));
	}
	IntervalTree<bam1_t*> r_it(r_iv), l_it(l_iv);
	std::unordered_set<bam1_t*> used_sm_reads;
	for (std::shared_ptr<insertion_cluster_t> rc : r_clusters) {
		int end = rc->end, start = end - stats.max_is;
		auto sm_reads = r_it.findContained(start, end);
		for (auto i : sm_reads) {
			bam1_t* r = i.value;
			if (!used_sm_reads.count(r)) {
				rc->add_semi_mapped_reads(bam_dup1(r));
				used_sm_reads.insert(r);
			}
		}
	}
	for (std::shared_ptr<insertion_cluster_t> rc : l_clusters) {
		int start, end;
		start = rc->start; end = start + stats.max_is;
		auto sm_reads = l_it.findContained(start, end);
		for (auto i : sm_reads) {
			bam1_t* r = i.value;
			if (!used_sm_reads.count(r)) {
				rc->add_semi_mapped_reads(bam_dup1(r));
				used_sm_reads.insert(r);
			}
		}
	}
	for (bam1_t* r : r_semi_mapped_pairs) {
		bam_destroy1(r);
	}
	for (bam1_t* r : l_semi_mapped_pairs) {
		bam_destroy1(r);
	}
}

int dfs(int curr_idx, std::vector<std::vector<int> >& dependencies, std::vector<int>& dfs_label, int& label) {
    dfs_label[curr_idx] = label;
    for (int next_idx : dependencies[curr_idx]) {
        if (dfs_label[next_idx] == -1) {
            dfs(next_idx, dependencies, dfs_label, label);
        }
    }
    return 0;
}

void find_insertions(int id, int contig_id, int comp_id, std::vector<cc_v_distance_t> cc_v_distances, bam_hdr_t* hdr,
                     std::shared_ptr<std::unordered_map<std::string, std::string>> mateseqs,
                     std::shared_ptr<std::unordered_map<std::string, std::string>> matequals) {

    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Aligner permissive_aligner(2, 2, 4, 1, false);
    StripedSmithWaterman::Aligner aligner_to_base(1, 4, 6, 1, true);
    StripedSmithWaterman::Aligner harsh_aligner(1, 4, 100, 1, true);

    for (cc_v_distance_t& cc_v_distance : cc_v_distances) {

        std::shared_ptr<insertion_cluster_t> c1 = cc_v_distance.c1;
        std::shared_ptr<insertion_cluster_t> c2 = cc_v_distance.c2;

        if (c1->cluster->used || c2->cluster->used) continue;

        // remap clusters
		bool success = find_insertion_from_cluster_pair(c1, c2, contig_id, hdr, *mateseqs.get(), *matequals.get(),
			aligner, permissive_aligner, aligner_to_base, harsh_aligner, stats);
        if (success) {
            c1->cluster->used = true; c2->cluster->used = true;
        }
    }
    std::unordered_set<std::shared_ptr<insertion_cluster_t> > clusters;
    for (cc_v_distance_t& cc_v_distance : cc_v_distances) {
        cc_v_distance.c1->deallocate_reads();
        cc_v_distance.c2->deallocate_reads();
        clusters.insert(cc_v_distance.c1);
        clusters.insert(cc_v_distance.c2);
    }
}

void find_contig_insertions(int contig_id, ctpl::thread_pool& thread_pool, std::vector<std::future<void> >& futures) {

    std::string r_dc_fname = workdir + "/workspace/fwd-stable/" + std::to_string(contig_id) + ".noremap.bam";
    std::string l_dc_fname = workdir + "/workspace/rev-stable/" + std::to_string(contig_id) + ".noremap.bam";
    if (!file_exists(r_dc_fname) || !file_exists(l_dc_fname)) return;

    open_samFile_t* r_dc_file = new open_samFile_t(r_dc_fname.c_str(), true);
    open_samFile_t* l_dc_file = new open_samFile_t(l_dc_fname.c_str(), true);

    std::string clip_consensus_fname = workdir + "/workspace/consensuses/" + std::to_string(contig_id) + ".txt";
    std::vector<std::shared_ptr<consensus_t>> rc_consensuses, lc_consensuses;
    if (file_exists(clip_consensus_fname)) {
        std::ifstream clipped_fin(clip_consensus_fname);
        std::string line;
        while (std::getline(clipped_fin, line)) {
            std::shared_ptr<consensus_t> consensus = std::make_shared<consensus_t>(line);
            if (consensus->is_hsr) continue;
            if (consensus->left_clipped) {
                consensus->start += consensus->lowq_prefix;
                if (consensus->breakpoint < consensus->start) {
                    consensus->breakpoint = consensus->start;
                }
                consensus->clip_len = consensus->breakpoint - consensus->start;
                consensus->sequence = consensus->sequence.substr(consensus->lowq_prefix);
                lc_consensuses.push_back(consensus);
            } else {
                consensus->end -= consensus->lowq_suffix;
                if (consensus->breakpoint > consensus->end) {
                    consensus->breakpoint = consensus->end;
                }
                consensus->clip_len = consensus->end - consensus->breakpoint;
                consensus->sequence = consensus->sequence.substr(0, consensus->sequence.length()-consensus->lowq_suffix);
                rc_consensuses.push_back(consensus);
            }
        }
    }

	std::vector<std::shared_ptr<insertion_cluster_t>> r_clusters = cluster_reads(r_dc_file, contig_id, rc_consensuses);
	std::vector<std::shared_ptr<insertion_cluster_t>> l_clusters = cluster_reads(l_dc_file, contig_id, lc_consensuses);

	std::string clipped_fname = workdir + "/workspace/sr/" + std::to_string(contig_id) + ".bam";
	add_semi_mapped_pairs(clipped_fname, contig_id, r_clusters, l_clusters);

    std::vector<cc_v_distance_t> pq;
    std::multimap<int, std::shared_ptr<insertion_cluster_t> > l_clusters_map;
    for (std::shared_ptr<insertion_cluster_t> l_cluster : l_clusters) {
        hts_pos_t cl_start = l_cluster->start;
        l_clusters_map.insert({cl_start, l_cluster});
    }
    for (std::shared_ptr<insertion_cluster_t> r_cluster : r_clusters) {
        hts_pos_t cl_end = r_cluster->end;
        auto begin = l_clusters_map.lower_bound(cl_end - stats.max_is);
        auto end = l_clusters_map.upper_bound(cl_end + stats.max_is);

        for (auto it = begin; it != end; it++) {
            std::shared_ptr<insertion_cluster_t> l_cluster = it->second;
            pq.push_back(cc_v_distance_t(r_cluster, l_cluster, r_cluster->size()*l_cluster->size()));
        }
    }
    std::sort(pq.begin(), pq.end(), [](const cc_v_distance_t& ccd1, const cc_v_distance_t& ccd2) {
        return ccd1.distance > ccd2.distance;
    });

    std::unordered_set<std::string> retained_read_names;
    for (std::shared_ptr<insertion_cluster_t> c : r_clusters) {
        for (std::shared_ptr<bam1_t> r : c->cluster->reads) {
            std::string qname = bam_get_qname(r);
            if (is_samechr(r.get())) {
                retained_read_names.insert(qname + "_1");
                retained_read_names.insert(qname + "_2");
            } else {
                retained_read_names.insert(qname);
            }
        }
    }
    for (std::shared_ptr<insertion_cluster_t> c : l_clusters) {
        for (std::shared_ptr<bam1_t> r : c->cluster->reads) {
            std::string qname = bam_get_qname(r);
            if (is_samechr(r.get())) {
                retained_read_names.insert(qname + "_1");
                retained_read_names.insert(qname + "_2");
            } else {
                retained_read_names.insert(qname);
            }
        }
    }

    auto mateseqs = std::make_shared<std::unordered_map<std::string, std::string> >();
    auto matequals = std::make_shared<std::unordered_map<std::string, std::string> >();
    std::ifstream mateseqs_fin(workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt");
    std::string qname, seq, qual; int mapq;
    while (mateseqs_fin >> qname >> seq >> qual >> mapq) {
        if (retained_read_names.count(qname)) {
            (*mateseqs)[qname] = seq;
            (*matequals)[qname] = qual;
        }
    }
    mateseqs_fin.close();

    // create dependencies graph for pq
    std::vector<std::vector<int> > dependencies(pq.size());
    std::unordered_map<std::shared_ptr<insertion_cluster_t>, std::vector<int> > observed_pos_in_pq_for_clusters;
    for (int curr_idx = 0; curr_idx < pq.size(); curr_idx++) {
        cc_v_distance_t cc_v_distance = pq[curr_idx];
        std::shared_ptr<insertion_cluster_t> c1 = cc_v_distance.c1;
        std::shared_ptr<insertion_cluster_t> c2 = cc_v_distance.c2;
        for (int prev_idx : observed_pos_in_pq_for_clusters[c1]) {
            dependencies[curr_idx].push_back(prev_idx);
            dependencies[prev_idx].push_back(curr_idx);
        }
        for (int prev_idx : observed_pos_in_pq_for_clusters[c2]) {
            dependencies[curr_idx].push_back(prev_idx);
            dependencies[prev_idx].push_back(curr_idx);
        }
        observed_pos_in_pq_for_clusters[c1].push_back(curr_idx);
        observed_pos_in_pq_for_clusters[c2].push_back(curr_idx);
    }

    // find connected components
    std::vector<int> dfs_labels(pq.size(), -1);
    int curr_label = 0;
    for (int i = 0; i < pq.size(); i++) {
        if (dfs_labels[i] == -1) {
            dfs(i, dependencies, dfs_labels, curr_label);
            curr_label++;
        }
    }
    
    std::vector<std::vector<cc_v_distance_t> > cc_v_distances_components(curr_label);
    for (int i = 0; i < pq.size(); i++) {
        cc_v_distance_t cc_v_distance = pq[i];
        int label = dfs_labels[i];
        cc_v_distances_components[label].push_back(cc_v_distance);
    }
    pq.clear();

    bam_hdr_t* header = bam_hdr_dup(l_dc_file->header);

    int MIN_COMP_SIZE = 1000;
    std::vector<std::vector<cc_v_distance_t> > batches(1);
    for (int comp_id = 0; comp_id < cc_v_distances_components.size(); comp_id++) {
        std::vector<cc_v_distance_t>& cc_v_distances_component = cc_v_distances_components[comp_id];
        std::sort(cc_v_distances_component.begin(), cc_v_distances_component.end(), [](const cc_v_distance_t& ccd1, const cc_v_distance_t& ccd2) {
            return ccd1.distance > ccd2.distance;
        });

        std::vector<cc_v_distance_t>& curr_batch = batches.back();
        curr_batch.insert(curr_batch.end(), cc_v_distances_component.begin(), cc_v_distances_component.end());
        if (curr_batch.size() > MIN_COMP_SIZE || comp_id == cc_v_distances_components.size()-1) {
            std::future<void> future = thread_pool.push(find_insertions, contig_id, comp_id, curr_batch, header,
                mateseqs, matequals);
            futures.push_back(std::move(future));
            batches.push_back(std::vector<cc_v_distance_t>());
        }
    }

    delete l_dc_file;
    delete r_dc_file;
}

int main(int argc, char* argv[]) {

    workdir = std::string(argv[1]);
    std::string workspace = workdir + "/workspace";
    std::string reference_fname  = argv[2];
    std::string sample_name  = argv[3];

    std::string full_cmd_fname = workdir + "/cmd.txt";
    std::ifstream full_cmd_fin(full_cmd_fname);
    std::string full_cmd_str;
    std::getline(full_cmd_fin, full_cmd_str);

    config.parse(workdir + "/config.txt");
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    contigs.read_fasta_into_map(reference_fname);
    contig_map.load(workdir);

    assembly_failed_no_seq.open(workdir + "/intermediate_results/assembly_failed.no_seq.sv");
    assembly_failed_cycle_writer.open(workdir + "/intermediate_results/assembly_failed.w_cycle.sv");
    assembly_failed_too_many_reads_writer.open(workdir + "/intermediate_results/assembly_failed.too_many_reads.sv");

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void> > futures;

    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        find_contig_insertions(contig_id, thread_pool, futures);
    }
    thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }

	bcf1_t* bcf_entry = bcf_init();

	// out transurveyor insertions
	std::string out_vcf_fname = workdir + "/intermediate_results/assembled_ins.vcf.gz";
    bcf_hdr_t* out_vcf_header = generate_vcf_header(contigs, sample_name, config, full_cmd_str);
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (out_vcf_file == NULL) {
		throw std::runtime_error("Unable to open file " + out_vcf_fname + ".");
	}
	if (bcf_hdr_write(out_vcf_file, out_vcf_header) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

    std::sort(insertions.begin(), insertions.end(), [](const std::shared_ptr<sv_t>& sv1, const std::shared_ptr<sv_t>& sv2) {
        int chr1 = contig_map.get_id(sv1->chr);
        int chr2 = contig_map.get_id(sv2->chr);
        if (chr1 != chr2) return chr1 < chr2;
        return sv_output_order(sv1, sv2);
    });

	int a_id = 0, t_id = 0;
	for (const std::shared_ptr<sv_t>& insertion : insertions) {
        if (insertion->ins_seq.length() < config.min_assembled_ins_size) {
            continue;
        }

        if (insertion->source == "REFERENCE_GUIDED_ASSEMBLY") insertion->id = "T_INS_" + std::to_string(t_id++);
        else insertion->id = "A_INS_" + std::to_string(a_id++);
        sv2bcf(out_vcf_header, bcf_entry, insertion.get(), contigs.get_seq(insertion->chr));
		if (bcf_write(out_vcf_file, out_vcf_header, bcf_entry) != 0) {
			throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
		}
	}

    bcf_close(out_vcf_file);
}
