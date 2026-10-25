#include <fstream>
#include <memory>
#include <string>
#include <mutex>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "../libs/cptl_stl.h"
#include "../libs/IntervalTree.h"
#include "sam_utils.h"
#include "sw_utils.h"
#include "types.h"
#include "utils.h"
#include "vcf_utils.h"
#include "clustering_utils.h"

std::string workdir;

config_t config;
stats_t stats;
chr_seqs_map_t chr_seqs;

std::unordered_map<std::string, std::vector<std::shared_ptr<deletion_t>>> deletions_by_chr, del_sr_entries_by_chr;
std::unordered_map<std::string, std::vector<std::shared_ptr<duplication_t>>> duplications_by_chr, dup_sr_entries_by_chr;
std::unordered_map<std::string, std::vector<std::shared_ptr<sv_t>>> other_svs_entries_by_chr;
std::mutex maps_mtx;

struct cc_pair {
	std::shared_ptr<cluster_t> c1, c2;
	hts_pos_t dist;

	cc_pair(std::shared_ptr<cluster_t> c1, std::shared_ptr<cluster_t> c2) : c1(c1), c2(c2), dist(cluster_t::distance(c1, c2)) {}

	bool compatible() {
		return dist <= stats.max_is;
	}
};
bool operator < (const cc_pair& cc1, const cc_pair& cc2) { return cc1.dist > cc2.dist; }

void cluster_clusters(std::vector<std::shared_ptr<cluster_t>>& clusters, int min_cluster_size, int max_cluster_size) {
	std::vector<std::shared_ptr<bam1_t>> reads;
	cluster_clusters(clusters, reads, stats.max_is, max_cluster_size, false);

	clusters.erase(std::remove_if(clusters.begin(), clusters.end(), [min_cluster_size](std::shared_ptr<cluster_t> c) {
		return c == NULL || c->count < min_cluster_size;
	}), clusters.end());

	std::sort(clusters.begin(), clusters.end());
}

std::vector<std::shared_ptr<cluster_t>> cluster_dps_support(int contig_id, std::string contig_name, std::string dp_dir) {

	std::string dp_fname = dp_dir + "/" + std::to_string(contig_id) + ".bam";
	if (!file_exists(dp_fname)) return std::vector<std::shared_ptr<cluster_t>>();

	std::vector<std::shared_ptr<cluster_t>> clusters;
	open_samFile_t* dp_bam_file = open_samFile(dp_fname, true);
	hts_itr_t* iter = sam_itr_querys(dp_bam_file->idx, dp_bam_file->header, contig_name.c_str());
	bam1_t* read = bam_init1();
	while (sam_itr_next(dp_bam_file->file, iter, read) >= 0) {
		std::string qname = bam_get_qname(read);
		std::shared_ptr<cluster_t> cluster = std::make_shared<cluster_t>(read, config.high_confidence_mapq);
		clusters.push_back(cluster);
	}
	close_samFile(dp_bam_file);

	int min_cluster_size = std::max(3, int(stats.get_median_depth(contig_name)+5)/10);
	int max_cluster_size = (stats.get_max_depth(contig_name) * stats.max_is)/stats.read_len;
	if (!clusters.empty()) cluster_clusters(clusters, min_cluster_size, max_cluster_size);

	// set ra_furthermost_seq
	std::unordered_map<std::string, std::shared_ptr<cluster_t>> mateseqs_to_retrieve;
	for (std::shared_ptr<cluster_t> c : clusters) {
		mateseqs_to_retrieve[c->ra_furthermost_seq] = c;
	}

	std::string qname, seq;
	int64_t nm;
	std::ifstream mateseqs_fin(dp_dir + std::to_string(contig_id) + ".txt");
	mateseqs_fin.clear();
	mateseqs_fin.seekg(0);
	while (mateseqs_fin >> qname >> seq >> nm) {
		if (!mateseqs_to_retrieve.count(qname)) continue;
		mateseqs_to_retrieve[qname]->ra_furthermost_seq = seq;
	}
	mateseqs_to_retrieve.clear();

	return clusters;
}

void cluster_lp_dps(int contig_id, std::string contig_name, std::vector<std::shared_ptr<deletion_t>>& deletions) {

	std::vector<std::shared_ptr<cluster_t>> clusters = cluster_dps_support(contig_id, contig_name, workdir + "/workspace/long-pairs");
	if (clusters.empty()) return;

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	for (std::shared_ptr<cluster_t> c : clusters) {
		if (c->used) continue;

		std::shared_ptr<consensus_t> rc_consensus = std::make_shared<consensus_t>(false, 0, c->la_end, 0, c->la_furthermost_seq, 0, 0, 0, 0, 0, 0);
		std::shared_ptr<consensus_t> lc_consensus = std::make_shared<consensus_t>(false, 0, c->ra_start, 0, c->ra_furthermost_seq, 0, 0, 0, 0, 0, 0);
		std::vector<std::shared_ptr<sv_t>> svs = detect_svs(contig_name, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), rc_consensus, lc_consensus, aligner, stats.read_len/3, config.min_clip_len, 0.0);

		std::shared_ptr<deletion_t> del = NULL;
		if (!svs.empty() && svs[0]->svtype() == "DEL") { // first check if I can obtain precise coordinates for the deletion
			del = std::dynamic_pointer_cast<deletion_t>(svs[0]);
			del->original_range = std::to_string(c->la_end) + "-" + std::to_string(c->ra_start);
			del->remapped = true;
		} else if (c->la_end < c->ra_start) {
			auto left_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(c->la_start, c->la_end, c->la_end-c->la_start, 0);
			auto right_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(c->ra_start, c->ra_end, c->ra_end-c->ra_start, 0);
			del = std::make_shared<deletion_t>(contig_name, c->la_end, c->ra_start, "", nullptr, nullptr, left_anchor_aln, right_anchor_aln);
			del->imprecise = true;
		}

		if (del && -del->svlen() >= config.min_sv_size) {
			del->source = "DP";
			deletions.push_back(del);
		}
	}
}

void cluster_ow_dps(int contig_id, std::string contig_name, std::vector<std::shared_ptr<duplication_t>>& duplications) {

	std::vector<std::shared_ptr<cluster_t>> clusters = cluster_dps_support(contig_id, contig_name, workdir + "/workspace/outward-pairs");
	if (clusters.empty()) return;
	
	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	for (std::shared_ptr<cluster_t> c : clusters) {
		if (c->used) continue;

		std::shared_ptr<duplication_t> dup = nullptr;
		if (c->la_start < c->ra_end) {
			auto left_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(c->la_start, c->la_end, c->la_end-c->la_start, 0);
			auto right_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(c->ra_start, c->ra_end, c->ra_end-c->ra_start, 0);
			dup = std::make_shared<duplication_t>(contig_name, c->la_start, c->ra_end, "", nullptr, nullptr, left_anchor_aln, right_anchor_aln);
			dup->imprecise = true;
		}
		
		if (dup && dup->svlen() >= config.min_sv_size) {
			dup->source = "DP";
			duplications.push_back(dup);
		}
	}
}

void cluster_dps(int id, int contig_id, std::string contig_name, std::string bam_fname) {

	std::vector<std::shared_ptr<deletion_t>> deletions;
	std::vector<std::shared_ptr<duplication_t>> duplications;
	cluster_lp_dps(contig_id, contig_name, deletions);
	cluster_ow_dps(contig_id, contig_name, duplications);

	maps_mtx.lock();
	deletions_by_chr[contig_name] = deletions;
	duplications_by_chr[contig_name] = duplications;
	maps_mtx.unlock();
}

void merge_sr_dp_dels(int id, int contig_id, std::string contig_name) {
	maps_mtx.lock();
	std::vector<std::shared_ptr<deletion_t>>& deletions = deletions_by_chr[contig_name];
	std::vector<std::shared_ptr<deletion_t>>& sv_sr_entries = del_sr_entries_by_chr[contig_name];
	maps_mtx.unlock();

	if (deletions.empty() || sv_sr_entries.empty()) return;

	std::vector<Interval<std::shared_ptr<sv_t>>> del_intervals;
	for (const auto& sr_entry : sv_sr_entries) {
		del_intervals.push_back(Interval<std::shared_ptr<sv_t>>(sr_entry->start, sr_entry->end, sr_entry));
	}
	IntervalTree<std::shared_ptr<sv_t>> del_tree = IntervalTree<std::shared_ptr<sv_t>>(del_intervals);

	for (int i = 0; i < deletions.size(); i++) {
		std::shared_ptr<deletion_t> del = deletions[i];
		std::vector<Interval<std::shared_ptr<sv_t>>> ivs = del_tree.findContained(del->left_anchor_aln->start, del->right_anchor_aln->end);

		std::vector<sv_t*> compatible_dels;
		for (Interval<std::shared_ptr<sv_t>> iv : ivs) {
			sv_t* corr_sr_del = iv.value.get();
			if (corr_sr_del->start-del->left_anchor_aln->start <= stats.max_is && del->right_anchor_aln->end-corr_sr_del->end <= stats.max_is &&
				corr_sr_del->start >= del->start-stats.read_len/2 && corr_sr_del->end <= del->end+stats.read_len/2) {
				compatible_dels.push_back(corr_sr_del);
			}
		}
		if (compatible_dels.size() == 1) {
			deletions[i] = nullptr;
		}
	}
	deletions.erase(std::remove(deletions.begin(), deletions.end(), nullptr), deletions.end());
}

void merge_sr_dp_dups(int id, int contig_id, std::string contig_name) {
	maps_mtx.lock();
	std::vector<std::shared_ptr<duplication_t>>& duplications = duplications_by_chr[contig_name];
	std::vector<std::shared_ptr<duplication_t>>& sv_sr_entries = dup_sr_entries_by_chr[contig_name];
	maps_mtx.unlock();

	if (duplications.empty() || sv_sr_entries.empty()) return;

	std::vector<Interval<std::shared_ptr<sv_t>>> dup_intervals;
	for (const auto& sr_entry : sv_sr_entries) {
		dup_intervals.push_back(Interval<std::shared_ptr<sv_t>>(sr_entry->start, sr_entry->end, sr_entry));
	}
	IntervalTree<std::shared_ptr<sv_t>> dup_tree = IntervalTree<std::shared_ptr<sv_t>>(dup_intervals);

	for (int i = 0; i < duplications.size(); i++) {
		std::shared_ptr<duplication_t> dup = duplications[i];
		std::vector<Interval<std::shared_ptr<sv_t>>> ivs = dup_tree.findContained(dup->start-stats.max_is, dup->end+stats.max_is);

		std::vector<sv_t*> compatible_dups;
		for (Interval<std::shared_ptr<sv_t>> iv : ivs) {
			sv_t* corr_sr_dup = iv.value.get();
			if (dup->left_anchor_aln->end-corr_sr_dup->start <= stats.max_is && corr_sr_dup->end-dup->right_anchor_aln->start <= stats.max_is &&
				corr_sr_dup->start <= dup->start+stats.read_len/2 && corr_sr_dup->end >= dup->end-stats.read_len/2) {
				compatible_dups.push_back(corr_sr_dup);
			}
		}
		if (compatible_dups.size() == 1) {
			duplications[i] = NULL;
		}
	}
	duplications.erase(std::remove(duplications.begin(), duplications.end(), nullptr), duplications.end());
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];

    workdir = argv[2];
    std::string workspace = workdir + "/workspace";
    std::string reference_fname = argv[3];
    std::string sample_name = argv[4];

    contig_map_t contig_map(workdir);
    config.parse(workdir + "/config.txt");

	stats.parse(workdir + "/stats.txt", config.per_contig_stats);

	chr_seqs.read_fasta_into_map(reference_fname);

	ctpl::thread_pool thread_pool1(config.threads);
	std::vector<std::future<void> > futures;
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool1.push(cluster_dps, contig_id, contig_name, bam_fname);
		futures.push_back(std::move(future));
	}
	thread_pool1.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		try {
			futures[i].get();
		} catch (char const* s) {
			std::cout << s << std::endl;
		}
	}
	futures.clear();

	std::string sr_vcf_fname = workdir + "/intermediate_results/sr.dedup.vcf.gz";
	htsFile* sr_vcf_file = bcf_open(sr_vcf_fname.c_str(), "r");
	bcf_hdr_t* hdr = bcf_hdr_read(sr_vcf_file);

	bcf1_t* bcf_entry = bcf_init();
	while (bcf_read(sr_vcf_file, hdr, bcf_entry) == 0) {
		std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, bcf_entry);
		if (sv->svtype() == "DEL" && sv->is_pass()) {
			del_sr_entries_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<deletion_t>(sv));
		} else if (sv->svtype() == "DUP" && sv->is_pass()) {
			dup_sr_entries_by_chr[sv->chr].push_back(std::dynamic_pointer_cast<duplication_t>(sv));
		} else {
			other_svs_entries_by_chr[sv->chr].push_back(sv);
		}
	}
	bcf_close(sr_vcf_file);

	// merge SR and DISC
	ctpl::thread_pool thread_pool2(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool2.push(merge_sr_dp_dels, contig_id, contig_name);
		futures.push_back(std::move(future));
		future = thread_pool2.push(merge_sr_dp_dups, contig_id, contig_name);
		futures.push_back(std::move(future));
	}
	thread_pool2.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		try {
			futures[i].get();
		} catch (char const* s) {
			std::cout << s << std::endl;
		}
	}
	futures.clear();

	std::string dp_vcf_fname = workdir + "/intermediate_results/dp.vcf.gz";
	htsFile* dp_vcf_file = bcf_open(dp_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(dp_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + dp_vcf_fname + ".");
	}

	// read dp inversions dp_invs.vcf.gz
	std::string dp_inv_vcf_fname = workdir + "/intermediate_results/dp_invs.vcf.gz";
	htsFile* dp_inv_vcf_file = bcf_open(dp_inv_vcf_fname.c_str(), "r");
	bcf_hdr_t* hdr_inv = bcf_hdr_read(dp_inv_vcf_file);

	std::unordered_map<std::string, std::vector<std::shared_ptr<sv_t>>> dp_svs_by_chr;
	while (bcf_read(dp_inv_vcf_file, hdr_inv, bcf_entry) == 0) {
		std::shared_ptr<sv_t> sv = bcf_to_sv(hdr_inv, bcf_entry);
		dp_svs_by_chr[sv->chr].push_back(sv);
	}

	int del_id = 0, dup_id = 0;
    for (std::string& contig_name : chr_seqs.ordered_contigs) {
		auto& deletions = deletions_by_chr[contig_name];
		for (const auto& del : deletions) del->id = "DEL_DP_" + std::to_string(del_id++);
		auto& duplications = duplications_by_chr[contig_name];
		for (const auto& dup : duplications) dup->id = "DUP_DP_" + std::to_string(dup_id++);

		dp_svs_by_chr[contig_name].insert(dp_svs_by_chr[contig_name].end(), deletions.begin(), deletions.end());
		dp_svs_by_chr[contig_name].insert(dp_svs_by_chr[contig_name].end(), duplications.begin(), duplications.end());
		std::sort(dp_svs_by_chr[contig_name].begin(), dp_svs_by_chr[contig_name].end(), [](const std::shared_ptr<sv_t>& sv1, const std::shared_ptr<sv_t>& sv2) {
			return sv1->start < sv2->start;
		});
		for (const auto& dp_sv : dp_svs_by_chr[contig_name]) {
			sv2bcf(hdr, bcf_entry, dp_sv.get(), chr_seqs.get_seq(contig_name));
			if (bcf_write(dp_vcf_file, hdr, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + dp_vcf_fname + ".");
			}
		}
    }

	bcf_close(dp_vcf_file);

	tbx_index_build(dp_vcf_fname.c_str(), 0, &tbx_conf_vcf);

	std::string merged_vcf_fname = workdir + "/intermediate_results/sr_dp.vcf.gz";
	htsFile* merged_vcf_file = bcf_open(merged_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(merged_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + merged_vcf_fname + ".");
	}

	for (std::string& contig_name : chr_seqs.ordered_contigs) {
		std::vector<std::shared_ptr<sv_t>> all_entries;
		all_entries.insert(all_entries.end(), dp_svs_by_chr[contig_name].begin(), dp_svs_by_chr[contig_name].end());
		all_entries.insert(all_entries.end(), del_sr_entries_by_chr[contig_name].begin(), del_sr_entries_by_chr[contig_name].end());
		all_entries.insert(all_entries.end(), dup_sr_entries_by_chr[contig_name].begin(), dup_sr_entries_by_chr[contig_name].end());
		all_entries.insert(all_entries.end(), other_svs_entries_by_chr[contig_name].begin(), other_svs_entries_by_chr[contig_name].end());
		std::sort(all_entries.begin(), all_entries.end(), [](const std::shared_ptr<sv_t>& sv1, const std::shared_ptr<sv_t>& sv2) {
			return sv1->start < sv2->start;
		});
		for (const auto& sv : all_entries) {
			sv2bcf(hdr, bcf_entry, sv.get(), chr_seqs.get_seq(contig_name));
			if (bcf_write(merged_vcf_file, hdr, bcf_entry) != 0) {
				throw std::runtime_error("Failed to write to " + merged_vcf_fname + ".");
			}
		}
	}

	bcf_close(merged_vcf_file);

	tbx_index_build(merged_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
