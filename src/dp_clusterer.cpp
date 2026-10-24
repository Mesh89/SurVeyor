#include <fstream>
#include <memory>
#include <string>
#include <mutex>
#include <map>
#include <set>
#include <tuple>
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
#include "coarse_coverage.h"

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
}

std::vector<std::shared_ptr<cluster_t>> cluster_dps_support(int contig_id, std::string contig_name, std::string dp_dir) {

	std::string dp_fname = dp_dir + "/" + std::to_string(contig_id) + ".bam";
	if (!file_exists(dp_fname)) return std::vector<std::shared_ptr<cluster_t>>();

	std::vector<std::shared_ptr<cluster_t>> clusters;
	open_samFile_t* dp_bam_file = new open_samFile_t(dp_fname, true);
	hts_itr_t* iter = sam_itr_querys(dp_bam_file->idx, dp_bam_file->header, contig_name.c_str());
	bam1_t* read = bam_init1();
	while (sam_itr_next(dp_bam_file->file, iter, read) >= 0) {
		std::string qname = bam_get_qname(read);
		std::shared_ptr<cluster_t> cluster = std::make_shared<cluster_t>(read, config.high_confidence_mapq);
		clusters.push_back(cluster);
	}
	hts_itr_destroy(iter);
	bam_destroy1(read);
	delete dp_bam_file;

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
	std::ifstream mateseqs_fin(dp_dir + "/" + std::to_string(contig_id) + ".txt");
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

		std::shared_ptr<consensus_t> rc_consensus = std::make_shared<consensus_t>(false, 0, c->la_end, 0, c->la_furthermost_seq, std::string(c->la_furthermost_seq.length(), '!'), 0, 0, 0, 0, 0, 0);
		std::shared_ptr<consensus_t> lc_consensus = std::make_shared<consensus_t>(false, 0, c->ra_start, 0, c->ra_furthermost_seq, std::string(c->ra_furthermost_seq.length(), '!'), 0, 0, 0, 0, 0, 0);
		std::vector<std::shared_ptr<sv_t>> svs = detect_svs(contig_name, chr_seqs.get_seq(contig_name), chr_seqs.get_len(contig_name), rc_consensus, lc_consensus, 
			aligner, stats.read_len/3, stats, config);

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

		if (del && del->svsize() >= config.min_sv_size) {
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
			dup->cn_unresolved = true;
		}

		if (dup && dup->svsize() >= config.min_sv_size) {
			dup->source = "DP";
			duplications.push_back(dup);
		}
	}
}

void cluster_lp_task(int id, int contig_id, std::string contig_name) {

	std::vector<std::shared_ptr<deletion_t>> deletions;
	cluster_lp_dps(contig_id, contig_name, deletions);

	std::lock_guard<std::mutex> lock(maps_mtx);
	deletions_by_chr[contig_name] = deletions;
}

void cluster_ow_task(int id, int contig_id, std::string contig_name) {

	std::vector<std::shared_ptr<duplication_t>> duplications;
	cluster_ow_dps(contig_id, contig_name, duplications);

	std::lock_guard<std::mutex> lock(maps_mtx);
	duplications_by_chr[contig_name] = duplications;
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

void refine_dp_dels(std::vector<std::shared_ptr<deletion_t>>& deletions, std::string contig_name, open_samFile_t& bam_file, const coarse_coverage_track_t& coverage) {
	char* contig_seq = chr_seqs.get_seq(contig_name);
	hts_pos_t contig_len = chr_seqs.get_len(contig_name);
	uint64_t max_depth = stats.get_max_depth(contig_name);
	int tid = bam_name2id(bam_file.header, contig_name.c_str());
	if (tid < 0) return;
	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	StripedSmithWaterman::Filter filter(false, false, 0, 0);
	for (auto& del : deletions) {
		if (!del->imprecise) continue;
		hts_pos_t la_start = del->left_anchor_aln->start, la_end = del->start;
		hts_pos_t ra_start = del->end, ra_end = del->right_anchor_aln->end;
		hts_pos_t left_bp_start = std::max(hts_pos_t(0), la_end-stats.read_len/2);
		hts_pos_t left_bp_end = std::min(contig_len-1, la_start+stats.max_is);
		hts_pos_t right_bp_start = std::max(hts_pos_t(0), ra_end-stats.max_is);
		hts_pos_t right_bp_end = std::min(contig_len-1, ra_start+stats.read_len/2);
		if (left_bp_start > left_bp_end || right_bp_start > right_bp_end) continue;

		hts_pos_t lh_start = std::max(hts_pos_t(0), left_bp_start-stats.read_len);
		hts_pos_t lh_end = std::min(contig_len, left_bp_end+stats.read_len+1);
		hts_pos_t rh_start = std::max(hts_pos_t(0), right_bp_start-stats.read_len);
		hts_pos_t rh_end = std::min(contig_len, right_bp_end+stats.read_len);
		// A 1 kb window above 4×max_depth indicates a collapsed repeat.
		if (coverage.exceeds(lh_start, lh_end, 4*max_depth) || coverage.exceeds(rh_start, rh_end, 4*max_depth)) continue;

		// Brace-quote contig names so colons are not parsed as coordinates.
		std::string left_region = "{" + contig_name + "}:" + std::to_string(left_bp_start+1) + "-" + std::to_string(left_bp_end+1);
		std::string right_region = "{" + contig_name + "}:" + std::to_string(right_bp_start+2) + "-" + std::to_string(right_bp_end+2);
		char* regions[2] = {&left_region[0], &right_region[0]};
		hts_itr_t* iter = sam_itr_regarray(bam_file.idx, bam_file.header, regions, 2);
		if (!iter) throw std::runtime_error("Failed to create refinement iterator for " + contig_name);
		std::vector<std::shared_ptr<bam1_t>> left_reads, right_reads;
		bam1_t* read = bam_init1();
		int ret;
		while ((ret = sam_itr_next(bam_file.file, iter, read)) >= 0) {
			if (!is_primary(read) || is_unmapped(read) || read->core.qual < config.min_stable_mapq) continue;
			int left_clip = get_left_clip_size(read), right_clip = get_right_clip_size(read);
			// At the left breakpoint, require a right clip and no substantial left clip.
			if (left_clip < config.min_clip_len && right_clip >= config.min_clip_len) {
				hts_pos_t bp = bam_endpos(read)-1;
				if (bp < left_bp_start || bp > left_bp_end) continue;
				left_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
			}
			// At the right breakpoint, require a left clip and no substantial right clip.
			if (left_clip >= config.min_clip_len && right_clip < config.min_clip_len) {
				hts_pos_t bp = read->core.pos-1;
				if (bp < right_bp_start || bp > right_bp_end) continue;
				right_reads.push_back(std::shared_ptr<bam1_t>(bam_dup1(read), bam_destroy1));
			}
		}
		bam_destroy1(read);
		hts_itr_destroy(iter);
		if (ret < -1) throw std::runtime_error("Failed to read refinement BAM for " + contig_name);

		struct junction_read_t {
			std::string name;
			std::shared_ptr<deletion_t> del;
		};
		using junction_reads_t = std::vector<junction_read_t>;
		std::map<std::tuple<hts_pos_t, hts_pos_t, std::string>, junction_reads_t> junctions;
		for (const auto* reads : {&left_reads, &right_reads}) {
			bool left_side = reads == &left_reads;
			for (const auto& read : *reads) {
				std::string seq = get_sequence(read.get());
				std::string qual = get_qual_ascii(read.get());

				// find best alignment score for the split read against the left and right haplotype sequences
				StripedSmithWaterman::Alignment lh_aln, rh_aln;
				aligner.Align(seq.c_str(), contig_seq+lh_start, lh_end-lh_start, filter, &lh_aln, 0);
				aligner.Align(seq.c_str(), contig_seq+rh_start, rh_end-rh_start, filter, &rh_aln, 0);
				int unsplit_score = std::max(lh_aln.sw_score, rh_aln.sw_score);

				// find the best variant supported by the split read
				hts_pair_pos_t split_lh = {lh_start, lh_end}, split_rh = {rh_start, rh_end};
				if (left_side) split_lh = {read->core.pos, bam_endpos(read.get())};
				else split_rh = {read->core.pos, bam_endpos(read.get())};
				std::vector<std::shared_ptr<sv_t>> svs = detect_svs_from_junction(contig_name, contig_seq, seq, qual, split_lh.beg, split_lh.end, split_rh.beg, split_rh.end, 
					aligner, 0, 0, stats, config, contig_len, stats.max_is, true);
				if (svs.empty() || svs[0]->svtype() != "DEL") continue;
				if (svs[0]->left_anchor_aln->best_score < config.min_clip_len || svs[0]->right_anchor_aln->best_score < config.min_clip_len) continue;
				int split_score = svs[0]->left_anchor_aln->best_score + svs[0]->right_anchor_aln->best_score;
				if (split_score - unsplit_score < config.min_clip_len) continue;

				// left-align
				auto candidate = std::make_shared<deletion_t>(*std::static_pointer_cast<deletion_t>(svs[0]));
				// Temporarily writes AUX SNP bases into contig_seq; safe because refinement runs one task per contig after all other pools finish.
				haplotype_normalization::left_align(candidate, contig_seq);
				if (candidate->start < left_bp_start || candidate->start > left_bp_end) continue;
				if (candidate->end < right_bp_start || candidate->end > right_bp_end) continue;
				if (candidate->svsize() < config.min_sv_size) continue;

				auto junction_key = std::make_tuple(candidate->start, candidate->end, candidate->ins_seq);
				auto unshifted_del = std::static_pointer_cast<deletion_t>(svs[0]);
				junctions[junction_key].push_back({bam_get_qname(read.get()), unshifted_del});
			}
		}

		// A junction is resolvable if exactly one junction supports it
		int supported_junctions = 0;
		const junction_reads_t* supported_junction = nullptr;
		for (const auto& junction : junctions) {
			std::set<std::string> read_names;
			for (const auto& support : junction.second) read_names.insert(support.name);
			if (read_names.size() >= 2) {
				supported_junctions++;
				supported_junction = &junction.second;
			}
		}
		if (supported_junctions != 1) continue;

		const auto& representative = supported_junction->front().del;
		auto refined = std::make_shared<deletion_t>(*del);
		refined->start = representative->start;
		refined->end = representative->end;
		refined->ins_seq = representative->ins_seq;
		refined->left_anchor_aln = representative->left_anchor_aln;
		refined->right_anchor_aln = representative->right_anchor_aln;
		refined->mh_len = representative->mh_len;
		refined->junction_remap_ref_beg = representative->junction_remap_ref_beg;
		refined->junction_remap_ref_end = representative->junction_remap_ref_end;
		refined->aux_snps.clear();
		refined->aux_indels.clear();
		if (refined->original_range.empty()) refined->original_range = std::to_string(del->start) + "-" + std::to_string(del->end);
		refined->imprecise = false;
		refined->remapped = true;
		del = refined;
	}
}

void refine_dp_dels_task(int id, int contig_id, std::string contig_name) {
	auto it = deletions_by_chr.find(contig_name);
	if (it == deletions_by_chr.end()) return;
	auto& deletions = it->second;
	if (std::none_of(deletions.begin(), deletions.end(), [](const std::shared_ptr<deletion_t>& del) { return del->imprecise; })) return;

	std::string sr_bam_fname = workdir + "/workspace/sr/" + std::to_string(contig_id) + ".bam";
	if (!file_exists(sr_bam_fname)) return;
	open_samFile_t bam_file(sr_bam_fname);
	std::string coverage_fname = workdir + "/workspace/coverage/" + std::to_string(contig_id) + ".bin";
	coarse_coverage_track_t coverage = read_coarse_coverage(coverage_fname, chr_seqs.get_len(contig_name), config.coverage_window_size);
	if (coverage.maxima.empty()) throw std::runtime_error("Missing coarse coverage track " + coverage_fname);
	refine_dp_dels(deletions, contig_name, bam_file, coverage);
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

	chr_seqs.read_fasta_into_map(reference_fname, true, config.threads);

	ctpl::thread_pool thread_pool1(config.threads);
	std::vector<std::future<void> > futures;
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool1.push(cluster_lp_task, contig_id, contig_name);
		futures.push_back(std::move(future));
		future = thread_pool1.push(cluster_ow_task, contig_id, contig_name);
		futures.push_back(std::move(future));
	}
	thread_pool1.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

	std::string sr_vcf_fname = workdir + "/intermediate_results/sr.norm.dedup.vcf.gz";
	htsFile* sr_vcf_file = bcf_open(sr_vcf_fname.c_str(), "r");
	bcf_hdr_t* hdr = bcf_hdr_read(sr_vcf_file);

	bcf1_t* bcf_entry = bcf_init();
	while (bcf_read(sr_vcf_file, hdr, bcf_entry) == 0) {
		std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, bcf_entry);
		if (sv == nullptr) {
			throw std::runtime_error("Unexpected unsupported variant in internal VCF " + sr_vcf_fname + ": " + std::string(bcf_entry->d.id ? bcf_entry->d.id : "<no-id>"));
		}
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
		futures[i].get();
	}
	futures.clear();

	ctpl::thread_pool thread_pool3(config.threads);
	for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		std::future<void> future = thread_pool3.push(refine_dp_dels_task, contig_id, contig_name);
		futures.push_back(std::move(future));
	}
	thread_pool3.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		futures[i].get();
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
		if (sv == nullptr) {
			throw std::runtime_error("Unexpected unsupported variant in internal VCF " + dp_inv_vcf_fname + ": " + std::string(bcf_entry->d.id ? bcf_entry->d.id : "<no-id>"));
		}
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
		std::sort(dp_svs_by_chr[contig_name].begin(), dp_svs_by_chr[contig_name].end(), sv_output_order);
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
		std::sort(all_entries.begin(), all_entries.end(), sv_output_order);
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
