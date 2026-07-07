#include <cstddef>
#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <stdexcept>
#include <mutex>
#include <atomic>
#include <limits>

#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/kseq.h"
#include "../libs/cptl_stl.h"
#include "../libs/cxxopts.h"
#include "common.h"
#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"

std::mutex mtx;

const std::string VERSION = "0.12";

std::unordered_map<std::string, int> sample2id;
std::vector<std::string> sample_names;
chr_seqs_map_t chr_seqs;
std::mutex samples_mtx;

bool keep_all_called = false;

struct sv_w_samplename_t { // minimal SV representation for memory efficiency
	std::string id, chr, svtype;
	hts_pos_t start, end, svlen;
	std::string ins_seq, ins_seq_inferred;
	bool imprecise = false, incomplete_ins_seq = false;
    std::string sample;
	float epr;
	std::vector<int> gt;

    sv_w_samplename_t() {}
    sv_w_samplename_t(std::shared_ptr<sv_t> sv, const std::string& sample) : sample(sample), chr(sv->chr), 
		svtype(sv->svtype()), svlen(sv->svlen()), ins_seq(sv->ins_seq), ins_seq_inferred(sv->inferred_ins_seq),
		start(sv->start), end(sv->end), imprecise(sv->imprecise), incomplete_ins_seq(sv->incomplete_ins_seq()), epr(sv->sample_info.epr) {
		gt = sv->sample_info.gt;
		if (incomplete_ins_seq) svlen = 0;
	}

	std::string unique_key() {
        return chr + ":" + std::to_string(start) + ":" + std::to_string(end) + ":" + svtype + ":" + ins_seq;
    }

	int ploidy() {
		return gt.size();
	}
};

int distance(sv_w_samplename_t& sv1, sv_w_samplename_t& sv2) {
    if (sv1.chr != sv2.chr) return INT32_MAX;
	if ((sv1.svtype == "DUP" && sv2.svtype == "INS") || (sv1.svtype == "INS" && sv2.svtype == "DUP")) {
		return std::min(abs(sv1.start-sv2.start), abs(sv1.end-sv2.end));
	}
    return std::max(abs(sv1.start-sv2.start), abs(sv1.end-sv2.end));
}

double overlap(sv_w_samplename_t& sv1, sv_w_samplename_t& sv2) {
	if (sv1.svtype == "INS" && sv2.svtype == "INS") return 1.0;
	if (sv1.end == sv1.start || sv2.end == sv2.start) return 1.0;
	return overlap(sv1.start, sv1.end, sv2.start, sv2.end)/double(std::min(sv1.end-sv1.start, sv2.end-sv2.start));
}

double len_ratio(hts_pos_t svlen1, hts_pos_t svlen2) {
	svlen1 = abs(svlen1);
	svlen2 = abs(svlen2);
	if (svlen1 == 0 && svlen2 == 0) return 1.0;
	else return double(std::min(svlen1, svlen2)) / std::max(svlen1, svlen2);
}

double len_ratio(sv_w_samplename_t& sv1, sv_w_samplename_t& sv2) {
	if (sv1.svtype == "DEL" && sv2.svtype == "DEL") {
		return len_ratio(sv1.svlen, sv2.svlen);
	} else if (sv1.svtype == "DUP" && sv2.svtype == "DUP") return 1.0;
	else if (sv1.svtype == "INS" && sv2.svtype == "INS") {
		if (sv1.incomplete_ins_seq || sv2.incomplete_ins_seq) {
			return 1.0;
		} else {
			return std::min(double(sv1.svlen), double(sv2.svlen)) / std::max(double(sv1.svlen), double(sv2.svlen));
		}
	} else if ((sv1.svtype == "INS" && sv2.svtype == "DUP") || (sv1.svtype == "DUP" && sv2.svtype == "INS")) return 1.0;
	else if (sv1.svtype == "INV" && sv2.svtype == "INV") {
		double sv1_len = sv1.end - sv1.start;
		double sv2_len = sv2.end - sv2.start;
		if (sv1_len == 0 && sv2_len == 0) return 1.0;
		else return std::min(sv1_len, sv2_len) / std::max(sv1_len, sv2_len);
	} else {
		return 0.0; // not compatible
	}
}
	

struct idx_size_t {
    int idx, size;

    idx_size_t() : idx(-1), size(0) {}
    idx_size_t(int idx, int size) : idx(idx), size(size) {}
};
bool operator < (const idx_size_t& s1, const idx_size_t& s2) {
    return s1.size < s2.size;
};
bool operator > (const idx_size_t& s1, const idx_size_t& s2) {
    return s1.size > s2.size;
};

std::unordered_map<std::string, std::vector<sv_w_samplename_t> > svs_by_chr;
std::unordered_map<std::string, std::vector<std::shared_ptr<bcf1_t>> > clustered_svs_by_chr;
int max_prec_dist, max_imprec_dist, max_dist;
double min_prec_frac_overlap, min_imprec_frac_overlap;
int max_prec_len_diff, max_imprec_len_diff;
double min_prec_len_ratio, min_imprec_len_ratio, min_ins_seq_sim;
bool overlap_for_ins;

std::unordered_set<std::string> called_by;
std::mutex called_by_mtx;
std::atomic<int> ignored_unsupported_svs(0);

int aligned_matching_bases(StripedSmithWaterman::Alignment& alignment) {
	int matches = 0;
	for (uint32_t cigar_op : alignment.cigar) {
		if (cigar_int_to_op(cigar_op) == '=') {
			matches += cigar_int_to_len(cigar_op);
		}
	}
	return matches;
}

double ins_seq_similarity(const std::string& seq1, const std::string& seq2) {
	size_t total_len = seq1.length() + seq2.length();
	if (total_len == 0) return 1.0;

	const std::string& query = seq1.length() < seq2.length() ? seq1 : seq2;
	const std::string& ref = seq1.length() < seq2.length() ? seq2 : seq1;
	static thread_local StripedSmithWaterman::Aligner aligner(1, 1, 1, 1, false);
	StripedSmithWaterman::Filter filter(true, true, 0, std::numeric_limits<uint16_t>::max());
	StripedSmithWaterman::Alignment alignment;
	if (!aligner.Align(query.c_str(), ref.c_str(), ref.length(), filter, &alignment, 0)) return 0.0;
	return 2.0 * aligned_matching_bases(alignment) / total_len;
}

bool ins_seq_similar_enough(sv_w_samplename_t& sv1, sv_w_samplename_t& sv2) {
	if (min_ins_seq_sim == 0.0 || sv1.svtype != "INS" || sv2.svtype != "INS") return true;
	if (sv1.incomplete_ins_seq || sv2.incomplete_ins_seq || sv1.ins_seq.empty() || sv2.ins_seq.empty()) return true;
	if (sv1.ins_seq.length() > UINT16_MAX || sv2.ins_seq.length() > UINT16_MAX) return true;
	return ins_seq_similarity(sv1.ins_seq, sv2.ins_seq) >= min_ins_seq_sim;
}

bool is_compatible(sv_w_samplename_t& sv1, sv_w_samplename_t& sv2) {
	if (sv1.svtype != sv2.svtype) return false;

	bool distance_ok, overlap_ok, len_diff_ok, len_ratio_ok, ins_seq_sim_ok;
	
	if (sv1.imprecise || sv2.imprecise) {
		distance_ok = distance(sv1, sv2) <= max_imprec_dist;
		overlap_ok = overlap(sv1, sv2) >= min_imprec_frac_overlap;
		len_diff_ok = abs(sv1.svlen-sv2.svlen) <= max_imprec_len_diff;
		len_ratio_ok = len_ratio(sv1, sv2) >= min_imprec_len_ratio;
	} else {
		distance_ok = distance(sv1, sv2) <= max_prec_dist;
		overlap_ok = overlap(sv1, sv2) >= min_prec_frac_overlap;
		len_diff_ok = abs(sv1.svlen-sv2.svlen) <= max_prec_len_diff;
		len_ratio_ok = len_ratio(sv1, sv2) >= min_prec_len_ratio;
	}

	ins_seq_sim_ok = ins_seq_similar_enough(sv1, sv2);

	return distance_ok && overlap_ok && len_diff_ok && len_ratio_ok && ins_seq_sim_ok;
}

bool can_join_clique(std::vector<int>& neighbor_clique, std::vector<sv_w_samplename_t>& svs, int curr_idx) {
	for (int clique_elem_idx : neighbor_clique) {
		if (!is_compatible(svs[curr_idx], svs[clique_elem_idx])) {
			return false;
		}
	}
	return true;
}

std::vector<std::vector<int>> compute_minimal_clique_cover(int start, int end, std::vector<sv_w_samplename_t>& svs,
		std::vector<std::vector<int>>& compatibility_list) {
    std::vector<idx_size_t> idx_by_neighborhood_size;
    for (int i = start; i < end; i++) {
        idx_by_neighborhood_size.emplace_back(i, compatibility_list[i].size());
    }
    std::sort(idx_by_neighborhood_size.begin(), idx_by_neighborhood_size.end(), std::greater<idx_size_t>());

    std::vector<std::vector<int>> cliques;
    int cliques_counter = 0;
    std::vector<int> cliques_idx(end-start, -1);
    for (idx_size_t& i : idx_by_neighborhood_size) {
    	int curr_idx = i.idx;

        // find cliqued neighbors
        std::vector<idx_size_t> cliqued_neighbors;
        for (int adj_idx : compatibility_list[curr_idx]) {
            int clique_idx = cliques_idx[adj_idx-start];
            if (clique_idx != -1 && is_compatible(svs[curr_idx], svs[adj_idx])) {
                cliqued_neighbors.emplace_back(adj_idx, cliques[clique_idx].size());
            }
        }

        // check if I can join any neighboring cliques
        std::sort(cliqued_neighbors.begin(), cliqued_neighbors.end(), std::greater<idx_size_t>());
        bool used = false;
        for (idx_size_t& neighbor_idx : cliqued_neighbors) {
            int neighboor_clique_idx = cliques_idx[neighbor_idx.idx-start];
            std::vector<int>& neighbor_clique = cliques[neighboor_clique_idx];
            if (can_join_clique(neighbor_clique, svs, curr_idx)) {
                neighbor_clique.push_back(curr_idx);
				// note that we do not need to update cliques_idx here. Only one node in the clique (the one which started it)
				// will have its cliques_idx set. This way, we do not test can_join_clique once per node in the clique
                used = true;
                break;
            }
        }

        if (!used) {
            cliques.push_back(std::vector<int>());
            cliques[cliques_counter].push_back(curr_idx);
            cliques_idx[curr_idx-start] = cliques_counter;
            cliques_counter++;
        }
    }

    std::vector<idx_size_t> clique_idxs_by_asc_size;
    for (int i = 0; i < cliques.size(); i++) {
    	clique_idxs_by_asc_size.emplace_back(i, cliques[i].size());
    }
    std::sort(clique_idxs_by_asc_size.begin(), clique_idxs_by_asc_size.end(), std::less<idx_size_t>());

    for (int i = 0; i < clique_idxs_by_asc_size.size(); i++) {
    	int donor_clique_idx = clique_idxs_by_asc_size[i].idx;
    	for (int j = 0; j < cliques[donor_clique_idx].size(); j++) {
    		int sv_idx = cliques[donor_clique_idx][j];
			// check if it can join a larger clique
    		for (int k = clique_idxs_by_asc_size.size()-1; k > i; k--) {
    			int recipient_clique_idx = clique_idxs_by_asc_size[k].idx;
				if (can_join_clique(cliques[recipient_clique_idx], svs, sv_idx)) {
					cliques[recipient_clique_idx].push_back(sv_idx);
					cliques[donor_clique_idx][j] = -1;
					break;
				}
    		}
    	}
    	std::vector<int>& clique = cliques[donor_clique_idx];
    	clique.erase(std::remove(clique.begin(), clique.end(), -1), clique.end());
    }

    return cliques;
}

sv_w_samplename_t choose_sv(std::vector<sv_w_samplename_t>& svs) {
    std::unordered_map<std::string, int> counts;
	std::unordered_map<std::string, float> eprs;
    for (sv_w_samplename_t& sv : svs) {
		if (!sv.imprecise && !sv.incomplete_ins_seq) {
			counts[sv.unique_key()]++;
        }
		eprs[sv.unique_key()] = sv.epr;
    }
    if (counts.empty()) { // if no precise accepts imprecise
        for (sv_w_samplename_t& sv : svs) {
            counts[sv.unique_key()]++;
        }
    }

    int max_freq = 0;
    std::string str;
    for (auto& e : counts) {
        if (e.second > max_freq || (e.second == max_freq && eprs[e.first] > eprs[str])) {
			// mtx.lock();
			// if (e.second == max_freq) std::cout << "Choosen " << e.first << " instead of " << str << " because of EPR " << eprs[e.first] << " vs " << eprs[str] << std::endl;
			// mtx.unlock();
            max_freq = e.second;
            str = e.first;
        }
    }

    for (sv_w_samplename_t& sv : svs) {
        if (sv.unique_key() == str) return sv;
    }
    std::cerr << "SHOULD NOT BE HERE." << std::endl;
    exit(1);
}

std::atomic<int> glob_cluster_id(0);
void print_cliques(std::vector<std::vector<int>>& cliques, std::vector<sv_w_samplename_t>& svs, bcf_hdr_t* out_hdr,
		chr_seqs_map_t& chr_seqs, std::vector<std::shared_ptr<bcf1_t>>& append_svs) {
	int n_samples = bcf_hdr_nsamples(out_hdr);
	const char** coos = new const char*[n_samples];
	const char** sizes = new const char*[n_samples]; // could be int, but I did not find an easy way to set a variable number of ints in the format
	const char** prec = new const char*[n_samples];
	const char** incomplete_ass = new const char*[n_samples];

	std::vector<std::pair<sv_w_samplename_t, int> > sv_order;
	for (int i = 0; i < cliques.size(); i++) {
		auto& clique = cliques[i];
		if (clique.empty()) continue;

		std::vector<sv_w_samplename_t> clique_svs;
		for (int idx : clique) {
			clique_svs.push_back(svs[idx]);
		}
		sv_w_samplename_t chosen_sv = choose_sv(clique_svs);
		sv_order.push_back({chosen_sv, i});
	}
	std::sort(sv_order.begin(), sv_order.end(),
			[out_hdr](const std::pair<sv_w_samplename_t, int>& sv1, const std::pair<sv_w_samplename_t, int>& sv2) {
		int cid1 = bcf_hdr_name2id(out_hdr, sv1.first.chr.c_str());
		int cid2 = bcf_hdr_name2id(out_hdr, sv2.first.chr.c_str());
		return std::tie(cid1, sv1.first.start) < std::tie(cid2, sv2.first.start);
	});

	std::vector<bcf1_t*> vcf_svs;

	bcf1_t* vcf_sv = bcf_init();
	for (auto& clique_sv_idx : sv_order) {
		int clique_idx = clique_sv_idx.second;
		std::vector<int>& clique = cliques[clique_idx];

		std::set<std::string> unique_samples;
		std::vector<sv_w_samplename_t> clique_svs;
		for (int idx : clique) {
			unique_samples.insert(svs[idx].sample);
			clique_svs.push_back(svs[idx]);
		}

		sv_w_samplename_t& chosen_sv = clique_sv_idx.first;

		bcf_clear(vcf_sv);

		int cluster_id = glob_cluster_id++;

		// set basic info
		vcf_sv->rid = bcf_hdr_name2id(out_hdr, chosen_sv.chr.c_str());
		vcf_sv->pos = chosen_sv.start;
		std::string id = "CLUSTER_" + std::to_string(cluster_id);
		bcf_update_id(out_hdr, vcf_sv, id.c_str());
		char* chr_seq = chr_seqs.get_seq(chosen_sv.chr);
		if (!chosen_sv.incomplete_ins_seq) {
			update_alleles_info(out_hdr, vcf_sv, chr_seq, chosen_sv.start, chosen_sv.end,
				chosen_sv.svtype, chosen_sv.svlen, chosen_sv.ins_seq);
		} else {
			std::string alleles = std::string(1, chr_seq[chosen_sv.start]) + ",<" + chosen_sv.svtype + ">";
			bcf_update_alleles_str(out_hdr, vcf_sv, alleles.c_str());
		}

		// set INFO
		int int_conv = chosen_sv.end+1;
		bcf_update_info_int32(out_hdr, vcf_sv, "END", &int_conv, 1);
		bcf_update_info_string(out_hdr, vcf_sv, "SVTYPE", chosen_sv.svtype.c_str());
		int_conv = unique_samples.size();
		bcf_update_info_int32(out_hdr, vcf_sv, "N_SAMPLES", &int_conv, 1);
		int_conv = clique.size();
		bcf_update_info_int32(out_hdr, vcf_sv, "N_SVS", &int_conv, 1);
		bcf_update_info_flag(out_hdr, vcf_sv, "IMPRECISE", "", chosen_sv.imprecise);
		bcf_update_info_flag(out_hdr, vcf_sv, "INCOMPLETE_ASSEMBLY", "", chosen_sv.incomplete_ins_seq);

		if (!chosen_sv.ins_seq.empty()) {
			bcf_update_info_string(out_hdr, vcf_sv, "SVINSSEQ", chosen_sv.ins_seq.c_str());
		}
		if (!chosen_sv.ins_seq_inferred.empty()) {
			bcf_update_info_string(out_hdr, vcf_sv, "SVINSSEQ_INFERRED", chosen_sv.ins_seq_inferred.c_str());
		}
		if (!chosen_sv.incomplete_ins_seq && chosen_sv.svtype != "DUP") {
			bcf_update_info_int32(out_hdr, vcf_sv, "SVLEN", &chosen_sv.svlen, 1);
		} else if (chosen_sv.svtype == "DUP") {
			bcf_update_info_int32(out_hdr, vcf_sv, "SVLEN", NULL, 0);
		}

		int sv_ploidy = 0;
		for (sv_w_samplename_t& sv : clique_svs) {
			if (sv_ploidy < sv.ploidy()) sv_ploidy = sv.ploidy();
		}
		for (sv_w_samplename_t& sv : clique_svs) if (sv.ploidy() != sv_ploidy) { // all SVs being clustered must have same ploidy
			mtx.lock();
			std::cerr << clique_svs[0].sample << "," << clique_svs[0].id << " and " << sv.sample << "," << sv.id << " have different ploidy." << std::endl;
			mtx.unlock();
		}

		// set GT
		int* gts = new int[n_samples*sv_ploidy];
		std::fill(gts, gts+n_samples*sv_ploidy, bcf_gt_unphased(0));
		for (sv_w_samplename_t& sv : clique_svs) {
			for (int i = 0; i < sv_ploidy; i++) {
				gts[sample2id[sv.sample]*sv_ploidy+i] = sv.gt[i];
			}
		}
		bcf_update_genotypes(out_hdr, vcf_sv, gts, n_samples*sv_ploidy);
		delete[] gts;

		// set CO (original coordinates)
		std::vector<std::string> coos_str(n_samples);
		std::vector<std::string> sizes_str(n_samples);
		std::vector<std::string> prec_str(n_samples);
		std::vector<std::string> incomplete_ass_str(n_samples);
		for (sv_w_samplename_t& sv : clique_svs) {
			if (!coos_str[sample2id[sv.sample]].empty()) {
				coos_str[sample2id[sv.sample]] += ",";
				sizes_str[sample2id[sv.sample]] += ",";
				prec_str[sample2id[sv.sample]] += ",";
				incomplete_ass_str[sample2id[sv.sample]] += ",";
			}
			coos_str[sample2id[sv.sample]] += std::to_string(sv.start+1) + "-" + std::to_string(sv.end+1); // report 1-based
			sizes_str[sample2id[sv.sample]] += std::to_string(sv.svlen);
			prec_str[sample2id[sv.sample]] += (sv.imprecise ? "I" : "P");
			incomplete_ass_str[sample2id[sv.sample]] += (sv.incomplete_ins_seq ? "T" : "F");
		}
		for (int i = 0; i < n_samples; i++) {
			if (coos_str[i].empty()) {
				coos_str[i] = ".";
				sizes_str[i] = ".";
				prec_str[i] = ".";
				incomplete_ass_str[i] = ".";
			}
			coos[i] = coos_str[i].c_str();
			sizes[i] = sizes_str[i].c_str();
			prec[i] = prec_str[i].c_str();
			incomplete_ass[i] = incomplete_ass_str[i].c_str();
		}
		bcf_update_format_string(out_hdr, vcf_sv, "CO", coos, n_samples);
		bcf_update_format_string(out_hdr, vcf_sv, "LN", sizes, n_samples);
		bcf_update_format_string(out_hdr, vcf_sv, "IP", prec, n_samples);
		bcf_update_format_string(out_hdr, vcf_sv, "IC", incomplete_ass, n_samples);

		append_svs.push_back(std::shared_ptr<bcf1_t>(bcf_dup(vcf_sv), bcf_destroy));

		cluster_id++;
    }
    bcf_destroy(vcf_sv);

    delete[] coos;
    delete[] sizes;
    delete[] prec;
	delete[] incomplete_ass;
}

bcf_hdr_t* generate_clustered_vcf_header(std::string command, std::unordered_set<std::string>& called_by) {

	bcf_hdr_t* out_hdr = bcf_hdr_init("w");

	// add contigs
	for (std::string& contig_name : chr_seqs.ordered_contigs) {
		bcf_hrec_t* hrec = generate_contig_hrec();
		int r1 = bcf_hrec_set_val(hrec, 0, contig_name.c_str(), contig_name.length(), false);
		std::string len_str = std::to_string(chr_seqs.get_len(contig_name));
		int r2 = bcf_hrec_set_val(hrec, 1, len_str.c_str(), len_str.length(), false);
		if (r1 || r2) {
			throw std::runtime_error("Failed to create contig to VCF header.");
		}
		bcf_hdr_add_hrec(out_hdr, hrec);
	}

	// add INFO tags
	int len;
	const char* svtype_tag = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the indel (DEL or DUP).\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, svtype_tag, &len));

	const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, end_tag, &len));

	const char* n_samples_tag = "##INFO=<ID=N_SAMPLES,Number=1,Type=Integer,Description=\"Number of samples having at least one SV in this cluster.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, n_samples_tag, &len));

	const char* n_svs_tag = "##INFO=<ID=N_SVS,Number=1,Type=Integer,Description=\"Number of SVs clustered. This may be different from N_SAMPLES as more than "
			"one SV per sample may be clustered.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, n_svs_tag, &len));

	const char* insseq_tag = "##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description=\"Inserted sequence.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, insseq_tag, &len));

	const char* infsvinsseq_tag = "##INFO=<ID=SVINSSEQ_INFERRED,Number=1,Type=String,Description=\"Inferred insertion sequence. When the inserted sequence is too long "
		"to be fully assembled but SurVeyor suspects it to be a transposition, it uses the reference to infer the content of the insertion. Not guaranteed to be accurate. \">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, infsvinsseq_tag, &len));

	const char* svlen_tag = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the resulting SV.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, svlen_tag, &len));

	const char* precise_tag = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Breakpoints and/or inserted sequence are imprecise.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, precise_tag, &len));

	const char* incomplete_ass_tag = "##INFO=<ID=INCOMPLETE_ASSEMBLY,Number=0,Type=Flag,Description=\"Inserted sequence is too long and only part of it could be assembled.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, incomplete_ass_tag, &len));

	// add FORMAT tags
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, gt_tag, &len));

	const char* og_sv_tag = "##FORMAT=<ID=CO,Number=.,Type=String,Description=\"Coordinates of the original SV(s).\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, og_sv_tag, &len));

	const char* ln_sv_tag = "##FORMAT=<ID=LN,Number=.,Type=String,Description=\"Len(s) of the original SV(s). 0 means unavailable, e.g. for incomplete inserted sequences.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, ln_sv_tag, &len));

	const char* pr_sv_tag = "##FORMAT=<ID=IP,Number=.,Type=String,Description=\"Whether the original SV(s) were precise (P) or imprecise (I). "
			"Anything without an explicit IMPRECISE tag is considered precise.\">";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, pr_sv_tag, &len));

	const char* ic_sv_tag = "##FORMAT=<ID=IC,Number=.,Type=String,Description=\"Whether the original insertion had an incomplete assembly. "
				"T is TRUE and F is FALSE.\">";

	std::string cmd_tag = "##SurVeyorCommand=" + command;
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, cmd_tag.c_str(), &len));

	auto now = std::chrono::system_clock::now();
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	std::string version_tag = "##SurVeyorVersion=" + VERSION + "; Date=" + std::ctime(&now_time);
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, version_tag.c_str(), &len));

	for (std::string s : called_by) {
		std::string called_by_tag = "##calledBy=" + s;
		bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, called_by_tag.c_str(), &len));
	}

	std::stringstream clustered_by_ss;
	clustered_by_ss << "##clusteredBy=SurVeyor " << VERSION << "; ";
	clustered_by_ss << "max-dist-precise: " << max_prec_dist << "; ";
	clustered_by_ss << "max-dist-imprecise: " << max_imprec_dist << "; ";
	clustered_by_ss << "min-overlap-precise: " << min_prec_frac_overlap << "; ";
	clustered_by_ss << "min-overlap-imprecise: " << min_imprec_frac_overlap << "; ";
	clustered_by_ss << "max-len-diff-precise: " << max_prec_len_diff << "; ";
	clustered_by_ss << "max-len-diff-imprecise: " << max_imprec_len_diff << "; ";
	clustered_by_ss << "overlap-for-ins: " << (overlap_for_ins ? "true" : "false") << "; ";
	clustered_by_ss << "min-ins-seq-sim: " << min_ins_seq_sim << "; ";
	bcf_hdr_add_hrec(out_hdr, bcf_hdr_parse_line(out_hdr, clustered_by_ss.str().c_str(), &len));

	int i = 0;
	for (std::string& sample_name : sample_names) {
		sample2id[sample_name] = i++; // bcf_hdr_nsamples(out_hdr) does not work unless we sync the header
		bcf_hdr_add_sample(out_hdr, sample_name.c_str());
	}

	return out_hdr;
}

void cluster_contig(int id, std::string contig_name, bcf_hdr_t* out_hdr) {
	mtx.lock();
	std::vector<sv_w_samplename_t>& svs = svs_by_chr[contig_name];
	std::vector<std::shared_ptr<bcf1_t>>& vcf_clustered_svs = clustered_svs_by_chr[contig_name];
	mtx.unlock();

    std::sort(svs.begin(), svs.end(), [](const sv_w_samplename_t& sv1, const sv_w_samplename_t& sv2) {
        return sv1.start < sv2.start;
    });

    std::vector<std::vector<int>> compatibility_list(svs.size());

	int i_prev = 0;
	for (int i = 0; i < svs.size(); i++) {
		sv_w_samplename_t& sv1 = svs[i];
		for (int j = i+1; j < svs.size(); j++) {
			sv_w_samplename_t& sv2 = svs[j];
			if (sv2.start-sv1.start > max_dist) break;
			if (is_compatible(sv1, sv2)) {
				compatibility_list[i].push_back(j);
				compatibility_list[j].push_back(i);
			}
		}
		compatibility_list[i].shrink_to_fit();

		if (i > 0 && sv1.start-svs[i-1].start > max_dist) {
			std::vector<std::vector<int>> cliques = compute_minimal_clique_cover(i_prev, i, svs, compatibility_list);
			print_cliques(cliques, svs, out_hdr, chr_seqs, vcf_clustered_svs);

			for (int j = i_prev; j < i; j++) {
				compatibility_list[j].clear();
				compatibility_list[j].shrink_to_fit();
			}
			i_prev = i;
		}
	}
	std::vector<std::vector<int>> cliques = compute_minimal_clique_cover(i_prev, svs.size(), svs, compatibility_list);
	print_cliques(cliques, svs, out_hdr, chr_seqs, vcf_clustered_svs);
}

void read_svs(int id, std::string sample_sv_fpath, std::string sample_name) {
	mtx.lock();
	std::cout << "Reading SVs from " << sample_sv_fpath << std::endl;
	mtx.unlock();

	htsFile* sample_sv_file = bcf_open(sample_sv_fpath.c_str(), "r");
	bcf_hdr_t* vcf_header = bcf_hdr_read(sample_sv_file);
	if (vcf_header == NULL) {
		throw std::runtime_error("Failed to read the VCF header of " + sample_sv_fpath + ".");
	}

	for (int i = 0; i < vcf_header->nhrec; i++) {
		if (strcmp(vcf_header->hrec[i]->key, "calledBy") == 0) {
			called_by_mtx.lock();
			called_by.insert(vcf_header->hrec[i]->value);
			called_by_mtx.unlock();
		}
	}

	bcf1_t* vcf_record = bcf_init();
	std::unordered_map<std::string, std::vector<sv_w_samplename_t> > local_svs_by_chr;
	while (bcf_read(sample_sv_file, vcf_header, vcf_record) == 0) {
		if (!keep_all_called && count_alt_alleles(vcf_header, vcf_record) < 1) continue; // skip records without ALT alleles
		
		auto sv = std::shared_ptr<sv_t>(bcf_to_sv(vcf_header, vcf_record));
		if (sv == nullptr) {
			ignored_unsupported_svs++;
			continue;
		} else if (sv->svtype() != "DEL" && sv->svtype() != "INS" && sv->svtype() != "DUP" && sv->svtype() != "INV") {
			ignored_unsupported_svs++;
			continue;
		}

		if (sv->is_pass()) {	
			sv_w_samplename_t sv_w_sample(sv, sample_name);
			std::string seqname = bcf_seqname(vcf_header, vcf_record);
			local_svs_by_chr[seqname].push_back(sv_w_sample);
		}
	}
	bcf_destroy(vcf_record);
	bcf_hdr_destroy(vcf_header);
	hts_close(sample_sv_file);
	
	samples_mtx.lock();
	sample_names.push_back(sample_name);
	samples_mtx.unlock();

	mtx.lock();
	for (auto& e : local_svs_by_chr) {
		std::string chr_name = e.first;
		std::vector<sv_w_samplename_t>& svs = e.second;
		svs_by_chr[chr_name].insert(svs_by_chr[chr_name].end(), svs.begin(), svs.end());
	}
	mtx.unlock();
}

int main(int argc, char* argv[]) {

	cxxopts::Options options("SurVClusterer", "Clusters SVs across different samples.");

	options.add_options()
			("filelist", "Filelist containing a line for each sample being clustered. The line must contain the name of the sample and the "
					"path of its VCF file, space or tab separated.", cxxopts::value<std::string>())
			("reference", "Reference genome in fasta format", cxxopts::value<std::string>())
			("o,out-prefix", "Files out_prefix.vcf.gz will be produced.", cxxopts::value<std::string>()->default_value("clustered"))
			("d,max-dist-precise", "Maximum distance allowed between the breakpoints of two precise variants.",
					cxxopts::value<int>()->default_value("100"))
			("D,max-dist-imprecise", "Maximum distance allowed between the breakpoints when at least one of the variants is imprecise.",
					cxxopts::value<int>()->default_value("500"))
			("x,min-overlap-precise", "Minimum overlap required between two precise variants, expressed as the fraction of the shortest variant"
					" that is covered by the longest variant.", cxxopts::value<double>()->default_value("0.8"))
			("X,min-overlap-imprecise", "Minimum overlap required between two variants when at least one is imprecise, "
					"expressed as the fraction of the shortest variant that is covered by the longest variant.",
					cxxopts::value<double>()->default_value("0.5"))
			("s,max-len-diff-precise", "Maximum length difference allowed between two precise variants.",
					cxxopts::value<int>()->default_value("100"))
			("S,max-len-diff-imprecise", "Maximum length difference allowed when at least one variant is imprecise.",
					cxxopts::value<int>()->default_value("500"))
			("l,min-len-ratio-precise", "Minimum length ratio allowed (smallest variant length / largest variant length) between two precise variants.", cxxopts::value<double>()->default_value("0.8"))
			("L,min-len-ratio-imprecise", "Minimum length ratio allowed (smallest variant length / largest variant length) when at least one variant is imprecise.", cxxopts::value<double>()->default_value("0.5"))
			("q,min-ins-seq-sim", "Minimum inserted-sequence similarity required between two complete insertions. Similarity is 2 * aligned matching bases / (len(S1) + len(S2)); 0 disables sequence comparison.",
					cxxopts::value<double>()->default_value("0"))
			("i,overlap-for-ins", "Require overlap for insertions.", cxxopts::value<bool>()->default_value("false"))
			("k,keep-all-called", "Keep all SVs, even if they have no alt alleles.", cxxopts::value<bool>()->default_value("false"))
			("t,threads", "Maximum number of threads used.", cxxopts::value<int>()->default_value("1"))
			("v,version", "Print the version number and exit.")
			("h,help", "Print usage");

	options.parse_positional({"filelist", "reference"});
	options.positional_help("filelist reference");
	options.show_positional_help();
	auto parsed_args = options.parse(argc, argv);

	if (parsed_args.count("version")) {
			std::cout << "SurVClusterer v" << VERSION << std::endl;
			exit(0);
	} else if (parsed_args.count("help") || !parsed_args.count("filelist") || !parsed_args.count("reference")) {
		std::cout << options.help() << std::endl;
		exit(0);
	}

	std::string full_cmd;
	for (int i = 0; i < argc; i++) {
		if (i > 0) full_cmd += " ";
		full_cmd += argv[i];
	}

    std::ifstream file_list(parsed_args["filelist"].as<std::string>());
    std::string reference_fname = parsed_args["reference"].as<std::string>();
    max_prec_dist = parsed_args["max-dist-precise"].as<int>();
    max_imprec_dist = parsed_args["max-dist-imprecise"].as<int>();
    max_dist = std::max(max_prec_dist, max_imprec_dist);
    min_prec_frac_overlap = parsed_args["min-overlap-precise"].as<double>();
    min_imprec_frac_overlap = parsed_args["min-overlap-imprecise"].as<double>();
    max_prec_len_diff = parsed_args["max-len-diff-precise"].as<int>();
    max_imprec_len_diff = parsed_args["max-len-diff-imprecise"].as<int>();
	min_prec_len_ratio = parsed_args["min-len-ratio-precise"].as<double>();
	min_imprec_len_ratio = parsed_args["min-len-ratio-imprecise"].as<double>();
	min_ins_seq_sim = parsed_args["min-ins-seq-sim"].as<double>();
	if (min_ins_seq_sim < 0.0 || min_ins_seq_sim > 1.0) {
		std::cerr << "--min-ins-seq-sim must be between 0.0 and 1.0." << std::endl;
		exit(1);
	}
    overlap_for_ins = parsed_args["overlap-for-ins"].as<bool>();
    std::string out_prefix = parsed_args["out-prefix"].as<std::string>();
    int n_threads = parsed_args["threads"].as<int>();
	keep_all_called = parsed_args.count("keep-all-called");

    std::string out_vcf_fname = out_prefix + ".vcf.gz";
    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");

    chr_seqs.read_fasta_into_map(reference_fname);

    std::string sample_name, sample_sv_fpath;
	ctpl::thread_pool thread_pool_read_svs(n_threads);
	std::vector<std::future<void> > futures;
    while (file_list >> sample_name >> sample_sv_fpath) {
		std::future<void> future = thread_pool_read_svs.push(read_svs, sample_sv_fpath, sample_name);
		futures.push_back(std::move(future));
    }
	thread_pool_read_svs.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

    std::cout << "Finished reading SVs." << std::endl;

	std::cout << sample_names.size() << " samples were read." << std::endl;

    bcf_hdr_t* out_hdr = generate_clustered_vcf_header(full_cmd, called_by);
    if (bcf_hdr_write(out_vcf_file, out_hdr) != 0) {
    	throw std::runtime_error("Could not write header to " + std::string(out_vcf_file->fn));
    }

    int nseqs;
    const char** contig_names = bcf_hdr_seqnames(out_hdr, &nseqs);
	ctpl::thread_pool thread_pool_cc(n_threads);
    for (int i = 0; i < nseqs; i++) {
		std::future<void> future = thread_pool_cc.push(cluster_contig, contig_names[i], out_hdr);
		futures.push_back(std::move(future));
    }
	thread_pool_cc.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		futures[i].get();
	}

    for (int i = 0; i < nseqs; i++) {
    	for (std::shared_ptr<bcf1_t> b : clustered_svs_by_chr[contig_names[i]]) {
			if (bcf_write(out_vcf_file, out_hdr, b.get()) != 0) {
				throw std::runtime_error("Failed to write record to " + std::string(out_vcf_file->fn));
			}
    	}
    }

    bcf_hdr_destroy(out_hdr);
    hts_close(out_vcf_file);

    std::cout << ignored_unsupported_svs << " unsupported SVs were ignored." << std::endl;
    std::cout << "Finished." << std::endl;
}
