#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include <mutex>

#include "../libs/cxxopts.h"
#include "../libs/IntervalTree.h"
#include "../libs/ssw_cpp.h"
#include "../libs/cptl_stl.h"
#include "../libs/ssw.h"
#include "common.h"
#include "../src/sw_utils.h"
#include "htslib/vcf.h"

chr_seqs_map_t chr_seqs;
int max_prec_dist, max_imprec_dist, max_repeat_dist;
double min_prec_frac_overlap, min_imprec_frac_overlap;
int max_prec_len_diff, max_imprec_len_diff;
double min_prec_len_ratio, min_imprec_len_ratio;
std::unordered_set<std::string> bdup_ids, cdup_ids;

std::vector<std::shared_ptr<sv_t>> benchmark_svs;
std::unordered_map<std::string, std::vector<std::shared_ptr<sv_t>>> called_dels_by_chr, called_inss_by_chr, called_invs_by_chr;
std::unordered_map<std::string, IntervalTree<repeat_t>*> reps_i;
std::unordered_map<std::string, IntervalTree<std::shared_ptr<sv_t>>*> dels_csv_i, inss_csv_i, invs_csv_i;

StripedSmithWaterman::Aligner edit_distance_aligner(1, 0, 1, 1, false);
bool ignore_seq = false;

std::mutex mtx;

bool compatible_gts(sv_t* sv1, sv_t* sv2) {
	int min_sv1_ac = sv1->allele_count(1), min_sv2_ac = sv2->allele_count(1);
	int max_sv1_ac = min_sv1_ac + sv1->missing_alleles(), max_sv2_ac = min_sv2_ac + sv2->missing_alleles();
	bool ac_iv_overlap = !(max_sv1_ac < min_sv2_ac || max_sv2_ac < min_sv1_ac);
	return ac_iv_overlap;
}

std::string get_bsv_type(sv_t* sv) {
	if (bdup_ids.count(sv->id)) return "DUP";
	else return sv->svtype();
}
std::string get_csv_type(sv_t* sv) {
	if (cdup_ids.count(sv->id)) return "DUP";
	else return sv->svtype();
}

int len_diff(sv_t* sv1, sv_t* sv2) {
	if (sv1->svtype() == "DEL" && sv2->svtype() == "DEL") return abs(sv1->svlen()-sv2->svlen());
	else if (sv1->svtype() == "DUP" && sv2->svtype() == "DUP") return 0;
	else if (sv1->svtype() == "INS" && sv2->svtype() == "INS") {
		if (sv1->incomplete_ins_seq() && sv2->incomplete_ins_seq()) {
			return 0;
		} else if (sv1->incomplete_ins_seq() || sv2->incomplete_ins_seq()) {
			sv_t* sv_full = sv1->incomplete_ins_seq() ? sv2 : sv1;
			sv_t* sv_incmpl = sv1->incomplete_ins_seq() ? sv1 : sv2;
			return std::max(0, (int) (sv_incmpl->ins_seq.length() - sv_full->svlen()));
		} else {
			return abs((int) (sv1->svlen()-sv2->svlen()));
		}
	} else if ((sv1->svtype() == "INS" && sv2->svtype() == "DUP") || (sv1->svtype() == "DUP" && sv2->svtype() == "INS")) return 0;
	else if (sv1->svtype() == "INV" && sv2->svtype() == "INV") return abs(sv1->svlen()-sv2->svlen());
	else return INT32_MAX;
}

double len_ratio(hts_pos_t svlen1, hts_pos_t svlen2) {
	svlen1 = abs(svlen1);
	svlen2 = abs(svlen2);
	if (svlen1 == 0 && svlen2 == 0) return 1.0;
	else return double(std::min(svlen1, svlen2)) / std::max(svlen1, svlen2);
}

double len_ratio(sv_t* sv1, sv_t* sv2) {
	if (sv1->svtype() == "DEL" && sv2->svtype() == "DEL") {
		return len_ratio(sv1->svlen(), sv2->svlen());
	} else if (sv1->svtype() == "DUP" && sv2->svtype() == "DUP") return 1.0;
	else if (sv1->svtype() == "INS" && sv2->svtype() == "INS") {
		if (sv1->incomplete_ins_seq() && sv2->incomplete_ins_seq()) { // if both are incomplete, we cannot compare lengths
			return 1.0;
		} else if (sv1->incomplete_ins_seq() || sv2->incomplete_ins_seq()) {
			// if one is incomplete, we know a lower bound for its length - if it is longer than the full one, we can compute an upper bound for the ratio
			sv_t* sv_full = sv1->incomplete_ins_seq() ? sv2 : sv1;
			sv_t* sv_incmpl = sv1->incomplete_ins_seq() ? sv1 : sv2;
			if (sv_incmpl->svlen() > sv_full->svlen()) {
				return double(sv_full->svlen()) / sv_incmpl->svlen();
			} else {
				return 1.0;
			}
		} else {
			return len_ratio(sv1->svlen(), sv2->svlen());
		}
	} else if ((sv1->svtype() == "INS" && sv2->svtype() == "DUP") || (sv1->svtype() == "DUP" && sv2->svtype() == "INS")) {
		return 1.0;
	} else if (sv1->svtype() == "INV" && sv2->svtype() == "INV") {
		return len_ratio(sv1->svlen(), sv2->svlen());
	} else {
		return 0.0; // not compatible
	}
}

char* generate_alt_allele(sv_t* sv, hts_pos_t start, hts_pos_t end) {
	start = std::min(sv->start, start);
	end = std::max(sv->end, end);
	int len = sv->start - start + sv->ins_seq.length() + end - sv->end;
	char* alt = new char[len+1];
	char* chr_seq = chr_seqs.get_seq(sv->chr);
	strncpy(alt, chr_seq+start+1, sv->start-start);
	strncpy(alt+sv->start-start, sv->ins_seq.c_str(), sv->ins_seq.length());
	strncpy(alt+sv->start-start+sv->ins_seq.length(), chr_seq+sv->end+1, end-sv->end);
	alt[len] = '\0';
	return alt;
}

bool alt_allele_match(sv_t* sv1, sv_t* sv2, int max_score_loss) {
	hts_pos_t start = std::min(sv1->start, sv2->start) - 100;
	if (start < 0) start = 0;
	hts_pos_t end = std::max(sv1->end, sv2->end) + 100;
	if (end >= chr_seqs.get_len(sv1->chr)) end = chr_seqs.get_len(sv1->chr)-1;
	char* alt1 = generate_alt_allele(sv1, start, end);
	char* alt2 = generate_alt_allele(sv2, start, end);
	int strlen1 = strlen(alt1);
	int strlen2 = strlen(alt2);
	StripedSmithWaterman::Alignment alignment;
	StripedSmithWaterman::Filter filter;
	edit_distance_aligner.Align(alt1, alt2, strlen2, filter, &alignment, 0);
	int left_clip = std::max(alignment.ref_begin, alignment.query_begin);
	int right_clip = std::max(strlen1 - (alignment.query_end + 1), strlen2 - (alignment.ref_end + 1));
	int score_loss = std::min(strlen1, strlen2) - alignment.sw_score + left_clip + right_clip;
	delete [] alt1;
	delete [] alt2;
	return score_loss <= max_score_loss;
}

bool is_compatible_ivals(sv_t* sv1, sv_t* sv2, bool repeat_mode) {
	bool imprecise_mode = sv1->imprecise || sv2->imprecise;
	int max_dist = imprecise_mode ? max_imprec_dist : max_prec_dist;
	if (repeat_mode) max_dist = max_repeat_dist;
	double min_frac_overlap = imprecise_mode ? min_imprec_frac_overlap : min_prec_frac_overlap;
	int max_len_diff = imprecise_mode ? max_imprec_len_diff : max_prec_len_diff;
	double min_len_ratio = imprecise_mode ? min_imprec_len_ratio : min_prec_len_ratio;
	return distance(sv1, sv2) <= max_dist &&
			(repeat_mode || overlap(sv1, sv2) >= min_frac_overlap) &&
			len_diff(sv1, sv2) <= max_len_diff &&
			len_ratio(sv1, sv2) >= min_len_ratio;
}

bool is_compatible_del_del(sv_t* sv1, sv_t* sv2, bool repeat_mode) {
	if (distance(sv1, sv2) > max_repeat_dist) return false;
	bool imprecise_mode = sv1->imprecise || sv2->imprecise;
	int max_dist = imprecise_mode ? max_imprec_dist : max_prec_dist;
	double min_len_ratio = imprecise_mode ? min_imprec_len_ratio : min_prec_len_ratio;
	int max_len_diff = imprecise_mode ? max_imprec_len_diff : max_prec_len_diff;
	int max_svlen = std::max(std::abs(sv1->svlen()), std::abs(sv2->svlen()));
	int max_score_loss = std::min({int(max_svlen * (1-min_len_ratio)), max_dist, max_len_diff});
	if (sv1->svsize() <= 50) { // for small deletions, only check that they produce similar alt alleles
		return alt_allele_match(sv1, sv2, max_score_loss);
	} else { // for large deletions, perform interval-based comparison first since alt allele generation is expensive
		return is_compatible_ivals(sv1, sv2, repeat_mode) &&
				alt_allele_match(sv1, sv2, max_score_loss);
	}
}
bool is_compatible_dup_dup(sv_t* sv1, sv_t* sv2, bool repeat_mode) {
	return is_compatible_ivals(sv1, sv2, repeat_mode);
}

bool check_ins_dup_seq(sv_t* ins_sv, sv_t* dup_sv, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Alignment& alignment) {
	if (ignore_seq) return true;

	int ins_seq_len = ins_sv->ins_seq.length();
	int max_len_diff = (ins_sv->imprecise || dup_sv->imprecise) ? max_imprec_len_diff : max_prec_len_diff;
	double min_len_ratio = (ins_sv->imprecise || dup_sv->imprecise) ? min_imprec_len_ratio : min_prec_len_ratio;
	if (dup_sv->svlen() > ins_seq_len+max_len_diff || double(ins_seq_len)/dup_sv->svlen() < min_len_ratio) return false;
	if (dup_sv->svlen() <= 0) return false;

	char* dup_seq = new char[dup_sv->end-dup_sv->start+1];
	strncpy(dup_seq, chr_seqs.get_seq(dup_sv->chr)+dup_sv->start, dup_sv->end-dup_sv->start);
	dup_seq[dup_sv->end-dup_sv->start] = '\0';
	std::string ext_dup_seq = dup_seq + dup_sv->ins_seq;
	int add_len = dup_sv->end - dup_sv->start + dup_sv->ins_seq.length();
	while (ext_dup_seq.length() <= ins_seq_len && ext_dup_seq.length() < UINT16_MAX) {
		ext_dup_seq += dup_seq + dup_sv->ins_seq;
	}
	delete[] dup_seq;

	if (ins_sv->start < dup_sv->start) {
		size_t rotation = dup_sv->start - ins_sv->start;
		rotation %= dup_sv->end-dup_sv->start;
		ext_dup_seq = ext_dup_seq.substr(dup_sv->end-dup_sv->start-rotation) + ext_dup_seq.substr(0, dup_sv->end-dup_sv->start-rotation);
	} else if (ins_sv->start > dup_sv->start) {
		size_t rotation = ins_sv->start - dup_sv->start;
		rotation %= dup_sv->end-dup_sv->start;
		ext_dup_seq = ext_dup_seq.substr(rotation) + ext_dup_seq.substr(0, rotation);
	}

	// ssw score is a uint16_t, so we cannot compare strings longer than that
	StripedSmithWaterman::Filter filter;
	std::string& smaller_seq = ins_seq_len < ext_dup_seq.length() ? ins_sv->ins_seq : ext_dup_seq;
	std::string& larger_seq = ins_seq_len < ext_dup_seq.length() ? ext_dup_seq : ins_sv->ins_seq;
	if (larger_seq.length() > UINT16_MAX) {
		std::string smaller_seq_prefix = smaller_seq.substr(0, UINT16_MAX-10);
		std::string smaller_seq_suffix = smaller_seq.substr(std::max(0, (int) smaller_seq.length()-(UINT16_MAX-10)));
		std::string larger_seq_prefix = larger_seq.substr(0, UINT16_MAX-10);
		std::string larger_seq_suffix = larger_seq.substr(std::max(0, (int) larger_seq.length()-(UINT16_MAX-10)));
		aligner.Align(smaller_seq_prefix.data(), larger_seq_prefix.data(), larger_seq_prefix.length(), filter, &alignment, 0);
		if (alignment.query_end-alignment.query_begin < smaller_seq_prefix.length()*0.8) return false;
		aligner.Align(smaller_seq_suffix.data(), larger_seq_suffix.data(), larger_seq_suffix.length(), filter, &alignment, 0);
		if (alignment.query_end-alignment.query_begin < smaller_seq_suffix.length()*0.8) return false;
		return true;
	}

	aligner.Align(smaller_seq.data(), larger_seq.data(), larger_seq.length(), filter, &alignment, 0);
	return alignment.query_end-alignment.query_begin >= smaller_seq.length()*0.8;
}
bool is_compatible_ins_dup(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner, bool repeat_mode) {
	StripedSmithWaterman::Alignment alignment;
	int max_dist = (sv1->imprecise || sv2->imprecise) ? max_imprec_dist : max_prec_dist;
	if (repeat_mode) max_dist = max_repeat_dist;
	return distance(sv1, sv2) <= max_dist && check_ins_dup_seq(sv1, sv2, aligner, alignment);
}

bool check_cmpl_cmpl_seq(std::string& cmpl_seq1, std::string& cmpl_seq2, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Alignment& alignment) {
	std::string& query = cmpl_seq1.length() < cmpl_seq2.length() ? cmpl_seq1 : cmpl_seq2;
	std::string& ref = cmpl_seq1.length() < cmpl_seq2.length() ? cmpl_seq2 : cmpl_seq1;
	StripedSmithWaterman::Filter filter;
	aligner.Align(query.data(), ref.data(), ref.length(), filter, &alignment, 0);
	return alignment.query_end-alignment.query_begin >= query.length()*0.8;
}

bool check_cmpl_incmpl_seq(std::string& cmpl_seq, std::string& incmpl_seq, bool incmpl_is_imprecise, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Alignment& alignment) {
	int dash_pos = incmpl_seq.find('-');
	std::string left_seq = incmpl_seq.substr(0, dash_pos);
	std::string right_seq = incmpl_seq.substr(dash_pos+1);
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment left_aln, right_aln;
	int max_len_diff = incmpl_is_imprecise ? max_imprec_len_diff : max_prec_len_diff;
	left_aln.Clear();
	right_aln.Clear();
	aligner.Align(left_seq.data(), cmpl_seq.data(), cmpl_seq.length(), filter, &left_aln, 0);
	aligner.Align(right_seq.data(), cmpl_seq.data(), cmpl_seq.length(), filter, &right_aln, 0);
	if (left_aln.query_end-left_aln.query_begin+right_aln.query_end-right_aln.query_begin < incmpl_seq.length()*0.8 ||
		left_aln.ref_begin > max_len_diff || cmpl_seq.length()-right_aln.ref_end > max_len_diff) return false;
	// set fields required to compute match score
	alignment.sw_score = left_aln.sw_score + right_aln.sw_score;
	alignment.query_end = cmpl_seq.length();
	return true;
}

bool check_incmpl_incmpl_seq(std::string& incmpl_seq1, std::string& incmpl_seq2, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Alignment& alignment) {
	// TODO
	return true;
}

bool check_ins_ins_seq(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner, StripedSmithWaterman::Alignment& alignment) {
	if (ignore_seq) return true;

	int max_len_diff = (sv1->imprecise || sv2->imprecise) ? max_imprec_len_diff : max_prec_len_diff;
	double min_len_ratio = (sv1->imprecise || sv2->imprecise) ? min_imprec_len_ratio : min_prec_len_ratio;
	if (len_diff(sv1, sv2) > max_len_diff || len_ratio(sv1, sv2) < min_len_ratio) return false;

	sv_t* l_sv = sv1->start < sv2->start ? sv1 : sv2;
	sv_t* r_sv = sv1->start < sv2->start ? sv2 : sv1;
	char* extra_seq = new char[r_sv->start-l_sv->start+1];
	strncpy(extra_seq, chr_seqs.get_seq(sv1->chr)+l_sv->start, r_sv->start-l_sv->start);
	extra_seq[r_sv->start-l_sv->start] = '\0';

	std::string lsv_seq = l_sv->ins_seq + extra_seq;
	std::string rsv_seq = extra_seq + r_sv->ins_seq;
	delete[] extra_seq;

	if (lsv_seq.length() > UINT16_MAX || rsv_seq.length() > UINT16_MAX) return true; // TODO: ssw score is a uint16_t, so we cannot compare strings longer than that

	if (!sv1->incomplete_ins_seq() && !sv2->incomplete_ins_seq()) {
		return check_cmpl_cmpl_seq(lsv_seq, rsv_seq, aligner, alignment);
	} else if (l_sv->incomplete_ins_seq() && !r_sv->incomplete_ins_seq()) {
		return check_cmpl_incmpl_seq(rsv_seq, lsv_seq, r_sv->imprecise, aligner, alignment);
	} else if (!l_sv->incomplete_ins_seq() && r_sv->incomplete_ins_seq()) {
		return check_cmpl_incmpl_seq(lsv_seq, rsv_seq, l_sv->imprecise, aligner, alignment);
	} else {
		return check_incmpl_incmpl_seq(lsv_seq, rsv_seq, aligner, alignment);
	}
}
bool is_compatible_ins_ins(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Alignment alignment;
	int max_dist = (sv1->imprecise || sv2->imprecise) ? max_imprec_dist : max_prec_dist;
	if (distance(sv1, sv2) > max_dist) return false;
	if (!ignore_seq && !sv1->incomplete_ins_seq() && !sv2->incomplete_ins_seq() && sv1->svlen() < 1000 && sv2->svlen() < 1000) {
		int max_len_diff = (sv1->imprecise || sv2->imprecise) ? max_imprec_len_diff : max_prec_len_diff;
		double min_len_ratio = (sv1->imprecise || sv2->imprecise) ? min_imprec_len_ratio : min_prec_len_ratio;
		if (len_diff(sv1, sv2) > max_len_diff || len_ratio(sv1, sv2) < min_len_ratio) return false;

		int max_svlen = std::max(std::abs(sv1->svlen()), std::abs(sv2->svlen()));
		int max_score_loss = max_svlen * (1-min_len_ratio);
		return alt_allele_match(sv1, sv2, max_score_loss);
	} else {
		return check_ins_ins_seq(sv1, sv2, aligner, alignment);
	}
}

bool check_ins_seq(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner) {
	StripedSmithWaterman::Alignment alignment;
	if (sv1->svtype() == "DUP" && sv2->svtype() == "DUP") {
		return true;
	} else if (sv1->svtype() == "DUP" && sv2->svtype() == "INS") {
		return check_ins_dup_seq(sv2, sv1, aligner, alignment);
	} else if (sv1->svtype() == "INS" && sv2->svtype() == "DUP") {
		return check_ins_dup_seq(sv1, sv2, aligner, alignment);
	} else {
		return check_ins_ins_seq(sv1, sv2, aligner, alignment);
	}
}

bool is_compatible_inv_inv(sv_t* sv1, sv_t* sv2) {
	return is_compatible_ivals(sv1, sv2, false);
}

bool is_compatible(sv_t* sv1, sv_t* sv2, StripedSmithWaterman::Aligner& aligner, bool repeat_mode) {
	if (sv1->svtype() == "DEL" && sv2->svtype() == "DEL") {
		return is_compatible_del_del(sv1, sv2, repeat_mode);
	} else if (sv1->svtype() == "DUP" && sv2->svtype() == "DUP") {
		return is_compatible_dup_dup(sv1, sv2, repeat_mode);
	} else if (sv1->svtype() == "DUP" && sv2->svtype() == "INS") {
		return is_compatible_ins_dup(sv2, sv1, aligner, repeat_mode);
	} else if (sv1->svtype() == "INS" && sv2->svtype() == "DUP") {
		return is_compatible_ins_dup(sv1, sv2, aligner, repeat_mode);
	} else if (sv1->svtype() == "INS" && sv2->svtype() == "INS") {
		return is_compatible_ins_ins(sv1, sv2, aligner);
	} else if (sv1->svtype() == "INV" && sv2->svtype() == "INV") {
		return is_compatible_inv_inv(sv1, sv2);
	} else {
		return false;
	}
}

bool is_exact_match(sv_t* sv1, sv_t* sv2) {
	if (sv1 == NULL || sv2 == NULL) return false;

	if (sv1->svtype() != sv2->svtype()) return false;
	if (sv1->chr != sv2->chr || sv1->start != sv2->start || sv1->end != sv2->end) return false;
	if (sv1->ins_seq != sv2->ins_seq) return false;
	return true;
}

struct sv_match_t {
	std::shared_ptr<sv_t> b_sv, c_sv;
	double score = 0;
	bool rep, exact, valid = true;

	sv_match_t(std::shared_ptr<sv_t> b_sv, std::shared_ptr<sv_t> c_sv, bool rep, StripedSmithWaterman::Aligner& aligner) : b_sv(b_sv), c_sv(c_sv) {
		if (b_sv == NULL || c_sv == NULL) {
			this->score = 0;
			this->rep = rep;
			this->exact = false;
			return;
		}

		this->exact = is_exact_match(b_sv.get(), c_sv.get());

		StripedSmithWaterman::Alignment alignment;
		alignment.Clear();
		int len_diff = 0;
		if (b_sv->svtype() == "INS" && c_sv->svtype() == "INS") {
			len_diff = abs(b_sv->svlen()-c_sv->svlen());
			check_ins_ins_seq(b_sv.get(), c_sv.get(), aligner, alignment);
		} else if (b_sv->svtype() == "INS" && c_sv->svtype() == "DUP") {
			check_ins_dup_seq(b_sv.get(), c_sv.get(), aligner, alignment);
		} else if (b_sv->svtype() == "DUP" && c_sv->svtype() == "INS") {
			check_ins_dup_seq(c_sv.get(), b_sv.get(), aligner, alignment);
		} else {
			len_diff = abs(b_sv->svlen()-c_sv->svlen());
		}

		int right_clip = get_right_clip_size(alignment);
		double dist_log = log(distance(b_sv.get(), c_sv.get())+1);
		double aln_score = alignment.sw_score/std::max(1.0, double(alignment.query_end+right_clip+1)) * 100;
		this->score = -len_diff - dist_log + aln_score;
		this->rep = rep;
	}
};
std::vector<sv_match_t> matches;

void find_match(int id, int start_idx, int end_idx) {

	StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	int query_dist = std::max({max_prec_dist, max_imprec_dist, max_repeat_dist});

	for (int i = start_idx; i < end_idx; i++) {
		std::shared_ptr<sv_t> bsv = benchmark_svs[i];

		IntervalTree<std::shared_ptr<sv_t>>* called_svs_chr_type = nullptr;
		if (bsv->svtype() == "DEL") {
			auto it = dels_csv_i.find(bsv->chr);
			if (it != dels_csv_i.end()) called_svs_chr_type = it->second;
		} else if (bsv->svtype() == "INS" || bsv->svtype() == "DUP") {
			auto it = inss_csv_i.find(bsv->chr);
			if (it != inss_csv_i.end()) called_svs_chr_type = it->second;
		} else if (bsv->svtype() == "INV") {
			auto it = invs_csv_i.find(bsv->chr);
			if (it != invs_csv_i.end()) called_svs_chr_type = it->second;
		} else {
			continue;
		}
		if (called_svs_chr_type == nullptr) continue;

		std::vector<repeat_t> reps_containing_bsv;
		auto it = reps_i.find(bsv->chr);
		if (it != reps_i.end() && it->second != nullptr) {
			std::vector<Interval<repeat_t>> intervals_temp = it->second->findOverlapping(bsv->start, bsv->end);
			for (auto &iv : intervals_temp) {
				repeat_t rep = iv.value;
				if (rep.contains(bsv.get())) {
					reps_containing_bsv.push_back(rep);
				}
			}
		}

		hts_pos_t query_start = std::max(hts_pos_t(0), bsv->start-query_dist);
		hts_pos_t query_end = bsv->end + query_dist;
		std::vector<Interval<std::shared_ptr<sv_t>>> candidate_intervals = called_svs_chr_type->findOverlapping(query_start, query_end);
		for (const Interval<std::shared_ptr<sv_t>>& candidate_interval : candidate_intervals) {
			const std::shared_ptr<sv_t>& csv = candidate_interval.value;
			if (is_compatible(bsv.get(), csv.get(), aligner, false)) {
				sv_match_t match(bsv, csv, false, aligner);
				mtx.lock();
				matches.push_back(match);
				mtx.unlock();
				continue;
			}

			if (is_compatible(bsv.get(), csv.get(), aligner, true)) {
				bool same_tr = false;
				for (repeat_t& rep : reps_containing_bsv) {
					if (rep.intersects(csv.get())) {
						same_tr = true;
						break;
					}
				}

				if (same_tr) {
					sv_match_t match(bsv, csv, true, aligner);
					mtx.lock();
					matches.push_back(match);
					mtx.unlock();
				}
			}
		}
	}
}

std::unordered_set<std::string> find_dup_ids(std::vector<std::shared_ptr<sv_t>>& svs) {
	std::unordered_set<std::string> seen_ids, dup_ids;
	for (const std::shared_ptr<sv_t>& sv : svs) {
		if (seen_ids.count(sv->id)) {
			dup_ids.insert(sv->id);
		} else {
			seen_ids.insert(sv->id);
		}
	}
	return dup_ids;
}

bool is_known_hom_alt(sv_t* sv) {
	return sv->sample_info.gt.size() > 0 &&
		sv->missing_alleles() == 0 &&
		sv->allele_count(1) == sv->sample_info.gt.size();
}

bool is_known_het_alt(sv_t* sv) {
	return sv->sample_info.gt.size() == 2 &&
		sv->missing_alleles() == 0 &&
		sv->allele_count(1) == 1;
}

bool is_known_het_alt(const std::vector<int>& gt) {
	if (gt.size() != 2) return false;
	int alt_alleles = 0;
	for (int allele : gt) {
		if (bcf_gt_is_missing(allele)) return false;
		alt_alleles += (bcf_gt_allele(allele) == 1);
	}
	return alt_alleles == 1;
}

std::vector<snp_t> read_het_benchmark_snps(const std::string& filename) {
	std::vector<snp_t> snps;
	htsFile* file = bcf_open(filename.c_str(), "r");
	if (file == NULL) {
		throw std::runtime_error("Error: could not open file " + filename);
	}

	bcf_hdr_t* hdr = bcf_hdr_read(file);
	if (hdr == NULL) {
		bcf_close(file);
		throw std::runtime_error("Error: could not read header from file " + filename);
	}

	bcf1_t* record = bcf_init();
	while (bcf_read(file, hdr, record) == 0) {
		std::string svtype = get_sv_type(hdr, record);
		if (svtype != "SNP" && svtype != "SNV") continue;
		if (!is_known_het_alt(get_bcf_gt(hdr, record))) continue;

		bcf_unpack(record, BCF_UN_STR);
		snps.push_back(snp_t(bcf_seqname_safe(hdr, record), record->pos, record->d.allele[1][0]));
	}

	bcf_destroy(record);
	bcf_hdr_destroy(hdr);
	bcf_close(file);
	return snps;
}

std::unordered_set<std::string> find_vars_overlapping_other_nonref(const std::vector<std::shared_ptr<sv_t>>& svs,
														   const std::vector<snp_t>& benchmark_het_snps) {
	std::unordered_map<std::string, std::vector<Interval<std::string, hts_pos_t>>> nonref_ivals_by_chr;
	for (const std::shared_ptr<sv_t>& sv : svs) {
		if (sv.get()->allele_count(1) == 0) continue;

		nonref_ivals_by_chr[sv->chr].push_back(Interval<std::string, hts_pos_t>(sv->start, sv->end, sv->id));
	}
	for (const snp_t& snp : benchmark_het_snps) {
		nonref_ivals_by_chr[snp.chr].push_back(Interval<std::string, hts_pos_t>(snp.pos, snp.pos, snp.unique_key()));
	}

	std::unordered_map<std::string, IntervalTree<std::string, hts_pos_t>*> nonref_itrees;
	for (auto& chr_windows : nonref_ivals_by_chr) {
		nonref_itrees[chr_windows.first] = new IntervalTree<std::string, hts_pos_t>(chr_windows.second);
	}

	std::unordered_set<std::string> het_with_other_nonref;
	for (const std::shared_ptr<sv_t>& sv : svs) {
		if (!is_known_het_alt(sv.get())) continue;

		auto nonref_it = nonref_itrees.find(sv->chr);
		if (nonref_it == nonref_itrees.end()) continue;

		for (const Interval<std::string, hts_pos_t>& iv : nonref_it->second->findOverlapping(sv->start, sv->end)) {
			if (iv.value != sv->id) {
				het_with_other_nonref.insert(sv->id);
				break;
			}
		}
	}

	for (auto& chr_tree : nonref_itrees) delete chr_tree.second;
	return het_with_other_nonref;
}

std::unordered_map<std::string, int> make_benchmark_match_capacity(const std::vector<std::shared_ptr<sv_t>>& svs, std::vector<snp_t>& benchmark_het_snps, int aux_radius) {

	std::unordered_map<std::string, int> benchmark_match_capacity;
	for (const std::shared_ptr<sv_t>& sv : svs) {
		benchmark_match_capacity[sv->id] = 1;
	}
	if (aux_radius <= 0) return benchmark_match_capacity;

	std::unordered_map<std::string, std::vector<Interval<bool, hts_pos_t>>> het_ivals_by_chr;
	for (const std::shared_ptr<sv_t>& sv : svs) {
		if (!is_known_het_alt(sv.get())) continue;

		het_ivals_by_chr[sv->chr].push_back(Interval<bool, hts_pos_t>(sv->start-aux_radius, sv->end+aux_radius, true));
	}
	for (const snp_t& snp : benchmark_het_snps) {
		het_ivals_by_chr[snp.chr].push_back(Interval<bool, hts_pos_t>(std::max(hts_pos_t(0), snp.pos-aux_radius), snp.pos+aux_radius, true));
	}

	std::unordered_map<std::string, IntervalTree<bool, hts_pos_t>*> het_itrees;
	for (auto& chr_windows : het_ivals_by_chr) {
		het_itrees[chr_windows.first] = new IntervalTree<bool, hts_pos_t>(chr_windows.second);
	}

	for (const std::shared_ptr<sv_t>& sv : svs) {
		if (!is_known_hom_alt(sv.get())) continue;

		auto het_it = het_itrees.find(sv->chr);
		if (het_it != het_itrees.end() && !het_it->second->findOverlapping(sv->start, sv->end).empty()) {
			benchmark_match_capacity[sv->id] = 2;
		}
	}

	for (auto& chr_tree : het_itrees) delete chr_tree.second;
	return benchmark_match_capacity;
}

int main(int argc, char* argv[]) {

	cxxopts::Options options("compare", "Given a VCF file with benchmark deletions and one with the called ones, reports for "
											"each benchmark deletion whether it was called.");

	options.add_options()
		("benchmark_file", "VCF file with the benchmark calls.", cxxopts::value<std::string>())
		("called_file", "VCF file with the calls to benchmark.", cxxopts::value<std::string>())
		("R,reference", "Reference file in fasta format.", cxxopts::value<std::string>())
		("T,tandem-repeats", "Tandem repeat file in BED format.", cxxopts::value<std::string>())
		("d,max_dist_precise", "Maximum distance allowed between the breakpoints of two precise variants.",
				cxxopts::value<int>()->default_value("100"))
		("D,max_dist_imprecise", "Maximum distance allowed between the breakpoints when at least one of the variants is imprecise.",
				cxxopts::value<int>()->default_value("500"))
		("v,min_overlap_precise", "Minimum overlap required between two precise variants, expressed as the fraction of the shortest variant"
				" that is covered by the longest variant.", cxxopts::value<double>()->default_value("0.8"))
		("V,min_overlap_imprecise", "Minimum overlap required between two variants when at least one is imprecise, "
				"expressed as the fraction of the shortest variant that is covered by the longest variant.",
				cxxopts::value<double>()->default_value("0.5"))
		("s,max_len_diff_precise", "Maximum length difference allowed between two precise variants.",
				cxxopts::value<int>()->default_value("100"))
		("S,max_len_diff_imprecise", "Maximum length difference allowed when at least one variant is imprecise.",
				cxxopts::value<int>()->default_value("500"))
		("l,min_len_ratio_precise", "Minimum length ratio allowed (smallest variant length / largest variant length) between two precise variants.", cxxopts::value<double>()->default_value("0.8"))
		("L,min_len_ratio_imprecise", "Minimum length ratio allowed (smallest variant length / largest variant length) when at least one variant is imprecise.", cxxopts::value<double>()->default_value("0.5"))
		("max-repeat-dist", "Maximum distance between two variant in the same tandem repeat.", cxxopts::value<int>()->default_value("1000"))
		("hom-alt-split-radius", "For called-to-benchmark GT labels, give 1/1 benchmark variants match capacity 2 when a 0/1 benchmark variant is within this distance. Disabled when <= 0.", cxxopts::value<int>()->default_value("0"))
		("r,report", "Print report only", cxxopts::value<bool>()->default_value("false"))
		("f,fps", "Print false positive SVs to file.", cxxopts::value<std::string>())
		("tp", "Print true positive called SVs to file.", cxxopts::value<std::string>())
		("fn", "Print false negative benchmark SVs to file.", cxxopts::value<std::string>())
		("I,ignore-seq", "Only compare insertions by position.", cxxopts::value<bool>()->default_value("false"))
		("a,all-imprecise", "Treat all deletions as imprecise.", cxxopts::value<bool>()->default_value("false"))
		("bdup-ids", "ID of benchmark SVs to be considered duplicatons. Note this only affects the final report, not how SVs are compared.", cxxopts::value<std::string>())
		("cdup-ids", "ID of called SVs to be considered duplicatons. Note this only affects the final report, not how SVs are compared.", cxxopts::value<std::string>())
		("wrong-gts", "Print pairs of matching SVs that have discordant genotypes.", cxxopts::value<std::string>())
		("keep-all-benchmark", "Keep all variants in the benchmark file, even if no alternative allele", cxxopts::value<bool>()->default_value("false"))
		("keep-all-called", "Keep all variants in the called file, even if no alternative allele", cxxopts::value<bool>()->default_value("false"))
		("c,called-to-benchmark-gts", "For each called SV, report benchmark GT and exact-match status.", cxxopts::value<std::string>())
		("e,exclusive", "SV cannot be used in multiple matches.", cxxopts::value<bool>()->default_value("false"))
		("ignore-ft", "Ignore the FT field, if present. Without this option, only SVs with FT=PASS or absent are considered.", cxxopts::value<bool>()->default_value("false"))
		("t,threads", "Number of threads to use.", cxxopts::value<int>()->default_value("1"))
		("h,help", "Print usage");

	options.parse_positional({"benchmark_file", "called_file"});
	options.positional_help("benchmark_file called_file");
	options.show_positional_help();
	auto parsed_args = options.parse(argc, argv);

	if (parsed_args.count("help") || !parsed_args.count("benchmark_file") || !parsed_args.count("called_file")) {
		std::cout << options.help() << std::endl;
		exit(0);
	}

    std::string benchmark_fname = parsed_args["benchmark_file"].as<std::string>();
    benchmark_svs = read_sv_list(benchmark_fname.c_str());
    std::string called_fname = parsed_args["called_file"].as<std::string>();
	std::vector<std::shared_ptr<sv_t>> called_svs = read_sv_list(called_fname.c_str());
	max_prec_dist = parsed_args["max_dist_precise"].as<int>();
	max_imprec_dist = parsed_args["max_dist_imprecise"].as<int>();
	max_repeat_dist = parsed_args["max-repeat-dist"].as<int>();
	min_prec_frac_overlap = parsed_args["min_overlap_precise"].as<double>();
	min_imprec_frac_overlap = parsed_args["min_overlap_imprecise"].as<double>();
	max_prec_len_diff = parsed_args["max_len_diff_precise"].as<int>();
	max_imprec_len_diff = parsed_args["max_len_diff_imprecise"].as<int>();
	min_prec_len_ratio = parsed_args["min_len_ratio_precise"].as<double>();
	min_imprec_len_ratio = parsed_args["min_len_ratio_imprecise"].as<double>();
    bool report = parsed_args["report"].as<bool>();
	bool exclusive = parsed_args["exclusive"].as<bool>();
	int threads = parsed_args["threads"].as<int>();
	ignore_seq = parsed_args["ignore-seq"].as<bool>();
	int hom_alt_split_radius = parsed_args["hom-alt-split-radius"].as<int>();
    if (parsed_args["all-imprecise"].as<bool>()) {
    	max_prec_dist = max_imprec_dist;
    	min_prec_frac_overlap = min_imprec_frac_overlap;
    	max_prec_len_diff = max_imprec_len_diff;
    }

	if (parsed_args.count("bdup-ids")) {
		std::string line;
		std::ifstream dup_f(parsed_args["bdup-ids"].as<std::string>());
		while (getline(dup_f, line)) {
			bdup_ids.insert(line);
		}
	}
	if (parsed_args.count("cdup-ids")) {
		std::string line;
		std::ifstream dup_f(parsed_args["cdup-ids"].as<std::string>());
		while (getline(dup_f, line)) {
			cdup_ids.insert(line);
		}
	}

	if (!parsed_args.count("reference") && !ignore_seq) {
		std::cerr << "If no reference is provided, the --ignore-seq flag must be provided." << std::endl;
		exit(0);
	}

	htsFile* fps_fout = NULL;
	bcf_hdr_t* fps_hdr = NULL;
	if (parsed_args.count("fps")) {
		std::string fps_fname = parsed_args["fps"].as<std::string>();
		fps_fout = hts_open(fps_fname.c_str(), "wz");
		// read header from called file
		htsFile* called_file = bcf_open(called_fname.c_str(), "r");
		fps_hdr = bcf_hdr_read(called_file);
		fps_hdr = bcf_hdr_dup(fps_hdr);
		if (bcf_hdr_write(fps_fout, fps_hdr) != 0) {
			std::cerr << "Error writing header to file " << fps_fname << std::endl;
			exit(1);
		}
		bcf_close(called_file);
	}

	htsFile* tp_fout = NULL;
	bcf_hdr_t* tp_hdr = NULL;
	if (parsed_args.count("tp")) {
		std::string tp_fname = parsed_args["tp"].as<std::string>();
		tp_fout = hts_open(tp_fname.c_str(), "wz");
		// read header from called file
		htsFile* called_file = bcf_open(called_fname.c_str(), "r");
		tp_hdr = bcf_hdr_read(called_file);
		tp_hdr = bcf_hdr_dup(tp_hdr);
		if (bcf_hdr_write(tp_fout, tp_hdr) != 0) {
			std::cerr << "Error writing header to file " << tp_fname << std::endl;
			exit(1);
		}
		bcf_close(called_file);
	}

	htsFile* fn_fout = NULL;
	bcf_hdr_t* fn_hdr = NULL;
	if (parsed_args.count("fn")) {
		std::string fn_fname = parsed_args["fn"].as<std::string>();
		fn_fout = hts_open(fn_fname.c_str(), "wz");
		// read header from benchmark file
		htsFile* benchmark_file = bcf_open(benchmark_fname.c_str(), "r");
		fn_hdr = bcf_hdr_read(benchmark_file);
		fn_hdr = bcf_hdr_dup(fn_hdr);
		if (bcf_hdr_write(fn_fout, fn_hdr) != 0) {
			std::cerr << "Error writing header to file " << fn_fname << std::endl;
			exit(1);
		}
		bcf_close(benchmark_file);
	}

	std::ofstream called_to_benchmark_gts_fout;
	if (parsed_args.count("called-to-benchmark-gts")) {
		std::string called_to_benchmark_gts_fname = parsed_args["called-to-benchmark-gts"].as<std::string>();
		called_to_benchmark_gts_fout.open(called_to_benchmark_gts_fname);
	}

	std::ofstream wrong_gts_fout;
	if (parsed_args.count("wrong-gts")) {
		std::string wrong_gts_fname = parsed_args["wrong-gts"].as<std::string>();
		wrong_gts_fout.open(wrong_gts_fname);
	}

	auto is_unsupported_func = [](const std::shared_ptr<sv_t>& sv) {return sv->svtype() != "DEL" && sv->svtype() != "INS" && sv->svtype() != "DUP" && sv->svtype() != "INV";};
	// erase and count elements from benchmark_svs that are not supported
	int n_unsupported = std::count_if(benchmark_svs.begin(), benchmark_svs.end(), is_unsupported_func);
	if (n_unsupported > 0) {
		std::cerr << "Warning: only SVTYPE=DEL, DUP, INS and INV are supported. " << n_unsupported << " unsupported variants in benchmark file." << std::endl;
	}
	benchmark_svs.erase(std::remove_if(benchmark_svs.begin(), benchmark_svs.end(), is_unsupported_func), benchmark_svs.end());

	// erase and count elements from called_svs that are not supported
	n_unsupported = std::count_if(called_svs.begin(), called_svs.end(), is_unsupported_func);
	if (n_unsupported > 0) {
		std::cerr << "Warning: only SVTYPE=DEL, DUP, INS and INV are supported. " << n_unsupported << " unsupported variants in called file." << std::endl;
	}
	called_svs.erase(std::remove_if(called_svs.begin(), called_svs.end(), is_unsupported_func), called_svs.end());

	if (parsed_args["keep-all-benchmark"].as<bool>()) {
		std::cerr << "Keeping all variants in the benchmark file, even if they have no ALT alleles." << std::endl;
	} else {
		int n_ac_0 = std::count_if(benchmark_svs.begin(), benchmark_svs.end(), [](const std::shared_ptr<sv_t>& sv) {return sv->allele_count(1)+sv->missing_alleles() == 0;});
		if (n_ac_0 > 0) {
			std::cerr << "Excluded " << n_ac_0 << " variants in benchmark file that have no ALT alleles." << std::endl;
		}
		benchmark_svs.erase(std::remove_if(benchmark_svs.begin(), benchmark_svs.end(), [](const std::shared_ptr<sv_t>& sv) {return sv->allele_count(1)+sv->missing_alleles() == 0;}), benchmark_svs.end());
	}

	if (parsed_args["keep-all-called"].as<bool>()) {
		std::cerr << "Keeping all variants in the called file, even if they have no ALT alleles." << std::endl;
	} else {
		int n_ac_0 = std::count_if(called_svs.begin(), called_svs.end(), [](const std::shared_ptr<sv_t>& sv) {return sv->allele_count(1) == 0;});
		if (n_ac_0 > 0) {
			std::cerr << "Excluded " << n_ac_0 << " variants in called file that have no ALT alleles." << std::endl;
		}
		called_svs.erase(std::remove_if(called_svs.begin(), called_svs.end(), [](const std::shared_ptr<sv_t>& sv) {return sv->allele_count(1) == 0;}), called_svs.end());
	}

	// erase and count elements from benchmark_svs that are not PASS
	if (!parsed_args["ignore-ft"].as<bool>()) {
		auto is_not_pass = [](const std::shared_ptr<sv_t>& sv) {return !sv->sample_info.is_pass();};
		int n_not_pass_b = std::count_if(benchmark_svs.begin(), benchmark_svs.end(), is_not_pass);
		int n_not_pass_c = std::count_if(called_svs.begin(), called_svs.end(), is_not_pass);
		if (n_not_pass_b > 0) {
			std::cerr << "Excluded " << n_not_pass_b << " variants in benchmark file that are not PASS." << std::endl;
		}
		if (n_not_pass_c > 0) {
			std::cerr << "Excluded " << n_not_pass_c << " variants in called file that are not PASS." << std::endl;
		}
		benchmark_svs.erase(std::remove_if(benchmark_svs.begin(), benchmark_svs.end(), is_not_pass), benchmark_svs.end());
		called_svs.erase(std::remove_if(called_svs.begin(), called_svs.end(), is_not_pass), called_svs.end());

	}

	std::unordered_set<std::string> bdup_ids = find_dup_ids(benchmark_svs);
	std::unordered_set<std::string> cdup_ids = find_dup_ids(called_svs);
	if (!bdup_ids.empty()) {
		std::cerr << "Warning: found duplicate IDs in benchmark file, assigning new IDs" << std::endl;
		for (std::shared_ptr<sv_t>& sv : benchmark_svs) {
			if (bdup_ids.count(sv->id)) {
				sv->id = sv->unique_key();
			}
		}
	}
	if (!cdup_ids.empty()) {
		std::cerr << "Warning: found duplicate IDs in called file, assigning new IDs" << std::endl;
		for (std::shared_ptr<sv_t>& sv : called_svs) {
			if (cdup_ids.count(sv->id)) {
				sv->id = sv->unique_key();
			}
		}
	}

	std::vector<snp_t> benchmark_het_snps;
	if (hom_alt_split_radius > 0) {
		benchmark_het_snps = read_het_benchmark_snps(benchmark_fname);
	}
	std::unordered_set<std::string> vars_overlapping_other_nonref = find_vars_overlapping_other_nonref(benchmark_svs, benchmark_het_snps);
	std::unordered_map<std::string, int> benchmark_match_capacity = make_benchmark_match_capacity(benchmark_svs, benchmark_het_snps, hom_alt_split_radius);

	if (parsed_args.count("reference")) {
		std::string ref_fname = parsed_args["reference"].as<std::string>();
		chr_seqs.read_fasta_into_map(ref_fname);
	}

    if (parsed_args.count("tandem-repeats")) {
		std::string line;
		std::unordered_map<std::string, std::vector<repeat_t>> reps;
		std::unordered_map<std::string, std::vector<Interval<repeat_t>>> reps_iv;
		std::ifstream rep_f(parsed_args["tandem-repeats"].as<std::string>());
		if (!rep_f.is_open()) {
			std::cerr << "Error opening tandem repeat file " << parsed_args["tandem-repeats"].as<std::string>() << std::endl;
			exit(1);
		}
			while (getline(rep_f, line)) {
			if (line.empty() || line[0] == '#') continue;

			std::stringstream ss(line);
			repeat_t r;
			if (!(ss >> r.chr >> r.start >> r.end)) {
				std::cerr << "Error parsing line in tandem repeat file: " << line << std::endl;
				exit(1);
			}

			reps[r.chr].push_back(r);
			reps_iv[r.chr].push_back(Interval<repeat_t>(r.start, r.end, r));
		}

		for (auto it = reps.begin(); it != reps.end(); it++) {
			reps_i[it->first] = new IntervalTree<repeat_t>(reps_iv[it->first]);
		}
    }

	std::unordered_map<std::string, std::vector<Interval<std::shared_ptr<sv_t>>>> dels_csv_iv, inss_csv_iv, invs_csv_iv;
	for (const std::shared_ptr<sv_t>& sv : called_svs) {
		if (sv->svtype() == "DEL") {
			called_dels_by_chr[sv->chr].push_back(sv);
			dels_csv_iv[sv->chr].push_back(Interval<std::shared_ptr<sv_t>>(sv->start, sv->end, sv));
		} else if (sv->svtype() == "INS" || sv->svtype() == "DUP") {
			called_inss_by_chr[sv->chr].push_back(sv);
			inss_csv_iv[sv->chr].push_back(Interval<std::shared_ptr<sv_t>>(sv->start, sv->end, sv));
		} else if (sv->svtype() == "INV") {
			called_invs_by_chr[sv->chr].push_back(sv);
			invs_csv_iv[sv->chr].push_back(Interval<std::shared_ptr<sv_t>>(sv->start, sv->end, sv));
		}
	}
	for (auto it = dels_csv_iv.begin(); it != dels_csv_iv.end(); it++) {
		dels_csv_i[it->first] = new IntervalTree<std::shared_ptr<sv_t>>(it->second);
	}
	for (auto it = inss_csv_iv.begin(); it != inss_csv_iv.end(); it++) {
		inss_csv_i[it->first] = new IntervalTree<std::shared_ptr<sv_t>>(it->second);
	}
	for (auto it = invs_csv_iv.begin(); it != invs_csv_iv.end(); it++) {
		invs_csv_i[it->first] = new IntervalTree<std::shared_ptr<sv_t>>(it->second);
	}

	ctpl::thread_pool thread_pool(threads);
    std::vector<std::future<void> > futures;
	int BLOCK_SIZE = 1;
	for (int i = 0; i < benchmark_svs.size(); i+=BLOCK_SIZE) {
		std::shared_ptr<sv_t> bsv = benchmark_svs[i];
		std::future<void> future = thread_pool.push(find_match, i, std::min(i+BLOCK_SIZE, (int) benchmark_svs.size()));
		futures.push_back(std::move(future));
	}
	thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
		futures[i].get();
    }
    futures.clear();

	std::set<std::string> b_tps, c_tps, b_gt_tps, c_gt_tps, c_unknown;
	std::unordered_map<std::string, sv_match_t*> unknown_match_by_called;

	auto match_for_called_cmp = [](const sv_match_t& a, const sv_match_t& b) {
		if (a.exact != b.exact) return a.exact;
		if (a.b_sv->missing_alleles() != b.b_sv->missing_alleles())
			return a.b_sv->missing_alleles() < b.b_sv->missing_alleles();
		if (a.rep != b.rep) return !a.rep;
		if (a.score != b.score) return a.score > b.score;
		return false;
	};

	// sort matches by
	// 1. exact matches first
	// 2. missing gt in the benchmark (less missing alleles first)
	// 3. rep (false first)
	// 4. score in descending order
	// 5. consistent ALT reads, in descending order
	// 6. ALT reads, in descending order
	// 7. extended ALT consensus score, in descending order
	// 8. benchmark id, called id (to make sort stable)
	std::vector<sv_match_t> accepted_matches;
	std::sort(matches.begin(), matches.end(), [&](sv_match_t& a, sv_match_t& b) {
		if (match_for_called_cmp(a, b)) return true;
		if (match_for_called_cmp(b, a)) return false;
		int a_consistent_reads = a.c_sv->sample_info.alt_bp1.reads_info.consistent_reads() + a.c_sv->sample_info.alt_bp2.reads_info.consistent_reads();
		int b_consistent_reads = b.c_sv->sample_info.alt_bp1.reads_info.consistent_reads() + b.c_sv->sample_info.alt_bp2.reads_info.consistent_reads();
		if (a_consistent_reads != b_consistent_reads) {
			return a_consistent_reads > b_consistent_reads;
		}
		int a_alt_reads = a.c_sv->sample_info.alt_bp1.reads_info.reads + a.c_sv->sample_info.alt_bp2.reads_info.reads;
		int b_alt_reads = b.c_sv->sample_info.alt_bp1.reads_info.reads + b.c_sv->sample_info.alt_bp2.reads_info.reads;
		if (a_alt_reads != b_alt_reads) {
			return a_alt_reads > b_alt_reads;
		}
		int a_ext_alt_consensus = a.c_sv->sample_info.ext_alt_consensus1_to_alt_score + a.c_sv->sample_info.ext_alt_consensus2_to_alt_score;
		int b_ext_alt_consensus = b.c_sv->sample_info.ext_alt_consensus1_to_alt_score + b.c_sv->sample_info.ext_alt_consensus2_to_alt_score;
		if (a_ext_alt_consensus != b_ext_alt_consensus) {
			return a_ext_alt_consensus > b_ext_alt_consensus;
		}
		return std::tie(a.b_sv->id, a.c_sv->id) < std::tie(b.b_sv->id, b.c_sv->id); // to make sort stable
	});

	// for each called SV, we only keep optimal matches, as defined by match_for_called_cmp
	std::vector<sv_match_t> optimal_matches_for_called;
	std::unordered_map<std::string, sv_match_t*> optimal_match_by_called; // we need one example of optimal match for each called SV
	for (sv_match_t& match : matches) {
		std::string c_id = match.c_sv->id;
		if (!optimal_match_by_called.count(c_id)) {
			optimal_match_by_called[c_id] = &match;
			optimal_matches_for_called.push_back(match);
		} else {
			sv_match_t& optimal_match = *optimal_match_by_called[c_id];
			if (!match_for_called_cmp(optimal_match, match)) {
				optimal_matches_for_called.push_back(match);
			}
		}
	}
	matches.swap(optimal_matches_for_called);

	// count tps (both in benchmark and called) from matches
	// in exclusive mode, each sv can only be used in one match
	// for matches where the benchmark sv has missing alleles, since the benchmark does not tell us whether the variant
	// is present or not, we do count the called SV it as neither a tp nor a fp, based on this match
	// (if the called SV is matched to another benchmark SV with known genotype, it will be counted based on that match)
	std::unordered_map<std::string, int> benchmark_match_count;
	std::unordered_set<std::string> used_c_sv_ids;
	for (sv_match_t& match : matches) {
		int b_capacity = benchmark_match_capacity.at(match.b_sv->id);
		if (exclusive && (benchmark_match_count[match.b_sv->id] >= b_capacity || used_c_sv_ids.count(match.c_sv->id))) continue;

		if (b_capacity > 1 && benchmark_match_count[match.b_sv->id] > 0 &&
			match.c_sv->sample_info.alt_bp1.reads_info.reads + match.c_sv->sample_info.alt_bp2.reads_info.reads == 0) {
			match.valid = false;
			continue;
		}

		if (match.b_sv->allele_count(1) == 0 && match.b_sv->missing_alleles() > 0) {
			c_unknown.insert(match.c_sv->id);
			unknown_match_by_called[match.c_sv->id] = &match;
			benchmark_match_count[match.b_sv->id]++;
			used_c_sv_ids.insert(match.c_sv->id);
			continue;
		}

		b_tps.insert(match.b_sv->id);
		benchmark_match_count[match.b_sv->id]++;
		c_tps.insert(match.c_sv->id);
		used_c_sv_ids.insert(match.c_sv->id);
		accepted_matches.push_back(match);
		if (compatible_gts(match.b_sv.get(), match.c_sv.get())) {
			b_gt_tps.insert(match.b_sv->id);
			c_gt_tps.insert(match.c_sv->id);
		} else if (wrong_gts_fout.is_open()) {
			wrong_gts_fout << match.b_sv->id << " " << match.c_sv->id << std::endl;
		}
	}

	StripedSmithWaterman::Aligner aligner(1,4,6,1,false);
	for (const std::shared_ptr<sv_t>& bsv : benchmark_svs) {
		if (!b_tps.count(bsv->id)) {
			accepted_matches.push_back(sv_match_t(bsv, nullptr, false, aligner));
		}
	}

	std::sort(accepted_matches.begin(), accepted_matches.end(), [](sv_match_t& a, sv_match_t& b) {
		return a.b_sv->id < b.b_sv->id;
	});
	if (!report) {
		for (sv_match_t& match : accepted_matches) {
			std::cout << match.b_sv->id << " " << (match.c_sv == NULL ? "NONE" : match.c_sv->id);
			std::cout << (match.rep ? " REP" : "") << std::endl;
		}
	}
	if (called_to_benchmark_gts_fout.is_open()) {
		std::unordered_set<std::string> written_ids;
		auto called_to_benchmark_gt = [&](sv_match_t& match) {
			int b_capacity = benchmark_match_capacity.at(match.b_sv->id);
			if (b_capacity > 1 && benchmark_match_count[match.b_sv->id] > 1) {
				return std::string("1/2");
			}
			std::string gt = match.b_sv->print_gt();
			if (gt == "0/1" && vars_overlapping_other_nonref.count(match.b_sv->id)) {
				return std::string("1/2");
			}
			return gt;
		};

		// Selected non-missing benchmark matches.
		for (sv_match_t& match : accepted_matches) {
			if (match.c_sv != NULL) {
				called_to_benchmark_gts_fout << match.c_sv->id << " " << called_to_benchmark_gt(match) << " " << match.exact << std::endl;
				written_ids.insert(match.c_sv->id);
			}
		}
		std::set<std::string> c_unchoosen_ids;
		std::unordered_map<std::string, sv_match_t*> unchoosen_match_by_called;
		for (sv_match_t& match : matches) {
			if (match.valid && match.c_sv != NULL && !used_c_sv_ids.count(match.c_sv->id)) {
				std::string c_id = match.c_sv->id;
				if (!c_unchoosen_ids.count(c_id)) {
					c_unchoosen_ids.insert(c_id);
					unchoosen_match_by_called[c_id] = &match;
				}
			}
		}

		// these SVs matched to a benchmark SV, but another match for the benchmark SV was chosen instead
		// report non-exact matches with the benchmark GT, but keep exact unselected matches ignored
		for (std::string id : c_unchoosen_ids) {
			sv_match_t* match = unchoosen_match_by_called[id];
			called_to_benchmark_gts_fout << id << " " << "./. " << match->exact << std::endl;
			written_ids.insert(id);
		}

		// Not yet written SVs are either selected missing-GT matches or negatives.
		for (const std::shared_ptr<sv_t>& csv : called_svs) {
			if (written_ids.count(csv->id)) continue;
			if (c_unknown.count(csv->id)) {
				called_to_benchmark_gts_fout << csv->id << " " << "./. " << unknown_match_by_called[csv->id]->exact << std::endl;
			} else {
				called_to_benchmark_gts_fout << csv->id << " " << "0/0 0" << std::endl;
			}
		}
	}

	// count tp and fn benchmark calls, by sv type
	std::unordered_map<std::string, int> n_benchmark_tp, n_benchmark_fn, n_benchmark_gt_tp, n_benchmark_gt_fn;
	for (const std::shared_ptr<sv_t>& bsv : benchmark_svs) {
		if (bsv->allele_count(1) == 0) continue;
		if (b_tps.count(bsv->id)) {
			n_benchmark_tp[get_bsv_type(bsv.get())]++;
		} else {
			n_benchmark_fn[get_bsv_type(bsv.get())]++;
		}
		if (b_gt_tps.count(bsv->id)) {
			n_benchmark_gt_tp[get_bsv_type(bsv.get())]++;
		} else {
			n_benchmark_gt_fn[get_bsv_type(bsv.get())]++;
		}
	}

	std::unordered_map<std::string, int> n_called_tp, n_called_fp, n_called_gt_tp, n_called_gt_fp;
	for (const std::shared_ptr<sv_t>& csv : called_svs) {
		if (c_tps.count(csv->id)) {
			n_called_tp[get_csv_type(csv.get())]++;
		} else {
			if (c_unknown.count(csv->id)) continue;
			n_called_fp[get_csv_type(csv.get())]++;
		}
		if (c_gt_tps.count(csv->id)) {
			n_called_gt_tp[get_csv_type(csv.get())]++;
		} else {
			n_called_gt_fp[get_csv_type(csv.get())]++;
		}
	}

	if (report) {
		std::cout.precision(2);
		for (std::string svtype : {"DEL", "DUP", "INS", "INV"}) {
			int tp = n_benchmark_tp[svtype], fn = n_benchmark_fn[svtype];
			std::cout << svtype << " SENSITIVITY: " << tp << "/" << (tp+fn) << " = " << tp/std::max(1.0, double(tp+fn)) << " ";

			tp = n_benchmark_gt_tp[svtype]; fn = n_benchmark_gt_fn[svtype];
			std::cout << "(considering GT: " << tp << "/" << (tp+fn) << " = " << tp/std::max(1.0, double(tp+fn)) << ")" << std::endl;
			
			tp = n_called_tp[svtype]; int fp = n_called_fp[svtype];
			std::cout << svtype << " PRECISION: " << tp << "/" << (tp+fp) << " = " << tp/std::max(1.0, double(tp+fp)) << " ";

			tp = n_called_gt_tp[svtype]; fp = n_called_gt_fp[svtype];
			std::cout << "(considering GT: " << tp << "/" << (tp+fp) << " = " << tp/std::max(1.0, double(tp+fp)) << ")" << std::endl;
		}
	}

	if (fps_fout != NULL) {
		for (const std::shared_ptr<sv_t>& csv : called_svs) {
			if (!c_tps.count(csv->id) && !c_unknown.count(csv->id)) {
				if (bcf_write(fps_fout, fps_hdr, csv->vcf_entry) != 0) {
					std::cerr << "Error writing variant to file." << std::endl;
					exit(1);
				}
			}
		}
		bcf_close(fps_fout);
		bcf_hdr_destroy(fps_hdr);
	}
	if (tp_fout != NULL) {
		for (const std::shared_ptr<sv_t>& csv : called_svs) {
			if (c_tps.count(csv->id)) {
				if (bcf_write(tp_fout, tp_hdr, csv->vcf_entry) != 0) {
					std::cerr << "Error writing variant to file." << std::endl;
					exit(1);
				}
			}
		}
		bcf_close(tp_fout);
		bcf_hdr_destroy(tp_hdr);
	}
	if (fn_fout != NULL) {
		for (const std::shared_ptr<sv_t>& bsv : benchmark_svs) {
			if (bsv->allele_count(1) > 0 && !b_tps.count(bsv->id)) {
				if (bcf_write(fn_fout, fn_hdr, bsv->vcf_entry) != 0) {
					std::cerr << "Error writing variant to file." << std::endl;
					exit(1);
				}
			}
		}
		bcf_close(fn_fout);
		bcf_hdr_destroy(fn_hdr);
	}
}
