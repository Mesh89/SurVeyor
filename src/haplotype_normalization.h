#ifndef HAPLOTYPE_NORMALIZATION_H
#define HAPLOTYPE_NORMALIZATION_H

#include <map>
#include "var_utils.h"
#include "../libs/ssw.h"
#include "../libs/ssw_cpp.h"

namespace haplotype_normalization {

struct normalization_context_t {
	char* chr_seq;
	hts_pos_t chr_len, extend;
	normalization_context_t(char* chr_seq, hts_pos_t chr_len, hts_pos_t extend) : chr_seq(chr_seq), chr_len(chr_len), extend(extend) {}
};

inline StripedSmithWaterman::Aligner& normalization_aligner() {
	static thread_local StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
	return aligner;
}

inline bool normalization_clipped(const StripedSmithWaterman::Alignment& aln, bool left) {
	if (aln.cigar.empty()) return false;
	uint32_t op = left ? aln.cigar.front() : aln.cigar.back();
	return cigar_int_to_op(op) == 'S' && cigar_int_to_len(op) > 0;
}

// Optional, in-memory provenance. Missing entries mean that a base cannot be
// traced to the original consensus; they never change normalization decisions.
struct normalization_provenance_t {
	const std::string* query = nullptr;
	std::map<hts_pos_t, hts_pos_t> reference;
	std::unordered_map<const sv_t*, std::vector<hts_pos_t>> inserted;
	std::map<std::pair<hts_pos_t, char>, hts_pos_t> snps;

	void add_reference(hts_pos_t ref_pos, hts_pos_t query_pos) {
		auto result = reference.emplace(ref_pos, query_pos);
		if (!result.second && result.first->second != query_pos) result.first->second = -1;
	}

	void add_insertion(const sv_t* sv, hts_pos_t query_begin) {
		auto& positions = inserted[sv];
		positions.clear();
		for (size_t i = 0; i < sv->ins_seq.length(); i++) positions.push_back(query_begin+i);
	}

	hts_pos_t insertion_origin(const sv_t* sv, size_t i) const {
		auto it = inserted.find(sv);
		return it == inserted.end() || i >= it->second.size() ? -1 : it->second[i];
	}

	void add_snp(hts_pos_t pos, char base, hts_pos_t query_pos) {
		auto result = snps.emplace(std::make_pair(pos, base), query_pos);
		if (!result.second && result.first->second != query_pos) result.first->second = -1;
	}

	hts_pos_t snp_origin(const snp_t& snp) const {
		auto it = snps.find({snp.pos, snp.alt_base});
		return it == snps.end() ? -1 : it->second;
	}

	hts_pos_t reference_origin(hts_pos_t pos, char base) const {
		auto snp = snps.find({pos, base});
		auto ref = reference.find(pos);
		hts_pos_t origin = snp != snps.end() ? snp->second : (ref != reference.end() ? ref->second : -1);
		if (query != nullptr && (origin < 0 || origin >= hts_pos_t(query->length()) || toupper((*query)[origin]) != toupper(base))) return -1;
		return origin;
	}
};

inline std::vector<hts_pos_t> query_origin_slice(const std::vector<hts_pos_t>* positions, size_t begin, size_t length) {
	std::vector<hts_pos_t> result(length, -1);
	if (positions != nullptr && begin <= positions->size() && length <= positions->size()-begin) std::copy(positions->begin()+begin, positions->begin()+begin+length, result.begin());
	return result;
}

inline std::vector<hts_pos_t> insertion_origin_slice(normalization_provenance_t* provenance, const sv_t* sv, size_t begin, size_t length) {
	if (provenance == nullptr) return {};
	auto it = provenance->inserted.find(sv);
	return query_origin_slice(it == provenance->inserted.end() ? nullptr : &it->second, begin, length);
}

// Use the same flank builders as normalization, including their trimming and
// edit ordering. Their ALT offsets locate inserted bases without sequence searches.
inline std::vector<hts_pos_t> flank_query_positions(const char* sequence, const std::vector<allele_base_mapping_t>& mapping,
	const std::vector<allele_edit_t>& edits, const std::vector<std::shared_ptr<sv_t>>& indels, const normalization_provenance_t& provenance) {
	std::vector<hts_pos_t> positions(mapping.size(), -1);
	for (size_t i = 0; i < mapping.size(); i++) {
		if (mapping[i].pos >= 0 && !mapping[i].reverse) positions[i] = provenance.reference_origin(mapping[i].pos, sequence[i]);
	}
	for (const auto& edit : edits) {
		if (edit.kind != allele_edit_kind_t::INDEL || edit.alt_begin == edit.alt_end) continue;
		if (edit.alt_begin < 0 || edit.alt_end > int(positions.size())) continue;
		const sv_t* source = nullptr;
		bool ambiguous = false;
		for (const auto& indel : indels) {
			if (indel->start+1 != edit.ref_begin || indel->end+1 != edit.ref_end || indel->ins_seq.length() != size_t(edit.alt_end-edit.alt_begin)) continue;
			if (source != nullptr) ambiguous = true;
			source = indel.get();
		}
		if (source == nullptr || ambiguous) continue;
		for (int i = edit.alt_begin; i < edit.alt_end; i++) positions[i] = provenance.insertion_origin(source, i-edit.alt_begin);
	}
	return positions;
}

// Find the most similar substring of 'text' to 'word' using a simple sliding window approach,
// returning the starting index of the most similar substring.
// Prefer, if possible, strings that are either the prefix or the suffix of 'text',
// as this can lead to a more parsimonious atomization
int find_most_similar_substring(char* text, int text_len, char* word) {
	int word_len = strlen(word);
	int best_start = -1;
	int best_score = -1;

	for (int i = 0; i <= text_len - word_len; i++) {
		int score = 0;
		for (int j = 0; j < word_len; j++) {
			if (toupper(text[i + j]) == toupper(word[j])) score++;
		}

		bool is_prefix_or_suffix = (i == 0 || i == text_len - word_len);
		if (score > best_score || (score == best_score && is_prefix_or_suffix)) {
			best_score = score;
			best_start = i;
		}
	}

	return best_start;
}

std::shared_ptr<sv_t> atomize_del(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {
	if (sv->ins_seq.empty()) return sv;

	char* chr_seq = context.chr_seq;

	if (sv->ins_seq.length() == 1) { // DEL + SNP. Mark last base as SNP and shorten deletion by 1 bp
		hts_pos_t snp_pos = sv->end;
		char ref_base = chr_seq[snp_pos];
		char alt_base = sv->ins_seq[0];
		if (toupper(ref_base) != toupper(alt_base) && is_genomic_base(alt_base)) {
			snp_t snp(snp_pos, alt_base);
			sv->aux_snps.push_back(snp);
			if (provenance != nullptr) provenance->add_snp(snp_pos, alt_base, provenance->insertion_origin(sv.get(), 0));
		}

		sv->end--;
		sv->ins_seq = "";
	} else {
		// transform DEL + INS into (possibly two) DELs + SNPs (as few as possible)
		int most_similar_pos = find_most_similar_substring(chr_seq + sv->start+1, sv->end - sv->start, (char*) sv->ins_seq.c_str());
		if (most_similar_pos == -1) {
			// inserted sequence is longer than the deleted sequence (probably due to a malformed variant)
			return sv;
		}

		// Create two new SVs for the split deletions
		hts_pos_t sv1_start = sv->start;
		hts_pos_t sv1_end = sv->start + most_similar_pos;
		hts_pos_t sv2_start = sv->start + most_similar_pos + sv->ins_seq.length();
		hts_pos_t sv2_end = sv->end;
		if (sv1_end-sv1_start >= sv2_end-sv2_start) {
			sv->end = sv1_end;
			if (sv2_end > sv2_start) {
				std::shared_ptr<sv_t> sv2 = std::make_shared<deletion_t>(sv->chr, sv2_start, sv2_end, "", nullptr, nullptr, nullptr, nullptr);
				sv->aux_indels.push_back(sv2);
			}
		} else {
			sv->start = sv2_start;
			if (sv1_end > sv1_start) {
				std::shared_ptr<sv_t> sv1 = std::make_shared<deletion_t>(sv->chr, sv1_start, sv1_end, "", nullptr, nullptr, nullptr, nullptr);
				sv->aux_indels.push_back(sv1);
			}
		}

		// Create SNPs for the mismatched bases
		for (size_t i = 0; i < sv->ins_seq.length(); i++) {
			char ref_base = chr_seq[sv1_start+1 + most_similar_pos + i];
			char alt_base = sv->ins_seq[i];
			if (toupper(ref_base) != toupper(alt_base) && is_genomic_base(alt_base)) {
				hts_pos_t snp_pos = sv1_start+1 + most_similar_pos + i;
				sv->aux_snps.push_back(snp_t(snp_pos, alt_base));
				if (provenance != nullptr) provenance->add_snp(snp_pos, alt_base, provenance->insertion_origin(sv.get(), i));
			}
		}

		std::sort(sv->aux_snps.begin(), sv->aux_snps.end(),
			[](const snp_t& a, const snp_t& b) {
				return a.pos < b.pos;
			});

		sv->ins_seq = "";
	}

	if (sv->start >= sv->end) {
		return nullptr;
	}
	return sv;
}

std::shared_ptr<sv_t> atomize_ins(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {
	if (sv->start == sv->end) return sv;
	std::vector<hts_pos_t> original_positions = insertion_origin_slice(provenance, sv.get(), 0, sv->ins_seq.length());

	char* chr_seq = context.chr_seq;

	char* deleted_seq = new char[sv->end - sv->start + 1];
	strncpy(deleted_seq, chr_seq + sv->start+1, sv->end - sv->start);
	deleted_seq[sv->end - sv->start] = '\0';

	// Find where the replaced reference sequence best fits within the inserted sequence.
	int most_similar_pos = find_most_similar_substring((char*) sv->ins_seq.c_str(), sv->ins_seq.length(), deleted_seq);
	if (most_similar_pos == -1) {
		// Deleted sequence is longer than the inserted sequence; leave the variant unchanged for now.
		delete[] deleted_seq;
		return sv;
	}

	int deleted_seq_len = sv->end - sv->start;
	hts_pos_t left_ins_pos = sv->start;
	hts_pos_t right_ins_pos = sv->end;
	hts_pos_t snp_start = sv->start+1;
	std::string left_ins = sv->ins_seq.substr(0, most_similar_pos);
	std::string matched_ins = sv->ins_seq.substr(most_similar_pos, deleted_seq_len);
	std::string right_ins = sv->ins_seq.substr(most_similar_pos + deleted_seq_len);

	if (left_ins.length() >= right_ins.length()) {
		sv->end = left_ins_pos;
		sv->ins_seq = left_ins;
		if (provenance != nullptr) provenance->inserted[sv.get()] = query_origin_slice(&original_positions, 0, left_ins.length());
		if (!right_ins.empty()) {
			std::shared_ptr<sv_t> right_ins_sv = std::make_shared<insertion_t>(sv->chr, right_ins_pos, right_ins_pos, right_ins, nullptr, nullptr, nullptr, nullptr);
			sv->aux_indels.push_back(right_ins_sv);
			if (provenance != nullptr) provenance->inserted[right_ins_sv.get()] = query_origin_slice(&original_positions, most_similar_pos+deleted_seq_len, right_ins.length());
		}
	} else {
		sv->start = right_ins_pos;
		sv->ins_seq = right_ins;
		if (provenance != nullptr) provenance->inserted[sv.get()] = query_origin_slice(&original_positions, most_similar_pos+deleted_seq_len, right_ins.length());
		if (!left_ins.empty()) {
			std::shared_ptr<sv_t> left_ins_sv = std::make_shared<insertion_t>(sv->chr, left_ins_pos, left_ins_pos, left_ins, nullptr, nullptr, nullptr, nullptr);
			sv->aux_indels.push_back(left_ins_sv);
			if (provenance != nullptr) provenance->inserted[left_ins_sv.get()] = query_origin_slice(&original_positions, 0, left_ins.length());
		}
	}

	// Create SNPs for the mismatched bases in the reference-like middle.
	for (int i = 0; i < deleted_seq_len; i++) {
		char ref_base = chr_seq[snp_start + i];
		char alt_base = matched_ins[i];
		if (toupper(ref_base) != toupper(alt_base)) {
			hts_pos_t snp_pos = snp_start + i;
			sv->aux_snps.push_back(snp_t(snp_pos, alt_base));
			if (provenance != nullptr) provenance->add_snp(snp_pos, alt_base, original_positions[most_similar_pos+i]);
		}
	}

	std::sort(sv->aux_snps.begin(), sv->aux_snps.end(),
		[](const snp_t& a, const snp_t& b) {
			return a.pos < b.pos;
		});

	if (sv->ins_seq.empty() && sv->aux_indels.empty()) {
		delete[] deleted_seq;
		return nullptr;
	}

	delete[] deleted_seq;
	return sv;
}

bool is_complex_indel(const std::shared_ptr<sv_t>& sv) {
	return sv->svtype() == "RPL" || (sv->svtype() == "DEL" && !sv->ins_seq.empty()) ||
		(sv->svtype() == "INS" && sv->start != sv->end);
}

std::shared_ptr<sv_t> replace_main_indel(std::shared_ptr<sv_t> old_main, std::shared_ptr<sv_t> new_main) {
	auto old_del = std::dynamic_pointer_cast<deletion_t>(old_main);
	auto new_del = std::dynamic_pointer_cast<deletion_t>(new_main);
	if (old_del && new_del) {
		new_del->remapped = old_del->remapped;
		new_del->original_range = old_del->original_range;
	}
	std::string chr = new_main->chr, ins_seq = new_main->ins_seq;
	hts_pos_t start = new_main->start, end = new_main->end;
	bcf1_t* vcf_entry = old_main->vcf_entry;
	old_main->vcf_entry = nullptr;
	static_cast<sv_t&>(*new_main) = static_cast<const sv_t&>(*old_main);
	new_main->chr = chr;
	new_main->start = start;
	new_main->end = end;
	new_main->ins_seq = ins_seq;
	new_main->vcf_entry = vcf_entry;
	return new_main;
}

std::shared_ptr<sv_t> merge_adjacent_main_indel(std::shared_ptr<sv_t> main_sv, normalization_provenance_t* provenance = nullptr) {
	if (is_complex_indel(main_sv) || (main_sv->svtype() != "DEL" && main_sv->svtype() != "INS")) return main_sv;
	for (auto it = main_sv->aux_indels.begin(); it != main_sv->aux_indels.end(); ++it) {
		if ((*it)->svtype() == main_sv->svtype() || (*it)->incomplete_ins_seq()) continue;
		auto del = main_sv->svtype() == "DEL" ? main_sv : *it;
		auto ins = main_sv->svtype() == "INS" ? main_sv : *it;
		if (del->svtype() != "DEL" || ins->svtype() != "INS") continue;
		if (ins->start != del->start && ins->start != del->end) continue;

		std::shared_ptr<sv_t> complex;
		if (del->svsize() > ins->svsize()) complex = std::make_shared<deletion_t>(del->chr, del->start, del->end, ins->ins_seq, nullptr, nullptr, nullptr, nullptr);
		else if (ins->svsize() > del->svsize()) complex = std::make_shared<insertion_t>(del->chr, del->start, del->end, ins->ins_seq, nullptr, nullptr, nullptr, nullptr);
		else complex = std::make_shared<replacement_t>(del->chr, del->start, del->end, ins->ins_seq, nullptr, nullptr, nullptr, nullptr);
		if (provenance != nullptr) provenance->inserted[complex.get()] = insertion_origin_slice(provenance, ins.get(), 0, ins->ins_seq.length());
		main_sv->aux_indels.erase(it);
		return replace_main_indel(main_sv, complex);
	}
	return main_sv;
}

std::shared_ptr<sv_t> atomize_limited(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {
	if (sv->svtype() == "DEL" && sv->svsize() <= 70) return atomize_del(sv, context, provenance);
	if (sv->svtype() == "INS" && !sv->incomplete_ins_seq() && sv->ins_seq.length() <= 70) return atomize_ins(sv, context, provenance);
	return sv;
}

std::shared_ptr<sv_t> atomize_haplotype(std::shared_ptr<sv_t> main_sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {
	main_sv = atomize_limited(main_sv, context, provenance);
	if (main_sv == nullptr) return nullptr;
	std::vector<std::shared_ptr<sv_t>> input_aux, output_aux;
	input_aux.swap(main_sv->aux_indels);
	for (auto& aux : input_aux) {
		aux = atomize_limited(aux, context, provenance);
		if (aux == nullptr) continue;
		main_sv->aux_snps.insert(main_sv->aux_snps.end(), aux->aux_snps.begin(), aux->aux_snps.end());
		output_aux.push_back(aux);
		output_aux.insert(output_aux.end(), aux->aux_indels.begin(), aux->aux_indels.end());
		aux->aux_indels.clear();
		aux->aux_snps.clear();
	}
	main_sv->aux_indels = std::move(output_aux);
	for (auto& aux : main_sv->aux_indels) aux->hpid = main_sv->hpid;
	return main_sv;
}

std::shared_ptr<sv_t> simplify_del(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {

	char* chr_seq = context.chr_seq;

	// try to shorten deletion if ins_seq
	int start_is = 0;
	while (start_is < sv->ins_seq.length() && sv->start < sv->end &&
		   toupper(sv->ins_seq[start_is]) == toupper(chr_seq[sv->start+1])) {
		if (provenance != nullptr) provenance->add_reference(sv->start+1, provenance->insertion_origin(sv.get(), start_is));
		start_is++;
		sv->start++;
	}
	if (provenance != nullptr) provenance->inserted[sv.get()] = insertion_origin_slice(provenance, sv.get(), start_is, sv->ins_seq.length()-start_is);
	sv->ins_seq = sv->ins_seq.substr(start_is);

	int end_is = sv->ins_seq.length();
	while (end_is > 0 && sv->start < sv->end &&
		   toupper(sv->ins_seq[end_is-1]) == toupper(chr_seq[sv->end])) {
		if (provenance != nullptr) provenance->add_reference(sv->end, provenance->insertion_origin(sv.get(), end_is-1));
		end_is--;
		sv->end--;
	}
	if (provenance != nullptr) provenance->inserted[sv.get()] = insertion_origin_slice(provenance, sv.get(), 0, end_is);
	sv->ins_seq = sv->ins_seq.substr(0, end_is);

	if (sv->start >= sv->end) {
		return nullptr;
	}
	return sv;
}

std::shared_ptr<sv_t> simplify_ins(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {

	char* chr_seq = context.chr_seq;

	// try to shorten insertion if start != end
	int start_is = 0;
	while (start_is < sv->ins_seq.length() && sv->start < sv->end && toupper(sv->ins_seq[start_is]) == toupper(chr_seq[sv->start+1])) {
		if (provenance != nullptr) provenance->add_reference(sv->start+1, provenance->insertion_origin(sv.get(), start_is));
		start_is++;
		sv->start++;
	}
	if (provenance != nullptr) provenance->inserted[sv.get()] = insertion_origin_slice(provenance, sv.get(), start_is, sv->ins_seq.length()-start_is);
	sv->ins_seq = sv->ins_seq.substr(start_is);

	int end_is = sv->ins_seq.length();
	while (end_is > 0 && sv->start < sv->end && toupper(sv->ins_seq[end_is-1]) == toupper(chr_seq[sv->end])) {
		if (provenance != nullptr) provenance->add_reference(sv->end, provenance->insertion_origin(sv.get(), end_is-1));
		end_is--;
		sv->end--;
	}
	if (provenance != nullptr) provenance->inserted[sv.get()] = insertion_origin_slice(provenance, sv.get(), 0, end_is);
	sv->ins_seq = sv->ins_seq.substr(0, end_is);

	if (sv->ins_seq.empty()) {
		return nullptr;
	}
	return sv;
}

std::shared_ptr<sv_t> simplify(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {
	std::string svtype = sv->svtype();
	if (svtype == "DEL") {
		return simplify_del(sv, context, provenance);
	} else if (svtype == "INS") {
		return simplify_ins(sv, context, provenance);
	} else {
		return sv;
	}
}

std::vector<std::shared_ptr<sv_t>> vars_from_alignment(StripedSmithWaterman::Alignment& aln, std::string chr, hts_pos_t ref_start, char* alt_seq, std::vector<snp_t>& snps, normalization_provenance_t* provenance = nullptr, const std::vector<hts_pos_t>* alt_query_positions = nullptr) {
	std::vector<std::shared_ptr<sv_t>> vars;
	hts_pos_t current_pos = ref_start + aln.ref_begin, query_pos = aln.query_begin;
	for (int i = 0; i < aln.cigar.size(); i++) {
		uint32_t c = aln.cigar[i];
		int op_length = cigar_int_to_len(c);
		char op = cigar_int_to_op(c);
		char next_op = i+1 < aln.cigar.size() ? cigar_int_to_op(aln.cigar[i+1]) : '\0';
		if ((op == 'D' && next_op == 'I') || (op == 'I' && next_op == 'D')) {
			int next_length = cigar_int_to_len(aln.cigar[++i]);
			int del_length = op == 'D' ? op_length : next_length;
			int ins_length = op == 'I' ? op_length : next_length;
			std::string ins_seq(alt_seq+query_pos, ins_length);
			if (del_length > ins_length) vars.push_back(std::make_shared<deletion_t>(chr, current_pos-1, current_pos+del_length-1, ins_seq, nullptr, nullptr, nullptr, nullptr));
			else if (ins_length > del_length) vars.push_back(std::make_shared<insertion_t>(chr, current_pos-1, current_pos+del_length-1, ins_seq, nullptr, nullptr, nullptr, nullptr));
			else vars.push_back(std::make_shared<replacement_t>(chr, current_pos-1, current_pos+del_length-1, ins_seq, nullptr, nullptr, nullptr, nullptr));
			if (provenance != nullptr) provenance->inserted[vars.back().get()] = query_origin_slice(alt_query_positions, query_pos, ins_length);
			current_pos += del_length;
			query_pos += ins_length;
			continue;
		}
		if (op == 'D') {
			vars.push_back(std::make_shared<deletion_t>(chr, current_pos-1, current_pos+op_length-1, "", nullptr, nullptr, nullptr, nullptr));
		} else if (op == 'I') {
			vars.push_back(std::make_shared<insertion_t>(chr, current_pos-1, current_pos-1, std::string(alt_seq+query_pos, op_length), nullptr, nullptr, nullptr, nullptr));
			if (provenance != nullptr) provenance->inserted[vars.back().get()] = query_origin_slice(alt_query_positions, query_pos, op_length);
		} else if (op == 'X') {
			for (int j = 0; j < op_length; j++) {
				if (is_genomic_base(alt_seq[query_pos+j])) {
					snps.push_back(snp_t(current_pos+j, alt_seq[query_pos+j]));
					if (provenance != nullptr) provenance->add_snp(current_pos+j, alt_seq[query_pos+j], query_origin_slice(alt_query_positions, query_pos+j, 1)[0]);
				}
			}
		}

		if (op != 'I' && op != 'S') current_pos += op_length;
		if (op != 'D') query_pos += op_length;
	}
	return vars;
}

bool realign_haplotype_atoms(std::shared_ptr<sv_t> sv, std::vector<std::shared_ptr<sv_t>>& indels, std::vector<snp_t>& snps, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr, const std::vector<hts_pos_t>* alt_query_positions = nullptr) {
	char* chr_seq = context.chr_seq;
	hts_pos_t chr_len = context.chr_len;
	hts_pos_t ins_start = sv->start, ins_end = sv->end;

	hts_pos_t extend = context.extend;
	std::vector<std::shared_ptr<sv_t>> aux_indels = sv->aux_indels;
	std::vector<snp_t> aux_snps = sv->aux_snps;
	std::vector<allele_edit_t> lf_edits, rf_edits;
	std::vector<allele_base_mapping_t> lf_mapping, rf_mapping;
	char* lf_seq = generate_haplotype_left(chr_seq, ins_start, extend, aux_indels, aux_snps, provenance ? &lf_edits : nullptr, provenance ? &lf_mapping : nullptr);
	char* rf_seq = generate_haplotype_right(chr_seq, chr_len, ins_end+1, extend, aux_indels, aux_snps, provenance ? &rf_edits : nullptr, provenance ? &rf_mapping : nullptr);

	hts_pos_t lf_len = strlen(lf_seq), rf_len = strlen(rf_seq), alt_len = lf_len + sv->ins_seq.length() + rf_len;
	char* putative_alt_allele = new char[alt_len + 1];
	memcpy(putative_alt_allele, lf_seq, lf_len);
	memcpy(putative_alt_allele + lf_len, sv->ins_seq.c_str(), sv->ins_seq.length());
	memcpy(putative_alt_allele + lf_len + sv->ins_seq.length(), rf_seq, rf_len);
	putative_alt_allele[alt_len] = '\0';
	std::vector<hts_pos_t> reconstructed_positions;
	if (provenance != nullptr && alt_query_positions == nullptr) {
		reconstructed_positions = flank_query_positions(lf_seq, lf_mapping, lf_edits, aux_indels, *provenance);
		std::vector<hts_pos_t> inserted_positions = insertion_origin_slice(provenance, sv.get(), 0, sv->ins_seq.length());
		reconstructed_positions.insert(reconstructed_positions.end(), inserted_positions.begin(), inserted_positions.end());
		std::vector<hts_pos_t> right_positions = flank_query_positions(rf_seq, rf_mapping, rf_edits, aux_indels, *provenance);
		reconstructed_positions.insert(reconstructed_positions.end(), right_positions.begin(), right_positions.end());
		alt_query_positions = &reconstructed_positions;
	}

	hts_pos_t ref_start = std::max(hts_pos_t(0), ins_start - extend + 1);
	hts_pos_t ref_end = std::min(chr_len, ins_end + extend + 1);
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	// keep extending the reference window until we find an alignment without soft clipping at either end (up to 10 times)
	for (int i = 0; i <= 10; i++) {
		normalization_aligner().Align(putative_alt_allele, chr_seq + ref_start, ref_end-ref_start, filter, &aln, 0);
		bool left_clipped = normalization_clipped(aln, true);
		bool right_clipped = normalization_clipped(aln, false);
		if ((!left_clipped && !right_clipped) || i == 10) break; // avoid expanding in the last loop

		bool changed = false;
		if (left_clipped && ref_start > 0) {
			ref_start = std::max(hts_pos_t(0), ref_start-50);
			changed = true;
		}
		if (right_clipped && ref_end < chr_len) {
			ref_end = std::min(chr_len, ref_end+50);
			changed = true;
		}
		if (!changed) break;
	}

	bool success = !normalization_clipped(aln, true) && !normalization_clipped(aln, false);
	if (success) {
		if (provenance != nullptr) provenance->snps.clear();
		indels = vars_from_alignment(aln, sv->chr, ref_start, putative_alt_allele, snps, provenance, alt_query_positions);
	}

	delete[] lf_seq;
	delete[] rf_seq;
	delete[] putative_alt_allele;

	return success;
}

std::shared_ptr<sv_t> update_main_var_from_realigned_vars(std::shared_ptr<sv_t> var, std::vector<std::shared_ptr<sv_t>>& vars, std::vector<snp_t>& snps) {
	if (vars.empty()) return nullptr;

	std::shared_ptr<sv_t> primary = vars.front();
	std::string aligned_chr = primary->chr;
	hts_pos_t aligned_start = primary->start;
	hts_pos_t aligned_end = primary->end;
	std::string aligned_ins_seq = primary->ins_seq;

	// Preserve type-specific deletion state when realignment keeps a deletion as
	// the primary variant. There is no corresponding state on insertion_t or RPL.
	std::shared_ptr<deletion_t> original_del = std::dynamic_pointer_cast<deletion_t>(var);
	std::shared_ptr<deletion_t> primary_del = std::dynamic_pointer_cast<deletion_t>(primary);
	if (original_del && primary_del) {
		primary_del->remapped = original_del->remapped;
		primary_del->original_range = original_del->original_range;
	}

	// Copy all common evidence and annotations without changing the dynamic type
	// selected by the alignment. vcf_entry ownership must be transferred because
	// sv_t destroys it in its destructor.
	bcf1_t* original_vcf_entry = var->vcf_entry;
	var->vcf_entry = nullptr;
	static_cast<sv_t&>(*primary) = static_cast<const sv_t&>(*var);
	primary->vcf_entry = original_vcf_entry;

	primary->chr = aligned_chr;
	primary->start = aligned_start;
	primary->end = aligned_end;
	primary->ins_seq = std::move(aligned_ins_seq);
	primary->aux_indels.assign(vars.begin() + 1, vars.end());
	primary->aux_snps = std::move(snps);
	for (const auto& aux_indel : primary->aux_indels) aux_indel->hpid = primary->hpid;

	return primary;
}

std::shared_ptr<sv_t> realign_indel_haplotype(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {
	bool has_replaced_ref = sv->start != sv->end && !sv->ins_seq.empty();
	bool has_aux_indels = !sv->aux_indels.empty();
	if (sv->incomplete_ins_seq() || (!has_replaced_ref && !has_aux_indels) || sv->svsize() > 100) return sv;

	std::vector<std::shared_ptr<sv_t>> realigned_svs;
	std::vector<snp_t> snps;
	if (realign_haplotype_atoms(sv, realigned_svs, snps, context, provenance)) {
		sv = update_main_var_from_realigned_vars(sv, realigned_svs, snps);
	}
	return sv;
}

std::shared_ptr<sv_t> realign(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {
	std::string svtype = sv->svtype();
	if (svtype == "INS" || svtype == "DEL" || svtype == "RPL") {
		return realign_indel_haplotype(sv, context, provenance);
	} else {
		return sv;
	}
}


std::shared_ptr<sv_t> normalize_haplotype(std::shared_ptr<sv_t> sv, const normalization_context_t& context, normalization_provenance_t* provenance = nullptr) {
	sv = merge_adjacent_main_indel(sv, provenance);
	if (sv != nullptr) sv = simplify(sv, context, provenance);
	if (sv != nullptr) {
		sv = realign(sv, context, provenance);
		if (sv != nullptr) sv = atomize_haplotype(sv, context, provenance);
	}
	return sv;
}

} // namespace haplotype_normalization

#endif
