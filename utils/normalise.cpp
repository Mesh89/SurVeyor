#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <memory>
#include <mutex>
#include <functional>

#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "../src/sam_utils.h"
#include "../src/sw_utils.h"
#include "../src/var_utils.h"
#include "../src/vcf_utils.h"
#include "../libs/cptl_stl.h"

chr_seqs_map_t chr_seqs;
bcf_hdr_t* hdr;
thread_local StripedSmithWaterman::Aligner realign_aligner(1, 4, 6, 1, false);
int normalise_max_is = 1000;

std::vector<bcf1_t*> normalised_vcf_records;
std::mutex mtx;

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

std::shared_ptr<sv_t> atomize_del(std::shared_ptr<sv_t> sv) {
	if (sv->ins_seq.empty()) return sv;

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	if (sv->ins_seq.length() == 1) { // DEL + SNP. Mark last base as SNP and shorten deletion by 1 bp
		hts_pos_t snp_pos = sv->end;
		char ref_base = chr_seq[snp_pos];
		char alt_base = sv->ins_seq[0];
		if (toupper(ref_base) != toupper(alt_base) && is_genomic_base(alt_base)) {
			snp_t snp(snp_pos, alt_base);
			sv->aux_snps.push_back(snp);
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

std::shared_ptr<sv_t> atomize_ins(std::shared_ptr<sv_t> sv) {
	if (sv->start == sv->end) return sv;

	char* chr_seq = chr_seqs.get_seq(sv->chr);

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
		if (!right_ins.empty()) {
			std::shared_ptr<sv_t> right_ins_sv = std::make_shared<insertion_t>(sv->chr, right_ins_pos, right_ins_pos, right_ins, nullptr, nullptr, nullptr, nullptr);
			sv->aux_indels.push_back(right_ins_sv);
		}
	} else {
		sv->start = right_ins_pos;
		sv->ins_seq = right_ins;
		if (!left_ins.empty()) {
			std::shared_ptr<sv_t> left_ins_sv = std::make_shared<insertion_t>(sv->chr, left_ins_pos, left_ins_pos, left_ins, nullptr, nullptr, nullptr, nullptr);
			sv->aux_indels.push_back(left_ins_sv);
		}
	}

	// Create SNPs for the mismatched bases in the reference-like middle.
	for (int i = 0; i < deleted_seq_len; i++) {
		char ref_base = chr_seq[snp_start + i];
		char alt_base = matched_ins[i];
		if (toupper(ref_base) != toupper(alt_base)) {
			hts_pos_t snp_pos = snp_start + i;
			sv->aux_snps.push_back(snp_t(snp_pos, alt_base));
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

std::shared_ptr<sv_t> atomize(int id, std::shared_ptr<sv_t> sv) {
	std::string svtype = sv->svtype();
	if (svtype == "DEL" && sv->svsize() <= 70) { 
		// only atomize small deletions
		// splitting large deletions can create problems to the current genotyping algorithm
		// especially when calculating features like discordant pairs or read depth
		return atomize_del(sv);
	} else if (svtype == "INS" && !sv->incomplete_ins_seq() && sv->ins_seq.length() <= 70) {
		return atomize_ins(sv);
	} else {
		return sv;
	}
}

hts_pos_t max_indel_size(const std::shared_ptr<sv_t>& sv) {
	hts_pos_t max_size = sv->svsize();
	for (const auto& aux_indel : sv->aux_indels) {
		max_size = std::max(max_size, aux_indel->svsize());
	}
	return max_size;
}

std::shared_ptr<sv_t> simplify_del(std::shared_ptr<sv_t> sv) {

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	// try to shorten deletion if ins_seq
	int start_is = 0;
	while (start_is < sv->ins_seq.length() && sv->start < sv->end &&
		   toupper(sv->ins_seq[start_is]) == toupper(chr_seq[sv->start+1])) {
		start_is++;
		sv->start++;
	}
	sv->ins_seq = sv->ins_seq.substr(start_is);

	int end_is = sv->ins_seq.length();
	while (end_is > 0 && sv->start < sv->end &&
		   toupper(sv->ins_seq[end_is-1]) == toupper(chr_seq[sv->end])) {
		end_is--;
		sv->end--;
	}
	sv->ins_seq = sv->ins_seq.substr(0, end_is);

	if (sv->start >= sv->end) {
		return nullptr;
	}
	return sv;
}

std::shared_ptr<sv_t> simplify_ins(std::shared_ptr<sv_t> sv) {

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	// try to shorten insertion if start != end
	int start_is = 0;
	while (start_is < sv->ins_seq.length() && sv->start < sv->end && toupper(sv->ins_seq[start_is]) == toupper(chr_seq[sv->start+1])) {
		start_is++;
		sv->start++;
	}
	sv->ins_seq = sv->ins_seq.substr(start_is);

	int end_is = sv->ins_seq.length();
	while (end_is > 0 && sv->start < sv->end && toupper(sv->ins_seq[end_is-1]) == toupper(chr_seq[sv->end])) {
		end_is--;
		sv->end--;
	}
	sv->ins_seq = sv->ins_seq.substr(0, end_is);

	if (sv->ins_seq.empty()) {
		return nullptr;
	}
	return sv;
}

std::shared_ptr<sv_t> simplify(std::shared_ptr<sv_t> sv) {
	std::string svtype = sv->svtype();
	if (svtype == "DEL") {
		return simplify_del(sv);
	} else if (svtype == "INS") {
		return simplify_ins(sv);
	} else {
		return sv;
	}
}

std::shared_ptr<sv_t> realign(std::shared_ptr<sv_t> sv);

void simplify_and_realign_svs(int id, std::vector<std::shared_ptr<sv_t>>& svs, int start, int end) {
	for (int i = start; i < end; i++) {
		svs[i] = simplify(svs[i]);
		if (svs[i] != nullptr) {
			svs[i] = realign(svs[i]);
		}
	}
}

void left_align_del(std::shared_ptr<sv_t> sv) {

	if (!sv->ins_seq.empty()) return;

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	for (snp_t& snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	hts_pos_t limit = 0; // we cannot left-align past this position
	for (std::shared_ptr<sv_t> indel : sv->aux_indels) {
		if (indel->end <= sv->start) {
			limit = std::max(limit, indel->end); 
		}
	}

	while (sv->start > limit && toupper(chr_seq[sv->start]) == toupper(chr_seq[sv->end])) {
		sv->start--;
		sv->end--;
	}

	for (snp_t& snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	// delete snps that fall within the deleted region
	sv->aux_snps.erase(std::remove_if(sv->aux_snps.begin(), sv->aux_snps.end(),
		[sv](const snp_t& snp) { return snp.pos > sv->start && snp.pos <= sv->end; }), sv->aux_snps.end());
}

void left_align_dup(std::shared_ptr<sv_t> sv) {

	if (!sv->ins_seq.empty()) return;
	if (!sv->aux_indels.empty()) return; // this may be complicated, let us skip for now

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	hts_pos_t limit = 0;
	for (snp_t& snp : sv->aux_snps) {
		if (snp.pos <= sv->start) {
			limit = snp.pos; // note that snps are sorted by position
		}
	}

	for (snp_t& snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	while (sv->start > limit && toupper(chr_seq[sv->start]) == toupper(chr_seq[sv->end])) {
		sv->start--;
		sv->end--;
	}

	for (snp_t& snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}
}

void left_align_ins(std::shared_ptr<sv_t> sv) {

	if (sv->ins_seq.empty()) return;
	if (!sv->aux_indels.empty()) return; // this may be complicated, let us skip for now

	char* chr_seq = chr_seqs.get_seq(sv->chr);

	for (snp_t& snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	if (sv->start == sv->end) {
		while (sv->start > 0 && toupper(chr_seq[sv->start]) == toupper(sv->ins_seq[sv->ins_seq.length()-1])) {
			for (int i = sv->ins_seq.length()-1; i >= 1; i--) {
				sv->ins_seq[i] = sv->ins_seq[i-1];
			}
			sv->ins_seq[0] = toupper(chr_seq[sv->start]);
			sv->start--;
			sv->end--;
		}
	}

	for (snp_t& snp : sv->aux_snps) {
		std::swap(chr_seq[snp.pos], snp.alt_base);
	}

	// delete snps that fall within the deleted region
	sv->aux_snps.erase(std::remove_if(sv->aux_snps.begin(), sv->aux_snps.end(),
		[sv](const snp_t& snp) { return snp.pos > sv->start && snp.pos <= sv->end; }), sv->aux_snps.end());
}


void left_align(std::shared_ptr<sv_t> sv) {
	std::string svtype = sv->svtype();
	if (svtype == "DEL") {
		left_align_del(sv);
	} else if (svtype == "DUP") {
		left_align_dup(sv);
	} else if (svtype == "INS") {
		left_align_ins(sv);
	}
}

void canonicalize_aux(std::shared_ptr<sv_t> sv) {
	std::sort(sv->aux_snps.begin(), sv->aux_snps.end(),
		[](const snp_t& a, const snp_t& b) {
			return std::tie(a.pos, a.alt_base) < std::tie(b.pos, b.alt_base);
		});
	std::sort(sv->aux_indels.begin(), sv->aux_indels.end(),
		[](const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) {
			return std::make_tuple(a->start, a->end, a->svtype(), a->ins_seq) <
				   std::make_tuple(b->start, b->end, b->svtype(), b->ins_seq);
		});
}

std::vector<std::shared_ptr<sv_t>> vars_from_alignment(StripedSmithWaterman::Alignment& aln, std::string chr, hts_pos_t ref_start, char* alt_seq, std::vector<snp_t>& snps) {
	std::vector<std::shared_ptr<sv_t>> vars;
	hts_pos_t current_pos = ref_start + aln.ref_begin, query_pos = aln.query_begin;
	for (int i = 0; i < aln.cigar.size(); i++) {
		uint32_t c = aln.cigar[i];
		int op_length = cigar_int_to_len(c);
		char op = cigar_int_to_op(c);
		if (op == 'D') {
			vars.push_back(std::make_shared<deletion_t>(chr, current_pos-1, current_pos+op_length-1, "", nullptr, nullptr, nullptr, nullptr));
		} else if (op == 'I') {
			vars.push_back(std::make_shared<insertion_t>(chr, current_pos-1, current_pos-1, std::string(alt_seq+query_pos, op_length), nullptr, nullptr, nullptr, nullptr));
		} else if (op == 'X') {
			for (int j = 0; j < op_length; j++) {
				if (is_genomic_base(alt_seq[query_pos+j])) {
					snps.push_back(snp_t(current_pos+j, alt_seq[query_pos+j]));
				}
			}
		}

		if (op != 'I' && op != 'S') current_pos += op_length;
		if (op != 'D') query_pos += op_length;
	}
	return vars;
}

std::shared_ptr<sv_t> update_main_var_from_realigned_vars(std::shared_ptr<sv_t> var, std::vector<std::shared_ptr<sv_t>>& vars, std::vector<snp_t>& snps) {
	if (vars.empty()) return nullptr;

	if (vars[0]->svtype() != var->svtype()) return var;

	var->start = vars[0]->start;
	var->end = vars[0]->end;
	var->ins_seq = vars[0]->ins_seq;
	var->aux_indels.clear();
	for (int i = 1; i < vars.size(); i++) var->aux_indels.push_back(vars[i]);
	var->aux_snps = snps;
	return var;
}

std::shared_ptr<sv_t> realign_indel_haplotype(std::shared_ptr<sv_t> sv) {
	bool has_replaced_ref = sv->start != sv->end && !sv->ins_seq.empty();
	bool has_aux_indels = !sv->aux_indels.empty();
	if (sv->incomplete_ins_seq() || (!has_replaced_ref && !has_aux_indels) || sv->svsize() > 100) return sv;

	char* chr_seq = chr_seqs.get_seq(sv->chr);
	hts_pos_t chr_len = chr_seqs.get_len(sv->chr);
	hts_pos_t ins_start = sv->start, ins_end = sv->end;

	hts_pos_t extend = normalise_max_is;
	std::vector<std::shared_ptr<sv_t>> aux_indels = sv->aux_indels;
	std::vector<snp_t> aux_snps = sv->aux_snps;
	char* lf_seq = generate_haplotype_left(chr_seq, ins_start, extend, aux_indels, aux_snps);
	char* rf_seq = generate_haplotype_right(chr_seq, chr_len, ins_end+1, extend, aux_indels, aux_snps);

	hts_pos_t lf_len = strlen(lf_seq), rf_len = strlen(rf_seq), alt_len = lf_len + sv->ins_seq.length() + rf_len;
	char* putative_alt_allele = new char[alt_len + 1];
	memcpy(putative_alt_allele, lf_seq, lf_len);
	memcpy(putative_alt_allele + lf_len, sv->ins_seq.c_str(), sv->ins_seq.length());
	memcpy(putative_alt_allele + lf_len + sv->ins_seq.length(), rf_seq, rf_len);
	putative_alt_allele[alt_len] = '\0';

	hts_pos_t ref_start = std::max(hts_pos_t(0), ins_start - extend + 1);
	hts_pos_t ref_end = std::min(chr_len, ins_end + extend + 1);
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment aln;
	// keep extending the reference window until we find an alignment without soft clipping at either end (up to 10 times)
	for (int i = 0; i <= 10; i++) {
		realign_aligner.Align(putative_alt_allele, chr_seq + ref_start, ref_end-ref_start, filter, &aln, 0);
		bool left_clipped = get_left_clip_size(aln) > 0;
		bool right_clipped = get_right_clip_size(aln) > 0;
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

	if (!is_left_clipped(aln) && !is_right_clipped(aln)) {
		std::vector<snp_t> snps;
		std::vector<std::shared_ptr<sv_t>> realigned_svs = vars_from_alignment(aln, sv->chr, ref_start, putative_alt_allele, snps);
		sv = update_main_var_from_realigned_vars(sv, realigned_svs, snps);
	}

	delete[] lf_seq;
	delete[] rf_seq;

	delete[] putative_alt_allele;

	return sv;
}

std::shared_ptr<sv_t> realign(std::shared_ptr<sv_t> sv) {
	std::string svtype = sv->svtype();
	if (svtype == "INS") {
		return realign_indel_haplotype(sv);
	} else if (svtype == "DEL") {
		return realign_indel_haplotype(sv);
	} else {
		return sv;
	}
}

int main(int argc, char* argv[]) {

	std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];
	std::string reference_fname = argv[3];

	int n_threads = 1;
	if (argc >= 5) {
		n_threads = std::max(1, std::stoi(argv[4]));
	}
	int min_indel_size = 0;
	if (argc >= 6) {
		min_indel_size = std::stoi(argv[5]);
	}
	if (argc >= 7) {
		normalise_max_is = std::stoi(argv[6]);
	}

	chr_seqs.read_fasta_into_map(reference_fname);

	htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	hdr = bcf_hdr_read(in_vcf_file);
	bcf1_t* vcf_record = bcf_init();

	std::vector<std::shared_ptr<sv_t>> input_svs;
	while (bcf_read(in_vcf_file, hdr, vcf_record) == 0) {
		std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, vcf_record);
		if (sv == nullptr) continue;

		sv->vcf_entry = bcf_dup(vcf_record);
		input_svs.push_back(sv);
	}

	const int block_size = 1000;
	ctpl::thread_pool simplify_realign_thread_pool(n_threads);
	std::vector<std::future<void> > futures;
	for (int i = 0; i < input_svs.size(); i += block_size) {
		int start = i, end = std::min(i+block_size, (int) input_svs.size());
		std::future<void> future = simplify_realign_thread_pool.push(simplify_and_realign_svs, std::ref(input_svs), start, end);
		futures.push_back(std::move(future));
	}
	simplify_realign_thread_pool.stop(true);
	for (int i = 0; i < futures.size(); i++) {
		futures[i].get();
	}
	futures.clear();

	std::vector<std::shared_ptr<sv_t>> svs;
	for (std::shared_ptr<sv_t> sv : input_svs) {
		if (sv == nullptr) continue;

		// TODO: it's possible that since we perform realignment, atomization is not really necessary. Verify
		sv = atomize(0, sv);
		if (sv == nullptr) continue;

		if ((sv->svtype() == "DEL" || sv->svtype() == "INS") && max_indel_size(sv) < min_indel_size) continue;
		
		left_align(sv);
		canonicalize_aux(sv);
		if (sv->source == "READ" || sv->source == "HP") {
			sv->junction_remap_ref_beg = sv->start;
			sv->junction_remap_ref_end = sv->end + 1;
		}
		svs.push_back(sv);
	}

	for (const auto& sv : svs) {
		bcf1_t* vcf_record_norm = bcf_init();
		sv2bcf(hdr, vcf_record_norm, sv.get(), chr_seqs.get_seq(sv->chr));
		copy_all_fmt(hdr, sv->vcf_entry, vcf_record_norm);
		normalised_vcf_records.push_back(vcf_record_norm);
	}

	std::stable_sort(normalised_vcf_records.begin(), normalised_vcf_records.end(),
			[](const bcf1_t* b1, const bcf1_t* b2) { return std::tie(b1->rid, b1->pos) < std::tie(b2->rid, b2->pos); });

	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, hdr) != 0) {
		throw std::runtime_error("Failed to write VCF header to " +  out_vcf_fname);
	}
	for (bcf1_t* vcf_record_norm : normalised_vcf_records) {
		if (bcf_write(out_vcf_file, hdr, vcf_record_norm) != 0) {
			throw std::runtime_error("Failed to write VCF record to " +  out_vcf_fname);
		}
	}

	hts_close(in_vcf_file);
	hts_close(out_vcf_file);

	tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf);
}
