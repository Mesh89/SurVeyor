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
#include "../src/haplotype_normalization.h"
#include "../src/vcf_utils.h"
#include "../libs/cptl_stl.h"

chr_seqs_map_t chr_seqs;
bcf_hdr_t* hdr;
int normalise_max_is = 1000;

std::vector<bcf1_t*> normalised_vcf_records;
std::mutex mtx;

hts_pos_t max_indel_size(const std::shared_ptr<sv_t>& sv) {
	hts_pos_t max_size = sv->svsize();
	for (const auto& aux_indel : sv->aux_indels) {
		max_size = std::max(max_size, aux_indel->svsize());
	}
	return max_size;
}

void simplify_and_realign_svs(int id, std::vector<std::shared_ptr<sv_t>>& svs, int start, int end) {
	for (int i = start; i < end; i++) {
		haplotype_normalization::normalization_context_t context(chr_seqs.get_seq(svs[i]->chr), chr_seqs.get_len(svs[i]->chr), normalise_max_is);
		svs[i] = haplotype_normalization::normalize_haplotype(svs[i], context);
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

		if ((sv->svtype() == "DEL" || sv->svtype() == "INS") && max_indel_size(sv) < min_indel_size) continue;

		left_align(sv);
		canonicalize_aux(sv);
		if (sv->source == "READ") {
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
