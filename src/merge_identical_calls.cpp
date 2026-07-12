#include <memory>
#include <unordered_map>
#include <unordered_set>

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "DynamicIntervalTree.h"
#include "types.h"
#include "vcf_utils.h"


int priority(sv_t* sv) {
	if (sv->source == "DE_NOVO_ASSEMBLY" || sv->source == "REFERENCE_GUIDED_ASSEMBLY") {
		if (sv->imprecise || sv->incomplete_ins_seq()) return 7;
		else return 1;
	} else if (sv->source == "2SR") return 2;
	else if (sv->source == "2SR-2SR") return 2;
	else if (sv->source == "HSR-SR" || sv->source == "SR-HSR") return 3;
	else if (sv->source == "2HSR") return 4;
	else if (sv->source == "2SR-DP" || sv->source == "DP-2SR") return 5;
	else if (sv->source == "1SR_LC" || sv->source == "1SR_RC") return 5;
	else if (sv->source == "1HSR_LC" || sv->source == "1HSR_RC") return 6;
	else if (sv->source == "DP-DP") return 7;
	else if (sv->source == "DP") return 7;
	else if (sv->source == "READ") return 8;
	throw std::runtime_error("Unknown source: " + sv->source);
}

bool has_junction_remap_ref_range(std::shared_ptr<sv_t> sv) {
	return sv->junction_remap_ref_beg != HTS_POS_MIN && sv->junction_remap_ref_end != HTS_POS_MIN;
}

bool junction_remap_ref_range_contains(std::shared_ptr<sv_t> outer, std::shared_ptr<sv_t> inner) {
	return has_junction_remap_ref_range(outer) && has_junction_remap_ref_range(inner) &&
		outer->chr == inner->chr &&
		outer->junction_remap_ref_beg <= inner->junction_remap_ref_beg &&
		inner->junction_remap_ref_end <= outer->junction_remap_ref_end;
}

bool sv_is_fully_contained_in_ref_range(std::shared_ptr<sv_t> sv, hts_pos_t ref_beg, hts_pos_t ref_end) {
	hts_pos_t sv_beg = std::min(sv->start, sv->end);
	hts_pos_t sv_end = std::max(sv->start, sv->end);
	if (sv_beg == sv_end) {
		return ref_beg <= sv_beg && sv_beg < ref_end;
	}
	return ref_beg <= sv_beg && sv_end <= ref_end;
}

bool snp_is_contained_in_ref_range(const snp_t& snp, hts_pos_t ref_beg, hts_pos_t ref_end) {
	return ref_beg <= snp.pos && snp.pos < ref_end;
}

bool junction_remap_ref_range_merge_condition(std::shared_ptr<sv_t> contained_sv, std::shared_ptr<sv_t> containing_sv) {
	if (!junction_remap_ref_range_contains(containing_sv, contained_sv)) return false;

	std::unordered_set<std::string> contained_variants;
	contained_variants.insert("SV:" + contained_sv->unique_key(false));
	for (const auto& aux_indel : contained_sv->aux_indels) {
		contained_variants.insert("SV:" + aux_indel->unique_key(false));
	}
	for (const auto& aux_snp : contained_sv->aux_snps) {
		contained_variants.insert(aux_snp.unique_key());
	}

	std::unordered_set<std::string> containing_variants_in_contained_range;
	if (sv_is_fully_contained_in_ref_range(containing_sv, contained_sv->junction_remap_ref_beg, contained_sv->junction_remap_ref_end)) {
		containing_variants_in_contained_range.insert("SV:" + containing_sv->unique_key(false));
	}
	for (const auto& aux_indel : containing_sv->aux_indels) {
		if (sv_is_fully_contained_in_ref_range(aux_indel, contained_sv->junction_remap_ref_beg, contained_sv->junction_remap_ref_end)) {
			containing_variants_in_contained_range.insert("SV:" + aux_indel->unique_key(false));
		}
	}
	for (const auto& aux_snp : containing_sv->aux_snps) {
		if (snp_is_contained_in_ref_range(aux_snp, contained_sv->junction_remap_ref_beg, contained_sv->junction_remap_ref_end)) {
			containing_variants_in_contained_range.insert(aux_snp.unique_key());
		}
	}

	return contained_variants == containing_variants_in_contained_range;
}

int main(int argc, char* argv[]) {

	std::string in_vcf_fname = argv[1];
	std::string out_vcf_fname = argv[2];
	std::string reference_fname = argv[3];

	htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
	bcf_hdr_t* in_vcf_hdr = bcf_hdr_read(in_vcf_file);
	bcf1_t* b = bcf_init();
	std::vector<std::shared_ptr<sv_t>> svs;
	while (bcf_read(in_vcf_file, in_vcf_hdr, b) == 0) {
		std::shared_ptr<sv_t> sv = bcf_to_sv(in_vcf_hdr, b);
		if (sv == nullptr) {
			throw std::runtime_error("Unexpected unsupported variant in internal VCF " + in_vcf_fname + ": " + std::string(b->d.id ? b->d.id : "<no-id>"));
		}
		svs.push_back(sv);
	}

	std::stable_sort(svs.begin(), svs.end(), [](const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) {
		hts_pos_t a_jrr_len = a->junction_remap_ref_end - a->junction_remap_ref_beg;
		hts_pos_t b_jrr_len = b->junction_remap_ref_end - b->junction_remap_ref_beg;
		if (a_jrr_len != b_jrr_len) return a_jrr_len > b_jrr_len;
		return priority(a.get()) < priority(b.get());
	});

	std::unordered_map<std::string, std::unique_ptr<DynamicIntervalTree<hts_pos_t, std::shared_ptr<sv_t>>>> svs_by_jrr;
	std::vector<std::shared_ptr<sv_t>> surviving_svs;
	std::unordered_map<std::string, std::shared_ptr<sv_t>> surviving_sv_keys;
	for (std::shared_ptr<sv_t> sv : svs) {
		if (surviving_sv_keys.count(sv->unique_key()) > 0) { // already seen an identical SV, skip this one
			std::shared_ptr<sv_t> existing_sv = surviving_sv_keys[sv->unique_key()];
			if (!existing_sv->is_pass() && sv->is_pass()) {
				existing_sv->sample_info.filters.swap(sv->sample_info.filters);
			}
			continue;
		}

		if (has_junction_remap_ref_range(sv) && sv->junction_remap_ref_end > sv->junction_remap_ref_beg) {
			auto tree_it = svs_by_jrr.find(sv->chr);

			bool skip = false;
			if (tree_it != svs_by_jrr.end()) {
				std::vector<std::shared_ptr<sv_t>> larger_overlapping_svs = tree_it->second->query(sv->junction_remap_ref_beg, sv->junction_remap_ref_end - 1);
				for (const auto& lo_sv : larger_overlapping_svs) {
					if (junction_remap_ref_range_merge_condition(sv, lo_sv)) {
						if (!lo_sv->is_pass() && sv->is_pass()) {
							lo_sv->sample_info.filters.swap(sv->sample_info.filters);
						}
						skip = true;
						break;
					}
				}
			}
			if (skip) continue;

			if (svs_by_jrr.count(sv->chr) == 0) {
				svs_by_jrr[sv->chr] = std::unique_ptr<DynamicIntervalTree<hts_pos_t, std::shared_ptr<sv_t>>>(
					new DynamicIntervalTree<hts_pos_t, std::shared_ptr<sv_t>>());
			}
			svs_by_jrr[sv->chr]->insert(sv->junction_remap_ref_beg, sv->junction_remap_ref_end - 1, sv);
		}
		surviving_svs.push_back(sv);
		surviving_sv_keys[sv->unique_key()] = sv;
	}

	// sv_output_order orders variants within a contig but does not compare
	// contig names. Include the VCF header RID first so chromosome blocks remain
	// continuous and downstream outputs can be indexed.
	std::sort(surviving_svs.begin(), surviving_svs.end(), [in_vcf_hdr](
			const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) {
		int a_rid = bcf_hdr_name2id(in_vcf_hdr, a->chr.c_str());
		int b_rid = bcf_hdr_name2id(in_vcf_hdr, b->chr.c_str());
		if (a_rid != b_rid) return a_rid < b_rid;
		return sv_output_order(a, b);
	});

	// output
	htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
	if (bcf_hdr_write(out_vcf_file, in_vcf_hdr) != 0) {
		throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
	}

	chr_seqs_map_t chr_seqs;
	chr_seqs.read_fasta_into_map(reference_fname);
	for (std::shared_ptr<sv_t> sv : surviving_svs) {
		sv2bcf(in_vcf_hdr, b, sv.get(), chr_seqs.get_seq(sv->chr), false);
		if (bcf_write(out_vcf_file, in_vcf_hdr, b) != 0) {
			throw std::runtime_error("Failed to write to " + out_vcf_fname + ".");
		}
	}

	bcf_close(out_vcf_file);
}
