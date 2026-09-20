#ifndef HAPLOTYPE_QUALITY_H
#define HAPLOTYPE_QUALITY_H

#include "haplotype_normalization.h"

using junction_query_origins_t = haplotype_normalization::normalization_provenance_t;

bool passes_aux_snp_quality(hts_pos_t query_pos, size_t query_len, const std::string& qual, int lowq_prefix, int lowq_suffix) {
	return query_pos >= lowq_prefix && query_pos < hts_pos_t(query_len)-lowq_suffix && query_pos >= 0 && query_pos < hts_pos_t(qual.length()) && (uint8_t) qual[query_pos] - 33 >= 40;
}

std::shared_ptr<sv_t> copy_discovery_indel(const std::shared_ptr<sv_t>& sv, junction_query_origins_t& origins) {
	std::shared_ptr<sv_t> copy;
	if (sv->svtype() == "DEL") copy = std::make_shared<deletion_t>(*std::static_pointer_cast<deletion_t>(sv));
	else if (sv->svtype() == "INS") copy = std::make_shared<insertion_t>(*std::static_pointer_cast<insertion_t>(sv));
	else if (sv->svtype() == "RPL") copy = std::make_shared<replacement_t>(*std::static_pointer_cast<replacement_t>(sv));
	else return nullptr;
	origins.inserted[copy.get()] = haplotype_normalization::insertion_origin_slice(&origins, sv.get(), 0, sv->ins_seq.length());
	return copy;
}

std::shared_ptr<sv_t> filter_decomposed_aux_snps(std::shared_ptr<sv_t> sv, const haplotype_normalization::normalization_context_t& context,
	const std::string& query, const std::string& qual, int lowq_prefix, int lowq_suffix, const junction_query_origins_t& query_origins) {
	if (context.extend <= 0) return sv;
	if (sv->aux_indels.empty() && (sv->start == sv->end || sv->ins_seq.empty())) return sv;
	if (lowq_prefix == 0 && lowq_suffix == 0 && qual.length() >= query.length() && std::all_of(qual.begin(), qual.end(), [](char q) { return (uint8_t) q - 33 >= 40; })) return sv;

	// Validate the haplotype before invoking the existing flank builders. Discovery
	// does not guess how overlapping, incomplete, or nested AUX edits should compose.
	std::vector<std::shared_ptr<sv_t>> edits{sv};
	edits.insert(edits.end(), sv->aux_indels.begin(), sv->aux_indels.end());
	for (const auto& edit : edits) {
		if (edit->vcf_entry != nullptr || edit->chr != sv->chr || edit->start < 0 || edit->end < edit->start || edit->end >= context.chr_len || edit->incomplete_ins_seq()) return sv;
		if (edit->svtype() != "INS" && edit->svtype() != "DEL" && edit->svtype() != "RPL") return sv;
		if (edit != sv && (!edit->aux_indels.empty() || !edit->aux_snps.empty() || haplotype_normalization::is_complex_indel(edit))) return sv;
	}
	std::sort(edits.begin(), edits.end(), aux_indel_haplotype_order);
	hts_pos_t cursor = edits.front()->start+1;
	for (size_t i = 0; i < edits.size(); i++) {
		const auto& edit = edits[i];
		if (edit->start+1 < cursor || (i && edit->start == edits[i-1]->start && edit->end == edits[i-1]->end)) return sv;
		cursor = edit->end+1;
		for (const auto& snp : sv->aux_snps) if (edit->start < snp.pos && snp.pos <= edit->end) return sv;
	}
	for (const auto& snp : sv->aux_snps) if (snp.pos < 0 || snp.pos >= context.chr_len) return sv;

	// Normalization may replace the primary variant or modify AUX edits. Use an
	// independent copy so no-failure and unresolved-provenance paths are exact no-ops.
	junction_query_origins_t provenance = query_origins;
	provenance.query = &query;
	std::shared_ptr<sv_t> candidate = copy_discovery_indel(sv, provenance);
	candidate->aux_indels.clear();
	for (const auto& aux : sv->aux_indels) candidate->aux_indels.push_back(copy_discovery_indel(aux, provenance));
	candidate = haplotype_normalization::normalize_haplotype(candidate, context, &provenance);
	if (candidate == nullptr) return sv;
	std::vector<snp_t> retained_snps;
	bool failed = false;
	for (const auto& snp : candidate->aux_snps) {
		hts_pos_t origin = provenance.snp_origin(snp);
		if (origin < 0 || origin >= hts_pos_t(query.length()) || toupper(query[origin]) != toupper(snp.alt_base)) return sv;
		if (passes_aux_snp_quality(origin, query.length(), qual, lowq_prefix, lowq_suffix)) retained_snps.push_back(snp);
		else failed = true;
	}
	if (!failed) return sv;
	candidate->aux_snps = std::move(retained_snps);
	return candidate;
}

#endif
