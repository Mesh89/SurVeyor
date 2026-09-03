#ifndef VAR_UTILS_H
#define VAR_UTILS_H

#include <algorithm>
#include <cstring>
#include "types.h"

const int MIN_HP_LEN_FOR_HP_GENOTYPING = 3;

inline bool is_homopolymer_indel(sv_t* sv, char* chr_seq) {
	if (sv->svtype() == "DEL") {
		// for deletions, we need
		// 1) empty inserted sequence
		// 2) all deleted bases are the same, and A/C/G/T
		if (!sv->ins_seq.empty() || sv->start >= sv->end) return false;

		char target_base = toupper(chr_seq[sv->start+1]);
		if (target_base != 'A' && target_base != 'C' && target_base != 'G' && target_base != 'T') {
			return false;
		}

		for (int i = sv->start+2; i <= sv->end; i++) {
			if (toupper(chr_seq[i]) != target_base) {
				return false;
			}
		}
		return true;
	} else if (sv->svtype() == "INS") {
		// for insertions, we need
		// 1) start == end
		// 2) all inserted bases are the same, and A/C/G/T
		if (sv->start != sv->end || sv->ins_seq.empty()) return false;
		
		char target_base = toupper(sv->ins_seq[0]);
		if (target_base != 'A' && target_base != 'C' && target_base != 'G' && target_base != 'T') {
			return false;
		}

		for (char c : sv->ins_seq) {
			if (toupper(c) != target_base) {
				return false;
			}
		}
		return true;
	}
	return false;
}

inline char get_homopolymer_base(sv_t* sv, char* chr_seq) {
	if (!is_homopolymer_indel(sv, chr_seq)) {
		return '\0';
	}
	if (sv->svtype() == "DEL") {
		return toupper(chr_seq[sv->start+1]);
	} else if (sv->svtype() == "INS") {
		return toupper(sv->ins_seq[0]);
	}
	return '\0';
}

// The returned range is 0-based, half-open [beg, end)
inline hts_pair_pos_t find_ref_hp_range_for_indel(sv_t* sv, char* chr_seq, hts_pos_t chr_len) {
	char hp_base = get_homopolymer_base(sv, chr_seq);
	if (hp_base == '\0') {
		throw std::runtime_error("SV " + sv->id + " is not a homopolymer indel.");
	}

	hts_pair_pos_t hp_range;
	hp_range.beg = sv->start;
	while (hp_range.beg >= 0 && toupper(chr_seq[hp_range.beg]) == hp_base) {
		hp_range.beg--;
	}
	hp_range.beg++; // move back to the first base of the homopolymer run

	hp_range.end = sv->end + 1;
	while (hp_range.end < chr_len && toupper(chr_seq[hp_range.end]) == hp_base) {
		hp_range.end++;
	}

	return hp_range;
}

inline bool should_genotype_as_hp_indel(sv_t* sv, char* chr_seq, hts_pos_t chr_len) {
	if (!is_homopolymer_indel(sv, chr_seq)) {
		return false;
	}
	hts_pair_pos_t hp_range = find_ref_hp_range_for_indel(sv, chr_seq, chr_len);
	hts_pos_t ref_hp_len = hp_range.end - hp_range.beg;
	if (sv->svtype() == "INS") {
		return ref_hp_len + sv->ins_seq.length() >= MIN_HP_LEN_FOR_HP_GENOTYPING;
	}
	return ref_hp_len >= MIN_HP_LEN_FOR_HP_GENOTYPING;
}

inline bool aux_indel_haplotype_order(const std::shared_ptr<sv_t>& a, const std::shared_ptr<sv_t>& b) {
    if (a->start != b->start) return a->start < b->start;
    if (a->end != b->end) return a->end < b->end;
    return a->ins_seq < b->ins_seq;
}

inline char* generate_haplotype_left(char* chrom_seq, hts_pos_t hap_end, hts_pos_t hap_len, 
    std::vector<std::shared_ptr<sv_t>>& aux_indels, std::vector<snp_t>& aux_snps, std::vector<allele_edit_t>* edits = nullptr) {
    
    std::sort(aux_indels.begin(), aux_indels.end(), aux_indel_haplotype_order);
    std::sort(aux_snps.begin(), aux_snps.end(), [](snp_t& a, snp_t& b) {
        return a.pos < b.pos;
    });

    // Find last SNP / indel strictly to the left of hap_end (i.e., with pos / end < hap_end)
    int curr_snp_idx = aux_snps.size()-1, curr_indel_idx = aux_indels.size()-1;
    while (curr_snp_idx >= 0 && aux_snps[curr_snp_idx].pos >= hap_end) curr_snp_idx--;
    while (curr_indel_idx >= 0 && aux_indels[curr_indel_idx]->end >= hap_end) curr_indel_idx--;
        
    // Output buffer
    char* hap_seq = new char[hap_len + 1];
    std::fill(hap_seq, hap_seq + hap_len, 'N');
    hap_seq[hap_len] = '\0';
    size_t first_edit_idx = edits == nullptr ? 0 : edits->size();
    
    hts_pos_t remaining_hap_len = hap_len;
    while (hap_end > 0 && remaining_hap_len > 0) {
        const hts_pos_t next_snp_pos = (curr_snp_idx >= 0) ? aux_snps[curr_snp_idx].pos : -10;
        const hts_pos_t next_indel_end = (curr_indel_idx >= 0) ? aux_indels[curr_indel_idx]->end : -10;
        
        hts_pos_t copy_start = hap_end - remaining_hap_len + 1;
        if (copy_start < 0) copy_start = 0;
        copy_start = max(copy_start, next_snp_pos+1, next_indel_end+1);
        hts_pos_t copy_len = hap_end - copy_start + 1;

        strncpy(hap_seq + remaining_hap_len - copy_len, chrom_seq + copy_start, copy_len);
        remaining_hap_len -= copy_len;
        hap_end -= copy_len;

        if (next_indel_end+1 == copy_start) { // we insert the indel
            std::shared_ptr<sv_t> indel = aux_indels[curr_indel_idx];
            if (indel->svtype() == "DEL") {
                if (edits != nullptr) edits->push_back({allele_edit_kind_t::INDEL, indel->start+1, indel->end+1, int(remaining_hap_len), int(remaining_hap_len), int(indel->end-indel->start)});
                hap_end = indel->start;
            } else if (indel->svtype() == "INS") {
                hts_pos_t copied = 0;
                hts_pos_t inserted_len = 0;
                for (char base : indel->ins_seq) if (base != '-') inserted_len++;
                for (auto it = indel->ins_seq.rbegin(); it != indel->ins_seq.rend() && copied < remaining_hap_len; ++it) {
                    if (*it == '-') continue;
                    hap_seq[remaining_hap_len - copied - 1] = *it;
                    copied++;
                }
                if (edits != nullptr && copied == inserted_len) edits->push_back({allele_edit_kind_t::INDEL, indel->start+1, indel->end+1, int(remaining_hap_len-copied), int(remaining_hap_len), int(indel->end-indel->start+inserted_len)});
                remaining_hap_len -= copied;
            }
            curr_indel_idx--;
        } else if (next_snp_pos+1 == copy_start && remaining_hap_len > 0) { // we insert the SNP
            hap_seq[remaining_hap_len-1] = aux_snps[curr_snp_idx].alt_base;
            if (edits != nullptr) edits->push_back({allele_edit_kind_t::SNP, aux_snps[curr_snp_idx].pos, aux_snps[curr_snp_idx].pos+1, int(remaining_hap_len-1), int(remaining_hap_len), 1});
            remaining_hap_len--;
            hap_end--;
            curr_snp_idx--;
        }
    }
    if (remaining_hap_len > 0) {
        // move the sequence to the start
        for (hts_pos_t i = 0; i < hap_len-remaining_hap_len; i++) {
            hap_seq[i] = hap_seq[i+remaining_hap_len];
        }
        hap_seq[hap_len-remaining_hap_len] = '\0';
        if (edits != nullptr) for (size_t i = first_edit_idx; i < edits->size(); i++) {
            (*edits)[i].alt_begin -= remaining_hap_len;
            (*edits)[i].alt_end -= remaining_hap_len;
        }
    }
    return hap_seq;
}

inline char* generate_haplotype_right(char* chrom_seq, hts_pos_t chrom_len, hts_pos_t hap_start, hts_pos_t hap_len,
    std::vector<std::shared_ptr<sv_t>>& aux_indels, std::vector<snp_t>& aux_snps, std::vector<allele_edit_t>* edits = nullptr) {

    // Note that aux_indels coordinates are in VCF format

    std::sort(aux_indels.begin(), aux_indels.end(), aux_indel_haplotype_order);
    std::sort(aux_snps.begin(), aux_snps.end(), [](const snp_t& a, const snp_t& b) {
        return a.pos < b.pos;
    });

    // Find first SNP / indel strictly to the right of hap_start (mirror of left's >= hap_end skip)
    int curr_snp_idx = 0, curr_indel_idx = 0;
    while (curr_snp_idx < (int)aux_snps.size() && aux_snps[curr_snp_idx].pos < hap_start) curr_snp_idx++;
    while (curr_indel_idx < (int)aux_indels.size() && aux_indels[curr_indel_idx]->start < hap_start-1) curr_indel_idx++;

    // Output buffer
    char* hap_seq = new char[hap_len + 1];
    std::fill(hap_seq, hap_seq + hap_len, 'N');
    hap_seq[hap_len] = '\0';

    hts_pos_t remaining_hap_len = hap_len;
    hts_pos_t out_pos = 0; // next write position in hap_seq (left-to-right)

    while (remaining_hap_len > 0 && hap_start < chrom_len) {
        const hts_pos_t next_snp_pos =
            (curr_snp_idx < aux_snps.size()) ? aux_snps[curr_snp_idx].pos : chrom_len + 10; // make sure we don't trigger SNP/indel logic when at the end of chromosome
        const hts_pos_t next_indel_start =
            (curr_indel_idx < aux_indels.size()) ? aux_indels[curr_indel_idx]->start : chrom_len + 10;

        // Copy reference until just before the next event (SNP or indel anchor).
        // copy_end is exclusive, i.e., the first base we do NOT copy
        hts_pos_t copy_end = hap_start + remaining_hap_len;
        copy_end = std::min(copy_end, next_snp_pos);
        copy_end = std::min(copy_end, next_indel_start+1); // indel start is the position *BEFORE* the indel
        copy_end = std::min(copy_end, chrom_len);
        hts_pos_t copy_len = copy_end - hap_start;

        strncpy(hap_seq + out_pos, chrom_seq + hap_start, copy_len);
        out_pos += copy_len;
        remaining_hap_len -= copy_len;
        hap_start = copy_end;

        if (hap_start == next_indel_start+1) {
            std::shared_ptr<sv_t> indel = aux_indels[curr_indel_idx];
            if (indel->svtype() == "DEL") {
                if (edits != nullptr) edits->push_back({allele_edit_kind_t::INDEL, indel->start+1, indel->end+1, int(out_pos), int(out_pos), int(indel->end-indel->start)});
                hap_start = indel->end + 1;
            } else if (indel->svtype() == "INS") {
                hts_pos_t copied = 0;
                hts_pos_t inserted_len = 0;
                for (char base : indel->ins_seq) if (base != '-') inserted_len++;
                for (auto it = indel->ins_seq.begin(); it != indel->ins_seq.end() && copied < remaining_hap_len; ++it) {
                    if (*it == '-') continue;
                    hap_seq[out_pos + copied] = *it;
                    copied++;
                }
                if (edits != nullptr && copied == inserted_len) edits->push_back({allele_edit_kind_t::INDEL, indel->start+1, indel->end+1, int(out_pos), int(out_pos+copied), int(indel->end-indel->start+inserted_len)});
                out_pos += copied;
                remaining_hap_len -= copied;
            }
            curr_indel_idx++;
        } else if (hap_start == next_snp_pos && remaining_hap_len > 0) {
            hap_seq[out_pos] = aux_snps[curr_snp_idx].alt_base;
            if (edits != nullptr) edits->push_back({allele_edit_kind_t::SNP, aux_snps[curr_snp_idx].pos, aux_snps[curr_snp_idx].pos+1, int(out_pos), int(out_pos+1), 1});
            out_pos++;
            remaining_hap_len--;
            hap_start++;       // consume this reference base
            curr_snp_idx++;
        } 
   }

    // If we couldn't fill the requested length (ran off chromosome), truncate to actual length
    hap_seq[out_pos] = '\0';
    return hap_seq;
}

#endif // VAR_UTILS_H
