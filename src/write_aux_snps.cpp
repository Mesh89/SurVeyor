#include <algorithm>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <fstream>

#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include "DynamicIntervalTree.h"
#include "SegTree.h"
#include "types.h"
#include "vcf_utils.h"

const hts_pos_t MIN_SV_SIZE = 50;
std::ofstream merge_identical_log;

void index_vcf(const std::string& vcf_fname) {
    if (tbx_index_build(vcf_fname.c_str(), 0, &tbx_conf_vcf) != 0) {
        throw std::runtime_error("Failed to index " + vcf_fname + ".");
    }
}

void copy_hp_info(bcf_hdr_t* hdr, bcf1_t* src, bcf1_t* dest) {
    int hp_genotyped = 0;
    int hp_genotyped_len = 0;
    if (bcf_get_info_flag(hdr, src, "HP_GENOTYPED", &hp_genotyped, &hp_genotyped_len) > 0) {
        bcf_update_info_flag(hdr, dest, "HP_GENOTYPED", "", 1);
    }

    int32_t* hp_ref_range = NULL;
    int n_hp_ref_range = 0;
    if (bcf_get_info_int32(hdr, src, "HP_REF_RANGE", &hp_ref_range, &n_hp_ref_range) > 0) {
        bcf_update_info_int32(hdr, dest, "HP_REF_RANGE", hp_ref_range, n_hp_ref_range);
    }
    free(hp_ref_range);
}

bool is_literal_allele_record(bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele < 2) return false;
    const char* alt = record->d.allele[1];
    return alt[0] != '<' && strchr(alt, '[') == NULL && strchr(alt, ']') == NULL;
}

bool same_record_identity(bcf_hdr_t* hdr, bcf1_t* b1, bcf1_t* b2) {
    bcf_unpack(b1, BCF_UN_STR);
    bcf_unpack(b2, BCF_UN_STR);

    if (b1->rid != b2->rid || b1->pos != b2->pos) return false;

    bool b1_is_literal = is_literal_allele_record(b1);
    bool b2_is_literal = is_literal_allele_record(b2);
    if (b1_is_literal != b2_is_literal) return false;

    if (b1_is_literal) {
        return strcmp(b1->d.allele[0], b2->d.allele[0]) == 0 &&
               strcmp(b1->d.allele[1], b2->d.allele[1]) == 0;
    }

    return get_sv_end(hdr, b1) == get_sv_end(hdr, b2) && get_ins_seq(hdr, b1) == get_ins_seq(hdr, b2) 
        && get_sv_type(hdr, b1) == get_sv_type(hdr, b2);
}

bool is_small_variant_by_svsize(bcf_hdr_t* hdr, bcf1_t* record) {
    std::string svtype = get_sv_type(hdr, record);
    if (svtype == "BND") return false;

    if (svtype == "DEL" || svtype == "DUP" || svtype == "INS" || svtype == "INV") {
        std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, record);
        if (sv == nullptr) return false;
        if (svtype == "INS" && sv->incomplete_ins_seq()) return false;
        return sv->svsize() < MIN_SV_SIZE;
    }

    int len = get_sv_len(hdr, record);
    return len != bcf_int32_missing && len > -MIN_SV_SIZE && len < MIN_SV_SIZE;
}

void divide_variants_by_svsize(bcf_hdr_t* hdr, const std::vector<bcf1_t*>& vcf_records,
                              std::vector<bcf1_t*>& small_variants, std::vector<bcf1_t*>& svs) {
    small_variants.clear();
    svs.clear();

    for (bcf1_t* record : vcf_records) {
        if (is_small_variant_by_svsize(hdr, record)) {
            small_variants.push_back(record);
        } else {
            svs.push_back(record);
        }
    }
}

std::string bcf_record_unique_key(bcf_hdr_t* hdr, bcf1_t* record) {
    return std::string(bcf_seqname_safe(hdr, record)) + ":" +
        std::to_string(record->pos) + ":" +
        std::to_string(get_sv_end(hdr, record)) + ":" +
        get_sv_type(hdr, record) + ":" +
        get_ins_seq(hdr, record);
}

int duplicate_priority(bcf_hdr_t* hdr, bcf1_t* record,
                       const std::unordered_set<bcf1_t*>& aux_records) {
    if (aux_records.count(record) > 0) return 1;

    std::string source = get_sv_info_str(hdr, record, "SOURCE");
    if (source == "READ") return 0;

    return 2;
}

void remove_lower_priority_duplicates(bcf_hdr_t* hdr, std::vector<bcf1_t*>& vcf_records,
                                        const std::unordered_set<bcf1_t*>& aux_records) {
    std::unordered_map<std::string, int> max_priority_by_key;
    for (bcf1_t* record : vcf_records) {
        if (count_alt_alleles(hdr, record) == 0) {
            continue;
        }
        std::string key = bcf_record_unique_key(hdr, record);
        int priority = duplicate_priority(hdr, record, aux_records);
        max_priority_by_key[key] = std::max(max_priority_by_key[key], priority);
    }

    for (bcf1_t*& record : vcf_records) {
        if (count_alt_alleles(hdr, record) == 0) {
            continue;
        }
        std::string key = bcf_record_unique_key(hdr, record);
        if (duplicate_priority(hdr, record, aux_records) < max_priority_by_key[key]) {
            bcf_destroy(record);
            record = nullptr;
        }
    }

    vcf_records.erase(std::remove(vcf_records.begin(), vcf_records.end(), nullptr), vcf_records.end());
}

std::string get_record_id(bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    return record->d.id == NULL ? "." : std::string(record->d.id);
}

bool is_expanded_aux_indel_record(bcf1_t* record) {
    return get_record_id(record).find(".AUX_INDEL.") != std::string::npos;
}

std::string get_gt_str(bcf_hdr_t* hdr, bcf1_t* record) {
    int* gt = nullptr;
    int ngt = 0;
    if (bcf_get_genotypes(hdr, record, &gt, &ngt) < 0 || ngt == 0) {
        free(gt);
        return ".";
    }

    std::stringstream ss;
    for (int i = 0; i < ngt; i++) {
        if (i > 0) ss << "/";
        if (bcf_gt_is_missing(gt[i])) {
            ss << ".";
        } else {
            ss << bcf_gt_allele(gt[i]);
        }
    }
    free(gt);
    return ss.str();
}

void log_merge_identical_operation(bcf_hdr_t* hdr, const char* operation,
                                   bcf1_t* first, bcf1_t* second,
                                   const std::string& detail = "") {
    if (!merge_identical_log.is_open()) return;

    merge_identical_log
        << operation
        << "\tfirst=" << get_record_id(first)
        << "\tsecond=" << get_record_id(second)
        << "\tfirst_gt=" << get_gt_str(hdr, first)
        << "\tsecond_gt=" << get_gt_str(hdr, second);
    if (!detail.empty()) {
        merge_identical_log << "\t" << detail;
    }
    merge_identical_log << "\n";
}

int get_required_hpid(bcf_hdr_t* hdr, bcf1_t* record) {
    int32_t* hpid = NULL;
    int hpid_len = 0;
    if (bcf_get_info_int32(hdr, record, "HPID", &hpid, &hpid_len) <= 0 ||
            hpid_len == 0 ||
            hpid[0] == bcf_int32_missing ||
            hpid[0] == bcf_int32_vector_end) {
        free(hpid);
        throw std::runtime_error("Missing HPID for " + get_record_id(record) + ".");
    }
    int value = hpid[0];
    free(hpid);
    return value;
}

void set_gt(bcf_hdr_t* hdr, bcf1_t* record, int allele1, int allele2) {
    int gt[2] = { bcf_gt_unphased(allele1), bcf_gt_unphased(allele2) };
    bcf_update_genotypes(hdr, record, gt, 2);
}

bool is_hom_alt(bcf_hdr_t* hdr, bcf1_t* record) {
    int* gt = nullptr;
    int ngt = 0;
    if (bcf_get_genotypes(hdr, record, &gt, &ngt) < 0 || ngt < 2) {
        free(gt);
        return false;
    }
    int allele1 = bcf_gt_allele(gt[0]);
    int allele2 = bcf_gt_allele(gt[1]);
    free(gt);
    return allele1 > 0 && allele1 == allele2;
}

int get_format_int32_or_default(bcf_hdr_t* hdr, bcf1_t* record, const char* tag, int default_value = 0) {
    int32_t* values = nullptr;
    int n_values = 0;
    int value = default_value;
    if (bcf_get_format_int32(hdr, record, tag, &values, &n_values) > 0 &&
            values[0] != bcf_int32_missing && values[0] != bcf_int32_vector_end) {
        value = values[0];
    }
    free(values);
    return value;
}

bool ref_support_exceeds_other_alt_support(bcf_hdr_t* hdr, bcf1_t* record) {
    int rr = get_format_int32_or_default(hdr, record, "RR1") + get_format_int32_or_default(hdr, record, "RR2");
    int orc = get_format_int32_or_default(hdr, record, "OR1") + get_format_int32_or_default(hdr, record, "OR2");
    return rr > orc;
}

void concat_record_ids(bcf_hdr_t* hdr, bcf1_t* record, bcf1_t* other_record) {
    std::string id = get_record_id(record);
    std::string other_id = get_record_id(other_record);
    if (id == ".") {
        id = other_id;
    } else if (other_id != ".") {
        id += ";" + other_id;
    }
    bcf_update_id(hdr, record, id.c_str());
}

void make_record_multiallelic(bcf_hdr_t* hdr, bcf1_t* record, bcf1_t* second_allele_record) {
    bcf_unpack(record, BCF_UN_STR);
    bcf_unpack(second_allele_record, BCF_UN_STR);
    std::string ref = record->d.allele[0];
    std::string alt1 = record->d.allele[1];
    std::string alt2 = second_allele_record->d.allele[1];

    std::string ref2 = second_allele_record->d.allele[0];
    if (ref2.length() > ref.length() && ref2.compare(0, ref.length(), ref) == 0) {
        alt1 += ref2.substr(ref.length());
        ref = ref2;
    } else if (ref.length() > ref2.length() && ref.compare(0, ref2.length(), ref2) == 0) {
        alt2 += ref.substr(ref2.length());
    }

    std::string alleles = ref + "," + alt1 + "," + alt2;
    bcf_update_alleles_str(hdr, record, alleles.c_str());
    concat_record_ids(hdr, record, second_allele_record);
    set_gt(hdr, record, 1, 2);
}

void add_star_allele(bcf_hdr_t* hdr, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    if (record->n_allele != 2 || count_alt_alleles(hdr, record) != 1) {
        return;
    }

    std::string alleles = record->d.allele[0];
    for (int i = 1; i < record->n_allele; i++) {
        if (strcmp(record->d.allele[i], "*") == 0) {
            return;
        }
        alleles += ",";
        alleles += record->d.allele[i];
    }

    int star_allele = record->n_allele;
    alleles += ",*";
    bcf_update_alleles_str(hdr, record, alleles.c_str());
    set_gt(hdr, record, 1, star_allele);
}

void add_star_alleles_for_overlapping(bcf_hdr_t* hdr, std::vector<bcf1_t*>& small_variants) {
    std::unordered_set<bcf1_t*> read_records_to_remove;
    int curr_rid = -1;
    hts_pos_t max_del_end = -1;
    hts_pos_t max_read_hom_del_end = -1;
    bcf1_t* max_read_hom_del_record = nullptr;
    hts_pos_t max_non_read_hom_del_end = -1;
    bcf1_t* max_non_read_hom_del_record = nullptr;

    for (bcf1_t* record : small_variants) {
        if (record->rid != curr_rid) {
            curr_rid = record->rid;
            max_del_end = -1;
            max_read_hom_del_end = -1;
            max_read_hom_del_record = nullptr;
            max_non_read_hom_del_end = -1;
            max_non_read_hom_del_record = nullptr;
        }

        std::string svtype = get_sv_type(hdr, record);
        std::string source = get_sv_info_str(hdr, record, "SOURCE");
        int alt_alleles = count_alt_alleles(hdr, record);
        if (is_hom_alt(hdr, record)) {
            if (source == "READ" && record->pos <= max_non_read_hom_del_end &&
                    max_non_read_hom_del_record != nullptr) {
                read_records_to_remove.insert(record);
            } else if (source != "READ" && record->pos <= max_read_hom_del_end &&
                    max_read_hom_del_record != nullptr) {
                read_records_to_remove.insert(max_read_hom_del_record);
            }
        }
        if (alt_alleles == 1 && record->pos <= max_del_end) {
            add_star_allele(hdr, record);
        }
        if (svtype == "DEL" && alt_alleles > 0) {
            hts_pos_t del_end = get_sv_end(hdr, record);
            max_del_end = std::max(max_del_end, del_end);
            if (is_hom_alt(hdr, record) && source == "READ" && del_end > max_read_hom_del_end) {
                max_read_hom_del_end = del_end;
                max_read_hom_del_record = record;
            } else if (is_hom_alt(hdr, record) && source != "READ" && del_end > max_non_read_hom_del_end) {
                max_non_read_hom_del_end = del_end;
                max_non_read_hom_del_record = record;
            }
        }
    }

    for (bcf1_t*& record : small_variants) {
        if (read_records_to_remove.count(record) > 0) {
            record = nullptr;
        }
    }
    small_variants.erase(std::remove(small_variants.begin(), small_variants.end(), nullptr), small_variants.end());
}

// Merge duplicate records for the same allele into the higher-priority record.
// Two het records are promoted to 1/1, unless the first record's
// classifier requires the other allele to be reference. In those exceptions,
// the first genotype is left unchanged
void merge_variant_records(bcf_hdr_t* hdr, bcf1_t* first, bcf1_t* second) {
    int first_alt_alleles = count_alt_alleles(hdr, first);
    int second_alt_alleles = count_alt_alleles(hdr, second);
    log_merge_identical_operation(hdr, "start", first, second,
        "first_ac=" + std::to_string(first_alt_alleles) +
        "\tsecond_ac=" + std::to_string(second_alt_alleles));
    if (first_alt_alleles == 0 || second_alt_alleles == 0) {
        log_merge_identical_operation(hdr, "skip_zero_ac", first, second);
        return;
    }

    // From here onwards, both records are assumed positive

    // If both records are identical, and EREFA does not force a ref allele to exist, merge them into a 1/1
    if (same_record_identity(hdr, first, second)) {
        int first_erefa = get_format_int32_or_default(hdr, first, "EREFA");
        if (first_erefa == 1) {
            log_merge_identical_operation(hdr, "skip_erefa", first, second,
                "EREFA=" + std::to_string(first_erefa));
                return;
            }

            if (first_alt_alleles == 1) concat_record_ids(hdr, first, second);
            log_merge_identical_operation(hdr, "concat_ids", first, second);
            set_gt(hdr, first, 1, 1);
            log_merge_identical_operation(hdr, "set_first_1_1", first, second);
            set_gt(hdr, second, 0, 0);
            log_merge_identical_operation(hdr, "set_second_0_0", first, second);
    }
    // If the first record is a 1/1, but the probability of the second record to exist (EPR) is higher
    // than the probability of the first record to be homozygous (HOPR), accept both as 0/1
    else if (first_alt_alleles >= 2 && !same_record_identity(hdr, first, second)) {
        float second_epr = get_sv_epr(hdr, second);
        float first_hopr = get_sv_hopr(hdr, first);
        std::string scores = "second_epr=" + std::to_string(second_epr) +
            "\tfirst_hopr=" + std::to_string(first_hopr);
        if (second_epr > first_hopr) {
            log_merge_identical_operation(hdr, "replace_one_alt", first, second, scores);
            set_gt(hdr, first, 0, 1);
            set_gt(hdr, second, 0, 1);
            log_merge_identical_operation(hdr, "set_for_multiallelic", first, second, scores);
        } else {
            set_gt(hdr, second, 0, 0);
            log_merge_identical_operation(hdr, "keep_first_multiallelic", first, second, scores);
        }
    } else {
        log_merge_identical_operation(hdr, "no_action", first, second);
    }

    // bool two_hets = first_alt_alleles == 1 && second_alt_alleles == 1;
    // if (two_hets && ref_support_exceeds_other_alt_support(hdr, second)) return;
}

// The group is already sorted by decreasing EPR. Keep at most two ALT alleles at
// each exact position for INS/DEL/SNP groups. The best call is always kept. If the first two
// calls are identical hets, merge them into a single 1/1 call unless the second
// call has more REF support than competing ALT support. Otherwise, if the
// best call is already 1/1, let the second call replace one allele only when its
// EPR is stronger than the best call's homozygous probability. All later ALT
// calls are suppressed. Accepted second alleles are folded into the first record
// as ALT2, and a 0/0 record ends the group because no later record can add an
// ALT allele.
void apply_multiallelic_logic_to_group(bcf_hdr_t* hdr, std::vector<bcf1_t*>& variant_group,
                                       std::vector<bcf1_t*>& records_to_remove) {
    
    if (variant_group.empty()) return;

    bcf1_t* first = variant_group[0];
    int first_alt_alleles = count_alt_alleles(hdr, first);
    if (first_alt_alleles == 0) return;

    if (variant_group.size() == 1) return;

    bcf1_t* second = variant_group[1];
    int second_alt_alleles = count_alt_alleles(hdr, second);
    if (second_alt_alleles == 0) return;

    size_t suppress_from = 2;

    if (first_alt_alleles > 0 && second_alt_alleles > 0 && same_record_identity(hdr, first, second)) {
        merge_variant_records(hdr, first, second);
        records_to_remove.push_back(second);
    } else if (first_alt_alleles >= 2) {
        float second_epr = get_sv_epr(hdr, second);
        float first_hopr = get_sv_hopr(hdr, first);
        if (second_epr > first_hopr) {
            set_gt(hdr, first, 0, 1);
            if (second_alt_alleles > 1) {
                set_gt(hdr, second, 0, 1);
            }
            make_record_multiallelic(hdr, first, second);
            records_to_remove.push_back(second);
        } else {
            set_gt(hdr, second, 0, 0);
            records_to_remove.push_back(second);
            suppress_from = 2;
        }
    } else {
        if (second_alt_alleles > 1) {
            set_gt(hdr, second, 0, 1);
        }
        make_record_multiallelic(hdr, first, second);
        records_to_remove.push_back(second);
    }

    for (size_t i = suppress_from; i < variant_group.size(); i++) {
        bcf1_t* record = variant_group[i];
        int alt_alleles = count_alt_alleles(hdr, record);
        if (alt_alleles == 0) {
            return;
        }

        set_gt(hdr, record, 0, 0);
        records_to_remove.push_back(record);
    }
}

void make_multiallelic(bcf_hdr_t* hdr, std::vector<bcf1_t*>& small_variants) {
    int curr_rid = -1;
    hts_pos_t curr_pos = -1;
    std::vector<bcf1_t*> variant_group;
    std::vector<bcf1_t*> records_to_remove;

    for (bcf1_t* record : small_variants) {
        if (!variant_group.empty() && (record->rid != curr_rid || record->pos != curr_pos)) {
            apply_multiallelic_logic_to_group(hdr, variant_group, records_to_remove);
            variant_group.clear();
        }

        if (variant_group.empty()) {
            curr_rid = record->rid;
            curr_pos = record->pos;
        }

        std::string svtype = get_sv_type(hdr, record);
        if (svtype == "INS" || svtype == "DEL" || svtype == "SNP") {
            variant_group.push_back(record);
        }
    }

    apply_multiallelic_logic_to_group(hdr, variant_group, records_to_remove);

    std::unordered_set<bcf1_t*> remove_set(records_to_remove.begin(), records_to_remove.end());
    for (bcf1_t*& record : small_variants) {
        if (remove_set.count(record) > 0) {
            record = nullptr;
        }
    }
    small_variants.erase(std::remove(small_variants.begin(), small_variants.end(), nullptr), small_variants.end());
}

void make_multiallelic_new(bcf_hdr_t* hdr, std::vector<bcf1_t*>& small_variants) {
    std::unordered_map<std::string, std::vector<std::pair<bcf1_t*, float>>> variants_by_chr;
    std::unordered_map<std::string, hts_pos_t> max_end_by_chr;

    for (bcf1_t* record : small_variants) {
        if (count_alt_alleles(hdr, record) == 0) continue;

        std::string svtype = get_sv_type(hdr, record);
        if (svtype != "INS" && svtype != "DEL" && svtype != "SNP") continue;

        std::string chr = bcf_seqname_safe(hdr, record);
        variants_by_chr[chr].emplace_back(record, get_sv_epr(hdr, record));
        max_end_by_chr[chr] = std::max(max_end_by_chr[chr], (hts_pos_t) get_sv_end(hdr, record));
    }

    for (auto& chr_variants : variants_by_chr) {
        const std::string& chr = chr_variants.first;
        std::vector<std::pair<bcf1_t*, float>>& variants = chr_variants.second;

        std::stable_sort(variants.begin(), variants.end(),
            [](const std::pair<bcf1_t*, float>& a, const std::pair<bcf1_t*, float>& b) {
                return a.second > b.second;
            });

        SegTree seg_tree(max_end_by_chr[chr] + 1);
        DynamicIntervalTree<hts_pos_t, bcf1_t*> interval_tree;
        for (const auto& variant_epr : variants) {
            bcf1_t* record = variant_epr.first;
            if (count_alt_alleles(hdr, record) == 0) continue;

            hts_pos_t start = record->pos;
            hts_pos_t end = get_sv_end(hdr, record);

            if (seg_tree.any_ge(start, end, 2)) {
                set_gt(hdr, record, 0, 0);
                continue;
            }

            std::vector<bcf1_t*> overlapping_variants;
            overlapping_variants = interval_tree.query(start, end);
            if (overlapping_variants.size() == 1) {
                merge_variant_records(hdr, overlapping_variants[0], record);
            }

            seg_tree.add(start, end, 1);
            interval_tree.insert(start, end, record);
        }
    }
}

int main(int argc, char* argv[]) {
    std::string in_vcf_fname = argv[1];
    std::string out_smvars_prefix = argv[2];
    std::string out_stvars_prefix = argv[3];
    std::string reference_fname = argv[4];
    std::string out_smvars_vcf_fname = out_smvars_prefix + ".vcf.gz";
    std::string out_stvars_vcf_fname = out_stvars_prefix + ".vcf.gz";
    std::string merge_identical_log_fname = out_smvars_prefix + ".merge-identical.log";

    merge_identical_log.open(merge_identical_log_fname);
    if (!merge_identical_log) {
        throw std::runtime_error("Failed to open " + merge_identical_log_fname + " for writing.");
    }
    merge_identical_log << "operation\tfirst\tsecond\tfirst_gt\tsecond_gt\tdetails\n";

    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(reference_fname);

    // count the number of records for each HPID in the input VCF
    // greater than 1 means that the parent record was already expanded in expand_aux_haplotypes.cpp
    // so we should not expand it again here, i.e., create aux INDEL records for it
    std::unordered_map<int, int> hpid_counts;
    htsFile* count_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    bcf_hdr_t* count_vcf_hdr = bcf_hdr_read(count_vcf_file);
    bcf1_t* count_record = bcf_init();
    while (bcf_read(count_vcf_file, count_vcf_hdr, count_record) == 0) {
        hpid_counts[get_required_hpid(count_vcf_hdr, count_record)]++;
    }
    bcf_destroy(count_record);
    bcf_hdr_destroy(count_vcf_hdr);
    hts_close(count_vcf_file);

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    bcf_hdr_t* in_vcf_hdr = bcf_hdr_read(in_vcf_file);

    std::vector<bcf1_t*> vcf_records;
    std::unordered_set<bcf1_t*> aux_records;
	bcf1_t* r = bcf_init();
	while (bcf_read(in_vcf_file, in_vcf_hdr, r) == 0) {
        bcf1_t* b = bcf_dup(r);
		vcf_records.push_back(b);

        char* s_data = NULL;
        int len = 0;
        bcf_get_info_string(in_vcf_hdr, b, "AUX_SNPS", (void**) &s_data, &len);
        if (len > 0 && !is_expanded_aux_indel_record(b)) {
            bcf_unpack(b, BCF_UN_ALL);
            int parent_hpid = get_required_hpid(in_vcf_hdr, b);
            std::istringstream ss(s_data);
            std::string snp_str;
            int i = 0;
            while (std::getline(ss, snp_str, ',')) {
                snp_t snp(snp_str);
                std::string contig_name = bcf_seqname(in_vcf_hdr, b);
                char* chr_seq = chr_seqs.get_seq(contig_name);
                std::string id = std::string(b->d.id) + ".SNP." + std::to_string(i++);
                std::vector<int> gt = get_bcf_gt(in_vcf_hdr, b);
                bcf1_t* snp_record = generate_snp(in_vcf_hdr, bcf_seqname(in_vcf_hdr, b), snp.pos,
                    chr_seq[snp.pos], snp.alt_base, id, gt);
                bcf_update_info_int32(in_vcf_hdr, snp_record, "HPID", &parent_hpid, 1);
                copy_all_fmt(in_vcf_hdr, b, snp_record);
                vcf_records.push_back(snp_record);
                aux_records.insert(snp_record);
            }
            // remove INFO/AUX_SNPS from the original record
            bcf_update_info_string(in_vcf_hdr, b, "AUX_SNPS", NULL);
        }
        free(s_data);

        s_data = NULL;
        len = 0;
        bcf_get_info_string(in_vcf_hdr, b, "AUX_INDELS", (void**) &s_data, &len);
        if (len > 0 && !is_expanded_aux_indel_record(b)) {
            bcf_unpack(b, BCF_UN_ALL);
            int parent_hpid = get_required_hpid(in_vcf_hdr, b);
            std::istringstream ss(s_data);;
            std::string indel_str;
            int i = 0;
            while (hpid_counts[parent_hpid] == 1 && std::getline(ss, indel_str, ',')) {
                std::stringstream indel_ss(indel_str);
                std::string start_str, end_str, ins_seq;
                std::getline(indel_ss, start_str, ':');
                std::getline(indel_ss, end_str, ':');
                std::getline(indel_ss, ins_seq, ':');
                std::string contig_name = bcf_seqname(in_vcf_hdr, b);
                hts_pos_t start = std::stoll(start_str)-1;
                hts_pos_t end = std::stoll(end_str)-1;
                std::shared_ptr<sv_t> indel;
                if (start == end) {
                    // insertion
                    indel = std::make_shared<insertion_t>(contig_name, start, end, ins_seq, nullptr, nullptr, nullptr, nullptr);
                } else {
                    // deletion
                    indel = std::make_shared<deletion_t>(contig_name, start, end, ins_seq, nullptr, nullptr, nullptr, nullptr);
                }
                indel->id = std::string(b->d.id) + ".INDEL." + std::to_string(i++);
                indel->hpid = parent_hpid;
                std::vector<int> gt = get_bcf_gt(in_vcf_hdr, b);
                indel->sample_info.gt = gt;
                bcf1_t* indel_record = bcf_init();
                sv2bcf(in_vcf_hdr, indel_record, indel.get(), chr_seqs.get_seq(contig_name));
                copy_all_fmt(in_vcf_hdr, b, indel_record);
                copy_hp_info(in_vcf_hdr, b, indel_record);
                vcf_records.push_back(indel_record);
                aux_records.insert(indel_record);
            }
            // remove INFO/AUX_INDELS from the original record
            bcf_update_info_string(in_vcf_hdr, b, "AUX_INDELS", NULL);
        }
        free(s_data);
    }

    // sort by pos and then by non-ascending EPR
    std::stable_sort(vcf_records.begin(), vcf_records.end(),
        [in_vcf_hdr](bcf1_t* b1, bcf1_t* b2) {
            if (std::tie(b1->rid, b1->pos) != std::tie(b2->rid, b2->pos)) {
                return std::tie(b1->rid, b1->pos) < std::tie(b2->rid, b2->pos);
            }
            return get_sv_epr(in_vcf_hdr, b1) > get_sv_epr(in_vcf_hdr, b2);
        });

    remove_lower_priority_duplicates(in_vcf_hdr, vcf_records, aux_records);

    std::vector<bcf1_t*> small_variants, svs;
    divide_variants_by_svsize(in_vcf_hdr, vcf_records, small_variants, svs);

    // make_multiallelic(in_vcf_hdr, small_variants);
    // add_star_alleles_for_overlapping(in_vcf_hdr, small_variants);
    make_multiallelic_new(in_vcf_hdr, small_variants);

    htsFile* out_smvars_vcf_file = bcf_open(out_smvars_vcf_fname.c_str(), "wz");
    if (bcf_hdr_write(out_smvars_vcf_file, in_vcf_hdr) != 0) {
        throw std::runtime_error("Failed to write the VCF header to " + out_smvars_vcf_fname + ".");
    }
    for (bcf1_t* small_variant : small_variants) {
        if (bcf_write(out_smvars_vcf_file, in_vcf_hdr, small_variant) != 0) {
            throw std::runtime_error("Failed to write VCF record to " + out_smvars_vcf_fname);
        }
    }
    hts_close(out_smvars_vcf_file);
    index_vcf(out_smvars_vcf_fname);

    htsFile* out_stvars_vcf_file = bcf_open(out_stvars_vcf_fname.c_str(), "wz");
    if (bcf_hdr_write(out_stvars_vcf_file, in_vcf_hdr) != 0) {
        throw std::runtime_error("Failed to write the VCF header to " + out_stvars_vcf_fname + ".");
    }
    for (bcf1_t* sv : svs) {
        if (bcf_write(out_stvars_vcf_file, in_vcf_hdr, sv) != 0) {
            throw std::runtime_error("Failed to write VCF record to " + out_stvars_vcf_fname);
        }
    }
    hts_close(out_stvars_vcf_file);
    index_vcf(out_stvars_vcf_fname);

    for (bcf1_t* vcf_record : vcf_records) {
        bcf_destroy(vcf_record);
    }
    bcf_destroy(r);

    bcf_hdr_destroy(in_vcf_hdr);
    hts_close(in_vcf_file);
}
