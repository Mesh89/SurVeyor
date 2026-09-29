#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "htslib/vcf.h"

#include "vcf_utils.h"

enum class gt_state_t {
    UNKNOWN,
    HOM_REF,
    NON_REF
};

struct allele_gt_status_t {
    bool has_hom_ref = false;
    bool has_non_ref = false;
};

gt_state_t get_gt_state(bcf_hdr_t* hdr, bcf1_t* record) {
    int32_t* gt = nullptr;
    int ngt = 0;
    int ret = bcf_get_genotypes(hdr, record, &gt, &ngt);
    if (ret <= 0) {
        free(gt);
        return gt_state_t::UNKNOWN;
    }

    bool saw_called_allele = false;
    bool saw_alt_allele = false;
    for (int i = 0; i < ret; i++) {
        if (gt[i] == bcf_int32_vector_end) break;
        if (bcf_gt_is_missing(gt[i])) {
            free(gt);
            return gt_state_t::UNKNOWN;
        }
        saw_called_allele = true;
        if (bcf_gt_allele(gt[i]) > 0) saw_alt_allele = true;
    }
    free(gt);

    if (!saw_called_allele) return gt_state_t::UNKNOWN;
    return saw_alt_allele ? gt_state_t::NON_REF : gt_state_t::HOM_REF;
}

void append_indel_atom_keys(std::vector<std::string>& keys, const std::shared_ptr<sv_t>& sv) {
    if (sv->svtype() != "INS" && sv->svtype() != "DEL") return;

    // Match expand_aux_haplotypes.cpp: a replacement is represented by a
    // deletion atom and an insertion atom in AUX_INDELS.
    if (sv->start != sv->end && !sv->ins_seq.empty()) {
        deletion_t deletion(sv->chr, sv->start, sv->end, "", nullptr, nullptr, nullptr, nullptr);
        insertion_t insertion(sv->chr, sv->start, sv->start, sv->ins_seq, nullptr, nullptr, nullptr, nullptr);
        keys.push_back(deletion.unique_key(false));
        keys.push_back(insertion.unique_key(false));
        return;
    }

    keys.push_back(sv->unique_key(false));
}

std::unordered_map<std::string, allele_gt_status_t> load_indel_gt_statuses(
        const std::string& input_vcf) {
    htsFile* input = bcf_open(input_vcf.c_str(), "r");
    if (input == nullptr) {
        throw std::runtime_error("Unable to open input VCF " + input_vcf + ".");
    }
    bcf_hdr_t* hdr = bcf_hdr_read(input);
    if (hdr == nullptr) {
        hts_close(input);
        throw std::runtime_error("Unable to read header from " + input_vcf + ".");
    }
    if (bcf_hdr_nsamples(hdr) != 1) {
        bcf_hdr_destroy(hdr);
        hts_close(input);
        throw std::runtime_error("prune_hom_ref_aux_indels requires exactly one VCF sample.");
    }

    std::unordered_map<std::string, allele_gt_status_t> statuses;
    bcf1_t* record = bcf_init();
    while (bcf_read(input, hdr, record) == 0) {
        gt_state_t gt_state = get_gt_state(hdr, record);
        if (gt_state == gt_state_t::UNKNOWN) continue;

        std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, record);
        if (sv == nullptr || (sv->svtype() != "INS" && sv->svtype() != "DEL")) continue;

        std::vector<std::string> atom_keys;
        append_indel_atom_keys(atom_keys, sv);
        for (const std::string& key : atom_keys) {
            allele_gt_status_t& status = statuses[key];
            if (gt_state == gt_state_t::HOM_REF) status.has_hom_ref = true;
            if (gt_state == gt_state_t::NON_REF) status.has_non_ref = true;
        }
    }

    bcf_destroy(record);
    bcf_hdr_destroy(hdr);
    hts_close(input);
    return statuses;
}

std::string aux_indels_string(const std::vector<std::shared_ptr<sv_t>>& aux_indels) {
    std::string value;
    for (const auto& indel : aux_indels) {
        if (!value.empty()) value += ",";
        value += std::to_string(indel->start + 1) + ":";
        value += std::to_string(indel->end + 1) + ":";
        value += indel->ins_seq;
    }
    return value;
}

size_t prune_aux_indels(std::shared_ptr<sv_t>& sv, const std::unordered_map<std::string, allele_gt_status_t>& statuses) {
    std::vector<std::shared_ptr<sv_t>> retained;
    retained.reserve(sv->aux_indels.size());
    size_t removed = 0;
    for (const auto& aux_indel : sv->aux_indels) {
        const std::string key = aux_indel->unique_key(false);
        auto status_it = statuses.find(key);
        bool remove = status_it != statuses.end() &&
            status_it->second.has_hom_ref && !status_it->second.has_non_ref;
        if (remove) {
            removed++;
        } else {
            retained.push_back(aux_indel);
        }
    }
    if (removed > 0) sv->aux_indels = std::move(retained);
    return removed;
}

bool has_dup_suffix(const std::string& id) {
    return id.size() > 4 && id.compare(id.size() - 4, 4, "_DUP") == 0;
}

std::unordered_map<std::string, int> load_post_pruning_insertion_hpids(
        const std::string& input_vcf, const std::unordered_map<std::string, allele_gt_status_t>& statuses) {
    htsFile* input = bcf_open(input_vcf.c_str(), "r");
    if (input == nullptr) {
        throw std::runtime_error("Unable to open input VCF " + input_vcf + ".");
    }
    bcf_hdr_t* hdr = bcf_hdr_read(input);
    if (hdr == nullptr) {
        hts_close(input);
        throw std::runtime_error("Unable to read header from " + input_vcf + ".");
    }

    std::unordered_map<std::string, int> insertion_hpids;
    bcf1_t* record = bcf_init();
    while (bcf_read(input, hdr, record) == 0) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, record);
        if (sv == nullptr || sv->svtype() != "INS" || has_dup_suffix(sv->id)) continue;

        size_t removed = prune_aux_indels(sv, statuses);
        int post_pruning_hpid = removed > 0 ? make_hpid(*sv) : sv->hpid;
        auto result = insertion_hpids.emplace(sv->id, post_pruning_hpid);
        if (!result.second && result.first->second != post_pruning_hpid) {
            bcf_destroy(record);
            bcf_hdr_destroy(hdr);
            hts_close(input);
            throw std::runtime_error("Insertion ID " + sv->id + " has conflicting post-pruning HPIDs.");
        }
    }

    bcf_destroy(record);
    bcf_hdr_destroy(hdr);
    hts_close(input);
    return insertion_hpids;
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: prune_hom_ref_aux_indels <input.vcf.gz> <output.vcf.gz> [--positive-only]\n";
        return 1;
    }

    const std::string input_vcf = argv[1];
    const std::string output_vcf = argv[2];
    bool positive_only = false;
    if (argc == 4) {
        if (std::string(argv[3]) != "--positive-only") {
            std::cerr << "Error: unknown option " << argv[3] << ".\n";
            return 1;
        }
        positive_only = true;
    }

    try {
        const auto statuses = load_indel_gt_statuses(input_vcf);
        const auto insertion_hpids = load_post_pruning_insertion_hpids(input_vcf, statuses);

        htsFile* input = bcf_open(input_vcf.c_str(), "r");
        if (input == nullptr) {
            throw std::runtime_error("Unable to open input VCF " + input_vcf + ".");
        }
        bcf_hdr_t* hdr = bcf_hdr_read(input);
        if (hdr == nullptr) {
            hts_close(input);
            throw std::runtime_error("Unable to read header from " + input_vcf + ".");
        }

        htsFile* output = bcf_open(output_vcf.c_str(), "wz");
        if (output == nullptr) {
            bcf_hdr_destroy(hdr);
            hts_close(input);
            throw std::runtime_error("Unable to open output VCF " + output_vcf + ".");
        }
        if (bcf_hdr_write(output, hdr) != 0) {
            hts_close(output);
            bcf_hdr_destroy(hdr);
            hts_close(input);
            throw std::runtime_error("Unable to write header to " + output_vcf + ".");
        }

        size_t records_changed = 0;
        size_t aux_indels_removed = 0;
        size_t dup_hpids_updated = 0;
        size_t non_positive_records_filtered = 0;
        bcf1_t* record = bcf_init();
        while (bcf_read(input, hdr, record) == 0) {
            std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, record);
            if (sv != nullptr) {
                size_t removed = prune_aux_indels(sv, statuses);
                if (removed > 0) {
                    records_changed++;
                    aux_indels_removed += removed;

                    std::string aux_value = aux_indels_string(sv->aux_indels);
                    if (aux_value.empty()) {
                        bcf_update_info_string(hdr, record, "AUX_INDELS", nullptr);
                    } else {
                        bcf_update_info_string(hdr, record, "AUX_INDELS", aux_value.c_str());
                    }

                    sv->hpid = make_hpid(*sv);
                    bcf_update_info_int32(hdr, record, "HPID", &sv->hpid, 1);
                }

                if (sv->svtype() == "DUP" && has_dup_suffix(sv->id)) {
                    std::string origin_id = sv->id.substr(0, sv->id.size() - 4);
                    auto origin_it = insertion_hpids.find(origin_id);
                    if (origin_it == insertion_hpids.end()) {
                        bcf_destroy(record);
                        hts_close(output);
                        bcf_hdr_destroy(hdr);
                        hts_close(input);
                        throw std::runtime_error("Unable to find origin insertion " + origin_id + " for duplication " + sv->id + ".");
                    }
                    if (sv->hpid != origin_it->second) dup_hpids_updated++;
                    sv->hpid = origin_it->second;
                    bcf_update_info_int32(hdr, record, "HPID", &sv->hpid, 1);
                }
            }

            if (positive_only && get_gt_state(hdr, record) != gt_state_t::NON_REF) {
                non_positive_records_filtered++;
                continue;
            }

            if (bcf_write(output, hdr, record) != 0) {
                bcf_destroy(record);
                hts_close(output);
                bcf_hdr_destroy(hdr);
                hts_close(input);
                throw std::runtime_error("Unable to write record to " + output_vcf + ".");
            }
        }

        bcf_destroy(record);
        hts_close(output);
        bcf_hdr_destroy(hdr);
        hts_close(input);

        std::cerr << "records_changed\t" << records_changed << "\n";
        std::cerr << "aux_indels_removed\t" << aux_indels_removed << "\n";
        std::cerr << "dup_hpids_updated\t" << dup_hpids_updated << "\n";
        if (positive_only) {
            std::cerr << "non_positive_records_filtered\t" << non_positive_records_filtered << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
