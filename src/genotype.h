#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <algorithm>
#include <cstring>
#include <fstream>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "../libs/ssw_cpp.h"
#include "extend_1sr_consensus.h"
#include "htslib/sam.h"
#include "sam_utils.h"
#include "types.h"
#include "var_utils.h"
#include "vcf_utils.h"

constexpr double MIN_EPR = 0.05;

// query and reference must be NUL-terminated. In exact-match mode, strstr uses
// the terminator to bound its search; reference_len is still used to validate
// that the returned match lies within the requested reference interval.
inline StripedSmithWaterman::Alignment align_fast(StripedSmithWaterman::Aligner& aligner, const char* query, const char* reference, int reference_len, 
    const StripedSmithWaterman::Filter& filter, bool require_exact_match = false) {

    StripedSmithWaterman::Alignment alignment;
    alignment.Clear();

    if (!require_exact_match) {
        aligner.Align(query, reference, reference_len, filter, &alignment, 0);
        return alignment;
    }

    if (query == nullptr || reference == nullptr || reference_len <= 0) {
        return alignment;
    }

    const int query_len = std::strlen(query);
    if (query_len <= 0 || query_len > reference_len) {
        return alignment;
    }

    // Queries and alleles are uppercase. Reject ambiguous query bases because
    // the genotyping SSW uses Nasmatch=false.
    if (std::strspn(query, "ACGT") != static_cast<size_t>(query_len)) {
        return alignment;
    }

    // Exact-match references must be NUL-terminated; constructed alleles are.
    const char* match = std::strstr(reference, query);
    if (match == nullptr || match - reference + query_len > reference_len) {
        return alignment;
    }

    const int ref_begin = match - reference;
    alignment.Clear();
    alignment.sw_score = query_len; // Genotyping aligners use match_score = 1.
    alignment.sw_score_next_best = 0;
    alignment.ref_end = ref_begin + query_len - 1;
    alignment.query_end = query_len - 1;
    alignment.ref_end_next_best = -1;
    alignment.mismatches = 0;

    if (filter.report_begin_position || filter.report_cigar) {
        alignment.ref_begin = ref_begin;
        alignment.query_begin = 0;
    } else {
        alignment.ref_begin = -1;
        alignment.query_begin = -1;
    }

    if (filter.report_cigar && alignment.sw_score >= filter.score_filter &&
            query_len - 1 <= filter.distance_filter) {
        alignment.cigar.push_back(bam_cigar_gen(query_len, BAM_CEQUAL));
        alignment.cigar_string = std::to_string(query_len) + "=";
    }

    return alignment;
}

std::string remove_svid_dup_suffix(const std::string& sv_id) {
    if (sv_id.size() > 4 && sv_id.substr(sv_id.size()-4) == "_DUP") {
        return sv_id.substr(0, sv_id.size()-4);
    }
    return sv_id;
}

struct bp_support_read_t {
    std::string read_name;
    int64_t mapq, mate_mapq;
    std::string seq;
    bool mate_is_reverse;
    hts_pos_t mate_pos;
    hts_pos_t mate_endpos;
    bool is_first_in_pair;

    bp_support_read_t() : mapq(0), mate_mapq(0), mate_is_reverse(false), mate_pos(-1), mate_endpos(-1), is_first_in_pair(false) {}

    explicit bp_support_read_t(bam1_t* read) :
        read_name(bam_get_qname(read)),
        mapq(read->core.qual),
        mate_mapq(get_mq(read)),
        seq(get_sequence(read)),
        mate_is_reverse(bam_is_mrev(read)),
        mate_pos(read->core.mpos),
        mate_endpos(get_mate_endpos(read)),
        is_first_in_pair(is_first_read(read))
     {}
};


std::string read_name_with_suffix(const bp_support_read_t& read) {
    return read.read_name + (read.is_first_in_pair ? "/1" : "/2");
}
std::string read_name_with_suffix(bam1_t* read) {
    return std::string(bam_get_qname(read)) + ((is_first_read(read)) ? "/1" : "/2");
}

struct evidence_logger_t {

    std::ofstream alt_reads_to_sv_associations, alt_pairs_to_sv_associations;
    std::mutex mtx;

    evidence_logger_t(const std::string& alt_reads_association_fname, const std::string& alt_pairs_association_fname) {
        alt_reads_to_sv_associations.open(alt_reads_association_fname);
        if (!alt_reads_to_sv_associations) {
            throw std::runtime_error("Unable to open file " + alt_reads_association_fname + ".");
        }
        alt_pairs_to_sv_associations.open(alt_pairs_association_fname);
        if (!alt_pairs_to_sv_associations) {
            throw std::runtime_error("Unable to open file " + alt_pairs_association_fname + ".");
        }
    }

    void log_pair_association(const std::string& sv_id, bam1_t* pair) {
        std::lock_guard<std::mutex> lock(mtx);
        alt_pairs_to_sv_associations << sv_id << " " << bam_get_qname(pair) << std::endl;
    }

    void log_reads_associations(std::string sv_id, int bp_n, std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<int>& scores) {
        std::lock_guard<std::mutex> lock(mtx);
        for (size_t i = 0; i < reads.size(); i++) {
            alt_reads_to_sv_associations << sv_id << " " << bp_n << " " << read_name_with_suffix(reads[i].get()) << " " << scores[i] << std::endl;
        }
    }
    void log_reads_associations(std::string sv_id, int bp_n, std::vector<bp_support_read_t>& reads, std::vector<int>& scores) {
        std::lock_guard<std::mutex> lock(mtx);
        for (size_t i = 0; i < reads.size(); i++) {
            alt_reads_to_sv_associations << sv_id << " " << bp_n << " " << read_name_with_suffix(reads[i]) << " " << scores[i] << std::endl;
        }
    }

    evidence_logger_t(const evidence_logger_t&) = delete;
    evidence_logger_t& operator=(const evidence_logger_t&) = delete;
};

struct evidence_map_t {
    std::unordered_map<std::string, int> read_to_hpid_map;

    // OAR/ORR consistency joins two observations that can arrive in either order on different genotyping threads:
    // a variant classifies an assigned-away read as ALT or REF, and the assigned variant later finds it consistent.
    std::mutex other_read_support_mtx;
    std::unordered_map<std::string, sv_t::other_read_info_t> assigned_consistent_reads;
    std::unordered_map<std::string, std::vector<std::pair<sv_t::sample_info_t*, int>>> oar_targets_by_read;
    std::unordered_map<std::string, std::vector<std::pair<sv_t::sample_info_t*, int>>> orr_targets_by_read;

    evidence_map_t() {}

    void load(const std::string& alt_reads_association_fname, const std::string& vcf_fname, config_t& config) {
        
        htsFile* vcf_file = bcf_open(vcf_fname.c_str(), "r");
        if (vcf_file == NULL) {
            throw std::runtime_error("Unable to open file " + vcf_fname + ".");
        }

        bcf_hdr_t* vcf_header = bcf_hdr_read(vcf_file);
        if (vcf_header == NULL) {
            throw std::runtime_error("Failed to read the VCF header.");
        }

        std::unordered_map<std::string, float> sv_epr_map;
        std::unordered_map<std::string, int> sv_hpid_map;

        bcf1_t* vcf_record = bcf_init();
        while (bcf_read(vcf_file, vcf_header, vcf_record) == 0) {
            bcf_unpack(vcf_record, BCF_UN_ALL);

            std::string id = vcf_record->d.id;
            id = remove_svid_dup_suffix(id);
            float epr = get_sv_epr(vcf_header, vcf_record);
            sv_epr_map[id] = epr;

            int* hpid_data = NULL;
            int hpid_len = 0;
            if (bcf_get_info_int32(vcf_header, vcf_record, "HPID", &hpid_data, &hpid_len) <= 0 || hpid_len == 0) {
                free(hpid_data);
                throw std::runtime_error("Missing HPID for SV " + id + ".");
            }
            sv_hpid_map[id] = hpid_data[0];
            free(hpid_data);
        }
        hts_close(vcf_file);
        bcf_hdr_destroy(vcf_header);
        bcf_destroy(vcf_record);

        std::ifstream alt_reads_association_fin(alt_reads_association_fname);
        if (!alt_reads_association_fin) {
            throw std::runtime_error("Unable to open file " + alt_reads_association_fname + ".");
        }
        std::string sv_id, read_name;
        int bp, score;

        struct best_assoc_t {
            bool passes_min_epr;
            int score;
            bool unique;
            int hpid;
            std::set<int> hpids;

            best_assoc_t() : passes_min_epr(false), score(0), unique(true), hpid(0), hpids() {}
            best_assoc_t(bool passes_min_epr, int score, int hpid) :
                passes_min_epr(passes_min_epr), score(score), unique(true), hpid(hpid), hpids{hpid} {}
        };

        // fields are: whether the best association passes MIN_EPR, score of
        // the best association found within that EPR class, whether the
        // association is unique or not, hpid of that association, and all hpids
        // with that best EPR class and score
        // note that unique means that there is a single best association to a HPID
        std::unordered_map<std::string, best_assoc_t> read_to_best_assoc_map;

        // For each read, find the best association. Furthermore, flag reads that have a "best" association to multiple SVs
        while (alt_reads_association_fin >> sv_id >> bp >> read_name >> score) {
            sv_id = remove_svid_dup_suffix(sv_id);
            bool passes_min_epr = sv_epr_map[sv_id] >= MIN_EPR;
            int hpid = sv_hpid_map[sv_id];
            if (!read_to_best_assoc_map.count(read_name)) {
                read_to_best_assoc_map[read_name] = best_assoc_t(passes_min_epr, score, hpid);
            } else {
                best_assoc_t& curr_best_assoc = read_to_best_assoc_map[read_name];
                if (passes_min_epr > curr_best_assoc.passes_min_epr || 
                    (passes_min_epr == curr_best_assoc.passes_min_epr && score > curr_best_assoc.score)) {
                    curr_best_assoc.passes_min_epr = passes_min_epr;
                    curr_best_assoc.score = score;
                    curr_best_assoc.hpid = hpid;
                    curr_best_assoc.unique = true;
                    curr_best_assoc.hpids.clear();
                    curr_best_assoc.hpids.insert(hpid);
                } else if (passes_min_epr == curr_best_assoc.passes_min_epr && score == curr_best_assoc.score && hpid != curr_best_assoc.hpid) {
                    curr_best_assoc.unique = false;
                    curr_best_assoc.hpids.insert(hpid);
                }
            }
        }


        alt_reads_association_fin.clear();
        alt_reads_association_fin.seekg(0, std::ios::beg);

        // Now, for each SV, we calculate U, the number of reads with a unique best
        // association to that SV.
        // Then, we assign reads as follows:
        // - we assign each read with a unique best association to that SV
        // - when a read has multiple equally good best associations, we keep the best two
        //   in terms of U, breaking U ties by EPR. Then, we assign the read to one of
        //   the two randomly, with probability proportional to the number U of each SV
        //   (e.g., if SV1 has U=6 and SV2 has U=4, then the read is assigned to SV1 with probability 0.6 and to SV2 with probability 0.4)

        using sv_with_bp_t = std::pair<std::string, int>;
        std::unordered_map<int, int> hpid_to_U_map;
        for (const auto& kv : read_to_best_assoc_map) {
            const best_assoc_t& best_assoc = kv.second;
            if (best_assoc.unique) {
                hpid_to_U_map[best_assoc.hpid] += 1;
            }
        }
        std::unordered_map<std::string, std::vector<sv_with_bp_t>> read_to_multiple_svs; // only for reads with multiple best associations
        while (alt_reads_association_fin >> sv_id >> bp >> read_name >> score) {
            sv_id = remove_svid_dup_suffix(sv_id);
            const best_assoc_t& best_assoc = read_to_best_assoc_map[read_name];
            bool passes_min_epr = sv_epr_map[sv_id] >= MIN_EPR;
            if (passes_min_epr == best_assoc.passes_min_epr && score == best_assoc.score) {
                if (best_assoc.unique) {
                    read_to_hpid_map[read_name] = sv_hpid_map[sv_id];
                } else {
                    read_to_multiple_svs[read_name].push_back({sv_id, bp});
                }
            }
        }

        // Now, assign non-uniquely assigned reads
        std::mt19937 gen(config.seed);
        std::vector<std::string> read_names_with_multiple_svs;
        read_names_with_multiple_svs.reserve(read_to_multiple_svs.size());
        for (const auto& kv : read_to_multiple_svs) {
            read_names_with_multiple_svs.push_back(kv.first);
        }
        std::sort(read_names_with_multiple_svs.begin(), read_names_with_multiple_svs.end());
        for (const std::string& read_name : read_names_with_multiple_svs) {
            std::vector<sv_with_bp_t>& sv_w_bps = read_to_multiple_svs[read_name];
            std::sort(sv_w_bps.begin(), sv_w_bps.end());
            // Find top two SVs, prioritizing variants above MIN_EPR, then U, then EPR.
            std::vector<std::pair<std::tuple<bool, int, float>, sv_with_bp_t>> sv_U_vec;
            for (const auto& sv_w_bp : sv_w_bps) {
                int U = hpid_to_U_map[sv_hpid_map[sv_w_bp.first]];
                float epr = sv_epr_map[sv_w_bp.first];
                bool passes_min_epr = epr >= MIN_EPR;
                sv_U_vec.push_back({{passes_min_epr, U, epr}, sv_w_bp});
            }
            std::sort(sv_U_vec.begin(), sv_U_vec.end(), std::greater<std::pair<std::tuple<bool, int, float>, sv_with_bp_t>>());
            sv_with_bp_t sv1 = sv_U_vec[0].second;

            // pick best SV that has a different HPID than sv1 (if any)
            size_t sv2_idx = 1;
            while (sv2_idx < sv_U_vec.size() && sv_hpid_map[sv_U_vec[sv2_idx].second.first] == sv_hpid_map[sv1.first]) {
                sv2_idx++;
            }
            if (sv2_idx == sv_U_vec.size()) {
                read_to_hpid_map[read_name] = sv_hpid_map[sv1.first];
                continue;
            }
            sv_with_bp_t sv2 = sv_U_vec[sv2_idx].second;

            int U1 = std::get<1>(sv_U_vec[0].first);
            int U2 = std::get<1>(sv_U_vec[sv2_idx].first);
            bool sv1_passes_min_epr = std::get<0>(sv_U_vec[0].first);
            bool sv2_passes_min_epr = std::get<0>(sv_U_vec[sv2_idx].first);

            if (U1 < 3) U1 = 0; // we require a minimum number of 3 uniquely assigned reads
            if (U2 < 3) U2 = 0; // we require a minimum number of 3 uniquely assigned reads
            int total_U = U1 + U2;
            if (sv1_passes_min_epr != sv2_passes_min_epr) {
                read_to_hpid_map[read_name] = sv_hpid_map[sv1.first];
            } else if (total_U == 0) {
                // assign to the highest-ranked SV when both candidates have no usable U
                read_to_hpid_map[read_name] = sv_hpid_map[sv1.first];
            } else {
                std::uniform_int_distribution<> dis(1, total_U);
                int r = dis(gen);
                if (r <= U1) {
                    read_to_hpid_map[read_name] = sv_hpid_map[sv1.first];
                } else {
                    read_to_hpid_map[read_name] = sv_hpid_map[sv2.first];
                }
            }

        }
    }

    bool is_read_assigned_to_different_sv(bam1_t* read, sv_t* sv) {
        std::string read_name = read_name_with_suffix(read);
        if (!read_to_hpid_map.count(read_name)) return false;
        return read_to_hpid_map[read_name] != sv->hpid;
    }

    // Check assignment for reads represented without their original BAM record, preserving the /1 or /2 suffix in the lookup key.
    bool is_read_assigned_to_different_sv(const bp_support_read_t& read, sv_t* sv) {
        std::string read_name = read_name_with_suffix(read);
        if (!read_to_hpid_map.count(read_name)) return false;
        return read_to_hpid_map[read_name] != sv->hpid;
    }

    bool is_read_assigned_to_this_sv(const bp_support_read_t& read, sv_t* sv) {
        std::string read_name = read_name_with_suffix(read);
        if (!read_to_hpid_map.count(read_name)) return false;
        return read_to_hpid_map[read_name] == sv->hpid;
    }

    // Explicitly insert one OAR*C read, or upgrade its flags if it was already inserted. The caller must hold other_read_support_mtx.
    void insert_or_update_oar_consistent_read(sv_t::sample_info_t& sample_info, int bp_n, const std::string& read_name, const sv_t::other_read_info_t& assigned_read_info) {
        if (sample_info.too_deep) return;

        std::unordered_map<std::string, sv_t::other_read_info_t>* consistent_reads;
        if (bp_n == 1) {
            consistent_reads = &sample_info.oar_bp1_consistent_reads;
        } else if (bp_n == 2) {
            consistent_reads = &sample_info.oar_bp2_consistent_reads;
        } else {
            throw std::runtime_error("Invalid OAR breakpoint number " + std::to_string(bp_n) + ".");
        }

        sv_t::other_read_info_t read_info = assigned_read_info;
        auto existing_read = consistent_reads->find(read_name);
        if (existing_read == consistent_reads->end()) {
            consistent_reads->emplace(read_name, read_info);
        } else {
            existing_read->second.hq |= read_info.hq;
            existing_read->second.exact |= read_info.exact;
        }
    }

    // Explicitly insert one ORR*C read, or upgrade its flags if it was already inserted. The caller must hold other_read_support_mtx.
    void insert_or_update_orr_consistent_read(sv_t::sample_info_t& sample_info, int bp_n, const std::string& read_name, const sv_t::other_read_info_t& assigned_read_info) {
        if (sample_info.too_deep) return;

        std::unordered_map<std::string, sv_t::other_read_info_t>* consistent_reads;
        if (bp_n == 1) {
            consistent_reads = &sample_info.orr_bp1_consistent_reads;
        } else if (bp_n == 2) {
            consistent_reads = &sample_info.orr_bp2_consistent_reads;
        } else {
            throw std::runtime_error("Invalid ORR breakpoint number " + std::to_string(bp_n) + ".");
        }

        sv_t::other_read_info_t read_info = assigned_read_info;
        auto existing_read = consistent_reads->find(read_name);
        if (existing_read == consistent_reads->end()) {
            consistent_reads->emplace(read_name, read_info);
        } else {
            existing_read->second.hq |= read_info.hq;
            existing_read->second.exact |= read_info.exact;
        }
    }

    // Remember that an assigned-away read supports this variant's ALT breakpoint. If the assigned
    // variant has already reported consistency, populate OAR*C immediately; otherwise record_assigned_read_consistency does it later.
    void register_oar_support(sv_t::sample_info_t& sample_info, int bp_n, const std::string& read_name) {
        std::lock_guard<std::mutex> lock(other_read_support_mtx);

        auto existing_targets = oar_targets_by_read.find(read_name);
        if (existing_targets == oar_targets_by_read.end()) {
            std::vector<std::pair<sv_t::sample_info_t*, int>> targets = {{&sample_info, bp_n}};
            oar_targets_by_read.emplace(read_name, targets);
        } else {
            bool already_registered = false;
            for (const auto& target : existing_targets->second) {
                if (target.first == &sample_info && target.second == bp_n) {
                    already_registered = true;
                    break;
                }
            }
            if (!already_registered) existing_targets->second.push_back({&sample_info, bp_n});
        }

        auto assigned_read_it = assigned_consistent_reads.find(read_name);
        if (assigned_read_it != assigned_consistent_reads.end()) insert_or_update_oar_consistent_read(sample_info, bp_n, read_name, assigned_read_it->second);
    }

    // Convenience overloads preserve the /1 or /2 suffix used as the read key.
    void register_oar_support(sv_t::sample_info_t& sample_info, int bp_n, bam1_t* read) {
        register_oar_support(sample_info, bp_n, read_name_with_suffix(read));
    }

    void register_oar_support(sv_t::sample_info_t& sample_info, int bp_n, const bp_support_read_t& read) {
        register_oar_support(sample_info, bp_n, read_name_with_suffix(read));
    }

    // Remember that an assigned-away read supports this variant's REF breakpoint. If the assigned
    // variant has already reported consistency, populate ORR*C immediately; otherwise record_assigned_read_consistency does it later.
    void register_orr_support(sv_t::sample_info_t& sample_info, int bp_n, const std::string& read_name) {
        std::lock_guard<std::mutex> lock(other_read_support_mtx);

        // Insert the read into the list of targets for this read name, if not already present.
        auto existing_targets = orr_targets_by_read.find(read_name);
        if (existing_targets == orr_targets_by_read.end()) {
            std::vector<std::pair<sv_t::sample_info_t*, int>> targets = {{&sample_info, bp_n}};
            orr_targets_by_read.emplace(read_name, targets);
        } else {
            bool already_registered = false;
            for (const auto& target : existing_targets->second) {
                if (target.first == &sample_info && target.second == bp_n) {
                    already_registered = true;
                    break;
                }
            }
            if (!already_registered) existing_targets->second.push_back({&sample_info, bp_n});
        }

        // If the assigned variant has already reported consistency, insert or update ORR*C immediately.
        auto assigned_read_it = assigned_consistent_reads.find(read_name);
        if (assigned_read_it != assigned_consistent_reads.end()) insert_or_update_orr_consistent_read(sample_info, bp_n, read_name, assigned_read_it->second);
    }

    // Convenience overloads preserve the /1 or /2 suffix used as the read key.
    void register_orr_support(sv_t::sample_info_t& sample_info, int bp_n, bam1_t* read) {
        register_orr_support(sample_info, bp_n, read_name_with_suffix(read));
    }

    void register_orr_support(sv_t::sample_info_t& sample_info, int bp_n, const bp_support_read_t& read) {
        register_orr_support(sample_info, bp_n, read_name_with_suffix(read));
    }

    // Record consistency with the variant that received this read, then propagate it to every variant
    // where the read was classified as ALT or REF. hq and exact are OR-upgraded across duplicate records for CHQ and E.
    void record_assigned_read_consistency(const std::string& read_name, bool hq, bool exact) {
        std::lock_guard<std::mutex> lock(other_read_support_mtx);
        sv_t::other_read_info_t assigned_read_info;
        assigned_read_info.hq = hq;
        assigned_read_info.exact = exact;
        auto existing_read = assigned_consistent_reads.find(read_name);
        if (existing_read == assigned_consistent_reads.end()) {
            assigned_consistent_reads.emplace(read_name, assigned_read_info);
        } else {
            existing_read->second.hq |= assigned_read_info.hq;
            existing_read->second.exact |= assigned_read_info.exact;
            assigned_read_info = existing_read->second;
        }
        auto oar_targets_it = oar_targets_by_read.find(read_name);
        if (oar_targets_it != oar_targets_by_read.end()) {
            for (const auto& target : oar_targets_it->second) insert_or_update_oar_consistent_read(*target.first, target.second, read_name, assigned_read_info);
        }
        auto orr_targets_it = orr_targets_by_read.find(read_name);
        if (orr_targets_it != orr_targets_by_read.end()) {
            for (const auto& target : orr_targets_it->second) insert_or_update_orr_consistent_read(*target.first, target.second, read_name, assigned_read_info);
        }
    }

    // Convenience overloads preserve the /1 or /2 suffix used as the read key.
    void record_assigned_read_consistency(bam1_t* read, bool hq, bool exact) {
        record_assigned_read_consistency(read_name_with_suffix(read), hq, exact);
    }

    void record_assigned_read_consistency(const bp_support_read_t& read, bool hq, bool exact) {
        record_assigned_read_consistency(read_name_with_suffix(read), hq, exact);
    }

    // Mark a variant as too deep and clear all OAR/ORR counts and consistency information while
    // holding the join mutex. Later consistency observations will see too_deep and leave the fields empty.
    void clear_other_read_support_for_too_deep(sv_t::sample_info_t& sample_info) {
        std::lock_guard<std::mutex> lock(other_read_support_mtx);
        sample_info.too_deep = true;
        sample_info.oar_bp1_reads = 0;
        sample_info.oar_bp2_reads = 0;
        sample_info.orr_bp1_reads = 0;
        sample_info.orr_bp2_reads = 0;
        sample_info.oar_bp1_consistent_reads.clear();
        sample_info.oar_bp2_consistent_reads.clear();
        sample_info.orr_bp1_consistent_reads.clear();
        sample_info.orr_bp2_consistent_reads.clear();
    }

};

std::vector<std::string> gen_consensus_seqs(std::string ref_seq, std::vector<std::string>& seqs);
std::vector<std::shared_ptr<bam1_t>> gen_consensus_and_find_consistent_seqs_subset(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, 
    std::vector<bool> revcomp_read, std::string& consensus_seq, double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read);
std::vector<std::shared_ptr<bam1_t>> find_seqs_consistent_with_ref_seq(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, 
    double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<bp_support_read_t>& consistent_reads,
    std::vector<bool>& is_exact_read, double consistent_avg_score, double consistent_stddev_score);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<std::shared_ptr<bam1_t>>& consistent_reads, 
    std::vector<bool>& is_exact_read, double consistent_avg_score, double consistent_stddev_score);

void read_mates(int contig_id);
void release_mates(int contig_id);

IntervalTree<ext_read_t*> get_candidate_reads_for_extension_itree(std::string contig_name, hts_pos_t contig_len, std::vector<hts_pair_pos_t> target_ivals, open_samFile_t* bam_file,
                                                                  std::vector<ext_read_t*>& candidate_reads_for_extension);


// Given a sequence alt_seq, a series of sequences ref_seqs and a read length read_len,
// return all the positions in alt_seq where a read of length read_len that is not present 
// in any of the ref_seqs starts.
std::vector<hts_pos_t> get_diff_reads_expected_positions(std::vector<char*>& ref_seqs, std::vector<hts_pos_t>& ref_lens, char* alt_seq, hts_pos_t alt_len, int read_len) {
    std::vector<hts_pos_t> positions;
    if (read_len > alt_len || read_len <= 0) {
        return positions;
    }
    for (hts_pos_t i = 0; i <= alt_len - read_len; i++) {
        char* read_begin = alt_seq + i;
        char* read_end = read_begin + read_len;
        bool found_in_ref = false;
        for (size_t j = 0; j < ref_seqs.size(); j++) {
            char* ref_begin = ref_seqs[j];
            char* ref_end = ref_begin + ref_lens[j];
            auto it = std::search(ref_begin, ref_end, read_begin, read_end);
            if (it != ref_end) {
                found_in_ref = true;
                break;
            }
        }
        if (!found_in_ref) {
            positions.push_back(i);
        }
    }
    return positions;
}

std::vector<int> get_consistent_reads_start_positions(std::vector<std::shared_ptr<bam1_t>>& consistent_reads,
    std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<int>& start_positions) {
    std::vector<int> consistent_positions;
    std::unordered_set<bam1_t*> consistent_set;
    for (auto& r : consistent_reads) {
        consistent_set.insert(r.get());
    }

    for (int i = 0; i < reads.size(); i++) {
        if (consistent_set.count(reads[i].get())) {
            consistent_positions.push_back(start_positions[i]);
        }
    }
    return consistent_positions;
}

// Given a vector of observed positions (with duplicates allowed) and the number of distinct valid positions
// observable, compute the occupancy ratio as the number of distinct observed positions divided 
// by expected number of observed positions under uniform distribution of observations across valid positions.
// We can assume that all positions in 'positions' are valid.
// If U is the number of unique observed positions, E[U] =  k * (1 - (1 - 1/k)^n), where k is valid_positions and n is positions.size()
double occ_ratio(std::vector<int>& positions, int valid_positions) {
    if (valid_positions <= 0) return sv_t::sample_info_t::NOT_COMPUTED;

    int n = positions.size();
    if (n == 0) return 1.0;

    double k = valid_positions;

    std::unordered_set<int> unique_positions(positions.begin(), positions.end());
    double U = unique_positions.size();
    double EU = k * (1.0 - std::exp((double)n * std::log1p(-1.0 / k)));

    if (EU <= 0.0) return sv_t::sample_info_t::NOT_COMPUTED;
    return U / EU;
}

#endif // GENOTYPE_H
