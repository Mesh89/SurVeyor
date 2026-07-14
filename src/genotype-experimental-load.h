#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <fstream>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "extend_1sr_consensus.h"
#include "htslib/sam.h"
#include "sam_utils.h"
#include "types.h"
#include "vcf_utils.h"

std::string read_name_with_suffix(bam1_t* read) {
    std::string read_name = bam_get_qname(read);
    read_name += (is_first_read(read) ? "/1" : "/2");
    return read_name;
}

std::string remove_svid_dup_suffix(const std::string& sv_id) {
    if (sv_id.size() > 4 && sv_id.substr(sv_id.size()-4) == "_DUP") {
        return sv_id.substr(0, sv_id.size()-4);
    }
    return sv_id;
}

struct evidence_logger_t {

    std::ofstream alt_reads_to_sv_associations, alt_pairs_to_sv_associations;
    std::mutex mtx;

    evidence_logger_t(const std::string& workdir) {
        alt_reads_to_sv_associations.open(workdir + "/alt_reads_to_sv_associations.txt");
        alt_pairs_to_sv_associations.open(workdir + "/alt_pairs_to_sv_associations.txt");
    }

    void log_pair_association(const std::string& sv_id, bam1_t* pair) {
        std::lock_guard<std::mutex> lock(mtx);
        alt_pairs_to_sv_associations << sv_id << " " << bam_get_qname(pair) << std::endl;
    }

    void log_reads_associations(std::string sv_id, int bp_n, std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<int>& scores) {
        std::lock_guard<std::mutex> lock(mtx);
        for (size_t i = 0; i < reads.size(); i++) {
            bam1_t* read = reads[i].get();
            alt_reads_to_sv_associations << sv_id << " " << bp_n << " " << read_name_with_suffix(read) << " " << scores[i] << std::endl;
        }
    }

    evidence_logger_t(const evidence_logger_t&) = delete;
    evidence_logger_t& operator=(const evidence_logger_t&) = delete;
};

struct evidence_map_t {
    std::unordered_map<std::string, std::string> read_to_sv_map;
    std::unordered_map<std::string, int> read_to_bp_map;
    std::unordered_map<std::string, std::vector<std::pair<std::string, int>>> read_to_non_chosen_svs_map;

    evidence_map_t() {}

    void load(std::string workdir, std::string vcf_fname, config_t& config) {
        
        htsFile* vcf_file = bcf_open(vcf_fname.c_str(), "r");
        if (vcf_file == NULL) {
            throw std::runtime_error("Unable to open file " + vcf_fname + ".");
        }

        bcf_hdr_t* vcf_header = bcf_hdr_read(vcf_file);
        if (vcf_header == NULL) {
            throw std::runtime_error("Failed to read the VCF header.");
        }

        std::unordered_map<std::string, float> sv_epr_map;

        bcf1_t* vcf_record = bcf_init();
        while (bcf_read(vcf_file, vcf_header, vcf_record) == 0) {
            bcf_unpack(vcf_record, BCF_UN_ALL);

            std::string id = vcf_record->d.id;
            id = remove_svid_dup_suffix(id);
            float epr = get_sv_epr(vcf_header, vcf_record);
            sv_epr_map[id] = epr;
        }
        hts_close(vcf_file);
        bcf_hdr_destroy(vcf_header);
        bcf_destroy(vcf_record);

        std::string alt_reads_association_fname = workdir + "/alt_reads_to_sv_associations.txt";
        std::ifstream alt_reads_association_fin(alt_reads_association_fname);
        std::string sv_id, read_name;
        int bp, score;

        // tuple fields are: score of the first best association found, sv_id of that association, 
        // whether the association is unique or not
        std::unordered_map<std::string, std::tuple<int, std::string, bool>> read_to_best_assoc_map;

        // For each read, find the best association. Furthermore, flag reads that have a "best" association to multiple SVs
        while (alt_reads_association_fin >> sv_id >> bp >> read_name >> score) {
            sv_id = remove_svid_dup_suffix(sv_id);
            if (!read_to_best_assoc_map.count(read_name)) {
                read_to_best_assoc_map[read_name] = {score, sv_id, true};
            } else {
                auto& curr_best_assoc = read_to_best_assoc_map[read_name];
                int& best_score = std::get<0>(curr_best_assoc);
                std::string& best_sv_id = std::get<1>(curr_best_assoc);
                bool& unique = std::get<2>(curr_best_assoc);
                if (score > best_score) {
                    best_score = score;
                    best_sv_id = sv_id;
                    unique = true;
                } else if (score == best_score && sv_id != best_sv_id) {
                    unique = false;
                }
            }
        }


        alt_reads_association_fin.clear();
        alt_reads_association_fin.seekg(0, std::ios::beg);

        // Now, for each SV, we calculate a score S that is the sum of scores of reads that have A (not necessarily unique) 
        // best association to that SV
        // Then, we assign reads as follows: 
        // - we assign each read with a unique best association to that SV (let us call their number U)
        // - when a read has multiple equally good best associations, we keep the best two in terms of S. Then,
        //   we assign the read to one of the two randomly, with probability proportional to the number U of each SV
        //   (e.g., if SV1 has U=6 and SV2 has U=4, then the read is assigned to SV1 with probability 0.6 and to SV2 with probability 0.4)

        using sv_with_bp_t = std::pair<std::string, int>;
        std::unordered_map<std::string, int> sv_to_S_map;
        std::unordered_map<std::string, int> sv_to_U_map;
        std::unordered_map<std::string, std::vector<sv_with_bp_t>> read_to_multiple_svs; // only for reads with multiple best associations
        while (alt_reads_association_fin >> sv_id >> bp >> read_name >> score) {
            sv_id = remove_svid_dup_suffix(sv_id);
            auto& best_assoc = read_to_best_assoc_map[read_name];
            int best_score = std::get<0>(best_assoc);
            if (score == best_score) {
                sv_to_S_map[sv_id] += score;
            }

            if (score == best_score) {
                bool unique = std::get<2>(best_assoc);
                if (unique) {
                    sv_to_U_map[sv_id] += 1;
                    read_to_sv_map[read_name] = sv_id;
                    read_to_bp_map[read_name] = bp;
                } else {
                    read_to_multiple_svs[read_name].push_back({sv_id, bp});
                }
            } else {
                read_to_non_chosen_svs_map[read_name].push_back({sv_id, bp});
            }
        }

        // Now, assign non-uniquely assigned reads
        std::mt19937 gen(config.seed);
        for (auto& kv : read_to_multiple_svs) {
            std::string read_name = kv.first;
            std::vector<sv_with_bp_t>& sv_w_bps = kv.second;
            // find top two svs in terms of S
            std::vector<std::pair<int, sv_with_bp_t>> sv_S_vec;
            for (const auto& sv_w_bp : sv_w_bps) {
                sv_S_vec.push_back({sv_to_S_map[sv_w_bp.first], sv_w_bp});
            }
            std::sort(sv_S_vec.begin(), sv_S_vec.end(), std::greater<std::pair<int, sv_with_bp_t>>());
            sv_with_bp_t sv1 = sv_S_vec[0].second;
            sv_with_bp_t sv2 = sv_S_vec[1].second;

            int U1 = sv_to_U_map[sv1.first];
            int U2 = sv_to_U_map[sv2.first];
            if (U1 < 3) U1 = 0; // we require a minimum number of 3 uniquely assigned reads
            if (U2 < 3) U2 = 0;
            int total_U = U1 + U2;
            if (total_U == 0) {
                // assign to the highest EPR SV
                float epr1 = sv_epr_map[sv1.first];
                float epr2 = sv_epr_map[sv2.first];
                if (epr1 >= epr2) {
                    read_to_sv_map[read_name] = sv1.first;
                    read_to_bp_map[read_name] = sv1.second;
                } else {
                    read_to_sv_map[read_name] = sv2.first;
                    read_to_bp_map[read_name] = sv2.second;
                }
            } else {
                std::uniform_int_distribution<> dis(1, total_U);
                int r = dis(gen);
                if (r <= U1) {
                    read_to_sv_map[read_name] = sv1.first;
                    read_to_bp_map[read_name] = sv1.second;
                } else {
                    read_to_sv_map[read_name] = sv2.first;
                    read_to_bp_map[read_name] = sv2.second;
                }
            }

            for (const auto& sv_w_bp : sv_w_bps) {
                if (sv_w_bp.first != read_to_sv_map[read_name]) {
                    read_to_non_chosen_svs_map[read_name].push_back(sv_w_bp);
                }
            }
        }
    }

    bool is_read_assigned_to_different_sv(bam1_t* read, std::string sv_id) {
        std::string read_name = read_name_with_suffix(read);
        if (!read_to_sv_map.count(read_name)) return false;
        sv_id = remove_svid_dup_suffix(sv_id);
        return read_to_sv_map[read_name] != sv_id;
    }

    std::vector<std::pair<std::string, int>> get_non_chosen_svs_for_read(bam1_t* read) {
        std::string read_name = read_name_with_suffix(read);
        if (!read_to_non_chosen_svs_map.count(read_name)) return {};
        return read_to_non_chosen_svs_map[read_name];
    }
};

std::mutex orc_mtx;

void increase_orc_supp(sv_t::sample_info_t& sample_info, int bp_n, bool hq) {
    std::lock_guard<std::mutex> lock(orc_mtx);
    if (bp_n == 1) {
        sample_info.assigned_to_other_sv_bp1_consistent++;
        if (hq) {
            sample_info.assigned_to_other_sv_bp1_consistent_highmq++;
        }
    } else if (bp_n == 2) {
        sample_info.assigned_to_other_sv_bp2_consistent++;
        if (hq) {
            sample_info.assigned_to_other_sv_bp2_consistent_highmq++;
        }
    }
}

void increase_orc(std::unordered_map<std::string, std::shared_ptr<sv_t>>& sv_map, std::string sv_id, int bp_n, bool hq) {
    if (!sv_map.count(sv_id)) return;
    increase_orc_supp(sv_map[sv_id]->sample_info, bp_n, hq);

    sv_id += "_DUP";
    if (!sv_map.count(sv_id)) return;
    increase_orc_supp(sv_map[sv_id]->sample_info, bp_n, hq);
}

std::vector<std::string> gen_consensus_seqs(std::string ref_seq, std::vector<std::string>& seqs);
std::vector<std::shared_ptr<bam1_t>> gen_consensus_and_find_consistent_seqs_subset(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, 
    std::vector<bool> revcomp_read, std::string& consensus_seq, double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read);
std::vector<std::shared_ptr<bam1_t>> find_seqs_consistent_with_ref_seq(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads, 
    double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, int n_reads, std::vector<std::shared_ptr<bam1_t>>& consistent_reads, 
    std::vector<bool>& is_exact_read, double consistent_avg_score, double consistent_stddev_score);

void read_mates(int contig_id);
void release_mates(int contig_id);

IntervalTree<ext_read_t*> get_candidate_reads_for_extension_itree(std::string contig_name, hts_pos_t contig_len, std::vector<hts_pair_pos_t> target_ivals, open_samFile_t* bam_file,
                                                                  std::vector<ext_read_t*>& candidate_reads_for_extension);


char* generate_haplotype_left(char* chrom_seq, hts_pos_t hap_end, hts_pos_t hap_len, 
    std::vector<std::shared_ptr<sv_t>>& aux_indels, std::vector<snp_t>& aux_snps);
char* generate_haplotype_right(char* chrom_seq, hts_pos_t chrom_len, hts_pos_t hap_start, hts_pos_t hap_len,
    std::vector<std::shared_ptr<sv_t>>& aux_indels, std::vector<snp_t>& aux_snps);

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
    if (n == 0) return sv_t::sample_info_t::NOT_COMPUTED;

    double k = valid_positions;

    std::unordered_set<int> unique_positions(positions.begin(), positions.end());
    double U = unique_positions.size();
    double EU = k * (1.0 - std::exp((double)n * std::log1p(-1.0 / k)));

    if (EU <= 0.0) return sv_t::sample_info_t::NOT_COMPUTED;
    return U / EU;
}

#endif // GENOTYPE_H
