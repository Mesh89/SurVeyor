#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <algorithm>
#include <cstring>
#include <fstream>
#include <list>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../libs/ssw_cpp.h"
#include "extend_1sr_consensus.h"
#include "htslib/sam.h"
#include "sam_utils.h"
#include "types.h"
#include "var_utils.h"
#include "vcf_utils.h"

constexpr double MIN_EPR = 0.05;
constexpr hts_pos_t GENOTYPE_CONSENSUS_EXTENSION = 500;

struct alignment_targets_t {
    char* alt_seq = NULL;
    int alt_len = 0;
    std::vector<char*> ref_seqs;
    std::vector<int> ref_lens;

    // The left flank is alt_seq[0:left_flank_end], while the right flank starts at right_flank_start.
    int left_flank_end = 0, right_flank_start = 0;
    int right_flank_end_offset = 1;
    bool has_left_split = true, has_right_split = true;

    char* left_independent_ref_seq = NULL;
    char* right_independent_ref_seq = NULL;
    int left_independent_ref_len = 0, right_independent_ref_len = 0;
};

inline int consensus_alignment_score(const StripedSmithWaterman::Alignment& alignment) {
    return alignment.query_end - alignment.query_begin - alignment.mismatches;
}

inline consensus_alignment_metrics_t score_consensus_alignment(const std::string& consensus_seq, const alignment_targets_t& targets, StripedSmithWaterman::Aligner& aligner) {
    consensus_alignment_metrics_t metrics;
    metrics.length = consensus_seq.length();

    if (consensus_seq.empty() || targets.alt_seq == NULL || targets.alt_len <= 0) return metrics;

    StripedSmithWaterman::Filter with_pos_and_cigar(true, true, 0, 32767);
    StripedSmithWaterman::Filter score_only(false, false, 0, 32767);

    StripedSmithWaterman::Alignment alt_alignment;
    alt_alignment.Clear();
    aligner.Align(consensus_seq.c_str(), targets.alt_seq, targets.alt_len, with_pos_and_cigar, &alt_alignment, 0);
    metrics.alt_score = consensus_alignment_score(alt_alignment);
    metrics.alt_ref_begin = alt_alignment.ref_begin;
    metrics.alt_ref_end = alt_alignment.ref_end;

    for (int i = 0; i < targets.ref_seqs.size() && i < targets.ref_lens.size(); i++) {
        if (targets.ref_seqs[i] == NULL || targets.ref_lens[i] <= 0) continue;
        StripedSmithWaterman::Alignment ref_alignment;
        ref_alignment.Clear();
        aligner.Align(consensus_seq.c_str(), targets.ref_seqs[i], targets.ref_lens[i], with_pos_and_cigar, &ref_alignment, 0);
        metrics.ref_score = std::max(metrics.ref_score, consensus_alignment_score(ref_alignment));
    }

    int left_ref_len = targets.has_left_split ? std::max(0, targets.left_flank_end-alt_alignment.ref_begin) : 0;
    int right_ref_len = targets.has_right_split ? std::max(0, alt_alignment.ref_end+targets.right_flank_end_offset-targets.right_flank_start) : 0;
    auto left = targets.has_left_split ? find_aln_prefix_score(alt_alignment.cigar, left_ref_len, 1, -4, -6, -1, true) : std::make_pair(0, 0);
    auto right = targets.has_right_split ? find_aln_suffix_score(alt_alignment.cigar, right_ref_len, 1, -4, -6, -1, true) : std::make_pair(0, 0);
    metrics.split_ref_lengths = {{left_ref_len, right_ref_len}};
    metrics.split_sizes = {{left.second, right.second}};
    metrics.split_scores = {{left.first, right.first}};

    if (left.second > 0 && targets.left_independent_ref_seq != NULL && targets.left_independent_ref_len > 0) {
        StripedSmithWaterman::Alignment alignment;
        alignment.Clear();
        const std::string query = consensus_seq.substr(0, left.second);
        aligner.Align(query.c_str(), targets.left_independent_ref_seq, targets.left_independent_ref_len, score_only, &alignment, 0);
        metrics.independent_ref_scores[0] = alignment.sw_score;
    }

    if (right.second > 0 && targets.right_independent_ref_seq != NULL && targets.right_independent_ref_len > 0) {
        StripedSmithWaterman::Alignment alignment;
        alignment.Clear();
        const std::string query = consensus_seq.substr(consensus_seq.length()-right.second);
        aligner.Align(query.c_str(), targets.right_independent_ref_seq, targets.right_independent_ref_len, score_only, &alignment, 0);
        metrics.independent_ref_scores[1] = alignment.sw_score;
    }

    return metrics;
}

// query and reference must be NUL-terminated. strstr uses the terminator to
// bound its search; reference_len is still used to validate that the returned
// match lies within the requested reference interval.
inline StripedSmithWaterman::Alignment align_fast(StripedSmithWaterman::Aligner& aligner, const char* query, const char* reference, int reference_len, 
    const StripedSmithWaterman::Filter& filter, bool require_exact_match = false) {

    StripedSmithWaterman::Alignment alignment;
    alignment.Clear();

    if (query == nullptr || reference == nullptr || reference_len <= 0) {
        return alignment;
    }

    const int query_len = std::strlen(query);
    if (query_len > 0 && query_len <= reference_len &&
            std::strspn(query, "ACGT") == static_cast<size_t>(query_len)) {
        const char* match = std::strstr(reference, query);
        if (match != nullptr && match - reference + query_len <= reference_len) {
            const int ref_begin = match - reference;
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
    }

    if (!require_exact_match) {
        aligner.Align(query, reference, reference_len, filter, &alignment, 0);
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

    static const size_t MAX_OPEN_READ_ASSOCIATION_FILES = 64;
    std::string alt_reads_association_dir;
    std::unordered_map<std::string, std::unique_ptr<std::ofstream>> alt_reads_to_sv_associations;
    std::list<std::string> alt_reads_association_lru;
    std::unordered_map<std::string, std::list<std::string>::iterator> alt_reads_association_lru_iters;
    std::ofstream alt_pairs_to_sv_associations;
    std::mutex mtx;

    evidence_logger_t(const std::string& alt_reads_association_dir, const std::string& alt_pairs_association_fname) : alt_reads_association_dir(alt_reads_association_dir) {
        alt_pairs_to_sv_associations.open(alt_pairs_association_fname);
        if (!alt_pairs_to_sv_associations) {
            throw std::runtime_error("Unable to open file " + alt_pairs_association_fname + ".");
        }
    }

    std::ofstream& get_reads_association_file(const std::string& chr, const std::string& suffix) {
        std::string file_key = chr + suffix;
        auto file_it = alt_reads_to_sv_associations.find(file_key);
        if (file_it != alt_reads_to_sv_associations.end()) {
            alt_reads_association_lru.erase(alt_reads_association_lru_iters[file_key]);
            alt_reads_association_lru.push_front(file_key);
            alt_reads_association_lru_iters[file_key] = alt_reads_association_lru.begin();
            return *file_it->second;
        }
        if (alt_reads_to_sv_associations.size() == MAX_OPEN_READ_ASSOCIATION_FILES) {
            const std::string& evicted_chr = alt_reads_association_lru.back();
            alt_reads_to_sv_associations.erase(evicted_chr);
            alt_reads_association_lru_iters.erase(evicted_chr);
            alt_reads_association_lru.pop_back();
        }
        std::unique_ptr<std::ofstream> file(new std::ofstream(alt_reads_association_dir + "/" + file_key + ".txt", std::ios::app));
        if (!*file) throw std::runtime_error("Unable to open read association file for " + file_key + ".");
        std::ofstream& file_ref = *file;
        alt_reads_to_sv_associations.emplace(file_key, std::move(file));
        alt_reads_association_lru.push_front(file_key);
        alt_reads_association_lru_iters[file_key] = alt_reads_association_lru.begin();
        return file_ref;
    }

    void log_pair_association(const std::string& sv_id, bam1_t* pair) {
        std::lock_guard<std::mutex> lock(mtx);
        alt_pairs_to_sv_associations << sv_id << " " << bam_get_qname(pair) << std::endl;
    }

    void log_reads_associations(const std::string& chr, std::string sv_id, int bp_n, std::vector<std::shared_ptr<bam1_t>>& reads, std::vector<int>& scores, std::vector<int>& positions, int alt_idx = -1) {
        if (reads.empty()) return;
        std::lock_guard<std::mutex> lock(mtx);
        std::ofstream& alt_reads_to_sv_associations = get_reads_association_file(chr, ".alt");
        for (size_t i = 0; i < reads.size(); i++) {
            alt_reads_to_sv_associations << sv_id << " " << bp_n << " " << read_name_with_suffix(reads[i].get()) << " " << scores[i] << " " << positions[i] << " " << alt_idx << '\n';
        }
    }
    void log_reads_associations(const std::string& chr, std::string sv_id, int bp_n, std::vector<bp_support_read_t>& reads, std::vector<int>& scores, std::vector<int>& positions, int alt_idx = -1) {
        if (reads.empty()) return;
        std::lock_guard<std::mutex> lock(mtx);
        std::ofstream& alt_reads_to_sv_associations = get_reads_association_file(chr, ".alt");
        for (size_t i = 0; i < reads.size(); i++) {
            alt_reads_to_sv_associations << sv_id << " " << bp_n << " " << read_name_with_suffix(reads[i]) << " " << scores[i] << " " << positions[i] << " " << alt_idx << '\n';
        }
    }

    void log_ref_reads_associations(const std::string& chr, const std::string& sv_id, int bp_n, const std::vector<std::shared_ptr<bam1_t>>& reads) {
        if (reads.empty()) return;
        std::lock_guard<std::mutex> lock(mtx);
        std::ofstream& ref_reads_to_sv_associations = get_reads_association_file(chr, ".ref");
        for (const std::shared_ptr<bam1_t>& read : reads) ref_reads_to_sv_associations << read_name_with_suffix(read.get()) << " " << sv_id << " " << bp_n << '\n';
    }

    void log_ref_reads_associations(const std::string& chr, const std::string& sv_id, int bp_n, const std::vector<bp_support_read_t>& reads) {
        if (reads.empty()) return;
        std::lock_guard<std::mutex> lock(mtx);
        std::ofstream& ref_reads_to_sv_associations = get_reads_association_file(chr, ".ref");
        for (const bp_support_read_t& read : reads) ref_reads_to_sv_associations << read_name_with_suffix(read) << " " << sv_id << " " << bp_n << '\n';
    }

    void log_er_reads_associations(const std::string& chr, const std::string& sv_id, const std::vector<std::string>& read_names) {
        if (read_names.empty()) return;
        std::lock_guard<std::mutex> lock(mtx);
        std::ofstream& er_reads_to_sv_associations = get_reads_association_file(chr, ".er");
        for (const std::string& read_name : read_names) er_reads_to_sv_associations << read_name << " " << sv_id << '\n';
    }

    evidence_logger_t(const evidence_logger_t&) = delete;
    evidence_logger_t& operator=(const evidence_logger_t&) = delete;
};

struct evidence_map_t {
    struct read_alt_association_t {
        sv_t* sv;
        int bp;
        int pos;
        int alt_idx;

        read_alt_association_t(sv_t* sv, int bp, int pos, int alt_idx) : sv(sv), bp(bp), pos(pos), alt_idx(alt_idx) {}
    };

    struct read_alt_associations_t {
        std::vector<read_alt_association_t> associations;
        bool oar_recorded = false;
        bool orr_recorded = false;
        bool consistent = false;
        bool hq = false;
        bool exact = false;
    };

    struct read_ref_association_t {
        sv_t* sv;
        int bp;

        read_ref_association_t(sv_t* sv, int bp) : sv(sv), bp(bp) {}
    };

    struct read_assignment_t {
        int hpid;
        // Singleton HPIDs derive their sole VID from vid_idxs_by_hpid. Multi-VID
        // HPIDs reference one shared, interned winning set.
        uint32_t vid_set_idx;

        read_assignment_t() : hpid(0), vid_set_idx(std::numeric_limits<uint32_t>::max()) {}
        read_assignment_t(int hpid, uint32_t vid_set_idx) : hpid(hpid), vid_set_idx(vid_set_idx) {}
    };

    std::unordered_map<std::string, read_assignment_t> read_assignments;
    std::unordered_map<std::string, read_alt_associations_t> read_alt_associations;
    std::unordered_map<std::string, std::vector<read_ref_association_t>> read_ref_associations;
    std::unordered_map<std::string, std::vector<sv_t*>> read_er_associations;
    std::unordered_map<std::string, uint32_t> sv_id_to_idx;
    std::vector<std::string> sv_id_by_idx;
    std::unordered_map<int, std::vector<uint32_t>> vid_idxs_by_hpid;
    std::vector<std::vector<uint32_t>> interned_vid_sets;

    std::mutex other_read_support_mtx;

    evidence_map_t() {}

    void load(const std::string& alt_reads_association_fname, const std::string& ref_reads_association_fname, const std::string& er_reads_association_fname, const std::vector<sv_t*>& svs, config_t& config) {
        read_assignments.clear();
        read_alt_associations.clear();
        read_ref_associations.clear();
        read_er_associations.clear();
        sv_id_to_idx.clear();
        sv_id_by_idx.clear();
        vid_idxs_by_hpid.clear();
        interned_vid_sets.clear();
        std::unordered_map<std::string, float> sv_epr_map;
        std::unordered_map<std::string, int> sv_hpid_map;
        std::unordered_map<std::string, sv_t*> sv_by_exact_id;
        for (sv_t* sv : svs) {
            std::string id = remove_svid_dup_suffix(sv->id);
            sv_epr_map[id] = sv->sample_info.epr;
            sv_hpid_map[id] = sv->hpid;
            sv_by_exact_id[sv->id] = sv;
        }

        for (const auto& entry : sv_hpid_map) {
            uint32_t vid_idx = sv_id_by_idx.size();
            sv_id_to_idx.emplace(entry.first, vid_idx);
            sv_id_by_idx.push_back(entry.first);
            vid_idxs_by_hpid[entry.second].push_back(vid_idx);
        }
        for (auto& entry : vid_idxs_by_hpid) {
            std::sort(entry.second.begin(), entry.second.end());
        }

        std::map<std::vector<uint32_t>, uint32_t> vid_set_to_idx;
        auto intern_vid_set = [&](std::vector<uint32_t> vid_idxs) {
            std::sort(vid_idxs.begin(), vid_idxs.end());
            vid_idxs.erase(std::unique(vid_idxs.begin(), vid_idxs.end()), vid_idxs.end());
            auto existing = vid_set_to_idx.find(vid_idxs);
            if (existing != vid_set_to_idx.end()) return existing->second;
            uint32_t set_idx = interned_vid_sets.size();
            interned_vid_sets.push_back(vid_idxs);
            vid_set_to_idx.emplace(std::move(vid_idxs), set_idx);
            return set_idx;
        };

        std::string sv_id, read_name;
        int bp, score, pos, alt_idx;

        if (!ref_reads_association_fname.empty()) {
            std::ifstream ref_reads_association_fin(ref_reads_association_fname);
            if (!ref_reads_association_fin) throw std::runtime_error("Unable to open file " + ref_reads_association_fname + ".");
            while (ref_reads_association_fin >> read_name >> sv_id >> bp) read_ref_associations[read_name].push_back(read_ref_association_t(sv_by_exact_id.at(sv_id), bp));
        }

        if (!er_reads_association_fname.empty()) {
            std::ifstream er_reads_association_fin(er_reads_association_fname);
            if (!er_reads_association_fin) throw std::runtime_error("Unable to open file " + er_reads_association_fname + ".");
            while (er_reads_association_fin >> read_name >> sv_id) read_er_associations[read_name].push_back(sv_by_exact_id.at(sv_id));
        }

        if (alt_reads_association_fname.empty()) return;
        std::ifstream alt_reads_association_fin(alt_reads_association_fname);
        if (!alt_reads_association_fin) throw std::runtime_error("Unable to open file " + alt_reads_association_fname + ".");

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
        while (alt_reads_association_fin >> sv_id >> bp >> read_name >> score >> pos >> alt_idx) {
            std::string exact_sv_id = sv_id;
            read_alt_associations[read_name].associations.push_back(read_alt_association_t(sv_by_exact_id.at(exact_sv_id), bp, pos, alt_idx));
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
        // Only multi-VID HPIDs need a temporary per-read VID collection. For a
        // singleton HPID the VID is derived from the HPID without per-read storage.
        std::unordered_map<std::string, std::vector<uint32_t>> unique_read_vid_idxs;
        while (alt_reads_association_fin >> sv_id >> bp >> read_name >> score >> pos >> alt_idx) {
            sv_id = remove_svid_dup_suffix(sv_id);
            const best_assoc_t& best_assoc = read_to_best_assoc_map[read_name];
            bool passes_min_epr = sv_epr_map[sv_id] >= MIN_EPR;
            if (passes_min_epr == best_assoc.passes_min_epr && score == best_assoc.score) {
                if (best_assoc.unique) {
                    int hpid = sv_hpid_map[sv_id];
                    read_assignments[read_name] = read_assignment_t(hpid, std::numeric_limits<uint32_t>::max());
                    if (vid_idxs_by_hpid[hpid].size() > 1) {
                        unique_read_vid_idxs[read_name].push_back(sv_id_to_idx.at(sv_id));
                    }
                } else {
                    read_to_multiple_svs[read_name].push_back({sv_id, bp});
                }
            }
        }
        for (auto& entry : unique_read_vid_idxs) {
            read_assignments[entry.first].vid_set_idx = intern_vid_set(std::move(entry.second));
        }
        unique_read_vid_idxs.clear();

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
            int assigned_hpid;
            if (sv2_idx == sv_U_vec.size()) {
                assigned_hpid = sv_hpid_map[sv1.first];
            } else {
                sv_with_bp_t sv2 = sv_U_vec[sv2_idx].second;

                int U1 = std::get<1>(sv_U_vec[0].first);
                int U2 = std::get<1>(sv_U_vec[sv2_idx].first);
                bool sv1_passes_min_epr = std::get<0>(sv_U_vec[0].first);
                bool sv2_passes_min_epr = std::get<0>(sv_U_vec[sv2_idx].first);

                if (U1 < 3) U1 = 0; // we require a minimum number of 3 uniquely assigned reads
                if (U2 < 3) U2 = 0; // we require a minimum number of 3 uniquely assigned reads
                int total_U = U1 + U2;
                if (sv1_passes_min_epr != sv2_passes_min_epr) {
                    assigned_hpid = sv_hpid_map[sv1.first];
                } else if (total_U == 0) {
                    // assign to the highest-ranked SV when both candidates have no usable U
                    assigned_hpid = sv_hpid_map[sv1.first];
                } else {
                    std::uniform_int_distribution<> dis(1, total_U);
                    int r = dis(gen);
                    if (r <= U1) {
                        assigned_hpid = sv_hpid_map[sv1.first];
                    } else {
                        assigned_hpid = sv_hpid_map[sv2.first];
                    }
                }
            }

            uint32_t assigned_vid_set_idx = std::numeric_limits<uint32_t>::max();
            std::vector<uint32_t> assigned_vid_idxs;
            for (const auto& sv_w_bp : sv_w_bps) {
                if (sv_hpid_map[sv_w_bp.first] == assigned_hpid) {
                    assigned_vid_idxs.push_back(sv_id_to_idx.at(sv_w_bp.first));
                }
            }
            if (vid_idxs_by_hpid[assigned_hpid].size() > 1) {
                assigned_vid_set_idx = intern_vid_set(std::move(assigned_vid_idxs));
            }
            read_assignments[read_name] = read_assignment_t(assigned_hpid, assigned_vid_set_idx);
        }

    }

    const std::vector<uint32_t>* get_assigned_vid_idxs(const std::string& read_name) const {
        auto assignment_it = read_assignments.find(read_name);
        if (assignment_it == read_assignments.end()) return NULL;
        const read_assignment_t& assignment = assignment_it->second;
        if (assignment.vid_set_idx == std::numeric_limits<uint32_t>::max()) {
            auto hpid_vids_it = vid_idxs_by_hpid.find(assignment.hpid);
            return hpid_vids_it == vid_idxs_by_hpid.end() ? NULL : &hpid_vids_it->second;
        }
        return &interned_vid_sets[assignment.vid_set_idx];
    }

    bool is_read_assigned_to_different_sv(bam1_t* read, sv_t* sv) {
        std::string read_name = read_name_with_suffix(read);
        auto assignment_it = read_assignments.find(read_name);
        return assignment_it != read_assignments.end() && assignment_it->second.hpid != sv->hpid;
    }

    // Check assignment for reads represented without their original BAM record, preserving the /1 or /2 suffix in the lookup key.
    bool is_read_assigned_to_different_sv(const bp_support_read_t& read, sv_t* sv) {
        std::string read_name = read_name_with_suffix(read);
        auto assignment_it = read_assignments.find(read_name);
        return assignment_it != read_assignments.end() && assignment_it->second.hpid != sv->hpid;
    }

    bool is_read_assigned_to_this_sv(const bp_support_read_t& read, sv_t* sv) {
        std::string read_name = read_name_with_suffix(read);
        auto assignment_it = read_assignments.find(read_name);
        return assignment_it != read_assignments.end() && assignment_it->second.hpid == sv->hpid;
    }

    const std::vector<read_alt_association_t>& get_read_alt_associations(const std::string& read_name) const {
        static const std::vector<read_alt_association_t> no_associations;
        auto association_it = read_alt_associations.find(read_name);
        return association_it == read_alt_associations.end() ? no_associations : association_it->second.associations;
    }

    const std::vector<read_alt_association_t>& get_read_alt_associations(bam1_t* read) const {
        return get_read_alt_associations(read_name_with_suffix(read));
    }

    const std::vector<read_alt_association_t>& get_read_alt_associations(const bp_support_read_t& read) const {
        return get_read_alt_associations(read_name_with_suffix(read));
    }

    const std::vector<read_ref_association_t>& get_read_ref_associations(const std::string& read_name) const {
        static const std::vector<read_ref_association_t> no_associations;
        auto association_it = read_ref_associations.find(read_name);
        return association_it == read_ref_associations.end() ? no_associations : association_it->second;
    }

    const std::vector<read_ref_association_t>& get_read_ref_associations(bam1_t* read) const {
        return get_read_ref_associations(read_name_with_suffix(read));
    }

    const std::vector<read_ref_association_t>& get_read_ref_associations(const bp_support_read_t& read) const {
        return get_read_ref_associations(read_name_with_suffix(read));
    }

    bool is_read_ref(const std::string& read_name, sv_t* sv, int bp) const {
        for (const read_ref_association_t& association : get_read_ref_associations(read_name)) {
            if (association.sv == sv && association.bp == bp) return true;
        }
        return false;
    }

    bool is_read_ref(bam1_t* read, sv_t* sv, int bp) const {
        return is_read_ref(read_name_with_suffix(read), sv, bp);
    }

    bool is_read_ref(const bp_support_read_t& read, sv_t* sv, int bp) const {
        return is_read_ref(read_name_with_suffix(read), sv, bp);
    }

    bool is_read_er(const std::string& read_name, sv_t* sv) const {
        auto association_it = read_er_associations.find(read_name);
        if (association_it == read_er_associations.end()) return false;
        return std::find(association_it->second.begin(), association_it->second.end(), sv) != association_it->second.end();
    }

    bool is_read_er(bam1_t* read, sv_t* sv) const {
        return is_read_er(read_name_with_suffix(read), sv);
    }

    bool is_read_er(const bp_support_read_t& read, sv_t* sv) const {
        return is_read_er(read_name_with_suffix(read), sv);
    }

    int get_read_alt_pos(const std::string& read_name, sv_t* sv, int bp) const {
        for (const read_alt_association_t& association : get_read_alt_associations(read_name)) {
            if (association.sv == sv && association.bp == bp) return association.pos;
        }
        return -1;
    }

    int get_read_alt_pos(bam1_t* read, sv_t* sv, int bp) const {
        return get_read_alt_pos(read_name_with_suffix(read), sv, bp);
    }

    int get_read_alt_pos(const bp_support_read_t& read, sv_t* sv, int bp) const {
        return get_read_alt_pos(read_name_with_suffix(read), sv, bp);
    }

    const std::string& get_sv_id(uint32_t vid_idx) const {
        return sv_id_by_idx[vid_idx];
    }

    // Once a variant has classified an owned ALT read, propagate its complete OAR contribution to every other ALT association.
    void record_assigned_alt_read(sv_t* assigned_sv, const std::string& read_name, bool consistent, bool hq, bool exact) {
        auto assignment_it = read_assignments.find(read_name);
        if (assignment_it == read_assignments.end() || assignment_it->second.hpid != assigned_sv->hpid) return;
        const std::vector<uint32_t>* supporting_vid_idxs = get_assigned_vid_idxs(read_name);
        std::lock_guard<std::mutex> lock(other_read_support_mtx);
        auto associations_it = read_alt_associations.find(read_name);
        if (associations_it == read_alt_associations.end()) return;
        bool first_observation = !associations_it->second.oar_recorded;
        bool first_orr_observation = !associations_it->second.orr_recorded;
        bool newly_consistent = consistent && !associations_it->second.consistent;
        bool newly_hq = consistent && hq && !associations_it->second.hq;
        bool newly_exact = consistent && exact && !associations_it->second.exact;
        associations_it->second.oar_recorded = true;
        associations_it->second.orr_recorded = true;
        associations_it->second.consistent |= consistent;
        associations_it->second.hq |= consistent && hq;
        associations_it->second.exact |= consistent && exact;
        std::set<std::pair<sv_t*, int>> registered_targets;
        for (const read_alt_association_t& association : associations_it->second.associations) {
            if (association.sv->hpid == assigned_sv->hpid || association.bp < 1 || association.bp > 2) continue;
            if (!registered_targets.insert({association.sv, association.bp}).second) continue;
            sv_t::sample_info_t& target_sample_info = association.sv->sample_info;
            if (target_sample_info.too_deep) continue;
            if (first_observation) {
                std::unordered_map<std::string, int>* reads_by_vid;
                if (association.bp == 1) {
                    target_sample_info.oar_bp1_reads++;
                    target_sample_info.oar_bp1_reads_by_hpid[assignment_it->second.hpid]++;
                    reads_by_vid = &target_sample_info.oar_bp1_reads_by_vid;
                } else {
                    target_sample_info.oar_bp2_reads++;
                    target_sample_info.oar_bp2_reads_by_hpid[assignment_it->second.hpid]++;
                    reads_by_vid = &target_sample_info.oar_bp2_reads_by_vid;
                }
                if (supporting_vid_idxs) for (uint32_t vid_idx : *supporting_vid_idxs) (*reads_by_vid)[get_sv_id(vid_idx)]++;
            }
            if (association.bp == 1) {
                if (newly_consistent) target_sample_info.oar_bp1_consistent_reads++;
                if (newly_hq) target_sample_info.oar_bp1_consistent_hq_reads++;
                if (newly_exact) target_sample_info.oar_bp1_exact_reads++;
            } else {
                if (newly_consistent) target_sample_info.oar_bp2_consistent_reads++;
                if (newly_hq) target_sample_info.oar_bp2_consistent_hq_reads++;
                if (newly_exact) target_sample_info.oar_bp2_exact_reads++;
            }
        }
        std::set<std::pair<sv_t*, int>> registered_ref_targets;
        for (const read_ref_association_t& association : get_read_ref_associations(read_name)) {
            if (association.sv->hpid == assigned_sv->hpid || association.bp < 1 || association.bp > 2) continue;
            if (!registered_ref_targets.insert({association.sv, association.bp}).second) continue;
            sv_t::sample_info_t& target_sample_info = association.sv->sample_info;
            if (target_sample_info.too_deep) continue;
            if (association.bp == 1) {
                if (first_orr_observation) target_sample_info.orr_bp1_reads++;
                if (newly_consistent) target_sample_info.orr_bp1_consistent_reads++;
                if (newly_hq) target_sample_info.orr_bp1_consistent_hq_reads++;
                if (newly_exact) target_sample_info.orr_bp1_exact_reads++;
            } else {
                if (first_orr_observation) target_sample_info.orr_bp2_reads++;
                if (newly_consistent) target_sample_info.orr_bp2_consistent_reads++;
                if (newly_hq) target_sample_info.orr_bp2_consistent_hq_reads++;
                if (newly_exact) target_sample_info.orr_bp2_exact_reads++;
            }
        }
    }

    void record_assigned_alt_read(sv_t* assigned_sv, bam1_t* read, bool consistent, bool hq, bool exact) {
        record_assigned_alt_read(assigned_sv, read_name_with_suffix(read), consistent, hq, exact);
    }

    void record_assigned_alt_read(sv_t* assigned_sv, const bp_support_read_t& read, bool consistent, bool hq, bool exact) {
        record_assigned_alt_read(assigned_sv, read_name_with_suffix(read), consistent, hq, exact);
    }

    // Mark a variant as too deep and clear all OAR/ORR counts and consistency information while
    // holding the join mutex. Later consistency observations will see too_deep and leave the fields empty.
    void clear_other_read_support_for_too_deep(sv_t::sample_info_t& sample_info) {
        std::lock_guard<std::mutex> lock(other_read_support_mtx);
        sample_info.too_deep = true;
        sample_info.oar_bp1_reads = 0;
        sample_info.oar_bp2_reads = 0;
        sample_info.oar_bp1_consistent_reads = 0;
        sample_info.oar_bp2_consistent_reads = 0;
        sample_info.oar_bp1_consistent_hq_reads = 0;
        sample_info.oar_bp2_consistent_hq_reads = 0;
        sample_info.oar_bp1_exact_reads = 0;
        sample_info.oar_bp2_exact_reads = 0;
        sample_info.orr_bp1_reads = 0;
        sample_info.orr_bp2_reads = 0;
        sample_info.orr_bp1_consistent_reads = 0;
        sample_info.orr_bp2_consistent_reads = 0;
        sample_info.orr_bp1_consistent_hq_reads = 0;
        sample_info.orr_bp2_consistent_hq_reads = 0;
        sample_info.orr_bp1_exact_reads = 0;
        sample_info.orr_bp2_exact_reads = 0;
        sample_info.oar_bp1_reads_by_hpid.clear();
        sample_info.oar_bp2_reads_by_hpid.clear();
        sample_info.oar_bp1_reads_by_vid.clear();
        sample_info.oar_bp2_reads_by_vid.clear();
    }

};

std::vector<std::string> gen_consensus_seqs(std::string ref_seq, std::vector<std::string>& seqs);
std::vector<bool> gen_consensus_and_classify_seqs(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads,
    std::vector<bool> revcomp_read, std::string& consensus_seq, double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read);
std::vector<bool> classify_seqs_with_ref_seq(std::string ref_seq, std::vector<std::shared_ptr<bam1_t>>& reads,
    const std::vector<bool>& is_eligible_read, double& avg_score, double& stddev_score, std::vector<bool>& is_exact_read);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, std::vector<bp_support_read_t>& reads,
    std::vector<bool>& is_consistent_read, std::vector<bool>& is_exact_read,
    double consistent_avg_score, double consistent_stddev_score);

void set_bp_consensus_info(sv_t::bp_reads_info_t& bp_reads_info, std::vector<std::shared_ptr<bam1_t>>& reads,
    std::vector<bool>& is_consistent_read, std::vector<bool>& is_exact_read,
    double consistent_avg_score, double consistent_stddev_score);

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

// Return the span covered by all distinguishable reads on the ALT allele.
// Positions are read starts, so the rightmost read contributes read_len bases.
int get_max_feasible_alt_len(const std::vector<hts_pos_t>& distinguishable_positions, int read_len) {
    if (distinguishable_positions.empty() || read_len <= 0) {
        return sv_t::sample_info_t::NOT_COMPUTED;
    }
    auto bounds = std::minmax_element(distinguishable_positions.begin(), distinguishable_positions.end());
    return *bounds.second - *bounds.first + read_len;
}

std::vector<int> get_consistent_reads_start_positions(const std::vector<bool>& is_consistent_read,
    std::vector<int>& start_positions) {
    std::vector<int> consistent_positions;
    if (is_consistent_read.size() != start_positions.size()) {
        throw std::runtime_error("Read consistency and start-position vector size mismatch.");
    }
    for (int i = 0; i < is_consistent_read.size(); i++) {
        if (is_consistent_read[i]) {
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
