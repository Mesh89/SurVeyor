#ifndef AUX_CONTEXT_ASSIGNMENT_H
#define AUX_CONTEXT_ASSIGNMENT_H

#include <array>
#include <cmath>

#include "genotype_dels.h"

std::vector<std::vector<uint8_t>> get_genotyping_consensus_qualities(const std::string& allele_seq, const std::vector<std::shared_ptr<bam1_t>>& reads, const std::vector<bool>& revcomp_read, std::vector<std::string>& seqs, std::vector<hts_pos_t>& read_start_offsets, std::vector<std::vector<uint8_t>>* original_quals_storage);

// TODO: Replace this temporary pass-1 AUX context validator when proper SNP genotyping is available.
inline void validate_aux_snp_contexts(evidence_map_t& evidence, const std::string& chr, const std::vector<sv_t*>& svs, chr_seqs_map_t& chr_seqs, stats_t& stats, config_t& config, const std::string& bam_fname, const std::string& reference_fname, StripedSmithWaterman::Aligner& aligner) {
    using group_key_t = std::pair<std::string, std::vector<std::string>>;
    std::map<group_key_t, std::vector<sv_t*>> candidates;
    for (sv_t* sv : svs) {
        if (sv->svtype() != "DEL") continue;
        std::vector<std::string> aux_indels;
        for (const auto& indel : sv->aux_indels) aux_indels.push_back(indel->unique_key(false));
        std::sort(aux_indels.begin(), aux_indels.end());
        candidates[{sv->unique_key(false), aux_indels}].push_back(sv);
    }

    struct pooled_read_t {
        std::vector<size_t> contexts;
        size_t best_context;
    };
    struct context_group_t {
        std::vector<sv_t*> contexts;
        std::vector<std::string> templates;
        std::vector<int> positions;
        std::unordered_map<std::string, pooled_read_t> reads;
        hts_pos_t alt_start, ref_bp1_end, ref_bp2_start, alt_end;
    };
    std::vector<context_group_t> groups;
    std::unordered_map<sv_t*, std::pair<size_t, size_t>> context_indexes;
    for (auto& entry : candidates) {
        auto& contexts = entry.second;
        if (contexts.size() < 2 || std::none_of(contexts.begin(), contexts.end(), [](sv_t* sv) { return !sv->aux_snps.empty(); })) continue;
        std::sort(contexts.begin(), contexts.end(), [](sv_t* a, sv_t* b) { return a->id < b->id; });
        context_group_t group;
        group.contexts = contexts;
        bool equal_lengths = true;
        for (sv_t* sv : contexts) {
            // Native template construction sorts AUX vectors; keep candidate records untouched.
            deletion_t context(sv->chr, sv->start, sv->end, sv->ins_seq, nullptr, nullptr, nullptr, nullptr);
            context.id = sv->id;
            context.aux_indels = sv->aux_indels;
            context.aux_snps = sv->aux_snps;
            auto sequences = build_del_read_alignment_sequences(&context, chr_seqs.get_seq(chr), chr_seqs.get_len(chr), sv->start+1, sv->end+1, stats);
            group.templates.push_back(sequences.alt_seq.get());
            if (group.templates.back().length() != group.templates.front().length()) equal_lengths = false;
            group.alt_start = sequences.targets.alt_start;
            group.ref_bp1_end = sequences.targets.ref_bp1_end;
            group.ref_bp2_start = sequences.targets.ref_bp2_start;
            group.alt_end = sequences.targets.alt_end;
        }
        if (!equal_lengths) continue;
        for (int pos = 0; pos < group.templates[0].length(); pos++) {
            if (std::any_of(group.templates.begin()+1, group.templates.end(), [&](const std::string& seq) { return seq[pos] != group.templates[0][pos]; })) group.positions.push_back(pos);
        }
        if (group.positions.empty()) continue;
        for (size_t i = 0; i < contexts.size(); i++) context_indexes[contexts[i]] = {groups.size(), i};
        groups.push_back(std::move(group));
    }
    if (groups.empty()) return;

    // Use the original global maxima for eligibility, before excluding any context.
    for (const auto& entry : evidence.read_alt_associations) {
        std::map<size_t, std::map<size_t, int>> scores_by_group;
        for (const auto& association : entry.second.associations) {
            auto index = context_indexes.find(association.sv);
            if (index == context_indexes.end()) continue;
            auto& scores = scores_by_group[index->second.first];
            auto score = scores.emplace(index->second.second, association.score);
            score.first->second = std::max(score.first->second, association.score);
        }
        for (const auto& group_scores : scores_by_group) {
            if (group_scores.second.size() < 2) continue;
            pooled_read_t read;
            int best_score = std::numeric_limits<int>::min();
            for (const auto& score : group_scores.second) {
                read.contexts.push_back(score.first);
                if (score.second > best_score) {
                    best_score = score.second;
                    read.best_context = score.first; // Contexts are sorted by candidate ID, retaining the first score tie.
                }
            }
            if (best_score == entry.second.max_score) groups[group_scores.first].reads.emplace(entry.first, std::move(read));
        }
    }

    std::unique_ptr<open_samFile_t> bam_file;
    std::unordered_set<sv_t*> excluded_contexts;
    const std::string bases = "ACGT";
    StripedSmithWaterman::Filter filter(true, true, 0, 32767);
    for (const auto& group : groups) {
        if (group.reads.empty()) continue;
        if (!bam_file) {
            bam_file.reset(new open_samFile_t(bam_fname));
            hts_set_fai_filename(bam_file->file, fai_path(reference_fname.c_str()));
        }
        std::string left = chr+":"+std::to_string(group.alt_start)+"-"+std::to_string(group.ref_bp1_end);
        std::string right = chr+":"+std::to_string(group.ref_bp2_start)+"-"+std::to_string(group.alt_end);
        char* regions[2] = {&left[0], &right[0]};
        hts_itr_t* iter = sam_itr_regarray(bam_file->idx, bam_file->header, regions, 2);
        std::shared_ptr<bam1_t> read(bam_init1(), bam_destroy1);
        std::vector<std::array<double, 10>> log_likelihoods(group.positions.size(), std::array<double, 10>{});
        std::vector<bool> observed(group.positions.size(), false);
        while (sam_itr_next(bam_file->file, iter, read.get()) >= 0) {
            if (is_unmapped(read.get()) || !is_primary(read.get())) continue;
            auto pooled = group.reads.find(read_name_with_suffix(read.get()));
            if (pooled == group.reads.end()) continue;
            std::string seq = get_sequence(read.get());
            std::vector<uint8_t> qualities(seq.length(), 255);
            for (size_t context : pooled->second.contexts) {
                std::vector<std::string> seqs;
                std::vector<hts_pos_t> offsets;
                auto calibrated = get_genotyping_consensus_qualities(group.templates[context], {read}, {false}, seqs, offsets, nullptr);
                for (size_t i = 0; i < qualities.size(); i++) qualities[i] = std::min(qualities[i], calibrated[0][i]);
            }
            const std::string& target = group.templates[pooled->second.best_context];
            StripedSmithWaterman::Alignment alignment;
            aligner.Align(seq.c_str(), target.c_str(), target.length(), filter, &alignment, 0);
            int ref_pos = alignment.ref_begin, query_pos = alignment.query_begin;
            for (uint32_t encoded_op : alignment.cigar) {
                int op = bam_cigar_op(encoded_op), len = bam_cigar_oplen(encoded_op);
                if (op == BAM_CSOFT_CLIP) continue; // query_begin already accounts for leading clipping.
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    auto position = std::lower_bound(group.positions.begin(), group.positions.end(), ref_pos);
                    for (; position != group.positions.end() && *position < ref_pos+len; position++) {
                        int query_index = query_pos+*position-ref_pos;
                        size_t base = bases.find(seq[query_index]);
                        if (base == std::string::npos) continue;
                        double error = std::min(0.75, std::pow(10.0, -qualities[query_index]/10.0));
                        if (error > config.max_seq_error) continue;
                        size_t p = position-group.positions.begin();
                        observed[p] = true;
                        double probabilities[4];
                        for (int a = 0; a < 4; a++) probabilities[a] = a == base ? 1-error : error/3;
                        int genotype = 0;
                        for (int a = 0; a < 4; a++) {
                            for (int b = a; b < 4; b++) log_likelihoods[p][genotype++] += std::log((probabilities[a]+probabilities[b])/2);
                        }
                    }
                }
                if (bam_cigar_type(op)&1) query_pos += len;
                if (bam_cigar_type(op)&2) ref_pos += len;
            }
        }
        hts_itr_destroy(iter);

        std::vector<unsigned> allowed(group.positions.size(), 0);
        for (size_t p = 0; p < group.positions.size(); p++) {
            if (!observed[p]) continue;
            double best = *std::max_element(log_likelihoods[p].begin(), log_likelihoods[p].end());
            int genotype = 0;
            for (int a = 0; a < 4; a++) {
                for (int b = a; b < 4; b++) {
                    if (best-log_likelihoods[p][genotype++] <= 1e-9) allowed[p] |= (1u<<a) | (1u<<b);
                }
            }
        }
        std::vector<bool> survives(group.contexts.size(), true);
        for (size_t i = 0; i < group.contexts.size(); i++) {
            for (size_t p = 0; p < group.positions.size(); p++) {
                if (!observed[p]) continue;
                size_t base = bases.find(group.templates[i][group.positions[p]]);
                if (base == std::string::npos || !(allowed[p] & (1u<<base))) survives[i] = false;
            }
        }
        if (std::any_of(survives.begin(), survives.end(), [](bool value) { return value; })) {
            for (size_t i = 0; i < group.contexts.size(); i++) if (!survives[i]) excluded_contexts.insert(group.contexts[i]);
        }
    }
    for (auto& entry : evidence.read_alt_associations) {
        for (auto& association : entry.second.associations) if (excluded_contexts.count(association.sv)) association.cached_owner = false;
    }
}

#endif
