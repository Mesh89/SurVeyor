#ifndef QUALITY_READ_SCORING_H
#define QUALITY_READ_SCORING_H

#include <array>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "../libs/ssw_selected.h"
#include "hp_mismatch_rate_thresholds.h"

inline hp_tail_quality_model_t& quality_read_scoring_model() {
    static hp_tail_quality_model_t model;
    return model;
}

inline int quality_read_base_code(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

struct quality_read_scoring_tables_t {
    enum { alphabet_size = 89 };
    std::array<int, 21> match, mismatch;
    std::array<char, 84> query_letters;
    std::array<int8_t, 256> translation;
    std::array<int8_t, alphabet_size*alphabet_size> matrix{};

    quality_read_scoring_tables_t() {
        double e37 = std::pow(10.0, -3.7);
        for (int c = 0; c <= 20; c++) {
            double e = std::pow(10.0, -2*c/10.0);
            match[c] = e >= 0.75 ? 0 : std::lround(10*std::log(4*(1-e))/std::log(4*(1-e37)));
            mismatch[c] = e >= 0.75 ? 0 : std::lround(-40*std::log(4*e/3)/std::log(4*e37/3));
            for (int b = 0; b < 4; b++) {
                int q = 5+4*c+b;
                for (int t = 0; t < 4; t++) matrix[t*alphabet_size+q] = matrix[q*alphabet_size+t] = t == b ? match[c] : mismatch[c];
            }
        }
        translation.fill(4);
        for (int b = 0; b < 4; b++) {
            translation["ACGT"[b]] = b;
            translation["acgt"[b]] = b;
        }
        int k = 0;
        for (int ch = 0x21; ch <= 0x7e; ch++) {
            if (std::strchr("ACGTNacgtn", ch)) continue;
            query_letters[k] = ch;
            translation[ch] = 5+k;
            k++;
        }
    }
};

inline const quality_read_scoring_tables_t& quality_read_scoring_tables() {
    static const quality_read_scoring_tables_t tables;
    return tables;
}

struct quality_read_encoding_t {
    std::string query;
    std::vector<int> classes; // BAM-stored positions; -1 for ambiguous bases.
    int match_sum = 0;
};

// reversed describes seq relative to the BAM-stored sequence, including palindromes.
inline quality_read_encoding_t encode_quality_read(bam1_t* read, const std::string& seq, bool reversed, char* contig_seq, hts_pos_t contig_len, const config_t& config) {
    std::string stored_seq = get_sequence(read);
    thread_local hp_tail_quality_table_t quality_cache;
    std::vector<uint8_t> quals = recalibrate_clip_read_qualities(read, config, quality_read_scoring_model(), quality_cache, nullptr, contig_seq, contig_len);
    const quality_read_scoring_tables_t& tables = quality_read_scoring_tables();
    quality_read_encoding_t encoded;
    encoded.classes.resize(seq.length(), -1);
    encoded.query.resize(seq.length(), 'N');
    for (int i = 0; i < seq.length(); i++) {
        if (quality_read_base_code(stored_seq[i]) == 4) continue;
        int c = std::min(20, (int)(quals[i]/2.0+0.5));
        encoded.classes[i] = c;
        encoded.match_sum += tables.match[c];
    }
    for (int i = 0; i < seq.length(); i++) {
        int c = encoded.classes[reversed ? seq.length()-i-1 : i];
        int b = quality_read_base_code(seq[i]);
        if (c >= 0 && b < 4) encoded.query[i] = tables.query_letters[4*c+b];
    }
    return encoded;
}

struct quality_read_aligner_t {
    StripedSmithWaterman::Aligner aligner;
    std::string target;

    quality_read_aligner_t() : aligner(quality_read_scoring_tables().matrix.data(), quality_read_scoring_tables_t::alphabet_size, quality_read_scoring_tables().translation.data(), 256) {
        aligner.SetGapPenalty(60, 10);
        aligner.SetTargetAlphabetSize(5);
    }
};

inline void quality_align(const quality_read_encoding_t& encoded, const char* target, int target_len, const StripedSmithWaterman::Filter& filter, StripedSmithWaterman::Alignment* alignment) {
    thread_local quality_read_aligner_t context;
    context.target.assign(target, target_len);
    // Keep target codes in 0..4: IUPAC bases and flank separators can be query letters.
    for (char& base : context.target) if (quality_read_base_code(base) == 4) base = 'N';
    context.aligner.Align(encoded.query.c_str(), context.target.c_str(), target_len, filter, alignment, 0);
}

inline StripedSmithWaterman::Alignment quality_align_fast(const char* raw_seq, const quality_read_encoding_t& encoded, const char* target, int target_len, const StripedSmithWaterman::Filter& filter, bool require_exact_match = false) {
    StripedSmithWaterman::Alignment alignment;
    alignment.Clear();
    if (raw_seq == nullptr || target == nullptr || target_len <= 0) return alignment;
    int query_len = std::strlen(raw_seq);
    if (query_len > 0 && query_len <= target_len && std::strspn(raw_seq, "ACGT") == static_cast<size_t>(query_len)) {
        const char* match = std::strstr(target, raw_seq);
        if (match != nullptr && match-target+query_len <= target_len) {
            int ref_begin = match-target;
            alignment.sw_score = encoded.match_sum;
            alignment.sw_score_next_best = 0;
            alignment.ref_end = ref_begin+query_len-1;
            alignment.query_end = query_len-1;
            alignment.ref_end_next_best = -1;
            alignment.mismatches = 0;
            if (filter.report_begin_position || filter.report_cigar) {
                alignment.ref_begin = ref_begin;
                alignment.query_begin = 0;
            } else {
                alignment.ref_begin = -1;
                alignment.query_begin = -1;
            }
            if (filter.report_cigar && alignment.sw_score >= filter.score_filter && query_len-1 <= filter.distance_filter) {
                alignment.cigar.push_back(bam_cigar_gen(query_len, BAM_CEQUAL));
                alignment.cigar_string = std::to_string(query_len)+"=";
            }
            return alignment;
        }
    }
    if (!require_exact_match) quality_align(encoded, target, target_len, filter, &alignment);
    return alignment;
}

// Walk the original CIGAR in BAM orientation, using the same availability rules
// as score_original_bam_placement and the recalibrated classes computed once above.
inline bool quality_score_original_bam_placement(const bam1_t* original_read, const char* contig_seq, hts_pos_t contig_len, const quality_read_encoding_t& encoded, int64_t& score) {
    if (!original_read || !contig_seq || (original_read->core.flag & BAM_FUNMAP) || original_read->core.tid < 0 || original_read->core.pos < 0 || original_read->core.pos >= contig_len || original_read->core.n_cigar == 0 || original_read->core.l_qseq <= 0) return false;
    const quality_read_scoring_tables_t& tables = quality_read_scoring_tables();
    const uint32_t* cigar = bam_get_cigar(original_read);
    const uint8_t* seq = bam_get_seq(original_read);
    hts_pos_t ref_pos = original_read->core.pos;
    int64_t query_pos = 0, placement_score = 0;
    bool has_aligned_bases = false, has_alignment_ops = false, trailing_clip = false;
    int previous_op = -1;
    for (uint32_t i = 0; i < original_read->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int64_t len = bam_cigar_oplen(cigar[i]);
        if (len == 0) return false;
        if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            if (has_alignment_ops) trailing_clip = true;
            if (op == BAM_CSOFT_CLIP) {
                if (len > original_read->core.l_qseq-query_pos) return false;
                query_pos += len;
            }
        } else {
            if (trailing_clip) return false;
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                if (len > original_read->core.l_qseq-query_pos || len > contig_len-ref_pos) return false;
                for (int64_t j = 0; j < len; j++) {
                    int read_base = bam_seqi(seq, query_pos+j);
                    if (read_base == 0) return false; // '=' omits the actual query base.
                    int c = encoded.classes[query_pos+j];
                    int query_code = quality_read_base_code(seq_nt16_str[read_base]);
                    int ref_code = quality_read_base_code(contig_seq[ref_pos+j]);
                    if (c >= 0 && query_code < 4 && ref_code < 4) placement_score += query_code == ref_code ? tables.match[c] : tables.mismatch[c];
                }
                query_pos += len;
                ref_pos += len;
                has_aligned_bases = true;
            } else if (op == BAM_CINS || op == BAM_CDEL) {
                if (op == BAM_CINS) {
                    if (len > original_read->core.l_qseq-query_pos) return false;
                    query_pos += len;
                } else {
                    if (len > contig_len-ref_pos) return false;
                    ref_pos += len;
                }
                placement_score -= previous_op == op ? 10*len : 60+10*(len-1);
            } else {
                return false; // N, P and B have no genotyping affine-gap meaning.
            }
            has_alignment_ops = true;
        }
        previous_op = op;
    }
    if (!has_aligned_bases || query_pos != original_read->core.l_qseq) return false;
    score = placement_score;
    return true;
}

inline bool quality_passes_original_bam_alt_veto(const bam1_t* original_read, const char* contig_seq, hts_pos_t contig_len, const quality_read_encoding_t& encoded, int alt_score) {
    int64_t original_score;
    return !quality_score_original_bam_placement(original_read, contig_seq, contig_len, encoded, original_score) || alt_score > original_score;
}

#endif
