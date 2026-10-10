#ifndef ORIGINAL_BAM_SCORE_H
#define ORIGINAL_BAM_SCORE_H

#include <cstdint>
#include "htslib/sam.h"

// The genotyping SSW translation treats non-ACGT bases as N (and U as A).
// With Nasmatch=false and NvsNasmatch=false, even N against N is a mismatch.
inline int genotyping_base_code(char base) {
    switch (base) {
        case 'A': case 'a': case 'U': case 'u': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

// Score only the original BAM placement, in BAM sequence orientation. A null
// original_read explicitly represents synthetic evidence without that placement.
// Unsupported CIGARs or incomplete sequence/reference leave the score unavailable.
inline bool score_original_bam_placement(const bam1_t* original_read, const char* contig_seq, hts_pos_t contig_len, int64_t& score) {
    if (!original_read || !contig_seq || (original_read->core.flag & BAM_FUNMAP) || original_read->core.tid < 0 || original_read->core.pos < 0 || original_read->core.pos >= contig_len || original_read->core.n_cigar == 0 || original_read->core.l_qseq <= 0) return false;
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
                    int query_code = genotyping_base_code(seq_nt16_str[read_base]);
                    int ref_code = genotyping_base_code(contig_seq[ref_pos+j]);
                    placement_score += query_code < 4 && query_code == ref_code ? 1 : -4;
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
                // Adjacent operations of the same type belong to one gap run.
                placement_score -= previous_op == op ? len : 6+len-1;
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

// Additional ALT-only veto after the caller requires ALT > local_REF.
// Rescore the original CIGAR: +1 match, -4 mismatch (including ambiguous bases),
// -(6+length-1) per gap run, and 0 for clips; no AS tag or alignment search.
// A reliable score requires ALT > original_BAM (ties fail), giving
// ALT > max(local_REF, original_BAM). Unavailable placements, including synthetic
// reads passed as nullptr, preserve existing ALT eligibility. On false, callers
// discard ALT support without changing REF/ER classification or depth accounting.
inline bool passes_original_bam_alt_veto(const bam1_t* original_read, const char* contig_seq, hts_pos_t contig_len, int alt_score) {
    int64_t original_score;
    return !score_original_bam_placement(original_read, contig_seq, contig_len, original_score) || alt_score > original_score;
}

#endif
