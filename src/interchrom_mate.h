#ifndef INTERCHROM_MATE_H
#define INTERCHROM_MATE_H

#include <map>
#include <mutex>
#include "htslib/bgzf.h"
#include <tuple>
#include "genotype.h"
#include "quality_read_scoring.h"

inline bool mapped_interchrom_pair(const bam1_t* b) {
    return (b->core.flag & BAM_FPAIRED) && !(b->core.flag & (BAM_FUNMAP | BAM_FMUNMAP))
        && b->core.tid >= 0 && b->core.mtid >= 0 && b->core.tid != b->core.mtid;
}

struct mate_alt_context_t {
    std::string sequence;
    hts_pos_t left = 0;
};

struct interchrom_mate_check_t {
    bam_pool_t readers;
    const chr_seqs_map_t& reference;
    struct scored_mate_t {
        bool found = false, reverse = false;
        std::string seq[2];
        quality_read_encoding_t encoded[2];
        int64_t bam_score = -1;
    };
    using mate_key_t = std::tuple<int, hts_pos_t, bool, std::string, bool, std::string, hts_pos_t>;
    std::vector<std::map<mate_key_t, scored_mate_t>> cache;
    interchrom_mate_check_t(int threads, const std::string& bam, const std::string& ref,
                            const chr_seqs_map_t& reference)
        : readers(threads, bam, ref), reference(reference), cache(threads) {}

    // File offsets are shared only between workers reading this same sorted BAM.
    std::mutex offsets_mutex;
    std::map<std::pair<int, hts_pos_t>, int64_t> offsets;

    int64_t preceding_offset(int tid, hts_pos_t pos, hts_pos_t begin) {
        std::lock_guard<std::mutex> lock(offsets_mutex);
        auto it = offsets.lower_bound({tid, pos});
        if (it == offsets.begin()) return -1;
        --it;
        if (it->first.first != tid || it->first.second < begin) return -1;
        return it->second;
    }

    void remember_offset(const bam1_t* read, int64_t after_record) {
        std::lock_guard<std::mutex> lock(offsets_mutex);
        auto& offset = offsets[{read->core.tid, read->core.pos}];
        offset = std::max(offset, after_record);
    }

    std::shared_ptr<bam1_t> mate(int worker, bam1_t* stable) {
        open_samFile_t* f = readers.get_bam_reader(worker);
        const int tid = stable->core.mtid;
        const hts_pos_t pos = stable->core.mpos;
        BGZF* bgzf = f->file->format.format == bam && f->file->is_bgzf ? f->file->fp.bgzf : nullptr;
        // On the first lookup in a 16 kb interval, retain positions along the way.
        const hts_pos_t begin = bgzf ? pos - pos % 16384 : pos;
        const int64_t offset = bgzf ? preceding_offset(tid, pos, begin) : -1;
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> it(nullptr, hts_itr_destroy);
        if (offset >= 0) {
            // Resume after a strictly earlier start. All records tied at pos are still ahead.
            if (bgzf_seek(bgzf, offset, SEEK_SET) < 0) throw std::runtime_error("Failed to seek to mapped mate");
        } else {
            it.reset(sam_itr_queryi(f->idx, tid, begin, pos + 1));
            if (!it) throw std::runtime_error("Invalid mapped mate interval");
        }
        int64_t last_bucket = -1;
        std::shared_ptr<bam1_t> m(bam_init1(), bam_destroy1);
        const int expected = (stable->core.flag & BAM_FREAD1) ? BAM_FREAD2 : BAM_FREAD1;
        int status;
        while ((status = it ? sam_itr_next(f->file, it.get(), m.get()) : sam_read1(f->file, f->header, m.get())) >= 0) {
            if (m->core.tid != tid || m->core.pos > pos) break;
            if (bgzf) {
                const int64_t after = bgzf_tell(bgzf);
                const int64_t bucket = after >> 32; // The upper bits select 64 KiB of compressed file.
                if (bucket != last_bucket) remember_offset(m.get(), after);
                last_bucket = bucket;
            }
            if (!is_primary(m.get()) || is_unmapped(m.get()) || m->core.pos != stable->core.mpos ||
                (m->core.flag & (BAM_FREAD1 | BAM_FREAD2)) != expected || strcmp(bam_get_qname(m.get()), bam_get_qname(stable))) continue;
            auto a = bam_aux_get(stable, "RG"), b = bam_aux_get(m.get(), "RG");
            if ((!a && !b) || (a && b && !strcmp(bam_aux2Z(a), bam_aux2Z(b)))) return m;
        }
        if (status < -1) throw std::runtime_error("Error reading mapped mate");
        return {};
    }

    const scored_mate_t& score_mate(int worker, bam1_t* stable, hts_pos_t pad, config_t& config) {
        auto group = bam_aux_get(stable, "RG");
        mate_key_t key(stable->core.mtid, stable->core.mpos, bool(stable->core.flag & BAM_FREAD1),
                  bam_get_qname(stable), group != nullptr, group ? bam_aux2Z(group) : "", pad);
        auto& entries = cache[worker];
        auto hit = entries.find(key);
        if (hit != entries.end()) return hit->second;
        // Each worker revisits the same mates while it processes nearby candidates.
        if (entries.size() >= 512) entries.clear();
        auto& result = entries[std::move(key)];
        auto m = mate(worker, stable);
        if (!m || !m->core.l_qseq) return result;
        result.found = true;
        result.reverse = bam_is_rev(m.get());
        auto reader = readers.get_bam_reader(worker);
        auto ctg = reference.seqs.at(reader->header->target_name[m->core.tid]);
        hts_pos_t from = std::max<hts_pos_t>(0, m->core.pos - pad);
        hts_pos_t to = std::min<hts_pos_t>(ctg->len, bam_endpos(m.get()) + pad);
        const std::string window(ctg->seq + from, to - from);
        const StripedSmithWaterman::Filter scores_only(false, false, 0, 32767);
        result.seq[0] = result.seq[1] = get_sequence(m.get());
        rc(result.seq[1]);
        for (int rev = 0; rev < 2; ++rev) {
            auto& enc = result.encoded[rev];
            enc = encode_quality_read(m.get(), result.seq[rev], rev, ctg->seq, ctg->len, config);
            if (!rev) quality_score_original_bam_placement(m.get(), ctg->seq, ctg->len, enc, result.bam_score);
            // Every substitution or gap scores at most as well as a full-length perfect match.
            if (result.bam_score >= enc.match_sum) continue;
            auto fit = quality_align_fast(result.seq[rev].c_str(), enc, window.c_str(), window.size(), scores_only);
            result.bam_score = std::max<int64_t>(result.bam_score, fit.sw_score);
        }
        return result;
    }

    // Each anchor is the stable end's projected outer 5' base, relative to the left ALT junction.
    bool supports(int worker, sv_t* sv, bam1_t* stable, std::initializer_list<hts_pos_t> anchors,
                  int max_is, int read_len, char* chr_seq, hts_pos_t chr_len, config_t& config,
                  mate_alt_context_t& context) {
        const bool insertion = sv->svtype() == "INS";
        if (sv->incomplete_ins_seq()) return insertion;
        const hts_pos_t pad = max_is + read_len;
        const auto& m = score_mate(worker, stable, pad, config);
        if (!m.found) return insertion;
        if (context.sequence.empty()) {
            std::unique_ptr<char[]> left(generate_haplotype_left(
                chr_seq, sv->start, pad, sv->aux_indels, sv->aux_snps, nullptr, nullptr, sv->id, sv->chr));
            std::unique_ptr<char[]> right(generate_haplotype_right(
                chr_seq, chr_len, sv->end + 1, pad, sv->aux_indels, sv->aux_snps, nullptr, nullptr, sv->id, sv->chr));
            context.sequence = std::string(left.get()) + sv->ins_seq + right.get();
            context.left = strlen(left.get());
            to_uppercase(&context.sequence[0]);
        }
        // Only the inward-facing orientation can complete this fragment.
        const int rev = m.reverse == bool(bam_is_rev(stable));
        const StripedSmithWaterman::Filter scores_only(false, false, 0, 32767);
        int alt_best = -1;
        for (hts_pos_t offset : anchors) {
            if (offset == HTS_POS_MAX) continue;
            hts_pos_t outer = context.left + offset;
            hts_pos_t begin = std::max<hts_pos_t>(0, bam_is_rev(stable) ? outer - max_is + 1 : outer);
            hts_pos_t end = std::min<hts_pos_t>(context.sequence.size(), bam_is_rev(stable) ? outer + 1 : outer + max_is);
            if (begin >= end) continue;
            const std::string window = context.sequence.substr(begin, end - begin);
            auto hit = quality_align_fast(m.seq[rev].c_str(), m.encoded[rev], window.c_str(), window.size(), scores_only);
            alt_best = std::max(alt_best, int(hit.sw_score));
        }
        return alt_best > 0 && alt_best >= m.bam_score;
    }
};

inline hts_pos_t mate_anchor(bam1_t* stable, const StripedSmithWaterman::Alignment& aln, hts_pos_t junction_relative_offset) {
    hts_pos_t outer = bam_is_rev(stable) ? aln.ref_end + stable->core.l_qseq - 1 - aln.query_end
                                         : aln.ref_begin - aln.query_begin;
    return junction_relative_offset + outer;
}
#endif
