#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "htslib/sam.h"
#include "htslib/tbx.h"
#include "htslib/vcf.h"

#include "../libs/cptl_stl.h"
#include "../libs/ssw_cpp.h"
#include "utils.h"
#include "sam_utils.h"
#include "hp_read_info.h"
#include "vcf_utils.h"

const int MIN_REF_HP_LEN = 5;
const int MIN_SUPPORTING_READS = 3;
const hts_pos_t CHUNK_SIZE = 1000000;

using mate_map_t = std::unordered_map<std::string, std::pair<std::string, int>>;
using max_3p_error_rates_t = std::map<int, double>;

struct hp_run_t {
    hts_pos_t beg, end;
    char base;
    int usable_reads = 0;
    std::unordered_map<int, int> hp_len_counts;

    hp_run_t(hts_pos_t beg, hts_pos_t end, char base) : beg(beg), end(end), base(base) {}
};

max_3p_error_rates_t read_max_3p_error_rates(std::string fname) {
    std::ifstream fin(fname);
    if (!fin) throw std::runtime_error("Unable to open " + fname + ".");

    std::string hp_len_header, max_error_rate_header;
    fin >> hp_len_header >> max_error_rate_header;
    if (hp_len_header != "HP_LEN" || max_error_rate_header.empty()) {
        throw std::runtime_error("Invalid header in " + fname + ".");
    }

    max_3p_error_rates_t max_3p_error_rates;
    int hp_len;
    double max_error_rate;
    while (fin >> hp_len >> max_error_rate) {
        if (hp_len < 0 || max_error_rate < 0 || max_error_rate > 1 || !max_3p_error_rates.emplace(hp_len, max_error_rate).second) {
            throw std::runtime_error("Invalid entry in " + fname + ".");
        }
    }
    if (!fin.eof() || max_3p_error_rates.empty()) throw std::runtime_error("Invalid data in " + fname + ".");
    return max_3p_error_rates;
}

double get_max_3p_error_rate(const max_3p_error_rates_t& max_3p_error_rates, int hp_len) {
    auto it = max_3p_error_rates.upper_bound(hp_len);
    if (it != max_3p_error_rates.begin()) it--;
    return it->second;
}

bool is_usable_hp_read(const hp_read_info_t& hp_read_info, int min_clip_len, double max_seq_error,
    const max_3p_error_rates_t& max_3p_error_rates) {

    if (hp_read_info.hp_deletion_extends_outside_hp || hp_read_info.hp_insertion_has_non_hp_bases) return false;
    if (hp_read_info.tail_5p_len < min_clip_len || hp_read_info.tail_3p_len < min_clip_len) return false;
    if (double(hp_read_info.tail_5p_mismatches) / hp_read_info.tail_5p_len > max_seq_error) return false;
    return double(hp_read_info.tail_3p_mismatches) / hp_read_info.tail_3p_len <= get_max_3p_error_rate(max_3p_error_rates, hp_read_info.hp_len);
}

std::vector<hp_run_t> find_hp_runs(char* contig_seq, hts_pos_t contig_len, hts_pos_t chunk_beg, hts_pos_t chunk_end) {

    std::vector<hp_run_t> hp_runs;
    hts_pos_t beg = chunk_beg;
    if (beg > 0 && beg < contig_len && contig_seq[beg] == contig_seq[beg - 1]) {
        char base = contig_seq[beg];
        while (beg < contig_len && contig_seq[beg] == base) beg++;
    }

    for (; beg < chunk_end;) {
        char base = contig_seq[beg];
        hts_pos_t end = beg + 1;
        while (end < contig_len && contig_seq[end] == base) end++;

        if ((base == 'A' || base == 'C' || base == 'G' || base == 'T') && end - beg >= MIN_REF_HP_LEN) {
            hp_runs.emplace_back(beg, end, base);
        }
        beg = end;
    }
    return hp_runs;
}

std::string workdir;
std::vector<mate_map_t> mateseqs_w_mapq;
std::vector<int> active_threads_per_chr;
std::vector<std::mutex> mutex_per_chr;

void read_mates(int contig_id) {
    std::lock_guard<std::mutex> lock(mutex_per_chr[contig_id]);
    if (active_threads_per_chr[contig_id] == 0) {
        std::string fname = workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt";
        std::ifstream fin(fname);
        std::string qname, read_seq, qual;
        int mapq;
        while (fin >> qname >> read_seq >> qual >> mapq) {
            mateseqs_w_mapq[contig_id][qname] = {read_seq, mapq};
        }
    }
    active_threads_per_chr[contig_id]++;
}

void release_mates(int contig_id) {
    std::lock_guard<std::mutex> lock(mutex_per_chr[contig_id]);
    active_threads_per_chr[contig_id]--;
    if (active_threads_per_chr[contig_id] == 0) {
        mateseqs_w_mapq[contig_id].clear();
    }
}

void add_rescued_hp_read(hp_run_t& hp_run,
    const std::pair<std::string, int>& mate,
    bool anchor_is_reverse, char* contig_seq, hts_pos_t contig_len,
    int read_len, int min_clip_len, double max_seq_error, StripedSmithWaterman::Aligner& aligner,
    const StripedSmithWaterman::Filter& filter, const max_3p_error_rates_t& max_3p_error_rates) {

    hts_pos_t extend = std::max(0, read_len - 1);
    hts_pos_t allele_beg = std::max((hts_pos_t) 0, hp_run.beg - extend);
    hts_pos_t allele_end = std::min(contig_len, hp_run.end + extend);
    std::string ref_allele(contig_seq + allele_beg, allele_end - allele_beg);
    hts_pair_pos_t allele_hp_range = {hp_run.beg - allele_beg, hp_run.end - allele_beg};

    std::string mate_seq = mate.first;
    if (!anchor_is_reverse) rc(mate_seq);
    StripedSmithWaterman::Alignment aln;
    if (!aligner.Align(mate_seq.c_str(), ref_allele.c_str(), ref_allele.length(), filter, &aln, 0)) {
        return;
    }
    bool aln_as_rev = !anchor_is_reverse;
    int left_clip = aln.query_begin;
    int right_clip = mate_seq.length() - aln.query_end - 1;
    if (aln.ref_begin < 0 || aln.cigar.empty() || (!aln_as_rev && left_clip > 0) || (aln_as_rev && right_clip > 0)) {
        return;
    }

    hp_read_info_t hp_read_info = calculate_hp_read_info(aln, mate_seq,
        allele_hp_range, hp_run.base, &ref_allele[0], ref_allele.length(), aln_as_rev, bp_support_read_t(), 0);
    if (!is_usable_hp_read(hp_read_info, min_clip_len, max_seq_error, max_3p_error_rates)) return;
    hp_run.usable_reads++;
    hp_run.hp_len_counts[hp_read_info.hp_len]++;
}


std::vector<std::shared_ptr<sv_t>> find_hp_indels_for_chunk(int id, size_t contig_id, std::string contig_name,
    char* contig_seq, hts_pos_t contig_len, hts_pos_t chunk_beg, hts_pos_t chunk_end,
    const config_t* config, const stats_t* stats, bam_pool_t* bam_pool,
    StripedSmithWaterman::Aligner& aligner, const StripedSmithWaterman::Filter& filter,
    const max_3p_error_rates_t& max_3p_error_rates) {

    std::vector<std::shared_ptr<sv_t>> hp_indels;
    if (contig_len == 0) return hp_indels;

    std::vector<hp_run_t> hp_runs = find_hp_runs(contig_seq, contig_len, chunk_beg, chunk_end);
    if (hp_runs.empty()) return hp_indels;

    std::vector<hts_pos_t> hp_run_ends;
    std::vector<hts_pos_t> hp_run_begs;
    hp_run_ends.reserve(hp_runs.size());
    hp_run_begs.reserve(hp_runs.size());
    for (const hp_run_t& hp_run : hp_runs) {
        hp_run_begs.push_back(hp_run.beg);
        hp_run_ends.push_back(hp_run.end);
    }

    open_samFile_t* alignment_file = bam_pool->get_bam_reader(id);
    read_mates(contig_id);
    try {
        const mate_map_t& mateseqs_w_mapq_chr = mateseqs_w_mapq[contig_id];
        int min_anchor_mapq = config->high_confidence_mapq;
        hts_pos_t query_beg = std::max((hts_pos_t) 0, hp_runs.front().beg - stats->max_is);
        hts_pos_t query_end = std::min(contig_len, hp_runs.back().end + stats->max_is);
        std::stringstream region_ss;
        region_ss << contig_name << ":" << query_beg + 1 << "-" << query_end;
        std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iter(sam_itr_querys(alignment_file->idx, alignment_file->header, region_ss.str().c_str()), &hts_itr_destroy);
        if (!iter) throw std::runtime_error("Unable to query alignments for " + contig_name + ".");

        std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read(bam_init1(), &bam_destroy1);
        if (!read) throw std::runtime_error("Unable to allocate an alignment record.");

        int read_status;
        while ((read_status = sam_itr_next(alignment_file->file, iter.get(), read.get())) >= 0) {
            if (is_unmapped(read.get()) || !is_primary(read.get()) || read->core.l_qseq <= 0) continue;

            hts_pos_t read_beg = read->core.pos;
            hts_pos_t read_end = bam_endpos(read.get());
            size_t hp_idx = std::upper_bound(hp_run_ends.begin(), hp_run_ends.end(), read_beg) - hp_run_ends.begin();
            for (; hp_idx < hp_runs.size() && hp_runs[hp_idx].beg < read_end; hp_idx++) {
                hp_run_t& hp_run = hp_runs[hp_idx];

                hts_pair_pos_t hp_range = {hp_run.beg, hp_run.end};
                hp_read_info_t hp_read_info = calculate_hp_read_info(read.get(), hp_range, hp_run.base, contig_seq, contig_len, 0);
                if (!is_usable_hp_read(hp_read_info, config->min_clip_len, config->max_seq_error, max_3p_error_rates)) continue;
                hp_run.usable_reads++;
                hp_run.hp_len_counts[hp_read_info.hp_len]++;
            }

            if (read->core.qual < min_anchor_mapq || !is_dc_pair(read.get()) || mateseqs_w_mapq_chr.empty()) continue;
            std::string mate_qname = get_mate_lookup_qname(read.get());
            auto mate_it = mateseqs_w_mapq_chr.find(mate_qname);
            if (mate_it == mateseqs_w_mapq_chr.end()) continue;

            if (!bam_is_rev(read.get())) {
                hts_pos_t min_hp_beg = read->core.pos + stats->read_len / 2;
                hts_pos_t max_hp_beg = read->core.pos + stats->max_is;
                size_t rescue_hp_idx = std::lower_bound(hp_run_begs.begin(), hp_run_begs.end(), min_hp_beg) - hp_run_begs.begin();
                for (; rescue_hp_idx < hp_runs.size() && hp_runs[rescue_hp_idx].beg <= max_hp_beg; rescue_hp_idx++) {
                    add_rescued_hp_read(hp_runs[rescue_hp_idx], mate_it->second, false,
                        contig_seq, contig_len, stats->read_len, config->min_clip_len, config->max_seq_error,
                        aligner, filter, max_3p_error_rates);
                }
            } else {
                hts_pos_t min_hp_end = read_end - stats->max_is;
                hts_pos_t max_hp_end = read_end - stats->read_len / 2;
                size_t rescue_hp_idx = std::lower_bound(hp_run_ends.begin(), hp_run_ends.end(), min_hp_end) - hp_run_ends.begin();
                for (; rescue_hp_idx < hp_runs.size() && hp_runs[rescue_hp_idx].end <= max_hp_end; rescue_hp_idx++) {
                    add_rescued_hp_read(hp_runs[rescue_hp_idx], mate_it->second, true,
                        contig_seq, contig_len, stats->read_len, config->min_clip_len, config->max_seq_error,
                        aligner, filter, max_3p_error_rates);
                }
            }
        }
        if (read_status < -1) throw std::runtime_error("Error while reading alignments for " + contig_name + ".");

        for (const hp_run_t& hp_run : hp_runs) {
            int ref_hp_len = hp_run.end - hp_run.beg;
            for (const auto& hp_len_count : hp_run.hp_len_counts) {
                int alt_hp_len = hp_len_count.first;
                int supporting_reads = hp_len_count.second;
                if (hp_run.beg == 0 || alt_hp_len == ref_hp_len || supporting_reads < MIN_SUPPORTING_READS) continue;

                int hp_len_diff = alt_hp_len - ref_hp_len;
                if (std::abs(hp_len_diff) < config->min_sv_size) continue;
                hts_pos_t anchor_pos = hp_run.beg - 1;
                std::shared_ptr<sv_t> hp_indel;
                if (hp_len_diff > 0) {
                    hp_indel = std::make_shared<insertion_t>(contig_name, anchor_pos, anchor_pos,
                        std::string(hp_len_diff, hp_run.base), nullptr, nullptr, nullptr, nullptr);
                } else {
                    hp_indel = std::make_shared<deletion_t>(contig_name, anchor_pos,
                        anchor_pos - hp_len_diff, "", nullptr, nullptr, nullptr, nullptr);
                }
                hp_indel->source = "HP";
                hp_indel->junction_remap_ref_beg = hp_indel->start;
                hp_indel->junction_remap_ref_end = hp_indel->end + 1;
                hp_indel->hp_ref_beg = hp_run.beg;
                hp_indel->hp_ref_end = hp_run.end;
                hp_indels.push_back(hp_indel);
            }
        }
    } catch (...) {
        release_mates(contig_id);
        throw;
    }
    release_mates(contig_id);
    return hp_indels;
}

int main(int argc, char* argv[]) {
    if (argc != 5 && argc != 6) {
        std::cerr << "Usage: find_hp_indels <workdir> <reference.fa> <alignments.bam|cram> <output.vcf.gz> [max-3p-error-rates.tsv]\n";
        return 1;
    }

    workdir = argv[1];
    std::string reference_fname = argv[2];
    std::string alignment_fname = argv[3];
    std::string out_vcf_fname = argv[4];
    std::string max_3p_error_rates_fname = argc == 6 ? argv[5] : "resources/max-3p-error-rate-by-hp-len.tsv";

    config_t config;
    config.parse(workdir + "/config.txt");

    stats_t stats;
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);

    max_3p_error_rates_t max_3p_error_rates = read_max_3p_error_rates(max_3p_error_rates_fname);

    contig_map_t contig_map(workdir);

    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(reference_fname);

    bam_pool_t bam_pool(config.threads, alignment_fname, reference_fname);
    mateseqs_w_mapq.resize(contig_map.size());
    active_threads_per_chr = std::vector<int>(contig_map.size());
    mutex_per_chr = std::vector<std::mutex>(contig_map.size());
    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    StripedSmithWaterman::Filter filter;

    std::vector<std::future<std::vector<std::shared_ptr<sv_t>>>> futures;
    ctpl::thread_pool thread_pool(config.threads);
    for (size_t contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        hts_pos_t contig_len = chr_seqs.get_len(contig_name);
        for (hts_pos_t chunk_beg = 0; chunk_beg < contig_len; chunk_beg += CHUNK_SIZE) {
            futures.push_back(thread_pool.push(find_hp_indels_for_chunk, contig_id, contig_name, chr_seqs.get_seq(contig_name), contig_len,
                chunk_beg, std::min(contig_len, chunk_beg + CHUNK_SIZE), &config, &stats, &bam_pool, std::ref(aligner), std::cref(filter),
                std::cref(max_3p_error_rates)));
        }
    }
    thread_pool.stop(true);

    std::vector<std::shared_ptr<sv_t>> hp_indels;
    for (std::future<std::vector<std::shared_ptr<sv_t>>>& future : futures) {
        std::vector<std::shared_ptr<sv_t>> chunk_hp_indels = future.get();
        std::sort(chunk_hp_indels.begin(), chunk_hp_indels.end(), sv_output_order);
        hp_indels.insert(hp_indels.end(), chunk_hp_indels.begin(), chunk_hp_indels.end());
    }

    std::string command;
    for (int i = 0; i < argc; i++) {
        if (!command.empty()) command += " ";
        command += argv[i];
    }
    std::unique_ptr<bcf_hdr_t, decltype(&bcf_hdr_destroy)> hdr( generate_vcf_header(chr_seqs, "", config, command), &bcf_hdr_destroy);
    std::unique_ptr<htsFile, decltype(&hts_close)> out(bcf_open(out_vcf_fname.c_str(), "wz"), &hts_close);
    if (!out) throw std::runtime_error("Unable to open " + out_vcf_fname + " for writing.");
    if (hts_set_threads(out.get(), std::max(1, std::min(config.threads, 4))) != 0) {
        throw std::runtime_error("Unable to enable threaded VCF compression.");
    }
    if (bcf_hdr_write(out.get(), hdr.get()) != 0) {
        throw std::runtime_error("Unable to write the VCF header to " + out_vcf_fname + ".");
    }
    std::unique_ptr<bcf1_t, decltype(&bcf_destroy)> record(bcf_init(), &bcf_destroy);
    if (!record) throw std::runtime_error("Unable to allocate a VCF record.");
    for (size_t i = 0; i < hp_indels.size(); i++) {
        std::shared_ptr<sv_t>& hp_indel = hp_indels[i];
        hp_indel->id = "HP_" + hp_indel->svtype() + "_" + std::to_string(i + 1);
        hp_indel->sample_info.gt.clear();
        sv2bcf(hdr.get(), record.get(), hp_indel.get(), chr_seqs.get_seq(hp_indel->chr));
        int hp_ref_range[] = {(int) hp_indel->hp_ref_beg + 1, (int) hp_indel->hp_ref_end + 1};
        bcf_update_info_int32(hdr.get(), record.get(), "HP_REF_RANGE", hp_ref_range, 2);
        if (bcf_write(out.get(), hdr.get(), record.get()) != 0) {
            throw std::runtime_error("Unable to write VCF record " + hp_indel->id + ".");
        }
    }
    out.reset();

    if (tbx_index_build(out_vcf_fname.c_str(), 0, &tbx_conf_vcf) != 0) {
        throw std::runtime_error("Unable to index " + out_vcf_fname + ".");
    }
}
