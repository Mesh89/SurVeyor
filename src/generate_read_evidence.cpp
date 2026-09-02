#include <fstream>
#include <future>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "htslib/vcf.h"
#include "../libs/cptl_stl.h"

#include "genotype.h"
#include "genotype_dels.h"
#include "genotype_dups.h"
#include "genotype_hp_indels.h"
#include "genotype_inss.h"

config_t config;
const bool USE_HP_SPECIFIC_PATH = false;

std::vector<ext_mate_map_t> mates_by_contig;
std::vector<std::mutex> mates_mutex_by_contig;
std::vector<int> pending_mates_users_by_contig;
std::vector<bool> mates_loaded_by_contig;
std::vector<bool> mates_scheduling_complete_by_contig;

ext_mate_map_t load_mates(const std::string& workdir, int contig_id) {
    ext_mate_map_t mates;
    std::ifstream fin(workdir + "/workspace/mateseqs/" + std::to_string(contig_id) + ".txt");
    std::string qname, read_seq, qual;
    int mapq;
    while (fin >> qname >> read_seq >> qual >> mapq) mates[qname] = {read_seq, qual, mapq};
    return mates;
}

void register_mates_user(int contig_id) {
    std::lock_guard<std::mutex> lock(mates_mutex_by_contig[contig_id]);
    pending_mates_users_by_contig[contig_id]++;
}

ext_mate_map_t* acquire_mates(const std::string& workdir, int contig_id) {
    std::lock_guard<std::mutex> lock(mates_mutex_by_contig[contig_id]);
    if (!mates_loaded_by_contig[contig_id]) {
        mates_by_contig[contig_id] = load_mates(workdir, contig_id);
        mates_loaded_by_contig[contig_id] = true;
    }
    return &mates_by_contig[contig_id];
}

void release_mates(int contig_id) {
    std::lock_guard<std::mutex> lock(mates_mutex_by_contig[contig_id]);
    pending_mates_users_by_contig[contig_id]--;
    if (pending_mates_users_by_contig[contig_id] == 0 && mates_scheduling_complete_by_contig[contig_id]) {
        mates_by_contig[contig_id].clear();
        mates_loaded_by_contig[contig_id] = false;
    }
}

void complete_mates_scheduling(int contig_id) {
    std::lock_guard<std::mutex> lock(mates_mutex_by_contig[contig_id]);
    mates_scheduling_complete_by_contig[contig_id] = true;
    if (pending_mates_users_by_contig[contig_id] == 0) {
        mates_by_contig[contig_id].clear();
        mates_loaded_by_contig[contig_id] = false;
    }
}

void write_too_deep_variants(const std::string& association_dir, const std::string& chr, const std::vector<sv_t*>& svs) {
    std::ofstream fout;
    for (sv_t* sv : svs) {
        if (!sv->sample_info.too_deep) continue;
        if (!fout.is_open()) fout.open(association_dir + "/" + chr + ".td.txt");
        fout << sv->id << '\n';
    }
}

void generate_read_evidence_block(int id, int contig_id, const std::string& chr, char* contig_seq, hts_pos_t contig_len, std::vector<sv_t*> svs, stats_t& stats, config_t& config, bam_pool_t* bam_pool, const std::string& workdir, evidence_logger_t* evidence_logger) {
    open_samFile_t* bam_file = bam_pool->get_bam_reader(id);
    StripedSmithWaterman::Aligner aligner(1, 4, 6, 1, false);
    std::unordered_map<std::string, std::vector<sv_t*>> hp_indels_by_range;

    for (sv_t* sv : svs) {
        if (USE_HP_SPECIFIC_PATH && should_genotype_as_hp_indel(sv, contig_seq, contig_len)) {
            hts_pair_pos_t ref_hp_range = find_ref_hp_range_for_indel(sv, contig_seq, contig_len);
            hp_indels_by_range[std::to_string(ref_hp_range.beg) + ":" + std::to_string(ref_hp_range.end)].push_back(sv);
        } else if (sv->svtype() == "DEL") {
            deletion_t* del = static_cast<deletion_t*>(sv);
            write_aligned_del_read_evidence(del, bam_file, contig_seq, contig_len, stats, config, aligner, *evidence_logger);
        } else if (sv->svtype() == "DUP") {
            duplication_t* dup = static_cast<duplication_t*>(sv);
            if (dup->svlen() <= stats.read_len-2*config.min_clip_len) write_aligned_small_dup_read_evidence(dup, bam_file, contig_seq, contig_len, stats, config, aligner, *evidence_logger);
            else write_aligned_large_dup_read_evidence(dup, bam_file, contig_seq, contig_len, stats, config, aligner, *evidence_logger);
        } else if (sv->svtype() == "INS") {
            insertion_t* ins = static_cast<insertion_t*>(sv);
            write_aligned_ins_read_evidence(ins, bam_file, contig_seq, contig_len, stats, config, aligner, *evidence_logger);
        }
    }

    if (hp_indels_by_range.empty()) return;
    ext_mate_map_t* mates = acquire_mates(workdir, contig_id);
    for (auto& hp_indels : hp_indels_by_range) {
        hts_pair_pos_t ref_hp_range = find_ref_hp_range_for_indel(hp_indels.second[0], contig_seq, contig_len);
        write_aligned_hp_indels_group_read_evidence(hp_indels.second, ref_hp_range, bam_file, contig_seq, contig_len, stats, config, aligner, *mates, *evidence_logger);
    }
    release_mates(contig_id);
}

int main(int argc, char* argv[]) {
    if (argc != 5) throw std::runtime_error("Usage: generate_read_evidence <input.vcf.gz> <alignment.bam> <reference.fa> <workdir>");

    std::string in_vcf_fname = argv[1];
    std::string bam_fname = argv[2];
    std::string reference_fname = argv[3];
    std::string workdir = argv[4];
    std::string association_dir = workdir + "/reads_to_sv_associations";

    contig_map_t contig_map;
    contig_map.load(workdir);
    mates_by_contig.resize(contig_map.size());
    mates_mutex_by_contig = std::vector<std::mutex>(contig_map.size());
    pending_mates_users_by_contig = std::vector<int>(contig_map.size());
    mates_loaded_by_contig = std::vector<bool>(contig_map.size());
    mates_scheduling_complete_by_contig = std::vector<bool>(contig_map.size());
    config.parse(workdir + "/config.txt");
    stats_t stats;
    stats.parse(workdir + "/stats.txt", config.per_contig_stats);
    chr_seqs_map_t chr_seqs;
    chr_seqs.read_fasta_into_map(reference_fname);
    bam_pool_t bam_pool(config.threads, bam_fname, reference_fname);
    evidence_logger_t evidence_logger(association_dir);

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    if (!in_vcf_file) throw std::runtime_error("Unable to open file " + in_vcf_fname + ".");
    bcf_hdr_t* in_vcf_header = bcf_hdr_read(in_vcf_file);
    if (!in_vcf_header) throw std::runtime_error("Failed to read the VCF header.");

    std::vector<std::vector<std::shared_ptr<sv_t>>> svs_by_contig(contig_map.size());
    bcf1_t* vcf_record = bcf_init();
    while (bcf_read(in_vcf_file, in_vcf_header, vcf_record) == 0) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(in_vcf_header, vcf_record);
        if (!sv || sv->start > sv->end) continue;
        svs_by_contig[contig_map.get_id(sv->chr)].push_back(sv);
    }

    ctpl::thread_pool thread_pool(config.threads);
    std::vector<std::future<void>> futures;
    const int BLOCK_SIZE = 20;
    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        if (svs_by_contig[contig_id].empty()) continue;
        std::string chr = contig_map.get_name(contig_id);
        char* contig_seq = chr_seqs.get_seq(chr);
        hts_pos_t contig_len = chr_seqs.get_len(chr);
        std::vector<sv_t*> hp_indels, dels, dups, inss;
        for (const std::shared_ptr<sv_t>& sv : svs_by_contig[contig_id]) {
            if (USE_HP_SPECIFIC_PATH && should_genotype_as_hp_indel(sv.get(), contig_seq, contig_len)) hp_indels.push_back(sv.get());
            else if (sv->svtype() == "DEL") dels.push_back(sv.get());
            else if (sv->svtype() == "DUP") dups.push_back(sv.get());
            else if (sv->svtype() == "INS") inss.push_back(sv.get());
        }

        auto schedule_block = [&](std::vector<sv_t*>& block) {
            std::future<void> future = thread_pool.push(generate_read_evidence_block, contig_id, chr, contig_seq, contig_len, block, std::ref(stats), std::ref(config), &bam_pool, workdir, &evidence_logger);
            futures.push_back(std::move(future));
            block.clear();
        };

        std::vector<hts_pair_pos_t> ref_hp_ranges;
        for (sv_t* hp_indel : hp_indels) ref_hp_ranges.push_back(find_ref_hp_range_for_indel(hp_indel, contig_seq, contig_len));
        std::vector<sv_t*> block;
        for (int i = 0; i < hp_indels.size(); i++) {
            block.push_back(hp_indels[i]);
            bool last = i == hp_indels.size()-1;
            bool next_range_differs = !last && (ref_hp_ranges[i].beg != ref_hp_ranges[i+1].beg || ref_hp_ranges[i].end != ref_hp_ranges[i+1].end);
            if (last || (block.size() >= BLOCK_SIZE && next_range_differs)) {
                register_mates_user(contig_id);
                schedule_block(block);
            }
        }
        complete_mates_scheduling(contig_id);
        for (std::vector<sv_t*>* variants : {&dels, &dups, &inss}) {
            for (sv_t* sv : *variants) {
                block.push_back(sv);
                if (block.size() == BLOCK_SIZE) schedule_block(block);
            }
            if (!block.empty()) schedule_block(block);
        }
    }
    thread_pool.stop(true);
    for (std::future<void>& future : futures) future.get();

    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        if (svs_by_contig[contig_id].empty()) continue;
        std::vector<sv_t*> svs;
        for (const std::shared_ptr<sv_t>& sv : svs_by_contig[contig_id]) svs.push_back(sv.get());
        write_too_deep_variants(association_dir, contig_map.get_name(contig_id), svs);
    }

    bcf_destroy(vcf_record);
    bcf_hdr_destroy(in_vcf_header);
    bcf_close(in_vcf_file);
    return 0;
}
