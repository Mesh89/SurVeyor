#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <memory>
#include <cerrno>
#include <sys/stat.h>

#include "htslib/faidx.h"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "sam_utils.h"
#include "utils.h"
#include "coarse_coverage.h"
#include "../libs/cptl_stl.h"

const int MIN_RND_POS = 1000;

std::string workdir;
std::string workspace;
std::string reference_fname;

config_t config;
stats_t stats;
chr_seqs_map_t chr_seqs;
htsThreadPool hts_pool = {NULL, 0};

std::mutex mtx;
std::mutex* mtx_contig;
std::vector<std::vector<std::string> > mate_seqs;

std::vector<hts_pos_t> general_isize_dist;
std::unordered_map<std::string, int> min_depth_by_contig, max_depth_by_contig, median_depth_by_contig;
std::vector<uint32_t> depths;
uint64_t qual_counts[256];

// dist_between_end_and_rnd[i]: Distribution of read pairs count with one read < x and the other between x and x+i, for many random x.
std::vector<std::vector<uint32_t> > dist_between_end_and_rnd; 
std::vector<std::pair<uint64_t, uint32_t> > partial_sums;

std::vector<uint32_t> isize_counts;
std::vector<std::vector<uint32_t> > isizes_count_geq_i;

std::ofstream open_mateseqs_fout(int contig_id) {
    std::string fname = std::to_string(contig_id) + ".txt";
    return std::ofstream(workspace + "/mateseqs/" + fname, std::ios_base::app);
}

void get_base_stats(int id, std::string bam_fname, std::string reference_fname, std::vector<std::pair<std::string, hts_pos_t> > rnd_positions) {
    open_samFile_t bam_file(bam_fname);
    if (hts_set_fai_filename(bam_file.file, fai_path(reference_fname.c_str())) != 0) {
        throw std::runtime_error("Failed to read reference " + reference_fname);
    }

    std::vector<hts_pos_t> local_dist;
    bam1_t* read = bam_init1();
    int read_len = 0;

    for (auto& rnd_pos : rnd_positions) {
        std::stringstream ss;
        ss << rnd_pos.first << ":" << rnd_pos.second << "-" << rnd_pos.second+1000;
        hts_itr_t* iter = sam_itr_querys(bam_file.idx, bam_file.header, ss.str().c_str());
        if (!iter) continue;

        int added = 0;
        while (sam_itr_next(bam_file.file, iter, read) >= 0) {
            if (!(read->core.flag & BAM_FPROPER_PAIR) || !is_primary(read)) continue;
            if (read->core.isize < 0 || read->core.isize > 20000) continue;

            added++;
            local_dist.push_back(read->core.isize);
            read_len = std::max(read_len, read->core.l_qseq);

            if (added > 100) break;
        }
        hts_itr_destroy(iter);

        if (local_dist.size() > 100000) break;
    }

    std::lock_guard<std::mutex> lock(mtx);
    general_isize_dist.insert(general_isize_dist.end(), local_dist.begin(), local_dist.end());
    stats.read_len = std::max(stats.read_len, read_len);

    bam_destroy1(read);
}

void categorize(int id, int contig_id, std::string contig_name, std::string bam_fname, std::string reference_fname, std::vector<hts_pos_t> rnd_positions) {
    std::sort(rnd_positions.begin(), rnd_positions.end());

    open_samFile_t bam_file(bam_fname);
    if (hts_set_fai_filename(bam_file.file, fai_path(reference_fname.c_str())) != 0) {
        throw std::runtime_error("Failed to read reference " + reference_fname);
    }
    if (hts_pool.pool && hts_set_opt(bam_file.file, HTS_OPT_THREAD_POOL, &hts_pool) != 0) {
        throw std::runtime_error("Failed to set thread pool for " + bam_fname);
    }

    std::unique_ptr<hts_itr_t, decltype(&hts_itr_destroy)> iter(sam_itr_querys(bam_file.idx, bam_file.header, contig_name.c_str()), &hts_itr_destroy);
    int bam_contig_id = bam_name2id(bam_file.header, contig_name.c_str());
    if (bam_contig_id < 0) throw std::runtime_error("Contig " + contig_name + " is missing from the alignment header");
    coarse_coverage_builder_t coverage(bam_file.header->target_len[bam_contig_id], config.coverage_window_size);
    std::string coverage_dir = workspace + "/coverage";
    if (::mkdir(coverage_dir.c_str(), 0755) != 0 && errno != EEXIST) throw std::runtime_error("Failed to create " + coverage_dir);
    std::string coverage_fname = coverage_dir + "/" + std::to_string(contig_id) + ".bin";
    if (!iter) { // no reads
        write_coarse_coverage(coverage_fname, coverage.maxima);
    	return;
    }

    // HTSlib retains these filename pointers until sam_idx_save finishes.
    const std::string sr_index_fname = workspace + "/sr/" + std::to_string(contig_id) + ".bam.bai";
    const std::string hsr_index_fname = workspace + "/hsr/" + std::to_string(contig_id) + ".bam.bai";
    std::unique_ptr<samFile, decltype(&hts_close)> sr_writer(nullptr, &hts_close);
    std::unique_ptr<samFile, decltype(&hts_close)> hsr_writer(nullptr, &hts_close);
    std::unique_ptr<samFile, decltype(&hts_close)> rdc_writer(nullptr, &hts_close);
    std::unique_ptr<samFile, decltype(&hts_close)> ldc_writer(nullptr , &hts_close);
    std::unique_ptr<samFile, decltype(&hts_close)> lp_writer(nullptr, &hts_close);
    std::unique_ptr<samFile, decltype(&hts_close)> ow_writer(nullptr, &hts_close);
    std::unique_ptr<samFile, decltype(&hts_close)> ss_writer(nullptr, &hts_close);
    std::ofstream lp_mateseqs_fout; 
    std::ofstream ow_mateseqs_fout;
    std::ofstream ss_mateseqs_fout;

    int curr_pos = 0;
    std::vector<uint32_t> local_depths(rnd_positions.size());
    std::vector<std::vector<uint32_t> > rnd_positions_dist_between_end_and_rnd(rnd_positions.size());
    uint64_t contig_qual_counts[256];
	std::fill(contig_qual_counts, contig_qual_counts+256, uint64_t(0));
    std::vector<uint32_t> local_isize_counts(stats.max_is+1);
    std::vector<std::vector<uint32_t> > isize_dist_by_rndpos(rnd_positions.size());
    uint64_t sum_is = 0;
    uint32_t n_is = 0;
    
    std::unique_ptr<bam1_t, decltype(&bam_destroy1)> read_ptr(bam_init1(), &bam_destroy1);
    bam1_t* read = read_ptr.get();
    int read_status;
    while ((read_status = sam_itr_next(bam_file.file, iter.get(), read)) >= 0) {
        if (!is_primary(read)) continue;

        if (is_dc_pair(read)) {
            if (is_stable_end(read, config)) {
                if (bam_is_rev(read)) {
                    if (!ldc_writer) ldc_writer.reset(open_writer_mt(workspace + "/rev-stable/" + std::to_string(contig_id) + ".noremap.bam", bam_file.header, &hts_pool));
                    int ok = sam_write1(ldc_writer.get(), bam_file.header, read);
                    if (ok < 0) throw std::runtime_error("Failed to write to " + std::string(ldc_writer->fn));
                } else {
                    if (!rdc_writer) rdc_writer.reset(open_writer_mt(workspace + "/fwd-stable/" + std::to_string(contig_id) + ".noremap.bam", bam_file.header, &hts_pool));
                    int ok = sam_write1(rdc_writer.get(), bam_file.header, read);
                    if (ok < 0) throw std::runtime_error("Failed to write to " + std::string(rdc_writer->fn));
                }
            }
            if (!is_mate_unmapped(read)) { // save read seq for remapping
                std::string qname = bam_get_qname(read), read_seq = get_sequence(read, true), qual_ascii = get_qual_ascii(read, true);
                if (is_samechr(read)) {
                    if (read->core.flag & BAM_FREAD1) qname += "_1";
                    else qname += "_2";
                }

                std::lock_guard<std::mutex> lock(mtx_contig[read->core.mtid]);
                mate_seqs[read->core.mtid].push_back(qname + " " + read_seq + " " + qual_ascii + " " + std::to_string(read->core.qual));
                if (mate_seqs[read->core.mtid].size() > 1000) {
                    std::ofstream mate_seqs_fout = open_mateseqs_fout(read->core.mtid);
                    for (std::string& mate_seq : mate_seqs[read->core.mtid]) {
                        mate_seqs_fout << mate_seq << std::endl;
                    }
                    mate_seqs[read->core.mtid].clear();
                }
            }
        }
        if (is_samechr(read)) {
        	if (is_long(read, stats.max_is)) {
                if (read->core.isize > 0) {
                    if (!lp_writer) lp_writer.reset(open_writer_mt(workspace + "/long-pairs/" + std::to_string(contig_id) + ".bam", bam_file.header, &hts_pool));

                    int ok = sam_write1(lp_writer.get(), bam_file.header, read);
                    if (ok < 0) throw std::runtime_error("Failed to write to " + std::string(lp_writer->fn));
                } else {
                    if (!lp_mateseqs_fout.is_open()) lp_mateseqs_fout.open(workspace + "/long-pairs/" + std::to_string(contig_id) + ".txt");
                    lp_mateseqs_fout << bam_get_qname(read) << " " << get_sequence(read) << " " << std::to_string(get_nm(read)) << "\n";
                }
            } else if (is_outward(read)) {
                if (read->core.isize > 0) {
                    if (!ow_writer) ow_writer.reset(open_writer_mt(workspace + "/outward-pairs/" + std::to_string(contig_id) + ".bam", bam_file.header, &hts_pool));

                    int ok = sam_write1(ow_writer.get(), bam_file.header, read);
                    if (ok < 0) throw std::runtime_error("Failed to write to " + std::string(ow_writer->fn));
                } else {
                    if (!ow_mateseqs_fout.is_open()) ow_mateseqs_fout.open(workspace + "/outward-pairs/" + std::to_string(contig_id) + ".txt");
                    ow_mateseqs_fout << bam_get_qname(read) << " " << get_sequence(read) << " " << std::to_string(get_nm(read)) << "\n";
                }
            } else if (is_samestr(read)) {
                if (read->core.pos < read->core.mpos) {
                    if (!ss_writer) ss_writer.reset(open_writer_mt(workspace + "/same-strand/" + std::to_string(contig_id) + ".bam", bam_file.header, &hts_pool));

                    int ok = sam_write1(ss_writer.get(), bam_file.header, read);
                    if (ok < 0) throw std::runtime_error("Failed to write to " + std::string(ss_writer->fn));
                } else {
                    if (!ss_mateseqs_fout.is_open()) ss_mateseqs_fout.open(workspace + "/same-strand/" + std::to_string(contig_id) + ".txt");
                    ss_mateseqs_fout << bam_get_qname(read) << " " << get_sequence(read) << " " << std::to_string(get_nm(read)) << "\n";
                }
        	}
        }

        if (is_unmapped(read)) continue;

        coverage.add_alignment(read->core.pos, bam_endpos(read));

        while (curr_pos < rnd_positions.size() && read->core.pos > rnd_positions[curr_pos]) curr_pos++;

        // sample depth
        bool sampled = false;
		hts_pos_t read_endpos = bam_endpos(read);
		for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < read_endpos; i++) {
			if (read->core.pos <= rnd_positions[i] && rnd_positions[i] <= read_endpos) {
                sampled = true;
				local_depths[i]++;
			}
		}
        if (sampled && read->core.qual >= config.high_confidence_mapq) contig_qual_counts[int(avg_qual(read)+0.5)]++;
        
        // sample pairs crossing
		hts_pos_t pair_startpos = read->core.pos + read->core.l_qseq/2, pair_endpos = get_mate_endpos(read) - read->core.l_qseq/2;
		if (read->core.isize > 0 && read->core.isize <= stats.max_is && is_samechr(read) && !is_samestr(read)
			&& !is_left_clipped(read, config.min_clip_len) && !is_right_clipped(read, config.min_clip_len)) {
			for (int i = curr_pos; i < rnd_positions.size() && rnd_positions[i] < pair_endpos; i++) {
				 if (pair_startpos <= rnd_positions[i] && rnd_positions[i] <= pair_endpos) {
                    hts_pos_t dist = std::min(pair_endpos-rnd_positions[i], (hts_pos_t) stats.max_is);
                    if (dist >= rnd_positions_dist_between_end_and_rnd[i].size()) rnd_positions_dist_between_end_and_rnd[i].resize(dist+1);
					rnd_positions_dist_between_end_and_rnd[i][dist]++;
                    local_isize_counts[read->core.isize]++;
                    isize_dist_by_rndpos[i].push_back(read->core.isize);
                    sum_is += read->core.isize;
                    n_is++;
				}
			}
		}

		if (is_dc_pair(read) && !is_stable_end(read, config)) continue;
		if (has_sequencing_3prime_poly_g_clip(read)) continue;
        if (is_left_clipped(read, config.min_clip_len) || is_right_clipped(read, config.min_clip_len)) {
			if (!sr_writer) {
				sr_writer.reset(open_writer_mt(workspace + "/sr/" + std::to_string(contig_id) + ".bam", bam_file.header, &hts_pool));
				if (sam_idx_init(sr_writer.get(), bam_file.header, 0, sr_index_fname.c_str()) < 0) throw std::runtime_error("Failed to initialize index for " + std::string(sr_writer->fn));
			}

			int ok = sam_write1(sr_writer.get(), bam_file.header, read);
			if (ok < 0) throw std::runtime_error("Failed to write to " + std::string(sr_writer->fn));
		} else if (is_hidden_split_read(read, config)) {
            if (!hsr_writer) {
                hsr_writer.reset(open_writer_mt(workspace + "/hsr/" + std::to_string(contig_id) + ".bam", bam_file.header, &hts_pool));
                if (sam_idx_init(hsr_writer.get(), bam_file.header, 0, hsr_index_fname.c_str()) < 0) throw std::runtime_error("Failed to initialize index for " + std::string(hsr_writer->fn));
            }

            int ok = sam_write1(hsr_writer.get(), bam_file.header, read);
            if (ok < 0) throw std::runtime_error("Failed to write to " + std::string(hsr_writer->fn));
        } 
    }

    if (read_status < -1) throw std::runtime_error("Failed to read " + bam_fname + " on " + contig_name);

    for (auto writer : {&sr_writer, &hsr_writer}) {
        if (!*writer) continue;
        std::string fname = (*writer)->fn;
        int index_status = sam_idx_save(writer->get());
        int close_status = sam_close(writer->release());
        if (index_status < 0) throw std::runtime_error("Failed to save index for " + fname);
        if (close_status < 0) throw std::runtime_error("Failed to close " + fname);
    }

    lp_mateseqs_fout.close();
    ow_mateseqs_fout.close();
    ss_mateseqs_fout.close();

    coverage.finish();
    write_coarse_coverage(coverage_fname, coverage.maxima);

    local_depths.erase(std::remove(local_depths.begin(), local_depths.end(), 0), local_depths.end());

    std::lock_guard<std::mutex> lock(mtx);
    if (local_depths.size() >= MIN_RND_POS) {
		std::sort(local_depths.begin(), local_depths.end());
		min_depth_by_contig[contig_name] = local_depths[local_depths.size()/100];
		median_depth_by_contig[contig_name] = local_depths[local_depths.size()/2];
		max_depth_by_contig[contig_name] = local_depths[local_depths.size()-local_depths.size()/100];
    }
    depths.insert(depths.end(), local_depths.begin(), local_depths.end());

    partial_sums.push_back({sum_is, n_is});
	for (int i = 0; i < 256; i++) qual_counts[i] += contig_qual_counts[i];

    for (int i = 0; i < rnd_positions.size(); i++) {
		std::vector<uint32_t> local_isizes_count_geq_i(stats.max_is+1); // position i contains the count of isizes that are >= than i
		for (uint32_t isize : isize_dist_by_rndpos[i]) local_isizes_count_geq_i[isize]++;
		for (int j = stats.max_is-1; j >= 0; j--) {
			local_isizes_count_geq_i[j] += local_isizes_count_geq_i[j+1];
		}
		if (local_isizes_count_geq_i[0] == 0) continue;
		for (int j = 0; j <= stats.max_is; j++) {
            if (isizes_count_geq_i[j].size() < local_isizes_count_geq_i[j]+1) isizes_count_geq_i[j].resize(local_isizes_count_geq_i[j]+1);
            isizes_count_geq_i[j][local_isizes_count_geq_i[j]]++;
		}

        rnd_positions_dist_between_end_and_rnd[i].resize(stats.max_is+1);
        dist_between_end_and_rnd[0].resize(1);
        dist_between_end_and_rnd[0][0]++;
        for (int j = 1; j <= stats.max_is; j++) {
            rnd_positions_dist_between_end_and_rnd[i][j] += rnd_positions_dist_between_end_and_rnd[i][j-1];
            if (dist_between_end_and_rnd[j].size() <= rnd_positions_dist_between_end_and_rnd[i][j]) dist_between_end_and_rnd[j].resize(rnd_positions_dist_between_end_and_rnd[i][j]+1);
            dist_between_end_and_rnd[j][rnd_positions_dist_between_end_and_rnd[i][j]]++;
        }
        rnd_positions_dist_between_end_and_rnd[i].clear();
    }

    for (int i = 0; i <= stats.max_is; i++) {
    	isize_counts[i] += local_isize_counts[i];
    }
}

void find_1_perc(std::vector<uint32_t>& v, uint32_t& min, uint32_t& max) {
    uint64_t sum = std::accumulate(v.begin(), v.end(), uint64_t(0));
    uint64_t _1_perc_count = sum/100;
    min = 0;
    max = v.size()-1;
    uint64_t curr_sum = 0;
    for (int i = 0; i < v.size(); i++) {
        curr_sum += v[i];
        if (curr_sum >= _1_perc_count) {
            min = i;
            break;
        }
    }
    curr_sum = 0;
    for (int i = v.size()-1; i >= 0; i--) {
        curr_sum += v[i];
        if (curr_sum >= _1_perc_count) {
            max = i;
            break;
        }
    }
}

uint32_t find_median(std::vector<uint32_t>& v) {
    uint64_t sum = std::accumulate(v.begin(), v.end(), uint64_t(0));
    uint64_t curr_sum = 0;
    for (int i = 0; i < v.size(); i++) {
        curr_sum += v[i];
        if (curr_sum >= sum/2) return i;
    }
    return 0;
}

int main(int argc, char* argv[]) {

    std::string bam_fname = argv[1];
    workdir = argv[2];
    workspace = workdir + "/workspace";
    reference_fname = argv[3];

    config.parse(workdir + "/config.txt");

    contig_map_t contig_map(workdir);
    chr_seqs.read_lens_into_map(reference_fname);

    open_samFile_t bam_file(bam_fname.c_str());
	if (hts_set_fai_filename(bam_file.file, fai_path(reference_fname.c_str())) != 0) {
		throw std::runtime_error("Failed to read reference " + reference_fname);
	}

    mtx_contig = new std::mutex[contig_map.size()];
    mate_seqs.resize(contig_map.size());

    // generate random positions
	std::vector<std::pair<std::string, hts_pos_t> > rnd_positions;
    std::unordered_map<std::string, std::vector<hts_pos_t> > rnd_pos_map;
	random_pos_generator_t random_pos_generator(bam_file.header, config.seed, config.sampling_regions);
    int n_rand_pos = std::max(hts_pos_t(1000), random_pos_generator.reference_len/1000);
	for (int i = 0; i < n_rand_pos; i++) {
        std::pair<std::string, hts_pos_t> rnd_pos = random_pos_generator.get_random_pos();
		rnd_positions.push_back(rnd_pos);
        rnd_pos_map[rnd_pos.first].push_back(rnd_pos.second);
	}

    std::mt19937 rng(config.seed);
    std::shuffle(rnd_positions.begin(), rnd_positions.end(), rng);    

    ctpl::thread_pool base_stats_thread_pool(config.threads);
    std::vector<std::future<void> > futures;
    int RND_POS_BLOCK_SIZE = std::max(1, (int) rnd_positions.size()/(config.threads*10));
    for (int i = 0; i < rnd_positions.size(); i+=RND_POS_BLOCK_SIZE) {
        std::vector<std::pair<std::string, hts_pos_t> > rnd_pos_block(rnd_positions.begin()+i, rnd_positions.begin()+std::min(i+RND_POS_BLOCK_SIZE, (int) rnd_positions.size()));
        std::future<void> future = base_stats_thread_pool.push(get_base_stats, bam_fname, reference_fname, rnd_pos_block);
        futures.push_back(std::move(future));
    }
    base_stats_thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    futures.clear();

    // Calculate base stats. These are about the general insert size distribution, not the distribution of insert sizes crossing sampled positions
    int mean_is = mean(general_isize_dist);
    int stddev_is = stddev(general_isize_dist);
    general_isize_dist.erase(std::remove_if(general_isize_dist.begin(), general_isize_dist.end(), [mean_is, stddev_is](int x) { return std::abs(x-mean_is) > 5*stddev_is; }), general_isize_dist.end());

    if (general_isize_dist.empty()) {
        throw std::runtime_error("Failed to estimate reliable insert size statistics: no read pairs found \
            crossing sampled positions. Consider using the --sampling-regions option.");
    }

    mean_is = mean(general_isize_dist);
    stddev_is = stddev(general_isize_dist);

    uint64_t lower_stddev_is = 0;
    int n_vals = 0;
    for (int x : general_isize_dist) {
        if (x < mean_is) {
            lower_stddev_is += (mean_is-x)*(mean_is-x);
            n_vals++;
        }
    }
    lower_stddev_is = std::sqrt(lower_stddev_is / std::max(1, n_vals));

    uint64_t higher_stddev_is = 0;
    n_vals = 0;
    for (int x : general_isize_dist) {
        if (x > mean_is) {
            higher_stddev_is += (x-mean_is)*(x-mean_is);
            n_vals++;
        }
    }
    higher_stddev_is = std::sqrt(higher_stddev_is / std::max(1, n_vals));

    stats.min_is = mean_is - 3*lower_stddev_is;
    stats.max_is = mean_is + 3.5*higher_stddev_is;

    isize_counts.resize(stats.max_is+1);
    isizes_count_geq_i.resize(stats.max_is+1);
    dist_between_end_and_rnd.resize(stats.max_is+1);

    // Categorize reads
    int hts_threads = std::min(8, config.threads / 2);
    int n_workers = config.threads - hts_threads;
    if (hts_threads > 0) {
        hts_pool.pool = hts_tpool_init(hts_threads);
        if (hts_pool.pool == NULL) throw std::runtime_error("Failed to initialize HTSlib thread pool");
    }
    ctpl::thread_pool categorize_thread_pool(n_workers);
    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
        std::string contig_name = contig_map.get_name(contig_id);
        std::future<void> future = categorize_thread_pool.push(categorize, contig_id, contig_name, bam_fname, reference_fname, rnd_pos_map[contig_name]);
        futures.push_back(std::move(future));
    }
    categorize_thread_pool.stop(true);
    for (int i = 0; i < futures.size(); i++) {
        futures[i].get();
    }
    if (hts_pool.pool) hts_tpool_destroy(hts_pool.pool);
    hts_pool.pool = NULL;

    if (depths.size() < MIN_RND_POS) {
        throw std::runtime_error(
            "Failed to estimate reliable depth statistics: only "
            + std::to_string(depths.size())
            + " informative sampled positions found, need at least "
            + std::to_string(MIN_RND_POS)
            + "; the sample is probably extremely sparse. Consider using the --sampling-regions option."
        );
    }

    for (int i = 0; i < contig_map.size(); i++) {
    	if (mate_seqs[i].empty()) continue;
		std::ofstream mate_seqs_fout = open_mateseqs_fout(i);
		for (std::string& mate_seq : mate_seqs[i]) {
			mate_seqs_fout << mate_seq << std::endl;
		}
		mate_seqs_fout.close();
	}

    // get 1-percentile of base qualities
    uint64_t tot_quals = std::accumulate(qual_counts, qual_counts+256, uint64_t(0));
    uint64_t _1_perc_count = tot_quals/100;
    int min_avg_base_qual = 0;
    for (int i = 0; i < 256; i++) {
    	if (_1_perc_count <= qual_counts[i]) {
    		min_avg_base_qual = i;
    		break;
    	}
    	_1_perc_count -= qual_counts[i];
    }

    uint64_t sum_is = 0;
    uint32_t n_is = 0;
    for (std::pair<uint64_t, uint32_t>& p : partial_sums) {
        sum_is += p.first;
        n_is += p.second;
    }

    if (n_is == 0) {
        throw std::runtime_error("Failed to estimate reliable insert size statistics: no read pairs\
             found crossing sampled positions. Consider using the --sampling-regions option.");
    }

    std::ofstream stats_out(workdir + "/stats.txt");
    stats_out << "read_len . " << stats.read_len << std::endl;
    stats_out << "min_is . " << stats.min_is << std::endl;
    stats_out << "max_is . " << stats.max_is << std::endl;
	std::sort(depths.begin(), depths.end());
	stats_out << "min_depth . " << depths[depths.size()/100] << std::endl;
	stats_out << "median_depth . " << depths[depths.size()/2] << std::endl;
	stats_out << "max_depth . " << depths[depths.size()-depths.size()/100] << std::endl;
    stats_out << "min_avg_base_qual . " << min_avg_base_qual << std::endl;
    for (int contig_id = 0; contig_id < contig_map.size(); contig_id++) {
		std::string contig_name = contig_map.get_name(contig_id);
		if (!min_depth_by_contig.count(contig_name)) continue;
    	stats_out << "min_depth " << contig_name << " " << min_depth_by_contig[contig_name] << std::endl;
    	stats_out << "median_depth " << contig_name << " " << median_depth_by_contig[contig_name] << std::endl;
    	stats_out << "max_depth " << contig_name << " " << max_depth_by_contig[contig_name] << std::endl;
    }
    stats_out << "pop_avg_crossing_is . " << sum_is / n_is << std::endl;

    std::ofstream crossing_isizes_dist_fout(workdir + "/crossing_isizes.txt");
    for (int i = 0; i <= stats.max_is; i++) {
    	crossing_isizes_dist_fout << i << " " << isize_counts[i] << std::endl;
    }
    crossing_isizes_dist_fout.close();

	for (int i = 0; i <= stats.max_is; i++) {
        uint32_t min_disc_pairs_by_insertion_size = 0, max_disc_pairs_by_insertion_size = 0;
        find_1_perc(dist_between_end_and_rnd[i], min_disc_pairs_by_insertion_size, max_disc_pairs_by_insertion_size);
        stats_out << "min_disc_pairs_by_insertion_size " << i << " " << min_disc_pairs_by_insertion_size << std::endl;
        stats_out << "max_disc_pairs_by_insertion_size " << i << " " << max_disc_pairs_by_insertion_size << std::endl;

        uint32_t min_pairs_crossing_gap = 0, max_pairs_crossing_gap = 0;
        find_1_perc(isizes_count_geq_i[i], min_pairs_crossing_gap, max_pairs_crossing_gap);
        uint32_t median_pairs_crossing_gap = find_median(isizes_count_geq_i[i]);
        stats_out << "min_pairs_crossing_gap " << i << " " << min_pairs_crossing_gap << std::endl;
		stats_out << "median_pairs_crossing_gap " << i << " " << median_pairs_crossing_gap << std::endl;
        stats_out << "max_pairs_crossing_gap " << i << " " << max_pairs_crossing_gap << std::endl;
	}
	stats_out.close();
}
