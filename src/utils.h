#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <queue>
#include <numeric>
#include <cmath>
#include <random>
#include <unistd.h>
#include <algorithm>
#include <cstring>

#include <htslib/sam.h>
#include "htslib/hts.h"
#include "htslib/kseq.h"
KSEQ_INIT(int, read)

struct config_t {

    int threads, seed;
    int min_clip_len, min_stable_mapq, min_diff_hsr;
    int high_confidence_mapq;
    int min_sv_size, max_trans_size;
    double max_seq_error;
    bool per_contig_stats;
    bool training_mode;
    std::string sampling_regions, version;

    const int flanking_size = 5000, indel_tested_region_size = 10000;
    const int min_assembled_ins_size = 50;

    void parse(std::string config_file) {
        std::unordered_map<std::string, std::string> config_params;
        std::ifstream fin(config_file);
        std::string name, value;
        while (fin >> name >> value) {
            config_params[name] = value;
        }
        fin.close();

        threads = std::stoi(config_params["threads"]);
        seed = std::stoi(config_params["seed"]);
        min_clip_len = std::stoi(config_params["min_clip_len"]);
        min_stable_mapq = std::stoi(config_params["min_stable_mapq"]);
        high_confidence_mapq = std::stoi(config_params["high_confidence_mapq"]);
        min_diff_hsr = std::stoi(config_params["min_diff_hsr"]);
        min_sv_size = std::stoi(config_params["min_sv_size"]);
        if (config_params.count("max_trans_size")) max_trans_size = std::stoi(config_params["max_trans_size"]);
        max_seq_error = std::stod(config_params["max_seq_error"]);
        per_contig_stats = std::stoi(config_params["per_contig_stats"]);
        training_mode = std::stoi(config_params["training_mode"]);
        sampling_regions = config_params["sampling_regions"];
        version = config_params["version"];
    }
};

struct stats_t {

    int min_is, max_is;
    int read_len;
    std::unordered_map<std::string, int> min_depths, median_depths, max_depths;
    std::unordered_map<int, int> min_disc_pairs_by_insertion_size, max_disc_pairs_by_insertion_size;
    int min_avg_base_qual;
    int pop_avg_crossing_is = 0;
    bool per_contig_stats = false;

    void parse(std::string stats_file, bool per_contig_stats) {
        std::unordered_map<std::string, std::unordered_map<std::string, std::string>> stats_params;
        std::ifstream fin(stats_file);
        std::string stat_name, stat_subname, value;
        while (fin >> stat_name >> stat_subname >> value) {
            stats_params[stat_name][stat_subname] = value;
            if (stat_name == "min_is" && stat_subname == ".") min_is = std::stoi(value);
            if (stat_name == "max_is" && stat_subname == ".") max_is = std::stoi(value);
            if (stat_name == "read_len" && stat_subname == ".") read_len = std::stoi(value);
            if (stat_name == "min_depth") min_depths[stat_subname] = std::stoi(value);
            if (stat_name == "median_depth") median_depths[stat_subname] = std::stoi(value);
            if (stat_name == "max_depth") max_depths[stat_subname] = std::stoi(value);
            if (stat_name == "min_disc_pairs_by_insertion_size") min_disc_pairs_by_insertion_size[std::stoi(stat_subname)] = std::stoi(value);
            if (stat_name == "max_disc_pairs_by_insertion_size") max_disc_pairs_by_insertion_size[std::stoi(stat_subname)] = std::stoi(value);
            if (stat_name == "min_avg_base_qual") min_avg_base_qual = std::stoi(value);
            if (stat_name == "pop_avg_crossing_is" && stat_subname == ".") pop_avg_crossing_is = std::stoi(value);
        }
        fin.close();
        this->per_contig_stats = per_contig_stats;
    }

    int get_min_depth(std::string contig_name) {
        if (per_contig_stats && min_depths.count(contig_name))
            return min_depths[contig_name];
        else return min_depths["."];
    }
    int get_median_depth(std::string contig_name) {
        if (per_contig_stats && median_depths.count(contig_name))
            return median_depths[contig_name];
        else return median_depths["."];
    }
    int get_max_depth(std::string contig_name) {
        if (per_contig_stats && max_depths.count(contig_name))
            return max_depths[contig_name];
        else return max_depths["."];
    }

    int get_min_disc_pairs_by_insertion_size(int is) {
        if (min_disc_pairs_by_insertion_size.count(is))
            return min_disc_pairs_by_insertion_size[is];
        else return min_disc_pairs_by_insertion_size[max_is];
    }
    int get_max_disc_pairs_by_insertion_size(int is) {
        if (max_disc_pairs_by_insertion_size.count(is))
            return max_disc_pairs_by_insertion_size[is];
        else return max_disc_pairs_by_insertion_size[max_is];
    }
};

struct contig_map_t {

    std::unordered_map<std::string, size_t> name_to_id;
    std::vector<std::string> id_to_name;

    contig_map_t() {}
    contig_map_t(std::string workdir) {
        load(workdir);
    }

    void load(std::string workdir) {
    	std::ifstream fin(workdir + "/contig_map");
		std::string name;
		int id = 0;
		while (fin >> name) {
			name_to_id[name] = id;
			id_to_name.push_back(name);
			id++;
		}
    }

    size_t size() {return id_to_name.size();}
    std::string get_name(size_t id) {return id_to_name[id];};
    size_t get_id(std::string& name) {return name_to_id[name];};
};

template<typename T>
double mean(std::vector<T>& v) {
    if (v.empty()) return 0;
    return double(std::accumulate(v.begin(), v.end(), 0.0))/v.size();
}

template<typename T>
double stddev(std::vector<T>& v) {
    if (v.empty()) return 0;
    double m = mean(v);
    double sum = 0;
    for (T& e : v) {
        sum += (e - m)*(e - m);
    }
    return std::sqrt(sum/v.size());
}

template<typename T>
inline T max(T a, T b, T c) { return std::max(a, std::max(b,c)); }

template<typename T>
inline T max(T a, T b, T c, T d) { return std::max(std::max(a,b), std::max(c,d)); }

bool file_exists(std::string& fname) {
	return std::ifstream(fname).good();
}

char toupper(char c) {
    if (c >= 'a' && c <= 'z') {
        return c - ('a' - 'A');
    }
    return c;
}

void to_uppercase(char* s) {
    while (*s) {
        *s = toupper(*s);
        s++;
    }
}

struct chr_seqs_map_t {
    struct chr_seq_t {
        char* seq;
        hts_pos_t len;

        chr_seq_t(char* seq, hts_pos_t len) : seq(seq), len(len) {}
        ~chr_seq_t() {delete[] seq;}
    };

    std::unordered_map<std::string, chr_seq_t*> seqs;
    std::vector<std::string> ordered_contigs;

    void read_fasta_into_map(std::string& reference_fname, bool uppercase_reference = true) {
        ordered_contigs.clear();
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            char* chr_seq = new char[seq->seq.l + 1];
            strncpy(chr_seq, seq->seq.s, seq->seq.l);
            chr_seq[seq->seq.l] = '\0';
            if (uppercase_reference) to_uppercase(chr_seq);
            if (seqs.count(seq_name)) delete seqs[seq_name];
            seqs[seq_name] = new chr_seq_t(chr_seq, seq->seq.l);
            ordered_contigs.push_back(seq_name);
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    void read_lens_into_map(std::string& reference_fname) {
        ordered_contigs.clear();
        FILE* fasta = fopen(reference_fname.c_str(), "r");
        kseq_t* seq = kseq_init(fileno(fasta));
        while (kseq_read(seq) >= 0) {
            std::string seq_name = seq->name.s;
            if (seqs.count(seq_name)) delete seqs[seq_name];
            seqs[seq_name] = new chr_seq_t(NULL, seq->seq.l);
            ordered_contigs.push_back(seq_name);
        }
        kseq_destroy(seq);
        fclose(fasta);
    }

    char* get_seq(std::string seq_name) {
        if (!seqs.count(seq_name)) {
            throw std::runtime_error("Sequence " + seq_name + " not found in reference.");
        }
        return seqs[seq_name]->seq;
    }

    hts_pos_t get_len(std::string seq_name) {
        if (!seqs.count(seq_name)) {
            throw std::runtime_error("Sequence " + seq_name + " not found in reference.");
        }
        return seqs[seq_name]->len;
    }

    void clear() {
        for (auto& e : seqs) {
            delete e.second;
            e.second = NULL;
        }
    }

    ~chr_seqs_map_t() {
        clear();
    }
};

struct base_frequencies_t {
    int a = 0, c = 0, g = 0, t = 0;

    base_frequencies_t() {}
    base_frequencies_t(int a, int c, int g, int t) : a(a), c(c), g(g), t(t) {}

    int max() {
        return std::max(std::max(a, c), std::max(g, t));
    }
    int tot() {
        return a + c + g + t;
    }
    double max_freq() {
        return max()/double(std::max(1, tot()));
    }

    bool empty() {
        return a + c + g + t == 0;
    }
};

// Map byte -> {0,1,2,3} for A,C,G,T, else 4.
alignas(64)
static const uint8_t nt_map[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

base_frequencies_t get_base_frequencies(const char* seq, int len) {
    int counts[5] = {0,0,0,0,0};
    for (int i = 0; i < len; i++) {
        counts[nt_map[(uint8_t)seq[i]]]++;
    }
    // ignore counts[4] which is for non-ACGT
    return base_frequencies_t(counts[0], counts[1], counts[2], counts[3]);
}

std::tuple<base_frequencies_t, base_frequencies_t, base_frequencies_t> get_left_flanking_base_frequencies_50_100_500(const char* contig_seq, hts_pos_t pos) {
    int counts[5] = {0,0,0,0,0};

    hts_pos_t start50 = std::max(hts_pos_t(0), pos - 50);
    for (hts_pos_t i = start50; i < pos; i++) {
        counts[nt_map[(uint8_t)contig_seq[i]]]++;
    }
    base_frequencies_t bf50 = base_frequencies_t(counts[0], counts[1], counts[2], counts[3]);

    hts_pos_t start100 = std::max(hts_pos_t(0), pos - 100);
    for (hts_pos_t i = start100; i < start50; i++) {
        counts[nt_map[(uint8_t)contig_seq[i]]]++;
    }
    base_frequencies_t bf100 = base_frequencies_t(counts[0], counts[1], counts[2], counts[3]);

    hts_pos_t start500 = std::max(hts_pos_t(0), pos - 500);
    for (hts_pos_t i = start500; i < start100; i++) {
        counts[nt_map[(uint8_t)contig_seq[i]]]++;
    }
    base_frequencies_t bf500 = base_frequencies_t(counts[0], counts[1], counts[2], counts[3]);

    return std::make_tuple(bf50, bf100, bf500);
}

std::tuple<base_frequencies_t, base_frequencies_t, base_frequencies_t> get_right_flanking_base_frequencies_50_100_500(const char* contig_seq, hts_pos_t pos, hts_pos_t contig_len) {
    int counts[5] = {0,0,0,0,0};

    hts_pos_t end = std::min(contig_len, pos + 50);
    for (hts_pos_t i = pos; i < end; i++) {
        counts[nt_map[(uint8_t)contig_seq[i]]]++;
    }
    base_frequencies_t bf50 = base_frequencies_t(counts[0], counts[1], counts[2], counts[3]);

    hts_pos_t end100 = std::min(contig_len, pos + 100);
    for (hts_pos_t i = end; i < end100; i++) {
        counts[nt_map[(uint8_t)contig_seq[i]]]++;
    }
    base_frequencies_t bf100 = base_frequencies_t(counts[0], counts[1], counts[2], counts[3]);

    hts_pos_t end500 = std::min(contig_len, pos + 500);
    for (hts_pos_t i = end100; i < end500; i++) {
        counts[nt_map[(uint8_t)contig_seq[i]]]++;
    }
    base_frequencies_t bf500 = base_frequencies_t(counts[0], counts[1], counts[2], counts[3]);

    return std::make_tuple(bf50, bf100, bf500);
}

bool is_homopolymer(const char* seq, int len, double threshold = 0.8) {
    return get_base_frequencies(seq, len).max_freq() >= threshold;
}
bool is_homopolymer(std::string seq, double threshold = 0.8) {
	return is_homopolymer(seq.data(), seq.length(), threshold);
}

// hp_prefix must be allocated and of size (at least) len
void is_homopolymer_prefix(const char* seq, int len, bool* hp_prefix) {
    int max_idx = 0;
    int counts[5] = {0,0,0,0,0};
    for (int i = 0; i < len; i++) {
        int curr_idx = nt_map[(uint8_t)seq[i]];
        counts[curr_idx]++;
        if (counts[curr_idx] > counts[max_idx]) {
            max_idx = curr_idx;
        }
        hp_prefix[i] = (counts[max_idx] >= (i+1)*0.8);
    }
}
void is_homopolymer_suffix(const char* seq, int len, bool* hp_suffix) {
    int max_idx = 0;
    int counts[5] = {0,0,0,0,0};
    for (int i = len-1; i >= 0; i--) {
        int curr_idx = nt_map[(uint8_t)seq[i]];
        counts[curr_idx]++;
        if (counts[curr_idx] > counts[max_idx]) {
            max_idx = curr_idx;
        }
        hp_suffix[i] = (counts[max_idx] >= (len - i)*0.8);
    }
}

int find_char_in_str(std::string s, char c, int start = 0, int stop = -1) {
    if (start < 0) start = 0;
    if (stop == -1 || stop > s.length()) stop = s.length();
    for (int i = start; i < stop; i++) {
        if (s[i] == c) return i;
    }
    return -1;
}

int64_t overlap(hts_pos_t s1, hts_pos_t e1, hts_pos_t s2, hts_pos_t e2) {
    int64_t overlap = std::min(e1, e2) - std::max(s1, s2);
    return std::max(int64_t(0), overlap);
}

template<typename T>
T bring_within_range(T val, T min_val, T max_val) {
    if (val < min_val) return min_val;
    if (val > max_val) return max_val;
    return val;
}

struct edge_t {
	int next, score, overlap;

	edge_t() : next(0), score(0), overlap(0) {}
	edge_t(int next, int score, int overlap) : next(next), score(score), overlap(overlap) {}
};

std::vector<int> find_rev_topological_order(int n, std::vector<int> out_edges, std::vector<std::vector<edge_t> >& l_adj_rev) {

	std::queue<int> sinks;
	for (int i = 0; i < n; i++) {
		if (!out_edges[i]) sinks.push(i);
	}

	std::vector<int> rev_topological_order;
	while (!sinks.empty()) {
		int s = sinks.front();
		sinks.pop();
		rev_topological_order.push_back(s);
		for (edge_t& e : l_adj_rev[s]) {
			out_edges[e.next]--;
			if (out_edges[e.next] == 0) sinks.push(e.next);
		}
	}
	return rev_topological_order;
}

struct random_pos_generator_t {

    const int NO_SAMPLING_PADDING = 1000;

    hts_pos_t reference_len = 0;
    std::mt19937 rng;
    std::uniform_int_distribution<hts_pos_t> dist;

    struct region_t {
        std::string chr;
        hts_pos_t start, end;
        region_t(std::string chr, hts_pos_t start, hts_pos_t end) : chr(chr), start(start), end(end) {}
    };
    std::vector<region_t> regions;

    random_pos_generator_t(bam_hdr_t* hdr, int seed, std::string sampling_regions_fname = "") {
        rng.seed(seed);

        if (sampling_regions_fname.empty()) {
            for (int i = 0; i < hdr->n_targets; i++) {
                std::string chr = hdr->target_name[i];
                if (hdr->target_len[i] <= 2*NO_SAMPLING_PADDING) continue;
                regions.push_back(region_t(chr, NO_SAMPLING_PADDING, hdr->target_len[i] - NO_SAMPLING_PADDING));
            }
        } else {
            std::ifstream fin(sampling_regions_fname);
            // check if file exists and is not empty
            if (!fin.good()) {
                std::string error_msg = "Error: sampling regions file " + sampling_regions_fname + " does not exist or is not readable\n";
                throw std::runtime_error(error_msg);
            }

            std::string line;
            while (std::getline(fin, line)) {
                std::istringstream iss(line);

                if (line.empty() || line[0] == '#') continue; // skip empty lines and comments

                // tokenize iss imitating >>, do not assume format
                std::vector<std::string> tokens;
                std::string token;
                while (iss >> token) {
                    tokens.push_back(token);
                }

                if (tokens.size() == 3) {
                    std::string chr = tokens[0];
                    hts_pos_t start = std::stoll(tokens[1]);
                    hts_pos_t end = std::stoll(tokens[2]);
                    regions.push_back(region_t(chr, start, end));
                } else {
                    std::string error_msg = "Error: invalid sampling region format in " + sampling_regions_fname + "\n";
                    throw std::runtime_error(error_msg);
                }
            }
        }
        regions.erase(std::remove_if(regions.begin(), regions.end(), [](region_t& r) {return r.start >= r.end;}), regions.end());

        for (region_t& r : regions) {
            reference_len += r.end - r.start;
        }

        if (reference_len == 0) {
            throw std::runtime_error("Error: the total length of sampling regions is 0.");
        }
        dist = std::uniform_int_distribution<hts_pos_t>(0, reference_len-1);
    }

    std::pair<std::string, hts_pos_t> get_random_pos() {
        hts_pos_t random_pos = dist(rng);

        for (region_t& r : regions) {
            if (random_pos < r.end - r.start) {
                return {r.chr, r.start + random_pos};
            } else {
                random_pos -= r.end - r.start;
            }
        }
        throw std::runtime_error("Error: random_pos_generator_t::get_random_pos() failed\n");
    }
};

bool is_genomic_string(std::string s) {
    for (char c : s) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
            c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return false;
        }
    }
    return true;
}

char* concat2(const char* s1, const char* s2, size_t len1, size_t len2) {
    char* result = new char[len1 + len2 + 1];
    memcpy(result, s1, len1);
    memcpy(result + len1, s2, len2);
    result[len1 + len2] = '\0';
    return result;
}

char* concat3(const char* s1, const char* s2, const char* s3, size_t len1, size_t len2, size_t len3) {
    char* result = new char[len1 + len2 + len3 + 1];
    memcpy(result, s1, len1);
    memcpy(result + len1, s2, len2);
    memcpy(result + len1 + len2, s3, len3);
    result[len1 + len2 + len3] = '\0';
    return result;
}

#endif
