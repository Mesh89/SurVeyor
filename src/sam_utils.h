#ifndef SAM_UTILS_H
#define SAM_UTILS_H

#include <iostream>
#include <algorithm>
#include <mutex>
#include <memory>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "utils.h"

bool is_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FUNMAP;
}
bool is_mate_unmapped(bam1_t* r) {
    return r->core.flag & BAM_FMUNMAP;
}
bool is_primary(bam1_t* r) {
    return !(r->core.flag & BAM_FSECONDARY) && !(r->core.flag & BAM_FSUPPLEMENTARY);
}
bool is_samechr(bam1_t* r) {
    return r->core.tid == r->core.mtid && !is_unmapped(r) && !is_mate_unmapped(r);
}
bool is_samestr(bam1_t* r) {
    return bam_is_rev(r) == bam_is_mrev(r);
}
bool is_dc_pair(bam1_t* r) {
    return !is_samechr(r) || std::abs(r->core.isize) > 100000 || is_unmapped(r) != is_mate_unmapped(r);
}
bool is_outward(bam1_t* r) {
	return is_samechr(r) && !is_samestr(r) && ((!bam_is_rev(r) && r->core.isize < 0) || (bam_is_rev(r) && r->core.isize > 0));
}
bool is_short(bam1_t* r, int min_is) {
    return is_samechr(r) && !is_samestr(r) && ((!bam_is_rev(r) && r->core.isize < min_is) || (bam_is_rev(r) && r->core.isize > -min_is));
}
bool is_long(bam1_t* r, int max_is) {
	return is_samechr(r) && !is_samestr(r) && ((!bam_is_rev(r) && r->core.isize > max_is) || (bam_is_rev(r) && r->core.isize < -max_is));
}
bool is_proper_pair(bam1_t* r, int min_is, int max_is) {
	return is_primary(r) && is_samechr(r) && !is_samestr(r) && !is_dc_pair(r) && !is_outward(r) && !is_short(r, min_is) && !is_long(r, max_is);
}
bool is_first_read(bam1_t* r) {
    return r->core.flag & BAM_FREAD1;
}

std::string get_mate_lookup_qname(bam1_t* read) {
    std::string qname = bam_get_qname(read);
    if (is_samechr(read)) {
        if (read->core.flag & BAM_FREAD1) qname += "_2";
        else qname += "_1";
    }
    return qname;
}

int get_endpoint(bam1_t* r) {
    return bam_is_rev(r) ? r->core.pos : bam_endpos(r);
}

int get_left_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[0]) == 'S' ? bam_cigar_oplen(cigar[0]): 0;
}
int get_right_clip_size(bam1_t* r) {
    uint32_t* cigar = bam_get_cigar(r);
    return bam_cigar_opchr(cigar[r->core.n_cigar-1]) == 'S' ? bam_cigar_oplen(cigar[r->core.n_cigar-1]): 0;
}
bool is_left_clipped(bam1_t* r, int min_clip_len) {
	if (is_unmapped(r)) return false;
    return get_left_clip_size(r) >= min_clip_len;
}
bool is_right_clipped(bam1_t* r, int min_clip_len) {
	if (is_unmapped(r)) return false;
    return get_right_clip_size(r) >= min_clip_len;
}
bool is_clipped(bam1_t* r, int min_clip_len) {
    return is_left_clipped(r, min_clip_len) || is_right_clipped(r, min_clip_len);
}

hts_pos_t get_unclipped_start(bam1_t* r) {
    return r->core.pos - get_left_clip_size(r);
}
hts_pos_t get_unclipped_end(bam1_t* r) {
    return bam_endpos(r) + get_right_clip_size(r);
}

int64_t get_mq(bam1_t* r) {
    if (is_unmapped(r) || is_mate_unmapped(r)) {
        return 0;
    }
    uint8_t* mq = bam_aux_get(r, "MQ");
    if (mq == NULL) {
        throw std::runtime_error("Read " + std::string(bam_get_qname(r)) + " does not have the MQ tag.");
    }
    return bam_aux2i(mq);
}

char* get_mc(bam1_t* r) {
    static char empty_mc[] = "";
    if (is_unmapped(r) || is_mate_unmapped(r)) {
        return empty_mc;
    }
    uint8_t* mc = bam_aux_get(r, "MC");
    if (mc == NULL) {
        throw std::runtime_error("Read " + std::string(bam_get_qname(r)) + " does not have the MC tag.");
    }
    return bam_aux2Z(mc);
}

int64_t get_nm(bam1_t* r) {
    uint8_t* nm = bam_aux_get(r, "NM");
    if (nm == NULL) {
        throw std::runtime_error("Read " + std::string(bam_get_qname(r)) + " does not have the NM tag.");
    }
    return bam_aux2i(nm);
}

struct indel_summary_t {
    int ops = 0;
    int dels = 0;
    int inss = 0;
};

indel_summary_t get_indel_summary(bam1_t* r) {
    indel_summary_t summary;
    uint32_t* cigar = bam_get_cigar(r);
    for (uint32_t i = 0; i < r->core.n_cigar; i++) {
        char op = bam_cigar_opchr(cigar[i]);
        if (op == 'D' || op == 'I') {
            int len = bam_cigar_oplen(cigar[i]);
            summary.ops++;
            if (op == 'D') summary.dels += len;
            else summary.inss += len;
        }
    }
    return summary;
}

char* get_md(bam1_t* r) {
    uint8_t* md = bam_aux_get(r, "MD");
    if (md == NULL) {
        throw std::runtime_error("Read " + std::string(bam_get_qname(r)) + " does not have the MD tag.");
    }
    return bam_aux2Z(md);
}

int get_mate_left_clip_size(bam1_t* r) {
    if (is_mate_unmapped(r)) return 0;

    char* mc = get_mc(r);
    int i = 0, left_clip = 0;
    while (mc[i] >= '0' && mc[i] <= '9') {
        left_clip = left_clip * 10 + (mc[i] - '0');
        i++;
    }
    return (mc[i] == 'S') ? left_clip : 0;
}

int get_mate_right_clip_size(bam1_t* r) {
    if (is_mate_unmapped(r)) return 0;

    char* mc = get_mc(r);
    int len = strlen(mc), i = len - 1;
    if (mc[i] != 'S') return 0;
    i--;
    int right_clip = 0, place = 1;
    while (i >= 0 && mc[i] >= '0' && mc[i] <= '9') {
        right_clip += (mc[i] - '0') * place;
        place *= 10;
        i--;
    }
    return right_clip;
}

bool is_mate_left_clipped(bam1_t* r) {
    if (is_mate_unmapped(r)) return false;

    char* mc_tag_str = get_mc(r);
    int i = 0;
    while (mc_tag_str[i] >= '0' && mc_tag_str[i] <= '9') i++;
    return mc_tag_str[i] == 'S';
}
bool is_mate_right_clipped(bam1_t* r) {
    if (is_mate_unmapped(r)) return false;

    char* mc_tag_str = get_mc(r);
    int i = strlen(mc_tag_str)-1;
    return mc_tag_str[i] == 'S';
}
bool is_mate_clipped(bam1_t* r) {
	return is_mate_left_clipped(r) || is_mate_right_clipped(r);
}

bool is_hidden_split_read(bam1_t* r, config_t config) {
    if (!is_samechr(r)) return false;
	if (is_left_clipped(r, config.min_clip_len) || is_right_clipped(r, config.min_clip_len)) return false;

    indel_summary_t indels = get_indel_summary(r);
    int mismatches = get_nm(r) - indels.dels - indels.inss;

    return mismatches + indels.ops >= config.min_diff_hsr || abs(indels.dels-indels.inss) >= config.min_sv_size;
}


hts_pos_t get_mate_endpos(bam1_t* r) {
    char* mc = get_mc(r);
    int i = 0, mclen = strlen(mc);

    int len = 0;
    hts_pos_t pos = r->core.mpos;
    while (i < mclen) {
        if (mc[i] >= '0' && mc[i] <= '9') {
            len = (len*10) + (mc[i]-'0');
        } else {
            if (mc[i] != 'I' && mc[i] != 'S') {
                pos += len;
            }
            len = 0;
        }
        i++;
    }
    return pos-1;
}

hts_pos_t get_mate_unclipped_start(bam1_t* r) {
    return r->core.mpos - get_mate_left_clip_size(r);
}
hts_pos_t get_mate_unclipped_end(bam1_t* r) {
    return get_mate_endpos(r) + get_mate_right_clip_size(r);
}

bool can_align_to(bam1_t* read, int max_is, int32_t contig_id, hts_pos_t coordinate) {
    if (is_mate_unmapped(read) || contig_id != read->core.mtid || max_is < 0 || read->core.l_qseq <= 0) return false;

    if (bam_is_mrev(read)) {
        hts_pos_t mate_endpos = get_mate_endpos(read);
        return mate_endpos - max_is <= coordinate && coordinate <= mate_endpos;
    } else {
        return read->core.mpos <= coordinate && coordinate <= read->core.mpos + max_is;
    }
}

bool has_no_indels(bam1_t* read) {
    uint32_t *cigar = bam_get_cigar(read);
    return read->core.n_cigar == 1 && bam_cigar_op(cigar[0]) == BAM_CMATCH && bam_cigar_oplen(cigar[0]) == read->core.l_qseq;
}

bool is_perfectly_aligned(bam1_t* read) {
    // All M
    if (has_no_indels(read)) {
        return get_nm(read) == 0; // Whether it is perfectly aligned based on NM
    }
    return false; // Not perfectly aligned or cannot confirm alignment without NM or MD
}

void rc(std::string& read) {
    int len = read.length();
    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
		char c = std::toupper(read[i]);
		if (c == 'A') read[i] = 'T';
		else if (c == 'C') read[i] = 'G';
		else if (c == 'G') read[i] = 'C';
		else if (c == 'T') read[i] = 'A';
		else c = 'N';
	}
}
void rc(char* read) {
    int len = strlen(read);

    for (int i = 0; i < len/2; i++) {
        std::swap(read[i], read[len-i-1]);
    }
    for (int i = 0; i < len; i++) {
    	char c = std::toupper(read[i]);
        if (c == 'A') read[i] = 'T';
        else if (c == 'C') read[i] = 'G';
        else if (c == 'G') read[i] = 'C';
        else if (c == 'T') read[i] = 'A';
        else c = 'N';
    }
}

alignas(16) 
static const char nucl2chr[16] = {
    0, 'A', 'C', 0, 'G', 0, 0, 0, 'T', 0, 0, 0, 0, 0, 0, 'N'
};
char get_base(const uint8_t* seq, int i) {
    return nucl2chr[bam_seqi(seq, i)];
}
std::string get_sequence(bam1_t* r, bool fastq_seq = false) {
    char seq[100000];
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (fastq_seq && bam_is_rev(r)) rc(seq);
    return std::string(seq);
}
std::string get_qual_ascii(bam1_t* r, bool fastq_seq = false) {
	uint8_t* qual = bam_get_qual(r);
	std::string qual_ascii(r->core.l_qseq, ' ');
	for (int i = 0; i < r->core.l_qseq; i++) {
		qual_ascii[i] = char(33 + qual[i]);
	}
	if (fastq_seq && bam_is_rev(r)) qual_ascii = std::string(qual_ascii.rbegin(), qual_ascii.rend());
	return qual_ascii;
}

std::string get_cigar_string(bam1_t* r) {
    std::string cigar_str;
    uint32_t* cigar = bam_get_cigar(r);
    for (uint32_t i = 0; i < r->core.n_cigar; i++) {
        cigar_str += std::to_string(bam_cigar_oplen(cigar[i])) + bam_cigar_opchr(cigar[i]);
    }
    return cigar_str;
}

double avg_qual(bam1_t* read) {
	double avg_qual = 0;
	uint8_t* qual = bam_get_qual(read);
	for (int i = 0; i < read->core.l_qseq; i++) {
		avg_qual += qual[i];
	}
	return avg_qual/read->core.l_qseq;
}
double avg_qual(std::string& qual_ascii, int offset = 33) {
	double tot = 0;
	for (char c : qual_ascii) tot += int(c)-offset;
	return tot/qual_ascii.length();
}

int get_aligned_portion_len(bam1_t* read) {
    return read->core.l_qseq - get_left_clip_size(read) - get_right_clip_size(read);
}

int64_t get_AS_tag(bam1_t* read) {
    uint8_t* aux_get = bam_aux_get(read, "AS");
    return bam_aux2i(aux_get);
}

void copy_sequence(bam1_t* r, char* seq, bool fastq_seq = false) { // assumes seq is long enough
    const uint8_t* bam_seq = bam_get_seq(r);
    for (int i = 0; i < r->core.l_qseq; i++) {
        seq[i] = get_base(bam_seq, i);
    }
    seq[r->core.l_qseq] = '\0';
    if (fastq_seq && bam_is_rev(r)) rc(seq);
}

samFile* open_writer(std::string filename, bam_hdr_t* header) {
    samFile* writer = sam_open(filename.c_str(), "wb");
    if (writer == NULL) {
        throw std::runtime_error("Unable to open " + filename);
    }
    if (sam_hdr_write(writer, header) != 0) {
        sam_close(writer);
        throw std::runtime_error("Could not write file " + filename);
    }
    return writer;
}

void write_and_index_file(std::vector<bam1_t*>& reads, std::string path, bam_hdr_t* header) {
    samFile* file = open_writer(path, header);
    if (file == NULL) {
        throw std::runtime_error("Unable to open " + path);
    }

    // write reads
    std::sort(reads.begin(), reads.end(), [](bam1_t *r1, bam1_t *r2) { return r1->core.pos < r2->core.pos; });
    for (bam1_t* r : reads) {
        int ok = sam_write1(file, header, r);
        if (ok < 0) throw std::runtime_error("Unable to write to " + path);
    }

    sam_close(file);

    file = sam_open(path.c_str(), "r");

    int code = sam_index_build(path.c_str(), 0);
    if (code != 0) {
        throw std::runtime_error("Cannot index " + path);
    }

    sam_close(file);
}


struct open_samFile_t {
    samFile* file = nullptr;
    bam_hdr_t* header = nullptr;
    hts_idx_t* idx = nullptr;

    open_samFile_t() = default;
    open_samFile_t(std::string fname_str, bool index_file = false) {
        std::unique_ptr<samFile, decltype(&hts_close)> file_ptr(sam_open(fname_str.c_str(), "r"), &hts_close);
        if (!file_ptr) {
            throw std::runtime_error("Could not open " + fname_str);
        }

        if (index_file) {
            int code = sam_index_build(fname_str.c_str(), 0);
            if (code != 0) {
                throw std::runtime_error("Cannot index " + fname_str);
            }
        }

        std::unique_ptr<hts_idx_t, decltype(&hts_idx_destroy)> idx(sam_index_load(file_ptr.get(), fname_str.c_str()), &hts_idx_destroy);
        if (!idx) {
            throw std::runtime_error("Unable to open index for " + fname_str);
        }

        std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> header_ptr(sam_hdr_read(file_ptr.get()), &bam_hdr_destroy);
        if (!header_ptr) {
            throw std::runtime_error("Unable to open header for " + fname_str);
        }

        file = file_ptr.release();
        header = header_ptr.release();
        this->idx = idx.release();
    }

    open_samFile_t(const open_samFile_t&) = delete;
    open_samFile_t& operator=(const open_samFile_t&) = delete;

    open_samFile_t(open_samFile_t&&) = delete;
    open_samFile_t& operator=(open_samFile_t&&) = delete;

    ~open_samFile_t() {
        if (idx) hts_idx_destroy(idx);
        if (header) bam_hdr_destroy(header);
        if (file) sam_close(file);
    }
};

struct bam_pool_t {
    std::mutex mtx;
    std::vector<open_samFile_t*> pool;
    std::string bam_fname, reference_fname;

    bam_pool_t(int size, std::string bam_fname, std::string reference_fname) : bam_fname(bam_fname), reference_fname(reference_fname) {
        for (int i = 0; i < size; i++) {
            open_samFile_t* o = new open_samFile_t(bam_fname.c_str());
            hts_set_fai_filename(o->file, fai_path(reference_fname.c_str()));
            pool.push_back(o);
        }
    }

    open_samFile_t* get_bam_reader(int i) {
        return pool[i];
    }

    ~bam_pool_t() {
        for (open_samFile_t* o : pool) {
            delete o;
        }
    }
};


#endif
