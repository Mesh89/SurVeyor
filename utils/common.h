#ifndef SVCOMPARE_COMMON_H
#define SVCOMPARE_COMMON_H

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

#include "../src/vcf_utils.h"
#include "../src/utils.h"

bool ends_with(const char* str, const char* suffix) {
	int l_string = strlen(str);
	int l_suffix = strlen(suffix);
    if (l_suffix > l_string) return false;
	return strcmp(str+(l_string-l_suffix), suffix) == 0;
}

int get_optional_format_int32(bcf_hdr_t* hdr, bcf1_t* line, const char* tag, int fallback = 0) {
    int32_t* vals = NULL;
    int nvals = 0;
    int value = fallback;
    int ret = bcf_get_format_int32(hdr, line, tag, &vals, &nvals);
    if (ret > 0 && nvals > 0 &&
        vals[0] != bcf_int32_missing && vals[0] != bcf_int32_vector_end) {
        value = vals[0];
    }
    free(vals);
    return value;
}

void load_compare_tiebreak_fields(bcf_hdr_t* hdr, bcf1_t* line, sv_t* sv) {
    sv->sample_info.alt_bp1.reads_info.reads = get_optional_format_int32(hdr, line, "AR1");
    sv->sample_info.alt_bp1.reads_info.consistent_fwd = get_optional_format_int32(hdr, line, "AR1CF");
    sv->sample_info.alt_bp1.reads_info.consistent_rev = get_optional_format_int32(hdr, line, "AR1CR");

    sv->sample_info.alt_bp2.reads_info.reads = get_optional_format_int32(hdr, line, "AR2");
    sv->sample_info.alt_bp2.reads_info.consistent_fwd = get_optional_format_int32(hdr, line, "AR2CF");
    sv->sample_info.alt_bp2.reads_info.consistent_rev = get_optional_format_int32(hdr, line, "AR2CR");

    sv->sample_info.ext_alt_consensus1_to_alt_score = get_optional_format_int32(hdr, line, "XAAS");
    sv->sample_info.ext_alt_consensus2_to_alt_score = get_optional_format_int32(hdr, line, "XAAS2");
}

std::vector<std::shared_ptr<sv_t>> read_sv_list(const char* filename) {
	std::vector<std::shared_ptr<sv_t>> svs;
    htsFile* file = bcf_open(filename, "r");
    if (file == NULL) {
        throw std::runtime_error("Error: could not open file " + std::string(filename));
    }
    bcf_hdr_t* hdr = bcf_hdr_read(file);
    if (hdr == NULL) {
        bcf_close(file);
        throw std::runtime_error("Error: could not read header from file " + std::string(filename));
    }
    bcf1_t* line = bcf_init();
    while (bcf_read(file, hdr, line) == 0) {
        std::shared_ptr<sv_t> sv = bcf_to_sv(hdr, line);
        if (sv != nullptr) {
            load_compare_tiebreak_fields(hdr, line, sv.get());
            sv->vcf_entry = bcf_dup(line);
            svs.push_back(sv);
        }
    }
    bcf_destroy(line);
    bcf_hdr_destroy(hdr);
    bcf_close(file);
    return svs;
}

struct repeat_t {
    std::string chr;
    int start, end;

    repeat_t() : start(0), end(0) {}

    bool contains(sv_t* sv) {
        return sv->chr == chr && start-10 <= sv->start && sv->end <= end+10;
    }
    bool intersects(sv_t* sv) {
        return sv->chr == chr && sv->start <= end && sv->end >= start;
    }
};
bool operator == (const repeat_t& r1, const repeat_t& r2) {
    return r1.chr == r2.chr && r1.start == r2.start && r1.end == r2.end;
}

double overlap_bps(sv_t* sv1, sv_t* sv2) {
    overlap(sv1->start, sv1->end, sv2->start, sv2->end);
	return std::max(hts_pos_t(0), std::min(sv1->end, sv2->end)-std::max(sv1->start, sv2->start));
}

double overlap(sv_t* sv1, sv_t* sv2) {
	if (sv1->svtype() == "INS" && sv2->svtype() == "INS") return 1.0;
	if (sv1->end == sv1->start || sv2->end == sv2->start) return 1.0;
    return overlap_bps(sv1, sv2)/double(std::min(sv1->end-sv1->start, sv2->end-sv2->start));
}

int distance(sv_t* sv1, sv_t* sv2) {
    if (sv1->chr != sv2->chr) return INT32_MAX;
	if ((sv1->svtype() == "DUP" && sv2->svtype() == "INS") || (sv1->svtype() == "INS" && sv2->svtype() == "DUP")) {
		return std::min(abs(sv1->start-sv2->start), abs(sv1->end-sv2->end));
	}
    return std::max(abs(sv1->start-sv2->start), abs(sv1->end-sv2->end));
}

#endif //SVCOMPARE_COMMON_H
