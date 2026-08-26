#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include "../libs/cptl_stl.h"

const double NAN_VALUE = std::numeric_limits<double>::quiet_NaN();
const int GENOTYPE_CONSENSUS_EXTENSION = 500;

using stats_t = std::unordered_map<std::string, std::unordered_map<std::string, int>>;
using features_t = std::unordered_map<std::string, double>;

struct alt_read_metrics_t {
    double ar1, ar2, ar1c, ar2c, ar1chq, ar2chq, ar1e, ar2e;
    bool hp_genotyped;
    std::string chrom;
    int64_t start, stop;
    double rr1, rr2, rr1c, rr2c, rr1e, rr2e;
    std::vector<std::string> oar1_vids, oar2_vids;
    double xaas_xars_diff_to_len;
    bool has_extension_evidence;
    double aas_ars_diff_to_len;
    bool has_assembly_evidence;

    double ar(int bp_idx) const { return bp_idx == 0 ? ar1 : ar2; }
    double arc(int bp_idx) const { return bp_idx == 0 ? ar1c : ar2c; }
    double archq(int bp_idx) const { return bp_idx == 0 ? ar1chq : ar2chq; }
    double are(int bp_idx) const { return bp_idx == 0 ? ar1e : ar2e; }
    const std::vector<std::string>& oar_vids(int bp_idx) const { return bp_idx == 0 ? oar1_vids : oar2_vids; }
};

struct selected_alt_read_metrics_t {
    double oar, oarc, oarchq, oare, hp_genotyped, min_archq, min_ar_over_nar, min_arc_over_narc, min_are_over_nare, xaas_xars_diff_to_len, aas_ars_diff_to_len;
};

struct model_data_t {
    std::vector<std::string> feature_names;
    std::vector<uint64_t> variant_ids;
    std::vector<uint64_t> record_keys;
    std::vector<double> values;
};

std::vector<std::string> split(const std::string& str, char delim) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        size_t end = str.find(delim, start);
        fields.push_back(str.substr(start, end == std::string::npos ? end : end-start));
        if (end == std::string::npos) break;
        start = end+1;
    }
    return fields;
}

stats_t load_stats(const std::string& fname) {
    std::ifstream in(fname);
    if (!in) throw std::runtime_error("Failed to open stats file " + fname + ".");
    stats_t stats;
    std::string stat, key;
    int value;
    while (in >> stat >> key >> value) stats[stat][key] = value;
    return stats;
}

int get_stat(const stats_t& stats, const std::string& stat, const std::string& chrom) {
    auto stat_it = stats.find(stat);
    if (stat_it == stats.end()) throw std::runtime_error("Missing statistic " + stat + ".");
    auto chrom_it = stat_it->second.find(chrom);
    if (chrom_it != stat_it->second.end()) return chrom_it->second;
    auto default_it = stat_it->second.find(".");
    if (default_it == stat_it->second.end()) throw std::runtime_error("Missing statistic " + stat + " for " + chrom + ".");
    return default_it->second;
}

bool has_format(bcf_hdr_t* hdr, bcf1_t* record, const char* key) {
    return bcf_get_fmt(hdr, record, key) != NULL;
}

bool has_info(bcf_hdr_t* hdr, bcf1_t* record, const char* key) {
    return bcf_get_info(hdr, record, key) != NULL;
}

std::vector<double> get_format_numbers(bcf_hdr_t* hdr, bcf1_t* record, const char* key, const std::vector<double>& defaults) {
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, key);
    if (id < 0 || !has_format(hdr, record, key)) return defaults;
    int type = bcf_hdr_id2type(hdr, BCF_HL_FMT, id);
    std::vector<double> values;
    if (type == BCF_HT_INT) {
        int32_t* data = NULL;
        int ndata = 0;
        int n = bcf_get_format_int32(hdr, record, key, &data, &ndata);
        for (int i = 0; i < n && data[i] != bcf_int32_vector_end; i++) values.push_back(data[i] == bcf_int32_missing ? (i < defaults.size() ? defaults[i] : defaults.back()) : data[i]);
        free(data);
    } else if (type == BCF_HT_REAL) {
        float* data = NULL;
        int ndata = 0;
        int n = bcf_get_format_float(hdr, record, key, &data, &ndata);
        for (int i = 0; i < n && !bcf_float_is_vector_end(data[i]); i++) values.push_back(bcf_float_is_missing(data[i]) ? (i < defaults.size() ? defaults[i] : defaults.back()) : data[i]);
        free(data);
    }
    return values.empty() ? defaults : values;
}

double get_format_number(bcf_hdr_t* hdr, bcf1_t* record, const char* key, double default_value, double norm_factor = 1.0) {
    return get_format_numbers(hdr, record, key, {default_value})[0]/norm_factor;
}

std::vector<double> get_format_numbers(bcf_hdr_t* hdr, bcf1_t* record, const char* key, std::initializer_list<double> defaults) {
    std::vector<double> default_values(defaults);
    std::vector<double> values = get_format_numbers(hdr, record, key, default_values);
    if (values.size() < default_values.size()) values.resize(default_values.size(), default_values.back());
    return values;
}

std::vector<double> get_info_numbers(bcf_hdr_t* hdr, bcf1_t* record, const char* key) {
    int id = bcf_hdr_id2int(hdr, BCF_DT_ID, key);
    if (id < 0 || !has_info(hdr, record, key)) return {};
    int type = bcf_hdr_id2type(hdr, BCF_HL_INFO, id);
    std::vector<double> values;
    if (type == BCF_HT_INT) {
        int32_t* data = NULL;
        int ndata = 0;
        int n = bcf_get_info_int32(hdr, record, key, &data, &ndata);
        for (int i = 0; i < n && data[i] != bcf_int32_vector_end; i++) values.push_back(data[i] == bcf_int32_missing ? NAN_VALUE : data[i]);
        free(data);
    } else if (type == BCF_HT_REAL) {
        float* data = NULL;
        int ndata = 0;
        int n = bcf_get_info_float(hdr, record, key, &data, &ndata);
        for (int i = 0; i < n && !bcf_float_is_vector_end(data[i]); i++) values.push_back(bcf_float_is_missing(data[i]) ? NAN_VALUE : data[i]);
        free(data);
    }
    return values;
}

std::string get_info_string(bcf_hdr_t* hdr, bcf1_t* record, const char* key, const std::string& default_value = "") {
    char* data = NULL;
    int ndata = 0;
    int n = bcf_get_info_string(hdr, record, key, &data, &ndata);
    std::string value = n > 0 && data != NULL ? data : default_value;
    free(data);
    return value;
}

std::string get_format_string(bcf_hdr_t* hdr, bcf1_t* record, const char* key, const std::string& default_value = "") {
    char** data = NULL;
    int ndata = 0;
    int n = bcf_get_format_string(hdr, record, key, &data, &ndata);
    std::string value = n > 0 && data != NULL && data[0] != NULL ? data[0] : default_value;
    if (data != NULL) free(data[0]);
    free(data);
    return value;
}

std::vector<std::string> get_info_strings(bcf_hdr_t* hdr, bcf1_t* record, const char* key) {
    if (!has_info(hdr, record, key)) return {};
    return split(get_info_string(hdr, record, key), ',');
}

std::string normalize_sv_id(const std::string& id) {
    return id.size() >= 4 && id.substr(id.size()-4) == "_DUP" ? id.substr(0, id.size()-4) : id;
}

std::vector<std::string> get_oar_vids(bcf_hdr_t* hdr, bcf1_t* record, const char* key) {
    std::vector<std::string> vids;
    if (!has_format(hdr, record, key)) return vids;
    for (const std::string& value : split(get_format_string(hdr, record, key), ',')) if (!value.empty() && value != ".") vids.push_back(normalize_sv_id(value));
    return vids;
}

std::string get_svtype(bcf_hdr_t* hdr, bcf1_t* record) {
    return get_info_string(hdr, record, "SVTYPE");
}

int64_t record_pos(bcf1_t* record) { return record->pos+1; }
int64_t record_stop(bcf1_t* record) { return record->pos+record->rlen; }

std::string get_svinsseq(bcf_hdr_t* hdr, bcf1_t* record) {
    bcf_unpack(record, BCF_UN_STR);
    std::string alt = record->n_allele > 1 ? record->d.allele[1] : "";
    if (alt.find('<') == std::string::npos) return alt.empty() ? "" : alt.substr(1);
    if (has_info(hdr, record, "SVINSSEQ")) {
        std::string svinsseq = get_info_string(hdr, record, "SVINSSEQ");
        size_t comma = svinsseq.find(',');
        return comma == std::string::npos ? svinsseq : svinsseq.substr(0, comma);
    }
    if (has_info(hdr, record, "LEFT_SVINSSEQ") || has_info(hdr, record, "RIGHT_SVINSSEQ")) {
        std::string left = get_info_string(hdr, record, "LEFT_SVINSSEQ"), right = get_info_string(hdr, record, "RIGHT_SVINSSEQ");
        size_t left_comma = left.find(','), right_comma = right.find(',');
        if (left_comma != std::string::npos) left.resize(left_comma);
        if (right_comma != std::string::npos) right.resize(right_comma);
        return left+"-"+right;
    }
    return "";
}

int64_t get_svlen(bcf_hdr_t* hdr, bcf1_t* record) {
    std::string svtype = get_svtype(hdr, record);
    std::string svinsseq = get_svinsseq(hdr, record);
    int64_t svlen;
    if (svtype == "INS" || svtype == "INS_TO_DUP" || svtype == "DEL") svlen = svinsseq.size()-(record_stop(record)-record_pos(record));
    else if (svtype == "DUP") svlen = record_stop(record)-record_pos(record)+svinsseq.size();
    else throw std::runtime_error("Unexpected SVTYPE " + svtype + ".");
    for (const std::string& indel : get_info_strings(hdr, record, "AUX_INDELS")) {
        std::vector<std::string> fields = split(indel, ':');
        if (fields.size() < 3) throw std::runtime_error("Malformed AUX_INDELS value " + indel + ".");
        svlen -= std::stoll(fields[1])-std::stoll(fields[0]), svlen += fields[2].size();
    }
    return svlen;
}

int64_t get_edit_distance(bcf_hdr_t* hdr, bcf1_t* record, int64_t svinslen) {
    int64_t edit_distance = record_stop(record)-record_pos(record)+svinslen;
    edit_distance += get_info_strings(hdr, record, "AUX_SNPS").size();
    for (const std::string& indel : get_info_strings(hdr, record, "AUX_INDELS")) {
        std::vector<std::string> fields = split(indel, ':');
        if (fields.size() < 3) throw std::runtime_error("Malformed AUX_INDELS value " + indel + ".");
        edit_distance += std::stoll(fields[1])-std::stoll(fields[0])+fields[2].size();
    }
    return edit_distance;
}

bool gt_as_homopolymer(bcf_hdr_t* hdr, bcf1_t* record) { return has_info(hdr, record, "HP_GENOTYPED"); }
bool has_assembly_evidence(bcf_hdr_t* hdr, bcf1_t* record) { return has_format(hdr, record, "AL") || has_format(hdr, record, "AL2"); }
bool has_extension_evidence(bcf_hdr_t* hdr, bcf1_t* record) { return has_format(hdr, record, "XAL") || has_format(hdr, record, "XAL2"); }

std::string get_model_name(bcf_hdr_t* hdr, bcf1_t* record, int max_is, int read_len) {
    if (gt_as_homopolymer(hdr, record)) return "HP";
    std::string svtype = get_svtype(hdr, record);
    if (svtype == "DUP" && has_info(hdr, record, "INS_TO_DUP")) {
        svtype = "INS_TO_DUP";
        if (get_svlen(hdr, record) > read_len-30) svtype += "_LARGE";
    }
    if (svtype == "DEL") {
        if (std::abs(get_svlen(hdr, record)) >= max_is) svtype += "_LARGE";
        if (!has_format(hdr, record, "XAL")) svtype += "_NOEXL";
    } else if (svtype == "DUP" && get_svlen(hdr, record) > read_len-30) {
        svtype += "_LARGE";
        if (!has_format(hdr, record, "XAL")) svtype += "_NOEXL";
    }
    return svtype;
}

double piecewise_normalise(double value, double minv, double maxv, double quantisation = 0.025) {
    bool neg = value < 0;
    value = std::abs(value);
    double ret_val;
    if (value <= minv) ret_val = value/std::max(1.0, minv)*0.25;
    else if (value > maxv) ret_val = value/std::max(1.0, maxv)*0.25+0.75;
    else ret_val = 0.25+(value-minv)/(maxv-minv)*0.75;
    if (neg) ret_val = -ret_val;
    return std::nearbyint(ret_val/quantisation)*quantisation;
}

std::pair<double, double> normalise_mate_coverage_spans(const std::vector<double>& spans, double reads, double read_len, double quantisation = 0.025) {
    double left = spans[0]/std::max(1.0, reads*read_len), right = spans[1]/std::max(1.0, reads*read_len);
    return {std::nearbyint(left/quantisation)*quantisation, std::nearbyint(right/quantisation)*quantisation};
}

double calculate_z_score(double mean1, double stddev1, double n1, double mean2, double stddev2, double n2) {
    if (std::isnan(mean1) || std::isnan(mean2) || n1 == 0 || n2 == 0) return NAN_VALUE;
    double std_error = std::sqrt(stddev1*stddev1/n1+stddev2*stddev2/n2);
    if (std_error == 0) std_error = 1;
    return (mean1-mean2)/std_error;
}

double consensus_alt_ref_score_diff_to_len(bcf_hdr_t* hdr, bcf1_t* record, const std::string& prefix) {
    double aas1 = get_format_number(hdr, record, (prefix+"AAS").c_str(), 0), aas2 = get_format_number(hdr, record, (prefix+"AAS2").c_str(), 0);
    double ars1 = get_format_number(hdr, record, (prefix+"ARS").c_str(), 0), ars2 = get_format_number(hdr, record, (prefix+"ARS2").c_str(), 0);
    return (aas1-ars1+aas2-ars2)/std::max<int64_t>(1, get_edit_distance(hdr, record, get_svinsseq(hdr, record).size()));
}

void add_consensus_alignment_features(features_t& features, bcf_hdr_t* hdr, bcf1_t* record, const std::string& prefix, int read_len, int64_t edit_distance) {
    double al1 = get_format_number(hdr, record, (prefix+"AL").c_str(), NAN_VALUE), al2 = get_format_number(hdr, record, (prefix+"AL2").c_str(), NAN_VALUE);
    double normalised_al1, normalised_al2;
    if (prefix.empty()) {
        std::vector<double> factors = get_format_numbers(hdr, record, "MFAL", {(double)read_len, (double)read_len});
        normalised_al1 = std::isfinite(al1) ? al1/factors[0] : NAN_VALUE;
        normalised_al2 = std::isfinite(al2) ? al2/factors[1] : NAN_VALUE;
    } else {
        double unextended_al1 = get_format_number(hdr, record, "AL", NAN_VALUE), unextended_al2 = get_format_number(hdr, record, "AL2", NAN_VALUE);
        normalised_al1 = std::isfinite(al1) && std::isfinite(unextended_al1) ? (al1-unextended_al1)/(2.0*GENOTYPE_CONSENSUS_EXTENSION) : NAN_VALUE;
        normalised_al2 = std::isfinite(al2) && std::isfinite(unextended_al2) ? (al2-unextended_al2)/(2.0*GENOTYPE_CONSENSUS_EXTENSION) : NAN_VALUE;
    }
    features[prefix+"AL1"] = normalised_al1;
    features[prefix+"AL2"] = normalised_al2;
    features[prefix+"AL"] = std::isfinite(normalised_al1) || std::isfinite(normalised_al2) ? (std::isfinite(normalised_al1) ? normalised_al1 : 0)+(std::isfinite(normalised_al2) ? normalised_al2 : 0) : NAN_VALUE;
    double aas1 = get_format_number(hdr, record, (prefix+"AAS").c_str(), 0), aas2 = get_format_number(hdr, record, (prefix+"AAS2").c_str(), 0);
    double ars1 = get_format_number(hdr, record, (prefix+"ARS").c_str(), 0), ars2 = get_format_number(hdr, record, (prefix+"ARS2").c_str(), 0);
    features[prefix+"AAS_"+prefix+"ARS_DIFF_TO_LEN"] = (aas1-ars1+aas2-ars2)/std::max<int64_t>(1, edit_distance);
    std::vector<double> ass1 = get_format_numbers(hdr, record, (prefix+"ASS").c_str(), {NAN_VALUE, NAN_VALUE});
    std::vector<double> ass2 = get_format_numbers(hdr, record, (prefix+"ASS2").c_str(), {NAN_VALUE, NAN_VALUE});
    features[prefix+"ASS1_LEFT_RATIO"] = ass1[0]/std::max(1.0, al1); features[prefix+"ASS1_RIGHT_RATIO"] = ass1[1]/std::max(1.0, al1);
    features[prefix+"ASS2_LEFT_RATIO"] = ass2[0]/std::max(1.0, al2); features[prefix+"ASS2_RIGHT_RATIO"] = ass2[1]/std::max(1.0, al2);
    std::vector<double> assc1 = get_format_numbers(hdr, record, (prefix+"ASSC").c_str(), {NAN_VALUE, NAN_VALUE});
    std::vector<double> assc2 = get_format_numbers(hdr, record, (prefix+"ASSC2").c_str(), {NAN_VALUE, NAN_VALUE});
    std::vector<double> asscia1 = get_format_numbers(hdr, record, (prefix+"ASSCIA").c_str(), {NAN_VALUE, NAN_VALUE});
    std::vector<double> asscia2 = get_format_numbers(hdr, record, (prefix+"ASSC2IA").c_str(), {NAN_VALUE, NAN_VALUE});
    features[prefix+"ASSC1_IA_RATIO"] = (assc1[0]+assc1[1])/std::max(1.0, asscia1[0]+asscia1[1]);
    features[prefix+"ASSC2_IA_RATIO"] = (assc2[0]+assc2[1])/std::max(1.0, asscia2[0]+asscia2[1]);
    features[prefix+"ASSC1_IA_DIFF"] = (asscia1[0]+asscia1[1]-assc1[0]-assc1[1])/std::max(1.0, ass1[0]+ass1[1]);
    features[prefix+"ASSC2_IA_DIFF"] = (asscia2[0]+asscia2[1]-assc2[0]-assc2[1])/std::max(1.0, ass2[0]+ass2[1]);
}

uint64_t fnv1a(const std::string& value) {
    uint64_t hash = 14695981039346656037ULL;
    for (unsigned char c : value) hash = (hash^c)*1099511628211ULL;
    return hash;
}

uint64_t generate_id(bcf_hdr_t* hdr, bcf1_t* record, const std::string& model_name) {
    std::ostringstream key;
    key << bcf_seqname(hdr, record) << ':' << record_pos(record) << '-' << record_stop(record) << ':' << get_svtype(hdr, record) << ':' << get_svlen(hdr, record) << ':' << get_svinsseq(hdr, record) << ':' << get_info_string(hdr, record, "AUX_SNPS") << ':' << get_info_string(hdr, record, "AUX_INDELS") << ':' << (record->d.id == NULL ? "." : record->d.id) << ':' << model_name;
    return fnv1a(key.str());
}

selected_alt_read_metrics_t select_alt_read_metrics_for_oar_vids(bcf_hdr_t* hdr, bcf1_t* record, const char* key, int bp_idx, const std::unordered_map<std::string, std::vector<alt_read_metrics_t>>& alt_reads_by_vid) {
    typedef std::pair<const alt_read_metrics_t*, std::string> candidate_t;
    std::vector<candidate_t> candidates;
    for (const std::string& vid : get_oar_vids(hdr, record, key)) {
        auto it = alt_reads_by_vid.find(vid);
        if (it != alt_reads_by_vid.end()) for (const alt_read_metrics_t& values : it->second) candidates.push_back({&values, vid});
    }
    if (candidates.empty()) return {0, 0, 0, 0, 0, 0, 0, 0, 0, NAN_VALUE, NAN_VALUE};
    auto better = [bp_idx](const candidate_t& a, const candidate_t& b) {
        if (a.first->ar(bp_idx) != b.first->ar(bp_idx)) return a.first->ar(bp_idx) > b.first->ar(bp_idx);
        if (a.first->arc(bp_idx) != b.first->arc(bp_idx)) return a.first->arc(bp_idx) > b.first->arc(bp_idx);
        if (a.first->are(bp_idx) != b.first->are(bp_idx)) return a.first->are(bp_idx) > b.first->are(bp_idx);
        return a.second < b.second;
    };
    const alt_read_metrics_t* selected = std::min_element(candidates.begin(), candidates.end(), better)->first;
    const alt_read_metrics_t* selected_oars[2] = {NULL, NULL};
    for (int oar_bp_idx = 0; oar_bp_idx < 2; oar_bp_idx++) {
        std::vector<candidate_t> oar_candidates;
        for (const std::string& vid : selected->oar_vids(oar_bp_idx)) {
            auto it = alt_reads_by_vid.find(vid);
            if (it != alt_reads_by_vid.end()) for (const alt_read_metrics_t& values : it->second) oar_candidates.push_back({&values, vid});
        }
        if (!oar_candidates.empty()) selected_oars[oar_bp_idx] = std::min_element(oar_candidates.begin(), oar_candidates.end(), [oar_bp_idx](const candidate_t& a, const candidate_t& b) {
            if (a.first->ar(oar_bp_idx) != b.first->ar(oar_bp_idx)) return a.first->ar(oar_bp_idx) > b.first->ar(oar_bp_idx);
            if (a.first->arc(oar_bp_idx) != b.first->arc(oar_bp_idx)) return a.first->arc(oar_bp_idx) > b.first->arc(oar_bp_idx);
            if (a.first->are(oar_bp_idx) != b.first->are(oar_bp_idx)) return a.first->are(oar_bp_idx) > b.first->are(oar_bp_idx);
            return a.second < b.second;
        })->first;
    }
    double min_ar_over_nar = std::min(selected->ar1/std::max(1.0, selected->ar1+(selected_oars[0] ? selected_oars[0]->ar1 : 0)+selected->rr1), selected->ar2/std::max(1.0, selected->ar2+(selected_oars[1] ? selected_oars[1]->ar2 : 0)+selected->rr2));
    double min_arc_over_narc = std::min(selected->ar1c/std::max(1.0, selected->ar1c+(selected_oars[0] ? selected_oars[0]->ar1c : 0)+selected->rr1c), selected->ar2c/std::max(1.0, selected->ar2c+(selected_oars[1] ? selected_oars[1]->ar2c : 0)+selected->rr2c));
    double min_are_over_nare = std::min(selected->ar1e/std::max(1.0, selected->ar1e+(selected_oars[0] ? selected_oars[0]->ar1e : 0)+selected->rr1e), selected->ar2e/std::max(1.0, selected->ar2e+(selected_oars[1] ? selected_oars[1]->ar2e : 0)+selected->rr2e));
    return {selected->ar(bp_idx), selected->arc(bp_idx), selected->archq(bp_idx), selected->are(bp_idx), (double)selected->hp_genotyped, std::min(selected->ar1chq, selected->ar2chq), min_ar_over_nar, min_arc_over_narc, min_are_over_nare, selected->has_extension_evidence ? selected->xaas_xars_diff_to_len : NAN_VALUE, selected->has_assembly_evidence ? selected->aas_ars_diff_to_len : NAN_VALUE};
}

void add_base_count_features(features_t& features, const std::vector<double>& counts, const std::string& prefix, const std::string& max_name) {
    double total = 0;
    for (double count : counts) total += count;
    double denom = std::max(1.0, total);
    features[prefix+"A_RATIO"] = counts[0]/denom;
    features[prefix+"C_RATIO"] = counts[1]/denom;
    features[prefix+"G_RATIO"] = counts[2]/denom;
    features[prefix+"T_RATIO"] = counts[3]/denom;
    features[max_name] = std::max(std::max(features[prefix+"A_RATIO"], features[prefix+"C_RATIO"]), std::max(features[prefix+"G_RATIO"], features[prefix+"T_RATIO"]));
}

void add_read_pair_features(features_t& features, bcf_hdr_t* hdr, bcf1_t* record, const std::string& prefix, double min_disc_pairs, double max_disc_pairs, double& reads, double& nma1, double& nma2, double& nms1, double& nms2) {
    reads = get_format_number(hdr, record, prefix.c_str(), 0);
    std::vector<double> hq = get_format_numbers(hdr, record, (prefix+"HQ").c_str(), {0, 0});
    std::vector<double> nma = get_format_numbers(hdr, record, (prefix+"NMA").c_str(), {NAN_VALUE, NAN_VALUE});
    std::vector<double> nms = get_format_numbers(hdr, record, (prefix+"NMS").c_str(), {NAN_VALUE, NAN_VALUE});
    nma1 = nma[0]; nma2 = nma[1]; nms1 = nms[0]; nms2 = nms[1];
    features[prefix] = piecewise_normalise(reads, min_disc_pairs, max_disc_pairs);
    features[prefix+"HQ_1"] = piecewise_normalise(hq[0], min_disc_pairs, max_disc_pairs);
    features[prefix+"HQ_2"] = piecewise_normalise(hq[1], min_disc_pairs, max_disc_pairs);
    features[prefix+"HQ_1_RATIO"] = hq[0]/std::max(1.0, reads);
    features[prefix+"HQ_2_RATIO"] = hq[1]/std::max(1.0, reads);
    std::vector<double> min_mq = get_format_numbers(hdr, record, (prefix+"mQ").c_str(), {NAN_VALUE, NAN_VALUE});
    std::vector<double> max_mq = get_format_numbers(hdr, record, (prefix+"MQ").c_str(), {NAN_VALUE, NAN_VALUE});
    features[prefix+"mQ_1"] = min_mq[0]; features[prefix+"mQ_2"] = min_mq[1];
    features[prefix+"MQ_1"] = max_mq[0]; features[prefix+"MQ_2"] = max_mq[1];
}

std::vector<double> record_to_features(bcf_hdr_t* hdr, bcf1_t* record, const stats_t& stats, const std::vector<std::string>& feature_names, const std::unordered_map<std::string, std::vector<alt_read_metrics_t>>& alt_reads_by_vid) {
    std::string chrom = bcf_seqname(hdr, record);
    double min_depth = get_stat(stats, "min_depth", chrom), max_depth = get_stat(stats, "max_depth", chrom);
    int max_is = get_stat(stats, "max_is", "."), read_len = get_stat(stats, "read_len", ".");
    double min_pairs_crossing_point = stats.at("min_pairs_crossing_gap").at("0"), max_pairs_crossing_point = stats.at("max_pairs_crossing_gap").at("0");
    std::string model_name = get_model_name(hdr, record, max_is, read_len), svtype = get_svtype(hdr, record), source = get_info_string(hdr, record, "SOURCE");
    features_t features;
    features["START_STOP_DIST"] = record_stop(record)-record_pos(record);
    int64_t signed_svlen = get_svlen(hdr, record);
    double svlen = std::abs(signed_svlen);
    features["SVLEN"] = std::log1p(svlen);
    std::string svinsseq = get_svinsseq(hdr, record);
    int64_t svinslen = svinsseq.size(), edit_distance = get_edit_distance(hdr, record, svinslen);
    features["SVINSLEN"] = svinslen; features["EDIT_DISTANCE"] = edit_distance;
    features["INS_SEQ_COV_PREFIX_LEN"] = 1; features["INS_SEQ_COV_SUFFIX_LEN"] = 1;
    size_t dash = svinsseq.find('-');
    if (dash != std::string::npos && !svinsseq.empty()) {
        features["INS_SEQ_COV_PREFIX_LEN"] = (double)dash/svinsseq.size();
        features["INS_SEQ_COV_SUFFIX_LEN"] = (double)(svinsseq.size()-dash-1)/svinsseq.size();
    }
    std::vector<double> earf = get_format_numbers(hdr, record, "EARF", {NAN_VALUE, NAN_VALUE});
    double exp_alt_reads_freq1 = earf[0], exp_alt_reads_freq2 = earf[1];
    features["EXP_ALT_READS_FREQ1"] = std::nearbyint(exp_alt_reads_freq1*100)/100;
    features["EXP_ALT_READS_FREQ2"] = std::nearbyint(exp_alt_reads_freq2*100)/100;
    if (model_name == "HP") {
        std::vector<double> hp_range = get_info_numbers(hdr, record, "HP_REF_RANGE");
        features["HP_REF_LEN"] = hp_range[1]-hp_range[0];
        features["HP_ALT_LEN"] = features["HP_REF_LEN"]+signed_svlen;
    } else features["HP_REF_LEN"] = features["HP_ALT_LEN"] = NAN_VALUE;

    add_base_count_features(features, get_info_numbers(hdr, record, "LEFT_ANCHOR_BASE_COUNT"), "LEFT_ANCHOR_", "MAX_LEFT_ANCHOR_BASE_RATIO");
    add_base_count_features(features, get_info_numbers(hdr, record, "LEFT_FLANKING_BASE_COUNT_50"), "LEFT_FLANKING_", "MAX_LEFT_FLANKING_BASE_RATIO_50");
    features["LEFT_FLANKING_A_RATIO_50"] = features["LEFT_FLANKING_A_RATIO"]; features["LEFT_FLANKING_C_RATIO_50"] = features["LEFT_FLANKING_C_RATIO"]; features["LEFT_FLANKING_G_RATIO_50"] = features["LEFT_FLANKING_G_RATIO"]; features["LEFT_FLANKING_T_RATIO_50"] = features["LEFT_FLANKING_T_RATIO"];
    add_base_count_features(features, get_info_numbers(hdr, record, "LEFT_FLANKING_BASE_COUNT_100"), "LEFT_FLANKING_100_", "MAX_LEFT_FLANKING_BASE_RATIO_100");
    features["LEFT_FLANKING_A_RATIO_100"] = features["LEFT_FLANKING_100_A_RATIO"]; features["LEFT_FLANKING_C_RATIO_100"] = features["LEFT_FLANKING_100_C_RATIO"]; features["LEFT_FLANKING_G_RATIO_100"] = features["LEFT_FLANKING_100_G_RATIO"]; features["LEFT_FLANKING_T_RATIO_100"] = features["LEFT_FLANKING_100_T_RATIO"];
    add_base_count_features(features, get_info_numbers(hdr, record, "LEFT_FLANKING_BASE_COUNT_500"), "LEFT_FLANKING_500_", "MAX_LEFT_FLANKING_BASE_RATIO_500");
    features["LEFT_FLANKING_A_RATIO_500"] = features["LEFT_FLANKING_500_A_RATIO"]; features["LEFT_FLANKING_C_RATIO_500"] = features["LEFT_FLANKING_500_C_RATIO"]; features["LEFT_FLANKING_G_RATIO_500"] = features["LEFT_FLANKING_500_G_RATIO"]; features["LEFT_FLANKING_T_RATIO_500"] = features["LEFT_FLANKING_500_T_RATIO"];
    add_base_count_features(features, get_info_numbers(hdr, record, "RIGHT_ANCHOR_BASE_COUNT"), "RIGHT_ANCHOR_", "MAX_RIGHT_ANCHOR_BASE_RATIO");
    add_base_count_features(features, get_info_numbers(hdr, record, "RIGHT_FLANKING_BASE_COUNT_50"), "RIGHT_FLANKING_", "MAX_RIGHT_FLANKING_BASE_RATIO_50");
    features["RIGHT_FLANKING_A_RATIO_50"] = features["RIGHT_FLANKING_A_RATIO"]; features["RIGHT_FLANKING_C_RATIO_50"] = features["RIGHT_FLANKING_C_RATIO"]; features["RIGHT_FLANKING_G_RATIO_50"] = features["RIGHT_FLANKING_G_RATIO"]; features["RIGHT_FLANKING_T_RATIO_50"] = features["RIGHT_FLANKING_T_RATIO"];
    add_base_count_features(features, get_info_numbers(hdr, record, "RIGHT_FLANKING_BASE_COUNT_100"), "RIGHT_FLANKING_100_", "MAX_RIGHT_FLANKING_BASE_RATIO_100");
    features["RIGHT_FLANKING_A_RATIO_100"] = features["RIGHT_FLANKING_100_A_RATIO"]; features["RIGHT_FLANKING_C_RATIO_100"] = features["RIGHT_FLANKING_100_C_RATIO"]; features["RIGHT_FLANKING_G_RATIO_100"] = features["RIGHT_FLANKING_100_G_RATIO"]; features["RIGHT_FLANKING_T_RATIO_100"] = features["RIGHT_FLANKING_100_T_RATIO"];
    add_base_count_features(features, get_info_numbers(hdr, record, "RIGHT_FLANKING_BASE_COUNT_500"), "RIGHT_FLANKING_500_", "MAX_RIGHT_FLANKING_BASE_RATIO_500");
    features["RIGHT_FLANKING_A_RATIO_500"] = features["RIGHT_FLANKING_500_A_RATIO"]; features["RIGHT_FLANKING_C_RATIO_500"] = features["RIGHT_FLANKING_500_C_RATIO"]; features["RIGHT_FLANKING_G_RATIO_500"] = features["RIGHT_FLANKING_500_G_RATIO"]; features["RIGHT_FLANKING_T_RATIO_500"] = features["RIGHT_FLANKING_500_T_RATIO"];
    add_base_count_features(features, get_info_numbers(hdr, record, "SV_REF_PREFIX_BASE_COUNT"), "SV_REF_PREFIX_", "MAX_SV_REF_PREFIX_BASE_RATIO");
    add_base_count_features(features, get_info_numbers(hdr, record, "SV_REF_SUFFIX_BASE_COUNT"), "SV_REF_SUFFIX_", "MAX_SV_REF_SUFFIX_BASE_RATIO");
    add_base_count_features(features, get_info_numbers(hdr, record, "INS_PREFIX_BASE_COUNT"), "INS_PREFIX_", "MAX_INS_PREFIX_BASE_COUNT_RATIO");
    add_base_count_features(features, get_info_numbers(hdr, record, "INS_SUFFIX_BASE_COUNT"), "INS_SUFFIX_", "MAX_INS_SUFFIX_BASE_COUNT_RATIO");
    features["TD"] = get_format_number(hdr, record, "TD", 0);

    double ar1 = get_format_number(hdr, record, "AR1", 0), ar1c = get_format_number(hdr, record, "AR1C", 0), ar1chq = get_format_number(hdr, record, "AR1CHQ", 0);
    double ar1ce = get_format_number(hdr, record, "AR1CE", 0), ar1cehq = get_format_number(hdr, record, "AR1CEHQ", 0), ar1e = get_format_number(hdr, record, "AR1E", 0), ar1ehq = get_format_number(hdr, record, "AR1EHQ", 0);
    double ar1cas = get_format_number(hdr, record, "AR1CAS", NAN_VALUE), ar1css = get_format_number(hdr, record, "AR1CSS", NAN_VALUE);
    double ar1_adj = exp_alt_reads_freq1 > 0 ? ar1/exp_alt_reads_freq1 : ar1, ar1c_adj = exp_alt_reads_freq1 > 0 ? ar1c/exp_alt_reads_freq1 : ar1c;
    features["AR1"] = piecewise_normalise(ar1, min_depth, max_depth); features["AR1C"] = piecewise_normalise(ar1c, min_depth, max_depth);
    features["AR1_ADJ"] = piecewise_normalise(ar1_adj, min_depth, max_depth); features["AR1C_ADJ"] = piecewise_normalise(ar1c_adj, min_depth, max_depth);
    features["AR1C_RATIO"] = ar1c/std::max(1.0, ar1); features["AR1CmQ"] = get_format_number(hdr, record, "AR1CmQ", NAN_VALUE); features["AR1CMQ"] = get_format_number(hdr, record, "AR1CMQ", NAN_VALUE);
    features["AR1CHQ"] = piecewise_normalise(ar1chq, min_depth, max_depth); features["AR1C_HQ_RATIO"] = ar1chq/std::max(1.0, ar1c);
    features["AR1E"] = piecewise_normalise(ar1e, min_depth, max_depth); features["AR1EmQ"] = get_format_number(hdr, record, "AR1EmQ", NAN_VALUE); features["AR1EMQ"] = get_format_number(hdr, record, "AR1EMQ", NAN_VALUE);
    features["AR1E_HQ_RATIO"] = ar1ehq/std::max(1.0, ar1e); features["AR1E_RATIO"] = ar1e/std::max(1.0, ar1);
    features["AR1CE"] = piecewise_normalise(ar1ce, min_depth, max_depth); features["AR1CEmQ"] = get_format_number(hdr, record, "AR1CEmQ", NAN_VALUE); features["AR1CEMQ"] = get_format_number(hdr, record, "AR1CEMQ", NAN_VALUE);
    features["AR1CE_HQ_RATIO"] = ar1cehq/std::max(1.0, ar1ce); features["AR1CE_RATIO"] = ar1ce/std::max(1.0, ar1c);
    std::pair<double, double> spans = normalise_mate_coverage_spans(get_format_numbers(hdr, record, "AR1CMSPAN", {0, 0}), ar1c, read_len);
    features["AR1CMSPAN_1"] = spans.first; features["AR1CMSPAN_2"] = spans.second;
    spans = normalise_mate_coverage_spans(get_format_numbers(hdr, record, "AR1CMHQSPAN", {0, 0}), ar1chq, read_len);
    features["AR1CMHQSPAN_1"] = spans.first; features["AR1CMHQSPAN_2"] = spans.second;
    features["AR1HPMODE"] = get_format_number(hdr, record, "AR1HPMODE", NAN_VALUE); features["AR1CHPMODE"] = get_format_number(hdr, record, "AR1CHPMODE", NAN_VALUE);
    features["AR1HPMODE_AR1CHPMODE_DIFF"] = features["AR1HPMODE"]-features["AR1CHPMODE"];
    features["AR1HPMODE_ALTLEN_DIFF"] = features["AR1HPMODE"]-features["HP_ALT_LEN"]; features["AR1CHPMODE_ALTLEN_DIFF"] = features["AR1CHPMODE"]-features["HP_ALT_LEN"];
    for (const char* key : {"AR1CHPIQR", "AR1CHPmQ", "AR1CHPMQ", "AR1CHPAQ", "AR1CHPSQ", "AR1HP5PMR", "AR1HP3PMR"}) features[key] = get_format_number(hdr, record, key, NAN_VALUE);

    bool has_ar2 = has_format(hdr, record, "AR2");
    double ar2 = get_format_number(hdr, record, "AR2", 0), ar2c = get_format_number(hdr, record, "AR2C", 0), ar2chq = get_format_number(hdr, record, "AR2CHQ", 0);
    double ar2ce = get_format_number(hdr, record, "AR2CE", 0), ar2cehq = get_format_number(hdr, record, "AR2CEHQ", 0), ar2e = get_format_number(hdr, record, "AR2E", 0), ar2ehq = get_format_number(hdr, record, "AR2EHQ", 0);
    double ar2cas = get_format_number(hdr, record, "AR2CAS", NAN_VALUE), ar2css = get_format_number(hdr, record, "AR2CSS", NAN_VALUE);
    if (!has_ar2) ar2 = ar1, ar2c = ar1c, ar2chq = ar1chq, ar2ce = ar1ce, ar2cehq = ar1cehq, ar2e = ar1e, ar2ehq = ar1ehq, ar2cas = ar1cas, ar2css = ar1css;
    double ar2_adj = has_ar2 && exp_alt_reads_freq2 > 0 ? ar2/exp_alt_reads_freq2 : (has_ar2 ? ar2 : ar1_adj), ar2c_adj = has_ar2 && exp_alt_reads_freq2 > 0 ? ar2c/exp_alt_reads_freq2 : (has_ar2 ? ar2c : ar1c_adj);
    features["AR2"] = piecewise_normalise(ar2, min_depth, max_depth); features["AR2C"] = piecewise_normalise(ar2c, min_depth, max_depth);
    features["AR2_ADJ"] = piecewise_normalise(ar2_adj, min_depth, max_depth); features["AR2C_ADJ"] = piecewise_normalise(ar2c_adj, min_depth, max_depth);
    features["AR2C_RATIO"] = ar2c/std::max(1.0, ar2); features["AR2CmQ"] = has_ar2 ? get_format_number(hdr, record, "AR2CmQ", NAN_VALUE) : features["AR1CmQ"]; features["AR2CMQ"] = has_ar2 ? get_format_number(hdr, record, "AR2CMQ", NAN_VALUE) : features["AR1CMQ"];
    features["AR2CHQ"] = piecewise_normalise(ar2chq, min_depth, max_depth); features["AR2C_HQ_RATIO"] = ar2chq/std::max(1.0, ar2c);
    features["AR2E"] = piecewise_normalise(ar2e, min_depth, max_depth); features["AR2EmQ"] = has_ar2 ? get_format_number(hdr, record, "AR2EmQ", NAN_VALUE) : features["AR1EmQ"]; features["AR2EMQ"] = has_ar2 ? get_format_number(hdr, record, "AR2EMQ", NAN_VALUE) : features["AR1EMQ"];
    features["AR2E_HQ_RATIO"] = ar2ehq/std::max(1.0, ar2e); features["AR2E_RATIO"] = ar2e/std::max(1.0, ar2);
    features["AR2CE"] = piecewise_normalise(ar2ce, min_depth, max_depth); features["AR2CEmQ"] = has_ar2 ? get_format_number(hdr, record, "AR2CEmQ", NAN_VALUE) : features["AR1CEmQ"]; features["AR2CEMQ"] = has_ar2 ? get_format_number(hdr, record, "AR2CEMQ", NAN_VALUE) : features["AR1CEMQ"];
    features["AR2CE_HQ_RATIO"] = ar2cehq/std::max(1.0, ar2ce); features["AR2CE_RATIO"] = ar2ce/std::max(1.0, ar2c);
    spans = has_ar2 ? normalise_mate_coverage_spans(get_format_numbers(hdr, record, "AR2CMSPAN", {0, 0}), ar2c, read_len) : std::make_pair(features["AR1CMSPAN_1"], features["AR1CMSPAN_2"]);
    features["AR2CMSPAN_1"] = spans.first; features["AR2CMSPAN_2"] = spans.second;
    spans = has_ar2 ? normalise_mate_coverage_spans(get_format_numbers(hdr, record, "AR2CMHQSPAN", {0, 0}), ar2chq, read_len) : std::make_pair(features["AR1CMHQSPAN_1"], features["AR1CMHQSPAN_2"]);
    features["AR2CMHQSPAN_1"] = spans.first; features["AR2CMHQSPAN_2"] = spans.second;

    double ar1cf_count = get_format_number(hdr, record, "AR1CF", 0), ar1cr_count = get_format_number(hdr, record, "AR1CR", 0), ar2cf_count = get_format_number(hdr, record, "AR2CF", 0), ar2cr_count = get_format_number(hdr, record, "AR2CR", 0);
    double ar1cf = ar1cf_count/std::max(1.0, ar1c), ar1cr = ar1cr_count/std::max(1.0, ar1c), ar2cf = has_ar2 ? ar2cf_count/std::max(1.0, ar2c) : ar1cf, ar2cr = has_ar2 ? ar2cr_count/std::max(1.0, ar2c) : ar1cr;
    features["ARCF"] = ar1cf+ar2cf; features["ARCR"] = ar1cr+ar2cr; features["MAXARCD"] = std::max(features["ARCF"], features["ARCR"]);
    double ar1cef = get_format_number(hdr, record, "AR1CEF", 0, std::max(1.0, ar1ce)), ar1cer = get_format_number(hdr, record, "AR1CER", 0, std::max(1.0, ar1ce));
    double ar2cef = has_ar2 ? get_format_number(hdr, record, "AR2CEF", 0, std::max(1.0, ar2ce)) : ar1cef, ar2cer = has_ar2 ? get_format_number(hdr, record, "AR2CER", 0, std::max(1.0, ar2ce)) : ar1cer;
    features["ARCEF"] = ar1cef+ar2cef; features["ARCER"] = ar1cer+ar2cer; features["MAXARCED"] = std::max(features["ARCEF"], features["ARCER"]);
    double ar1ef = get_format_number(hdr, record, "AR1EF", 0, std::max(1.0, ar1e)), ar1er = get_format_number(hdr, record, "AR1ER", 0, std::max(1.0, ar1e));
    double ar2ef = has_ar2 ? get_format_number(hdr, record, "AR2EF", 0, std::max(1.0, ar2e)) : ar1ef, ar2er = has_ar2 ? get_format_number(hdr, record, "AR2ER", 0, std::max(1.0, ar2e)) : ar1er;
    features["AREF"] = ar1ef+ar2ef; features["ARER"] = ar1er+ar2er; features["MAXARED"] = std::max(features["AREF"], features["ARER"]);

    selected_alt_read_metrics_t other1 = select_alt_read_metrics_for_oar_vids(hdr, record, "OAR1VID", 0, alt_reads_by_vid), other2 = select_alt_read_metrics_for_oar_vids(hdr, record, "OAR2VID", 1, alt_reads_by_vid);
    double oar1 = other1.oar, oar1c = other1.oarc, oar1chq = other1.oarchq, oar1e = other1.oare, oar1all = get_format_number(hdr, record, "OAR1ALL", 0);
    double oar2 = other2.oar, oar2c = other2.oarc, oar2chq = other2.oarchq, oar2e = other2.oare, oar2all = get_format_number(hdr, record, "OAR2ALL", 0);
    if (!has_ar2) other2 = other1, oar2 = oar1, oar2c = oar1c, oar2chq = oar1chq, oar2e = oar1e, oar2all = oar1all;
    features["OTHER_HP_GENOTYPED"] = std::max(other1.hp_genotyped, other2.hp_genotyped);
    for (int bp = 1; bp <= 2; bp++) {
        selected_alt_read_metrics_t other = bp == 1 ? other1 : other2;
        double oar = bp == 1 ? oar1 : oar2, oarc = bp == 1 ? oar1c : oar2c, oarchq = bp == 1 ? oar1chq : oar2chq, oare = bp == 1 ? oar1e : oar2e, oarall = bp == 1 ? oar1all : oar2all;
        std::string n = std::to_string(bp);
        features["OAR"+n] = piecewise_normalise(oar, min_depth, max_depth); features["OAR"+n+"ALL"] = piecewise_normalise(oarall, min_depth, max_depth); features["OAR"+n+"_ALL_RATIO"] = oar/std::max(1.0, oarall);
        features["OTHER"+n+"_MIN_ARCHQ"] = piecewise_normalise(other.min_archq, min_depth, max_depth); features["OTHER"+n+"_MIN_AR_OVER_NAR"] = other.min_ar_over_nar; features["OTHER"+n+"_MIN_ARC_OVER_NARC"] = other.min_arc_over_narc; features["OTHER"+n+"_MIN_ARE_OVER_NARE"] = other.min_are_over_nare;
        features["OTHER"+n+"_AAS_ARS_DIFF_TO_LEN"] = other.aas_ars_diff_to_len; features["OTHER"+n+"_XAAS_XARS_DIFF_TO_LEN"] = other.xaas_xars_diff_to_len;
        features["OAR"+n+"C"] = piecewise_normalise(oarc, min_depth, max_depth); features["OAR"+n+"CHQ"] = piecewise_normalise(oarchq, min_depth, max_depth); features["OAR"+n+"C_HQ_RATIO"] = oarchq/std::max(1.0, oarc); features["OAR"+n+"E"] = piecewise_normalise(oare, min_depth, max_depth);
    }

    double orr1 = get_format_number(hdr, record, "ORR1", 0), orr1c = get_format_number(hdr, record, "ORR1C", 0), orr1e = get_format_number(hdr, record, "ORR1E", 0);
    double orr2 = get_format_number(hdr, record, "ORR2", 0), orr2c = get_format_number(hdr, record, "ORR2C", 0), orr2e = get_format_number(hdr, record, "ORR2E", 0);
    if (!has_format(hdr, record, "ORR2")) orr2 = orr1, orr2c = orr1c, orr2e = orr1e;
    double rr1 = get_format_number(hdr, record, "RR1", 0), rr1c = get_format_number(hdr, record, "RR1C", 0), rr1chq = get_format_number(hdr, record, "RR1CHQ", 0), rr1e = get_format_number(hdr, record, "RR1E", 0), rr1ehq = get_format_number(hdr, record, "RR1EHQ", 0);
    double rr1ce = get_format_number(hdr, record, "RR1CE", 0), rr1cehq = get_format_number(hdr, record, "RR1CEHQ", 0), rr1cas = get_format_number(hdr, record, "RR1CAS", NAN_VALUE), rr1css = get_format_number(hdr, record, "RR1CSS", NAN_VALUE);
    double rr2 = get_format_number(hdr, record, "RR2", 0), rr2c = get_format_number(hdr, record, "RR2C", 0), rr2chq = get_format_number(hdr, record, "RR2CHQ", 0), rr2e = get_format_number(hdr, record, "RR2E", 0), rr2ehq = get_format_number(hdr, record, "RR2EHQ", 0);
    double rr2ce = get_format_number(hdr, record, "RR2CE", 0), rr2cehq = get_format_number(hdr, record, "RR2CEHQ", 0), rr2cas = get_format_number(hdr, record, "RR2CAS", NAN_VALUE), rr2css = get_format_number(hdr, record, "RR2CSS", NAN_VALUE);
    bool has_rr2 = has_format(hdr, record, "RR2");
    if (!has_rr2) rr2 = rr1, rr2c = rr1c, rr2chq = rr1chq, rr2e = rr1e, rr2ehq = rr1ehq, rr2ce = rr1ce, rr2cehq = rr1cehq, rr2cas = rr1cas, rr2css = rr1css;
    for (int bp = 1; bp <= 2; bp++) {
        bool first = bp == 1;
        double rr = first ? rr1 : rr2, rrc = first ? rr1c : rr2c, rrchq = first ? rr1chq : rr2chq, rre = first ? rr1e : rr2e, rrehq = first ? rr1ehq : rr2ehq, rrce = first ? rr1ce : rr2ce, rrcehq = first ? rr1cehq : rr2cehq;
        std::string n = std::to_string(bp), prefix = "RR"+n;
        features[prefix] = piecewise_normalise(rr, min_depth, max_depth); features[prefix+"C"] = piecewise_normalise(rrc, min_depth, max_depth);
        features[prefix+"CmQ"] = !first && !has_rr2 ? features["RR1CmQ"] : get_format_number(hdr, record, (prefix+"CmQ").c_str(), NAN_VALUE); features[prefix+"CMQ"] = !first && !has_rr2 ? features["RR1CMQ"] : get_format_number(hdr, record, (prefix+"CMQ").c_str(), NAN_VALUE);
        features[prefix+"C_HQ_RATIO"] = rrchq/std::max(1.0, rrc); features[prefix+"E"] = piecewise_normalise(rre, min_depth, max_depth);
        features[prefix+"EmQ"] = !first && !has_rr2 ? features["RR1EmQ"] : get_format_number(hdr, record, (prefix+"EmQ").c_str(), NAN_VALUE); features[prefix+"EMQ"] = !first && !has_rr2 ? features["RR1EMQ"] : get_format_number(hdr, record, (prefix+"EMQ").c_str(), NAN_VALUE);
        features[prefix+"E_HQ_RATIO"] = rrehq/std::max(1.0, rre); features[prefix+"E_RATIO"] = rre/std::max(1.0, rr); features[prefix+"CE"] = piecewise_normalise(rrce, min_depth, max_depth);
        features[prefix+"CEmQ"] = !first && !has_rr2 ? features["RR1CEmQ"] : get_format_number(hdr, record, (prefix+"CEmQ").c_str(), NAN_VALUE); features[prefix+"CEMQ"] = !first && !has_rr2 ? features["RR1CEMQ"] : get_format_number(hdr, record, (prefix+"CEMQ").c_str(), NAN_VALUE);
        features[prefix+"CE_HQ_RATIO"] = rrcehq/std::max(1.0, rrce); features[prefix+"CE_RATIO"] = rrce/std::max(1.0, rrc);
        spans = !first && !has_rr2 ? std::make_pair(features["RR1CMSPAN_1"], features["RR1CMSPAN_2"]) : normalise_mate_coverage_spans(get_format_numbers(hdr, record, (prefix+"CMSPAN").c_str(), {0, 0}), rrc, read_len);
        features[prefix+"CMSPAN_1"] = spans.first; features[prefix+"CMSPAN_2"] = spans.second;
        spans = !first && !has_rr2 ? std::make_pair(features["RR1CMHQSPAN_1"], features["RR1CMHQSPAN_2"]) : normalise_mate_coverage_spans(get_format_numbers(hdr, record, (prefix+"CMHQSPAN").c_str(), {0, 0}), rrchq, read_len);
        features[prefix+"CMHQSPAN_1"] = spans.first; features[prefix+"CMHQSPAN_2"] = spans.second;
    }
    features["RR1HPMODE"] = get_format_number(hdr, record, "RR1HPMODE", NAN_VALUE); features["RR1CHPMODE"] = get_format_number(hdr, record, "RR1CHPMODE", NAN_VALUE);
    features["RR1HPMODE_RR1CHPMODE_DIFF"] = features["RR1HPMODE"]-features["RR1CHPMODE"];
    features["RR1HPMODE_REFLEN_DIFF"] = features["RR1HPMODE"]-features["HP_REF_LEN"]; features["RR1CHPMODE_REFLEN_DIFF"] = features["RR1CHPMODE"]-features["HP_REF_LEN"];
    for (const char* key : {"RR1CHPIQR", "RR1CHPmQ", "RR1CHPMQ", "RR1CHPAQ", "RR1CHPSQ", "RR1HP5PMR", "RR1HP3PMR"}) features[key] = get_format_number(hdr, record, key, NAN_VALUE);
    double rr1cf = get_format_number(hdr, record, "RR1CF", 0), rr1cr = get_format_number(hdr, record, "RR1CR", 0), rr2cf = has_rr2 ? get_format_number(hdr, record, "RR2CF", 0) : rr1cf, rr2cr = has_rr2 ? get_format_number(hdr, record, "RR2CR", 0) : rr1cr;
    double rr1cf_ratio = rr1cf/std::max(1.0, rr1cf+rr1cr), rr1cr_ratio = rr1cr/std::max(1.0, rr1cf+rr1cr), rr2cf_ratio = rr2cf/std::max(1.0, rr2cf+rr2cr), rr2cr_ratio = rr2cr/std::max(1.0, rr2cf+rr2cr);
    features["RRCF"] = rr1cf_ratio+rr2cf_ratio; features["RRCR"] = rr1cr_ratio+rr2cr_ratio; features["MAXRRCD"] = std::max(features["RRCF"], features["RRCR"]);
    double rr1cef = get_format_number(hdr, record, "RR1CEF", 0, std::max(1.0, rr1ce)), rr1cer = get_format_number(hdr, record, "RR1CER", 0, std::max(1.0, rr1ce));
    double rr2cef = has_rr2 ? get_format_number(hdr, record, "RR2CEF", 0, std::max(1.0, rr2ce)) : rr1cef, rr2cer = has_rr2 ? get_format_number(hdr, record, "RR2CER", 0, std::max(1.0, rr2ce)) : rr1cer;
    features["RRCEF"] = rr1cef+rr2cef; features["RRCER"] = rr1cer+rr2cer; features["MAXRRCED"] = std::max(features["RRCEF"], features["RRCER"]);
    double rr1ef = get_format_number(hdr, record, "RR1EF", 0, std::max(1.0, rr1e)), rr1er = get_format_number(hdr, record, "RR1ER", 0, std::max(1.0, rr1e));
    double rr2ef = has_rr2 ? get_format_number(hdr, record, "RR2EF", 0, std::max(1.0, rr2e)) : rr1ef, rr2er = has_rr2 ? get_format_number(hdr, record, "RR2ER", 0, std::max(1.0, rr2e)) : rr1er;
    features["RREF"] = rr1ef+rr2ef; features["RRER"] = rr1er+rr2er; features["MAXRRED"] = std::max(features["RREF"], features["RRER"]);

    double nar1 = rr1+oar1, nar1c = rr1c+oar1c, nar1chq = rr1chq+oar1chq, nar1ce = rr1ce+oar1e, nar1e = rr1e+oar1e;
    double nar2 = rr2+oar2, nar2c = rr2c+oar2c, nar2chq = rr2chq+oar2chq, nar2ce = rr2ce+oar2e, nar2e = rr2e+oar2e;
    for (int bp = 1; bp <= 2; bp++) {
        std::string n = std::to_string(bp); double nar = bp == 1 ? nar1 : nar2, narc = bp == 1 ? nar1c : nar2c, narchq = bp == 1 ? nar1chq : nar2chq, narce = bp == 1 ? nar1ce : nar2ce, nare = bp == 1 ? nar1e : nar2e;
        features["NAR"+n] = piecewise_normalise(nar, min_depth, max_depth); features["NAR"+n+"C"] = piecewise_normalise(narc, min_depth, max_depth); features["NAR"+n+"CHQ"] = piecewise_normalise(narchq, min_depth, max_depth);
        features["NAR"+n+"C_HQ_RATIO"] = narchq/std::max(1.0, narc); features["NAR"+n+"CE"] = piecewise_normalise(narce, min_depth, max_depth); features["NAR"+n+"E"] = piecewise_normalise(nare, min_depth, max_depth);
    }
    double er = get_format_number(hdr, record, "ER", 0), erhq = get_format_number(hdr, record, "ERHQ", 0), er1_total = er+ar1+rr1, er2_total = er+ar2+rr2;
    features["ER1_DEVIATION"] = er1_total > 0 ? er/er1_total-(1-exp_alt_reads_freq1) : NAN_VALUE; features["ER2_DEVIATION"] = er2_total > 0 ? er/er2_total-(1-exp_alt_reads_freq2) : NAN_VALUE; features["ER_HQ_RATIO"] = er > 0 ? erhq/er : NAN_VALUE;
    features["AR1_RR1_CAS_Z_SCORE"] = calculate_z_score(ar1cas, ar1css, ar1c, rr1cas, rr1css, rr1c); features["AR2_RR2_CAS_Z_SCORE"] = calculate_z_score(ar2cas, ar2css, ar2c, rr2cas, rr2css, rr2c);

#define ADD_RATIO(NAME, NUM, DEN) features[NAME] = (NUM)/std::max(1.0, (DEN))
    ADD_RATIO("AR1_OVER_RR1", ar1, ar1+rr1); ADD_RATIO("AR2_OVER_RR2", ar2, ar2+rr2); ADD_RATIO("AR1C_OVER_RR1C", ar1c, ar1c+rr1c); ADD_RATIO("AR2C_OVER_RR2C", ar2c, ar2c+rr2c);
    ADD_RATIO("AR1CE_OVER_RR1CE", ar1ce, ar1ce+rr1ce); ADD_RATIO("AR2CE_OVER_RR2CE", ar2ce, ar2ce+rr2ce); ADD_RATIO("AR1E_OVER_RR1E", ar1e, ar1e+rr1e); ADD_RATIO("AR2E_OVER_RR2E", ar2e, ar2e+rr2e);
    ADD_RATIO("AR1_OVER_OAR1", ar1, ar1+oar1); ADD_RATIO("AR2_OVER_OAR2", ar2, ar2+oar2); ADD_RATIO("AR1C_OVER_OAR1C", ar1c, ar1c+oar1c); ADD_RATIO("AR2C_OVER_OAR2C", ar2c, ar2c+oar2c);
    ADD_RATIO("AR1CE_OVER_OAR1E", ar1ce, ar1ce+oar1e); ADD_RATIO("AR2CE_OVER_OAR2E", ar2ce, ar2ce+oar2e); ADD_RATIO("AR1E_OVER_OAR1E", ar1e, ar1e+oar1e); ADD_RATIO("AR2E_OVER_OAR2E", ar2e, ar2e+oar2e);
    ADD_RATIO("OAR1_OVER_NAR1", oar1, nar1); ADD_RATIO("OAR2_OVER_NAR2", oar2, nar2); ADD_RATIO("OAR1C_OVER_NAR1C", oar1c, nar1c); ADD_RATIO("OAR2C_OVER_NAR2C", oar2c, nar2c);
    ADD_RATIO("OAR1E_OVER_NAR1CE", oar1e, nar1ce); ADD_RATIO("OAR2E_OVER_NAR2CE", oar2e, nar2ce); ADD_RATIO("OAR1E_OVER_NAR1E", oar1e, nar1e); ADD_RATIO("OAR2E_OVER_NAR2E", oar2e, nar2e);
    ADD_RATIO("OAR1_OVER_TOTAL1", oar1, nar1+ar1); ADD_RATIO("OAR2_OVER_TOTAL2", oar2, nar2+ar2); ADD_RATIO("OAR1C_OVER_TOTAL1C", oar1c, nar1c+ar1c); ADD_RATIO("OAR2C_OVER_TOTAL2C", oar2c, nar2c+ar2c);
    ADD_RATIO("OAR1E_OVER_TOTAL1CE", oar1e, nar1ce+ar1ce); ADD_RATIO("OAR2E_OVER_TOTAL2CE", oar2e, nar2ce+ar2ce); ADD_RATIO("OAR1E_OVER_TOTAL1E", oar1e, nar1e+ar1e); ADD_RATIO("OAR2E_OVER_TOTAL2E", oar2e, nar2e+ar2e);
    ADD_RATIO("ORR1_RATIO", orr1, rr1); ADD_RATIO("ORR2_RATIO", orr2, rr2); ADD_RATIO("ORR1C_RATIO", orr1c, rr1); ADD_RATIO("ORR2C_RATIO", orr2c, rr2); ADD_RATIO("ORR1E_RATIO", orr1e, rr1); ADD_RATIO("ORR2E_RATIO", orr2e, rr2);
    ADD_RATIO("AR1_OVER_NAR1", ar1, ar1+nar1); ADD_RATIO("AR2_OVER_NAR2", ar2, ar2+nar2); ADD_RATIO("AR1C_OVER_NAR1C", ar1c, ar1c+nar1c); ADD_RATIO("AR2C_OVER_NAR2C", ar2c, ar2c+nar2c);
    ADD_RATIO("AR1CE_OVER_NAR1CE", ar1ce, ar1ce+nar1ce); ADD_RATIO("AR2CE_OVER_NAR2CE", ar2ce, ar2ce+nar2ce); ADD_RATIO("AR1E_OVER_NAR1E", ar1e, ar1e+nar1e); ADD_RATIO("AR2E_OVER_NAR2E", ar2e, ar2e+nar2e);
#undef ADD_RATIO

    std::vector<double> md = get_format_numbers(hdr, record, "MD", {0, 0, 0, 0});
    features["MDLF"] = piecewise_normalise(md[0], min_depth, max_depth); features["MDSP"] = piecewise_normalise(md[1], min_depth, max_depth); features["MDSF"] = piecewise_normalise(md[2], min_depth, max_depth); features["MDRF"] = piecewise_normalise(md[3], min_depth, max_depth);
    features["MDSP_OVER_MDLF"] = md[1]/std::max(1.0, md[0]); features["MDSF_OVER_MDRF"] = md[2]/std::max(1.0, md[3]);
    std::vector<double> mdhq = get_format_numbers(hdr, record, "MDHQ", {0, 0, 0, 0});
    features["MDLFHQ"] = piecewise_normalise(mdhq[0], min_depth, max_depth); features["MDSPHQ"] = piecewise_normalise(mdhq[1], min_depth, max_depth); features["MDSFHQ"] = piecewise_normalise(mdhq[2], min_depth, max_depth); features["MDRFHQ"] = piecewise_normalise(mdhq[3], min_depth, max_depth);
    features["MDSP_OVER_MDLF_HQ"] = piecewise_normalise(mdhq[1]-mdhq[0], min_depth, max_depth); features["MDSF_OVER_MDRF_HQ"] = piecewise_normalise(mdhq[2]-mdhq[3], min_depth, max_depth);
    std::vector<double> clmd = get_format_numbers(hdr, record, "CLMD", {0, 0}), clmdhq = get_format_numbers(hdr, record, "CLMDHQ", {0, 0});
    features["MDLC"] = piecewise_normalise(clmd[0], min_depth, max_depth); features["MDRC"] = piecewise_normalise(clmd[1], min_depth, max_depth); features["MDLCHQ"] = piecewise_normalise(clmdhq[0], min_depth, max_depth); features["MDRCHQ"] = piecewise_normalise(clmdhq[1], min_depth, max_depth);
    features["KS_PVAL"] = get_format_number(hdr, record, "KSPVAL", NAN_VALUE); features["SIZE_NORM"] = NAN_VALUE;
    if (has_format(hdr, record, "MAXSIZE")) {
        double min_size = get_format_number(hdr, record, "MINSIZE", 0), max_size = get_format_number(hdr, record, "MAXSIZE", 0);
        features["SIZE_NORM"] = max_size == min_size ? svlen/2-min_size : (svlen/2-min_size)/(max_size-min_size);
    }

    double min_disc_pairs, max_disc_pairs;
    if (svtype == "DEL") {
        std::string key = std::to_string((int)std::max(0.0, max_is-svlen));
        min_disc_pairs = stats.at("min_pairs_crossing_gap").at(key); max_disc_pairs = stats.at("max_pairs_crossing_gap").at(key);
    } else if (svtype == "INS" && (source == "DE_NOVO_ASSEMBLY" || source == "REFERENCE_GUIDED_ASSEMBLY")) {
        std::string key = std::to_string((int)std::min<double>(max_is, svinslen));
        min_disc_pairs = stats.at("min_disc_pairs_by_insertion_size").at(key); max_disc_pairs = stats.at("max_disc_pairs_by_insertion_size").at(key);
    } else min_disc_pairs = min_pairs_crossing_point, max_disc_pairs = max_pairs_crossing_point;

    double asp1, asp1nma1, asp1nma2, asp1nms1, asp1nms2, asp2, asp2nma1, asp2nma2, asp2nms1, asp2nms2;
    double rsp1, rsp1nma1, rsp1nma2, rsp1nms1, rsp1nms2, rsp2, rsp2nma1, rsp2nma2, rsp2nms1, rsp2nms2;
    double nsp1, nsp1nma1, nsp1nma2, nsp1nms1, nsp1nms2, nsp2, nsp2nma1, nsp2nma2, nsp2nms1, nsp2nms2;
    add_read_pair_features(features, hdr, record, "ASP1", min_disc_pairs, max_disc_pairs, asp1, asp1nma1, asp1nma2, asp1nms1, asp1nms2);
    add_read_pair_features(features, hdr, record, "ASP2", min_disc_pairs, max_disc_pairs, asp2, asp2nma1, asp2nma2, asp2nms1, asp2nms2);
    add_read_pair_features(features, hdr, record, "RSP1", min_disc_pairs, max_disc_pairs, rsp1, rsp1nma1, rsp1nma2, rsp1nms1, rsp1nms2);
    add_read_pair_features(features, hdr, record, "RSP2", min_disc_pairs, max_disc_pairs, rsp2, rsp2nma1, rsp2nma2, rsp2nms1, rsp2nms2);
    add_read_pair_features(features, hdr, record, "NSP1", min_disc_pairs, max_disc_pairs, nsp1, nsp1nma1, nsp1nma2, nsp1nms1, nsp1nms2);
    add_read_pair_features(features, hdr, record, "NSP2", min_disc_pairs, max_disc_pairs, nsp2, nsp2nma1, nsp2nma2, nsp2nms1, nsp2nms2);
    features["ASP1_ASP2_RATIO"] = std::max(asp1, asp2)/std::max(1.0, asp1+asp2);
    std::vector<double> asp1span = get_format_numbers(hdr, record, "ASP1SPAN", {0, 0}), asp2span = get_format_numbers(hdr, record, "ASP2SPAN", {0, 0});
    features["ASP1SPAN_1"] = asp1span[0]/max_is; features["ASP2SPAN_2"] = asp2span[1]/max_is;
    features["ASP1SPAN_2"] = asp1span[1]/(svtype == "INS" ? std::max<double>(1, std::max<int64_t>(max_is, svinslen)) : max_is);
    features["ASP2SPAN_1"] = asp2span[0]/(svtype == "INS" ? std::max<double>(1, std::max<int64_t>(max_is, svinslen)) : max_is);

    if (!has_format(hdr, record, "ASP2")) asp2 = asp1, asp2nma1 = asp1nma1, asp2nma2 = asp1nma2, asp2nms1 = asp1nms1, asp2nms2 = asp1nms2;
    if (!has_format(hdr, record, "RSP2")) rsp2 = rsp1, rsp2nma1 = rsp1nma1, rsp2nma2 = rsp1nma2, rsp2nms1 = rsp1nms1, rsp2nms2 = rsp1nms2;
    if (!has_format(hdr, record, "NSP2")) nsp2 = nsp1, nsp2nma1 = nsp1nma1, nsp2nma2 = nsp1nma2, nsp2nms1 = nsp1nms1, nsp2nms2 = nsp1nms2;
    features["ASP1_OVER_RSP1"] = asp1/std::max(1.0, asp1+rsp1); features["ASP2_OVER_RSP2"] = asp2/std::max(1.0, asp2+rsp2);
#define ADD_Z(NAME, M1, S1, N1, M2, S2, N2) features[NAME] = calculate_z_score(M1, S1, N1, M2, S2, N2)
    ADD_Z("ASP1_RSP1_1_NM_Z_SCORE", asp1nma1, asp1nms1, asp1, rsp1nma1, rsp1nms1, rsp1); ADD_Z("ASP1_RSP1_2_NM_Z_SCORE", asp1nma2, asp1nms2, asp1, rsp1nma2, rsp1nms2, rsp1);
    ADD_Z("ASP2_RSP2_1_NM_Z_SCORE", asp2nma1, asp2nms1, asp2, rsp2nma1, rsp2nms1, rsp2); ADD_Z("ASP2_RSP2_2_NM_Z_SCORE", asp2nma2, asp2nms2, asp2, rsp2nma2, rsp2nms2, rsp2);
    ADD_Z("ASP1_NSP1_1_NM_Z_SCORE", asp1nma1, asp1nms1, asp1, nsp1nma1, nsp1nms1, nsp1); ADD_Z("ASP1_NSP1_2_NM_Z_SCORE", asp1nma2, asp1nms2, asp1, nsp1nma2, nsp1nms2, nsp1);
    ADD_Z("ASP2_NSP2_1_NM_Z_SCORE", asp2nma1, asp2nms1, asp2, nsp2nma1, nsp2nms1, nsp2); ADD_Z("ASP2_NSP2_2_NM_Z_SCORE", asp2nma2, asp2nms2, asp2, nsp2nma2, nsp2nms2, nsp2);

    double ssp1, ssp1nma1, ssp1nma2, ssp1nms1, ssp1nms2, ssp2, ssp2nma1, ssp2nma2, ssp2nms1, ssp2nms2;
    add_read_pair_features(features, hdr, record, "SSP1", min_disc_pairs, max_disc_pairs, ssp1, ssp1nma1, ssp1nma2, ssp1nms1, ssp1nms2);
    add_read_pair_features(features, hdr, record, "SSP2", min_disc_pairs, max_disc_pairs, ssp2, ssp2nma1, ssp2nma2, ssp2nms1, ssp2nms2);
    if (svtype == "DEL" || svtype == "INS") ssp1nma2 = NAN_VALUE, ssp2nma1 = NAN_VALUE;
    else if (svtype == "DUP") ssp1nma1 = NAN_VALUE, ssp2nma2 = NAN_VALUE;
    ADD_Z("SSP1_RSP1_1_NM_Z_SCORE", ssp1nma1, ssp1nms1, ssp1, rsp1nma1, rsp1nms1, rsp1); ADD_Z("SSP1_RSP1_2_NM_Z_SCORE", ssp1nma2, ssp1nms2, ssp1, rsp1nma2, rsp1nms2, rsp1);
    ADD_Z("SSP2_RSP2_1_NM_Z_SCORE", ssp2nma1, ssp2nms1, ssp2, rsp2nma1, rsp2nms1, rsp2); ADD_Z("SSP2_RSP2_2_NM_Z_SCORE", ssp2nma2, ssp2nms2, ssp2, rsp2nma2, rsp2nms2, rsp2);
#undef ADD_Z
    add_consensus_alignment_features(features, hdr, record, "", read_len, edit_distance);
    add_consensus_alignment_features(features, hdr, record, "X", read_len, edit_distance);

    std::vector<double> values;
    values.reserve(feature_names.size());
    for (const std::string& feature_name : feature_names) {
        auto it = features.find(feature_name);
        if (it == features.end()) throw std::runtime_error("Feature '" + feature_name + "' required for model " + model_name + " is not produced by extract_features.");
        values.push_back(it->second);
    }
    return values;
}

std::map<std::string, model_data_t> load_model_features(const std::string& model_stage_dir) {
    DIR* dir = opendir(model_stage_dir.c_str());
    if (dir == NULL) throw std::runtime_error("Failed to open model directory " + model_stage_dir + ".");
    std::map<std::string, model_data_t> models;
    dirent* entry;
    while ((entry = readdir(dir)) != NULL) {
        std::string fname = entry->d_name;
        if (fname.size() < 9 || fname.substr(fname.size()-9) != ".features") continue;
        std::string model_name = fname.substr(0, fname.size()-9);
        std::ifstream in(model_stage_dir+"/"+fname);
        if (!in) { closedir(dir); throw std::runtime_error("Failed to open " + model_stage_dir+"/"+fname + "."); }
        std::string feature_name;
        while (std::getline(in, feature_name)) if (!feature_name.empty()) models[model_name].feature_names.push_back(feature_name);
    }
    closedir(dir);
    if (models.empty()) throw std::runtime_error("No .features files found in " + model_stage_dir + ".");
    return models;
}

bcf_srs_t* open_region_reader(const std::string& vcf_fname, const std::string& chrom) {
    bcf_srs_t* reader = bcf_sr_init();
    bcf_sr_set_opt(reader, BCF_SR_REQUIRE_IDX);
    std::string region = "{"+chrom+"}";
    if (bcf_sr_set_regions(reader, region.c_str(), 0) < 0 || !bcf_sr_add_reader(reader, vcf_fname.c_str())) {
        std::string error = bcf_sr_strerror(reader->errnum);
        bcf_sr_destroy(reader);
        throw std::runtime_error("Failed to open " + vcf_fname + " for chromosome " + chrom + ": " + error + ".");
    }
    return reader;
}

std::unordered_map<std::string, std::vector<alt_read_metrics_t>> load_alt_read_metrics(const std::string& vcf_fname, const std::string& chrom) {
    bcf_srs_t* reader = open_region_reader(vcf_fname, chrom);
    bcf_hdr_t* hdr = reader->readers[0].header;
    std::unordered_map<std::string, std::vector<alt_read_metrics_t>> alt_reads_by_vid;
    while (bcf_sr_next_line(reader)) {
        bcf1_t* record = bcf_sr_get_line(reader, 0);
        bcf_unpack(record, BCF_UN_ALL);
        double ar1 = get_format_number(hdr, record, "AR1", 0), ar2 = get_format_number(hdr, record, "AR2", ar1);
        double rr1 = get_format_number(hdr, record, "RR1", 0), rr2 = get_format_number(hdr, record, "RR2", rr1), rr1c = get_format_number(hdr, record, "RR1C", 0), rr2c = get_format_number(hdr, record, "RR2C", rr1c), rr1e = get_format_number(hdr, record, "RR1E", 0), rr2e = get_format_number(hdr, record, "RR2E", rr1e);
        double ar1c = get_format_number(hdr, record, "AR1C", 0), ar2c = get_format_number(hdr, record, "AR2C", ar1c), ar1chq = get_format_number(hdr, record, "AR1CHQ", 0), ar2chq = get_format_number(hdr, record, "AR2CHQ", ar1chq), ar1e = get_format_number(hdr, record, "AR1E", 0), ar2e = get_format_number(hdr, record, "AR2E", ar1e);
        std::vector<std::string> oar1_vids = get_oar_vids(hdr, record, "OAR1VID"), oar2_vids = has_format(hdr, record, "AR2") ? get_oar_vids(hdr, record, "OAR2VID") : oar1_vids;
        std::string id = record->d.id == NULL ? "." : record->d.id;
        alt_reads_by_vid[normalize_sv_id(id)].push_back({ar1, ar2, ar1c, ar2c, ar1chq, ar2chq, ar1e, ar2e, gt_as_homopolymer(hdr, record), bcf_seqname(hdr, record), record->pos, record_stop(record), rr1, rr2, rr1c, rr2c, rr1e, rr2e, oar1_vids, oar2_vids, consensus_alt_ref_score_diff_to_len(hdr, record, "X"), has_extension_evidence(hdr, record), consensus_alt_ref_score_diff_to_len(hdr, record, ""), has_assembly_evidence(hdr, record)});
    }
    bcf_sr_destroy(reader);
    return alt_reads_by_vid;
}

std::map<std::string, model_data_t> process_chrom(const std::string& vcf_fname, const std::string& chrom, const stats_t& stats, const std::map<std::string, model_data_t>& model_definitions) {
    std::unordered_map<std::string, std::vector<alt_read_metrics_t>> alt_reads_by_vid = load_alt_read_metrics(vcf_fname, chrom);
    std::map<std::string, model_data_t> models = model_definitions;
    bcf_srs_t* reader = open_region_reader(vcf_fname, chrom);
    bcf_hdr_t* hdr = reader->readers[0].header;
    uint64_t record_idx = 0;
    while (bcf_sr_next_line(reader)) {
        bcf1_t* record = bcf_sr_get_line(reader, 0);
        bcf_unpack(record, BCF_UN_ALL);
        uint64_t record_key = (uint64_t(uint32_t(record->rid))<<32)|record_idx++;
        if (get_svtype(hdr, record).find("INV") == 0) continue;
        std::string model_name = get_model_name(hdr, record, get_stat(stats, "max_is", chrom), get_stat(stats, "read_len", chrom));
        auto model_it = models.find(model_name);
        if (model_it == models.end()) { bcf_sr_destroy(reader); throw std::runtime_error("No .features file found for model " + model_name + "."); }
        std::vector<double> values = record_to_features(hdr, record, stats, model_it->second.feature_names, alt_reads_by_vid);
        model_it->second.variant_ids.push_back(generate_id(hdr, record, model_name));
        model_it->second.record_keys.push_back(record_key);
        model_it->second.values.insert(model_it->second.values.end(), values.begin(), values.end());
    }
    bcf_sr_destroy(reader);
    return models;
}

void extract_features(const std::string& vcf_fname, const std::string& stats_fname, const std::string& model_stage_dir, const std::string& output_fname, int n_threads) {
    stats_t stats = load_stats(stats_fname);
    const std::map<std::string, model_data_t> model_definitions = load_model_features(model_stage_dir);
    std::map<std::string, model_data_t> models = model_definitions;
    tbx_t* index = tbx_index_load(vcf_fname.c_str());
    if (index == NULL) throw std::runtime_error("Failed to load the index for " + vcf_fname + ".");
    std::vector<std::string> chroms;
    int nseq = 0;
    const char** seqnames = tbx_seqnames(index, &nseq);
    for (int rid = 0; rid < nseq; rid++) chroms.push_back(seqnames[rid]);
    free(seqnames);
    tbx_destroy(index);
    ctpl::thread_pool pool(n_threads);
    std::vector<std::future<std::map<std::string, model_data_t>>> futures;
    for (const std::string& chrom : chroms) futures.push_back(pool.push([&, chrom](int) { return process_chrom(vcf_fname, chrom, stats, model_definitions); }));
    for (std::future<std::map<std::string, model_data_t>>& future : futures) {
        std::map<std::string, model_data_t> chrom_models = future.get();
        for (auto& model_entry : chrom_models) {
            model_data_t& dest = models.at(model_entry.first);
            model_data_t& src = model_entry.second;
            dest.variant_ids.insert(dest.variant_ids.end(), src.variant_ids.begin(), src.variant_ids.end());
            dest.record_keys.insert(dest.record_keys.end(), src.record_keys.begin(), src.record_keys.end());
            dest.values.insert(dest.values.end(), src.values.begin(), src.values.end());
        }
    }

    std::ofstream out(output_fname, std::ios::binary);
    if (!out) throw std::runtime_error("Failed to create " + output_fname + ".");
    const char magic[8] = {'S', 'V', 'F', 'E', 'A', 'T', '2', '\0'};
    out.write(magic, sizeof(magic));
    uint32_t n_models = models.size();
    out.write((char*)&n_models, sizeof(n_models));
    for (const auto& model_entry : models) {
        const std::string& model_name = model_entry.first;
        const model_data_t& model = model_entry.second;
        uint32_t name_len = model_name.size(), n_features = model.feature_names.size();
        uint64_t n_variants = model.variant_ids.size();
        out.write((char*)&name_len, sizeof(name_len)); out.write(model_name.data(), name_len);
        out.write((char*)&n_variants, sizeof(n_variants)); out.write((char*)&n_features, sizeof(n_features));
        if (n_variants > 0) out.write((char*)model.variant_ids.data(), n_variants*sizeof(uint64_t));
        if (n_variants > 0) out.write((char*)model.record_keys.data(), n_variants*sizeof(uint64_t));
        if (!model.values.empty()) out.write((char*)model.values.data(), model.values.size()*sizeof(double));
    }
    if (!out) throw std::runtime_error("Failed while writing " + output_fname + ".");
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: extract_features input.vcf.gz stats.txt model_stage_dir output.bin threads\n";
        return 1;
    }
    try {
        int n_threads = std::stoi(argv[5]);
        if (n_threads < 1) throw std::runtime_error("The number of threads must be positive.");
        extract_features(argv[1], argv[2], argv[3], argv[4], n_threads);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
