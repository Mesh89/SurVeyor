#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <htslib/vcf.h>

struct prediction_t {
    int8_t gt;
    bool force_missing_ref;
    float epr, hopr, expr;
    int8_t erefa;
};

template<class T> void read_value(std::ifstream& in, T& value) {
    in.read((char*)&value, sizeof(value));
    if (!in) throw std::runtime_error("Malformed prediction bundle.");
}

std::unordered_map<uint64_t, prediction_t> load_predictions(const std::string& predictions_fname) {
    std::ifstream in(predictions_fname, std::ios::binary);
    if (!in) throw std::runtime_error("Failed to open " + predictions_fname + ".");
    char magic[8];
    in.read(magic, sizeof(magic));
    const char expected_magic[8] = {'S', 'V', 'P', 'R', 'E', 'D', '1', '\0'};
    if (!in || memcmp(magic, expected_magic, sizeof(magic)) != 0) throw std::runtime_error("Invalid prediction bundle.");
    uint64_t n_predictions;
    read_value(in, n_predictions);
    std::unordered_map<uint64_t, prediction_t> predictions;
    predictions.reserve(n_predictions);
    for (uint64_t i = 0; i < n_predictions; i++) {
        uint64_t id;
        int8_t force_missing_ref;
        prediction_t prediction;
        read_value(in, id); read_value(in, prediction.gt); read_value(in, force_missing_ref); prediction.force_missing_ref = force_missing_ref;
        read_value(in, prediction.epr); read_value(in, prediction.hopr); read_value(in, prediction.expr); read_value(in, prediction.erefa);
        predictions[id] = prediction;
    }
    return predictions;
}

void add_format_header(bcf_hdr_t* hdr, const char* id, const char* type, const char* description) {
    int header_id = bcf_hdr_id2int(hdr, BCF_DT_ID, id);
    if (header_id >= 0 && bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, header_id)) return;
    std::string line = "##FORMAT=<ID=" + std::string(id) + ",Number=1,Type=" + type + ",Description=\"" + description + "\">";
    if (bcf_hdr_append(hdr, line.c_str()) != 0) throw std::runtime_error("Failed to add FORMAT/" + std::string(id) + " to the VCF header.");
}

void update_header(bcf_hdr_t* hdr, const std::string& training_set_sha256) {
    add_format_header(hdr, "EPR", "Float", "Probability of the SV existing in the sample, according to the ML model.");
    add_format_header(hdr, "HOPR", "Float", "Probability of an existing SV to be homozygous, according to the ML model.");
    add_format_header(hdr, "EXPR", "Float", "Probability of the SV to be represented exactly, according to the ML model.");
    add_format_header(hdr, "EREFA", "Integer", "Whether the EREFA-stage classifier requires the other allele to be reference.");
    if (training_set_sha256 != ".") {
        bcf_hdr_remove(hdr, BCF_HL_GEN, "SurVeyorTrainingSetSHA256");
        std::string line = "##SurVeyorTrainingSetSHA256=" + training_set_sha256;
        if (bcf_hdr_append(hdr, line.c_str()) != 0) throw std::runtime_error("Failed to add the training-set hash to the VCF header.");
    }
    if (bcf_hdr_sync(hdr) != 0) throw std::runtime_error("Failed to synchronize the VCF header.");
}

void set_genotype(bcf_hdr_t* hdr, bcf1_t* record, int gt, bool force_missing_ref) {
    int32_t values[2];
    values[0] = force_missing_ref ? bcf_gt_missing : bcf_gt_unphased(gt/2);
    values[1] = bcf_gt_unphased(gt >= 1 ? 1 : 0);
    if (bcf_update_genotypes(hdr, record, values, 2) != 0) throw std::runtime_error("Failed to update FORMAT/GT.");
}

void set_format_float(bcf_hdr_t* hdr, bcf1_t* record, const char* key, float value) {
    if (std::isnan(value)) bcf_float_set_missing(value);
    if (bcf_update_format_float(hdr, record, key, &value, 1) != 0) throw std::runtime_error("Failed to update FORMAT/" + std::string(key) + ".");
}

void set_format_int32(bcf_hdr_t* hdr, bcf1_t* record, const char* key, int value) {
    int32_t int_value = value < 0 ? bcf_int32_missing : value;
    if (bcf_update_format_int32(hdr, record, key, &int_value, 1) != 0) throw std::runtime_error("Failed to update FORMAT/" + std::string(key) + ".");
}

void clear_ml_fields(bcf_hdr_t* hdr, bcf1_t* record) {
    if (bcf_update_filter(hdr, record, NULL, 0) != 0) throw std::runtime_error("Failed to clear the VCF filters.");
    set_format_float(hdr, record, "EPR", NAN);
    set_format_float(hdr, record, "HOPR", NAN);
    set_format_float(hdr, record, "EXPR", NAN);
    set_format_int32(hdr, record, "EREFA", -1);
}

void apply_prediction(bcf_hdr_t* hdr, bcf1_t* record, const prediction_t& prediction) {
    set_genotype(hdr, record, prediction.gt, prediction.force_missing_ref);
    set_format_float(hdr, record, "EPR", prediction.epr);
    set_format_float(hdr, record, "HOPR", prediction.hopr);
    set_format_float(hdr, record, "EXPR", prediction.expr);
    set_format_int32(hdr, record, "EREFA", prediction.erefa);
}

void write_classifier_vcf(const std::string& input_fname, const std::string& predictions_fname, const std::string& output_fname, const std::string& training_set_sha256, int n_threads) {
    std::unordered_map<uint64_t, prediction_t> predictions = load_predictions(predictions_fname);
    
    htsFile* in = bcf_open(input_fname.c_str(), "r");
    if (in == NULL) throw std::runtime_error("Failed to open " + input_fname + ".");
    if (hts_set_threads(in, n_threads) != 0) { hts_close(in); throw std::runtime_error("Failed to enable multithreaded decompression for " + input_fname + "."); }
    
    bcf_hdr_t* hdr = bcf_hdr_read(in);
    if (hdr == NULL) { hts_close(in); throw std::runtime_error("Failed to read the header from " + input_fname + "."); }
    if (bcf_hdr_nsamples(hdr) != 1) { bcf_hdr_destroy(hdr); hts_close(in); throw std::runtime_error("The classifier VCF writer requires exactly one sample."); }
    update_header(hdr, training_set_sha256);
    
    htsFile* out = bcf_open(output_fname.c_str(), "wz");
    if (out == NULL) { bcf_hdr_destroy(hdr); hts_close(in); throw std::runtime_error("Failed to create " + output_fname + "."); }
    if (hts_set_threads(out, n_threads) != 0) { hts_close(out); bcf_hdr_destroy(hdr); hts_close(in); throw std::runtime_error("Failed to enable multithreaded compression for " + output_fname + "."); }
    if (bcf_hdr_write(out, hdr) != 0) { hts_close(out); bcf_hdr_destroy(hdr); hts_close(in); throw std::runtime_error("Failed to write the VCF header."); }
    
    int32_t current_rid = -1;
    uint64_t record_idx = 0;
    bcf1_t* record = bcf_init();
    if (record == NULL) { hts_close(out); bcf_hdr_destroy(hdr); hts_close(in); throw std::runtime_error("Failed to allocate a VCF record."); }
    while (bcf_read(in, hdr, record) == 0) {
        if (record->rid != current_rid) { current_rid = record->rid; record_idx = 0; }
        uint64_t record_key = (uint64_t(uint32_t(record->rid))<<32)|record_idx++;
        auto prediction = predictions.find(record_key);
        if (prediction != predictions.end()) { clear_ml_fields(hdr, record); apply_prediction(hdr, record, prediction->second); }
        if (bcf_write(out, hdr, record) != 0) { bcf_destroy(record); hts_close(out); bcf_hdr_destroy(hdr); hts_close(in); throw std::runtime_error("Failed to write a VCF record."); }
    }
    bcf_destroy(record);
    int in_close_result = hts_close(in), out_close_result = hts_close(out);
    bcf_hdr_destroy(hdr);

    if (in_close_result != 0) throw std::runtime_error("Failed to finalize reading " + input_fname + ".");
    if (out_close_result != 0) throw std::runtime_error("Failed to finalize " + output_fname + ".");
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: write_classifier_vcf input.vcf.gz predictions.bin output.vcf.gz training_set_sha256 threads\n";
        return 1;
    }
    try {
        int n_threads = std::stoi(argv[5]);
        if (n_threads < 1) throw std::runtime_error("The number of threads must be positive.");
        write_classifier_vcf(argv[1], argv[2], argv[3], argv[4], n_threads);
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    return 0;
}
