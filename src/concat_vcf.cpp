#include <htslib/vcf.h>
#include <htslib/tbx.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>

bcf_hdr_t* read_bcf(std::string vcf_fname, std::vector<bcf1_t*>& svs) {
    htsFile* vcf = bcf_open(vcf_fname.c_str(), "r");
    if (vcf == NULL) {
        throw std::runtime_error("Unable to open file " + vcf_fname + ".");
    }
    bcf_hdr_t* hdr = bcf_hdr_read(vcf);
    if (hdr == NULL) {
        throw std::runtime_error("Unable to read header from file " + vcf_fname + ".");
    }

    bcf1_t* bcf = bcf_init();
    while (bcf_read(vcf, hdr, bcf) == 0) {
        svs.push_back(bcf_dup(bcf));
    }

    bcf_close(vcf);
    bcf_destroy(bcf);
    return hdr;
}

int main(int argc, char* argv[]) {
    std::vector<bcf1_t*> svs;
    bcf_hdr_t* hdr = read_bcf(argv[1], svs);
    read_bcf(argv[2], svs);

    std::stable_sort(svs.begin(), svs.end(), [](bcf1_t* a, bcf1_t* b) {
        return std::tie(a->rid, a->pos) < std::tie(b->rid, b->pos);
    });

    std::string out_fname = argv[3];
    htsFile* out_vcf = bcf_open(out_fname.c_str(), "wz");
    if (out_vcf == NULL) {
        throw std::runtime_error("Unable to open file " + out_fname + ".");
    }

    if (bcf_hdr_write(out_vcf, hdr) != 0) {
        throw std::runtime_error("Unable to write header to file " + out_fname + ".");
    }
    for (bcf1_t* bcf : svs) {
        if (bcf_write(out_vcf, hdr, bcf) != 0) {
            throw std::runtime_error("Unable to write record to file " + out_fname + ".");
        }
    }

    bcf_close(out_vcf);
    bcf_hdr_destroy(hdr);
    for (bcf1_t* bcf : svs) {
        bcf_destroy(bcf);
    }

    tbx_index_build(out_fname.c_str(), 0, &tbx_conf_vcf);
}
