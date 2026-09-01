#include <string>
#include <htslib/sam.h>
#include <htslib/tbx.h>
#include <htslib/vcf.h>
#include <unordered_map>

#include "htslib/hts.h"
#include "vcf_utils.h"
#include "SegTree.h"

const SegTree::T max_alt_alleles = 2; // maximum number of alt alleles allowed

int main(int argc, char** argv) {
    
    std::string in_vcf_fname = argv[1];
    std::string out_vcf_fname = argv[2];

    htsFile* in_vcf_file = bcf_open(in_vcf_fname.c_str(), "r");
    if (!in_vcf_file) {
        throw std::runtime_error("Failed to open input VCF file: " + in_vcf_fname);
    }
    bcf_hdr_t* hdr = bcf_hdr_read(in_vcf_file);
    if (!hdr) {
        bcf_close(in_vcf_file);
        throw std::runtime_error("Failed to read header from VCF file: " + in_vcf_fname);
    }

    // build interval trees for het and hom SVs
    bcf1_t* b = bcf_init();
    std::vector<bcf1_t*> all_svs;
    std::unordered_map<std::string, std::vector<std::pair<bcf1_t*, double>>> svs_with_epr;
    std::unordered_map<std::string, hts_pos_t> chr_lengths;
    while (bcf_read(in_vcf_file, hdr, b) == 0) {
        bcf1_t* sv = bcf_dup(b);
        all_svs.push_back(sv);

        if (get_sv_type(hdr, sv) != "DEL" && get_sv_type(hdr, sv) != "INS") continue;

        std::string chr = bcf_seqname_safe(hdr, sv);
        float* epr = NULL;
        int nvalues = 0, len = 0;
        bcf_get_format_float(hdr, sv, "EPR", &epr, &nvalues);
        if (nvalues <= 0) continue;

        if (count_alt_alleles(hdr, sv) > 0) {
            svs_with_epr[chr].emplace_back(sv, epr[0]);
            chr_lengths[chr] = std::max(chr_lengths[chr], (hts_pos_t) get_sv_end(hdr, sv));
        }
        free(epr);
    }

    // mark as 0/0 deletions that overlap a better (i.e., higher EPR) 1/1 deletion
    for (auto& p1 : svs_with_epr) {
        const std::string& chr = p1.first;
        std::vector<std::pair<bcf1_t*, double>>& svs = p1.second;

        // sort svs by descending epr
        std::sort(svs.begin(), svs.end(), [](const std::pair<bcf1_t*, double>& a, const std::pair<bcf1_t*, double>& b) {
            return a.second > b.second;
        });

        SegTree seg_tree_del(chr_lengths[chr]), seg_tree_ins(chr_lengths[chr]);
        for (const auto& p2 : svs) {
            bcf1_t* sv = p2.first;
            float epr = p2.second;

            std::string svtype = get_sv_type(hdr, sv);
            SegTree& seg_tree = (svtype == "DEL") ? seg_tree_del : seg_tree_ins;
            hts_pos_t range_start = sv->pos, range_end = get_sv_end(hdr, sv);
            if (svtype == "INS") {
                range_start -= 100;
                range_end += 100;
            }

            int n_alt_alleles = count_alt_alleles(hdr, sv);
            int max_alt_alleles = 2 - n_alt_alleles;
            if (seg_tree.any_ge(range_start, range_end, max_alt_alleles+1)) {
                int* gt = (int*)malloc(2 * sizeof(int));
                gt[0] = bcf_gt_unphased(0);
                gt[1] = bcf_gt_unphased(0);
                bcf_update_genotypes(hdr, sv, gt, 2);
                free(gt);
            } else {
                // update segment tree with the new deletion
                seg_tree.add(range_start, range_end, n_alt_alleles);
            }
        }
    }

    htsFile* out_vcf_file = bcf_open(out_vcf_fname.c_str(), "wz");
    if (!out_vcf_file) {
        bcf_hdr_destroy(hdr);
        bcf_close(in_vcf_file);
        throw std::runtime_error("Failed to open output VCF file: " + out_vcf_fname);
    }
    if (bcf_hdr_write(out_vcf_file, hdr) != 0) {
        bcf_hdr_destroy(hdr);
        bcf_close(in_vcf_file);
        bcf_close(out_vcf_file);
        throw std::runtime_error("Failed to write the VCF header to " + out_vcf_fname + ".");
    }
    for (bcf1_t* sv : all_svs) {
        if (bcf_write1(out_vcf_file, hdr, sv) != 0) {
            bcf_hdr_destroy(hdr);
            bcf_close(in_vcf_file);
            bcf_close(out_vcf_file);
            throw std::runtime_error("Failed to write SV to " + out_vcf_fname + ".");
        }
    }
    bcf_hdr_destroy(hdr);
    bcf_close(in_vcf_file);
    bcf_close(out_vcf_file);
    for (bcf1_t* sv : all_svs) {
        bcf_destroy(sv);
    }
}
