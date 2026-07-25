#ifndef VCF_UTILS_H
#define VCF_UTILS_H

#include <cstddef>
#include <iostream>
#include <chrono>
#include <ctime>
#include <htslib/vcf.h>
#include <memory>
#include <sstream>
#include "htslib/hts.h"
#include "types.h"
#include "utils.h"

bcf_hrec_t* generate_contig_hrec() {
	bcf_hrec_t* contig_hrec = (bcf_hrec_t*) malloc(sizeof(bcf_hrec_t));
	contig_hrec->type = BCF_HL_CTG;
	contig_hrec->key = strdup("contig");
	contig_hrec->value = NULL;
	contig_hrec->keys = contig_hrec->vals = NULL;
	contig_hrec->nkeys = 0;
	int r1 = bcf_hrec_add_key(contig_hrec, "ID", 2);
	int r2 = bcf_hrec_add_key(contig_hrec, "length", 6);
	if (r1 || r2) {
		throw std::runtime_error("Failed to create contig to VCF header.");
	}
	return contig_hrec;
}

// Define a struct to hold all properties for each format tag
typedef struct {
    const char* suffix;       // Tag suffix (e.g., "C", "CF", etc.)
	int n;					  // Number of values	
    const char* type;         // Data type ("Integer" or "Float")
    const char* desc_format;  // Description format string
} format_tag_def_t;

void add_headers(bcf_hdr_t* hdr, const char prefix, int bp_n, const char* bp_name, const char* allele_type, const format_tag_def_t formats[], int n_formats) {
	char tag_id[20];
    char tag_str[1024];
    char desc_str[200];
    int len;

    for (size_t i = 0; i < n_formats; i++) {
        const format_tag_def_t* fmt = &formats[i];

        snprintf(tag_id, sizeof(tag_id), fmt->suffix, prefix, bp_n, fmt->suffix);
        snprintf(desc_str, sizeof(desc_str), fmt->desc_format, bp_name, allele_type);
        snprintf(tag_str, sizeof(tag_str), 
                "##FORMAT=<ID=%s,Number=%d,Type=%s,Description=\"%s\">",
                tag_id, fmt->n, fmt->type, desc_str);

        bcf_hdr_remove(hdr, BCF_HL_FMT, tag_id);
        bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, tag_str, &len));
    }
}

void add_read_support_headers(bcf_hdr_t* hdr, const char prefix, int bp_n, const char* bp_name, const char* allele_type) {
    // Array of tag definitions
    const format_tag_def_t formats[] = {
        {
            "%cR%d", 1, "Integer",
            prefix == 'R'
                ? "Number of reads supporting %s in the %s allele, including reads that also support another variant."
                : "Number of reads supporting %s in the %s allele."
        },
        {
            "%cR%dC", 1, "Integer", 
            "Number of consistent reads supporting %s in the %s allele."
        },
        {
            "%cR%dCF", 1, "Integer",
            "Number of consistent forward reads supporting %s in the %s allele."
        },
        {
            "%cR%dCR", 1, "Integer",
            "Number of consistent reverse reads supporting %s in the %s allele."
        },
        {
            "%cR%dE", 1, "Integer",
            "Number of consistent exact (i.e., perfectly matching) reads supporting %s in the %s allele."
        },
        {
            "%cR%dEHQ", 1, "Integer",
            "Number of high-quality consistent exact (i.e., perfectly matching) reads supporting %s in the %s allele."
        },
        {
            "%cR%dEmQ", 1, "Integer",
            "Minimum mate mapping quality of consistent exact (i.e., perfectly matching) reads supporting %s in the %s allele."
        },
        {
            "%cR%dEMQ", 1, "Integer",
            "Maximum mate mapping quality of consistent exact (i.e., perfectly matching) reads supporting %s in the %s allele."
        },
        {
            "%cR%dEF", 1, "Integer",
            "Number of consistent exact (i.e., perfectly matching) forward reads supporting %s in the %s allele."
        },
        {
            "%cR%dER", 1, "Integer",
            "Number of consistent exact (i.e., perfectly matching) reverse reads supporting %s in the %s allele."
        },
        {
            "%cR%dCAS", 1, "Float",
            "Average aln score of consistent reads supporting %s of the SV to the %s allele consensus."
        },
		{
			"%cR%dCSS", 1, "Float",
			"Standard deviation of aln score of consistent reads supporting %s of the SV to the %s allele consensus."
		},
        {
            "%cR%dCHQ", 1, "Integer",
            "Number of high-quality consistent reads supporting %s in the %s allele."
        },
		{
            "%cR%dCmQ", 1, "Integer",
            "Minimum mate mapping quality of consistent reads supporting %s in the %s allele."
        },
        {
            "%cR%dCMQ", 1, "Integer",
            "Maximum mate mapping quality of consistent reads supporting %s in the %s allele."
        },
		{
            "%cR%dCAQ", 1, "Float",
            "Average mate mapping quality of consistent reads supporting %s in the %s allele."
        },
		{
			"%cR%dCSQ", 1, "Float",
			"Standard deviation of mate mapping quality of consistent reads supporting %s in the %s allele."
		},
		{
			"%cR%dCMSPAN", 2, "Integer",
			"Number of base pairs to the left and right, respectively, of the %s in the %s allele covered by mates of consistent supporting reads it."
		},
		{
			"%cR%dCMHQSPAN", 2, "Integer",
			"Number of base pairs to the left and right, respectively, of the %s in the %s allele covered by high-quality mates of consistent supporting reads it."
		}
    };

    add_headers(hdr, prefix, bp_n, bp_name, allele_type, formats, sizeof(formats)/sizeof(formats[0]));
}

void add_pairs_support_headers(bcf_hdr_t* hdr, const char prefix, int bp_n, const char* bp_name, const char* allele_type) {
	const format_tag_def_t formats[] = {
		{
			"%cSP%d", 1, "Integer",
			"Number of pairs supporting %s in the %s allele."
		},
		{
			"%cSP%dHQ", 2, "Integer",
			"Number of high-quality positive and negative reads, respectively, within pairs supporting %s in the %s allele."
		},
		{
			"%cSP%dmQ", 2, "Integer",
			"Minimum mapping quality of positive and negative reads, respectively, within pairs supporting %s in the %s allele."
		},
		{
			"%cSP%dMQ", 2, "Integer",
			"Maximum mapping quality of positive and negative reads, respectively, within pairs supporting %s in the %s allele."
		},
		{
			"%cSP%dAQ", 2, "Float",
			"Average mapping quality of positive and negative reads, respectively, within pairs supporting %s in the %s allele."
		},
		{
			"%cSP%dSQ", 2, "Float",
			"Standard deviation of mapping quality of positive and negative reads, respectively, within pairs supporting %s in the %s allele."
		},
		{
			"%cSP%dSPAN", 2, "Integer",
			"Number of base pairs to the left and right, respectively, of the %s in the %s allele covered by reads within its supporting pairs."
		},
		{
			"%cSP%dNMA", 2, "Float",
			"Average NM (non-matches) of positive and negative reads, respectively, within pairs supporting %s in the %s allele."
		},
		{
			"%cSP%dNMS", 2, "Float",
			"Standard deviation of NM (non-matches) of positive and negative reads, respectively, within pairs supporting %s in the %s allele."
		}
	};

	add_headers(hdr, prefix, bp_n, bp_name, allele_type, formats, sizeof(formats)/sizeof(formats[0]));
}

void add_stray_pairs_headers(bcf_hdr_t* hdr, int bp_number) {
	const format_tag_def_t formats[] = {
		{
			"%cSP%d", 1, "Integer",
			"Number of pairs around the %s in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dHQ", 2, "Integer",
			"Number of high-quality pairs around the %s in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dmQ", 2, "Integer",
			"Minimum mapping quality of pairs around the %s in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dMQ", 2, "Integer",
			"Maximum mapping quality of pairs around the %s in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dAQ", 2, "Float",
			"Average mapping quality of pairs around the %s in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dSQ", 2, "Float",
			"Standard deviation of mapping quality of pairs around the %s in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dNMA", 2, "Float",
			"Average NM (non-matches) of pairs around the %s in the %s allele that are discordant and yet do not support the SV."
		},
		{
			"%cSP%dNMS", 2, "Float",
			"Standard deviation of NM (non-matches) of pairs around the %s in the %s allele that are discordant and yet do not support the SV."
		}
	};

	std::string bp_name = "breakpoint " + std::to_string(bp_number);
	add_headers(hdr, 'S', bp_number, bp_name.c_str(), "reference", formats, sizeof(formats)/sizeof(formats[0]));
}

void add_fmt_tags(bcf_hdr_t* hdr) {
	int len = 0;

	bcf_hdr_remove(hdr, BCF_HL_INFO, "END");
	const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, end_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "SVLEN");
	const char* svlen_tag = "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, svlen_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "SVINSLEN");
	const char* svinslen_tag = "##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description=\"Length of the inserted sequence.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, svinslen_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "HPID");
	const char* hpid_tag = "##INFO=<ID=HPID,Number=1,Type=Integer,Description=\"Identifier of the local haplotype represented by this call.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, hpid_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "DUP_CN_UNRESOLVED");
	const char* dup_cn_unresolved_tag = "##INFO=<ID=DUP_CN_UNRESOLVED,Number=0,Type=Flag,Description=\"For this DUP record, the number of additional tandem copies is unresolved; END and SVLEN describe the duplicated reference interval, not necessarily the total inserted allele length.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, dup_cn_unresolved_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "LEFT_ANCHOR_BASE_COUNT");
	const char* labc_tag = "##INFO=<ID=LEFT_ANCHOR_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the left anchor region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,labc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "LEFT_FLANKING_BASE_COUNT_50");
	const char* labc50_tag = "##INFO=<ID=LEFT_FLANKING_BASE_COUNT_50,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the first 50 bp of the left flanking region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,labc50_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "LEFT_FLANKING_BASE_COUNT_100");
	const char* labc100_tag = "##INFO=<ID=LEFT_FLANKING_BASE_COUNT_100,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the first 100 bp of the left flanking region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,labc100_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "LEFT_FLANKING_BASE_COUNT_500");
	const char* labc500_tag = "##INFO=<ID=LEFT_FLANKING_BASE_COUNT_500,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the first 500 bp of the left flanking region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,labc500_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "RIGHT_ANCHOR_BASE_COUNT");
	const char* rabc_tag = "##INFO=<ID=RIGHT_ANCHOR_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the right anchor region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rabc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "RIGHT_FLANKING_BASE_COUNT_50");
	const char* rabc50_tag = "##INFO=<ID=RIGHT_FLANKING_BASE_COUNT_50,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the first 50 bp of the right flanking region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rabc50_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "RIGHT_FLANKING_BASE_COUNT_100");
	const char* rabc100_tag = "##INFO=<ID=RIGHT_FLANKING_BASE_COUNT_100,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the first 100 bp of the right flanking region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rabc100_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "RIGHT_FLANKING_BASE_COUNT_500");
	const char* rabc500_tag = "##INFO=<ID=RIGHT_FLANKING_BASE_COUNT_500,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the first 500 bp of the right flanking region of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,rabc500_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "SV_REF_PREFIX_BASE_COUNT");
	const char* svrefpbc_tag = "##INFO=<ID=SV_REF_PREFIX_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the first 5000 bp of the reference region of the SV (if the SV affects less than 5000 bp, the whole SV is considered).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,svrefpbc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "SV_REF_SUFFIX_BASE_COUNT");
	const char* svrefsbc_tag = "##INFO=<ID=SV_REF_SUFFIX_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the last 5000 bp of the reference region of the SV (if the SV affects less than 5000 bp, the whole SV is considered).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,svrefsbc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "INS_PREFIX_BASE_COUNT");
	const char* pbc_tag = "##INFO=<ID=INS_PREFIX_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the prefix of the inserted sequence (for incomplete assemblies) or in the full inserted sequence.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,pbc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "INS_SUFFIX_BASE_COUNT");
	const char* sbc_tag = "##INFO=<ID=INS_SUFFIX_BASE_COUNT,Number=4,Type=Integer,Description=\"Number of As, Cs, Gs, and Ts in the suffix of the inserted sequence (for incomplete assemblies) or in the full inserted sequence. For insertion not marked with INCOMPLETE_ASSEMBLY, this will be identical to PREFIX_BASE_COUNT.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,sbc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "MH_LEN");
	const char* mh_len = "##INFO=<ID=MH_LEN,Number=1,Type=Integer,Description=\"Length of the microhomology.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr,mh_len, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "HP_GENOTYPED");
	const char* hp_genotyped_tag = "##INFO=<ID=HP_GENOTYPED,Number=0,Type=Flag,Description=\"This variant was genotyped using the homopolymer-specific genotyping path.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, hp_genotyped_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_INFO, "HP_REF_RANGE");
	const char* hp_ref_range_tag = "##INFO=<ID=HP_REF_RANGE,Number=2,Type=Integer,Description=\"1-based half-open [start,end) coordinates of the reference homopolymer run used in homopolymer-specific genotyping.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, hp_ref_range_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "GT");
	const char* gt_tag = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, gt_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "FT");
	const char* ft_tag = "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Filter. PASS indicates a reliable call. Any other value means the call is not reliable.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ft_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "MD");
	const char* md_tag = "##FORMAT=<ID=MD,Number=4,Type=Integer,Description=\"Median depth in the left flanking, prefix, suffix and right flanking regions of the SV.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, md_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "MDHQ");
	const char* mdhq_tag = "##FORMAT=<ID=MDHQ,Number=4,Type=Integer,Description=\"Median depth in the left flanking, prefix, suffix and right flanking regions of the SV. Only high MAPQ reads are considered.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, mdhq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "CLMD");
	const char* clmd_tag = "##FORMAT=<ID=CLMD,Number=2,Type=Integer,Description=\"Median depth in the left and right flanking regions, considering only the regions where the SV evidence (DP and SR) lies.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, clmd_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "CLMDHQ");
	const char* clmdhq_tag = "##FORMAT=<ID=CLMDHQ,Number=2,Type=Integer,Description=\"Median depth in the left and right flanking regions, considering only the regions where the SV evidence (DP and SR) lies. Only high MAPQ reads are considered.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, clmdhq_tag, &len));

	std::string bp1_name = "breakpoint 1";
	std::string bp2_name = "breakpoint 2";
	add_read_support_headers(hdr, 'A', 1, bp1_name.c_str(), "alternate");
	add_pairs_support_headers(hdr, 'A', 1, bp1_name.c_str(), "alternate");

	add_read_support_headers(hdr, 'A', 2, bp2_name.c_str(), "alternate");
	add_pairs_support_headers(hdr, 'A', 2, bp2_name.c_str(), "alternate");

	add_read_support_headers(hdr, 'R', 1, bp1_name.c_str(), "reference");
	add_pairs_support_headers(hdr, 'R', 1, bp1_name.c_str(), "reference");

	add_read_support_headers(hdr, 'R', 2, bp2_name.c_str(), "reference");
	add_pairs_support_headers(hdr, 'R', 2, bp2_name.c_str(), "reference");

	add_pairs_support_headers(hdr, 'N', 1, bp1_name.c_str(), "neutral");
	add_pairs_support_headers(hdr, 'N', 2, bp2_name.c_str(), "neutral");

	add_stray_pairs_headers(hdr, 1);
	add_stray_pairs_headers(hdr, 2);

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ER");
	const char* er_tag = "##FORMAT=<ID=ER,Number=1,Type=Integer,Description=\"Number of reads supporting equally well reference and alternate allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, er_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ERHQ");
	const char* erhq_tag = "##FORMAT=<ID=ERHQ,Number=1,Type=Integer,Description=\"Number of high-quality reads supporting equally well reference and alternate allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, erhq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR1");
	const char* oar1_tag = "##FORMAT=<ID=OAR1,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 1 of the ALT allele of this SV, but assigned to another variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar1_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR1HPID");
	const char* oar1hpid_tag = "##FORMAT=<ID=OAR1HPID,Number=1,Type=Integer,Description=\"HPID receiving the most reads counted in OAR1. Ties are resolved in favor of the smallest HPID.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar1hpid_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR1MAX");
	const char* oar1max_tag = "##FORMAT=<ID=OAR1MAX,Number=1,Type=Integer,Description=\"Number of OAR1 reads assigned to the HPID reported in OAR1HPID.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar1max_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR1C");
	const char* oar1c_tag = "##FORMAT=<ID=OAR1C,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 1 of the ALT allele of this SV that are assigned to another variant and are consistent with that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar1c_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR1CHQ");
	const char* oar1chq_tag = "##FORMAT=<ID=OAR1CHQ,Number=1,Type=Integer,Description=\"Number of high-quality reads supporting breakpoint 1 of the ALT allele of this SV that are assigned to another variant and are consistent with that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar1chq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR1E");
	const char* oar1e_tag = "##FORMAT=<ID=OAR1E,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 1 of the ALT allele of this SV that are assigned to another variant and are exact, consistent reads for that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar1e_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR2");
	const char* oar2_tag = "##FORMAT=<ID=OAR2,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 2 of the ALT allele of this SV, but assigned to another variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR2HPID");
	const char* oar2hpid_tag = "##FORMAT=<ID=OAR2HPID,Number=1,Type=Integer,Description=\"HPID receiving the most reads counted in OAR2. Ties are resolved in favor of the smallest HPID.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar2hpid_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR2MAX");
	const char* oar2max_tag = "##FORMAT=<ID=OAR2MAX,Number=1,Type=Integer,Description=\"Number of OAR2 reads assigned to the HPID reported in OAR2HPID.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar2max_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR2C");
	const char* oar2c_tag = "##FORMAT=<ID=OAR2C,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 2 of the ALT allele of this SV that are assigned to another variant and are consistent with that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar2c_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR2CHQ");
	const char* oar2chq_tag = "##FORMAT=<ID=OAR2CHQ,Number=1,Type=Integer,Description=\"Number of high-quality reads supporting breakpoint 2 of the ALT allele of this SV that are assigned to another variant and are consistent with that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar2chq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "OAR2E");
	const char* oar2e_tag = "##FORMAT=<ID=OAR2E,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 2 of the ALT allele of this SV that are assigned to another variant and are exact, consistent reads for that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, oar2e_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ORR1");
	const char* orr1_tag = "##FORMAT=<ID=ORR1,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 1 of the REF allele of this SV that also support another variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, orr1_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ORR1C");
	const char* orr1c_tag = "##FORMAT=<ID=ORR1C,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 1 of the REF allele of this SV that also support another variant and are consistent with that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, orr1c_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ORR1CHQ");
	const char* orr1chq_tag = "##FORMAT=<ID=ORR1CHQ,Number=1,Type=Integer,Description=\"Number of high-quality reads supporting breakpoint 1 of the REF allele of this SV that also support another variant and are consistent with that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, orr1chq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ORR1E");
	const char* orr1e_tag = "##FORMAT=<ID=ORR1E,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 1 of the REF allele of this SV that also support another variant and are exact, consistent reads for that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, orr1e_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ORR2");
	const char* orr2_tag = "##FORMAT=<ID=ORR2,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 2 of the REF allele of this SV that also support another variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, orr2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ORR2C");
	const char* orr2c_tag = "##FORMAT=<ID=ORR2C,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 2 of the REF allele of this SV that also support another variant and are consistent with that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, orr2c_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ORR2CHQ");
	const char* orr2chq_tag = "##FORMAT=<ID=ORR2CHQ,Number=1,Type=Integer,Description=\"Number of high-quality reads supporting breakpoint 2 of the REF allele of this SV that also support another variant and are consistent with that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, orr2chq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "ORR2E");
	const char* orr2e_tag = "##FORMAT=<ID=ORR2E,Number=1,Type=Integer,Description=\"Number of reads supporting breakpoint 2 of the REF allele of this SV that also support another variant and are exact, consistent reads for that variant.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, orr2e_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "TD");
	const char* nc_tag = "##FORMAT=<ID=TD,Number=1,Type=Integer,Description=\"The variant region is too deep to be genotyped reliably.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, nc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "MINSIZE");
	const char* minsize_tag = "##FORMAT=<ID=MINSIZE,Number=1,Type=Integer,Description=\"Minimum size of the event calculated based on insert size distribution."
			"Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, minsize_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "MAXSIZE");
	const char* maxsize_tag = "##FORMAT=<ID=MAXSIZE,Number=1,Type=Integer,Description=\"Maximum size of the event calculated based on insert size distribution."
			"Note that this is calculated on the assumption of HOM_ALT events, and should be doubled to accommodate HET events.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, maxsize_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "KSPVAL");
	const char* kspval_tag = "##FORMAT=<ID=KSPVAL,Number=1,Type=Float,Description=\"p-value of the KS test between the insert size distibution of pairs surrounding the putative SV "
			"and the global distribution (insert size of pairs randomly sampled). If the two are significantly different, the likelihood of the SV being real is increased. \">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, kspval_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AXR");
	const char* axr_tag = "##FORMAT=<ID=AXR,Number=2,Type=Integer,Description=\"Number of reads used to extend the alternative allele consensus to the left and the right, respectively.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, axr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AXRHQ");
	const char* axrhq_tag = "##FORMAT=<ID=AXRHQ,Number=2,Type=Integer,Description=\"Number of high-quality reads used to extend the alternative allele consensus to the left and the right, respectively.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, axrhq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXL");
	const char* exl_tag = "##FORMAT=<ID=EXL,Number=1,Type=Integer,Description=\"Length of the extended alternative allele consensus.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exl_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXL2");
	const char* exl2_tag = "##FORMAT=<ID=EXL2,Number=1,Type=Integer,Description=\"Length of the extended alternative allele consensus for the second breakpoint (only for insertions).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exl2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXAS");
	const char* exas_tag = "##FORMAT=<ID=EXAS,Number=1,Type=Integer,Description=\"Edit distance between the extended alternative allele consensus and the original alternative allele (the reference with the variants applied).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exas_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXAS2");
	const char* exas2_tag = "##FORMAT=<ID=EXAS2,Number=1,Type=Integer,Description=\"Edit distance between the extended alternative allele consensus for the second breakpoint and the original alternative allele (the reference with the SV applied).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exas2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXRS");
	const char* exrs_tag = "##FORMAT=<ID=EXRS,Number=1,Type=Integer,Description=\"Score of the alignment between the extended alternative allele consensus and the reference.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exrs_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXRS2");
	const char* exrs2_tag = "##FORMAT=<ID=EXRS2,Number=1,Type=Integer,Description=\"Score of the alignment between the extended alternative allele consensus for the second breakpoint and the reference.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exrs2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSS");
	const char* exss_tag = "##FORMAT=<ID=EXSS,Number=2,Type=Integer,Description=\"How much of the extended alternative allele aligns to the left and right, respectively, of the breakpoint.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exss_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSS2");
	const char* exss2_tag = "##FORMAT=<ID=EXSS2,Number=2,Type=Integer,Description=\"How much of the extended alternative allele aligns to the left and right, respectively, of the second breakpoint (only for insertions).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exss2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSSC");
	const char* exssc_tag = "##FORMAT=<ID=EXSSC,Number=2,Type=Integer,Description=\"Score of the left half and right half of the extended alternative allele consensus.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exssc_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSSCIA");
	const char* exsscia_tag = "##FORMAT=<ID=EXSSCIA,Number=2,Type=Integer,Description=\"Score of the left half and right half of the extended alternative allele consensus for the first breakpoint, when allowed to map to the reference independently.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exsscia_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSSC2");
	const char* exssc2_tag = "##FORMAT=<ID=EXSSC2,Number=2,Type=Integer,Description=\"Score of the left half and right half of the extended alternative allele consensus for the second breakpoint (only for insertions).\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exssc2_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXSSC2IA");
	const char* exssc2ia_tag = "##FORMAT=<ID=EXSSC2IA,Number=2,Type=Integer,Description=\"Score of the left half and right half of the extended alternative allele consensus for the second breakpoint, when allowed to map to the reference independently.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, exssc2ia_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1C_OCCR");
	const char* ar1c_occr_tag = "##FORMAT=<ID=AR1C_OCCR,Number=1,Type=Float,Description=\"Occupancy ratio of starts of reads supporting breakpoint 1 of the ALT allele, normalised by expected occupancy ratio.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1c_occr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR2C_OCCR");
	const char* ar2c_occr_tag = "##FORMAT=<ID=AR2C_OCCR,Number=1,Type=Float,Description=\"Occupancy ratio of starts of reads supporting breakpoint 2 of the ALT allele, normalised by expected occupancy ratio.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar2c_occr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1HPMODE");
	const char* ar1m_tag = "##FORMAT=<ID=AR1HPMODE,Number=1,Type=Integer,Description=\"Mode of homopolymer length among all reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1m_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1CHPMODE");
	const char* ar1cm_tag = "##FORMAT=<ID=AR1CHPMODE,Number=1,Type=Integer,Description=\"Mode of homopolymer length among consistent reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1cm_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1CHPIQR");
	const char* ar1hpiqr_tag = "##FORMAT=<ID=AR1CHPIQR,Number=1,Type=Float,Description=\"Interquartile range of observed homopolymer lengths among good homopolymer reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1hpiqr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1CHPmQ");
	const char* ar1hpminmq_tag = "##FORMAT=<ID=AR1CHPmQ,Number=1,Type=Integer,Description=\"Minimum read mapping quality among good homopolymer reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1hpminmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1CHPMQ");
	const char* ar1hpmaxmq_tag = "##FORMAT=<ID=AR1CHPMQ,Number=1,Type=Integer,Description=\"Maximum read mapping quality among good homopolymer reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1hpmaxmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1CHPAQ");
	const char* ar1hpavgmq_tag = "##FORMAT=<ID=AR1CHPAQ,Number=1,Type=Float,Description=\"Average read mapping quality among good homopolymer reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1hpavgmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1CHPSQ");
	const char* ar1hpstdmq_tag = "##FORMAT=<ID=AR1CHPSQ,Number=1,Type=Float,Description=\"Standard deviation of read mapping quality among good homopolymer reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1hpstdmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1HP5PMR");
	const char* ar1hp5pmr_tag = "##FORMAT=<ID=AR1HP5PMR,Number=1,Type=Float,Description=\"Average 5' tail mismatch rate among homopolymer reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1hp5pmr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "AR1HP3PMR");
	const char* ar1hp3pmr_tag = "##FORMAT=<ID=AR1HP3PMR,Number=1,Type=Float,Description=\"Average 3' tail mismatch rate among homopolymer reads assigned to breakpoint 1 of the ALT allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, ar1hp3pmr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1HPMODE");
	const char* rr1m_tag = "##FORMAT=<ID=RR1HPMODE,Number=1,Type=Integer,Description=\"Mode of homopolymer length among all reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1m_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1CHPMODE");
	const char* rr1cm_tag = "##FORMAT=<ID=RR1CHPMODE,Number=1,Type=Integer,Description=\"Mode of homopolymer length among consistent reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1cm_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1CHPIQR");
	const char* rr1hpiqr_tag = "##FORMAT=<ID=RR1CHPIQR,Number=1,Type=Float,Description=\"Interquartile range of observed homopolymer lengths among good homopolymer reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1hpiqr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1CHPmQ");
	const char* rr1hpminmq_tag = "##FORMAT=<ID=RR1CHPmQ,Number=1,Type=Integer,Description=\"Minimum read mapping quality among good homopolymer reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1hpminmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1CHPMQ");
	const char* rr1hpmaxmq_tag = "##FORMAT=<ID=RR1CHPMQ,Number=1,Type=Integer,Description=\"Maximum read mapping quality among good homopolymer reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1hpmaxmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1CHPAQ");
	const char* rr1hpavgmq_tag = "##FORMAT=<ID=RR1CHPAQ,Number=1,Type=Float,Description=\"Average read mapping quality among good homopolymer reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1hpavgmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1CHPSQ");
	const char* rr1hpstdmq_tag = "##FORMAT=<ID=RR1CHPSQ,Number=1,Type=Float,Description=\"Standard deviation of read mapping quality among good homopolymer reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1hpstdmq_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1HP5PMR");
	const char* rr1hp5pmr_tag = "##FORMAT=<ID=RR1HP5PMR,Number=1,Type=Float,Description=\"Average 5' tail mismatch rate among homopolymer reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1hp5pmr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "RR1HP3PMR");
	const char* rr1hp3pmr_tag = "##FORMAT=<ID=RR1HP3PMR,Number=1,Type=Float,Description=\"Average 3' tail mismatch rate among homopolymer reads assigned to breakpoint 1 of the REF allele.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, rr1hp3pmr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EPR");
	const char* epr_tag = "##FORMAT=<ID=EPR,Number=1,Type=Float,Description=\"Probability of the SV existing in the sample, according to the ML model.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, epr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "HOPR");
	const char* hopr_tag = "##FORMAT=<ID=HOPR,Number=1,Type=Float,Description=\"Probability of the SV to be homozygous, according to the ML model.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, hopr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EXPR");
	const char* expr_tag = "##FORMAT=<ID=EXPR,Number=1,Type=Float,Description=\"Probability of the SV to be represented exactly, according to the ML model.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, expr_tag, &len));

	bcf_hdr_remove(hdr, BCF_HL_FMT, "EREFA");
	const char* erefa_tag = "##FORMAT=<ID=EREFA,Number=1,Type=Integer,Description=\"Whether the GT-stage classifier requires the other allele to be reference.\">";
	bcf_hdr_add_hrec(hdr, bcf_hdr_parse_line(hdr, erefa_tag, &len));
}

bcf_hdr_t* generate_vcf_header(chr_seqs_map_t& contigs, std::string sample_name, config_t config, std::string command) {
	bcf_hdr_t* header = bcf_hdr_init("w");

	// add contigs
	for (std::string contig_name : contigs.ordered_contigs) {
		bcf_hrec_t* hrec = generate_contig_hrec();
		int r1 = bcf_hrec_set_val(hrec, 0, contig_name.c_str(), contig_name.length(), false);
		std::string len_str = std::to_string(contigs.get_len(contig_name));
		int r2 = bcf_hrec_set_val(hrec, 1, len_str.c_str(), len_str.length(), false);
		if (r1 || r2) {
			throw std::runtime_error("Failed to create contig to VCF header.");
		}
		bcf_hdr_add_hrec(header, hrec);
	}

	int len;

	// add FILTER tags
	const char* remap_boundary_flt_tag = "##FILTER=<ID=REMAP_BOUNDARY_FILTER,Description=\"One of the breakpoints is incompatible with mate locations.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, remap_boundary_flt_tag, &len));

	const char* failed_to_ext_flt_tag = "##FILTER=<ID=FAILED_TO_EXTEND,Description=\"No reads can extend the consensus.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, failed_to_ext_flt_tag, &len));

	const char* anomalous_fl_depth_flt_tag = "##FILTER=<ID=ANOMALOUS_FLANKING_DEPTH,Description=\"The insertion region has anomalous depth.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_fl_depth_flt_tag, &len));

	const char* anomalous_sc_flt_tag = "##FILTER=<ID=ANOMALOUS_SC_NUMBER,Description=\"The number of soft-clipped reads supporting this call is too large.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, anomalous_sc_flt_tag, &len));

	const char* not_enough_disc_pairs_flt_tag = "##FILTER=<ID=NOT_ENOUGH_DISC_PAIRS,Description=\"Not enough discordant pairs supporting the event.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, not_enough_disc_pairs_flt_tag, &len));

	const char* low_mapq_disc_pairs_flt_tag = "##FILTER=<ID=LOW_MAPQ_DISC_PAIRS,Description=\"Discordant pairs supporting the SV have low MAPQ.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_mapq_disc_pairs_flt_tag, &len));

	const char* weak_support_flt_tag = "##FILTER=<ID=WEAK_SUPPORT,Description=\"Remapped breakpoint has low support from local reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, weak_support_flt_tag, &len));

	const char* low_ptn_ratio_flt_tag = "##FILTER=<ID=LOW_PTN_RATIO,Description=\"Low positive-to-negative ratio.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_ptn_ratio_flt_tag, &len));

	const char* short_anchor_flt_tag = "##FILTER=<ID=SHORT_ANCHOR,Description=\"Read pairs supporting the SV are limited to a very small portion of the genome.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, short_anchor_flt_tag, &len));

	const char* low_qual_flt_tag = "##FILTER=<ID=LOW_ALT_CONSENSUS_SCORE,Description=\"Assembled alternative consensus aligns better to REF than to ALT.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, low_qual_flt_tag, &len));

	// add INFO tags
	const char* svtype_tag = "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svtype_tag, &len));

	const char* end_tag = "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, end_tag, &len));

	const char* invpos_tag = "##INFO=<ID=INVPOS,Number=2,Type=Integer,Description=\"The start and end positions of the inverted sequence. This is only relevant for inversions.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, invpos_tag, &len));

	const char* svlen_tag = "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svlen_tag, &len));

	const char* svinslen_tag = "##INFO=<ID=SVINSLEN,Number=1,Type=Integer,Description=\"Length of the inserted sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svinslen_tag, &len));

	const char* hpid_tag = "##INFO=<ID=HPID,Number=1,Type=Integer,Description=\"Identifier of the local haplotype represented by this call.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, hpid_tag, &len));

	const char* junction_remap_ref_range_tag = "##INFO=<ID=JUNCTION_REMAP_REF_RANGE,Number=2,Type=Integer,Description=\"1-based half-open [start,end) coordinates of the reference interval covered when the putative alternate junction allele was remapped to the reference.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, junction_remap_ref_range_tag, &len));

	const char* svinsseq_tag = "##INFO=<ID=SVINSSEQ,Number=1,Type=String,Description=\"Inserted sequence.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, svinsseq_tag, &len));

	const char* left_svinsseq_tag = "##INFO=<ID=LEFT_SVINSSEQ,Number=1,Type=String,Description=\"Known left side of insertion for an insertion of unknown length.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, left_svinsseq_tag, &len));

	const char* right_svinsseq_tag = "##INFO=<ID=RIGHT_SVINSSEQ,Number=1,Type=String,Description=\"Known right side of insertion for an insertion of unknown length.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, right_svinsseq_tag, &len));

	const char* infsvinsseq_tag = "##INFO=<ID=SVINSSEQ_INFERRED,Number=1,Type=String,Description=\"Inferred insertion sequence. When the inserted sequence is too long "
		"to be fully assembled but SurVeyor suspects it to be a transposition, it uses the reference to infer the content of the insertion. Not guaranteed to be accurate. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, infsvinsseq_tag, &len));

	const char* source_tag = "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source of the SV.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, source_tag, &len));

	const char* est_size_tag = "##INFO=<ID=EST_SIZE,Number=1,Type=Integer,Description=\"Estimated size of the imprecise event. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, est_size_tag, &len));

	const char* og_range_tag = "##INFO=<ID=ORIGINAL_RANGE,Number=1,Type=String,Description=\"Unadjusted imprecise range predicted by discordant pairs. \">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, og_range_tag, &len));

	const char* imprecise_tag = "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"The reported boundaries are not precise.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, imprecise_tag, &len));

	const char* ia_tag = "##INFO=<ID=INCOMPLETE_ASSEMBLY,Number=0,Type=Flag,Description=\"The inserted sequence is too long, and it could not be fully assembled using short reads.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ia_tag, &len));

	const char* aux_indels_tag = "##INFO=<ID=AUX_INDELS,Number=1,Type=String,Description=\"Auxiliary indels supporting the variant.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, aux_indels_tag, &len));

	const char* aux_snps_tag = "##INFO=<ID=AUX_SNPS,Number=1,Type=String,Description=\"Auxiliary SNPs supporting the variant.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, aux_snps_tag, &len));

	const char* exp_alt_freq_tag = "##INFO=<ID=EXP_ALT_READS_FREQ,Number=2,Type=Float,Description=\"Expected frequency of reads supporting the ALT bp1 and bp2 allele.\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, exp_alt_freq_tag, &len));

	// add ALT
	const char* del_alt_tag = "##ALT=<ID=DEL,Description=\"Deletion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, del_alt_tag, &len));

	const char* dup_alt_tag = "##ALT=<ID=DUP,Description=\"Tandem Duplication\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, dup_alt_tag, &len));

	const char* ins_alt_tag = "##ALT=<ID=INS,Description=\"Insertion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, ins_alt_tag, &len));

	const char* inv_alt_tag = "##ALT=<ID=INV,Description=\"Invertion\">";
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, inv_alt_tag, &len));

	add_fmt_tags(header);

	std::string cmd_tag = "##SurVeyorCommand=" + command;
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, cmd_tag.c_str(), &len));

	auto now = std::chrono::system_clock::now();
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	std::string version_tag = "##SurVeyorVersion=" + config.version + "; Date=" + std::ctime(&now_time);
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, version_tag.c_str(), &len));

	std::stringstream called_by_ss;
	called_by_ss << "##calledBy=SurVeyor " << config.version << "; ";
	called_by_ss << "seed: " << config.seed << "; ";
	called_by_ss << "min-sv-size: " << config.min_sv_size << "; ";
	called_by_ss << "min-clip-len: " << config.min_clip_len << "; ";
	called_by_ss << "max-seq-error: " << config.max_seq_error << "; ";
	called_by_ss << "max-trans-size: " << config.max_trans_size << "; ";
	called_by_ss << "min-stable-mapq: " << config.min_stable_mapq << "; ";
	called_by_ss << "min-diff-hsr: " << config.min_diff_hsr << "; ";
	called_by_ss << "sampling-regions: " << (config.sampling_regions.empty() ? "no" : config.sampling_regions) << "; ";
	called_by_ss << "per-contig-stats: " << (config.per_contig_stats ? "true" : "false") << "; ";
	std::string called_by = called_by_ss.str();
	bcf_hdr_add_hrec(header, bcf_hdr_parse_line(header, called_by.c_str(), &len));

	if (!sample_name.empty()) {
		bcf_hdr_add_sample(header, sample_name.c_str());
	}
	
	return header;
}

void update_alleles_info(bcf_hdr_t* hdr, bcf1_t* bcf_entry, char* chr_seq, hts_pos_t start, hts_pos_t end, std::string svtype, int svlen, std::string ins_seq) {
	
	std::string ref_allele = std::string(1, chr_seq[start]);
	std::string alt_allele = std::string(1, chr_seq[start]);
	if (svtype == "BND") {
		throw std::runtime_error("BND not supported in update_alleles_info.");
	} else if (svtype == "DEL" && end-start <= 100) {
		ref_allele += std::string(chr_seq + start + 1, end-start);
		alt_allele += ins_seq;
	} else if (svtype == "INS" && end-start <= 100 && !ins_seq.empty()) {
		ref_allele += std::string(chr_seq + start + 1, end-start);
		alt_allele += ins_seq;
	} else {
		alt_allele = "<" + svtype + ">";
	}
	std::string alleles = ref_allele + "," + alt_allele;
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	int int_conv = end+1;
	bcf_update_info_int32(hdr, bcf_entry, "END", &int_conv, 1);

	if (svtype != "BND" && svtype != "DUP") {
		int_conv = svlen;
		bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
		if (!ins_seq.empty()) {
			int_conv = ins_seq.length();
			bcf_update_info_int32(hdr, bcf_entry, "SVINSLEN", &int_conv, 1);
		}
	} else if (svtype == "DUP") {
		bcf_update_info_int32(hdr, bcf_entry, "SVLEN", NULL, 0);
	}
}

void sv2bcf(bcf_hdr_t* hdr, bcf1_t* bcf_entry, sv_t* sv, char* chr_seq, bool for_gt = false) {
	bcf_clear(bcf_entry);
	
	bcf_entry->rid = bcf_hdr_name2id(hdr, sv->chr.c_str());
	bcf_entry->pos = sv->start;
	bcf_update_id(hdr, bcf_entry, sv->id.c_str());

	std::string ref_allele = std::string(1, chr_seq[sv->start]);
	std::string alt_allele = std::string(1, chr_seq[sv->start]);
	if (sv->svtype() == "BND") {
		breakend_t* bnd = (breakend_t*) sv;
		if (bnd->left_facing) {
			alt_allele = std::string("[") + bcf_seqname(hdr, bcf_entry) + ":" + std::to_string(bnd->end+1) + "[" + ref_allele;
		} else {
			alt_allele = std::string("]") + bcf_seqname(hdr, bcf_entry) + ":" + std::to_string(bnd->end+1) + "]" + ref_allele; 
		}
	} else if (sv->svtype() == "DEL" && sv->end-sv->start <= 100) {
		ref_allele += std::string(chr_seq + sv->start + 1, sv->end-sv->start);
		alt_allele += sv->ins_seq;
	} else if (sv->svtype() == "INS" && sv->end-sv->start <= 100 && !sv->incomplete_ins_seq()) {
		ref_allele += std::string(chr_seq + sv->start + 1, sv->end-sv->start);
		alt_allele += sv->ins_seq;
	} else {
		alt_allele = "<" + sv->svtype() + ">";
		if (!sv->ins_seq.empty()) {
			if (sv->incomplete_ins_seq()) {
				std::string left_ins_seq = sv->ins_seq.substr(0, sv->known_seq_prefix_len());
				std::string right_ins_seq = sv->ins_seq.substr(sv->ins_seq.length()-sv->known_seq_suffix_len(), sv->known_seq_suffix_len());
				if (!left_ins_seq.empty()) {
					bcf_update_info_string(hdr, bcf_entry, "LEFT_SVINSSEQ", left_ins_seq.c_str());
				}
				if (!right_ins_seq.empty()) {
					bcf_update_info_string(hdr, bcf_entry, "RIGHT_SVINSSEQ", right_ins_seq.c_str());
				}
			} else {
				bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ", sv->ins_seq.c_str());
			}
		}
	}
	std::string alleles = ref_allele + "," + alt_allele;
	bcf_update_alleles_str(hdr, bcf_entry, alleles.c_str());

	int n_samples = bcf_hdr_nsamples(hdr);

	// add FILTERs
	if (n_samples) {
		std::string ft;
		for (std::string& filter : sv->sample_info.filters) {
			if (!ft.empty()) ft += ";";
			ft += filter;
		}
		if (ft.empty()) {
			ft = "PASS";
		}
		const char* ft_cstr[1];
		ft_cstr[0] = (char*) ft.c_str();
		bcf_update_format_string(hdr, bcf_entry, "FT", ft_cstr, 1);
	}

	// add INFO
	int int_conv = sv->end+1;
	bcf_update_info_int32(hdr, bcf_entry, "END", &int_conv, 1);
	
	bcf_update_info_string(hdr, bcf_entry, "SVTYPE", sv->svtype().c_str());
	if (sv->svtype() == "DUP") {
		bcf_update_info_int32(hdr, bcf_entry, "SVLEN", NULL, 0);
	} else if (!sv->incomplete_ins_seq()) {
		int_conv = sv->svlen();
		bcf_update_info_int32(hdr, bcf_entry, "SVLEN", &int_conv, 1);
		if (!sv->ins_seq.empty()) {
			int_conv = sv->ins_seq.length();
			bcf_update_info_int32(hdr, bcf_entry, "SVINSLEN", &int_conv, 1);
		}
	}
	if (!sv->source.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SOURCE", sv->source.c_str());
	}
	if (sv->hpid == 0) sv->hpid = make_hpid(*sv);
	bcf_update_info_int32(hdr, bcf_entry, "HPID", &sv->hpid, 1);
	if (sv->junction_remap_ref_beg != HTS_POS_MIN && sv->junction_remap_ref_end != HTS_POS_MIN) {
		int junction_remap_ref_range[] = {(int) sv->junction_remap_ref_beg + 1, (int) sv->junction_remap_ref_end + 1};
		bcf_update_info_int32(hdr, bcf_entry, "JUNCTION_REMAP_REF_RANGE", junction_remap_ref_range, 2);
	}
	if (!sv->inferred_ins_seq.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "SVINSSEQ_INFERRED", sv->inferred_ins_seq.c_str());
	}

	bcf_update_info_flag(hdr, bcf_entry, "IMPRECISE", "", sv->imprecise);
	bcf_update_genotypes(hdr, bcf_entry, sv->sample_info.gt.data(), sv->sample_info.gt.size());

	if (sv->incomplete_ins_seq()) {
		bcf_update_info_flag(hdr, bcf_entry, "INCOMPLETE_ASSEMBLY", "", 1);
	}

	// auxiliary SNPs
	std::string aux_snps_str;
	for (const auto& snp : sv->aux_snps) {
		if (!aux_snps_str.empty()) {
			aux_snps_str += ",";
		}
		aux_snps_str += std::to_string(snp.pos+1) + ":" + snp.alt_base;
	}
	if (!aux_snps_str.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "AUX_SNPS", aux_snps_str.c_str());
	}

	// auxiliary indels
	std::string aux_indels_str;
	for (const auto& indel : sv->aux_indels) {
		if (!aux_indels_str.empty()) {
			aux_indels_str += ",";
		}
		aux_indels_str += std::to_string(indel->start+1) + ":" + std::to_string(indel->end+1) + ":" + indel->ins_seq;
	}
	if (!aux_indels_str.empty()) {
		bcf_update_info_string(hdr, bcf_entry, "AUX_INDELS", aux_indels_str.c_str());
	}

	if (sv->svtype() == "DEL") {
		if (!for_gt) {
			deletion_t* del = (deletion_t*) sv;
			if (!del->original_range.empty()) {
				bcf_update_info_string(hdr, bcf_entry, "ORIGINAL_RANGE", del->original_range.c_str());
			}
			if (del->imprecise && del->estimated_size != deletion_t::SIZE_NOT_COMPUTED) {
				bcf_update_info_int32(hdr, bcf_entry, "EST_SIZE", &(del->estimated_size), 1);
			}
		}
	} else if (sv->svtype() == "DUP") {
		duplication_t* dup = static_cast<duplication_t*>(sv);
		bcf_update_info_flag(hdr, bcf_entry, "DUP_CN_UNRESOLVED", "", dup->cn_unresolved);
	} else if (sv->svtype() == "INV") {
		inversion_t* inv = (inversion_t*) sv;
		if (inv->inv_start != inv->start || inv->inv_end != inv->end) {
			int invpos[2] = {(int) inv->inv_start+1, (int) inv->inv_end+1};
			bcf_update_info_int32(hdr, bcf_entry, "INVPOS", invpos, 2);
		}
	}
}

std::string get_sv_type(bcf_hdr_t* hdr, bcf1_t* sv) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, "SVTYPE", &data, &len) < 0) {
		// try to determine svtype from ALT
		bcf_unpack(sv, BCF_UN_STR);
		if (sv->d.allele[1][0] == '<') {
			std::string alt = sv->d.allele[1];
			return alt.substr(1, alt.length()-2);
		} else {
			int ref_len = strlen(sv->d.allele[0]);
			int alt_len = strlen(sv->d.allele[1]);
			if (alt_len < ref_len) return "DEL";
			if (alt_len > ref_len) return "INS";
			return ref_len == 1 ? "SNP" : "MNP";
		}
		std::cerr << "Failed to determine SVTYPE for sv " << std::string(sv->d.id) << std::endl;
		return "";
    }
    std::string svtype = data;
    free(data);
    return svtype;
}

int get_sv_len(bcf_hdr_t* hdr, bcf1_t* sv) {
	int* data = NULL;
	int size = 0;
	bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
	if (size > 0) {
		int len = data[0];
		free(data);
		return len;
	}

	if (sv->d.allele[1][0] != '<') {
		return strlen(sv->d.allele[1])-strlen(sv->d.allele[0]);
	}

	return bcf_int32_missing;
}

int get_sv_end(bcf_hdr_t* hdr, bcf1_t* sv) {
    int* data = NULL;
    int size = 0;
    bcf_get_info_int32(hdr, sv, "END", &data, &size);
    if (size > 0) {
        int end = data[0];
        free(data);
        return end-1; // return 0-based
    }

	if (sv->d.allele[1][0] != '<') {
		return sv->pos + strlen(sv->d.allele[0]) - 1;
	}

	if (get_sv_type(hdr, sv) == "INS") return sv->pos;

	int svlen = get_sv_len(hdr, sv);
	if (svlen != bcf_int32_missing) {
		return sv->pos + abs(svlen);
	}

    throw std::runtime_error("SV " + std::string(sv->d.id) + "has no END or SVLEN annotation.");
}

std::string get_sv_info_str(bcf_hdr_t* hdr, bcf1_t* sv, std::string info) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, info.c_str(), &data, &len) < 0) {
		return "";
    }
    std::string svtype = data;
    free(data);
    return svtype;
}

std::string get_ins_seq(bcf_hdr_t* hdr, bcf1_t* sv) {
	
	bcf_unpack(sv, BCF_UN_INFO);

	// LEFT_SVINSSEQ + "-" + RIGHT_SVINSSEQ, if they exist
	std::string left_ins_seq = get_sv_info_str(hdr, sv, "LEFT_SVINSSEQ");
	std::string right_ins_seq = get_sv_info_str(hdr, sv, "RIGHT_SVINSSEQ");
	if (!left_ins_seq.empty() || !right_ins_seq.empty()) {
		return left_ins_seq + "-" + right_ins_seq;
	}

	// priority to the ALT allele, if it is not symbolic and longer than just the padding base
	char c = toupper(sv->d.allele[1][0]);
	if ((c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') && strlen(sv->d.allele[1]) > 1) {
		// return everything after the padding base
		return sv->d.allele[1] + 1;
	}

	// otherwise, look for SVINSSEQ (compliant with Manta)
	std::string ins_seq = get_sv_info_str(hdr, sv, "SVINSSEQ");
	if (!ins_seq.empty()) {
		return ins_seq;
	}

	return "";
}

void copy_all_info(bcf_hdr_t* hdr, bcf1_t* src, bcf1_t* dest) {
    bcf_unpack(src, BCF_UN_INFO);

	for (int i = 0; i < src->n_info; ++i) {
        bcf_info_t* info = &src->d.info[i];
        
        // Get the string tag (key) from the header using the internal ID
        const char *key = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);
        if (!key) continue;

        switch (info->type) {
            case BCF_HT_INT: {
                int32_t *vals = NULL; 
                int nvals = 0;
                if (bcf_get_info_int32(hdr, src, key, &vals, &nvals) > 0) {
                    bcf_update_info_int32(hdr, dest, key, vals, nvals);
                }
                free(vals);
                break;
            }
            case BCF_HT_REAL: {
                float *vals = NULL; 
                int nvals = 0;
                if (bcf_get_info_float(hdr, src, key, &vals, &nvals) > 0) {
                    bcf_update_info_float(hdr, dest, key, vals, nvals);
                }
                free(vals);
                break;
            }
            case BCF_HT_STR: {
                char *vals = NULL; 
                int nvals = 0;
                if (bcf_get_info_string(hdr, src, key, &vals, &nvals) > 0) {
                    // bcf_update_info_string expects a single null-terminated string
                    bcf_update_info_string(hdr, dest, key, vals);
                }
                free(vals);
                break;
            }
            case BCF_HT_FLAG: {
                bcf_update_info_flag(hdr, dest, key, NULL, 1);
                break;
            }
            default:
                break;
        }
    }
}

void copy_all_fmt(bcf_hdr_t* hdr, bcf1_t* src, bcf1_t* dest) {
    bcf_unpack(src, BCF_UN_FMT);

    int i;
    for (i = 0; i < src->n_fmt; i++) {
        bcf_fmt_t* fmt = &src->d.fmt[i];
        const char* key = bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id);
        int res = 0;

        if (fmt->type == BCF_BT_FLOAT) {
            float *vals = NULL;
            int n_vals = 0;
            // Read values (HTSlib handles memory allocation)
            if (bcf_get_format_float(hdr, src, key, &vals, &n_vals) >= 0) {
                bcf_update_format_float(hdr, dest, key, vals, n_vals);
                free(vals);
            }
        } else if (fmt->type == BCF_BT_CHAR) {
            char** vals = NULL; // Array of strings
            int n_vals = 0;
            // Note: bcf_get_format_string allocates an array of pointers
            if (bcf_get_format_string(hdr, src, key, &vals, &n_vals) >= 0) {
                bcf_update_format_string(hdr, dest, key, (const char**)vals, bcf_hdr_nsamples(hdr));
                free(vals[0]); 
                free(vals);
            }
        } else {
            // Handle all Integer types (INT8, INT16, INT32) via INT32 interface
            int32_t* vals = NULL;
            int n_vals = 0;
            if (bcf_get_format_int32(hdr, src, key, &vals, &n_vals) >= 0) {
                // Special handling for Genotype (GT) fields
                if (strcmp(key, "GT") == 0) {
                     bcf_update_genotypes(hdr, dest, vals, n_vals);
                } else {
                     bcf_update_format_int32(hdr, dest, key, vals, n_vals);
                }
                free(vals);
            }
        }
    }
}

float get_sv_epr(bcf_hdr_t *hdr, bcf1_t *b) {
    int ngt = 0;
    float *epr = NULL;
    if (bcf_get_format_float(hdr, b, "EPR", &epr, &ngt) > 0) {
        float value = epr[0];
        free(epr);
        return value;
    }
    free(epr);
    return -1.0;
}

float get_sv_hopr(bcf_hdr_t *hdr, bcf1_t *b) {
    int ngt = 0;
    float *hopr = NULL;
    if (bcf_get_format_float(hdr, b, "HOPR", &hopr, &ngt) > 0) {
        float value = hopr[0];
        free(hopr);
        return value;
    }
    free(hopr);
    return -1.0;
}

std::vector<int> get_bcf_gt(bcf_hdr_t *hdr, bcf1_t *b) {
	int ngt = 0;
	int* gt_arr = NULL;
	std::vector<int> gt_vector;
	if (bcf_get_genotypes(hdr, b, &gt_arr, &ngt) > 0) {
		for (int i = 0; i < ngt; i++) {
			gt_vector.push_back(gt_arr[i]);
		}
	}
	free(gt_arr);
	return gt_vector;
}

std::shared_ptr<sv_t> bcf_to_sv(bcf_hdr_t* hdr, bcf1_t* b) {

	int* data = NULL;
	int len = 0;

	std::shared_ptr<consensus_t> rc_consensus = nullptr, lc_consensus = nullptr;

	data = NULL;
	len = 0;
	int imprecise = bcf_get_info_flag(hdr, b, "IMPRECISE", &data, &len);
	free(data);

	std::string svtype = get_sv_type(hdr, b);
	hts_pos_t start = b->pos, end = get_sv_end(hdr, b);
	hts_pos_t left_anchor_mapping_start, left_anchor_mapping_end, right_anchor_mapping_start, right_anchor_mapping_end;
	if (svtype == "DUP") {
		left_anchor_mapping_start = start, left_anchor_mapping_end = start + 150;
		right_anchor_mapping_start = end - 150, right_anchor_mapping_end = end;
	} else if (svtype == "INV") {
		left_anchor_mapping_start = start - 150, left_anchor_mapping_end = start;
		right_anchor_mapping_start = start, right_anchor_mapping_end = start + 150;
	} else {
		left_anchor_mapping_start = start - 150, left_anchor_mapping_end = start;
		right_anchor_mapping_start = end, right_anchor_mapping_end = end + 150;
	}

	auto left_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(std::max(hts_pos_t(0), left_anchor_mapping_start), left_anchor_mapping_end, 0, 0);
	auto right_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(std::max(hts_pos_t(0), right_anchor_mapping_start), right_anchor_mapping_end, 0, 0);

	std::shared_ptr<sv_t> sv;
	if (svtype == "DEL") {
		sv = std::make_shared<deletion_t>(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln);
	} else if (svtype == "DUP") {
		auto dup = std::make_shared<duplication_t>(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln);
		data = NULL;
		len = 0;
		dup->cn_unresolved = bcf_get_info_flag(hdr, b, "DUP_CN_UNRESOLVED", &data, &len) > 0;
		free(data);
		sv = dup;
	} else if (svtype == "INS") {
		sv = std::make_shared<insertion_t>(bcf_seqname_safe(hdr, b), b->pos, end, get_ins_seq(hdr, b), rc_consensus, lc_consensus, left_anchor_aln, right_anchor_aln);
	} else if (svtype == "INV") {
		auto lbp_left_anchor_aln = left_anchor_aln, lbp_right_anchor_aln = right_anchor_aln;
		hts_pos_t rbp_left_split_mapping_start, rbp_left_split_mapping_end, rbp_right_split_mapping_start, rbp_right_split_mapping_end;
		rbp_left_split_mapping_start = end - 150, rbp_left_split_mapping_end = end;
		rbp_right_split_mapping_start = end, rbp_right_split_mapping_end = end + 150;
		auto rbp_left_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(std::max(hts_pos_t(0), rbp_left_split_mapping_start), rbp_left_split_mapping_end, 0, 0);
		auto rbp_right_anchor_aln = std::make_shared<sv_t::anchor_aln_t>(std::max(hts_pos_t(0), rbp_right_split_mapping_start), rbp_right_split_mapping_end, 0, 0);
		sv = std::make_shared<inversion_t>(bcf_seqname_safe(hdr, b), b->pos, get_sv_end(hdr, b), get_ins_seq(hdr, b), rc_consensus, lc_consensus, lbp_left_anchor_aln, lbp_right_anchor_aln, rbp_left_anchor_aln, rbp_right_anchor_aln);

		data = NULL;
		bcf_get_info_int32(hdr, b, "INVPOS", &data, &len);
		inversion_t* inv = (inversion_t*) sv.get();
		if (len > 1) {
			inv->inv_start = data[0]-1;
			inv->inv_end = data[1]-1;
		}
		free(data);
	} else if (svtype == "BND") {
		std::string alt = b->d.allele[1];
		size_t colon_pos = alt.find(':');
		if (colon_pos == std::string::npos) {
			std::cerr << "Invalid BND ALT format for sv " << b->d.id << ": " << alt << std::endl;
			return NULL;
		}
		size_t pos_start = colon_pos + 1, pos_end = alt.find_last_of("[]");
		if (pos_end == std::string::npos || pos_end <= pos_start) {
			std::cerr << "Invalid BND ALT format for sv " << b->d.id << ": " << alt << std::endl;
			return NULL;
		}
    	hts_pos_t pos2 = std::stoi(alt.substr(pos_start, pos_end-pos_start));
		bool left_facing = false;
		if (alt[pos_end] == '[') {
			left_facing = true;
		}
			sv = std::make_shared<breakend_t>(bcf_seqname_safe(hdr, b), b->pos, pos2-1, get_ins_seq(hdr, b), left_anchor_aln, right_anchor_aln, left_facing);
	} else {
		return NULL;
	}
	sv->imprecise = (imprecise == 1);

	char* s_data = NULL;
	len = 0;
	int res = bcf_get_info_string(hdr, b, "SVINSSEQ_INFERRED", (void**) &s_data, &len);
	if (len > 0) {
		sv->inferred_ins_seq = s_data;
	}
	free(s_data);

	s_data = NULL;
	len = 0;
	bcf_get_info_string(hdr, b, "AUX_SNPS", (void**) &s_data, &len);
	if (len > 0) {
		std::istringstream ss(s_data);
		std::string snp_str;
		while (std::getline(ss, snp_str, ',')) {
			sv->aux_snps.push_back(snp_t(snp_str));
		}
	}
	free(s_data);

	s_data = NULL;
	len = 0;
	bcf_get_info_string(hdr, b, "AUX_INDELS", (void**) &s_data, &len);
	if (len > 0) {
		std::istringstream ss(s_data);
		std::string indel_str;
		while (std::getline(ss, indel_str, ',')) {
			// parse start:end:ins_seq
			std::istringstream indel_ss(indel_str);
			std::string start_str, end_str, ins_seq;
			std::getline(indel_ss, start_str, ':');
			std::getline(indel_ss, end_str, ':');
			std::getline(indel_ss, ins_seq, ':');
			hts_pos_t start = std::stoll(start_str);
			hts_pos_t end = std::stoll(end_str);
			if (start == end) {
				sv->aux_indels.push_back(std::make_shared<insertion_t>(bcf_seqname_safe(hdr, b), start-1, end-1, ins_seq, nullptr, nullptr, nullptr, nullptr));
			} else {
				sv->aux_indels.push_back(std::make_shared<deletion_t>(bcf_seqname_safe(hdr, b), start-1, end-1, ins_seq, nullptr, nullptr, nullptr, nullptr));
			}
		}
	}
	free(s_data);

	sv->id = b->d.id;
	sv->source = get_sv_info_str(hdr, b, "SOURCE");
	data = NULL;
	len = 0;
	if (bcf_get_info_int32(hdr, b, "HPID", &data, &len) > 0 && len > 0) {
		sv->hpid = data[0];
	}
	free(data);

	data = NULL;
	len = 0;
	if (bcf_get_info_int32(hdr, b, "JUNCTION_REMAP_REF_RANGE", &data, &len) > 0 && len > 1) {
		sv->junction_remap_ref_beg = data[0] - 1;
		sv->junction_remap_ref_end = data[1] - 1;
	}
	free(data);
	if (sv->hpid == 0) sv->hpid = make_hpid(*sv);
	for (auto& aux_indel : sv->aux_indels) {
		aux_indel->hpid = sv->hpid;
	}

	char** ss_data = NULL;
	len = 0;
	bcf_get_format_string(hdr, b, "FT", &ss_data, &len);
	std::string ft = "PASS";
	if (len > 0) {
		ft = ss_data[0];
		free(ss_data[0]);
	}
	free(ss_data);

	std::istringstream tokenStream(ft);
	std::string token;
	while (std::getline(tokenStream, token, ';')) {
		std::string filter = token;
		if (filter != "PASS") sv->sample_info.filters.push_back(token);
	}

	sv->sample_info.gt = get_bcf_gt(hdr, b);
	std::sort(sv->sample_info.gt.begin(), sv->sample_info.gt.end());

	float* f_data = NULL;
	len = 0;
	if (bcf_get_format_float(hdr, b, "EPR", &f_data, &len) > 0 && len > 0) {
		sv->sample_info.epr = f_data[0];
		free(f_data);
	}

	return sv;
}

int find_sample_index(bcf_hdr_t* hdr, std::string& sample_name) {
    for (int i = 0; i < bcf_hdr_nsamples(hdr); i++) {
        if (sample_name == hdr->samples[i]) {
            return i;
        }
    }
    return -1;  // Sample not found
}

// returns (via parameter) the imap to be used with bcf_subset
bcf_hdr_t* bcf_subset_header(bcf_hdr_t* in_hdr, std::string sample_name, int*& imap) {
	imap = new int[1];
    bcf_hdr_t* out_hdr;
    int sample_idx = find_sample_index(in_hdr, sample_name);
    if (sample_idx >= 0) {
        char* samples[1] = { strdup(sample_name.c_str()) };
        out_hdr = bcf_hdr_subset(in_hdr, 1, samples, imap);
    } else {
        out_hdr = bcf_hdr_subset(in_hdr, 0, NULL, NULL);
        bcf_hdr_add_sample(out_hdr, sample_name.c_str());
        imap[0] = -1;
    }
	if (bcf_hdr_sync(out_hdr) < 0) {
		throw std::runtime_error("Failed to sync header.");
	}
	return out_hdr;
}

int count_alt_alleles(bcf_hdr_t* hdr, bcf1_t* sv) {
    int* gt = nullptr;
    int ngt = 0;
    if (bcf_get_genotypes(hdr, sv, &gt, &ngt) < 0 || ngt < 2) {
        free(gt);
        return 0;
    }
    int count = 0;
    for (int i = 0; i < ngt; i++) {
        if (bcf_gt_allele(gt[i]) > 0) count++;
    }
    free(gt);
    return count;
}

bcf1_t* generate_snp(bcf_hdr_t* hdr, std::string chr, hts_pos_t pos, char ref_base, char alt_base, 
	std::string sv_id, std::vector<int>& gt) {

	bcf1_t* vcf_record = bcf_init();
	vcf_record->rid = bcf_hdr_name2id(hdr, chr.c_str());
	vcf_record->pos = pos; // 0-based position
	vcf_record->d.id = strdup(sv_id.c_str());

	std::string ref_allele(1, ref_base);
	std::string alt_allele(1, alt_base);
	const char* alleles[] = {ref_allele.c_str(), alt_allele.c_str()};
	bcf_update_alleles(hdr, vcf_record, alleles, 2);

	bcf_update_genotypes(hdr, vcf_record, gt.data(), gt.size());

	return vcf_record;

}

#endif /* VCF_UTILS_H */
