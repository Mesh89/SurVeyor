import pysam, argparse
import features
import xgboost as xgb
import os
import numpy as np

# DUPs with a copy number model probability at least this high are reported as having an unresolved copy number
CNU_THRESHOLD = 0.5

class Classifier:
    def write_vcf(vcf_reader, vcf_header, svid_to_gt, svid_to_epr, svid_to_hopr, svid_to_cnupr, out_vcf_fname, stats_fname):
        stats = features.load_stats(stats_fname)

        vcf_writer = pysam.VariantFile(out_vcf_fname, 'wz', header=vcf_header)
        for record in vcf_reader.fetch():
            if features.Features.get_svtype(record) != "INV":
                record.filter.clear()
                record_id = features.Features.generate_id(record)
                max_is, read_len = features.get_stat(stats, 'max_is', record.chrom), features.get_stat(stats, 'read_len', record.chrom)
                model_name = features.Features.get_model_name(record, max_is, read_len)
                if record_id in svid_to_gt:
                    gt = (svid_to_gt[record_id]//2, 1 if svid_to_gt[record_id] >= 1 else 0)
                    if model_name in ("DUP_LARGE", "DUP_LARGE_IMPRECISE"):# or \
                        # for large duplications, and for likely multi-allelic duplications, only output ./1 genotypes
                        if gt[1] == 1:
                            gt = (None, 1)
                    record.samples[0]['GT'] = gt
                    record.samples[0]['EPR'] = float(svid_to_epr[record_id])
                    if svid_to_gt[record_id] >= 1:
                        record.samples[0]['HOPR'] = float(svid_to_hopr[record_id])
                else:
                    record.samples[0]['GT'] = (None, None)
                # generate_id is not unique (an INS_TO_DUP record can be identical to a DUP record), so check the model too
                if record_id in svid_to_cnupr and model_name in features.CNU_MODELS:
                    # rounded, so that CNU agrees with the CNUPR printed in the VCF (floats are printed with 6 significant digits)
                    cnupr = round(float(svid_to_cnupr[record_id]), 4)
                    record.samples[0]['CNUPR'] = cnupr
                    record.samples[0]['CNU'] = 1 if cnupr >= CNU_THRESHOLD else 0
                vcf_writer.write(record)
            else:
                vcf_writer.write(record)
        vcf_writer.close()

    def run_classifier(in_vcf, out_vcf, stats_fname, model_dir):
        test_data, _, _, test_variant_ids = \
            features.parse_vcf(in_vcf, stats_fname, "XXX", tolerate_no_gts = True)

        svid_to_gt = dict()
        svid_to_epr, svid_to_hopr = dict(), dict()
        for model_name in test_data:
            model_file = os.path.join(model_dir, "yes_or_no", model_name + '.ubj')

            if model_name.startswith("INV"):
                continue

            classifier = xgb.XGBClassifier()
            classifier.load_model(model_file)
            predictions = classifier.predict(test_data[model_name])
            eprs = classifier.predict_proba(test_data[model_name])
            for i in range(len(predictions)):
                svid_to_gt[test_variant_ids[model_name][i]] = predictions[i]
                svid_to_epr[test_variant_ids[model_name][i]] = eprs[i][1]

            if len(test_data[model_name]) == 0:
                continue

            positive_mask = (predictions == 1)
            positive_data = test_data[model_name][positive_mask]
            positive_variant_ids = test_variant_ids[model_name][positive_mask]

            if len(positive_data) == 0:
                continue

            model_file = os.path.join(model_dir, "gts", model_name + '.ubj')
            classifier.load_model(model_file)
            predictions = classifier.predict(positive_data)
            hoprs = classifier.predict_proba(positive_data)
            for i in range(len(predictions)):
                svid_to_gt[positive_variant_ids[i]] = predictions[i] + 1
                svid_to_hopr[positive_variant_ids[i]] = hoprs[i][1]

        # copy number model, for DUPs
        svid_to_cnupr = dict()
        classifier = xgb.XGBClassifier()
        classifier.load_model(os.path.join(model_dir, "cnu", "DUP.ubj"))
        for model_name in features.CNU_MODELS:
            if model_name not in test_data or len(test_data[model_name]) == 0:
                continue
            cnuprs = classifier.predict_proba(test_data[model_name])[:, 1]
            for i in range(len(cnuprs)):
                svid_to_cnupr[test_variant_ids[model_name][i]] = cnuprs[i]

        # write the predictions to a VCF file
        vcf_reader = pysam.VariantFile(in_vcf)
        header = vcf_reader.header
        header.add_line('##FORMAT=<ID=EPR,Number=1,Type=Float,Description="Probability of the SV existing in the sample, according to the ML model.">')
        header.add_line('##FORMAT=<ID=HOPR,Number=1,Type=Float,Description="Probability of an existing SV to be homozygous, according to the ML model.">')
        header.add_line('##FORMAT=<ID=CNUPR,Number=1,Type=Float,Description="Probability that the number of additional tandem copies of this DUP is unresolved, i.e. that the inserted allele is longer than one copy of the duplicated reference interval, according to the ML model.">')
        header.add_line('##FORMAT=<ID=CNU,Number=1,Type=Integer,Description="1 if the number of additional tandem copies of this DUP is unresolved (CNUPR >= 0.5): END and SVLEN describe the duplicated reference interval, but the inserted allele is likely longer. 0 otherwise.">')
        Classifier.write_vcf(vcf_reader, header, svid_to_gt, svid_to_epr, svid_to_hopr, svid_to_cnupr, out_vcf, stats_fname)

if __name__ == "__main__":
    cmd_parser = argparse.ArgumentParser(description='Classify SVs using a built ML model.')
    cmd_parser.add_argument('in_vcf', help='Input VCF file.')
    cmd_parser.add_argument('out_vcf', help='Output VCF file.')
    cmd_parser.add_argument('stats', help='Stats of the test VCF file.')
    cmd_parser.add_argument('model_dir', help='Directory containing the trained model.')
    cmd_args = cmd_parser.parse_args()
    Classifier.run_classifier(cmd_args.in_vcf, cmd_args.out_vcf, cmd_args.stats, cmd_args.model_dir)
