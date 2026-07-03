import pysam, argparse
import features
import xgboost as xgb
import os
import numpy as np
import timeit

class Classifier:
    def load_features_files(model_stage_dir):
        feature_names_by_model = dict()
        if not os.path.isdir(model_stage_dir):
            return feature_names_by_model

        for fname in os.listdir(model_stage_dir):
            if not fname.endswith(".features"):
                continue
            model_name = fname[:-len(".features")]
            features_fname = os.path.join(model_stage_dir, fname)
            with open(features_fname) as f:
                feature_names_by_model[model_name] = [line.strip() for line in f if line.strip()]
        return feature_names_by_model

    def write_vcf(vcf_reader, vcf_header, svid_to_gt, svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa, out_vcf_fname, stats_fname):
        stats = features.load_stats(stats_fname)

        vcf_writer = pysam.VariantFile(out_vcf_fname, 'wz', header=vcf_header)
        for record in vcf_reader.fetch():
            if not features.Features.skips_ml_genotyping(record):
                record.filter.clear()
                record.samples[0]['EPR'] = None
                record.samples[0]['HOPR'] = None
                record.samples[0]['EXPR'] = None
                record.samples[0]['EREFA'] = None
                record_id = features.Features.generate_id(record)
                if record_id in svid_to_gt:
                    gt = (svid_to_gt[record_id]//2, 1 if svid_to_gt[record_id] >= 1 else 0)
                    max_is, read_len = features.get_stat(stats, 'max_is', record.chrom), features.get_stat(stats, 'read_len', record.chrom)
                    if features.Features.get_model_name(record, max_is, read_len) in ("DUP_LARGE", "DUP_LARGE_NOEXL", "INS_TO_DUP_LARGE"):# or \
                        # for large duplications, and for likely multi-allelic duplications, only output ./1 genotypes
                        if gt[1] == 1:
                            gt = (None, 1)
                    record.samples[0]['GT'] = gt
                    record.samples[0]['EPR'] = float(svid_to_epr[record_id])
                    if svid_to_gt[record_id] >= 1:
                        record.samples[0]['HOPR'] = float(svid_to_hopr[record_id])
                    if record_id in svid_to_expr:
                        record.samples[0]['EXPR'] = float(svid_to_expr[record_id])
                    if record_id in svid_to_erefa:
                        record.samples[0]['EREFA'] = int(svid_to_erefa[record_id])
                else:
                    record.samples[0]['GT'] = (None, None)
                vcf_writer.write(record)
            else:
                vcf_writer.write(record)
        vcf_writer.close()

    def run_classifier(in_vcf, out_vcf, stats_fname, model_dir, threads=1):
        cmd = "Classifier.run_classifier %s %s %s %s --threads %d" % (in_vcf, out_vcf, stats_fname, model_dir, threads)
        start_time = timeit.default_timer()
        print("Executing:", cmd)

        parse_start_time = timeit.default_timer()
        feature_names_by_model = Classifier.load_features_files(os.path.join(model_dir, "yes_or_no"))
        test_data, _, test_variant_ids, _ = \
            features.parse_vcf(in_vcf, stats_fname, "XXX", ignore_gts = True, feature_names_by_model = feature_names_by_model)
        parse_elapsed = timeit.default_timer() - parse_start_time
        print("Feature parsing was run in %.2f seconds" % parse_elapsed)

        svid_to_gt = dict()
        svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa = dict(), dict(), dict(), dict()
        for model_name in test_data:
            model_start_time = timeit.default_timer()
            model_file = os.path.join(model_dir, "yes_or_no", model_name + '.ubj')

            if model_name.startswith("INV"):
                continue

            classifier = xgb.XGBClassifier(n_jobs=threads)
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
                svid_to_gt[positive_variant_ids[i]] = 2 if predictions[i] == 1 else 1
                svid_to_hopr[positive_variant_ids[i]] = hoprs[i][1]

            het_mask = np.array([svid_to_gt[variant_id] == 1 for variant_id in positive_variant_ids])
            het_data = positive_data[het_mask]
            het_variant_ids = positive_variant_ids[het_mask]

            if len(het_data) > 0:
                for variant_id in het_variant_ids:
                    svid_to_erefa[variant_id] = 0

                model_file = os.path.join(model_dir, "eref", model_name + ".ubj")
                if os.path.exists(model_file):
                    classifier.load_model(model_file)
                    erefa_predictions = classifier.predict(het_data)
                    for i in range(len(erefa_predictions)):
                        svid_to_erefa[het_variant_ids[i]] = int(erefa_predictions[i])

            model_file = os.path.join(model_dir, "exact", model_name + ".ubj")
            if os.path.exists(model_file):
                classifier.load_model(model_file)
                exprs = classifier.predict_proba(positive_data)
                for i in range(len(positive_variant_ids)):
                    svid_to_expr[positive_variant_ids[i]] = exprs[i][1]

        # write the predictions to a VCF file
        write_start_time = timeit.default_timer()
        vcf_reader = pysam.VariantFile(in_vcf)
        header = vcf_reader.header
        if 'EPR' not in header.formats:
            header.add_line('##FORMAT=<ID=EPR,Number=1,Type=Float,Description="Probability of the SV existing in the sample, according to the ML model.">')
        if 'HOPR' not in header.formats:
            header.add_line('##FORMAT=<ID=HOPR,Number=1,Type=Float,Description="Probability of an existing SV to be homozygous, according to the ML model.">')
        if 'EXPR' not in header.formats:
            header.add_line('##FORMAT=<ID=EXPR,Number=1,Type=Float,Description="Probability of the SV to be represented exactly, according to the ML model.">')
        if 'EREFA' not in header.formats:
            header.add_line('##FORMAT=<ID=EREFA,Number=1,Type=Integer,Description="Whether the EREFA-stage classifier requires the other allele to be reference.">')
        Classifier.write_vcf(vcf_reader, header, svid_to_gt, svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa, out_vcf, stats_fname)
        write_elapsed = timeit.default_timer() - write_start_time
        print("VCF writing was run in %.2f seconds" % write_elapsed)

if __name__ == "__main__":
    cmd_parser = argparse.ArgumentParser(description='Classify SVs using a built ML model.')
    cmd_parser.add_argument('in_vcf', help='Input VCF file.')
    cmd_parser.add_argument('out_vcf', help='Output VCF file.')
    cmd_parser.add_argument('stats', help='Stats of the test VCF file.')
    cmd_parser.add_argument('model_dir', help='Directory containing the trained model.')
    cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use during prediction.')
    cmd_args = cmd_parser.parse_args()
    Classifier.run_classifier(cmd_args.in_vcf, cmd_args.out_vcf, cmd_args.stats, cmd_args.model_dir, threads=cmd_args.threads)
