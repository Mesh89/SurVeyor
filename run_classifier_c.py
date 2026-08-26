import argparse
import os
import struct
import subprocess
import tempfile
import timeit

import numpy as np
import pysam
import xgboost as xgb

import features
from training_manifest import load_training_set_sha256, set_training_set_vcf_header


class Classifier:
    def get_extractor_path():
        surveyor_path = os.path.dirname(os.path.realpath(__file__))
        extractor_path = os.path.join(surveyor_path, "bin", "extract_features")
        if os.path.isfile(extractor_path):
            return extractor_path
        build_extractor_path = os.path.join(surveyor_path, "build", "bin", "extract_features")
        if os.path.isfile(build_extractor_path):
            return build_extractor_path
        raise RuntimeError("extract_features was not found. Build SurVeyor before running run_classifier_c.py.")

    def load_feature_bundle(feature_bundle_fname):
        test_data, test_variant_ids = dict(), dict()
        with open(feature_bundle_fname, "rb") as feature_bundle:
            if feature_bundle.read(8) != b"SVFEAT1\0":
                raise RuntimeError("Invalid C++ feature bundle.")
            n_models = struct.unpack("<I", feature_bundle.read(4))[0]
            for _ in range(n_models):
                model_name_len = struct.unpack("<I", feature_bundle.read(4))[0]
                model_name = feature_bundle.read(model_name_len).decode("utf-8")
                n_variants = struct.unpack("<Q", feature_bundle.read(8))[0]
                n_features = struct.unpack("<I", feature_bundle.read(4))[0]
                variant_ids = np.fromfile(feature_bundle, dtype="<u8", count=n_variants)
                feature_values = np.fromfile(feature_bundle, dtype="<f8", count=n_variants*n_features).reshape((n_variants, n_features))
                if n_variants > 0:
                    test_data[model_name] = feature_values
                    test_variant_ids[model_name] = variant_ids
        return test_data, test_variant_ids

    def get_info_string(record, key):
        if key not in record.info:
            return ""
        value = record.info[key]
        if isinstance(value, (list, tuple)):
            return ",".join(str(x) for x in value)
        return str(value)

    def generate_id(record, model_name):
        canonical_id = "%s:%d-%d:%s:%d:%s:%s:%s:%s:%s" % (
            record.chrom, record.pos, record.stop, features.Features.get_svtype(record), features.Features.get_svlen(record),
            features.Features.get_svinsseq(record), Classifier.get_info_string(record, "AUX_SNPS"), Classifier.get_info_string(record, "AUX_INDELS"),
            record.id if record.id is not None else ".", model_name)
        variant_id = 14695981039346656037
        for byte in canonical_id.encode("utf-8"):
            variant_id = ((variant_id ^ byte)*1099511628211) & 0xffffffffffffffff
        return variant_id

    def write_vcf(vcf_reader, vcf_header, svid_to_gt, svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa, out_vcf_fname, stats_fname, threads):
        stats = features.load_stats(stats_fname)
        vcf_writer = pysam.VariantFile(out_vcf_fname, 'wz', header=vcf_header, threads=threads)
        for record in vcf_reader.fetch():
            if not features.Features.skips_ml_genotyping(record):
                record.filter.clear()
                record.samples[0]['EPR'] = None
                record.samples[0]['HOPR'] = None
                record.samples[0]['EXPR'] = None
                record.samples[0]['EREFA'] = None
                max_is, read_len = features.get_stat(stats, 'max_is', record.chrom), features.get_stat(stats, 'read_len', record.chrom)
                model_name = features.Features.get_model_name(record, max_is, read_len)
                record_id = Classifier.generate_id(record, model_name)
                if record_id in svid_to_gt:
                    gt = (svid_to_gt[record_id]//2, 1 if svid_to_gt[record_id] >= 1 else 0)
                    if model_name in ("DUP_LARGE", "DUP_LARGE_NOEXL", "INS_TO_DUP_LARGE"):
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
        print("Executing:", cmd)
        training_set_sha256 = load_training_set_sha256(model_dir)
        if training_set_sha256 is None:
            print("Warning: model directory has no training_manifest.json; the output VCF will not identify its training set.")

        parse_start_time = timeit.default_timer()
        output_dir = os.path.dirname(os.path.abspath(out_vcf))
        with tempfile.TemporaryDirectory(prefix=".surveyor-features-", dir=output_dir) as temp_dir:
            feature_bundle_fname = os.path.join(temp_dir, "features.bin")
            subprocess.run([Classifier.get_extractor_path(), in_vcf, stats_fname, os.path.join(model_dir, "yes_or_no"), feature_bundle_fname, str(threads)], check=True)
            test_data, test_variant_ids = Classifier.load_feature_bundle(feature_bundle_fname)
        print("Feature parsing was run in %.2f seconds" % (timeit.default_timer()-parse_start_time))

        classification_start_time = timeit.default_timer()
        svid_to_gt = dict()
        svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa = dict(), dict(), dict(), dict()
        for model_name in test_data:
            model_file = os.path.join(model_dir, "yes_or_no", model_name+'.ubj')
            classifier = xgb.XGBClassifier(n_jobs=threads)
            classifier.load_model(model_file)
            predictions = classifier.predict(test_data[model_name])
            eprs = classifier.predict_proba(test_data[model_name])
            for i in range(len(predictions)):
                svid_to_gt[test_variant_ids[model_name][i]] = predictions[i]
                svid_to_epr[test_variant_ids[model_name][i]] = eprs[i][1]

            positive_mask = predictions == 1
            positive_data = test_data[model_name][positive_mask]
            positive_variant_ids = test_variant_ids[model_name][positive_mask]
            if len(positive_data) == 0:
                continue
            classifier.load_model(os.path.join(model_dir, "gts", model_name+'.ubj'))
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
                model_file = os.path.join(model_dir, "eref", model_name+".ubj")
                if os.path.exists(model_file):
                    classifier.load_model(model_file)
                    erefa_predictions = classifier.predict(het_data)
                    for i in range(len(erefa_predictions)):
                        svid_to_erefa[het_variant_ids[i]] = int(erefa_predictions[i])

            model_file = os.path.join(model_dir, "exact", model_name+".ubj")
            if os.path.exists(model_file):
                classifier.load_model(model_file)
                exprs = classifier.predict_proba(positive_data)
                for i in range(len(positive_variant_ids)):
                    svid_to_expr[positive_variant_ids[i]] = exprs[i][1]
        print("Classification was run in %.2f seconds" % (timeit.default_timer()-classification_start_time))

        write_start_time = timeit.default_timer()
        vcf_reader = pysam.VariantFile(in_vcf, threads=threads)
        header = vcf_reader.header
        if 'EPR' not in header.formats:
            header.add_line('##FORMAT=<ID=EPR,Number=1,Type=Float,Description="Probability of the SV existing in the sample, according to the ML model.">')
        if 'HOPR' not in header.formats:
            header.add_line('##FORMAT=<ID=HOPR,Number=1,Type=Float,Description="Probability of an existing SV to be homozygous, according to the ML model.">')
        if 'EXPR' not in header.formats:
            header.add_line('##FORMAT=<ID=EXPR,Number=1,Type=Float,Description="Probability of the SV to be represented exactly, according to the ML model.">')
        if 'EREFA' not in header.formats:
            header.add_line('##FORMAT=<ID=EREFA,Number=1,Type=Integer,Description="Whether the EREFA-stage classifier requires the other allele to be reference.">')
        if training_set_sha256 is not None:
            set_training_set_vcf_header(header, training_set_sha256)
        Classifier.write_vcf(vcf_reader, header, svid_to_gt, svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa, out_vcf, stats_fname, threads)
        print("VCF writing was run in %.2f seconds" % (timeit.default_timer()-write_start_time))


if __name__ == "__main__":
    cmd_parser = argparse.ArgumentParser(description='Classify SVs using C++ feature extraction and the existing XGBoost models.')
    cmd_parser.add_argument('in_vcf', help='Input VCF file.')
    cmd_parser.add_argument('out_vcf', help='Output VCF file.')
    cmd_parser.add_argument('stats', help='Stats of the test VCF file.')
    cmd_parser.add_argument('model_dir', help='Directory containing the trained model.')
    cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use during prediction.')
    cmd_args = cmd_parser.parse_args()
    Classifier.run_classifier(cmd_args.in_vcf, cmd_args.out_vcf, cmd_args.stats, cmd_args.model_dir, threads=cmd_args.threads)
