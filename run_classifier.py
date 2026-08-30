import argparse
import os
import struct
import subprocess
import tempfile
import timeit

import numpy as np
import xgboost as xgb

from training_manifest import load_training_set_sha256


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

    def get_writer_path():
        surveyor_path = os.path.dirname(os.path.realpath(__file__))
        writer_path = os.path.join(surveyor_path, "bin", "write_classifier_vcf")
        if os.path.isfile(writer_path):
            return writer_path
        build_writer_path = os.path.join(surveyor_path, "build", "bin", "write_classifier_vcf")
        if os.path.isfile(build_writer_path):
            return build_writer_path
        raise RuntimeError("write_classifier_vcf was not found. Build SurVeyor before running run_classifier_c.py.")

    def load_feature_bundle(feature_bundle_fname):
        test_data, test_variant_ids, test_record_keys = dict(), dict(), dict()
        with open(feature_bundle_fname, "rb") as feature_bundle:
            if feature_bundle.read(8) != b"SVFEAT2\0":
                raise RuntimeError("Invalid C++ feature bundle.")
            n_models = struct.unpack("<I", feature_bundle.read(4))[0]
            for _ in range(n_models):
                model_name_len = struct.unpack("<I", feature_bundle.read(4))[0]
                model_name = feature_bundle.read(model_name_len).decode("utf-8")
                n_variants = struct.unpack("<Q", feature_bundle.read(8))[0]
                n_features = struct.unpack("<I", feature_bundle.read(4))[0]
                variant_ids = np.fromfile(feature_bundle, dtype="<u8", count=n_variants)
                record_keys = np.fromfile(feature_bundle, dtype="<u8", count=n_variants)
                feature_values = np.fromfile(feature_bundle, dtype="<f8", count=n_variants*n_features).reshape((n_variants, n_features))
                if n_variants > 0:
                    test_data[model_name] = feature_values
                    test_variant_ids[model_name] = variant_ids
                    test_record_keys[model_name] = record_keys
        return test_data, test_variant_ids, test_record_keys

    def write_prediction_bundle(predictions_fname, svid_to_gt, svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa, svid_to_record_key, svid_to_model_name):
        with open(predictions_fname, "wb") as predictions_file:
            predictions_file.write(b"SVPRED1\0")
            predictions_file.write(struct.pack("<Q", len(svid_to_gt)))
            for variant_id, gt in svid_to_gt.items():
                force_missing_ref = int(gt >= 1 and svid_to_model_name[variant_id] in ("DUP_LARGE", "DUP_LARGE_NOEXL", "INS_TO_DUP_LARGE"))
                predictions_file.write(struct.pack("<Qbbfffb", svid_to_record_key[variant_id], int(gt), force_missing_ref, float(svid_to_epr[variant_id]), float(svid_to_hopr.get(variant_id, np.nan)), float(svid_to_expr.get(variant_id, np.nan)), int(svid_to_erefa.get(variant_id, -1))))

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
            test_data, test_variant_ids, test_record_keys = Classifier.load_feature_bundle(feature_bundle_fname)
        print("Feature parsing was run in %.2f seconds" % (timeit.default_timer()-parse_start_time))

        classification_start_time = timeit.default_timer()
        svid_to_gt = dict()
        svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa = dict(), dict(), dict(), dict()
        svid_to_record_key, svid_to_model_name = dict(), dict()
        for model_name in test_data:
            model_file = os.path.join(model_dir, "yes_or_no", model_name+'.ubj')
            classifier = xgb.XGBClassifier(n_jobs=threads)
            classifier.load_model(model_file)
            predictions = classifier.predict(test_data[model_name])
            eprs = classifier.predict_proba(test_data[model_name])
            for i in range(len(predictions)):
                svid_to_gt[test_variant_ids[model_name][i]] = predictions[i]
                svid_to_epr[test_variant_ids[model_name][i]] = eprs[i][1]
                svid_to_record_key[test_variant_ids[model_name][i]] = test_record_keys[model_name][i]
                svid_to_model_name[test_variant_ids[model_name][i]] = model_name

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
        output_dir = os.path.dirname(os.path.abspath(out_vcf))
        with tempfile.TemporaryDirectory(prefix=".surveyor-predictions-", dir=output_dir) as temp_dir:
            predictions_fname = os.path.join(temp_dir, "predictions.bin")
            Classifier.write_prediction_bundle(predictions_fname, svid_to_gt, svid_to_epr, svid_to_hopr, svid_to_expr, svid_to_erefa, svid_to_record_key, svid_to_model_name)
            subprocess.run([Classifier.get_writer_path(), in_vcf, predictions_fname, out_vcf, training_set_sha256 if training_set_sha256 is not None else ".", str(threads)], check=True)
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
