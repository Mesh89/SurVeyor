import argparse, time
import numpy as np
import features
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import xgboost as xgb
from training_manifest import build_training_manifest, write_training_manifest

cmd_parser = argparse.ArgumentParser(description='Train ML model.')
cmd_parser.add_argument('training_prefixes', help='Prefix of the training VCF and FP files.')
cmd_parser.add_argument('outdir')
cmd_parser.add_argument('--model_name', default='ALL', help='Restrict to this model.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use for training.')
cmd_parser.add_argument('--cross-species', action='store_true', help='Use cross-species model.')
cmd_args = cmd_parser.parse_args()

training_data, training_gts, training_exact = defaultdict(list), defaultdict(list), defaultdict(list)

def write_features_file(model_fname, features_names):
    features_fname = os.path.splitext(model_fname)[0] + ".features"
    with open(features_fname, "w") as f:
        for feature_name in features_names:
            f.write(feature_name + "\n")

def remove_model_artifacts(model_fname, importance_fname = None):
    artifact_fnames = [
        model_fname,
        os.path.splitext(model_fname)[0] + ".features",
    ]
    if importance_fname is not None:
        artifact_fnames.append(importance_fname)

    for artifact_fname in artifact_fnames:
        if os.path.exists(artifact_fname):
            os.remove(artifact_fname)

def make_classifier():
    if cmd_args.cross_species:
        return xgb.XGBClassifier(n_estimators=50, max_depth=7, min_child_weight=42, learning_rate=0.1, n_jobs=cmd_args.threads, random_state=42, tree_method='hist')
    return xgb.XGBClassifier(n_estimators=1000, max_depth=8, min_child_weight=16, learning_rate=0.1, n_jobs=cmd_args.threads, random_state=42, tree_method='hist')

def compute_keep_indices(X):
    keep_indices = []
    for i in range(X.shape[1]):
        col = X[:, i]
        finite = np.isfinite(col)
        if finite.sum() == 0:
            continue

        finite_vals = col[finite]
        values_constant = np.all(finite_vals == finite_vals[0])
        missingness_constant = np.all(finite) or np.all(~finite)
        if not (values_constant and missingness_constant):
            keep_indices.append(i)

    return np.array(keep_indices, dtype=np.int32)

def process_vcf(training_prefix, restrict_to_model_name = None):
    gt_labels = features.read_gt_labels(training_prefix + ".gts")
    vcf_training_data, vcf_training_gts, _, vcf_training_exact = \
        features.parse_vcf(training_prefix + ".vcf.gz", training_prefix + ".stats", training_prefix + ".gts", 
                           ignore_gts = False, restrict_to_model_name = restrict_to_model_name, gt_labels = gt_labels)
    if restrict_to_model_name in (None, "ALL", "INS_TO_DUP", "INS_TO_DUP_LARGE"):
        ins_to_dup_vcf_training_data, ins_to_dup_vcf_training_gts, _, ins_to_dup_vcf_training_exact = \
            features.parse_vcf(training_prefix + ".INS_TO_DUP.vcf.gz", training_prefix + ".stats",
                                training_prefix + ".gts", ignore_gts = False, restrict_to_model_name = restrict_to_model_name,
                                gt_labels = gt_labels)

        for model_name in ("INS_TO_DUP", "INS_TO_DUP_LARGE"):
            if model_name in ins_to_dup_vcf_training_data:
                vcf_training_data[model_name] = ins_to_dup_vcf_training_data[model_name]
                vcf_training_gts[model_name] = ins_to_dup_vcf_training_gts[model_name]
                vcf_training_exact[model_name] = ins_to_dup_vcf_training_exact[model_name]

    return vcf_training_data, vcf_training_gts, vcf_training_exact

if __name__ == '__main__':
    import multiprocessing
    multiprocessing.set_start_method("spawn", force=True)  # optional but recommended

    training_prefixes = cmd_args.training_prefixes.split(",")
    restrict_to_model_name = None if cmd_args.model_name == "ALL" else cmd_args.model_name
    training_manifest = build_training_manifest(
        training_prefixes,
        model_name=cmd_args.model_name,
        cross_species=cmd_args.cross_species,
    )

    processed_samples = [None] * len(training_prefixes)
    with ProcessPoolExecutor(max_workers=cmd_args.threads) as executor:
        future_to_sample = {
            executor.submit(process_vcf, prefix, restrict_to_model_name): (sample_index, prefix)
            for sample_index, prefix in enumerate(training_prefixes)
        }
        for future in as_completed(future_to_sample):
            sample_index, prefix = future_to_sample[future]
            try:
                processed_samples[sample_index] = future.result()
            except Exception as e:
                raise RuntimeError(f"Failed while processing training prefix: {prefix}") from e

    for vcf_training_data, vcf_training_gts, vcf_training_exact in processed_samples:
        for model in vcf_training_data:
            training_data[model].append(vcf_training_data[model])
            training_gts[model].append(vcf_training_gts[model])
            training_exact[model].append(vcf_training_exact[model])

    for model in training_data:
        training_data[model] = np.concatenate(training_data[model])
        training_gts[model] = np.concatenate(training_gts[model])
        training_exact[model] = np.concatenate(training_exact[model])

    yes_or_no_outdir = os.path.join(cmd_args.outdir, "yes_or_no")
    os.makedirs(yes_or_no_outdir, exist_ok=True)

    gts_outdir = os.path.join(cmd_args.outdir, "gts")
    os.makedirs(gts_outdir, exist_ok=True)

    eref_outdir = os.path.join(cmd_args.outdir, "eref")
    os.makedirs(eref_outdir, exist_ok=True)

    exact_outdir = os.path.join(cmd_args.outdir, "exact")
    os.makedirs(exact_outdir, exist_ok=True)

    manifest_path = write_training_manifest(cmd_args.outdir, training_manifest)
    print(f"Wrote training manifest to {manifest_path}")

    for model_name in training_data:
        if cmd_args.model_name != "ALL" and model_name != cmd_args.model_name:
            continue

        start_time = time.time()

        all_training_data = training_data[model_name]

        features_names = features.Features.get_feature_names(model_name)
        keep_indices = compute_keep_indices(all_training_data)
        if len(keep_indices) == 0:
            raise RuntimeError(f"Model {model_name} has no usable features after pruning.")
        features_names = [features_names[i] for i in keep_indices]

        known_gt_mask = training_gts[model_name] != "./."
        training_data[model_name] = training_data[model_name][known_gt_mask]
        training_gts[model_name] = training_gts[model_name][known_gt_mask]
        training_exact[model_name] = training_exact[model_name][known_gt_mask]
        training_data[model_name] = training_data[model_name][:, keep_indices]
        if len(training_data[model_name]) == 0:
            raise RuntimeError(f"Model {model_name} has no non-missing training examples.")

        end_time = time.time()
        print(f"Preprocessing for model {model_name} took {end_time - start_time} seconds")

        training_labels = np.array([0 if x == "0/0" else 1 for x in training_gts[model_name]])
        unique_labels = np.unique(training_labels)
        if len(unique_labels) == 1:
            raise RuntimeError(f"Only one label ({unique_labels[0]}) present in yes/no training data for model {model_name}. Cannot train the first-stage classifier.")

        start_time = time.time()
        classifier = make_classifier()
        classifier.fit(training_data[model_name], training_labels)

        importances = classifier.feature_importances_
        indices = np.argsort(importances)[::-1]
        with open(os.path.join(yes_or_no_outdir, model_name + ".importance.txt"), 'w') as f:
            for i in range(len(features_names)):
                f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

        model_fname = os.path.join(yes_or_no_outdir, model_name + ".ubj")
        classifier.save_model(model_fname)
        write_features_file(model_fname, features_names)

        known_positive_mask = features.gt_is_known_positive_array(training_gts[model_name])
        positive_training_data = training_data[model_name][known_positive_mask]
        positive_training_labels = np.array([features.gt_stage_label(gt) for gt in training_gts[model_name][known_positive_mask]], dtype=np.int32)

        if len(positive_training_data) == 0:
            raise RuntimeError(f"No known positive training examples found for model {model_name}. Cannot train the GT-stage classifier.")

        unique_labels = np.unique(positive_training_labels)
        if len(unique_labels) == 1:
            raise RuntimeError(f"Only one label ({unique_labels[0]}) present in positive training data for model {model_name}. Cannot train the GT-stage classifier.")

        unique, counts = np.unique(positive_training_labels, return_counts=True)

        classifier = make_classifier()
        classifier.fit(positive_training_data, positive_training_labels)

        model_fname = os.path.join(gts_outdir, model_name + ".ubj")
        importance_fname = os.path.join(gts_outdir, model_name + ".importance.txt")
        remove_model_artifacts(model_fname, importance_fname)

        importances = classifier.feature_importances_
        indices = np.argsort(importances)[::-1]
        with open(importance_fname, 'w') as f:
            for i in range(len(features_names)):
                f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

        classifier.save_model(model_fname)
        write_features_file(model_fname, features_names)

        known_het_mask = features.gt_is_known_het_array(training_gts[model_name])
        eref_training_data = training_data[model_name][known_het_mask]
        eref_training_labels = np.array(
            [features.gt_eref_label(gt) for gt in training_gts[model_name][known_het_mask]],
            dtype=np.int32
        )

        eref_stage_classes = np.unique(eref_training_labels)
        model_fname = os.path.join(eref_outdir, model_name + ".ubj")
        importance_fname = os.path.join(eref_outdir, model_name + ".importance.txt")
        remove_model_artifacts(model_fname, importance_fname)

        if len(eref_stage_classes) >= 2:
            write_features_file(model_fname, features_names)
            classifier = make_classifier()
            classifier.fit(eref_training_data, eref_training_labels)

            importances = classifier.feature_importances_
            indices = np.argsort(importances)[::-1]
            with open(importance_fname, 'w') as f:
                for i in range(len(features_names)):
                    f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

            classifier.save_model(model_fname)

        if model_name == "HP" or model_name.startswith("DEL") or model_name.startswith("INS"):
            exact_start_time = time.time()
            positive_mask = features.gt_has_alt_array(training_gts[model_name])
            exact_training_data = training_data[model_name][positive_mask]
            exact_training_labels = training_exact[model_name][positive_mask].astype(int)

            if len(exact_training_data) == 0:
                raise RuntimeError(f"No positive {model_name} examples found. Cannot train exact-stage classifier.")

            unique_labels = np.unique(exact_training_labels)
            if len(unique_labels) == 1:
                raise RuntimeError(
                    f"Only one exact label ({unique_labels[0]}) present in {model_name} training data. "
                    "Cannot train the exact-stage classifier."
                )

            n_exact = np.sum(exact_training_labels == 1)
            n_inexact = np.sum(exact_training_labels == 0)
	            
            classifier = make_classifier()
            classifier.fit(exact_training_data, exact_training_labels)

            importances = classifier.feature_importances_
            indices = np.argsort(importances)[::-1]
            with open(os.path.join(exact_outdir, model_name + ".importance.txt"), "w") as f:
                for i in range(len(features_names)):
                    f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

            model_fname = os.path.join(exact_outdir, model_name + ".ubj")
            classifier.save_model(model_fname)
            write_features_file(model_fname, features_names)
            exact_end_time = time.time()

        end_time = time.time()
        print(f"Training for model {model_name} took {end_time - start_time} seconds")
