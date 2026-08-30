import argparse, time
import numpy as np
import os
import pysam
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import xgboost as xgb
from run_classifier import Classifier
from training_manifest import build_training_manifest, write_training_manifest

cmd_parser = argparse.ArgumentParser(description='Train ML model.')
cmd_parser.add_argument('training_prefixes', help='Prefix of the training VCF and FP files.')
cmd_parser.add_argument('outdir')
cmd_parser.add_argument('--model_name', default='ALL', help='Restrict to this model.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use for training.')
cmd_parser.add_argument('--cross-species', action='store_true', help='Use cross-species model.')
cmd_args = cmd_parser.parse_args()

training_data, training_gts, training_exact = defaultdict(list), defaultdict(list), defaultdict(list)

MODEL_NAMES = (
    "HP", "DEL", "DEL_NOEXL", "DEL_LARGE", "DEL_LARGE_NOEXL", "DUP", "DUP_LARGE", "DUP_LARGE_NOEXL",
    "INS", "INS_TO_DUP", "INS_TO_DUP_LARGE",
)

def log(message):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)

def gt_alleles(gt):
    return gt.replace("|", "/").split("/")

def gt_has_alt(gt):
    return "1" in gt_alleles(gt)

def gt_is_known_positive(gt):
    alleles = gt_alleles(gt)
    return "1" in alleles and "." not in alleles

def gt_is_hom_alt(gt):
    alleles = gt_alleles(gt)
    return len(alleles) > 0 and all(allele == "1" for allele in alleles)

def gt_is_known_het(gt):
    alleles = gt_alleles(gt)
    return "." not in alleles and "1" in alleles and not gt_is_hom_alt(gt)

def gt_stage_label(gt):
    if gt == "0/1" or gt == "1/2":
        return 0
    if gt == "1/1":
        return 1
    raise RuntimeError(f"Unsupported GT-stage label: {gt}")

def gt_eref_label(gt):
    if gt == "0/1":
        return 1
    if gt == "1/2":
        return 0
    raise RuntimeError(f"Unsupported EREFA-stage label: {gt}")

def gt_has_alt_array(gts):
    return np.array([gt_has_alt(gt) for gt in gts])

def gt_is_known_positive_array(gts):
    return np.array([gt_is_known_positive(gt) for gt in gts])

def gt_is_known_het_array(gts):
    return np.array([gt_is_known_het(gt) for gt in gts])

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

def write_extractor_feature_definitions(feature_dir, feature_names):
    for model_name in MODEL_NAMES:
        with open(os.path.join(feature_dir, model_name + ".features"), "w") as feature_file:
            for feature_name in feature_names:
                feature_file.write(feature_name + "\n")

def load_record_ids_by_key(vcf_fname):
    record_ids = dict()
    record_idxs = defaultdict(int)
    with pysam.VariantFile(vcf_fname) as vcf_reader:
        for record in vcf_reader.fetch():
            record_idx = record_idxs[record.rid]
            record_idxs[record.rid] += 1
            record_ids[(record.rid << 32) | record_idx] = record.id
    return record_ids

def extract_vcf_features(vcf_fname, stats_fname, gt_labels, feature_dir, temp_dir, extractor_threads, restrict_to_model_name = None):
    start_time = time.time()
    log(f"Extracting features from {vcf_fname} with {extractor_threads} thread(s)")
    feature_bundle_fname = os.path.join(temp_dir, os.path.basename(vcf_fname) + ".features.bin")
    command = [Classifier.get_extractor_path(), vcf_fname, stats_fname, feature_dir, feature_bundle_fname, str(extractor_threads)]
    extractor = subprocess.Popen(command)
    while True:
        try:
            returncode = extractor.wait(timeout=60)
            break
        except subprocess.TimeoutExpired:
            log(f"Still extracting features from {vcf_fname} ({time.time()-start_time:.0f}s elapsed)")
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, command)
    vcf_training_data, _, record_keys_by_model = Classifier.load_feature_bundle(feature_bundle_fname)
    record_ids = load_record_ids_by_key(vcf_fname)
    vcf_training_gts, vcf_training_exact = dict(), dict()
    for model_name, record_keys in record_keys_by_model.items():
        if restrict_to_model_name not in (None, "ALL") and model_name != restrict_to_model_name:
            continue
        model_gts, model_exact = [], []
        for record_key in record_keys:
            record_key = int(record_key)
            if record_key not in record_ids:
                raise RuntimeError(f"Could not match C++ feature row {record_key} to a record in {vcf_fname}")
            record_id = record_ids[record_key]
            label = gt_labels.get(record_id)
            if label is None:
                raise RuntimeError(f"Missing GT label for record {record_id} in {vcf_fname}")
            gt, exact = label
            model_gts.append(gt)
            model_exact.append(exact)
        vcf_training_gts[model_name] = np.array(model_gts)
        vcf_training_exact[model_name] = np.array(model_exact)
    if restrict_to_model_name not in (None, "ALL"):
        vcf_training_data = {model_name: values for model_name, values in vcf_training_data.items() if model_name == restrict_to_model_name}
    n_variants = sum(len(model_data) for model_data in vcf_training_data.values())
    log(f"Extracted {n_variants} variants from {vcf_fname} across {len(vcf_training_data)} model(s) in {time.time()-start_time:.1f}s")
    return vcf_training_data, vcf_training_gts, vcf_training_exact

def select_gt(gt1, gt2):
    if gt1 == "./." and gt2 != "./.":
        return gt2
    if gt1 != "./." and gt2 == "./.":
        return gt1
    elif gt1 == gt2:
        return gt1
    else:
        return "./."

def read_gt_labels(file_path):
    if not os.path.exists(file_path):
        raise RuntimeError(f"Genotype labels file not found: {file_path}")
    gt_labels = dict()
    with open(file_path, 'r') as file:
        for line in file:
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 3:
                raise RuntimeError(f"Malformed genotype labels line in {file_path}: {' '.join(fields)}")
            id, gt, exact = fields[0], fields[1], int(fields[2])
            if id not in gt_labels:
                gt_labels[id] = (gt, exact)
            else:
                prev_gt, prev_exact = gt_labels[id]
                gt_labels[id] = (
                    select_gt(prev_gt, gt),
                    max(prev_exact, exact),
                )
    return gt_labels

def process_vcf(training_prefix, feature_dir, temp_root, extractor_threads, restrict_to_model_name = None):
    start_time = time.time()
    log(f"Starting training prefix {training_prefix}")
    gt_labels = read_gt_labels(training_prefix + ".gts")
    log(f"Loaded {len(gt_labels)} genotype labels for {training_prefix}")
    with tempfile.TemporaryDirectory(prefix=".surveyor-training-features-", dir=temp_root) as temp_dir:
        vcf_training_data, vcf_training_gts, vcf_training_exact = extract_vcf_features(training_prefix + ".vcf.gz", training_prefix + ".stats", gt_labels, feature_dir, temp_dir, extractor_threads, restrict_to_model_name)
        if restrict_to_model_name in (None, "ALL", "INS_TO_DUP", "INS_TO_DUP_LARGE"):
            ins_to_dup_vcf_training_data, ins_to_dup_vcf_training_gts, ins_to_dup_vcf_training_exact = extract_vcf_features(training_prefix + ".INS_TO_DUP.vcf.gz", training_prefix + ".stats", gt_labels, feature_dir, temp_dir, extractor_threads, restrict_to_model_name)

            for model_name in ("INS_TO_DUP", "INS_TO_DUP_LARGE"):
                if model_name in ins_to_dup_vcf_training_data:
                    vcf_training_data[model_name] = ins_to_dup_vcf_training_data[model_name]
                    vcf_training_gts[model_name] = ins_to_dup_vcf_training_gts[model_name]
                    vcf_training_exact[model_name] = ins_to_dup_vcf_training_exact[model_name]

    n_variants = sum(len(model_data) for model_data in vcf_training_data.values())
    log(f"Processed training prefix {training_prefix}: {n_variants} variants in {time.time()-start_time:.1f}s")
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

    os.makedirs(cmd_args.outdir, exist_ok=True)
    default_feature_names = Classifier.get_default_feature_names()
    log(f"Loaded {len(default_feature_names)} default feature names")
    n_workers = min(cmd_args.threads, len(training_prefixes))
    extractor_threads = max(1, cmd_args.threads // n_workers)
    log(f"Processing {len(training_prefixes)} training prefixes with {n_workers} worker(s) and {extractor_threads} extractor thread(s) per worker")
    processed_samples = [None] * len(training_prefixes)
    with tempfile.TemporaryDirectory(prefix=".surveyor-training-feature-defs-", dir=cmd_args.outdir) as feature_dir:
        write_extractor_feature_definitions(feature_dir, default_feature_names)
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_sample = {
                executor.submit(process_vcf, prefix, feature_dir, cmd_args.outdir, extractor_threads, restrict_to_model_name): (sample_index, prefix)
                for sample_index, prefix in enumerate(training_prefixes)
            }
            log(f"Submitted all {len(future_to_sample)} training prefixes")
            n_completed = 0
            for future in as_completed(future_to_sample):
                sample_index, prefix = future_to_sample[future]
                try:
                    processed_samples[sample_index] = future.result()
                    n_completed += 1
                    log(f"Collected training prefix {prefix} ({n_completed}/{len(training_prefixes)})")
                except Exception as e:
                    raise RuntimeError(f"Failed while processing training prefix: {prefix}") from e
            del future
            future_to_sample.clear()
            log("All training prefixes collected; shutting down worker processes")

    log("Worker processes shut down; organizing extracted arrays by model")
    for sample_index in range(len(processed_samples)):
        vcf_training_data, vcf_training_gts, vcf_training_exact = processed_samples[sample_index]
        for model in vcf_training_data:
            training_data[model].append(vcf_training_data[model])
            training_gts[model].append(vcf_training_gts[model])
            training_exact[model].append(vcf_training_exact[model])
        processed_samples[sample_index] = None
    del processed_samples, vcf_training_data, vcf_training_gts, vcf_training_exact

    for model in training_data:
        n_variants = sum(len(sample_data) for sample_data in training_data[model])
        input_size_gib = sum(sample_data.nbytes for sample_data in training_data[model])/2**30
        concatenate_start_time = time.time()
        log(f"Aggregating model {model}: {n_variants} variants from {len(training_data[model])} samples ({input_size_gib:.2f} GiB of feature data)")
        training_data[model] = np.concatenate(training_data[model])
        training_gts[model] = np.concatenate(training_gts[model])
        training_exact[model] = np.concatenate(training_exact[model])
        log(f"Aggregated model {model} in {time.time()-concatenate_start_time:.1f}s")

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
        log(f"Preparing model {model_name}: {all_training_data.shape[0]} variants and {all_training_data.shape[1]} raw features")

        features_names = default_feature_names
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
        log(f"Training yes/no classifier for {model_name} on {len(training_data[model_name])} variants and {len(features_names)} features")
        classifier = make_classifier()
        classifier.fit(training_data[model_name], training_labels)
        log(f"Finished yes/no classifier for {model_name} in {time.time()-start_time:.1f}s")

        importances = classifier.feature_importances_
        indices = np.argsort(importances)[::-1]
        with open(os.path.join(yes_or_no_outdir, model_name + ".importance.txt"), 'w') as f:
            for i in range(len(features_names)):
                f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

        model_fname = os.path.join(yes_or_no_outdir, model_name + ".ubj")
        classifier.save_model(model_fname)
        write_features_file(model_fname, features_names)

        known_positive_mask = gt_is_known_positive_array(training_gts[model_name])
        positive_training_data = training_data[model_name][known_positive_mask]
        positive_training_labels = np.array([gt_stage_label(gt) for gt in training_gts[model_name][known_positive_mask]], dtype=np.int32)

        if len(positive_training_data) == 0:
            raise RuntimeError(f"No known positive training examples found for model {model_name}. Cannot train the GT-stage classifier.")

        unique_labels = np.unique(positive_training_labels)
        if len(unique_labels) == 1:
            raise RuntimeError(f"Only one label ({unique_labels[0]}) present in positive training data for model {model_name}. Cannot train the GT-stage classifier.")

        unique, counts = np.unique(positive_training_labels, return_counts=True)

        gt_start_time = time.time()
        log(f"Training GT classifier for {model_name} on {len(positive_training_data)} positive variants")
        classifier = make_classifier()
        classifier.fit(positive_training_data, positive_training_labels)
        log(f"Finished GT classifier for {model_name} in {time.time()-gt_start_time:.1f}s")

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

        known_het_mask = gt_is_known_het_array(training_gts[model_name])
        eref_training_data = training_data[model_name][known_het_mask]
        eref_training_labels = np.array(
            [gt_eref_label(gt) for gt in training_gts[model_name][known_het_mask]],
            dtype=np.int32
        )

        eref_stage_classes = np.unique(eref_training_labels)
        model_fname = os.path.join(eref_outdir, model_name + ".ubj")
        importance_fname = os.path.join(eref_outdir, model_name + ".importance.txt")
        remove_model_artifacts(model_fname, importance_fname)

        if len(eref_stage_classes) >= 2:
            write_features_file(model_fname, features_names)
            eref_start_time = time.time()
            log(f"Training EREFA classifier for {model_name} on {len(eref_training_data)} heterozygous variants")
            classifier = make_classifier()
            classifier.fit(eref_training_data, eref_training_labels)
            log(f"Finished EREFA classifier for {model_name} in {time.time()-eref_start_time:.1f}s")

            importances = classifier.feature_importances_
            indices = np.argsort(importances)[::-1]
            with open(importance_fname, 'w') as f:
                for i in range(len(features_names)):
                    f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

            classifier.save_model(model_fname)

        if model_name == "HP" or model_name.startswith("DEL") or model_name.startswith("INS"):
            exact_start_time = time.time()
            positive_mask = gt_has_alt_array(training_gts[model_name])
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
	            
            log(f"Training exact classifier for {model_name} on {len(exact_training_data)} positive variants ({n_exact} exact, {n_inexact} inexact)")
            classifier = make_classifier()
            classifier.fit(exact_training_data, exact_training_labels)
            log(f"Finished exact classifier for {model_name} in {time.time()-exact_start_time:.1f}s")

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
