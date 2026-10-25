import argparse, time
import numpy as np
import features
import os
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
import xgboost as xgb

cmd_parser = argparse.ArgumentParser(description='Train ML model.')
cmd_parser.add_argument('training_prefixes', help='Prefix of the training VCF and FP files.')
cmd_parser.add_argument('outdir')
cmd_parser.add_argument('--model_name', default='ALL', help='Restrict to this model. CNU trains only the DUP copy number model.')
cmd_parser.add_argument('--threads', type=int, default=1, help='Number of threads to use for training.')
cmd_parser.add_argument('--cross-species', action='store_true', help='Use cross-species model.')
cmd_args = cmd_parser.parse_args()

training_data, training_gts = defaultdict(list), defaultdict(list)
training_cn_unresolved = defaultdict(list)

# models of the DUP records (INS_TO_DUP excluded) used to train the copy number model
CNU_MODELS = ("DUP", "DUP_LARGE", "DUP_LARGE_NOEXL")

def process_vcf(training_prefix):
    vcf_training_data, vcf_training_gts, vcf_training_cn_unresolved, _ = \
        features.parse_vcf(training_prefix + ".vcf.gz", training_prefix + ".stats", training_prefix + ".gts", 
                           tolerate_no_gts = False)
    ins_to_dup_vcf_training_data, ins_to_dup_vcf_training_gts, _, _ = \
        features.parse_vcf(training_prefix + ".INS_TO_DUP.vcf.gz", training_prefix + ".stats",
                            training_prefix + ".gts", tolerate_no_gts = False)
    
    vcf_training_data["INS_TO_DUP"] = ins_to_dup_vcf_training_data["INS_TO_DUP"]
    vcf_training_gts["INS_TO_DUP"] = ins_to_dup_vcf_training_gts["INS_TO_DUP"]
    vcf_training_data["INS_TO_DUP_LARGE"] = ins_to_dup_vcf_training_data["INS_TO_DUP_LARGE"]
    vcf_training_gts["INS_TO_DUP_LARGE"] = ins_to_dup_vcf_training_gts["INS_TO_DUP_LARGE"]
    return vcf_training_data, vcf_training_gts, vcf_training_cn_unresolved

if __name__ == '__main__':
    import multiprocessing
    multiprocessing.set_start_method("spawn", force=True)  # optional but recommended

    training_prefixes = cmd_args.training_prefixes.split(",")
    with ProcessPoolExecutor(max_workers=cmd_args.threads) as executor:
        futures = [executor.submit(process_vcf, prefix) for prefix in training_prefixes]
        # concatenate in the order of the prefixes, so that training is reproducible
        for future in futures:
            vcf_training_data, vcf_training_gts, vcf_training_cn_unresolved = future.result()
            for model in vcf_training_data:
                training_data[model].append(vcf_training_data[model])
                training_gts[model].append(vcf_training_gts[model])
                if model in CNU_MODELS:
                    training_cn_unresolved[model].append(vcf_training_cn_unresolved[model])

    for model in training_data:
        training_data[model] = np.concatenate(training_data[model])
        training_gts[model] = np.concatenate(training_gts[model])
    for model in training_cn_unresolved:
        training_cn_unresolved[model] = np.concatenate(training_cn_unresolved[model])

    yes_or_no_outdir = os.path.join(cmd_args.outdir, "yes_or_no")
    os.makedirs(yes_or_no_outdir, exist_ok=True)

    gts_outdir = os.path.join(cmd_args.outdir, "gts")
    os.makedirs(gts_outdir, exist_ok=True)

    for model_name in training_data:
        if cmd_args.model_name != "ALL" and model_name != cmd_args.model_name:
            continue

        training_labels = np.array([0 if x == "0/0" else 1 for x in training_gts[model_name]])
        if cmd_args.cross_species:
            classifier = xgb.XGBClassifier(n_estimators=50, max_depth=7, min_child_weight=42, learning_rate=0.1, n_jobs=cmd_args.threads, random_state=42, tree_method='hist')
        else:
            classifier = xgb.XGBClassifier(n_estimators=1000, max_depth=8, min_child_weight=16, learning_rate=0.1, n_jobs=cmd_args.threads, random_state=42, tree_method='hist')

        start_time = time.time()
        classifier.fit(training_data[model_name], training_labels)

        features_names = features.Features.get_feature_names(model_name)
        importances = classifier.feature_importances_
        indices = np.argsort(importances)[::-1]
        with open(os.path.join(yes_or_no_outdir, model_name + ".importance.txt"), 'w') as f:
            for i in range(len(features_names)):
                f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

        model_fname = os.path.join(yes_or_no_outdir, model_name + ".ubj")
        classifier.save_model(model_fname)

        positive_training_data = training_data[model_name][(training_gts[model_name] == "0/1") | (training_gts[model_name] == "1/1")]
        positive_training_labels = np.array([1 if x == "1/1" else 0 for x in training_gts[model_name] if x == "0/1" or x == "1/1"])

        unique_labels = np.unique(positive_training_labels)
        if len(unique_labels) == 1:
            print(f"Only one label present in positive training data for model {model_name}.")
            first_entry = positive_training_data[0:1]
            first_label = positive_training_labels[0]
            opposite_label = 1 - first_label  
            positive_training_data = np.concatenate([positive_training_data, first_entry], axis=0)
            positive_training_labels = np.concatenate([positive_training_labels, np.array([opposite_label])])

        unique, counts = np.unique(positive_training_labels, return_counts=True)

        classifier.fit(positive_training_data, positive_training_labels)

        importances = classifier.feature_importances_
        indices = np.argsort(importances)[::-1]
        with open(os.path.join(gts_outdir, model_name + ".importance.txt"), 'w') as f:
            for i in range(len(features_names)):
                f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))

        model_fname = os.path.join(gts_outdir, model_name + ".ubj")
        classifier.save_model(model_fname)

        end_time = time.time()
        print(f"Training for model {model_name} took {end_time - start_time} seconds")

    if cmd_args.model_name in ("ALL", "CNU"):
        # copy number model: DUPs present in the benchmark, labelled 1 if their copy number is unresolved (the benchmark
        # insertion is longer than a single copy of the DUP) and 0 otherwise
        cnu_training_data, cnu_training_labels = [], []
        for model_name in CNU_MODELS:
            if model_name not in training_cn_unresolved:
                continue
            mask = np.isin(training_gts[model_name], ("0/1", "1/1", "./1")) & (training_cn_unresolved[model_name] >= 0)
            cnu_training_data.append(training_data[model_name][mask])
            cnu_training_labels.append(training_cn_unresolved[model_name][mask])
        cnu_training_data, cnu_training_labels = np.concatenate(cnu_training_data), np.concatenate(cnu_training_labels)

        start_time = time.time()
        classifier = xgb.XGBClassifier(n_estimators=600, max_depth=6, learning_rate=0.05, subsample=0.8, colsample_bytree=0.6, n_jobs=cmd_args.threads, random_state=0, tree_method='hist')
        classifier.fit(cnu_training_data, cnu_training_labels)

        cnu_outdir = os.path.join(cmd_args.outdir, "cnu")
        os.makedirs(cnu_outdir, exist_ok=True)
        features_names = features.Features.get_feature_names("DUP")
        importances = classifier.feature_importances_
        indices = np.argsort(importances)[::-1]
        with open(os.path.join(cnu_outdir, "DUP.importance.txt"), 'w') as f:
            for i in range(len(features_names)):
                f.write("%d. %s (%f)\n" % (i + 1, features_names[indices[i]], importances[indices[i]]))
        classifier.save_model(os.path.join(cnu_outdir, "DUP.ubj"))

        end_time = time.time()
        print(f"Training for the copy number model took {end_time - start_time} seconds ({len(cnu_training_labels)} DUPs, {int(cnu_training_labels.sum())} with unresolved copy number)")
