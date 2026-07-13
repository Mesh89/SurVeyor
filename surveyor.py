import sys, os, argparse, hashlib, pysam, timeit, shutil
from run_classifier import Classifier

VERSION = "0.13"

def valid_min_sv_size(arg):
    try:
        sv_size = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Value must be an integer.")
    if sv_size < 1:
        raise argparse.ArgumentTypeError("Value must be at least 1.")
    return sv_size

parser = argparse.ArgumentParser(description='SurVeyor, an SV caller.')
parser.add_argument('--version', action='version', version="SurVeyor v%s" % VERSION, help='Print version number.')

subparsers = parser.add_subparsers(dest='command', help='Commands', required=True)

common_parser = argparse.ArgumentParser(add_help=False)
common_parser.add_argument('--threads', type=int, default=1, help='Number of threads to be used.')

call_genotype_shared_options_parser = argparse.ArgumentParser(add_help=False)
call_genotype_shared_options_parser.add_argument('--seed', type=int, default=0, help='Seed for random sampling of genomic positions.')
call_genotype_shared_options_parser.add_argument('--max-seq-error', type=float, default=0.04, help='Max sequencing error admissible on the platform used.')
call_genotype_shared_options_parser.add_argument('--min-sv-size', type=valid_min_sv_size, default=50, help='Min SV size.')
call_genotype_shared_options_parser.add_argument('--min-clip-len', type=int, default=10, help='Min length for a clip to be used.')
call_genotype_shared_options_parser.add_argument('--sampling-regions', help='File in BED format containing regions to be used to estimate statistics such as depth.')
call_genotype_shared_options_parser.add_argument('--min-stable-mapq', type=int, default=20, help='Minimum MAPQ for a stable read.')
call_genotype_shared_options_parser.add_argument('--high-confidence-mapq', type=int, default=60, help='MAPQ threshold above which a read is considered high-confidence.')
call_genotype_shared_options_parser.add_argument('--skip-bam-ref-validation', action='store_true',
                        help='Skip BAM/reference compatibility checks. Use this at your own risk.')
call_genotype_shared_options_parser.add_argument('--force-workdir', action='store_true',
                        help='Allow execution in a non-empty workdir. Use this with caution.')
call_genotype_shared_options_parser.add_argument('--per-contig-stats', action='store_true',
                        help='Depth statistics are computed separately for each contig. Useful when one or more of the target contigs are expected to have '
                        'dramatically different depth than others. Otherwise, it is not recommended to use this option.')
call_genotype_shared_options_parser.add_argument('--two-pass', action='store_true', help='Activate two-pass genotyping.')
call_genotype_shared_options_parser.add_argument('--min-diff-hsr', type=int, default=3, help='Minimum number of differences with the reference \
                        (considered as number of insertions, deletions and mismatches) for a read to be considered a hidden split read.')
call_genotype_shared_options_parser.add_argument('--tr-bed', help='BED file with tandem repetitive regions. If provided, it will be used for a more aggrestive duplicate removal.')
call_genotype_shared_options_parser.add_argument('--haploid_contigs', '--haploid-contigs', default='',
                        help='Comma-separated contigs to treat as haploid (for example, chrX,chrY). '
                             'Retained positive variants on these contigs are written as 1/1 and may not overlap.')

call_parser = subparsers.add_parser('call', parents=[common_parser, call_genotype_shared_options_parser], help='Call SVs denovo.')
call_parser.add_argument('bam_file', help='Input bam file.')
call_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
call_parser.add_argument('reference', help='Reference genome in FASTA format.')
call_parser.add_argument('--generate-training-data', action='store_true', help='Generate data needed to train a genotyping ML model.')
call_parser.add_argument('--ml-model', help='Path to the ML model to be used for filtering and genotyping.')
call_parser.add_argument('--samplename', default='', help='Name of the sample to be used in the VCF output.'
                                                         'If not provided, the basename of the bam/cram file will be used, up until the first \'.\'')
call_parser.add_argument('--max-trans-size', type=int, default=10000, help='Maximum size of the transpositions which '
                                                                          'SurVeyor will predict when only one side is available.')

genotype_parser = subparsers.add_parser('genotype', parents=[common_parser, call_genotype_shared_options_parser], help='Genotype SVs.')
genotype_parser.add_argument('in_vcf_file', help='Input VCF file.')
genotype_parser.add_argument('bam_file', help='Input bam file.')
genotype_parser.add_argument('workdir', help='Working directory for Surveyor to use.')
genotype_parser.add_argument('reference', help='Reference genome in FASTA format.')
genotype_parser.add_argument('--generate-training-data', action='store_true', help='Generate data needed to train a genotyping ML model.')
genotype_parser.add_argument('--ml-model', help='Path to the ML model to be used for genotyping.')
genotype_parser.add_argument('--samplename', default='', help='Name of the sample to be used in the VCF output.'
                                                         'If not provided, the basename of the bam/cram file will be used,'
                                                         'up until the first \'.\'')

generate_training_data_parser = subparsers.add_parser('generate-training-data', parents=[common_parser], help='Generate training data for ML model starting from the workdir of a call run.')
generate_training_data_parser.add_argument('benchmark_vcf', help='Benchmark VCF file used to assign labels to the training data.')
generate_training_data_parser.add_argument('workdir', help='Working directory used by the call command.')
generate_training_data_parser.add_argument('outdir', help='Output directory to store the training data.')
generate_training_data_parser.add_argument('simple_repeat_bed', help='BED file with simple repeat regions.')
generate_training_data_parser.add_argument('reference', help='Reference genome in FASTA format.')
generate_training_data_parser.add_argument('samplename', help='Name of the sample used in the call command.')
generate_training_data_parser.add_argument('--unreliable-gts', help='File with list of ID of variants with unreliable genotypes, to be set to ./1 if found')
generate_training_data_parser.add_argument('--restrict-to-bed', help='Restrict the training data VCF to indels fully contained in an interval from this BED file.')
generate_training_data_parser.add_argument('--use-reassigned-training-data', action='store_true',
                                           help='Use reassigned training data if present. Requires training-data.reassigned.vcf.gz '
                                                'and training-data.reassigned.INS_TO_DUP.vcf.gz in the workdir.')

cmd_args = parser.parse_args()

def run_cmd(cmd, error_msg=None):
    start_time = timeit.default_timer()
    print("Executing:", cmd, flush=True)
    return_code = os.system(cmd)
    if return_code != 0:
        if error_msg:
            print(error_msg, flush=True)
        else:
            print("Error executing:", cmd, flush=True)
            print("Return code:", return_code, flush=True)
        exit(1)
    elapsed = timeit.default_timer() - start_time
    print(cmd, "was run in %.2f seconds" % elapsed, flush=True)

SURVEYOR_PATH = os.path.dirname(os.path.realpath(__file__))

def mkdir(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def rmdir(dirname):
    if os.path.exists(dirname):
        shutil.rmtree(dirname)

def mkdir_clean(dirname):
    rmdir(dirname)
    mkdir(dirname)

def fail(message):
    print(message, flush=True)
    exit(1)

def get_max_is_from_stats(stats_fname):
    with open(stats_fname) as stats_file:
        for line in stats_file:
            tokens = line.strip().split()
            if len(tokens) >= 3 and tokens[0] == "max_is" and tokens[1] == ".":
                return int(tokens[2])
    fail("Error: max_is not found in %s." % stats_fname)

def get_stat_from_stats(stats_fname, stat_name, subname="."):
    with open(stats_fname) as stats_file:
        for line in stats_file:
            tokens = line.strip().split()
            if len(tokens) >= 3 and tokens[0] == stat_name and tokens[1] == subname:
                return int(tokens[2])
    fail("Error: %s %s not found in %s." % (stat_name, subname, stats_fname))

def get_config_value(config_fname, key):
    if os.path.exists(config_fname):
        with open(config_fname) as config_file:
            for line in config_file:
                tokens = line.strip().split()
                if len(tokens) >= 2 and tokens[0] == key:
                    return tokens[1]
    fail("Error: %s not found in %s." % (key, config_fname))

def validate_workdir_for_execution(workdir, force_workdir):
    if not os.path.exists(workdir):
        return
    if not os.path.isdir(workdir):
        fail("Error: workdir %s exists but is not a directory." % workdir)
    if os.listdir(workdir) and not force_workdir:
        fail("Error: workdir %s already exists and is not empty. Use a new workdir, empty it first, or rerun with --force-workdir to force the use of a non-empty workdir." % workdir)

def compute_reference_contig_md5(reference_fasta, contig_name, chunk_size=1000000):
    contig_len = reference_fasta.get_reference_length(contig_name)
    md5 = hashlib.md5()
    for start in range(0, contig_len, chunk_size):
        end = min(contig_len, start + chunk_size)
        md5.update(reference_fasta.fetch(contig_name, start, end).upper().encode())
    return md5.hexdigest()

def parse_bed_regions(bed_fname):
    regions = []
    try:
        with open(bed_fname) as bed_file:
            for line_no, line in enumerate(bed_file, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue

                tokens = stripped.split()
                if len(tokens) < 3:
                    fail("Error: invalid format in BED file %s at line %d. Expected at least 3 BED fields: 'contig start end'." %
                         (bed_fname, line_no))

                try:
                    start = int(tokens[1])
                    end = int(tokens[2])
                except ValueError:
                    fail("Error: invalid coordinates in BED file %s at line %d." %
                         (bed_fname, line_no))

                if start < 0 or end < 0:
                    fail("Error: BED file %s line %d has a negative coordinate." %
                         (bed_fname, line_no))
                if start >= end:
                    fail("Error: BED file %s line %d has start >= end." %
                         (bed_fname, line_no))

                regions.append((line_no, tokens[0], start, end))
    except OSError as exc:
        fail("Error: could not read BED file %s: %s." %
             (bed_fname, exc))
    return regions

# This function assumes that the input VCF is sorted by coordinate
def restrict_vcf_to_regions_sweep(in_vcf_fname, out_vcf_fname, regions):
    regions_by_contig = {}
    for _, contig_name, start, end in regions:
        regions_by_contig.setdefault(contig_name, []).append((start, end))

    for contig_name, contig_regions in regions_by_contig.items():
        contig_regions.sort()
        merged_regions = []
        for start, end in contig_regions:
            if merged_regions and start < merged_regions[-1][1]:
                merged_regions[-1] = (merged_regions[-1][0], max(merged_regions[-1][1], end))
            else:
                merged_regions.append((start, end))
        regions_by_contig[contig_name] = merged_regions

    with pysam.VariantFile(in_vcf_fname) as in_vcf, \
         pysam.VariantFile(out_vcf_fname, 'wz', header=in_vcf.header) as out_vcf:
        current_contig = None
        current_regions = []
        current_region_idx = 0

        for record in in_vcf:
            if record.contig != current_contig:
                current_contig = record.contig
                current_regions = regions_by_contig.get(current_contig, [])
                current_region_idx = 0

            while current_region_idx < len(current_regions) and current_regions[current_region_idx][1] <= record.start:
                current_region_idx += 1

            if current_region_idx == len(current_regions):
                continue

            region_start, region_end = current_regions[current_region_idx]
            if region_start <= record.start and record.stop <= region_end:
                out_vcf.write(record)

def validate_bam_reference_compatibility(bam_fname, reference_fname, sampling_regions_fname=None):
    skip_hint = " You can bypass this check with --skip-bam-ref-validation. This is generally discouraged. Only use it if you are sure this will not cause problems downstream."
    bam_file = pysam.AlignmentFile(bam_fname, reference_filename=reference_fname)
    reference_fasta = pysam.FastaFile(reference_fname)
    try:
        bam_contig_lengths = dict(zip(bam_file.references, bam_file.lengths))
        bam_sq_records = {
            sq["SN"]: sq
            for sq in bam_file.header.to_dict().get("SQ", [])
            if "SN" in sq
        }
        reference_contig_lengths = {
            contig_name: reference_fasta.get_reference_length(contig_name)
            for contig_name in reference_fasta.references
        }
        missing_contigs = [
            contig_name for contig_name in bam_contig_lengths
            if contig_name not in reference_contig_lengths
        ]
        if missing_contigs:
            fail("Error: the reference FASTA is missing contig(s) present in the BAM header: %s" %
                 ", ".join(missing_contigs[:10]) + (" ..." if len(missing_contigs) > 10 else "") + "." + skip_hint)

        length_mismatches = []
        for contig_name, bam_len in bam_contig_lengths.items():
            reference_len = reference_contig_lengths[contig_name]
            if bam_len != reference_len:
                length_mismatches.append((contig_name, bam_len, reference_len))

        if length_mismatches:
            formatted = ", ".join(
                "%s (BAM=%d, FASTA=%d)" % (contig_name, bam_len, reference_len)
                for contig_name, bam_len, reference_len in length_mismatches[:10]
            )
            if len(length_mismatches) > 10:
                formatted += " ..."
            fail("Error: the BAM header and reference FASTA disagree on contig length(s): %s.%s" % (formatted, skip_hint))

        md5_mismatches = []
        md5_cache = {}
        for contig_name in bam_contig_lengths:
            bam_md5 = bam_sq_records.get(contig_name, {}).get("M5")
            if not bam_md5:
                continue

            if contig_name not in md5_cache:
                md5_cache[contig_name] = compute_reference_contig_md5(reference_fasta, contig_name)
            if md5_cache[contig_name] != bam_md5.lower():
                md5_mismatches.append((contig_name, bam_md5, md5_cache[contig_name]))

        if md5_mismatches:
            formatted = ", ".join(
                "%s (BAM=%s, FASTA=%s)" % (contig_name, bam_md5, reference_md5)
                for contig_name, bam_md5, reference_md5 in md5_mismatches[:10]
            )
            if len(md5_mismatches) > 10:
                formatted += " ..."
            fail("Error: the BAM header and reference FASTA disagree on contig MD5 checksum(s): %s.%s" % (formatted, skip_hint))

        if sampling_regions_fname:
            sampling_regions = parse_bed_regions(sampling_regions_fname)
            for line_no, contig_name, start, end in sampling_regions:
                if contig_name not in bam_contig_lengths:
                    fail("Error: sampling regions file %s line %d references contig %s, which is missing from the BAM header." %
                         (sampling_regions_fname, line_no, contig_name))
                if contig_name not in reference_contig_lengths:
                    fail("Error: sampling regions file %s line %d references contig %s, which is missing from the reference FASTA." %
                         (sampling_regions_fname, line_no, contig_name))
                contig_len = reference_contig_lengths[contig_name]
                if end > contig_len:
                    fail("Error: sampling regions file %s line %d extends past the end of contig %s (end=%d, contig_len=%d)." %
                         (sampling_regions_fname, line_no, contig_name, end, contig_len))
    finally:
        bam_file.close()
        reference_fasta.close()

if cmd_args.command in ('call', 'genotype') and not cmd_args.skip_bam_ref_validation:
    validate_bam_reference_compatibility(cmd_args.bam_file, cmd_args.reference, cmd_args.sampling_regions)

if cmd_args.command in ('call', 'genotype'):
    validate_workdir_for_execution(cmd_args.workdir, cmd_args.force_workdir)

mkdir(cmd_args.workdir)

if cmd_args.command in ('call', 'genotype'):
    with open(cmd_args.workdir + "/cmd.txt", "w") as full_cmd_out:
        print(" ".join(sys.argv[:]), file=full_cmd_out)

mkdir(cmd_args.workdir + "/intermediate_results")

def reads_categorizer(workdir):
    mkdir(workdir)
    with open(workdir + "/config.txt", "w") as config_file:
        config_file.write("threads %d\n" % cmd_args.threads)
        config_file.write("seed %d\n" % cmd_args.seed)
        config_file.write("min_sv_size %s\n" % cmd_args.min_sv_size)
        config_file.write("min_clip_len %s\n" % cmd_args.min_clip_len)
        config_file.write("max_seq_error %s\n" % cmd_args.max_seq_error)
        if cmd_args.sampling_regions:
            config_file.write("sampling_regions %s\n" % cmd_args.sampling_regions)
        config_file.write("per_contig_stats %d\n" % cmd_args.per_contig_stats)
        config_file.write("training_mode %d\n" % cmd_args.generate_training_data)
        config_file.write("version %s\n" % VERSION)
        config_file.write("min_diff_hsr %s\n" % cmd_args.min_diff_hsr)
        config_file.write("min_stable_mapq %d\n" % cmd_args.min_stable_mapq)
        config_file.write("high_confidence_mapq %d\n" % cmd_args.high_confidence_mapq)
        if cmd_args.haploid_contigs:
            config_file.write("haploid_contigs %s\n" % cmd_args.haploid_contigs)

    mkdir(workdir + "/workspace")
    mkdir_clean(workdir + "/workspace/sr")
    mkdir_clean(workdir + "/workspace/hsr")
    mkdir_clean(workdir + "/workspace/fwd-stable")
    mkdir_clean(workdir + "/workspace/rev-stable")
    mkdir_clean(workdir + "/workspace/long-pairs")
    mkdir_clean(workdir + "/workspace/outward-pairs")
    mkdir_clean(workdir + "/workspace/same-strand")
    mkdir_clean(workdir + "/workspace/mateseqs")

    with open("%s/contig_map" % workdir, "w") as contig_map:
        bam_file = pysam.AlignmentFile(cmd_args.bam_file, reference_filename=cmd_args.reference)
        for k in bam_file.references:
            contig_map.write("%s\n" % (k))

    read_categorizer_cmd = SURVEYOR_PATH + "/bin/reads_categorizer %s %s %s" % (cmd_args.bam_file, workdir, cmd_args.reference)
    run_cmd(read_categorizer_cmd)

def deduplicate_vcf(vcf_fname, deduped_vcf_fname):
    if cmd_args.tr_bed:
        compare_cmd = SURVEYOR_PATH + "/bin/compare %s %s -R %s -T %s > %s/compare.txt" % (vcf_fname, vcf_fname, cmd_args.reference, cmd_args.tr_bed, cmd_args.workdir)
    else:
        compare_cmd = SURVEYOR_PATH + "/bin/compare %s %s -R %s > %s/compare.txt" % (vcf_fname, vcf_fname, cmd_args.reference, cmd_args.workdir)
    run_cmd(compare_cmd)

    with open(cmd_args.workdir + "/compare.txt") as compare_file:
        epr_vals, imprecise_vals = {}, {}
        with pysam.VariantFile(vcf_fname) as vcf:
            for record in vcf:
                epr_vals[record.id] = record.samples[0].get('EPR', 0)
                imprecise_vals[record.id] = 'IMPRECISE' in record.info

        removed_ids = set()
        for line in compare_file:
            id1, id2 = line.strip().split()[:2]
            if id1 >= id2 or id2 == "NONE": # each pair will be output twice, e.g. A B and B A. Let us process it once
                continue

            if imprecise_vals[id1] and not imprecise_vals[id2]:
                removed_ids.add(id1)
            elif imprecise_vals[id2] and not imprecise_vals[id1]:
                removed_ids.add(id2)
            elif epr_vals[id1] > epr_vals[id2]:
                removed_ids.add(id2)
            elif epr_vals[id2] > epr_vals[id1]:
                removed_ids.add(id1)

        vcf_deduped = deduped_vcf_fname
        with pysam.VariantFile(vcf_fname) as vcf, pysam.VariantFile(vcf_deduped, 'w', header=vcf.header) as out_vcf:
            for record in vcf:
                if record.id not in removed_ids:
                    out_vcf.write(record)

def separate_ins_to_dup(in_vcf_fname, ins_to_dup_vcf_fname, remaining_vcf_fname):
    with pysam.VariantFile(in_vcf_fname) as in_vcf, \
         pysam.VariantFile(ins_to_dup_vcf_fname, 'w', header=in_vcf.header) as ins_to_dup_vcf, \
         pysam.VariantFile(remaining_vcf_fname, 'w', header=in_vcf.header) as remaining_vcf:
        for record in in_vcf:
            if "INS_TO_DUP" in record.info:
                record.id = record.id[:-4] # remove the _DUP suffix
                ins_to_dup_vcf.write(record)
            else:
                remaining_vcf.write(record)


def call_candidate_variants(bam_fname, workdir, reference_fname, sample_name):
    reads_categorizer(workdir)
    max_is = get_max_is_from_stats(workdir + "/stats.txt")
    with open(workdir + "/config.txt", "a") as config_file:
        config_file.write("max_trans_size %d\n" % cmd_args.max_trans_size)

    mkdir_clean(workdir + "/workspace/consensuses")
    clip_consensus_builder_cmd = SURVEYOR_PATH + "/bin/clip_consensus_builder %s %s %s" % (workdir, reference_fname, sample_name)
    run_cmd(clip_consensus_builder_cmd)

    find_svs_from_sr_consensuses_cmd = SURVEYOR_PATH + "/bin/find_svs_from_sr_consensuses %s %s %s %s" % (bam_fname, workdir, reference_fname, sample_name)
    run_cmd(find_svs_from_sr_consensuses_cmd)

    normalise_cmd = SURVEYOR_PATH + "/bin/normalise %s/intermediate_results/sr.vcf.gz %s/intermediate_results/sr.norm.vcf.gz %s %d %d %d" % (workdir, workdir, reference_fname, cmd_args.threads, cmd_args.min_sv_size, max_is)
    run_cmd(normalise_cmd)

    merge_identical_calls_cmd = SURVEYOR_PATH + "/bin/merge_identical_calls %s/intermediate_results/sr.norm.vcf.gz %s/intermediate_results/sr.norm.dedup.vcf.gz %s" % (workdir, workdir, reference_fname)
    run_cmd(merge_identical_calls_cmd)

    dp_clusterer = SURVEYOR_PATH + "/bin/dp_clusterer %s %s %s %s" % (bam_fname, workdir, reference_fname, sample_name)
    run_cmd(dp_clusterer)

    ins_assembler_cmd = SURVEYOR_PATH + "/bin/insertions_assembler %s %s %s" % (workdir, reference_fname, sample_name)
    run_cmd(ins_assembler_cmd)

    concat_cmd = SURVEYOR_PATH + "/bin/concat_vcf %s/intermediate_results/sr_dp.vcf.gz %s/intermediate_results/assembled_ins.vcf.gz %s/intermediate_results/out.vcf.gz" % (workdir, workdir, workdir)
    run_cmd(concat_cmd)

    normalise_cmd = SURVEYOR_PATH + "/bin/normalise %s/intermediate_results/out.vcf.gz %s/intermediate_results/out.norm.vcf.gz %s %d %d %d" % (workdir, workdir, reference_fname, cmd_args.threads, cmd_args.min_sv_size, max_is)
    run_cmd(normalise_cmd)
    
    expand_aux_haplotypes_cmd = SURVEYOR_PATH + "/bin/expand_aux_haplotypes %s/intermediate_results/out.norm.vcf.gz %s/intermediate_results/out.norm.hap.vcf.gz %s" % (workdir, workdir, reference_fname)
    run_cmd(expand_aux_haplotypes_cmd)

    merge_identical_calls_cmd = SURVEYOR_PATH + "/bin/merge_identical_calls %s/intermediate_results/out.norm.hap.vcf.gz %s/intermediate_results/calls-raw.vcf.gz %s" % (workdir, workdir, reference_fname)
    run_cmd(merge_identical_calls_cmd)


def genotype_variants(bam_fname, workdir, reference_fname, sample_name, ml_model, n_iters, generate_training_data):
    insertions_to_duplications_cmd = SURVEYOR_PATH + "/bin/insertions_to_duplications %s/intermediate_results/calls-raw.vcf.gz %s/intermediate_results/calls-for-genotyping.vcf.gz %s %s" % (workdir, workdir, reference_fname, workdir)
    run_cmd(insertions_to_duplications_cmd)

    genotype_cmd = SURVEYOR_PATH + "/bin/genotype %s/intermediate_results/calls-for-genotyping.vcf.gz %s/intermediate_results/calls-with-fmt.vcf.gz %s %s %s %s" % (workdir, workdir, bam_fname, reference_fname, workdir, sample_name)
    run_cmd(genotype_cmd)

    if generate_training_data:
        separate_ins_to_dup(workdir + "/intermediate_results/calls-with-fmt.vcf.gz", workdir + "/training-data.INS_TO_DUP.vcf.gz", workdir + "/training-data.vcf.gz")

    if not ml_model:
        print("No model provided. Skipping the classification step.")
        return

    Classifier.run_classifier(workdir + "/intermediate_results/calls-with-fmt.vcf.gz", workdir + "/intermediate_results/calls-with-gt.vcf.gz", workdir + "/stats.txt", ml_model, threads=cmd_args.threads)

    reconcile_vcf_gt_cmd = SURVEYOR_PATH + "/bin/reconcile_vcf_gt %s %s %s %s" % (workdir + "/intermediate_results/calls-raw.vcf.gz", workdir + "/intermediate_results/calls-with-gt.vcf.gz", workdir + "/intermediate_results/calls-with-gt.reconciled.vcf.gz", sample_name)
    run_cmd(reconcile_vcf_gt_cmd)

    write_aux_snps_cmd = SURVEYOR_PATH + "/bin/write_aux_snps %s/intermediate_results/calls-with-gt.reconciled.vcf.gz %s/calls-genotyped.smvars %s/calls-genotyped.stvars %s %s" % (workdir, workdir, workdir, reference_fname, workdir)
    run_cmd(write_aux_snps_cmd)

    if cmd_args.two_pass:
        cp_cmd = "cp %s/intermediate_results/calls-with-gt.vcf.gz %s/intermediate_results/calls-with-gt.iter0.vcf.gz" % (workdir, workdir)
        run_cmd(cp_cmd)

        for i in range(n_iters, n_iters+1):
            prev_iter_gt_file = workdir + "/intermediate_results/calls-with-gt.iter%d.vcf.gz" % (i-1)
            next_iter_fmt_file = workdir + "/intermediate_results/calls-with-fmt.iter%d.vcf.gz" % i
            genotype_cmd = SURVEYOR_PATH + "/bin/genotype %s %s %s %s %s %s --reassign-evidence" % (prev_iter_gt_file, next_iter_fmt_file, bam_fname, reference_fname, workdir, sample_name)
            run_cmd(genotype_cmd)

            next_iter_gt_file = workdir + "/intermediate_results/calls-with-gt.iter%d.vcf.gz" % i
            Classifier.run_classifier(next_iter_fmt_file, next_iter_gt_file, workdir + "/stats.txt", ml_model, threads=cmd_args.threads)

        final_iter_gt_file = next_iter_gt_file
        final_iter_fmt_file = next_iter_fmt_file
        if generate_training_data and os.path.exists(final_iter_fmt_file):
            separate_ins_to_dup(final_iter_fmt_file, workdir + "/training-data.reassigned.INS_TO_DUP.vcf.gz", workdir + "/training-data.reassigned.vcf.gz")

        reconciled_file = final_iter_gt_file.replace(".vcf.gz", ".reconciled.vcf.gz")
        reconcile_vcf_gt_cmd = SURVEYOR_PATH + "/bin/reconcile_vcf_gt %s %s %s %s" % (workdir + "/intermediate_results/calls-raw.vcf.gz", final_iter_gt_file, reconciled_file, sample_name)
        run_cmd(reconcile_vcf_gt_cmd)

        write_aux_snps_cmd = SURVEYOR_PATH + "/bin/write_aux_snps %s %s/calls-genotyped.reassigned.smvars %s/calls-genotyped.reassigned.stvars %s %s" % (reconciled_file, workdir, workdir, reference_fname, workdir)
        run_cmd(write_aux_snps_cmd)



if cmd_args.samplename:
    sample_name = cmd_args.samplename
else:
    sample_name = os.path.basename(cmd_args.bam_file).split(".")[0]


if cmd_args.command == 'call':

    if not cmd_args.ml_model:
        print("No model provided. Genotyped variants will not be output.", flush=True)

        if not cmd_args.generate_training_data:
            print("At least one of --ml-model or --generate-training-data must be provided.", flush=True)
            exit(1)
        if cmd_args.two_pass:
            print("Error: --two-pass requires --ml-model to be provided.", flush=True)
            exit(1)

    call_candidate_variants(cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name)

    n_iters = 1
    genotype_variants(cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name, cmd_args.ml_model, n_iters, cmd_args.generate_training_data)

    if not cmd_args.ml_model:
        exit(0)

    deduplicate_vcf(cmd_args.workdir + "/calls-genotyped.smvars.vcf.gz", cmd_args.workdir + "/calls-genotyped.smvars.deduped.vcf.gz")
    deduplicate_vcf(cmd_args.workdir + "/calls-genotyped.stvars.vcf.gz", cmd_args.workdir + "/calls-genotyped.stvars.deduped.vcf.gz")

elif cmd_args.command == 'genotype':

    if not cmd_args.ml_model:
        print("No model provided. Genotyped variants will not be output.", flush=True)

        if not cmd_args.generate_training_data:
            print("At least one of --ml-model or --generate-training-data must be provided.", flush=True)
            exit(1)
        if cmd_args.two_pass:
            print("Error: --two-pass requires --ml-model to be provided.", flush=True)
            exit(1)

    check_duplicate_ids_cmd = SURVEYOR_PATH + "/bin/check_duplicate_ids %s" % cmd_args.in_vcf_file
    run_cmd(check_duplicate_ids_cmd, "Error: Duplicate IDs found in the input VCF file. Please remove duplicates before running the genotype command.")

    reads_categorizer(cmd_args.workdir)

    cp_cmd = "cp %s %s/intermediate_results/calls-raw.vcf.gz" % (cmd_args.in_vcf_file, cmd_args.workdir)
    run_cmd(cp_cmd)

    n_iters = 1
    genotype_variants(cmd_args.bam_file, cmd_args.workdir, cmd_args.reference, sample_name, cmd_args.ml_model, n_iters, cmd_args.generate_training_data)

    if not cmd_args.ml_model:
        exit(0)

    deduplicate_vcf(cmd_args.workdir + "/calls-genotyped.smvars.vcf.gz", cmd_args.workdir + "/calls-genotyped.smvars.deduped.vcf.gz")
    deduplicate_vcf(cmd_args.workdir + "/calls-genotyped.stvars.vcf.gz", cmd_args.workdir + "/calls-genotyped.stvars.deduped.vcf.gz")

elif cmd_args.command == 'generate-training-data':

    original_training_data_vcf = cmd_args.workdir + "/training-data.vcf.gz"
    original_ins_to_dup_training_data_vcf = cmd_args.workdir + "/training-data.INS_TO_DUP.vcf.gz"
    reassigned_training_data_vcf = cmd_args.workdir + "/training-data.reassigned.vcf.gz"
    reassigned_ins_to_dup_training_data_vcf = cmd_args.workdir + "/training-data.reassigned.INS_TO_DUP.vcf.gz"

    if not os.path.exists(original_training_data_vcf):
        print("Error: training-data.vcf.gz not found in the workdir %s. Please make sure to run the call command with --generate-training-data first." % cmd_args.workdir, flush=True)
        exit(1)

    source_training_data_vcf = original_training_data_vcf
    source_ins_to_dup_training_data_vcf = original_ins_to_dup_training_data_vcf
    if cmd_args.use_reassigned_training_data:
        if not os.path.exists(reassigned_training_data_vcf) or not os.path.exists(reassigned_ins_to_dup_training_data_vcf):
            print("Error: reassigned training data not found in the workdir %s. Please make sure to run call or genotype with --generate-training-data, --ml-model, and --two-pass first." % cmd_args.workdir, flush=True)
            exit(1)
        source_training_data_vcf = reassigned_training_data_vcf
        source_ins_to_dup_training_data_vcf = reassigned_ins_to_dup_training_data_vcf

    restrict_to_bed_regions = None
    if cmd_args.restrict_to_bed:
        restrict_to_bed_regions = parse_bed_regions(cmd_args.restrict_to_bed)

    mkdir(cmd_args.outdir)
    sample_training_data_vcf = os.path.join(cmd_args.outdir, cmd_args.samplename + ".vcf.gz")
    sample_ins_to_dup_training_data_vcf = os.path.join(cmd_args.outdir, cmd_args.samplename + ".INS_TO_DUP.vcf.gz")
    if restrict_to_bed_regions:
        restrict_vcf_to_regions_sweep(source_training_data_vcf, sample_training_data_vcf, restrict_to_bed_regions)
        restrict_vcf_to_regions_sweep(source_ins_to_dup_training_data_vcf, sample_ins_to_dup_training_data_vcf, restrict_to_bed_regions)
    else:
        shutil.copyfile(source_training_data_vcf, sample_training_data_vcf)
        shutil.copyfile(source_ins_to_dup_training_data_vcf, sample_ins_to_dup_training_data_vcf)
    shutil.copyfile(cmd_args.workdir + "/stats.txt", os.path.join(cmd_args.outdir, cmd_args.samplename + ".stats"))

    # if training-data.reassigned.vcf.gz and training-data.reassigned.INS_TO_DUP.vcf.gz are present,
    # mark unreliable genotypes in the benchmark variants where FMT/AR1 or FMT/AR2 is different between
    # the non-reassigned and reassigned training data
    unreliable_cids = set()
    if not cmd_args.use_reassigned_training_data and \
        os.path.exists(reassigned_training_data_vcf) and \
        os.path.exists(reassigned_ins_to_dup_training_data_vcf):
        with pysam.VariantFile(original_training_data_vcf) as original_vcf, \
             pysam.VariantFile(reassigned_training_data_vcf) as reassigned_vcf:
            for orig_record, reassigned_record in zip(original_vcf, reassigned_vcf):
                if orig_record.id != reassigned_record.id:
                    print("Error: Variant IDs do not match between training-data.vcf.gz and training-data.reassigned.vcf.gz", flush=True)
                    print("Found %s and %s" % (orig_record.id, reassigned_record.id), flush=True)
                    exit(1)
                if 'AR1' in orig_record.samples[0]:
                    o_a1 = orig_record.samples[0]['AR1']
                    r_a1 = reassigned_record.samples[0]['AR1']
                    if o_a1 > 0 and r_a1/o_a1 <= 0.9:
                        unreliable_cids.add(orig_record.id)
                if 'AR2' in orig_record.samples[0] and orig_record.samples[0]['AR2'] != reassigned_record.samples[0]['AR2']:
                    o_a2 = orig_record.samples[0]['AR2']
                    r_a2 = reassigned_record.samples[0]['AR2']
                    if o_a2 > 0 and r_a2/o_a2 <= 0.9:
                        unreliable_cids.add(orig_record.id)

    def has_alt_allele(gt):
        """Return True if GT tuple has at least one ALT allele (allele index > 0)."""
        if gt is None:
            return False
        for a in gt:
            if a is not None and a > 0:
                return True
        return False

    updated_benchmark_vcf_path = os.path.join(cmd_args.outdir, cmd_args.samplename + "-benchmark.vcf.gz")
    if cmd_args.unreliable_gts:
        unreliable_bids = set()
        with open(cmd_args.unreliable_gts) as unreliable_file:
            for line in unreliable_file:
                unreliable_bids.add(line.strip())

        with pysam.VariantFile(cmd_args.benchmark_vcf) as benchmark_vcf, \
             pysam.VariantFile(updated_benchmark_vcf_path, 'wz', header=benchmark_vcf.header) as updated_benchmark_vcf:
            for record in benchmark_vcf:
                # if the variant is in the unreliable list and has a non-ref genotype, set it to ./1
                if record.id in unreliable_bids and has_alt_allele(record.samples[0]['GT']):
                    record.samples[0]['GT'] = (None, 1)
                updated_benchmark_vcf.write(record)
    else:
        shutil.copyfile(cmd_args.benchmark_vcf, updated_benchmark_vcf_path)

    read_len = get_stat_from_stats(cmd_args.workdir + "/stats.txt", "read_len")
    min_clip_len = int(get_config_value(cmd_args.workdir + "/config.txt", "min_clip_len"))
    hom_alt_split_radius = max(0, read_len - 2 * min_clip_len)

    compare_cmd = SURVEYOR_PATH + "/bin/compare %s %s -T %s -R %s --report -c %s -e -t %d --keep-all-called --hom-alt-split-radius %d" % (
        updated_benchmark_vcf_path,
        sample_training_data_vcf,
        cmd_args.simple_repeat_bed,
        cmd_args.reference,
        os.path.join(cmd_args.outdir, cmd_args.samplename + ".gts.tmp"),
        cmd_args.threads,
        hom_alt_split_radius
    )
    run_cmd(compare_cmd)

    os.remove(updated_benchmark_vcf_path)

    final_gts_path = os.path.join(cmd_args.outdir, cmd_args.samplename + ".gts")
    with open(os.path.join(cmd_args.outdir, cmd_args.samplename + ".gts.tmp")) as tmp_gts_file, \
         open(final_gts_path, 'w') as final_gts_file:
        for line in tmp_gts_file:
            fields = line.strip().split()
            id, gt, perfect = fields[0], fields[1], fields[2]
            if id in unreliable_cids and "1" in gt:
                gt = "./."
                perfect = "0"
            final_gts_file.write("%s %s %s\n" % (id, gt, perfect))

    os.remove(os.path.join(cmd_args.outdir, cmd_args.samplename + ".gts.tmp"))
