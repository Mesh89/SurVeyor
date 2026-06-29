# SurVeyor
Structural variant discovery and genotyping pipeline from paired-end NGS data. SurVeyor can detect and genotype deletions, tandem duplications, insertions and inversions that are 50 base pairs or longer using WGS paired-end sequencing data. Furthermore, when multiple samples from a cohort are available, it is able to leverage the shared information to increase sample-level recall.

## Installation

**Please download the source code from the latest release, since an appropriate trained model is provided.**
**We discourage cloning the repository, as it may not be compatible with the latest trained model.**

In order to compile the code, the following are required:
- A C and a C++ compiler are required. If the GCC suite is used, version 4.9.3 or above is required.
- CMake (3.5 or above)

Downloaded the latest release from https://github.com/Mesh89/SurVeyor/releases , uncompress it and enter the directory.
The following commands should be sufficient

```
./build_htslib.sh
cmake -DCMAKE_BUILD_TYPE=Release . && make
```

If you are compiling on the same platform as where you will execute it, you can use -DNATIVE=ON to create faster executables
```
cmake -DCMAKE_BUILD_TYPE=Release -DNATIVE=ON . && make
```

Python 3 is necessary to run SurVeyor. Libraries NumPy (http://www.numpy.org/), PySam (https://github.com/pysam-developers/pysam), scikit-learn (https://scikit-learn.org/stable/) and xgboost (https://xgboost.readthedocs.io/en/stable/) are also necessary.

Please also download trained-model.zip from the release you are using, and uncompress it in a location of your convenience.

## Demo

You can verify that the software is running correctly using the provided demo.
Please run (assuming trained-model/ is the location of the uncompressed model)
```
mkdir workdir
python3 surveyor.py call demo/reads.bam workdir/ demo/ref.fa --ml-model trained-model/
```
The software should finish in a few seconds. Then, 
```
bcftools view -H workdir/calls-genotyped.vcf.gz --min-ac=1
```
should output 4 variants: two insertions, one deletion and one duplication. The expected variants are in demo/result.txt

## Data preparation

SurVeyor needs a BAM/CRAM file, a (possibly empty) working directory and reference genome in FASTA format. For the genotype command, it will also require an input VCF.

The BAM/CRAM file must be coordinate-sorted and indexed. Furthermore, the MD, MC and the MQ tag must be present for all primary alignments, when applicable.

Recent versions of BWA MEM will automatically add all the necessary tags. If you used an old version, a different aligner, or if for any reason your BAM/CRAM does not have the required tags, you can use the samtools command fixmate (for MC and MQ) and callmd (for MD) to add them.

## Denovo calling of structural variants

The basic command to call structural variants 
```bash
python surveyor.py call --threads N_THREADS BAM_FILE WORKDIR REFERENCE_FASTA --ml-model TRAINED_MODEL
```
TRAINED_MODEL should be the location of the uncompressed trained-model.zip

For other parameters, please see the help with
```bash
python surveyor.py call -h
```

## Genotyping a set of structural variants

The basic command to genotype a set of structural variants in VCF format (IN_VCF) is
```bash
python surveyor.py genotype --threads N_THREADS IN_VCF BAM_FILE WORKDIR REFERENCE_FASTA --ml-model TRAINED_MODEL
```
The genotyped variants will be reported in genotyped.vcf.gz (version with all of the original calls) and genotyped.deduped.vcf.gz (version with duplicate calls removed).

For other parameters, please see the help with
```
python surveyor.py genotype -h
```

## Calling structural variants on a cohort

SurVeyor can perform cohort-aware SV detection when a cohort of multiple samples is available. All samples must be aligned to the same reference genome.

Cohort-aware SV detection consists of three steps:
1) First, denovo calling must be performed for each sample individually
2) Next, calls from all samples must be clustered. This produces a non-redundant catalogue of SVs within the cohort
3) Finally, the catalogue must be genotyped on the sample(s) of interest

Let us assume we have $n$ samples, called $S_1$ to $S_n$.

### Denovo calling for all samples

First, we perform denovo calling for each sample, i.e., for each $i = 1,\dots,n$:
```bash
python surveyor.py call --threads N_THREADS BAM_FILE_S_i WORKDIR_S_i REFERENCE_FASTA --ml-model TRAINED_MODEL
```
where `BAM_FILE_S_i` is the BAM or CRAM file relative to sample $S_i$.

### Clustering individual calls into a unified catalogue

A text file `FILELIST` should be produced where each line contains a sample name and the corresponding SurVeyor VCF:

For our example, `FILELIST` should look like this
```text
SAMPLE_1 WORKDIR_S_1/calls-genotyped.vcf.gz
SAMPLE_2 WORKDIR_S_2/calls-genotyped.vcf.gz
...
SAMPLE_n WORKDIR_S_n/calls-genotyped.vcf.gz
```

Then, run the command
```bash
./bin/cluster FILELIST REFERENCE_FASTA -o OUT_PREFIX -t N_THREADS
```
A file `OUT_PREFIX.vcf.gz` will be produced.

For the sake of our example, let `OUT_PREFIX=catalogue`, meaning the output of the clustering step is `catalogue.vcf.gz`

### Re-genotyping the catalogue on all samples

We re-genotype `catalogue.vcf.gz` produced by the clustering step on each sample individually.

If you are genotyping a catalogue generated from many samples and you notice that the cohort-aware calls tend to contain many duplicate SVs, add the following flag to the genotype command
```
--tr-bed SIMPLE_REPEATS_BED
```
where `SIMPLE_REPEATS_BED` is a list of repetitive regions for the reference, in BED format. We recommend using the simpleRepeats table from UCSC Table browser. For convenience, the list for GRCh38 is provided in `simple_repeats_bed/`.

Concretely, for each $i = 1,\dots,n$:
```bash
python surveyor.py genotype --threads N_THREADS catalogue.vcf.gz BAM_FILE_S_i WORKDIR_S_i-regt REFERENCE_FASTA --ml-model TRAINED_MODEL
```
where `BAM_FILE_S_i` is the BAM or CRAM file relative to sample $S_i$.

For each $i$, final calls for sample $S_i$ will be in `WORKDIR_S_i-regt/genotyped.deduped.vcf.gz`

### Genotyping a precomputed cohort catalogue

If a clustered cohort catalogue is already available, users do not need to repeat denovo calling and clustering. The clustered VCF can be used directly as input to the `genotype` command.

For example, if `catalogue.vcf.gz` is a precomputed clustered cohort catalogue, run:
```bash
python surveyor.py genotype --threads N_THREADS catalogue.vcf.gz BAM_FILE WORKDIR REFERENCE_FASTA --ml-model TRAINED_MODEL
```

The genotyped calls will be written to:
```text
WORKDIR/genotyped.vcf.gz
WORKDIR/genotyped.deduped.vcf.gz
```

As mentioned previously, if the catalogue was generated from many samples and you notice an excessive number of duplicate calls, we recommend using the option 
```
--tr-bed SIMPLE_REPEATS_BED
```
where `SIMPLE_REPEATS_BED` is a list of repetitive regions for the reference in BED format, and using the `genotyped.deduped.vcf.gz` file.

### Tandem duplications with unresolved copy number (to appear in 0.13)
Some tandem duplications reported by SurVeyor have an unresolved copy number. In these cases, SurVeyor can identify the reference interval involved in the duplication, but it cannot determine from short-read evidence how additional copies are present. Such records are annotated with the INFO flag DUP_CN_UNRESOLVED.

For example:

```bash
SVTYPE=DUP;END=100500;SVLEN=500;DUP_CN_UNRESOLVED
````

When DUP_CN_UNRESOLVED is present, END and SVLEN describe the duplicated reference interval. They should not be interpreted as the total inserted allele length. Rather, it should be intepreted as the reference interval being expanded in the alternative allele. The event may correspond to a larger inserted sequence containing one or more copies of the reported reference interval.

When DUP_CN_UNRESOLVED is absent from a DUP record, SurVeyor found evidence consistent with a single additional tandem copy of the reported interval.

This annotation is intended to make tandem duplication calls easier to interpret, especially in repetitive regions where short reads may identify the duplicated sequence but not resolve the number of copied units.

## Citation

Currently, there is no manuscript for SurVeyor.

If you use deletions or tandem duplications, please cite: 
Rajaby, R., Sung, WK. SurVIndel2: improving copy number variant calling from next-generation sequencing using hidden split reads. Nat Commun 15, 10473 (2024). https://doi.org/10.1038/s41467-024-53087-7

If you use insertions, please cite: Rajaby, R., Liu, DX., Au, C.H. et al. 
INSurVeyor: improving insertion calling from short read sequencing data. Nat Commun 14, 3243 (2023). https://doi.org/10.1038/s41467-023-38870-2
