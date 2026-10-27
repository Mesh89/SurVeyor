# SurVeyor
An SV caller from paired-end NGS data.

## Installation

**Please download the source code from the latest release, since an appropriate trained model is provided.**
**We discourage cloning the repository, as it may not be compatible with the latest trained model.**

In order to compile the code, the following are required:
- A C and a C++ compiler are required. If the GCC suite is used, version 4.9.3 or above is required.
- CMake (2.8 or above)

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
```
python surveyor.py call --threads N_THREADS BAM_FILE WORKDIR REFERENCE_FASTA --ml-model TRAINED_MODEL
```
TRAINED_MODEL should be the location of the uncompressed trained-model.zip

For other parameters, please see the help with
```
python surveyor.py call -h
```

## Genotyping a set of structural variants

The basic command to genotype a set of structural variants in VCF format (IN_VCF) is
```
python surveyor.py genotype --threads N_THREADS IN_VCF BAM_FILE WORKDIR REFERENCE_FASTA --ml-model TRAINED_MODEL
```
The genotyped variants will be reported in genotyped.vcf.gz (version with all of the original calls) and genotyped.deduped.vcf.gz (version with duplicate calls removed).

For other parameters, please see the help with
```
python surveyor.py genotype -h
```

## Tandem duplications with unresolved copy number (since 0.12.3)
Some tandem duplications reported by SurVeyor have an unresolved copy number. In these cases, SurVeyor can identify the reference interval involved in the duplication, but it cannot determine from short-read evidence how many additional copies are present. For every DUP record and every sample, SurVeyor reports this in the FORMAT field CNU, both when calling and when genotyping.

For example:

```
FORMAT  GT:...:CNU
SAMPLE  0/1:...:1
```

When CNU is 1, END and SVLEN describe the duplicated reference interval. They should not be interpreted as the total inserted allele length. Rather, the reference interval is expanded in the alternative allele, which may be more than twice (or less than half) as long, for example because it contains several copies of the reported interval.

When CNU is 0, the evidence is consistent with an inserted allele between half and twice the length of the reported interval, i.e., roughly one additional tandem copy.

Since CNU is predicted separately for each sample, it is also available when genotyping a catalogue of SVs discovered in other samples with the genotype command above. CNU says nothing about whether the duplication is present in the sample; that is given by GT.

CNU requires the trained model released with version 0.12.3 or later. It replaces the INFO flag DUP_CN_UNRESOLVED of version 0.12.2.

## Citation

Currently, there is no manuscript for SurVeyor.

If you use deletions or tandem duplications, please cite: 
Rajaby, R., Sung, WK. SurVIndel2: improving copy number variant calling from next-generation sequencing using hidden split reads. Nat Commun 15, 10473 (2024). https://doi.org/10.1038/s41467-024-53087-7

If you use insertions, please cite: Rajaby, R., Liu, DX., Au, C.H. et al. 
INSurVeyor: improving insertion calling from short read sequencing data. Nat Commun 14, 3243 (2023). https://doi.org/10.1038/s41467-023-38870-2
