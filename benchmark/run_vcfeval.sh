#!/bin/bash
set -euo pipefail

# Usage: ./run_vcfeval.sh TRUTH_VCF QUERY_VCF REF_FA OUTDIR [REGIONS_BED] [--no-rtg-all] [--rtg RTG_PATH]
#          [--query-include EXPR]
# Example:
#   ./run_vcfeval.sh GIAB_truth.vcf.gz mycalls.vcf.gz hg38.fa results giab_confident.bed
# --query-include keeps only query records passing the bcftools expression (e.g. 'REF!=ALT').

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 TRUTH_VCF QUERY_VCF REF_FA OUTDIR [REGIONS_BED] [--no-rtg-all] [--rtg RTG_PATH] [--query-include EXPR]" >&2
  exit 1
fi

TRUTH_VCF=$1
QUERY_VCF=$2
REF=$3
REF_SDF=${REF%.fa}.sdf  # Assuming .fa extension; adjust if needed
OUTDIR=$4
OUT_PREFIX="out"
REGIONS_BED=""         # optional BED / BED.GZ of regions to evaluate
RUN_RTG_ALL=true
RTG=rtg               # defaults to rtg on PATH
QUERY_INCLUDE=""      # optional bcftools expression applied to the query first

shift 4
while [[ $# -gt 0 ]]; do
  ARG=$1
  case "$ARG" in
    --rtg)
      if [[ $# -lt 2 || -z "$2" || "$2" == --* ]]; then
        echo "[ERROR] --rtg requires an executable path or command name" >&2
        exit 1
      fi
      RTG=$2
      shift
      ;;
    --no-rtg-all)
      RUN_RTG_ALL=false
      ;;
    --query-include)
      if [[ $# -lt 2 || -z "$2" ]]; then
        echo "[ERROR] --query-include requires a bcftools expression" >&2
        exit 1
      fi
      QUERY_INCLUDE=$2
      shift
      ;;
    --*)
      echo "[ERROR] Unknown option: $ARG" >&2
      exit 1
      ;;
    *)
      if [[ -z "$REGIONS_BED" ]]; then
        REGIONS_BED=$ARG
      else
        echo "[ERROR] Unknown argument: $ARG" >&2
        exit 1
      fi
      ;;
  esac
  shift
done

mkdir -p "$OUTDIR"
VCFEVAL_TMPDIR="$(cd "$OUTDIR" && pwd)/tmp"
mkdir -p "$VCFEVAL_TMPDIR"
export TMPDIR="$VCFEVAL_TMPDIR"
echo "[INFO] Temporary files: $VCFEVAL_TMPDIR"

SCRIPT_DIRNAME="$(dirname "$(realpath "$0")")"
SET_GT="$SCRIPT_DIRNAME/../bin/set_gt"
NORM_ARGS=(-f "$REF" -m -any --check-ref w --multi-overlaps 0)

wait_for_jobs() {
  local PID
  local STATUS=0
  for PID in "$@"; do
    wait "$PID" || {
      STATUS=$?
      echo "[ERROR] Background job $PID failed with exit status $STATUS" >&2
    }
  done
  return "$STATUS"
}

#  1. Original Files (Preserves original GT, excludes 0/0)
TRUTH_ORIG="$OUTDIR/truth.orig.vcf.gz"
QUERY_ORIG="$OUTDIR/query.orig.vcf.gz"

#  2. Standard Normalized Files (Preserves original GT, excludes 0/0)
TRUTH_NORM="$OUTDIR/truth.norm.vcf.gz"
QUERY_NORM="$OUTDIR/query.norm.vcf.gz"

# 3. Forced Heterozygous Files (POSITIVE variants GT forced to ./1 for Presence/Absence check)
TRUTH_HET="$OUTDIR/truth.het.vcf.gz"
QUERY_HET="$OUTDIR/query.het.vcf.gz"

# 4. Forced Heterozygous Files (ALL variants)
QUERY_NORM_ALL="$OUTDIR/query.norm.all.vcf.gz"
QUERY_HET_ALL="$OUTDIR/query.het.all.vcf.gz"

# Optional BED handling
if [[ -n "${REGIONS_BED}" ]]; then
  if [[ ! -f "$REGIONS_BED" ]]; then
    echo "[ERROR] REGIONS_BED not found: $REGIONS_BED" >&2
    exit 1
  fi
  echo "[INFO] Restricting benchmarking to regions in: $REGIONS_BED"
fi

# --- PREPARATION FUNCTIONS ---
# Every input is processed in chunks of whole contigs, in parallel, and the sorted chunks are
# concatenated in header contig order. Contigs never interact in filtering, normalization or
# set_gt, and bcftools sort is stable, so the result equals one pass over the whole file.

shopt -s nullglob
CHUNK_DIR="$VCFEVAL_TMPDIR/chunks"
# Separate from CHUNK_DIR: htslib also looks for a regions file's index under its name
# with the extension replaced, so truth.0001 would pick up an index named truth.tbi.
INPUT_DIR="$VCFEVAL_TMPDIR/inputs"
rm -rf "$CHUNK_DIR" "$INPUT_DIR"
mkdir -p "$CHUNK_DIR" "$INPUT_DIR"

# vcfeval reports only records inside REGIONS_BED, so other contigs (thousands of decoys,
# each a vcfeval job) are dropped up front.
TARGET_CONTIGS="$CHUNK_DIR/contigs.txt"
if [[ -n "$REGIONS_BED" ]]; then
  zcat -f "$REGIONS_BED" | awk '$1 !~ /^(#|track|browser)/ && NF { print $1 }' | sort -u > "$TARGET_CONTIGS"
else
  : > "$TARGET_CONTIGS"
fi

# Print a path for VCF that region queries can read, indexing an unindexed input into
# INPUT_DIR and compressing it first if it is plain text.
indexed_input() {
  local VCF=$1 NAME=$2
  if [[ -f "$VCF.tbi" || -f "$VCF.csi" ]]; then
    echo "$VCF"
  elif bcftools index -t --threads 3 -o "$INPUT_DIR/$NAME.tbi" "$VCF" 2> /dev/null; then
    echo "$VCF##idx##$INPUT_DIR/$NAME.tbi"
  else
    # Chained: errexit does not apply inside the command substitution that calls this.
    bcftools view --no-version --threads 3 -Oz -o "$INPUT_DIR/$NAME.vcf.gz" "$VCF" &&
      bcftools index -t --threads 3 "$INPUT_DIR/$NAME.vcf.gz" &&
      echo "$INPUT_DIR/$NAME.vcf.gz"
  fi
}

# Group the contigs with records in VCF, in header order (the order bcftools sort uses) and
# restricted to TARGET_CONTIGS when a BED is given, into regions files PREFIX.NNNN of
# similar record counts.
write_chunks() {
  local VCF=$1 PREFIX=$2
  bcftools view -h "$VCF" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' > "$PREFIX.contigs"
  bcftools index -s "$VCF" > "$PREFIX.stats"
  LC_ALL=C awk -F'\t' -v prefix="$PREFIX" -v restrict="${REGIONS_BED:+1}" -v n=32 '
    FILENAME == ARGV[1] { keep[$1] = 1; next }
    FILENAME == ARGV[2] { if (!($1 in rank)) rank[$1] = ++nr; next }
    !restrict || ($1 in keep) {
      if (!($1 in rank)) rank[$1] = ++nr
      name[rank[$1]] = $1; count[$1] = $3; total += $3
    }
    END {
      k = 1; acc = 0
      for (i = 1; i <= nr; i++) {
        if (!(i in name)) continue
        file = sprintf("%s.%04d", prefix, k)
        print name[i] "\t1\t2147483647" > file
        acc += count[name[i]]
        if (acc * n >= total) { close(file); k++; acc = 0 }
      }
    }' "$TARGET_CONTIGS" "$PREFIX.contigs" "$PREFIX.stats"
}

# Filter one chunk with the given bcftools view options, normalize it if NORMALIZE is true,
# and sort it into OUT. Messages go to OUT.log, shown only on failure.
prepare_part() {
  local IN=$1 REGIONS=$2 OUT=$3 NORMALIZE=$4
  shift 4
  {
    if $NORMALIZE; then
      bcftools view --no-version "$@" -R "$REGIONS" -Ou "$IN" | \
        bcftools norm --no-version "${NORM_ARGS[@]}" -Ou | \
        bcftools sort -T "$VCFEVAL_TMPDIR/bcftools.XXXXXX" -Oz -o "$OUT"
    else
      bcftools view --no-version "$@" -R "$REGIONS" -Ou "$IN" | \
        bcftools sort -T "$VCFEVAL_TMPDIR/bcftools.XXXXXX" -Oz -o "$OUT"
    fi
  } 2> "$OUT.log" || {
    local STATUS=$?
    cat "$OUT.log" >&2
    return "$STATUS"
  }
}

# Prepare one normalized chunk, then its ./1 version keeping at most DEDUP copies per
# exact allele (0 keeps all). Remaining arguments are bcftools view options.
prepare_het_part() {
  local IN=$1 REGIONS=$2 NORM_OUT=$3 HET_OUT=$4 DEDUP=$5
  shift 5
  prepare_part "$IN" "$REGIONS" "$NORM_OUT" true "$@"
  "$SET_GT" "$NORM_OUT" "$HET_OUT" --deduplicate "$DEDUP"
}

# Concatenate the chunk parts of one output in contig order and index it. Without parts
# (no records on the evaluated contigs) the output keeps only the header of SOURCE.
concat_parts() {
  local OUT=$1 PART_DIR=$2 SOURCE=$3
  local -a PARTS=("$PART_DIR"/*.vcf.gz)
  if [[ ${#PARTS[@]} -eq 0 ]]; then
    bcftools view --no-version -h -Oz -o "$OUT" "$SOURCE"
  else
    bcftools concat --naive -o "$OUT" "${PARTS[@]}" 2> "$PART_DIR/concat.log" || {
      local STATUS=$?
      cat "$PART_DIR/concat.log" >&2
      return "$STATUS"
    }
  fi
  bcftools index -t --threads 3 -f "$OUT"
}

# --- HELPER FUNCTIONS ---

# Helper to safely compute ratio (Returns 6 decimal float or -)
ratio() {
    local num=$1
    local den=$2
    if (( den == 0 )); then
        echo "-"
    else
        awk -v n="$num" -v d="$den" 'BEGIN { printf "%.4f", n/d }'
    fi
}

# Helper to calculate F1 score (Returns 6 decimal float or -)
calc_f1() {
    local s=$1
    local p=$2
    # Check if inputs are numbers (not "-") and greater than 0
    awk -v s="$s" -v p="$p" 'BEGIN { 
        if (s=="-" || p=="-" || s==0 || p==0) print "-"; 
        else printf "%.4f", 2*s*p/(s+p) 
    }'
}

print_table_header() {
    echo ""
    printf "%-12s | %-4s | %-8s | %-8s | %-8s | %-7s | %-7s | %-7s | %-7s\n" \
           "Bin" "Type" "Recall" "Prec" "F1" "TP_bs" "FN" "TP_cl" "FP"
    echo "-----------------------------------------------------------------------------------------"
}

# --- ANALYSIS FUNCTIONS ---

# Size bins: "MIN MAX"
BINS=("1 5" "6 10" "11 20" "21 30" "31 40" "41 50" "1 50")

# Count records per size bin in one decompression pass. Prints "DEL INS" per bin, then the
# SNP/MNP count. Each test matches bcftools -i on strlen(REF)/strlen(ALT): a multi-allelic
# record passes a comparison when any ALT does, and the comparisons are combined after that.
count_bins() {
    local FILE=$1
    local BIN_LIST
    BIN_LIST=$(IFS=,; echo "${BINS[*]}")
    bgzip -dc "$FILE" | cut -f1,4,5 | LC_ALL=C awk -F'\t' -v bins="$BIN_LIST" '
        BEGIN {
            nb = split(bins, b, ",")
            for (i = 1; i <= nb; i++) { split(b[i], mm, " "); lo[i] = mm[1]; hi[i] = mm[2]; del[i] = 0; ins[i] = 0 }
            snp = 0
        }
        substr($1, 1, 1) == "#" { next }
        {
            rl = length($2); n = split($3, alt, ",")
            shorter = 0; longer = 0; equal = 0; amin = -1; amax = -1
            for (j = 1; j <= n; j++) {
                al = length(alt[j])
                if (al < rl) shorter = 1; else if (al > rl) longer = 1; else equal = 1
                if (amin < 0 || al < amin) amin = al
                if (al > amax) amax = al
            }
            snp += equal
            for (i = 1; i <= nb; i++) {
                if (shorter && rl - 1 >= lo[i] && rl - 1 <= hi[i]) del[i]++
                if (longer && amax - 1 >= lo[i] && amin - 1 <= hi[i]) ins[i]++
            }
        }
        END { for (i = 1; i <= nb; i++) print del[i], ins[i]; print snp }'
}

# Print the DEL, INS and ALL rows of one bin from its counts
analyze_bin() {
    local MIN=$1 MAX=$2
    local TP_DEL_SENS=$3 TP_INS_SENS=$4 FN_DEL_SENS=$5 FN_INS_SENS=$6
    local TP_DEL_PREC=$7 TP_INS_PREC=$8 FP_DEL_PREC=$9 FP_INS_PREC=${10}

    # --- SENSITIVITY (TP-base & FN) ---
    TP_SENS=$(($TP_DEL_SENS + $TP_INS_SENS))
    FN_SENS=$(($FN_DEL_SENS + $FN_INS_SENS))

    SENS_DEL=$(ratio "$TP_DEL_SENS" "$((TP_DEL_SENS + FN_DEL_SENS))")
    SENS_INS=$(ratio "$TP_INS_SENS" "$((TP_INS_SENS + FN_INS_SENS))")
    SENS_TOT=$(ratio "$TP_SENS" "$((TP_SENS + FN_SENS))")

    # --- PRECISION (TP-call & FP) ---
    TP_PREC=$(($TP_DEL_PREC + $TP_INS_PREC))
    FP_PREC=$(($FP_DEL_PREC + $FP_INS_PREC))

    PREC_DEL=$(ratio "$TP_DEL_PREC" "$((TP_DEL_PREC + FP_DEL_PREC))")
    PREC_INS=$(ratio "$TP_INS_PREC" "$((TP_INS_PREC + FP_INS_PREC))")
    PREC_TOT=$(ratio "$TP_PREC" "$((TP_PREC + FP_PREC))")

    # --- CALCULATE F1 ---
    F1_DEL=$(calc_f1 "$SENS_DEL" "$PREC_DEL")
    F1_INS=$(calc_f1 "$SENS_INS" "$PREC_INS")
    F1_TOT=$(calc_f1 "$SENS_TOT" "$PREC_TOT")

    # --- PRINT TABLE ROWS ---
    # Format: Bin | Type | Recall | Prec | F1 | TP_bs | FN | TP_cl | FP
    
    printf "%-12s | %-4s | %-8s | %-8s | %-8s | %-7s | %-7s | %-7s | %-7s\n" \
           "${MIN}-${MAX}bp" "DEL" "$SENS_DEL" "$PREC_DEL" "$F1_DEL" "$TP_DEL_SENS" "$FN_DEL_SENS" "$TP_DEL_PREC" "$FP_DEL_PREC"

    printf "%-12s | %-4s | %-8s | %-8s | %-8s | %-7s | %-7s | %-7s | %-7s\n" \
           "${MIN}-${MAX}bp" "INS" "$SENS_INS" "$PREC_INS" "$F1_INS" "$TP_INS_SENS" "$FN_INS_SENS" "$TP_INS_PREC" "$FP_INS_PREC"

    printf "%-12s | %-4s | %-8s | %-8s | %-8s | %-7s | %-7s | %-7s | %-7s\n" \
           "${MIN}-${MAX}bp" "ALL" "$SENS_TOT" "$PREC_TOT" "$F1_TOT" "$TP_SENS" "$FN_SENS" "$TP_PREC" "$FP_PREC"
}

# Print the SNP/MNP row (REF and ALT of equal length) from its counts
run_snp_metrics() {
    local TP_SENS=$1 FN_SENS=$2 TP_PREC=$3 FP_PREC=$4
    local SENS PREC F1

    SENS=$(ratio "$TP_SENS" "$((TP_SENS + FN_SENS))")
    PREC=$(ratio "$TP_PREC" "$((TP_PREC + FP_PREC))")
    F1=$(calc_f1 "$SENS" "$PREC")

    printf "%-12s | %-4s | %-8s | %-8s | %-8s | %-7s | %-7s | %-7s | %-7s\n" \
           "SNP/MNP" "SNP" "$SENS" "$PREC" "$F1" "$TP_SENS" "$FN_SENS" "$TP_PREC" "$FP_PREC"
}

run_metrics() {
    local DIR=$1
    local LABEL=$2

    echo "=========================================================="
    echo "Running analysis for: $LABEL"
    echo "Directory: $DIR"
    echo "=========================================================="

    if [ ! -d "$DIR" ]; then
        echo "[ERROR] Directory $DIR does not exist. Skipping."
        return
    fi

    # Count all four vcfeval outputs concurrently, once each
    local COUNT_TMPDIR
    COUNT_TMPDIR=$(mktemp -d "$VCFEVAL_TMPDIR/bin-metrics.XXXXXX")
    local -a COUNT_PIDS=()
    local NAME STATUS
    for NAME in tp-baseline fn tp fp; do
        count_bins "$DIR/$NAME.vcf.gz" > "$COUNT_TMPDIR/$NAME" &
        COUNT_PIDS+=("$!")
    done
    wait_for_jobs "${COUNT_PIDS[@]}" || {
        STATUS=$?
        rm -rf "$COUNT_TMPDIR"
        return "$STATUS"
    }
    local -a TPB FN TP FP
    mapfile -t TPB < "$COUNT_TMPDIR/tp-baseline"
    mapfile -t FN < "$COUNT_TMPDIR/fn"
    mapfile -t TP < "$COUNT_TMPDIR/tp"
    mapfile -t FP < "$COUNT_TMPDIR/fp"
    rm -rf "$COUNT_TMPDIR"

    # Print Table Header
    print_table_header

    # Each count line is "DEL INS"; word splitting passes both numbers
    local I MIN MAX
    for I in "${!BINS[@]}"; do
        read -r MIN MAX <<< "${BINS[$I]}"
        analyze_bin "$MIN" "$MAX" ${TPB[$I]} ${FN[$I]} ${TP[$I]} ${FP[$I]}
    done

    echo "-----------------------------------------------------------------------------------------"
    # Run SNPs at the bottom of the table
    local N=${#BINS[@]}
    run_snp_metrics "${TPB[$N]}" "${FN[$N]}" "${TP[$N]}" "${FP[$N]}"
    echo ""
}

# --- PREPARATION ---

echo "[INFO] Preparing truth and query VCFs..."
TRUTH_IN=$(indexed_input "$TRUTH_VCF" truth)
QUERY_IN=$(indexed_input "$QUERY_VCF" query)
write_chunks "$TRUTH_IN" "$CHUNK_DIR/truth"
write_chunks "$QUERY_IN" "$CHUNK_DIR/query"

PART_DIR="$VCFEVAL_TMPDIR/parts"
rm -rf "$PART_DIR"
OUTPUTS=(truth.orig truth.norm truth.het query.orig query.norm query.het)
if $RUN_RTG_ALL; then
  OUTPUTS+=(query.norm.all query.het.all)
fi
for NAME in "${OUTPUTS[@]}"; do
  mkdir -p "$PART_DIR/$NAME"
done

PIDS=()
for REGIONS in "$CHUNK_DIR"/truth.[0-9]*; do
  PART=$(basename "$REGIONS").vcf.gz
  prepare_part "$TRUTH_IN" "$REGIONS" "$PART_DIR/truth.orig/$PART" false --min-ac=1 &
  PIDS+=("$!")
  prepare_het_part "$TRUTH_IN" "$REGIONS" "$PART_DIR/truth.norm/$PART" "$PART_DIR/truth.het/$PART" 0 --min-ac=1 &
  PIDS+=("$!")
done
for REGIONS in "$CHUNK_DIR"/query.[0-9]*; do
  PART=$(basename "$REGIONS").vcf.gz
  prepare_part "$QUERY_IN" "$REGIONS" "$PART_DIR/query.orig/$PART" false --min-ac=1 ${QUERY_INCLUDE:+-i "$QUERY_INCLUDE"} &
  PIDS+=("$!")
  prepare_het_part "$QUERY_IN" "$REGIONS" "$PART_DIR/query.norm/$PART" "$PART_DIR/query.het/$PART" 0 --min-ac=1 ${QUERY_INCLUDE:+-i "$QUERY_INCLUDE"} &
  PIDS+=("$!")
  if $RUN_RTG_ALL; then
    # Retain at most two copies per exact normalized allele in the ALL query.
    prepare_het_part "$QUERY_IN" "$REGIONS" "$PART_DIR/query.norm.all/$PART" "$PART_DIR/query.het.all/$PART" 2 ${QUERY_INCLUDE:+-i "$QUERY_INCLUDE"} &
    PIDS+=("$!")
  fi
done
wait_for_jobs "${PIDS[@]}"

PIDS=()
for NAME in "${OUTPUTS[@]}"; do
  if [[ $NAME == truth.* ]]; then
    SOURCE=$TRUTH_IN
  else
    SOURCE=$QUERY_IN
  fi
  concat_parts "$OUTDIR/$NAME.vcf.gz" "$PART_DIR/$NAME" "$SOURCE" &
  PIDS+=("$!")
done
wait_for_jobs "${PIDS[@]}"
rm -rf "$CHUNK_DIR" "$INPUT_DIR" "$PART_DIR"

# --- EXECUTION ---

echo "[INFO] Running RTG vcfeval..."

# Evaluate one comparison, then tabulate its size bins into results.txt
run_comparison() {
  local BASELINE=$1 CALLS=$2 DIR=$3 LABEL=$4
  rm -rf "$DIR"
  "$RTG" RTG_MEM=64g vcfeval -b "$BASELINE" -c "$CALLS" -o "$DIR/" -t "$REF_SDF" ${REGIONS_BED:+-e "$REGIONS_BED"} --ref-overlap
  run_metrics "$DIR" "$LABEL" > "$DIR/results.txt"
}

PIDS=()
run_comparison "$TRUTH_ORIG" "$QUERY_ORIG" "$OUTDIR/rtg" "GT AWARE (Standard)" &
PIDS+=("$!")
run_comparison "$TRUTH_NORM" "$QUERY_NORM" "$OUTDIR/rtg-bi" "BI-ALLELIC (Ref Overlap)" &
PIDS+=("$!")
run_comparison "$TRUTH_HET" "$QUERY_HET" "$OUTDIR/rtg-nogt" "NON-GT AWARE (Squash)" &
PIDS+=("$!")
if $RUN_RTG_ALL; then
  run_comparison "$TRUTH_HET" "$QUERY_HET_ALL" "$OUTDIR/rtg-all" "ALL VARIANTS (Forced HET)" &
  PIDS+=("$!")
fi
wait_for_jobs "${PIDS[@]}"

cat "$OUTDIR/rtg-nogt/results.txt"
cat "$OUTDIR/rtg-bi/results.txt"
cat "$OUTDIR/rtg/results.txt"
if $RUN_RTG_ALL; then
  cat "$OUTDIR/rtg-all/results.txt"
fi

echo "[INFO] All comparisons done."
