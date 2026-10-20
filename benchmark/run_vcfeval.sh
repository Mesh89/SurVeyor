#!/bin/bash
set -euo pipefail

# Usage: ./run_vcfeval.sh TRUTH_VCF QUERY_VCF REF_FA OUTDIR [REGIONS_BED] [--no-rtg-all] [--rtg RTG_PATH]
# Example:
#   ./run_vcfeval.sh GIAB_truth.vcf.gz mycalls.vcf.gz hg38.fa results giab_confident.bed

TRUTH_VCF=$1
QUERY_VCF=$2
REF=$3
REF_SDF=${REF%.fa}.sdf  # Assuming .fa extension; adjust if needed
OUTDIR=$4
OUT_PREFIX="out"
REGIONS_BED=""         # optional BED / BED.GZ of regions to evaluate
RUN_RTG_ALL=true
RTG=rtg               # defaults to rtg on PATH

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

PIDS=()
bcftools view --min-ac=1 "$TRUTH_VCF" | \
  bcftools sort -T "$VCFEVAL_TMPDIR/bcftools.XXXXXX" -Oz -o "$TRUTH_ORIG" &
PIDS+=("$!")

bcftools view --min-ac=1 "$QUERY_VCF" | \
  bcftools sort -T "$VCFEVAL_TMPDIR/bcftools.XXXXXX" -Oz -o "$QUERY_ORIG" &
PIDS+=("$!")

echo "[INFO] Normalizing truth VCF..."
bcftools view --min-ac=1 "$TRUTH_VCF" | \
  bcftools norm -f "$REF" -m -any --check-ref w --multi-overlaps 0 | \
  bcftools sort -T "$VCFEVAL_TMPDIR/bcftools.XXXXXX" -Oz -o "$TRUTH_NORM" &
PIDS+=("$!")

echo "[INFO] Normalizing query VCF..."
bcftools view --min-ac=1 "$QUERY_VCF" | \
  bcftools norm -f "$REF" -m -any --check-ref w --multi-overlaps 0 | \
  bcftools sort -T "$VCFEVAL_TMPDIR/bcftools.XXXXXX" -Oz -o "$QUERY_NORM" &
PIDS+=("$!")

if $RUN_RTG_ALL; then
  bcftools view "$QUERY_VCF" | \
    bcftools norm -f "$REF" -m -any --check-ref w --multi-overlaps 0 | \
    bcftools sort -T "$VCFEVAL_TMPDIR/bcftools.XXXXXX" -Oz -o "$QUERY_NORM_ALL" &
  PIDS+=("$!")
fi
wait_for_jobs "${PIDS[@]}"

PIDS=()
tabix -f -p vcf "$TRUTH_ORIG" &
PIDS+=("$!")
tabix -f -p vcf "$QUERY_ORIG" &
PIDS+=("$!")
tabix -f -p vcf "$TRUTH_NORM" &
PIDS+=("$!")
tabix -f -p vcf "$QUERY_NORM" &
PIDS+=("$!")
if $RUN_RTG_ALL; then
  tabix -f -p vcf "$QUERY_NORM_ALL" &
  PIDS+=("$!")
fi
wait_for_jobs "${PIDS[@]}"

PIDS=()
"$SCRIPT_DIRNAME/../bin/set_gt" "$TRUTH_NORM" "$TRUTH_HET" &
PIDS+=("$!")
"$SCRIPT_DIRNAME/../bin/set_gt" "$QUERY_NORM" "$QUERY_HET" &
PIDS+=("$!")
if $RUN_RTG_ALL; then
  "$SCRIPT_DIRNAME/../bin/set_gt" "$QUERY_NORM_ALL" "$QUERY_HET_ALL" --deduplicate 2 &
  PIDS+=("$!")
fi
wait_for_jobs "${PIDS[@]}"

PIDS=()
tabix -f -p vcf "$TRUTH_HET" &
PIDS+=("$!")
tabix -f -p vcf "$QUERY_HET" &
PIDS+=("$!")
if $RUN_RTG_ALL; then
  tabix -f -p vcf "$QUERY_HET_ALL" &
  PIDS+=("$!")
fi
wait_for_jobs "${PIDS[@]}"

# Optional BED handling
if [[ -n "${REGIONS_BED}" ]]; then
  if [[ ! -f "$REGIONS_BED" ]]; then
    echo "[ERROR] REGIONS_BED not found: $REGIONS_BED" >&2
    exit 1
  fi
  echo "[INFO] Restricting benchmarking to regions in: $REGIONS_BED"
fi

echo "[INFO] Running RTG vcfeval..."

PIDS=()
rm -rf "$OUTDIR/rtg/"
"$RTG" RTG_MEM=64g vcfeval -b $TRUTH_ORIG -c $QUERY_ORIG -o $OUTDIR/rtg/ -t $REF_SDF ${REGIONS_BED:+-e "$REGIONS_BED"} --ref-overlap &
PIDS+=("$!")

rm -rf "$OUTDIR/rtg-bi/"
"$RTG" RTG_MEM=64g vcfeval -b $TRUTH_NORM -c $QUERY_NORM -o $OUTDIR/rtg-bi/ -t $REF_SDF ${REGIONS_BED:+-e "$REGIONS_BED"} --ref-overlap &
PIDS+=("$!")

rm -rf "$OUTDIR/rtg-nogt/"
"$RTG" RTG_MEM=64g vcfeval -b $TRUTH_HET -c $QUERY_HET -o $OUTDIR/rtg-nogt/ -t $REF_SDF ${REGIONS_BED:+-e "$REGIONS_BED"} --ref-overlap &
PIDS+=("$!")

if $RUN_RTG_ALL; then
  rm -rf "$OUTDIR/rtg-all/"
  "$RTG" RTG_MEM=64g vcfeval -b $TRUTH_HET -c $QUERY_HET_ALL -o $OUTDIR/rtg-all/ -t $REF_SDF ${REGIONS_BED:+-e "$REGIONS_BED"} --ref-overlap &
  PIDS+=("$!")
fi
wait_for_jobs "${PIDS[@]}"

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

# Function to analyze a specific bin size
analyze_bin() {
    local DIR=$1
    local MIN=$2
    local MAX=$3

    # --- SENSITIVITY (TP-base & FN) ---
    
    # Deletions
    TP_DEL_SENS=$(bcftools view -H -i "strlen(REF)>strlen(ALT) && (strlen(REF)-1>=$MIN && strlen(REF)-1<=$MAX)" "$DIR/tp-baseline.vcf.gz" 2>/dev/null | wc -l)
    FN_DEL_SENS=$(bcftools view -H -i "strlen(REF)>strlen(ALT) && (strlen(REF)-1>=$MIN && strlen(REF)-1<=$MAX)" "$DIR/fn.vcf.gz" 2>/dev/null | wc -l)

    # Insertions
    TP_INS_SENS=$(bcftools view -H -i "strlen(ALT)>strlen(REF) && (strlen(ALT)-1>=$MIN && strlen(ALT)-1<=$MAX)" "$DIR/tp-baseline.vcf.gz" 2>/dev/null | wc -l)
    FN_INS_SENS=$(bcftools view -H -i "strlen(ALT)>strlen(REF) && (strlen(ALT)-1>=$MIN && strlen(ALT)-1<=$MAX)" "$DIR/fn.vcf.gz" 2>/dev/null | wc -l)

    TP_SENS=$(($TP_DEL_SENS + $TP_INS_SENS))
    FN_SENS=$(($FN_DEL_SENS + $FN_INS_SENS))

    SENS_DEL=$(ratio "$TP_DEL_SENS" "$((TP_DEL_SENS + FN_DEL_SENS))")
    SENS_INS=$(ratio "$TP_INS_SENS" "$((TP_INS_SENS + FN_INS_SENS))")
    SENS_TOT=$(ratio "$TP_SENS" "$((TP_SENS + FN_SENS))")

    # --- PRECISION (TP-call & FP) ---

    # Deletions
    TP_DEL_PREC=$(bcftools view -H -i "strlen(REF)>strlen(ALT) && (strlen(REF)-1>=$MIN && strlen(REF)-1<=$MAX)" "$DIR/tp.vcf.gz" 2>/dev/null | wc -l)
    FP_DEL_PREC=$(bcftools view -H -i "strlen(REF)>strlen(ALT) && (strlen(REF)-1>=$MIN && strlen(REF)-1<=$MAX)" "$DIR/fp.vcf.gz" 2>/dev/null | wc -l)

    # Insertions
    TP_INS_PREC=$(bcftools view -H -i "strlen(ALT)>strlen(REF) && (strlen(ALT)-1>=$MIN && strlen(ALT)-1<=$MAX)" "$DIR/tp.vcf.gz" 2>/dev/null | wc -l)
    FP_INS_PREC=$(bcftools view -H -i "strlen(ALT)>strlen(REF) && (strlen(ALT)-1>=$MIN && strlen(ALT)-1<=$MAX)" "$DIR/fp.vcf.gz" 2>/dev/null | wc -l)

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
           
    # Optional separator between bins
    # echo "-----------------------------------------------------------------------------------------"
}

run_snp_metrics() {
    local DIR=$1
    local TP_SENS FN_SENS SENS TP_PREC FP_PREC PREC F1
    
    # Filter: Lengths are equal
    local FILTER="strlen(REF)=strlen(ALT)"

    # SENSITIVITY
    TP_SENS=$(bcftools view -H -i "$FILTER" "$DIR/tp-baseline.vcf.gz" 2>/dev/null | wc -l)
    FN_SENS=$(bcftools view -H -i "$FILTER" "$DIR/fn.vcf.gz" 2>/dev/null | wc -l)
    SENS=$(ratio "$TP_SENS" "$((TP_SENS + FN_SENS))")

    # PRECISION
    TP_PREC=$(bcftools view -H -i "$FILTER" "$DIR/tp.vcf.gz" 2>/dev/null | wc -l)
    FP_PREC=$(bcftools view -H -i "$FILTER" "$DIR/fp.vcf.gz" 2>/dev/null | wc -l)
    PREC=$(ratio "$TP_PREC" "$((TP_PREC + FP_PREC))")

    # F1
    F1=$(calc_f1 "$SENS" "$PREC")

    # Print Row
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

    # Print Table Header
    print_table_header

    # Define bins here: "MIN MAX"
    local BINS=("1 5" "6 10" "11 20" "21 30" "31 40" "41 50" "1 50")
    local BIN_TMPDIR
    BIN_TMPDIR=$(mktemp -d "$VCFEVAL_TMPDIR/bin-metrics.XXXXXX")
    local -a BIN_PIDS=()
    local -a BIN_OUTPUTS=()
    local BIN
    local MIN
    local MAX
    local OUTFILE
    local STATUS

    for BIN in "${BINS[@]}"; do
        read -r MIN MAX <<< "$BIN"
        OUTFILE="${BIN_TMPDIR}/bin-${MIN}-${MAX}.txt"
        BIN_OUTPUTS+=("$OUTFILE")
        analyze_bin "$DIR" "$MIN" "$MAX" > "$OUTFILE" &
        BIN_PIDS+=("$!")
    done

    wait_for_jobs "${BIN_PIDS[@]}" || {
        STATUS=$?
        rm -rf "$BIN_TMPDIR"
        return "$STATUS"
    }
    for OUTFILE in "${BIN_OUTPUTS[@]}"; do
        cat "$OUTFILE"
    done

    rm -rf "$BIN_TMPDIR"
    
    echo "-----------------------------------------------------------------------------------------"
    # Run SNPs at the bottom of the table
    run_snp_metrics "$DIR"
    echo ""
}

# --- EXECUTION ---

PIDS=()
run_metrics "$OUTDIR/rtg-nogt" "NON-GT AWARE (Squash)" > $OUTDIR/rtg-nogt/results.txt &
PIDS+=("$!")
run_metrics "$OUTDIR/rtg-bi" "BI-ALLELIC (Ref Overlap)" > $OUTDIR/rtg-bi/results.txt &
PIDS+=("$!")
run_metrics "$OUTDIR/rtg" "GT AWARE (Standard)" > $OUTDIR/rtg/results.txt &
PIDS+=("$!")
if $RUN_RTG_ALL; then
  run_metrics "$OUTDIR/rtg-all" "ALL VARIANTS (Forced HET)" > $OUTDIR/rtg-all/results.txt &
  PIDS+=("$!")
fi
wait_for_jobs "${PIDS[@]}"

cat $OUTDIR/rtg-nogt/results.txt
cat $OUTDIR/rtg-bi/results.txt
cat $OUTDIR/rtg/results.txt
if $RUN_RTG_ALL; then
  cat $OUTDIR/rtg-all/results.txt
fi

echo "[INFO] All comparisons done."
