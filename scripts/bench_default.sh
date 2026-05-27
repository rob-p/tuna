#!/usr/bin/env bash
# =============================================================================
# bench_compare.sh — k-mer counter benchmark: tuna · KMC · FastK · KFC
# Scenario : count + text dump, k=31, single-genome inputs
# Datasets : 1000 E. coli genomes, 10 human genomes (from a file-of-files)
# Threads  : 1, 4, 8, 16   —   RAM budget: 256 GB
#
# Each genome is counted individually by all 4 tools.
# To avoid systematic OS page-cache bias (first tool always sees cold data),
# the tool order rotates:
#   E. coli : rotates every 250 files  (4 groups × 250 = 1000)
#   human   : rotates every 2 files    (5 groups × 2 = 10)
# With 4 tools the rotation cycles in 4 steps, so each tool starts first
# equally often across the dataset.
#
# Notes on tool differences:
#   tuna  : single command, writes TSV directly (no dump step)
#   KMC   : count → binary, then kmc_dump → TSV  (both timed together)
#           -fm is required for assembled FASTA (multi-sequence per file)
#   FastK : count → .ktab binary, then Tabex -A → TSV (both timed together)
#           .fna extension is not recognised; symlinks to .fa are created once
#           at setup (not timed)
#   KFC   : count → .kff binary, then kfc dump → TSV (both timed together)
#           no RAM budget flag exposed
#
# Output: $WORKDIR/results.log
#   one "=== ... ===" header + /usr/bin/time -v block per run
#   bring this file back locally and grep/parse for wall time and peak RSS
# =============================================================================
set -uo pipefail

# =============================================================================
# CONFIGURE — fill in before submitting to the cluster
# (values already set in the environment take precedence)
# =============================================================================
: "${TUNA:=""}"       # /path/to/tuna
: "${KMC:=""}"        # /path/to/kmc
: "${KMC_DUMP:=""}"   # /path/to/kmc_dump
: "${FASTK:=""}"      # /path/to/FastK
: "${TABEX:=""}"      # /path/to/Tabex
: "${KFC:=""}"        # /path/to/kfc

WORKDIR="/WORKS/vlevallois/expes_tuna"      
ECOLI_FOF="/WORKS/vlevallois/data/dataset_genome_ecoli/fof_1000.list"    # absolute path, one genome path per line, 1000 E. coli genomes
HUMAN_FOF="/WORKS/vlevallois/data/dataset_genome_human/fof_10.list"      # absolute path, one genome path per line, 10 human genomes

K=31
RAM_GB=256
THREADS_LIST=(1 4 8 16)

# Tool rotation order — do not change (drives the cycle logic below)
TOOLS=(tuna kmc fastk kfc)

# =============================================================================
# Sanity checks
# =============================================================================
err=0
for var in TUNA KMC KMC_DUMP FASTK TABEX KFC WORKDIR ECOLI_FOF HUMAN_FOF; do
    [[ -z "${!var}" ]] && { echo "[error] $var is not set"; err=1; }
done
for var in TUNA KMC KMC_DUMP FASTK TABEX KFC; do
    [[ -n "${!var}" && ! -x "${!var}" ]] && { echo "[error] ${!var} is not executable"; err=1; }
done
for var in ECOLI_FOF HUMAN_FOF; do
    [[ -n "${!var}" && ! -f "${!var}" ]] && { echo "[error] ${!var}: file not found"; err=1; }
done
[[ "$err" -eq 1 ]] && exit 1

ecoli_count=$(wc -l < "$ECOLI_FOF")
human_count=$(wc -l < "$HUMAN_FOF")
[[ "$ecoli_count" -lt 1000 ]] && echo "[warn] ECOLI_FOF has $ecoli_count lines (expected 1000)"
[[ "$human_count"  -lt 10   ]] && echo "[warn] HUMAN_FOF has $human_count lines (expected 10)"

# =============================================================================
# Directory setup
# =============================================================================
mkdir -p "$WORKDIR"/fastk_links
for tool in tuna kmc fastk kfc; do
    mkdir -p "$WORKDIR/$tool/runs"
done

RESULTS="$WORKDIR/results.log"
{
    echo "# bench_compare.sh — $(date)"
    echo "# TUNA=$TUNA"
    echo "# KMC=$KMC  KMC_DUMP=$KMC_DUMP"
    echo "# FASTK=$FASTK  TABEX=$TABEX"
    echo "# KFC=$KFC"
    echo "# WORKDIR=$WORKDIR"
    echo "# ECOLI_FOF=$ECOLI_FOF  ($ecoli_count genomes)"
    echo "# HUMAN_FOF=$HUMAN_FOF  ($human_count genomes)"
    echo "# K=$K  RAM=${RAM_GB}GB  THREADS=${THREADS_LIST[*]}"
    echo "# Tool rotation: E. coli every 250 files, human every 2 files"
} > "$RESULTS"

# =============================================================================
# FastK: create .fa / .fa.gz symlinks once for all input files
# FastK determines file type from extension; .fna is not recognised.
# =============================================================================
fastk_link_path() {
    # Returns the path of the FastK-compatible symlink for a given input file.
    local f="$1"
    local base
    base=$(basename "$f")
    case "$base" in
        *.fna.gz) echo "$WORKDIR/fastk_links/${base%.fna.gz}.fa.gz" ;;
        *.fna)    echo "$WORKDIR/fastk_links/${base%.fna}.fa"       ;;
        *)        echo "$WORKDIR/fastk_links/$base"                  ;;
    esac
}

echo "[setup] Creating FastK-compatible symlinks..."
while read -r f; do
    link=$(fastk_link_path "$f")
    ln -sf "$f" "$link"
done < <(sort -u "$ECOLI_FOF" "$HUMAN_FOF")

# =============================================================================
# Helper: write a labelled header to results.log
# =============================================================================
log_header() {
    { echo ""; echo "=== $* DATE=$(date +%Y-%m-%dT%H:%M:%S) ==="; } >> "$RESULTS"
}

# =============================================================================
# Per-tool run functions (single genome, single thread count)
# Each function creates an outdir, runs the full pipeline under /usr/bin/time,
# appends stdout+stderr to results.log, then removes the outdir.
# =============================================================================

_run_tuna() {
    local genome_file="$1" outdir="$2" threads="$3"
    # tuna writes TSV directly — one command, no separate dump step
    /usr/bin/time -v \
        "$TUNA" \
            -k "$K" \
            -t "$threads" \
            -ram "$RAM_GB" \
            -hp \
            -w "$outdir/work" \
            "$genome_file" \
            "$outdir/out.tsv" \
        >> "$RESULTS" 2>&1
}

_run_kmc() {
    local genome_file="$1" outdir="$2" threads="$3"
    mkdir -p "$outdir/tmp"
    # -fm : multi-FASTA (required for assembled genomes with multiple contigs)
    # -ci1: include k-mers with count >= 1
    # -hp : hide percentage progress
    /usr/bin/time -v bash -c "
        set -e
        \"$KMC\" \
            -k${K} -ci1 -fm -m${RAM_GB} -hp -t${threads} \
            \"$genome_file\" \
            \"$outdir/out\" \
            \"$outdir/tmp\" \
        && \"$KMC_DUMP\" \"$outdir/out\" \"$outdir/out.tsv\"
    " >> "$RESULTS" 2>&1
}

_run_fastk() {
    local genome_file="$1" outdir="$2" threads="$3"
    local link
    link=$(fastk_link_path "$genome_file")
    mkdir -p "$outdir/tmp"
    # -t1: report k-mers with count >= 1
    # -N : output prefix
    # -P : temp directory for sorting
    /usr/bin/time -v bash -c "
        set -e
        \"$FASTK\" \
            -k${K} -t1 -T${threads} -M${RAM_GB} \
            -N\"$outdir/out\" \
            -P\"$outdir/tmp\" \
            \"$link\" \
        && \"$TABEX\" -A \"$outdir/out\" > \"$outdir/out.tsv\"
    " >> "$RESULTS" 2>&1
}

_run_kfc() {
    local genome_file="$1" outdir="$2" threads="$3"
    # -t 1 (build): superkmer solidity threshold = 1 (keep singletons)
    # -t 1 (dump) : output k-mers with count >= 1
    # KFC has no RAM budget flag
    /usr/bin/time -v bash -c "
        set -e
        \"$KFC\" build \
            -k ${K} \
            -t 1 \
            -T ${threads} \
            -i \"$genome_file\" \
            -o \"$outdir/out.kff\" \
        && \"$KFC\" dump \
            -i \"$outdir/out.kff\" \
            -o \"$outdir/out.tsv\" \
            -t 1 \
            -T ${threads}
    " >> "$RESULTS" 2>&1
}

# =============================================================================
# Dispatcher: run one tool on one genome and clean up
# =============================================================================
run_one() {
    local tool="$1" dataset="$2" gen_idx="$3" genome_file="$4" threads="$5"
    local outdir="$WORKDIR/$tool/runs/${dataset}_g${gen_idx}_t${threads}"
    mkdir -p "$outdir"

    log_header "TOOL=$tool DS=$dataset GEN=$gen_idx FILE=$(basename "$genome_file") T=$threads"

    case "$tool" in
        tuna)  _run_tuna  "$genome_file" "$outdir" "$threads" ;;
        kmc)   _run_kmc   "$genome_file" "$outdir" "$threads" ;;
        fastk) _run_fastk "$genome_file" "$outdir" "$threads" ;;
        kfc)   _run_kfc   "$genome_file" "$outdir" "$threads" ;;
    esac || echo "[warn] $tool failed: $dataset gen $gen_idx t=$threads" | tee -a "$RESULTS"

    rm -rf "$outdir"
}

# =============================================================================
# Main benchmark loop
#
# Outer loop : thread counts (1, 4, 8, 16)
# Inner loop : genomes from FOF, one at a time
#   For each genome, all 4 tools run in the current rotated order.
#   Tool order rotates every CYCLE_SIZE files:
#     E. coli : CYCLE_SIZE=250  → groups 1-250, 251-500, 501-750, 751-1000
#     human   : CYCLE_SIZE=2    → pairs 1-2, 3-4, 5-6, 7-8, 9-10
#
# Rotation:
#   group 0 → [tuna, kmc, fastk, kfc]
#   group 1 → [kmc, fastk, kfc, tuna]
#   group 2 → [fastk, kfc, tuna, kmc]
#   group 3 → [kfc, tuna, kmc, fastk]
# =============================================================================
echo "[bench] Starting benchmark — $(date)"
n_tools=${#TOOLS[@]}

for threads in "${THREADS_LIST[@]}"; do
    echo "[bench] ── threads=$threads ───────────────────────────────"

    # ── E. coli: 1000 genomes, rotate tool order every 250 ──
    echo "[bench]   E. coli (1000 genomes, cycle=250)"
    gen_idx=0
    while IFS= read -r genome_file; do
        gen_idx=$(( gen_idx + 1 ))
        rotation=$(( (gen_idx - 1) / 250 % n_tools ))
        [[ $(( (gen_idx - 1) % 50 )) -eq 0 ]] && \
            echo "[bench]     ecoli gen $gen_idx/1000 rotation=$rotation t=$threads — $(date +%H:%M:%S)"
        for (( i=0; i<n_tools; i++ )); do
            tool_idx=$(( (i + rotation) % n_tools ))
            run_one "${TOOLS[$tool_idx]}" ecoli "$gen_idx" "$genome_file" "$threads"
        done
    done < "$ECOLI_FOF"

    # ── human: 10 genomes, rotate tool order every 2 ──
    echo "[bench]   human (10 genomes, cycle=2)"
    gen_idx=0
    while IFS= read -r genome_file; do
        gen_idx=$(( gen_idx + 1 ))
        rotation=$(( (gen_idx - 1) / 2 % n_tools ))
        echo "[bench]     human gen $gen_idx/10 rotation=$rotation t=$threads — $(date +%H:%M:%S)"
        for (( i=0; i<n_tools; i++ )); do
            tool_idx=$(( (i + rotation) % n_tools ))
            run_one "${TOOLS[$tool_idx]}" human "$gen_idx" "$genome_file" "$threads"
        done
    done < "$HUMAN_FOF"

done

echo "" >> "$RESULTS"
echo "# Done — $(date)" >> "$RESULTS"
echo "[bench] Done — results in $RESULTS"
