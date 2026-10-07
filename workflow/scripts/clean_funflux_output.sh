#!/usr/bin/env bash
# clean_funflux_output.sh
#
# Remove bulky intermediate files from a completed FunFlux, FunFluxL,
# FunFluxL+, or Funnotator output directory while keeping final result files.
# The script is dry-run by default and deletes only when called with --run.
#
# Usage:
#   clean_funflux_output.sh [OUTPUT_DIR]
#   clean_funflux_output.sh --target OUTPUT_DIR
#   clean_funflux_output.sh --run [OUTPUT_DIR]
#   clean_funflux_output.sh --run --target OUTPUT_DIR
#
# If OUTPUT_DIR is omitted, the current directory is used.

set -euo pipefail

DO_RUN=0
TARGET_DIR=""

usage() {
  sed -n '1,22p' "$0"
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run)
      DO_RUN=1
      shift
      ;;
    --dry-run)
      DO_RUN=0
      shift
      ;;
    -t|--target|--output-dir)
      [[ $# -ge 2 ]] || die "$1 requires a directory argument."
      TARGET_DIR="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      break
      ;;
    -*)
      die "Unknown option: $1"
      ;;
    *)
      [[ -z "$TARGET_DIR" ]] || die "Only one output directory can be specified."
      TARGET_DIR="$1"
      shift
      ;;
  esac
done

[[ $# -eq 0 ]] || die "Unexpected extra argument(s): $*"

TARGET_DIR="${TARGET_DIR:-.}"
[[ -d "$TARGET_DIR" ]] || die "Target is not a directory: $TARGET_DIR"

TARGET_DIR="$(cd "$TARGET_DIR" && pwd -P)"
[[ "$TARGET_DIR" != "/" ]] || die "Refusing to operate on /"
[[ "$TARGET_DIR" != "$HOME" ]] || die "Refusing to operate on HOME: $TARGET_DIR"

cd "$TARGET_DIR"

has_dirs() {
  local d
  for d in "$@"; do
    [[ -d "$d" ]] || return 1
  done
  return 0
}

detect_workflow() {
  if has_dirs 01.read_processing 02.short-read_assembly 03.read_selection 04.long-read_assembly 05.assembly_curation 06.polishing 07.genome_qc 08.annotation logs 09.report; then
    echo "FunFluxL+"
  elif has_dirs 01.pre-processing 02.short-read-assembly 02.assembly 03.short-read-selection 03.post-processing 04.annotation logs report; then
    echo "FunFluxL+_legacy"
  elif has_dirs 01.pre-processing 02.assembly 03.post-processing 04.annotation logs report; then
    echo "FunFlux/FunFluxL"
  elif has_dirs 01.pre-processing 02.annotation logs report; then
    echo "Funnotator"
  else
    die "Target does not look like a supported FunFlux/FunFluxL/FunFluxL+/Funnotator output directory."
  fi
}

WORKFLOW="$(detect_workflow)"
TARGETS=()

case "$WORKFLOW" in
  FunFluxL+)
    TARGETS+=(
      # Filtered long reads and SPAdes guide-assembly intermediates
      "01.read_processing/long/*/filtered.fastq"
      "02.short-read_assembly/*/spades/K*"
      "02.short-read_assembly/*/spades/misc"
      "02.short-read_assembly/*/spades/pipeline_state"
      "02.short-read_assembly/*/spades/tmp"

      # Assembler-specific intermediates. Keep the standardized primary FASTA,
      # GFA, and provenance files; only disposable Flye work directories and
      # Hifiasm overlap/error-correction caches are selected here.
      "04.long-read_assembly/*/flye/00-assembly"
      "04.long-read_assembly/*/flye/10-consensus"
      "04.long-read_assembly/*/flye/20-repeat"
      "04.long-read_assembly/*/flye/30-contigger"
      "04.long-read_assembly/*/flye/40-polishing"
      "04.long-read_assembly/*/hifiasm/*.bin"

      # BUSCO and annotation intermediates
      "07.genome_qc/completeness_evaluation/*/tmp"
      "07.genome_qc/completeness_evaluation/*/busco_tmp"
      "08.annotation/*/funannotate/predict_misc"
      "08.annotation/*/funannotate/predict_results"
      "08.annotation/*/funannotate/annotate_misc"
      "08.annotation/antismash_db"
      "08.annotation/*/iprscan"
      "08.annotation/*/repeatmasking"
    )
    ;;
  FunFluxL+_legacy|FunFlux/FunFluxL)
    TARGETS+=(
      # FunFluxL filtered reads
      "01.pre-processing/*.fastq"
      "01.pre-processing/*.fq"
      "01.pre-processing/*.fastq.gz"
      "01.pre-processing/*.fq.gz"

      # SPAdes intermediates
      "02.assembly/*/K*"
      "02.assembly/*/misc"
      "02.assembly/*/pipeline_state"
      "02.assembly/*/tmp"

      # FunFluxL+ SPAdes guide-assembly intermediates
      "02.short-read-assembly/*/spades/K*"
      "02.short-read-assembly/*/spades/misc"
      "02.short-read-assembly/*/spades/pipeline_state"
      "02.short-read-assembly/*/spades/tmp"

      # Flye intermediates
      "02.assembly/*/00-assembly"
      "02.assembly/*/10-consensus"
      "02.assembly/*/20-repeat"
      "02.assembly/*/30-contigger"
      "02.assembly/*/40-polishing"

      # Hifiasm cache files can be large and are rebuilt automatically from the
      # retained primary GFA plus the original reads when the rule is rerun.
      "02.assembly/*/hifiasm/*.bin"

      # BUSCO temporary directories
      "03.post-processing/completeness_evaluation/*/tmp"
      "03.post-processing/completeness_evaluation/*/busco_tmp"

      # Bulky funannotate internals; final annotate_results and logfiles are kept.
      "04.annotation/*/funannotate/predict_misc"
      "04.annotation/*/funannotate/predict_results"
      "04.annotation/*/funannotate/annotate_misc"

      # Large standalone annotation/intermediate working directories
      "04.annotation/antismash_db"
      "04.annotation/*/iprscan"
      "04.annotation/*/repeatmasking"
    )
    ;;
  Funnotator)
    TARGETS+=(
      # BUSCO temporary directories
      "01.pre-processing/completeness_evaluation/*/tmp"
      "01.pre-processing/completeness_evaluation/*/busco_tmp"

      # Bulky funannotate internals; final annotate_results and logfiles are kept.
      "02.annotation/*/funannotate/predict_misc"
      "02.annotation/*/funannotate/predict_results"
      "02.annotation/*/funannotate/annotate_misc"

      # Large standalone annotation/intermediate working directories
      "02.annotation/antismash_db"
      "02.annotation/*/iprscan"
      "02.annotation/*/repeatmasking"
    )
    ;;
esac

shopt -s nullglob dotglob

to_delete=()
for pattern in "${TARGETS[@]}"; do
  matches=( $pattern )
  if (( ${#matches[@]} )); then
    to_delete+=( "${matches[@]}" )
  fi
done

echo "== Fungal BioFlux output cleanup =="
echo "Workflow: $WORKFLOW"
echo "Target: $TARGET_DIR"
if [[ "$DO_RUN" -eq 1 ]]; then
  echo "Mode: RUN (will delete)"
else
  echo "Mode: DRY-RUN (no deletion)"
fi
echo

if (( ${#to_delete[@]} == 0 )); then
  echo "Nothing to clean; no matching files or directories were found."
  exit 0
fi

echo "Targets:"
for path in "${to_delete[@]}"; do
  echo "  - $path"
done
echo

if [[ "$DO_RUN" -eq 1 ]]; then
  rm -rf -- "${to_delete[@]}"
  echo "Cleanup finished."
else
  echo "Dry-run only. Re-run with --run to delete."
fi
