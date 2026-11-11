#!/usr/bin/env bash
# clean_funflux_output.sh
#
# Frees the space of a FINISHED FunFlux output folder: deletes the bulky working files that
# no result needs and keeps every result (the genomes, the annotation, the reports, the
# tables, the logs, the PROVENANCE files). Run it once you have checked the results and
# want to archive or share the folder. It is a dry run by default: it lists what it would
# delete and deletes nothing. Only --run deletes.
#
# Do not run FunFlux on a cleaned folder again. Some deleted files (the read mappings, the
# downloaded databases and, with --include-reads, the reads) are files that later steps
# read: Snakemake would make them again and repeat every step after them, which can be
# most of each sample.
#
# Usage:
#   clean_funflux_output.sh [OUTPUT_DIR]                        list what would be deleted
#   clean_funflux_output.sh --run [OUTPUT_DIR]                  delete it
#   clean_funflux_output.sh --run --include-reads [OUTPUT_DIR]  delete the kept reads too
# Without OUTPUT_DIR the current folder is used; -t, --target or --output-dir OUTPUT_DIR
# names it too. --dry-run cancels an earlier --run; -h or --help prints this text.
#
# FunFlux v2 (every mode) is recognised by its folders 00.resources, 02.assembly,
# 04.annotation and 06.report, and is cleaned only when every sample folder in 02.assembly/
# has 04.annotation/{sample}/.cleanup_complete.txt, the marker written once that sample's
# annotation is done. These show that each sample reached the end once; they cannot show
# whether a later run on the folder (a sample added, a setting changed) is unfinished or
# still running. Make sure no FunFlux run is using the folder before you clean it.
#
# What is deleted, and why no result needs it:
#   02.assembly/{sample}/
#     spades/K*, misc, pipeline_state, tmp    SPAdes' working folders (also guide/spades/ in
#                                             hybrid); the contigs, graphs and log stay
#     hifiasm/*.bin                           hifiasm's caches of corrected reads and overlaps;
#                                             the graphs stay. Snakemake empties hifiasm/
#                                             before any new assembly, so they are never reused
#     flye/00-assembly ... 40-polishing       Flye's stage folders; assembly.fasta,
#                                             assembly_info.txt and the graphs stay
#     contaminants/{sample}_map.bam, .bai     the reads mapped to the draft (also in
#                                             guide/contaminants/ in hybrid); BlobTools,
#                                             Qualimap and contig_depth have used them, and
#                                             their results stay
#     contaminants/contigs_decontam.fasta     the index files Medaka's aligner writes next to
#       .fai, .*.mmi                          the draft it polishes
#     medaka/calls_to_draft.bam, .bai,        Medaka's read alignment and the per-base scores
#       consensus_probs.hdf                   of its network; consensus.fasta stays
#     eval/busco/tmp                          BUSCO's copy of the lineage proteins
#     eval/busco/logs/miniprot_align_*_out.log  miniprot's alignments of the lineage proteins
#                                             (about 1 GB), which BUSCO turned into
#                                             full_table.tsv; busco.log stays
#     eval/busco/run_*/miniprot_output,       miniprot's index of the genome and its working
#       auto_lineage/run_*/miniprot_output    files
#   04.annotation/{sample}/
#     funannotate/predict_misc/genemark, EVM  GeneMark-ES's and EVidenceModeler's working folders
#     funannotate/predict_misc/tbl2asn/ and   tbl2asn's copies of the genome and its annotation,
#       annotate_misc/tbl2asn/ (and 1/ in     the same files as in predict_results/ and
#       each): genome.fsa, .gbf, .sqn, .tbl   annotate_results/. Its validation reports (*.val,
#                                             genome.ecn, discrepency.report.txt) stay: those
#                                             of the final annotation exist only there
#     funannotate/predict_misc/               the protein evidence funannotate copied from its
#       proteins.combined.fa                  database. predict_misc/ab_initio_parameters/,
#                                             which predict_results/*.parameters.json points
#                                             to, stays
#     funannotate/annotate_misc/iprscan.xml   funannotate's copy of the InterProScan XML (a few
#                                             GB); iprscan/interproscan_out.xml stays
#     iprscan/tmp                             InterProScan's scratch folder
#     contigs_sort.fasta.mmi, .fai            index files that tools run on the annotated genome
#                                             leave next to it (FunFlux writes none)
#   00.resources/ (one per run)
#     antismash/databases                     the antiSMASH database FunFlux downloaded (9 GB).
#                                             With directories.antismash_db set, the folder
#                                             holds only links to your copy and is kept, as
#                                             deleting it would free nothing
#     unite/unite_its_sintax.fasta,           the UNITE reference and its filtered copy;
#       unite_its_sintax.identified.fasta     PROVENANCE.txt (address, date, md5) stays. With
#                                             directories.unite_db set, the reference is a link
#                                             to your copy and is kept
#     busco/lineages                          the BUSCO lineage data FunFlux downloaded when
#                                             parameters.busco.lineage names a lineage (a few
#                                             hundred MB); busco/file_versions.tsv, the list of
#                                             the release used, stays
#     signalp6/package                        SignalP 6 with its model weights (1.6 GB). It is
#                                             licensed by DTU for your own use and must not be
#                                             in a copy you share. Without it FunFlux cannot run
#                                             SignalP on this folder again: delete
#                                             00.resources/signalp6/ as a whole before a new run,
#                                             and SignalP 6 is installed again from your archive
#
# --include-reads (off by default) also deletes the reads FunFlux keeps:
#   01.reads/{sample}/long/{sample}_filtered.fastq               longread, hybrid
#   01.reads/{sample}/selected/{sample}_selected_R1.fastq, _R2    hybrid
# They are the largest files of a read-based run, but they are not your raw data: FunFlux
# reads the raw reads where they are and never copies them, so these can be made again.
# They are kept unless you ask, because they are the exact reads the genome was assembled
# and polished from (assembly_provenance.tsv names the filtered file): keep them if you may
# want to polish again or map the reads to the genome.
#
# Never deleted: the genomes, the annotation files (GFF3, GenBank, proteins), the reports,
# the tables, the logs, the PROVENANCE files, the QC folders, and anything outside the
# output folder (the funannotate database; .snakemake/ in the folder FunFlux was started
# from). This holds even when a folder in the output folder is a link to another place,
# such as a 00.resources/ that several runs share: paths that lie there are skipped, with a
# note.
#
# The v1 layouts (FunFlux, FunFluxL, FunFluxL+, Funnotator) are still recognised and cleaned
# with their own lists, unchanged; --include-reads does not apply to them.

set -euo pipefail

DO_RUN=0
INCLUDE_READS=0
TARGET_DIR=""

# --help prints the header above: every line from the second one down to the line before
# "set -euo pipefail". The end is found by that text, not by a line number, so the help
# stays whole when the header changes.
usage() {
  sed -n '2,/^set -euo pipefail/p' "$0" | sed '$d'
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
    --include-reads)
      INCLUDE_READS=1
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

# True when the folder holds at least one entry that is not a link. A downloaded antiSMASH
# database has real folders in it; one built from directories.antismash_db has only links
# to your copy. The marker .snakemake_timestamp, which Snakemake leaves in a folder output,
# does not count. A folder that does not exist holds nothing.
folder_has_real_entries() {
  local entry
  for entry in "$1"/*; do
    # Without a match the pattern comes back as it is written, which does not exist.
    [[ -e "$entry" || -L "$entry" ]] || continue
    if [[ "$(basename -- "$entry")" == ".snakemake_timestamp" ]]; then
      continue
    fi
    if [[ ! -L "$entry" ]]; then
      return 0
    fi
  done
  return 1
}

# Which workflow wrote the folder, told by its numbered folders. FunFlux v2 is tested
# first. Its four marker folders exist in every mode; 01.reads (absent in contigs mode) and
# 05.amr (only with the AMR screen on) are not markers for that reason. The v1 FunFlux and
# FunFluxL folders have 02.assembly and 04.annotation too, but no v1 layout has 06.report,
# which keeps the two apart.
detect_workflow() {
  if has_dirs 00.resources 02.assembly 04.annotation 06.report; then
    echo "FunFlux v2"
  elif has_dirs 01.read_processing 02.short-read_assembly 03.read_selection 04.long-read_assembly 05.assembly_curation 06.polishing 07.genome_qc 08.annotation logs 09.report; then
    echo "FunFluxL+"
  elif has_dirs 01.pre-processing 02.short-read-assembly 02.assembly 03.short-read-selection 03.post-processing 04.annotation logs report; then
    echo "FunFluxL+_legacy"
  elif has_dirs 01.pre-processing 02.assembly 03.post-processing 04.annotation logs report; then
    echo "FunFlux/FunFluxL"
  elif has_dirs 01.pre-processing 02.annotation logs report; then
    echo "Funnotator"
  else
    die "Target does not look like a supported output directory: a finished FunFlux v2 run (with 06.report/), or a FunFlux/FunFluxL/FunFluxL+/Funnotator v1 run."
  fi
}

WORKFLOW="$(detect_workflow)"

# A FunFlux v2 folder is cleaned only when every sample in it has finished. 06.report alone
# does not show it: MultiQC does not wait for every step of a sample, and a later run on a
# finished folder (a sample added) finds 06.report already there. So each sample folder in
# 02.assembly/ must have the marker rule cleanup writes once that sample's annotation is
# done, 04.annotation/{sample}/.cleanup_complete.txt. A sample without it is still
# running, has failed, or was added by a run that has not finished, and its working files
# may still be needed: the script stops before listing anything.
if [[ "$WORKFLOW" == "FunFlux v2" ]]; then
  unfinished_samples=()
  for sample_dir in 02.assembly/*/; do
    # If 02.assembly/ is empty, the pattern comes back as it is written, which is no folder.
    [[ -d "$sample_dir" ]] || continue
    sample="$(basename -- "$sample_dir")"
    if [[ ! -e "04.annotation/$sample/.cleanup_complete.txt" ]]; then
      unfinished_samples+=( "$sample" )
    fi
  done
  if (( ${#unfinished_samples[@]} )); then
    die "Not every sample of this FunFlux v2 run has finished: 04.annotation/<sample>/.cleanup_complete.txt is missing for: ${unfinished_samples[*]}. Clean the folder once the run is complete."
  fi
fi

TARGETS=()

# The reads a FunFlux v2 run keeps (longread and hybrid modes), deleted only with
# --include-reads; the header says why. Listed once here, because the note near the end
# looks for the same files.
V2_READ_PATTERNS=(
  "01.reads/*/long/*_filtered.fastq"
  "01.reads/*/selected/*_selected_R1.fastq"
  "01.reads/*/selected/*_selected_R2.fastq"
)

case "$WORKFLOW" in
  "FunFlux v2")
    # The header explains each group. Every pattern is relative to the output folder and
    # "*" stands for the sample name; a pattern that matches nothing (a mode or step that
    # did not run) is simply dropped.
    TARGETS+=(
      # The assemblers' working files. Snakemake treats spades/, hifiasm/ and flye/ as
      # whole-folder outputs and empties them before any new assembly, so nothing in them
      # is ever reused.
      "02.assembly/*/spades/K*"
      "02.assembly/*/spades/misc"
      "02.assembly/*/spades/pipeline_state"
      "02.assembly/*/spades/tmp"
      "02.assembly/*/guide/spades/K*"
      "02.assembly/*/guide/spades/misc"
      "02.assembly/*/guide/spades/pipeline_state"
      "02.assembly/*/guide/spades/tmp"
      "02.assembly/*/hifiasm/*.bin"
      "02.assembly/*/flye/00-assembly"
      "02.assembly/*/flye/10-consensus"
      "02.assembly/*/flye/20-repeat"
      "02.assembly/*/flye/30-contigger"
      "02.assembly/*/flye/40-polishing"

      # The reads mapped to the draft: BlobTools took its coverage from them and Qualimap
      # its report; both results stay. (In contigs mode the BAM is temporary and already
      # gone.)
      "02.assembly/*/contaminants/*_map.bam"
      "02.assembly/*/contaminants/*_map.bam.bai"
      "02.assembly/*/guide/contaminants/*_map.bam"
      "02.assembly/*/guide/contaminants/*_map.bam.bai"

      # Medaka (ONT samples): the samtools and minimap2 indexes its aligner writes next to
      # the draft, then its own alignment and network scores. consensus.fasta stays.
      "02.assembly/*/contaminants/contigs_decontam.fasta.fai"
      "02.assembly/*/contaminants/contigs_decontam.fasta.*.mmi"
      "02.assembly/*/medaka/calls_to_draft.bam"
      "02.assembly/*/medaka/calls_to_draft.bam.bai"
      "02.assembly/*/medaka/consensus_probs.hdf"

      # BUSCO's working files. The miniprot "log" is miniprot's alignment output, which
      # BUSCO links into miniprot_output/; the summaries, full_table.tsv,
      # missing_busco_list.tsv and busco.log stay. The auto_lineage pattern comes first, so
      # that the folder BUSCO also links as run_<lineage> is listed under its real name.
      "02.assembly/*/eval/busco/tmp"
      "02.assembly/*/eval/busco/logs/miniprot_align_*_out.log"
      "02.assembly/*/eval/busco/auto_lineage/run_*/miniprot_output"
      "02.assembly/*/eval/busco/run_*/miniprot_output"

      # funannotate's working files and copies. Not predict_misc/ as a whole: the
      # trained gene-finder models in predict_misc/ab_initio_parameters/ are what
      # predict_results/*.parameters.json points to.
      "04.annotation/*/funannotate/predict_misc/genemark"
      "04.annotation/*/funannotate/predict_misc/EVM"
      "04.annotation/*/funannotate/predict_misc/proteins.combined.fa"
      "04.annotation/*/funannotate/annotate_misc/iprscan.xml"

      # tbl2asn's copies of the genome (.fsa), its feature table (.tbl), GenBank file (.gbf)
      # and submission file (.sqn), twice in each folder (tbl2asn/ and tbl2asn/1/): the same
      # files as in predict_results/ and annotate_results/. Not the folders as a whole:
      # the validation reports next to them (*.val, genome.ecn, discrepency.report.txt)
      # stay, because annotate_results/ has no copy of the final annotation's reports.
      "04.annotation/*/funannotate/predict_misc/tbl2asn/genome.fsa"
      "04.annotation/*/funannotate/predict_misc/tbl2asn/genome.tbl"
      "04.annotation/*/funannotate/predict_misc/tbl2asn/genome.gbf"
      "04.annotation/*/funannotate/predict_misc/tbl2asn/genome.sqn"
      "04.annotation/*/funannotate/predict_misc/tbl2asn/1/genome.fsa"
      "04.annotation/*/funannotate/predict_misc/tbl2asn/1/genome.tbl"
      "04.annotation/*/funannotate/predict_misc/tbl2asn/1/genome.gbf"
      "04.annotation/*/funannotate/predict_misc/tbl2asn/1/genome.sqn"
      "04.annotation/*/funannotate/annotate_misc/tbl2asn/genome.fsa"
      "04.annotation/*/funannotate/annotate_misc/tbl2asn/genome.tbl"
      "04.annotation/*/funannotate/annotate_misc/tbl2asn/genome.gbf"
      "04.annotation/*/funannotate/annotate_misc/tbl2asn/genome.sqn"
      "04.annotation/*/funannotate/annotate_misc/tbl2asn/1/genome.fsa"
      "04.annotation/*/funannotate/annotate_misc/tbl2asn/1/genome.tbl"
      "04.annotation/*/funannotate/annotate_misc/tbl2asn/1/genome.gbf"
      "04.annotation/*/funannotate/annotate_misc/tbl2asn/1/genome.sqn"

      # InterProScan's scratch folder, and the index files a tool left next to the
      # annotated genome.
      "04.annotation/*/iprscan/tmp"
      "04.annotation/*/contigs_sort.fasta.mmi"
      "04.annotation/*/contigs_sort.fasta.fai"
    )

    # Run-level material that can be fetched or installed again. The antiSMASH database
    # and the UNITE reference are listed only when FunFlux downloaded them into this
    # folder. With directories.antismash_db set, databases/ holds only links to your copy,
    # and with directories.unite_db set the reference is a link to yours: deleting those
    # would free nothing and only make a later run repeat antiSMASH or the ITS
    # classification. (rm -rf would remove only the links, never your copy.)
    if folder_has_real_entries "00.resources/antismash/databases"; then
      TARGETS+=( "00.resources/antismash/databases" )
    fi
    if [[ ! -L "00.resources/unite/unite_its_sintax.fasta" ]]; then
      TARGETS+=( "00.resources/unite/unite_its_sintax.fasta" )
    fi
    # The BUSCO lineage data exist only when parameters.busco.lineage named a lineage:
    # rule busco_data downloaded them once for the run (in auto mode each BUSCO job
    # downloads its own data and deletes it when done). Rule busco_data would download
    # them again, possibly a newer release of the lineage.
    TARGETS+=(
      "00.resources/unite/unite_its_sintax.identified.fasta"
      "00.resources/busco/lineages"
      "00.resources/signalp6/package"
    )

    if [[ "$INCLUDE_READS" -eq 1 ]]; then
      TARGETS+=( "${V2_READ_PATTERNS[@]}" )
    fi
    ;;
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

# nullglob: a pattern whose wildcard matches nothing gives nothing, not the pattern itself.
# dotglob: "*" also matches names that start with a dot.
shopt -s nullglob dotglob

# Where a path really sits on disk: its folder with every link resolved, plus its own name.
# Two patterns can reach one folder under two names (BUSCO links run_<lineage> to
# auto_lineage/run_<lineage>); comparing these places lists that folder once.
physical_path() {
  local parent
  parent="$(cd -- "$(dirname -- "$1")" && pwd -P)"
  echo "$parent/$(basename -- "$1")"
}

# Turn the patterns into the paths that exist, each once. The list printed below is the
# list rm is given.
to_delete=()
skipped_outside=()
declare -A already_listed=()
for pattern in "${TARGETS[@]}"; do
  # $pattern is unquoted on purpose: that is what makes the shell expand the wildcards into
  # the paths they match. The patterns contain no spaces.
  # shellcheck disable=SC2206
  matches=( $pattern )
  # Bash older than 4.4 (still common on computing clusters) stops with "unbound variable"
  # when an empty list is read under set -u, so a pattern that matched nothing is passed
  # over here.
  if (( ${#matches[@]} == 0 )); then
    continue
  fi
  for match in "${matches[@]}"; do
    # nullglob does not touch a pattern without a wildcard (00.resources/signalp6/package),
    # which comes back as itself whether or not it exists. Keep only what is there; -L lets
    # a broken link count too.
    [[ -e "$match" || -L "$match" ]] || continue
    place="$(physical_path "$match")"
    # Never delete outside the output folder. A folder on the way can be a link to another
    # place (00.resources/ linked to a folder several runs share, a sample folder on
    # another disk), and the shell follows such links when it expands a pattern, so the
    # real place is checked here. A match that is itself a link passes: rm -rf removes the
    # link, never what it points to.
    if [[ "$place" != "$TARGET_DIR"/* ]]; then
      skipped_outside+=( "$match (really $place)" )
      continue
    fi
    [[ -z "${already_listed[$place]:-}" ]] || continue
    already_listed[$place]=1
    to_delete+=( "$match" )
  done
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
else
  echo "Targets:"
  for path in "${to_delete[@]}"; do
    echo "  - $path"
  done
  echo
fi

# Paths that matched a pattern but really lie outside the output folder (see the loop
# above): named here, so that the user sees why they are not in the list.
if (( ${#skipped_outside[@]} )); then
  echo "NOTE: skipped, because they lie outside the output folder (a folder on the way is a link):"
  for path in "${skipped_outside[@]}"; do
    echo "  - $path"
  done
  echo
fi

# Notes for FunFlux v2: the kept reads that were left out (the largest files of a
# read-based run, easy to miss), the licence of the SignalP 6 package, and the two
# warnings of the header about runs on this folder.
if [[ "$WORKFLOW" == "FunFlux v2" ]]; then
  if [[ "$INCLUDE_READS" -eq 0 ]]; then
    kept_reads=()
    for pattern in "${V2_READ_PATTERNS[@]}"; do
      # Unquoted on purpose, as above.
      # shellcheck disable=SC2206
      matches=( $pattern )
      # An empty list is skipped for older bash, as above.
      if (( ${#matches[@]} )); then
        kept_reads+=( "${matches[@]}" )
      fi
    done
    if (( ${#kept_reads[@]} )); then
      echo "NOTE: ${#kept_reads[@]} read file(s) in 01.reads/ are kept. Add --include-reads to delete them too."
      echo
    fi
  fi
  if [[ -d "00.resources/signalp6/package" ]]; then
    echo "NOTE: 00.resources/signalp6/package is SignalP 6, licensed by DTU for your own use; it must"
    echo "      not be in a copy you share. Without it FunFlux cannot run SignalP on this folder again:"
    echo "      delete 00.resources/signalp6/ as a whole before a new run."
    echo
  fi
  echo "NOTE: Check that no FunFlux run is using this folder. Once it is cleaned, do not run FunFlux"
  echo "      on it again (to add a sample or change a setting): the deleted mappings and databases"
  echo "      would be made again, and most of each sample's steps repeated."
  echo
elif [[ "$INCLUDE_READS" -eq 1 ]]; then
  echo "NOTE: --include-reads applies to FunFlux v2 folders only and was ignored."
  echo
fi

if (( ${#to_delete[@]} == 0 )); then
  exit 0
fi

if [[ "$DO_RUN" -eq 1 ]]; then
  rm -rf -- "${to_delete[@]}"
  echo "Cleanup finished."
else
  echo "Dry-run only. Re-run with --run to delete."
fi
