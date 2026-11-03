# FunFlux v2.0.0 — contigs front end: make a genome you already have usable
# (rules/contigs/20_assembly.smk).
#
# This mode annotates genomes that are ALREADY assembled: the user hands over FASTA
# files, not reads, so there is no read QC and no assembler. The front end is one rule
# that tidies the FASTA (filter_contigs) and, when the contamination and mitochondrion
# screen is switched on (parameters.contigs.screen: true; see docs/modes/contigs.md),
# one rule that passes the screened genome on to the mitochondrion split (finalize_assembly).
# The Snakefile includes this file in contigs mode only, right after 00_common.smk.
#
# What happens to a sample in this mode. {sample} is a wildcard: a placeholder that
# Snakemake fills in with each sample's name, so every rule below runs once per sample.
# {sample}.fasta is the user's own file in input.contigs_dir (.fa and .fna are accepted
# too); every file made from it in the two chains below lies under 02.assembly/{sample}/
# in the output folder.
#
#   screen: false (the default)
#     {sample}.fasta --filter_contigs--> contigs_final.fasta  (FINAL_CONTIGS)
#     That file goes straight to QUAST, BUSCO, ITS extraction and funannotate. No other
#     genome file is made: no contigs_all.fasta, no contaminants/ or mitochondria/
#     folder.
#
#   screen: true
#     {sample}.fasta --filter_contigs--> contigs_filt.fasta  (DRAFT_CONTIGS)
#       --> shared/10_decontam.smk: coverage track, BLAST, BlobTools, taxonomy selector
#       --> contigs_decontam.fasta  (DECONTAM_CONTIGS)
#       --finalize_assembly--> contigs_all.fasta  (CONTIGS_ALL)
#       --> shared/15_mitochondria.smk splits it into contigs_final.fasta (the nuclear
#           genome) and contigs_mitochondrial.fasta.
#     This is the same chain the read-based modes run after their assembler.
#
# Everything referenced here is defined once in 00_common.smk: CONTIGS_DIR,
# CONTIGS_EXTN, DRAFT_CONTIGS, DECONTAM_CONTIGS, CONTIGS_ALL, CONTIGS_SCREEN, LOGS,
# SPADES_HEADER_REGEX and the three awk programs FASTA_LIN_CMD, FASTA_SEL_CMD and
# FASTA_HEAD_CMD.

# The user's own FASTA for one sample: <input.contigs_dir>/<sample>.<extension>. The
# extension is the one 00_common.smk found and accepted for this batch (every sample of
# a batch must use the same one). No rule makes this file; sample discovery found it by
# listing the folder, so a contigs run enters the workflow here. Only this file uses the
# pattern, so it is not kept in 00_common.smk.
_INPUT_FASTA = os.path.join(CONTIGS_DIR, "{sample}." + CONTIGS_EXTN)


# ───────────── filter_contigs — tidy the supplied FASTA ─────────────
# A genome can come from anywhere, so the rule LOOKS at the first header and picks one of
# two treatments.
#
#   SPAdes-style header (">NODE_12_length_45678_cov_23.4"): the length and the read
#   coverage of every contig are written in its name, so the cut-offs of the Illumina
#   mode are applied here too: keep contigs of at least 500 bp and at least 2x coverage.
#   The headers are left exactly as they were; each sequence is written on one line.
#   The test is strict on purpose. SPADES_HEADER_REGEX in 00_common.smk, a regular
#   expression (a text pattern), asks for the whole NODE_<n>_length_<n>_cov_ start: a
#   header such as ">scaffold_1 length_12345" must not be taken for a SPAdes one, because
#   the numbers read out of it would be wrong.
#
#   Any other header (NCBI, Flye, a colleague's pipeline): those numbers do not exist,
#   and inventing a filter would silently delete real sequence. So only the headers are
#   cut at their first space or tab, which every later step needs (BLAST, BlobTools and
#   the mitochondrion screen all use the first word as the contig name), and NOTHING is
#   removed. The 500 bp floor is applied anyway, later, by `funannotate clean` in the
#   annotation (rule funannotate_preprocess).
#
# Takes in:    {sample}.<ext> from input.contigs_dir (see _INPUT_FASTA above).
# Does:        decides between the two treatments from the first header, writes the
#              tidied FASTA, and counts the records before and after.
# Produces:    DRAFT_CONTIGS. Without the screen that is contigs_final.fasta, the file
#              QUAST, BUSCO, ITS extraction and funannotate read; with the screen it
#              is contigs_filt.fasta, the draft the screen examines.
# Consumed by: with the screen, the screen rules that read the draft: map_contigs,
#              blast_contigs, blob_json and select_contigs (shared/10_decontam.smk) and
#              screen_mitochondria (shared/15_mitochondria.smk). Without it, the rules
#              that read the final genome: genome_assembly_evaluation and completeness
#              (QUAST and BUSCO, shared/20_qc.smk), its_windows (shared/30_taxonomy.smk)
#              and funannotate_preprocess (shared/40_annotation.smk).
# Log:         which of the two treatments ran, and how many records were kept out of
#              how many. That is the record of whether anything was filtered at all.
#
# The rule stops with an error when no contig is left, with the counts in the message:
# an empty genome would otherwise travel on and fail hours later, inside funannotate.
#
# (v1 differences: Funnotator removed every contig of 500 bp or less from any FASTA and
# did not trim headers. FunFlux applied the 500 bp and 2x cut-offs to its own SPAdes
# assembly only. The strict header test and the empty-result stop are new, and the rule
# now has a log.)
#
# `|| true` after the grep commands is required. grep ends with status 1 when it finds
# nothing, and Snakemake runs every shell command with `set -e` (stop at the first
# command that fails), so without it the rule would die on every non-SPAdes genome,
# which is exactly the case the second treatment is for. The record counts need it for
# the same reason: grep -c also ends with status 1 when the count is zero.
rule filter_contigs:
    input:
        contigs = _INPUT_FASTA,
    output:
        contigs = DRAFT_CONTIGS,
    log:
        LOGS + "/filter_contigs_{sample}.log"
    shell:
        """
        # Count the records (lines starting with ">") before filtering, for the log.
        total=$(grep -c '^>' {input.contigs} || true)

        # Look at the first header only: SPAdes writes the same style for every
        # record, so one line is enough to decide.
        spades_header=$(head -n 1 {input.contigs} | grep -E {SPADES_HEADER_REGEX:q} || true)

        if [ -n "$spades_header" ]; then
            echo "SPAdes-style header detected for {wildcards.sample}: removing contigs shorter than 500 bp or below 2x coverage; headers kept unchanged." > {log}
            # Join each sequence onto one line, then keep a record only if its header,
            # split at "_", has a length (field 4) of at least 500 and a coverage
            # (field 6) of at least 2.
            awk {FASTA_LIN_CMD:q} {input.contigs} | \
            awk -F"_" {FASTA_SEL_CMD:q} > {output.contigs}
        else
            echo "Header of {wildcards.sample} is not SPAdes-style: trimming headers to their first word; NO length or coverage filtering applied." > {log}
            # Cut every header at its first space or tab; copy sequence lines unchanged.
            awk {FASTA_HEAD_CMD:q} {input.contigs} > {output.contigs}
        fi

        # Count the records that are left and log both counts; stop if none is left.
        kept=$(grep -c '^>' {output.contigs} || true)
        echo "Kept $kept of $total records." >> {log}

        if [ "$kept" -eq 0 ]; then
            echo "[FunFlux] filter_contigs: no contig of {wildcards.sample} is left (kept $kept of $total records). Check the input FASTA, and see {log}." | tee -a {log} >&2
            exit 1
        fi
        """


# ───────────── finalize_assembly — name the screened genome ─────────────
# Defined only when the screen is on: the `if` below is decided once, when Snakemake
# reads the workflow, before any job runs. With the screen off, filter_contigs has
# already written the final genome and nothing remains to do.
#
# Why a rule that only copies a file: every mode ends its assembly chain with the same
# file name, contigs_all.fasta, so the mitochondrion split and everything after it read
# one name whichever way the genome arrived. In the longread and hybrid modes this is the
# step that picks which version becomes contigs_all.fasta, the polished one where
# polishing ran (FINALIZE_SOURCE in 00_common.smk); here there is nothing to pick, so the
# screened genome is simply copied.
#
# Takes in:    DECONTAM_CONTIGS, the genome that rule select_contigs
#              (shared/10_decontam.smk) wrote after applying the taxonomy selector
#              (every contig, when parameters.decontamination.mode is off, the default).
# Does:        copies it and records the source in the log.
# Produces:    CONTIGS_ALL, the complete screened genome, mitochondrion included.
# Consumed by: rule separate_mitochondria in shared/15_mitochondria.smk, which splits it
#              into the nuclear genome (contigs_final.fasta) and the mitochondrial
#              contigs.
if CONTIGS_SCREEN:

    rule finalize_assembly:
        input:
            contigs = DECONTAM_CONTIGS,
        output:
            contigs = CONTIGS_ALL,
        log:
            LOGS + "/finalize_assembly_{sample}.log"
        shell:
            """
            echo "contigs_all.fasta of {wildcards.sample} is a copy of {input.contigs}" > {log}
            cp {input.contigs} {output.contigs}
            """
