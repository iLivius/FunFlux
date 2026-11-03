# FunFlux v2.0.0 — illumina front end, assembly: SPAdes, the contig filter and the
# hand-off to the mitochondrion split (rules/illumina/20_assembly.smk).
#
# The trimmed, PhiX-free read pairs of rules/illumina/10_reads.smk are assembled into
# contigs with SPAdes; the contigs that are too short or too thinly covered to be trusted
# are removed; the rest is the DRAFT that the contamination and mitochondrion screen
# examines. The screen itself is shared by every mode (shared/10_decontam.smk and
# shared/15_mitochondria.smk). The last rule of this file, finalize_assembly, takes what
# the screen kept and gives it the name every mode's assembly chain ends with, so that
# the mitochondrion split and everything after it read one file name whichever way the
# genome was made.
#
# The Snakefile includes every file in rules/illumina/ in illumina mode only, right after
# 00_common.smk, in name order: rules/illumina/10_reads.smk first, then this file.
# rules/hybrid/20_assembly.smk holds copies of the first two rules below for the hybrid
# mode's Illumina "guide" assembly (its header says so): a change made here must be made
# there too.
#
# What happens to a sample. {sample} is a wildcard: a placeholder that Snakemake fills in
# with each sample's name, so every rule below runs once per sample. Everything lies under
# 02.assembly/{sample}/ in the output folder.
#
#   {sample}_trim_R1.fastq, _R2  (TRIM_R1/R2, from rule trim_adapters)
#     --illumina_assembly (SPAdes --isolate)--> spades/contigs.fasta   (and the whole spades/ tree)
#     --filter_contigs (>= 500 bp, >= 2x)------> contigs_filt.fasta    (DRAFT_CONTIGS)
#       --> shared/10_decontam.smk: bowtie2 coverage track, BLAST, BlobTools, taxonomy selector
#       --> contaminants/contigs_decontam.fasta  (DECONTAM_CONTIGS)
#     --finalize_assembly (a copy)-------------> contigs_all.fasta     (CONTIGS_ALL)
#       --> shared/15_mitochondria.smk splits it into contigs_final.fasta (the nuclear
#           genome; QUAST, BUSCO, ITS extraction and funannotate read this) and
#           contigs_mitochondrial.fasta.
#
# The contamination screen always runs in this mode (SCREEN is always true for Illumina
# reads), so all three rules below always exist.
#
# Everything referenced here is defined once in 00_common.smk: TRIM_R1, TRIM_R2,
# SPADES_DIR, SPADES_KMERS, DRAFT_CONTIGS, DECONTAM_CONTIGS, CONTIGS_ALL, FASTA_LIN_CMD,
# FASTA_SEL_CMD, LOGS, CPUS and RAM.
#
# conda: Snakemake resolves a rule's environment path relative to the rule file, not to
# the working directory, so the files in workflow/envs/ are written "../../envs/<tool>.yaml".
# Only SPAdes gets an environment: the filter is awk and the hand-off is cp, both from the
# computer.
#
# CPUs: a rule that uses several CPUs declares Snakemake's own `threads:` and passes
# {threads} to its tool. That is what makes `--cores N` actually limit the run.

# The unfiltered contigs SPAdes writes inside its own folder. Two rules of this file use
# the path (illumina_assembly writes it, filter_contigs reads it) and no other file does,
# so it is defined here and not in 00_common.smk. SPAdes also writes scaffolds.fasta
# (contigs joined across gaps with runs of N where the read pairs say they belong
# together); FunFlux uses the contigs and never reads the scaffolds, so that nothing
# downstream has to deal with stretches of N.
_SPADES_CONTIGS = SPADES_DIR + "/contigs.fasta"


# ───────────── illumina_assembly — SPAdes, the de novo assembly ─────────────
# Takes in:    TRIM_R1 and TRIM_R2, the PhiX-free, adapter-trimmed read pairs from rule
#              trim_adapters (rules/illumina/10_reads.smk).
# Does:        spades.py assembles them, with
#                --isolate        SPAdes' setting for a single organism sequenced to high,
#                                 even coverage (a pure culture). It skips SPAdes' own read
#                                 error correction step (the same as --only-assembler) and
#                                 rules out the --careful mismatch correction: with deep,
#                                 even coverage the correct k-mers outnumber the error
#                                 k-mers in the graph anyway, so on isolate data those steps
#                                 cost hours and change little. The reads arrive trimmed by
#                                 fastp.
#                -k               the k-mer sizes SPAdes builds its graphs with, one graph
#                                 per size, each refining the last: small k joins regions of
#                                 low coverage, large k (up to 127) resolves repeats that a
#                                 short k cannot tell apart. The list is
#                                 parameters.spades_kmers (SPADES_KMERS, checked by
#                                 00_common.smk; default 21,33,55,77,99,127, the list of
#                                 the v1 workflows; `auto` lets SPAdes choose from the read
#                                 length, section 9g says when). It is passed through `params:` on
#                                 purpose: a changed list re-runs the assembly, as it should,
#                                 since the result depends on it.
#                --pe1-1/--pe1-2  the two mates of paired-end library 1 (the only library).
#                -t               the threads; OMP_NUM_THREADS tells the parts of SPAdes
#                                 that use OpenMP (a way for a program to run on several
#                                 CPUs) the same number, so they do not start one thread
#                                 per CPU of the computer whatever --cores says.
#                -m               the memory limit in GB; see "Memory" below.
# Produces:    the whole SPAdes working folder, SPADES_DIR (02.assembly/{sample}/spades/),
#              declared as a directory() output so that Snakemake treats the folder as the
#              result, and _SPADES_CONTIGS inside it, the file the next rule needs. The
#              folder holds one sub-folder per k-mer size (K21, K33, ...), the scaffolds,
#              the assembly graph and SPAdes' own spades.log; it is kept, as the record of
#              how the assembly was made. (Snakemake accepts a file output inside the same
#              rule's directory() output.)
# Consumed by: filter_contigs (contigs.fasta only).
# Log:         SPAdes' console output (the stages and their timings); the detailed log is
#              spades/spades.log.
#
# Memory. SPAdes is told how much memory it may use (-m, in GB) and stops with a clear
# message if it needs more; the value is resources.ram_gb (RAM) from the configuration.
# It is passed as a Snakemake `resources:` entry (mem_gb), not as a param: a memory cap
# changes nothing in the result, so a changed ram_gb must not re-run an assembly of several
# hours, and `snakemake --resources mem_gb=N` can limit how many assemblies run at once.
# N must be at least resources.ram_gb: Snakemake scales a job's resource down to the
# global limit, so a smaller N also lowers the -m cap SPAdes is given. On a cluster the job
# asks the scheduler for the same memory (mem_mb, the same number in MB) and for the time
# limit of resources.runtime_hours (runtime, in minutes); a single computer ignores both
# (00_common.smk, section 4).
#
# CAUTION, about where the next file goes: Snakemake deletes a directory() output before it
# re-runs the rule that makes it, so a re-run of SPAdes empties spades/. That is why
# filter_contigs writes contigs_filt.fasta NEXT TO spades/, never inside it (see
# DRAFT_CONTIGS in 00_common.smk).
#
# Threads: all CPUS. The assembly is the longest step of this mode.
rule illumina_assembly:
    input:
        r1 = TRIM_R1,
        r2 = TRIM_R2,
    output:
        assembly_dir = directory(SPADES_DIR),
        contigs = _SPADES_CONTIGS,
    params:
        kmers = SPADES_KMERS,
    threads: CPUS
    resources:
        mem_gb = RAM,
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/spades.yaml"
    log:
        LOGS + "/illumina_assembly_{sample}.log"
    shell:
        """
        OMP_NUM_THREADS={threads} \
        spades.py \
            -k {params.kmers} \
            --isolate \
            --pe1-1 {input.r1} \
            --pe1-2 {input.r2} \
            -o {output.assembly_dir} \
            -t {threads} \
            -m {resources.mem_gb} > {log} 2>&1
        """


# ───────────── filter_contigs — drop the short and the thinly covered contigs ─────────────
# SPAdes writes every path of its final graph as a contig, including the debris every
# assembly leaves: fragments under 500 bp (mostly unresolved repeat ends and read-error
# paths) and contigs at a coverage far below the genome's (sequencing errors, traces of
# a contaminant, leftover PhiX). Both numbers stand in the SPAdes header of every contig,
# ">NODE_<n>_length_<L>_cov_<C>", where cov is the k-mer coverage of the largest k, so the
# filter needs nothing but the headers. Contigs of at least 500 bp AND at least 2x are
# kept, headers unchanged; the cut-offs are fixed, not configuration keys (docs/modes/illumina.md
# says so).
#
# Takes in:    _SPADES_CONTIGS, the unfiltered contigs.fasta from illumina_assembly.
# Does:        two awk programs from 00_common.smk, in a pipe (the contigs mode's filter uses
#              the same two). FASTA_LIN_CMD joins each wrapped sequence onto one line, so
#              that FASTA_SEL_CMD can take a header and the line after it as one record and
#              apply the two cut-offs to the header's fields.
# Produces:    DRAFT_CONTIGS, 02.assembly/{sample}/contigs_filt.fasta, one line per
#              sequence, next to spades/ (never inside it: see the CAUTION above). Kept.
# Consumed by: the screen rules that read the draft: index_contigs, blast_contigs,
#              blob_json and select_contigs (shared/10_decontam.smk) and screen_mitochondria
#              (shared/15_mitochondria.smk).
# Log:         how many contigs were kept out of how many.
#
# The rule stops with an error when no contig is left, with the counts in the message: an
# empty draft would otherwise go on through the screen and fail hours later, inside BLAST or
# funannotate. A draft with nothing left means the reads or the assembly went wrong, not
# the filter.
#
# `|| true` after the grep commands is required. grep -c ends with status 1 when the count
# is zero, and Snakemake runs every shell command with `set -e` (stop at the first command
# that fails), so without it the rule would die on an empty file instead of reporting it.
rule filter_contigs:
    input:
        contigs = _SPADES_CONTIGS,
    output:
        contigs = DRAFT_CONTIGS,
    log:
        LOGS + "/filter_contigs_{sample}.log"
    shell:
        """
        # Count the records (lines starting with ">") before filtering, for the log.
        total=$(grep -c '^>' {input.contigs} || true)
        echo "Removing SPAdes contigs of {wildcards.sample} shorter than 500 bp or below 2x coverage; headers kept unchanged." > {log}

        # Join each sequence onto one line, then keep a record only if its header, split
        # at "_", has a length (field 4) of at least 500 and a coverage (field 6) of at
        # least 2.
        awk {FASTA_LIN_CMD:q} {input.contigs} | \
        awk -F"_" {FASTA_SEL_CMD:q} > {output.contigs}

        # Count the records that are left and log both counts; stop if none is left.
        kept=$(grep -c '^>' {output.contigs} || true)
        echo "Kept $kept of $total records." >> {log}

        if [ "$kept" -eq 0 ]; then
            echo "[FunFlux] filter_contigs: no contig of {wildcards.sample} is left (kept $kept of $total records). Check the SPAdes log (spades/spades.log next to this file) and the reads, and see {log}." | tee -a {log} >&2
            exit 1
        fi
        """


# ───────────── finalize_assembly — name the screened genome ─────────────
# Why a rule that only copies a file: every mode ends its assembly chain with the same file
# name, contigs_all.fasta, so the mitochondrion split and everything after it read one name
# whichever way the genome arrived. In the longread and hybrid modes this is the step that
# picks which version becomes contigs_all.fasta, the polished one where polishing ran; here
# there is no polishing (the Illumina reads already built the contigs), so the screened
# genome is simply copied.
#
# Takes in:    DECONTAM_CONTIGS, the genome that rule select_contigs (shared/10_decontam.smk)
#              wrote after applying the taxonomy selector (every contig, when
#              parameters.decontamination.mode is off, the default), with the confirmed
#              mitochondrial contigs kept whatever the genus filter says.
# Does:        copies it and records the source in the log.
# Produces:    CONTIGS_ALL, the complete screened genome, mitochondrion included. Kept; a
#              delivered file.
# Consumed by: rule separate_mitochondria in shared/15_mitochondria.smk, which splits it
#              into the nuclear genome (contigs_final.fasta) and the mitochondrial contigs.
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
