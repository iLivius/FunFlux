# FunFlux v2.0.0 — hybrid front end, step 2: the Illumina GUIDE assembly, SPAdes and the
# contig filter (rules/hybrid/20_assembly.smk).
#
# Identical to rules/illumina/20_assembly.smk: keep in sync. The two rules below,
# illumina_assembly and filter_contigs, are copies of the first two rules of that file,
# comments and commands alike. The only text that differs is two output names and their
# folders: SPAdes writes into GUIDE_SPADES_DIR (02.assembly/{sample}/guide/spades/) instead
# of SPADES_DIR, and the filter writes GUIDE_DRAFT (02.assembly/{sample}/guide/contigs_filt.fasta)
# instead of DRAFT_CONTIGS. That file's third rule, finalize_assembly, has no copy here: in
# hybrid mode the finished genome is named at the end of the long-read chain
# (rules/hybrid/50_polish.smk). Whoever edits one file edits the other.
#
# Why an Illumina assembly that is not delivered. In hybrid mode the genome is built from
# the long reads (rules/hybrid/40_long_assembly.smk). The Illumina reads are assembled all
# the same, into the GUIDE: a short-read assembly of the same isolate that exists for
# two reasons. (1) Its contamination screen (rules/hybrid/25_guide.smk) decides which
# Illumina pairs belong to the isolate; only those pairs then score the ONT reads in
# Filtlong and correct the final assembly with Polypolish. (2) QUAST reports the guide's
# statistics next to the final genome's (label {sample}_guide), so the report shows what
# the long reads added. Nothing after the screen and QUAST reads the guide.
#
# The Snakefile includes every file in rules/hybrid/ in hybrid mode only, in name order;
# this file comes after 10_reads.smk and before 25_guide.smk.
#
# What happens to a sample. {sample} is a wildcard: a placeholder that Snakemake fills in
# with each sample's name, so both rules run once per sample. Everything lies under
# 02.assembly/{sample}/guide/ in the output folder.
#
#   {sample}_trim_R1.fastq, _R2  (TRIM_R1/R2, from rule trim_adapters in rules/hybrid/10_reads.smk)
#     --illumina_assembly (SPAdes --isolate)--> guide/spades/contigs.fasta  (and the whole spades/ tree)
#     --filter_contigs (>= 500 bp, >= 2x)------> guide/contigs_filt.fasta   (GUIDE_DRAFT)
#       --> rules/hybrid/25_guide.smk: bowtie2 coverage track, Qualimap, BLAST, BlobTools,
#           mitochondrion screen, taxonomy selector --> guide/contaminants/contigs_decontam.fasta
#           (GUIDE_DECONTAM), then the read selection --> the selected pairs (SEL_R1, SEL_R2)
#
# Where the rule comments below name the illumina mode's neighbours, read the hybrid
# ones: rule trim_adapters is in rules/hybrid/10_reads.smk; the readers of the filtered
# guide are index_guide_contigs, blast_guide_contigs, blob_guide_json,
# screen_guide_mitochondria and select_guide_contigs in rules/hybrid/25_guide.smk (not the
# shared screen rules, which screen the long-read assembly); and the longest step of this
# mode is the long-read assembler, not SPAdes. The SPAdes settings (--isolate, the k-mer
# list parameters.spades_kmers, the memory limit resources.ram_gb) and the filter's cut-offs
# are the same as in illumina mode, and docs/modes/hybrid.md says so.
#
# Everything referenced here is defined once in 00_common.smk: TRIM_R1, TRIM_R2,
# GUIDE_SPADES_DIR, SPADES_KMERS, GUIDE_DRAFT, FASTA_LIN_CMD, FASTA_SEL_CMD, LOGS, CPUS and
# RAM.
#
# conda: Snakemake resolves a rule's environment path relative to the rule file, not to
# the working directory, so the files in workflow/envs/ are written "../../envs/<tool>.yaml".
# Only SPAdes gets an environment: the filter is awk, from the computer.
#
# CPUs: a rule that uses several CPUs declares Snakemake's own `threads:` and passes
# {threads} to its tool. That is what makes `--cores N` actually limit the run.

# The unfiltered contigs SPAdes writes inside its own folder. Two rules of this file use
# the path (illumina_assembly writes it, filter_contigs reads it) and no other file does,
# so it is defined here and not in 00_common.smk. SPAdes also writes scaffolds.fasta
# (contigs joined across gaps with runs of N where the read pairs say they belong
# together); FunFlux uses the contigs and never reads the scaffolds, so that nothing
# downstream has to deal with stretches of N.
_SPADES_CONTIGS = GUIDE_SPADES_DIR + "/contigs.fasta"


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
# Produces:    the whole SPAdes working folder, GUIDE_SPADES_DIR (02.assembly/{sample}/guide/spades/),
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
# GUIDE_DRAFT in 00_common.smk).
#
# Threads: all CPUS. The assembly is the longest step of this mode.
rule illumina_assembly:
    input:
        r1 = TRIM_R1,
        r2 = TRIM_R2,
    output:
        assembly_dir = directory(GUIDE_SPADES_DIR),
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
# Produces:    GUIDE_DRAFT, 02.assembly/{sample}/guide/contigs_filt.fasta, one line per
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
        contigs = GUIDE_DRAFT,
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
