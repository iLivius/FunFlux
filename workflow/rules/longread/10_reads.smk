# FunFlux v2.1.0 — longread front end, step 1: read filtering and read QC
# (rules/longread/10_reads.smk).
#
# This mode assembles a genome from long reads alone: Oxford Nanopore (ONT) or PacBio
# HiFi, decided per sample by the file name ({sample}_ont.<ext> or {sample}_pacbio.<ext>
# in input.longread_dir; 00_common.smk, section 12). Both technologies take the same
# three steps here, with the same settings. Before them, ONT reads alone are trimmed of the
# sequencing adapter (rules/shared/05_ont_adapters.smk), so Filtlong reads the trimmed file
# for an ONT sample. After that, the technology matters at the assembler
# (rules/longread/20_assembly.smk) and in polishing (rules/longread/30_polish.smk).
# The Snakefile includes this file in longread mode only, right after 00_common.smk.
#
# What happens to a sample's reads. {sample} is a wildcard: a placeholder that Snakemake
# fills in with each sample's name, so every rule below runs once per sample. The raw
# file is the user's own; everything made from it lies under 01.reads/{sample}/long/ in
# the output folder.
#
#   {sample}_ont.fastq.gz (or _pacbio; plain .fastq is accepted too)
#     [ONT only: --find_ont_adapter, --trim_ont_adapters--> {sample}_trimmed.fastq.gz
#      (TRIMMED_LONG, temporary; rules/shared/05_ont_adapters.smk): the ONT adapter and the
#      barcode after it removed from the read ends]
#     --filter_long_reads--> {sample}_filtered.fastq  (FILTERED_LONG, kept)
#       --> rule long_read_assembly (20_assembly.smk): the reads the genome is built from
#       --> rule map_contigs (shared/10_decontam.smk): the reads aligned back to the draft,
#           the coverage track of the contamination screen and the BAM Qualimap evaluates
#       --> rule long_read_consensus (30_polish.smk): the reads Medaka polishes with (ONT)
#     --raw_long_read_qc------> raw_qc/       (NANOPLOT_RAW_DIR)      --> MultiQC
#     --filtered_long_read_qc-> filtered_qc/  (NANOPLOT_FILTERED_DIR) --> MultiQC
#
# Why filter at all. A long-read run holds reads from a few hundred bases to over a
# hundred kilobases, of very uneven quality. Short reads add little to an assembly
# (the assembler needs reads that span repeats), and the worst reads add errors. Filtlong
# drops the reads under a minimum length, ranks the rest by quality and length, and keeps
# the best share of the bases. The two NanoPlot reports, one before and one after, show
# what adapter trimming (ONT samples) and the filter did together: the read-length and
# quality distributions as delivered and as handed to the assembler. The settings live
# in parameters.longread.filtlong (00_common.smk, section 11) and are sized for ONT
# reads; for PacBio HiFi, which is already accurate, they are provisional
# (docs/modes/longread.md).
#
# The hybrid mode reuses this step: rules/hybrid/30_long_reads.smk holds a copy of the
# two NanoPlot rules and a Filtlong rule of its own, in which the Illumina reads of the
# sample guide the filter (docs/methods/filtlong-guidance.md).
#
# Everything referenced here is defined once in 00_common.smk: FILTERED_LONG,
# NANOPLOT_RAW_DIR, NANOPLOT_FILTERED_DIR, FILTLONG_MIN_LENGTH, FILTLONG_KEEP_PERCENT,
# FILTLONG_LENGTH_WEIGHT, FILTLONG_TARGET_ARG, LOGS and CPUS, and the functions
# raw_long_fastq() (the sample's raw long-read file, found by sample discovery),
# filtlong_input() (the file Filtlong starts from: the trimmed reads of an ONT sample, the raw
# file otherwise) and capped_cpus().
#
# conda: environment paths are relative to THIS file (workflow/rules/longread/), so
# "../../envs/filtlong.yaml" is workflow/envs/filtlong.yaml.
#
# CPUs: a rule that uses several CPUs declares Snakemake's own `threads:` and passes
# {threads} to its tool. That is what makes `--cores N` actually limit the run. Filtlong
# runs on one CPU and declares none.
#
# In the shell commands below, a doubled brace ({{ or }}) stands for a single brace in the
# shell. Snakemake reads single braces as its own placeholders, such as {input.reads}, and
# turns a doubled brace into a single one.


# ───────────── filter_long_reads — Filtlong: keep the best reads ─────────────
# Takes in:    filtlong_input: for an ONT sample, the reads fastplong trimmed (TRIMMED_LONG,
#              rule trim_ont_adapters in rules/shared/05_ont_adapters.smk); for a PacBio
#              sample, or with parameters.longread.ont_adapter off, the sample's RAW long-read
#              FASTQ, which no rule makes: sample discovery found it by listing
#              input.longread_dir. Plain or gzip-compressed; Filtlong reads both.
# Does:        filtlong with the settings of parameters.longread.filtlong:
#                --min_length    reads shorter than this (bp) are dropped outright
#                                (default 1000).
#                --keep_percent  of the remaining BASES, keep the best this-many percent
#                                (default 90): reads are ranked, the worst are removed
#                                until that share is left.
#                --length_weight how much a read's length counts, relative to its mean
#                                quality, in that ranking (default 10, ten times
#                                Filtlong's own default of 1: for assembly, a long read
#                                is worth more than a slightly cleaner short one).
#                --target_bases  an absolute cap on the bases kept, for a very deep run
#                                (params target_bases: "--target_bases <size>", or
#                                nothing when no cap is set; 00_common.smk builds the
#                                text as FILTLONG_TARGET_ARG). It is written without
#                                quotes on purpose: an empty value then leaves nothing
#                                behind, no bare flag and no empty argument.
# Produces:    FILTERED_LONG, {sample}_filtered.fastq, uncompressed and KEPT: the assembler,
#              the read mapping and Medaka all read it, and a re-run of any of them must
#              not have to filter the reads again.
# Consumed by: rule long_read_assembly (20_assembly.smk), rule map_contigs
#              (shared/10_decontam.smk), rule long_read_consensus (30_polish.smk) and rule
#              filtered_long_read_qc below. Rule all asks for it too.
# Log:         Filtlong's own report: how many reads and bases came in, how many were kept.
#
# v1 -> v2: FunFluxL 0.1 ran Filtlong with fixed values (--min_length 1000 --keep_percent 90
# --target_bases 500000000). Here the values come from parameters.<mode>.filtlong,
# --length_weight is passed too, and --target_bases is added only when a cap is set:
# 00_common.smk prepares that text once (FILTLONG_TARGET_ARG) and the command just writes
# it. The output moved from 01.pre-processing/{sample}_long_filt.fastq to the sample's
# long-read folder, and the log is named after the rule.
rule filter_long_reads:
    input:
        reads = filtlong_input,
    output:
        reads = FILTERED_LONG,
    params:
        min_length = FILTLONG_MIN_LENGTH,
        keep_percent = FILTLONG_KEEP_PERCENT,
        length_weight = FILTLONG_LENGTH_WEIGHT,
        target_bases = FILTLONG_TARGET_ARG,
    # Asked of a cluster scheduler: a time limit only, as the step reads every read but
    # holds little in memory (00_common.smk, section 4).
    resources:
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/filtlong.yaml"
    log:
        LOGS + "/filter_long_reads_{sample}.log"
    shell:
        """
        filtlong \
            --min_length {params.min_length} \
            --keep_percent {params.keep_percent} \
            --length_weight {params.length_weight} \
            {params.target_bases} \
            {input.reads} > {output.reads} 2>{log}
        """


# ───────────── raw_long_read_qc — NanoPlot on the reads as delivered ─────────────
# Takes in:    the sample's RAW long-read FASTQ (raw_long_fastq), as delivered. Filtlong
#              starts from the same file for a PacBio HiFi sample, or with ont_adapter:
#              off; for an ONT sample it starts from the reads fastplong trimmed.
# Does:        NanoPlot reads every record and draws the read-length and read-quality
#              distributions (histograms, length against quality, --loglength adds the
#              log-scaled length plots), and writes a NanoStats table of the numbers:
#              read count, total bases, mean and median length and quality, N50, and the
#              share of reads above each quality cut-off.
#              --no_static leaves out the static (PNG) copies of the plots: NanoPlot draws
#              those with a Chrome or Chromium browser installed on the computer (through
#              the Kaleido library), and on a computer or cluster node without one it
#              writes no PNG and only a warning. Without the option the output therefore
#              differed from machine to machine; with it, every machine writes the same
#              interactive HTML plots. MultiQC reads only the NanoStats table.
#              --prefix "{sample}_" puts the sample name in front of every file name,
#              so the table is {sample}_NanoStats.txt. MultiQC (rule multiqc in
#              shared/90_report.smk) takes the sample name from that file name; without
#              the prefix every sample's table would be called NanoStats.txt and the rows
#              could not be told apart.
# Produces:    NANOPLOT_RAW_DIR (01.reads/{sample}/long/raw_qc/), declared as a directory()
#              output, so Snakemake treats the whole folder as the result and starts a
#              re-run from an empty folder; NanoPlot names the files in it (the NanoStats
#              table, an HTML report and the plots as HTML pages).
# Consumed by: MultiQC (rule multiqc in shared/90_report.smk), which reads the NanoStats
#              table. Rule all asks for the folder too, so the job list shows this step.
# Log:         NanoPlot's own messages.
#
# NanoPlot is only a report. No later step needs its result, so a NanoPlot crash must not
# stop a batch that may spend days on assembly and annotation. The output folder is
# created first (so the rule's declared output exists whatever happens), NanoPlot is
# allowed to fail (the part after || runs only then), and a [FunFlux] line in the log says
# what happened. The job still counts as successful, and MultiQC simply has no NanoPlot
# row for that sample. The same guard protects QUAST and Qualimap (shared/20_qc.smk).
#
# v1 -> v2: the NanoPlot command is FunFluxL's plus --no_static (see above; FunFluxL also
# wrote the plots as PNG); the CPU count is now Snakemake's own threads: (capped at 8),
# which --cores enforces, where FunFluxL used a private "resources: cpus" that --cores did
# not count. The failure guard is new. The folder
# moved from 01.pre-processing/read_qc/{sample}/raw_read_qc to the sample's long-read
# folder, and the log is named after the rule.
rule raw_long_read_qc:
    input:
        reads = raw_long_fastq,
    output:
        qc_dir = directory(NANOPLOT_RAW_DIR),
    threads: capped_cpus(8)
    conda:
        "../../envs/nanoplot.yaml"
    log:
        LOGS + "/raw_long_read_qc_{sample}.log"
    shell:
        """
        mkdir -p {output.qc_dir}

        NanoPlot \
            --fastq {input.reads} \
            --threads {threads} \
            --loglength \
            --no_static \
            --prefix "{wildcards.sample}_" \
            --outdir {output.qc_dir} > {log} 2>&1 || {{
            echo "" >> {log}
            echo "[FunFlux] raw_long_read_qc: NanoPlot exited with an error for sample {wildcards.sample}. NanoPlot is only a report that no other step reads, so the run continues; the raw read QC will be missing from the MultiQC page for this sample. See the NanoPlot messages above." >> {log}
        }}
        """


# ───────────── filtered_long_read_qc — NanoPlot on the reads Filtlong kept ─────────────
# The same report as raw_long_read_qc, on the Filtlong output. Side by side, the two
# folders show what the filter removed: how the length distribution moved up, how much
# of the low-quality tail is gone, and how many bases are left for the assembler. If the
# filtered set looks too thin (too few bases for the expected genome size), the Filtlong
# settings in parameters.<mode>.filtlong (the block named after the mode, longread or
# hybrid) are the place to loosen.
#
# Takes in:    FILTERED_LONG, from rule filter_long_reads above.
# Does:        NanoPlot, exactly as above, --prefix "{sample}_" included.
# Produces:    NANOPLOT_FILTERED_DIR (01.reads/{sample}/long/filtered_qc/), a directory()
#              output.
# Consumed by: MultiQC (rule multiqc in shared/90_report.smk); rule all asks for it.
# Log:         NanoPlot's own messages.
#
# Non-fatal, like raw_long_read_qc: a NanoPlot failure leaves an empty folder and a
# [FunFlux] line in the log, and the run goes on.
#
# v1 -> v2: as for raw_long_read_qc. The folder moved from
# 01.pre-processing/read_qc/{sample}/filt_read_qc to the sample's long-read folder.
rule filtered_long_read_qc:
    input:
        reads = FILTERED_LONG,
    output:
        qc_dir = directory(NANOPLOT_FILTERED_DIR),
    threads: capped_cpus(8)
    conda:
        "../../envs/nanoplot.yaml"
    log:
        LOGS + "/filtered_long_read_qc_{sample}.log"
    shell:
        """
        mkdir -p {output.qc_dir}

        NanoPlot \
            --fastq {input.reads} \
            --threads {threads} \
            --loglength \
            --no_static \
            --prefix "{wildcards.sample}_" \
            --outdir {output.qc_dir} > {log} 2>&1 || {{
            echo "" >> {log}
            echo "[FunFlux] filtered_long_read_qc: NanoPlot exited with an error for sample {wildcards.sample}. NanoPlot is only a report that no other step reads, so the run continues; the filtered read QC will be missing from the MultiQC page for this sample. See the NanoPlot messages above." >> {log}
        }}
        """
