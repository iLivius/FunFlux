# FunFlux v2.1.0 — hybrid front end, step 4: long-read filtering, guided by the Illumina
# pairs when parameters.hybrid.short_read_guidance is on, and long-read QC
# (rules/hybrid/30_long_reads.smk).
#
# Three rules, the long-read front end of this mode:
#   filter_long_reads     : Filtlong keeps the best long reads. For an ONT sample the
#                           sample's selected Illumina pairs are the yardstick (this file's
#                           own version of the rule; see below).
#   raw_long_read_qc      : NanoPlot on the reads as delivered.
#   filtered_long_read_qc : NanoPlot on the reads Filtlong kept.
# The two NanoPlot rules are copies (Identical to rules/longread/10_reads.smk: keep in
# sync), comments and commands alike; a change there must be made here too, and the other
# way round. The Filtlong rule is NOT a copy: it is the longread mode's command with the
# guidance options in front, and it is explained in full here.
#
# The Snakefile includes every file in rules/hybrid/ in hybrid mode only, in name order;
# this file comes after 25_guide.smk (which makes the selected pairs this file needs) and
# before 40_long_assembly.smk.
#
# What happens to a sample's long reads. {sample} is a wildcard: a placeholder that
# Snakemake fills in with each sample's name, so every rule below runs once per sample.
# The raw file is the user's own ({sample}_ont.<ext> or {sample}_pacbio.<ext> in
# input.longread_dir; the suffix decides the technology); everything made from it lies
# under 01.reads/{sample}/long/ in the output folder.
#
#   {sample}_ont.fastq.gz
#     --find_ont_adapter, --trim_ont_adapters--> {sample}_trimmed.fastq.gz (TRIMMED_LONG,
#       temporary; rules/shared/05_ont_adapters.smk): the ONT adapter and the barcode after it
#       removed from the read ends (unless parameters.hybrid.ont_adapter is off)
#   {sample}_trimmed.fastq.gz  + {sample}_selected_R1.fastq, _R2  (SEL_R1/R2, rule map_sel_contigs
#                                                                 in rules/hybrid/25_guide.smk)
#     --filter_long_reads (Filtlong -1 -2 [--trim --split N] ...)--> {sample}_filtered.fastq  (FILTERED_LONG, kept)
#   {sample}_pacbio.fastq.gz  (PacBio HiFi: never guided)
#     --filter_long_reads (Filtlong, the longread mode's command)--> {sample}_filtered.fastq  (FILTERED_LONG, kept)
#       --> rule long_read_assembly (40_long_assembly.smk): the reads the genome is built from
#       --> rule map_contigs (shared/10_decontam.smk): the reads aligned back to the draft,
#           the coverage track of the contamination screen and the BAM Qualimap evaluates
#       --> rule long_read_consensus (50_polish.smk): the reads Medaka polishes with (ONT)
#     --raw_long_read_qc------> raw_qc/       (NANOPLOT_RAW_DIR)      --> MultiQC
#     --filtered_long_read_qc-> filtered_qc/  (NANOPLOT_FILTERED_DIR) --> MultiQC
#
# Why the Illumina pairs are brought in. Filtlong ranks long reads by length and quality
# and keeps the best share of the bases. Without a reference, "quality" is the read's own
# Phred scores, which for ONT reads are a rough guide. With the Illumina pairs as a
# reference (-1 and -2), Filtlong does not read the Phred scores at all: it counts how many
# of each read's 16-mers (stretches of 16 bases) also occur in the Illumina reads, which
# are far more accurate. A stretch of a long read that no Illumina read supports is either
# a sequencing error, a chimeric join, or a region the Illumina library did not cover; it
# scores low. --trim removes such bases at the two ends of a read and --split N cuts the
# read wherever N or more consecutive bases have no support, keeping the pieces as
# separate reads. Only pairs the guide screen kept are used (SEL_R1/R2), so a contaminant
# cannot vouch for its own ONT reads. docs/methods/filtlong-guidance.md explains the
# settings. The guidance is OFF by default: the Illumina reads cover repeats worst, and the
# long reads that span repeats hold an assembly together, so scoring them low or cutting
# them there can break the assembly apart (BacFlux measured this on five bacterial
# genomes). With it on, the split length, 1000, is FunFluxL+'s.
#
# The comments of the two NanoPlot copies below name the Filtlong settings by the block
# of the active mode (parameters.<mode>.filtlong: here parameters.hybrid.filtlong) and
# rule filter_long_reads as "above", which is true in both files.
#
# Everything referenced here is defined once in 00_common.smk: FILTERED_LONG,
# NANOPLOT_RAW_DIR, NANOPLOT_FILTERED_DIR, FILTLONG_MIN_LENGTH, FILTLONG_KEEP_PERCENT,
# FILTLONG_LENGTH_WEIGHT, FILTLONG_TARGET_ARG, LOGS and CPUS, and the functions
# raw_long_fastq() (the sample's raw long-read file), filtlong_input() (the file Filtlong
# starts from: the trimmed reads of an ONT sample, the raw file otherwise),
# filtlong_guidance() (the sample's guidance options, or nothing), filtlong_guide_reads()
# (the selected pairs those options name, or nothing) and capped_cpus().
#
# conda: environment paths are relative to THIS file (workflow/rules/hybrid/), so
# "../../envs/filtlong.yaml" is workflow/envs/filtlong.yaml.
#
# CPUs: a rule that uses several CPUs declares Snakemake's own `threads:` and passes
# {threads} to its tool. That is what makes `--cores N` actually limit the run. Filtlong
# runs on one CPU and declares none.
#
# In the shell commands below, a doubled brace ({{ or }}) stands for a single brace in the
# shell. Snakemake reads single braces as its own placeholders, such as {input.reads}, and
# turns a doubled brace into a single one.


# ───────────── filter_long_reads — Filtlong, with the Illumina pairs as the yardstick ─────────────
# Takes in:    filtlong_input: for an ONT sample, the reads fastplong trimmed (TRIMMED_LONG,
#              rule trim_ont_adapters in rules/shared/05_ont_adapters.smk); for a PacBio
#              sample, or with parameters.hybrid.ont_adapter off, the sample's RAW long-read
#              FASTQ, found by sample discovery in input.longread_dir. Plain or
#              gzip-compressed; Filtlong reads both. And, for a guided ONT sample, the sample's own
#              selected Illumina pairs SEL_R1 and SEL_R2 (rule map_sel_contigs in
#              rules/hybrid/25_guide.smk), through filtlong_guide_reads() in 00_common.smk:
#              the files the -1 and -2 below name, listed as inputs so that Filtlong waits
#              for them. A PacBio sample, or any sample with guidance off, gets an empty
#              list, which Snakemake accepts as "nothing to wait for": its filtering and
#              assembly start without waiting for the guide screen they do not use.
# Does:        filtlong with the same four settings as in longread mode, here from
#              parameters.hybrid.filtlong:
#                --min_length    reads shorter than this (bp) are dropped outright
#                                (default 1000).
#                --keep_percent  of the remaining BASES, keep the best this-many percent
#                                (default 90).
#                --length_weight how much a read's length counts, relative to its quality,
#                                in the ranking (default 10, ten times Filtlong's own
#                                default: for assembly, a long read is worth more than a
#                                slightly cleaner short one).
#                --target_bases  an absolute cap on the bases kept (params target_bases,
#                                "--target_bases <size>" or nothing; 00_common.smk builds
#                                the text as FILTLONG_TARGET_ARG).
#              and, FIRST on the command line, the sample's guidance options (params
#              guidance, from filtlong_guidance() in 00_common.smk), which are one of:
#                ONT, guidance on, filtlong_split N (default 1000):
#                    -1 <SEL_R1> -2 <SEL_R2> --trim --split N
#                ONT, guidance on, filtlong_split null:
#                    -1 <SEL_R1> -2 <SEL_R2>      scored against the pairs, never cut.
#                    --trim goes with --split: with filtlong_split null the reads are
#                    neither split nor trimmed, so that a read is either cut where the
#                    Illumina reads stop supporting it or left whole, not something in
#                    between.
#                PacBio HiFi, or guidance off:
#                    nothing: the command is the longread mode's.
#              Both texts are written UNQUOTED on purpose: an empty value then leaves
#              nothing behind (no bare flag, no empty argument), and a filled one expands
#              into its several options.
# Produces:    FILTERED_LONG, {sample}_filtered.fastq, uncompressed and KEPT: the assembler,
#              the read mapping and Medaka all read it, and a re-run of any of them must
#              not have to filter the reads again.
# Consumed by: rule long_read_assembly (40_long_assembly.smk), rule map_contigs
#              (shared/10_decontam.smk), rule long_read_consensus (50_polish.smk) and rule
#              filtered_long_read_qc below. Rule all asks for it too.
# Log:         Filtlong's own report: how many reads and bases came in, how many were kept,
#              and with -1/-2 how many Illumina 16-mers it loaded.
#
# PacBio HiFi reads are never guided, whatever the configuration says: they are already
# accurate (each is the consensus of many passes over one molecule), and cutting them
# where the Illumina reads leave a gap would throw away the long, accurate reads HiFi is
# used for. The dry run prints, per sample, whether Filtlong is guided.
#
# v1 -> v2: FunFluxL+ passed -1 -2 --trim --split to EVERY sample, PacBio included, and had
# no switch for it. Guidance is now a setting (parameters.hybrid.short_read_guidance,
# parameters.hybrid.filtlong_split) and applies to ONT samples only; the options arrive as
# one text prepared per sample instead of being typed into the command, so --trim --split
# now stand before --min_length instead of after it (Filtlong does not care about the
# order), and --target_bases reaches the command as in the longread rule (FunFluxL+ built
# it in the shell as a bash array). For a guided ONT sample the options and their values
# are FunFluxL+'s. The output moved from 01.read_processing/long/{sample}/filtered.fastq to
# the sample's long-read folder, and the log is named after the rule (FunFluxL+:
# filter_reads_{sample}.log).
rule filter_long_reads:
    input:
        reads = filtlong_input,
        guide_reads = filtlong_guide_reads,
    output:
        reads = FILTERED_LONG,
    params:
        guidance = filtlong_guidance,
        min_length = FILTLONG_MIN_LENGTH,
        keep_percent = FILTLONG_KEEP_PERCENT,
        length_weight = FILTLONG_LENGTH_WEIGHT,
        target_bases = FILTLONG_TARGET_ARG,
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/filtlong.yaml"
    log:
        LOGS + "/filter_long_reads_{sample}.log"
    shell:
        """
        filtlong \
            {params.guidance} \
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
