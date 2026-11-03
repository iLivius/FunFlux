# FunFlux v2.0.0 — hybrid front end, step 1: Illumina read preparation, PhiX removal and
# trimming (rules/hybrid/10_reads.smk).
#
# Identical to rules/illumina/10_reads.smk: keep in sync. Everything after this header,
# from the comment above _RAW_R1 to the end of the file, is a byte-for-byte copy of that
# file, comments and commands alike; only this header is the hybrid mode's own. Whoever
# edits one file edits the other; a diff of the two from that line on must be empty.
#
# Why the Illumina reads are prepared in hybrid mode, when the genome is built from the
# long reads. The Illumina pairs are not assembled into the delivered genome; they build
# the GUIDE assembly (SPAdes, rules/hybrid/20_assembly.smk), whose contamination screen
# tells which pairs belong to the isolate (rules/hybrid/25_guide.smk). Those selected
# pairs then score the ONT reads in Filtlong (rules/hybrid/30_long_reads.smk) and correct
# the finished long-read assembly with Polypolish (rules/hybrid/50_polish.smk). The PhiX
# removal and the trimming are the same as in illumina mode: a PhiX contig in the guide
# would draw Illumina pairs to it, and adapter bases would fragment the guide just as they
# would fragment a delivered assembly.
#
# The Snakefile includes every file in rules/hybrid/ in hybrid mode only, right after
# 00_common.smk, in name order: this file first, then 20_assembly.smk (the guide
# assembly), 25_guide.smk (the guide screen and the read selection), 30_long_reads.smk
# (Filtlong and NanoPlot), 40_long_assembly.smk (hifiasm or Flye) and 50_polish.smk
# (Medaka, Polypolish, finalize_assembly).
#
# What happens to a sample's Illumina reads here. {sample} is a wildcard: a placeholder
# that Snakemake fills in with each sample's name. The PhiX genome and its index are made
# ONCE per run, in 00.resources/phix/; everything per sample lies under
# 01.reads/{sample}/illumina/ in the output folder.
#
#   links.phix_link --download_phix--> phix.fna.gz                 (PHIX_FASTA, kept)
#                   --build_phix-----> phix.1.bt2 ... phix.rev.2.bt2 (PHIX_INDEX_PREFIX + the six endings, kept)
#                                           |
#   {sample}_R1.<ext>, {sample}_R2.<ext>    |          (the user's raw pairs in input.illumina_dir)
#     --map_phix (bowtie2 --local, --un-conc)--> {sample}_nophix_R1.fastq, _R2  (NOPHIX_R1/R2, temporary)
#     --trim_adapters (fastp)-----------------> {sample}_trim_R1.fastq, _R2    (TRIM_R1/R2, temporary)
#                                           +-> {sample}_fastp.html, .json      (FASTP_HTML/FASTP_JSON, kept)
#
# Where the trimmed pairs go in THIS mode. The rule comments below are the illumina
# mode's. Where they name the readers of the trimmed pairs as rule illumina_assembly in
# rules/illumina/20_assembly.smk and rule map_contigs in shared/10_decontam.smk, read: rule
# illumina_assembly in rules/hybrid/20_assembly.smk (the guide assembly), and rules
# map_guide_contigs (the coverage track of the guide screen) and map_sel_contigs (the read
# selection) in rules/hybrid/25_guide.smk. Snakemake deletes the trimmed pairs once all
# three have finished. In hybrid mode the shared rule map_contigs aligns the filtered LONG
# reads to the long-read assembly and never sees the Illumina pairs. The fastp JSON goes
# into the MultiQC page (rule multiqc in shared/90_report.smk) as the row "fastp read QC"
# and is one of the files rule all asks for, exactly as in illumina mode.
#
# Everything referenced here is defined once in 00_common.smk: ILLUMINA_DIR, SHORT_EXTN,
# PHIXURL, PHIX_FASTA, PHIX_INDEX_PREFIX, READS_ILLUMINA_DIR, NOPHIX_R1, NOPHIX_R2, TRIM_R1,
# TRIM_R2, FASTP_JSON, FASTP_HTML, LOGS and CPUS. So is the function capped_cpus().
#
# conda: Snakemake resolves a rule's environment path relative to the rule file, not to
# the working directory, so the files in workflow/envs/ are written "../../envs/<tool>.yaml".
# bowtie2 (two rules) and fastp get an environment, and download_phix the small download
# environment (envs/download.yaml, see the rule).
#
# CPUs: a rule that uses several CPUs declares Snakemake's own `threads:` and passes
# {threads} to its tool. That is what makes `--cores N` actually limit the run.

# The user's own read pair for one sample: <input.illumina_dir>/<sample>_R1.<extension>
# and the matching _R2 file. The extension is the one 00_common.smk found and accepted for
# this batch (fastq, fq, fastq.gz or fq.gz; every sample of a batch must use the same
# one). No rule makes these files; sample discovery found them by listing the folder, so
# an illumina run enters the workflow here. Only rule map_phix uses the two patterns, so
# they are not kept in 00_common.smk.
_RAW_R1 = os.path.join(ILLUMINA_DIR, "{sample}_R1." + SHORT_EXTN)
_RAW_R2 = os.path.join(ILLUMINA_DIR, "{sample}_R2." + SHORT_EXTN)

# The file-name pattern bowtie2 is given for the pairs that are NOT PhiX (--un-conc, see
# rule map_phix). bowtie2 replaces the % with 1 for the first mates and 2 for the second,
# so this pattern MUST give exactly NOPHIX_R1 and NOPHIX_R2 (00_common.smk), the two files
# the rule declares as its outputs; otherwise bowtie2 writes files Snakemake does not
# expect and the rule fails for "missing output files". Only map_phix uses it.
_NOPHIX_PATTERN = READS_ILLUMINA_DIR + "/{sample}_nophix_R%.fastq"


# ───────────── download_phix — fetch the PhiX174 genome ─────────────
# Once per run, not per sample: the rule has no {sample} in its outputs, so Snakemake runs
# it a single time and every map_phix job waits for it.
#
# Takes in:    nothing from other rules; the download link links.phix_link from the
#              configuration (PHIXURL; the template points at NCBI's assembly of the
#              5.4 kb PhiX174 genome as a gzip-compressed FASTA of about 2 kB).
# Does:        downloads the file with curl. --fail turns an HTTP error (a moved or missing
#              file) into a failed command instead of a saved error page, and a failed
#              command stops the rule (Snakemake runs every shell block with
#              `set -euo pipefail`). curl's own messages go to the log.
# Produces:    PHIX_FASTA, 00.resources/phix/phix.fna.gz, still compressed: bowtie2-build
#              reads gzip directly, so nothing unpacks it. Kept, as all run-level material
#              in 00.resources/ is, so a resumed run or a second batch into the same output
#              folder does not download it again.
# Consumed by: build_phix.
#
# Environment: envs/download.yaml, which supplies curl (a minimal Linux installation may have
# none), as for the other download rules (funannotate_db in shared/40_annotation.smk,
# tiberius_model in shared/42_tiberius.smk, ...).
#
# Where it runs: always on the computer that runs Snakemake itself (localrule: True). On a
# single computer that changes nothing; on a cluster it keeps the download off the compute
# nodes, which often have no internet.
rule download_phix:
    localrule: True
    output:
        phix = PHIX_FASTA,
    params:
        link = PHIXURL,
    conda:
        "../../envs/download.yaml"
    log:
        LOGS + "/download_phix.log"
    shell:
        """
        curl -L --fail -o "{output.phix}" "{params.link}" > {log} 2>&1
        """


# ───────────── build_phix — the bowtie2 index of the PhiX genome ─────────────
# Once per run, like download_phix. bowtie2 cannot align reads to a plain FASTA; it first
# needs its own index of the reference. The index is six files that bowtie2 names itself
# from the prefix it is given (<prefix>.1.bt2 to <prefix>.4.bt2 and <prefix>.rev.1.bt2,
# <prefix>.rev.2.bt2), so the rule declares the six names and hands the prefix to
# bowtie2-build.
#
# Takes in:    PHIX_FASTA from download_phix (bowtie2-build reads the gzip file as it is).
# Does:        bowtie2-build writes the index under PHIX_INDEX_PREFIX, next to the genome.
# Produces:    the six .bt2 files in 00.resources/phix/, kept (run-level material; they are
#              a few kilobytes).
# Consumed by: map_phix.
#
# Threads: one. The PhiX genome is 5.4 kb; indexing it takes a second.
rule build_phix:
    input:
        phix = PHIX_FASTA,
    output:
        idx1 = PHIX_INDEX_PREFIX + ".1.bt2",
        idx2 = PHIX_INDEX_PREFIX + ".2.bt2",
        idx3 = PHIX_INDEX_PREFIX + ".3.bt2",
        idx4 = PHIX_INDEX_PREFIX + ".4.bt2",
        ridx1 = PHIX_INDEX_PREFIX + ".rev.1.bt2",
        ridx2 = PHIX_INDEX_PREFIX + ".rev.2.bt2",
    params:
        basename = PHIX_INDEX_PREFIX,
    conda:
        "../../envs/bowtie.yaml"
    log:
        LOGS + "/build_phix.log"
    shell:
        """
        bowtie2-build {input.phix} {params.basename} > {log} 2>&1
        """


# ───────────── map_phix — remove the PhiX read pairs ─────────────
# Every read pair of the sample is aligned to the PhiX genome, and the pairs that do NOT
# align are written out: they are the library without the sequencing control. A pair that
# aligns concordantly is PhiX and is simply not written anywhere.
#
# Takes in:    the six index files from build_phix, and the sample's raw read pair
#              (_RAW_R1, _RAW_R2, the files in input.illumina_dir; bowtie2 reads gzip
#              FASTQ directly). bowtie2 is given the index PREFIX (params db), not the
#              files; the files are listed as inputs so that Snakemake builds the index
#              first.
# Does:        bowtie2 --local aligns the pairs. Local alignment may leave the ends of a
#              read unaligned (soft-clipped), so a PhiX read that still carries adapter
#              sequence at its end, or a low-quality tail, is recognised as PhiX all the
#              same; in the default end-to-end mode the adapter bases would count as
#              mismatches and could let a phage read through. --un-conc writes the pairs
#              that did not align concordantly (both mates aligned, facing each other at
#              a plausible distance) to two files, one per mate, named from the pattern
#              _NOPHIX_PATTERN. The alignments themselves are not kept: -S /dev/null
#              throws the SAM text away, because nothing reads it and for a deep library
#              it would be tens of gigabytes.
# Produces:    NOPHIX_R1 and NOPHIX_R2, the PhiX-free pairs, uncompressed, temporary:
#              Snakemake deletes them as soon as trim_adapters has read them. (Each is
#              about as large as the raw file unpacked.)
# Consumed by: trim_adapters.
# Log:         bowtie2's alignment summary. "overall alignment rate" there is the share of
#              the library that was PhiX.
#
# Threads: all CPUS (bowtie2 --threads). The whole library is aligned, which is the
# longest of the four steps of this file.
rule map_phix:
    input:
        idx1 = PHIX_INDEX_PREFIX + ".1.bt2",
        idx2 = PHIX_INDEX_PREFIX + ".2.bt2",
        idx3 = PHIX_INDEX_PREFIX + ".3.bt2",
        idx4 = PHIX_INDEX_PREFIX + ".4.bt2",
        ridx1 = PHIX_INDEX_PREFIX + ".rev.1.bt2",
        ridx2 = PHIX_INDEX_PREFIX + ".rev.2.bt2",
        r1 = _RAW_R1,
        r2 = _RAW_R2,
    output:
        r1 = temp(NOPHIX_R1),
        r2 = temp(NOPHIX_R2),
    params:
        db = PHIX_INDEX_PREFIX,
        unaligned = _NOPHIX_PATTERN,
    threads: CPUS
    # Asked of a cluster scheduler: a time limit only, as the step reads every read but
    # holds little in memory (00_common.smk, section 4).
    resources:
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/bowtie.yaml"
    log:
        LOGS + "/map_phix_{sample}.log"
    shell:
        """
        bowtie2 \
            -x {params.db} \
            -1 {input.r1} \
            -2 {input.r2} \
            --threads {threads} \
            --un-conc {params.unaligned} \
            -S /dev/null \
            --local > {log} 2>&1
        """


# ───────────── trim_adapters — fastp: adapters off, bad ends off, read QC report ─────────────
# The last step before assembly. Illumina reads can end in adapter sequence (when the DNA
# fragment was shorter than the read), and base quality drops towards the 3' end of a
# read; fastp removes both. It also writes the read-quality report of the sample: read
# counts, lengths, quality and GC content before and after, duplication rate, insert-size
# estimate, adapter content.
#
# Takes in:    NOPHIX_R1 and NOPHIX_R2 from map_phix.
# Does:        fastp on the pair, with these settings (all fixed; there is no config key):
#                --detect_adapter_for_pe  also work out the adapter sequence from the reads
#                                         themselves for paired-end data. For pairs, fastp
#                                         otherwise finds adapters only where the two mates
#                                         overlap; this adds the per-read detection it
#                                         always does for single-end reads.
#                --cut_front              slide a window from the 5' end of each read and
#                                         cut while the window's mean quality is low.
#                --cut_right              slide a window from 5' to 3' and, at the first
#                                         window of low mean quality, cut it and everything
#                                         after it.
#                                         Window size and quality threshold are fastp's
#                                         defaults: 4 bases, mean quality 20.
#                --length_required 75     a read left shorter than 75 bases is dropped,
#                                         with its mate: SPAdes' larger k-mers (up to 127)
#                                         need reads well longer than that.
#                --verbose                progress messages, into the log.
# Produces:    TRIM_R1 and TRIM_R2, the trimmed pairs, temporary: Snakemake deletes them
#              once both readers below have finished. FASTP_HTML and FASTP_JSON, the report
#              as a page for people and as numbers for MultiQC; both are kept in
#              01.reads/{sample}/illumina/ and the JSON is one of the files rule all asks
#              for.
# Consumed by: the trimmed pairs by rule illumina_assembly (rules/illumina/20_assembly.smk)
#              and rule map_contigs (shared/10_decontam.smk); the JSON by rule multiqc
#              (shared/90_report.smk).
# Log:         fastp's summary (reads in, reads out, how many failed each filter), which
#              fastp writes to the error stream.
#
# Threads: up to 16 (capped_cpus). fastp gains little beyond that, and the cap leaves CPUs
# for other samples' jobs.
rule trim_adapters:
    input:
        r1 = NOPHIX_R1,
        r2 = NOPHIX_R2,
    output:
        r1 = temp(TRIM_R1),
        r2 = temp(TRIM_R2),
        html = FASTP_HTML,
        json = FASTP_JSON,
    threads: capped_cpus(16)
    conda:
        "../../envs/fastp.yaml"
    log:
        LOGS + "/trim_adapters_{sample}.log"
    shell:
        """
        fastp \
            --detect_adapter_for_pe \
            --length_required 75 \
            --cut_front \
            --cut_right \
            --thread {threads} \
            --verbose \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -j {output.json} \
            -h {output.html} > {log} 2>&1
        """
