# FunFlux v2.1.0 — ONT adapter trimming (rules/shared/05_ont_adapters.smk): remove the Oxford
# Nanopore sequencing adapter, and the barcode that follows it, from the ends of ONT reads before
# Filtlong and the assembler see them.
#
# Why. An ONT read starts with the adapter that brought the DNA molecule to the pore, and in a
# barcoded library the barcode block follows it (native barcoding: AAGGTTAA, the 24 bp barcode,
# CAGCACCT). Dorado removes both only when it is run by hand (`dorado basecaller`, kit14 kits);
# reads basecalled live in MinKNOW, or with `dorado duplex`, usually keep them. The assemblers do
# not remove them either: in the hifiasm assembly of the public Botrytis cinerea B05.10 ONT reads
# (SRR40741277, about 70 % of them untrimmed), every contig end with a telomere carried about
# 80 bp of adapter and barcode beyond it. docs/methods/adapter-trimming.md gives the evidence.
#
# Three rules, once per ONT sample. The chain:
#
#   raw ONT FASTQ (the user's file, raw_long_fastq)
#     --find_ont_adapter-->  adapters/{sample}_adapter_scan.tsv      (ONT_ADAPTER_SCAN)
#                            adapters/{sample}_adapter_choice.tsv    (ONT_ADAPTER_CHOICE)
#                            adapters/{sample}_fastplong_options.txt (ONT_FASTPLONG_OPTIONS)
#                            adapters/{sample}_adapter_sequence.fasta (ONT_ADAPTER_SEQUENCE)
#                            (and reads ONT's published barcodes, ONT_BARCODE_LIST)
#     --trim_ont_adapters--> {sample}_trimmed.fastq.gz (TRIMMED_LONG, temporary)
#                            adapters/{sample}_fastplong.json, .html (FASTPLONG_JSON, FASTPLONG_HTML)
#       --> rule filter_long_reads (Filtlong; longread/10_reads.smk, hybrid/30_long_reads.smk)
#           reads the trimmed reads instead of the raw ones (filtlong_input in 00_common.smk)
#           and writes {sample}_filtered.fastq (FILTERED_LONG, kept)
#             --check_ont_trimming--> adapters/{sample}_adapter_check.tsv  (ONT_ADAPTER_CHECK)
#   The sequence FASTA also goes to rule blast_adapters (25_contig_summary.smk), so that the
#   contig adapter screen searches the assembly for exactly the sequences that were trimmed; the
#   choice and check files go to rule screen_adapters, which reports the numbers.
#
# What decides the sequence (parameters.<mode>.ont_adapter, read as ONT_ADAPTER_SETTING):
#   detect (the default) : find_ont_adapter.py works the sequence out from all the reads: the
#                          adapter, followed by the block that follows it in most reads (a
#                          barcode and its flanks), or nothing when the reads are already
#                          trimmed. It compares the reads with the barcodes ONT publishes
#                          (ONT_BARCODE_LIST), names the barcode it finds, and in a file that
#                          holds several barcodes trims each of them. The script's opening text
#                          explains the method.
#   a DNA sequence       : that sequence is trimmed as it is (the scan still runs and reports
#                          how many reads carry the adapter).
#   off                  : none of these rules exists, and Filtlong reads the raw reads, as in
#                          FunFlux v2.0.0.
# fastplong's own adapter guess is never used: it looks only at the first reads of the file and
# picks the most common sequence at the read ends, which in fungal reads can be the telomere
# repeat. On B05.10 it did, and cut the telomeres from the reads: the assembly kept telomeres at
# 2 contig ends instead of 34.
#
# PacBio HiFi reads are not trimmed. The PacBio software that makes HiFi reads moves reads with
# adapter sequence into a separate file of failed reads; the contig adapter screen still checks
# their assemblies.
#
# The Snakefile includes this file with the other shared files. The rules exist only in longread
# and hybrid mode AND only when at least one sample has ONT reads and trimming is not switched
# off: the `if` below tests HAS_LONG_READS and ONT_TRIM_SAMPLES (00_common.smk, section 13).
#
# Everything referenced here is defined once in 00_common.smk: HAS_LONG_READS, ONT_TRIM_SAMPLES,
# ONT_ADAPTER_SETTING, ONT_BARCODE_LIST, ONT_ADAPTER_SCAN, ONT_ADAPTER_CHOICE, ONT_FASTPLONG_OPTIONS,
# ONT_ADAPTER_SEQUENCE, TRIMMED_LONG, FILTERED_LONG, FASTPLONG_JSON, FASTPLONG_HTML, ONT_ADAPTER_CHECK,
# ONT_ADAPTER_SCRIPT, LONG_RUNTIME_MINUTES and LOGS, and the functions raw_long_fastq() (the
# sample's raw long-read file) and capped_cpus().
#
# conda: environment paths are written relative to THIS file (workflow/rules/shared/), so
# "../../envs/fastplong.yaml" goes up two folders, to workflow/, and then into workflow/envs/.
# The two script rules run in envs/python.yaml (a pinned Python 3.12): the script uses only
# Python's standard library.


# The rules below exist only when an ONT sample of this run is trimmed.
if HAS_LONG_READS and ONT_TRIM_SAMPLES:

    # The samples these rules may run for, written as a regular expression (the names joined by
    # |, "or"; re.escape() makes a dot in a name match only a dot). The rule-level wildcard
    # constraint narrows the workflow-wide one (00_common.smk, section 14), so the trimmed reads
    # of a PacBio sample can never be asked for, not even by mistake.
    _ONT_TRIM_SAMPLE_CONSTRAINT = "|".join(re.escape(sample) for sample in ONT_TRIM_SAMPLES)

    # ───────────── find_ont_adapter — work out what fastplong removes ─────────────
    # Takes in:    the sample's RAW ONT reads (raw_long_fastq; plain or .gz), and ONT_BARCODE_LIST,
    #              the native and rapid barcodes ONT publishes (shipped in workflow/data/adapters/).
    # Does:        python find_ont_adapter.py scan (ONT_ADAPTER_SCRIPT) with the setting
    #              (detect, or the given sequence). The script reads every read once; on the
    #              5.8 Gb B05.10 file that takes about two minutes. The setting reaches the script
    #              through params, so a changed ont_adapter re-runs the scan and, through its
    #              files, the trimming and everything after it. The barcode list goes through
    #              params as well, not as an input: a FunFlux update that adds barcodes to the
    #              list must not re-assemble every ONT sample of an existing output folder. The
    #              last line of the command copies the script's summary line, and its NOTE if
    #              there is one, from the log to the console.
    # Produces:    ONT_ADAPTER_SCAN       one row per 10,000 reads: how many start with the
    #                                     adapter (a batch of trimmed reads in an untrimmed file
    #                                     shows up here);
    #              ONT_ADAPTER_CHOICE     the decision and the numbers behind it;
    #              ONT_FASTPLONG_OPTIONS  the fastplong options, one line ("-s <the most common
    #                                     sequence>", with --disable_chimera_splitting when a
    #                                     sequence to trim is shorter than 50 bp, or "-A" when
    #                                     nothing is to be trimmed);
    #              ONT_ADAPTER_SEQUENCE   every sequence to trim, each in both orientations (FASTA),
    #                                     or an empty file; fastplong reads it with -a.
    #              All kept: they are small, and they document what was done to the reads.
    # Consumed by: trim_ont_adapters (the options and the sequences), check_ont_trimming (the
    #              choice), and blast_adapters (the sequences) and screen_adapters (the choice) in
    #              25_contig_summary.smk.
    rule find_ont_adapter:
        input:
            reads = raw_long_fastq,
        output:
            scan = ONT_ADAPTER_SCAN,
            choice = ONT_ADAPTER_CHOICE,
            options = ONT_FASTPLONG_OPTIONS,
            sequence = ONT_ADAPTER_SEQUENCE,
        params:
            script = ONT_ADAPTER_SCRIPT,
            setting = ONT_ADAPTER_SETTING,
            barcode_list = ONT_BARCODE_LIST,
        wildcard_constraints:
            sample = _ONT_TRIM_SAMPLE_CONSTRAINT,
        resources:
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/python.yaml"
        log:
            LOGS + "/find_ont_adapter_{sample}.log"
        shell:
            """
            set -euo pipefail

            python {params.script:q} scan \
                --reads {input.reads:q} \
                --setting {params.setting:q} \
                --sample {wildcards.sample:q} \
                --barcode-list {params.barcode_list:q} \
                --out-scan {output.scan:q} \
                --out-choice {output.choice:q} \
                --out-options {output.options:q} \
                --out-sequence {output.sequence:q} > {log} 2>&1

            grep -F '[FunFlux]' {log} >&2 || true
            """


    # ───────────── trim_ont_adapters — fastplong removes the sequence ─────────────
    # Takes in:    the sample's RAW ONT reads (raw_long_fastq), and ONT_FASTPLONG_OPTIONS and
    #              ONT_ADAPTER_SEQUENCE from find_ont_adapter.
    # Does:        fastplong with the options of the scan, read from the file with $(cat ...),
    #              the sequence file with -a, and these fixed options:
    #                -Q -L                  fastplong's quality and length filters off: Filtlong,
    #                                       the next step, is the one read filter of the workflow.
    #                --trimming_extension 0 cut exactly the sequence found. fastplong's default
    #                                       (10) also cuts the next 10 bases, which after a whole
    #                                       adapter + barcode are genome: on B05.10 it shortened
    #                                       every read that starts in a telomere.
    #                -z 1                   fast compression: the file is temporary.
    #              With -s, fastplong searches the first and last 200 bp of every read for the
    #              sequence (at the end, for its reverse complement), allowing a quarter of its
    #              bases to differ, and removes it with everything outside it. -a adds every
    #              sequence of the file, searched at both ends as written: this is how a file
    #              that holds several barcodes gets all of them trimmed (the file carries each
    #              sequence in both orientations, because -a does not add the reverse complement).
    #              With one sequence, -a repeats -s and changes nothing (checked on B05.10 reads).
    #              The file is given as an argument of its own, not through the options line, so
    #              that its path is quoted like every other path. With "-A" fastplong copies the
    #              reads unchanged, and the empty file is not read. The command first makes sure
    #              that the options file is not empty: without -s (or -A) fastplong would guess the
    #              adapter itself, which is exactly what must never happen.
    # Produces:    TRIMMED_LONG, {sample}_trimmed.fastq.gz, TEMPORARY: Snakemake deletes it once
    #              Filtlong has read it. Filtlong is its only reader, so that nothing but a change
    #              to the trimming itself can make Snakemake write it, and Filtlong and the
    #              assembly run, again. FASTPLONG_JSON and FASTPLONG_HTML,
    #              fastplong's reports (reads and bases in and out, reads trimmed, reads split),
    #              kept for people; MultiQC is not given them, because its fastp module would show
    #              them as a fastp run.
    # Consumed by: filter_long_reads (through filtlong_input in 00_common.smk).
    #
    # Threads: fastplong uses at most 16 worker threads (capped_cpus(16)); on B05.10 it took
    # about 4 minutes with 16.
    rule trim_ont_adapters:
        input:
            reads = raw_long_fastq,
            options = ONT_FASTPLONG_OPTIONS,
            sequence = ONT_ADAPTER_SEQUENCE,
        output:
            reads = temp(TRIMMED_LONG),
            json = FASTPLONG_JSON,
            html = FASTPLONG_HTML,
        wildcard_constraints:
            sample = _ONT_TRIM_SAMPLE_CONSTRAINT,
        threads: capped_cpus(16)
        resources:
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/fastplong.yaml"
        log:
            LOGS + "/trim_ont_adapters_{sample}.log"
        shell:
            """
            set -euo pipefail

            if [ ! -s {input.options:q} ]; then
                echo "[FunFlux] trim_ont_adapters {wildcards.sample}: the options file {input.options} is empty; fastplong is not started, because without them it would guess the adapter itself." | tee {log} >&2
                exit 1
            fi

            fastplong \
                -i {input.reads:q} \
                -o {output.reads:q} \
                -Q -L \
                --trimming_extension 0 \
                -z 1 \
                -w {threads} \
                -j {output.json:q} \
                -h {output.html:q} \
                -a {input.sequence:q} \
                $(cat {input.options:q}) > {log} 2>&1
            """


    # ───────────── check_ont_trimming — count the adapters left after Filtlong ─────────────
    # Takes in:    FILTERED_LONG, the reads Filtlong kept from the trimmed reads (rule
    #              filter_long_reads), and ONT_ADAPTER_CHOICE from find_ont_adapter (to know which
    #              sequence was trimmed). The kept Filtlong output is read, not the temporary
    #              trimmed reads: a re-run of this small check (after a FunFlux update, or when its
    #              file was deleted) then never makes Snakemake trim, filter and assemble again.
    #              These are also the reads the assembler gets.
    # Does:        python find_ont_adapter.py check: counts the reads that still start with the
    #              adapter and, when a barcode block was trimmed, with the end of the trimmed
    #              sequence. More than 1 % of reads still starting with the adapter prints a NOTE
    #              on the console; it does not stop the run.
    # Produces:    ONT_ADAPTER_CHECK, "name<TAB>value" lines, kept.
    # Consumed by: screen_adapters (25_contig_summary.smk), which puts the numbers into the
    #              Adapters table of the MultiQC report. Rule all asks for it too.
    rule check_ont_trimming:
        input:
            reads = FILTERED_LONG,
            choice = ONT_ADAPTER_CHOICE,
        output:
            check = ONT_ADAPTER_CHECK,
        params:
            script = ONT_ADAPTER_SCRIPT,
        wildcard_constraints:
            sample = _ONT_TRIM_SAMPLE_CONSTRAINT,
        resources:
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/python.yaml"
        log:
            LOGS + "/check_ont_trimming_{sample}.log"
        shell:
            """
            set -euo pipefail

            python {params.script:q} check \
                --reads {input.reads:q} \
                --choice {input.choice:q} \
                --sample {wildcards.sample:q} \
                --out {output.check:q} > {log} 2>&1

            grep -F '[FunFlux]' {log} >&2 || true
            """
