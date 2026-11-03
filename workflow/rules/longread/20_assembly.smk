# FunFlux v2.0.0 — longread front end, step 2: the assembler
# (rules/longread/20_assembly.smk).
#
# Two rules. long_read_assembly builds the draft genome from the filtered long reads. It
# comes in two versions, hifiasm (the default) and Flye, and only ONE of them exists in a
# run: the `if LONG_READ_ASSEMBLER == "hifiasm": ... else: ...` below is decided once, when
# Snakemake reads the workflow, from parameters.longread.assembler (00_common.smk, section
# 11). One assembler serves the whole batch; the technology of each sample (ONT or PacBio
# HiFi, from the file-name suffix) only changes the assembler's read-type flag.
# assembly_evidence then turns the assembler's own per-contig file into the evidence table.
# The Snakefile includes this file in longread mode only, after 10_reads.smk.
#
# Whichever assembler runs, the two rules end with the SAME three files directly in the
# sample's assembly folder, so that nothing after this point has to know which assembler
# was used:
#
#   {sample}_filtered.fastq  (FILTERED_LONG, from rule filter_long_reads)
#     --long_read_assembly--> 02.assembly/{sample}/
#          hifiasm/ or flye/           the assembler's own files (HIFIASM_DIR / FLYE_DIR)
#          assembly.fasta              DRAFT_CONTIGS: the draft, one record per contig
#            --> shared/10_decontam.smk: coverage track, BLAST, BlobTools, taxonomy selector
#            --> shared/15_mitochondria.smk: the mitochondrion screen
#          assembly_provenance.tsv     ASSEMBLY_PROVENANCE: tool, version, command
#            --> nobody; a delivered record, rule all asks for it
#   the assembler's own per-contig file (in hifiasm/ or flye/)
#     --assembly_evidence--> 02.assembly/{sample}/
#          assembly_evidence.tsv       ASSEMBLY_EVIDENCE: per-contig depth and circularity
#            --> rule screen_mitochondria (shared/15_mitochondria.smk) and rule
#                contig_summary (shared/25_contig_summary.smk)
#
# Why the two tables. The mitochondrial genome of a fungus is a small circle present in
# many copies per cell, so its contig is many times deeper than the nuclear contigs and
# the assembler often reports it as circular. The mitochondrion screen weighs that next to
# the BLAST hits. Both assemblers know depth and circularity, but write them in their own
# formats (Flye in assembly_info.txt, hifiasm as tags in its graph), so rule
# assembly_evidence converts them into one common table with mitochondria.py
# (evidence-flye / evidence-hifiasm; the script explains the columns). The provenance table
# records which tool and version made the draft and with what command, because the choice
# of assembler is a per-run setting that would otherwise be visible only in the config.
#
# Before the assembler starts, an ONT sample that will be polished with Medaka waits for
# rule check_medaka_model (shared/12_medaka_check.smk): medaka_gate() returns that rule's
# output file for such a sample and nothing for any other, and the rule lists it as an
# input. So a wrong Medaka model name stops the run here, in seconds, instead of after the
# assembly has run for hours. The file is not read by the assembler; it is only waited for.
#
# The hybrid mode runs the same assembler: rules/hybrid/40_long_assembly.smk is a copy of
# this file.
#
# Everything referenced here is defined once in 00_common.smk: LONG_READ_ASSEMBLER,
# HIFIASM_COMMON_ARGS, FILTERED_LONG, HIFIASM_DIR, FLYE_DIR, DRAFT_CONTIGS,
# ASSEMBLY_EVIDENCE, ASSEMBLY_PROVENANCE, MITOCHONDRIA_SCRIPT, LOGS and CPUS, and the
# per-sample functions hifiasm_tech_flag(), flye_read_flag(), read_tech() and
# medaka_gate(). Two file names only this file uses are defined just below, and the three
# names that set up the evidence step just above rule assembly_evidence.
#
# conda: environment paths are relative to THIS file (workflow/rules/longread/), so
# "../../envs/hifiasm.yaml" is workflow/envs/hifiasm.yaml. The assembler runs in its own
# environment (hifiasm.yaml or flye.yaml) and starts no Python. The evidence table is made
# by a Python script that uses only the standard library, in rule assembly_evidence, which
# runs in envs/python.yaml, a pinned Python 3.12 and nothing else; the comment above that
# rule says why it is a rule of its own.
#
# CPUs: the assembler is the most CPU-hungry step of the run and gets the whole budget,
# `threads: CPUS` (the config's resources.threads), which --cores enforces. The evidence
# step needs one.
#
# In the shell commands below, a doubled brace ({{ or }}) stands for a single brace in the
# shell. Snakemake reads single braces as its own placeholders, such as {input.reads}, and
# turns a doubled brace into a single one. The {{ ... }} > file blocks group several
# printf lines so that their output goes to one file.


# The assembler's own files that the two rules of this file read. hifiasm writes its
# assembly as graphs (GFA files), not FASTA; with the --primary option it writes one primary
# and one alternate assembly, named <prefix>.p_ctg.gfa and <prefix>.a_ctg.gfa, and the prefix
# the rule passes (-o) is HIFIASM_DIR/hifiasm: long_read_assembly reads the contigs from the
# primary graph, and assembly_evidence reads their depth and circularity from it. Flye writes
# assembly.fasta and, next to it, assembly_info.txt, its per-contig table (length, depth,
# circular yes/no, repeat yes/no, alternative-haplotype group), which assembly_evidence
# reads. Only this file reads these two paths, so they are kept here rather than in
# 00_common.smk, as contigs/20_assembly.smk does with _INPUT_FASTA.
_HIFIASM_PRIMARY_GFA = HIFIASM_DIR + "/hifiasm.p_ctg.gfa"
_FLYE_INFO           = FLYE_DIR + "/assembly_info.txt"


if LONG_READ_ASSEMBLER == "hifiasm":

    # ───────────── long_read_assembly (hifiasm) ─────────────
    # hifiasm builds a string graph from read overlaps and, by default, tries to separate
    # the two haplotypes of a diploid genome. Most fungal isolates are haploid, so the
    # options 00_common.smk prepares (HIFIASM_COMMON_ARGS) start with --primary, which
    # asks for one primary assembly plus one alternate set instead of two phased
    # haplotypes, and -f 0, which switches off hifiasm's k-mer Bloom filter (a memory table
    # sized for large genomes; 16 GB at the default setting, and hifiasm's manual advises
    # 0 for small genomes). The optional settings of parameters.<mode>.hifiasm (genome
    # size, duplicate purging, homozygous coverage, telomere motif; <mode> is the block
    # named after the mode, longread or hybrid) follow in the same text. ONT reads get
    # --ont on top (params tech_flag, from hifiasm_tech_flag()): hifiasm then expects
    # noisier, simplex reads instead of HiFi reads, and corrects them first. PacBio HiFi
    # samples get nothing, which is hifiasm's native HiFi mode.
    #
    # Takes in:    FILTERED_LONG (rule filter_long_reads), and, for an ONT sample that
    #              uses Medaka, the model-check file (medaka_gate; see the head of this
    #              file). The per-run options and the technology flag are params, so a
    #              changed hifiasm setting re-runs the assembly.
    # Does:        1. hifiasm into HIFIASM_DIR with the prefix "hifiasm". -i tells hifiasm
    #                 to ignore the error-corrected reads and overlaps it saved in an
    #                 earlier run (hifiasm.ec.bin, hifiasm.ovlp.*.bin): a Snakemake re-run
    #                 after a changed input or setting must never reuse that stale cache.
    #                 The price is that hifiasm cannot resume its own work; keep -i.
    #                 Everything hifiasm prints goes to the log.
    #              2. awk turns the primary graph into FASTA: every "S" line of a GFA file
    #                 is one segment (here one contig) with its name in field 2 and its
    #                 sequence in field 3, so the program writes ">name" and the sequence,
    #                 one line each, for every S line whose sequence is not "*" (a graph
    #                 can carry segments without sequence), and counts them. END ... exit 1
    #                 fails the rule when no contig was written: an empty assembly would
    #                 otherwise travel on and fail much later. hifiasm contig names end in
    #                 l (linear) or c (circular), for example ptg000001l; rule
    #                 assembly_evidence (below) reads that letter.
    #              3. the provenance table, written by the rule itself with printf: one
    #                 "field<TAB>value" line each for the assembler, its version (hifiasm
    #                 --version), the sample's technology (params read_type: ont or
    #                 pacbio), the input file and the whole command line.
    # Produces:    HIFIASM_DIR (directory() output: a re-run starts from an empty folder,
    #              which is also why -i has nothing stale to find), the primary graph
    #              inside it (named as an output of its own because rule assembly_evidence
    #              reads it), and two of the three standardised files: DRAFT_CONTIGS and
    #              ASSEMBLY_PROVENANCE. The third, ASSEMBLY_EVIDENCE, is made by rule
    #              assembly_evidence.
    # Consumed by: see the head of this file.
    # Log:         hifiasm's whole output (progress, k-mer histogram, coverage estimate).
    #              The awk and printf lines after it are not redirected: they print nothing
    #              unless they fail, and then their error appears in Snakemake's own output,
    #              not in this log.
    #
    # In the provenance line the technology flag is written inside plain double quotes
    # ("{params.tech_flag}"), not with the usual {...:q} quoting: Snakemake writes an
    # EMPTY value with :q as nothing at all, and a missing argument would shift every
    # later %s of the printf by one place for a PacBio sample. Inside the quotes an empty
    # flag is an empty argument, and the line simply shows a double space there.
    #
    # Notes on the command. (1) The options arrive in two pieces, the per-run text and the
    # per-sample flag, so --ont stands last; hifiasm does not care about the order, and the
    # provenance line shows it. (2) The CPU count is Snakemake's own threads:, which --cores
    # enforces. (3) No assembly_linearized.fasta copy is written: every reader of the draft
    # (BLAST, BlobTools, the mitochondrion screen, the selector) reads multi-line FASTA, and
    # the awk above writes one line per record anyway. (4) The log is named after the rule,
    # and the explanation of -i sits in this comment, not inside the command, so that
    # rewording it never re-runs the assembly. (5) The evidence table is written by rule
    # assembly_evidence (below).
    rule long_read_assembly:
        input:
            reads = FILTERED_LONG,
            model_check = medaka_gate,
        output:
            native_dir = directory(HIFIASM_DIR),
            primary_gfa = _HIFIASM_PRIMARY_GFA,
            contigs = DRAFT_CONTIGS,
            provenance = ASSEMBLY_PROVENANCE,
        params:
            common_args = HIFIASM_COMMON_ARGS,
            tech_flag = hifiasm_tech_flag,
            read_type = read_tech,
        threads: CPUS
        # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
        resources:
            mem_mb = HEAVY_MEM_MB,
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/hifiasm.yaml"
        log:
            LOGS + "/long_read_assembly_{sample}.log"
        shell:
            """
            mkdir -p {output.native_dir}
            hifiasm \
                -i \
                {params.common_args} {params.tech_flag} \
                -t {threads} \
                -o {output.native_dir}/hifiasm \
                {input.reads} > {log} 2>&1

            awk '$1 == "S" && $3 != "*" {{print ">" $2; print $3; count++}} END {{if (count == 0) exit 1}}' \
                {output.primary_gfa} > {output.contigs}
            {{
                printf 'field\tvalue\n'
                printf 'assembler\thifiasm\n'
                printf 'version\t'
                hifiasm --version
                printf 'read_type\t%s\n' {params.read_type:q}
                printf 'input_fastq\t%s\n' {input.reads:q}
                printf 'command\thifiasm -i %s %s -t %s -o %s/hifiasm %s\n' \
                    {params.common_args:q} "{params.tech_flag}" {threads} {output.native_dir:q} {input.reads:q}
            }} > {output.provenance}
            """

else:

    # ───────────── long_read_assembly (Flye) ─────────────
    # Flye assembles from approximate read overlaps (a repeat graph), then polishes its
    # own contigs with the reads for a number of rounds. The read-type flag (params
    # read_flag, from flye_read_flag()) tells Flye what error profile to expect:
    # --pacbio-hifi for every PacBio HiFi sample; for ONT samples --nano-hq (reads from
    # recent chemistry; the default) or --nano-raw (older, noisier reads), chosen in
    # 00_common.smk, section 11, from parameters.<mode>.flye_input_mode and the Medaka
    # model. --iterations 5 asks for five rounds of Flye's own polishing instead of its
    # default one; a FunFluxL setting, kept. Flye's result is not byte-identical from run
    # to run, so two runs of the same reads may differ in small ways.
    #
    # Takes in:    FILTERED_LONG (rule filter_long_reads), and, for an ONT sample that
    #              uses Medaka, the model-check file (medaka_gate; see the head of this
    #              file).
    # Does:        1. flye into FLYE_DIR (--out-dir). Flye writes assembly.fasta,
    #                 assembly_info.txt, its graph and its working files there.
    #                 Everything Flye prints goes to the log.
    #              2. cp copies assembly.fasta up into the sample's assembly folder, where
    #                 every later step expects the draft.
    #              3. the provenance table, written by the rule itself with printf: one
    #                 "field<TAB>value" line each for the assembler, its version (flye
    #                 --version), the sample's technology (params read_type), the input
    #                 file and the whole command line.
    # Produces:    FLYE_DIR (directory() output: a re-run starts from an empty folder) with
    #              assembly_info.txt inside it (named as an output of its own because rule
    #              assembly_evidence reads it), and two of the three standardised files:
    #              DRAFT_CONTIGS and ASSEMBLY_PROVENANCE. The third, ASSEMBLY_EVIDENCE, is
    #              made by rule assembly_evidence.
    # Consumed by: see the head of this file.
    # Log:         Flye's whole output (its stage messages and its final statistics).
    #
    # Why Flye gets a folder of its own. A directory() output is DELETED by Snakemake
    # before its rule re-runs. In FunFluxL that output was the whole sample folder (Flye's
    # --out-dir), which in this layout also holds contaminants/, mitochondria/, medaka/
    # and eval/: a re-run of the assembler would have wiped the screen, the polishing and
    # the QC of the sample. So Flye works in flye/ and the draft is copied up one level.
    #
    # v1 -> v2: the flye command is FunFluxL's with one change of place, --out-dir FLYE_DIR
    # instead of the sample folder (see above), and the cp line that follows from it. The
    # CPU count is Snakemake's own threads:, which --cores enforces, where FunFluxL used a
    # private "resources: cpus". The assembly_linearized.fasta copy is no longer written
    # (every reader of the draft reads multi-line FASTA). The log is named after the rule,
    # and the evidence table for the mitochondrion screen is written by rule
    # assembly_evidence (below).
    rule long_read_assembly:
        input:
            reads = FILTERED_LONG,
            model_check = medaka_gate,
        output:
            native_dir = directory(FLYE_DIR),
            native_info = _FLYE_INFO,
            contigs = DRAFT_CONTIGS,
            provenance = ASSEMBLY_PROVENANCE,
        params:
            read_flag = flye_read_flag,
            iterations = 5,
            read_type = read_tech,
        threads: CPUS
        # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
        resources:
            mem_mb = HEAVY_MEM_MB,
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/flye.yaml"
        log:
            LOGS + "/long_read_assembly_{sample}.log"
        shell:
            """
            mkdir -p {output.native_dir}
            flye \
                {params.read_flag} {input.reads} \
                --out-dir {output.native_dir} \
                --threads {threads} \
                --iterations {params.iterations} > {log} 2>&1

            cp {output.native_dir}/assembly.fasta {output.contigs}
            {{
                printf 'field\tvalue\n'
                printf 'assembler\tflye\n'
                printf 'version\t'
                flye --version
                printf 'read_type\t%s\n' {params.read_type:q}
                printf 'input_fastq\t%s\n' {input.reads:q}
                printf 'command\tflye %s --out-dir %s --threads %s --iterations %s\n' \
                    {params.read_flag:q} {output.native_dir:q} {threads} {params.iterations:q}
            }} > {output.provenance}
            """


# ───────────── assembly_evidence — the assembler's depth and circularity, in one table ─────────────
# The evidence table is read from the assembler's own per-contig file, and which file that is
# depends on the assembler. The choice is made once, when Snakemake reads the workflow, from
# LONG_READ_ASSEMBLER, as for long_read_assembly above: three names say which subcommand of
# mitochondria.py to run, which option hands it the file, and which file.
if LONG_READ_ASSEMBLER == "hifiasm":
    _EVIDENCE_SUBCOMMAND     = "evidence-hifiasm"
    _EVIDENCE_SOURCE_OPTION  = "--gfa"
    _EVIDENCE_SOURCE         = _HIFIASM_PRIMARY_GFA
else:
    _EVIDENCE_SUBCOMMAND     = "evidence-flye"
    _EVIDENCE_SOURCE_OPTION  = "--flye-info"
    _EVIDENCE_SOURCE         = _FLYE_INFO

# Takes in:    the assembler's own per-contig file, written by rule long_read_assembly above:
#              hifiasm's primary graph (_HIFIASM_PRIMARY_GFA) or Flye's assembly_info.txt
#              (_FLYE_INFO), whichever assembler the run uses.
# Does:        python mitochondria.py (MITOCHONDRIA_SCRIPT) with one of two subcommands:
#                evidence-hifiasm --gfa <graph>   reads the depth tag (rd:i:<n>) and the name
#                                                 of every contig in the graph; hifiasm names
#                                                 end in c (circular) or l (linear).
#                evidence-flye --flye-info <file> reads the columns #seq_name, cov., circ.,
#                                                 repeat and alt_group of Flye's table.
# Produces:    ASSEMBLY_EVIDENCE (assembly_evidence.tsv next to assembly.fasta): one row per
#              contig with its read depth, whether it is circular and, for Flye, whether it is
#              a repeat or an alternative copy of another contig. The script explains the
#              columns.
# Consumed by: rule screen_mitochondria (shared/15_mitochondria.smk), which weighs the depth
#              and the circularity next to the BLAST hits, and rule contig_summary
#              (shared/25_contig_summary.smk).
# Log:         the script's one summary line (how many contigs, and for hifiasm how many
#              carried a depth tag), or its error message.
#
# Why a rule of its own, and not two more lines at the end of long_read_assembly:
#   1. A failed table never costs the assembly. When any command of a job fails, Snakemake
#      deletes every output of that job, so a Python error after hours of assembly also
#      deleted the assembly. Now only this small job fails, the assembly stays on disk, and
#      the next run of the workflow repeats only this step.
#   2. Every FunFlux script runs with a pinned Python. This rule runs in envs/python.yaml
#      (Python 3.12 and nothing else). Inside the assembler's job, `python` came from the
#      assembler's environment: the hifiasm environment holds none, so it was whatever
#      `python` the computer had, which may be missing, old, or different on each cluster
#      node, and the Flye environment holds one only because Flye needs it.
# One CPU is enough: the script reads one file once, line by line.
#
# v1 -> v2: new rule. FunFluxL ran the same command at the end of its assembly rule.
rule assembly_evidence:
    input:
        source = _EVIDENCE_SOURCE,
    output:
        evidence = ASSEMBLY_EVIDENCE,
    params:
        script = MITOCHONDRIA_SCRIPT,
        subcommand = _EVIDENCE_SUBCOMMAND,
        source_option = _EVIDENCE_SOURCE_OPTION,
    threads: 1
    conda:
        "../../envs/python.yaml"
    log:
        LOGS + "/assembly_evidence_{sample}.log"
    shell:
        """
        python {params.script:q} {params.subcommand} \
            {params.source_option} {input.source:q} \
            --output {output.evidence:q} > {log} 2>&1
        """
