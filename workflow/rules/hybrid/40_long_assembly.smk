# FunFlux v2.1.0 — hybrid front end, step 5: the long-read assembler
# (rules/hybrid/40_long_assembly.smk).
#
# Identical to rules/longread/20_assembly.smk: keep in sync. Everything below this header
# comment, from the two rule-local names onwards, is a copy of that file, comments and
# commands alike; only this header is the hybrid mode's own. A change there must be made
# here too, and the other way round. To check, diff the two files from the comment above
# _HIFIASM_PRIMARY_GFA to the end.
#
# Two rules. long_read_assembly comes in two versions, hifiasm (the default) and Flye; only
# ONE exists in a run, chosen once from parameters.hybrid.assembler when Snakemake reads the
# workflow. It builds the draft genome from the filtered long reads. assembly_evidence then
# turns the assembler's own per-contig file into the evidence table. Together they end with
# the same three files whichever assembler ran: assembly.fasta (DRAFT_CONTIGS), the
# per-contig evidence table for the mitochondrion screen (ASSEMBLY_EVIDENCE) and the
# provenance table (ASSEMBLY_PROVENANCE). The rule comments below explain both versions;
# why the two tables exist is in the header of rules/longread/20_assembly.smk.
#
# What is different in hybrid mode is not the assembler but what went in and what comes
# after. The filtered reads of an ONT sample were chosen with the Illumina pairs as the
# yardstick (rule filter_long_reads in rules/hybrid/30_long_reads.smk, when
# parameters.hybrid.short_read_guidance is on): Filtlong scored each ONT read by how much
# of it the Illumina reads support, and cut the reads where they do not. For a PacBio HiFi
# sample the reads are filtered as in longread mode. The draft then goes through the
# shared screen (shared/10_decontam.smk, shared/15_mitochondria.smk), Medaka for ONT
# samples and, new in this mode, Polypolish with the selected Illumina pairs (both in
# rules/hybrid/50_polish.smk) before it is named contigs_all.fasta. The assembler itself
# never sees an Illumina read.
#
# The Snakefile includes every file in rules/hybrid/ in hybrid mode only, in name order;
# this file comes after 30_long_reads.smk and before 50_polish.smk.
#
# Where the rule comments below name the longread mode's neighbours and settings, read the
# hybrid ones: rule filter_long_reads is in rules/hybrid/30_long_reads.smk, not
# 10_reads.smk; the settings are parameters.hybrid.assembler, parameters.hybrid.hifiasm and
# parameters.hybrid.flye_input_mode (00_common.smk reads the block named after the mode,
# with the same keys and the same checks); and the technology flags come from the sample's
# file name in input.longread_dir, exactly as in longread mode. The Medaka model check the
# assembler waits for (rule check_medaka_model in shared/12_medaka_check.smk) reads the RAW
# ONT reads, which matters here: with the Illumina guidance on, the filtered reads exist only
# after the whole guide arm has finished, and the check must answer before that.
#
# Everything referenced here is defined once in 00_common.smk: LONG_READ_ASSEMBLER,
# HIFIASM_COMMON_ARGS, FILTERED_LONG, HIFIASM_DIR, FLYE_DIR, DRAFT_CONTIGS,
# ASSEMBLY_EVIDENCE, ASSEMBLY_PROVENANCE, MITOCHONDRIA_SCRIPT, LOGS and CPUS, and the
# per-sample functions hifiasm_tech_flag(), flye_read_flag(), read_tech() and
# medaka_gate(). Two file names only this file uses are defined just below, and the three
# names that set up the evidence step just above rule assembly_evidence.
#
# conda: environment paths are relative to THIS file (workflow/rules/hybrid/), so
# "../../envs/hifiasm.yaml" is workflow/envs/hifiasm.yaml. Which environment runs the
# evidence script, why the assembler gets `threads: CPUS`, and what a doubled brace in a
# shell command means: see the header of rules/longread/20_assembly.smk; all three hold
# here unchanged.


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
