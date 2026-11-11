# FunFlux v2.1.0 — mitochondrion screen (rules/shared/15_mitochondria.smk): find the
# mitochondrial genome among the contigs and set it apart from the nuclear genome.
#
# Why. A fungal cell carries two genomes: the nuclear genome (the chromosomes, tens of
# megabases) and a small mitochondrial genome (usually some tens of kb) that exists in many
# copies per cell. An assembler cannot tell them apart and
# delivers the mitochondrial genome as one more ordinary contig. Left in the assembly, that
# contig would go through the nuclear gene predictor, which is wrong: fungal mitochondria
# read the genetic code differently from the nucleus (code 4 in the filamentous fungi: UGA
# means tryptophan, not "stop"), so the predictor would call false genes on it. This file
# finds the mitochondrial contig and writes it to a file of its own, so that QUAST, BUSCO,
# the ITS search and funannotate only ever see the nuclear genome.
#
# Two rules, at two different places in the workflow:
#   screen_mitochondria    runs EARLY, on the draft assembly, right after the BLAST (rule
#                          blast_contigs in shared/10_decontam.smk). It reads that BLAST
#                          table, decides which contigs are mitochondrial, and writes a
#                          report plus the list of confirmed contig IDs.
#   separate_mitochondria  runs LATE, when the genome is finished (decontaminated and, in
#                          the long-read modes, polished). It reads that genome
#                          (CONTIGS_ALL) and the list, and writes the nuclear genome
#                          (FINAL_CONTIGS) and the mitochondrial contigs (MITO_CONTIGS).
#
# The order matters. The list is made on the draft, before the contig selector runs, because
# the selector (select_contigs in shared/10_decontam.smk) has to read it: it keeps the listed
# contigs whatever the genus filter says. If the screen waited until after the selector, a
# genus filter could already have thrown the mitochondrion away.
#
# The chain, step by step:
#   1. screen_mitochondria    BLASTOUT + DRAFT_CONTIGS + nodes.dmp (the NCBI taxonomy)
#                             -> MITO_SCREEN_TSV (the report, for people) + CONFIRMED_MITO_IDS
#   2. select_contigs         (shared/10_decontam.smk) reads CONFIRMED_MITO_IDS and keeps those
#                             contigs whatever the genus filter says
#                             -> DECONTAM_CONTIGS
#   3. polishing              longread and hybrid only (Medaka for ONT samples, Polypolish in
#                             hybrid); then, in every mode, finalize_assembly -> CONTIGS_ALL
#                             (the finished genome, mitochondrion still in it)
#   4. separate_mitochondria  CONTIGS_ALL + CONFIRMED_MITO_IDS
#                             -> FINAL_CONTIGS (the nuclear genome) + MITO_CONTIGS (the
#                                mitochondrial contigs; an empty file if there are none)
#
# How a contig is called mitochondrial. The script,
# workflow/scripts/15_mitochondria/mitochondria.py, never uses circularity or read depth
# alone. It starts from the BLAST table and looks for a "direct organelle hit": a hit whose
# reference title says the reference IS a mitochondrial genome (for example "Fusarium
# graminearum mitochondrion, complete sequence"). A nuclear gene that is merely described as
# "mitochondrial ribosomal protein" does not count. For each contig it takes the direct hit
# that covers the largest percentage of the contig (the BLAST column qcovs). Then:
#   * Without assembler evidence (contigs and illumina mode, where BLAST is the only
#     information): the contig is confirmed when that hit covers at least
#     parameters.mitochondria.min_direct_blast_qcov percent of it (35 by default).
#   * With assembler evidence (the long-read modes): a hit of that size is confirmed only
#     if the assembler agrees as well, that is, the contig is much deeper than the typical
#     nuclear contig (at least parameters.mitochondria.min_assembler_coverage_ratio times
#     the median depth, 3 by default) or it is a clean circle (closed by the assembler, not
#     a repeat and not an alternative copy of another contig). A hit that covers at least
#     parameters.mitochondria.very_strong_qcov percent of the contig (70 by default) is
#     confirmed on its own.
#   * In every mode, a contig whose best direct hit comes from an animal or a plant (human
#     mitochondrial DNA from the library, for example) is not confirmed: it is reported as
#     foreign_mitochondrial and stays in the assembly (the species check, explained at rule
#     screen_mitochondria).
# Every contig gets a row in the report with its class and the reason, so each call can be
# checked afterwards.
#
# parameters.mitochondria.mode decides what the screen does with its findings:
#   off                 no contig is examined (every row says not_screened) and the list of
#                       confirmed IDs is empty.
#   report              the classes are written, but the list of confirmed IDs is left empty,
#                       so nothing is protected and nothing is split off. The assembly goes
#                       on exactly as it was.
#   separate_confirmed  (default) the confirmed contigs are written to the list.
# Both rules exist whenever the screen is on, whatever the mode. With off or report they
# simply produce an empty list and an empty MITO_CONTIGS, so a finished sample has the same
# set of files in every case.
#
# Both rules are inside `if SCREEN:`. SCREEN (00_common.smk) is true in the three read-based
# modes and, in contigs mode, only when parameters.contigs.screen is true. The `if` is
# decided once, when Snakemake reads the workflow, before any job runs. When SCREEN is false
# this file defines no rule, and in contigs mode rule filter_contigs writes the final genome
# (FINAL_CONTIGS) directly.
#
# These names come from 00_common.smk, where each is defined once and never worked out
# again: SCREEN, HAS_LONG_READS, BLASTOUT, DRAFT_CONTIGS, ASSEMBLY_EVIDENCE, MITO_SCREEN_TSV,
# CONFIRMED_MITO_IDS, CONTIGS_ALL, FINAL_CONTIGS, MITO_CONTIGS, MITOCHONDRIA_SCRIPT, MITO_MODE,
# MITO_MIN_QCOV, MITO_MIN_RATIO, MITO_VERY_STRONG_QCOV, BLASTDB (the folder that holds
# nodes.dmp) and LOGS.
#
# conda: the path is written relative to THIS file (workflow/rules/shared/), so
# "../../envs/python.yaml" climbs shared/ -> rules/ -> workflow/ and then goes into envs/.
# The script uses only the Python standard library and calls no other program, so both rules
# run in envs/python.yaml, which holds a pinned Python 3.12 and nothing else. Activating that
# environment makes `python` that Python, instead of whatever `python` the computer has (a
# computer with only `python3`, or an old Python, would stop the job).


# Both rules below exist only when the contamination and mitochondrion screen is on.
if SCREEN:

    # ── The assembler evidence, or none: decided once for the whole run ──
    # In the long-read modes rule assembly_evidence (longread/20_assembly.smk and
    # hybrid/40_long_assembly.smk) writes ASSEMBLY_EVIDENCE next to the assembly, and the
    # screen reads it: it is listed as an input (so Snakemake waits for it) and its path is
    # handed to the script. Contigs and illumina mode have no such
    # table: the input is an empty list, which Snakemake accepts as "nothing to wait for",
    # and the script receives an empty value (see the rule for how that is written).
    if HAS_LONG_READS:
        _EVIDENCE_INPUT = ASSEMBLY_EVIDENCE
        _EVIDENCE_PATH = ASSEMBLY_EVIDENCE
    else:
        _EVIDENCE_INPUT = []
        _EVIDENCE_PATH = ""

    # ───────────── screen_mitochondria: decide which contigs are mitochondrial ─────────────
    # Takes in:    BLASTOUT, the BLAST table of the draft (from blast_contigs in
    #              shared/10_decontam.smk), and DRAFT_CONTIGS, the contigs that BLAST was run
    #              on. The FASTA gives the contig lengths and makes sure that every contig,
    #              with or without a hit, gets a row in the report. In the long-read modes
    #              also ASSEMBLY_EVIDENCE, the assembler's table (see below). And nodes.dmp,
    #              the family tree of the NCBI taxonomy, from the BLAST database folder (the
    #              file BlobTools reads too), for the species check below.
    # Does:        python mitochondria.py screen (MITOCHONDRIA_SCRIPT). The settings come from
    #              the configuration through `params:`, so that Snakemake re-runs the rule
    #              when one changes:
    #                --mode                         parameters.mitochondria.mode
    #                --min-direct-blast-qcov        parameters.mitochondria.min_direct_blast_qcov
    #                --min-assembler-coverage-ratio parameters.mitochondria.min_assembler_coverage_ratio
    #                --very-strong-qcov             parameters.mitochondria.very_strong_qcov
    #                --assembly-evidence            the assembler's table, or empty (see below)
    #                --taxonomy-nodes               nodes.dmp, for the species check
    # Produces:    MITO_SCREEN_TSV: the report, one row per contig, with its class and the
    #              reason (confirmed_mitochondrial, mitochondrial_candidate,
    #              foreign_mitochondrial, not_mitochondrial, ...), the best organelle hit and,
    #              in the long-read modes, the assembler's numbers. In contigs mode it is
    #              also a delivered file: rule all asks for it.
    #              CONFIRMED_MITO_IDS: the confirmed contigs, one name per line (empty unless
    #              the mode is separate_confirmed).
    # Consumed by: select_contigs in shared/10_decontam.smk (protects the listed contigs from
    #              the genus filter) and separate_mitochondria below.
    #
    # The species check. The script confirms a contig from the wording of its best BLAST
    # match ("... mitochondrion, complete genome"), which says that the match IS a
    # mitochondrial genome but not whose. A few human reads in an Illumina library can
    # assemble into small contigs of human mitochondrial DNA that pass every other check.
    # So the script looks the organism of the best match up in the taxonomy: an animal or a
    # plant, and the contig is reported as foreign_mitochondrial, neither confirmed nor
    # split off. It stays in the assembly, where the genus filter of the contamination
    # screen treats it like any other contig (with the default,
    # parameters.decontamination.mode: off, it is kept), and the contig summary gives it the
    # verdict foreign_mitochondrion. FunFlux 2.0.0 confirmed these animal and plant contigs
    # and wrote them to contigs_mitochondrial.fasta with the fungus's own mitochondrion. A
    # best match that the taxonomy places in none of the three kingdoms (fungi, animals,
    # plants) changes nothing: an environmental record filed under Eukaryota only, or the
    # mitochondrion of an oomycete, which NCBI does not count among the fungi.
    #
    # The assembler evidence. In the long-read modes rule assembly_evidence writes a small
    # table next to the assembly, made from the assembler's own output (hifiasm or Flye):
    # the read depth of every contig, whether it is circular and (Flye only) whether it is
    # a repeat or an alternative copy of another contig (ASSEMBLY_EVIDENCE in 00_common.smk).
    # There the table is an input of this rule and its path goes to --assembly-evidence.
    # Contigs and illumina mode have no such table (no long-read assembler; SPAdes writes
    # the depth into its headers but not circularity), so BLAST is the only information and
    # --assembly-evidence is given an EMPTY value. The option is written with an explicit
    # pair of double quotes so that the script always receives a value, an empty one in
    # those modes: the script's argument parser rejects a bare flag without a value. With no
    # assembler evidence the script decides on --min-direct-blast-qcov alone, so
    # --min-assembler-coverage-ratio and --very-strong-qcov change nothing there; they are
    # passed anyway so that the command has the same shape in every mode.
    #
    # Notes on the command. The ratio option is --min-assembler-coverage-ratio (the script
    # also accepts its older name, --min-flye-coverage-ratio). --very-strong-qcov passes the
    # cut-off of 70 percent, which is an option of the script rather than a fixed number.
    # --assembly-evidence "" is passed in the modes without assembler evidence. The script is
    # version 0.2 of the mitochondrial screen with three small changes made for FunFlux v2
    # and the species check added in FunFlux 2.1.0. Its report has the columns assembler_*,
    # three informational columns at the end (best_mito_slen, length_ratio,
    # concatemer_flag) and, last, best_mito_kingdom, the answer of the species check
    # (Fungi, Metazoa, Viridiplantae, or empty). The files go to the sample's mitochondria folder,
    # 02.assembly/<sample>/mitochondria/ (MITO_DIR). The script runs with the pinned Python of
    # envs/python.yaml.
    rule screen_mitochondria:
        input:
            blast = BLASTOUT,
            contigs = DRAFT_CONTIGS,
            assembly_evidence = _EVIDENCE_INPUT,
            nodes = os.path.join(BLASTDB, "nodes.dmp"),
        output:
            report = MITO_SCREEN_TSV,
            confirmed_ids = CONFIRMED_MITO_IDS,
        params:
            script = MITOCHONDRIA_SCRIPT,
            assembly_evidence = _EVIDENCE_PATH,
            mode = MITO_MODE,
            minimum_qcov = MITO_MIN_QCOV,
            minimum_coverage_ratio = MITO_MIN_RATIO,
            very_strong_qcov = MITO_VERY_STRONG_QCOV,
        conda:
            "../../envs/python.yaml"
        log:
            LOGS + "/screen_mitochondria_{sample}.log"
        shell:
            """
            python {params.script:q} screen \
              --blast {input.blast:q} \
              --contigs {input.contigs:q} \
              --assembly-evidence "{params.assembly_evidence}" \
              --taxonomy-nodes {input.nodes:q} \
              --report {output.report:q} \
              --confirmed-ids {output.confirmed_ids:q} \
              --mode {params.mode:q} \
              --min-direct-blast-qcov {params.minimum_qcov:q} \
              --min-assembler-coverage-ratio {params.minimum_coverage_ratio:q} \
              --very-strong-qcov {params.very_strong_qcov:q} > {log} 2>&1
            """


    # ───────────── separate_mitochondria: nuclear genome and mitochondrial contigs ─────────────
    # Takes in:    CONTIGS_ALL, the complete finished genome with the mitochondrion still in it,
    #              written by the mode's own rule finalize_assembly (in contigs/20_assembly.smk,
    #              illumina/20_assembly.smk, longread/30_polish.smk or hybrid/50_polish.smk,
    #              after any polishing), and CONFIRMED_MITO_IDS from screen_mitochondria.
    # Does:        python mitochondria.py split. A contig whose name is on the list goes to the
    #              mitochondrial file, every other contig to the nuclear file. Order and
    #              headers are not changed; the sequence is written 80 bases per line. The list
    #              holds contig NAMES (the first word of the FASTA header), so every step
    #              between the screen and the split must keep them. The script stops with a
    #              message, before it writes anything, if a listed name is missing from
    #              CONTIGS_ALL (the names changed on the way) or if nothing would be left for
    #              the nuclear genome.
    # Produces:    FINAL_CONTIGS: the NUCLEAR genome. Every later step that looks at the
    #              assembly (assembly statistics, ITS search, gene prediction) reads this file
    #              and no other.
    #              MITO_CONTIGS: the confirmed mitochondrial contig(s), kept as a delivered
    #              file; an empty file when none was confirmed (or the mode is off or report).
    # Consumed by: FINAL_CONTIGS by rule genome_assembly_evaluation and rule completeness
    #              (QUAST and BUSCO, shared/20_qc.smk), rule its_windows (shared/30_taxonomy.smk)
    #              and rule funannotate_preprocess (shared/40_annotation.smk). With the AMR
    #              screen on, rule amr_chroquetas (shared/50_amr.smk) reads both files:
    #              FINAL_CONTIGS for the nuclear targets, MITO_CONTIGS for cytochrome b, which
    #              the mitochondrion encodes. MITO_CONTIGS is a delivered file; in contigs mode
    #              rule all asks for it.
    #
    # Notes: the input is CONTIGS_ALL, and the two outputs sit in the sample's assembly folder,
    # contigs_final.fasta (the nuclear genome) and contigs_mitochondrial.fasta. The script runs
    # with the pinned Python of envs/python.yaml.
    rule separate_mitochondria:
        input:
            contigs = CONTIGS_ALL,
            confirmed_ids = CONFIRMED_MITO_IDS,
        output:
            nuclear = FINAL_CONTIGS,
            mitochondrial = MITO_CONTIGS,
        params:
            script = MITOCHONDRIA_SCRIPT,
        conda:
            "../../envs/python.yaml"
        log:
            LOGS + "/separate_mitochondria_{sample}.log"
        shell:
            """
            python {params.script:q} split \
              --fasta {input.contigs:q} \
              --confirmed-ids {input.confirmed_ids:q} \
              --nuclear-fasta {output.nuclear:q} \
              --mitochondrial-fasta {output.mitochondrial:q} > {log} 2>&1
            """
