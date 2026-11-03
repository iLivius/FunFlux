# FunFlux v2.0.0 — hybrid front end, step 6: polishing and the finished assembly
# (rules/hybrid/50_polish.smk).
#
# Three rules close the assembly chain of the hybrid mode:
#   long_read_consensus   : Medaka corrects an Oxford Nanopore (ONT) assembly with the
#                           reads it was built from. ONT samples only, and only when
#                           Medaka is on. Identical to rules/longread/30_polish.smk: keep
#                           in sync (one word differs; see below).
#   short_read_correction : Polypolish corrects the assembly a second time, with the
#                           sample's selected Illumina pairs. Every sample. This mode's own
#                           rule, and the one that defines the delivered genome.
#   finalize_assembly     : names the finished genome, contigs_all.fasta: in this mode
#                           always the Polypolish result.
# The Snakefile includes every file in rules/hybrid/ in hybrid mode only, in name order;
# this file comes last, after 40_long_assembly.smk.
#
# The one word that differs in the Medaka copy: the [FunFlux] line printed when Medaka
# fails names the setting that switches Medaka off, parameters.hybrid.medaka_model here
# and parameters.longread.medaka_model in the longread rule. Everything else in that
# block, comments and commands, stays the same. Its comments name the neighbours by role,
# not by file, so that they hold in both modes: here, rule filter_long_reads is in
# rules/hybrid/30_long_reads.smk, and the consensus is read by rule short_read_correction
# below (polish_source() points at it), not by finalize_assembly.
#
# Why two polishing steps. Nanopore basecalling leaves systematic errors in the reads
# (mostly in homopolymer runs, where the same base repeats) that the assembler cannot
# average out, because every read makes the same mistake; Medaka corrects them with a
# neural network trained for the basecaller that made the reads. Illumina reads make
# different, and far fewer, errors, so a second pass with them (Polypolish) catches what
# Medaka left, and corrects a PacBio HiFi assembly too, where Medaka never runs. The
# delivered genome of this mode is therefore the long-read assembly, screened, polished
# with the long reads where that applies, and corrected with the Illumina reads. Medaka
# keeps the contig names; Polypolish appends a word to each header, which the rule
# removes, so the mitochondrion split finds the names the screen used.
#
# Where polishing sits. It runs on the SCREENED draft, not on the raw assembly: the
# contamination screen has already dropped the contigs the taxonomy selector rejected
# (rule select_contigs, shared/10_decontam.smk), so no read of a removed contig can pull
# the consensus of the kept ones. The mitochondrion is still in the file at this point; it
# is split off after polishing (rule separate_mitochondria, shared/15_mitochondria.smk).
#
#   contigs_decontam.fasta  (DECONTAM_CONTIGS, from rule select_contigs)
#     ONT sample with Medaka:
#       + {sample}_filtered.fastq (FILTERED_LONG) + {sample}_medaka_model.txt (MEDAKA_MODEL_RESOLVED)
#       --long_read_consensus--> medaka/consensus.fasta  (MEDAKA_CONSENSUS)  = polish_source()
#     every other sample (PacBio HiFi, or Medaka off):
#       contigs_decontam.fasta as it is                                     = polish_source()
#   polish_source() + {sample}_selected_R1.fastq, _R2  (SEL_R1/R2, rule map_sel_contigs in
#                                                       rules/hybrid/25_guide.smk)
#     --short_read_correction (Polypolish)--> polypolish/polished.fasta  (POLYPOLISH_CONTIGS)
#     --finalize_assembly--------------------> contigs_all.fasta          (CONTIGS_ALL)
#   contigs_all.fasta --> rule separate_mitochondria --> contigs_final.fasta + contigs_mitochondrial.fasta
#
# Who gets Medaka is decided once, when Snakemake reads the workflow: MEDAKA_SAMPLES
# (00_common.smk, section 13) holds the ONT samples of the batch when
# parameters.hybrid.medaka_model is a model name or auto, and nothing when it is false or
# the batch has no ONT sample. The whole Medaka rule sits inside `if MEDAKA_SAMPLES:`, so a
# PacBio-only batch, or a run with Medaka off, has no Medaka rule at all, and a rule-level
# wildcard constraint keeps the rule to those samples. polish_source() (00_common.smk)
# resolves the choice per sample: the Medaka consensus for the samples in MEDAKA_SAMPLES,
# the screened draft for the others. finalize_source() is POLYPOLISH_CONTIGS for every
# hybrid sample.
#
# Everything referenced here is defined once in 00_common.smk: MEDAKA_SAMPLES,
# FILTERED_LONG, DECONTAM_CONTIGS, MEDAKA_MODEL_RESOLVED, MEDAKA_DIR, MEDAKA_CONSENSUS,
# SEL_R1, SEL_R2, POLYPOLISH_DIR, POLYPOLISH_CONTIGS, CONTIGS_ALL, FASTA_HEAD_CMD,
# FASTA_LIN_CMD, LOGS and CPUS, and the functions capped_cpus(), polish_source() and
# finalize_source(). The wildcard constraint only the Medaka rule uses is defined just
# below the `if`; the Polypolish working files are named just above that rule.
#
# conda: environment paths are relative to THIS file (workflow/rules/hybrid/), so
# "../../envs/medaka.yaml" is workflow/envs/medaka.yaml. Medaka and Polypolish (with bwa)
# get an environment; finalize_assembly only copies a file and needs none.
#
# CPUs: a rule that uses several CPUs declares Snakemake's own `threads:` and passes
# {threads} to its tool. That is what makes `--cores N` actually limit the run.
#
# In the shell commands below, a doubled brace ({{ or }}) stands for a single brace in the
# shell. Snakemake reads single braces as its own placeholders, such as {input.reads}, and
# turns a doubled brace into a single one.


# The Medaka rule exists only when a sample of this run is polished with Medaka (see above).
if MEDAKA_SAMPLES:

    # The samples the Medaka rule may run for, written as a regular expression (the names
    # joined by |, "or"; re.escape() makes a dot in a name match only a dot). The
    # rule-level wildcard constraint below narrows the workflow-wide one (00_common.smk,
    # section 14), so a PacBio sample's consensus can never be asked for, not even by
    # mistake. 12_medaka_check.smk builds the same expression under a name of its own:
    # it is a shared file and must not depend on a name that only this front end defines.
    _CONSENSUS_SAMPLE_CONSTRAINT = "|".join(re.escape(sample) for sample in MEDAKA_SAMPLES)

    # ───────────── long_read_consensus — Medaka polishes the screened draft ─────────────
    # Takes in:    FILTERED_LONG, the Filtlong output (rule filter_long_reads of this
    #              mode's front end): the same reads the draft was assembled from;
    #              DECONTAM_CONTIGS, the screened draft (rule select_contigs,
    #              shared/10_decontam.smk); and MEDAKA_MODEL_RESOLVED, the one-line file
    #              rule check_medaka_model wrote before the assembler ran (shared/
    #              12_medaka_check.smk): the model name the user set, after the check, or
    #              the model Medaka's auto mode worked out from the read headers (then a
    #              path to the model file; medaka_consensus -m accepts both).
    # Does:        reads that line into the shell variable model, then medaka_consensus:
    #              -i the reads, -d the draft, -x (rebuild the draft's alignment index;
    #              see below), -t the CPUs, -m the model, -o the output folder. Medaka
    #              aligns the reads to the draft (minimap2), runs the model over the
    #              pile-up and writes the corrected sequence. Everything it prints goes to
    #              the log.
    # Produces:    MEDAKA_DIR (02.assembly/{sample}/medaka/), declared as a directory()
    #              output, so Snakemake starts a re-run from an empty folder and stale
    #              intermediate files of an earlier run (calls_to_draft.bam, the
    #              consensus_probs.hdf of the network) cannot be mistaken for fresh ones;
    #              and MEDAKA_CONSENSUS, consensus.fasta inside it, the polished assembly,
    #              declared as well so that the rule only counts as done when the FASTA
    #              exists. Both are kept.
    # Consumed by: the next step of this mode's assembly chain, which 00_common.smk points
    #              at per sample: rule finalize_assembly through finalize_source() in
    #              longread mode, rule short_read_correction (Polypolish) through
    #              polish_source() in hybrid mode. Rule all asks for the consensus too.
    # Log:         Medaka's own messages, and on failure one [FunFlux] line.
    #
    # Why -x. Medaka writes two index files of the draft (contigs_decontam.fasta.fai and
    # .map-ont.mmi) next to the draft, not in Medaka's folder, and reuses them when they
    # are there; a re-run of the screen (a changed setting, a new assembly) rewrites the
    # draft but leaves the old index files behind. -x rebuilds them on every run (a few
    # seconds). Without it, Medaka would align the reads to the contigs of the earlier
    # draft: a contig can drop out of the polished genome, or get another contig's
    # sequence under its name, and Medaka reports no error.
    #
    # The failure guard (the part after ||) runs only when Medaka fails. It appends one
    # [FunFlux] line to the log, prints the same line on the console (tee -a writes to the
    # log and to standard output, which >&2 turns into standard error) and stops the run:
    # the polished assembly is a real result, not a report, so unlike the QC rules this one
    # must not be skipped.
    #
    # CPUs: capped at 12, the FunFluxL value. Polishing is one job per ONT sample, and the
    # cap leaves CPUs for the other samples' jobs that run at the same time.
    #
    # v1 -> v2: the medaka_consensus command is FunFluxL's plus -x (see above), in one
    # rule instead of two. FunFluxL had one rule for a named model, which first checked
    # the name against `medaka tools list_models` here, after the assembly, and one for
    # auto mode, which passed no -m and let medaka_consensus work the model out itself,
    # late and with a terse error when the read headers carried no basecaller tag. Both
    # decisions now happen before the assembler, in rule check_medaka_model, and this rule
    # reads the answer from its file: with a named model the command line is FunFluxL's
    # plus -x; in auto mode -m now names the model file the check found. The CPU count is
    # Snakemake's own threads: (capped at 12, as FunFluxL's private "resources: cpus"
    # was). The folder moved from 03.post-processing/consensus/{sample}/ to the sample's
    # assembly folder, and the log is named after the rule (FunFluxL:
    # consensus_long_{sample}.log).
    rule long_read_consensus:
        input:
            reads = FILTERED_LONG,
            contigs = DECONTAM_CONTIGS,
            model = MEDAKA_MODEL_RESOLVED,
        output:
            consensus_dir = directory(MEDAKA_DIR),
            consensus = MEDAKA_CONSENSUS,
        wildcard_constraints:
            sample = _CONSENSUS_SAMPLE_CONSTRAINT,
        threads: capped_cpus(12)
        # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
        resources:
            mem_mb = HEAVY_MEM_MB,
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/medaka.yaml"
        log:
            LOGS + "/long_read_consensus_{sample}.log"
        shell:
            """
            model=$(cat {input.model})

            medaka_consensus \
                -i {input.reads} \
                -d {input.contigs} \
                -x \
                -t {threads} \
                -m "$model" \
                -o {output.consensus_dir} > {log} 2>&1 || {{
                echo "[FunFlux] long_read_consensus: Medaka failed for sample {wildcards.sample} with model '$model' (the line in {input.model}); the Medaka messages are in {log}. The model itself passed rule check_medaka_model, so look there for another cause (memory, a damaged read file). To run without polishing, set parameters.hybrid.medaka_model to false." | tee -a {log} >&2
                exit 1
            }}
            """


# ── Polypolish's working files: only rule short_read_correction uses them ──
# All in POLYPOLISH_DIR (02.assembly/{sample}/polypolish/), the folder that also holds the
# result, polished.fasta (POLYPOLISH_CONTIGS, 00_common.smk). The draft is COPIED into this
# folder before it is indexed: bwa index writes its five index files (.amb, .ann, .bwt,
# .pac, .sa) next to the FASTA it indexes, and no rule should write next to another rule's
# output (the Medaka consensus or the screened draft). The one exception is Medaka itself:
# it leaves its index files next to the screened draft and rebuilds them on every run (-x
# in long_read_consensus). Two alignments, one per mate file, then the two filtered
# alignments, then Polypolish's own output before the header cut. All six are temporary
# (temp() in the rule): Snakemake deletes them when the job is done.
_POLYPOLISH_DRAFT        = POLYPOLISH_DIR + "/draft.fasta"
_POLYPOLISH_SAM_1        = POLYPOLISH_DIR + "/align_1.sam"
_POLYPOLISH_SAM_2        = POLYPOLISH_DIR + "/align_2.sam"
_POLYPOLISH_FILT_SAM_1   = POLYPOLISH_DIR + "/filtered_align_1.sam"
_POLYPOLISH_FILT_SAM_2   = POLYPOLISH_DIR + "/filtered_align_2.sam"
_POLYPOLISH_RAW_POLISHED = POLYPOLISH_DIR + "/polished_raw.fasta"


# ───────────── short_read_correction — Polypolish corrects the assembly with the Illumina pairs ─────────────
# The last polishing step, and the one that defines the delivered genome of this mode. The
# selected Illumina pairs are aligned to the assembly and Polypolish corrects each position
# where the aligned reads clearly agree on a different base from the assembly's. Polypolish
# is written for exactly this: it wants the two mate files aligned SEPARATELY and every
# read reported at ALL the places it could come from (bwa mem -a), which lets it polish
# repeats too (a read from a repeat counts for every copy of it), and it is conservative:
# a position where the reads do not agree with each other is left as it was.
#
# Takes in:    polish_source(): the draft to correct, per sample. MEDAKA_CONSENSUS (rule
#              long_read_consensus above) for an ONT sample with Medaka, otherwise
#              DECONTAM_CONTIGS, the screened long-read assembly (rule select_contigs,
#              shared/10_decontam.smk). Snakemake calls the function once per sample, with
#              the sample's name, so each job gets its own source file. And SEL_R1, SEL_R2,
#              the Illumina pairs selected on the screened guide (rule map_sel_contigs,
#              rules/hybrid/25_guide.smk): only the reads the guide screen kept correct the
#              genome.
# Does:        1. cp copies the draft into the rule's own folder (see the note above the
#                 working files).
#              2. bwa index builds bwa's index of the copy; then bwa mem -a twice, once per
#                 mate file (-t: the threads), each into its own SAM. -a reports every
#                 alignment of a read, not only its best one.
#              3. polypolish filter reads both SAMs and removes the alignments whose insert
#                 size (the distance between the two mates of a pair) is out of the
#                 ordinary: a mate that landed far from its partner is in the wrong place.
#                 It needs both files to know which alignments make a pair.
#              4. polypolish polish writes the corrected FASTA to standard output and, in
#                 doing so, appends the word "polypolish" to every header.
#              5. two awk programs from 00_common.smk, in a pipe: FASTA_HEAD_CMD cuts every
#                 header at its first blank or tab, so the appended word goes and the
#                 contig NAME, the first word, stays (every later step identifies a contig
#                 by it: the mitochondrion split, QUAST, BUSCO, the ITS search,
#                 funannotate); then FASTA_LIN_CMD writes each sequence on one line, as
#                 v1's polished.fasta was written.
#              6. rm -f removes the five bwa index files, which are not declared outputs
#                 (they would otherwise stay behind next to the deleted draft copy).
# Produces:    POLYPOLISH_CONTIGS (polypolish/polished.fasta), kept; rule all asks for it.
#              The six working files above are temp(). Disk: the two SAMs hold EVERY
#              alignment of EVERY selected read as text; for a deep Illumina library they
#              are tens of gigabytes each while the job runs.
# Consumed by: finalize_assembly below (finalize_source() points here for every hybrid
#              sample).
# Log:         bwa's messages, then Polypolish's (how many alignments the filter kept, how
#              many positions were changed).
#
# Polypolish runs for every hybrid sample, PacBio HiFi included: a HiFi assembly is
# accurate already, so it changes little there, but the step is the same for the whole
# batch and the log says how much it changed.
#
# Threads: up to 24 (capped_cpus), for bwa mem (-t). Polypolish itself runs on one CPU.
#
# v1 -> v2: the cp, bwa and polypolish commands are FunFluxL+'s (rule short_read_correction),
# unchanged. The one change is the header cut: FunFluxL+ wrote the result with FASTA_LIN_CMD
# alone and left " polypolish" in every header, which the mitochondrion split ignored (it
# reads the first word) but QUAST, BUSCO and the ITS search showed; FASTA_HEAD_CMD now runs
# first. The working files moved from 06.polishing/polypolish/{sample}/ to the sample's
# polypolish folder, and the CPU count is Snakemake's own threads: (capped at 24, as
# FunFluxL+'s private "resources: cpus" was), which --cores enforces.
rule short_read_correction:
    input:
        contigs = polish_source,
        r1 = SEL_R1,
        r2 = SEL_R2,
    output:
        draft_copy = temp(_POLYPOLISH_DRAFT),
        sam_1 = temp(_POLYPOLISH_SAM_1),
        sam_2 = temp(_POLYPOLISH_SAM_2),
        filt_sam_1 = temp(_POLYPOLISH_FILT_SAM_1),
        filt_sam_2 = temp(_POLYPOLISH_FILT_SAM_2),
        raw_polished = temp(_POLYPOLISH_RAW_POLISHED),
        polished = POLYPOLISH_CONTIGS,
    threads: capped_cpus(24)
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/polypolish.yaml"
    log:
        LOGS + "/short_read_correction_{sample}.log"
    shell:
        """
        cp {input.contigs} {output.draft_copy}
        bwa index {output.draft_copy} > {log} 2>&1
        bwa mem \
            -t {threads} \
            -a {output.draft_copy} \
            {input.r1} > {output.sam_1} 2>> {log}
        bwa mem \
            -t {threads} \
            -a {output.draft_copy} \
            {input.r2} > {output.sam_2} 2>> {log}
        polypolish filter \
            --in1 {output.sam_1} \
            --in2 {output.sam_2} \
            --out1 {output.filt_sam_1} \
            --out2 {output.filt_sam_2} >> {log} 2>&1
        polypolish polish \
            {output.draft_copy} \
            {output.filt_sam_1} \
            {output.filt_sam_2} > {output.raw_polished} 2>> {log}
        awk {FASTA_HEAD_CMD:q} {output.raw_polished} | \
        awk {FASTA_LIN_CMD:q} > {output.polished}
        rm -f \
            {output.draft_copy}.amb \
            {output.draft_copy}.ann \
            {output.draft_copy}.bwt \
            {output.draft_copy}.pac \
            {output.draft_copy}.sa
        """


# ───────────── finalize_assembly — name the finished genome ─────────────
# Why a rule that only copies a file: every mode ends its assembly chain with the same
# file name, contigs_all.fasta, so the mitochondrion split and everything after it read
# one name whichever way the genome arrived. In this mode the source is always the
# Polypolish result; the per-sample Medaka decision was resolved one step earlier, in
# polish_source(). The dry run prints that decision per sample, and a sample polished
# with Medaka has a medaka/ folder.
#
# Takes in:    finalize_source(): POLYPOLISH_CONTIGS (rule short_read_correction above) for
#              every hybrid sample (00_common.smk, section 13). Snakemake calls the function
#              once per sample, with the sample's name.
# Does:        copies it and records the source in the log.
# Produces:    CONTIGS_ALL, the complete screened and polished genome, mitochondrion
#              included.
# Consumed by: rule separate_mitochondria in shared/15_mitochondria.smk, which splits it
#              into the nuclear genome (contigs_final.fasta) and the mitochondrial
#              contigs.
#
# v1 -> v2: new rule. FunFluxL+ had no file for the finished genome; a helper function
# handed polished.fasta to the mitochondrion split at the moment it asked for its input.
# Now the file exists under one name in every mode and the log says where it came from.
# The copy is the same as in the other modes' finalize_assembly.
rule finalize_assembly:
    input:
        contigs = finalize_source,
    output:
        contigs = CONTIGS_ALL,
    log:
        LOGS + "/finalize_assembly_{sample}.log"
    shell:
        """
        echo "contigs_all.fasta of {wildcards.sample} is a copy of {input.contigs}" > {log}
        cp {input.contigs} {output.contigs}
        """
