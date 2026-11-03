# FunFlux v2.0.0 — longread front end, step 3: polishing and the finished assembly
# (rules/longread/30_polish.smk).
#
# Two rules close the assembly chain of the longread mode:
#   long_read_consensus : Medaka corrects an Oxford Nanopore (ONT) assembly with the reads
#                         it was built from. ONT samples only, and only when Medaka is on.
#   finalize_assembly   : names the finished genome, contigs_all.fasta, for every sample:
#                         the Medaka consensus where Medaka ran, the screened draft as it
#                         is everywhere else.
# The Snakefile includes this file in longread mode only, after 20_assembly.smk.
#
# Why polish, and why only ONT. Nanopore basecalling leaves systematic errors in the reads
# (mostly in homopolymer runs, where the same base repeats) that the assembler cannot
# average out, because every read makes the same mistake. Medaka reads the pile-up of all
# reads on every position with a neural network trained for the basecaller that made the
# reads, and corrects the sequence. A wrong model, one trained for another chemistry or
# basecaller version, polishes badly, which is why the model is checked before the
# assembler starts (rule check_medaka_model, shared/12_medaka_check.smk) and why its answer
# is written down and read here. PacBio HiFi reads are already accurate (consensus of many
# passes over the same molecule) and never go through Medaka, whatever the config says.
#
# Where polishing sits. It runs on the SCREENED draft, not on the raw assembly: the
# contamination screen has already dropped the contigs the taxonomy selector rejected
# (rule select_contigs, shared/10_decontam.smk), so no read of a removed contig can pull
# the consensus of the kept ones. The mitochondrion is still in the file at this point;
# it is split off after polishing (rule separate_mitochondria, shared/15_mitochondria.smk).
#
#   contigs_decontam.fasta  (DECONTAM_CONTIGS, from rule select_contigs)
#     ONT sample with Medaka:
#       + {sample}_filtered.fastq (FILTERED_LONG) + {sample}_medaka_model.txt (MEDAKA_MODEL_RESOLVED)
#       --long_read_consensus--> medaka/consensus.fasta  (MEDAKA_CONSENSUS)
#       --finalize_assembly----> contigs_all.fasta       (CONTIGS_ALL)
#     every other sample (PacBio HiFi, or Medaka off):
#       --finalize_assembly----> contigs_all.fasta       (CONTIGS_ALL)
#   contigs_all.fasta --> rule separate_mitochondria --> contigs_final.fasta + contigs_mitochondrial.fasta
#
# Who gets Medaka is decided once, when Snakemake reads the workflow: MEDAKA_SAMPLES
# (00_common.smk, section 13) holds the ONT samples of the batch when
# parameters.longread.medaka_model is a model name or auto, and nothing when it is false
# or the batch has no ONT sample. The whole Medaka rule sits inside `if MEDAKA_SAMPLES:`,
# so a PacBio-only batch, or a run with Medaka off, has no Medaka rule at all, and a
# rule-level wildcard constraint keeps the rule to those samples. finalize_source()
# (00_common.smk) resolves the choice per sample: the Medaka consensus for the samples in
# MEDAKA_SAMPLES, the screened draft for the others.
#
# The hybrid mode shares this step: rules/hybrid/50_polish.smk holds a copy of the Medaka
# rule, then Polypolish, which corrects the assembly a second time with the sample's
# Illumina reads, and its own finalize_assembly, which takes the Polypolish result.
#
# Everything referenced here is defined once in 00_common.smk: MEDAKA_SAMPLES,
# FILTERED_LONG, DECONTAM_CONTIGS, MEDAKA_MODEL_RESOLVED, MEDAKA_DIR, MEDAKA_CONSENSUS,
# CONTIGS_ALL, LOGS and CPUS, and the functions capped_cpus() and finalize_source(). The
# wildcard constraint only this file uses is defined just below the `if`.
#
# conda: environment paths are relative to THIS file (workflow/rules/longread/), so
# "../../envs/medaka.yaml" is workflow/envs/medaka.yaml. finalize_assembly only copies a
# file and needs no environment.
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
                echo "[FunFlux] long_read_consensus: Medaka failed for sample {wildcards.sample} with model '$model' (the line in {input.model}); the Medaka messages are in {log}. The model itself passed rule check_medaka_model, so look there for another cause (memory, a damaged read file). To run without polishing, set parameters.longread.medaka_model to false." | tee -a {log} >&2
                exit 1
            }}
            """


# ───────────── finalize_assembly — name the finished genome ─────────────
# Why a rule that only copies a file: every mode ends its assembly chain with the same
# file name, contigs_all.fasta, so the mitochondrion split and everything after it read
# one name whichever way the genome arrived. In this mode it is also where the per-sample
# Medaka decision is resolved, through finalize_source() (Takes in, below), and the log of
# each sample names the source.
#
# Takes in:    finalize_source(): MEDAKA_CONSENSUS (rule long_read_consensus above) for a
#              sample in MEDAKA_SAMPLES, otherwise DECONTAM_CONTIGS, the genome that rule
#              select_contigs (shared/10_decontam.smk) wrote after applying the taxonomy
#              selector (every contig, when parameters.decontamination.mode is off, the
#              default). Snakemake calls the function once per sample, with the sample's
#              name, so each job gets its own source file.
# Does:        copies it and records the source in the log.
# Produces:    CONTIGS_ALL, the complete screened and polished genome, mitochondrion
#              included.
# Consumed by: rule separate_mitochondria in shared/15_mitochondria.smk, which splits it
#              into the nuclear genome (contigs_final.fasta) and the mitochondrial
#              contigs.
#
# v1 -> v2: new rule. FunFluxL had no file for the finished genome; a helper function
# picked the consensus or the screened draft for each sample at the moment the
# mitochondrion split asked for its input, and nothing recorded the choice. Now the file
# exists under one name in every mode and the log says where it came from. The copy is
# the same as in contigs/20_assembly.smk.
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
