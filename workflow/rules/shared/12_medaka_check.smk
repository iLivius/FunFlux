# FunFlux v2.0.0 — Medaka model check (rules/shared/12_medaka_check.smk): make sure the
# Medaka model exists, or work one out from the reads, BEFORE the assembler starts.
#
# Why. Medaka polishes an Oxford Nanopore (ONT) assembly with the reads it was built
# from: it looks at the reads piled up on every position and corrects the systematic
# basecalling errors that the assembler could not. To do that it needs a MODEL that
# matches how the reads were basecalled (flow cell chemistry, basecaller version,
# accuracy mode). Polishing is one of the last steps of a long-read run, because it
# needs the finished, screened assembly. A model name that this Medaka version does
# not have (a typo, or a model dropped in a newer release) would therefore only come
# to light after the assembler had run for an hour or more. This file runs the check
# first, on the raw reads, so that a wrong model stops the run within seconds, and
# the answer is written down once so that the polishing step cannot disagree with it.
#
# One rule, check_medaka_model, once per ONT sample that uses Medaka. The chain:
#
#   raw ONT FASTQ (the user's file, raw_long_fastq)
#     --check_medaka_model--> {sample}_medaka_model.txt  (MEDAKA_MODEL_RESOLVED)
#       --> rule long_read_assembly lists it as an input (medaka_gate in 00_common.smk),
#           so the assembler of an ONT sample waits for the check and never starts on
#           a model that would fail at the end
#       --> rule long_read_consensus (the Medaka step of the longread and hybrid front
#           ends) reads the one line in it with cat and passes it to medaka_consensus -m
#
# What the script does with parameters.<mode>.medaka_model (00_common.smk reads it as
# MEDAKA_MODEL; docs/methods/medaka-model.md explains the choices):
#   a model name   : checked against the models this Medaka version knows
#                    (medaka tools list_models). Unknown: the run stops, and the log
#                    holds a table of the closest available names. Known: its weights
#                    are downloaded now if Medaka does not hold them yet (see "Where it
#                    runs" below). A path to a model file is accepted as it is.
#   auto, or empty : the model is worked out from the basecaller tag in the headers of
#                    the first reads (medaka tools resolve_model --auto_model consensus),
#                    which also downloads its weights if Medaka does not hold them yet.
#   false (or an off word) : Medaka is off; MEDAKA_SAMPLES is then empty and this file
#                    defines no rule at all.
#   medaka_model_fallback_auto: true : a named model that does not exist falls back to
#                    the auto answer, with a loud warning, instead of stopping the run.
#
# The Snakefile includes this file with the other shared files. The rule exists only in
# longread and hybrid mode AND only when at least one sample needs Medaka: the `if`
# below tests HAS_LONG_READS and MEDAKA_SAMPLES (the ONT samples of the batch when
# Medaka is on; PacBio HiFi reads are never polished with Medaka). It is decided once,
# when Snakemake reads the workflow, before any job runs. A PacBio-only batch, or a
# batch with Medaka switched off, therefore has no check rule and no model file.
#
# Everything referenced here is defined once in 00_common.smk: HAS_LONG_READS,
# MEDAKA_SAMPLES, MEDAKA_MODEL, MEDAKA_FALLBACK_AUTO, MEDAKA_MODEL_RESOLVED,
# MEDAKA_CHECK_SCRIPT and LOGS, and the functions raw_long_fastq() (the sample's raw
# long-read file) and has_explicit_medaka_model() (True only for a concrete model name).
#
# conda: environment paths are written relative to THIS file (workflow/rules/shared/),
# so "../../envs/medaka.yaml" goes up two folders, to workflow/, and then into
# workflow/envs/. The script needs the `medaka` program on its PATH, which is why it
# runs inside the Medaka environment although it is plain Python.
#
# CPUs: the check reads a few read headers and asks Medaka for its model list; it needs
# one thread and declares none.
#
# Where it runs: always on the computer that runs Snakemake itself (localrule: True). Medaka
# 2.2.2 ships the weights of only four models; for any other model it downloads them from
# GitHub the first time it needs them, and stores them in its conda environment (or in
# ~/.medaka/data when that folder cannot be written). The check makes that first time happen
# here, where there is internet, so the polishing job (rule long_read_consensus), which on a
# cluster runs on a compute node that may have no internet, finds the weights stored. On a
# single computer this changes nothing. The script explains the download.


# The rule below exists only when a sample of this run is polished with Medaka.
if HAS_LONG_READS and MEDAKA_SAMPLES:

    # ── The three values the rule hands to the script, worked out once ──
    # The --model value: the model name as the user wrote it, or EMPTY for every setting
    # that means "work it out" (auto, true, an empty value). has_explicit_medaka_model()
    # makes that distinction; MEDAKA_MODEL itself holds the raw setting, which is the
    # text "auto" when the key is missing, so the raw value cannot be passed as it is.
    if has_explicit_medaka_model(MEDAKA_MODEL):
        _MEDAKA_MODEL_ARG = str(MEDAKA_MODEL).strip()
    else:
        _MEDAKA_MODEL_ARG = ""

    # The --fallback-to-auto value: the script wants the words true or false.
    if MEDAKA_FALLBACK_AUTO:
        _MEDAKA_FALLBACK_ARG = "true"
    else:
        _MEDAKA_FALLBACK_ARG = "false"

    # The samples this rule may run for: only the ONT samples that use Medaka, written
    # as a regular expression (the names joined by |, "or"; re.escape() makes a dot in
    # a name match only a dot). The rule-level wildcard constraint below narrows the
    # workflow-wide one (00_common.smk, section 14), so a PacBio sample's model file can
    # never be asked for, not even by mistake.
    _MEDAKA_SAMPLE_CONSTRAINT = "|".join(re.escape(sample) for sample in MEDAKA_SAMPLES)

    # ───────────── check_medaka_model — check or work out the Medaka model ─────────────
    # Takes in:    the sample's RAW long-read FASTQ (raw_long_fastq; plain or .gz), the
    #              file the user gave in input.longread_dir. The raw file, not the Filtlong
    #              output: in hybrid mode with short_read_guidance: true the filtered file
    #              exists only hours into the run (guided Filtlong waits for the Illumina
    #              guide assembly), and this check has to
    #              answer before the assembler starts. Only the basecaller tag in the read
    #              headers is read, which filtering does not change.
    # Does:        python medaka_model_check.py (MEDAKA_CHECK_SCRIPT) with the model
    #              setting, the reads, the fallback switch and the output path; see the
    #              head of this file for what the script does with each setting. The
    #              setting reaches the script through `params:`, so a changed
    #              medaka_model re-runs the check (and, through the model file, the
    #              polishing). The `|| { cat ...; exit 1; }` tail runs only when the
    #              script fails: it copies the log, with the script's table of suggested
    #              models, to the console, so the fix is visible without opening the log.
    # Produces:    MEDAKA_MODEL_RESOLVED, one line: the model name, or the path of the
    #              model file that Medaka's auto mode found. Kept (it is small and it
    #              documents which model polished the sample).
    # Consumed by: rule long_read_assembly (as a gate: it waits for the file) and rule
    #              long_read_consensus (reads the line and passes it to Medaka).
    #
    # --model is written with an explicit pair of double quotes so that the script
    # always receives a value, an empty one in auto mode; the usual {...:q} quoting is
    # kept for the paths, which never are empty.
    #
    # v1 -> v2: new rule. FunFluxL checked a named model inside the Medaka rule itself,
    # after the assembly had run, and in auto mode let medaka_consensus work the model
    # out on its own, so a mistake surfaced hours late and the choice was never written
    # down.
    rule check_medaka_model:
        localrule: True
        input:
            reads = raw_long_fastq,
        output:
            resolved = MEDAKA_MODEL_RESOLVED,
        params:
            script = MEDAKA_CHECK_SCRIPT,
            model = _MEDAKA_MODEL_ARG,
            fallback = _MEDAKA_FALLBACK_ARG,
        wildcard_constraints:
            sample = _MEDAKA_SAMPLE_CONSTRAINT,
        conda:
            "../../envs/medaka.yaml"
        log:
            LOGS + "/check_medaka_model_{sample}.log"
        shell:
            """
            python {params.script:q} \
              --model "{params.model}" \
              --reads {input.reads:q} \
              --fallback-to-auto {params.fallback} \
              --out {output.resolved:q} > {log} 2>&1 || {{ cat {log} >&2; exit 1; }}
            """
