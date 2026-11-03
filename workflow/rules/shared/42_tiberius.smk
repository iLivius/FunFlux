# FunFlux v2.0.0 — Tiberius, a gene finder beside those of funannotate (rules/shared/42_tiberius.smk)
#
# Tiberius finds the protein-coding genes of a genome with a deep neural network joined to a
# hidden Markov model (the kind of model AUGUSTUS and GeneMark use to describe the order of
# exons, introns and the DNA between genes). Like the other ab initio gene finders, it reads
# only the DNA sequence. Unlike them it is not trained on the genome being annotated: its
# authors trained one "fungi" model on many fungal genomes, and FunFlux uses that model as it
# is. Tiberius is on by default (parameters.tiberius.mode; see the section "Tiberius" of
# docs/analysis/gene-prediction.md).
#
# How its genes reach the final gene set: funannotate predict hands the models of every gene
# finder, and the protein alignments, to EVidenceModeler (EVM). EVM builds one consensus gene
# per locus: every source votes for the exons, introns and gene ends it supports, each vote
# multiplied by the weight of its source, and the gene structure with the highest total wins.
# Tiberius's models enter as one more source, "other_pred1" (a source of the kind EVM calls
# OTHER_PREDICTION), with weight 1 (TIBERIUS_EVM_WEIGHT): the same as GeneMark-ES, SNAP,
# GlimmerHMM and the protein alignments; AUGUSTUS gets 1, and 2 for its models that the
# evidence supports well. GeneMark-ES is independent of this: with both switched on, both go
# into EVM. The models Tiberius predicts are not the final genes; EVM's choice is.
#
# Where it sits: in the annotation stage, between repeat masking and gene prediction.
# "{sample}" in the file names is a Snakemake wildcard: a rule with it runs once for every
# sample, with the sample name filled in.
#
#   links.tiberius_model_link -> tiberius_model (once per run) -> 00.resources/tiberius/model/
#   CONTIGS_MASK (from repeat_masking) + the model
#       -> tiberius_prediction -> 04.annotation/{sample}/tiberius/tiberius.gtf      (TIBERIUS_GTF)
#       -> tiberius_gff3       -> 04.annotation/{sample}/tiberius/tiberius.evm.gff3 (TIBERIUS_GFF3)
#       -> funannotate_prediction (40_annotation.smk), as --other_gff <file>:1
#
# Only defined when Tiberius is switched on (TIBERIUS_ENABLED in 00_common.smk, section 9e).
# With Tiberius off the three rules do not exist, and funannotate_prediction neither waits for
# nor receives anything from them: tiberius_gate and tiberius_arguments (00_common.smk,
# section 12) then return nothing.
#
# Names used here, all defined in 00_common.smk: TIBERIUS_ENABLED, TIBERIUS_MODEL_URL,
# TIBERIUS_MODEL_SHA256, TIBERIUS_MODEL_ARCHIVE, TIBERIUS_MODEL_DIR, TIBERIUS_MODEL_CONFIG,
# TIBERIUS_MODEL_WEIGHTS, TIBERIUS_SEQ_LEN, CONTIGS_MASK, TIBERIUS_GTF, TIBERIUS_GFF3, CPUS and
# LOGS.
#
# conda: the paths are written relative to THIS file (workflow/rules/shared/), so
# "../../envs/tiberius.yaml" lands on workflow/envs/tiberius.yaml. The model download runs in
# envs/download.yaml, Tiberius in envs/tiberius.yaml and the conversion to GFF3 in
# envs/funannotate.yaml.
#
# New in v2: no v1 workflow ran Tiberius.


if TIBERIUS_ENABLED:

    # ── tiberius_model — download, check and unpack the Tiberius fungi model ──
    # Takes in:    nothing from other rules; the address links.tiberius_model_link
    #              (TIBERIUS_MODEL_URL) and the checksum parameters.tiberius.model_sha256
    #              (TIBERIUS_MODEL_SHA256).
    # Does:        downloads the archive (about 30 MB) from the Tiberius authors' server with
    #              curl, checks its SHA-256 checksum, and only then unpacks it. The address has
    #              no version number, so the authors can replace the file at the same address;
    #              a replaced file fails the check and the run stops with a message, instead of
    #              giving different gene models without a word. The archive holds one folder
    #              (fungi/) with three files: model_config.json (the settings of the network),
    #              model_layers.json (its layers) and weights.h5 (the trained weights).
    #              tar --strip-components=1 drops that folder level, so the three files land
    #              directly in TIBERIUS_MODEL_DIR. The archive is then deleted; after a failed
    #              check it is kept, so that it can be looked at.
    # Produces:    TIBERIUS_MODEL_CONFIG and TIBERIUS_MODEL_WEIGHTS, the two files Snakemake
    #              checks for, in TIBERIUS_MODEL_DIR (00.resources/tiberius/model/). They stay,
    #              so a later run into the same output folder does not download them again.
    # Consumed by: tiberius_prediction, for every sample.
    #
    # Environment: envs/download.yaml, for curl and tar (a minimal Linux installation may
    # lack either), as in rule funannotate_db; sha256sum comes from the computer (GNU
    # coreutils, part of every Linux base system).
    # Where it runs: always on the computer that runs Snakemake itself (localrule: True). On
    # a single computer that changes nothing; on a cluster it keeps the download off the
    # compute nodes, which often have no internet.
    #
    # How the shell works:
    #   * curl -L follows a redirect; --fail turns an HTTP error into a failed command instead
    #     of a saved error page, and a failed command stops the rule (Snakemake runs every
    #     shell block with `set -euo pipefail`).
    #   * sha256sum -c - reads a line "<checksum>  <file>" (two spaces) from its input, here
    #     the checksum from the configuration and the downloaded archive, and writes OK or
    #     FAILED for that file. If the check fails, the "then" branch prints a message that
    #     says what to do and ends the rule with `exit 1`, so tar is never reached.
    rule tiberius_model:
        localrule: True
        output:
            model_config = TIBERIUS_MODEL_CONFIG,
            model_weights = TIBERIUS_MODEL_WEIGHTS,
        params:
            url = TIBERIUS_MODEL_URL,
            sha256 = TIBERIUS_MODEL_SHA256,
            archive = TIBERIUS_MODEL_ARCHIVE,
            model_dir = TIBERIUS_MODEL_DIR,
        conda:
            "../../envs/download.yaml"
        log:
            LOGS + "/tiberius_model.log"
        shell:
            """
            mkdir -p "{params.model_dir}"

            echo "[tiberius_model] downloading {params.url}" > {log}
            curl -L --fail -o "{params.archive}" "{params.url}" >> {log} 2>&1

            echo "[tiberius_model] checking the archive against parameters.tiberius.model_sha256 ({params.sha256})" >> {log}
            if ! echo "{params.sha256}  {params.archive}" | sha256sum -c - >> {log} 2>&1; then
                echo "[FunFlux] tiberius_model: the file downloaded from {params.url} does not match parameters.tiberius.model_sha256, so it was NOT unpacked (details in {log}; the file is kept as {params.archive}). The Tiberius authors may have replaced the model at that address. Either set parameters.tiberius.mode: off, or check the new model and put its checksum (sha256sum {params.archive}) into parameters.tiberius.model_sha256." | tee -a {log} >&2
                exit 1
            fi

            echo "[tiberius_model] unpacking into {params.model_dir}" >> {log}
            tar -xzf "{params.archive}" --strip-components=1 -C "{params.model_dir}" >> {log} 2>&1
            rm -f "{params.archive}"
            echo "[tiberius_model] done: $(ls {params.model_dir} | tr '\\n' ' ')" >> {log}
            """


    # ── tiberius_prediction — Tiberius's gene models of one sample ──
    # Takes in:    CONTIGS_MASK from repeat_masking: the genome that funannotate predict will
    #              receive, so the scaffold names match (funannotate places Tiberius's genes on
    #              its own copy of the genome by those names). The model files from
    #              tiberius_model.
    # Does:        runs Tiberius 1.1.8 on the CPU. The bioconda package has no `tiberius`
    #              command, so the program is started as the Python module tiberius.main
    #              (python -I: Python ignores PYTHONPATH and the user's own packages, so only
    #              the environment's packages load). Its options:
    #                --model DIR       the unpacked fungi model. (Tiberius could also fetch a
    #                                  model itself with --model_cfg, but from the same
    #                                  unversioned address, without a checksum, into the conda
    #                                  environment.)
    #                --no_softmasking  the fungi network reads only the bases (A, C, G, T, N),
    #                                  not which ones are repeats (lowercase). Without this flag
    #                                  Tiberius 1.1.8 hands the network the masking as well, and
    #                                  stops with an input-shape error. So Tiberius does not use
    #                                  the repeat masking of CONTIGS_MASK; the other gene finders do.
    #                --seq_len         TIBERIUS_SEQ_LEN, 100080: the length of the pieces the genome
    #                                  is cut into (00_common.smk, section 9e, says why it is
    #                                  given).
    #                --out FILE        the GTF to write.
    #              Threads: Tiberius has no option for them, so the shell sets them through
    #              environment variables: TF_NUM_INTRAOP_THREADS, the threads TensorFlow (the
    #              library that runs the network) uses inside one operation, and the thread counts
    #              of the maths libraries under it (OMP_, MKL_ and OPENBLAS_NUM_THREADS for OpenMP,
    #              MKL and OpenBLAS) are set to {threads}; TF_NUM_INTEROP_THREADS=2 lets at most two
    #              TensorFlow operations run side by side. CUDA_VISIBLE_DEVICES="" tells TensorFlow
    #              not to look for a graphics card: the environment holds the CPU build of
    #              TensorFlow. The shell holds no comments on purpose: Snakemake compares the
    #              shell text with that of the last run, and a reworded comment there would run
    #              Tiberius, and the gene prediction after it, again for every sample.
    # Produces:    TIBERIUS_GTF, Tiberius's gene models in GTF, its only output format: for each
    #              gene one transcript, with its exons, coding parts (CDS), start and stop codon.
    # Consumed by: tiberius_gff3.
    #
    # Resources: about 16 GB of memory for a fungal genome of about 40 Mb, and roughly half an
    # hour on 24 to 32 CPUs. On a single computer Snakemake starts it by thread count alone (it
    # ignores the memory declared below unless --resources limits it), so switch Tiberius off
    # on a computer without that much free memory. On a cluster the job asks the scheduler for
    # resources.ram_gb of memory, which must then be at least that.
    # It takes all of CPUS, as funannotate_prediction does, because the prediction of the same
    # sample waits for it anyway.
    #
    # Repeatability: TensorFlow's CPU maths may add numbers in a different order from run to run
    # (it warns about "slightly different numerical results" at start-up), which could change a
    # score in its last digits and, rarely, a gene model. Two runs on the same genome, one
    # with 25 threads through this rule and one with 32, gave the same GTF byte for byte.
    rule tiberius_prediction:
        input:
            genome = CONTIGS_MASK,
            model_config = TIBERIUS_MODEL_CONFIG,
            model_weights = TIBERIUS_MODEL_WEIGHTS,
        output:
            gtf = TIBERIUS_GTF,
        params:
            model_dir = TIBERIUS_MODEL_DIR,
            seq_len = TIBERIUS_SEQ_LEN,
        threads: CPUS
        # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
        resources:
            mem_mb = HEAVY_MEM_MB,
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/tiberius.yaml"
        log:
            LOGS + "/tiberius_prediction_{sample}.log"
        shell:
            """
            export TF_NUM_INTRAOP_THREADS={threads}
            export TF_NUM_INTEROP_THREADS=2
            export OMP_NUM_THREADS={threads}
            export MKL_NUM_THREADS={threads}
            export OPENBLAS_NUM_THREADS={threads}
            export CUDA_VISIBLE_DEVICES=""

            python -I -m tiberius.main \
                --genome {input.genome} \
                --model {params.model_dir} \
                --no_softmasking \
                --seq_len {params.seq_len} \
                --out {output.gtf} > {log} 2>&1
            """


    # ── tiberius_gff3 — Tiberius's GTF as the GFF3 that funannotate hands to EVM ──
    # Takes in:    TIBERIUS_GTF from tiberius_prediction.
    # Does:        converts it with augustus_GTF_to_EVM_GFF3.pl, a converter that ships with
    #              EVidenceModeler 1.1.1 inside FunFlux's funannotate environment (activating
    #              that environment sets EVM_HOME, the folder of EVidenceModeler). Tiberius
    #              writes its GTF the way AUGUSTUS does, which is the format this converter
    #              reads. Two general converters, gffread and AGAT, were tried as well; their
    #              files did not pass funannotate's check of --other_gff files.
    #              The converter keeps the genes in a Perl hash (a lookup table), and Perl lists
    #              a hash in an order that changes from run to run. PERL_HASH_SEED=0 and
    #              PERL_PERTURB_KEYS=0 fix that order, so the same GTF always gives the same
    #              GFF3, byte for byte.
    # Produces:    TIBERIUS_GFF3: gene, mRNA, exon and CDS lines for every Tiberius gene.
    # Consumed by: funannotate_prediction, through tiberius_gate and tiberius_arguments in
    #              00_common.smk (--other_gff <TIBERIUS_GFF3>:1).
    rule tiberius_gff3:
        input:
            gtf = TIBERIUS_GTF,
        output:
            gff3 = TIBERIUS_GFF3,
        conda:
            "../../envs/funannotate.yaml"
        log:
            LOGS + "/tiberius_gff3_{sample}.log"
        shell:
            """
            PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0 \
                perl "$EVM_HOME/EvmUtils/misc/augustus_GTF_to_EVM_GFF3.pl" {input.gtf} \
                > {output.gff3} 2> {log}
            """
