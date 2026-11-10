# FunFlux v2.1.0 — SignalP 6: which of the predicted proteins are secreted? (optional)
# (rules/shared/45_signalp6.smk)
#
# Many proteins that matter for how a fungus lives — digestive enzymes, toxins, effectors
# that a pathogen delivers to its host, cell-wall proteins — are secreted: they start with a
# short stretch, the signal peptide, that sends them into the secretory pathway, the route
# that secreted proteins take out of the cell. SignalP 6 recognises signal peptides in
# protein sequences. funannotate annotate uses its result to build the secretome, the list
# of secreted proteins of the genome (annotations.secretome.txt).
#
# SignalP 6 is licensed software of the Technical University of Denmark (DTU): the official
# archive is handed out after the licence terms have been accepted on DTU's website. So the
# workflow never downloads it and contains no DTU address. A licensed user downloads the
# archive by hand and names it in files.signalp6_archive; the workflow installs exactly that
# archive into a conda environment of its own.
#
# Everything in this file is defined only when SignalP 6 is switched on
# (parameters.signalp6.mode is fast, slow or slow-sequential; SIGNALP6_ENABLED in 00_common.smk).
# With the mode off (the default) this file defines no rule at all, and funannotate_annotation
# waits for no SignalP result. The three modes, as SignalP's own README describes them:
#   fast             a smaller model that comes close to the full one, with a fraction of the
#                    memory and time
#   slow             the full model, run in parallel; needs more than 14 GB of RAM
#   slow-sequential  the full model, run one part after another: the memory of fast, but
#                    6 times slower
#
# Where it sits, and the chain. SignalP 6 comes after gene prediction and before the final
# annotation; it runs side by side with InterProScan, eggNOG and antiSMASH (40_annotation.smk).
#
#   install_signalp6    : once per run. Installs the archive you supplied and writes
#                         SIGNALP6_MANIFEST (installation.json), the record that says "installed".
#   signalp6_annotation : once per sample. Runs SignalP 6 on the proteins that
#                         funannotate_prediction predicted, and writes SIGNALP6_RESULTS.
#                         funannotate_annotation then hands that table to funannotate annotate.
#
# Names used here, all defined once in 00_common.smk: SIGNALP6_ENABLED, SIGNALP6_ARCHIVE,
# SIGNALP6_ARCHIVE_SHA256, SIGNALP6_MODE, SIGNALP6_TORCH_NUM_THREADS, SIGNALP6_BATCH_SIZE,
# SIGNALP6_MANIFEST, SIGNALP6_RESULTS, SIGNALP6_INSTALL_SCRIPT, FUNANNOTATE_DIR,
# funannotate_proteins and LOGS.
#
# conda: the path is written relative to THIS file (workflow/rules/shared/), so
# "../../envs/signalp6.yaml" lands on workflow/envs/signalp6.yaml. Both rules use that one
# environment. It holds the Python packages SignalP 6 needs (PyTorch, NumPy, ...) but not
# SignalP itself; the installer adds that from the archive.

if SIGNALP6_ENABLED:

    # ── install_signalp6 — install SignalP 6 from the archive you supplied ──
    # Takes in:    nothing from other rules; the official DTU archive, files.signalp6_archive
    #              (SIGNALP6_ARCHIVE; an absolute path, checked when Snakemake reads the workflow).
    #              The settings are the mode (SIGNALP6_MODE: fast, slow or slow-sequential, which
    #              decides which model the archive must contain) and, optionally, the checksum you
    #              expect (parameters.signalp6.archive_sha256).
    # Does:        runs scripts/45_signalp6/install_signalp6.py. In order, the script:
    #                * computes the SHA-256 checksum of the archive, and stops if you gave one
    #                  and it does not match;
    #                * unpacks the archive into a temporary folder inside 00.resources/signalp6/
    #                  (removed at the end), refusing any member that could write outside it
    #                  (absolute or ".." paths, links);
    #                * checks that the archive holds the model for the chosen mode;
    #                * installs the Python package of the archive, and only that (pip with
    #                  --no-deps, no download), into the conda environment of this rule;
    #                * reads the SignalP version (kept for the manifest, not printed), checks that
    #                  the signalp6 command now exists, and runs `signalp6 --help` once;
    #                * moves the unpacked package, which holds the large model weights, to
    #                  00.resources/signalp6/package;
    #                * writes the manifest, first under a temporary name and then renamed, so a
    #                  failed job never leaves half a manifest behind.
    #              matplotlib's cache goes to 00.resources/signalp6/.matplotlib, not to your home
    #              folder. "script:" means that Snakemake runs the Python file itself and gives it
    #              an object called snakemake that holds this rule's input, output and params.
    # Produces:    SIGNALP6_MANIFEST (00.resources/signalp6/installation.json): the archive's name
    #              and checksum, the SignalP version, the mode, where the model folder is (written
    #              relative to the manifest, so the output folder can be moved) and conda_prefix,
    #              the folder of the conda environment the package went into.
    # Consumed by: signalp6_annotation, which reads the model folder from the manifest.
    #
    # Why the manifest records conda_prefix. The package was installed INTO the conda environment,
    # but the manifest sits in the output folder and outlives that environment. Snakemake names an
    # environment's folder after a fingerprint of its recipe, and builds a fresh one when the recipe
    # changes, when --conda-prefix points elsewhere or when .snakemake/conda is deleted. The new
    # environment has no SignalP in it, yet the old manifest is still there, so without a check the
    # annotation rule would fail late. signalp6_annotation therefore checks first, and if SignalP
    # is gone it tells you to delete the manifest, so that this rule runs again in the current
    # environment. conda_prefix says which environment the package was installed into; it is for
    # diagnosis only, and nothing reads it to find files.
    #
    # The log. Snakemake does not copy the output of a `script:` rule into its log file, so the
    # script does it itself: everything it prints, including the output of pip and of
    # `signalp6 --help`, goes to the screen and to this log.
    rule install_signalp6:
        input:
            archive = SIGNALP6_ARCHIVE,
        output:
            manifest = SIGNALP6_MANIFEST,
        params:
            mode = SIGNALP6_MODE,
            archive_sha256 = SIGNALP6_ARCHIVE_SHA256,
        conda:
            "../../envs/signalp6.yaml"
        log:
            LOGS + "/install_signalp6.log"
        script:
            SIGNALP6_INSTALL_SCRIPT


    # ── signalp6_annotation — predict the secreted proteins ──
    # Takes in:    FUNANNOTATE_DIR from rule funannotate_prediction, listed so that this rule waits
    #              for the prediction. The proteins are this file inside that folder:
    #                predict_results/<Genus_species>_<sample>.proteins.fa
    #              (found through funannotate_proteins in 00_common.smk, a function that Snakemake
    #              calls for each job). Also SIGNALP6_MANIFEST from rule install_signalp6, which
    #              says where the SignalP model is.
    # Does:        in this order (the log is emptied first, and every check that fails writes its
    #              reason to the log and stops the job):
    #                0. NEW: checks that SignalP 6 is really installed in the environment this job
    #                   runs in (the signalp6 command exists and the signalp package can be
    #                   imported). If not, it stops at once and names the manifest file to delete.
    #                1. checks that funannotate produced a protein file that is not empty;
    #                2. makes a scratch folder next to the result folder, removed again if anything
    #                   goes wrong, and keeps matplotlib's cache inside it, so that nothing is
    #                   written into your home folder;
    #                3. reads the model folder from the manifest and checks that it exists;
    #                4. writes a copy of the proteins whose FASTA headers hold only the first word,
    #                   the transcript ID. funannotate writes, for example, ">SAMPLE1_000001-T1
    #                   SAMPLE1_000001": the transcript ID, a space, then the gene ID. SignalP copies
    #                   the whole header into the first column of its table, and funannotate takes
    #                   that column, as it is, for the transcript ID, so the column must hold the
    #                   ID alone. The copy is refused if an ID is empty or occurs twice (SignalP
    #                   would quietly rename a repeated ID by adding _2, _3, ...);
    #                5. runs signalp6 on that copy:
    #                     --organism eukarya     only the standard signal peptides of eukaryotes
    #                                            (type Sec/SPI)
    #                     --format none          only the summary files: no file and no plot per
    #                                            protein, which would slow a large job down
    #                     --mode                 fast, slow or slow-sequential, as configured
    #                     --torch_num_threads    the job's threads, for PyTorch (the neural-network
    #                                            library SignalP is written with)
    #                     --bsize                parameters.signalp6.batch_size, the number of
    #                                            proteins the model takes at a time (10 by default)
    #                     --model_dir            the model folder from step 3;
    #                6. checks the result: not empty; holding a line that starts with
    #                   "# SignalP-6" (funannotate reads that line to tell the SignalP version, and
    #                   without it would read the table as an old SignalP 4 one); and at least one
    #                   protein row, every row with at least five tab-separated columns (what
    #                   funannotate reads);
    #                7. only then replaces the result folder by the scratch folder, so that a failed
    #                   or interrupted run cannot leave half a result that looks finished.
    # Produces:    SIGNALP6_RESULTS (prediction_results.txt), in the signalp6 folder of the sample
    #              (04.annotation/{sample}/signalp6/; {sample} is the Snakemake wildcard: the rule
    #              runs once per sample). The folder also holds SignalP's other summary files
    #              (output.gff3, region_output.gff3, processed_entries.fasta, output.json), the
    #              ID-only protein copy and the matplotlib cache.
    # Consumed by: funannotate_annotation, which passes it to funannotate annotate (--signalp).
    #
    # Notes:
    #   * The first step (step 0) checks the installation, because of the rebuilt-environment
    #     problem described under install_signalp6.
    #   * threads: the number of PyTorch threads, SIGNALP6_TORCH_NUM_THREADS, from
    #     parameters.signalp6.torch_num_threads (8 by default). 00_common.smk lowers it to
    #     resources.threads when it is higher, and Snakemake cuts it down to --cores if that is
    #     smaller.
    #   * --write_procs is not passed. It is the number of processes that write the files per
    #     protein, and SignalP 6 reads it only when it writes such files, which --format none
    #     switches off.
    rule signalp6_annotation:
        input:
            funannotate_dir = FUNANNOTATE_DIR,
            installation = SIGNALP6_MANIFEST,
        output:
            signalp_results = SIGNALP6_RESULTS,
        params:
            protein_file = funannotate_proteins,
            mode = SIGNALP6_MODE,
            batch_size = SIGNALP6_BATCH_SIZE,
        threads: SIGNALP6_TORCH_NUM_THREADS
        # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
        resources:
            mem_mb = HEAVY_MEM_MB,
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/signalp6.yaml"
        log:
            LOGS + "/signalp6_annotation_{sample}.log"
        shell:
            """
            set -euo pipefail

            : > "{log}"

            # Stop at once if SignalP 6 is not installed in this environment (see install_signalp6).
            if ! command -v signalp6 > /dev/null 2>&1 || ! python -c "import signalp" >> "{log}" 2>&1; then
                echo "[FunFlux] signalp6_annotation: SignalP 6 is not installed in the conda environment of this rule (the environment was probably rebuilt after SignalP 6 was installed). Delete the file {input.installation} and run the workflow again, so that install_signalp6 installs SignalP 6 into the current environment." | tee -a "{log}" >&2
                exit 1
            fi

            proteins=$(readlink -f "{params.protein_file}")
            if [ ! -s "$proteins" ]; then
                echo "Funannotate did not produce the expected non-empty protein FASTA: $proteins" >> "{log}"
                exit 1
            fi
            output_dir=$(dirname "{output.signalp_results}")
            output_dir=$(readlink -m "$output_dir")
            parent_dir=$(dirname "$output_dir")
            mkdir -p "$parent_dir"
            temporary_dir=$(mktemp -d "$parent_dir/.signalp6-{wildcards.sample}.XXXXXX")
            trap 'rm -rf -- "$temporary_dir"' EXIT
            export MPLCONFIGDIR="$temporary_dir/.matplotlib"
            mkdir -p "$MPLCONFIGDIR"

            installation=$(readlink -f "{input.installation}")
            model_dir=$(python -c 'import json, os, sys; manifest = os.path.abspath(sys.argv[1]); data = json.load(open(manifest)); print(os.path.abspath(os.path.join(os.path.dirname(manifest), data["model_dir"])))' "$installation")
            if [ ! -d "$model_dir" ]; then
                echo "SignalP installation manifest points to a missing model directory: $model_dir" >> "{log}"
                exit 1
            fi

            # funannotate writes each protein header as the transcript ID, a space and the
            # gene ID (>SAMPLE1_000001-T1 SAMPLE1_000001). SignalP copies the whole header into the
            # first column of its table, and funannotate annotate takes that column as the
            # transcript ID. So SignalP gets a copy of the proteins whose headers hold the
            # transcript ID alone, and funannotate can read SignalP's table as it is.
            signalp_proteins="$temporary_dir/proteins.for_signalp.faa"
            if ! awk '
                /^>/ {{
                    identifier = substr($0, 2)
                    sub(/[[:space:]].*$/, "", identifier)
                    if (identifier == "" || seen[identifier]++) exit 1
                    print ">" identifier
                    next
                }}
                {{ print }}
            ' "$proteins" > "$signalp_proteins"; then
                echo "Could not create a unique ID-only protein FASTA for SignalP." >> "{log}"
                exit 1
            fi
            if [ ! -s "$signalp_proteins" ]; then
                echo "The ID-only protein FASTA for SignalP is empty." >> "{log}"
                exit 1
            fi

            echo "Running SignalP 6 from the dedicated workflow environment on ID-normalized Funannotate proteins" >> "{log}"
            signalp6 \
                --fastafile "$signalp_proteins" \
                --organism eukarya \
                --output_dir "$temporary_dir" \
                --format none \
                --mode "{params.mode}" \
                --torch_num_threads {threads} \
                --bsize {params.batch_size} \
                --model_dir "$model_dir" >> "{log}" 2>&1

            results="$temporary_dir/prediction_results.txt"
            if [ ! -s "$results" ]; then
                echo "SignalP 6 did not create a non-empty prediction_results.txt file." >> {log}
                exit 1
            fi
            if ! grep -q '^# SignalP-6' "$results"; then
                echo "SignalP output is not a SignalP 6 prediction summary." >> {log}
                exit 1
            fi
            if ! awk -F '\t' '
                /^#/ {{ next }}
                NF < 5 {{ invalid = 1; next }}
                {{ records += 1 }}
                END {{ exit(invalid || records == 0) }}
            ' "$results"; then
                echo "SignalP 6 summary is empty or does not have the columns Funannotate parses." >> {log}
                exit 1
            fi

            # Only now, with the result checked, replace the sample's signalp6 folder (it holds
            # nothing but this rule's files) by the scratch folder, so that a failed or
            # interrupted run cannot leave half a result that looks finished.
            rm -rf -- "$output_dir"
            mv "$temporary_dir" "$output_dir"
            trap - EXIT
            """
