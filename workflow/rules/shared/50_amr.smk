# FunFlux v2.1.0 — Stage 50: antifungal resistance mutations, with ChroQueTas (optional)
# (rules/shared/50_amr.smk)
#
# A fungus becomes resistant to an antifungal drug or a fungicide mostly through a point
# mutation in the protein the drug acts on, its target: beta-tubulin for the benzimidazoles
# (carbendazim, benomyl), cytochrome b for the QoI fungicides (the strobilurins), Cyp51 for
# the azoles, succinate dehydrogenase for the SDHIs, the Fks glucan synthase for the
# echinocandins, and more. One changed amino acid at the right place, such as E198A in
# beta-tubulin, is enough. Acquired resistance genes, which a bacterial AMR screen looks
# for, play almost no part, so a gene catalogue would find nothing here. What is needed is
# the target protein of the genome, compared position by position with the mutations known
# to cause resistance.
#
# That is what ChroQueTas does (Bédard et al. 2025, Nature Microbiology 10: 2338-2352). Its
# database, the curated part of FungAMR, holds for each species (57 in the version FunFlux
# pins) the reference protein of each target and the single mutations FungAMR records at
# known positions, with the drugs and the evidence for resistance (a positive level) or
# against it (a negative level). For each target of the species, ChroQueTas
#   1. finds the gene in the genome with miniprot, which aligns the reference protein to the
#      DNA and handles the introns, and translates it;
#   2. keeps the hit only if BLASTP finds it at more than 90 % identity over more than
#      75 % of the reference (below that, the hit is more likely a related protein or a
#      broken copy than the target itself);
#   3. aligns it to the reference with MAFFT and reads every FungAMR position: the same
#      amino acid ("No mutation"), the change FungAMR lists ("FungAMR mutation", with the
#      drugs and the evidence) or another change at that position ("New mutation").
# Not detected: changes outside the protein (promoter insertions, overexpression), gene
# duplication or loss, aneuploidy, and resistance that needs several mutations together.
#
# Up to two runs per sample, because of the genetic code. ChroQueTas translates with one
# genetic code for a whole run, and a fungus has two: the nucleus reads the standard code
# (or, in the yeasts of the CTG clade, code 12, which ChroQueTas picks for those species
# itself), and the mitochondrion its own (code 4 in most fungi, where TGA codes for
# tryptophan instead of a stop; code 3 in the Saccharomycetaceae). The one mitochondrially
# encoded target of the database is cytochrome b, the target of the QoI fungicides. So rule
# amr_chroquetas runs ChroQueTas
#   * once on the nuclear genome (FINAL_CONTIGS), with the code ChroQueTas chooses for the
#     species; every target but cytochrome b is read from this run;
#   * once on the mitochondrial genome, with the code of parameters.amr.mito_genetic_code
#     (auto: 3 or 4 by species), only for a species whose database folder holds
#     cytochrome b; cytochrome b is read from this run. Which file holds the mitochondrion
#     depends on the run (00_common.smk, section 9i): MITO_CONTIGS when the screen split it
#     off, otherwise FINAL_CONTIGS, among whose contigs it then lies.
# The wrapper script does the choosing and the reading (chroquetas_amr.py, in
# scripts/50_amr/; its opening text has the details). Each target ends with one of four
# statuses: evaluated (its positions were compared), not evaluated (no position was
# compared; for example, the protein was found with insertions against the reference, and
# a common cause is an intron the protein aligner could not place, such as an intron of a
# mitochondrial cytochrome b gene), not found (not in the genome at more than 90 %
# identity and more than 75 % coverage), or not screened (the species is not in the
# database: the whole sample gets that one row, and a NOTE names the database species of
# the same genus, if any).
#
# The species. ChroQueTas needs the species as its database folder names it
# (Botrytis_cinerea). The wrapper takes column 2 of the annotation table
# (files.annotation_params, through species_name in 00_common.smk) and finds the folder
# itself: the name with underscores for spaces and without "f. sp.", then the first two
# words, then a short list of other names of the same species whose folder carries a
# different name (Candida glabrata -> Nakaseomyces_glabratus, for example). ChroQueTas's
# own species check only looks for the text anywhere in the folder names, so it would
# accept "Botrytis" and fail later; the wrapper's exact match comes first.
#
# Licence, and why the step is off by default. ChroQueTas and its database are licensed
# CC BY-NC-ND 4.0: non-commercial use only, and no modified copies may be passed on. They
# are not on bioconda, and FunFlux does not ship them. Rule chroquetas_db downloads the
# unchanged files at run time from one pinned commit of the authors' GitHub repository
# (links.chroquetas_link) and checks them against two checksums. Everything in this file is
# defined only when parameters.amr.chroquetas is on (AMR_ENABLED in 00_common.smk, section
# 9i); switching it on is the user's acceptance of the licence. With the step off (the
# default) this file defines no rule, no 05.amr folder is made and MultiQC gets no AMR
# table.
#
# Where it sits, and the chain. The AMR screen needs only the finished genome, so its jobs
# can run as soon as the assembly (and the mitochondrion split) is done, side by side with
# QC and annotation. "{sample}" in the file names is a Snakemake wildcard: a rule with it
# runs once for every sample, with the sample name filled in.
#
#   links.chroquetas_link -> chroquetas_db (once per run) -> 00.resources/chroquetas/
#                              ChroQueTas.sh, FungAMR_db/<species>/, species.txt,
#                              PROVENANCE.txt, LICENSE
#   FINAL_CONTIGS (+ MITO_CONTIGS where the screen split it off) + the files above
#       -> amr_chroquetas -> 05.amr/{sample}/chroquetas/{nuclear,mitochondrial}/
#                                            (ChroQueTas's own output, untouched)
#                            05.amr/{sample}/amr_positions.tsv
#                            05.amr/{sample}/amr_summary.tsv      (rule all asks for it)
#                            05.amr/{sample}/amr_summary_mqc.tsv  (rule multiqc reads it)
#
# Names used here, all defined once in 00_common.smk: AMR_ENABLED, CHROQUETAS_URL,
# CHROQUETAS_SCRIPT_SHA256, CHROQUETAS_DB_SHA256, CHROQUETAS_SCRIPT,
# CHROQUETAS_DB_ARCHIVE, CHROQUETAS_DB_DIR, CHROQUETAS_SPECIES_LIST, CHROQUETAS_PROVENANCE,
# CHROQUETAS_LICENSE, FINAL_CONTIGS, AMR_MITO_CONTIGS, AMR_OPTIONAL_INPUTS,
# AMR_MITO_GENETIC_CODE, AMR_CHROQUETAS_DIR, AMR_POSITIONS, AMR_SUMMARY, AMR_SUMMARY_MQC,
# CHROQUETAS_AMR_SCRIPT, LOGS, and the functions species_name() and capped_cpus(). This
# file defines no name of its own.
#
# conda: the path is written relative to THIS file (workflow/rules/shared/), so
# "../../envs/chroquetas.yaml" climbs shared/ -> rules/ -> workflow/ and then goes into
# envs/. The environment holds ChroQueTas's dependencies (miniprot, BLAST, MAFFT, SeqKit,
# EMBOSS) at the versions its authors use, and a pinned Python for the wrapper; it does not
# hold ChroQueTas itself, which rule chroquetas_db downloads. chroquetas_db runs in
# envs/download.yaml for its curl and tar, like rules telomere_db and tiberius_model.
#
# New in v2. There is no v1 counterpart.

if AMR_ENABLED:

    # ── chroquetas_db — download and check ChroQueTas and the FungAMR database ──
    # Takes in:    nothing from other rules; the address links.chroquetas_link
    #              (CHROQUETAS_URL, one pinned commit of the ChroQueTas repository) and the
    #              two checksums parameters.amr.chroquetas_script_sha256 and
    #              chroquetas_db_sha256 (CHROQUETAS_SCRIPT_SHA256, CHROQUETAS_DB_SHA256).
    # Does:        downloads three files with curl: bin/ChroQueTas.sh (the program, one bash
    #              script of about 30 KB), FungAMR_db/db-v20250811.tgz (the database, about
    #              140 KB) and LICENSE. It checks the SHA-256 checksum of the script and of
    #              the archive, and only then unpacks the archive. The address names one
    #              commit, so the files at it should never change; a file that fails the
    #              check is a changed or cut-off download, and the run stops with a message
    #              instead of screening with a program or a database nobody has looked at.
    #              The failed file is renamed with a .failed suffix and kept, so that it can
    #              be looked at (Snakemake deletes the outputs of a failed job, but not a
    #              file under another name). The archive holds one folder, db-v0.6.0/, with
    #              one folder per species inside; tar --strip-components=1 drops the outer
    #              folder, so the species folders land directly in CHROQUETAS_DB_DIR, the
    #              folder ChroQueTas's -f option must name. The archive is then deleted.
    # Produces:    in 00.resources/chroquetas/:
    #                CHROQUETAS_SCRIPT        ChroQueTas.sh, unchanged.
    #                CHROQUETAS_LICENSE       LICENSE, the licence text, unchanged.
    #                CHROQUETAS_DB_DIR        FungAMR_db/, one folder per species, each with
    #                                         a <target>.faa (the reference protein) and a
    #                                         <target>.txt (its resistance positions) per
    #                                         target. Not a declared output (Snakemake checks
    #                                         files, and the species list below stands for
    #                                         it); rule amr_chroquetas gets it through params.
    #                CHROQUETAS_SPECIES_LIST  species.txt, the species folder names, one per
    #                                         line, sorted: the record that the database is
    #                                         unpacked, and a list a user can read.
    #                CHROQUETAS_PROVENANCE    PROVENANCE.txt, one item per line: the three
    #                                         addresses, the date in UTC, the three
    #                                         checksums, the number of species and the
    #                                         licence.
    #              They stay, so a later run into the same output folder does not download
    #              them again.
    # Consumed by: amr_chroquetas, for every sample.
    #
    # Environment: envs/download.yaml, for curl and tar (a minimal Linux installation may
    # lack either); sha256sum, find, sort, chmod, date and wc come from the computer (GNU
    # coreutils and findutils, part of every Linux base system).
    # Where it runs: always on the computer that runs Snakemake itself (localrule: True). On
    # a single computer that changes nothing; on a cluster it keeps the download off the
    # compute nodes, which often have no internet.
    #
    # How the shell works:
    #   * curl -L follows a redirect; --fail turns an HTTP error into a failed command
    #     instead of a saved error page, and a failed command stops the rule (the shell runs
    #     with `set -euo pipefail`).
    #   * sha256sum -c - reads a line "<checksum>  <file>" (two spaces) from its input, here
    #     the checksum from the configuration and the downloaded file, and writes OK or
    #     FAILED for that file. If the check fails, the "then" branch renames the file,
    #     prints a message that says what to do and ends the rule with `exit 1`, so nothing
    #     is unpacked and the species list and the provenance file are never written.
    #   * The database folder is emptied before the archive is unpacked into it, so that a
    #     species dropped from a newer database cannot stay behind from an older one.
    #   * The archive stores its files readable by their owner only; chmod -R a+rX makes
    #     the folders and files readable (and the folders enterable) like the other files
    #     of the run, so that a co-worker who can read the output folder can read these too.
    #     The contents stay as they are.
    #   * find lists the folders directly inside the database folder (-mindepth 1
    #     -maxdepth 1: one level down, not the folder itself; -type d: folders only) and
    #     prints the name of each, one per line (-printf '%f\n'; the \n is written \\n in
    #     the command, because in the rule's command text, a Python string, a single \n
    #     would already be a line break and never reach find); sort puts them in a
    #     fixed order (LC_ALL=C: plain character order, the same on every computer). An
    #     empty list (test -s: "the file has content") means the archive held no species
    #     folder, and the rule stops.
    #   * wc -l counts the lines of the species list for the provenance file and the log.
    rule chroquetas_db:
        localrule: True
        output:
            script = CHROQUETAS_SCRIPT,
            license_text = CHROQUETAS_LICENSE,
            species_list = CHROQUETAS_SPECIES_LIST,
            provenance = CHROQUETAS_PROVENANCE,
        params:
            script_url = CHROQUETAS_URL + "/bin/ChroQueTas.sh",
            database_url = CHROQUETAS_URL + "/FungAMR_db/db-v20250811.tgz",
            license_url = CHROQUETAS_URL + "/LICENSE",
            script_sha256 = CHROQUETAS_SCRIPT_SHA256,
            database_sha256 = CHROQUETAS_DB_SHA256,
            archive = CHROQUETAS_DB_ARCHIVE,
            database_dir = CHROQUETAS_DB_DIR,
        conda:
            "../../envs/download.yaml"
        log:
            LOGS + "/chroquetas_db.log"
        shell:
            """
            set -euo pipefail
            mkdir -p "$(dirname {output.script})"

            echo "[chroquetas_db] downloading {params.script_url}" > {log}
            curl -L --fail -o "{output.script}" "{params.script_url}" >> {log} 2>&1
            echo "[chroquetas_db] downloading {params.database_url}" >> {log}
            curl -L --fail -o "{params.archive}" "{params.database_url}" >> {log} 2>&1
            echo "[chroquetas_db] downloading {params.license_url}" >> {log}
            curl -L --fail -o "{output.license_text}" "{params.license_url}" >> {log} 2>&1

            echo "[chroquetas_db] checking ChroQueTas.sh against parameters.amr.chroquetas_script_sha256 ({params.script_sha256})" >> {log}
            if ! echo "{params.script_sha256}  {output.script}" | sha256sum -c - >> {log} 2>&1; then
                mv "{output.script}" "{output.script}.failed"
                echo "[FunFlux] chroquetas_db: the file downloaded from {params.script_url} does not match parameters.amr.chroquetas_script_sha256, so it was NOT used (details in {log}; the file is kept as {output.script}.failed). The address is pinned to one commit of the ChroQueTas repository, so a mismatch means a changed or cut-off download. Look at the kept file; to accept it, put its checksum (sha256sum {output.script}.failed) into parameters.amr.chroquetas_script_sha256." | tee -a {log} >&2
                exit 1
            fi

            echo "[chroquetas_db] checking the database archive against parameters.amr.chroquetas_db_sha256 ({params.database_sha256})" >> {log}
            if ! echo "{params.database_sha256}  {params.archive}" | sha256sum -c - >> {log} 2>&1; then
                mv "{params.archive}" "{params.archive}.failed"
                echo "[FunFlux] chroquetas_db: the file downloaded from {params.database_url} does not match parameters.amr.chroquetas_db_sha256, so it was NOT unpacked (details in {log}; the file is kept as {params.archive}.failed). The address is pinned to one commit of the ChroQueTas repository, so a mismatch means a changed or cut-off download. Look at the kept file; to accept it, put its checksum (sha256sum {params.archive}.failed) into parameters.amr.chroquetas_db_sha256." | tee -a {log} >&2
                exit 1
            fi

            echo "[chroquetas_db] unpacking into {params.database_dir}" >> {log}
            rm -rf "{params.database_dir}"
            mkdir -p "{params.database_dir}"
            tar -xzf "{params.archive}" --strip-components=1 -C "{params.database_dir}" >> {log} 2>&1
            chmod -R a+rX "{params.database_dir}"
            rm -f "{params.archive}"

            find "{params.database_dir}" -mindepth 1 -maxdepth 1 -type d -printf '%f\\n' | LC_ALL=C sort > "{output.species_list}"
            if [ ! -s "{output.species_list}" ]; then
                echo "[FunFlux] chroquetas_db: the database archive from {params.database_url} holds no species folder (details in {log}). Check links.chroquetas_link: it must name a commit of the ChroQueTas repository whose FungAMR_db/db-v20250811.tgz holds one folder per species." | tee -a {log} >&2
                exit 1
            fi
            species_count=$(wc -l < "{output.species_list}")

            echo "script_url: {params.script_url}" > {output.provenance}
            echo "database_url: {params.database_url}" >> {output.provenance}
            echo "license_url: {params.license_url}" >> {output.provenance}
            echo "downloaded_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> {output.provenance}
            echo "script_sha256: {params.script_sha256}" >> {output.provenance}
            echo "database_sha256: {params.database_sha256}" >> {output.provenance}
            echo "license_sha256: $(sha256sum < "{output.license_text}" | cut -d ' ' -f 1)" >> {output.provenance}
            echo "species: $species_count" >> {output.provenance}
            echo "licence: CC BY-NC-ND 4.0 (Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International; text in LICENSE)" >> {output.provenance}
            echo "[chroquetas_db] done: $species_count species" >> {log}
            """


    # ── amr_chroquetas — the resistance mutations of one sample ──
    # Takes in:    FINAL_CONTIGS, the nuclear genome (rule separate_mitochondria, or rule
    #              filter_contigs in contigs mode without the screen); AMR_OPTIONAL_INPUTS,
    #              which holds MITO_CONTIGS where the screen split the mitochondrion off and
    #              nothing otherwise (00_common.smk, section 9i), so the job waits for it;
    #              CHROQUETAS_SCRIPT and CHROQUETAS_SPECIES_LIST from rule chroquetas_db (the
    #              species list stands for the unpacked database, whose folder,
    #              CHROQUETAS_DB_DIR, comes in through params).
    # Does:        python chroquetas_amr.py (CHROQUETAS_AMR_SCRIPT), once. The script finds
    #              the species' database folder, runs ChroQueTas (as `bash ChroQueTas.sh
    #              ...`) on the nuclear genome and, when the species has cytochrome b, on the
    #              mitochondrial genome with the mitochondrial code, takes each target's
    #              result from the run of its compartment and writes the three tables
    #              (method at the top of this file). Through params come the species
    #              (column 2 of the annotation table), the mitochondrial genetic code (auto,
    #              or a number), the database folder, and the mitochondrial contigs as a
    #              path: in the runs where AMR_MITO_CONTIGS is empty, the :q leaves
    #              --mito-contigs without a value, which the script reads as "not given" and
    #              then takes the mitochondrion from FINAL_CONTIGS. The script judges each
    #              ChroQueTas run by the files it wrote, not by its exit status, which is no
    #              guide (ChroQueTas carries on after a failed tool, and ends with status 1
    #              when it finds no target at all); a run that left neither its summary file
    #              nor the log line that every target was missing stops the job, with the
    #              last lines of ChroQueTas's log. A species not in the database is not an
    #              error: the sample gets one "not screened" row.
    #              Everything the script prints goes to the log; the shell then copies the
    #              lines that start with [FunFlux] (the summary line and, if any, the NOTE)
    #              to the console (grep -F looks for the fixed text; the "|| true" is needed
    #              because grep ends with status 1 when it finds nothing, which under the
    #              strict shell would fail the job).
    # Produces:    AMR_CHROQUETAS_DIR  05.amr/{sample}/chroquetas/, with nuclear/ and
    #                                  mitochondrial/: ChroQueTas's own files, untouched
    #                                  (a directory() output: Snakemake deletes it before
    #                                  the job runs again, so ChroQueTas, which refuses an
    #                                  output folder that exists, always starts clean).
    #              AMR_POSITIONS       05.amr/{sample}/amr_positions.tsv: every position
    #                                  checked (sample, target, compartment, fragment,
    #                                  position, reference_aa, found_aa, call, verdict,
    #                                  resistance_to, no_resistance_to,
    #                                  known_resistance_changes, fungamr_drugs). A change
    #                                  FungAMR lists is not always a resistance change:
    #                                  the verdict says which way the evidence points
    #                                  (resistance, conflicting, unknown, no resistance,
    #                                  wild type, or not compared).
    #              AMR_SUMMARY         05.amr/{sample}/amr_summary.tsv: one row per target
    #                                  (sample, species, target, compartment, status,
    #                                  positions_checked, verdict, fungamr_mutations,
    #                                  new_mutations, resistance_to, no_resistance_to,
    #                                  note).
    #              AMR_SUMMARY_MQC     05.amr/{sample}/amr_summary_mqc.tsv: the summary with
    #                                  MultiQC's custom-content header lines, the rows named
    #                                  "<sample>: <target>", so that the tables of all
    #                                  samples merge into one in the report.
    # Consumed by: rule multiqc (90_report.smk) reads the _mqc copy; rule all asks for the
    #              summary. The positions table and ChroQueTas's folder are for people.
    #
    # Threads: 4, or resources.threads if that is smaller (capped_cpus); ChroQueTas hands
    # them to miniprot and MAFFT. A run takes seconds: about 4 s and 0.6 GB for the nuclear
    # genome of Botrytis cinerea B05.10, about 1 s for its mitochondrion.
    rule amr_chroquetas:
        input:
            genome = FINAL_CONTIGS,
            optional_files = AMR_OPTIONAL_INPUTS,
            chroquetas_script = CHROQUETAS_SCRIPT,
            species_list = CHROQUETAS_SPECIES_LIST,
        output:
            positions = AMR_POSITIONS,
            summary = AMR_SUMMARY,
            mqc = AMR_SUMMARY_MQC,
            chroquetas_dir = directory(AMR_CHROQUETAS_DIR),
        params:
            script = CHROQUETAS_AMR_SCRIPT,
            database = CHROQUETAS_DB_DIR,
            mito_contigs = AMR_MITO_CONTIGS,
            species = species_name,
            mito_genetic_code = AMR_MITO_GENETIC_CODE,
        threads: capped_cpus(4)
        conda:
            "../../envs/chroquetas.yaml"
        log:
            LOGS + "/amr_chroquetas_{sample}.log"
        shell:
            """
            set -euo pipefail

            python {params.script:q} \
                --chroquetas-script {input.chroquetas_script:q} \
                --database {params.database:q} \
                --genome {input.genome:q} \
                --mito-contigs {params.mito_contigs:q} \
                --species {params.species:q} \
                --mito-genetic-code {params.mito_genetic_code:q} \
                --threads {threads} \
                --sample {wildcards.sample:q} \
                --work-dir {output.chroquetas_dir:q} \
                --out-summary {output.summary:q} \
                --out-positions {output.positions:q} \
                --out-mqc {output.mqc:q} > {log} 2>&1

            # Copy the script's summary line and its NOTE, if there is one, to the console.
            grep -F '[FunFlux]' {log} >&2 || true
            """
