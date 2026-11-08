# FunFlux v2.0.0 — stage 04, annotation (rules/shared/40_annotation.smk): from the finished
# nuclear genome to a functionally annotated one.
#
# The genome that arrives here (FINAL_CONTIGS, 02.assembly/{sample}/contigs_final.fasta) is
# the NUCLEAR genome. The earlier stages have already removed contaminants and split the
# mitochondrion off, because a fungal mitochondrion uses a different genetic code (code 4)
# and must not go through a nuclear gene predictor. In this stage that genome is
#   1. tidied      : very short contigs, and contigs that are near-copies of part of another
#                    contig, are dropped; the rest are sorted by length and renamed scaffold_1,
#                    scaffold_2, ... (funannotate_preprocess);
#   2. soft-masked : repeats are written in lowercase, so that repeat sequence does not
#                    give rise to false genes (repeat_masking);
#   3. searched for genes (funannotate_prediction): gene finders that learn from the genome
#                    itself, plus protein evidence from a related fungus and from UniProt;
#   4. described   : the predicted proteins are scanned for protein domains (InterProScan) and
#                    assigned to orthologous groups and functions (eggNOG-mapper), the genes
#                    are scanned for biosynthetic gene clusters of secondary metabolites
#                    (antiSMASH) and, if switched on, the proteins for secretion signals
#                    (SignalP 6, in shared/45_signalp6.smk);
#   5. merged      : funannotate annotate puts all of this into the final GenBank, GFF3 and
#                    table files (funannotate_annotation), and the scratch files are removed
#                    (cleanup).
#
# Nothing in this file depends on the mode (illumina, longread, hybrid or contigs): every
# mode ends its assembly stages with the same FINAL_CONTIGS file, so all four modes run
# exactly these rules.
#
# The chain for one sample ({sample} is a Snakemake "wildcard": every rule below runs once
# for each sample name, with the name filled in):
#
#   FINAL_CONTIGS
#     -> funannotate_preprocess : CONTIGS_CLEAN, CONTIGS_SORT and the record CONTIGS_REMOVED
#     -> repeat_masking         : CONTIGS_MASK
#     -> funannotate_prediction : FUNANNOTATE_DIR (genes, proteins, GenBank file). Before it starts:
#                                 augustus_model_check (below) records whether AUGUSTUS uses a
#                                 stored model or is trained on the genome, and stops on a model name
#                                 the database does not hold; when GeneMark is switched on,
#                                 genemark_check (shared/41_genemark.smk) checks that the GeneMark
#                                 licence key still works; when Tiberius is switched on (the
#                                 default), tiberius_prediction and tiberius_gff3
#                                 (shared/42_tiberius.smk) predict genes with Tiberius, which the
#                                 prediction then adds to its evidence.
#     -> four analyses that all read the prediction and can run side by side:
#          iprscan_annotation   : IPRSCAN_XML         protein domains (InterProScan)
#          eggnog_annotation    : EGGNOG_DIR          orthologous groups and functions (eggNOG-mapper)
#          antismash_annotation : ANTISMASH_DIR       secondary-metabolite gene clusters (antiSMASH)
#          signalp6_annotation  : SIGNALP6_RESULTS    secreted proteins (SignalP 6; it lives in
#                                                     shared/45_signalp6.smk and exists only when
#                                                     SignalP 6 is switched on)
#     -> funannotate_annotation : ANNOTATION_COMPLETE (the merged annotation, in FUNANNOTATE_DIR)
#     -> cleanup                : CLEANUP_COMPLETE
#
# No later stage reads the files of this one. rule all (workflow/Snakefile) asks for the
# main ones directly, through all_targets() in 00_common.smk, section 15, with
# CLEANUP_COMPLETE last. The MultiQC report (90_report.smk) uses only the QC results of the
# earlier stages.
#
# Once per run, not per sample (the reference databases the chain needs):
#
#   funannotate_db            : downloads, verifies and unpacks the funannotate database, if
#                               directories.funannotate_db does not hold one yet.
#   normalize_funannotate_db  : repairs the absolute paths stored inside that database
#                               (FUNANNOTATE_DB_MARKER says it is ready for use).
#   antismash_db              : downloads the antiSMASH database into ANTISMASH_DB_DIR, or
#   antismash_db_local        : builds ANTISMASH_DB_DIR as a set of links to a copy you already
#                               hold (directories.antismash_db). Only one of the two exists in a run.
#
# Every name used here is defined once in 00_common.smk, which workflow/Snakefile includes
# before this file; no path, stage number or setting is worked out again here:
#   paths    : FINAL_CONTIGS, FUNANNOTATE_DB_MARKER, ANNOTATION_DIR, CONTIGS_CLEAN, CONTIGS_SORT,
#              CONTIGS_MASK, CONTIGS_REMOVED, REPEATMASK_DIR, FUNANNOTATE_DIR, IPRSCAN_DIR,
#              IPRSCAN_XML, EGGNOG_DIR, ANTISMASH_DIR, ANTISMASH_DB_DIR, SIGNALP6_RESULTS,
#              ANNOTATION_COMPLETE, CLEANUP_COMPLETE, AUGUSTUS_MODEL_RECORD, LOGS, NORMALIZE_DB_SCRIPT,
#              AUGUSTUS_CHECK_SCRIPT
#   settings : FUNANNOTATEDB, FUNANNOTATEID, FUNANNOTATEURL, FUNANNOTATESHAURL, EGGNOGDB,
#              ANTISMASHDB, IPRSCAN, MASKING_METHOD, REPEATMODELER_QUICK, REPEATMODELER_LTRSTRUCT,
#              REPEATMODELER_SEED, PREDICTION_BUSCO_DB, GENEMARK_ENABLED, GENEMARKDIR,
#              SIGNALP6_ENABLED, CPUS, and capped_cpus(N), which is the smaller of N and CPUS
#   functions: species_name, protein_evidence, augustus_model, augustus_optimize_for_sample,
#              augustus_arguments, funannotate_proteins, funannotate_gbk, genemark_gate,
#              tiberius_gate, tiberius_arguments. Each is
#              a small "named function": Snakemake calls
#              it once per job, hands it that job's wildcards (so it knows the sample), and uses
#              what it returns, for example the species name of that sample.
# The only names defined in this file are _FUNANNOTATE_DB_INFO (the leading underscore marks a
# name that belongs to this file alone) and the two small functions signalp6_results_if_enabled
# and signalp6_option.
#
# conda: paths are written relative to THIS file. Snakemake resolves a rule's conda
# environment against the file that defines the rule, and this one sits in
# workflow/rules/shared/, so "../../envs/x.yaml" climbs shared/ -> rules/ -> workflow/ and
# lands on workflow/envs/x.yaml. funannotate_db runs in envs/download.yaml, which supplies curl
# and tar (with gzip); its sha256sum comes from the computer. Two rules use no conda
# environment: antismash_db_local and cleanup need only standard shell commands (mkdir, ln,
# rm). normalize_funannotate_db and augustus_model_check
# run a Python script that uses only the standard library, in envs/python.yaml, which holds a
# pinned Python 3.12 and nothing else: activating that environment makes `python` that Python,
# instead of whatever `python` the computer has (a computer with only `python3`, or an old
# Python, would stop the job).
#
# Threads: every rule that tells its tool how many CPUs to use declares Snakemake's built-in
# `threads:` and writes {threads} in the shell. Most declare capped_cpus(N), where N is the
# number of threads the tool can use well and CPUS is the whole budget set in the config;
# funannotate_prediction takes all of CPUS. Snakemake counts declared threads against the
# number given with --cores: it never runs jobs that together ask for more, and it cuts a
# single larger request down to that number. (A custom `resources: cpus`, as v1 used, is not
# counted.)
#
# Shell blocks are Snakemake format strings: {name} is replaced by Snakemake, so a literal
# curly brace for awk or bash has to be written doubled, {{ like this }}. Snakemake starts
# every shell block with `set -euo pipefail`, so a command that fails stops the job, unless
# the shell tests for the failure itself (with if or ||).


# The index file of a funannotate database. Its presence is how the workflow knows that a
# database is in place: funannotate_db produces it, and normalize_funannotate_db reads it.
# Defined once here because both rules name it. (00_common.smk, section 9d, looks for the
# same file to decide whether the database is already on disk.)
_FUNANNOTATE_DB_INFO = os.path.join(FUNANNOTATEDB, "funannotate-db-info.txt")


# ──────────── funannotate_db — download, verify and unpack the database ────────────
# funannotate needs a large reference database (UniProt, Pfam, dbCAN, MEROPS, MIBiG,
# InterPro, and the trained gene-finder settings of about 166 fungal species). A frozen
# snapshot of it is published on Zenodo, together with a checksum file, and
# directories.funannotate_db is where it lives. The folder is shared between runs, so this
# rule runs only when it holds no database yet (Snakemake asks for funannotate-db-info.txt,
# the index file of a real snapshot). On a machine that already has the database it does not
# run, except after a change to its shell (see the last bullet of "How the shell works").
#
# Takes in:    nothing from the workflow. It fetches the archive (FUNANNOTATEURL, from
#              links.funannotate_link) and its SHA-256 checksum file (FUNANNOTATESHAURL, the
#              .sha256 file of the same Zenodo record) unless they are already in the folder.
# Does:        downloads what is missing (both files when the checksum file is not there,
#              otherwise only the archive, if that is missing), CHECKS the archive against the
#              checksum, and only then unpacks it into the database folder. A SHA-256 checksum
#              is a fingerprint of a file: if even one byte of the archive differs from the
#              published one, the fingerprint differs.
# Produces:    the unpacked database. Its index file, funannotate-db-info.txt, is the file
#              Snakemake watches. The downloaded .tar.gz (about 2.1 GB for the October 2026
#              snapshot) and .sha256 stay in the folder.
# Consumed by: normalize_funannotate_db, which repairs the paths inside the unpacked database.
#
# Environment: envs/download.yaml, for curl and tar (a minimal Linux installation may lack
# either); sha256sum comes from the computer (GNU coreutils, part of every Linux base system).
# Where it runs: always on the computer that runs Snakemake itself (localrule: True). On a
# single computer that changes nothing; on a cluster it keeps the download off the compute
# nodes, which often have no internet. The archive is about 2.1 GB, and unpacking it takes a
# few minutes of that computer's time.
#
# v1 -> v2:
#   * NEW: the archive is checked against its checksum BEFORE it is unpacked. v1 downloaded
#     the checksum file and then never looked at it, so a truncated or damaged archive was
#     unpacked without a word. Now `sha256sum -c` runs first; on a mismatch the rule stops
#     with a message, and nothing is unpacked.
#   * The log is LOGS/funannotate_db.log. The v1 message: and priority: lines are not carried over.
#   * Every step now writes a line to the log that says what it found and what it does (v1
#     wrote "missing, downloading tar.gz only" without saying what was missing). The first
#     line starts the log with > ; every later line is added with >> (v1 started it with > in
#     one branch only).
#
# How the shell works:
#   * curl -L follows a redirect if the server answers with one. --fail turns an HTTP error (a
#     missing file, for example) into a failed command instead of a saved error page, and a
#     failed command stops the rule (see `set -euo pipefail` in the file header).
#   * The checksum file names the archive without any folder ("<hash>  <name>.tar.gz"), so the
#     check has to run from inside the database folder. The round brackets ( ... ) run those
#     two commands in a sub-shell, a private copy of the shell, so the cd does not carry over
#     to the lines after it. If the check fails, the "then" branch prints a message that names
#     both files and ends the rule with `exit 1`, so tar is never reached and the rule counts
#     as failed.
#   * The archive holds one top folder (funannotate_db_2026.10.08/ in the October 2026
#     snapshot); tar --strip-components=1 drops it, so the files land directly in the
#     database folder.
#   * Change this shell, or its environment (envs/download.yaml and its pin list), only when
#     needed. Snakemake re-runs a rule whose shell text or environment changed (comment lines
#     included) when its records, in the .snakemake folder of the folder it was started from,
#     show that this rule made the output. A re-run of this one unpacks the archive again over
#     the shared, already repaired database, and every step after it runs again too.
rule funannotate_db:
    localrule: True
    output:
        db_info = _FUNANNOTATE_DB_INFO,
    params:
        funannotate_tar_url = FUNANNOTATEURL,
        funannotate_sha_url = FUNANNOTATESHAURL,
        funannotate_tar = os.path.join(FUNANNOTATEDB, f"{FUNANNOTATEID}.tar.gz"),
        funannotate_sha = os.path.join(FUNANNOTATEDB, f"{FUNANNOTATEID}.sha256"),
        funannotate_sha_name = f"{FUNANNOTATEID}.sha256",
        funannotate_db = FUNANNOTATEDB,
    conda:
        "../../envs/download.yaml"
    log:
        LOGS + "/funannotate_db.log"
    shell:
        """
        mkdir -p "{params.funannotate_db}"

        # The first line starts the log (>); every later line is added to it (>>).
        echo "[funannotate_db] database folder: {params.funannotate_db}" > {log}

        # No checksum file yet: download the checksum file and the archive, in that order.
        if [ ! -f "{params.funannotate_sha}" ]; then
            echo "[funannotate_db] checksum file {params.funannotate_sha} not found: downloading it from {params.funannotate_sha_url}" >> {log}
            curl -L --fail -o "{params.funannotate_sha}" "{params.funannotate_sha_url}" >> {log} 2>&1

            echo "[funannotate_db] downloading the archive {params.funannotate_tar} from {params.funannotate_tar_url}" >> {log}
            curl -L --fail -o "{params.funannotate_tar}" "{params.funannotate_tar_url}" >> {log} 2>&1
        else
            # The checksum file is already there, so it is not downloaded again. The archive is
            # downloaded only if it is missing, for example when an earlier run stopped after
            # the checksum file and before any of the archive arrived. An archive cut off
            # half-way is not missing: it fails the checksum test below and stops the rule.
            echo "[funannotate_db] checksum file {params.funannotate_sha} found: not downloaded again" >> {log}

            if [ ! -f "{params.funannotate_tar}" ]; then
                echo "[funannotate_db] archive {params.funannotate_tar} not found: downloading it from {params.funannotate_tar_url}" >> {log}
                curl -L --fail -o "{params.funannotate_tar}" "{params.funannotate_tar_url}" >> {log} 2>&1
            else
                echo "[funannotate_db] archive {params.funannotate_tar} found: not downloaded again" >> {log}
            fi
        fi

        # Check the archive against the checksum file before anything is unpacked. sha256sum
        # writes its verdict to the log: the archive name followed by OK, or by FAILED and a
        # warning line.
        echo "[funannotate_db] checking the archive against the checksum file" >> {log}
        if ! ( cd "{params.funannotate_db}" && sha256sum -c "{params.funannotate_sha_name}" ) >> {log} 2>&1; then
            echo "[FunFlux] funannotate_db: {params.funannotate_tar} does not match its checksum file {params.funannotate_sha}, so it was NOT unpacked (details in {log}). Delete both files and run the workflow again to download them afresh." | tee -a {log} >&2
            exit 1
        fi

        # The archive is intact: unpack it into the database folder.
        echo "[funannotate_db] unpacking the archive into {params.funannotate_db}" >> {log}
        tar -xzf "{params.funannotate_tar}" --strip-components=1 -C "{params.funannotate_db}" >> {log} 2>&1
        echo "[funannotate_db] done" >> {log}
        """


# ──────── normalize_funannotate_db — repair the paths inside the database ────────
# The funannotate snapshot was built on someone else's computer and remembers where its
# files were, as absolute paths, in two kinds of small text files: the database index
# (funannotate-db-info.txt) and one info.json per trained species. On this computer those
# paths point nowhere, and funannotate cannot find its own files. This rule rewrites them so
# that they point into directories.funannotate_db.
#
# Takes in:    funannotate-db-info.txt in the database folder: made by funannotate_db, or
#              already there from an earlier run.
# Does:        runs workflow/scripts/40_annotation/normalize_funannotate_db_paths.py on the
#              database folder. The script touches only those small metadata files (the big
#              data files are only checked for existence), leaves a file alone when its paths
#              are already right, and stops with an error if a corrected path would point at a
#              file that is not there, instead of writing a path that cannot work. Running it
#              twice changes nothing the second time.
# Produces:    FUNANNOTATE_DB_MARKER, a small JSON file with the counts (how many files were
#              examined and how many paths were rewritten) and the database path. Snakemake
#              needs only its existence: it means "the database is ready to use".
# Consumed by: funannotate_preprocess, augustus_model_check and funannotate_prediction, which
#              list the marker as an input so that none can start before the database has
#              been repaired.
#
# CAUTION: this rule rewrites files INSIDE the shared database folder, in place, and keeps no
# copy of the old content. There is no lock, so two runs that both start on a database that
# has not been repaired yet can write the same temporary file together and damage it. Start
# one run, let it get past this step, and only then start another one that shares the
# database. Once the database is repaired this step only reads it.
#
# v1 -> v2:
#   * New rule. The unpublished FunFluxL+ did this in a Python block inside its Snakefile (a
#     run: block); here it is a script with its own tests, and the logic is unchanged. The
#     rule runs it with the pinned Python of envs/python.yaml (the script needs nothing beyond
#     the Python standard library).
rule normalize_funannotate_db:
    input:
        db_info = _FUNANNOTATE_DB_INFO,
    output:
        marker = FUNANNOTATE_DB_MARKER,
    params:
        funannotate_db = FUNANNOTATEDB,
    conda:
        "../../envs/python.yaml"
    log:
        LOGS + "/normalize_funannotate_db.log"
    shell:
        # {NORMALIZE_DB_SCRIPT:q} inserts the path of the script, quoted for the shell. When the
        # script succeeds it writes its own one-line summary into the log. When it fails, 2>> {log}
        # adds its error message to the log (the script then writes neither the marker nor the
        # summary line, so Snakemake sees the rule fail).
        """
        python {NORMALIZE_DB_SCRIPT:q} \
            --database "{params.funannotate_db}" \
            --summary-json {output.marker} \
            --log {log} 2>> {log}
        """


# ──────────── funannotate_preprocess — tidy and sort the genome ────────────
# Two funannotate steps run back to back on the nuclear genome:
#   funannotate clean : drops contigs shorter than 500 bp (funannotate's default minimum);
#                       contigs that hold any character other than the DNA letters A, C, G, T,
#                       N and the ambiguity codes R, Y, S, W, K, M, B, D, H, V (upper or lower
#                       case); contigs with fewer than four different letters (a run of N, for
#                       example); and contigs that are near-copies of another contig (more than
#                       95% identical over more than 95% of their length, found with minimap2;
#                       only contigs shorter than the assembly's N50 are tested). Such
#                       duplicates are typical of a provider's assembly that still carries
#                       haplotigs (the second copy of a heterozygous region). --cpus sets how
#                       many minimap2 searches run at the same time.
#   funannotate sort  : sorts the remaining contigs from longest to shortest and renames them
#                       scaffold_1, scaffold_2, ... because AUGUSTUS can have problems with
#                       complicated contig names. --minlen 0 keeps everything clean kept, and
#                       it is needed: without --minlen, funannotate sort 1.8.17 stops with a
#                       Python error (TypeError).
# A third step is new in v2 and only reads: it lists the contigs that clean dropped.
#
# Takes in:    FINAL_CONTIGS, the nuclear genome. In contigs mode without the screen,
#              filter_contigs writes it directly; with the screen, separate_mitochondria
#              (shared/15_mitochondria.smk) splits it off the complete assembly; the read-based
#              modes end the same way. Also FUNANNOTATE_DB_MARKER from normalize_funannotate_db.
#              The marker is only waited for, never read: funannotate takes its files from the
#              database folder, and that folder must have been repaired first.
# Does:        funannotate clean, funannotate sort, then compares the contig names before and
#              after clean (the first word of each FASTA header is the name).
# Produces:    CONTIGS_CLEAN   the genome after clean, with the original contig names;
#              CONTIGS_SORT    the genome that is annotated from here on, sorted and renamed
#                              scaffold_N (the QUAST, BUSCO and ITS results of the earlier stages
#                              use the original names; the annotated_as column of the contig
#                              summary, written by rule contig_summary in 25_contig_summary.smk,
#                              links the two names);
#              CONTIGS_REMOVED one contig name per line: every contig present in FINAL_CONTIGS but
#                              absent from CONTIGS_CLEAN, in the order of the input FASTA. An EMPTY
#                              file means clean dropped nothing, so a contig never disappears silently.
# Consumed by: repeat_masking reads CONTIGS_SORT. CONTIGS_CLEAN and CONTIGS_REMOVED are records
#              for the reader; no rule reads them.
#
# v1 -> v2:
#   * NEW output CONTIGS_REMOVED (see above).
#   * The input is FINAL_CONTIGS, the nuclear genome (v1 FunFlux read contigs_sel.fasta, which
#     still held the mitochondrial genome).
#   * threads: capped_cpus(2), passed to funannotate clean as --cpus {threads}. v1 declared
#     `resources: cpus = CPUS` but passed no CPU number, so funannotate clean used its built-in
#     default, --cpus 2 (two minimap2 searches at a time); v2 passes the same 2 (1 when
#     resources.threads is 1), and Snakemake now counts the job as 2 threads, not 1. Each
#     minimap2 search can use more than one CPU (minimap2 starts 3 threads unless told
#     otherwise), so the job uses a little more than it declares: somewhat more than 2 CPUs
#     on average (the figure varies with the load on the
#     computer). The cleaned genome does not depend on --cpus: each contig is tested on its
#     own, and the output is written in the order of the input. funannotate sort has no CPU
#     option and runs on one CPU.
#
# Where it runs: funannotate clean works in a scratch folder called clean_<random ID>, which
# it creates in the CURRENT directory and deletes when it finishes. In v2 the current directory
# would be the folder Snakemake was started from (v1 changed into the output folder first), so
# a clean that is killed half-way would leave that folder in your project. The shell therefore
# changes into the sample's annotation folder first (every path it uses is absolute), where the
# clean_* line of the cleanup rule looks for leftovers. This is new in v2.
#
# How the list is made: awk reads two files in one go. While it reads the first one, the
# cleaned genome (FILENAME == ARGV[1] means "the file being read is the first one on the
# command line"), it remembers the name of every contig in the array `kept`. For the second
# file, the input genome, it prints every name that is not in `kept`. substr($1, 2) is the
# first word of the header without the ">". Comparing the file name, not the line number, is
# what keeps this correct even if the cleaned file were empty.
rule funannotate_preprocess:
    input:
        contigs = FINAL_CONTIGS,
        db_marker = FUNANNOTATE_DB_MARKER,
    output:
        contigs_clean = CONTIGS_CLEAN,
        contigs_sort = CONTIGS_SORT,
        contigs_removed = CONTIGS_REMOVED,
    params:
        funannotate_db = FUNANNOTATEDB,
        annotation_dir = ANNOTATION_DIR,
    threads: capped_cpus(2)
    conda:
        "../../envs/funannotate.yaml"
    log:
        LOGS + "/funannotate_preprocess_{sample}.log"
    shell:
        """
        set -euo pipefail

        # funannotate finds its database through this variable.
        export FUNANNOTATE_DB="{params.funannotate_db}"

        # funannotate clean makes a scratch folder in the current directory (see the comment above).
        cd "{params.annotation_dir}"

        # Drop short, non-nucleotide and duplicated contigs; the contig names stay as they were.
        funannotate clean \
            -i {input.contigs} \
            -o {output.contigs_clean} \
            --cpus {threads} > {log} 2>&1

        # Sort the contigs from longest to shortest and rename them scaffold_1, scaffold_2, ...
        funannotate sort \
            -i {output.contigs_clean} \
            -o {output.contigs_sort} \
            --minlen 0 >> {log} 2>&1

        # List the contigs that funannotate clean dropped (an empty file if none).
        awk 'FILENAME == ARGV[1] {{ if (/^>/) kept[substr($1, 2)] = 1; next }}
             /^>/ {{ id = substr($1, 2); if (!(id in kept)) print id }}' \
            {output.contigs_clean} {input.contigs} > {output.contigs_removed}
        """


# ──────────────── repeat_masking — soft-mask the repeats ────────────────
# Genes are predicted on a genome whose repeats are SOFT-masked: written in lowercase, not
# replaced by N, so no sequence is lost, but the gene finders treat the lowercase stretches
# with caution. Without masking, the many copies of transposons would produce false genes.
# Which tool does the masking is a setting (parameters.masking_method; 00_common.smk accepts
# only these two values):
#   tantan                     : fast. Masks tandem repeats and low-complexity sequence.
#   repeatmodeler_repeatmasker : slow (hours). RepeatModeler first builds a library of the
#                                repeat families of THIS genome, then RepeatMasker uses that
#                                library to mask the genome, so interspersed repeats such as
#                                transposons are found too.
# Both branches are in the one shell block below, and the setting picks one when the job runs.
#
# Takes in:    CONTIGS_SORT from funannotate_preprocess.
# Does:        starts from an empty work folder (REPEATMASK_DIR), then masks with the chosen
#              method. The RepeatModeler branch stops with an error if RepeatModeler did not
#              produce a usable repeat library, so an unmasked genome is never passed on as masked.
# Produces:    CONTIGS_MASK, the soft-masked genome. REPEATMASK_DIR keeps the work files of the
#              RepeatModeler branch (the library and the RepeatMasker reports); the tantan branch
#              leaves it empty.
# Consumed by: funannotate_prediction (and genemark_check, when GeneMark is switched on).
#
# v1 -> v2:
#   * NEW: when parameters.repeatmodeler_seed is set, RepeatModeler is called with
#     `-srand <seed>`, and only then. RepeatModeler normally seeds its random choice of genome
#     samples from the clock, so two runs on the same genome give slightly different repeat
#     libraries; a fixed seed fixes that choice, but whether the masking then repeats exactly
#     has not been tested. Empty (the default) keeps the v1 behaviour.
#   * v1 began with `mkdir -p logs`, which makes a folder called logs inside the CURRENT directory
#     (the output folder in v1, but the folder Snakemake was started from in v2). It is dropped:
#     Snakemake creates the log folder itself.
#   * threads: capped_cpus(12) (v1: min(CPUS, 12)); the number is the same.
#   * Note carried over from v1, not changed: RepeatMasker's -pa option counts parallel
#     batches, and each batch runs its search engine on several CPUs of its own (RepeatMasker
#     4.2.3 documents 4 per RMBlast batch), so -pa 12 can keep up to 48 CPUs busy, more than
#     the {threads} Snakemake reserved.
rule repeat_masking:
    input:
        contigs_sort = CONTIGS_SORT,
    output:
        contigs_mask = CONTIGS_MASK,
    params:
        method = MASKING_METHOD,
        work_dir = REPEATMASK_DIR,
        quick_arg = "-quick" if REPEATMODELER_QUICK else "",
        ltrstruct_arg = "-LTRStruct" if REPEATMODELER_LTRSTRUCT else "",
        seed_arg = f"-srand {REPEATMODELER_SEED}" if REPEATMODELER_SEED != "" else "",
    threads: capped_cpus(12)
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/repeatmasking.yaml"
    log:
        LOGS + "/repeat_masking_{sample}.log"
    shell:
        """
        set -euo pipefail

        # Full paths, so that they still point at the right files after the cd below. In v2
        # every path is absolute already; the readlink -f calls are kept from v1, where the
        # paths were relative.
        input_fasta=$(readlink -f {input.contigs_sort})
        output_fasta=$(readlink -f {output.contigs_mask})
        log_file=$(readlink -f {log})
        work_dir="{params.work_dir}"

        # Start from an empty work folder: the files of an earlier attempt are removed.
        rm -rf "$work_dir"
        mkdir -p "$work_dir"

        if [ "{params.method}" = "tantan" ]; then
            # tantan writes the genome with its repeats in lowercase.
            tantan "$input_fasta" > "$output_fasta" 2> "$log_file"
        elif [ "{params.method}" = "repeatmodeler_repeatmasker" ]; then
            # Inside the work folder: BuildDatabase prepares the genome for RepeatModeler, which
            # then builds the library of the repeat families of this genome.
            cd "$work_dir"
            BuildDatabase -name {wildcards.sample}_rm "$input_fasta" > "$log_file" 2>&1
            RepeatModeler -database {wildcards.sample}_rm -threads {threads} {params.quick_arg} {params.ltrstruct_arg} {params.seed_arg} >> "$log_file" 2>&1
            # The library is <sample>_rm-families.fa. If that file is missing or empty, a file
            # called consensi.fa.classified in the work folder is used instead (kept from v1).
            families="{wildcards.sample}_rm-families.fa"
            if [ ! -s "$families" ] && [ -s "consensi.fa.classified" ]; then
                families="consensi.fa.classified"
            fi
            # No usable library: stop, so that an unmasked genome is never passed on as masked.
            if [ ! -s "$families" ]; then
                echo "RepeatModeler did not produce a usable repeat library." >> "$log_file"
                exit 1
            fi
            # Mask the genome with that library. -xsmall writes the repeats in lowercase, and
            # RepeatMasker names its result after the input file, with .masked added.
            RepeatMasker -pa {threads} -xsmall -lib "$families" -dir . "$input_fasta" >> "$log_file" 2>&1
            cp "$(basename "$input_fasta").masked" "$output_fasta"
        fi
        """


# ──────────────── augustus_model_check — which AUGUSTUS parameters? ────────────────
# AUGUSTUS, the gene finder that weighs most in the final gene set, either uses the stored
# parameters of the species named in column 4 of the annotation table, as they are, or is
# trained on this genome's BUSCO genes when column 4 is empty or "-" (00_common.smk, section 9d,
# explains the choice and why FunFlux checks it). This rule makes the choice visible and safe.
#
# Takes in:    FUNANNOTATE_DB_MARKER: the database is unpacked and its paths repaired. That
#              comes early in a run, long before the gene prediction, so a misspelt name stops
#              the run within minutes even on a first run, when the database was not yet on disk
#              while Snakemake read the workflow (when it was, 00_common.smk has already made the
#              same check then). The sample's species (column 2), model (column 4) and whether
#              the long training applies to it come in as params, through species_name,
#              augustus_model and augustus_optimize_for_sample; a change of any of them runs this
#              check, and the prediction after it, again for that sample only.
# Does:        runs workflow/scripts/40_annotation/check_augustus_model.py, which decides the way
#              funannotate does: a named model must be a folder of that exact name (capitals
#              count) in <database>/trained_species/, or the run stops with the closest names
#              the database holds (funannotate itself would silently train under the misspelt
#              name); with no name, a model stored earlier under the sample's own name
#              (<species>_<sample> in lower case) would be used instead of training, and the
#              record says so. The script's output goes to the log; when it stops, the log is
#              also printed to the terminal ("|| { ... }" below runs only when the script fails).
# Produces:    AUGUSTUS_MODEL_RECORD, three tab-separated lines: the sample, the mode
#              (pretrained = a stored model used as it is; trained = trained on this genome's
#              BUSCO genes, optimised or not) and the model name ("-" when trained).
# Consumed by: funannotate_prediction, which waits for it and copies it into its log.
#
# Like normalize_funannotate_db, it runs the script with the pinned Python of envs/python.yaml
# (the script needs only the Python standard library), on the machine that runs Snakemake
# (localrule), because it only reads a folder listing.
#
# CAUTION: Snakemake compares this shell command with the one of the last run. Any change to
# it, a reworded comment included, re-runs this check for every sample, and through the record
# the gene prediction and the whole annotation after it (hours per sample). So does any change
# to its environment, envs/python.yaml (a comment included) or its pin list. Change the script
# (check_augustus_model.py) instead when you can: Snakemake does not compare its content.
localrules: augustus_model_check

rule augustus_model_check:
    input:
        db_marker = FUNANNOTATE_DB_MARKER,
    output:
        record = AUGUSTUS_MODEL_RECORD,
    params:
        funannotate_db = FUNANNOTATEDB,
        species = species_name,
        model = augustus_model,
        optimize = augustus_optimize_for_sample,
    conda:
        "../../envs/python.yaml"
    log:
        LOGS + "/augustus_model_check_{sample}.log"
    shell:
        """
        python {AUGUSTUS_CHECK_SCRIPT:q} \
            --database "{params.funannotate_db}" \
            --sample {wildcards.sample} \
            --species {params.species:q} \
            --model "{params.model}" \
            --optimize {params.optimize} \
            --record {output.record} > {log} 2>&1 || {{ cat {log} >&2; exit 1; }}
        """


# ──────────────── funannotate_prediction — find the genes ────────────────
# funannotate predict finds the protein-coding genes of the genome. It runs several ab initio
# gene finders (programs that recognise genes from the sequence alone: AUGUSTUS, SNAP,
# GlimmerHMM and, if available, GeneMark-ES), aligns protein evidence to the genome (the
# proteins of a related fungus given in the annotation table, and UniProt/SwissProt from the
# funannotate database), and lets EVidence Modeler combine all of these into one consensus gene
# set. When Tiberius is on (the default), the gene models Tiberius predicted before
# (shared/42_tiberius.smk) join that combination as one more source, with weight 1. SNAP and GlimmerHMM are trained on this genome, on a set of conserved genes found with a
# BUSCO dataset of the funannotate database (parameters.prediction.busco_db); GeneMark-ES trains
# itself on the genome. AUGUSTUS is trained on the same BUSCO genes when column 4 of the
# annotation table is empty, and otherwise uses the stored model named there as it is, without
# training (augustus_model_check, above, records which). An exception to "trained on this
# genome": a model someone added to the database with `funannotate species` (named in column 4,
# or stored under this sample's own name) can also hold SNAP, GlimmerHMM and GeneMark-ES
# parameters, which funannotate then uses as they are instead of training. The 166 species of
# the database snapshot hold AUGUSTUS parameters only.
#
# Takes in:    CONTIGS_MASK from repeat_masking, FUNANNOTATE_DB_MARKER (waited for, like in
#              funannotate_preprocess), and AUGUSTUS_MODEL_RECORD from augustus_model_check, which
#              this rule copies into its log. When GeneMark is switched on, also the small file
#              genemark_check writes (through genemark_gate in 00_common.smk): it proves that
#              GeneMark-ES starts with the user's licence key, so an expired key stops the run
#              before the prediction, not partway through it. With GeneMark off nothing is
#              waited for. When Tiberius is on, also TIBERIUS_GFF3 from tiberius_gff3, through
#              tiberius_gate (with Tiberius off, nothing). The species name, the protein evidence
#              file and the AUGUSTUS options of this sample come from the annotation table
#              (files.annotation_params) through the named functions species_name,
#              protein_evidence and augustus_arguments:
#              "--augustus_species <name>" for a stored model; nothing for training on the
#              genome, or "--optimize_augustus" when parameters.prediction.augustus_optimize is
#              true. The Tiberius option comes from tiberius_arguments: "--other_gff
#              <TIBERIUS_GFF3>:1" when Tiberius is on, nothing when it is off.
# Does:        funannotate predict. Before it starts, PATH (the list of folders the shell
#              searches for programs) is cut down to the rule's conda environment plus the
#              system tools, and GENEMARK_PATH is cleared, so that a licensed GeneMark
#              installation that happens to exist on the computer cannot be found by accident;
#              GeneMark is used only when directories.genemark_dir names it and
#              parameters.genemark.mode is not off. The sample name is used as the locus tag
#              (the prefix of every gene ID) and as the strain.
# Produces:    FUNANNOTATE_DIR, declared as a directory(): funannotate chooses the file names
#              inside it, notably predict_results/<Genus_species>_<sample>.{proteins.fa,gbk,gff3,...}.
#              A directory() output means Snakemake treats the whole folder as the product, and
#              deletes the folder before it runs this rule again.
# Consumed by: iprscan_annotation, eggnog_annotation, antismash_annotation and, when
#              switched on, signalp6_annotation (each reads the proteins or the GenBank file),
#              and funannotate_annotation, which adds their results to the same folder.
#
# v1 -> v2:
#   * NEW: --busco_db <parameters.prediction.busco_db>. The default, dikarya, is what funannotate
#     used anyway, so the default result is unchanged.
#   * --cpus takes {threads} (all of CPUS, cut down to --cores if that is smaller); v1 used
#     `resources: cpus`, which Snakemake does not count against --cores.
#   * Species, protein evidence and AUGUSTUS model come from named functions, not lambdas
#     (unnamed one-line functions written inside the rule).
#   * NEW input genemark_check (through genemark_gate); see above.
#   * NEW: Tiberius's gene models, with --other_gff (on by default; parameters.tiberius.mode:
#     off gives the command without it).
#   * NEW: the AUGUSTUS model is optional (v1 required it, so AUGUSTUS was never trained on the
#     genome): column 4 empty or "-" trains AUGUSTUS on the BUSCO genes, optionally with
#     --optimize_augustus. A named model must exist in the database (augustus_model_check), and
#     the log starts with the record of what AUGUSTUS does. For a named model the
#     AUGUSTUS options of funannotate predict are the same as in v1.
#   * The shell adds to the log with >> without emptying it first, as in v1. Attempts still do
#     not mix: Snakemake deletes a job's log file before every attempt, a retry after a failure
#     included (checked in Snakemake 9.14.6), so the log holds only the latest attempt.
rule funannotate_prediction:
    input:
        contigs_mask = CONTIGS_MASK,
        db_marker = FUNANNOTATE_DB_MARKER,
        augustus_record = AUGUSTUS_MODEL_RECORD,
        genemark_check = genemark_gate,
        tiberius_gff3 = tiberius_gate,
    output:
        funannotate_dir = directory(FUNANNOTATE_DIR),
    params:
        genemark_dir = GENEMARKDIR,
        genemark_enabled = str(GENEMARK_ENABLED).lower(),
        funannotate_db = FUNANNOTATEDB,
        busco_db = PREDICTION_BUSCO_DB,
        species = species_name,
        proteins = protein_evidence,
        augustus = augustus_arguments,
        tiberius = tiberius_arguments,
    threads: CPUS
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/funannotate.yaml"
    log:
        LOGS + "/funannotate_prediction_{sample}.log"
    shell:
        """
        export FUNANNOTATE_DB="{params.funannotate_db}"
        # Keep only the conda environment and the system folders on PATH, and forget any
        # GENEMARK_PATH set outside, so that funannotate cannot pick up a licensed GeneMark
        # installation the config did not name. GeneMark is handed over only when it is on.
        export PATH="$CONDA_PREFIX/bin:/usr/bin:/bin"
        unset GENEMARK_PATH
        if [ "{params.genemark_enabled}" = "true" ]; then
            export GENEMARK_PATH="{params.genemark_dir}"
        fi

        # What AUGUSTUS does for this sample (stored model or training), first in the log.
        cat {input.augustus_record} >> {log}

        # Gene prediction. The sample name is the locus-tag prefix (--name) and the strain; the
        # AUGUSTUS options are --augustus_species <model>, nothing, or --optimize_augustus
        # (augustus_arguments); the Tiberius option is --other_gff <file>:1 or nothing
        # (tiberius_arguments); the second file after --protein_evidence is UniProt/SwissProt
        # from the funannotate database.
        funannotate predict \
            -i {input.contigs_mask} \
            -o {output.funannotate_dir} \
            -s {params.species:q} \
            --name {wildcards.sample} \
            --strain {wildcards.sample} \
            {params.augustus} \
            {params.tiberius} \
            --busco_db {params.busco_db} \
            --protein_evidence {params.proteins} {params.funannotate_db}/uniprot_sprot.fasta \
            --cpus {threads} >> {log} 2>&1
        """


# ──────────── iprscan_annotation — protein domains (InterProScan) ────────────
# InterProScan compares every predicted protein with the member databases of InterPro
# (Pfam and others) and reports which protein families and domains it contains. With
# -goterms it also reports the Gene Ontology terms linked to those matches, and with -pa the
# pathway annotation. funannotate annotate later turns this into the domain and GO labels of
# the final annotation.
#
# Takes in:    FUNANNOTATE_DIR from funannotate_prediction; the shell looks inside it for the
#              protein FASTA (predict_results/*.proteins.fa; update_results/ is also checked, for
#              a workflow that includes a funannotate update step). The InterProScan launcher
#              is the program named by files.iprscan.
# Does:        runs InterProScan directly, not through the `funannotate iprscan` wrapper (that
#              avoids problems the wrapper has when it assembles the XML from newer InterProScan
#              releases). It removes an old result file first, so a stale file can never be taken
#              for a fresh one, and afterwards checks that the XML parses and has the expected
#              top-level element (protein-matches), which catches an empty, cut-off or wrong file
#              before funannotate_annotation reads it. The checking Python program sits in the shell
#              as a "here-document": the lines between <<'PY' and PY are handed to `python -`.
# Produces:    IPRSCAN_XML. IPRSCAN_DIR/tmp is InterProScan's scratch folder; it is left in place
#              after the run (InterProScan normally removes its own files inside it when it
#              finishes, so the folder is empty or nearly so).
# Consumed by: funannotate_annotation.
#
# InterProScan runs in its own conda environment (workflow/envs/interproscan.yaml), which
# supplies Java 11 (openjdk 11), Python 3.11 and Perl 5.32: the launcher of InterProScan
# 5.78-109.0 refuses any Java older than 11, InterProScan runs `python3` for seven of its analyses
# (PANTHER, Gene3D and others; the XML check below uses the same Python), and `perl`, with
# thread support, for SUPERFAMILY and the PROSITE patterns. Kept apart, InterProScan does
# not depend on the Java that conda installs in the funannotate environment (v1 noted that
# environment could end up with Java 8; the one pinned for v2 has Java 25). InterProScan
# itself is not in conda: it is a large manual installation, called through files.iprscan.
#
# v1 -> v2:
#   * NEW: -T <IPRSCAN_DIR>/tmp. InterProScan writes its scratch files to a folder called temp
#     inside the CURRENT directory unless told otherwise. In v1 the current directory was the
#     output folder, so this was a top-level temp/ there; v2 does not change the working
#     directory (it is the folder you started Snakemake from), so the scratch is sent to a folder
#     of the sample explicitly. The rule creates it first.
#   * -cpu takes {threads}, with threads capped at 10 (v1: `resources: cpus = min(CPUS, 10)`).
#     InterProScan gains little from more, and very high values put needless load on a shared
#     computer.
#   * The search for the proteins, the removal of the old XML and the XML check are unchanged.
#     Within the check, chr(125) is the closing curly brace; writing it that way keeps Snakemake's
#     own brace-based string formatting from mistaking it for a placeholder.
#   * NEW: retries: 1. A failed InterProScan job is started once more before the run reports
#     it. InterProScan keeps an internal database of its progress in the scratch folder, and on
#     a busy computer it has been seen to lose one of that database's files and stop ("Could not
#     open table", a missing kvstore/*.sst file); the second attempt, with a fresh scratch
#     subfolder, then ran through. A real problem (a wrong launcher path, a broken installation)
#     fails twice and stops the run as before. The log shows the last attempt only.
rule iprscan_annotation:
    input:
        funannotate_dir = FUNANNOTATE_DIR,
    output:
        iprscan_out = IPRSCAN_XML,
    params:
        iprscan = IPRSCAN,
        scratch_dir = IPRSCAN_DIR + "/tmp",
    threads: capped_cpus(10)
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    retries: 1
    conda:
        "../../envs/interproscan.yaml"
    log:
        LOGS + "/iprscan_annotation_{sample}.log"
    shell:
        """
        # Stop at the first command that fails, at an unset variable, and when any command
        # of a pipe fails.
        set -euo pipefail

        # Full path of the prediction folder. In v2 every path is absolute already; the
        # readlink -f calls are kept from v1, where the paths were relative.
        funannotate_dir=$(readlink -f {input.funannotate_dir})

        # Create the InterProScan output folder before turning its path into a full one, so
        # that the rule also works on a brand-new output folder.
        iprscan_dir=$(dirname "{output.iprscan_out}")
        mkdir -p "$iprscan_dir"
        iprscan_dir=$(readlink -f "$iprscan_dir")
        iprscan_xml="$iprscan_dir/$(basename "{output.iprscan_out}")"

        # InterProScan's scratch folder (see -T below); it must exist before the run.
        mkdir -p "{params.scratch_dir}"

        # InterProScan works on the predicted proteins. funannotate predict writes them to
        # predict_results; update_results is also checked, for a workflow that adds a
        # funannotate update step. The first *.proteins.fa found is used.
        proteins=""
        for results_dir in "$funannotate_dir/predict_results" "$funannotate_dir/update_results"; do
            if [ -d "$results_dir" ]; then
                proteins=$(find "$results_dir" -maxdepth 1 -type f -name "*.proteins.fa" -print -quit)
                if [ -n "$proteins" ]; then
                    break
                fi
            fi
        done

        # Stop here if gene prediction did not produce a protein FASTA. Running
        # InterProScan without a real input would only create a misleading error.
        if [ -z "$proteins" ]; then
            echo "Could not find a *.proteins.fa file in $funannotate_dir/predict_results or update_results." > {log}
            exit 1
        fi

        # Remove any previous XML before the run. A stale file must never be
        # mistaken for a successful InterProScan result.
        rm -f "$iprscan_xml"
        echo "Running InterProScan directly on $proteins" > {log}

        # Run InterProScan directly, not through the funannotate iprscan wrapper, which has
        # problems putting the XML together with newer InterProScan releases. This writes the
        # XML funannotate annotate expects. -goterms and -pa add the GO terms and pathways.
        "{params.iprscan}" \
            -i "$proteins" \
            -f XML \
            -goterms \
            -pa \
            -cpu {threads} \
            -T "{params.scratch_dir}" \
            -o "$iprscan_xml" >> {log} 2>&1

        # Check the XML before funannotate_annotation reads it: an empty, cut-off or damaged
        # file, or XML of the wrong kind, stops the rule here. The lines between the two PY
        # markers are the Python program.
        python - "$iprscan_xml" >> {log} 2>&1 <<'PY'
import sys
import xml.etree.ElementTree as ET

# The shell passes the InterProScan XML path as the first Python argument.
path = sys.argv[1]
try:
    # Parse the XML and keep the root element for a type check.
    root = ET.parse(path).getroot()
except Exception as exc:
    sys.stderr.write("InterProScan XML validation failed for %s: %s\\n" % (path, exc))
    sys.exit(1)

# When the XML declares a namespace, ElementTree writes the root tag as the namespace address
# in curly braces followed by protein-matches. Keep only what follows the closing brace.
# chr(125) is that closing curly brace; spelling it this way keeps Snakemake's own
# brace-based string formatting from mistaking it for a placeholder.
tag = root.tag.split(chr(125))[-1]
if tag != "protein-matches":
    sys.stderr.write("InterProScan XML validation failed for %s: unexpected root tag '%s'\\n" % (path, root.tag))
    sys.exit(1)
PY
        """


# ──────────── antiSMASH reference database — one of two rules ────────────
# antiSMASH recognises biosynthetic gene clusters by matching curated cluster and domain
# models, and those models live in a reference database that has to be put in place once.
# Whichever of the two rules below exists in a run makes the same folder, ANTISMASH_DB_DIR
# (00.resources/antismash/databases), and every antiSMASH job waits for it. The choice is made
# while Snakemake reads the workflow, from directories.antismash_db, so a run contains exactly
# one of them and the other is not defined at all:
#   directories.antismash_db empty : antismash_db downloads the database.
#   directories.antismash_db set   : antismash_db_local links to the copy you already hold.
#
# The folder is persistent: it stays in 00.resources/ after the run and is reused by the next
# run on the same output folder. It is not a temp() output (one that Snakemake deletes once no
# later job needs it), and it has no marker file. The workflow NEVER writes into your own
# antiSMASH database: the local rule only fills ANTISMASH_DB_DIR with symbolic links to it (a
# symbolic link is a small file that points at another file or folder), so that Snakemake's
# habit of deleting a directory() output before re-running its rule can only ever delete
# links, never your database.
#
# v1 -> v2: v1 downloaded the database into a temp() folder, so it was deleted after every
# run and downloaded again on every resume. Now it is kept, and an existing copy can be reused.
if ANTISMASHDB:

    # ── antismash_db_local — a view of the antiSMASH database you already hold ──
    # Takes in:    directories.antismash_db (ANTISMASHDB): the folder that holds clusterblast/,
    #              pfam/, ... It is only read.
    # Does:        creates ANTISMASH_DB_DIR and puts one symbolic link in it for every entry of
    #              your folder (ln -sfn: make a link, replace an existing one, and do not
    #              follow a link that already points at a folder). Nothing is copied or rebuilt.
    # Produces:    ANTISMASH_DB_DIR, a folder of links that looks like the database.
    # Consumed by: antismash_annotation, as its --databases folder.
    # Conda:       none; mkdir, ln and basename are enough.
    rule antismash_db_local:
        input:
            source_db = ANTISMASHDB,
        output:
            antismash_db = directory(ANTISMASH_DB_DIR),
        log:
            LOGS + "/antismash_db_local.log"
        shell:
            """
            mkdir -p {output.antismash_db}
            {{
              echo "Building a local antiSMASH database view"
              echo "  source (read-only): {input.source_db}"
              echo "  view:               {output.antismash_db}"
            }} > {log}
            for entry in "{input.source_db}"/*; do
                ln -sfn "$entry" "{output.antismash_db}/$(basename "$entry")"
            done
            """

else:

    # ── antismash_db — download the antiSMASH database ──
    # Takes in:    nothing from the workflow; it downloads from the antiSMASH servers.
    # Does:        download-antismash-databases, a command of the antiSMASH package.
    # Produces:    ANTISMASH_DB_DIR, the database (a large download).
    # Consumed by: antismash_annotation, as its --databases folder.
    # Where it runs: always on the computer that runs Snakemake itself (localrule: True). On
    # a single computer that changes nothing; on a cluster it keeps the download off the
    # compute nodes, which often have no internet.
    rule antismash_db:
        localrule: True
        output:
            antismash_db = directory(ANTISMASH_DB_DIR),
        conda:
            "../../envs/antismash.yaml"
        log:
            LOGS + "/antismash_db.log"
        shell:
            """
            download-antismash-databases \
                --database-dir {output.antismash_db} > {log} 2>&1
            """


# ──────────── antismash_annotation — secondary-metabolite gene clusters ────────────
# Many fungal toxins, pigments, antibiotics and other bioactive compounds (secondary
# metabolites) are made by enzymes whose genes sit next to each other in the genome as a
# biosynthetic gene cluster. antiSMASH finds such clusters and reports each one with its type
# and its closest known relative.
#
# Takes in:    ANTISMASH_DB_DIR from antismash_db or antismash_db_local, and FUNANNOTATE_DIR from
#              funannotate_prediction. The rule lists the whole prediction folder as its input
#              and the shell reaches inside it for the one file it needs, the GenBank file of the
#              predicted genes (predict_results/<Genus_species>_<sample>.gbk, found through the
#              named function funannotate_gbk).
# Does:        runs antiSMASH with the fungal rule set (--taxon fungi) and with
#              --genefinding-tool none, which means "do not predict genes again": the genes
#              funannotate already predicted travel inside the GenBank file. The cluster
#              coordinates therefore sit on the same gene set as the rest of the annotation.
# Produces:    ANTISMASH_DIR, declared as a directory(): antiSMASH writes many files into it (an
#              HTML report, the cluster GenBank files, ...). The file antismash.gbk inside it is
#              the one funannotate_annotation reads, and its name comes from the output basename
#              "antismash" set in params below; change one and the other must follow.
# Consumed by: funannotate_annotation.
#
# v1 -> v2: threads: capped_cpus(12) (v1: `resources: cpus`, capped at 24 in FunFlux and at 12 in
# Funnotator, FunFluxL and FunFluxL+); --cpus takes {threads}. The GenBank file comes from a
# named function, not a lambda.
rule antismash_annotation:
    input:
        antismash_db = ANTISMASH_DB_DIR,
        funannotate_dir = FUNANNOTATE_DIR,
    output:
        antismash_dir = directory(ANTISMASH_DIR),
    params:
        basename = "antismash",
        genefinding_tool = "none",
        gbk_file = funannotate_gbk,
    threads: capped_cpus(12)
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/antismash.yaml"
    log:
        LOGS + "/antismash_annotation_{sample}.log"
    shell:
        """
        antismash \
            --taxon fungi \
            --output-dir {output.antismash_dir} \
            --output-basename {params.basename} \
            --databases {input.antismash_db} \
            --genefinding-tool {params.genefinding_tool} \
            --cpus {threads} \
            {params.gbk_file} > {log} 2>&1
        """


# ──────────── eggnog_annotation — orthologous groups and functions ────────────
# eggNOG-mapper looks for the closest relative of every predicted protein in the eggNOG
# database (a DIAMOND search, -m diamond) and passes on what is known about the orthologous
# group that relative belongs to: COG functional category, Gene Ontology terms, KEGG
# pathways, EC numbers, a preferred gene name and a description. It is the "what do these
# genes do" layer of the annotation.
#
# Takes in:    FUNANNOTATE_DIR from funannotate_prediction as the dependency; the proteins it
#              reads are predict_results/<Genus_species>_<sample>.proteins.fa inside that folder
#              (found through the named function funannotate_proteins). The eggNOG database is
#              the folder given as directories.eggnog_db (EGGNOGDB), prepared by the user.
# Does:        runs emapper.py with DIAMOND (--override replaces the tables of an earlier
#              attempt), then makes sure that the version stamp in the header of the result
#              table is the plain version number ("## emapper-2.1.13"). funannotate annotate
#              reads the version from that line to decide how to read the table's columns.
# Produces:    EGGNOG_DIR with <sample>.emapper.annotations and the other emapper tables, and a
#              scratch folder EGGNOG_DIR/eggnog_tmp inside it. The scratch folder is declared as
#              temp(directory(...)): Snakemake deletes a temp() output after the last rule that
#              needs it has finished; nobody needs this one, so it goes as soon as the job ends.
# Consumed by: funannotate_annotation reads <sample>.emapper.annotations.
#
# The version stamp. eggNOG-mapper 2.1.13 asks git for a version: first in the folder its own
# package is installed in (inside the conda environment) and, if that gives no answer, in the
# folder it was started in. git can answer whenever such a folder lies inside a git working
# copy, and its answer (a tag or a commit id) then ends up in the stamp, for example
# "## emapper-2.1.13-<commit id>". Two lines in the shell deal with this:
#   * the export of GIT_DIR points git at a folder that does not exist, so that git finds no
#     working copy, wherever Snakemake was started and wherever the conda environment sits;
#   * the sed line then rewrites whatever follows "## emapper-" up to the first blank into the
#     plain version number of the installed eggnog-mapper package, so a stamp that did get a
#     suffix (or a "v") is repaired too.
#
# v1 -> v2:
#   * threads: capped_cpus(12) (v1: `resources: cpus = min(CPUS, 24)` in FunFlux, 12 in
#     Funnotator, FunFluxL and FunFluxL+); --cpu takes {threads}. 12 and not 24 because
#     Snakemake reserves the threads a rule declares for as long as the job runs: asking for
#     24, eggNOG could not start until InterProScan (10 threads, for several hours) had
#     finished, which added its half hour to the end of every run; with 12 the two run side by
#     side.
#   * The database folder and the protein file come from EGGNOGDB and a named function. The two
#     outputs are declared exactly as in v1.
#   * The git handling is new. v1 ran in the output folder, which is normally not a git working
#     copy, and exported GIT_CEILING_DIRECTORIES="$PWD", which does not exclude the current folder
#     itself. v2 runs in the folder Snakemake was started from; when that is a git working copy
#     (the FunFlux repository, for example) the stamp of every eggNOG table would carry the
#     commit id of FunFlux instead of the plain version ("## emapper-2.1.13"). GIT_DIR
#     removes the cause, and the sed pattern, now written without the "v" (emapper 2.1.13 writes
#     "emapper-2.1.13", not "emapper-v2.1.13"), repairs a stamp that has a suffix anyway.
rule eggnog_annotation:
    input:
        funannotate_dir = FUNANNOTATE_DIR,
    output:
        temp_dir = temp(directory(EGGNOG_DIR + "/eggnog_tmp")),
        eggnog_dir = directory(EGGNOG_DIR),
    params:
        dmnd_db = EGGNOGDB,
        protein_file = funannotate_proteins,
    threads: capped_cpus(12)
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/eggnog-mapper.yaml"
    log:
        LOGS + "/eggnog_annotation_{sample}.log"
    shell:
        """
        mkdir -p {output.temp_dir} {output.eggnog_dir}

        # eggNOG-mapper can write git's answer (a commit id or a tag) into its version stamp when
        # its package or the current folder lies inside a git working copy. Point git at a
        # folder that does not exist, so that never happens (see the rule's comment).
        export GIT_DIR=/nonexistent

        # Search the predicted proteins against the eggNOG database with DIAMOND. The result
        # tables are named after the sample.
        emapper.py \
            -i {params.protein_file} \
            --output_dir {output.eggnog_dir} \
            --cpu {threads} \
            -m diamond \
            --data_dir {params.dmnd_db} \
            --output {wildcards.sample} \
            --temp_dir {output.temp_dir} \
            --override > {log} 2>&1

        # Rewrite the version stamp of the annotation table to the plain version number of the
        # installed eggnog-mapper package, the form funannotate annotate reads.
        emapper_version=$(python -c 'import importlib.metadata as m; print(m.version("eggnog-mapper"))')
        sed -i "s/^## emapper-[^[:space:]]*/## emapper-${{emapper_version}}/" \
            {output.eggnog_dir}/{wildcards.sample}.emapper.annotations
        """


# signalp6_results_if_enabled — the SignalP 6 input of funannotate_annotation, if there is one.
# Takes in:  the wildcards of one funannotate_annotation job. Snakemake calls the function once
#            per job and hands them over; wildcards.sample is the sample name.
# Returns:   the SIGNALP6_RESULTS file of that sample (the sample name filled into the path)
#            when SignalP 6 is switched on in the config (SIGNALP6_ENABLED, fixed when Snakemake
#            reads the workflow). Otherwise an empty list, which Snakemake reads as "no input
#            here", so nothing waits for a file that no rule will make (signalp6_annotation
#            exists only when SignalP 6 is on).
# Used by:   the input signalp_results of rule funannotate_annotation.
def signalp6_results_if_enabled(wildcards):
    if SIGNALP6_ENABLED:
        return SIGNALP6_RESULTS.format(sample=wildcards.sample)
    return []


# signalp6_option — the SignalP 6 option of funannotate annotate, if there is one.
# Takes in:  the wildcards of one funannotate_annotation job, as above.
# Returns:   ["--signalp", <SIGNALP6_RESULTS of that sample>] when SignalP 6 is switched on,
#            and an empty list otherwise.
# Used by:   the params signalp_option of rule funannotate_annotation. The shell writes it as
#            {params.signalp_option:q}: Snakemake puts each item of the list in quotes for the
#            shell and joins them with a space, and an empty list becomes no text at all, so
#            the command then simply has no --signalp option.
# Why not in the shell: the option used to be a bash list, filled only with SignalP 6 on.
# The shell runs with `set -u` (stop on an unset variable), and bash before 4.4 (CentOS and
# RHEL 7, Ubuntu 16.04) treats an EMPTY list as unset there, so with SignalP 6 off, the
# default, funannotate annotate never started on those systems.
def signalp6_option(wildcards):
    if SIGNALP6_ENABLED:
        return ["--signalp", SIGNALP6_RESULTS.format(sample=wildcards.sample)]
    return []


# ──────────── funannotate_annotation — merge everything into the annotation ────────────
# funannotate annotate takes the gene prediction and adds the functional evidence from the
# four analyses (InterProScan domains, eggNOG functions, antiSMASH clusters and, when
# switched on, SignalP 6 secretion signals) to its own database searches (Pfam, MEROPS
# proteases, dbCAN carbohydrate-active enzymes, UniProt and others). The final annotated
# genome (GenBank, GFF3, tables, protein and transcript FASTA files) is written into the
# annotate_results folder inside FUNANNOTATE_DIR, next to the earlier prediction results.
#
# Takes in:    FUNANNOTATE_DIR (funannotate_prediction), IPRSCAN_XML (iprscan_annotation),
#              EGGNOG_DIR (eggnog_annotation), ANTISMASH_DIR (antismash_annotation) and, only
#              when SignalP 6 is switched on, SIGNALP6_RESULTS (signalp6_annotation), which is
#              listed through the function signalp6_results_if_enabled above.
# Does:        in this order:
#                1. with SignalP 6 on, stops if its result file is missing or empty;
#                2. changes into a work folder, FUNANNOTATE_DIR/annotate_work, because funannotate
#                   annotate drops some of its own files into the folder it is started in: its
#                   log (funannotate-annotate.<random>.log, moved to logfiles/ at the end) and,
#                   at times, genome.fixedproducts;
#                3. cuts PATH down to the rule's conda environment plus the system tools.
#                   funannotate annotate runs SignalP or Phobius by itself when it finds them
#                   and has no results for them (eggNOG-mapper too, but this rule always hands
#                   over the eggNOG results with --eggnog), so an installation found elsewhere
#                   on the computer could make hidden predictions; the results to use are all
#                   handed over explicitly;
#                4. runs funannotate annotate;
#                5. with SignalP 6 on, stops with an error if funannotate did not turn the
#                   SignalP results into the secretome list (annotations.secretome.txt);
#                6. keeps a non-empty genome.fixedproducts in annotate_misc/, writes the marker,
#                   removes the leftover log files, and removes the work folder if it is empty.
# Produces:    the annotation files in FUNANNOTATE_DIR, and ANNOTATION_COMPLETE, an empty marker
#              file (a file whose only job is to show that a step finished). funannotate annotate
#              writes its results (annotate_results/, annotate_misc/) into the folder given with
#              -i (it would ignore -o; the rule passes none). That folder, FUNANNOTATE_DIR, is
#              already the directory() output of funannotate_prediction, and Snakemake refuses
#              an output of one rule inside another rule's directory() output, so the marker is
#              the only output this rule can declare. The marker is permanent.
# Consumed by: cleanup. The annotation files themselves are the end product you collect.
#
# v1 -> v2: threads: capped_cpus(12) (v1: `resources: cpus = min(CPUS, 24)` in FunFlux and
# Funnotator, 12 in FunFluxL and FunFluxL+); --cpus takes {threads}. SignalP 6 enters through
# the named function instead of a lambda. The commands of the shell are unchanged.
#
# The {params.signalp_option:q} in the shell is the SignalP 6 option: --signalp and the results
# path when SignalP 6 is on, nothing when it is off (signalp6_option, above, explains why it is
# built in Snakemake and not in the shell).
rule funannotate_annotation:
    input:
        funannotate_dir = FUNANNOTATE_DIR,
        iprscan_out = IPRSCAN_XML,
        eggnog_dir = EGGNOG_DIR,
        antismash_dir = ANTISMASH_DIR,
        signalp_results = signalp6_results_if_enabled,
    output:
        annotate_complete = ANNOTATION_COMPLETE,
    params:
        funannotate_db = FUNANNOTATEDB,
        signalp_enabled = str(SIGNALP6_ENABLED).lower(),
        signalp_option = signalp6_option,
    threads: capped_cpus(12)
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/funannotate.yaml"
    log:
        LOGS + "/funannotate_annotation_{sample}.log"
    shell:
        """
        set -euo pipefail

        # funannotate finds its database through this variable. Every rule sets it again,
        # because each job runs in a shell of its own.
        export FUNANNOTATE_DB="{params.funannotate_db}"

        # Full paths of the inputs, the log and the marker, so that they still work after
        # the cd into the work folder below. In v2 every path is absolute already; the
        # readlink -f calls are kept from v1, where the paths were relative.
        funannotate_dir=$(readlink -f {input.funannotate_dir})
        antismash_gbk=$(readlink -f {input.antismash_dir}/antismash.gbk)
        iprscan_xml=$(readlink -f {input.iprscan_out})
        eggnog_annotations=$(readlink -f {input.eggnog_dir}/{wildcards.sample}.emapper.annotations)
        log_file=$(readlink -f {log})
        complete_file=$(readlink -f {output.annotate_complete})
        annotate_work="$funannotate_dir/annotate_work"

        # SignalP 6 on: its result must exist and not be empty. It is passed on below as
        # --signalp followed by the file (params.signalp_option).
        if [ "{params.signalp_enabled}" = "true" ]; then
            signalp_results=$(readlink -f {input.signalp_results})
            if [ ! -s "$signalp_results" ]; then
                echo "Expected non-empty SignalP 6 results at $signalp_results" > "$log_file"
                exit 1
            fi
        fi

        # funannotate annotate drops some of its own files into the folder it is started in,
        # so start it in a work folder of its own.
        mkdir -p "$annotate_work"
        cd "$annotate_work"

        # Only the conda environment and the system folders stay on PATH. funannotate annotate
        # runs SignalP or Phobius by itself when it finds them, so a copy installed elsewhere
        # on the computer could add predictions the config did not ask for.
        export PATH="$CONDA_PREFIX/bin:/usr/bin:/bin"

        funannotate annotate \
            -i "$funannotate_dir" \
            --antismash "$antismash_gbk" \
            --iprscan "$iprscan_xml" \
            --eggnog "$eggnog_annotations" \
            {params.signalp_option:q} \
            --cpus {threads} > "$log_file" 2>&1

        # With SignalP 6 on, funannotate must have turned its results into the secretome list.
        if [ "{params.signalp_enabled}" = "true" ] && [ ! -e "$funannotate_dir/annotate_misc/annotations.secretome.txt" ]; then
            echo "Funannotate did not create annotations.secretome.txt from the supplied SignalP 6 results." >> "$log_file"
            exit 1
        fi

        # Keep a non-empty genome.fixedproducts with the other annotation files.
        if [ -s genome.fixedproducts ]; then
            mkdir -p "$funannotate_dir/annotate_misc"
            mv genome.fixedproducts "$funannotate_dir/annotate_misc/genome.fixedproducts"
        fi

        # The marker file: funannotate annotate wrote its results into the prediction folder,
        # which gives Snakemake no new file of its own to watch.
        touch "$complete_file"

        # Remove the *.log files left directly in the prediction folder and funannotate's own
        # log in the work folder, then the work folder itself if it is empty.
        rm -f "$funannotate_dir"/*.log
        rm -f "$annotate_work"/funannotate-annotate.*.log
        rmdir "$annotate_work" 2>/dev/null || true
        """


# ──────────── cleanup — remove the leftover scratch files ────────────
# This last rule carries over v1's clean-up of three kinds of leftovers in the sample's annotation
# folder: clean_* scratch folders (left by a funannotate clean that was stopped half-way, see
# funannotate_preprocess), *.log files and a genome.fixedproducts file. It runs once the
# annotation is complete.
#
# Takes in:    ANNOTATION_COMPLETE from funannotate_annotation, so that it runs only after the
#              annotation has finished.
# Does:        removes clean_* folders, *.log files and genome.fixedproducts from the sample's
#              annotation folder (ANNOTATION_DIR; passed in as params.annotation_dir, with the
#              sample name already filled in), and lists what it removed in the log.
# Produces:    CLEANUP_COMPLETE, an empty marker file. The marker is permanent.
# Consumed by: nothing; it is the last file of the annotation chain, and rule all asks for it.
#
# v1 -> v2:
#   * The paths in the three rm lines were relative to the output folder; they now start from
#     ANNOTATION_DIR, which is an absolute path.
#   * The rule now has a real log: rm -rf became rm -rfv ("verbose": it prints every item it
#     removes) with the output appended to the log. v1 declared a log but wrote nothing to it.
rule cleanup:
    input:
        annotate_complete = ANNOTATION_COMPLETE,
    output:
        cleanup_complete = CLEANUP_COMPLETE,
    params:
        annotation_dir = ANNOTATION_DIR,
    log:
        LOGS + "/cleanup_{sample}.log"
    shell:
        """
        # The annotation marker (.annotation_complete.txt) always exists when this job starts,
        # because Snakemake starts a job only when all its inputs exist; the test is kept from v1.
        if [ -f {input.annotate_complete} ]; then
            echo "Removing scratch files from {params.annotation_dir}" > {log}

            # Remove the scratch files; -v lists every removed item in the log. A pattern that
            # matches nothing removes nothing.
            rm -rfv {params.annotation_dir}/clean_* >> {log} 2>&1
            rm -rfv {params.annotation_dir}/*.log >> {log} 2>&1
            rm -rfv {params.annotation_dir}/genome.fixedproducts >> {log} 2>&1

            # The marker file that tells Snakemake this step finished.
            touch {output.cleanup_complete}
        fi
        """
