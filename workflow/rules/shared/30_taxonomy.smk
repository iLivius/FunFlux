# FunFlux v2.1.0 — Stage 30: the fungal ITS barcode, cut out of the genome and classified
# (rules/shared/30_taxonomy.smk)
#
# The ITS (internal transcribed spacer) is the standard DNA barcode of fungi. It sits in
# the ribosomal DNA (rDNA), which is present as a block of tens to hundreds of identical or
# nearly identical copies in a row (a "tandem array"). Estimated from the sequencing depth
# of 91 fungal genomes, the number ranges from 14 to 1,442, median 82 (Lofgren et al. 2019,
# Molecular Ecology 28: 721-730). One copy reads, in this order:
#
#       SSU  -  ITS1  -  5.8S  -  ITS2  -  LSU
#
# The 5.8S gene sits between the two spacers ITS1 and ITS2. ITSx (a program that finds
# these regions with statistical models called profile HMMs) cuts ITS1 and ITS2 out of
# a sequence, and SINTAX (a classifier that vsearch offers) compares them with UNITE,
# the curated reference database of fungal ITS sequences, and reports a taxonomy with a
# confidence for each rank (phylum, class, ..., species).
#
# Why this module does not simply hand the whole genome to ITSx (see
# docs/methods/its-long-read-contigs.md): on a long-read assembly the whole array
# can sit on one contig of millions of bases, often on the MINUS strand, i.e. reading along
# the reverse complement of the contig as it is stored. ITSx then reports "domains found in
# wrong order, may be chimeric" and writes nothing. The module first finds every 5.8S gene
# and cuts a short window around it, turned so that it reads SSU - ITS1 - 5.8S - ITS2 - LSU,
# and gives ITSx only those windows.
#
# The chain, in order:
#
#   unite_db_local OR unite_db : the UNITE reference, once per run (not per sample): a local
#                                copy that the user already holds is linked, otherwise the
#                                file is downloaded. Only one of the two rules exists in a
#                                given run, decided when the workflow is read.
#   unite_identified           : unite_identified.py keeps the UNITE records whose taxonomy
#                                names a genus, once per run; a record without one can stop the
#                                SINTAX result above the genus (why: the rule's comment). Not
#                                planned when parameters.its.unite_identified_only is false.
#   its_windows                : nhmmer finds the 5.8S genes with ITSx's own profiles, and
#                                its_windows.py cuts one oriented window per gene copy.
#   its_extraction             : ITSx on the windows -> ITS1, ITS2 and the whole region, one
#                                record per rDNA copy.
#   its_classification         : its_distinct.py collapses the identical copies, counts them
#                                and drops records of implausible length; vsearch --sintax
#                                classifies the distinct ITS1 and ITS2 sequences against UNITE.
#
# The same chain as files. "{sample}" in a file name is a Snakemake wildcard: a per-sample rule
# runs once for every sample, with the sample name filled in.
#
#   contigs_final.fasta ──its_windows──> its_windows.fasta, its_windows.tsv
#        ──its_extraction──> {sample}.ITS1.fasta, .ITS2.fasta, .full.fasta, .summary.txt
#        ──its_classification──> {sample}_ITS_copies.txt, {sample}_ITS1_taxonomy.txt,
#                                {sample}_ITS2_taxonomy.txt           (UNITE comes in at the last step)
#
# and, once per run, the UNITE reference that the last step reads:
#
#   unite_its_sintax.fasta ──unite_identified──> unite_its_sintax.identified.fasta, identified_filter.txt
#
# Where this file sits in the workflow. The chain starts from FINAL_CONTIGS
# (02.assembly/{sample}/contigs_final.fasta), the nuclear genome that the assembly stages end
# with, in every mode. (The contigs-mode predecessor of this workflow ran ITSx on the file
# exactly as supplied.) The results go to 03.taxonomy/{sample}/; only the UNITE files, shared
# by all samples, go to 00.resources/unite/. No rule outside this file reads the results: they
# are for people. rule all (the rule at the end of the Snakefile that lists the files a finished
# run must have) asks for the four ITSx files, the copies table and the two taxonomy files. The
# other files are made on the way: rule its_windows writes the two windows files because ITSx
# needs the windows, rule unite_db (or unite_db_local) writes the UNITE files because SINTAX
# needs UNITE, and rule unite_identified writes the filtered reference because SINTAX reads
# that one. The annotation stage (40_annotation.smk) starts from the same genome and runs
# alongside.
#
# Everything referenced here is defined once in 00_common.smk and not worked out again here:
# DIR_RESOURCES, UNITE_FASTA, UNITE_FASTA_IDENTIFIED, UNITE_IDENTIFIED_COUNTS, UNITE_REFERENCE,
# UNITEDB, UNITEURL, FINAL_CONTIGS, ITS_WINDOWS, ITS_WINDOW_TABLE, ITS_FULL, ITS_ITS1, ITS_ITS2,
# ITS_SUMMARY, ITS_COPIES, ITS1_TAXONOMY, ITS2_TAXONOMY, TAXONOMY_DIR, UNITE_IDENTIFIED_SCRIPT,
# ITS_WINDOWS_SCRIPT, ITS_DISTINCT_SCRIPT, UNITE_PROVENANCE, ITS_WINDOW_FLANK,
# ITS_MIN_PART_LENGTH, ITS_MAX_PART_LENGTH, ITS_MAX_FULL_LENGTH, ITS_TAXONOMY_CUTOFF, LOGS and
# capped_cpus() (capped_cpus(N) asks for N threads, or for resources.threads when that is
# lower). This file defines no constant of its own.
#
# conda: environment paths are written relative to THIS file (workflow/rules/shared/), so
# "../../envs/itsx.yaml" climbs shared/ -> rules/ -> workflow/ and then goes into envs/.
# The itsx and vsearch environments each hold their own Python 3.12, pinned like every other
# package. The three small Python scripts called below, its_windows.py (rule its_windows, itsx
# environment), unite_identified.py (rule unite_identified, vsearch environment) and
# its_distinct.py (rule its_classification, vsearch environment), use only the standard
# library and are started as `python`. Activating an environment puts its own programs
# first on the PATH (the list of folders the shell searches for a program), so `python` is the
# environment's Python. Before Python was added to these two environments (on 2026-10-08), the
# scripts ran with whatever `python` the computer has (for example /usr/bin/python), and a
# computer with only `python3` stopped the job.
#
# Two shell habits appear below. Snakemake runs every shell block in bash "strict mode"
# (set -euo pipefail; some rules also write it out): a command that fails stops the job at
# once, so a command that may fail harmlessly is followed by "|| true". And a line
# trap '<command>' EXIT runs <command> when the job's shell ends, after success or failure;
# here it deletes temporary folders.


# ═════════════ The UNITE reference database (once per run) ═════════════
# Both rules below make the same two files, UNITE_FASTA (00.resources/unite/
# unite_its_sintax.fasta, the UNITE sequences in the format SINTAX reads) and a short text
# file PROVENANCE.txt beside it, so that anybody can later see where the database came from
# and that it is the file that was used. Exactly one of the two rules is defined, chosen
# while Snakemake reads the workflow: it depends only on the configuration, never on the data.
# UNITEDB is the folder named by directories.unite_db, already checked in 00_common.smk, or
# empty when that key is not set; empty means "download it", so rule unite_db is defined.
if UNITEDB:

    # ───────────── unite_db_local — use the UNITE copy the user already holds ─────────────
    # The user pointed directories.unite_db at a folder that contains
    # unite_its_sintax.fasta (checked in 00_common.smk). UNITE is large (about 217 MB for the
    # release the config template downloads), and many projects share one copy, so nothing is
    # downloaded or copied.
    #
    # Takes in:    nothing from other rules; the folder named by directories.unite_db
    #              (UNITEDB). The workflow only ever READS it.
    # Does:        makes a symbolic link (a pointer, not a copy) at UNITE_FASTA to the file in
    #              that folder: the "local view". -s makes a symbolic link, -f replaces a link
    #              that is already there, -n treats an existing link to a folder as a plain
    #              name. The log gets ln's messages and one line saying what was linked. Then
    #              it writes PROVENANCE.txt, one item per line: the source folder, the date in
    #              UTC, and the md5 checksum and size in bytes of the FASTA (md5sum and wc read
    #              through the link, so the checksum identifies the exact file that was used).
    # Produces:    UNITE_FASTA (the link) and PROVENANCE.txt, both in 00.resources/unite/.
    # Consumed by: rule unite_identified (the FASTA); rule its_classification reads the FASTA
    #              directly only when parameters.its.unite_identified_only is false.
    #
    # New in v2: v1 always downloaded UNITE. No conda environment: ln, md5sum and wc come
    # from the computer.
    rule unite_db_local:
        output:
            unite_fasta = UNITE_FASTA,
            provenance = UNITE_PROVENANCE,
        params:
            unite_dir = UNITEDB,
        log:
            LOGS + "/unite_db.log"
        shell:
            """
            mkdir -p $(dirname {output.unite_fasta})
            ln -sfn {params.unite_dir:q}/unite_its_sintax.fasta {output.unite_fasta} > {log} 2>&1
            echo "Linked the local UNITE database {params.unite_dir}/unite_its_sintax.fasta as {output.unite_fasta}" >> {log}

            echo "source_directory: {params.unite_dir}" > {output.provenance}
            echo "linked_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> {output.provenance}
            echo "md5: $(md5sum {output.unite_fasta} | cut -d' ' -f1)" >> {output.provenance}
            echo "size_bytes: $(wc -c < {output.unite_fasta})" >> {output.provenance}
            """

else:

    # ───────────── unite_db — download the UNITE reference ─────────────
    # Takes in:    nothing from other rules; the download link links.unite_its_link from the
    #              configuration (UNITEURL). The link points to a gzip-compressed FASTA.
    # Does:        downloads it with curl (-L follows redirects, --fail makes an HTTP error
    #              stop the rule instead of saving the error page) and unpacks it on the fly
    #              with gzip into UNITE_FASTA, so the compressed file is never stored. The
    #              messages of curl and gzip go to the log. Then it writes PROVENANCE.txt, one
    #              item per line: the link, the date in UTC, and the md5 checksum and size in
    #              bytes of the unpacked FASTA.
    # Produces:    UNITE_FASTA and PROVENANCE.txt in 00.resources/unite/. They are kept, so
    #              a later run (or a second project using the same output folder) does not
    #              download the database again.
    # Consumed by: rule unite_identified (the FASTA); rule its_classification reads the FASTA
    #              directly only when parameters.its.unite_identified_only is false.
    #
    # Environment: envs/download.yaml, which supplies curl and gzip (a minimal Linux
    # installation may have no curl); md5sum, date and wc come from the computer (GNU
    # coreutils, part of every Linux base system).
    # Where it runs: always on the computer that runs Snakemake itself (localrule: True). On
    # a single computer that changes nothing; on a cluster it keeps the download off the
    # compute nodes, which often have no internet.
    #
    # v1 -> v2: the database now lives in 00.resources/unite/ (run-level material) instead
    # of inside the ITS folder of the first stage that needed it, and PROVENANCE.txt is new.
    # The download command itself is unchanged. v1 ran it in the vsearch environment, which
    # holds neither curl nor gzip, so both came from the computer.
    rule unite_db:
        localrule: True
        output:
            unite_fasta = UNITE_FASTA,
            provenance = UNITE_PROVENANCE,
        params:
            link = UNITEURL,
        conda:
            "../../envs/download.yaml"
        log:
            LOGS + "/unite_db.log"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.unite_fasta})
            curl -L --fail "{params.link}" 2> {log} | gzip -dc > {output.unite_fasta} 2>> {log}

            echo "source_url: {params.link}" > {output.provenance}
            echo "downloaded_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> {output.provenance}
            echo "md5: $(md5sum {output.unite_fasta} | cut -d' ' -f1)" >> {output.provenance}
            echo "size_bytes: $(wc -c < {output.unite_fasta})" >> {output.provenance}
            """


# ───────────── unite_identified — keep the UNITE records whose taxonomy names a genus ─────────────
# Why: every UNITE record carries its taxonomy in the header, rank by rank
# (">ACC|SH...;tax=d:Fungi,p:Ascomycota,...,g:Botrytis,s:Botrytis_cinerea;"). Many records are
# identified only down to some rank above the genus, and their taxonomy then simply stops
# there: "tax=d:Fungi;" alone, for example. In the release the config template downloads,
# about a third of the records name no genus. SINTAX works out its confidence values by
# classifying the query many times from random samples of its k-mers (short words of the
# sequence) and counting how often each name comes out; when a kingdom-only record is
# identical to the query, SINTAX follows it as often as the identified records, and the
# result stops at the kingdom even for a genome whose ITS matches well-identified species
# perfectly. Classified against the records that carry a genus, the same ITS gets its genus
# and species with the usual confidence values.
#
# Takes in:    UNITE_FASTA, the whole release, from rule unite_db or unite_db_local.
# Does:        unite_identified.py (scripts/30_taxonomy/unite_identified.py) reads the file
#              record by record and copies to the output every record whose tax= field holds
#              a rank "g:" followed by a name. A record whose taxonomy stops above the genus,
#              or whose header has no tax= field, is dropped. The script prints the two
#              counts on standard output, which the shell saves as identified_filter.txt, and
#              a one-line summary on standard error, which goes to the log; the counts are
#              then appended to the log too, so it tells the whole story on its own. The
#              script stops the job if no record at all carries a genus (then the file is
#              not a UNITE SINTAX release, or the download was empty). It takes a second or
#              two for the default release (about 217 MB).
# Produces:    UNITE_FASTA_IDENTIFIED (00.resources/unite/unite_its_sintax.identified.fasta)
#              and UNITE_IDENTIFIED_COUNTS (00.resources/unite/identified_filter.txt, two
#              lines: kept_with_genus: <n> and dropped_without_genus: <n>). Both are kept
#              with the UNITE files, so a later run on the same output folder reuses them.
# Consumed by: rule its_classification (the FASTA, as UNITE_REFERENCE) while
#              parameters.its.unite_identified_only is true, the default. With false,
#              UNITE_REFERENCE is UNITE_FASTA, nothing asks for this rule's files, and the
#              rule is defined but never planned. The counts file is for people.
#
# New in v2. The conda environment is vsearch's, for its Python (see the note at the top of
# this file).
rule unite_identified:
    input:
        unite_fasta = UNITE_FASTA,
    output:
        identified = UNITE_FASTA_IDENTIFIED,
        counts = UNITE_IDENTIFIED_COUNTS,
    params:
        script = UNITE_IDENTIFIED_SCRIPT,
    conda:
        "../../envs/vsearch.yaml"
    log:
        LOGS + "/unite_identified.log"
    shell:
        """
        set -euo pipefail

        # The counts (standard output) become identified_filter.txt; the summary line
        # (standard error) starts the log, and the counts are appended to it afterwards.
        python {params.script:q} \
            --fasta {input.unite_fasta} \
            --out-fasta {output.identified} > {output.counts} 2> {log}
        cat {output.counts} >> {log}
        """


# ═════════════ The ITS barcode of each sample (once per sample) ═════════════
# The three rules below run one after the other for every sample: its_windows finds the rDNA
# copies, its_extraction cuts ITS1 and ITS2 out of them, and its_classification classifies the
# distinct sequences against UNITE.

# ───────────── its_windows — find every rDNA copy and cut a window around it ─────────────
# Takes in:    FINAL_CONTIGS (contigs_final.fasta), the nuclear genome, from rule
#              separate_mitochondria when the contamination screen runs, otherwise (contigs
#              mode without the screen) from rule filter_contigs.
# Does:        four steps in one shell block, numbered (1) to (4) there too:
#              (1) limits OpenBLAS to one thread. nhmmer, the HMMER search program that this
#                  rule and ITSx both run, is built with OpenBLAS (bioconda's hmmer 3.4), a
#                  maths library that starts one worker thread per CPU core of the computer,
#                  whatever --cpu says: on a large server, hundreds of threads for one ITSx
#                  call, far more than the threads Snakemake booked for the job.
#                  OMP_NUM_THREADS=1 and OPENBLAS_NUM_THREADS=1 stop that, so the job uses only
#                  the threads it asked for.
#              (2) cuts ITSx's own 5.8S profiles out of the profile file that ITSx ships
#                  (ITSx_db/HMMs/F.hmm, the fungal profiles, inside the active conda
#                  environment). Of its 538 profiles, those with "58S" in their NAME line
#                  (358 with ITSx 1.1.3) are the anchors at the start and the end of the
#                  5.8S gene. One profile name occurs twice in the file, so hmmfetch (which
#                  fetches by name) cannot be used; instead a small awk program splits the
#                  file at the "//" line that ends each profile and keeps every profile whose
#                  NAME line has 58S. Two checks stop the job with a [FunFlux] message, in the
#                  log and on the screen: before awk, when F.hmm is missing or empty (the job
#                  then does not run in the itsx environment, or that environment does not
#                  hold an intact ITSx 1.1.3); after awk, when no NAME line has 58S (this ITSx
#                  version would then name its 5.8S profiles differently). The number of
#                  profiles kept goes to the log. grep -c exits with status 1 when it counts
#                  nothing, hence its "|| true".
#              (3) nhmmer searches the genome with those profiles, on the job's threads.
#                  -E 0.001 reports only hits with an E-value of at most 0.001 (the E-value is
#                  the number of hits this good that would be expected by chance in a search
#                  of this size; smaller is better), and --dna declares the alphabet as DNA.
#                  The hit table (--tblout) goes to a temporary folder that is deleted when
#                  the job ends, and nhmmer's long printed alignments are discarded.
#              (4) its_windows.py keeps the hits with an E-value of at most 1e-5, merges the
#                  hits on the same contig and strand into one site (one 5.8S gene) as long as
#                  each starts at most 1000 bases after the end of the site built so far, cuts
#                  parameters.its.window_flank bases on each side of the middle of the site
#                  (2000 by default, so a full window is 4001 bases: room for the end of the
#                  SSU, both spacers and the start of the LSU), turns the window around
#                  (reverse complement) if the gene is on the minus strand, and writes the
#                  windows and a table. The 1e-5 and the 1000 are the script's own defaults,
#                  which the rule keeps. Its one-line summary goes to the log.
#              A genome in which no 5.8S is recognised is NOT an error: the script writes an
#              empty FASTA and a table with only a header. The two later rules still finish:
#              the four ITSx files and the two taxonomy files are empty, and
#              {sample}_ITS_copies.txt holds its header and one line saying that no ITS region
#              was extracted.
# Produces:    ITS_WINDOWS (03.taxonomy/{sample}/its_windows.fasta; headers
#              contig:start-end_strand, one window per rDNA copy) and ITS_WINDOW_TABLE
#              (its_windows.tsv: the same windows with their evidence, one row each).
# Consumed by: rule its_extraction (the FASTA). The table is for people who want to see where
#              the rDNA copies lie in the genome.
#
# Threads: 8. nhmmer searches a 40 Mb genome in about a minute with 4 threads, and more gains
# little. The number matters beyond speed: Snakemake reserves the threads a rule declares for
# as long as the job runs, so a rule that asks for more than it can use keeps other jobs
# waiting. (With 24, this rule and the two below could not start while InterProScan held 10 of
# a run's 25 threads, which it does for hours.)
#
# New in v2 (docs/methods/its-long-read-contigs.md); v1 gave the whole genome to ITSx.
rule its_windows:
    input:
        contigs = FINAL_CONTIGS,
    output:
        windows = ITS_WINDOWS,
        table = ITS_WINDOW_TABLE,
    params:
        script = ITS_WINDOWS_SCRIPT,
        flank = ITS_WINDOW_FLANK,
    threads: capped_cpus(8)
    conda:
        "../../envs/itsx.yaml"
    log:
        LOGS + "/its_windows_{sample}.log"
    shell:
        """
        # (1) One OpenBLAS thread, so that nhmmer uses only the threads this job asked for.
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

        # A temporary folder for the anchor profiles and the hit table, deleted when the job ends.
        work_dir=$(mktemp -d)
        trap 'rm -rf "$work_dir"' EXIT

        # (2) F.hmm, ITSx's file of fungal profiles, inside the active itsx environment. Stop
        # with a clear message if it is missing or empty (-s: the file exists and is not
        # empty). Without this test, a missing file would make awk fail with an error that
        # reaches only Snakemake's screen, and an empty one would reach the anchor count
        # below, whose message blames the profile names.
        fungal_profiles="$CONDA_PREFIX/bin/ITSx_db/HMMs/F.hmm"
        if [ ! -s "$fungal_profiles" ]; then
            echo "[FunFlux] its_windows: ITSx's fungal profile file $fungal_profiles is missing or empty. Is the job running in its itsx conda environment (--sdm conda), with the pinned ITSx 1.1.3?" | tee -a {log} >&2
            exit 1
        fi

        # Copy every profile whose NAME line contains 58S into anchors_58S.hmm. awk adds
        # each line of F.hmm to record (ORS is awk's output line break); a NAME line with 58S
        # marks the profile as wanted; the "//" line that ends every profile prints the record
        # if it was wanted and starts the next one empty.
        awk '
            {{ record = record $0 ORS }}
            $1 == "NAME" && index($0, "58S") > 0 {{ keep = 1 }}
            $1 == "//" {{ if (keep) printf "%s", record; record = ""; keep = 0 }}
        ' "$fungal_profiles" > "$work_dir/anchors_58S.hmm"

        # Count the anchors for the log; none at all stops the job with a clear message.
        n_anchors=$(grep -c '^NAME' "$work_dir/anchors_58S.hmm" || true)
        echo "its_windows: $n_anchors 5.8S anchor profiles taken from ITSx's F.hmm" > {log}
        if [ "$n_anchors" -eq 0 ]; then
            echo "[FunFlux] its_windows: no 5.8S profile found in $fungal_profiles (no NAME line contains 58S). Does this ITSx version name its 5.8S profiles differently?" | tee -a {log} >&2
            exit 1
        fi

        # (3) Search the genome for the 5.8S anchors. Only the hit table is kept: the printed
        # alignments are thrown away, nhmmer's warnings go to the log.
        nhmmer --cpu {threads} --dna -E 0.001 --tblout "$work_dir/anchors_vs_genome.tbl" "$work_dir/anchors_58S.hmm" {input.contigs} > /dev/null 2>> {log}

        # (4) Merge the hits into one site per 5.8S gene and cut an oriented window around
        # each site: the windows FASTA (read by ITSx next) and the table.
        python {params.script:q} \
            --tblout "$work_dir/anchors_vs_genome.tbl" \
            --fasta {input.contigs} \
            --flank {params.flank} \
            --out-fasta {output.windows} \
            --out-table {output.table} >> {log} 2>&1
        """


# ───────────── its_extraction — ITSx: cut ITS1, ITS2 and the whole region out of the windows ─────────────
# Takes in:    ITS_WINDOWS from rule its_windows, one oriented window per rDNA copy.
# Does:        runs ITSx on the windows with the options of the v1 workflow (all but one, see
#              v1 -> v2 below):
#                -E 0.001      the E-value limit for a hit to count (ITSx's own default is 1e-5)
#                --complement T  also check the opposite strand (the windows are already turned
#                              to read SSU - ITS1 - 5.8S - ITS2 - LSU, so this rarely matters)
#                --nhmmer T    search with nhmmer instead of hmmsearch
#                --preserve T  keep the window names as the record names, instead of ITSx's own
#                --only_full F also report partial regions, not only complete ones
#                --graphical F --positions F --not_found F   do not write the picture, the table of
#                              positions and the list of sequences without ITS
#                --cpu         the job's threads
#              ITSx works in a temporary folder next to its output ({sample}.ITSx_tmp, given
#              with --temp), which is deleted when the job ends, even on failure (the trap line).
#              The -o prefix is TAXONOMY_DIR/{sample}, which makes ITSx write {sample}.ITS1.fasta,
#              {sample}.ITS2.fasta, {sample}.full.fasta and {sample}.summary.txt; its messages go
#              to the log. The OpenBLAS line is explained at rule its_windows.
#              If the windows file is EMPTY (the genome has no recognisable 5.8S), ITSx is not
#              run at all: the four output files are created empty and a [FunFlux] line in the
#              log says why. This branch is necessary. Given an empty file, ITSx 1.1.3 writes
#              only its summary file; with more than one thread it prints "FATAL ERROR" and still
#              exits with status 0, with one thread it exits with status 1. Either way the job
#              would fail: at once on status 1, or on status 0 because the other three outputs
#              are missing.
# Produces:    ITS_ITS1, ITS_ITS2 and ITS_FULL: one record per rDNA copy, named by its
#              window, plus ITS_SUMMARY (ITSx's own report). All four are empty for a genome
#              without a recognisable 5.8S. The "full" file is the whole ITS1 + 5.8S + ITS2
#              region.
# Consumed by: rule its_classification (ITS1, ITS2 and full). All four are also final results
#              that rule all asks for.
#
# v1 -> v2: (1) the input is the windows from rule its_windows and not the genome (see the
# top of this file); (2) the CPU count is Snakemake's threads:, 8 (the windows are a few
# kilobases each, so ITSx finishes in seconds; see the note on threads at rule its_windows),
# instead of a private "resources: cpus" of 24; (3) the OpenBLAS line and the empty-windows
# branch are new; (4) v1's "--plus T" is gone: ITSx 1.1.3 has no option of that name (it is not
# in its help, and its script compares each argument with a fixed list of option names that
# does not contain it), so ITSx skipped it and its value, and leaving it out changes nothing;
# (5) every other flag, and their order, are v1's.
rule its_extraction:
    input:
        windows = ITS_WINDOWS,
    output:
        full = ITS_FULL,
        its1 = ITS_ITS1,
        its2 = ITS_ITS2,
        summary = ITS_SUMMARY,
    params:
        basename = TAXONOMY_DIR + "/{sample}",
    threads: capped_cpus(8)
    conda:
        "../../envs/itsx.yaml"
    log:
        LOGS + "/its_extraction_{sample}.log"
    shell:
        """
        set -euo pipefail
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

        mkdir -p $(dirname "{params.basename}")

        # Windows found: run ITSx on them, in a fresh temporary folder that is removed when
        # the job ends.
        if [ -s {input.windows} ]; then
            temp_dir="{params.basename}.ITSx_tmp"
            rm -rf "$temp_dir"
            mkdir -p "$temp_dir"
            trap 'rm -rf "$temp_dir"' EXIT

            ITSx \
                -i {input.windows} \
                -o {params.basename} \
                --temp "$temp_dir" \
                --preserve T \
                --only_full F \
                --graphical F \
                --complement T \
                --positions F \
                --not_found F \
                -E 0.001 \
                --nhmmer T \
                --cpu {threads} > {log} 2>&1
        else
            # No window (no 5.8S in this genome): ITSx would fail on the empty file, so write
            # the four outputs empty and say why in the log.
            echo "[FunFlux] its_extraction: no 5.8S window was found in the genome of {wildcards.sample} (the windows file is empty), so ITSx was not run. The ITS files are written empty." > {log}
            touch {output.full} {output.its1} {output.its2} {output.summary}
        fi
        """


# ───────────── its_classification — collapse the copies, then classify with SINTAX ─────────────
# Takes in:    ITS_ITS1, ITS_ITS2 and ITS_FULL from rule its_extraction (one record per rDNA
#              copy), and UNITE_REFERENCE: the UNITE records whose taxonomy names a genus
#              (UNITE_FASTA_IDENTIFIED, from rule unite_identified) by default, or the whole
#              release (UNITE_FASTA, from rule unite_db or unite_db_local) when
#              parameters.its.unite_identified_only is false. 00_common.smk picks the one
#              when the workflow is read.
# Does:        (1) its_distinct.py (scripts/30_taxonomy/its_distinct.py) collapses the copies
#                  and checks their length. A long-read assembly can hold several copies, and
#                  ITSx returns one record per copy: five identical ITS2 records for an assembly
#                  with five copies. The script keeps one record per distinct sequence, counts
#                  its copies, and applies a LENGTH GUARD: ITS1 and ITS2 of fungi are 50 to 600
#                  bases long, the whole region at most parameters.its.max_full_length (1500
#                  bases). A record outside these limits is a mis-extraction (v1 could report
#                  a huge "ITS2" that SINTAX then classified as a non-fungus). Such records are
#                  not classified; they are listed with the reason. It writes the table
#                  {sample}_ITS_copies.txt (one row per distinct sequence: region, sequence
#                  name, length, number of copies, status ok or rejected_..., and the windows it
#                  came from) and one FASTA of the distinct, accepted sequences for ITS1 and one
#                  for ITS2, named <sample>_ITS1_<n> and <sample>_ITS2_<n>, so every line of the
#                  SINTAX output can be found in the table. The two FASTA files go to a
#                  temporary folder that is deleted when the job ends; the script's one-line
#                  summary goes to the log.
#              (2) vsearch --sintax classifies each of those FASTA files against the UNITE
#                  reference (--db, the file named by UNITE_REFERENCE, see Takes in):
#                  --tabbedout writes the result, --sintax_cutoff (parameters.its_taxonomy_cutoff,
#                  0.8 by default) is the confidence a rank needs to appear in the last column of
#                  the output, and --randseed 1 with --threads 1 makes the confidences repeatable
#                  (see v1 -> v2 below). If a FASTA is empty (ITSx found nothing, or every record
#                  was rejected), vsearch is not run: an empty taxonomy file is written and a line
#                  in the log says why, so the files always exist.
# Produces:    ITS_COPIES ({sample}_ITS_copies.txt), ITS1_TAXONOMY and ITS2_TAXONOMY
#              ({sample}_ITS1_taxonomy.txt, {sample}_ITS2_taxonomy.txt), in 03.taxonomy/{sample}/.
#              A taxonomy file has no header line. Each line has four tab-separated columns: the
#              record name, the taxonomy with the confidence of each rank, the strand, and the
#              taxonomy kept after the cutoff. An example line of SAMPLE1_ITS2_taxonomy.txt,
#              column by column:
#                1  SAMPLE1_ITS2_1
#                2  d:Fungi(1.00),p:Ascomycota(1.00),...,s:Fusarium_graminearum(1.00)
#                3  +
#                4  d:Fungi,p:Ascomycota,...,s:Fusarium_graminearum
# Consumed by: nothing in the workflow; these are final results, read by people (rule all asks
#              for them).
#
# v1 -> v2: v1 classified every record ITSx returned. v2 first collapses identical copies and
# applies the length guard, so the SINTAX input is one record per
# distinct sequence, with a new name. The vsearch options, their order and the branch for an
# empty FASTA are v1's, with one option added and one thread instead of 24:
#   * --randseed 1 is new. SINTAX gets its confidence values by classifying a query many times,
#     each time from a random sample of its k-mers (short words of the sequence). Left alone,
#     vsearch seeds that random choice from the clock, and the same sequence gets slightly
#     different confidences on every run (the genus confidence, for example, can move by a few
#     hundredths). With a fixed seed, repeated runs give identical results. The taxonomy itself
#     usually stays the same and only the confidence numbers wobble, which
#     would make every later comparison of two runs fail for no reason.
#   * threads: 1 (v1 gave vsearch 24). The seed alone is not enough with several threads: the
#     vsearch manual warns that "when using multiple threads, the --randseed option may not work
#     as intended, because sequences may be processed in a random order by different threads",
#     and asks for --threads 1 together with a fixed seed. One thread costs nothing here: each
#     vsearch call takes a few seconds, almost all of it spent reading and indexing UNITE, and
#     eight threads are no faster (vsearch 2.30.0).
rule its_classification:
    input:
        its1 = ITS_ITS1,
        its2 = ITS_ITS2,
        full = ITS_FULL,
        its_db = UNITE_REFERENCE,
    output:
        its1_tax = ITS1_TAXONOMY,
        its2_tax = ITS2_TAXONOMY,
        copies = ITS_COPIES,
    params:
        script = ITS_DISTINCT_SCRIPT,
        min_part_length = ITS_MIN_PART_LENGTH,
        max_part_length = ITS_MAX_PART_LENGTH,
        max_full_length = ITS_MAX_FULL_LENGTH,
        cutoff = ITS_TAXONOMY_CUTOFF,
    threads: 1
    conda:
        "../../envs/vsearch.yaml"
    log:
        LOGS + "/its_classification_{sample}.log"
    shell:
        """
        work_dir=$(mktemp -d)
        trap 'rm -rf "$work_dir"' EXIT

        # Collapse identical copies and apply the length guard. The distinct ITS1 and ITS2
        # sequences go to the temporary folder, the copies table to its final place.
        python {params.script:q} \
            --sample {wildcards.sample} \
            --its1 {input.its1} \
            --its2 {input.its2} \
            --full {input.full} \
            --min-part-length {params.min_part_length} \
            --max-part-length {params.max_part_length} \
            --max-full-length {params.max_full_length} \
            --out-its1 "$work_dir/ITS1.distinct.fasta" \
            --out-its2 "$work_dir/ITS2.distinct.fasta" \
            --out-table {output.copies} > {log} 2>&1

        # Classify the distinct ITS1 sequences if there is at least one; otherwise write an
        # empty taxonomy file and say why in the log.
        if [ -s "$work_dir/ITS1.distinct.fasta" ]; then
            vsearch \
                --sintax "$work_dir/ITS1.distinct.fasta" \
                --db {input.its_db} \
                --tabbedout {output.its1_tax} \
                --sintax_cutoff {params.cutoff} \
                --randseed 1 \
                --threads {threads} >> {log} 2>&1
        else
            echo "No ITS1 sequence to classify for sample {wildcards.sample} (ITSx found none, or none passed the length guard). Skipping classification." >> {log} 2>&1
            touch {output.its1_tax}
        fi

        # The same for ITS2.
        if [ -s "$work_dir/ITS2.distinct.fasta" ]; then
            vsearch \
                --sintax "$work_dir/ITS2.distinct.fasta" \
                --db {input.its_db} \
                --tabbedout {output.its2_tax} \
                --sintax_cutoff {params.cutoff} \
                --randseed 1 \
                --threads {threads} >> {log} 2>&1
        else
            echo "No ITS2 sequence to classify for sample {wildcards.sample} (ITSx found none, or none passed the length guard). Skipping classification." >> {log} 2>&1
            touch {output.its2_tax}
        fi
        """
