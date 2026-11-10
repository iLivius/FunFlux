# FunFlux v2.0.0 — Stage 25: the contig summary, one table that says what every contig is
# (rules/shared/25_contig_summary.smk)
#
# A finished assembly is a list of contigs, and the list alone does not say what each
# contig is. A long-read assembly of a filamentous fungus typically holds, side by side,
# whole chromosomes (one contig from telomere to telomere), chromosome arms (the assembler
# stopped at a repeat, often the rDNA array or the centromere), the mitochondrial genome (a
# small circle present in many copies per cell), pieces of the ribosomal DNA (rDNA) array
# (tens to hundreds of tandem copies of the SSU-ITS1-5.8S-ITS2-LSU unit, which no assembler
# resolves completely), alternative copies of a region (a second haplotype, or a stretch
# assembled twice), low-depth debris and, now and then, a contaminant that the contamination
# screen removed; a short-read assembly adds collapsed repeats (several genome copies of a
# transposon or satellite folded into one deep contig) and short gene-less pieces of the
# nuclear genome by the hundred. The earlier stages already produce the evidence for all of
# this, spread over half a dozen files. The rules in this file gather it into one row per
# contig and give every contig a verdict with its reasons: chromosome, chromosome_arm,
# mitochondrion, rDNA_array, low_depth, repeat, nuclear_contig, nuclear_fragment,
# contaminant or unclassified. The same file holds the contig adapter screen, which looks for
# sequencing adapters and barcodes left in the contigs. Six rules:
#
#   telomere_db      once per run: downloads the telomere repeat database of tidk and checks
#                    its checksum.
#   univec_illumina  once per run: downloads NCBI's UniVec, checks its checksum and keeps the
#                    Illumina adapter records.
#   telomere_units   once per sample: works out the telomere repeat unit of the assembly
#                    with tidk, scores every candidate unit, and calls the telomere at both
#                    ends of every contig for the chosen unit.
#   blast_adapters   once per sample: searches the draft contigs for adapter and barcode
#                    sequences with BLAST.
#   screen_adapters  once per sample: turns the hits into adapter loci at contig ends and
#                    inside contigs, and writes the sample's row of the Adapters table of the
#                    MultiQC report.
#   contig_summary   once per sample: joins the telomere call and the adapter loci with the
#                    rest of the run's evidence and writes the table.
#
# The scripts that do the work, telomere_units.py, adapter_screen.py and contig_summary.py
# (all in scripts/25_contig_summary/), explain the method and the rules in their opening
# text; this file says what goes in and what comes out, and how the telomere step works in
# outline.
#
# The evidence, where it comes from, and what it says about a contig:
#
#   telomeres             the telomere is a tandem array of a short repeat unit (TTAGGG in
#                         most filamentous fungi) that caps every chromosome end: an array
#                         at both ends of a contig means a whole chromosome, at one end a
#                         chromosome arm. The unit differs between fungal groups, so rule
#                         telomere_units works it out from the assembly itself (the next
#                         section) and hands rule contig_summary tidk's call for both ends
#                         of every contig (TELOMERE_ENDS) and the unit (TELOMERE_UNIT_FILE).
#   read depth            the BlobTools table of the contamination screen (BLOB_TABLE), the
#                         sample's reads mapped to the draft: the nuclear contigs share one
#                         depth, the mitochondrion is several times deeper, an alternative
#                         copy of a region or noise is much shallower. Read-based modes only:
#                         in contigs mode the screen maps the contigs to THEMSELVES, and the
#                         depth it reports reflects contig length and nothing else
#                         (10_decontam.smk explains it at its contigs-mode rule map_contigs),
#                         so the table is not given to the script there and the depth columns
#                         read NA.
#   assembler depth and   ASSEMBLY_EVIDENCE, the assembler's own per-contig table (longread
#   circularity           and hybrid mode): a circular contig at a high depth is the usual
#                         shape of a mitochondrion.
#   rDNA units            ITS_WINDOW_TABLE, one row per 5.8S gene found by rule its_windows
#                         (30_taxonomy.smk): a contig made of rDNA units is a piece of the
#                         rDNA array, not a chromosome.
#   BUSCO markers         BUSCO_DIR (rule completeness, 20_qc.smk): the complete markers that
#                         lie on each contig, which tell genuine nuclear sequence from rDNA
#                         and debris.
#   genes                 the funannotate GFF3 of the final annotation (the mRNA features),
#                         the same kind of evidence. funannotate renames the contigs
#                         scaffold_1, scaffold_2, ... in order of length, so the GFF3 cannot
#                         be joined to the draft contigs by name: the script joins it by
#                         sequence, through CONTIGS_SORT, the renamed genome that
#                         funannotate_preprocess wrote (40_annotation.smk).
#   adapter loci          ADAPTER_SCREEN (rule screen_adapters, below): adapter or barcode
#                         sequence within 100 bp of a contig end, or a whole adapter inside a
#                         contig. It fills the adapter column and changes no verdict.
#   the two screens       CONTIG_DECISIONS (every contig kept or removed, with the genus the
#                         BLAST hits gave it; rule select_contigs, 10_decontam.smk),
#                         MITO_SCREEN_TSV (confirmed mitochondrial or not; rule
#                         screen_mitochondria, 15_mitochondria.smk) and MITO_CONTIGS (the
#                         mitochondrial contigs the split wrote; rule separate_mitochondria).
#                         Only in runs with the screen.
#
# How the telomere repeat unit is found (rule telomere_units, script telomere_units.py).
# tidk (Brown et al. 2025) does the two heavy parts. `tidk explore` lists every short unit
# (5 to 30 bp here) that forms a tandem run somewhere in a sequence set, with its number of
# copies. `tidk ends` takes one unit and reports, for both ends of every sequence, whether
# an array of it sits there, how long the array is, how far from the end it lies and which
# strand it is on. FunFlux adds the choice between units:
#
#   1. Contig ends. The first and the last telomere_window bases of every contig are
#      written to a FASTA (a contig shorter than twice the window goes in whole), and tidk
#      explore lists the tandem units found there; it runs once more on the whole draft, so
#      that the candidates table can show how many copies of a unit lie at the ends and how
#      many anywhere. Most units it lists are microsatellites (AAAGAG, AGATGG, ...), some
#      with nearly as many copies at the ends as the telomere unit, so the copy number
#      cannot choose.
#   2. Candidates. Every unit found at the ends, every unit of the tidk database with strong
#      or weak support (from every kingdom), every unit the literature names for a fungus (a
#      table in the script) and, when telomere_motif names one, the configured unit. All of
#      them in tidk's canonical form (the alphabetically smallest of the unit's rotations and
#      of its reverse complement: TTAGGG becomes AACCCT), so that one array is one candidate
#      whichever way it was read. Three filters, applied to every candidate but the
#      configured one: a length of 5 to 30 bp (no telomere repeat shorter than 5 bp is
#      known, and 4-bp units are low-complexity artefacts); not a repeat of a shorter unit
#      (ACACAC is AC); and a G-rich strand, that is, a count of G different from the count
#      of C (every telomerase-made repeat has one; AAAT has none, and a unit with as many G
#      as C reads the same on both strands).
#   3. Score. tidk ends runs once per candidate on the whole draft, with its window and the
#      farthest an array may lie from the end both set to telomere_window and the shortest
#      array it accepts set to min_telomere_copies x the unit length. A telomere repeat
#      points one way at every chromosome end (its G-rich strand runs towards the end), so
#      tidk reports an array as present when it points the way a telomere points and as
#      wrong_strand when it points the other way. The same-array rule comes first: when the
#      start and the end of a short contig both report an array and the two coordinate
#      ranges overlap, they are one array seen from both ends, and the row with the larger
#      distance from its end gets the FunFlux status other_end instead. Then telomere_ends
#      is the number of present rows, opposite_ends the number of wrong_strand rows, the
#      score is their difference, and t2t_contigs the contigs with present at both ends. A
#      satellite that happens to sit at contig ends points either way with equal chance, so
#      its score stays near zero however often it occurs; a telomere unit scores close to
#      its number of ends.
#   4. Choice. A candidate is eligible with at least 2 telomere_ends and more than twice as
#      many telomere_ends as opposite_ends. The eligible candidates are ranked by score,
#      then known (literature or database) before unknown, then the shorter unit, then
#      alphabetically; the first is the telomere unit of the assembly, and none is written
#      when no candidate is eligible. When telomere_motif names a unit, that unit is the one
#      written to the unit file and the ends table, whatever its score; the ranking still
#      runs, and a NOTE says when the data favour another unit or give the configured one
#      no support.
#
# What the step writes, in 02.assembly/{sample}/eval/ ({sample} is a Snakemake wildcard:
# the per-sample rules run once for every sample, with the sample name filled in):
#   telomere_candidates.tsv  one row per candidate found at the ends or forming at least
#                            one array: the unit, its length, its copies at the ends and
#                            genome-wide, telomere_ends, opposite_ends, score, t2t_contigs,
#                            what the literature and the database say about it, and
#                            whether it was chosen.
#   telomere_unit.txt        one line: the chosen unit in canonical form, or none.
#   telomere_ends.tsv        tidk's per-end table for the chosen unit after the same-array
#                            rule: contig, contig length, side (start or end), status
#                            (present, absent, wrong_strand, not_terminal or other_end), the
#                            array's coordinates and length, its distance from the end, its
#                            strand. A header line alone when the unit is none.
# The script's summary line and its NOTE, if any, start with "[FunFlux]"; the rule copies
# them from the log to the Snakemake console, where a user looks first.
#
# The database. tidk's telomere repeat database (telo_database_v2.tsv, about 600 KB: one
# row per species whose chromosome ends have been examined, with its telomere unit, the
# support for it and how many of its chromosome ends carried it) serves two things: its
# strong and weak units are candidates in every run, and the candidates table says how
# many species of each kingdom carry a unit. Its repository,
# tolkit/a-telomeric-repeat-database, carries no licence file, so FunFlux does not ship the
# file: rule telomere_db downloads it once per output directory from links.telomere_db_link
# (an address pinned to one revision of the repository) into 00.resources/telomeres/ and
# checks it against parameters.contig_summary.telomere_db_sha256 before anything reads it.
#
# Which of the evidence files a run has depends on the mode, and 00_common.smk (section 9h)
# works that out once: the five files of the screen and the assembler are listed as inputs
# of rule contig_summary only where the mode makes them, and the matching command-line
# option gets the file's path there and nothing in the other modes. The script treats an
# option without a value as "not given" and writes NA into the columns that need it; the
# verdict rules then skip the evidence they lack. The other inputs exist in every mode and
# are always given: the ITS windows, the BUSCO folder, the annotated genome with its GFF3,
# and the two telomere files of rule telomere_units.
#
# Both per-sample rules read the DRAFT assembly (DRAFT_CONTIGS: contigs_filt.fasta in
# contigs and illumina mode, assembly.fasta in longread and hybrid mode, and
# contigs_final.fasta in contigs mode without the screen), not the final genome: the
# contigs the screen removed and the mitochondrion get a row too, so nothing disappears
# silently.
#
# Where this file sits in the workflow. telomere_db needs nothing from other rules and runs
# early, like the other database downloads. telomere_units needs only the draft and the
# database, so it can run as soon as the assembly exists. contig_summary waits for the
# finished annotation (ANNOTATION_COMPLETE, rule funannotate_annotation in
# 40_annotation.smk), so it is one of the last jobs of a sample. The files of both
# per-sample rules go to the sample's QC folder, 02.assembly/{sample}/eval/, next to the
# QUAST and BUSCO results. Rule multiqc (90_report.smk) lists the _mqc copy of the table as
# an input and hands it to MultiQC, so the table appears in the report; rule all asks for
# the table itself (00_common.smk, section 15). No other rule reads the other files: they
# are for people.
#
# Everything referenced here is defined once in 00_common.smk and not worked out again
# here: TELOMERE_DB, TELOMERE_DB_PROVENANCE, TELOMERE_DB_URL, TELOMERE_DB_SHA256,
# DRAFT_CONTIGS, TELOMERE_CANDIDATES, TELOMERE_UNIT_FILE, TELOMERE_ENDS,
# TELOMERE_UNITS_SCRIPT, ITS_WINDOW_TABLE, BUSCO_DIR, CONTIGS_SORT, ANNOTATION_COMPLETE,
# CONTIG_SUMMARY_OPTIONAL_INPUTS, CONTIG_SUMMARY_MITO_CONTIGS, CONTIG_SUMMARY_DECISIONS,
# CONTIG_SUMMARY_MITO_SCREEN, CONTIG_SUMMARY_BLOB_TABLE, CONTIG_SUMMARY_ASSEMBLY_EVIDENCE,
# CONTIG_SUMMARY, CONTIG_SUMMARY_MQC, CONTIG_SUMMARY_SCRIPT, CONTIG_SUMMARY_TELOMERE_MOTIF,
# CONTIG_SUMMARY_TELOMERE_WINDOW, CONTIG_SUMMARY_MIN_TELOMERE_COPIES,
# CONTIG_SUMMARY_MIN_RDNA_UNITS, UNIVEC_URL, UNIVEC_SHA256, UNIVEC_ILLUMINA,
# UNIVEC_PROVENANCE, ADAPTER_LIST_FILES, ADAPTER_QUERY, ADAPTER_BLAST, ADAPTER_SCREEN,
# ADAPTER_SCREEN_MQC, ADAPTER_SCREEN_SCRIPT and LOGS, plus the functions funannotate_gff3()
# (the GFF3 path of a sample) and trimmed_adapter_inputs(), trimmed_adapter_sequence(),
# trimmed_adapter_choice() and trimmed_adapter_check() (the read-level files of a trimmed ONT
# sample, or nothing). This file defines no name of its own.
#
# conda: environment paths are written relative to THIS file (workflow/rules/shared/), so
# "../../envs/tidk.yaml" climbs shared/ -> rules/ -> workflow/ and then goes into envs/.
# Five environments. blast_adapters runs in the BLAST environment of the contamination screen
# (envs/blast.yaml) and screen_adapters in envs/python.yaml (a pinned Python 3.12); univec_illumina
# downloads with envs/download.yaml. telomere_units runs in the tidk environment (envs/tidk.yaml), which
# holds tidk and a pinned Python 3.12, so the script and the tool it calls come from one
# place. contig_summary uses only Python's standard library and runs in the vsearch
# environment, which holds a pinned Python 3.12 (added for the ITS scripts of
# 30_taxonomy.smk, which run with it for the same reason: activating the environment makes
# `python` that Python, instead of whatever `python` the computer has). telomere_db runs in
# envs/download.yaml for its curl, like the other download rules (tiberius_model,
# chroquetas_db, ...); its sha256sum comes from the computer.


# ───────────── telomere_db — download and check the telomere repeat database ─────────────
# Takes in:    nothing from other rules; the address links.telomere_db_link
#              (TELOMERE_DB_URL) and the checksum parameters.contig_summary.telomere_db_sha256
#              (TELOMERE_DB_SHA256).
# Does:        downloads the tab-separated file (about 600 KB) with curl, checks its SHA-256
#              checksum, and only then writes the provenance file. The address names one
#              revision of the database's repository, so the file at it should never change;
#              a file that fails the check is a changed or cut-off download, and the run
#              stops with a message instead of choosing telomere units from a file nobody
#              has looked at. The failed file is renamed with a .failed suffix and kept, so
#              that it can be looked at (Snakemake deletes the outputs of a failed job, but
#              not a file under another name).
# Produces:    TELOMERE_DB, the database, and TELOMERE_DB_PROVENANCE, one item per line:
#              the address, the date in UTC, the checksum and the number of rows without
#              the header line. Both in 00.resources/telomeres/. They stay, so a later run
#              into the same output folder does not download the file again.
# Consumed by: telomere_units, for every sample.
#
# Environment: envs/download.yaml, for curl (a minimal Linux installation may have none);
# sha256sum, date and wc come from the computer (GNU coreutils, part of every Linux base
# system).
# Where it runs: always on the computer that runs Snakemake itself (localrule: True). On a
# single computer that changes nothing; on a cluster it keeps the download off the compute
# nodes, which often have no internet.
#
# How the shell works:
#   * curl -L follows a redirect; --fail turns an HTTP error into a failed command instead
#     of a saved error page, and a failed command stops the rule (the shell runs with
#     `set -euo pipefail`).
#   * sha256sum -c - reads a line "<checksum>  <file>" (two spaces) from its input, here
#     the checksum from the configuration and the downloaded file, and writes OK or FAILED
#     for that file. If the check fails, the "then" branch renames the file, prints a
#     message that says what to do and ends the rule with `exit 1`, so the provenance file
#     is never written.
#   * wc -l counts the lines of the file; the header line is taken off for the row count.
rule telomere_db:
    localrule: True
    output:
        db = TELOMERE_DB,
        provenance = TELOMERE_DB_PROVENANCE,
    params:
        url = TELOMERE_DB_URL,
        sha256 = TELOMERE_DB_SHA256,
    conda:
        "../../envs/download.yaml"
    log:
        LOGS + "/telomere_db.log"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.db})"

        echo "[telomere_db] downloading {params.url}" > {log}
        curl -L --fail -o "{output.db}" "{params.url}" >> {log} 2>&1

        echo "[telomere_db] checking the file against parameters.contig_summary.telomere_db_sha256 ({params.sha256})" >> {log}
        if ! echo "{params.sha256}  {output.db}" | sha256sum -c - >> {log} 2>&1; then
            mv "{output.db}" "{output.db}.failed"
            echo "[FunFlux] telomere_db: the file downloaded from {params.url} does not match parameters.contig_summary.telomere_db_sha256, so it was NOT used (details in {log}; the file is kept as {output.db}.failed). The address is pinned to one revision of the database, so a mismatch means a changed or cut-off download. Look at the kept file; to accept it, put its checksum (sha256sum {output.db}.failed) into parameters.contig_summary.telomere_db_sha256." | tee -a {log} >&2
            exit 1
        fi

        echo "source_url: {params.url}" > {output.provenance}
        echo "downloaded_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> {output.provenance}
        echo "sha256: {params.sha256}" >> {output.provenance}
        echo "rows: $(($(wc -l < "{output.db}") - 1))" >> {output.provenance}
        echo "[telomere_db] done: $(($(wc -l < "{output.db}") - 1)) rows" >> {log}
        """


# ───────────── univec_illumina — download UniVec and keep its Illumina records ─────────────
# Takes in:    nothing from other rules; the address links.univec_link (UNIVEC_URL, by default
#              build 10.1 in NCBI's archive) and the checksum
#              parameters.contig_summary.univec_sha256 (UNIVEC_SHA256).
# Does:        downloads UniVec (about 1.7 MB) with curl, checks its SHA-256 checksum, keeps
#              the records whose description names Illumina (151 in build 10.1: the TruSeq,
#              Nextera and other adapters and primers) and deletes the rest. A file that fails
#              the check is renamed with a .failed suffix and kept, and the run stops, as in
#              rule telomere_db.
# Produces:    UNIVEC_ILLUMINA, the Illumina records (FASTA), and UNIVEC_PROVENANCE, one item
#              per line: the address, the date in UTC, the checksum and the number of records.
#              Both in 00.resources/adapters/, kept.
# Consumed by: blast_adapters, for every sample.
#
# Why downloaded and not shipped: UniVec notes on every Illumina record that the
# oligonucleotide sequence is copyrighted by Illumina. The ONT, PacBio and MGI adapters are
# shipped in workflow/data/adapters/ (the README there names their sources).
#
# Environment and where it runs: as rule telomere_db (curl from envs/download.yaml; sha256sum,
# awk, grep and date from the computer; on the computer that runs Snakemake).
#
# How the shell works: the download goes to a file of its own (params download); after the
# check, awk copies a record when its header line (starting with >) contains "Illumina": the
# variable keep is set at every header line and every line is printed while it is true. awk also
# renames the record, from UniVec's gnl|uv|NGB00747.1:1-70 to Illumina_UniVec_NGB00747.1:1-70, so
# that the adapter screen's tables say at a glance that a hit is an Illumina adapter (the
# bracket expressions [|] stand for a literal |). The doubled braces are Snakemake's way of
# writing single braces in a shell command. grep -c ends with status 1 when it counts nothing,
# which the strict shell would take for a failure; "|| true" keeps the count of 0, and a NOTE then
# says that the screen searches for the shipped adapters only.
rule univec_illumina:
    localrule: True
    output:
        fasta = UNIVEC_ILLUMINA,
        provenance = UNIVEC_PROVENANCE,
    params:
        url = UNIVEC_URL,
        sha256 = UNIVEC_SHA256,
        download = UNIVEC_ILLUMINA + ".download",
    conda:
        "../../envs/download.yaml"
    log:
        LOGS + "/univec_illumina.log"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.fasta})"

        echo "[univec_illumina] downloading {params.url}" > {log}
        curl -L --fail -o "{params.download}" "{params.url}" >> {log} 2>&1

        echo "[univec_illumina] checking the file against parameters.contig_summary.univec_sha256 ({params.sha256})" >> {log}
        if ! echo "{params.sha256}  {params.download}" | sha256sum -c - >> {log} 2>&1; then
            mv "{params.download}" "{params.download}.failed"
            echo "[FunFlux] univec_illumina: the file downloaded from {params.url} does not match parameters.contig_summary.univec_sha256, so it was NOT used (details in {log}; the file is kept as {params.download}.failed). To accept it, put its checksum (sha256sum {params.download}.failed) into parameters.contig_summary.univec_sha256." | tee -a {log} >&2
            exit 1
        fi

        awk '/^>/ {{ keep = /Illumina/; if (keep) sub(/^>gnl[|]uv[|]/, ">Illumina_UniVec_") }} keep' "{params.download}" > "{output.fasta}"
        rm "{params.download}"
        records=$(grep -c '^>' "{output.fasta}" || true)
        if [ "$records" = 0 ]; then
            echo "[FunFlux] NOTE univec_illumina: no record of {params.url} names Illumina, so the adapter screen searches for the adapters shipped with FunFlux (ONT, PacBio, MGI) only." | tee -a {log} >&2
        fi

        echo "source_url: {params.url}" > {output.provenance}
        echo "downloaded_utc: $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> {output.provenance}
        echo "sha256: {params.sha256}" >> {output.provenance}
        echo "illumina_records: $records" >> {output.provenance}
        echo "[univec_illumina] done: $records Illumina records" >> {log}
        """


# ───────────── telomere_units — the telomere repeat unit of one assembly ─────────────
# Takes in:    DRAFT_CONTIGS, the draft assembly (made by the mode's front end: rule
#              filter_contigs in contigs and illumina mode, rule long_read_assembly in
#              longread and hybrid mode), and TELOMERE_DB from rule telomere_db.
# Does:        python telomere_units.py (TELOMERE_UNITS_SCRIPT), once: the method at the top
#              of this file. tidk itself is called by the script, several times (explore
#              twice, ends once per candidate), from the same conda environment. The three
#              settings of parameters.contig_summary that the search uses come in through
#              params: the motif (auto, or the unit to use as given), the window and the
#              copies. A changed setting re-runs this step, contig_summary and MultiQC, and
#              nothing else.
#              Everything the script prints goes to the log; the shell then copies the lines
#              that start with [FunFlux] (the summary of the choice and, if any, the NOTE)
#              to the console (the grep: -F looks for the fixed text, and the "|| true" is
#              needed because grep ends with status 1 when it finds nothing, which under the
#              strict shell would fail the job).
# Produces:    TELOMERE_CANDIDATES  02.assembly/{sample}/eval/telomere_candidates.tsv: every
#                                   candidate unit with its counts and score, the chosen
#                                   one marked (described at the top of this file).
#              TELOMERE_UNIT_FILE   02.assembly/{sample}/eval/telomere_unit.txt: the chosen
#                                   unit, or none.
#              TELOMERE_ENDS        02.assembly/{sample}/eval/telomere_ends.tsv: tidk's call
#                                   at both ends of every contig for that unit.
# Consumed by: rule contig_summary reads the unit file and the ends table. The candidates
#              table is for people.
#
# One thread: tidk explore takes a second or two on a 40 Mb genome, each tidk ends call
# well under a second, and the script runs for seconds to a minute in all. A failure stops
# the run, as for the other per-sample steps; the script exits with a readable message on a
# missing file, a malformed database or a tidk error, so the log says what to look at.
#
# New in v2. There is no v1 counterpart.
rule telomere_units:
    input:
        contigs = DRAFT_CONTIGS,
        telomere_db = TELOMERE_DB,
    output:
        candidates = TELOMERE_CANDIDATES,
        unit = TELOMERE_UNIT_FILE,
        ends = TELOMERE_ENDS,
    params:
        script = TELOMERE_UNITS_SCRIPT,
        telomere_motif = CONTIG_SUMMARY_TELOMERE_MOTIF,
        telomere_window = CONTIG_SUMMARY_TELOMERE_WINDOW,
        min_telomere_copies = CONTIG_SUMMARY_MIN_TELOMERE_COPIES,
    conda:
        "../../envs/tidk.yaml"
    log:
        LOGS + "/telomere_units_{sample}.log"
    shell:
        """
        set -euo pipefail

        python {params.script:q} \
            --contigs {input.contigs:q} \
            --telomere-db {input.telomere_db:q} \
            --telomere-motif {params.telomere_motif:q} \
            --telomere-window {params.telomere_window} \
            --min-telomere-copies {params.min_telomere_copies} \
            --sample {wildcards.sample:q} \
            --out-candidates {output.candidates:q} \
            --out-unit {output.unit:q} \
            --out-ends {output.ends:q} > {log} 2>&1

        # Copy the script's summary line and its NOTE, if there is one, to the console.
        grep -F '[FunFlux]' {log} >&2 || true
        """


# ───────────── blast_adapters — search the draft contigs for adapters and barcodes ─────────────
# Takes in:    DRAFT_CONTIGS, the draft assembly (see rule telomere_units); the shipped adapter
#              lists (ADAPTER_LIST_FILES: ONT adapters and barcodes, PacBio and MGI adapters);
#              UNIVEC_ILLUMINA from rule univec_illumina; and, for an ONT sample whose reads
#              were trimmed, ONT_ADAPTER_SEQUENCE, the very sequences fastplong removed (rule
#              find_ont_adapter, 05_ont_adapters.smk). That last file catches a barcode that no
#              published list holds; trimmed_adapter_sequence() gives it as a list of one file,
#              or an empty list for any other sample.
# Does:        joins the lists into one FASTA (ADAPTER_QUERY) with cat: an empty list of trimmed
#              sequences is written as nothing at all, so cat then reads only the other files.
#              It then runs blastn with the adapters
#              as queries and the contigs as subjects: -task blastn-short is BLAST's setting for
#              queries shorter than 50 bp (word size 7), and -evalue 1 keeps short matches,
#              which adapter_screen.py then judges by length, identity and position. The search
#              takes about half a minute on a 40 Mb genome.
# Produces:    ADAPTER_QUERY (temporary) and ADAPTER_BLAST, the hits in BLAST's tabular format
#              with the columns named after -outfmt 6 (kept).
# Consumed by: screen_adapters.
rule blast_adapters:
    input:
        contigs = DRAFT_CONTIGS,
        shipped = ADAPTER_LIST_FILES,
        illumina = UNIVEC_ILLUMINA,
        trimmed = trimmed_adapter_sequence,
    output:
        query = temp(ADAPTER_QUERY),
        hits = ADAPTER_BLAST,
    conda:
        "../../envs/blast.yaml"
    log:
        LOGS + "/blast_adapters_{sample}.log"
    shell:
        """
        set -euo pipefail

        cat {input.shipped:q} {input.illumina:q} {input.trimmed:q} > {output.query:q}

        blastn -task blastn-short \
            -query {output.query:q} \
            -subject {input.contigs:q} \
            -evalue 1 \
            -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen" \
            -out {output.hits:q} 2> {log}
        """


# ───────────── screen_adapters — adapter loci at contig ends and inside contigs ─────────────
# Takes in:    ADAPTER_BLAST from blast_adapters, DRAFT_CONTIGS (for the contig lengths) and,
#              for a trimmed ONT sample, the choice and check files of 05_ont_adapters.smk
#              (trimmed_adapter_choice() and trimmed_adapter_check() give their paths, or
#              nothing).
# Does:        python adapter_screen.py (ADAPTER_SCREEN_SCRIPT): keeps the hits of at least
#              90 % identity that are at least 20 bp long within 100 bp of a contig end, or at
#              least 28 bp long (a whole ONT adapter) inside a contig, merges overlapping hits
#              into loci, and writes the report. The script's opening text gives the reasons
#              for these rules. Its summary line is copied from the log to the console.
# Produces:    ADAPTER_SCREEN, one row per locus (a header line alone when there is none), and
#              ADAPTER_SCREEN_MQC, the sample's row of the Adapters table of the MultiQC
#              report. Both kept.
# Consumed by: contig_summary (the adapter column) and multiqc (90_report.smk). Rule all asks
#              for ADAPTER_SCREEN.
rule screen_adapters:
    input:
        hits = ADAPTER_BLAST,
        contigs = DRAFT_CONTIGS,
        read_files = trimmed_adapter_inputs,
    output:
        loci = ADAPTER_SCREEN,
        mqc = ADAPTER_SCREEN_MQC,
    params:
        script = ADAPTER_SCREEN_SCRIPT,
        read_choice = trimmed_adapter_choice,
        read_check = trimmed_adapter_check,
    conda:
        "../../envs/python.yaml"
    log:
        LOGS + "/screen_adapters_{sample}.log"
    shell:
        """
        set -euo pipefail

        python {params.script:q} \
            --blast {input.hits:q} \
            --contigs {input.contigs:q} \
            --sample {wildcards.sample:q} \
            --read-choice {params.read_choice:q} \
            --read-check {params.read_check:q} \
            --out {output.loci:q} \
            --out-mqc {output.mqc:q} > {log} 2>&1

        grep -F '[FunFlux]' {log} >&2 || true
        """


# ───────────── contig_summary — one row per draft contig, with a verdict ─────────────
# Takes in:    DRAFT_CONTIGS, the draft assembly (see rule telomere_units); TELOMERE_UNIT_FILE
#              and TELOMERE_ENDS (rule telomere_units: the unit, and tidk's call at both
#              ends of every contig); ITS_WINDOW_TABLE (rule its_windows); BUSCO_DIR (rule
#              completeness; a whole folder, the script finds the full_table.tsv of the
#              lineage BUSCO settled on); CONTIGS_SORT (rule funannotate_preprocess);
#              ANNOTATION_COMPLETE, the marker that the final annotation is written (rule
#              funannotate_annotation; the GFF3 itself lies inside FUNANNOTATE_DIR, a
#              folder another rule declares as its output, so the marker is the input and
#              the GFF3 path goes to the script through params); and
#              CONTIG_SUMMARY_OPTIONAL_INPUTS, the files of the screen and the assembler
#              that THIS mode makes (see the top of this file), so the job waits for them;
#              and ADAPTER_SCREEN (rule screen_adapters), for the adapter column.
# Does:        python contig_summary.py (CONTIG_SUMMARY_SCRIPT), once. The optional files
#              are passed through params: in the modes without one, its param is empty
#              and the :q leaves the option without a value, which the script reads as
#              "not given" (the way the policy options of rule select_contigs in
#              10_decontam.smk work). The two settings of parameters.contig_summary that
#              the verdicts use, the copies and the rDNA units, come in through params too,
#              so a changed setting re-runs this step and MultiQC and nothing else.
#              Everything the script prints (its one-line summary) goes to the log.
# Produces:    CONTIG_SUMMARY         02.assembly/{sample}/eval/contig_summary.tsv: one row
#                                     per draft contig (length, GC, the two depths,
#                                     circularity, the telomere copies at each end and
#                                     their distance from the contig end, the telomere
#                                     unit of the sample, rDNA units, BUSCO markers, genes,
#                                     genus, the screens' verdicts, the adapter loci
#                                     (none, at_end, internal), the scaffold name the
#                                     annotation gave the contig, the verdict and its
#                                     reason). A column a mode cannot fill reads NA.
#              CONTIG_SUMMARY_MQC     the same table with MultiQC's custom-content header
#                                     lines, the rows named "<sample>: <contig>", so that
#                                     the tables of all samples merge into one in the
#                                     report.
# Consumed by: rule multiqc (90_report.smk) reads the _mqc copy; rule all asks for the
#              table.
#
# One thread: the script runs for about ten seconds on a 40 Mb genome. A failure stops the
# run, as for the other per-sample steps; the script exits with a readable message on a
# malformed table, so the log says what to look at.
#
# New in v2. There is no v1 counterpart: the v1 workflows left the reader to compare the
# QUAST, BUSCO, BlobTools, mitochondrion and ITS results by hand.
rule contig_summary:
    input:
        contigs = DRAFT_CONTIGS,
        telomere_unit = TELOMERE_UNIT_FILE,
        telomere_ends = TELOMERE_ENDS,
        its_windows = ITS_WINDOW_TABLE,
        busco_dir = BUSCO_DIR,
        annotated_fasta = CONTIGS_SORT,
        annotation_complete = ANNOTATION_COMPLETE,
        optional_files = CONTIG_SUMMARY_OPTIONAL_INPUTS,
        adapter_screen = ADAPTER_SCREEN,
    output:
        table = CONTIG_SUMMARY,
        mqc = CONTIG_SUMMARY_MQC,
    params:
        script = CONTIG_SUMMARY_SCRIPT,
        mito_contigs = CONTIG_SUMMARY_MITO_CONTIGS,
        decisions = CONTIG_SUMMARY_DECISIONS,
        mito_screen = CONTIG_SUMMARY_MITO_SCREEN,
        blob_table = CONTIG_SUMMARY_BLOB_TABLE,
        assembly_evidence = CONTIG_SUMMARY_ASSEMBLY_EVIDENCE,
        gff3 = funannotate_gff3,
        min_telomere_copies = CONTIG_SUMMARY_MIN_TELOMERE_COPIES,
        min_rdna_units = CONTIG_SUMMARY_MIN_RDNA_UNITS,
    conda:
        "../../envs/vsearch.yaml"
    log:
        LOGS + "/contig_summary_{sample}.log"
    shell:
        """
        set -euo pipefail

        python {params.script:q} \
            --contigs {input.contigs:q} \
            --mito-contigs {params.mito_contigs:q} \
            --decisions {params.decisions:q} \
            --mito-screen {params.mito_screen:q} \
            --blob-table {params.blob_table:q} \
            --assembly-evidence {params.assembly_evidence:q} \
            --its-windows {input.its_windows:q} \
            --busco-dir {input.busco_dir:q} \
            --gff3 {params.gff3:q} \
            --annotated-fasta {input.annotated_fasta:q} \
            --telomere-unit-file {input.telomere_unit:q} \
            --telomere-ends {input.telomere_ends:q} \
            --adapter-screen {input.adapter_screen:q} \
            --min-telomere-copies {params.min_telomere_copies} \
            --min-rdna-units {params.min_rdna_units} \
            --sample {wildcards.sample:q} \
            --out {output.table:q} \
            --out-mqc {output.mqc:q} > {log} 2>&1
        """
