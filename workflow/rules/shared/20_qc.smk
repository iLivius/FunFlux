# FunFlux v2.0.0 — Stage 20: assembly quality control (rules/shared/20_qc.smk)
#
# Two checks of the finished NUCLEAR genome (FINAL_CONTIGS, contigs_final.fasta), one tool
# each:
#
#   QUAST — how good the ASSEMBLY is: how many contigs it took to represent the genome, the
#           largest contig, the total length, the GC content, and N50 / N90 (the contig
#           length such that half / 90 % of the assembly lies in contigs at least that long;
#           bigger is better). It measures how fragmented the assembly is, and fragmentation
#           limits every later step that depends on gene context. It says nothing about
#           whether the genes are all there.
#   BUSCO — how complete the GENOME is. BUSCO looks for a set of genes that nearly every
#           member of a lineage has exactly once ("single-copy orthologs"). A gene found
#           whole and once is "complete"; found in two or more copies, "duplicated"; found
#           in pieces, "fragmented"; not found, "missing". A well-assembled fungal genome is
#           typically more than 95 % complete.
#
# A third check, whether the reads map back sensibly, only makes sense in the modes that
# have reads:
#
#   Qualimap — how the sample's own reads fit the DRAFT they were assembled into: the share
#              of reads that map, the mean depth and how evenly it is spread, the insert
#              sizes of paired reads and the GC content. A low mapping rate or a very uneven
#              depth points at contamination or at an assembly that misses part of the
#              genome. The BAM it reads is the one the contamination screen made.
#
# QUAST and BUSCO exist in every run, with or without the screen; the Qualimap rule only
# when the mode has reads (`if HAS_READS:`, decided once when Snakemake reads the workflow).
# The three are independent of each other. Two other rules read their results: the contig
# summary (rule contig_summary in shared/25_contig_summary.smk) reads BUSCO's folder, and the
# MultiQC report (rule multiqc in shared/90_report.smk) shows all three side by side.
# Rule all asks for each result folder the mode makes (00_common.smk, section 15).
#
#   genome_assembly_evaluation : QUAST on contigs_final.fasta (in hybrid mode also on the
#                                screened Illumina guide assembly, in the same call)
#   busco_data                 : downloads the BUSCO lineage data once per run, only when
#                                parameters.busco.lineage names a lineage
#   completeness               : BUSCO on contigs_final.fasta
#   map_evaluation             : Qualimap on the reads mapped to the draft (read modes)
#
# Why contigs_final.fasta: it is the NUCLEAR genome, the same file that ITS extraction and
# funannotate read. When the contamination and mitochondrion screen is on and
# parameters.mitochondria.mode is separate_confirmed (the default), the mitochondrial genome
# has been split off before this point (rule separate_mitochondria); with report or off it
# stays in. Without the screen (contigs mode with parameters.contigs.screen: false, its
# default) contigs_final.fasta is the supplied FASTA after rule filter_contigs,
# mitochondrion included if it was there. The mitochondrion has a different genetic code,
# other genes and a very different copy number, so keeping it in would distort both the
# contig statistics and the BUSCO counts of a genome we are about to annotate as a nuclear
# genome.
#
# These names come from 00_common.smk, where each is defined once and never worked out
# again: FINAL_CONTIGS, GUIDE_DECONTAM, DECONTAM_BAM, QUAST_DIR, BUSCO_DIR, QUALIMAP_DIR,
# QC_LABELS, BUSCO_LINEAGE, BUSCO_LINEAGE_ARGS, BUSCO_DATA_DIR, BUSCO_DATA_READY,
# BUSCO_DATA_INPUT, BUSCO_DOWNLOAD_ARGS, IS_HYBRID, HAS_READS, RAM, HEAVY_MEM_MB,
# LONG_RUNTIME_MINUTES and LOGS. So is the function capped_cpus().
#
# conda: environment paths are written relative to THIS file (workflow/rules/shared/), so
# "../../envs/quast.yaml" goes up two folders, to workflow/, and then into workflow/envs/.
#
# In the shell commands below, a doubled brace ({{ or }}) stands for a single brace in the
# shell. Snakemake reads single braces as its own placeholders, such as {input.contigs}, and
# turns a doubled brace into a single one.


# ───────────── quast_labels: the name QUAST prints for each genome ─────────────
# Takes in:    the wildcards of one QUAST job (wildcards are the parts of a file name that
#              change from job to job; here only {sample}), and QC_LABELS from 00_common.smk:
#              one label per genome QUAST is given, each still containing the {sample}
#              placeholder. One genome gives ["{sample}"]. In hybrid mode QC_LABELS holds
#              two entries, ["{sample}_guide", "{sample}"], one per genome in _QUAST_GENOMES
#              below, in the same order.
# Does:        fills in this job's sample name and joins the labels with commas, which is
#              how QUAST wants its --labels list. In contigs mode the result is simply the
#              sample name.
# Produces:    the text for QUAST's --labels option.
# Called by:   Snakemake, once per job of rule genome_assembly_evaluation, because that
#              rule's params: labels names this function. The label is also the sample's
#              name in QUAST's report.tsv, which is how MultiQC finds the sample again (the
#              rename lines in shared/90_report.smk are written for it).
def quast_labels(wildcards):
    labels = []
    for label_template in QC_LABELS:
        labels.append(label_template.format(sample=wildcards.sample))
    return ",".join(labels)


# ── The genomes QUAST is given, in the order of QC_LABELS ──
# One in every mode but hybrid: the final nuclear genome. In hybrid mode two: the screened
# Illumina "guide" assembly first (GUIDE_DECONTAM, written by rule select_guide_contigs in
# rules/hybrid/25_guide.smk), then the final genome. The guide is the short-read view of the
# same isolate, so its statistics next to the final genome's show what the long reads added
# (fewer, longer contigs) and whether the two agree on the total length. The order matters:
# QUAST pairs its --labels with the input files in the order given, and QC_LABELS lists
# "{sample}_guide" first.
if IS_HYBRID:
    _QUAST_GENOMES = [GUIDE_DECONTAM, FINAL_CONTIGS]
else:
    _QUAST_GENOMES = [FINAL_CONTIGS]


# ───────────── genome_assembly_evaluation: QUAST statistics of the assembly ─────────────
# Takes in:    _QUAST_GENOMES: FINAL_CONTIGS (contigs_final.fasta), the nuclear genome that
#              the mode's front end made: rule filter_contigs in contigs mode without the
#              screen, rule separate_mitochondria when the screen is on. In hybrid mode the
#              screened guide assembly comes first (see above).
# Does:        creates the output folder, then runs QUAST on the genome(s). --labels gives
#              each assembly its name in the report (see quast_labels above): the sample's
#              name, and in hybrid mode {sample}_guide for the guide.
#              --no-icarus leaves out QUAST's interactive contig viewer: nobody opens it in
#              a batch run and it only adds files.
# Produces:    QUAST_DIR (02.assembly/{sample}/eval/quast/), declared as a directory()
#              output, so Snakemake treats the whole folder as the result; QUAST names the
#              files in it (report.tsv, report.html, ...).
# Consumed by: MultiQC (rule multiqc in shared/90_report.smk), which reads report.tsv.
#
# QUAST is only a report. No later step needs its result, so a QUAST crash must not stop a
# batch that may have spent days on assembly and annotation. The output folder is created
# first (so the rule's declared output exists whatever happens), QUAST is allowed to fail
# (the part after || runs only then), and a [FunFlux] line in the log says what happened.
# The job still counts as successful, and MultiQC simply has no QUAST table for that sample.
#
# v1 -> v2: (1) the CPU count is Snakemake's own threads: (capped at 24), which --cores
# enforces; v1 used a private "resources: cpus" that --cores did not count. (2) The
# labels come from QC_LABELS; for one genome that is the sample name, exactly as in v1.
# (3) The failure guard is new. docs/getting-started/from-v1.md names
# QUAST, Qualimap and NanoPlot as the QC rules that no longer stop the run when their tool
# fails; BUSCO is not on that list (see the note at rule completeness).
# (4) The log is named after the rule (v1: assembly_evaluation_{sample}.log). (5) In hybrid
# mode FunFluxL+ ran QUAST twice, once on the guide (label {sample}_illumina) and once on
# the final genome ({sample}_hybrid); one call with two labels replaces both, so the two
# assemblies stand side by side in one table.
rule genome_assembly_evaluation:
    input:
        genomes = _QUAST_GENOMES,
    output:
        quast_dir = directory(QUAST_DIR),
    params:
        labels = quast_labels,
    threads: capped_cpus(24)
    conda:
        "../../envs/quast.yaml"
    log:
        LOGS + "/genome_assembly_evaluation_{sample}.log"
    shell:
        """
        mkdir -p {output.quast_dir}

        quast \
            {input.genomes} \
            -o {output.quast_dir} \
            --labels {params.labels:q} \
            --no-icarus \
            -t {threads} > {log} 2>&1 || {{
            echo "" >> {log}
            echo "[FunFlux] genome_assembly_evaluation: QUAST exited with an error for sample {wildcards.sample}. QUAST is only a report that no other step reads, so the run continues; the assembly table will be missing from the MultiQC page for this sample. See the QUAST messages above." >> {log}
        }}
        """


# ───────────── busco_data: the BUSCO lineage data, once per run ─────────────
# Defined only when parameters.busco.lineage names a lineage (BUSCO_LINEAGE is not auto;
# decided once, when Snakemake reads the workflow). In auto mode each BUSCO job downloads
# what it needs itself (see rule completeness), because the lineage is known only once
# BUSCO has placed the genome.
#
# Takes in:    nothing from other rules; the lineage name (BUSCO_LINEAGE, for example
#              hypocreales_odb12).
# Does:        busco --download <lineage> --download_path BUSCO_DATA_DIR: BUSCO fetches the
#              newest release of that lineage from its server (busco-data.ezlab.org), checks
#              it against the checksum the server lists, and unpacks it into
#              BUSCO_DATA_DIR/lineages/<lineage>/. BUSCO writes its own messages to the log.
#              A name the server does not know makes BUSCO print an error but end without
#              one, so the shell then checks that the lineage folder's description file,
#              dataset.cfg, is there, and otherwise stops with a [FunFlux] message.
# Produces:    BUSCO_DATA_READY (00.resources/busco/lineages/<lineage>/dataset.cfg), the
#              file Snakemake watches, and the rest of the lineage folder next to it (a few
#              hundred MB unpacked), plus BUSCO's list of dataset versions, file_versions.tsv.
#              All of it stays, so a later run into the same output folder uses the same
#              release of the lineage and downloads nothing.
# Consumed by: completeness, for every sample (BUSCO_DATA_INPUT), which reads the folder
#              with --offline.
#
# Where it runs: always on the computer that runs Snakemake itself (localrule: True). On a
# single computer that changes nothing; on a cluster it keeps the download off the compute
# nodes, which often have no internet, and the BUSCO jobs on those nodes need none.
#
# New in v2: v1 downloaded the lineage data in every BUSCO job.
if BUSCO_LINEAGE != "auto":

    rule busco_data:
        localrule: True
        output:
            ready = BUSCO_DATA_READY,
        params:
            lineage = BUSCO_LINEAGE,
            data_dir = BUSCO_DATA_DIR,
        conda:
            "../../envs/busco.yaml"
        log:
            LOGS + "/busco_data.log"
        shell:
            """
            busco --download {params.lineage} --download_path {params.data_dir} > {log} 2>&1

            if [ ! -s {output.ready} ]; then
                echo "[FunFlux] busco_data: BUSCO did not download the lineage {params.lineage} (details in {log}). Check parameters.busco.lineage against the list that busco --list-datasets prints, or set it to auto." | tee -a {log} >&2
                exit 1
            fi
            """


# ───────────── completeness: BUSCO, how many of the expected genes are there ─────────────
# Takes in:    FINAL_CONTIGS (contigs_final.fasta), the same nuclear genome as QUAST.
# Does:        runs BUSCO in genome mode. The lineage it tests against is set by
#              parameters.busco.lineage, which 00_common.smk turns into BUSCO options
#              (BUSCO_LINEAGE_ARGS):
#                auto (the default)  --auto-lineage-euk: BUSCO places the genome itself,
#                                    first in a broad eukaryote lineage, then in the most
#                                    specific lineage it can reach (for a hypocrealean
#                                    fungus such as Fusarium, eukaryota_odb12 and then
#                                    hypocreales_odb12). It writes TWO summaries, a
#                                    "generic" one for the broad lineage and a "specific"
#                                    one for the precise lineage.
#                a lineage name      for example hypocreales_odb12, passed with -l: BUSCO
#                                    tests against that one lineage and writes only the
#                                    "specific" summary.
#              Where the lineage data (the gene sets) come from (BUSCO_DOWNLOAD_ARGS,
#              00_common.smk, section 9c):
#                auto                the job downloads them, for EVERY SAMPLE, into the
#                                    sample's own folder busco_tmp inside BUSCO_DIR, which
#                                    the shell deletes when BUSCO is done. That is a few
#                                    hundred MB of archives: the specific lineage, the broad
#                                    eukaryote lineage and the placement files, with which
#                                    the automatic lineage choice places the genome's marker
#                                    genes on a reference tree (mostly a reference protein
#                                    alignment; the tree and the marker lists are small).
#                                    BUSCO unpacks each archive and deletes it; unpacked, the
#                                    data can take more than a GB in busco_tmp while the job
#                                    runs. BUSCO 6.0.0 has no lock on its download folder and
#                                    silently updates a shared one, so two samples that
#                                    downloaded into one shared folder at the same time could
#                                    corrupt it. Never share that folder. The computer that
#                                    runs the job needs internet.
#                a lineage name      rule busco_data (below) has downloaded the lineage once,
#                                    into 00.resources/busco/, and BUSCO reads it with
#                                    --offline: it downloads nothing and writes nothing
#                                    there (it unpacks what it needs into the sample's own
#                                    result folder), so all samples can read one copy, and
#                                    the job needs no internet.
#              --opt-out-run-stats: without it, BUSCO 6.0.0 ends every run, a failed one too, by
#              looking up the country of the computer's IP address at ipinfo.io and sending
#              usage statistics to busco-data.ezlab.org: the tool versions, how BUSCO was
#              installed (conda, Docker or Singularity), the address it downloads from, the
#              options (without the input and output names; other paths replaced by
#              "user-specified"), the datasets, run times and results, the MD5 checksum and size
#              of the input FASTA, that country, a hash of the IP address and, for a failed run,
#              the error message (BUSCO replaces the paths in it only for configuration errors).
#              With the flag BUSCO does neither and only downloads its datasets (checked in
#              run_BUSCO.py of BUSCO 6.0.0; the log then lacks the closing line "Anonymous usage
#              data is gathered to improve the tool"). The results are the same either way.
#              After BUSCO has finished, the download folder busco_tmp is removed (auto mode
#              only; with a named lineage there is none), and the result folder is tidied in
#              two steps:
#                1. Two bulky sub-folders are removed wherever they sit in the result
#                   tree: hmmer_output (one search-result file per gene) and
#                   busco_sequences (the FASTA of every gene found). find looks for them
#                   (-prune: do not look inside a folder that is about to be deleted) and
#                   hands them to rm -rf ({{}} becomes {}, where find puts the folder
#                   names). The rest of BUSCO's own result tree stays: run_<lineage>/ with
#                   full_table.tsv (one line per expected gene: complete, duplicated,
#                   fragmented or missing, and where) and missing_busco_list.tsv,
#                   auto_lineage/, logs/ and tmp/.
#                2. BUSCO writes everything into a folder named after the sample, inside
#                   the output folder it was given. The content of that folder is moved up
#                   into BUSCO_DIR itself and the folder is removed, so that the summaries
#                   sit directly in 02.assembly/{sample}/eval/busco/.
# Produces:    BUSCO_DIR (02.assembly/{sample}/eval/busco/) with
#              short_summary.generic.<lineage>.<sample>.txt / .json and
#              short_summary.specific.<lineage>.<sample>.txt / .json (the lineage name is
#              part of the file name, which is how the two are told apart), plus the kept
#              result tree described above.
# Consumed by: MultiQC (rule multiqc in shared/90_report.smk), which reads the .txt
#              summaries; the kept tree is for the people who want to see which genes are
#              missing.
#
# Unlike QUAST, a BUSCO failure stops the run: this command has no failure guard, and
# only QUAST, Qualimap and NanoPlot may fail without stopping it
# (docs/getting-started/from-v1.md).
#
# v1 -> v2: (1) the lineage flag comes from BUSCO_LINEAGE_ARGS instead of being fixed to
# --auto-lineage-euk, so a lineage can be chosen in the configuration. (2) CHANGED: v1 deleted
# BUSCO's run_* folders and auto_lineage after the run, keeping only the summaries and
# BUSCO's logs/ and tmp/ folders; v2 keeps the result tree and removes just hmmer_output and
# busco_sequences. MultiQC therefore has to be told to ignore the kept folders, or it would
# report every sample several times (see the --ignore lines in shared/90_report.smk).
# (3) Threads are Snakemake's own, as for QUAST. (4) The tidying steps now write their
# messages to the log as well (the log name, completeness_{sample}.log, is the same as in
# v1). (5) --opt-out-run-stats is new (see above). (6) v1's --force is gone. BUSCO reads it in
# one place only: when its result folder for the sample (BUSCO_DIR/{sample}) already exists,
# --force deletes it, and without --force BUSCO refuses to start. That folder never exists when
# the job starts, because it lies inside BUSCO_DIR, a directory() output, which Snakemake
# deletes with everything in it before every job, so the flag never did anything here.
#
# (7) NEW: a named lineage is downloaded once per run by rule busco_data and read with
# --offline (see above). v1 had every sample download its data into busco_tmp, declared as
# a second, temporary output; auto mode still does, and the shell now removes that folder
# itself, because with a named lineage there is none.
#
# The one-lineage case: with a lineage given, BUSCO writes no generic summary, so MultiQC
# shows ONE BUSCO row per sample instead of two.
rule completeness:
    input:
        contigs = FINAL_CONTIGS,
        busco_data = BUSCO_DATA_INPUT,
    output:
        busco_dir = directory(BUSCO_DIR),
    params:
        lineage_args = BUSCO_LINEAGE_ARGS,
        download_args = BUSCO_DOWNLOAD_ARGS,
    threads: capped_cpus(24)
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/busco.yaml"
    log:
        LOGS + "/completeness_{sample}.log"
    shell:
        """
        # BUSCO writes its results into a folder named after the sample, inside BUSCO_DIR.
        # The download options (params.download_args) say where the lineage data are (see the
        # comment above).
        busco \
            --mode genome \
            -i {input.contigs} \
            --out {wildcards.sample} \
            --out_path {output.busco_dir} \
            {params.download_args} \
            {params.lineage_args} \
            --opt-out-run-stats \
            --cpu {threads} > {log} 2>&1

        # Remove the download folder of auto mode (with a named lineage there is none).
        rm -rf {output.busco_dir}/busco_tmp

        # Tidy step 1: remove the two bulky sub-folders, wherever they sit in the result tree.
        find {output.busco_dir}/{wildcards.sample} -type d -name hmmer_output -prune -exec rm -rf {{}} + >> {log} 2>&1
        find {output.busco_dir}/{wildcards.sample} -type d -name busco_sequences -prune -exec rm -rf {{}} + >> {log} 2>&1

        # Tidy step 2: move the results up into BUSCO_DIR and remove the sample folder.
        mv {output.busco_dir}/{wildcards.sample}/* {output.busco_dir} >> {log} 2>&1
        rm -rf {output.busco_dir}/{wildcards.sample}
        """


# ───────────── map_evaluation: Qualimap, how the reads map back to the draft ─────────────
# Defined only in the read-based modes: the `if` below is decided once, when Snakemake reads
# the workflow, before any job runs. Contigs mode has no reads, so there is nothing to map.
#
# Takes in:    DECONTAM_BAM and its index, the sample's reads aligned to the draft assembly
#              by rule map_contigs (shared/10_decontam.smk): the trimmed Illumina pairs in
#              illumina mode, the filtered long reads in longread and hybrid mode. The same
#              BAM gives BlobTools its coverage track, so the reads are mapped once. (In
#              hybrid mode the Illumina reads mapped to the guide assembly get a Qualimap
#              report of their own, from rule map_guide_evaluation in rules/hybrid/25_guide.smk.)
# Does:        qualimap bamqc: reads the whole BAM and writes the mapping rate, the depth
#              per contig and its histogram, the insert-size distribution (paired reads) and
#              the GC content of the mapped reads, as a page (-outformat html) with the
#              numbers in genome_results.txt and raw_data_qualimapReport/ beside it.
# Produces:    QUALIMAP_DIR (02.assembly/{sample}/eval/qualimap/), declared as a directory()
#              output, so Snakemake treats the whole folder as the result.
# Consumed by: MultiQC (rule multiqc in shared/90_report.smk), which reads
#              genome_results.txt and the histograms.
#
# Memory. Qualimap is a Java program, and Java uses only as much memory as its heap allows.
# --java-mem-size sets the heap: resources.ram_gb (RAM), capped at 64 GB, which is plenty
# for a fungal BAM. The value is a Snakemake `resources:` entry (java_mem), not a param: a
# memory cap changes nothing in the result, so a changed ram_gb must not re-run the report.
# `snakemake --resources java_mem=N` then limits how many reports run at once. N must be at
# least the heap value: Snakemake scales a job's resource down to the global limit, so a
# smaller N also lowers the heap Qualimap is given. On a cluster the job asks the scheduler
# for the heap plus 4 GB, but never for more than ram_gb (the mem_mb line in the rule says
# why). Threads: -nt, capped at 24.
#
# Headless Java. Qualimap draws its plots with Java's graphics library, which by default
# looks for a screen (the DISPLAY variable) and fails on a server without one. Two things
# make it draw without a screen: `unset DISPLAY` (so no half-set value from a login session
# can mislead it) and JAVA_TOOL_OPTIONS="-Djava.awt.headless=true", a variable the Java
# runtime itself reads and ADDS to its options. It must not be JAVA_OPTS: Qualimap's start
# script takes JAVA_OPTS, when set, as the complete list of Java options and then drops its
# own --java-mem-size, so the heap would silently fall back to Java's default. The Java
# runtime prints one line, "Picked up JAVA_TOOL_OPTIONS: ...", to the error stream; it lands
# in the log and is expected.
#
# Qualimap is only a report, so like QUAST it may fail without stopping the run: the output
# folder is created first, the part after || runs only when Qualimap fails and writes a
# [FunFlux] line to the log, and MultiQC then simply has no mapping table for the sample.
#
# v1 -> v2: the Qualimap command is v1 FunFlux's (FunFluxL and FunFluxL+ had the same one).
# (1) The CPU count is Snakemake's own threads: (capped at 24) and the heap a Snakemake
# resource, where v1 used its private resources: for both. (2) The headless Java settings
# are new: the v1 runs relied on the server's environment allowing a display. (3) The
# failure guard is new (docs/getting-started/from-v1.md names Qualimap among the QC rules
# that no longer stop the run). (4) The BAM comes from the sample's contaminants folder and
# the report goes to eval/qualimap/ (v1: 03.post-processing/mapping_evaluation/{sample}/).
# (5) The rule is shared by the three read-based modes; FunFluxL called it map_qc.
if HAS_READS:

    rule map_evaluation:
        input:
            bam = DECONTAM_BAM,
            bai = DECONTAM_BAM + ".bai",
        output:
            qualimap_dir = directory(QUALIMAP_DIR),
        threads: capped_cpus(24)
        resources:
            java_mem = min(RAM, 64),
            # Asked of a cluster scheduler: the Java heap plus 4 GB for Java itself, in MB, but
            # never more than resources.ram_gb (HEAVY_MEM_MB), and a time limit (00_common.smk,
            # section 4). The cap matters when ram_gb is below 68: heap plus 4 GB would then
            # be more than the ram_gb the user set to fit one compute node, and the scheduler
            # could never place the job. The heap is only Java's upper limit, which it grows
            # into as needed, not memory it takes at the start.
            mem_mb = min((min(RAM, 64) + 4) * 1024, HEAVY_MEM_MB),
            runtime = LONG_RUNTIME_MINUTES,
        conda:
            "../../envs/qualimap.yaml"
        log:
            LOGS + "/map_evaluation_{sample}.log"
        shell:
            """
            mkdir -p {output.qualimap_dir}

            unset DISPLAY
            JAVA_TOOL_OPTIONS="-Djava.awt.headless=true" \
            qualimap bamqc \
                -bam {input.bam} \
                --java-mem-size={resources.java_mem}G \
                -nt {threads} \
                -outdir {output.qualimap_dir} \
                -outformat html > {log} 2>&1 || {{
                echo "" >> {log}
                echo "[FunFlux] map_evaluation: Qualimap exited with an error for sample {wildcards.sample}. Qualimap is only a report that no other step reads, so the run continues; the mapping table will be missing from the MultiQC page for this sample. See the Qualimap messages above." >> {log}
            }}
            """
