# FunFlux v2.0.0 — hybrid front end, step 3: the guide screen and the read selection
# (rules/hybrid/25_guide.smk).
#
# What this file is for. In hybrid mode the Illumina reads are not assembled into the
# delivered genome; that is the long reads' job. The Illumina reads do two other things:
# they score the ONT reads in Filtlong (rules/hybrid/30_long_reads.smk) and they correct
# the finished long-read assembly with Polypolish (rules/hybrid/50_polish.smk). Both uses
# want ONLY the reads of the isolate. A read pair from a contaminant (a bacterium or a
# yeast that grew in the same flask) would score the contaminant's ONT reads as good reads
# and would polish the assembly with the wrong sequence. The GUIDE assembly, SPAdes on
# the Illumina pairs (rules/hybrid/20_assembly.smk), exists to tell the two apart: it is
# screened like any assembly (coverage track, BLAST, BlobTools, mitochondrion screen,
# taxonomy selector), and the Illumina pairs that align properly onto the contigs the
# screen KEPT are the "selected" pairs the two later steps use. The guide itself is not
# delivered: after this file only QUAST reads it (rule genome_assembly_evaluation in
# shared/20_qc.smk), to print its statistics next to the final genome's under the label
# {sample}_guide, so that the report shows what the long reads added.
#
# The Snakefile includes every file in rules/hybrid/ in hybrid mode only, in name order;
# this file comes after 20_assembly.smk and before 30_long_reads.smk.
#
# The ten rules, in the order the data flows. {sample} is a wildcard: a placeholder that
# Snakemake fills in with each sample's name, so every rule runs once per sample.
# Everything lies under 02.assembly/{sample}/guide/ in the output folder, except the
# selected pairs, which go to 01.reads/{sample}/selected/.
#
#   index_guide_contigs        bowtie2 index of the filtered guide (GUIDE_DRAFT)
#   map_guide_contigs          the trimmed Illumina pairs aligned to it: the coverage track
#                              (GUIDE_BAM, kept)
#   map_guide_evaluation       Qualimap on that BAM: the "mapping Illumina QC" row of the
#                              MultiQC page (GUIDE_QUALIMAP_DIR)
#   blast_guide_contigs        megablast of every guide contig against NCBI nt (GUIDE_BLASTOUT)
#   blob_guide_json            BlobTools joins the contigs, the coverage and the BLAST hits
#   blob_guide_table           ... and writes the per-contig taxonomy table
#   screen_guide_mitochondria  the mitochondrion screen of the guide, from the BLAST hits alone
#                              (GUIDE_MITO_SCREEN_TSV, GUIDE_MITO_IDS)
#   select_guide_contigs       the taxonomy selector: the screened guide (GUIDE_DECONTAM) plus
#                              the audit files, with the confirmed mitochondrial contigs
#                              kept whatever the genus filter says
#   index_selected_contigs     bowtie2 index of the screened guide
#   map_sel_contigs            the trimmed pairs aligned to the screened guide; the properly
#                              paired ones are written out as the selected pairs (SEL_R1,
#                              SEL_R2, kept)
#
#   GUIDE_DRAFT  (guide/contigs_filt.fasta, rule filter_contigs in rules/hybrid/20_assembly.smk)
#     |
#     +--> index_guide_contigs --> map_guide_contigs (+ TRIM_R1/R2) --> GUIDE_BAM
#     |                                 GUIDE_BAM --> map_guide_evaluation --> GUIDE_QUALIMAP_DIR --> MultiQC
#     |                                 GUIDE_BAM --> blob_guide_json (below)
#     +--> blast_guide_contigs --> GUIDE_BLASTOUT
#     |                                 GUIDE_BLASTOUT --> blob_guide_json --> blob_guide_table --> taxonomy table
#     |                                 GUIDE_BLASTOUT --> screen_guide_mitochondria --> GUIDE_MITO_IDS
#     +--> select_guide_contigs (taxonomy table + GUIDE_MITO_IDS + GUIDE_DRAFT)
#            --> GUIDE_DECONTAM  (guide/contaminants/contigs_decontam.fasta)
#                  --> index_selected_contigs --> map_sel_contigs (+ TRIM_R1/R2)
#                        --> SEL_R1, SEL_R2  (01.reads/{sample}/selected/{sample}_selected_R1.fastq, _R2)
#                              --> rule filter_long_reads (rules/hybrid/30_long_reads.smk), ONT samples
#                              --> rule short_read_correction (rules/hybrid/50_polish.smk), every sample
#                  --> rule genome_assembly_evaluation (shared/20_qc.smk), next to the final genome
#
# The same screen as the long-read assembly's, as a second set of rules. The first eight
# rules below are copies of the shared screen rules with the guide's files in place of the
# draft's: index_contigs and map_contigs (the illumina-mode versions), blast_contigs,
# blob_json, blob_table and select_contigs in shared/10_decontam.smk; screen_mitochondria
# in shared/15_mitochondria.smk; map_evaluation in shared/20_qc.smk. They are copies, not
# the shared rules themselves, because a hybrid run screens TWO assemblies of every sample,
# and Snakemake needs two rules with two sets of output paths for that: the long-read
# assembly (assembly.fasta) is screened by the shared rules, under 02.assembly/{sample}/
# contaminants/ and mitochondria/; the guide by these, under 02.assembly/{sample}/guide/.
# The two screens are independent and apply ONE policy: parameters.decontamination
# (DECONTAMINATION) and parameters.mitochondria (MITO_MODE and the cut-offs) for both. The
# guide screen decides which Illumina PAIRS are used; the long-read screen decides which
# CONTIGS are delivered. The long explanations of each step (why megablast, what the
# "bestsum" rule is, how a contig is called mitochondrial, how the selector decides and
# what its audit files hold) stand in the shared files and are not repeated here; each
# rule below says what is the same and what differs. Keep in sync: a change to a shared
# screen command must be made here too.
#
# With parameters.decontamination.mode off (the default, and the only setting validated for
# fungi) the selector keeps every contig: the screened guide is then the filtered guide
# re-wrapped at 80 characters per line, the selected pairs are the pairs properly paired on
# any filtered guide contig, and the mitochondrion protection changes nothing. The whole
# screen still runs, for the audit files and the Qualimap report.
#
# Everything referenced here is defined once in 00_common.smk: GUIDE_DRAFT,
# GUIDE_DECONTAM_DIR, GUIDE_BAM, GUIDE_BLASTOUT, GUIDE_DECONTAM, GUIDE_MITO_SCREEN_TSV,
# GUIDE_MITO_IDS, GUIDE_QUALIMAP_DIR, TRIM_R1, TRIM_R2, SEL_R1, SEL_R2, BLASTDB, NT_VERSION,
# BLAST_MAX_TARGET_SEQS, BLAST_MAX_HSPS, DECONTAMINATION, SELECT_TAXONOMY_SCRIPT,
# MITOCHONDRIA_SCRIPT, MITO_MODE, MITO_MIN_QCOV, MITO_MIN_RATIO, MITO_VERY_STRONG_QCOV, LOGS,
# CPUS, RAM, HEAVY_MEM_MB and LONG_RUNTIME_MINUTES, and the function capped_cpus(). The
# working files that only the rules of this file use are named just below, with an
# underscore in front, as the other rule files do with their own.
#
# conda: environment paths are written relative to THIS file (workflow/rules/hybrid/), so
# "../../envs/bowtie.yaml" goes up two folders, to workflow/, and then into workflow/envs/.
# The tools get their own environments (bowtie2 and samtools, Qualimap, BLAST+, BlobTools).
# The mitochondrion screen and the selector are Python scripts that use only the standard
# library; their two rules run them in envs/python.yaml, which holds a pinned Python 3.12 and
# nothing else, as the shared rules do. Activating that environment makes `python` that
# Python, instead of whatever `python` the computer has.
#
# CPUs: a rule that uses several CPUs declares Snakemake's own `threads:` and passes
# {threads} to its tool. That is what makes `--cores N` actually limit the run. Note that
# the guide BLAST is a second BLAST run of the same size as the long-read assembly's, each
# of up to 24 threads: with --cores 24 the two run one after the other.
#
# In the shell commands below, a doubled brace ({{ or }}) stands for a single brace in the
# shell. Snakemake reads single braces as its own placeholders, such as {input.bam}, and
# turns a doubled brace into a single one.


# ── The guide screen's working files: only this file uses them ──
# All in the guide's contaminants folder (GUIDE_DECONTAM_DIR, 02.assembly/{sample}/guide/
# contaminants/), the twin of the long-read screen's folder. The names follow the shared
# screen's (00_common.smk, section 5), so that the two folders read alike.
#
# The two bowtie2 index prefixes. bowtie2 adds the six .bt2 endings itself. "filtered" is
# the index of the filtered guide (GUIDE_DRAFT) for the coverage track, "selected" the
# index of the screened guide (GUIDE_DECONTAM) for the read selection; the two names are
# FunFluxL+'s, and they must differ because both indexes lie in the same folder.
_GUIDE_DRAFT_INDEX_PREFIX    = GUIDE_DECONTAM_DIR + "/filtered"
_GUIDE_SELECTED_INDEX_PREFIX = GUIDE_DECONTAM_DIR + "/selected"

# BlobTools' files. BLOB_PREFIX is handed to blobtools create, which adds the endings
# itself: the database (JSON) and the coverage file. CAUTION: BlobTools names the coverage
# file after the BAM's file name, so _GUIDE_BLOB_COV is DERIVED from GUIDE_BAM instead of
# being typed out again (exactly as 00_common.smk derives BLOB_COV from DECONTAM_BAM).
# Renaming the BAM would otherwise leave this rule waiting for a file BlobTools never
# writes. The table prefix is handed to blobtools view, which adds .blob.blobDB.table.txt.
_GUIDE_BLOB_PREFIX       = GUIDE_DECONTAM_DIR + "/blob"
_GUIDE_BLOB_JSON         = _GUIDE_BLOB_PREFIX + ".blobDB.json"
_GUIDE_BLOB_COV          = _GUIDE_BLOB_PREFIX + "." + os.path.basename(GUIDE_BAM) + ".cov"
_GUIDE_BLOB_TABLE_PREFIX = GUIDE_DECONTAM_DIR + "/bestscore"
_GUIDE_BLOB_TABLE        = _GUIDE_BLOB_TABLE_PREFIX + ".blob.blobDB.table.txt"

# The selector's three audit files (the same three the long-read screen writes: the genus
# shares, the IDs of the kept contigs, and the decision per contig with its reason).
_GUIDE_COMPOSITION      = GUIDE_DECONTAM_DIR + "/{sample}_composition.txt"
_GUIDE_CONTIG_LIST      = GUIDE_DECONTAM_DIR + "/contigs.list"
_GUIDE_CONTIG_DECISIONS = GUIDE_DECONTAM_DIR + "/contig_taxonomy_decisions.tsv"

# The alignments of the read selection (rule map_sel_contigs), SAM text, temporary.
_GUIDE_SELECTION_SAM = GUIDE_DECONTAM_DIR + "/selected_map.sam"


# ───────────── index_guide_contigs — the bowtie2 index of the filtered guide ─────────────
# The copy of rule index_contigs (illumina mode, shared/10_decontam.smk) for the guide;
# that rule's comment explains the six index files and why one thread is enough.
#
# Takes in:    GUIDE_DRAFT, the SPAdes guide after the length and coverage filter (rule
#              filter_contigs in rules/hybrid/20_assembly.smk): guide/contigs_filt.fasta.
# Does:        bowtie2-build -f writes the index under _GUIDE_DRAFT_INDEX_PREFIX, in the
#              guide's contaminants folder.
# Produces:    the six .bt2 files, temporary: only map_guide_contigs reads them.
# Consumed by: map_guide_contigs below.
#
# v1 -> v2: the command is unchanged (FunFluxL+ called the rule index_short_contigs). The
# six files moved from 03.read_selection/contamination/{sample}/ to the guide's
# contaminants folder.
rule index_guide_contigs:
    input:
        contigs = GUIDE_DRAFT,
    output:
        idx1 = temp(_GUIDE_DRAFT_INDEX_PREFIX + ".1.bt2"),
        idx2 = temp(_GUIDE_DRAFT_INDEX_PREFIX + ".2.bt2"),
        idx3 = temp(_GUIDE_DRAFT_INDEX_PREFIX + ".3.bt2"),
        idx4 = temp(_GUIDE_DRAFT_INDEX_PREFIX + ".4.bt2"),
        ridx1 = temp(_GUIDE_DRAFT_INDEX_PREFIX + ".rev.1.bt2"),
        ridx2 = temp(_GUIDE_DRAFT_INDEX_PREFIX + ".rev.2.bt2"),
    params:
        basename = _GUIDE_DRAFT_INDEX_PREFIX,
    conda:
        "../../envs/bowtie.yaml"
    log:
        LOGS + "/index_guide_contigs_{sample}.log"
    shell:
        """
        bowtie2-build -f {input.contigs} {params.basename} > {log} 2>&1
        """


# ───────────── map_guide_contigs — the trimmed pairs aligned to the filtered guide ─────────────
# The copy of rule map_contigs (illumina mode, shared/10_decontam.smk) for the guide: the
# coverage track of the guide screen, and the BAM Qualimap evaluates. It says how many of
# the sample's own Illumina pairs land on each guide contig. The contigs of the isolate
# share one read depth, while a contaminant that grew to another cell density sits at a
# depth of its own, which is what the BlobTools plot (GC against coverage) separates.
#
# Takes in:    the six index files from index_guide_contigs, and the trimmed, PhiX-free
#              read pairs TRIM_R1 and TRIM_R2 (rule trim_adapters in
#              rules/hybrid/10_reads.smk). bowtie2 is given the index PREFIX (params db),
#              not the files; the files are listed as inputs so that Snakemake builds the
#              index first.
# Does:        bowtie2 aligns the pairs in its default end-to-end mode (-p: the threads)
#              and writes the alignments as SAM text; bowtie2's own summary (how many
#              pairs there were and what share aligned) goes to the log. samtools view -b
#              -u turns the SAM into an uncompressed BAM for the pipe, samtools sort
#              orders the BAM by position, and samtools index writes the .bai index next
#              to it.
# Produces:    GUIDE_BAM (guide/contaminants/{sample}_map.bam) and its index, GUIDE_BAM +
#              ".bai". Both are KEPT, not temporary: Qualimap and BlobTools read the same
#              BAM, and when a later step is re-run (a changed BLAST depth re-runs
#              blob_guide_json, for instance) a temporary BAM would have to be rebuilt by
#              aligning tens of millions of reads again.
# Consumed by: blob_guide_json (the coverage) and map_guide_evaluation (the mapping
#              statistics).
#
# Threads: all CPUS, for bowtie2 (-p) and the three samtools commands (-@). The whole
# Illumina library is aligned; with map_sel_contigs below it is one of the two longest
# bowtie2 steps of the mode.
#
# v1 -> v2: bowtie2's command is FunFluxL+'s (rule map_short_contigs), unchanged; it has no
# -t (bowtie2's timings), which the illumina-mode map_contigs has from v1 FunFlux. The
# samtools commands are written as in the shared rule (samtools view -b -u for FunFluxL+'s
# -hbS; every command writes its messages to the log); the alignments are the same. The BAM
# moved from 03.read_selection/contamination/{sample}/filtered_map.bam to the guide's
# contaminants folder and carries the long-read screen's BAM name, {sample}_map.bam. The CPU
# count is Snakemake's own threads:, which --cores enforces, where FunFluxL+ used a private
# "resources: cpus".
rule map_guide_contigs:
    input:
        idx1 = _GUIDE_DRAFT_INDEX_PREFIX + ".1.bt2",
        idx2 = _GUIDE_DRAFT_INDEX_PREFIX + ".2.bt2",
        idx3 = _GUIDE_DRAFT_INDEX_PREFIX + ".3.bt2",
        idx4 = _GUIDE_DRAFT_INDEX_PREFIX + ".4.bt2",
        ridx1 = _GUIDE_DRAFT_INDEX_PREFIX + ".rev.1.bt2",
        ridx2 = _GUIDE_DRAFT_INDEX_PREFIX + ".rev.2.bt2",
        r1 = TRIM_R1,
        r2 = TRIM_R2,
    output:
        bam = GUIDE_BAM,
        bai = GUIDE_BAM + ".bai",
    params:
        db = _GUIDE_DRAFT_INDEX_PREFIX,
    threads: CPUS
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/bowtie.yaml"
    log:
        LOGS + "/map_guide_contigs_{sample}.log"
    shell:
        """
        bowtie2 \
          -x {params.db} \
          -1 {input.r1} \
          -2 {input.r2} \
          -p {threads} 2> {log} | \
        samtools view \
          -b \
          -u \
          -@ {threads} 2>> {log} | \
        samtools sort \
          -o {output.bam} \
          -@ {threads} 2>> {log}

        samtools index \
          {output.bam} \
          -@ {threads} 2>> {log}
        """


# ───────────── map_guide_evaluation — Qualimap, how the Illumina reads map to the guide ─────────────
# The copy of rule map_evaluation (shared/20_qc.smk) for the guide. The two Qualimap
# reports of a hybrid sample answer two questions: this one, how the Illumina reads fit
# the assembly made from them (mapping rate, insert sizes, depth per contig and how evenly
# it is spread); the shared one, how the long reads fit the long-read assembly. On the
# MultiQC page this report is the row "mapping Illumina QC", the shared one "mapping QC".
#
# Takes in:    GUIDE_BAM and its index (rule map_guide_contigs).
# Does:        qualimap bamqc: reads the whole BAM and writes the mapping rate, the depth
#              per contig and its histogram, the insert-size distribution and the GC
#              content of the mapped reads, as a page (-outformat html) with the numbers in
#              genome_results.txt and raw_data_qualimapReport/ beside it.
# Produces:    GUIDE_QUALIMAP_DIR (02.assembly/{sample}/guide/eval/qualimap/), declared as
#              a directory() output, so Snakemake treats the whole folder as the result.
#              Rule all asks for it.
# Consumed by: MultiQC (rule multiqc in shared/90_report.smk).
#
# Memory, Java and the failure guard: as in the shared rule, whose comment explains each.
# --java-mem-size is the Java heap (resources.ram_gb capped at 64 GB, a `resources:` entry
# so that a changed cap does not re-run the report); `unset DISPLAY` and JAVA_TOOL_OPTIONS
# let Qualimap draw its plots without a screen; a failing Qualimap writes a [FunFlux] line
# to the log and the run continues, with no "mapping Illumina QC" row for the sample.
#
# Threads: up to 24 (capped_cpus), as in the shared rule.
#
# v1 -> v2: the Qualimap command is FunFluxL+'s (rule short_mapping_qc). As for the
# shared rule: Snakemake's own threads: and a resource for the heap instead of private
# resources, the headless Java settings and the failure guard are new, and the report
# moved from 03.read_selection/mapping_evaluation/{sample}/ to the guide's eval folder.
rule map_guide_evaluation:
    input:
        bam = GUIDE_BAM,
        bai = GUIDE_BAM + ".bai",
    output:
        qualimap_dir = directory(GUIDE_QUALIMAP_DIR),
    threads: capped_cpus(24)
    resources:
        java_mem = min(RAM, 64),
        # Asked of a cluster scheduler: the Java heap plus 4 GB for Java itself, in MB, but
        # never more than resources.ram_gb (HEAVY_MEM_MB), and a time limit (00_common.smk,
        # section 4; rule map_evaluation explains the cap).
        mem_mb = min((min(RAM, 64) + 4) * 1024, HEAVY_MEM_MB),
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/qualimap.yaml"
    log:
        LOGS + "/map_guide_evaluation_{sample}.log"
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
            echo "[FunFlux] map_guide_evaluation: Qualimap exited with an error for sample {wildcards.sample}. Qualimap is only a report that no other step reads, so the run continues; the Illumina mapping table of the guide assembly will be missing from the MultiQC page for this sample. See the Qualimap messages above." >> {log}
        }}
        """


# ───────────── blast_guide_contigs — the taxonomy evidence for the guide ─────────────
# The copy of rule blast_contigs (shared/10_decontam.smk) for the guide. BLAST asks, for
# every guide contig, which sequences in the NCBI nucleotide database look like it and
# which organism each came from; BlobTools sums those matches into a taxonomic call per
# contig, and the mitochondrion screen reads the same matches to find a contig that is a
# mitochondrial genome.
#
# Takes in:    GUIDE_DRAFT, the filtered guide (rule filter_contigs in
#              rules/hybrid/20_assembly.smk). The screen works on the filtered guide on
#              purpose: it has to see every contig it is asked to judge, including the ones
#              the selector is about to drop. The database folder comes from
#              directories.blast_db (BLASTDB) and the database name from
#              parameters.nt_version (NT_VERSION).
# Does:        one `blastn -task megablast` run per sample, with the same database, the same
#              hit depth (parameters.blast_max_target_seqs, parameters.blast_max_hsps; the
#              shared rule explains why the two numbers matter, and
#              docs/methods/blast-hit-depth.md the choice) and the same 15-column -outfmt
#              as the shared rule. Do NOT change the column list: BlobTools and the
#              mitochondrion screen pick their columns by POSITION.
# Produces:    GUIDE_BLASTOUT (guide/contaminants/blastout), the tabular hits, kept.
# Consumed by: blob_guide_json and screen_guide_mitochondria below.
#
# Threads: up to 24 (capped_cpus). This is the second BLAST run of the hybrid mode, as
# large as the one on the long-read assembly; see the CPU note at the head of this file.
#
# v1 -> v2: FunFluxL+ (rule blast_short_contigs) ran with -max_target_seqs 50 -max_hsps 5
# and the database name typed into the rule; both depths and the name now come from the
# configuration, with the same defaults, and the CPU count is Snakemake's own threads:.
# The hits moved from 03.read_selection/contamination/{sample}/blastout to the guide's
# contaminants folder.
rule blast_guide_contigs:
    input:
        contigs = GUIDE_DRAFT,
    output:
        blast = GUIDE_BLASTOUT,
    params:
        dir = BLASTDB,
        db = os.path.join(BLASTDB, NT_VERSION),
        max_target_seqs = BLAST_MAX_TARGET_SEQS,
        max_hsps = BLAST_MAX_HSPS,
    threads: capped_cpus(24)
    # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
    resources:
        mem_mb = HEAVY_MEM_MB,
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/blast.yaml"
    log:
        LOGS + "/blast_guide_contigs_{sample}.log"
    shell:
        """
        BLASTDB={params.dir} \
        blastn \
          -task megablast \
          -query {input.contigs} \
          -db {params.db} \
          -outfmt '6 qseqid staxids bitscore pident evalue length qlen slen qcovs qcovhsp sskingdoms scomnames sscinames sblastnames stitle' \
          -num_threads {threads} \
          -evalue 1e-5 \
          -max_target_seqs {params.max_target_seqs} \
          -max_hsps {params.max_hsps} \
          -out {output.blast} > {log} 2>&1
        """


# ───────────── blob_guide_json — BlobTools: one database from guide + coverage + BLAST ─────────────
# The copy of rule blob_json (shared/10_decontam.smk) for the guide. blobtools create
# writes one database (a JSON file) that holds, for every contig, its length and GC, its
# Illumina read coverage from the BAM, and a taxonomy worked out from the BLAST hits.
#
# Takes in:    GUIDE_DRAFT, the same contigs that were aligned and BLASTed; GUIDE_BAM and
#              its index (map_guide_contigs), the coverage track; GUIDE_BLASTOUT
#              (blast_guide_contigs), the taxonomy evidence; and nodes.dmp and names.dmp
#              from the BLAST database folder, the NCBI taxonomy dump BlobTools needs to
#              turn a taxonomy ID into a name. The two dump files are listed as INPUTS so
#              that a missing file is reported when the run starts, not an hour into it.
# Does:        blobtools create, which writes its results under the prefix
#              _GUIDE_BLOB_PREFIX and adds the file endings itself.
# Produces:    _GUIDE_BLOB_JSON (the database) and _GUIDE_BLOB_COV (the coverage BlobTools
#              read from the BAM, named after the BAM's file name, see above). Both temp():
#              Snakemake deletes the JSON once blob_guide_table has read it, and the
#              coverage file, which no rule reads, as soon as this job has finished.
# Consumed by: blob_guide_table (the JSON only).
#
# v1 -> v2: the command is FunFluxL+'s (rule short_blob_json), unchanged. The files moved
# from 03.read_selection/contamination/{sample}/ to the guide's contaminants folder, and
# the coverage file's name follows the renamed BAM.
rule blob_guide_json:
    input:
        contigs = GUIDE_DRAFT,
        bam = GUIDE_BAM,
        bai = GUIDE_BAM + ".bai",
        blast = GUIDE_BLASTOUT,
        nodes = os.path.join(BLASTDB, "nodes.dmp"),
        names = os.path.join(BLASTDB, "names.dmp"),
    output:
        json = temp(_GUIDE_BLOB_JSON),
        cov = temp(_GUIDE_BLOB_COV),
    params:
        basename = _GUIDE_BLOB_PREFIX,
    conda:
        "../../envs/blobtools.yaml"
    log:
        LOGS + "/blob_guide_json_{sample}.log"
    shell:
        """
        blobtools create \
          -i {input.contigs} \
          -b {input.bam} \
          -t {input.blast} \
          --nodes {input.nodes} \
          --names {input.names} \
          -o {params.basename} > {log} 2>&1
        """


# ───────────── blob_guide_table — one taxonomic call per guide contig and rank ─────────────
# The copy of rule blob_table (shared/10_decontam.smk) for the guide. The "bestsum" rule
# adds up the bit scores of all hits per taxon, separately at each rank, and the taxon with
# the highest total wins that rank; the result is one row per contig that names its
# winner from superkingdom down to species, with the length, GC and coverage beside it.
#
# Takes in:    _GUIDE_BLOB_JSON, the database made by blob_guide_json.
# Does:        blobtools view with --taxrule bestsum, --rank all (every rank) and --hits
#              (also list which taxa contributed which score). BlobTools builds the file
#              name itself from the prefix _GUIDE_BLOB_TABLE_PREFIX.
# Produces:    _GUIDE_BLOB_TABLE (guide/contaminants/bestscore.blob.blobDB.table.txt),
#              kept: it is the evidence behind every keep-or-drop decision of the guide.
# Consumed by: select_guide_contigs, which reads the contig name and the genus column.
#
# v1 -> v2: the command is FunFluxL+'s (rule short_blob_table), unchanged. The table
# moved from 03.read_selection/contamination/{sample}/ to the guide's contaminants folder.
rule blob_guide_table:
    input:
        json = _GUIDE_BLOB_JSON,
    output:
        bestscore = _GUIDE_BLOB_TABLE,
    params:
        basename = _GUIDE_BLOB_TABLE_PREFIX,
    conda:
        "../../envs/blobtools.yaml"
    log:
        LOGS + "/blob_guide_table_{sample}.log"
    shell:
        """
        blobtools view \
          --input {input.json} \
          --out {params.basename} \
          --taxrule bestsum \
          --rank all \
          --hits > {log} 2>&1
        """


# ───────────── screen_guide_mitochondria — which guide contigs are mitochondrial ─────────────
# The copy of rule screen_mitochondria (shared/15_mitochondria.smk) for the guide, from the
# BLAST hits alone. Why the guide is screened for the mitochondrion at all: the guide's
# selector (next rule) applies the same genus filter as the long-read screen, and BLAST
# often gives a fungal mitochondrion a genus that is not the one most of the assembly
# gets. With a genus filter on (parameters.decontamination.mode include or exclude) the
# mitochondrial guide contig could be dropped from the screened guide; its Illumina pairs
# would then not be selected, would not score the ONT reads of the mitochondrion in
# Filtlong, and would not polish the mitochondrial contig of the final genome. This rule
# finds that contig, and select_guide_contigs keeps it whatever the filter says. With mode
# off (the default) nothing is dropped anyway, and the rule only writes its report.
#
# BLAST only. There is no assembler evidence table for a SPAdes assembly (SPAdes reports a
# k-mer coverage per contig but nothing about circularity, and mitochondria.py has no
# evidence-spades command), so --assembly-evidence is an empty value, written as "" for the
# reason the shared rule gives, and the script decides on --min-direct-blast-qcov alone.
# The other two cut-offs change nothing here; they are passed so that the command has the
# shared rule's shape.
#
# Takes in:    GUIDE_BLASTOUT (blast_guide_contigs) and GUIDE_DRAFT, the contigs BLAST was
#              run on (for the contig lengths, and so that every contig gets a row).
# Does:        python mitochondria.py screen (MITOCHONDRIA_SCRIPT) with the settings of
#              parameters.mitochondria through `params:`, as the shared rule.
# Produces:    GUIDE_MITO_SCREEN_TSV (guide/mitochondria/mitochondrial_screen.tsv), the
#              report, one row per guide contig with its class and the reason; and
#              GUIDE_MITO_IDS, the confirmed guide contigs, one name per line (empty unless
#              parameters.mitochondria.mode is separate_confirmed).
# Consumed by: select_guide_contigs (--always-keep), and nothing else. The guide's
#              mitochondrion is never split off, because the guide is not delivered; the
#              mitochondrial contig of the DELIVERED genome is found by the shared screen
#              of the long-read assembly, with the assembler's evidence.
#
# v1 -> v2: new rule. FunFluxL+ screened the guide for contaminants but not for the
# mitochondrion, and its guide selector had no protected list. With mode off the clean
# guide is the same as before.
rule screen_guide_mitochondria:
    input:
        blast = GUIDE_BLASTOUT,
        contigs = GUIDE_DRAFT,
    output:
        report = GUIDE_MITO_SCREEN_TSV,
        confirmed_ids = GUIDE_MITO_IDS,
    params:
        script = MITOCHONDRIA_SCRIPT,
        mode = MITO_MODE,
        minimum_qcov = MITO_MIN_QCOV,
        minimum_coverage_ratio = MITO_MIN_RATIO,
        very_strong_qcov = MITO_VERY_STRONG_QCOV,
    conda:
        "../../envs/python.yaml"
    log:
        LOGS + "/screen_guide_mitochondria_{sample}.log"
    shell:
        """
        python {params.script:q} screen \
          --blast {input.blast:q} \
          --contigs {input.contigs:q} \
          --assembly-evidence "" \
          --report {output.report:q} \
          --confirmed-ids {output.confirmed_ids:q} \
          --mode {params.mode:q} \
          --min-direct-blast-qcov {params.minimum_qcov:q} \
          --min-assembler-coverage-ratio {params.minimum_coverage_ratio:q} \
          --very-strong-qcov {params.very_strong_qcov:q} > {log} 2>&1
        """


# ───────────── select_guide_contigs — the screened guide and its audit trail ─────────────
# The copy of rule select_contigs (shared/10_decontam.smk) for the guide: where guide
# contigs are kept or dropped. The selector script reads the BlobTools table, takes each
# contig's genus and applies parameters.decontamination.mode (off: keep everything, the
# default; include: keep only the listed genera; exclude: drop only the listed genera),
# plus discard_no_hit for the contigs BLAST could not place. A contig on the protected
# list (GUIDE_MITO_IDS) is kept whatever the mode says. The shared rule's comment explains
# the policy, why there is no mode auto, the per-sample tables and when the script stops
# with an error (no contig kept, a protected ID missing from the FASTA, unusable input).
#
# Takes in:    _GUIDE_BLOB_TABLE (blob_guide_table), GUIDE_DRAFT (the sequences, read as
#              they are) and GUIDE_MITO_IDS (screen_guide_mitochondria).
# Does:        python select_contigs_by_taxonomy.py (SELECT_TAXONOMY_SCRIPT) with the
#              settings of parameters.decontamination (DECONTAMINATION, checked once by
#              00_common.smk and turned into plain text): the SAME policy as the long-read
#              screen, so one sample is screened by one rule on both sides. The :q in
#              {params.include_genera:q} and the others quotes a value for the shell and
#              writes an empty value as nothing, which the script reads as "not set".
# Produces:    _GUIDE_COMPOSITION, _GUIDE_CONTIG_LIST and _GUIDE_CONTIG_DECISIONS, the
#              three audit files (genus shares; the IDs kept; every contig with kept or
#              removed and the reason), kept, for people; and GUIDE_DECONTAM, the clean
#              guide (guide/contaminants/contigs_decontam.fasta), wrapped at 80 characters
#              per line, kept.
# Consumed by: GUIDE_DECONTAM by index_selected_contigs below (the read selection) and by
#              rule genome_assembly_evaluation in shared/20_qc.smk (QUAST, next to the final
#              genome). No rule reads the three audit files.
# Log:         the script's warnings (two or more genera above 5 % of the DNA, more than
#              20 % of the DNA removed, contigs in the table but not in the FASTA).
#
# v1 -> v2: the command is FunFluxL+'s (rule select_short_contigs) plus --always-keep, the
# protected list from the new screen_guide_mitochondria, and a log, which FunFluxL+ had
# not. The screened guide moved from 02.short-read_assembly/{sample}/contigs_selected.fasta to
# the guide's contaminants folder; the audit files from 03.read_selection/contamination/.
# The script runs with the pinned Python of envs/python.yaml, where FunFluxL+ used the
# Python that started Snakemake.
rule select_guide_contigs:
    input:
        bestscore = _GUIDE_BLOB_TABLE,
        contigs = GUIDE_DRAFT,
        confirmed_ids = GUIDE_MITO_IDS,
    output:
        abund = _GUIDE_COMPOSITION,
        list = _GUIDE_CONTIG_LIST,
        decisions = _GUIDE_CONTIG_DECISIONS,
        contigs = GUIDE_DECONTAM,
    params:
        selector = SELECT_TAXONOMY_SCRIPT,
        mode = DECONTAMINATION["mode"],
        include_genera = DECONTAMINATION["include_genera"],
        include_genera_by_sample = DECONTAMINATION["include_genera_by_sample"],
        exclude_genera = DECONTAMINATION["exclude_genera"],
        exclude_genera_file = DECONTAMINATION["exclude_genera_file"],
        sample_overrides = DECONTAMINATION["sample_overrides"],
        discard_no_hit = DECONTAMINATION["discard_no_hit"],
    conda:
        "../../envs/python.yaml"
    log:
        LOGS + "/select_guide_contigs_{sample}.log"
    shell:
        """
        python {params.selector:q} \
          --bestscore {input.bestscore:q} \
          --contigs {input.contigs:q} \
          --always-keep {input.confirmed_ids:q} \
          --sample {wildcards.sample:q} \
          --mode {params.mode:q} \
          --include-genera {params.include_genera:q} \
          --include-genera-by-sample {params.include_genera_by_sample:q} \
          --exclude-genera {params.exclude_genera:q} \
          --exclude-genera-file {params.exclude_genera_file:q} \
          --sample-overrides {params.sample_overrides:q} \
          --discard-no-hit {params.discard_no_hit:q} \
          --output-list {output.list:q} \
          --output-fasta {output.contigs:q} \
          --composition {output.abund:q} \
          --decisions {output.decisions:q} > {log} 2>&1
        """


# ───────────── index_selected_contigs — the bowtie2 index of the screened guide ─────────────
# A second bowtie2 index, of the CLEAN guide this time, for the read selection. The first
# index (index_guide_contigs) was of the filtered guide, before the screen; the pairs must
# now be aligned to the contigs the screen kept, and nothing else, so that a pair whose
# contig was dropped cannot be selected.
#
# Takes in:    GUIDE_DECONTAM (select_guide_contigs).
# Does:        bowtie2-build -f, under _GUIDE_SELECTED_INDEX_PREFIX, in the guide's
#              contaminants folder next to the first index.
# Produces:    the six .bt2 files, temporary: only map_sel_contigs reads them.
# Consumed by: map_sel_contigs below.
#
# Threads: one, as for index_guide_contigs.
#
# v1 -> v2: the command is unchanged (FunFluxL+ called the rule
# index_selected_short_contigs); the files moved to the guide's contaminants folder.
rule index_selected_contigs:
    input:
        contigs = GUIDE_DECONTAM,
    output:
        idx1 = temp(_GUIDE_SELECTED_INDEX_PREFIX + ".1.bt2"),
        idx2 = temp(_GUIDE_SELECTED_INDEX_PREFIX + ".2.bt2"),
        idx3 = temp(_GUIDE_SELECTED_INDEX_PREFIX + ".3.bt2"),
        idx4 = temp(_GUIDE_SELECTED_INDEX_PREFIX + ".4.bt2"),
        ridx1 = temp(_GUIDE_SELECTED_INDEX_PREFIX + ".rev.1.bt2"),
        ridx2 = temp(_GUIDE_SELECTED_INDEX_PREFIX + ".rev.2.bt2"),
    params:
        basename = _GUIDE_SELECTED_INDEX_PREFIX,
    conda:
        "../../envs/bowtie.yaml"
    log:
        LOGS + "/index_selected_contigs_{sample}.log"
    shell:
        """
        bowtie2-build -f {input.contigs} {params.basename} > {log} 2>&1
        """


# ───────────── map_sel_contigs — select the Illumina pairs that belong to the isolate ─────────────
# HOW the Illumina reads are selected and WHAT for. The trimmed pairs are aligned to the
# screened guide, and only the pairs that align PROPERLY (both mates aligned, facing each
# other, at a plausible distance apart) are written out as FASTQ again. The result is the
# clean Illumina read set: the pairs that belong to contigs the screen kept. It is used
# twice: as Filtlong's yardstick for the ONT reads (rule filter_long_reads in
# rules/hybrid/30_long_reads.smk) and as Polypolish's reads (rule short_read_correction in
# rules/hybrid/50_polish.smk). With parameters.decontamination.mode off these are the pairs
# properly paired on any filtered guide contig; with a genus filter on, the pairs of the
# removed contigs drop out.
#
# Takes in:    the six index files from index_selected_contigs, and the trimmed pairs
#              TRIM_R1 and TRIM_R2 (rule trim_adapters in rules/hybrid/10_reads.smk).
# Does:        1. bowtie2 aligns the pairs to the screened guide in its default end-to-end
#                 mode. --no-unal leaves the reads that did not align out of the SAM, which
#                 keeps it smaller; -S writes the alignments as SAM text to a temporary
#                 file. bowtie2's summary goes to the log.
#              2. samtools fastq turns the alignments back into reads. -f 0x2 keeps only
#                 the reads whose "properly paired" flag is set (0x2 is that flag's value);
#                 -1 and -2 name the two output files, one per mate; -0 /dev/null and -s
#                 /dev/null throw away the reads whose flags say neither first nor second
#                 mate, and the singletons (one mate without the other); -n leaves the read
#                 names as they are, without the /1 and /2 endings samtools would otherwise
#                 add. The two files then carry identical names in the same order, which
#                 Polypolish needs to match the two alignments of a pair.
# Produces:    SEL_R1 and SEL_R2 (01.reads/{sample}/selected/{sample}_selected_R1.fastq,
#              _R2), uncompressed, KEPT: Filtlong and Polypolish read them, and a resumed
#              run must not have to align the library again. Rule all asks for them. The
#              SAM, _GUIDE_SELECTION_SAM, is temp(): Snakemake deletes it when the job is
#              done. Disk: it holds every aligned read as text, tens of gigabytes for a deep
#              library while the job runs.
# Consumed by: rule filter_long_reads (ONT samples with guidance; rules/hybrid/
#              30_long_reads.smk) and rule short_read_correction (every sample;
#              rules/hybrid/50_polish.smk).
# Log:         bowtie2's alignment summary, then samtools' messages.
#
# Threads: all CPUS for bowtie2 (-p); samtools fastq runs on one CPU.
#
# The trimmed pairs are aligned to the guide twice in this file: in map_guide_contigs to
# the filtered guide, for the coverage track, and here to the screened guide, for the
# selection. With mode off the two guides hold the same contigs, so the second alignment
# repeats the first; that is how FunFluxL+ did it, kept as it is.
#
# v1 -> v2: the two commands are FunFluxL+'s (rule select_short_reads), unchanged. The
# selected pairs moved from 03.read_selection/contamination/{sample}/ to the sample's
# selected folder under 01.reads/, the SAM to the guide's contaminants folder, and the CPU
# count is Snakemake's own threads:, which --cores enforces.
rule map_sel_contigs:
    input:
        idx1 = _GUIDE_SELECTED_INDEX_PREFIX + ".1.bt2",
        idx2 = _GUIDE_SELECTED_INDEX_PREFIX + ".2.bt2",
        idx3 = _GUIDE_SELECTED_INDEX_PREFIX + ".3.bt2",
        idx4 = _GUIDE_SELECTED_INDEX_PREFIX + ".4.bt2",
        ridx1 = _GUIDE_SELECTED_INDEX_PREFIX + ".rev.1.bt2",
        ridx2 = _GUIDE_SELECTED_INDEX_PREFIX + ".rev.2.bt2",
        r1 = TRIM_R1,
        r2 = TRIM_R2,
    output:
        sam = temp(_GUIDE_SELECTION_SAM),
        r1 = SEL_R1,
        r2 = SEL_R2,
    params:
        db = _GUIDE_SELECTED_INDEX_PREFIX,
    threads: CPUS
    # Asked of a cluster scheduler: a time limit only, as the step reads every read but
    # holds little in memory (00_common.smk, section 4).
    resources:
        runtime = LONG_RUNTIME_MINUTES,
    conda:
        "../../envs/bowtie.yaml"
    log:
        LOGS + "/map_sel_contigs_{sample}.log"
    shell:
        """
        bowtie2 \
            -x {params.db} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
            --no-unal \
            -S {output.sam} > {log} 2>&1
        samtools fastq \
            -f 0x2 \
            -1 {output.r1} \
            -2 {output.r2} \
            -0 /dev/null \
            -s /dev/null \
            -n {output.sam} >> {log} 2>&1
        """
