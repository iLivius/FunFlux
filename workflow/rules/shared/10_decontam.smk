# FunFlux v2.1.0 — contamination screen of the draft assembly (rules/shared/10_decontam.smk)
#
# What it does: every contig of the draft assembly is searched with BLAST against NCBI nt
# (the core_nt database by default) to learn its taxonomy; BlobTools joins those hits with a
# coverage track into one table; and a Python selector applies the user's keep/drop policy
# and writes down a decision, with the reason, for every contig.
#
# Why. A culture that is called one fungus can still contain other organisms: bacteria or
# yeasts that grew in the same flask, DNA carried over from a neighbouring sample, adapter
# or vector debris. Each of them reaches the assembly as extra contigs. Two kinds of
# evidence per contig can set them apart from the isolate: its taxonomy (the BLAST hits)
# and its read depth (a contaminant is often covered at a different depth from the
# isolate's genome). BlobTools puts both into one table; the selector decides on the
# taxonomy alone. The screen runs BEFORE annotation: a contaminant contig left in the
# assembly would be handed to the fungal gene predictor, which calls false genes on it, and
# it would distort the assembly statistics (QUAST, BUSCO) and the ITS search.
#
# Nothing is removed by default. The fungal default is parameters.decontamination.mode:
# off, so every contig is kept, and every contig still gets a decision in
# contig_taxonomy_decisions.tsv. The reason: a fungus without a close reference genome in
# NCBI nt often gets the wrong genus from BLAST, so a person should read that table before
# any contig is thrown away.
#
# The rules, in the order the data flows (one line each):
#   map_contigs    : makes the coverage track, a BAM file (reads aligned to the draft) from
#                    which BlobTools counts how deeply each contig is covered. (In illumina
#                    mode rule index_contigs, the bowtie2 index of the draft, comes first.)
#   blast_contigs  : megablast of every draft contig against NCBI nt, the taxonomy evidence.
#   blob_json      : BlobTools joins the assembly, the coverage and the BLAST hits in one
#                    database file.
#   blob_table     : BlobTools writes that database out as a flat table: one row per
#                    contig, one taxonomic call per rank.
#   select_contigs : applies the policy, writes the kept contigs and the audit table.
#
# The mitochondrion screen has its own file (shared/15_mitochondria.smk), and its first rule,
# screen_mitochondria, runs between BLAST and the selector. It reads the same BLAST table,
# finds the contig(s) that match a mitochondrial genome and writes their IDs
# (CONFIRMED_MITO_IDS). select_contigs reads that list and keeps those contigs whatever the
# genus filter says. BLAST can give a fungal mitochondrion a genus that is not the
# organism's own, and even a correct genus can differ from the genus most of the
# assembly gets (see rule select_contigs). Either way a genus filter could throw the
# mitochondrion away before it could be set apart. That is why the mitochondrion screen runs
# on the DRAFT, before the selector.
#
#   DRAFT_CONTIGS
#     |
#     +--> map_contigs -----------+      (the coverage track: a BAM file)
#     |                           +--> blob_json --> blob_table --+    (BlobTools: one table,
#     +--> blast_contigs ---------+                               |     taxonomy per contig)
#               |                                                 +--> select_contigs --> DECONTAM_CONTIGS
#               +--> screen_mitochondria --> CONFIRMED_MITO_IDS --+           (and the audit files)
#                    (shared/15_mitochondria.smk)
#
# DRAFT_CONTIGS also goes straight into blob_json, screen_mitochondria and select_contigs.
#
# Where the screen sits in each mode. DRAFT_CONTIGS is what gets screened, and each mode's
# front end declares the rule that makes it: rule filter_contigs in contigs and illumina
# mode, rule long_read_assembly in longread and hybrid mode.
#
#   mode      screened (DRAFT_CONTIGS)                  DECONTAM_CONTIGS is read next by
#   --------  ----------------------------------------  --------------------------------------
#   contigs   your FASTA after filter_contigs           finalize_assembly, which copies it to
#             (contigs_filt.fasta); only when           CONTIGS_ALL
#             parameters.contigs.screen is true
#   illumina  the SPAdes contigs after the length and   finalize_assembly, which copies it to
#             coverage filter (contigs_filt.fasta)      CONTIGS_ALL
#   longread  the hifiasm or Flye assembly              Medaka polishes ONT samples (PacBio HiFi
#             (assembly.fasta)                          samples skip it), then finalize_assembly
#   hybrid    the hifiasm or Flye assembly              Medaka (ONT samples), then Polypolish
#             (assembly.fasta)                          with the Illumina reads, then
#                                                       finalize_assembly
#
# That is why select_contigs writes DECONTAM_CONTIGS and not the final genome: in the
# read-based modes the genome is only finished after polishing. In every mode the last step
# of the screen, rule separate_mitochondria in shared/15_mitochondria.smk, reads CONTIGS_ALL
# and writes the final nuclear genome (FINAL_CONTIGS). In hybrid mode the Illumina "guide"
# assembly is screened by its own copies of these rules, in rules/hybrid/25_guide.smk, not by
# this file.
#
# Everything is inside one `if SCREEN:` block. SCREEN (00_common.smk) is true in the three
# read-based modes and, in contigs mode, only when parameters.contigs.screen is true. The
# `if` is decided once, when Snakemake reads the workflow, before any job runs. When SCREEN
# is false this file defines no rule at all, so a contigs run without the screen has no
# BLAST rule that could ask for a database it was never given.
#
# These names come from 00_common.smk, where each is defined once and never worked out
# again: SCREEN, DECONTAM_READS, DRAFT_CONTIGS, TRIM_R1, TRIM_R2, FILTERED_LONG,
# DECONTAM_INDEX_PREFIX, DECONTAM_BAM, BLASTOUT, BLOB_PREFIX, BLOB_JSON, BLOB_COV,
# BLOB_TABLE_PREFIX, BLOB_TABLE, COMPOSITION, CONTIG_LIST, CONTIG_DECISIONS,
# DECONTAM_CONTIGS, CONFIRMED_MITO_IDS, BLASTDB, NT_VERSION, BLAST_MAX_TARGET_SEQS,
# BLAST_MAX_HSPS, DECONTAMINATION, SELECT_TAXONOMY_SCRIPT, LOGS and CPUS. So are the
# functions capped_cpus() and minimap_preset() (the minimap2 preset of a long-read sample).
#
# conda: environment paths are written relative to THIS file (workflow/rules/shared/), so
# "../../envs/blast.yaml" goes up two folders, to workflow/, and then into workflow/envs/.
# The tools get their own environments (bowtie2 and samtools, minimap2 and samtools, BLAST+,
# BlobTools). The selector is a Python script that uses only the standard library, and
# select_contigs runs it in envs/python.yaml, which holds a pinned Python 3.12 and nothing
# else: activating that environment makes `python` that Python, instead of whatever `python`
# the computer has (a computer with only `python3`, or an old Python, would stop the job).
#
# CPUs: a rule that uses several CPUs declares Snakemake's own `threads:` and passes
# {threads} to its tool. That is what makes `--cores N` actually limit the run.


# Every rule below exists only when the contamination and mitochondrion screen is on
# (SCREEN, see above).
if SCREEN:

    # ───────────────── Coverage track: map_contigs ─────────────────
    # blobtools create stops unless it is given a coverage source. Here that source is a BAM
    # file (reads aligned to the assembly), and map_contigs makes it. WHICH reads are aligned
    # depends on the mode, and so do the aligner and its settings, the conda environment,
    # whether an index rule comes first (illumina only) and whether the BAM is kept. So the
    # rule comes in three versions, each written out as a rule of its own, and only ONE of
    # them exists in any run. The choice is a plain `if` on DECONTAM_READS (00_common.smk),
    # decided when Snakemake reads the workflow; the mode cannot change during a run, so
    # neither can the choice.
    #   "none"      contigs mode: there are no reads, so the contigs are aligned to themselves.
    #   "illumina"  illumina mode: the trimmed read pairs, aligned with bowtie2 (two rules,
    #               index_contigs and map_contigs).
    #   "long"      longread and hybrid mode: the filtered long reads, aligned with minimap2.
    # Whichever version runs, the result has the same name (DECONTAM_BAM), so blob_json and
    # Qualimap read it without knowing which reads were aligned.
    if DECONTAM_READS == "none":

        # ── map_contigs (contigs mode): a deliberately FAKE coverage track ──
        # In this mode the user supplies finished assemblies, so there are no reads. BlobTools
        # still needs a BAM, so every contig is aligned to ITSELF with minimap2 (the draft is
        # both the reference and the query). The BAM exists only so that BlobTools runs; the
        # coverage it yields has no biological meaning. BlobTools estimates coverage as
        # (reads on the contig) x (mean read length) / (contig length). Here every contig is
        # ONE "read" that maps once, to itself, so the figure is roughly the mean contig
        # length divided by the contig's own length: well below 1 for the longest contigs
        # and far above 1 for the shortest ones. Ignore the coverage
        # column of the BlobTools table in this mode. Only the BLAST taxonomy does real
        # work, and the selector reads only the taxonomy columns anyway.
        #
        # Takes in:    DRAFT_CONTIGS, the supplied FASTA after rule filter_contigs
        #              (contigs/20_assembly.smk): contigs_filt.fasta.
        # Does:        minimap2 -a aligns the draft against itself, with no preset, and writes
        #              the alignments as SAM text. samtools view -b -u turns them into an
        #              uncompressed BAM for the pipe (samtools recognises SAM input by
        #              itself). samtools sort orders the BAM by position, and samtools index
        #              writes the .bai index next to it.
        # Produces:    DECONTAM_BAM and its index, DECONTAM_BAM + ".bai", both temporary
        #              (see below).
        # Consumed by: blob_json, and nothing else. There is no Qualimap report in contigs
        #              mode: rule map_evaluation (shared/20_qc.smk), the Qualimap step of the
        #              read-based modes, needs real reads.
        #
        # Why temporary: Snakemake deletes an output wrapped in temp() as soon as the last
        # job that needs it has finished. In the read-based modes the BAM is KEPT (see
        # DECONTAM_BAM in 00_common.smk), because Qualimap reads it as well and aligning tens
        # of millions of reads a second time would be expensive. Here nothing else reads it,
        # there are no reads to align again, and rebuilding it from the contigs takes a few
        # seconds of minimap2 for a fungal genome, so a kept copy would only add a large,
        # meaningless file to the output folder.
        #
        # Threads: the rule reserves all CPUS threads, and the three samtools commands use
        # them ({threads}). minimap2 runs without -t, exactly as in BacFlux, so it uses its
        # own default of 3 threads.
        #
        # v1 -> v2: new. Funnotator (the v1 contigs workflow) had no contamination screen.
        # The command is the contigs-mode map_contigs of BacFlux, the sister bacterial
        # workflow, with two changes: samtools view also writes its messages to the log
        # (2>> {log}), so that nothing in the pipe can print to the screen; and BacFlux's
        # samtools view -S is gone. samtools 1.21, the version pinned in
        # workflow/envs/minimap.yaml, ignores -S ("Ignored for compatibility with previous
        # samtools versions", its manual says), so dropping it changes no alignment.
        rule map_contigs:
            input:
                contigs = DRAFT_CONTIGS,
            output:
                bam = temp(DECONTAM_BAM),
                bai = temp(DECONTAM_BAM + ".bai"),
            threads: CPUS
            conda:
                "../../envs/minimap.yaml"
            log:
                LOGS + "/map_contigs_{sample}.log"
            shell:
                """
                minimap2 \
                  -a {input.contigs} \
                  {input.contigs} 2> {log} | \
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

    elif DECONTAM_READS == "illumina":

        # ── index_contigs (illumina mode): the bowtie2 index of the draft ──
        # bowtie2 cannot align reads to a plain FASTA; it first needs its own index of the
        # reference, built once per draft. The index is six files that bowtie2 names itself
        # from the prefix it is given (<prefix>.1.bt2 to <prefix>.4.bt2 and <prefix>.rev.1.bt2,
        # <prefix>.rev.2.bt2), so the rule declares the six names and hands the prefix to
        # bowtie2-build.
        #
        # Takes in:    DRAFT_CONTIGS, the SPAdes contigs after the length and coverage filter
        #              (rule filter_contigs in rules/illumina/20_assembly.smk): contigs_filt.fasta.
        # Does:        bowtie2-build -f (the input is FASTA) writes the index under
        #              DECONTAM_INDEX_PREFIX, in the contaminants folder, next to the BAM.
        # Produces:    the six .bt2 files, all temporary: only map_contigs reads them, and
        #              rebuilding them takes about a minute for a fungal genome.
        # Consumed by: map_contigs below.
        #
        # Threads: one. bowtie2-build could use several (--threads), but a 40 Mb genome
        # indexes in about a minute, so the rule declares no `threads:` and the run keeps its
        # CPUs for the mapping and the BLAST.
        #
        # v1 -> v2: the command is unchanged; the index now lies in the sample's contaminants
        # folder.
        rule index_contigs:
            input:
                contigs = DRAFT_CONTIGS,
            output:
                idx1 = temp(DECONTAM_INDEX_PREFIX + ".1.bt2"),
                idx2 = temp(DECONTAM_INDEX_PREFIX + ".2.bt2"),
                idx3 = temp(DECONTAM_INDEX_PREFIX + ".3.bt2"),
                idx4 = temp(DECONTAM_INDEX_PREFIX + ".4.bt2"),
                ridx1 = temp(DECONTAM_INDEX_PREFIX + ".rev.1.bt2"),
                ridx2 = temp(DECONTAM_INDEX_PREFIX + ".rev.2.bt2"),
            params:
                basename = DECONTAM_INDEX_PREFIX,
            conda:
                "../../envs/bowtie.yaml"
            log:
                LOGS + "/index_contigs_{sample}.log"
            shell:
                """
                bowtie2-build -f {input.contigs} {params.basename} > {log} 2>&1
                """

        # ── map_contigs (illumina mode): the trimmed read pairs aligned to the draft ──
        # The coverage track with real meaning: how many of the sample's own Illumina reads
        # land on each contig. The contigs of the isolate share one read depth, while a
        # contaminant that grew to another cell density sits at a depth of its own, which is
        # what the BlobTools plot (GC against coverage) separates.
        #
        # Takes in:    the six index files from index_contigs, and the trimmed, PhiX-free
        #              read pairs TRIM_R1 and TRIM_R2 (rule trim_adapters in
        #              rules/illumina/10_reads.smk). bowtie2 is given the index PREFIX
        #              (params db), not the files; the files are listed as inputs so that
        #              Snakemake builds the index first.
        # Does:        bowtie2 aligns the pairs in its default end-to-end mode (-p: the
        #              threads; -t prints the time each step took to the log) and writes the
        #              alignments as SAM text. samtools view -b -u turns them into an
        #              uncompressed BAM for the pipe, samtools sort orders the BAM by position,
        #              and samtools index writes the .bai index next to it.
        # Produces:    DECONTAM_BAM and its index, DECONTAM_BAM + ".bai". Both are KEPT, not
        #              temporary: Qualimap (rule map_evaluation in shared/20_qc.smk) and
        #              BlobTools read the same BAM, and when a later step is re-run (a changed
        #              BLAST depth re-runs blob_json, for instance) a temporary BAM would have
        #              to be rebuilt by aligning tens of millions of reads again.
        # Consumed by: blob_json (the coverage) and map_evaluation (the mapping statistics).
        #
        # Threads: all CPUS, for bowtie2 (-p) and the three samtools commands (-@).
        #
        # v1 -> v2: bowtie2's command is unchanged, and so is the shape of the pipe. The three
        # samtools commands are written as in the contigs-mode version above: samtools view -b
        # -u (uncompressed BAM; v1 wrote -hbS, whose -h is implied by BAM output and whose -S
        # is ignored by samtools 1.21), every samtools command writes its messages to the log
        # (2>> {log}), and the index prefix is no longer wrapped in temp(), which did nothing
        # on a param. The BAM moved from 03.post-processing/{sample}/ to the sample's
        # contaminants folder.
        rule map_contigs:
            input:
                idx1 = DECONTAM_INDEX_PREFIX + ".1.bt2",
                idx2 = DECONTAM_INDEX_PREFIX + ".2.bt2",
                idx3 = DECONTAM_INDEX_PREFIX + ".3.bt2",
                idx4 = DECONTAM_INDEX_PREFIX + ".4.bt2",
                ridx1 = DECONTAM_INDEX_PREFIX + ".rev.1.bt2",
                ridx2 = DECONTAM_INDEX_PREFIX + ".rev.2.bt2",
                r1 = TRIM_R1,
                r2 = TRIM_R2,
            output:
                bam = DECONTAM_BAM,
                bai = DECONTAM_BAM + ".bai",
            params:
                db = DECONTAM_INDEX_PREFIX,
            threads: CPUS
            # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
            resources:
                mem_mb = HEAVY_MEM_MB,
                runtime = LONG_RUNTIME_MINUTES,
            conda:
                "../../envs/bowtie.yaml"
            log:
                LOGS + "/map_contigs_{sample}.log"
            shell:
                """
                bowtie2 \
                  -x {params.db} \
                  -1 {input.r1} \
                  -2 {input.r2} \
                  -p {threads} \
                  -t 2> {log} | \
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

    elif DECONTAM_READS == "long":

        # ── map_contigs (longread and hybrid mode): the filtered long reads aligned to the draft ──
        # The coverage track: how deeply the sample's own long reads cover each contig of the
        # assembly they were assembled into. The contigs of the isolate share one depth; the
        # mitochondrion is many times deeper (many copies per cell), and a contaminant
        # usually sits at a depth of its own, which is what the BlobTools plot separates.
        #
        # Takes in:    FILTERED_LONG, the Filtlong output (rule filter_long_reads in the
        #              mode's front end), and DRAFT_CONTIGS, the hifiasm or Flye assembly
        #              (rule long_read_assembly): assembly.fasta.
        # Does:        minimap2 -a writes the alignments as SAM text; -x takes the preset of
        #              the sample's technology (params preset, from minimap_preset() in
        #              00_common.smk): map-ont for Nanopore reads, map-hifi for PacBio HiFi
        #              reads. The two have different error profiles, and minimap2 tunes its
        #              seeding and scoring to the preset, so in a mixed batch each sample gets
        #              its own. samtools then makes the sorted, indexed BAM exactly as in the
        #              two versions above.
        # Produces:    DECONTAM_BAM and its index, both KEPT (Qualimap reads the BAM too; the
        #              illumina version above says why a temporary BAM would not do).
        # Consumed by: blob_json (the coverage) and rule map_evaluation in shared/20_qc.smk
        #              (the mapping statistics).
        #
        # Threads: all CPUS, for minimap2 (-t) and the three samtools commands (-@). The
        # contigs-mode self-map above keeps minimap2's default of 3 threads on purpose, because
        # it is a seconds-long job; this one aligns a whole long-read set, so it gets every
        # thread.
        #
        # v1 -> v2: FunFluxL ran minimap2 WITHOUT -t, so minimap2 used its default 3 threads
        # while samtools got all CPUs; -t {threads} is new and changes no alignment, only the
        # time it takes. samtools view lost its -S (ignored by samtools 1.21) and gained
        # 2>> {log}; the rest of the command is unchanged, the preset was already per sample,
        # and the BAM was already kept. It moved from 03.post-processing/contaminants/ to the
        # sample's contaminants folder.
        rule map_contigs:
            input:
                reads = FILTERED_LONG,
                contigs = DRAFT_CONTIGS,
            output:
                bam = DECONTAM_BAM,
                bai = DECONTAM_BAM + ".bai",
            params:
                preset = minimap_preset,
            threads: CPUS
            # Asked of a cluster scheduler: memory and a time limit (00_common.smk, section 4).
            resources:
                mem_mb = HEAVY_MEM_MB,
                runtime = LONG_RUNTIME_MINUTES,
            conda:
                "../../envs/minimap.yaml"
            log:
                LOGS + "/map_contigs_{sample}.log"
            shell:
                """
                minimap2 \
                  -ax {params.preset} \
                  -t {threads} \
                  {input.contigs} \
                  {input.reads} 2> {log} | \
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


    # ───────────────── Taxonomy evidence: blast_contigs ─────────────────
    # BLAST asks, for every contig, which sequences in the NCBI nucleotide database (nt) look
    # like it, and says which organism each of them came from. BlobTools sums those matches
    # into a taxonomic call per contig; the mitochondrion screen reads the same matches to
    # find the contig that is a mitochondrial genome.
    #
    # Takes in:    DRAFT_CONTIGS, the draft assembly. The screen works on the draft on
    #              purpose: it has to see every contig it is asked to judge, including the
    #              ones the selector is about to drop.
    #              The database folder comes from directories.blast_db (BLASTDB) and the
    #              database name from parameters.nt_version (NT_VERSION, core_nt unless
    #              changed).
    # Does:        one `blastn -task megablast` run per sample. megablast is BLAST's mode for
    #              finding sequences that are highly similar to something already known. It
    #              is much faster than the more sensitive blastn task, and an assembly of an
    #              organism that is already in the database is expected to match it closely,
    #              which is the case megablast is built for.
    # Produces:    BLASTOUT: a plain-text table with 15 tab-separated columns, one line per
    #              alignment BLAST kept, with the subject title (stitle) last.
    # Consumed by: blob_json, and rule screen_mitochondria in shared/15_mitochondria.smk.
    #
    # Hit depth (docs/methods/blast-hit-depth.md explains the choice). Two numbers control how
    # much BLAST keeps per contig, and they are NOT two ways of saying "how many hits":
    #   -max_target_seqs  how many different database sequences ("subjects") to keep per
    #                     contig. BLAST applies it DURING the search, not as a final "best N"
    #                     filter, so a small value can lose better matches; BLAST itself
    #                     recommends 5 or more. FunFlux uses 50.
    #   -max_hsps         how many alignments to keep per subject (an HSP, a high-scoring
    #                     segment pair, is one local alignment between the contig and the
    #                     subject). This one decides the mitochondrion screen: BLAST works
    #                     out the column qcovs ("query coverage per subject", the percent of
    #                     the contig covered by that subject) from the alignments it KEPT.
    #                     With 1 it keeps only the best block, so a complete mitochondrial
    #                     genome that matches its closest reference in several blocks looks
    #                     only partly covered and can fall below the screen's threshold.
    #                     FunFlux uses 5.
    # Both come from the configuration (parameters.blast_max_target_seqs and
    # parameters.blast_max_hsps) and reach the command through `params:`, not as a global
    # name typed into the command. Snakemake compares params with the previous run and
    # re-runs BLAST when one changed; it would not notice a changed global inside the command
    # text, because it compares only the text as written.
    #
    # Do NOT change the -outfmt string (the list of columns). Both programs that read this
    # table pick columns by POSITION: BlobTools takes the contig name, the taxonomy ID and the
    # bit score (columns 1, 2 and 3), and the mitochondrion screen takes the contig name,
    # qlen, slen, qcovs and stitle (columns 1, 7, 8, 9 and 15). Moving or dropping a column
    # can break both without any error message.
    #
    # BLASTDB={params.dir} points BLAST at the folder that also holds its taxonomy files
    # (taxdb.btd and taxdb.bti); without them the organism-name columns come out as N/A.
    # -db names the database itself, inside that folder. -evalue 1e-5 ignores matches that
    # chance alone would produce more than once in 100,000 searches.
    #
    # v1 -> v2: v1 FunFlux ran with -max_target_seqs 10 -max_hsps 1 and the database name
    # core_nt typed into the rule. Both depths are now configurable and default to 50 and 5
    # (see above), and the database name is parameters.nt_version. CPU use is declared with
    # `threads:` (capped at 24, as before) instead of the custom `resources: cpus`, which
    # --cores did not count.
    rule blast_contigs:
        input:
            contigs = DRAFT_CONTIGS,
        output:
            blast = BLASTOUT,
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
            LOGS + "/blast_contigs_{sample}.log"
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


    # ───────────────── BlobTools: blob_json and blob_table ─────────────────
    # BlobTools (version 1.1.1 here) was written to draw "blob plots": contigs as bubbles, GC
    # against coverage, coloured by taxon, which is how a mixed culture shows up as separate
    # clouds. Here its job is the table behind such a plot. blob_json builds a database (a
    # JSON file) that holds, for every contig, its length and GC, its coverage from the BAM
    # and a taxonomy worked out from the BLAST hits. blob_table writes that database out as a
    # flat table.

    # ── blob_json: one database from assembly + coverage + BLAST ──
    # Takes in:    DRAFT_CONTIGS, the same contigs that were aligned and BLASTed;
    #              DECONTAM_BAM and its index (from map_contigs), the coverage track;
    #              BLASTOUT (from blast_contigs), the taxonomy evidence; and nodes.dmp and
    #              names.dmp from the BLAST database folder (directories.blast_db). Those two
    #              files belong to the NCBI taxonomy dump (the taxdump archive): nodes.dmp is
    #              the family tree of all taxa and names.dmp their names, which BlobTools
    #              needs to turn a taxonomy ID into "Fusarium" or "Ascomycota". They are listed
    #              as INPUTS, not hidden in the command, so that a missing file is reported
    #              when the run starts and not an hour into it. (As a result even a dry run
    #              needs the two files to exist; empty stand-ins are enough.)
    # Does:        blobtools create, which writes its results under the prefix BLOB_PREFIX
    #              and adds the file endings itself.
    # Produces:    BLOB_JSON (the database) and BLOB_COV (the coverage per contig that
    #              BlobTools read from the BAM, as a text file). BlobTools names the coverage
    #              file after the BAM's file name, which is why 00_common.smk builds BLOB_COV
    #              from DECONTAM_BAM instead of typing the name a second time. Both are
    #              temp(): Snakemake deletes BLOB_JSON once blob_table has read it, and
    #              BLOB_COV, which no rule reads, as soon as this job has finished.
    # Consumed by: blob_table (BLOB_JSON only).
    #
    # v1 -> v2: v1 FunFlux sent this rule's messages to the log of blob_table, and this rule
    # overwrote it first. Each rule now has its own log. Everything else in the command is
    # unchanged.
    rule blob_json:
        input:
            contigs = DRAFT_CONTIGS,
            bam = DECONTAM_BAM,
            bai = DECONTAM_BAM + ".bai",
            blast = BLASTOUT,
            nodes = os.path.join(BLASTDB, "nodes.dmp"),
            names = os.path.join(BLASTDB, "names.dmp"),
        output:
            json = temp(BLOB_JSON),
            cov = temp(BLOB_COV),
        params:
            basename = BLOB_PREFIX,
        conda:
            "../../envs/blobtools.yaml"
        log:
            LOGS + "/blob_json_{sample}.log"
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

    # ── blob_table: one taxonomic call per contig and rank ──
    # A contig usually has many BLAST hits, often to several different organisms. The
    # "bestsum" rule adds up the bit scores (BLAST's measure of how good a match is) of all
    # hits per taxon, separately at each rank, and the taxon with the highest total wins that
    # rank. The result is one row per contig that names its winner from superkingdom down to
    # species, with the length, GC and coverage beside it.
    #
    # Takes in:    BLOB_JSON, the database made by blob_json.
    # Does:        blobtools view with --taxrule bestsum (the rule above), --rank all (write
    #              the columns of every rank, not only one) and --hits (also list which taxa
    #              contributed which score). BlobTools builds the file name itself: the
    #              prefix BLOB_TABLE_PREFIX plus the ending .blob.blobDB.table.txt.
    # Produces:    BLOB_TABLE. It is kept, not temporary: it is the evidence behind every
    #              keep-or-drop decision, and worth re-reading.
    # Consumed by: select_contigs, which reads the contig name and the genus column.
    #
    # v1 -> v2: v1 appended this rule's messages (>>) to the log it shared with blob_json.
    # With a log of its own the rule now starts a fresh one (>), so a re-run does not pile
    # new messages onto old ones.
    rule blob_table:
        input:
            json = BLOB_JSON,
        output:
            bestscore = BLOB_TABLE,
        params:
            basename = BLOB_TABLE_PREFIX,
        conda:
            "../../envs/blobtools.yaml"
        log:
            LOGS + "/blob_table_{sample}.log"
        shell:
            """
            blobtools view \
              --input {input.json} \
              --out {params.basename} \
              --taxrule bestsum \
              --rank all \
              --hits > {log} 2>&1
            """


    # ───────────────── Contig selection and audit trail: select_contigs ─────────────────
    # Where contigs are actually kept or dropped. The selector script reads the BlobTools
    # table, takes each contig's genus, and applies the policy in
    # parameters.decontamination.mode:
    #   off      keep everything (the audit files are still written). The fungal default.
    #   include  keep only the listed genera.
    #   exclude  drop only the listed genera (exact names; upper and lower case are not
    #            told apart).
    # plus discard_no_hit, which, when true, drops the contigs that BLAST could not place at
    # all (in include and exclude mode; mode off keeps them). Whatever the mode, a contig
    # whose ID is on the protected list is kept (see below). In include mode genus names are
    # also matched through a few built-in synonyms and name prefixes that were written for
    # bacteria; the script says what they can do to a fungal name.
    #
    # There is no mode auto, BacFlux's default (keep the genus that holds the most contigs),
    # because BLAST genus labels are unreliable for poorly sequenced fungi
    # (docs/analysis/decontamination.md, "Why there is no auto mode"). Many short contigs
    # can get no hit at all, and when these no-hit contigs are the largest group, auto would
    # keep them and remove most of the genome. With discard_no_hit true, the genus with the
    # most contigs can be a mislabel (a better-sequenced relative of the isolate), and auto
    # would still remove every contig of the isolate that BLAST labelled with another
    # genus. 00_common.smk refuses mode auto, also in the sample_overrides table, when it
    # reads the configuration, before any job starts. The script refuses it again if it
    # still receives it (a run by hand, or a table edited while the run is going); its
    # docstring has the details.
    #
    # Takes in:    BLOB_TABLE, the per-contig taxonomy (from blob_table); DRAFT_CONTIGS, the
    #              sequences themselves, read directly (the script reads wrapped FASTA, so no
    #              one-line-per-record copy is needed); and CONFIRMED_MITO_IDS, the contigs
    #              the mitochondrion screen confirmed (rule screen_mitochondria in
    #              shared/15_mitochondria.smk). Those contigs are kept whatever their genus.
    #              A mitochondrion's genus can differ from the genus most of the assembly
    #              gets, and then a filter that keeps only that genus would drop it. This
    #              happens, for example, when the isolate's genus is poorly represented in
    #              nt: its mitochondrion can match a complete mitochondrial genome of its
    #              own, true genus, while most of the nuclear DNA is labelled with a
    #              better-sequenced relative. The audit table gives the reason
    #              protected_confirmed_mitochondrial. Without the screen's list these
    #              contigs would be exposed to the genus filter, and a contig lost here
    #              cannot be recovered by the later split.
    # Does:        python select_contigs_by_taxonomy.py (SELECT_TAXONOMY_SCRIPT) with the
    #              settings of parameters.decontamination, which 00_common.smk checks once and
    #              turns into plain text (DECONTAMINATION), so that every sample of a run
    #              follows one policy. Three of the settings are paths to files
    #              (include_genera_by_sample, exclude_genera_file, sample_overrides). The two
    #              per-sample ones let one awkward sample deviate without changing the policy
    #              for the rest of the batch. A setting left empty reaches the script as a
    #              bare flag: the :q in {params.include_genera:q} quotes a value for the
    #              shell when it holds spaces or special characters, but writes an empty
    #              value as nothing. The script reads a bare flag as "not set".
    # Produces:    four files, all in the sample's contaminants folder:
    #                COMPOSITION       one line per genus with its share of the assembly's
    #                                  DNA and of the contig count, over EVERY contig in the
    #                                  table, kept or dropped.
    #                CONTIG_LIST       the IDs of the kept contigs, one per line.
    #                CONTIG_DECISIONS  the audit table: for every contig its genus, whether
    #                                  it was kept or removed, and the reason. Every
    #                                  filtering decision in this workflow has to be
    #                                  traceable, and this is the trace.
    #                DECONTAM_CONTIGS  the kept sequences, wrapped at 80 characters per line.
    # Consumed by: DECONTAM_CONTIGS by the next step of the assembly chain, which depends on
    #              the mode (see the table at the top of this file). In contigs mode that is
    #              finalize_assembly, which copies it to CONTIGS_ALL. No rule reads the other
    #              three files; they are records for people.
    #
    # The script ends with an error, so that the rule counts as failed and Snakemake deletes
    # its outputs, when no contig is kept or when a protected ID is missing from the FASTA
    # (which would mean the list and the FASTA come from different assemblies). It also stops
    # on input it cannot use, for example a BlobTools table or FASTA with no records in it, a
    # sample_overrides or include_genera_by_sample table without its required columns, or
    # mode include with no genus given. 00_common.smk checks the last two when it reads the
    # configuration, before any job starts; the script checks them again, as it does auto.
    # The script prints two warnings that change no decision: two or more genera each hold
    # at least 5 % of the DNA, or more than 20 % of the DNA is being removed. On the error
    # stream it warns when the BlobTools table names contigs that are not in the FASTA. The
    # log catches all of it.
    #
    # v1 -> v2: v1 had no log, so those warnings only reached the screen and were lost in a
    # big batch; they now go to the log (open it to see them). The output moved from
    # contigs_sel.fasta to contaminants/contigs_decontam.fasta. Mode auto, which v1 offered,
    # is gone (see above). The script runs with the pinned Python of envs/python.yaml, where
    # v1 used the Python that started Snakemake. The other settings, the order of the flags
    # and the --always-keep flag are unchanged.
    rule select_contigs:
        input:
            bestscore = BLOB_TABLE,
            contigs = DRAFT_CONTIGS,
            confirmed_ids = CONFIRMED_MITO_IDS,
        output:
            abund = COMPOSITION,
            list = CONTIG_LIST,
            decisions = CONTIG_DECISIONS,
            contigs = DECONTAM_CONTIGS,
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
            LOGS + "/select_contigs_{sample}.log"
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
