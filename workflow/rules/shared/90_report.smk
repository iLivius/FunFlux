# FunFlux v2.0.0 — Stage 90: the report (rules/shared/90_report.smk)
#
# One rule: gather the quality-control results of the whole run into a single MultiQC page
# (06.report/multiqc_report.html), with sample names a person can read. MultiQC is a program
# that searches folders for the output of many bioinformatics tools, recognises each tool by
# its file names and contents, and draws their numbers side by side for all samples.
#
# The Snakefile includes this file LAST, after every file whose rules make a QC result. No
# other rule reads the report; rule all asks for it.
#
# What goes in depends on what the mode produced, so the input list is built from the
# capability flags of 00_common.smk (HAS_SHORT_READS, HAS_LONG_READS, HAS_READS, IS_HYBRID) and
# not from the mode name:
#
#   QC result                 illumina   longread   hybrid    contigs   made by
#   ------------------------  --------   --------   -------   -------   --------------------------
#   fastp (read QC)               x          -         x         -      rule trim_adapters
#   NanoPlot (read QC)            -          x         x         -      rules raw_long_read_qc and
#                                                                       filtered_long_read_qc
#   Qualimap (mapping QC)         x          x       x (2)       -      rule map_evaluation; in hybrid
#                                                                       mode also map_guide_evaluation
#   QUAST (assembly QC)           x          x       x (2)       x      rule genome_assembly_evaluation
#   BUSCO (completeness)          x          x         x         x      rule completeness
#   Contig summary (table)        x          x         x         x      rule contig_summary
#   Adapters (table)              x          x         x         x      rule screen_adapters
#   AMR mutations (table)        (3)        (3)       (3)       (3)     rule amr_chroquetas
#                                                      (2) = the Illumina "guide" assembly gets its own
#                                                            QUAST label and its own Qualimap report
#                                                      (3) = only when parameters.amr.chroquetas is on
#
# The contig summary (shared/25_contig_summary.smk) is not the output of a tool MultiQC knows:
# it is a table FunFlux writes itself, behind a few header lines that tell MultiQC what it is
# ("custom content", in MultiQC's words), one row per contig of the draft assembly, named
# "<sample>: <contig>". MultiQC shows the tables of all samples as one table. The Adapters table
# (shared/25_contig_summary.smk) is custom content too, one row per sample, named after the
# sample: the ONT adapter trimming of the reads and the adapter found on the contigs. The AMR table
# (shared/50_amr.smk) is custom content in the same way: one row per target protein of the
# species, named "<sample>: <target>", with the mutations found, a verdict on resistance, and
# the drugs with evidence for and against it.
#
# Four details here go wrong quietly, or with a misleading message, so each is explained where it
# is used:
#   1. MultiQC is started from the output folder with relative paths (see params in rule multiqc).
#      If wrong: the sample names are not rewritten.
#   2. The report and the data folder must have the names the workflow declares as outputs.
#      If wrong: MultiQC finishes, and Snakemake then reports the outputs as missing.
#   3. BUSCO's kept result tree contains copies of the summary files that must be ignored.
#      If wrong: extra BUSCO rows.
#   4. MultiQC's own shortening of sample names must be off (fn_clean_sample_names, explained in
#      the command of rule multiqc). If wrong: a sample whose name contains a file ending MultiQC
#      knows, such as SAMPLE1.clean or X.sorted, gets a wrong name, and two samples can merge.
#
# Everything referenced here is defined once in 00_common.smk and never worked out again: OUT (the
# output folder), DIR_REPORT, MULTIQC_REPORT, FASTP_JSON, NANOPLOT_RAW_DIR, NANOPLOT_FILTERED_DIR,
# QUALIMAP_DIR, GUIDE_QUALIMAP_DIR, QUAST_DIR, BUSCO_DIR, CONTIG_SUMMARY_MQC, ADAPTER_SCREEN_MQC, AMR_SUMMARY_MQC,
# AMR_ENABLED, SAMPLES, LOGS and the capability flags. A name that 00_common.smk does not have is
# defined below, starting with an underscore.
#
# conda: environment paths are written relative to THIS file (workflow/rules/shared/), so
# "../../envs/multiqc.yaml" climbs shared/ -> rules/ -> workflow/ and then goes into envs/.

# The two other outputs of rule multiqc, next to the page: the configuration file the rule writes
# for MultiQC (kept, so that the rename rules can be checked afterwards), and MultiQC's data
# folder, with the numbers behind the plots.
_MULTIQC_CONFIG = DIR_REPORT + "/multiqc_config.yaml"
_MULTIQC_DATA = DIR_REPORT + "/multiqc_data"


# ───────────── Which result files and folders MultiQC reads, in THIS mode ─────────────
# Takes in:    the capability flags, AMR_ENABLED and SAMPLES (the sample names found in the input
#              folder) from 00_common.smk.
# Does:        lists, for every sample, the files and folders that hold the QC results of the
#              active mode, in the order the table above has them: the read QC (fastp's JSON
#              report, NanoPlot's two folders), the mapping QC (the Qualimap folder, and in
#              hybrid mode the guide's Qualimap folder too), then QUAST_DIR, BUSCO_DIR, the
#              contig summary table, the Adapters table and, when the AMR screen is on, the AMR
#              table. MultiQC
#              finds the right files inside the folders. Each group is added only when the run
#              makes it: asking for a file that no rule makes would stop the workflow.
# Produces:    a list of absolute paths. It is used twice: as the input of rule multiqc, which
#              makes Snakemake build them first, and, turned into paths relative to the output
#              folder, as the paths handed to MultiQC (see _relative_to_output_root).
# Called by:   the line just below it, once, when Snakemake reads this file.
def multiqc_qc_inputs():
    paths = []

    # expand() is a Snakemake helper: it fills {sample} in with every name in SAMPLES, giving
    # one path per sample.
    if HAS_SHORT_READS:
        # fastp's report of the Illumina reads (rule trim_adapters): the JSON holds the numbers.
        paths += expand(FASTP_JSON, sample=SAMPLES)
    if HAS_LONG_READS:
        # NanoPlot on the raw and on the filtered long reads (rules raw_long_read_qc and
        # filtered_long_read_qc): MultiQC reads the NanoStats table in each folder.
        paths += expand(NANOPLOT_RAW_DIR, sample=SAMPLES)
        paths += expand(NANOPLOT_FILTERED_DIR, sample=SAMPLES)
    if HAS_READS:
        # Qualimap on the reads mapped to the draft (rule map_evaluation in shared/20_qc.smk).
        paths += expand(QUALIMAP_DIR, sample=SAMPLES)
    if IS_HYBRID:
        # Qualimap on the Illumina reads mapped to the guide assembly (rule
        # map_guide_evaluation in rules/hybrid/25_guide.smk).
        paths += expand(GUIDE_QUALIMAP_DIR, sample=SAMPLES)
    paths += expand(QUAST_DIR, sample=SAMPLES)
    paths += expand(BUSCO_DIR, sample=SAMPLES)
    # The contig summary of every sample (rule contig_summary in shared/25_contig_summary.smk),
    # in the copy written for MultiQC: the table behind header lines that name the section
    # and the plot type. MultiQC reads it as it reads a tool's output, so listing the file is
    # all that is needed. Its rows are named "<sample>: <contig>" by the script itself, and
    # MultiQC 1.33 keeps the row names written inside a custom-content table as they are:
    # -d (see rule multiqc) prefixes only the sample names it takes from FILE names. So the
    # rename lines below need no entry for this table.
    paths += expand(CONTIG_SUMMARY_MQC, sample=SAMPLES)
    # The Adapters table: one row per sample from the contig adapter screen (rule
    # screen_adapters in shared/25_contig_summary.smk), with the read-level numbers of the ONT
    # samples whose reads fastplong trimmed. Read the same way as the contig summary: the row
    # names are the sample names, written by the script, and the rows of all samples merge
    # into one table because their files share the table id. fastplong's own JSON report is not
    # listed: MultiQC's fastp module would read it as a fastp run.
    paths += expand(ADAPTER_SCREEN_MQC, sample=SAMPLES)
    if AMR_ENABLED:
        # The AMR table of every sample (rule amr_chroquetas in shared/50_amr.smk), in its copy
        # written for MultiQC. Read the same way as the contig summary just above: the script
        # names the rows "<sample>: <target>" itself, so no rename line is needed either.
        paths += expand(AMR_SUMMARY_MQC, sample=SAMPLES)
    return paths


_MULTIQC_INPUTS = multiqc_qc_inputs()


# Takes in a list of absolute paths (_MULTIQC_INPUTS) and returns the same paths relative to the
# output folder OUT, e.g. <OUT>/02.assembly/SAMPLE1/eval/quast -> 02.assembly/SAMPLE1/eval/quast.
# Called by: params qc_rel of rule multiqc. These short paths work because the rule changes into
# OUT before it starts MultiQC; the params of rule multiqc explain why short paths are needed.
def _relative_to_output_root(paths):
    relative_paths = []
    for path in paths:
        relative_paths.append(os.path.relpath(path, OUT))
    return relative_paths


# ───────────── The sample-name rewrite rules, one literal block per mode ─────────────
# Why they are needed. With -d (see rule multiqc), MultiQC names each result after the folders it
# was found in, followed by the tool's own name for the sample, all joined with " | ". Left alone,
# the page would show
#       02.assembly | SAMPLE1 | eval | quast | SAMPLE1
#       02.assembly | SAMPLE1 | eval | busco | short_summary.specific.hypocreales_odb12.SAMPLE1.txt
# The rules below, written into MultiQC's configuration file, turn these into
#       QUAST assembly QC | SAMPLE1
#       BUSCO completeness QC (hypocreales_odb12) | SAMPLE1
# One label scheme for every tool: "<what it is> | <sample>".
#
# How to read one line. Each line is a pair  'pattern': 'new name'.
#   * The pattern is a regular expression (a text-matching language). MultiQC uses a line only
#     when its pattern matches the WHOLE name (sample_names_replace_exact, set in rule multiqc).
#   * ^ marks the start of the name and $ its end. \. is a real dot and \| a real "|" (without
#     the backslash, both have a special meaning in a pattern).
#   * [^|]+ means "one or more characters other than |". The round brackets around it keep what
#     it matched: in ([^|]+), the sample folder, e.g. SAMPLE1.
#   * \1 inside the pattern means "the text the first round brackets kept, again". So a name is
#     only rewritten when it ends with the same sample as its folder: QUAST's label for the
#     assembly, or the end of BUSCO's file name. A name that does not fit is left as it is.
#   * In the new name, \1 puts the kept sample name back, and \2 puts back the text of the
#     second pair of brackets that keeps text: the BUSCO lineage. A pair that starts with ?:,
#     such as (?:generic|specific), keeps nothing and is not counted.
#   Worked example, the QUAST line: in "02.assembly | SAMPLE1 | eval | quast | SAMPLE1",
#   ([^|]+) keeps "SAMPLE1", \1 matches the SAMPLE1 at the end, and the new name is
#   "QUAST assembly QC | SAMPLE1".
#
#   QUAST    the name at the end is the label the QUAST rule gave the assembly in report.tsv: the
#            sample name (quast_labels() in shared/20_qc.smk). In hybrid mode the same report
#            has a second row for the guide assembly, labelled <sample>_guide, so the hybrid
#            block has a second QUAST line that ends in \1_guide$ and gives the label
#            "QUAST guide assembly QC".
#   BUSCO    one row per summary file. BUSCO writes a "generic" summary (the broad lineage it
#            tried first) and a "specific" one (the precise lineage it ended with), and puts the
#            lineage name in the file name, so the lineage in brackets tells the two rows apart.
#            With a lineage chosen in the configuration there is only the specific row.
#            In the pattern, (?:generic|specific) accepts either word, and ([^.|]+_odb[0-9]+)
#            keeps the lineage: a name without dots that ends in _odb and a number, e.g.
#            hypocreales_odb12. MultiQC's own shortening of names is off (fn_clean_sample_names,
#            see rule multiqc), so the name ends with the whole file name, which ends with the
#            sample and ".txt": \.\1\.txt$ at the end of the pattern matches exactly that.
#   fastp    MultiQC does not take fastp's sample name from a file name but from the read files
#            on the command line stored inside the JSON (-i and -I): it shortens each to its
#            file name and, with the shortening of names off, cannot strip the _R1/_R2 endings,
#            so it joins the two names into one long one. The pattern therefore keeps the sample
#            from the folder, 01.reads | <sample> | illumina |, and accepts whatever follows
#            (.+$). The loose ending is safe: the fastp JSON is the only file in that folder
#            that MultiQC recognises; the rest are FASTQ files and fastp's HTML page.
#   NanoPlot MultiQC recognises the NanoStats table by its content and names the row after the
#            file, <sample>_NanoStats.txt (NanoPlot is run with --prefix <sample>_, and the
#            shortening that would cut _NanoStats.txt off is switched off). One line per
#            folder: raw_qc for the raw reads, filtered_qc for the filtered reads.
#   Qualimap one Qualimap report gives several rows: the summary row is named after the BAM
#            file (<sample>_map.bam), the histogram rows (coverage, insert size, GC) after the
#            folder two levels above them, which here is the folder "qualimap" itself. The
#            pattern keeps the sample from the path and accepts any ending (.+$). In hybrid
#            mode the guide's report sits under <sample> | guide | eval | qualimap, so its line
#            is told apart by the extra guide folder and labelled "mapping Illumina QC" (the
#            mapped reads are the Illumina ones; the plain "mapping QC" row is the long reads).
#
# Each mode gets its OWN complete block of lines, written out in full (nothing is generated), so
# that a reader sees exactly the names a mode will get; the QUAST and BUSCO lines are the same in
# every block. The conditions below test the capability flags and not the mode name, but each
# branch is exactly one mode:
#    no reads = contigs,  IS_HYBRID = hybrid,  short reads only = illumina,  long reads only = longread.
# The lines are written for the full names: MultiQC's own shortening is off in every mode, so a
# tool's name for a sample keeps its file ending (BUSCO's .txt, NanoPlot's _NanoStats.txt).
# 06.report/multiqc_data/multiqc_sources.txt lists the name MultiQC gave each file: the place
# to look when a row keeps its raw folder name.
#
# The text of the chosen block is copied into multiqc_config.yaml exactly as it stands, so its
# "# ..." lines are comments of that YAML file, not of this one. r""" keeps the backslashes (in a
# plain Python string \1 would turn into an invisible control character), and the single quotes
# make YAML keep them as well. Write no braces in a pattern (such as [0-9]{2}): Snakemake would
# take them for a wildcard and stop (see the notes in rule multiqc).
#
# These lines, and the --ignore lines in rule multiqc, assume the folder layout of 00_common.smk
# (01.reads/<sample>/illumina, 01.reads/<sample>/long/raw_qc and filtered_qc,
# 02.assembly/<sample>/eval/quast, .../eval/busco, .../eval/qualimap and, in hybrid mode,
# 02.assembly/<sample>/guide/eval/qualimap). If that layout changes, change these lines with it:
# otherwise nothing fails, and MultiQC just shows the raw folder names.
_MULTIQC_RENAME_CONTIGS = r"""
  # assembly QC (QUAST)
  '^02\.assembly \| ([^|]+) \| eval \| quast \| \1$': 'QUAST assembly QC | \1'
  # genome completeness (BUSCO)
  '^02\.assembly \| ([^|]+) \| eval \| busco \| short_summary\.(?:generic|specific)\.([^.|]+_odb[0-9]+)\.\1\.txt$': 'BUSCO completeness QC (\2) | \1'
"""

_MULTIQC_RENAME_ILLUMINA = r"""
  # read QC (fastp)
  '^01\.reads \| ([^|]+) \| illumina \| .+$': 'fastp read QC | \1'
  # assembly QC (QUAST)
  '^02\.assembly \| ([^|]+) \| eval \| quast \| \1$': 'QUAST assembly QC | \1'
  # genome completeness (BUSCO)
  '^02\.assembly \| ([^|]+) \| eval \| busco \| short_summary\.(?:generic|specific)\.([^.|]+_odb[0-9]+)\.\1\.txt$': 'BUSCO completeness QC (\2) | \1'
  # mapping QC (Qualimap): the trimmed Illumina reads on the draft
  '^02\.assembly \| ([^|]+) \| eval \| qualimap \| .+$': 'mapping QC | \1'
"""

_MULTIQC_RENAME_LONGREAD = r"""
  # read QC (NanoPlot), raw and filtered long reads
  '^01\.reads \| ([^|]+) \| long \| raw_qc \| \1_NanoStats\.txt$': 'NanoPlot raw read QC | \1'
  '^01\.reads \| ([^|]+) \| long \| filtered_qc \| \1_NanoStats\.txt$': 'NanoPlot filtered read QC | \1'
  # assembly QC (QUAST)
  '^02\.assembly \| ([^|]+) \| eval \| quast \| \1$': 'QUAST assembly QC | \1'
  # genome completeness (BUSCO)
  '^02\.assembly \| ([^|]+) \| eval \| busco \| short_summary\.(?:generic|specific)\.([^.|]+_odb[0-9]+)\.\1\.txt$': 'BUSCO completeness QC (\2) | \1'
  # mapping QC (Qualimap): the filtered long reads on the draft
  '^02\.assembly \| ([^|]+) \| eval \| qualimap \| .+$': 'mapping QC | \1'
"""

_MULTIQC_RENAME_HYBRID = r"""
  # read QC (fastp), the Illumina reads
  '^01\.reads \| ([^|]+) \| illumina \| .+$': 'fastp read QC | \1'
  # read QC (NanoPlot), raw and filtered long reads
  '^01\.reads \| ([^|]+) \| long \| raw_qc \| \1_NanoStats\.txt$': 'NanoPlot raw read QC | \1'
  '^01\.reads \| ([^|]+) \| long \| filtered_qc \| \1_NanoStats\.txt$': 'NanoPlot filtered read QC | \1'
  # assembly QC (QUAST): the screened Illumina guide assembly and the final genome, one report
  '^02\.assembly \| ([^|]+) \| eval \| quast \| \1_guide$': 'QUAST guide assembly QC | \1'
  '^02\.assembly \| ([^|]+) \| eval \| quast \| \1$': 'QUAST assembly QC | \1'
  # genome completeness (BUSCO)
  '^02\.assembly \| ([^|]+) \| eval \| busco \| short_summary\.(?:generic|specific)\.([^.|]+_odb[0-9]+)\.\1\.txt$': 'BUSCO completeness QC (\2) | \1'
  # mapping QC (Qualimap): the Illumina reads on the guide, and the long reads on the draft
  '^02\.assembly \| ([^|]+) \| guide \| eval \| qualimap \| .+$': 'mapping Illumina QC | \1'
  '^02\.assembly \| ([^|]+) \| eval \| qualimap \| .+$': 'mapping QC | \1'
"""

# Pick this mode's block. _MULTIQC_RENAME_BLOCK is the text rule multiqc writes into the
# configuration file (params rename_block).
if not HAS_READS:
    # contigs mode: no read QC and no mapping QC, so the assembly QC is all there is.
    _MULTIQC_RENAME_BLOCK = _MULTIQC_RENAME_CONTIGS
elif IS_HYBRID:
    # hybrid mode: fastp, NanoPlot, both QUAST labels, BUSCO and both Qualimap reports.
    _MULTIQC_RENAME_BLOCK = _MULTIQC_RENAME_HYBRID
elif HAS_SHORT_READS:
    # illumina mode: fastp, QUAST, BUSCO and Qualimap.
    _MULTIQC_RENAME_BLOCK = _MULTIQC_RENAME_ILLUMINA
else:
    # longread mode (long reads only): NanoPlot, QUAST, BUSCO and Qualimap.
    _MULTIQC_RENAME_BLOCK = _MULTIQC_RENAME_LONGREAD

# NanoPlot writes a second NanoStats table, <sample>_NanoStats_post_filtering.txt, when it is
# asked to filter reads itself. FunFlux never asks it to, so the file should not exist, but
# MultiQC would count it as one more sample if it did, so in the modes with long reads rule
# multiqc tells MultiQC to ignore it. The other modes have no NanoPlot output, so the option
# is left out there.
if HAS_LONG_READS:
    _MULTIQC_NANOPLOT_IGNORE = "--ignore '*NanoStats_post_filtering.txt'"
else:
    _MULTIQC_NANOPLOT_IGNORE = ""


# ───────────── multiqc — one HTML report for the whole run ─────────────
# Takes in:    the QC results of every sample (see multiqc_qc_inputs): the read QC of the
#              mode's front end (fastp, NanoPlot), the Qualimap folder(s) of rule map_evaluation
#              (and rule map_guide_evaluation in hybrid mode), the QUAST and BUSCO folders of
#              rule genome_assembly_evaluation and rule completeness in shared/20_qc.smk, the
#              contig summary table of rule contig_summary and the Adapters table of rule
#              screen_adapters in shared/25_contig_summary.smk and, when the AMR screen is on,
#              the AMR table of rule amr_chroquetas in shared/50_amr.smk.
# Does:        (1) writes MultiQC's configuration file: the names of the two outputs, two options
#                  that keep the folder paths and the creation date out of the page header, the
#                  option that stops MultiQC from shortening sample names on its own, and the
#                  rename rules above;
#              (2) changes into the output folder and runs MultiQC on the QC paths, given
#                  relative to it (explained in params);
#              (3) tells MultiQC to skip the folders in which BUSCO keeps its own copies of the
#                  summaries, and in the long-read modes NanoPlot's second table (explained
#                  above the command).
# Produces:    06.report/multiqc_report.html (MULTIQC_REPORT), the configuration it was made with
#              (multiqc_config.yaml, kept) and multiqc_data/, the numbers behind the plots in
#              plain text and JSON, e.g. multiqc_general_stats.txt, one row per sample and tool.
# Consumed by: people. No rule reads it; rule all asks for the page.
#
# This rule needs one thread. A failing report stops the run, unlike QUAST: if MultiQC cannot
# produce the page, something is wrong with the QC results themselves.
#
# v1 -> v2: rewritten after BacFlux's report rule. (1) The input list comes from the capability
# flags. (2) Relative paths: v1 got short sample names because the whole workflow ran inside the
# output directory; v2 never changes directory, so this rule does it. (3) --force is new. (4) The
# configuration file is kept (v1 deleted it) and the data folder is declared as an output. (5) The
# rename rules are written out literally per mode (those of the old FunFluxL+ workflow never matched
# anything). (6) The BUSCO folders are ignored, because v2 keeps them. (7) MultiQC's own shortening
# of sample names is switched off. (8) The --ignore of NanoPlot's second table is kept from
# FunFluxL, in the long-read modes only.
rule multiqc:
    input:
        qc_inputs = _MULTIQC_INPUTS,
    output:
        report_html = MULTIQC_REPORT,
        multiqc_yaml = _MULTIQC_CONFIG,
        data_dir = directory(_MULTIQC_DATA),
    params:
        # Why relative paths: with -d, MultiQC puts the folder names in front of the sample
        # name, and a rename rule can only be written for names that look the same on every
        # computer. Absolute paths would put the whole installation path (/home/someone/project/
        # ...) in front. So input: keeps the absolute paths (the same names the QUAST and BUSCO
        # rules declare as outputs, which is how Snakemake links the jobs), the shell changes
        # into the output folder, and MultiQC gets the same folders as paths relative to it,
        # joined by spaces: 02.assembly/SAMPLE1/eval/quast and so on, so the names always begin
        # with 02.assembly. None of these paths can contain a blank or a "|", because
        # 00_common.smk refuses sample names with either (BAD_CHARS); the rename patterns rely
        # on that too.
        out_root = OUT,
        report_dir = DIR_REPORT,
        qc_rel = " ".join(_relative_to_output_root(_MULTIQC_INPUTS)),
        # Why the names of the outputs are written into the configuration: MultiQC calls its page
        # multiqc_report.html and its data folder multiqc_data by default, but the --filename
        # option changes BOTH (to <filename>.html and <filename>_data), and a personal MultiQC
        # configuration of the user can change them as well (a report title, for instance, gives
        # <title>_multiqc_report.html). Names set in a configuration file win over both, and this
        # run's file is read after the personal ones, so the outputs appear where this rule
        # declares them. (Checked in the MultiQC 1.33 code.)
        report_name = os.path.basename(MULTIQC_REPORT),
        data_name = os.path.basename(_MULTIQC_DATA),
        rename_block = _MULTIQC_RENAME_BLOCK,
        nanoplot_ignore = _MULTIQC_NANOPLOT_IGNORE,
    conda:
        "../../envs/multiqc.yaml"
    log:
        LOGS + "/multiqc.log"
    shell:
        # The command, step by step. The notes are kept here, outside the command text: a comment
        # line inside the here-document (see step 1) would be written into the configuration
        # file, one between the continued lines of the multiqc command would cut it in two, and
        # any change to the command text makes Snakemake re-run the rule.
        #
        # (1) Write MultiQC's configuration file. printf writes its first four lines:
        #     show_analysis_paths and show_analysis_time set to False (the page header then leaves
        #     out the absolute paths of the folders MultiQC read and the date the page was made),
        #     and the names of the page and the data folder (see params). cat then adds the rest
        #     from a here-document: the text between << 'EOF' and the line EOF, copied into the
        #     file as it stands (the quotes around EOF stop the shell from changing anything in it).
        #       fn_clean_sample_names: false       MultiQC does not shorten the names (see below)
        #       sample_names_replace_regex: true   the patterns are regular expressions
        #       sample_names_replace_exact: true   a pattern must match the whole name
        #     Why the shortening is off. By default MultiQC 1.33 shortens every name after -d has
        #     put the folder names in front, and before the rename rules see it: it cuts the whole
        #     name at the first place where any of 152 file endings occurs (.clean, .sorted, .trim,
        #     .fa, .contigs, ccs, ...), then trims 26 more (_trimmed, .r, _, -, ...) off its start
        #     and end. Sample names may contain dots (in contigs mode the file
        #     SAMPLE1.clean.fasta is the sample SAMPLE1.clean), so
        #     "02.assembly | SAMPLE1.clean | eval | quast | SAMPLE1.clean" would become
        #     "02.assembly | SAMPLE1": no rename rule fits it, and a sample SAMPLE1.sorted
        #     would get the same name and be merged with it. With the shortening off, a name is
        #     the folder names and the tool's own name in full, so BUSCO's ends in .txt (checked
        #     in the MultiQC 1.33 code, base_module.py, and by running MultiQC on a test tree).
        #     Those lines and the closing EOF start at the beginning of the line on purpose:
        #     Snakemake runs the command with its indentation, MultiQC needs sample_names_replace
        #     as a top-level key of the YAML file, and the shell recognises the closing EOF only at
        #     the start of a line. Indenting the YAML lines would write an invalid file; indenting
        #     the closing EOF would make the shell copy the rest of the command (cd, multiqc) into
        #     the file instead of running it, and Snakemake would then report the outputs missing.
        #     The rename lines arrive through {params.rename_block} and are written out as they
        #     are; ( ) | \ and $ mean nothing to Snakemake. Braces are the exception: Snakemake
        #     reads {...} in a params text as a wildcard, so a pattern such as [0-9]{2} would stop
        #     the workflow with a WildcardError, even with the braces doubled (tested with
        #     Snakemake 9.14.6). r""" keeps the backslashes of the command text (such as the \n of
        #     printf) as they are.
        #
        # (2) Change into the output folder, so that MultiQC can be given the QC folders as
        #     relative paths (see params).
        #
        # (3) Run MultiQC; everything it prints goes to the log.
        #     -d puts the names of the folders in front of each sample name, joined with " | ":
        #     the rename rules above are written for exactly those names.
        #     --force: without it, MultiQC does not overwrite an earlier report, data folder or
        #     plots folder in 06.report/; it writes multiqc_report_1.html beside them instead, and
        #     the rule then fails on a missing output although MultiQC finished well. Snakemake
        #     deletes the declared outputs before the job starts, but not files it does not know
        #     about, such as the multiqc_plots/ folder that a personal MultiQC configuration can
        #     switch on.
        #     The two --ignore lines: BUSCO's own result tree is kept, and inside it every lineage
        #     run has a copy of its summary (run_<lineage>/short_summary.txt, also under
        #     auto_lineage/). BUSCO also puts run_eukaryota_odb12 directly in the busco folder, as
        #     a link to auto_lineage/run_eukaryota_odb12, and MultiQC follows links, so it reads
        #     that copy twice. MultiQC would count each copy as one more sample: with BUSCO's
        #     auto-lineage, several BUSCO rows per sample instead of two. The patterns name the
        #     folders run_<lineage> and auto_lineage directly below a folder called busco, so
        #     MultiQC does not enter them. (The patterns '*/run_*/*' and '*/auto_lineage/*' do NOT
        #     work with the pinned MultiQC 1.33: MultiQC applies path patterns to folders only, so
        #     the files that lie directly inside run_<lineage> would still be read; checked on a
        #     BUSCO-shaped test tree, 4 rows instead of 2.)
        #     {params.nanoplot_ignore} is a third --ignore in the modes with long reads, for
        #     NanoPlot's optional second table (see _MULTIQC_NANOPLOT_IGNORE above), and nothing
        #     at all in the other modes.
        r"""
        mkdir -p "{params.report_dir}"

        printf "%s\n" "show_analysis_paths: False" "show_analysis_time: False" "output_fn_name: {params.report_name}" "data_dir_name: {params.data_name}" > "{output.multiqc_yaml}"

        cat >> "{output.multiqc_yaml}" << 'EOF'
fn_clean_sample_names: false
sample_names_replace_regex: true
sample_names_replace_exact: true

sample_names_replace:
{params.rename_block}
EOF

        cd "{params.out_root}"

        multiqc \
          --config "{output.multiqc_yaml}" \
          --force \
          -d \
          --ignore '*/busco/run_*' \
          --ignore '*/busco/auto_lineage' {params.nanoplot_ignore} \
          {params.qc_rel} \
          --outdir "{params.report_dir}" > "{log}" 2>&1
        """
