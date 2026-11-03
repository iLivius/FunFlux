# Report

The MultiQC report, `06.report/multiqc_report.html`, gathers the QC results of all samples of a
run into one page. MultiQC 1.33 makes it in rule `multiqc`, once per run, when the QC results,
the contig summary and, with the AMR screen on, the AMR table of every sample are ready. The
contig summary waits for the finished annotation, so the report is one of the last jobs of a run. This page describes what MultiQC
reads, how the rows are named, which numbers the report shows and which files rule `multiqc`
writes. Read it to find a sample's row in the report or the file behind a number.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **General Statistics table** | the MultiQC summary table; for an assembly it shows the N50, the length and BUSCO's complete single-copy and missing counts |
    | **draft assembly (the draft)** | the assembly before the contamination and mitochondrion screen removes any contig, such as `contigs_filt.fasta` or `assembly.fasta`; the contig summary describes it |
    | **guide assembly** | the SPAdes assembly of the Illumina reads; it is screened, and only the pairs that map to it go on |
    | **contig summary** | the table `02.assembly/{sample}/eval/contig_summary.tsv`: one row per contig of the draft, with its evidence (telomeres, depth, rDNA units, genes) and a verdict |
    | **verdict (contig summary)** | the label the contig summary gives each draft contig, one of ten such as `chromosome`, `chromosome_arm`, `rDNA_array` or `mitochondrion`, with its reasons |
    | **AMR screen** | the optional step, off by default (`parameters.amr.chroquetas`), in which ChroQueTas looks in each genome for the changes FungAMR lists for the species; results in `05.amr/` |
    | **target protein** | a protein whose FungAMR positions ChroQueTas checks, such as `Beta-tubulin` or `Cyp51` |
    | **verdict (AMR screen)** | the reading of FungAMR's evidence for a change or a target: `resistance`, `conflicting`, `no resistance`, `unknown`, `wild type`, or `not compared` for a position |

## What feeds it

QUAST, BUSCO and the contig summary feed the report in every mode. The read and mapping QC depend
on the mode, and the AMR table appears only with the AMR screen on:

| Tool | What MultiQC reads | Modes |
|---|---|---|
| fastp 1.0.1 | `01.reads/{sample}/illumina/{sample}_fastp.json` | `illumina`, `hybrid` |
| NanoPlot 1.46.2 | `01.reads/{sample}/long/raw_qc/`, `filtered_qc/` (`{sample}_NanoStats.txt`) | `longread`, `hybrid` |
| Qualimap 2.3 | `02.assembly/{sample}/eval/qualimap/`; in `hybrid` also `guide/eval/qualimap/` | `illumina`, `longread`, `hybrid` |
| QUAST 5.3.0 | `02.assembly/{sample}/eval/quast/`; in `hybrid` the report holds the guide assembly and the final genome, two rows | all |
| BUSCO 6.0.0 | `02.assembly/{sample}/eval/busco/` | all |
| FunFlux contig summary | `02.assembly/{sample}/eval/contig_summary_mqc.tsv`: one row per contig of the draft assembly, with its verdict | all |
| FunFlux AMR screen (ChroQueTas) | `05.amr/{sample}/amr_summary_mqc.tsv`: one row per target protein of the species, with the FungAMR mutations found and a verdict on resistance | all, with `parameters.amr.chroquetas: on` |

The ITS taxonomy and the annotation are not in the report. Their results are in `03.taxonomy/`
and `04.annotation/`.

QUAST, BUSCO and Qualimap themselves are described under [Assembly QC](assembly-qc.md). The pages
of the individual modes, reached from [Modes](../modes/index.md), describe fastp and NanoPlot. The
contig summary is described under [Contig summary](contig-summary.md), and the AMR table under
[Antifungal resistance mutations](amr.md).

!!! warning "A failed QUAST, Qualimap or NanoPlot run does not stop the workflow"

    The row of the sample concerned is then missing from the report. A failed MultiQC run, on the
    other hand, stops the workflow.

## Sample names

The report names each row after what it holds and the sample it belongs to. These are the names,
the files they come from and the modes that have them:

| Row | From | Modes |
|---|---|---|
| `fastp read QC | <sample>` | `01.reads/{sample}/illumina/{sample}_fastp.json` | `illumina`, `hybrid` |
| `NanoPlot raw read QC | <sample>`, `NanoPlot filtered read QC | <sample>` | `01.reads/{sample}/long/raw_qc/{sample}_NanoStats.txt`, `filtered_qc/{sample}_NanoStats.txt` | `longread`, `hybrid` |
| `mapping QC | <sample>` | `02.assembly/{sample}/eval/qualimap/`: the sample's reads on the draft | `illumina`, `longread`, `hybrid` |
| `mapping Illumina QC | <sample>` | `02.assembly/{sample}/guide/eval/qualimap/`: the Illumina pairs on the guide | `hybrid` |
| `QUAST guide assembly QC | <sample>` | `eval/quast/report.tsv`, the row labelled `{sample}_guide` | `hybrid` |
| `QUAST assembly QC | <sample>` | `eval/quast/report.tsv` | all |
| `BUSCO completeness QC (<lineage>) | <sample>` | `eval/busco/short_summary.generic.<lineage>.{sample}.txt` or `short_summary.specific.<lineage>.{sample}.txt` | all |
| `<sample>: <contig>` | `eval/contig_summary_mqc.tsv`, one row per contig of the draft assembly | all |
| `<sample>: <target>` | `05.amr/{sample}/amr_summary_mqc.tsv`, one row per target protein | all, AMR screen on |

These names were checked in finished runs of the three read-based modes on one public genome, and
of `longread` and `hybrid` with PacBio HiFi reads on a second. The runs are described on
[Validation on public genomes](../about/validation.md). A row whose name does not come out as
above shows the raw folder path instead. Its numbers are not affected.

With `parameters.busco.lineage: auto`, BUSCO writes a generic summary, for the broad lineage it
tried first, and a specific one, for the lineage it ended with. Every sample then gets two BUSCO
rows. With a named lineage, it gets one.

Sample names are shown in full. Otherwise MultiQC would cut a name at the first file ending it
knows, such as `.clean` or `.sorted`. A sample named `SAMPLE1.clean` would then get a wrong name,
and two such samples would share one row. The names do not depend on where the output folder lies.

## What the report shows

For QUAST and BUSCO, the General Statistics table of the report shows by default the N50, the
assembly length and the complete single-copy and missing BUSCO counts. Its numbers are in
`06.report/multiqc_data/multiqc_general_stats.txt`. That file has one line for each row name
listed above except the contig-summary and AMR rows, for example `QUAST assembly QC | SAMPLE1`.
Each line fills only the columns of its own tool. The QUAST and BUSCO columns are these:

| Column | Unit | In the General Statistics table of the report |
|---|---|---|
| `quast-N50` | kb | shown |
| `quast-Total_length` | Mb | shown, as the assembly length |
| `busco-complete` | BUSCO genes | hidden by default |
| `busco-complete_single_copy` | BUSCO genes | shown |
| `busco-complete_duplicated` | BUSCO genes | hidden by default |
| `busco-fragmented` | BUSCO genes | hidden by default |
| `busco-missing` | BUSCO genes | shown |

The complete, duplicated and fragmented counts can be switched on in the column settings of the
General Statistics table.

The sections of the report for QUAST, BUSCO and the two FunFlux tables hold the following:

| Section of the report | What it holds |
|---|---|
| QUAST | the full QUAST table and a plot of contig counts by length |
| BUSCO | one bar chart per dataset |
| **Contig summary** | `contig_summary_mqc.tsv` of every sample as one table, one row per contig of the draft assembly with all 20 columns of the file; the table's column settings show or hide columns, and `multiqc_data/multiqc_contig_summary.txt` holds the numbers. It is not part of the General Statistics table |
| **AMR mutations (ChroQueTas)**, with the AMR screen on | `amr_summary_mqc.tsv` of every sample as one table, one row per target protein of the species: its status, the positions compared, the verdict, the FungAMR mutations and new mutations found, and the drugs with evidence of resistance and of no resistance, with their levels ([The tables](amr.md#tables)) |

!!! note "MultiQC shows a table with more than 500 rows as a violin plot"

    The limit is MultiQC's `max_table_rows` setting. The contig summary has one row per contig of
    every sample, so for a batch of fragmented assemblies, read
    `02.assembly/{sample}/eval/contig_summary.tsv` itself, as described under
    [Contig summary](contig-summary.md#multiqc).

## Files

Rule `multiqc` writes the report, the configuration it was made with, the numbers behind it and a
log:

| Path | Content |
|---|---|
| `06.report/multiqc_report.html` | the report |
| `06.report/multiqc_config.yaml` | the configuration the report was made with: no absolute paths and no creation date and time at the top of the page (`show_analysis_paths: False`, `show_analysis_time: False`), sample names not shortened (`fn_clean_sample_names: false`), and the rules that turn folder paths into the sample names above (`sample_names_replace`) |
| `06.report/multiqc_data/` | the numbers behind the page as text and JSON: `multiqc_general_stats.txt` (one line per row name of the report), `multiqc_quast.txt`, `multiqc_busco.txt`, `multiqc_contig_summary.txt`, `multiqc_data.json`, MultiQC's own `multiqc.log`, and more |
| `logs/multiqc.log` | what MultiQC printed |

The report covers the samples of the current run only and is made again when they change.
`multiqc_data/` is then replaced as a whole.

To combine several runs in one report, run MultiQC yourself with
`--ignore '*/busco/run_*' --ignore '*/busco/auto_lineage'`. The kept BUSCO folders hold extra
copies of each BUSCO summary, and without these two options each copy shows as one more BUSCO row.

A personal MultiQC configuration cannot move the report or `multiqc_data/`, because
`multiqc_config.yaml` sets their names and is read last.
