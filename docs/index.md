---
hide:
  - toc
---

# Overview

FunFlux is a [Snakemake](https://snakemake.github.io/) workflow for fungal isolate genome
analysis. It accepts Illumina reads, Oxford Nanopore (ONT) or PacBio HiFi reads, short and long
reads from the same isolate, or an existing assembly. The workflow combines assembly,
contamination and mitochondrion screening, quality control, ITS-based taxonomy, repeat masking,
gene prediction and functional annotation.

FunFlux was developed for the genomic characterisation of fungal isolates. Its results describe
the assembled genome, its predicted genes and their functional annotation, including the
secondary-metabolite gene clusters found by antiSMASH. They are predictions from the genome
sequence, to be weighed with curation and experimental evidence.

An optional screen looks for the point mutations that the FungAMR database lists with evidence
for or against resistance to antifungal drugs and fungicides. The configuration key `mode` selects
the input mode. This site documents version **v2.0.0**.

FunFlux belongs to the [BioFlux](https://github.com/stars/iLivius/lists/bioflux) family of
pipelines.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **ONT** | Oxford Nanopore Technologies, a long-read sequencing platform |
    | **HiFi** | PacBio's accurate long reads, each the consensus of several passes over one molecule (CCS) |
    | **ITS** | internal transcribed spacer — the rDNA region used as the barcode of fungi |
    | **BUSCO** | a set of single-copy genes expected in every genome of a lineage, used to estimate completeness |
    | **soft-masking** | writing repeats in lower case, so that gene finders can treat them differently without losing the sequence |
    | **AMR** | antimicrobial resistance; in FunFlux, point mutations studied for resistance against antifungal drugs and fungicides |

## The workflow at a glance

![Detailed metro map of FunFlux v2.0.0. Blue (illumina), green (longread), purple (hybrid) and gold (contigs) lines run through read QC, assembly, the Illumina guide of hybrid, the contamination and mitochondrion screen and polishing to contigs_final.fasta, and then together through assembly QC, ITS taxonomy, gene prediction, functional annotation and the contig summary to the MultiQC report. Where the lines split into parallel tracks, the steps run side by side; open rings are steps the configuration can switch on or off. A dashed red branch, the optional resistance-mutation screen, reads the final genome and the mitochondrial contigs and ends in the report.](assets/funflux_metro.png#only-light)
![Detailed metro map of FunFlux v2.0.0. Blue (illumina), green (longread), purple (hybrid) and gold (contigs) lines run through read QC, assembly, the Illumina guide of hybrid, the contamination and mitochondrion screen and polishing to contigs_final.fasta, and then together through assembly QC, ITS taxonomy, gene prediction, functional annotation and the contig summary to the MultiQC report. Where the lines split into parallel tracks, the steps run side by side; open rings are steps the configuration can switch on or off. A dashed red branch, the optional resistance-mutation screen, reads the final genome and the mitochondrial contigs and ends in the report.](assets/funflux_metro_dark.png#only-dark)

Each coloured line is one mode: follow yours from its input to the report. The four lines take
their own routes to `contigs_final.fasta`, the genome every mode delivers, and from there run on
together through the analyses. In `hybrid` the purple line has two arms: the Illumina arm builds a
guide assembly with SPAdes, screens it and keeps the Illumina pairs that map to it, and Polypolish
uses those pairs to correct the long-read assembly. In `contigs` mode the solid gold line goes
straight to `contigs_final.fasta`; the dashed gold path is the contamination and mitochondrion
screen, which the configuration can switch on. Where the lines split into parallel tracks, the
steps run side by side. Analysis groups after the final genome may run in parallel. The order
along the line is a reading order, not the order in which the steps must run. An open ring marks a
step the configuration switches on or off; the two Medaka steps apply to ONT samples only. Thin
lines are extra inputs, and the dashed red line is the optional resistance-mutation screen.
Select the map to open it larger; select it again to see it at full size, and close it with × or Esc.

The stages, their tools and when they run:

| Stage | Tools | Runs |
|---|---|---|
| Assembly | depends on the mode ([below](#the-four-modes)) | `illumina`, `longread`, `hybrid`; `contigs` mode starts from your assembly |
| Contamination screen | BLAST against NCBI core_nt, BlobTools | always; in `contigs` mode only when the screen is switched on |
| Mitochondrion screen | BLAST, with the assembler's evidence in the long-read modes | with the contamination screen |
| Assembly QC | QUAST, BUSCO | always |
| Contig summary | tidk finds the telomeres; one table per sample gives every contig of the draft a verdict ([Contig summary](analysis/contig-summary.md)) | always |
| ITS taxonomy | ITSx, SINTAX against UNITE | always |
| Repeat masking | tantan, or RepeatModeler and RepeatMasker | always |
| Gene prediction | funannotate with its own gene finders, Tiberius (on by default) and GeneMark-ES (with a licence) | always |
| Functional annotation | InterProScan, eggNOG-mapper, antiSMASH, SignalP 6 | always; SignalP 6 is optional |
| Antifungal resistance mutations | ChroQueTas against the FungAMR database | only when switched on |
| Report | MultiQC, one report for all samples | always |

From the contamination screen on, every mode runs the same steps. Only polishing differs between
the read-based modes.

## The four modes

|  | `illumina` | `longread` | `hybrid` | `contigs` |
|---|---|---|---|---|
| **Input** | Illumina paired-end reads, `{sample}_R1.<ext>` and `{sample}_R2.<ext>` | ONT reads `{sample}_ont.<ext>` or PacBio HiFi reads `{sample}_pacbio.<ext>`; one batch may hold both | Illumina pairs and ONT or PacBio HiFi reads for every sample | one assembly per sample, `{sample}.fasta`, `.fa` or `.fna` |
| **Read QC** | PhiX removal (bowtie2), then fastp | ONT adapters removed (fastplong); NanoPlot before and after Filtlong | both of the above | none |
| **Assembly** | SPAdes | hifiasm (default) or Flye | an Illumina guide assembly (SPAdes) to select the isolate's read pairs, then hifiasm or Flye on the long reads | no assembler |
| **Contig filter** | ≥ 500 bp and ≥ 2× coverage, read off SPAdes' own headers | none | the same cut-offs, on the guide assembly only | the same cut-offs when the headers are SPAdes-style; otherwise headers are cut at the first space or tab and nothing is removed |
| **Adapter screen** | the draft contigs are searched for adapter and barcode sequence (BLAST+); reported, nothing removed | the same | the same | the same |
| **Polishing** | none | Medaka, ONT samples\* | Medaka (ONT samples\*), then Polypolish with the Illumina pairs | none |
| **Delivered genome** | the screened SPAdes contigs | the screened long-read assembly, polished with Medaka for ONT samples | the screened long-read assembly, corrected with the Illumina pairs | the input contigs, screened only with `parameters.contigs.screen: true` |

\* Medaka is on by default for ONT samples, and `medaka_model: false` switches it off. PacBio HiFi
samples are never polished with Medaka. `<ext>` is `fastq`, `fq`, `fastq.gz` or `fq.gz`, the same
for all files of one kind in a batch, and `_pacbio` files must hold HiFi (CCS) reads. CLR reads are not
supported.

`hybrid` requires more computation than either single-technology mode, because it assembles and
screens the Illumina reads as well. In `hybrid`, Filtlong can also score the ONT reads against
the selected Illumina pairs, but only when `short_read_guidance` is switched on. If you used
FunFlux 1.x, Funnotator or FunFluxL, [Coming from v1](getting-started/from-v1.md) lists what
changed. See [Modes](modes/index.md) and [Choosing a mode](getting-started/choosing-a-mode.md).

## Delivered genome files

Every sample ends with up to three FASTA files in `02.assembly/{sample}/`. Their content depends
on two settings: whether the contamination and mitochondrion screen runs, and
`parameters.mitochondria.mode`. The screen always runs in `illumina`, `longread` and `hybrid`, and
in `contigs` mode only with `parameters.contigs.screen: true`.

| File | Screen on, mitochondrion split off (`separate_confirmed`, the default) | Screen on, `parameters.mitochondria.mode` set to `off` or `report` | `contigs` mode without the screen |
|---|---|---|---|
| `contigs_all.fasta` | the assembly after the contamination screen and, in `longread` and `hybrid`, after polishing; mitochondrion included | the same | not written |
| `contigs_final.fasta` | the nuclear genome: `contigs_all.fasta` without the confirmed mitochondrial contigs | every contig | the input after the contig filter, including a mitochondrial genome if the input contains one |
| `contigs_mitochondrial.fasta` | the confirmed mitochondrial contigs; an empty file when none was confirmed | an empty file | not written |

!!! note "`contigs_final.fasta` is the genome that QC, ITS and annotation read"

    QUAST, BUSCO, the ITS extraction and funannotate run on `contigs_final.fasta`. It is the
    nuclear genome alone only when the screen runs with the default `separate_confirmed`.
    Otherwise nothing is split off.

The gene models refer to a renamed copy of that genome. `funannotate clean` and
`funannotate sort` turn `contigs_final.fasta` into `04.annotation/{sample}/contigs_sort.fasta`,
whose contigs are named `scaffold_1`, `scaffold_2`, ... from the longest down, and the gene
models use these names. QUAST, BUSCO and the ITS files keep the names of `contigs_final.fasta`.
`contigs_removed_by_clean.txt`, next to `contigs_sort.fasta`, lists the contigs
`funannotate clean` removed.

## One command, one configuration, one report

All four modes are launched the same way, from the repository root:

```bash
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml
```

`--sdm conda` builds the conda environment of each rule on first use. Besides Snakemake, you
obtain InterProScan yourself, and the optional GeneMark-ES and SignalP 6 if you use them
([Licensed and external tools](getting-started/external-tools.md)). `--cores` sets the CPU budget
of a local run and is used without `--jobs`. See [Running FunFlux](reference/running.md).

Pass the configuration file explicitly with `--configfile`. FunFlux ships a template,
`config/config.yaml`, but does not load it on its own: copy it, then set the mode, the input and
output paths, the database locations and the annotation table. A wrong mode, a missing input
directory, a bad sample name or a missing annotation row stops the dry run before any job starts.
See [Configuration](reference/configuration.md) and [Quick start](getting-started/quick-start.md).

Results are written to numbered folders that are the same in every mode: `00.resources` for
material shared by the samples of a run, `01.reads` (not in `contigs` mode), `02.assembly`,
`03.taxonomy`, `04.annotation`, `05.amr` (only with the resistance-mutation screen on),
`06.report` and `logs`. The funannotate database stays in `directories.funannotate_db`, outside
the output folder. Every run ends with a [MultiQC](https://multiqc.info/) report,
`06.report/multiqc_report.html`, which collects the QC of all samples. The annotation itself is
in `04.annotation/{sample}/funannotate/annotate_results/`. See
[Output files](reference/output.md).

!!! note "What has to be prepared in advance"

    FunFlux installs its own software, but three resources must be on disk before the first run:
    NCBI core_nt with its taxonomy files (not read in `contigs` mode without the screen), the
    eggNOG-mapper data, and an InterProScan installation. GeneMark-ES and SignalP 6 are optional
    and need a licence. Everything else is downloaded on first use, or can be pointed at a copy
    you already hold. See [Reference databases](getting-started/databases.md) and
    [Licensed and external tools](getting-started/external-tools.md).

## Where the resistance-mutation screen fits

Stage `05.amr` is disabled by default (`parameters.amr.chroquetas: off`). When enabled,
ChroQueTas searches each genome for the point mutations that the FungAMR database lists with
evidence for or against resistance to antifungal drugs and fungicides, and FunFlux summarises for each sample the
known changes found and the drugs FungAMR links them to. A sample of a species that FungAMR does
not cover is reported as not screened.

These calls describe known mutations in the genes and species that FungAMR covers. They do not
establish a resistance phenotype, and the absence of a call is not evidence of susceptibility.
Interpret them together with the evidence levels FungAMR gives each mutation. See
[Antifungal resistance mutations](analysis/amr.md).

## Where to go next

| Page | Read it for |
|---|---|
| [Choosing a mode](getting-started/choosing-a-mode.md) | which mode a dataset belongs in, and how to name the input files |
| [Installation](getting-started/installation.md) | Snakemake (9.14.6 or a later 9.x; tested with 9.14.6), the repository and the conda environments |
| [Reference databases](getting-started/databases.md) | what to download before the first run, and what FunFlux fetches |
| [Licensed and external tools](getting-started/external-tools.md) | InterProScan, which every run needs, and the optional GeneMark-ES and SignalP 6 |
| [Quick start](getting-started/quick-start.md) | the commands from a filled-in configuration to a finished run |
| [Configuration](reference/configuration.md) | every configuration key, the modes that read it, its default and allowed values |
| [Running FunFlux](reference/running.md) | the dry run, CPU cores and threads, resuming an interrupted run, running on a cluster |
| [Troubleshooting](troubleshooting.md) | what a `[FunFlux]` message or a failed step means, and what to do |
| [Validation on public genomes](about/validation.md) | what each read-based mode produced on *Botrytis cinerea* B05.10, what PacBio HiFi reads of *Cordyceps cateniannulata* GXU-8616 gave, and how to repeat the runs |

FunFlux is released under the MIT License; the tools and databases it calls carry their own
terms, listed on [Licensing](about/licensing.md).
