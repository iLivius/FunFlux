# FunFlux

```
__________             _______________
___  ____/___  ___________  ____/__  /___  _____  __
__  /_   _  / / /_  __ \_  /_   __  /_  / / /_  |/_/
_  __/   / /_/ /_  / / /  __/   _  / / /_/ /__>  <
/_/      \__,_/ /_/ /_//_/      /_/  \__,_/ /_/|_|

FunFlux v2.1.0
```

**Fungal genome assembly and annotation from Illumina short reads, ONT or PacBio HiFi long reads,
a hybrid of short and long reads, or an assembly you already have.**

<!-- Badges. The last one, "ci", shows the result of the automatic checks in .github/workflows/ci.yml (documentation build, unit tests, Snakemake dry runs) on the newest commit of main. GitHub updates it after every push; it does not run the workflow on real data. -->
[![Snakemake](https://img.shields.io/badge/snakemake-tested%209.14.6-brightgreen.svg)](https://snakemake.github.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13612159.svg)](https://doi.org/10.5281/zenodo.13612159)
[![CI](https://github.com/iLivius/FunFlux/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/iLivius/FunFlux/actions/workflows/ci.yml)

<!-- Workflow overview as a metro map, one coloured line per mode. GitHub shows the dark version to readers in dark mode and the light one otherwise. The detailed map is on the documentation home page (docs/index.md). -->
<p align="center"><picture><source media="(prefers-color-scheme: dark)" srcset="docs/assets/funflux_overview_dark.png"><img src="docs/assets/funflux_overview.png" alt="FunFlux v2.1.0 workflow overview as a metro map. One coloured line per mode (Illumina, ONT or PacBio HiFi long reads, hybrid, contigs) runs from the reads, or from an assembly, to one final genome file, contigs_final.fasta; ONT reads are trimmed of the sequencing adapter during read QC, and the hybrid line also builds the Illumina guide assembly. From there the four lines run on together through assembly QC, ITS taxonomy, gene prediction, functional annotation and the contig summary, which also screens the draft for adapter sequence, to a MultiQC report. A dashed red branch is the optional resistance-mutation screen; open rings are steps the configuration can switch on or off."></picture></p>

The detailed map, step by step, is on the [documentation home page](https://iLivius.github.io/FunFlux/).

## 📖 Documentation

**Full documentation: [iLivius.github.io/FunFlux](https://iLivius.github.io/FunFlux/)**

This README is deliberately short: it covers what FunFlux is and how to get a first run going.
Everything else — choosing a mode, installation, the reference databases and licensed tools,
every configuration key, the output files, running on a cluster, troubleshooting, and the
[validation on public genomes](https://iLivius.github.io/FunFlux/about/validation/) — lives on
the documentation site.

## Synopsis

FunFlux is a [Snakemake](https://snakemake.github.io/) workflow that takes a fungal isolate from
raw reads — or from an assembly somebody else made — through assembly, contamination and
mitochondrion screening, quality control, ITS taxonomy, repeat masking, gene prediction and
functional annotation, with one MultiQC report per run. One key in the configuration, `mode`,
decides which of four front ends runs:

```
               illumina        longread           hybrid                contigs
01.reads       PhiX removal    fastplong (ONT)    PhiX removal, fastp   -
               fastp           Filtlong           fastplong (ONT)
                               NanoPlot           Filtlong, NanoPlot
02.assembly    SPAdes          hifiasm or Flye    SPAdes guide,         contig filter
                                                  hifiasm or Flye
               ── contamination screen (BLAST+, BlobTools), mitochondrion screen * ──
               ── adapter screen of the draft contigs (BLAST+) ──
                               Medaka (ONT)       Medaka (ONT),
                                                  Polypolish
               ── 02.assembly/{sample}/contigs_final.fasta: the genome all later steps read ──
02.assembly    QUAST, BUSCO
03.taxonomy    ITS: 5.8S windows (nhmmer), ITSx, SINTAX against UNITE
04.annotation  funannotate clean and sort; tantan or RepeatModeler + RepeatMasker;
               Tiberius ‡ and funannotate predict (GeneMark-ES †); InterProScan,
               eggNOG-mapper, antiSMASH, SignalP 6 †; funannotate annotate
05.amr         resistance mutations: ChroQueTas against FungAMR †
06.report      MultiQC

* contigs mode: only with parameters.contigs.screen: true    † optional    ‡ on by default
```

FunFlux was developed for the genomic characterisation of fungal isolates. Its results describe
the assembled genome, its predicted genes and their functional annotation, including the
secondary-metabolite gene clusters found by antiSMASH. They are predictions from the genome
sequence, to be weighed with curation and experimental evidence.

Each mode produces `02.assembly/{sample}/contigs_final.fasta`, the genome that every later
analysis reads. When the mitochondrion screen runs with its default settings, the confirmed
mitochondrial contigs are split off into `contigs_mitochondrial.fasta`, and `contigs_final.fasta`
holds the nuclear genome alone. The screen always runs in `illumina`, `longread` and `hybrid`, and
in `contigs` mode only when it is switched on. Common processing improves consistency across datasets, while the
read type, the sequencing depth and the fragmentation of the assembly still affect the results.

**The resistance-mutation screen** (`05.amr`, disabled by default) searches each genome with
ChroQueTas for the point mutations that the FungAMR database lists with evidence for or against
resistance to antifungal drugs and fungicides. It reports known mutations in the genes and species the
database covers. It does not establish a resistance phenotype, and the absence of a call is not
evidence of susceptibility. See
[Antifungal resistance mutations](https://iLivius.github.io/FunFlux/analysis/amr/).

v2.0.0 brings together the earlier FunFlux 1.x, Funnotator and FunFluxL workflows as modes of one
workflow. The earlier versions retain their DOIs and remain citable for work performed with them.
A v1 configuration does not run in v2. If you are coming from one of them, start at
[Coming from v1](https://iLivius.github.io/FunFlux/getting-started/from-v1/).

`FunFlux` belongs to the [BioFlux](https://github.com/stars/iLivius/lists/bioflux) family of
pipelines.

## Quick start

```bash
# 1. Clone
git clone https://github.com/iLivius/FunFlux.git
cd FunFlux

# 2. Launcher environment: Snakemake 9.14.6 or a later 9.x release (FunFlux v2 was tested
#    with 9.14.6 and stops on anything older; 10 and later have not been tested)
conda create -n snakemake -c conda-forge -c bioconda "snakemake>=9.14.6,<10"
conda activate snakemake

# 3. Copy the template and edit the copy: mode, the input directory of that mode,
#    output_dir, the database paths, files.iprscan and files.annotation_params.
#    config/config_custom.yaml is git-ignored.
cp config/config.yaml config/config_custom.yaml

# 4. Dry run, then run. --cores takes the value of resources.threads (24 in the template).
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml -n
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml
```

Run both commands from the repository root: `workflow/Snakefile` is found there automatically.
Besides Snakemake, three tools are not installed by conda: InterProScan, which every run needs, and
the optional GeneMark-ES and SignalP 6, which need a licence from their authors. You obtain them
yourself ([Licensed and external tools](https://iLivius.github.io/FunFlux/getting-started/external-tools/)).
`--sdm conda` builds each rule's own conda environment on first use, so every other analysis tool
is provisioned for you. The databases are listed below.

For a local run, `--cores N` sets the CPU budget and should take the value of
`resources.threads`. Do not add `--jobs`/`-j`: for a local run it is another name for `--cores`, and giving both can let
the jobs use more CPUs than the machine has. See
[Running FunFlux](https://iLivius.github.io/FunFlux/reference/running/).

Always dry-run first. The configuration is checked while the plan is being built, so `-n` takes
seconds and stops on a wrong mode, a missing input directory, a bad sample name or a missing
annotation row before the first job starts.

Every sample needs a row in the tab-separated table named by `files.annotation_params`.
`config/annotation_parameters.tsv` is an example:

| Column | Content |
|---|---|
| 1 | sample name |
| 2 | species |
| 3 | a protein FASTA of related fungi (gene evidence) |
| 4, optional | a stored AUGUSTUS model to use instead of training AUGUSTUS on the genome |

**Three resources must be on disk before the first run**, and the two optional licensed tools
if you use them:

| Config key | What | Version | Size |
|---|---|---|---|
| `directories.blast_db` | [NCBI `core_nt`](https://ftp.ncbi.nlm.nih.gov/blast/db/) with `taxdb.btd`/`taxdb.bti` and the taxdump files `nodes.dmp`/`names.dmp` in the same folder; not read in `contigs` mode without the screen | — | 264 GB (build of 11 July 2026, index files included; the taxonomy files add 0.8 GB) |
| `directories.eggnog_db` | [eggNOG-mapper data](https://github.com/eggnogdb/eggnog-mapper/wiki) (`eggnog.db`, `eggnog_proteins.dmnd`) | eggNOG 5.0.2 | 48 GB |
| `files.iprscan` | [InterProScan](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html), the path of `interproscan.sh` | 5.78-109.0 | 50 GB after setup |
| `directories.genemark_dir` | [GeneMark-ES](https://genemark.bme.gatech.edu/license_download.cgi) with its licence key in `~/.gm_key`; optional | 4 | — |
| `files.signalp6_archive` | [SignalP 6](https://services.healthtech.dtu.dk/services/SignalP-6.0/) archive from DTU; optional, switched on with `parameters.signalp6.mode` | 6.0 | — |

The eggNOG version is fixed by eggnog-mapper 2.1.13, whose `download_eggnog_data.py` fetches
eggNOG 5.0.2, so data made for eggnog-mapper 3 (eggNOG 7) does not work with it.

Everything else — the funannotate database snapshot (Zenodo record 23265189, checked against its
SHA-256 before unpacking), UNITE, the antiSMASH databases, the Tiberius model, the telomere repeat
database, the Illumina adapter records of NCBI UniVec, the PhiX genome, the BUSCO lineage data,
the Medaka model weights and, with the resistance-mutation screen on, ChroQueTas with its FungAMR
database — is downloaded on the first run that needs it, most of it once per output folder, or all
at once with
`snakemake download_all`. `directories.unite_db` and `directories.antismash_db` take copies you
already hold. Download commands, disk sizes and every optional key are on
[Reference databases](https://iLivius.github.io/FunFlux/getting-started/databases/). Running on a
cluster whose compute nodes have no internet is described on
[On a cluster](https://iLivius.github.io/FunFlux/reference/running/#on-a-cluster).

## Output

Everything lands under `directories.output_dir`, in numbered folders that are the same in every
mode: `00.resources` for material shared by the samples of a run, then `01.reads` (not in
`contigs` mode) through `06.report`, and `logs`. The delivered genome is always
`02.assembly/{sample}/contigs_final.fasta`. The annotation — GenBank, GFF3, proteins and the
functional annotation as a table with one row per transcript — is in
`04.annotation/{sample}/funannotate/annotate_results/`, and the QC of all samples in
`06.report/multiqc_report.html`.

funannotate renames the contigs `scaffold_1`, `scaffold_2` and so on, while QUAST, BUSCO and the
ITS files keep the original names. The funannotate database is unpacked into `directories.funannotate_db`,
outside the output folder, and Snakemake keeps its records and the conda environments in
`.snakemake/` in the directory you launch from. One fungal genome in `contigs` mode takes several
hours and leaves in the order of 10 GB. What each stage holds, file by file, is on
[Output files](https://iLivius.github.io/FunFlux/reference/output/).

## Citation

> Antonielli, L., & Brader, G. (2026). *FunFlux: Integrated workflow for fungal genome assembly
> and annotation* (v2.1.0). Zenodo. <https://doi.org/10.5281/zenodo.13612159>

That DOI is the concept DOI and always resolves to the newest release; each tagged release also
gets its own version DOI, which is the one to cite when the exact code matters.

Most of the results FunFlux reports are produced by other people's tools and databases, so
**cite those too**. The list is on
[Citation and references](https://iLivius.github.io/FunFlux/about/citation/).

## Acknowledgements

This work was supported by the [RATION project](https://cordis.europa.eu/project/id/101084163)
(Risk AssessmenT InnOvatioN for low-risk pesticides), funded by the European Union's Horizon
Europe research and innovation programme under grant agreement No. 101084163.

Much of v2.0.0 — merging FunFlux 1.x, Funnotator and FunFluxL into one workflow, adding the hybrid
mode and writing the documentation — was done with [Claude Code](https://claude.com/claude-code).
Anthropic provided six months of Claude Max through their Open Source programme, and the scale of
the v2 rewrite would not have been realistic without it.

## Licence

FunFlux is released under the [MIT License](LICENSE). The repository ships code, download
addresses and two small files of adapter sequences (from Oxford Nanopore and NCBI UniVec), and no
databases: every database is either one you already hold or one the workflow fetches from its
publisher, under that publisher's terms. The tools and databases
FunFlux runs carry their own licences, listed with their conditions on
[Licensing](https://iLivius.github.io/FunFlux/about/licensing/). GeneMark-ES, SignalP 6, Medaka,
ChroQueTas with its FungAMR database (the resistance-mutation screen, off by default), and the
PROSITE and SMART databases that InterProScan searches in every run restrict commercial use.
