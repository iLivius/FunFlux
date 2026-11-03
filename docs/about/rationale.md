# Rationale

This page explains why FunFlux is built as it is. One workflow takes four kinds of input, and only
its front end depends on the mode. The page also explains why the mitochondrial genome is kept out
of the nuclear genome. Read it to understand what the `mode` key changes and what it leaves the
same.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **front end** | the steps that differ by mode, chosen by `mode`: read processing, assembly and contig filter before the screen, and polishing after it |
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **draft assembly (the draft)** | the assembly the screen examines, before any contig is removed: `contigs_filt.fasta` in `illumina` mode and in `contigs` mode with the screen, `assembly.fasta` in `longread` and `hybrid` |
    | **guide assembly** | in `hybrid` mode, the SPAdes assembly of the Illumina reads; it is screened, and only the pairs that map to it go on |
    | **coverage track** | the reads aligned back to the draft, `contaminants/{sample}_map.bam`, from which BlobTools takes each contig's coverage; in `contigs` mode, the contigs aligned to themselves |
    | **genus filter** | the setting `parameters.decontamination.mode: include` or `exclude`, under which the contig selector keeps or removes contigs by the genus BlobTools gave them |
    | **confirmed mitochondrial contig** | a contig the mitochondrion screen classes `confirmed_mitochondrial`; with the default `separate_confirmed`, it is kept whatever its genus and moved to `contigs_mitochondrial.fasta` |
    | **nuclear genome** | the delivered genome without the confirmed mitochondrial contigs, `contigs_final.fasta`; it holds every contig when nothing is split off (no screen, or `parameters.mitochondria.mode` `off` or `report`) |

## One workflow instead of four

FunFlux 2.0.0 joins four earlier workflows into one, and the `mode` key chooses which kind of
input it starts from. [Coming from v1](../getting-started/from-v1.md) lists what changed for a v1
user.

## What is shared and what depends on the mode

Only the front end differs between the modes:

| Part | Steps |
|---|---|
| Front end: differs by mode | read QC, read filtering, assembly and the contig filter before the screen, and polishing after it; in `hybrid` also the Illumina guide assembly, its own screen, and the selection of the Illumina read pairs for Polypolish and, only with `short_read_guidance: true` (off by default), for Filtlong on ONT samples |
| Shared: the same in every mode | the contamination and mitochondrion screen, the split into `contigs_final.fasta` and `contigs_mitochondrial.fasta`, QUAST and BUSCO, the contig summary, the ITS extraction and classification, repeat masking, gene prediction, functional annotation, the AMR screen and the MultiQC report |

From the screen on, every mode runs the same code, and only polishing differs between the
read-based modes. The stage folders, `02.assembly` to `06.report`, have the same names in every
mode. Where a shared step depends on the data, it follows what the run has. The coverage track of
the screen, for example, comes from the reads the run has, or from the contigs themselves in
`contigs` mode, as described on [Decontamination](../analysis/decontamination.md#coverage-track).

The long-read technology is not a mode of its own. It is a property of each sample, read from the
file name (`_ont` or `_pacbio`), so one `longread` or `hybrid` batch can hold both. See
[One technology per sample](../modes/longread.md#one-technology-per-sample).

## Why the mitochondrion is kept out of the nuclear genome

In every mode the assembly stage ends with `02.assembly/{sample}/contigs_final.fasta`, and QUAST,
BUSCO, the ITS extraction and funannotate read that file. The mitochondrial genome has to stay out
of these steps, because fungal mitochondria use a different genetic code. The mitochondria of most
fungi use code 4, in which UGA codes for tryptophan instead of stop, and those of
*Saccharomyces cerevisiae* and related yeasts use code 3. On a mitochondrial contig, the nuclear
gene predictors would call false genes.

`contigs_final.fasta` is the nuclear genome alone when the screen runs with the default
`parameters.mitochondria.mode: separate_confirmed`. Otherwise it holds every contig, a
mitochondrion included.

The steps that set the mitochondrion apart run in this order:

1. The mitochondrion screen runs on the draft assembly.
2. With the default `separate_confirmed`, the confirmed mitochondrial contigs are protected from
   the genus filter, because BLAST can assign a fungal mitochondrion to another genus.
3. Polishing comes next, in the long-read modes only. In `longread`, Medaka polishes the screened
   assembly of ONT samples unless `medaka_model: false` switches it off, and PacBio HiFi samples
   are not polished. In `hybrid`, Polypolish runs for every sample, after Medaka for ONT samples
   with Medaka on.
4. Only after that are the contigs split into `contigs_final.fasta` and
   `contigs_mitochondrial.fasta`, so the mitochondrion is polished together with the rest wherever
   a sample is polished. `contigs_all.fasta` keeps the two together.

See [Mitochondria](../analysis/mitochondria.md) for the mitochondrion screen.
