# Modes

The configuration key `mode` selects the front end, the steps that turn the input into an assembly.
The launch command does not change with the mode. This page sets the four modes side by side and
compares their steps, the contamination and mitochondrion screen as each mode runs it, and the
files each mode writes. Read it when you choose a mode, or to look up what a mode does.

From the screen on, every mode runs the same steps, and only polishing differs between the
read-based modes. Every mode produces `02.assembly/{sample}/contigs_final.fasta`. This file is the
nuclear genome alone when the screen runs with the default
`parameters.mitochondria.mode: separate_confirmed`. Otherwise nothing is split off, and any
mitochondrial contig stays in it.

From that file on, all four modes run QUAST, BUSCO, ITS extraction and classification, repeat
masking, gene prediction, functional annotation and the MultiQC report. Two of these steps take
extra inputs in some modes. In `hybrid` mode QUAST also reports the Illumina guide assembly, and
MultiQC adds the read-QC and mapping tables of the modes that have reads.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **front end** | the steps that differ by mode, chosen by `mode`: read processing, assembly and contig filter before the shared screen, and polishing after it |
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **draft assembly (the draft)** | the assembly the screen examines, before any contig is removed |
    | **contig selector** | rule `select_contigs`, the last step of the screen: it keeps or removes each contig by genus, as `parameters.decontamination` says; the default, `off`, keeps every contig |
    | **genus filter** | the setting `parameters.decontamination.mode: include` or `exclude`, under which the contig selector keeps or removes contigs by the genus BlobTools gave them |
    | **guide assembly** | in `hybrid` mode, the SPAdes assembly of the Illumina reads; it is screened, and only the pairs that map to it go on |
    | **assembler evidence** | what hifiasm or Flye report for each contig, above all read depth and circularity, written to `assembly_evidence.tsv` for the mitochondrion screen; `longread` and `hybrid` only |

## The four modes side by side

Three modes build an assembly from reads, and `contigs` mode starts from one. The table compares
the four modes step by step, and each column heading links to the page of its mode.

| | [`illumina`](illumina.md) | [`longread`](longread.md) | [`hybrid`](hybrid.md) | [`contigs`](contigs.md) |
|---|---|---|---|---|
| Status | run on one public genome ([Validation](../about/validation.md)) | run on ONT reads of one public genome and PacBio HiFi reads of another, with hifiasm and with Flye | run on ONT and Illumina reads of one public genome and PacBio HiFi and short reads of another, with hifiasm; Flye in dry runs only | run on the assembly of one public genome ([Validation](../about/validation.md#gene-finders)) |
| Input | Illumina paired-end reads | ONT or PacBio HiFi reads, technology set per sample | Illumina pairs plus ONT or PacBio HiFi reads of the same isolate | an assembled genome (FASTA) |
| Input folder and file names | `input.illumina_dir`: `{sample}_R1.<ext>`, `{sample}_R2.<ext>` | `input.longread_dir`: `{sample}_ont.<ext>` or `{sample}_pacbio.<ext>` | both folders (may be the same one) and both sets of names | `input.contigs_dir`: `{sample}.fasta`, `.fa` or `.fna` |
| Other required keys | `directories.blast_db`, `links.phix_link` | `directories.blast_db` | `directories.blast_db`, `links.phix_link` | `directories.blast_db`, only with the screen |
| Mode settings | `parameters.spades_kmers` | `parameters.longread` | `parameters.hybrid`, `parameters.spades_kmers` | `parameters.contigs` |
| Read processing | PhiX removal (bowtie2), fastp | Filtlong; NanoPlot before and after | PhiX removal, fastp; Filtlong (ONT reads scored against the Illumina pairs only if you switch it on); NanoPlot | none |
| Assembly | SPAdes `--isolate` | hifiasm (default) or Flye | SPAdes guide assembly, then hifiasm or Flye | none |
| Contig filter | at least 500 bp and 2x coverage | none | guide assembly only: at least 500 bp and 2x | SPAdes headers: at least 500 bp and 2x; other headers: none |
| Polishing | none | Medaka, ONT samples | Medaka (ONT samples), then Polypolish with the Illumina pairs | none |
| Mapping QC | Qualimap | Qualimap | Qualimap, guide and long-read assembly | none |
| `contigs_all.fasta` is | the screened contigs | the Medaka consensus (ONT) or the screened assembly (PacBio HiFi, or Medaka off) | the Polypolish output | the screened FASTA; written only with the screen |
| Not run | long-read steps | PhiX removal, fastp, SPAdes, Polypolish; Medaka for PacBio HiFi | Medaka and Filtlong guidance for PacBio HiFi | read processing, assembly, polishing, mapping QC; the screen unless `parameters.contigs.screen: true` |

`contigs` mode skips every read and assembly step and runs the screen only with
`parameters.contigs.screen: true`. `hybrid` mode combines the other two read-based modes. Its
Illumina steps build a guide assembly and select reads, its long-read steps build the genome that
is delivered, and Polypolish corrects that genome with the Illumina pairs.

In `longread` and `hybrid` mode the read technology is taken from each file name, so one batch can
hold ONT and PacBio HiFi samples. The flags that differ between the two technologies are listed on
[Long reads](longread.md#one-technology-per-sample).

## The contamination and mitochondrion screen

The screen checks every contig of the draft assembly. BLAST compares each contig with NCBI
core_nt, and BlobTools gives it a genus. The mitochondrion screen looks for the mitochondrial
genome among the contigs, and the contig selector keeps or removes each contig and records the
decision. With the default `parameters.mitochondria.mode: separate_confirmed`, the contig selector
keeps a confirmed mitochondrial contig whatever the genus filter says. Polishing comes after the
screen, on the contigs that were kept.

The screen is the same in every mode. Three things change with the mode: the draft it examines,
the reads behind the BlobTools coverage and the evidence for the mitochondrion.

| | `illumina` | `longread` | `hybrid` | `contigs` with the screen |
|---|---|---|---|---|
| Draft screened | `contigs_filt.fasta` | `assembly.fasta` | the guide's `contigs_filt.fasta`, then `assembly.fasta` | `contigs_filt.fasta` |
| Reads behind the BlobTools coverage | trimmed Illumina pairs, bowtie2 | filtered long reads, minimap2 | guide: Illumina pairs, bowtie2; long-read assembly: filtered long reads, minimap2 | none; the contigs are mapped to themselves and the coverage carries no information |
| Evidence for the mitochondrion | BLAST | BLAST plus the assembler's read depth and circularity | guide: BLAST; long-read assembly: BLAST plus assembler evidence | BLAST |

A long-read assembler adds read depth and circularity to BLAST when the mitochondrion is judged.
`hybrid` mode runs the screen twice, on the guide as in `illumina` and on the long-read assembly
as in `longread`. In `contigs` mode there are no reads, so the BlobTools coverage carries no
information and BLAST alone decides.

!!! note "With the default settings, the screen removes no contig"

    The screen is on in every read-based mode and off in `contigs` mode unless
    `parameters.contigs.screen: true`. With the default `parameters.decontamination.mode: off`, it
    removes no contig but records a decision for each one and splits off the confirmed
    mitochondrion. See [Decontamination](../analysis/decontamination.md) and
    [Mitochondria](../analysis/mitochondria.md).

## What `01.reads` and `02.assembly` contain

The table shows which folders and files each mode writes to `01.reads` and `02.assembly`.
`contigs_final.fasta`, QUAST, BUSCO, the telomere files and the contig summary are there in every
mode. The other rows depend on each mode's reads, assembler, screen and polisher.
Rows without a stage prefix are under `02.assembly/{sample}/`, and `01.reads/` is not created in
`contigs` mode.

| Path | `illumina` | `longread` | `hybrid` | `contigs` |
|---|---|---|---|---|
| `01.reads/{sample}/illumina/`: fastp reports (kept), PhiX-free and trimmed reads (temporary) | yes | | yes | |
| `01.reads/{sample}/long/`: `{sample}_filtered.fastq`, NanoPlot `raw_qc/` and `filtered_qc/` | | yes | yes | |
| `01.reads/{sample}/selected/`: Illumina pairs that map to the screened guide assembly | | | yes | |
| `02.assembly/{sample}/spades/` | yes | | | |
| `02.assembly/{sample}/guide/`: the Illumina guide assembly, its screen and Qualimap report | | | yes | |
| `hifiasm/` or `flye/`, `assembly.fasta`, `assembly_evidence.tsv`, `assembly_provenance.tsv` | | yes | yes | |
| `contigs_filt.fasta` | yes | | | with the screen |
| `{sample}_medaka_model.txt`, `medaka/consensus.fasta` | | ONT samples with Medaka | ONT samples with Medaka | |
| `polypolish/polished.fasta` | | | yes | |
| `contaminants/`, `mitochondria/` | yes | yes | yes | with the screen |
| `contigs_all.fasta`, `contigs_mitochondrial.fasta` | yes | yes | yes | with the screen |
| `contigs_final.fasta` | yes | yes | yes | yes |
| `eval/quast/`, `eval/busco/` | yes | yes | yes | yes |
| `eval/telomere_candidates.tsv`, `telomere_unit.txt`, `telomere_ends.tsv`: the telomere repeat unit of the assembly and the arrays at the contig ends | yes | yes | yes | yes |
| `eval/contig_summary.tsv`, `_mqc.tsv`: what each contig of the draft is | yes | yes | yes | yes |
| `eval/qualimap/` | yes | yes | yes | |

Every file is described on [Output files](../reference/output.md).

[Choosing a mode](../getting-started/choosing-a-mode.md) explains which mode fits a dataset.
