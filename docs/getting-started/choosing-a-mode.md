# Choosing a mode

Choose the mode from the data you have: Illumina paired-end reads, ONT or PacBio HiFi reads, both
kinds of reads for the same isolates, or an existing assembly. This page shows which mode fits
which data, how to name the input files, and what happens when the mode or a file name is wrong.

The `mode` key selects the front end, the steps that differ between modes. These are read
processing, assembly and the contig filter before the contamination and mitochondrion screen, and
polishing after it. From the screen on, every mode runs the same steps except polishing, which
differs between the read-based modes.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **front end** | the steps that differ by mode, chosen by `mode`: read processing, assembly and contig filter before the screen, and polishing after it |
    | **the screen** | the contamination and mitochondrion screen: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **batch** | the samples one run finds in the input folders of its mode; files of one kind share one extension |
    | **guide assembly** | in `hybrid`, the SPAdes assembly of the Illumina reads; it is screened, and only the pairs that map to it go on |
    | **dry run** | a Snakemake run with `-n`: FunFlux checks the configuration and the inputs, and Snakemake lists the jobs it plans without running any; it takes seconds |
    | **NOTE** | a line starting `[FunFlux] NOTE` that reports something to check without stopping the run; read them after every dry run |
    | **nuclear genome** | the delivered genome without the confirmed mitochondrial contigs, `contigs_final.fasta`; it holds every contig when nothing is split off (no screen, or `parameters.mitochondria.mode` `off` or `report`) |

## Which mode for which data

Each kind of data has its own mode, input key and file-name convention, and the input files must
follow that convention.

| Your data | `mode` | Input key | File names |
|---|---|---|---|
| Illumina paired-end reads | `illumina` | `input.illumina_dir` | `{sample}_R1.<ext>` and `{sample}_R2.<ext>` |
| ONT reads, PacBio HiFi reads, or a batch with both | `longread` | `input.longread_dir` | `{sample}_ont.<ext>` or `{sample}_pacbio.<ext>` |
| Illumina pairs plus ONT or PacBio HiFi reads of the same isolates | `hybrid` | `input.illumina_dir` and `input.longread_dir` (may be the same folder) | both conventions, every sample on both sides |
| An assembly in FASTA format | `contigs` | `input.contigs_dir` | `{sample}.fasta`, `{sample}.fa` or `{sample}.fna` |
| Some isolates with Illumina reads only, others with ONT or PacBio HiFi reads only | one run per mode: `illumina` for the first group, `longread` for the second, each with its own configuration file | each run its own input folder | as in those rows |

`<ext>` is `fastq`, `fq`, `fastq.gz` or `fq.gz`. Only the input key of the active mode is read, and
the others may stay empty. Isolates with different kinds of data go into separate runs, as the
last row shows. These runs may share one `output_dir` and one annotation parameter table when
their sample names differ, as described in
[Running several batches](../reference/running.md#running-several-batches).

The modes differ in the steps below, and the [Modes](../modes/index.md) pages describe each mode
in detail.

| | `illumina` | `longread` | `hybrid` | `contigs` |
|---|---|---|---|---|
| Reads | PhiX removal (bowtie2), adapter and quality trimming (fastp) | ONT adapter and barcode removal (fastplong); Filtlong; NanoPlot before and after | both; Filtlong as in `longread`, optionally guided by the Illumina pairs | none |
| Assembly | SPAdes `--isolate`; contigs of at least 500 bp and 2x coverage kept | hifiasm (default) or Flye | SPAdes guide assembly of the Illumina reads, then hifiasm or Flye of the long reads | your FASTA ([below](#contigs-from-another-assembler)) |
| Polishing | none | Medaka, ONT samples only | Medaka for ONT samples, then Polypolish with the Illumina reads | none |
| Mitochondrion | found by BLAST | found by BLAST, supported by read depth or circularity from the assembler; a BLAST hit covering at least 70 % confirms alone | as `longread`; the guide assembly is screened as well | only with `parameters.contigs.screen: true`, by BLAST |

`contigs` mode starts from a finished assembly and has no read steps. `hybrid` does the most work,
because it assembles twice. The Illumina reads build a guide assembly whose screen selects the
isolate's read pairs, and Polypolish corrects the final long-read assembly with these selected
pairs. If you switch on the optional Filtlong guidance, the same pairs also score the ONT reads.
`hybrid` needs both read types for every sample.

## File names

FunFlux finds the samples from their file names. A file that breaks its mode's convention stops
the run or is not found, so check the names before the first run.

| Item | Accepted | Not accepted |
|---|---|---|
| Sample name | the part of the file name before `_R1`/`_R2`, `_ont`/`_pacbio`, or, in `contigs` mode, before the last dot: `SAMPLE1_R1.fastq.gz` is sample `SAMPLE1`, `my.genome.fasta` is sample `my.genome`; underscores and dots are allowed | a space or one of the characters listed below the table stops the run |
| Extension | one extension for all files of one kind in a batch | a folder with both `.fastq.gz` and `.fq.gz` reads stops the run (`more than one file extension in this folder`) |
| Illumina reads | paired: one R1 and one R2 file per sample, and the name must end in `_R1.<ext>`; concatenate the files of several lanes or runs first | an R1 file without its R2 stops the run; a name with more after `_R1`, such as `SAMPLE1_S1_L001_R1_001.fastq.gz`, is not found, so rename or link the files |
| Long reads | one file per sample | an `_ont` and a `_pacbio` file for the same sample stop the run |
| `contigs` FASTA | uncompressed FASTA only | `.fasta.gz` stops the run; every file in the folder whose name contains a dot is read as a sample, so a `README.txt` there stops the run as well |

The characters that stop the run in a sample name are a space and any of
`* # @ % ^ / ! ? & : ; | < >`.

The sample name also becomes the locus-tag prefix of every gene that funannotate writes. NCBI
accepts a locus-tag prefix of 3 to 12 letters or digits that starts with a letter. A sample name
that does not fit still runs, and a NOTE then says that the gene IDs cannot be submitted to
GenBank as they are. Every sample needs a line in the
[annotation parameter table](../reference/configuration.md#annotation-parameter-table).

## Mixing ONT and PacBio HiFi samples

In `longread` and `hybrid` mode, the file-name suffix sets the technology of each sample.
`{sample}_ont.<ext>` is Oxford Nanopore and `{sample}_pacbio.<ext>` is PacBio HiFi (CCS). One batch
may hold both, but one sample must not have both files. Some settings are chosen once for the
whole run, and others for each sample.

| Setting | Chosen |
|---|---|
| The assembler (`parameters.<mode>.assembler`) | once for the whole run |
| The Filtlong settings | once for the whole run, for both technologies; their defaults are sized for ONT reads |
| ONT adapter trimming, the assembler's flags, the minimap2 preset, Medaka and, in `hybrid`, the Filtlong guidance when it is on | per sample, from the suffix; the dry run prints the decision for every sample ([One technology per sample](../modes/longread.md#one-technology-per-sample)) |
| Everything after polishing (`contigs_all.fasta`) | the same for both technologies |

One batch can therefore mix the two technologies, and only the steps up to polishing treat them
differently. In `hybrid`, PacBio samples get Filtlong's plain length and quality filter and no
Medaka, and Polypolish still corrects them with the Illumina reads.

## PacBio CLR reads are not supported

`_pacbio` means HiFi (CCS) reads, and no step is set up for continuous long reads (CLR). For a
`_pacbio` sample, hifiasm runs in its HiFi mode, Flye gets `--pacbio-hifi`, minimap2 gets
`map-hifi`, and Medaka is skipped. FunFlux does not measure read accuracy, so a `_pacbio` file of
CLR reads is assembled with these settings without an error. The only message is the NOTE of
every dry run with `_pacbio` samples:

```text
[FunFlux] NOTE: samples with a _pacbio file are assumed to hold PacBio HiFi (CCS) reads: GXU8616. Continuous long reads (CLR) are not supported by this workflow.
```

## Contigs from another assembler

In `contigs` mode, the contig filter looks at the first header of each FASTA, and only a FASTA
with SPAdes headers loses contigs at this step.

| First header | Contigs removed | Headers |
|---|---|---|
| SPAdes (`>NODE_<n>_length_<L>_cov_<C>`) | shorter than 500 bp or with a coverage below 2, as in `illumina` mode | stay as they are |
| Any other (Flye, hifiasm, a sequencing provider, a database download) | none | cut at the first space or tab |

`logs/filter_contigs_{sample}.log` names the branch taken and reports `Kept N of M records.`, and
the run stops when no contig is left.

Later, before gene prediction, `funannotate clean` removes contigs shorter than 500 bp and
duplicated contigs shorter than the N50, such as duplicated haplotypes in a provider assembly. Its
rules are given on [Repeat masking](../analysis/repeat-masking.md#funannotate-clean).
`04.annotation/{sample}/contigs_removed_by_clean.txt` lists the removed contigs, while QUAST,
BUSCO and the ITS steps still see every contig of `contigs_final.fasta`.

## When to switch on the screen in `contigs` mode

In `contigs` mode the screen runs only with `parameters.contigs.screen: true`, and it then needs
`directories.blast_db`, the folder that holds NCBI core_nt:

```yaml
directories:
  blast_db: /path/to/NCBI_core_nt

parameters:
  contigs:
    screen: true
```

The screen searches every contig with BLAST against NCBI core_nt, assigns it a genus with
BlobTools and tests it for a match to a mitochondrial genome. With the default
`parameters.mitochondria.mode: separate_confirmed`, the screen moves confirmed mitochondrial
contigs to `02.assembly/{sample}/contigs_mitochondrial.fasta`. `contigs_final.fasta` is then the
nuclear genome, and `contigs_all.fasta` holds both.

Switch the screen on when the FASTA may hold a mitochondrion or you want to look for
contamination, and core_nt is at hand.

| Screen | When |
|---|---|
| on | The assembly may contain the mitochondrial genome; an assembly from a sequencing provider usually does. Without the screen it is annotated with the nuclear gene predictor, although fungal mitochondria use a different genetic code (code 4 in most fungi, code 3 in *Saccharomyces cerevisiae* and related yeasts). |
| on | You want a BLAST taxonomy per contig to look for contamination: `02.assembly/{sample}/contaminants/contig_taxonomy_decisions.tsv` has a line per contig and `{sample}_composition.txt` the share of each genus. |
| off | The FASTA is already a nuclear genome, for example the `contigs_final.fasta` of an earlier FunFlux run that split the mitochondrion off. |
| off | core_nt is not available: the screen needs `directories.blast_db` with core_nt, `taxdb.btd`, `taxdb.bti`, `nodes.dmp` and `names.dmp`. |

With the default `parameters.decontamination.mode: off`, the screen removes no contig and only
splits off the confirmed mitochondrion. The screen adds 8 jobs per sample.

Without reads, the contigs are mapped against themselves only because BlobTools needs a BAM file.
That coverage carries no information, and a contig is confirmed as mitochondrial by BLAST alone,
when a mitochondrial genome covers at least `min_direct_blast_qcov` (35 %) of it, unless that
genome is an animal's or a plant's (class `foreign_mitochondrial`). See
[Decontamination](../analysis/decontamination.md) and [Mitochondria](../analysis/mitochondria.md).

## What a wrong mode does

The four file-name conventions do not overlap, so most wrong modes and folders stop the dry run.
The checks run in a fixed order and stop at the first problem.

A wrong `mode` value:

```text
[FunFlux] mode: must be one of illumina | longread | hybrid | contigs (got: 'ilumina')
```

The input key of the mode is not set, or is not a folder:

```text
[FunFlux] mode=illumina requires input.illumina_dir to be set in the config.
[FunFlux] input.illumina_dir: '/path/to/reads' is not an existing directory.
```

`illumina` on long reads, or on files named `_R1_001.fastq.gz`:

```text
[FunFlux] input.illumina_dir: no input files found. Expected names like {sample}_R1.fastq.gz and {sample}_R2.fastq.gz.
```

`longread` on Illumina reads:

```text
[FunFlux] input.longread_dir: no input files found. Expected names like {sample}_ont.fastq.gz or {sample}_pacbio.fastq.gz.
```

`contigs` on a folder of reads, or on compressed FASTA:

```text
[FunFlux] input.contigs_dir: file extension 'gz' is not recognised (expected one of: fasta, fa, fna).
```

`hybrid` with a sample that has reads of one kind only:

```text
[FunFlux] hybrid mode needs Illumina and long reads for every sample.
  No _ont or _pacbio file for: GXU8616
```

A sample with both long-read suffixes:

```text
[FunFlux] sample 'SAMPLE1' has both _ont and _pacbio files in /path/to/reads. Provide one long-read technology per sample.
```

Without `directories.blast_db`, a run with the screen (every read-based mode, and `contigs` mode
with the screen on) stops before the input folder is read:

```text
[FunFlux] directories.blast_db: required, but it is missing or empty in the config.
```

!!! warning "Three wrong choices do not stop the run"

    - **`illumina` or `longread` on a folder that holds both kinds of reads** runs with the files
      of its own convention and ignores the others. If both belong to the same isolates, the mode
      is `hybrid`.
    - **`contigs` without the screen** annotates a mitochondrial genome in the FASTA as nuclear.
      The NOTE printed on every such run says so.
    - **CLR reads named `_pacbio`** are assembled with HiFi settings.
