# Coming from v1

This page explains how to move from FunFlux 1.x, Funnotator, FunFluxL or FunFluxL+ to FunFlux
v2. Version 2.0.0 replaces these four workflows, and the `mode` key selects what used to be a
separate workflow.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **NOTE** | a line starting `[FunFlux] NOTE` that reports something to check without stopping the run; read them after every dry run |
    | **launch folder** | the folder Snakemake is started from, normally the FunFlux repository; relative paths are read from it, and `.snakemake/` is created there |
    | **funannotate database** | the frozen snapshot of funannotate's data, such as Swiss-Prot, Pfam, the BUSCO set and stored AUGUSTUS models, downloaded once from Zenodo into `directories.funannotate_db` |
    | **path repair** | rule `normalize_funannotate_db`, which rewrites the absolute paths stored in the funannotate database to your folder, in place and without a lock |
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **confirmed mitochondrial contig** | a contig the mitochondrion screen classes `confirmed_mitochondrial`; with the default `parameters.mitochondria.mode: separate_confirmed`, it is kept whatever its genus |
    | **guide assembly** | in `hybrid` mode, the SPAdes assembly of the Illumina reads; it is screened, and only the pairs that map to it go on |

| Was | Input | Now |
|---|---|---|
| FunFlux 1.x | Illumina paired-end reads | `mode: illumina` |
| Funnotator, shipped inside FunFlux 1.x as `workflow/Funnotator` | assembled genomes (FASTA) | `mode: contigs` |
| FunFluxL 0.1 | ONT or PacBio HiFi reads | `mode: longread` |
| FunFluxL+, never published | Illumina pairs plus ONT or PacBio HiFi reads | `mode: hybrid` |

The FunFlux repository continues as v2 under the same concept DOI, 10.5281/zenodo.13612159.
Cite results produced with a v1 workflow by the version you ran, as described on
[Citation and references](../about/citation.md).

The v1 workflows are no longer in the main branch, but their releases stay available. FunFlux
v1.1.0, with Funnotator inside it, is the tag `v1.1.0` of this repository:
`git clone --branch v1.1.0 https://github.com/iLivius/FunFlux.git FunFlux-v1.1.0` puts it in
a folder of its own. FunFluxL 0.1.0 is the tag `v0.1.0` of the
[FunFluxL repository](https://github.com/iLivius/FunFluxL).

Neither a v1 configuration nor a v1 output folder can be reused as it is. Make these changes, in
this order:

1. Start from the v2 template, add `mode` and the `input:` block, and move the long-read settings
   into `parameters.longread` or `parameters.hybrid`
   ([Why a v1 configuration stops](#why-a-v1-configuration-stops)).
2. Name the configuration with `--configfile`, give `--cores` the value of `resources.threads`,
   and, on a single computer, drop `--jobs` ([Command line](#command-line)).
3. Give v2 a new, empty `output_dir` ([Output folders](#output-folders)).
4. Point `directories.funannotate_db` at a new, empty folder, to get the newer snapshot of the
   funannotate database ([Behaviour changes](#behaviour-changes)).
5. Replace `decontamination.mode: auto` and `parameters.genus`, if you used them
   ([Behaviour changes](#behaviour-changes)).
6. Update scripts that read renamed files or columns ([Output folders](#output-folders),
   [Behaviour changes](#behaviour-changes)).

The "What to do" column of [Behaviour changes](#behaviour-changes) lists the rest.

## Why a v1 configuration stops

A v1 configuration has no `mode` key, so FunFlux v2 stops before any job starts:

```text
[FunFlux] mode: this configuration has no 'mode' key. FunFlux v2 needs one of illumina | longread | hybrid | contigs, and an 'input:' block that names the input directory for that mode. A configuration written for FunFlux, FunFluxL or FunFluxL+ v1 does not run unchanged: start from the config template shipped with FunFlux v2 (see 'Coming from v1' in the documentation).
```

Adding `mode:` is not enough. The input folder moved out of `directories`, and FunFluxL and
FunFluxL+ kept their long-read settings at the top of `parameters`, where v2 stops on them:

```text
[FunFlux] mode=illumina requires input.illumina_dir to be set in the config.
[FunFlux] parameters.flye_input_mode: this setting moved into the block named after your mode: parameters.longread.flye_input_mode or parameters.hybrid.flye_input_mode.
```

Start from the v2 template and copy your paths and settings into it:

```bash
cp config/config.yaml config/config_custom.yaml
```

The table lists the v1 keys that moved or were renamed, and some new keys. The other new keys are
named under [Behaviour changes](#behaviour-changes).

| v1 key | v2 key | Note |
|---|---|---|
| none | `mode` | required: `illumina`, `longread`, `hybrid` or `contigs` |
| `directories.input_dir` | `input.illumina_dir`, `input.longread_dir` or `input.contigs_dir` | `hybrid` reads two keys, which may name the same folder; v2 ignores `directories.input_dir` without a message |
| none in FunFluxL 0.1, which assembled with Flye; `parameters.long_read_assembler` in FunFluxL+ | `parameters.<mode>.assembler` | `hifiasm` (the default) or `flye` |
| `parameters.flye_input_mode` | `parameters.<mode>.flye_input_mode` | the value `pacbio-hifi` is gone: name the reads `{sample}_pacbio.<ext>` |
| `parameters.medaka_model` | `parameters.<mode>.medaka_model` | |
| none in FunFluxL 0.1, which ran Filtlong with `--min_length 1000 --keep_percent 90 --target_bases 500000000`; `parameters.filtlong.*` in FunFluxL+ | `parameters.<mode>.filtlong.*` | `min_length` 1000, `keep_percent` 90, `length_weight` 10, `target_bases` empty (no cap); see [Behaviour changes](#behaviour-changes) |
| `parameters.filtlong.split` (FunFluxL+) | `parameters.hybrid.filtlong_split` | a `split` key inside `filtlong` stops the run |
| `parameters.hifiasm.*` (FunFluxL+) | `parameters.<mode>.hifiasm.*` | |
| `parameters.genus` | `parameters.decontamination` | translated with a NOTE when there is no `decontamination` block: a genus becomes `mode: include`, an empty value `mode: off` (v1: `mode: auto`, which v2 does not have) |
| none (FunFlux 1.x searched `core_nt`) | `parameters.nt_version` | now read in every mode |
| none (FunFlux 1.x kept 10 database sequences per contig, with 1 alignment each) | `parameters.blast_max_target_seqs`, `parameters.blast_max_hsps` | defaults 10 and 5 |

On this page, `<mode>` is `longread` or `hybrid`. The long-read settings live in the block named
after the mode, and only that block is read. All other v1 keys keep their names and meanings,
apart from the exceptions listed under [Behaviour changes](#behaviour-changes), such as
`decontamination.mode: auto`, which is gone. In addition, `'null'` in quotes as `medaka_model`
means off, as described on [The Medaka model](../methods/medaka-model.md#settings).

## Command line

The table compares the v1 and v2 command lines. Relative paths and the conda environments now
belong to the launch folder, not to `output_dir`.

| | v1 | v2 |
|---|---|---|
| Command | `snakemake --configfile config/config.yaml --sdm conda --cores 12 --jobs 2`; Funnotator added `--snakefile workflow/Funnotator` | `snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml` |
| Configuration | `config/config.yaml` was loaded by default | must be named with `--configfile`; without it FunFlux stops with `[FunFlux] No configuration supplied.` |
| CPUs | most rules took their threads from `resources.threads` through `resources: cpus`, which `--cores` does not limit, so tools could run on more threads than `--cores` allowed | every CPU-bound rule declares `threads:`, so `--cores`, equal to `resources.threads`, limits the whole run; no `--jobs` |
| Working folder | Snakemake moved into `output_dir` (`workdir:`): relative paths in the config were read from there, and the conda environments were built in `<output_dir>/.snakemake/conda/` | Snakemake stays in the launch folder: relative paths are read from there, and `.snakemake/` with the conda environments is created there |

## Output folders

!!! warning "FunFlux v2 does not upgrade a v1 output folder in place"

    Give v2 a new `output_dir`. Do not resume a v1 run with v2, and do not point v2 at a v1 output
    folder. The old layout is not converted, and some paths exist in both layouts
    (`02.assembly/{sample}/`, `04.annotation/{sample}/`), so old and new files would mix.

The v2 layout is the same in every mode, and a mode that has nothing for a folder leaves it out:

```text
<output_dir>/
├── 00.resources/    material shared by the run: UNITE, antiSMASH databases, telomere repeat database, SignalP 6, ChroQueTas, PhiX
├── 01.reads/        read QC, filtered and selected reads (not in contigs mode)
├── 02.assembly/     assembly, screen, mitochondrion split, QUAST, BUSCO, Qualimap
├── 03.taxonomy/     ITS extraction and classification
├── 04.annotation/   masking, gene prediction, functional annotation
├── 05.amr/          antifungal resistance mutations (new; only with the AMR screen on)
├── 06.report/       MultiQC
└── logs/            one log per rule and sample
```

**FunFlux 1.x and FunFluxL**

| v1 | v2 |
|---|---|
| `01.pre-processing/` (PhiX index, fastp reports, trimmed reads) | `00.resources/phix/`, `01.reads/{sample}/illumina/` |
| `01.pre-processing/{sample}_long_filt.fastq`, `01.pre-processing/read_qc/` (FunFluxL) | `01.reads/{sample}/long/` |
| `02.assembly/{sample}/` | `02.assembly/{sample}/` |
| `03.post-processing/contaminants/{sample}/`, and the read-mapping BAM in `03.post-processing/{sample}/` (FunFlux 1.x) or `03.post-processing/contaminants/{sample}_map.bam` (FunFluxL) | `02.assembly/{sample}/contaminants/` |
| `03.post-processing/consensus/{sample}/` (FunFluxL, Medaka) | `02.assembly/{sample}/medaka/` |
| `03.post-processing/assembly_evaluation/{sample}/` (FunFlux 1.x) | `02.assembly/{sample}/eval/quast/`, now also in `longread` mode |
| `03.post-processing/completeness_evaluation/{sample}/` | `02.assembly/{sample}/eval/busco/` |
| `03.post-processing/mapping_evaluation/{sample}/` | `02.assembly/{sample}/eval/qualimap/` |
| `03.post-processing/ITS_extraction/{sample}/` | `03.taxonomy/{sample}/` |
| `03.post-processing/ITS_extraction/unite_its_sintax.fasta` | `00.resources/unite/` |
| `04.annotation/{sample}/` | `04.annotation/{sample}/` |
| `04.annotation/antismash_db/` | `00.resources/antismash/databases/` |
| `report/` | `06.report/` |

**FunFluxL+**

| FunFluxL+ | v2 |
|---|---|
| `01.read_processing/illumina/`, `01.read_processing/long/` | `00.resources/phix/`, `01.reads/{sample}/illumina/`, `01.reads/{sample}/long/` |
| `02.short-read_assembly/{sample}/` | `02.assembly/{sample}/guide/` |
| `03.read_selection/contamination/{sample}/` | `02.assembly/{sample}/guide/contaminants/`; the selected pairs in `01.reads/{sample}/selected/` |
| `03.read_selection/mapping_evaluation/{sample}/` | `02.assembly/{sample}/guide/eval/qualimap/` |
| `03.read_selection/assembly_evaluation/{sample}/` and `07.genome_qc/assembly_evaluation/{sample}/` | `02.assembly/{sample}/eval/quast/`, one report with the labels `{sample}_guide` and `{sample}` |
| `04.long-read_assembly/{sample}/` | `02.assembly/{sample}/` (`hifiasm/` or `flye/`, `assembly.fasta`) |
| `05.assembly_curation/contamination/{sample}/` | `02.assembly/{sample}/contaminants/` |
| `05.assembly_curation/mapping_evaluation/{sample}/` | `02.assembly/{sample}/eval/qualimap/` |
| `06.polishing/medaka/{sample}/`, `06.polishing/polypolish/{sample}/` | `02.assembly/{sample}/medaka/`, `02.assembly/{sample}/polypolish/` |
| `07.genome_qc/completeness_evaluation/{sample}/` | `02.assembly/{sample}/eval/busco/` |
| `07.genome_qc/ITS_extraction/{sample}/` | `03.taxonomy/{sample}/` |
| `07.genome_qc/ITS_extraction/unite_its_sintax.fasta` | `00.resources/unite/` |
| `08.annotation/{sample}/`, `08.annotation/antismash_db/` | `04.annotation/{sample}/`, `00.resources/antismash/databases/` |
| `09.report/` | `06.report/` |

**Funnotator**

| Funnotator | v2 |
|---|---|
| `01.pre-processing/pre_filtering/{sample}/contigs_filt.fasta` | `02.assembly/{sample}/contigs_final.fasta` |
| `01.pre-processing/assembly_evaluation/{sample}/` | `02.assembly/{sample}/eval/quast/` |
| `01.pre-processing/completeness_evaluation/{sample}/` | `02.assembly/{sample}/eval/busco/` |
| `01.pre-processing/ITS_extraction/{sample}/` | `03.taxonomy/{sample}/` |
| `01.pre-processing/ITS_extraction/unite_its_sintax.fasta` | `00.resources/unite/` |
| `02.annotation/{sample}/`, `02.annotation/antismash_db/` | `04.annotation/{sample}/`, `00.resources/antismash/databases/` |
| `report/` | `06.report/` |

The delivered assembly files have new names, and the mitochondrion now has a file of its own:

| v1 | v2, in `02.assembly/{sample}/` |
|---|---|
| `contigs_sel.fasta` (FunFlux 1.x), the Medaka `consensus.fasta` (FunFluxL), `polished.fasta` (FunFluxL+) | `contigs_all.fasta`: the complete assembly after the screen and polishing, mitochondrion included |
| none | `contigs_final.fasta`: the nuclear genome, read by QUAST, BUSCO, the ITS steps and funannotate |
| none | `contigs_mitochondrial.fasta`: the confirmed mitochondrial contigs; an empty file when there are none |

A new file, `04.annotation/{sample}/contigs_removed_by_clean.txt`, lists the contigs that
`funannotate clean` removed before gene prediction.

## Behaviour changes

The tables below group the changes by part of the workflow. The "Affects" column names the v1
workflows that a change concerns, and the "What to do" column names the action. "Nothing" means
that the change needs no action.

**Running and configuration**

| Change | Affects | What to do |
|---|---|---|
| One workflow with a `mode` key; the configuration must be named with `--configfile` | all | Start from the v2 template; add `mode` and the `input` block |
| Every CPU-bound rule declares `threads:`, so `--cores` limits the whole run | all | Give `--cores` the value of `resources.threads` (24 in the template, as in v1); drop `--jobs` on a single computer |
| FunFlux runs on a cluster with `--executor slurm --jobs N`. Every heavy step asks the scheduler for `resources.ram_gb` of memory and for the time limit of the new key `resources.runtime_hours` (default 48). Every download step runs where Snakemake runs, and `snakemake download_all` fetches all downloads in advance ([On a cluster](../reference/running.md#on-a-cluster)) | all | On a cluster, set the three `resources` keys from the limits of a node and a partition |
| curl, tar and gzip come from a conda environment (`download`), and InterProScan's Perl from the `interproscan` environment; v1 took them from the computer | all | Nothing |
| No `workdir:`: relative paths are read from the launch folder | all | Use absolute paths, or paths relative to the launch folder |
| Thread caps: eggNOG-mapper 12 (FunFlux 1.x: 24), nhmmer and ITSx 8, SINTAX 1, SignalP 6 at most `resources.threads` ([CPU cores and threads](../reference/running.md#cpu-cores-and-threads)) | all | Nothing |
| Conda environments built from pinned package lists of the tested tool versions (`workflow/envs/*.linux-64.pin.txt`), with the Python version pinned in every environment that contains Python | all | Nothing; the first v2 run builds new environments |
| New optional keys whose defaults keep the v1 behaviour: `spades_kmers`, `busco.lineage`, `prediction.busco_db`, `repeatmodeler_seed`, `directories.antismash_db`, `directories.unite_db`, `contigs.screen`, `medaka_model_fallback_auto` | all | Nothing |

**Inputs**

| Change | Affects | What to do |
|---|---|---|
| One input key per technology: `input.illumina_dir`, `input.longread_dir`, `input.contigs_dir` replace the v1 key `directories.input_dir` | all | Move the path |
| The suffixes stay `_ont` and `_pacbio`, and `flye_input_mode: pacbio-hifi` is removed. Illumina plus PacBio HiFi is allowed in `hybrid` | FunFluxL, FunFluxL+ | Name PacBio HiFi reads `{sample}_pacbio.<ext>` |
| Underscores are allowed in sample names (v1 rejected them). FunFlux prints a NOTE when a name is not a valid locus-tag prefix | all | For GenBank submission, use 3 to 12 letters or digits starting with a letter |

**Assembly and screen**

| Change | Affects | What to do |
|---|---|---|
| Long reads are assembled with hifiasm by default. FunFluxL 0.1 used Flye, which stays available | FunFluxL, FunFluxL+ | Set `parameters.<mode>.assembler: flye` to keep Flye |
| Since FunFlux 2.1.0, the ONT adapter and the barcode after it are trimmed from the reads of ONT samples with fastplong before Filtlong. By default (`parameters.<mode>.ont_adapter: detect`) FunFlux works out from the reads what to trim, and leaves reads that are already trimmed as they are. FunFluxL 0.1 and FunFluxL+ assembled the reads untrimmed ([Adapter trimming](../methods/adapter-trimming.md)) | FunFluxL, FunFluxL+, ONT samples | Nothing; `ont_adapter: off` assembles the reads untrimmed, as v1 did |
| Filtlong keeps every read that passes its length and quality filters: `target_bases` is empty (no cap) by default, and `length_weight: 10` favours longer reads. FunFluxL 0.1 kept at most 500 Mb (`--target_bases 500000000`) with Filtlong's default length weight of 1 | FunFluxL | To repeat FunFluxL 0.1's read selection, set `parameters.longread.filtlong.target_bases: 500000000` and `length_weight: 1` |
| Before the assembly starts, a named Medaka model is checked against Medaka's list of models, and `auto` is worked out from the basecaller tag of the first 100 reads of the raw FASTQ. A wrong model therefore fails in seconds instead of after the assembly. `auto` needs network access once and exactly one basecaller tag in those reads. New key `medaka_model_fallback_auto`. The assembler waits for the check's file, so a changed `medaka_model` also re-runs the assembly of the ONT samples | FunFluxL, FunFluxL+, ONT samples | For public or re-headered ONT reads, name the model or set `medaka_model: false` ([The Medaka model](../methods/medaka-model.md)) |
| The PhiX genome and its bowtie2 index are downloaded once per output folder and kept in `00.resources/phix/`. In v1, they were kept with the processed reads and deleted after the run. The PhiX-free pairs are called `{sample}_nophix_R1.fastq` and `_R2.fastq` | FunFlux 1.x, FunFluxL+ | Nothing |
| SPAdes writes into `02.assembly/{sample}/spades/` (in `hybrid`, `guide/spades/`), so a re-run of SPAdes no longer deletes the screen's files of the sample. In FunFlux 1.x the whole sample folder was SPAdes' output. Flye likewise works in `flye/`, and its assembly is copied up one level | FunFlux 1.x, FunFluxL, FunFluxL+ | Nothing |
| Filtlong guidance by the Illumina reads is off by default. Switched on (`short_read_guidance: true`, `filtlong_split: 1000`), it applies to ONT samples only. FunFluxL+ always guided and cut every sample, PacBio included | FunFluxL+ | For FunFluxL+'s behaviour on ONT samples, set `short_read_guidance: true`; move `filtlong.split` to `parameters.hybrid.filtlong_split` ([Filtlong guidance](../methods/filtlong-guidance.md)) |
| Polypolish output headers are cut at the first space or tab (v1: `>ptg000001l polypolish`) | FunFluxL+ | Nothing |
| `assembly_linearized.fasta` and `assembly_single_line.fasta` are no longer written | FunFluxL, FunFluxL+ | Update scripts that read them |
| Since FunFlux 2.1.0 the BLAST search keeps 10 database sequences per contig, with up to 5 alignments each, in every mode, and drops matches with an e-value above 1e-10 (FunFlux 2.0.0: 50, 5 and 1e-5; FunFlux 1.x: 10, 1 and 1e-5; FunFluxL and FunFluxL+: 50, 5 and 1e-5). FunFlux prints a NOTE when `blast_max_hsps` is below 5 and the mitochondrion screen is on | FunFlux 1.x, FunFluxL, FunFluxL+ | Nothing. To come close to earlier numbers, set `blast_max_hsps: 1` (FunFlux 1.x) or `blast_max_target_seqs: 50` (FunFluxL, FunFluxL+) ([BLAST hit depth](../methods/blast-hit-depth.md)) |
| `parameters.nt_version` is read in `illumina` mode too | FunFlux 1.x | Set it if your BLAST database has another name than `core_nt` |
| New: a mitochondrion screen and split in every mode with the screen. The screen classes each contig of the draft from its BLAST hits and, in `longread` and `hybrid`, from the assembler's read depth and circularity. With the default `parameters.mitochondria.mode: separate_confirmed`, the confirmed mitochondrial contigs are kept whatever their genus, on the `hybrid` guide assembly too. After polishing they are split off into `contigs_mitochondrial.fasta`, so that QUAST, BUSCO, the ITS steps and funannotate read only the nuclear genome, `contigs_final.fasta`. Since FunFlux 2.1.0, a contig whose best mitochondrial match is an animal's or a plant's, such as human mitochondrial DNA from the library, is not confirmed (`foreign_mitochondrial`): it stays in the assembly, where the genus filter treats it like any other contig. A confirmed mitochondrial contig missing from the BlobTools table is kept as `no-taxonomy`, and a confirmed ID missing from the FASTA stops the run ([Mitochondria](../analysis/mitochondria.md)) | all | Read `contigs_final.fasta` as the genome; `parameters.mitochondria.mode: report` or `off` keeps the mitochondrial contigs in it |
| `decontamination.mode: auto` is removed and stops the run before any job starts, also as a mode in `sample_overrides` when the screen runs. BLAST genus labels are unreliable for poorly sequenced fungi, so `auto` can remove most of the genome ([Configuration](../reference/configuration.md#decontamination)) | FunFlux 1.x, FunFluxL, FunFluxL+ | Use `off` (the default), read `contig_taxonomy_decisions.tsv`, then `exclude` the contaminant genera or `include` the genera BLAST assigned to the isolate's own contigs |
| With the screen, the decontamination files (`include_genera_by_sample`, `exclude_genera_file`, `sample_overrides`) and an `include` mode without a genus are checked before any job starts (v1: when `select_contigs` ran, after BLAST) | FunFlux 1.x, FunFluxL, FunFluxL+ | Nothing |
| `parameters.genus` is translated with a NOTE, and an empty `genus` now means `mode: off` (v1: `auto`, the genus with the most contigs). The template sets `discard_no_hit: false`, where the v1 template had `true` (the v1 code default was already `false`) | all | Replace `genus` with `parameters.decontamination`; set `discard_no_hit: true` if you relied on the v1 template |
| With SPAdes headers, the `contigs` mode filter keeps contigs of at least 500 bp and 2x coverage. Any other FASTA only has its headers cut at the first space or tab (Funnotator removed every contig of 500 bp or less). An empty result stops the run, and the screen is optional | Funnotator | Expect 500 bp contigs to be kept and no length filter on non-SPAdes FASTA |
| Delivered files renamed (table above) | all | Update downstream scripts |

**QC and ITS**

| Change | Affects | What to do |
|---|---|---|
| QUAST also runs in `longread` mode. In `hybrid`, one QUAST report holds the guide (`{sample}_guide`) and the final assembly (`{sample}`) | FunFluxL, FunFluxL+ | Nothing |
| BUSCO's result tree is now kept, except `hmmer_output/` and `busco_sequences/` (v1 deleted `run_*/` and `auto_lineage/`). With the `auto` lineage, that is up to about 2 GB per sample. With `auto`, the lineage data is still downloaded per sample into a temporary folder, while a named lineage is downloaded once per output folder into `00.resources/busco/` and read offline. New key `busco.lineage`. With a named lineage, MultiQC shows one BUSCO row per sample instead of two | all | Allow up to about 2 GB per sample, or name a lineage, for example `hypocreales_odb12`, which keeps the result tree of that lineage only |
| New: the contig summary, `02.assembly/{sample}/eval/contig_summary.tsv`. It has one row per contig of the draft assembly, with a verdict (`chromosome`, `chromosome_arm`, `mitochondrion`, `rDNA_array`, `low_depth`, `repeat`, `nuclear_contig`, `nuclear_fragment`, `contaminant`, `foreign_mitochondrion` or `unclassified`) and the evidence behind it, and is also shown as a table in the MultiQC report. Its `annotated_as` column links each contig to its `scaffold_N` name in the annotation ([Contig summary](../analysis/contig-summary.md)). The telomere repeat unit of each assembly is chosen from the data with tidk and checked against tidk's telomere repeat database, downloaded once per output folder. New keys `parameters.contig_summary.*`, new required link `links.telomere_db_link` | all | Copy the `links.telomere_db_link` line from the v2 template; the run stops without it |
| New since FunFlux 2.1.0: the contig adapter screen. Every draft assembly is searched with BLAST for ONT, PacBio, MGI and Illumina adapters and barcodes; the loci go to `02.assembly/{sample}/eval/adapter_screen.tsv`, to the `adapter` column of the contig summary and to the Adapters table of the MultiQC report. The Illumina records of NCBI UniVec are downloaded once per output folder. New optional link `links.univec_link` and key `parameters.contig_summary.univec_sha256` ([The contig adapter screen](../methods/adapter-trimming.md#the-contig-adapter-screen)) | all | Nothing; the screen only reports, and changes no genome file |
| ITS is cut out of `contigs_final.fasta` in one strand-oriented window per rDNA unit before ITSx runs (Funnotator ran ITSx on its raw input). New keys `its.window_flank` (2,000 bp on each side of the 5.8S gene) and `its.max_full_length` (1,500 bp). Identical ITS copies are collapsed and counted. ITS1 or ITS2 outside 50-600 bp is listed in `{sample}_ITS_copies.txt` and not classified, and a full ITS longer than `its.max_full_length` is marked `rejected_too_long` there (only ITS1 and ITS2 are classified, as in v1). Long-read assemblies, which could give empty or wrong ITS files in v1, now give the ITS | all | Update scripts that parse ITS record names: records are named after their window (`contig:start-end_strand`), classified copies `{sample}_ITS1_1`, `{sample}_ITS2_1`, ... ([ITS on long-read contigs](../methods/its-long-read-contigs.md)) |
| SINTAX runs on one thread with `--randseed 1`, so the ITS confidences repeat from run to run. In v1, SINTAX used several threads and no fixed seed | all | Nothing |
| SINTAX classifies against the UNITE records whose taxonomy names a genus (new rule `unite_identified`, new key `its.unite_identified_only`, default `true`). In v1, SINTAX used the whole release, in which about a third of the records stop above the genus, and one identical to the query can stop the result at `d:Fungi` | all | Nothing; `its.unite_identified_only: false` gives the v1 reference |
| QC-only rules (QUAST, Qualimap, NanoPlot) no longer stop the run when the tool fails | all | Check the log and the MultiQC page for a missing table |
| NanoPlot writes its plots as HTML pages only (`--no_static`). FunFluxL also wrote them as PNG, but the pinned NanoPlot 1.46.2 draws PNG only where a Chrome or Chromium browser is installed, so the output would differ between computers | FunFluxL, FunFluxL+ | Open the HTML plots; MultiQC reads the NanoStats table as before |
| Qualimap runs without a display, so it works on a server with no X session. Its Java heap is `min(ram_gb, 64)` GB, as in v1 | FunFlux 1.x, FunFluxL, FunFluxL+ | Nothing |
| BUSCO runs with `--opt-out-run-stats`, so it sends no usage statistics and makes no ipinfo.io lookup (v1 sent both). MultiQC no longer shortens sample names, so a sample named `SAMPLE1.clean` keeps its name | all | Nothing |

**Annotation and databases**

| Change | Affects | What to do |
|---|---|---|
| The antiSMASH database stays in `00.resources/antismash/databases/` between runs on the same output folder (v1 deleted it and downloaded it again on every resume), or is linked from `directories.antismash_db`. The same holds for UNITE, with `PROVENANCE.txt` | all | Optionally point `directories.antismash_db` and `directories.unite_db` at copies you have |
| The funannotate database archive is checked with `sha256sum -c` before it is unpacked (v1 downloaded the checksum file but did not check it). A mismatch stops the run | all | Delete the archive and its `.sha256` file and run again |
| The template links a newer snapshot of the funannotate database, `funannotate_db_1.8.17_2026-10-08` (Zenodo record 23265189), with newer Swiss-Prot, Pfam, InterPro, Gene Ontology and gene2product releases. The v1 templates link the 2025-12-15 snapshot (record 18271295). What changed is listed under [funannotate database snapshot](databases.md#funannotate-db) | all | Point `directories.funannotate_db` at a new, empty folder, because a folder that holds the v1 snapshot keeps being used |
| Since FunFlux 2.1.0, InterProScan runs with `-dp`: it computes every match itself and does not ask EBI's lookup service for matches already computed (the service did not answer InterProScan 5 in November 2026), so it needs no internet connection ([Functional annotation](../analysis/functional-annotation.md#interproscan)) | all | Nothing |
| With GeneMark-ES on, the new rule `genemark_check` runs GeneMark-ES on a 200 kb slice of each masked genome before gene prediction (about a minute). An expired licence key therefore stops the run before `funannotate predict` starts instead of partway through it | all, with GeneMark-ES | Renew the key when the message asks for it |
| The AUGUSTUS model (column 4 of the annotation parameter table, `files.annotation_params`) is optional. Empty or `-` makes funannotate train AUGUSTUS on the genome's BUSCO genes, which v1 never did, because v1 required a model. A named model must exist in the funannotate database, or the run stops with the closest names (v1 passed any name on, and funannotate silently trained AUGUSTUS under a misspelt one). The new rule `augustus_model_check` records the choice in `04.annotation/{sample}/augustus_model.txt` and at the top of the gene prediction log. `parameters.prediction.augustus_optimize` adds the long training | all | Nothing for a v1 table with a valid model name: AUGUSTUS gets the same option as in v1 |
| Tiberius, a gene finder with a model trained by its authors on many fungal genomes, runs by default before gene prediction, and `funannotate predict` adds its gene models to the evidence it combines (`--other_gff`, weight 1). It needs about 16 GB of memory and downloads its model (30 MB) once per output folder. New keys `parameters.tiberius.mode` and `model_sha256`, new link `links.tiberius_model_link` | all | Nothing; `parameters.tiberius.mode: off` gives gene prediction without it, as in v1 |
| New key `parameters.genemark.mode`: `auto`, the default, uses GeneMark-ES when `directories.genemark_dir` is set, and `off` never uses it. An empty value (`mode:` or `mode: ""`) means `auto`, as for the other mode keys | all | Nothing; leave `directories.genemark_dir` empty to run without GeneMark-ES |
| New, off by default: SignalP 6, which predicts signal peptides for funannotate annotate. It needs the licensed archive from DTU (`files.signalp6_archive`) and a `parameters.signalp6.mode` other than `off` ([Licensed and external tools](external-tools.md)) | all | Nothing |
| New: rule `normalize_funannotate_db`, the path repair, rewrites the absolute paths stored inside the funannotate database to your folder. It runs in place and without a lock, so let one run pass it before others start on the same database ([Reference databases](databases.md#what-the-workflow-does)) | all | Nothing |
| New, off by default: the AMR screen. ChroQueTas screens each genome for the point mutations that the FungAMR database lists with evidence for or against resistance, and writes `05.amr/{sample}/`. New keys `parameters.amr.*`, new link `links.chroquetas_link` ([Antifungal resistance mutations](../analysis/amr.md)) | all | Nothing; `parameters.amr.chroquetas: on` switches it on, and accepts the non-commercial licence of ChroQueTas |

**Logs**

| Change | Affects | What to do |
|---|---|---|
| One log per rule and sample in `logs/`; InterProScan's scratch folder is `04.annotation/{sample}/iprscan/tmp/` (no top-level `temp/`) | all | Nothing |

**Kept as in v1**

- `decontamination.mode: off` as the default
- BUSCO on the final genome only
- the filtered long reads, the selected pairs and the read-mapping BAM of the read-based modes,
  kept after the run
