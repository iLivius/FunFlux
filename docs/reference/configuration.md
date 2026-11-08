# Configuration

One YAML file describes a FunFlux run: its mode, its input and output folders, the databases, the
resources and the analysis settings. This page lists every key of that file, the modes that read
it, its default and its allowed values. The template is `config/config.yaml`. Copy it, fill in the
copy and pass it with `--configfile`:

```bash
cp config/config.yaml config/config_custom.yaml
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml
```

Without `--configfile`, the run stops at once with `[FunFlux] No configuration supplied.` Git
ignores `config/config_custom.yaml`, so your filled-in copy is never committed.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **YAML** | the text format of the configuration file, with keys indented under their block; unquoted `yes`, `on`, `true` read as true and `no`, `off`, `false` as false |
    | **NOTE** | a line starting `[FunFlux] NOTE` that reports something to check without stopping the run; read them after every dry run |
    | **launch folder** | the folder Snakemake is started from, normally the FunFlux repository; relative paths are read from it, and `.snakemake/` is created there |
    | **the screen** | the contamination and mitochondrion screen of the draft assembly: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector |
    | **batch** | the samples one run finds in its mode's input folder; each needs a row in the annotation parameter table |
    | **funannotate database** | the frozen snapshot of funannotate's data, such as Swiss-Prot, Pfam, the BUSCO set and stored AUGUSTUS models, downloaded once from Zenodo into `directories.funannotate_db` |
    | **checksum (SHA-256)** | a 64-character fingerprint of a file's content, printed by `sha256sum`; FunFlux checks several downloads against stored checksums and stops on a mismatch |
    | **heavy step, long step** | on a cluster, heavy steps ask the scheduler for `ram_gb` and `runtime_hours`, and long steps for `runtime_hours` only |

## How the file is read

FunFlux checks the configuration before any job starts, so a dry run (`snakemake -n`) finds
configuration errors in seconds. The first failing check stops the run, so repeat the dry run
until it passes. A line `[FunFlux] NOTE: ...` does not stop the run.
[Troubleshooting](../troubleshooting.md) lists the messages.

FunFlux treats a wrong, empty, missing or unknown entry as follows.

| In the file | What FunFlux does |
|---|---|
| a value outside the allowed values listed below | stops the run with a message that starts `[FunFlux] <key>:` |
| `nan` or `inf` (YAML `.nan`, `.inf`) where a number is expected | refuses it, as above |
| a missing optional key | uses the default given in the tables below. The template ships the same value unless a table says otherwise |
| a missing or empty required key | stops the run |
| a number key that is present but empty | stops the run, unless its table lists empty as allowed. Delete the line to get the default |
| an unknown or misplaced key under `parameters` | stops the run with a message that lists the known keys or names the right place |
| a misspelt optional key in the blocks `input`, `directories`, `files`, `links` or `resources` | ignores it without a message, because these blocks are not checked for unknown keys. With `antismash` instead of `antismash_db` in `directories`, for example, the antiSMASH database is downloaded |
| a misspelt top-level key | ignores it without a message, because top-level keys are not checked either. With a misspelt `parameter:` block, every key under `parameters` takes its default |

!!! warning "Outside `parameters`, a misspelt optional key is ignored without a message"

    A misspelt key or a refused value under `parameters` stops the run before any job starts. A
    misspelt optional key in another block is ignored, and its default is used. A misspelt
    `parameters` or `resources` block is ignored too, and all its keys take their defaults. Check
    those names yourself.

YAML reads unquoted `off`, `no` and `false` as false, and `on`, `yes` and `true` as true, before
FunFlux sees them. What a key then does depends on the key:

| Keys | Unquoted `off`, `no` or `false` | Empty (`key:` or `key: ""`) |
|---|---|---|
| `mode` of `decontamination`, `mitochondria`, `genemark`, `tiberius` and `signalp6`, and `amr.chroquetas` | read as `off`, so `decontamination.mode: off` works without quotes | the default of the key (`off` for `amr.chroquetas`) |
| `amr.mito_genetic_code` | refused, and so is true | `auto` |
| `nt_version`, `masking_method`, `busco.lineage`, `prediction.busco_db`, `assembler`, `flye_input_mode` | read as empty, so the default applies: `masking_method: off` still masks with tantan | the default |

In a key that has an `off` value, `off` needs no quotes. The six keys of the last row
have no such value, and an unquoted `off` gives their default. `tiberius.mode` and
`amr.chroquetas` also read an unquoted `on`, `yes` or `true` as `on`.

A relative path is read from the launch folder, not from the folder of the configuration file and
not from `output_dir`:

| Keys | A relative path |
|---|---|
| `output_dir`, the `directories.*_db` paths, `files.iprscan`, `files.annotation_params`, the protein files of the [annotation parameter table](#annotation-parameter-table) | is made absolute as FunFlux reads it |
| `directories.genemark_dir`, `files.signalp6_archive` | stops the run when GeneMark or SignalP 6 is on, because these must already be absolute |

Write absolute paths unless you know the launch folder in advance.

The Modes column of the tables below names the modes that read a key, and the other modes ignore
it. "screen" stands for the contamination and mitochondrion screen, which runs in `illumina`,
`longread` and `hybrid` mode, and in `contigs` mode with `parameters.contigs.screen: true`.

## What every run needs

These are the keys to fill in. The `links` come filled in, and every key under `parameters` has a
default.

| Key | Set it to |
|---|---|
| `mode` | `illumina`, `longread`, `hybrid` or `contigs` ([Choosing a mode](../getting-started/choosing-a-mode.md)) |
| `input.illumina_dir`, `input.longread_dir` or `input.contigs_dir` | the input folder of the mode, and in `hybrid` mode the first two |
| `directories.output_dir` | the folder for results and logs |
| `directories.eggnog_db`, `directories.funannotate_db` | the eggNOG-mapper data, and the folder for the funannotate database ([Reference databases](../getting-started/databases.md)) |
| `directories.blast_db` | the NCBI BLAST database folder, when the screen runs |
| `files.annotation_params` | the [annotation parameter table](#annotation-parameter-table) |
| `files.iprscan` | the InterProScan launcher ([Licensed and external tools](../getting-started/external-tools.md)) |
| `resources.threads` | the number you give to `--cores` |

## `mode`

The top-level key `mode` selects the input mode of the run.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `mode` | required; `illumina` in the template | `illumina`, `longread`, `hybrid`, `contigs` (case and surrounding spaces are ignored) | all | Selects the input type, how the reads are cleaned and assembled, and the input folder that is searched for samples. The contamination and mitochondrion screen and every later step except polishing are the same in all modes. |

A configuration without `mode` is taken for a FunFlux 1.x, FunFluxL or FunFluxL+ configuration, and
the run stops with a message that points to [Coming from v1](../getting-started/from-v1.md).

## `input`

The `input` block names the folders that FunFlux searches for samples. A run reads only the keys
of its mode.

```yaml
input:
  illumina_dir: ""             # illumina, hybrid
  longread_dir: ""             # longread, hybrid
  contigs_dir:  ""             # contigs
```

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `illumina_dir` | required in its modes | an existing folder holding `{sample}_R1.<ext>` and `{sample}_R2.<ext>` | `illumina`, `hybrid` | Paired-end reads. Samples are found from the R1 files, and an R1 file without its R2 stops the run. |
| `longread_dir` | required in its modes | an existing folder holding `{sample}_ont.<ext>` or `{sample}_pacbio.<ext>` | `longread`, `hybrid` | Long reads. The suffix sets the technology of each sample, `_ont` for Oxford Nanopore and `_pacbio` for PacBio HiFi (CCS). One batch may mix both, but a sample with both files stops the run. In `hybrid` mode it may be the same folder as `illumina_dir`. |
| `contigs_dir` | required in its mode | an existing folder holding `{sample}.fasta`, `{sample}.fa` or `{sample}.fna` | `contigs` | One assembly per sample. |

- Read files end in `fastq`, `fq`, `fastq.gz` or `fq.gz`. All files of one kind in a batch
  (Illumina reads, long reads, assemblies) must share one extension, and two extensions stop the
  run.
- A `_pacbio` file is taken to hold HiFi reads, and a NOTE says so. PacBio CLR reads are not
  supported.
- In `hybrid` mode every sample needs Illumina and long reads. The run stops and lists the samples
  that lack either.
- In `contigs_dir` every file whose name contains a dot is taken as a sample, and the last dot
  splits off the extension (`my.genome.fasta` is the sample `my.genome`). Keep nothing else there.
  A `README.txt` stops the run with the two-extensions error, and a FASTA in a subfolder becomes a
  sample name containing `/`, which is refused.

Sample names become file names, report labels and the locus-tag prefix of every gene
(`funannotate predict --name`):

| Sample name | Result |
|---|---|
| contains a space or one of `* # @ % ^ / ! ? & : ; | < >` | stops the run |
| contains underscores or dots | allowed |
| is not 3 to 12 letters or digits starting with a letter | runs, with a NOTE that its gene IDs cannot be submitted to GenBank as they are |

A name such as `SAMPLE1`, 3 to 12 letters or digits starting with a letter, gets neither
the stop nor the NOTE.

## `directories`

The `directories` block names the output folder and the folders of the databases and of
GeneMark-ES.

```yaml
directories:
  output_dir:     ""
  blast_db:       ""
  eggnog_db:      ""
  funannotate_db: ""
  genemark_dir:   ""
  antismash_db:   ""
  unite_db:       ""
```

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `output_dir` | required | any path | all | Root of every result and log, created by the first job. Later runs with the same `output_dir` reuse the databases and the SignalP 6 installation stored in `00.resources/`. Two things are written elsewhere: the funannotate database, in `directories.funannotate_db`, and Snakemake's `.snakemake/` folder, in the launch folder. The `.snakemake/` folder holds Snakemake's records of the run, the run logs and, unless `--conda-prefix` is given, the conda environments. |
| `blast_db` | required with the screen | a folder | screen | The NCBI nucleotide BLAST database. The folder holds the files named after `parameters.nt_version` (for `core_nt`, `core_nt.nal` and the `core_nt.NN.*` volumes), the BLAST taxonomy files `taxdb.btd` and `taxdb.bti`, and the NCBI taxonomy files `nodes.dmp` and `names.dmp`. FunFlux prints a NOTE when the folder does not exist. Without the two `.dmp` files the dry run stops with Snakemake's `MissingInputException in rule blob_json`, and without `taxdb.*` the organism names in the BLAST table read `N/A`. |
| `eggnog_db` | required | a folder | all | The eggNOG-mapper data (`eggnog.db`, `eggnog_proteins.dmnd`), passed to `emapper.py --data_dir`. FunFlux prints a NOTE when the folder does not exist. |
| `funannotate_db` | required | a folder, which may not exist yet | all | The funannotate database. When the folder holds no `funannotate-db-info.txt`, the archive named in `links.funannotate_link` is downloaded into it, checked against its `.sha256` file and unpacked. The absolute paths stored inside the database are rewritten to this location once per output folder. |
| `genemark_dir` | empty: no GeneMark-ES | absolute path of the folder holding an executable `gmes_petap.pl` | all | With `parameters.genemark.mode: auto`, funannotate predict runs GeneMark-ES from here. The licence key `~/.gm_key` of the user who runs Snakemake must then exist, or the run stops. |
| `antismash_db` | empty: download | a folder containing `clusterblast/` | all | An antiSMASH database you already hold. It is linked into `00.resources/antismash/databases/` and never written to. When empty, `download-antismash-databases` fills that folder (about 9.4 GB for antiSMASH 8.0.4). |
| `unite_db` | empty: download | a folder containing an uncompressed `unite_its_sintax.fasta` | all | A UNITE SINTAX reference you already hold. It is linked into `00.resources/unite/`, and `links.unite_its_link` is not used (a NOTE says so). |

`antismash_db` and `unite_db` are checked as soon as they are set. A missing folder, or a folder
without `clusterblast/` or `unite_its_sintax.fasta`, stops the run, dry runs included.

!!! warning "Start one run first when several runs share a funannotate database"

    Several runs can share one funannotate database, but nothing stops two of them from repairing
    its stored paths at the same time. Start one run, let it pass the step
    `normalize_funannotate_db`, then start the others, and reach the folder under the same path in
    every run ([Reference databases](../getting-started/databases.md#what-the-workflow-does)).

## `files`

The `files` block names three single files: the annotation parameter table, the InterProScan
launcher and the SignalP 6 archive.

```yaml
files:
  annotation_params: ""
  iprscan:           ""
  signalp6_archive:  ""
```

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `annotation_params` | required | path of a tab-separated table | all | The species name, protein evidence and AUGUSTUS model of each sample. See [Annotation parameter table](#annotation-parameter-table). |
| `iprscan` | required | path of `interproscan.sh` | all | The InterProScan 5 launcher, installed by hand (FunFlux v2 was tested with 5.78-109.0). It runs on the Java 11 of the workflow's `interproscan` environment. FunFlux prints a NOTE when the file does not exist. |
| `signalp6_archive` | empty | absolute path of a readable file | all, when `parameters.signalp6.mode` is not `off` | The SignalP 6 archive downloaded from DTU (for example `signalp-6.0i.fast.tar.gz`). It is installed once per output folder into `00.resources/signalp6/`, and the workflow never downloads SignalP. |

## `links`

The `links` block holds the download addresses of what the workflow fetches itself, and the
template fills them in.

```yaml
links:
  phix_link: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
  funannotate_link: https://zenodo.org/records/23265189/files/funannotate_db_1.8.17_2026-10-08.tar.gz?download=1
  unite_its_link: https://s3.hpc.ut.ee/plutof-public/original/338a1413-6039-4e00-b5cf-410346a1e366.gz
  tiberius_model_link: https://bioinf.uni-greifswald.de/bioinf/tiberius/models/fungi.tar.gz
  telomere_db_link: https://raw.githubusercontent.com/tolkit/a-telomeric-repeat-database/73c635c7b80e5ab369b4c7ec67d23db32dcaa1b8/data/telo_database_v2.tsv
  chroquetas_link: https://raw.githubusercontent.com/nmquijada/ChroQueTas/1aad051d6ad48f5e403636699d19095e655764f7
```

The workflow has no built-in links, so the run stops when a link that the mode needs is missing.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `phix_link` | template URL | a URL | `illumina`, `hybrid` (required there) | The PhiX174 genome (RefSeq GCF_000819615.1) that Illumina libraries are spiked with. Read pairs that align to it concordantly (both mates, facing each other) are removed; a pair with only one mate on PhiX is kept. |
| `funannotate_link` | template URL | a Zenodo URL of the form `https://zenodo.org/records/<id>/files/<archive>` | all, and required even when the database is already on disk | The archive of the funannotate database (Zenodo record 23265189, 2.1 GB), downloaded only while `directories.funannotate_db` holds no database. To switch to another version of the database, change the link and point `directories.funannotate_db` at a new, empty folder, for the reasons given under [Newer snapshots, and building your own](../getting-started/databases.md#build-snapshot). The checksum file, the archive name with `.sha256` in place of `.tar.gz`, is taken from the same record and checked with `sha256sum -c`. A mismatch stops the run before anything is unpacked. |
| `unite_its_link` | template URL | a URL of a gzip-compressed SINTAX FASTA | all, and required unless `directories.unite_db` is set | The UNITE reference for the ITS classification. The template points to UNITE's SINTAX release for eukaryotes of 19.02.2025 (266,589 sequences). It is downloaded once per output folder and unpacked on the fly into `00.resources/unite/unite_its_sintax.fasta`, and `PROVENANCE.txt` beside it records the URL, date, md5 and size. No checksum is checked. |
| `tiberius_model_link` | template URL | a URL of a `.tar.gz` archive holding one folder with the model files | all, and required unless `parameters.tiberius.mode` is `off` | The `fungi` model of the Tiberius gene finder, from its authors' server (about 30 MB). It is downloaded once per output folder, checked against [`parameters.tiberius.model_sha256`](#tiberius) and unpacked into `00.resources/tiberius/model/`. A mismatch stops the run before anything is unpacked. |
| `telomere_db_link` | template URL | the URL of the tab-separated database file | all, and required | The telomere repeat database of tidk, `telo_database_v2.tsv` (0.6 MB), with the known telomere repeat units that the contig summary checks its candidates against. It is downloaded once per output folder into `00.resources/telomeres/` and checked against [`parameters.contig_summary.telomere_db_sha256`](#contig_summary), and a mismatch stops the run. The template address names one commit of the repository `tolkit/a-telomeric-repeat-database`, so the file behind it does not change ([Reference databases](../getting-started/databases.md#telomere-db)). |
| `chroquetas_link` | template URL | the raw-file address of one commit of the ChroQueTas repository on GitHub, without a file name | all, required when [`parameters.amr.chroquetas`](#amr) is `on` and not read when it is `off` | ChroQueTas 1.0.1 and its FungAMR database. Rule `chroquetas_db` downloads `bin/ChroQueTas.sh`, `FungAMR_db/db-v20250811.tgz` and `LICENSE` from this address once per output folder into `00.resources/chroquetas/`. It checks the script and the archive against `parameters.amr.chroquetas_script_sha256` and `chroquetas_db_sha256`, and a mismatch stops the run before anything runs or is unpacked ([Reference databases](../getting-started/databases.md#chroquetas)). |

## `resources`

The `resources` block holds three numbers: the threads, the memory in GB and, for a cluster, the
time limit in hours.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `threads` | 24 | an integer ≥ 1 | all | The most threads any rule asks for (see below). On a single computer, give the same number to `--cores` and do not add `--jobs`/`-j`. For a local run `--jobs` is another name for `--cores`, and giving both can let the jobs use more CPUs than the machine has. On a cluster each job asks the scheduler for its rule's threads, so use at most the CPUs of one compute node (Slurm: `sinfo -o "%n %c"`). |
| `ram_gb` | 128 | an integer ≥ 1 | all | On a single computer, the memory limit SPAdes is given (`-m`, step `illumina_assembly`) and, capped at 64 GB, the Java heap of Qualimap (`map_evaluation`, and in `hybrid` also `map_guide_evaluation`). Snakemake reserves no memory for these steps unless `--resources` sets a limit, `mem_gb` for SPAdes and `java_mem` for Qualimap ([Running FunFlux](running.md#cpu-cores-and-threads)). On a cluster, it is also the memory every heavy step asks the scheduler for, in every mode, and Qualimap asks for its heap plus 4 GB, never more than `ram_gb`. Use at most the memory the scheduler offers one compute node, or no node can take those steps (Slurm: `sinfo -o "%n %m"`, in MB). See [How the settings reach Slurm](running.md#cluster-resources). A change re-runs nothing. |
| `runtime_hours` | 48 | an integer ≥ 1 | all | On a cluster, the time limit that every heavy and every long step asks the scheduler for. These are the steps that can run for more than an hour, such as the assemblers, BLAST, gene prediction and InterProScan, and a step still running at the limit is stopped. Use at most the time limit of the partition, or the scheduler refuses those steps or never starts them (Slurm: `sinfo -o "%P %l"`). A single computer ignores it. A change re-runs nothing. |

The rules ask for threads in four ways:

| Rules | Threads they ask for |
|---|---|
| `funannotate_prediction`, `tiberius_prediction`, the assemblers (`illumina_assembly`, `long_read_assembly`), the read mappings (`map_phix`, `map_contigs`, `map_guide_contigs`, `map_sel_contigs`) | `threads` |
| the other CPU-bound rules | `threads`, up to a fixed cap between 2 and 24 |
| SignalP 6 (`signalp6_annotation`) | `parameters.signalp6.torch_num_threads`, lowered to `threads` when it is larger |
| every other rule | 1 |

`threads` is the most any one job takes, so it is the number to give `--cores`. When `--cores`
is lower than a rule's number, Snakemake runs the job with `--cores` threads. The cap of every
rule, with the reasons, is listed on [CPU cores and threads](running.md#cpu-cores-and-threads). On a
cluster the three keys must each fit one compute node and one partition
([On a cluster](running.md#on-a-cluster)).

## `parameters`

The `parameters` block holds the analysis settings. FunFlux checks every key under `parameters` in
every mode, so a wrong value stops the run even in a mode that does not use the key. There are
three exceptions: `spades_kmers` is read only in `illumina` and `hybrid` mode, the `longread` and
`hybrid` blocks only in their own mode, and the files that the `decontamination` keys name only
when the screen runs. After a change of most settings, the next run re-runs the step that uses it
and every step after it.

### BLAST database and hit depth

These keys choose the BLAST database of the screen and how many hits BLAST keeps for each contig.

```yaml
parameters:
  nt_version: core_nt
  blast_max_target_seqs: 50
  blast_max_hsps: 5
```

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `nt_version` | `core_nt` | the name of the BLAST database in `directories.blast_db`, usually `core_nt` or `nt`, not checked | screen | The database BLAST searches, as `blastn -db <blast_db>/<nt_version>`. |
| `blast_max_target_seqs` | 50 | an integer ≥ 1, with a NOTE below 5 | screen | `-max_target_seqs`, the number of database sequences kept per contig. BLAST applies it during the search, not as a final best-N filter, so a small value can lose better matches. |
| `blast_max_hsps` | 5 | an integer ≥ 1, with a NOTE below 5 while the mitochondrion screen is on | screen | `-max_hsps`, the number of alignments kept per database sequence. The mitochondrion screen reads `qcovs`, the share of the contig covered by one database sequence, which BLAST computes from the kept alignments only. |

One BLAST table feeds both BlobTools (the genus of each contig) and the mitochondrion screen. With
`blast_max_hsps: 1` BLAST keeps only the single best local alignment per database sequence, so a
complete mitochondrial genome that aligns to its closest reference in several blocks can look only
half covered. Setting `blast_max_target_seqs: 10` and `blast_max_hsps: 1` gives the BLAST settings
of FunFlux 1.1.x. These values also change the hits whose bit scores BlobTools adds up, and so can change the genus it gives a contig. See
[BLAST hit depth](../methods/blast-hit-depth.md).

### `spades_kmers`

`spades_kmers` lists the k-mer sizes that SPAdes assembles the Illumina reads with.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `spades_kmers` | `"21,33,55,77,99,127"` | `auto`, or comma-separated text or a YAML list, with every k odd, from 1 to 127, and strictly ascending | `illumina`, `hybrid` | The SPAdes `-k` list. `auto` lets SPAdes choose from the read length ([Illumina mode](../modes/illumina.md)). |

### `decontamination`

These keys tell the contig selector which contigs to keep or remove by genus.

```yaml
parameters:
  decontamination:
    mode: off                  # off | include | exclude
    discard_no_hit: false
    include_genera:
    include_genera_by_sample:
    exclude_genera:
    exclude_genera_file:
    sample_overrides:
```

The contig selector, rule `select_contigs`, takes the genus BlobTools assigned to each contig (the
taxon with the highest summed BLAST bit score) and keeps or removes the contig. With every policy,
`off` included, it writes `contaminants/contig_taxonomy_decisions.tsv`, one line per contig with
its genus, `keep` or `remove`, and the reason. With `mitochondria.mode: separate_confirmed` (the
default), the contig selector always keeps the contigs that the mitochondrion screen confirmed (reason
`protected_confirmed_mitochondrial`). With `report` or `off`, nothing is protected.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `mode` | `off` | `off`, `include`, `exclude`, while `auto` stops the run (see the warning box below) | screen | `off` keeps every contig. `include` keeps only the listed genera. Name the genera that `contig_taxonomy_decisions.tsv` gives the isolate's own contigs, which for a poorly sequenced fungus are often not its true genus. When the screen runs, a sample in `include` mode without a genus stops the run before any job starts. `exclude` removes the listed genera. |
| `discard_no_hit` | `false` | `true`, `false` (also `yes`/`no`, `on`/`off`, `1`/`0`) | screen | With `exclude`, also removes the contigs without a BLAST hit (reason `discarded_no_hit`). `include` removes them anyway, because no-hit is not on its list, and the setting changes only the reason. No effect with `off`. |
| `include_genera` | empty | a genus, a YAML list, or names separated by `;` or `,` | screen | Genera that `include` keeps. |
| `include_genera_by_sample` | empty | path of a tab-separated table with the header columns `sample` and `genus`, in which a sample may have several rows | screen | Adds the sample's genera to `include_genera`. |
| `exclude_genera` | empty | as `include_genera` | screen | Genera that `exclude` removes. |
| `exclude_genera_file` | empty | path of a text file with one genus per line, in which lines starting with `#` are skipped | screen | Adds its genera to `exclude_genera`. |
| `sample_overrides` | empty | path of a tab-separated table with the columns `sample` and `mode`, and optionally `include_genera`, `exclude_genera`, `discard_no_hit` | screen | Replaces the named settings for single samples, and an empty cell keeps the run-wide value. A `mode` cell holds `off`, `include`, `exclude` or nothing, and a `discard_no_hit` cell a yes/no word or nothing. |

Genus names are compared without regard to case. `include` also accepts a few aliases, from a table
of bacterial genus splits and from the prefixes brady, meso, neo, sino, aeri, caldi and geo.
*Neocosmospora* and *Cosmospora*, for example, match each other. An alias can only keep more
contigs. `exclude` compares exact names.

When the screen runs, the three files are checked before any job starts. Each must exist and be
UTF-8 text:

| File | Must also hold |
|---|---|
| `include_genera_by_sample` | a header line with the columns `sample` and `genus` |
| `exclude_genera_file` | no further check |
| `sample_overrides` | a header line with the columns `sample` and `mode`. Every `mode` cell holds `off`, `include`, `exclude` or nothing, and the `discard_no_hit` cell of each sample of the batch a yes/no word or nothing |

Every sample of the batch that ends up in `include` mode needs at least one genus. The log
`logs/select_contigs_{sample}.log` warns when more than 20 % of the bases are removed, and when two
or more genera each hold at least 5 % of the bases.

The v1 key `parameters.genus` is still read when the `decontamination` block is missing or empty,
and a NOTE reports the translation:

| `decontamination` block | v1 key `parameters.genus` | Becomes |
|---|---|---|
| missing or empty | a genus | `mode: include` with that genus |
| missing or empty | empty | `mode: off`, so every contig is kept |
| not empty | any value | nothing, and `genus` is ignored without a message |

Move the genus into the `decontamination` block, because once that block holds anything, `genus`
has no effect.

!!! warning "There is no `auto` mode, because BLAST genus labels are unreliable for poorly sequenced fungi"

    The v1 workflows had a `mode: auto` that kept the genus BlobTools gave the most contigs.
    FunFlux v2 refuses it before any job starts and, when the screen runs, also in the `mode`
    column of `sample_overrides`.

    | Why `auto` fails | What happens |
    |---|---|
    | it counts contigs, not bases | when the contigs without a BLAST hit outnumber those of any genus, it keeps them and removes the rest of the genome |
    | the isolate's genus is poorly represented in core_nt | BLAST labels much of the isolate's own DNA with related genera, so `auto` keeps a wrong genus, and even `include` with the true genus can remove most of the genome |

    `auto` can therefore remove most of the genome. Run with `off` (the default), read
    `contig_taxonomy_decisions.tsv` and `{sample}_composition.txt`, then `exclude` the contaminant
    genera or `include` the genera BLAST gave the isolate's own contigs
    ([Why there is no auto mode](../analysis/decontamination.md#why-there-is-no-auto-mode)).

### `mitochondria`

These keys set how the mitochondrion screen recognises the mitochondrial genome, and what happens
to it.

```yaml
parameters:
  mitochondria:
    mode: separate_confirmed   # off | report | separate_confirmed
    min_direct_blast_qcov: 35
    min_assembler_coverage_ratio: 3.0
    very_strong_qcov: 70
```

A direct hit is a BLAST hit whose title names a mitochondrial genome, and `qcovs` is the
percentage of the contig that its best direct hit covers.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `mode` | `separate_confirmed` | `off`, `report`, `separate_confirmed` | screen | `separate_confirmed` writes `mitochondria/mitochondrial_screen.tsv`, protects the confirmed contigs from the contig selector, and writes them to `contigs_mitochondrial.fasta` instead of `contigs_final.fasta`, the nuclear genome that QC, ITS and gene prediction read. `report` writes the table only, so nothing is protected or split off. `off` classifies nothing. |
| `min_direct_blast_qcov` | 35 | a number above 0 and at most 100 | screen | Without assembler evidence (`illumina` and `contigs` mode, and the `hybrid` guide), a direct hit with `qcovs` at or above this value confirms the contig, and a weaker one makes it a `mitochondrial_candidate`. With assembler evidence, it is the minimum that depth or circularity must back up. |
| `min_assembler_coverage_ratio` | 3.0 | a number above 0 | `longread`, `hybrid` | Contig depth divided by the median contig depth, as reported by the assembler (Flye depth, hifiasm `rd:i`). A direct hit at or above `min_direct_blast_qcov`, plus a ratio at or above this value or a clean circle (circular, not a repeat, not an alternative copy), confirms the contig. |
| `very_strong_qcov` | 70 | a number above 0 and at most 100 | `longread`, `hybrid` | A direct hit at or above this value confirms the contig without depth or circularity. |

The assembler evidence comes from the hifiasm or Flye assembly of `longread` and `hybrid` mode.
`illumina` and `contigs` mode, and the guide screen of `hybrid`, give the screen none, so the last
two keys change nothing there.
Neither circularity nor high depth confirms a contig without a direct hit. See
[Mitochondria](../analysis/mitochondria.md).

### `busco`

The `busco` block chooses the lineage dataset of the BUSCO completeness check.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `lineage` | `auto` | `auto`, or a BUSCO dataset name ending in `_odb` and a number, such as `hypocreales_odb12` (case is ignored) | all | `auto` runs `busco --auto-lineage-euk`. BUSCO then reports a generic (eukaryota) and a specific lineage, and MultiQC shows two BUSCO rows per sample. A name runs `busco -l <name>`, which gives one row and the same dataset for every sample. |

Where the lineage data come from depends on the setting:

| Setting | Where the lineage data come from | Internet needed |
|---|---|---|
| `auto` | downloaded again for every sample in every run, into a temporary folder of that sample | in every BUSCO job |
| a name | downloaded once per output folder by rule `busco_data` into `00.resources/busco/`, and read by every BUSCO job with `--offline`. Later runs into the same folder reuse it | once, on the computer that runs Snakemake |

!!! note "A named lineage is the one choice for compute nodes without internet"

    A named lineage is downloaded once, on the computer that runs Snakemake, while `auto` needs
    internet in every BUSCO job.

The dry run checks only the form of the name. A name BUSCO does not offer stops `busco_data` with a
`[FunFlux] busco_data:` message. A dataset archive takes tens to hundreds of MB (274.6 MB for
`hypocreales_odb12`), and with `auto` the BUSCO folder of one sample holds several hundred MB. See
[Assembly QC](../analysis/assembly-qc.md).

### `contig_summary`

These keys steer the contig summary, and above all how it finds the telomeres.

```yaml
parameters:
  contig_summary:
    telomere_motif: auto
    telomere_db_sha256: 1d4f5c8c2cce6a76aaa8de01602317a8e0b8d0401711b96f3c7a826c385da699
    telomere_window: 25000
    min_telomere_copies: 5
    min_rdna_units: 3
```

After the annotation, rule `contig_summary` writes `02.assembly/{sample}/eval/contig_summary.tsv`.
It has one row per contig of the draft assembly, with the evidence the run produced about it and a
verdict. The evidence consists of the telomere arrays at the contig ends, read depth, the
assembler's depth and circularity flag, rDNA units, BUSCO markers, predicted genes and the calls of
the contamination and mitochondrion screen. The verdict is one of `chromosome`, `chromosome_arm`,
`mitochondrion`, `rDNA_array`, `low_depth`, `repeat`, `nuclear_contig`, `nuclear_fragment`,
`contaminant` or `unclassified`.

The telomere arrays are the one piece of evidence worked out from the sequence itself. Their repeat
unit differs between fungal groups, so rule `telomere_units` chooses it from the assembly with tidk.
It checks the unit against tidk's telomere repeat database, which rule `telomere_db` downloads once
per output folder ([Telomeres](../analysis/contig-summary.md#telomeres)). The first four keys steer
that choice.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `telomere_motif` | `auto` | `auto` (any case), or letters A, C, G and T only, 4 to 30 long, with lower case upper-cased | all | With `auto`, the telomere repeat unit is chosen from the data. Every tandem repeat unit found at the contig ends, every unit of the tidk database and every unit of the literature is scored by how it behaves at the contig ends. The best-supported one is the unit of the assembly. A unit you set is used as given, whatever the data say, and the data's own choice is reported next to it, with a NOTE when it differs. The known fungal units are listed under [Which fungi carry which repeat](../analysis/contig-summary.md#which-fungi-carry-which-repeat). Either strand and any rotation work, so `CCCTAA` gives the same result as `TTAGGG`. |
| `telomere_db_sha256` | the checksum of the database file at the template address | 64 hexadecimal characters (capitals are accepted) | all | The SHA-256 that the file downloaded from `links.telomere_db_link` must have, checked by rule `telomere_db` before it is used. A mismatch stops the run. The template address names one commit of the database repository, so the file behind it does not change. To use another version, change the link and this checksum together (`sha256sum <file>` prints it). |
| `telomere_window` | 25000 | an integer from 500 to 100000 | all | The terminal window, in bp. It is the stretch at each contig end in which the candidate units are collected, and how far from the end `tidk ends` accepts an array. The telomere array can sit some kb inside a contig end when the assembler ran past the telomere, so the window must cover that. The two offset columns of the table show where each array lies. |
| `min_telomere_copies` | 5 | an integer ≥ 2 | all | The copies of the unit an array needs to count as a telomere. `tidk ends` requires this many unit lengths (30 bp for a 6-bp unit), and a unit is collected at the contig ends only when it forms a tandem run of this many copies. |
| `min_rdna_units` | 3 | an integer ≥ 1 | all | The number of rDNA units (the ITS windows of `its_windows`) from which a contig can be an `rDNA_array`. |

A `parameters.contig_summary.telomere_motif` that is neither `auto` nor a plain unit stops the run
with
`[FunFlux] parameters.contig_summary.telomere_motif: must be auto (the unit is chosen from the data) or the telomere repeat unit of the species, 4 to 30 letters A, C, G and T, such as TTAGGG (got: <value>)`.
A malformed checksum stops it with
`[FunFlux] parameters.contig_summary.telomere_db_sha256: must be the SHA-256 checksum of the telomere database file, 64 hexadecimal characters (got: <value>)`.
The three numbers are checked as whole numbers against the limits above.

Every dry run prints the settings in one line:
`Contig summary: telomere motif auto (the repeat unit is chosen from the data), a telomere from 5 tandem copies within 25000 bp of a contig end, an rDNA array from 3 rDNA unit(s).`
With a unit set, the line reads `telomere motif TTAGGG (used as given; the data's own choice is reported next to it)`.

A change of these settings re-runs the following steps:

| Change | Re-runs |
|---|---|
| `telomere_motif`, `telomere_window`, `min_telomere_copies` | `telomere_units` and `contig_summary`, then MultiQC |
| `min_rdna_units` | `contig_summary` only |
| `links.telomere_db_link` or `telomere_db_sha256` | `telomere_db`, then `telomere_units` of every sample |

These keys re-run only the telomere and contig-summary steps, never the assembly or the
annotation. hifiasm has a separate key of its own, `parameters.longread.hifiasm.telomere_motif` or
`parameters.hybrid.hifiasm.telomere_motif`. It takes the unit as it reads at the 5' end of a
chromosome (`CCCTAA` for `TTAGGG`). See [Contig summary](../analysis/contig-summary.md).

### `its` and `its_taxonomy_cutoff`

These keys steer how the ITS is extracted from the genome and classified against UNITE.

```yaml
parameters:
  its:
    window_flank: 2000
    max_full_length: 1500
    unite_identified_only: true
  its_taxonomy_cutoff: 0.8
```

nhmmer finds every 5.8S gene in `contigs_final.fasta` with ITSx's own profiles. A window is cut
around each, and windows on the minus strand are reverse-complemented. ITSx extracts ITS1, 5.8S and
ITS2 from the windows, and SINTAX classifies ITS1 and ITS2 against UNITE.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `its.window_flank` | 2000 | an integer ≥ 200 | all | The bases kept on each side of the middle of each 5.8S site. 2,000 bp hold the end of the SSU, ITS1, ITS2 and the start of the LSU, which ITSx needs to recognise the unit. |
| `its.max_full_length` | 1500 | an integer ≥ 400 | all | The longest ITS1-5.8S-ITS2 record listed as `ok` in `{sample}_ITS_copies.txt`. A longer one is listed as `rejected_too_long` (1,500 bp is ITSx's own "suspiciously long" limit). Only ITS1 and ITS2 are classified, so this key changes that table only. |
| `its.unite_identified_only` | `true` | `true` or `false` | all | With `true`, rule `unite_identified` keeps the UNITE records whose taxonomy names a genus (`g:` followed by a name) in `00.resources/unite/unite_its_sintax.identified.fasta`, and SINTAX classifies against that file. About a third of the records of the default release stop above the genus, and one that is identical to the query can hold the result at `d:Fungi`. With `false`, SINTAX classifies against the whole release, and rule `unite_identified` does not run. See [The UNITE reference filter](../analysis/its-taxonomy.md#the-unite-reference-filter). |
| `its_taxonomy_cutoff` | 0.8 | a number above 0 and at most 1 | all | Passed as `vsearch --sintax_cutoff`. `{sample}_ITS1_taxonomy.txt` and `{sample}_ITS2_taxonomy.txt` hold every rank with its bootstrap confidence in column 2, and in column 4 only the ranks whose confidence reaches the cutoff. |

ITS1 and ITS2 sequences outside 50 to 600 bp are listed in `{sample}_ITS_copies.txt` and not
classified. That range is fixed. See [ITS taxonomy](../analysis/its-taxonomy.md) and
[ITS on long-read contigs](../methods/its-long-read-contigs.md).

### Repeat masking

Before gene prediction, the genome is soft-masked by tantan, or by RepeatModeler and RepeatMasker.

```yaml
parameters:
  masking_method: tantan       # tantan | repeatmodeler_repeatmasker
  repeatmodeler_quick: true
  repeatmodeler_ltrstruct: false
  repeatmodeler_seed:
```

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `masking_method` | `tantan` | `tantan`, `repeatmodeler_repeatmasker` | all | How `contigs_sort.fasta` is soft-masked (repeats in lower case) before gene prediction. `tantan` masks tandem repeats and low-complexity sequence. `repeatmodeler_repeatmasker` builds a repeat library from this genome with RepeatModeler and masks with RepeatMasker, which also finds interspersed repeats such as transposons. This method takes hours. |
| `repeatmodeler_quick` | `true` | `true`, `false` | all (RepeatModeler only) | `RepeatModeler -quick`, a shorter run that finds fewer repeat families. |
| `repeatmodeler_ltrstruct` | `false` | `true`, `false` | all (RepeatModeler only) | `RepeatModeler -LTRStruct`, which adds the LTR retrotransposon structural module. |
| `repeatmodeler_seed` | empty | empty, or an integer ≥ 0 | all (RepeatModeler only) | `RepeatModeler -srand <seed>`. When empty, RepeatModeler seeds itself from the clock, and two runs on the same genome build slightly different libraries. |

RepeatMasker needs no Dfam or RepBase library here, because it masks with the library
RepeatModeler built. See [Repeat masking](../analysis/repeat-masking.md).

### `prediction`

The `prediction` keys pass two settings to `funannotate predict`.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `busco_db` | `dikarya` | letters, digits and `_`, and once the database is on disk, a folder of that name must exist in it | all | `funannotate predict --busco_db`, the conserved-gene set that SNAP and GlimmerHMM are trained on, and AUGUSTUS too for the samples without a model in column 4 of the annotation table. The funannotate database as downloaded holds only `dikarya` (BUSCO `dikarya_odb9`, 1,312 genes). For a fungus outside the Dikarya, add a set such as `fungi` ([Other BUSCO sets for gene prediction](../getting-started/databases.md#busco-sets)) and set `busco_db: fungi`. |
| `augustus_optimize` | `false` | `true`, `false` | all | `funannotate predict --optimize_augustus`, for the samples whose column 4 of the annotation table is empty. After the quick training, `optimize_augustus.pl` tunes AUGUSTUS's meta-parameters, the settings the quick training keeps fixed, by 8-fold cross-validation: each eighth of the BUSCO test genes is predicted by a model trained without it. AUGUSTUS is then trained again (see [The AUGUSTUS model](../analysis/gene-prediction.md#augustus-model)). It adds about a quarter of an hour for a 40 Mb genome and uses at most 8 CPUs. It does not apply to samples that name a stored AUGUSTUS model, and a NOTE names them. Switching it re-runs gene prediction only for the samples whose column 4 is empty. |

### `genemark`

The `genemark` block decides whether funannotate predict uses GeneMark-ES.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `mode` | `auto` | `auto`, `off` (also `false`, `no`, `0`, `none`, and `null` in quotes). An empty value (`mode:`, `mode: ""` or an unquoted `null`) means `auto`, as for the other mode keys | all | With `auto`, funannotate predict uses GeneMark-ES when `directories.genemark_dir` is set. With `off` it never does, and funannotate runs without GeneMark. |

With GeneMark on, the run stops before any job starts when `genemark_dir` is not absolute, holds no
executable `gmes_petap.pl`, or `~/.gm_key` is missing. The step `genemark_check` also tests the
licence key before gene prediction
([The licence check before gene prediction](../getting-started/external-tools.md#the-licence-check-before-gene-prediction)).
FunFlux v2 was tested with GeneMark-ES on, and with `mode: off` on one public genome
([Validation on public genomes](../about/validation.md#gene-finders)).

### `tiberius`

The `tiberius` block switches the Tiberius gene finder on or off and pins its model.

```yaml
parameters:
  tiberius:
    mode: on                   # on | off
    model_sha256: 4a309a630574f395d82d747dcf64716704dab056db534e427ae0cec1b7410e20
```

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `mode` | `on` | `on` (also `true`, `yes`, `1`, and an empty value), `off` (also `false`, `no`, `0`, `none`, and `null` in quotes) | all | With `on`, Tiberius predicts genes (step `tiberius_prediction`) before `funannotate predict`, which adds them to its evidence (`--other_gff <file>:1`). With `off`, Tiberius does not run, its model is not downloaded, and `links.tiberius_model_link` is not read. Tiberius is independent of GeneMark-ES, and both can be on. |
| `model_sha256` | the checksum of the `fungi` model of 8 February 2026 | 64 hexadecimal characters (capitals are accepted) | Tiberius on | The SHA-256 that the file at `links.tiberius_model_link` must have, checked before it is unpacked. The address carries no version number. When its authors replace the model, the run stops until you accept the new file by putting its checksum (`sha256sum <archive>`) here. |

Tiberius needs about 16 GB of memory, so switch it off on a computer with less free memory. What it
does, and how its genes are combined with the others, is described under
[Tiberius](../analysis/gene-prediction.md#tiberius).

### `signalp6`

The `signalp6` block sets whether and how SignalP 6 runs.

```yaml
parameters:
  signalp6:
    mode: off                  # off | fast | slow | slow-sequential
    archive_sha256:
    torch_num_threads: 8
    batch_size: 10
```

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `mode` | `off` | `off`, `fast`, `slow`, `slow-sequential` (`off` also as empty, `false`, `no`, `0`, `none` or `null`) | all | Any value but `off` installs SignalP 6 from `files.signalp6_archive` once per output folder, predicts the signal peptides of the predicted proteins, and passes them to `funannotate annotate --signalp`. The archive must hold the model of the mode (table below). |
| `archive_sha256` | empty | 64 hexadecimal characters | SignalP on | Compared with the SHA-256 of the archive before installing, by the step that installs SignalP 6 (`install_signalp6`) and not by the dry run. `sha256sum <archive>` prints it. |
| `torch_num_threads` | 8, or `resources.threads` when that is smaller | an integer ≥ 1. A value above `resources.threads` is lowered to it, with a NOTE | SignalP on | The threads of `signalp6_annotation`, passed as `--torch_num_threads`. |
| `batch_size` | 10 | an integer ≥ 1 | SignalP on | `--bsize`, the number of proteins that SignalP 6 processes together. Raise it only after checking memory use. |

Each mode needs its own model in the archive:

| `mode` | The archive must hold |
|---|---|
| `fast` | `distilled_model_signalp6.pt` |
| `slow` | `ensemble_model_signalp6.pt` |
| `slow-sequential` | `sequential_models_signalp6` |

Choose the mode that matches the archive you downloaded, because a fast archive works only with
`fast`. The speed and memory of the three modes, from the SignalP 6 README, are given under
[SignalP 6](../getting-started/external-tools.md#signalp-6). FunFlux v2 was tested with `fast`.
See [Functional annotation](../analysis/functional-annotation.md).

### `amr`

The `amr` block switches the AMR screen (the resistance-mutation screen) on or off, and sets the
genetic code of its mitochondrial run and the checksums of its downloads.

```yaml
parameters:
  amr:
    chroquetas: off            # off | on
    mito_genetic_code: auto
    chroquetas_script_sha256: d8961a0d311be57d4a9939bd9e4d90dbc796589b3298dd1835359ed776abfce3
    chroquetas_db_sha256: 6f89ecb3117110aad85b8b6ca7b8a3d3715cb5984b6cb75450a6cda5633d16e7
```

In the AMR screen, rule `amr_chroquetas` runs ChroQueTas on every genome. For each sample it writes
`05.amr/{sample}/amr_summary.tsv`, with one row per target protein of the species, each with the
mutations of the FungAMR database found in it and a verdict on resistance. It also writes `amr_positions.tsv`,
every position compared. The species is column 2 of the
[annotation table](#annotation-parameter-table).

!!! note "Switching the AMR screen on accepts a non-commercial licence"

    ChroQueTas and its database are licensed CC BY-NC-ND 4.0, non-commercial use only. FunFlux
    downloads them at run time and does not ship them, and `on` is your acceptance of that licence
    ([Antifungal resistance mutations](../analysis/amr.md)).

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `chroquetas` | `off` | `off` (also empty, `false`, `no`, `0`, `none`, `null`), `on` (also `true`, `yes`, `1`) | all | With `on`, ChroQueTas and its database are downloaded once per output folder (step `chroquetas_db`), every sample is screened (step `amr_chroquetas`) and gets `05.amr/{sample}/amr_summary.tsv`, and MultiQC shows the AMR table. `links.chroquetas_link` is then required. With `off`, neither step runs, and the link is not read. |
| `mito_genetic_code` | `auto` | `auto` (any case, or an empty value), or the number of an NCBI genetic code, 1-6, 9-16 or 21-33 (7, 8 and 17-20 are not NCBI codes, and miniprot refuses them) | AMR on | The genetic code of the second ChroQueTas run, on the mitochondrial genome, where cytochrome b is encoded. That run is made only for the 15 species whose database entry includes cytochrome b. With `auto`, the code is chosen per sample: code 3 (yeast mitochondrial) for the four Saccharomycetaceae of the database (*Saccharomyces cerevisiae*, *S. paradoxus*, *Nakaseomyces glabratus*, *Kluyveromyces marxianus*), and code 4 (mold mitochondrial) for every other species. A number sets that code for every sample. The nuclear run always uses the code ChroQueTas chooses for the species ([Two runs](../analysis/amr.md#two-runs)). |
| `chroquetas_script_sha256` | the checksum of `bin/ChroQueTas.sh` at the template address | 64 hexadecimal characters (capitals are accepted) | AMR on | The SHA-256 that the downloaded script must have, checked before anything runs it. |
| `chroquetas_db_sha256` | the checksum of `FungAMR_db/db-v20250811.tgz` at the template address | 64 hexadecimal characters (capitals are accepted) | AMR on | The SHA-256 that the downloaded database archive must have, checked before it is unpacked. |

The template address names one commit, so the files behind it do not change. To use another
version, change the link and both checksums together (`sha256sum <file>` prints them). A wrong value
stops the run with `[FunFlux] parameters.amr.<key>: ...`, for example
`[FunFlux] parameters.amr.mito_genetic_code: must be auto or the number of an NCBI genetic code: 1-6, 9-16 or 21-33 (got: <value>)`.
Every dry run prints one line, `AMR screen: ChroQueTas on (licence CC BY-NC-ND 4.0: non-commercial use).`
or `AMR screen: off.`

A change re-runs the following steps:

| Change | Re-runs |
|---|---|
| `mito_genetic_code` | `amr_chroquetas` and MultiQC |
| `links.chroquetas_link`, `chroquetas_script_sha256` or `chroquetas_db_sha256` | `chroquetas_db`, then `amr_chroquetas` of every sample |
| the species, column 2 of the annotation table | also the gene prediction and everything after it, because funannotate reads the species too |

The `amr` keys and the link re-run only the AMR steps and the report, while a new species
re-runs the whole annotation.

### `longread` and `hybrid`

These two blocks hold the long-read settings: the assembler, Filtlong, Medaka and hifiasm, and in
`hybrid` mode also the Filtlong guidance.

```yaml
parameters:
  longread:
    assembler: hifiasm         # hifiasm | flye
    flye_input_mode: auto      # auto | nano-raw | nano-hq
    medaka_model: auto
    medaka_model_fallback_auto: false
    filtlong:
      min_length: 1000
      keep_percent: 90
      length_weight: 10
      target_bases:
    hifiasm:
      genome_size:
      bloom_filter_bits: 0
      purge_dups: auto
      homozygous_coverage: auto
      telomere_motif:
  hybrid:
    # the same keys as longread, plus:
    short_read_guidance: false
    filtlong_split: 1000
```

The long-read settings sit in the block named after the mode, `parameters.longread` in `longread`
mode and `parameters.hybrid` in `hybrid` mode. Each block is read only in its own mode. The same
keys at the top of `parameters`, where FunFluxL and FunFluxL+ kept their long-read settings, stop
the run in every mode.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `assembler` | `hifiasm` | `hifiasm`, `flye` | `longread`, `hybrid` | One assembler for the whole run. Flye does not repeat its assembly byte for byte from run to run. |
| `flye_input_mode` | `auto` | `auto`, `nano-raw`, `nano-hq`, while `pacbio-hifi` stops the run | `longread`, `hybrid` | The Flye read mode of ONT samples. `auto` gives `--nano-hq`, or `--nano-raw` when `medaka_model` names a model containing `fast`. PacBio samples always get `--pacbio-hifi`. Ignored with hifiasm, with a NOTE when it is set. |
| `medaka_model` | `auto` | `auto` (also empty or an unquoted `null`, `true`, `yes`, `on`, `1`), `false` for no Medaka (also `no`, `off`, `0`, `none`, `null` in quotes), or any other text, taken as a model name | `longread`, `hybrid` | Medaka polishing of ONT samples. PacBio samples never use Medaka. `auto` takes the model from the basecaller tag in the first 100 reads of the raw FASTQ (exactly one tag must be present) and needs network access once. The model is checked before the assembler starts. |
| `medaka_model_fallback_auto` | `false` | `true`, `false` | `longread`, `hybrid` | When a named model is not valid, try `auto` instead of stopping. |
| `filtlong.min_length` | 1000 | an integer ≥ 1 | `longread`, `hybrid` | Filtlong `--min_length`, in bp. |
| `filtlong.keep_percent` | 90 | a number above 0 and at most 100 | `longread`, `hybrid` | `--keep_percent`, the share of the bases kept, best reads first. |
| `filtlong.length_weight` | 10 | a number above 0 | `longread`, `hybrid` | `--length_weight`. Higher values favour longer reads. |
| `filtlong.target_bases` | empty | empty or 0 for no cap, or a size such as `1500000000` or `1.5g` | `longread`, `hybrid` | `--target_bases`, a cap on the bases kept, for very deep runs. With Flye and very deep ONT reads, set it to about 100 times the genome size ([Flye](../modes/longread.md#flye)). |
| `hifiasm.genome_size` | empty | empty or 0, or a size such as `40000000` or `40m` | `longread`, `hybrid` | `--hg-size`, the haploid genome size. When it is empty, hifiasm infers the depth itself. |
| `hifiasm.bloom_filter_bits` | 0 | an integer from 0 to 63 | `longread`, `hybrid` | `-f`, where 0 suits small fungal genomes. |
| `hifiasm.purge_dups` | `auto` | `auto` (also empty), or 0 to 3 | `longread`, `hybrid` | `-l`, how strongly duplicated haplotigs are purged. `auto` keeps hifiasm's own choice. |
| `hifiasm.homozygous_coverage` | `auto` | `auto` (also empty), or an integer ≥ 1 | `longread`, `hybrid` | `--hom-cov`, the read depth of the homozygous peak. |
| `hifiasm.telomere_motif` | empty | DNA or IUPAC letters, such as `CCCTAA` | `longread`, `hybrid` | `--telo-m`. Set it only when the organism's telomere repeat is known. |
| `short_read_guidance` | `false` | `true`, `false` | `hybrid` | With `true`, Filtlong scores each long read of an ONT sample against the selected pairs, the Illumina pairs properly paired on the screened guide assembly (`-1`, `-2`). With `false`, every sample gets the plain filter of `longread` mode. PacBio samples are never guided. |
| `filtlong_split` | 1000 when the key is missing | an integer ≥ 1, or empty or `null` to score but never cut | `hybrid` | With guidance, adds `--trim --split <N>`, so that a read is cut where the Illumina reads stop supporting it for N bp. Not used when `short_read_guidance` is false, and a number other than 1000 then gets a NOTE. |

PacBio samples are not guided, because cutting accurate HiFi reads where the Illumina coverage has
a gap would throw away their length. Filtlong's README advises against a short-read reference
unless you are confident in the quality of the short reads, since long reads that span a region the
short reads miss can score low. That is why the guidance is off by default.
`short_read_guidance: true` turns it on, and `filtlong_split: null` then keeps the scoring without
cutting.

The Filtlong defaults are sized for ONT reads and provisional for PacBio HiFi. hifiasm always runs
with `--primary`. Every dry run prints the decisions for each sample, and
[One technology per sample](../modes/longread.md#one-technology-per-sample) lists what differs
between ONT and PacBio HiFi samples.

See [Long reads](../modes/longread.md), [Hybrid](../modes/hybrid.md),
[The Medaka model](../methods/medaka-model.md) and
[Filtlong guidance in hybrid mode](../methods/filtlong-guidance.md).

### `contigs`

The `contigs` block has one key, which switches the screen on for supplied assemblies.

| Key | Default | Allowed values | Modes | What it does |
|---|---|---|---|---|
| `screen` | `false` | `true`, `false` | `contigs` (checked in every mode) | With `true`, the contamination and mitochondrion screen runs on the supplied FASTA. It needs `directories.blast_db`, and adds `contigs_filt.fasta`, the `contaminants/` and `mitochondria/` folders, `contigs_all.fasta` and `contigs_mitochondrial.fasta` to the outputs. With `false` there is no screen, and `contigs_final.fasta` is the FASTA after the contig filter. With SPAdes-style headers the filter keeps the contigs of at least 500 bp and 2x coverage, and other headers are only cut at the first space or tab. With `false`, a NOTE also warns that a mitochondrial genome in `contigs_final.fasta` would be annotated as nuclear. |

Assemblies from a sequencing provider usually still contain the mitochondrial genome. See
[Pre-assembled contigs](../modes/contigs.md).

## Annotation parameter table

The annotation parameter table gives FunFlux the species, the protein evidence and the AUGUSTUS
model of each sample. It is a tab-separated text file, named by `files.annotation_params`, with
three or four columns (the fourth is optional) and one row per sample.

| Column | Content | Used as |
|---|---|---|
| 1 `sample` | the sample name, as found in the input folder | the link between the row and the sample |
| 2 `species` | the species name, such as `Fusarium graminearum` | `funannotate predict -s`, the organism of the GenBank files and of the `.sqn` submission file. With spaces turned into `_` and `_{sample}` appended, it is the stem of funannotate's file names (`Fusarium_graminearum_SAMPLE1.gbk`). With the AMR screen on, it is also the species looked up in the FungAMR database ([The species](../analysis/amr.md#species)). |
| 3 `proteins` | path of a protein FASTA of related fungi, such as the RefSeq proteome of a species of the same genus ([how to download one](../getting-started/quick-start.md#6-write-the-annotation-parameter-table)) | gene evidence, as `funannotate predict --protein_evidence <file> <funannotate_db>/uniprot_sprot.fasta` |
| 4 `augustus_species` | optional: the name of a stored AUGUSTUS model, or empty, or `-` | a name gives `funannotate predict --augustus_species <name>`, and AUGUSTUS uses that model as it is, without training on this genome. Empty or `-`, or no fourth column at all, means that AUGUSTUS is trained on the genome's BUSCO genes, the long way with [`prediction.augustus_optimize`](#prediction) |

```text
#sample	species	proteins	augustus_species
SAMPLE1	Fusarium graminearum	/path/to/fusarium_proteins.faa	-
ARSEF3097	Beauveria bassiana	/path/to/beauveria_bassiana_proteins.faa	fusarium_graminearum
```

The first row trains AUGUSTUS on the genome (the default), and the second uses the stored
`fusarium_graminearum` model as it is.

FunFlux checks the table before any job starts:

- Lines starting with `#` and blank lines are skipped, and spaces around each field are removed.
  Column names are not read, so the position decides. The shipped example
  `config/annotation_parameters.tsv`, for instance, calls the fourth column `model`.
- Every data line, including a line for a sample outside the batch, has three or four tab-separated
  fields, and no sample name appears twice. Columns separated by spaces give, for example,
  `[FunFlux] files.annotation_params: line 2 has 1 tab-separated columns; expected 3 or 4 (sample, species, proteins and, optionally, augustus_species).`
- Every sample of the batch has a row, and its protein file exists and is not empty. A relative
  path is read from the launch folder. The protein files of samples outside the batch are not
  checked, so one table can serve several batches.

A named AUGUSTUS model must be one of the stored models, the folder names in
`<funannotate_db>/trained_species/`. The default funannotate database holds 166 of them, including
`fusarium_graminearum`, `aspergillus_nidulans`, `botrytis_cinerea` and `neurospora_crassa`. With
`FUNANNOTATE_DB` set, `funannotate species` lists them. Capitals count. When the name is checked
depends on whether the database is on disk yet, and every dry run prints the choice for each
sample:

| The funannotate database | The name is checked | The dry-run line, for example |
|---|---|---|
| is on disk | before any job starts | `AUGUSTUS: SAMPLE1 uses the stored model fusarium_graminearum.` |
| is not downloaded yet | by rule `augustus_model_check`, right after the database is unpacked, with its log in `logs/augustus_model_check_{sample}.log` | `SAMPLE1 names the stored model fusarium_graminearum (checked once the database is unpacked)`, and a sample with an empty column 4 gets the same `(checked once the database is unpacked)` |

A name that is not in the database stops the run, at the latest right after the database
is unpacked. The message gives the closest names, for example
`[FunFlux] files.annotation_params: sample 'SAMPLE1' names the AUGUSTUS model 'fusarium_graminaerum', which is not a trained species of the funannotate database at '<db>'. Did you mean 'fusarium_graminearum'? ...`.

With column 4 empty, AUGUSTUS is trained on the genome and needs at least 200 BUSCO models for it
([The AUGUSTUS model](../analysis/gene-prediction.md#augustus-model)). If the database holds a
model stored earlier under the sample's own name, funannotate uses that instead, and a NOTE says
so. Which of the two to choose is discussed on
[Gene prediction](../analysis/gene-prediction.md#augustus-model). A changed row re-runs gene
prediction and every annotation step of that sample.

## Coming from a v1 configuration

A v1 configuration stops at once, because it has no `mode`. Start from the v2 template and carry
your values over. [Coming from v1](../getting-started/from-v1.md) lists every key that was renamed,
moved or removed, and every v1 value that now stops the run.
