# Installation

This page lists what FunFlux needs on the computer and how to install it before the first run.
You install conda, Snakemake, InterProScan and one or two databases yourself. Snakemake installs
every other tool, each in its own conda environment, and the workflow downloads the other
databases.

The installation takes six steps, in this order:

| Step | What you do | Where |
|---|---|---|
| 1 | Clone FunFlux | [Get FunFlux](#get-funflux) |
| 2 | Install conda, then Snakemake in the launcher environment | [The launcher environment](#the-launcher-environment) |
| 3 | Check the commands FunFlux takes from the machine | [Commands outside the conda environments](#host-commands) |
| 4 | Download the eggNOG-mapper data and, for the contamination and mitochondrion screen, NCBI `core_nt` | [Reference databases](databases.md) |
| 5 | Install InterProScan and, optionally, GeneMark-ES and SignalP 6 | [Licensed and external tools](external-tools.md) |
| 6 | Optionally, build the tool environments before the first run, and on a cluster also fetch the downloads | [Building the environments in advance](#build-in-advance), [On a cluster](../reference/running.md#on-a-cluster) |

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **launcher environment** | the conda environment you create by hand for Snakemake itself; activate it before every run, and install no analysis tool into it |
    | **launch folder** | the folder Snakemake is started from, normally the FunFlux repository; relative paths are read from it, and `.snakemake/` is created there |
    | **rule** | one step of the workflow as Snakemake defines it, such as `select_contigs`; the console and the log files in `logs/` use its name |
    | **job** | one run of a rule, for one sample or once for the whole run; the dry run lists the jobs it plans |
    | **heavy step** | on a cluster, a step that asks the scheduler for both `ram_gb` and `runtime_hours`, such as the assemblers, BLAST and InterProScan |
    | **pin file** | the file `<name>.linux-64.pin.txt` beside each environment file: the exact package builds FunFlux v2 was tested with, installed without solving dependencies |
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against `core_nt`, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **AMR screen** | the optional step, off by default (`parameters.amr.chroquetas`), in which ChroQueTas looks in each genome for the changes FungAMR lists for the species; results in `05.amr/` |

## System requirements

Check that the computer meets these requirements before you install FunFlux:

| Item | Requirement |
|---|---|
| Operating system | Linux on x86-64. The environment pin files are for `linux-64`, and GeneMark-ES and InterProScan are distributed as 64-bit Linux builds. |
| conda | 24.7.1 or newer. Snakemake 9.14.6 stops with an error on older versions. |
| Snakemake | 9.14.6 or a later 9.x release. FunFlux v2 was developed and tested with 9.14.6 and stops at start-up on an older version. Version 10 or later has not been tested. |
| Commands | `bash` and the GNU base commands (`sha256sum`, `md5sum` and the rest of coreutils, `find`, `grep`, `sed`, `awk`) on the machine's own `PATH`, outside conda. With the AMR screen on, also a `grep` that knows `-P` (GNU grep). No `python`, `perl`, `curl`, `tar` or `gzip` is needed. See [Commands outside the conda environments](#host-commands). |
| CPU | the number in `resources.threads` (template: 24), given again as `--cores` |
| Memory | about 16 GB for Tiberius, on by default, and more than 14 GB for SignalP 6 in `slow` mode. See [Memory](#memory). |
| Disk | about 380 GB of databases. The conda environments take 20 to 25 GB for a `contigs` run, 25 GB with the screen and SignalP 6, and the read-based modes add 4 to 9 environments, not measured. The output takes in the order of 10 GB per sample in `contigs` mode, plus the reads and assembly files of the read-based modes. See [Disk](#disk). |
| Network | conda packages and several database downloads, on the computer that runs Snakemake; inside the jobs, only BUSCO with the automatic lineage. See [Network access](#network). |

### Memory { #memory }

With the default settings, plan for at least 16 GB of free memory, because of Tiberius. The table
gives what is known of the memory each step needs.

| Step | Memory |
|---|---|
| Tiberius, on by default | about 16 GB for a fungal genome of about 40 Mb; switch it off (`parameters.tiberius.mode: off`) on a computer with less free memory |
| SignalP 6 in `slow` mode | more than 14 GB ([SignalP 6](external-tools.md#signalp-6)) |
| SPAdes (`illumina`, `hybrid`) | `resources.ram_gb` (template: 128), given to SPAdes as its memory limit (`-m`) |
| Qualimap (`illumina`, `longread`, `hybrid`) | `resources.ram_gb`, capped at 64 GB, as its Java heap |
| The other steps | not measured |

On a single computer, `resources.ram_gb` limits only SPAdes and Qualimap. Snakemake itself starts
jobs by thread count alone and does not limit memory, unless `--resources mem_gb=N` is passed to
cap how many SPAdes jobs run at once, or `java_mem=N` for Qualimap. See
[CPU cores and threads](../reference/running.md#cpu-cores-and-threads). On a cluster, `ram_gb` is
also the memory every heavy step asks the scheduler for, as described under
[How the settings reach Slurm](../reference/running.md#cluster-resources).

### Run time { #run-time }

One fungal genome in `contigs` mode takes several hours, and InterProScan alone usually takes a
few hours for a fungal proteome. Each rule declares its threads, and Snakemake runs jobs side by
side as long as their threads fit into `--cores`. The steps that take longest are listed on
[Run time](../reference/running.md#run-time), and the threads each step asks for on
[CPU cores and threads](../reference/running.md#cpu-cores-and-threads).

### Disk { #disk }

Disk space goes to three parts: the databases, the conda environments and the output of each
sample. NCBI `core_nt` takes 264 GB of the about 380 GB of databases, and only the screen reads
it.

| What | Size | Where |
|---|---|---|
| NCBI `core_nt`, build of 11 July 2026; its taxonomy files add 0.8 GB | 264 GB | `directories.blast_db`; only for the screen |
| InterProScan 5.78-109.0, installed | 50 GB | `files.iprscan` |
| eggNOG 5.0.2 data | 48 GB | `directories.eggnog_db` |
| antiSMASH 8.0.4 databases | 9.4 GB | `00.resources/` of every output folder, or once in `directories.antismash_db` |
| funannotate database, unpacked, archive included | 8.1 GB | `directories.funannotate_db` |
| conda environments: the 18 of a `contigs` run with the screen, Tiberius and SignalP 6 | about 25 GB; the `download` environment among them takes 120 MB | `.snakemake/conda/`, or the folder given to `--conda-prefix` |
| conda environments added by the read-based modes: 4 in `illumina`, 4 or 5 in `longread`, 8 or 9 in `hybrid`, listed under [`--sdm conda`](#sdm-conda) | not measured | the same folder |
| conda environment of the AMR screen (`chroquetas`), when it is on | about 2 GB | the same folder |
| output of one sample, `contigs` mode | in the order of 10 GB | `directories.output_dir` |
| output of one sample, read-based modes: the kept reads (`01.reads/`), the assembler's working folder, the read-mapping BAM, the polishing files | whole sample: 13 GB in `illumina`, 20 to 31 GB in `longread` and 39 GB in `hybrid` on the 42.6 Mb genome of *Botrytis cinerea* B05.10; 16 to 17 GB in `longread` and 29 GB in `hybrid` on the 32.7 Mb genome of *Cordyceps cateniannulata* GXU-8616, with PacBio HiFi reads | `directories.output_dir` |

In the output of the read-based modes, the filtered long reads are kept uncompressed, and so are
the selected Illumina pairs of `hybrid`.

Sizes are as `du -h` prints them, in powers of 1024. The environments were measured with the
prefix on another file system than conda's package cache. On the same file system, conda
hard-links files from its package cache, and the environments add less to the disk.

Each measured run is listed under [Disk use per run](../reference/output.md#disk-use-per-run), and
[Disk use](../reference/output.md#disk-use) shows where the space of the output folder goes.

### Network access { #network }

Most downloads happen on the computer that runs Snakemake. The download steps are local rules, so
they run on that computer and never on a compute node of a cluster. Snakemake builds the conda
environments there too. Only three things reach the network from inside a job: BUSCO with the
automatic lineage, and two lookups that the results do not need.

| Host | What | When | Runs on |
|---|---|---|---|
| `conda.anaconda.org` | the packages listed in the pin files | building an environment | the computer that runs Snakemake |
| `zenodo.org` | the funannotate database and its checksum file | while `directories.funannotate_db` does not hold it (`funannotate_db`) | the computer that runs Snakemake |
| `s3.hpc.ut.ee` | UNITE | first run in an output folder, unless `directories.unite_db` is set (`unite_db`) | the computer that runs Snakemake |
| `dl.secondarymetabolites.org`, `ftp.ebi.ac.uk` | antiSMASH databases | first run in an output folder, unless `directories.antismash_db` is set (`antismash_db`) | the computer that runs Snakemake |
| `bioinf.uni-greifswald.de` | the Tiberius `fungi` model | first run in an output folder, unless `parameters.tiberius.mode` is `off` (`tiberius_model`) | the computer that runs Snakemake |
| `raw.githubusercontent.com` | the telomere repeat database of tidk; ChroQueTas and its FungAMR database | first run in an output folder (`telomere_db`; with the AMR screen on also `chroquetas_db`) | the computer that runs Snakemake |
| `ftp.ncbi.nlm.nih.gov` | PhiX genome | first run in an output folder, `illumina` and `hybrid` (`download_phix`) | the computer that runs Snakemake |
| chosen by Medaka (GitHub) | model weights | the model check of an ONT sample with Medaka, the first time a model is used (`check_medaka_model`), before the assembly | the computer that runs Snakemake |
| `busco-data.ezlab.org` | BUSCO lineage datasets | with a named lineage, once per output folder (`busco_data`) | the computer that runs Snakemake |
| `busco-data.ezlab.org` | BUSCO lineage datasets | with `parameters.busco.lineage: auto`, the default: every BUSCO job (`completeness`) | the job |
| `www.ebi.ac.uk` | InterProScan's lookup service, a shortcut: without it, InterProScan computes every match itself ([Functional annotation](../analysis/functional-annotation.md#interproscan)) | every InterProScan job | the job |
| `api.multiqc.info` | MultiQC's check for a newer version | the report (`multiqc`) | the job |

FunFlux runs BUSCO 6.0.0 with `--opt-out-run-stats`. Without it, BUSCO ends each job by looking up
the country of the machine's IP address online and sending usage statistics. See
[Assembly QC](../analysis/assembly-qc.md#busco).

When the jobs have no internet, BUSCO, InterProScan and MultiQC behave as follows.

- **BUSCO**, with `lineage: auto`, cannot download its lineage data, and the BUSCO job fails. Name
  a lineage instead, as described under
  [BUSCO on compute nodes without internet](../reference/running.md#cluster-busco).
- **InterProScan** prints
  `Pre-calculated match lookup service failed - analysis proceeding to run locally` and computes
  every match itself. The matches are the same, and the run is longer
  ([Functional annotation](../analysis/functional-annotation.md#interproscan)).
- **MultiQC** gives up its version check after about 2 seconds, and the report is the same.

!!! warning "Compute nodes without internet need a named BUSCO lineage"

    On a single computer with internet, nothing needs preparing. On a cluster whose compute nodes
    have no internet, name a BUSCO lineage. The rest is fetched on the login node.

Snakemake itself may have to run on a computer without internet, for example when the whole run
goes inside one batch job on a compute node. In that case, fetch everything first on a computer
with internet, with `snakemake download_all`, as described under
[Preparation, on a computer with internet](../reference/running.md#cluster-preparation).

#### Behind a proxy that re-signs HTTPS { #proxy }

FunFlux's downloads trust only the list of certificate authorities inside their conda
environments, not the computer's own. Behind an institute proxy or mirror that re-signs HTTPS with
an in-house certificate, they therefore fail with `SSL certificate problem` or
`CERTIFICATE_VERIFY_FAILED`.

`curl`, in the six download steps that use it, and Medaka read `CURL_CA_BUNDLE`, while BUSCO and
antiSMASH read `SSL_CERT_FILE`. Point both variables at the computer's list before you start
Snakemake:

| System | Commands |
|---|---|
| Debian, Ubuntu | `export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt` |
| Red Hat Enterprise Linux and its rebuilds (Rocky, Alma) | `export CURL_CA_BUNDLE=/etc/pki/tls/certs/ca-bundle.crt SSL_CERT_FILE=/etc/pki/tls/certs/ca-bundle.crt` |

## Get FunFlux

The FunFlux repository holds the workflow, the configuration template and the environment files.
Clone it:

```bash
git clone https://github.com/iLivius/FunFlux.git
cd FunFlux
```

The clone holds, among other files:

| Path | Content |
|---|---|
| `config/config.yaml` | the commented configuration template; copy it, for example to `config/config_custom.yaml`, which git ignores |
| `config/annotation_parameters.tsv` | an example annotation parameter table, the table `files.annotation_params` names |
| `workflow/envs/` | `<name>.yaml`, one conda environment per tool, each with its pin file `<name>.linux-64.pin.txt` (see [`--sdm conda`](#sdm-conda)) |

Start Snakemake from this folder, the launch folder. Snakemake finds `workflow/Snakefile` there
and keeps its run records in `.snakemake/`, which git ignores. The conda environments go there
too, unless `--conda-prefix` names another folder.

## The launcher environment

Snakemake runs in a conda environment of its own, the launcher environment. Create it and check
the version:

```bash
conda create -c conda-forge -c bioconda -n snakemake "snakemake>=9.14.6,<10"
conda activate snakemake
snakemake --version
```

The launcher environment holds Snakemake and its Python, which reads the configuration and runs
FunFlux's checks before any job starts. Activate it before every run. Install no analysis tool
into it, because each tool runs in its own environment and every FunFlux script runs with the
pinned Python of its rule's environment.

Snakemake prints a warning unless conda uses strict channel priority, set with
`conda config --set channel_priority strict`. To let Snakemake send the steps to a Slurm cluster,
add `snakemake-executor-plugin-slurm` to this environment, as described under
[Run A: Snakemake sends each step to Slurm](../reference/running.md#cluster-run-a).

FunFlux sets `PYTHONNOUSERSITE=1` for every job, so packages installed with `pip install --user`
cannot be used instead of the packages of a rule environment.

### Commands outside the conda environments { #host-commands }

A few commands come from the computer itself, not from a conda environment, and FunFlux does not
install them. Two kinds of rule take commands from outside a tool environment:

- **Rules without a conda environment** run standard shell commands only, and take them from the
  launcher environment and the machine. They are the contig filter (`filter_contigs`),
  `finalize_assembly`, the links to local databases (`unite_db_local`, `antismash_db_local`) and
  `cleanup`.
- **Rules with a conda environment** take their commands from the environment, and from the
  machine for every command the environment lacks. Most environments lack the standard commands
  of the shell lines around the tool, such as `awk`, `grep`, `sed`, `find`, `mkdir` and `rm`. The
  download steps take `sha256sum` and `md5sum` from the machine.

!!! warning "Every command in the next table must be installed on the machine itself, outside conda"

    When Snakemake activates a rule's environment, conda, with its default settings, takes the
    launcher environment off the `PATH`. The rules that run in a conda environment therefore find
    these commands only on the machine's own `PATH`.

| Command | Used for |
|---|---|
| `bash` | every rule: Snakemake runs the shell lines of each rule with bash |
| `sha256sum`, `md5sum` | checking the downloads: the funannotate database, the Tiberius model, the telomere repeat database and ChroQueTas (`sha256sum`), UNITE (`md5sum`) |
| `awk`, `grep`, `sed`, `find`, and basic commands such as `head`, `cat`, `cp`, `ln`, `mkdir`, `rm` | the contig filter, copies, links to local databases, and the shell lines around the tool in most rules (`awk` and `grep` in `its_windows` and `signalp6_annotation`, `grep` in `genemark_check`, `sed` in `eggnog_annotation`, `find` in `completeness` and `iprscan_annotation`; in the read-based modes `awk` turns hifiasm's graph into FASTA and cuts the Polypolish headers, `grep` counts the SPAdes contigs the filter kept, `cat` reads the Medaka model file) |
| `grep` with `-P` (GNU grep), `sed`, and basic commands such as `cut`, `paste`, `tr` | ChroQueTas (rule `amr_chroquetas`, AMR screen on), a bash script that runs inside the chroquetas environment; the environment holds the tools it calls, not these commands |

These commands come from conda instead, so the machine does not need them:

| Command | Comes from |
|---|---|
| `python` | the pinned Python of each rule's environment ([Which Python runs the helper scripts](../architecture.md#which-python-runs-the-helper-scripts)) |
| `curl`, `tar`, `gzip` | the `download` environment of the six download steps that use them (`download_phix`, `funannotate_db`, `unite_db`, `telomere_db`, `tiberius_model`, `chroquetas_db`) |
| `perl` | the `interproscan` environment, for InterProScan's SUPERFAMILY analysis and PROSITE patterns; the `funannotate` environment for GeneMark-ES, once its scripts start the first `perl` on the `PATH` ([GeneMark-ES](external-tools.md#genemark-es)) |

Check the machine's commands with no conda environment active:

```bash
conda deactivate     # repeat until no environment is active, base included
command -v bash sha256sum md5sum awk grep sed find
command -v cut paste tr && echo ok | grep -P 'o\w'   # AMR screen on: the last line reads ok
```

Install a missing command with the system package manager. Installing it into the launcher
environment does not reach the rules that run in a conda environment. FunFlux has been tested on
x86-64 Linux only.

## `--sdm conda`: one pinned environment per tool { #sdm-conda }

Every tool has its own environment file in `workflow/envs/`. With `--sdm conda`
(short for `--software-deployment-method conda`), Snakemake builds the environments of the planned
jobs before the first job starts and reuses them in later runs. A dry run lists the missing ones
as `Conda environment workflow/envs/<name>.yaml will be created.`

Beside each `<name>.yaml` is its pin file, `<name>.linux-64.pin.txt`: the exact package builds of
an environment FunFlux v2 was tested with, as written by `conda list --explicit`.

What Snakemake installs depends on whether it can use the pin file:

| Situation | What Snakemake installs |
|---|---|
| The pin file exists | that list, without solving dependencies; it logs `Using pinnings from <path>/<name>.linux-64.pin.txt.` |
| Installing the pin file fails | it prints a warning and solves `<name>.yaml` instead |
| A platform other than `linux-64` | there is no pin file, and Snakemake always solves the yaml |

!!! note "Only the pin file reproduces the tested builds"

    The yaml fixes the versions of the main tools and of Python, but few of their dependencies.
    When Snakemake solves it, the versions and builds of everything else can differ from the
    pinned ones. That includes tools in the table below, such as AUGUSTUS, DIAMOND and HMMER, and
    the OpenJDK and PyTorch patch releases.

The table lists each environment, the tools in it and the runs that build it.

| Environment | Tools, versions from the pin file | Built when |
|---|---|---|
| `funannotate` | funannotate 1.8.17 with AUGUSTUS 3.5.0, SNAP 2017_03_01, GlimmerHMM 3.0.4, EVidenceModeler 1.1.1, DIAMOND 2.2.4, tRNAscan-SE 2.0.13, Perl 5.32.1 | always |
| `interproscan` | OpenJDK 11.0.30, Python 3.11.15 and Perl 5.32.1, the Java, Python and Perl your InterProScan runs with | always |
| `eggnog-mapper` | eggNOG-mapper 2.1.13, DIAMOND 2.0.15 | always |
| `antismash` | antiSMASH 8.0.4 | always |
| `busco` | BUSCO 6.0.0 | always |
| `quast` | QUAST 5.3.0 | always |
| `itsx` | ITSx 1.1.3, HMMER 3.4, Python 3.12.13 | always |
| `vsearch` | VSEARCH 2.30.0, Python 3.12.15 | always |
| `tidk` | tidk 0.3.0, Python 3.12 | always |
| `download` | curl 8.22.0, tar 1.35, gzip 1.14, for the download steps; 120 MB | always |
| `python` | Python 3.12.15 and nothing else, for the FunFlux scripts whose rule has no tool environment: the contig selector, the mitochondrion screen and split, the assembler's evidence table, the repair of the funannotate database paths, the AUGUSTUS model check, the ONT adapter scan and the adapter screen report | always |
| `repeatmasking` | tantan 51, RepeatModeler 2.0.8, RepeatMasker 4.2.3 | always |
| `multiqc` | MultiQC 1.33 | always |
| `blast` | BLAST+ 2.16.0 | always (the contig adapter screen; also the contamination screen) |
| `blobtools` | BlobTools 1.1.1 | the screen |
| `minimap` | minimap2 2.30, samtools 1.21 | the screen in `contigs`, `longread` and `hybrid` mode (`illumina` maps with bowtie2) |
| `signalp6` | Python 3.10, PyTorch 1.13.1 (CPU); SignalP 6 itself comes from your archive | `parameters.signalp6.mode` not `off` |
| `tiberius` | Tiberius 1.1.8, TensorFlow 2.17.0 (CPU build), Python 3.12.15; about 2 GB | `parameters.tiberius.mode` `on`, the default |
| `bowtie` | bowtie2 2.5.4, samtools 1.21, Python 3.14.6 | `illumina`, `hybrid` (PhiX removal, read mapping, read selection) |
| `fastp` | fastp 1.0.1 | `illumina`, `hybrid` |
| `spades` | SPAdes 4.2.0, Python 3.14.6 | `illumina`, `hybrid` |
| `qualimap` | Qualimap 2.3, OpenJDK 25.0.2, Python 3.14.6 | `illumina`, `longread`, `hybrid` |
| `fastplong` | fastplong 0.7.1 | `longread`, `hybrid`: ONT samples, unless `ont_adapter: off` |
| `filtlong` | Filtlong 0.3.1 | `longread`, `hybrid` |
| `nanoplot` | NanoPlot 1.46.2, Python 3.13.14 | `longread`, `hybrid` |
| `hifiasm` | hifiasm 0.25.0 | `longread`, `hybrid` with `assembler: hifiasm`, the default |
| `flye` | Flye 2.9.6, Python 3.13.16 | `longread`, `hybrid` with `assembler: flye` |
| `medaka` | Medaka 2.2.2, Python 3.12.13 (with its own minimap2 2.31 and samtools 1.23.1) | `longread`, `hybrid`: ONT samples with Medaka (the model check and the polishing) |
| `polypolish` | Polypolish 0.6.1, BWA 0.7.19 | `hybrid` |
| `chroquetas` | miniprot 0.14, BLAST+ 2.16.0, MAFFT 7.525, SeqKit 2.1.0, EMBOSS 6.6.0 (the versions the ChroQueTas authors pin in their own conda package), Python 3.12.15; ChroQueTas itself is downloaded by rule `chroquetas_db`; about 2 GB | `parameters.amr.chroquetas: on` |

By default every mode builds the same 15 environments, the `always` rows and tiberius, most of
them for the steps after the assembly. The screen and the read-based front ends add theirs to
these 15, up to 27 in `hybrid`:

| Mode | Added | Environments |
|---|---|---|
| `contigs` with the screen | blobtools, minimap | 17 |
| `illumina` | blobtools, bowtie, fastp, spades, qualimap | 20 |
| `longread` | blobtools, minimap, filtlong, nanoplot, qualimap, hifiasm or flye; fastplong for ONT samples (unless `ont_adapter: off`); medaka for ONT samples with Medaka | 21 to 23 |
| `hybrid` | the `longread` set plus bowtie, fastp, spades, polypolish | 25 to 27 |

In every mode, SignalP 6 adds signalp6 and the AMR screen adds chroquetas, while switching
Tiberius off removes tiberius.

### Building the environments in advance { #build-in-advance }

The environments can be built before the first run, with one command:

```bash
snakemake --sdm conda --conda-create-envs-only --configfile config/config_custom.yaml
```

This builds every environment the configured run needs and exits without starting a job.
Snakemake plans the run to find them, so the configuration must pass FunFlux's checks (mode, input
folder, annotation table). With `--list-conda-envs` in place of `--conda-create-envs-only`,
Snakemake prints each environment file with the folder it is built in.

Next, run `snakemake download_all` to fetch, in the same way, every file the run downloads:
databases, models, a named BUSCO lineage and the Medaka weights. It is needed only when Snakemake
itself will run without internet, as in a batch job on a compute node. See
[Preparation, on a computer with internet](../reference/running.md#cluster-preparation).

## Sharing environments: `--conda-prefix` { #conda-prefix }

Without `--conda-prefix`, the environments are built in `.snakemake/conda/` of the launch folder,
so all runs started from one clone share them. `--conda-prefix` puts them in a folder of your
choice, which several clones, users or projects can share:

```bash
snakemake --sdm conda --conda-prefix /path/to/funflux_envs --cores 24 \
    --configfile config/config_custom.yaml
```

`export SNAKEMAKE_CONDA_PREFIX=/path/to/funflux_envs` before the command has the same effect.

- **Each environment** is a folder named `<hash>_`, a code Snakemake computes from the yaml file,
  the pin file and the path of the prefix. A FunFlux update that changes a yaml or pin file builds
  a new folder beside the old one.
- **Moving or renaming the prefix** rebuilds every environment, because conda environments
  contain absolute paths. In an output folder made from the same launch folder, Snakemake then
  re-runs every finished step that has an environment, and the steps after it, as described under
  [Conda environments](../reference/running.md#conda-environments). Choose the prefix once.
- **Old folders** stay until you delete them. To see which ones are in use, run
  `--list-conda-envs` with the same `--conda-prefix` and configuration. A configuration without
  the screen or SignalP 6 lists fewer. A `<hash>_` folder, with its `<hash>_.*` files, can be
  deleted when no configuration or clone using the prefix lists it. In a shared prefix, collect
  the lists of all configurations first. Do not use `--conda-cleanup-envs` for this, because in
  Snakemake 9.14.6 it deletes the environments of the current workflow.
- **With the screen or SignalP 6**, the prefix must be writable by whoever runs FunFlux. BlobTools
  writes its taxonomy index (`nodesDB.txt`) into its environment, and `install_signalp6` installs
  SignalP 6 into the `signalp6` environment.
- **A rebuilt `signalp6` environment** no longer contains SignalP 6. The next run stops in
  `signalp6_annotation` and asks you to delete `00.resources/signalp6/installation.json`. See
  [Troubleshooting](../troubleshooting.md#signalp-6).
- **With SignalP 6 on**, everyone who uses the prefix also uses the SignalP 6 installed in it,
  which is licensed to the person who installed it.
- **On a cluster**, the prefix must be on storage the compute nodes see, under the same path. See
  [Shared storage](../reference/running.md#cluster-storage).

Continue with [Reference databases](databases.md) and
[Licensed and external tools](external-tools.md), then [Quick start](quick-start.md).
