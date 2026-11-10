# Running FunFlux

This page describes how to launch FunFlux and how to manage a run. It covers the dry run, CPU
cores and memory, conda environments and run times. It also explains what to do when a run stops
or a step must run again, how to run several batches or use a cluster, and what to remove before
you share an output folder. Read it before the first run.

A run starts from the repository root, with a dry run first:

```bash
cd /path/to/FunFlux
cp config/config.yaml config/config_custom.yaml     # once; edit the copy

# check the configuration and list the jobs, without running anything
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml -n

# run
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml
```

The configuration must be named with `--configfile`, and without it the run stops with
`[FunFlux] No configuration supplied.` FunFlux has no default configuration, so nothing missing
from your copy is filled in from the template. `config/config_custom.yaml` is ignored by git, so a
`git pull` leaves it alone. Every key is described on [Configuration](configuration.md).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **launch folder** | the folder Snakemake is started from, normally the FunFlux repository; relative paths are read from it, and `.snakemake/` is created there |
    | **dry run** | a Snakemake run with `-n`: FunFlux checks the configuration and the inputs, and Snakemake lists the jobs it plans without running any |
    | **rule** | one step of the workflow as Snakemake defines it, such as `select_contigs`; the console and the log files in `logs/` use its name |
    | **job** | one run of a rule, for one sample or once for the whole run |
    | **lock** | the record, in `.snakemake/`, of the files a run in progress reads and writes; a killed run leaves it behind, and `--unlock` removes it |
    | **batch** | the samples one run finds in its mode's input folder; files of one kind share one extension |
    | **heavy step, long step, light step** | on a cluster, heavy steps ask the scheduler for `ram_gb` and `runtime_hours`, long steps for `runtime_hours` only, and light steps for neither |
    | **login node, compute node** | on a cluster, the computer you log in to, which usually has internet, and the machines the scheduler runs the jobs on, which may have none |

Launch from the repository root every time, because three things depend on the launch folder:

- `workflow/Snakefile` is found there without `--snakefile`.
- `.snakemake/` is created there. It holds Snakemake's record of how each output was made, its
  locks, the markers of unfinished jobs, its run logs and, by default, the conda environments.
- Relative paths in the configuration are read from there, by Snakemake and by every job.

FunFlux treats the paths in the configuration in three ways:

| Paths | A relative path is |
|---|---|
| `directories.output_dir`, the database folders, `files.iprscan`, `files.annotation_params`, the protein paths in the annotation parameter table | turned into an absolute path when FunFlux starts |
| the `input.*_dir` folders, the decontamination tables | used as given |
| `directories.genemark_dir`, `files.signalp6_archive` | not allowed: these must be absolute |

With absolute paths throughout, a configuration works from any launch folder.

## Dry run

A dry run shows whether a configuration is ready, without running anything. `-n` (`--dry-run`)
runs every FunFlux check, builds the job list and stops. FunFlux checks the configuration before
Snakemake plans any job, so a dry run takes seconds:

| Problem | What the dry run does |
|---|---|
| a wrong mode, a missing input folder, a bad sample name, a missing row in the annotation table, a mistyped key; when the contamination and mitochondrion screen runs, a decontamination table that is missing or malformed, or a sample in `include` mode without a genus | the first failed check prints a `[FunFlux]` message and stops ([Troubleshooting](../troubleshooting.md#stops-before-any-job-runs)) |
| a missing eggNOG or InterProScan path (these are only tested for existence) | prints a `[FunFlux] NOTE:`; the run fails later, in that rule |
| a missing `directories.blast_db`, with the screen on | prints the same NOTE, then stops with `MissingInputException in rule blob_json`: the screen needs `nodes.dmp` and `names.dmp` from that folder |

!!! warning "A missing eggNOG or InterProScan path does not stop the dry run"

    Most configuration mistakes stop the dry run. These two paths only get a NOTE, and the run
    fails when the rule that reads them starts.

A dry run first prints the FunFlux banner and the settings FunFlux read, with any
`[FunFlux] NOTE:` lines among them:

```text
Mode: contigs — Annotation of fungal genome assemblies you already have.
[FunFlux] NOTE: contigs mode without parameters.contigs.screen: your FASTA gets no contamination or mitochondrion screen. ...
Decontamination mode: off (discard_no_hit=false).
GeneMark-ES: off (funannotate runs its no-GeneMark branch).
Tiberius: on (its gene models join the others in EVidenceModeler, weight 1).
SignalP 6: off.
Contig summary: telomere motif auto (the repeat unit is chosen from the data), a telomere from 5 tandem copies within 25000 bp of a contig end, an rDNA array from 3 rDNA unit(s).
AMR screen: off.
Sample SAMPLE1 will be processed.
AUGUSTUS: SAMPLE1 is trained on its BUSCO genes (checked once the database is unpacked).
```

Then comes Snakemake's plan:

- `Conda environment workflow/envs/<name>.yaml will be created.` for each environment not built
  yet
- the `Job stats` table
- one block per job, with its `reason:`
- `Would remove temporary output <path>` for each temporary file
- last, `This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.`

When the output is redirected to a file, FunFlux's lines can land after Snakemake's.

Check the `Job stats`. With Tiberius on, and GeneMark-ES, SignalP 6 and the AMR screen (the
resistance-mutation screen) off, the first run of one sample plans:

| Mode | Per sample | Once per run | First run, one sample |
|---|---|---|---|
| `illumina` | 34 | 12 | 46 |
| `longread` | 34, plus `find_ont_adapter`, `trim_ont_adapters` and `check_ont_trimming` for an ONT sample (unless `ont_adapter` is `off`), and `check_medaka_model` and `long_read_consensus` for an ONT sample with Medaka | 10 | 44 (PacBio HiFi) or 49 (ONT with Medaka) |
| `hybrid` | 49, plus the same five for an ONT sample | 12 | 61 (PacBio HiFi) or 66 (ONT with Medaka) |
| `contigs` | 21; the screen adds 8 | 10 | 31; 39 with the screen |

`illumina` and `hybrid` add `download_phix` and `build_phix` to the jobs that run once per run.
The other once-per-run jobs are `normalize_funannotate_db`, `antismash_db` or
`antismash_db_local`, `unite_db` or `unite_db_local`, `telomere_db`, `univec_illumina`, `multiqc`
and `all`. Six
more are added only when their condition holds. The counts in the table include the first three,
which hold on a first run with the default settings; with the funannotate database already in
place, `funannotate_db` drops out and every count is one lower:

- `funannotate_db`, only while the database is missing
- `unite_identified`, unless `parameters.its.unite_identified_only` is `false`
- `tiberius_model`, with Tiberius on
- `install_signalp6`, with SignalP 6 on
- `chroquetas_db`, with the AMR screen on
- `busco_data`, for a named BUSCO lineage

The optional steps change these counts in every mode:

| Setting | Jobs per sample | Once per run |
|---|---|---|
| Tiberius off | 2 fewer (a `contigs` sample without the screen: 19 instead of 21) | `tiberius_model` left out |
| GeneMark-ES on | 1 more | |
| SignalP 6 on | 1 more | `install_signalp6` added |
| AMR screen on | 1 more (`amr_chroquetas`) | `chroquetas_db` added |
| `parameters.its.unite_identified_only: false` | | `unite_identified` left out |
| `parameters.busco.lineage` names a lineage | | `busco_data` added |
| `parameters.contigs.screen: true` (`contigs` mode only) | 8 more | |
| `parameters.<mode>.ont_adapter: off` (ONT samples in `longread` and `hybrid`) | 3 fewer | |

For example, the first run of one `contigs` sample without the screen plans 31 jobs, and 33 with
the AMR screen on. The two-sample examples of `config/config_test_examples.yaml` plan 83 jobs in
`longread` mode, with one ONT and one PacBio sample, and 115 in `hybrid` mode. With the AMR screen
on, they plan 86 and 118. On a finished output folder, a dry run prints
`Nothing to be done (all requested files are present and up to date).`

## CPU cores and threads

`--cores N` is the only CPU setting. Give it the same number as `resources.threads`.

Snakemake starts a job only when the threads it asks for fit into what is left of N, and keeps
them reserved until the job ends. Each step asks for what its tool can use. The assemblers, the
read mappings, Tiberius and funannotate predict take every thread, while the other multi-threaded
tools are capped, so that several run side by side. In the table, "up to 12" means 12, or
`resources.threads` when that is smaller:

| Threads asked for | Steps (rule names as in the console) |
|---|---|
| `resources.threads` | SPAdes (`illumina_assembly`), hifiasm or Flye (`long_read_assembly`), Tiberius (`tiberius_prediction`), funannotate predict (`funannotate_prediction`); the read mappings: PhiX removal (`map_phix`), the coverage track (`map_contigs`) and, in `hybrid`, the guide mapping and the read selection (`map_guide_contigs`, `map_sel_contigs`) |
| up to 24 | BLAST (`blast_contigs`, `blast_guide_contigs`), QUAST (`genome_assembly_evaluation`), BUSCO (`completeness`), Qualimap (`map_evaluation`, `map_guide_evaluation`), bwa mem before Polypolish (`short_read_correction`; Polypolish itself runs on one thread) |
| up to 16 | fastp (`trim_adapters`), fastplong (`trim_ont_adapters`) |
| up to 12 | Medaka (`long_read_consensus`), RepeatModeler and RepeatMasker (`repeat_masking`), antiSMASH (`antismash_annotation`), eggNOG-mapper (`eggnog_annotation`), funannotate annotate (`funannotate_annotation`) |
| up to 10 | InterProScan (`iprscan_annotation`) |
| up to 8 | nhmmer (`its_windows`), ITSx (`its_extraction`), NanoPlot (`raw_long_read_qc`, `filtered_long_read_qc`) |
| `parameters.signalp6.torch_num_threads` (default 8), at most `resources.threads` | SignalP 6 (`signalp6_annotation`) |
| up to 4 | ChroQueTas, which hands them to miniprot and MAFFT (`amr_chroquetas`) |
| up to 2 | funannotate clean (`funannotate_preprocess`) |
| 1 | every other step, among them SINTAX (`its_classification`), Filtlong (`filter_long_reads`), the bowtie2 indexes (`build_phix`, `index_contigs`, `index_guide_contigs`, `index_selected_contigs`), the Medaka model check (`check_medaka_model`) and the GeneMark-ES licence test (`genemark_check`) |

- With `--cores` equal to `resources.threads`, a step that asks for every thread runs alone.
  InterProScan holds its 10 threads for hours, and eggNOG-mapper (12) and the ITS steps (8, 8, 1)
  are capped so that they can run beside it.
- In `hybrid`, the guide BLAST and the long-read BLAST of a sample each ask for 24 threads, so
  with `--cores 24` they run one after the other.
- At the start, Snakemake prints `Provided cores: N` and
  `Rules claiming more threads will be scaled down.`, which means that a step asking for more
  than N gets N.
- Two steps use a little more than they ask for. RepeatMasker's `-pa` counts parallel batches,
  and each batch can run several RMBlast threads. `funannotate clean` runs 2 minimap2 searches at
  a time, and each starts 3 threads.
- A `parameters.signalp6.torch_num_threads` above `resources.threads` is lowered to
  `resources.threads`, with a NOTE.

!!! warning "Do not add `--jobs` on a single computer"

    For a local run, Snakemake treats `--jobs` (`-j`) as another name for `--cores`. Passing both
    can let the rules use more CPUs than the machine has, so use `--cores` alone. On a cluster it
    is the other way round: `--executor slurm` needs `--jobs` (see [On a cluster](#on-a-cluster)).

Snakemake does not limit memory on its own. `resources.ram_gb` reaches two steps:

| Step | Gets `resources.ram_gb` as | To cap how many run at once |
|---|---|---|
| SPAdes (`illumina_assembly`) | its memory limit, `-m` | `--resources mem_gb=N` |
| Qualimap (`map_evaluation`; in `hybrid` also `map_guide_evaluation`) | its Java heap: the smaller of `ram_gb` and 64 GB | `--resources java_mem=N` |

Snakemake passes the value on but limits nothing unless you ask it to. With
`--resources mem_gb=N`, it starts only as many SPAdes jobs as fit into N GB at the same time.
Changing `ram_gb` re-runs nothing. On a cluster, `ram_gb` is also the memory every heavy step asks
the scheduler for, as described under [How the settings reach Slurm](#cluster-resources).
InterProScan's Java heap is set in its own `interproscan.sh` (see
[Troubleshooting](../troubleshooting.md#interproscan)).

## Conda environments

`--sdm conda`, short for `--software-deployment-method conda`, makes Snakemake build one conda
environment per tool from `workflow/envs/<name>.yaml` the first time a rule needs it, and run the
rule inside it. The exact package builds come from `workflow/envs/<name>.linux-64.pin.txt` next to
each file, captured from tested environments. Besides Snakemake, only InterProScan and the
licensed tools are installed by hand, as described on
[Installation](../getting-started/installation.md) and
[Licensed and external tools](../getting-started/external-tools.md).

Building needs network access, once per environment. A `contigs` run builds 15 environments by
default and a `hybrid` run up to 27. The full list is on
[Installation](../getting-started/installation.md#sdm-conda).

The environments are stored in `.snakemake/conda/` of the launch folder. `--conda-prefix` stores
them elsewhere, so that several clones or launch folders share one set:

```bash
snakemake --sdm conda --conda-prefix /path/to/funflux_envs --cores 24 \
    --configfile config/config_custom.yaml
```

The environments are built once and reused. The table lists two options that act on them, and
what happens when they move or are deleted:

| Option or event | Effect |
|---|---|
| `--conda-create-envs-only` | builds the environments and exits ([Building the environments in advance](../getting-started/installation.md#build-in-advance)) |
| `--list-conda-envs` (with `--sdm conda --configfile ...`) | prints the folder of each environment, for example to run `funannotate setup` inside the funannotate one with `conda activate <folder>` |
| another `--conda-prefix`, or the prefix moved, for runs from the same launch folder | new environments, and the finished steps run again: Snakemake records where each step's environment was, and a step whose environment now sits elsewhere runs again (`Software environment definition has changed since last execution`), with every step after it |
| a deleted `.snakemake/conda/` or prefix folder | new environments in the same place; finished steps are kept |
| either of the two rows above, with SignalP 6 on | SignalP 6 was installed into the old signalp6 environment and has to be installed again: delete `00.resources/signalp6/installation.json` in the output folder ([Troubleshooting](../troubleshooting.md#signalp-6)) |

!!! warning "Keep one `--conda-prefix` for an output folder"

    A new place for the environments means building them, and installing SignalP 6, again. From
    the same launch folder, it also re-runs every finished step.

Every job runs with `PYTHONNOUSERSITE=1`, so packages installed with `pip install --user` cannot
be used instead of those of an environment.

## Run time

One fungal genome in `contigs` mode takes several hours:

| Step | Time for one fungal genome |
|---|---|
| InterProScan | usually a few hours for a fungal proteome |
| Tiberius | about half an hour for a 40 Mb genome, with about 16 GB of memory |
| GeneMark-ES licence test | about one minute |
| repeat masking with tantan | seconds |
| repeat masking with `masking_method: repeatmodeler_repeatmasker` | far longer than with tantan ([Repeat masking](../analysis/repeat-masking.md)) |
| most other steps | minutes |

Tiberius, gene prediction, InterProScan and funannotate annotate run one after the other, so
these steps set the length of a run. Because eggNOG-mapper is capped at 12 threads, it can run
beside InterProScan with `--cores` 22 or more. With fewer, one of the two waits for the other. The
read-based modes add read processing and assembly before this chain, and polishing where it runs
(`hybrid`, and ONT samples with Medaka in `longread`).

Whole runs from reads, one sample each, took 6.6 to 11.5 hours on 16 to 32 cores, on two public
genomes ([Validation on public genomes](../about/validation.md)):

| Genome | Mode and reads | Assembler | Cores | Run time (h) |
|---|---|---|---|---|
| *Botrytis cinerea* B05.10, 42.6 Mb | `illumina`, Illumina 79x | SPAdes | 24 | 7.9 |
| | `longread`, ONT 136x | hifiasm | 28 | 8.9 |
| | `longread`, ONT 136x | Flye | 20 | 8.4 |
| | `longread`, ONT 42x, with Medaka | hifiasm | 16 | 9.2 |
| | `hybrid`, ONT 136x plus Illumina 79x | hifiasm | 32 | 11.5 |
| *Cordyceps cateniannulata* GXU-8616, 32.7 Mb | `longread`, PacBio HiFi 85x | hifiasm | 24 | 6.6 |
| | `longread`, PacBio HiFi 85x | Flye | 20 | 7.3 |
| | `hybrid`, PacBio HiFi 85x plus MGISEQ short reads 2 x 150 bp, 151x | hifiasm | 24 | 7.8 |

The runs shared one server with other runs going on at the same time, so the times are indicative
only.

- **Run time** is that of the main Snakemake run, from the first to the last time stamp of its log.
- **The Flye run on GXU-8616** took two attempts, because InterProScan failed once and the run was
  started again by hand. Its run time, 7.3 h, is the sum of the two attempts, 3.7 and 3.6 h
  ([Validation](../about/validation.md#hifi-assembly)).

## Resuming after an interruption

Start the same command again from the same launch folder. Snakemake skips every job whose outputs
exist and are up to date and runs the rest. The launch folder matters, because its `.snakemake/`
records which outputs were being written when the run stopped.

| Snakemake prints | What to do |
|---|---|
| `Error: Directory cannot be locked. ...` | A killed run (kill signal, power loss) left its lock. If no other Snakemake process runs from this folder, run the same command with `--unlock`, then the command again. |
| `The files below seem to be incomplete. ...` | Jobs were running when the run stopped; their outputs may be half-written. Add `--rerun-incomplete` (`--ri`): those jobs start again from the beginning. |
| `Exiting because a job execution failed. Look below for error messages` | A job failed and Snakemake deleted its outputs (`Removing output files of failed job ...`). The block `Error in rule <rule>:` names the log; fix the cause and run again ([Troubleshooting](../troubleshooting.md)). |

Run the same command again, and add `--unlock` or `--rerun-incomplete` only when Snakemake prints
the matching message:

```bash
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml --unlock
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml --rerun-incomplete
```

- `--unlock` needs `--configfile` as well, because FunFlux reads the configuration before
  Snakemake removes the lock, and stops without one.
- `--keep-going` (`-k`) lets the jobs that do not depend on a failed job carry on, so the other
  samples of a batch finish. It covers failures of running jobs, not `[FunFlux]` stops.
- A resumed run does not download or install again what is already in `00.resources/`. That
  covers the UNITE reference, the antiSMASH database, the telomere repeat database, the SignalP 6
  installation, ChroQueTas and its database, the data of a named BUSCO lineage and the repaired
  funannotate database paths. The funannotate database is not downloaded again either: it is
  downloaded only while `funannotate-db-info.txt` is missing from `directories.funannotate_db`.
  With `parameters.busco.lineage: auto`, BUSCO downloads its lineage datasets again in every
  BUSCO job.
- Snakemake deletes the log of a job when it starts that job again, so copy the log of a failed
  job before you resume if you need it.

## What makes Snakemake re-run a rule

A dry run lists every job Snakemake plans to run again, each with its reason. With its default
`--rerun-triggers code input mtime params software-env`, Snakemake runs a job again when:

| What happened | The `reason:` a dry run prints |
|---|---|
| an output is missing | `Missing output files` |
| an input file is newer than the output; for an input of at most 1,000,000 bytes only if its content changed as well, so touching a small file triggers nothing | `Updated input files` |
| a job before it runs again; a re-run spreads to every later step of that sample | `Input files updated by another job` |
| a value FunFlux passes to the rule changed | `Params have changed since last execution` |
| the rule's command changed, for example after updating FunFlux | `Code has changed since last execution` |
| its environment file changed, or its environment now sits in another folder (another `--conda-prefix`) | `Software environment definition has changed since last execution` |
| its list of input files changed, for example `multiqc` when the samples of the batch change | `Set of input files has changed since last execution` |

A changed setting re-runs the first rule that reads it, and every later step follows. The
settings FunFlux passes to rules, and the first rule each one re-runs:

| Changed in the configuration | First rule that runs again |
|---|---|
| `directories.blast_db`, `parameters.nt_version`, `blast_max_target_seqs`, `blast_max_hsps` | `blast_contigs` |
| `parameters.decontamination.*` | `select_contigs` |
| `parameters.mitochondria.*` | `screen_mitochondria` |
| `parameters.busco.lineage` | `busco_data` when a named lineage is not in `00.resources/busco/` yet, then `completeness` |
| `parameters.contig_summary.telomere_motif`, `telomere_window`, `min_telomere_copies` | `telomere_units`, then `contig_summary` and `multiqc`; nothing else |
| `parameters.contig_summary.min_rdna_units` | `contig_summary`, then `multiqc`; nothing else |
| `links.telomere_db_link`, `parameters.contig_summary.telomere_db_sha256` | `telomere_db`, then `telomere_units` of every sample, `contig_summary` and `multiqc` |
| `links.univec_link`, `parameters.contig_summary.univec_sha256` | `univec_illumina`, then `blast_adapters` and `screen_adapters` of every sample, `contig_summary` and `multiqc` |
| `parameters.its.window_flank` | `its_windows` |
| `parameters.its.max_full_length`, `parameters.its_taxonomy_cutoff` | `its_classification` |
| `parameters.its.unite_identified_only` | `its_classification` (its UNITE input changes); switched on, `unite_identified` runs first if its files are missing |
| `directories.unite_db`, `links.unite_its_link` | `unite_db_local` or `unite_db` |
| `links.funannotate_link` | `funannotate_db`, if a run from this launch folder downloaded the database: the archive the new link names is downloaded unless present, checked and unpacked into `directories.funannotate_db` |
| `directories.funannotate_db` | `normalize_funannotate_db` |
| `directories.antismash_db` | `antismash_db_local` or `antismash_db` |
| `parameters.masking_method`, `repeatmodeler_quick`, `repeatmodeler_ltrstruct`, `repeatmodeler_seed` | `repeat_masking` |
| `directories.genemark_dir`, `parameters.genemark.mode` | `genemark_check` when GeneMark-ES is on after the change, otherwise `funannotate_prediction` |
| `parameters.prediction.busco_db`, column 3 (proteins) of a sample's row in the annotation table | `funannotate_prediction` |
| column 2 (species) or 4 (AUGUSTUS model) of a sample's row in the annotation table | `augustus_model_check`, then `funannotate_prediction`; for column 2 also `amr_chroquetas` (AMR screen on) |
| `parameters.prediction.augustus_optimize` | `augustus_model_check` of the samples whose column 4 is empty, then their `funannotate_prediction` |
| `parameters.tiberius.mode` | from `off` to `on`: `tiberius_model` (when the model is not there yet), then `tiberius_prediction`; from `on` to `off`: `funannotate_prediction` |
| `links.tiberius_model_link`, `parameters.tiberius.model_sha256` | `tiberius_model`, then `tiberius_prediction` of every sample |
| `files.iprscan` | `iprscan_annotation` |
| `directories.eggnog_db` | `eggnog_annotation` |
| `parameters.signalp6.mode` | `install_signalp6` when SignalP 6 is on after the change, otherwise `funannotate_annotation` |
| `parameters.signalp6.archive_sha256`, `files.signalp6_archive` | `install_signalp6` |
| `parameters.signalp6.batch_size` | `signalp6_annotation` |
| `parameters.amr.chroquetas` | from `off` to `on`: `chroquetas_db` (when its files are not there yet), then `amr_chroquetas` of every sample and `multiqc`; from `on` to `off`: `multiqc` only, and `05.amr/` stays on disk |
| `parameters.amr.mito_genetic_code` | `amr_chroquetas`, then `multiqc`; nothing else |
| `links.chroquetas_link`, `parameters.amr.chroquetas_script_sha256`, `parameters.amr.chroquetas_db_sha256` | `chroquetas_db`, then `amr_chroquetas` of every sample and `multiqc` |
| `links.phix_link` | `download_phix`, then `build_phix` and every `map_phix` |
| `parameters.spades_kmers` | `illumina_assembly` (in `hybrid`: the guide assembly) |
| `parameters.<mode>.ont_adapter` | `find_ont_adapter` of every ONT sample, then the trimming, Filtlong and the whole assembly; switched off, `filter_long_reads` reads the raw reads, which re-runs the same chain |
| `parameters.<mode>.filtlong.*`, `parameters.hybrid.short_read_guidance`, `parameters.hybrid.filtlong_split` | `filter_long_reads` |
| `parameters.<mode>.assembler`, `flye_input_mode`, `hifiasm.*` | `long_read_assembly` |
| `parameters.<mode>.medaka_model`, `medaka_model_fallback_auto` | `check_medaka_model`; `long_read_assembly` of every ONT sample waits for its file and runs again too. Switching Medaka off or on also changes which file `finalize_assembly` (`longread`) or `short_read_correction` (`hybrid`) copies or corrects |

A new BUSCO lineage, for example, re-runs BUSCO, the contig summary and the report, and a new
`directories.funannotate_db` re-runs the whole annotation of every sample.

Other changes re-run nothing, except where a changed file is itself an input of a step:

| Changed | Re-runs |
|---|---|
| `resources.threads`, `resources.ram_gb`, `resources.runtime_hours`, `--cores`, `parameters.signalp6.torch_num_threads` | nothing |
| the content of the databases and of the protein evidence files (only their paths are passed) | nothing, except for the files in the next row |
| a newer copy of `nodes.dmp` or `names.dmp` in `directories.blast_db`, `funannotate-db-info.txt`, the SignalP 6 archive, the folder named in `directories.antismash_db` (its own modification time, which changes when an entry is added or removed), the `unite_its_sintax.fasta` that `directories.unite_db` points to (Snakemake follows the link) | the steps that read it |

!!! note "Updating a database in place re-runs nothing, unless a file in the last row becomes newer"

    Only the paths of the databases and of the protein evidence files are passed to the rules. A
    change to their content therefore re-runs nothing, unless one of the files in the last row of
    the table becomes newer.

To redo a step anyway:

```bash
# this rule for every sample, and every step after it
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml \
    --forcerun iprscan_annotation
```

To redo a step for one sample only, delete that step's output, for example
`04.annotation/{sample}/iprscan/interproscan_out.xml`. A deleted final result is rebuilt on the
next run, with every step after it. A deleted intermediate file that no pending job needs is not
rebuilt. The final results are the files every run makes sure exist:

| Mode | Final results |
|---|---|
| every mode | `eval/quast/`, `eval/busco/`, `eval/contig_summary.tsv`, the seven ITS files, `contigs_clean.fasta`, `contigs_sort.fasta`, `contigs_mask.fasta`, `funannotate/`, `iprscan/interproscan_out.xml`, `eggnog/`, `antismash/`, `signalp6/prediction_results.txt` (SignalP 6 on) and `.cleanup_complete.txt` of every sample, `05.amr/{sample}/amr_summary.tsv` (AMR screen on), plus `00.resources/funannotate_db_paths.normalized.json` and `06.report/multiqc_report.html` |
| `contigs` | also `contigs_final.fasta`; with the screen also `contigs_mitochondrial.fasta` and `mitochondria/mitochondrial_screen.tsv` |
| `illumina` | also the three genome files, always (`contigs_final.fasta`, `contigs_mitochondrial.fasta`, `mitochondria/mitochondrial_screen.tsv`); `01.reads/{sample}/illumina/{sample}_fastp.json` and `eval/qualimap/` |
| `longread` | also the three genome files; `01.reads/{sample}/long/{sample}_filtered.fastq`, `raw_qc/`, `filtered_qc/`, `assembly_provenance.tsv`, `eval/qualimap/` and, for ONT samples with Medaka, `medaka/consensus.fasta` |
| `hybrid` | also the three genome files; all of the `longread` files plus the fastp JSON, the two selected read files, `polypolish/polished.fasta` and `guide/eval/qualimap/` |

`--list-changes params` (or `code`, `input`) lists the output files whose parameters, code or
inputs changed since they were made.

## Running several batches

A run processes every sample in the input folder of its mode, and those samples form one batch.
All files of one kind in a batch share one file extension, but in `hybrid` the Illumina files and
the long-read files may have different ones. The annotation table may hold rows for samples of
other batches, and FunFlux ignores the rows without a sample in the batch.

- **One output folder per batch** gives each folder its own `00.resources/`. The UNITE reference
  (207 MB), the antiSMASH database and a named BUSCO lineage are downloaded into each, and
  SignalP 6 is installed into each (1.6 GB). Point `directories.unite_db` and
  `directories.antismash_db` at shared copies to download them once, and give every run the same
  `--conda-prefix`.
- **Several batches can share one output folder** when their sample names differ. `00.resources/`
  is reused, and `06.report/multiqc_report.html` is rewritten with the samples of the latest run
  only. A sample name that comes back with another input file is processed again, and its earlier
  results are replaced.
- **Runs at the same time** each count only their own `--cores`, and the sum must fit the machine.
  From one launch folder, Snakemake refuses a second run (`Directory cannot be locked`) only when
  both would write the same file, or one would write a file the other reads. Runs from different
  launch folders are not checked against each other.
- **Runs that name the same `directories.funannotate_db`**: start one, let it pass
  `normalize_funannotate_db`, then start the others. Every run must reach the folder under the
  same real path (see [Reference databases](../getting-started/databases.md#what-the-workflow-does)).

## On a cluster { #on-a-cluster }

FunFlux runs on a Slurm cluster in one of two ways. In run A, Snakemake spreads the samples and
steps over many nodes and keeps every download on the login node. Run B is a single-computer run
inside one batch job:

| | A. Snakemake sends each step to Slurm | B. The whole run in one batch job |
|---|---|---|
| Where Snakemake runs | on the login node, or another node with internet, for the whole run | inside one batch job, on one compute node |
| Where the steps run | each step is a Slurm job of its own, on any compute node; the download steps and `augustus_model_check` run where Snakemake runs | all on that node |
| Command | `--executor slurm --jobs N` ([Run A](#cluster-run-a)) | `--cores N`, as on a single computer ([Run B](#cluster-run-b)) |
| Memory and time | each heavy step asks Slurm for `resources.ram_gb` and `resources.runtime_hours` | the batch job asks once, for the whole run |
| Internet on the compute nodes | not needed, except by BUSCO with `parameters.busco.lineage: auto` | not needed once the [preparation](#cluster-preparation) has run, except by BUSCO with `lineage: auto` |
| Needs | the Slurm executor plugin; a login node on which Snakemake may run for many hours | one node with the cores, memory and time of the whole run |
| Checked | with Snakemake's generic cluster executor and a submit command that only records each job; not yet on a Slurm cluster | as on a single computer |

Both start with the same preparation, and both need the same [shared storage](#cluster-storage).

### Preparation, on a computer with internet { #cluster-preparation }

The preparation builds the environments, fetches what FunFlux downloads itself and checks the
result, from the folder and with the `--conda-prefix` of the real run. Run these commands from the
FunFlux folder on a computer that has internet and sees the same storage as the compute nodes:
the login node or a data-transfer node.

1. Copy the template and fill it in. Put every path on [shared storage](#cluster-storage), and
   set `resources.threads`, `ram_gb` and `runtime_hours` from the limits of the cluster
   ([table below](#cluster-resources)). When the compute nodes have no internet, also name a
   [BUSCO lineage](#cluster-busco). The command that lists the lineage names runs in the busco
   environment, which step 2 builds, so you can name the lineage after step 2.

    ```bash
    conda activate snakemake
    cd /path/to/FunFlux
    cp config/config.yaml config/config_custom.yaml
    ```

2. Build the conda environments of every step of the run. To build them, conda first downloads
   every package into its package cache, by default inside the conda installation. If that sits
   in a home folder with a small quota, move the cache to project space first, as in the first
   line:

    ```bash
    export CONDA_PKGS_DIRS=/path/to/project_space/conda_pkgs
    snakemake --sdm conda --conda-prefix /path/to/project_space/funflux_envs \
        --conda-create-envs-only --configfile config/config_custom.yaml
    ```

3. Fetch everything FunFlux downloads itself, and nothing else. Add `-n` first to see the plan.

    ```bash
    snakemake download_all --sdm conda --conda-prefix /path/to/project_space/funflux_envs \
        --cores 4 --configfile config/config_custom.yaml
    ```

4. Check. The dry run must print
   `Nothing to be done (all requested files are present and up to date).`

    ```bash
    snakemake download_all -n --sdm conda --conda-prefix /path/to/project_space/funflux_envs \
        --configfile config/config_custom.yaml
    ```

    Two more lines can follow it, `1 jobs have missing provenance/metadata ...` and
    `Rules with missing metadata: funannotate_db`. Both are harmless: the funannotate database
    was already in its folder, downloaded by hand or by a run from another launch folder, so
    Snakemake holds no record of downloading it.

Four points about these commands, each with its reason:

| Point | Why |
|---|---|
| Write `download_all` right after `snakemake`. | `--configfile` and `--sdm` take one or more values. At the end of the line `download_all` is read as a second configuration file, and Snakemake stops with `No such file or directory: 'download_all'`. |
| Run steps 2 to 4 and the real run from the same folder, with the same configuration and the same `--conda-prefix`. | Relative paths in the configuration are read from the launch folder, and Snakemake keeps its records in `.snakemake/` there. Without `--conda-prefix` the environments are in that `.snakemake/conda/` too, and Medaka keeps the model weights the check downloads inside its environment. From another folder, the real run builds the environments again, without the weights. With another prefix it also runs every download step again: Snakemake records where each step's environment was (reason `Software environment definition has changed since last execution`), and step 4 plans the downloads instead of printing `Nothing to be done`. |
| Run steps 2 to 4 again after any change to the configuration or to the input folders. | A new sample, another Medaka model, the AMR screen switched on or a lineage named each add a download. |
| Run step 4 even when step 3 went well. | FunFlux checks `parameters.prediction.busco_db` against the funannotate database only once the database is on disk. A set the newly downloaded database lacks stops step 4 with the command that adds it, and that command needs internet too ([Troubleshooting](../troubleshooting.md#databases-and-links)). |

`download_all` covers every download of the run except BUSCO with `lineage: auto`, and with that
setting it prints a note that every BUSCO job will need internet. What it fetches depends on the
configuration:

| Download | Rule | When | Where it lands |
|---|---|---|---|
| PhiX genome | `download_phix` | `illumina`, `hybrid` | `00.resources/phix/` |
| funannotate database, then the repair of the paths stored in it | `funannotate_db`, `normalize_funannotate_db` | the download only while `directories.funannotate_db` holds no database | `directories.funannotate_db` |
| UNITE reference | `unite_db` | unless `directories.unite_db` names your copy | `00.resources/unite/` |
| telomere repeat database | `telomere_db` | always | `00.resources/telomeres/` |
| Illumina adapter records of NCBI UniVec | `univec_illumina` | always | `00.resources/adapters/` |
| Tiberius `fungi` model | `tiberius_model` | Tiberius on, the default | `00.resources/tiberius/` |
| ChroQueTas and its FungAMR database | `chroquetas_db` | AMR screen on | `00.resources/chroquetas/` |
| antiSMASH databases (9.4 GB) | `antismash_db` | unless `directories.antismash_db` names your copy | `00.resources/antismash/` |
| BUSCO lineage data | `busco_data` | `parameters.busco.lineage` names a lineage | `00.resources/busco/` |
| Medaka model weights | `check_medaka_model`, which reads the first reads of each sample, so the reads must be in place | ONT samples polished with Medaka | inside the `medaka` environment, or `~/.medaka/data` when that cannot be written |

The large databases and the licensed tools stay yours to install. You provide them yourself and
name them by path, as described on [Reference databases](../getting-started/databases.md) and
[Licensed and external tools](../getting-started/external-tools.md):

| What | Key | Needed |
|---|---|---|
| NCBI `core_nt` and its taxonomy files | `directories.blast_db` | with the screen on |
| eggNOG-mapper data | `directories.eggnog_db` | always |
| InterProScan | `files.iprscan` | always |
| GeneMark-ES and its key | `directories.genemark_dir`, `~/.gm_key` | GeneMark-ES on |
| SignalP 6 archive | `files.signalp6_archive` | SignalP 6 on |

### Run A: Snakemake sends each step to Slurm { #cluster-run-a }

Install the Slurm executor plugin once into the launcher environment, Snakemake's own conda
environment:

```bash
conda activate snakemake
conda install -c conda-forge -c bioconda snakemake-executor-plugin-slurm
```

Then, on the login node, start Snakemake inside `tmux` or `screen`, so that it keeps running when
you log out:

```bash
snakemake --executor slurm --jobs 20 --latency-wait 60 --sdm conda \
    --conda-prefix /path/to/project_space/funflux_envs \
    --default-resources slurm_account=<account> slurm_partition=<partition> mem_mb=8000 runtime=240 \
    --configfile config/config_custom.yaml
```

| Option | What it does |
|---|---|
| `--executor slurm` | submits every step to Slurm with `sbatch`, except the local ones: the download steps and `augustus_model_check` run on the login node |
| `--jobs 20` | the most jobs Snakemake keeps queued or running at once; it refuses to start without `--jobs` |
| `--latency-wait 60` | after a job ends, Snakemake waits up to 60 seconds for its files to appear on the login node. Shared storage can show a file written on a compute node only some seconds later, and the default wait, 5 seconds, is often too short |
| `--default-resources` | the values for every job that does not set its own: the account and partition for every job (leave them out if your cluster needs none), memory (MB) and time (minutes) for the light steps. A rule's own value wins over it |

`--cores` is not needed, because `resources.threads` already caps the CPUs of each job. Given
together with `--executor`, it lowers the threads of every job to at most that number, and so what
each job asks Slurm for.

Each Slurm job writes a log of its own, `<Slurm job id>.log`, under
`.snakemake/slurm_logs/rule_<rule>/` in the launch folder. For the steps that run per sample, the
log is one folder further down, in a folder named after the sample. When a job fails, the
`Error in rule` block names that file. The plugin deletes the logs of jobs that succeeded, and any
log older than 10 days. The messages of the tool itself are in `logs/` of the output folder, as on
a single computer.

#### How the settings reach Slurm { #cluster-resources }

What each job asks Slurm for comes from three settings of the configuration and from
`--default-resources`:

| Setting | Each job asks Slurm for | Must not exceed | Shown by |
|---|---|---|---|
| `resources.threads` | `--cpus-per-task`: the threads of its step, at most this number ([CPU cores and threads](#cpu-cores-and-threads)) | the CPUs of one compute node | `sinfo -o "%n %c"` |
| `resources.ram_gb` | `--mem`: `ram_gb` × 1024 MB for every heavy step (table below) | the memory Slurm offers one compute node | `sinfo -o "%n %m"`, in MB: divide by 1024 and round down; it is often a little less than the size the node is sold with |
| `resources.runtime_hours` | `-t`: `runtime_hours` × 60 minutes for every heavy and every long step | the time limit of the partition | `sinfo -o "%P %l"` |
| `--default-resources mem_mb=... runtime=...` | `--mem` and `-t` for the light steps | the same limits | the same commands |

The memory and time a step asks for depend on the kind of step:

| Steps | Memory asked for | Time asked for |
|---|---|---|
| the heavy steps: SPAdes (`illumina_assembly`), hifiasm or Flye (`long_read_assembly`), Medaka (`long_read_consensus`), Polypolish (`short_read_correction`), the read mappings with BAM sorting (`map_contigs`, `map_guide_contigs`), BLAST (`blast_contigs`, `blast_guide_contigs`), BUSCO (`completeness`), `repeat_masking`, Tiberius (`tiberius_prediction`), funannotate predict and annotate, InterProScan, eggNOG-mapper, antiSMASH, SignalP 6, and in `hybrid` Filtlong (`filter_long_reads`), guided or not | `ram_gb` | `runtime_hours` |
| Qualimap (`map_evaluation`, `map_guide_evaluation`) | its Java heap (`ram_gb`, at most 64 GB) plus 4 GB, never more than `ram_gb` | `runtime_hours` |
| the long steps that hold little in memory: the PhiX filter (`map_phix`), the read selection of `hybrid` (`map_sel_contigs`), Filtlong in `longread` (`filter_long_reads`), the ONT adapter scan, trimming and check (`find_ont_adapter`, `trim_ont_adapters`, `check_ont_trimming`) | nothing | `runtime_hours` |
| every other step (seconds to minutes, a few GB at most) | nothing | nothing |

Set `threads`, `ram_gb` and `runtime_hours` so that one job fits one node and one partition. The
largest memory peaks measured on the 42.6 Mb genome of B05.10 were 17 GB for hifiasm, about 16 GB
for Tiberius, about 15 GB for SPAdes and 11 GB for Flye. The other heavy steps were not measured.
A `ram_gb` well below the node's memory lets several heavy jobs share a node.

- A step that asks for nothing gets Slurm's defaults, often a GB or two per CPU and the
  partition's default time. The Slurm plugin then warns `No wall time information given` and
  `No job memory information ('mem_mb' or 'mem_mb_per_cpu') is given - submitting without`. If
  the defaults are smaller than a few GB or an hour, give
  `--default-resources mem_mb=8000 runtime=240` as above.
- Changing `threads`, `ram_gb` or `runtime_hours` re-runs nothing.

### Run B: the whole run in one batch job { #cluster-run-b }

Run B treats the allocated node like a single computer:

| Step | What to do |
|---|---|
| Request | one node with N cores and enough memory: about 16 GB for Tiberius, `resources.ram_gb` (default 128) for SPAdes in `illumina` and `hybrid` mode, and up to 64 GB for Qualimap in the read-based modes ([Memory](../getting-started/installation.md#memory)) |
| Configure | `resources.threads: N`, the number of cores you request |
| Prepare | when the node has no internet: the [preparation](#cluster-preparation) above, on the login node, with the folder and `--conda-prefix` of the batch script |
| Launch | inside the job, from the same folder, the usual command with `--cores N` ([script below](#a-batch-script)) |
| Time limit | allow for the run times [above](#run-time), 6.6 to 11.5 hours from reads, plus the building of the environments and the downloads when they were not prepared. After a job is stopped at its time limit, submit the same script again ([Resuming after an interruption](#resuming-after-an-interruption)); if Snakemake asks for `--unlock`, run the command once with it first |

On one node, Snakemake counts only threads, so `resources.ram_gb` is only the memory cap of SPAdes
and the heap of Qualimap, and `runtime_hours` is ignored.

#### A batch script { #a-batch-script }

The script below is for Slurm, and other schedulers take the same commands with their own
resource lines. Save it as `funflux.sbatch` and submit it with `sbatch funflux.sbatch`.

```bash
#!/bin/bash
#SBATCH --job-name=funflux
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24      # N: the same number as resources.threads
#SBATCH --mem=160G              # more than resources.ram_gb (Request row above)
#SBATCH --time=24:00:00         # Time limit row above

# Make "conda activate" work inside a batch job, then activate Snakemake.
eval "$(conda shell.bash hook)"
conda activate snakemake

cd /path/to/FunFlux
snakemake --sdm conda --conda-prefix /path/to/project_space/funflux_envs \
    --cores "$SLURM_CPUS_PER_TASK" --configfile config/config_custom.yaml
```

`$SLURM_CPUS_PER_TASK` holds the number of cores Slurm gave the job, so `--cores` always matches
the request.

### Shared storage { #cluster-storage }

Both runs need shared storage, which on most clusters means project space plus the home folder:

| Must be on storage every compute node sees, under the same path | Why |
|---|---|
| the FunFlux folder and the launch folder, with its `.snakemake/` | every job reads the workflow and Snakemake's records |
| the conda installation, with the launcher environment (Snakemake's conda environment) in it | in A, each job starts Snakemake again with the same Python, and Snakemake switches to each step's environment with that installation's `activate` script; in B, the batch script activates it |
| the `--conda-prefix` folder | every step runs in its environment from there |
| the output folder | every step reads and writes it |
| every path in the configuration: the input folders, the annotation table and its protein files, the BLAST, eggNOG, funannotate, antiSMASH and UNITE folders, the InterProScan installation, the GeneMark-ES folder, the SignalP 6 archive | the steps read them |
| your home folder: `~/.gm_key` with GeneMark-ES on, and `~/.medaka/data` when Medaka could not write into its environment | GeneMark-ES and Medaka look there, in the home folder of the account the job runs under |

- Put the environments on project space with one `--conda-prefix`, not on a home folder with a
  small quota. A `contigs` run needs 20 to 25 GB, and the read-based modes more
  ([Sharing environments](../getting-started/installation.md#conda-prefix)).
  `export CONDA_PKGS_DIRS=/path/to/project_space/conda_pkgs` moves conda's package cache there too.
- With SignalP 6 on, a prefix that other people use also gives them the SignalP 6 installation,
  which is licensed to you ([SignalP 6](../getting-started/external-tools.md#signalp-6)).
- A compute node can have other base commands than the login node, so run the check under
  [Commands outside the conda environments](../getting-started/installation.md#host-commands)
  inside a job on a compute node.

### BUSCO on compute nodes without internet { #cluster-busco }

BUSCO is the one step that can need internet on the compute nodes. Where its lineage data come
from depends on `parameters.busco.lineage`:

| `parameters.busco.lineage` | Where the lineage data come from | Needs internet |
|---|---|---|
| `auto`, the default | downloaded inside every BUSCO job, on the compute node, because BUSCO picks the lineage only after it has looked at the genome | on the compute nodes |
| a lineage, for example `hypocreales_odb12` | downloaded once per output folder by `busco_data` into `00.resources/busco/`, where Snakemake runs; every BUSCO job reads it with `--offline` | only where Snakemake runs |

!!! warning "Name a BUSCO lineage when the compute nodes have no internet"

    Choose one that holds every sample of the batch. The report then shows one BUSCO row per
    sample, without the generic eukaryota row (see
    [Assembly QC](../analysis/assembly-qc.md#busco)).

`busco --list-datasets` lists the names. Run it on the login node, in FunFlux's `busco`
environment, after step 2 of the [preparation](#cluster-preparation):

```bash
# Print the folder of each environment; the busco one is on the line workflow/envs/busco.yaml.
snakemake --sdm conda --conda-prefix /path/to/project_space/funflux_envs --list-conda-envs \
    --configfile config/config_custom.yaml
conda activate /path/to/project_space/funflux_envs/<the folder on that line>
busco --list-datasets     # needs internet; leaves a small busco_downloads/ folder here
rm -r busco_downloads
conda activate snakemake
```

The list prints each name with its version after a dot, for example `hypocreales_odb12.2`. Write
the name without the dot and the number after it, as in `lineage: hypocreales_odb12`. BUSCO then
downloads the newest version.

### The GeneMark-ES key { #cluster-genemark }

GeneMark-ES reads its licence key from `~/.gm_key` in the home folder of the account the job runs
under, so the key must be there on every compute node. FunFlux checks that the file exists before
any job starts. In run A every job reads the workflow again, so a compute node that does not see
the key stops that job with `[FunFlux] GeneMark is enabled but there is no licence key at <home>/.gm_key`.
On clusters whose home folder is shared, copying the key once is enough (see
[GeneMark-ES](../getting-started/external-tools.md#genemark-es)).

## Cleaning up

During the run, FunFlux removes the temporary files and the scratch folders of the rules. As the
last step of each sample, rule `cleanup` also removes the `clean_*` folders, `*.log` files and
`genome.fixedproducts` left in `04.annotation/{sample}/`, and lists them in
`logs/cleanup_{sample}.log`. The full list of what FunFlux removes is on
[Output files](output.md#temporary-files-and-folders-a-re-run-replaces).

Before you share or archive an output folder, see [Disk use](output.md#disk-use) for its largest
items. Deleting a final result makes the next run on that folder rebuild it, with every step
after it.

!!! warning "Remove the SignalP 6 package from a copy you share"

    `00.resources/signalp6/package/` (1.6 GB) holds DTU's licensed SignalP 6 model weights, and
    FunFlux does not remove it. In a folder you keep running, a missing package stops the next
    SignalP job with `SignalP installation manifest points to a missing model directory`.

The environments are in `.snakemake/conda/` or in the `--conda-prefix` folder. After deleting that
folder, they are built again on the next run, and SignalP 6 needs a new installation (see
[Conda environments](#conda-environments)).
