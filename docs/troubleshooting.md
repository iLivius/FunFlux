# Troubleshooting

This page lists the messages that stop a FunFlux run or leave a result missing, with the cause of
each and what to change. It follows a run from start to end: the checks of the configuration
before any job starts, the failures of single steps while the workflow runs, and the stops
particular to a cluster.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **dry run** | a Snakemake run with `-n`: FunFlux checks the configuration and the inputs, and Snakemake lists the jobs it plans without running any |
    | **rule** | one step of the workflow as Snakemake defines it, such as `select_contigs`; the console and the log files in `logs/` use its name |
    | **launch folder** | the folder Snakemake is started from, normally the FunFlux repository; relative paths are read from it, and `.snakemake/` is created there |
    | **lock** | the record, in `.snakemake/`, of the files a run in progress reads and writes; a killed run leaves it behind, and `--unlock` removes it |
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **contig selector** | rule `select_contigs`, the last step of the screen: it keeps or removes each contig by genus, as `parameters.decontamination` says; the default, `off`, keeps every contig |
    | **annotation parameter table** | the tab-separated file named by `files.annotation_params`, one row per sample: name, species, a protein FASTA of related fungi and, optionally, an AUGUSTUS model |
    | **heavy step, long step, light step** | on a cluster, heavy steps ask the scheduler for `ram_gb` and `runtime_hours`, long steps for `runtime_hours` only, and light steps for neither |

## Where to look

What you see, and when it appears, tells you where the cause is and which part of this page to
read.

| What you see | Where the cause is |
|---|---|
| A `[FunFlux] ...` line before Snakemake's `Building DAG of jobs...`, and the run ends | A configuration check failed and nothing ran. The message names the key. See [Stops before any job runs](#stops-before-any-job-runs). |
| A `[FunFlux] <step>: ...` line after jobs have started, for example `[FunFlux] funannotate_db: ...` | That step failed or, if it is a report step such as QUAST, it warns and the run goes on. Look the step up under [Failures while the workflow runs](#failures-while-the-workflow-runs). Its log is `logs/<step>.log` or `logs/<step>_{sample}.log`. |
| `Error in rule <rule>:` while the workflow runs | The block names the job's `log:`, the file that holds the tool's own messages. See [Failures while the workflow runs](#failures-while-the-workflow-runs). |
| `Error: Directory cannot be locked.` or `The files below seem to be incomplete.` | A previous run was killed. `Directory cannot be locked` has a second cause: another run from the same launch folder would write the same files, or one of the two runs would write a file the other reads. Use `--unlock` only when no other Snakemake process runs from that folder. See [Resuming after an interruption](reference/running.md#resuming-after-an-interruption). |
| A finished run with an empty or missing result | Some messages go only to the logs and never stop the run. See [Messages that do not stop the run](#messages-that-do-not-stop-the-run). |
| A message from Slurm, a job killed for memory or time, `SSL certificate problem`, or a step that works on one computer but not on the cluster | See [On a cluster](#on-a-cluster). |

Three kinds of log are kept:

| Log of | Where |
|---|---|
| every per-sample rule | `<output_dir>/logs/<rule>_{sample}.log` |
| the rules that run once for all samples | `<output_dir>/logs/`, without a sample name: `funannotate_db.log`, `normalize_funannotate_db.log`, `unite_db.log`, `unite_identified.log`, `antismash_db.log` or `antismash_db_local.log`, `tiberius_model.log`, `telomere_db.log`, `busco_data.log` (a named BUSCO lineage), `download_phix.log`, `build_phix.log`, `install_signalp6.log`, `chroquetas_db.log`, `multiqc.log` |
| Snakemake itself | `.snakemake/log/` in the launch folder, one file per run, named after its start time |

`--show-failed-logs` makes Snakemake print the log of every failed job at the end of the run.

## Stops before any job runs

FunFlux checks the whole configuration before any job starts, on every run and every dry run. A
failed check prints one `[FunFlux]` line and exits with status 1. Only the first failed check is
reported, so fix it and run `-n` again. Lines that start with `[FunFlux] NOTE:` do not stop the
run.

Placeholders in the messages below are in angle brackets. An empty cell means that the message
itself already says what is wrong or what to change.

### Configuration file and mode

These checks concern the configuration file itself and its `mode` key.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] No configuration supplied. Pass one explicitly, e.g.:` | No `--configfile` on the command line. FunFlux has no default configuration. | Add `--configfile config/config_custom.yaml`. |
| `[FunFlux] mode: this configuration has no 'mode' key. ...` | A v1 configuration (FunFlux 1.x, FunFluxL, FunFluxL+), or an empty `mode:`. | Start from the v2 template. [Coming from v1](getting-started/from-v1.md) lists what moved. |
| `[FunFlux] mode: must be one of illumina | longread | hybrid | contigs (got: '<value>')` | Another word in `mode`. | Use one of the four. |

### Keys and values

These messages name a key that is missing, unknown or in the wrong block, or whose value has the
wrong form or is out of range.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] <block>: must be a block of 'key: value' lines.` | `input`, `directories`, `files`, `links`, `resources`, `parameters` or a block inside `parameters` holds a single value or a list. | Indent the keys under the block name. |
| `[FunFlux] <section>.<key>: required, but it is missing or empty in the config.` | A required entry is empty: `directories.output_dir`, `eggnog_db`, `funannotate_db`, `blast_db` (when the screen runs), `files.iprscan`, `files.annotation_params`, `links.funannotate_link`, `links.phix_link` (in `illumina` and `hybrid` mode). The template leaves them empty on purpose. | Fill it in. |
| `[FunFlux] <key>: expected true or false (got <value>).` | A yes/no key holds another word. `parameters.contigs.screen` is checked in every mode, although only `contigs` mode uses it. | Write `true` or `false`. |
| `[FunFlux] <key>: expected a whole number (got <value>).`, `... must be at least N (got M).`, `... must be at most N (got M).` | The value is not an integer, or it is out of range. | The allowed values are listed on [Configuration](reference/configuration.md). |
| `[FunFlux] <key>: expected a number (got <value>).`, `... must be above 0 (got X).`, `... must be at most N (got X).` | The same for decimal keys, for example `min_direct_blast_qcov: 150`. `nan` and `inf`, and YAML's `.nan` and `.inf`, count as not a number. | As above. |
| `[FunFlux] parameters.<block>: unsupported key(s): <keys>. Supported keys: <list>.` | A misspelled key inside a block of `parameters`. | Use a key from the list. |
| `[FunFlux] parameters.<key>: unknown setting (a typo?). Known settings: ...` | A misspelled key directly under `parameters`. | Use a name from the list. |
| `[FunFlux] parameters.<key>: this setting moved into the block named after your mode: ...` | `long_read_assembler`, `flye_input_mode`, `medaka_model`, `filtlong` or `hifiasm` at the top of `parameters`, as in FunFluxL. | Move it into `parameters.longread` or `parameters.hybrid`. `long_read_assembler` is now called `assembler`. |

!!! warning "A misspelled optional key outside `parameters` is ignored without a message"

    Only `parameters` and its blocks are checked for unknown keys. A misspelled optional key in
    `directories`, `files`, `links`, `input` or `resources` is ignored, and the default applies.
    For example, `antismash:` instead of `antismash_db:` in `directories` means that the antiSMASH
    database is downloaded, and `resources.thread:` means 24 threads.

### Input files and sample names

These checks look for the input files of the active mode and take the sample names from the file
names.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] mode=<mode> requires input.<key> to be set in the config.` | The key for the input folder of the active mode is missing or empty. | Set `input.contigs_dir`, `input.illumina_dir` or `input.longread_dir` or, in `hybrid` mode, the last two. |
| `[FunFlux] input.<key>: '<value>' is not an existing directory.` | Wrong path. A relative path is read from the launch folder, not from `output_dir`. | Use an absolute path. |
| `[FunFlux] input.<key>: no input files found. Expected names like <example>.` | No file follows the naming pattern. | Name the files as the mode expects: `{sample}.fasta` (or `.fa`, `.fna`) for assemblies, `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz` for Illumina pairs, and `{sample}_ont.fastq.gz` or `{sample}_pacbio.fastq.gz` for long reads. |
| `[FunFlux] input.<key>: more than one file extension in this folder (<a>, <b>). ...` | Mixed extensions, or other files in the folder. In `contigs` mode every file whose name contains a dot is taken as a sample, so a `README.txt` gives `(fasta, txt)`. | Use one extension for all files of one kind. In `hybrid` mode the Illumina files and the long-read files may have different ones. Keep other files elsewhere. |
| `[FunFlux] input.<key>: file extension '<ext>' is not recognised (expected one of: ...).` | An extension FunFlux does not recognise, for example `.fas`. | Use `fasta`, `fa` or `fna` for assemblies, and `fastq`, `fq`, `fastq.gz` or `fq.gz` for reads. |
| `[FunFlux] sample name '<name>' contains characters that are not allowed (...). Rename the input files.` | The name holds a space or one of `* # @ % ^ / ! ? & : ; | < >`. A FASTA in a subfolder of `contigs_dir` becomes the sample `subfolder/name` and stops here. | Rename the files. Underscores and dots are allowed. Keep the input folder flat. |
| `[FunFlux] input.illumina_dir: sample '<name>' has an R1 file but no R2 mate: <path>` | The R2 file is missing or named differently. | Name the pair `{sample}_R1.<ext>` and `{sample}_R2.<ext>`. |
| `[FunFlux] sample '<name>' has both _ont and _pacbio files in <dir>. ...` | One sample has files of both long-read technologies. | Keep one long-read file per sample. |
| `[FunFlux] hybrid mode needs Illumina and long reads for every sample.` | Some samples have reads of one kind only. The message goes on with `No _ont or _pacbio file for:`, `No Illumina _R1 file for:` or both, each followed by the sample names. | Supply both kinds of reads, or run the sample in another batch. |

### Annotation parameter table

These messages concern the annotation parameter table, its rows and its columns.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] files.annotation_params: '<path>' not found.` | Wrong path in `files.annotation_params`. | |
| `[FunFlux] files.annotation_params: line N has K tab-separated columns; expected 3 or 4 (sample, species, proteins and, optionally, augustus_species).` | Spaces instead of tabs, or a missing or extra column. A species name with a space is still one column. | Separate the columns with tabs. |
| `[FunFlux] files.annotation_params: sample '<name>' names the AUGUSTUS model '<model>', which is not a trained species of the funannotate database at '<db>'. Did you mean '<model>'? Leave column 4 empty, or write '-', to train AUGUSTUS on this genome instead. ...` | Column 4 names a model that the database does not hold, for example through a typo, or with capitals where the folder name has none. The message appears at start-up when the database is already on disk. Otherwise it comes from rule `augustus_model_check`, right after the database is unpacked, and is also written to `logs/augustus_model_check_{sample}.log`. | Use one of the suggested names, which are the folder names in `<db>/trained_species/`, or leave column 4 empty to train AUGUSTUS on the genome. |
| `[FunFlux] files.annotation_params: sample '<name>' names the AUGUSTUS model '<model>', which is not a possible model name (letters, digits, '_', '-' and '.' only). ...` | Column 4 holds a path, a space or another character that no model name has. | Write the model's folder name only, or leave column 4 empty. |
| `[FunFlux] files.annotation_params: sample '<name>' appears more than once (line N).` | A duplicate row. | Keep one row per sample. |
| `[FunFlux] files.annotation_params: sample '<name>' has no row in '<path>'.` | Every sample of the batch needs a row. Rows for samples outside the batch are ignored. | Add the row. |
| `[FunFlux] files.annotation_params: protein file '<path>' of sample '<name>' does not exist.` (or `... is empty.`) | A relative path in the table is read from the launch folder. | Give an absolute path to a non-empty protein FASTA. |

A sample whose column 4 is empty can also stop later, during gene prediction, with
`Not enough gene models <N> to train Augustus`; see
[Other stops in the shared stages](#other-stops-in-the-shared-stages).

### Databases and links

These messages concern the reference databases and the download links of the configuration.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] directories.antismash_db: '<path>' is not a directory. ...` (or `... has no 'clusterblast'. ...`) | The key must name the antiSMASH databases folder itself, the one that holds `clusterblast/`, `pfam/` and the rest. | Point it there, or leave it empty to download the database. |
| `[FunFlux] directories.unite_db: '<path>' is not a directory. ...` (or `... has no 'unite_its_sintax.fasta'. ...`) | The folder must hold the UNITE SINTAX reference as an uncompressed `unite_its_sintax.fasta`. | Unpack or rename the file, or leave the key empty to download it. |
| `[FunFlux] links.unite_its_link: required, unless directories.unite_db names a local UNITE copy.` | Both are empty. | Set one of them. |
| `[FunFlux] links.funannotate_link: expected a Zenodo URL of the form https://zenodo.org/records/<id>/files/<archive>, got: <url>` | FunFlux finds the checksum file in the same Zenodo record, so the link must be a Zenodo record URL. | Use a Zenodo record URL, as the template does. |
| `[FunFlux] links.funannotate_link: cannot read a file name from the URL: <url>` | The URL ends without a file name. | As above. |
| `[FunFlux] links.telomere_db_link: required, but it is missing or empty in the config.` | The link to tidk's telomere repeat database is required in every mode. A configuration written before the contig summary existed lacks it. | Copy the line from `config/config.yaml`. |
| `[FunFlux] parameters.prediction.busco_db: '<name>' is not in the funannotate database at '<db>'. Add it with FunFlux's funannotate conda environment activated ('Other BUSCO sets for gene prediction' in the documentation shows how to find it): funannotate setup -i busco -b <name> -d <db> -w` | The funannotate database ships only `dikarya`. | Run the command shown, in the activated environment (see [Other BUSCO sets for gene prediction](getting-started/databases.md#busco-sets)), or use `dikarya`. |
| `[FunFlux] parameters.prediction.busco_db: expected a plain dataset name such as dikarya (got: '<value>')` | The name holds characters other than letters, digits and underscores. | |
| `MissingInputException in rule blob_json ...` (from Snakemake) | The screen runs, and `nodes.dmp` or `names.dmp` is missing from `directories.blast_db`. When the folder itself is missing, `[FunFlux] NOTE: directories.blast_db '<path>' does not exist on this machine; ...` comes first. | Correct the path, or unpack the NCBI taxonomy dump into that folder (see [Reference databases](getting-started/databases.md)). |

### Screen, masking and BUSCO settings

These messages concern the settings of the screen, of repeat masking, of BUSCO and of the contig
summary.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] parameters.decontamination.mode: must be one of off | include | exclude (got: <value>)` | Another word. | Use one of the three. |
| `[FunFlux] parameters.decontamination.mode: 'auto' is not available in FunFlux. BLAST genus labels are unreliable for poorly sequenced fungi, so auto can remove most of the genome. Use off (the default), read contig_taxonomy_decisions.tsv, then exclude the contaminant genera or include the genera that BLAST assigned to the isolate's own contigs.` | `mode: auto`, which FunFlux 1.x, FunFluxL and FunFluxL+ had, is removed, for the reasons given on [Configuration](reference/configuration.md#decontamination). | Use `off`. To filter, use `exclude` or `include` as the message says. |
| `[FunFlux] parameters.mitochondria.mode: must be one of off | report | separate_confirmed (got: <value>)` | Another word, or `true` or `yes`, which YAML reads as true. An unquoted `off` is accepted. | Use one of the three. |
| `[FunFlux] parameters.masking_method: must be tantan | repeatmodeler_repeatmasker (got: <value>)` | Another word. | Use one of the two. |
| `[FunFlux] parameters.busco.lineage: must be 'auto' or a BUSCO lineage name such as hypocreales_odb12 (got: <value>). ...` | A lineage name must end in `_odb` and a number. `busco --list-datasets` prints each name with its version after a dot, such as `hypocreales_odb12.2`, and the dot and the number after it are not part of the name. | Write `auto`, or a name such as `hypocreales_odb12`. |
| `[FunFlux] parameters.contig_summary.telomere_motif: must be auto (the unit is chosen from the data) or the telomere repeat unit of the species, 4 to 30 letters A, C, G and T, such as TTAGGG (got: <value>)` | A value that is neither `auto` nor a plain unit. It may hold a character other than A, C, G and T, for example an IUPAC code such as N, a bracket or a digit. It may also be shorter than 4 or longer than 30 letters. | Write `auto`, or the unit in plain letters. The known units are listed on [Contig summary](analysis/contig-summary.md#which-fungi-carry-which-repeat). |
| `[FunFlux] parameters.contig_summary.telomere_db_sha256: must be the SHA-256 checksum of the telomere database file, 64 hexadecimal characters (got: <value>)` | A shortened or mistyped checksum. | Copy the 64 characters from `config/config.yaml` or, for another version of the database, the ones `sha256sum` prints for that file (see [Telomere repeat database](getting-started/databases.md#telomere-db)). |
| `[FunFlux] parameters.decontamination.sample_overrides: sample '<sample>' has mode 'auto' in '<path>'. Mode auto is not available in FunFlux. BLAST genus labels are unreliable ...` (the rest as in the `auto` message above) | `auto` in the `mode` column of the override table. | Write `off`, `include`, `exclude` or nothing in that cell. |
| `[FunFlux] parameters.decontamination.sample_overrides: sample '<sample>' has mode '<mode>' in '<path>'; the mode column must be off, include, exclude or blank (blank keeps the run-wide mode).` | Another word in the `mode` column. Every row is checked, including the rows of samples outside the batch. | |
| `[FunFlux] parameters.decontamination.sample_overrides: sample '<sample>' has discard_no_hit '<value>' in '<path>'; the discard_no_hit column must be true or false (also yes/no, on/off, 1/0) or blank (blank keeps the run-wide setting).` | Another word in the `discard_no_hit` column of a sample of the batch. | |
| `[FunFlux] parameters.decontamination.<key>: '<path>' is not an existing file.` | `include_genera_by_sample`, `exclude_genera_file` or `sample_overrides` names a file that does not exist. A relative path is read from the launch folder. | Correct the path, or use an absolute one. |
| `[FunFlux] parameters.decontamination.<key>: '<path>' is not UTF-8 text. Save it again as plain text with UTF-8 encoding.` | The file has another encoding, for example a table that a spreadsheet program saved as UTF-16. | |
| `[FunFlux] parameters.decontamination.<key>: '<path>' must start with a tab-separated header line that has the columns sample and genus (missing: <columns>).` (`sample and mode` for `sample_overrides`) | The first line lacks a column name, or the columns are separated by spaces. | Start the table with the column names, separated by tabs. |
| `[FunFlux] parameters.decontamination: mode include keeps only the genera you name, but no genus is given for sample(s) <samples>. Name the genus in include_genera, in include_genera_by_sample (one row per sample) or in the include_genera column of sample_overrides.` | `include` is set, run-wide or in `sample_overrides`, for samples without a genus. | |

The two `decontamination.mode` messages can stop any run. The checks of the decontamination files
and of `include` without a genus run only when the screen runs, after the checks of the input
folder and of the annotation table.

### GeneMark-ES settings

These messages concern GeneMark-ES and its licence key.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] parameters.genemark.mode: must be auto | off (got: <value>)` | Another word. An empty value (`mode:` or `mode: ""`) means `auto`. | Write `auto` or `off`. |
| `[FunFlux] directories.genemark_dir: must be an absolute path when GeneMark is enabled.` | A relative path. | Use an absolute path. |
| `[FunFlux] directories.genemark_dir: '<dir>' has no executable gmes_petap.pl.` | Wrong folder, or the script lacks execute permission. | Point at the folder that holds `gmes_petap.pl`, or make it executable (`chmod u+x`). |
| `[FunFlux] GeneMark is enabled but there is no licence key at <home>/.gm_key. ...` | GeneMark-ES is on (`directories.genemark_dir` set, `genemark.mode: auto`), but the user running Snakemake has no `~/.gm_key`. | Install the key (see [Licensed and external tools](getting-started/external-tools.md)), or set `parameters.genemark.mode: off`. |

### Tiberius settings

These messages concern Tiberius, which is on by default, and the download of its model.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] parameters.tiberius.mode: must be on | off (got: <value>)` | Another word. An empty value means `on`. | Write `on` or `off`. |
| `[FunFlux] parameters.tiberius.model_sha256: expected the SHA-256 checksum of the model archive, 64 characters 0-9 and a-f (got: <value>)` | A shortened or mistyped checksum. | Copy the 64 characters `sha256sum` prints, or leave the key out for the default. |
| `[FunFlux] links.tiberius_model_link: required, but it is missing or empty in the config.` | Tiberius is on (the default) and the link is not in the configuration, for example in a configuration written before Tiberius was added. | Copy the line from `config/config.yaml`, or set `parameters.tiberius.mode: off`. |

### SignalP 6 settings

These messages concern SignalP 6 and the DTU archive it is installed from.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] parameters.signalp6.mode: must be off | fast | slow | slow-sequential (got: <value>)` | Another word. | Use one of the four. |
| `[FunFlux] files.signalp6_archive: required when parameters.signalp6.mode is not off.` | SignalP 6 is on without an archive. | Name the DTU archive, or set `mode: off`. |
| `[FunFlux] files.signalp6_archive: must be an absolute path to the official DTU archive.` | A relative path. | Use an absolute path. |
| `[FunFlux] files.signalp6_archive: '<path>' is not a readable file.` | Wrong path or permissions. | |

### AMR screen settings { #amr-screen-settings }

These messages concern the AMR screen (the resistance-mutation screen), which runs ChroQueTas
against the FungAMR database.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] parameters.amr.chroquetas: must be on | off (got: <value>)` | Another word. An empty value means `off`. | Write `on` or `off`. |
| `[FunFlux] parameters.amr.mito_genetic_code: must be auto or the number of an NCBI genetic code: 1-6, 9-16 or 21-33 (got: <value>)` | Another word, or a number that is not an NCBI genetic code (7, 8 and 17-20 are none, and miniprot refuses them). An unquoted `yes`, `no`, `on` or `off` fails too, because YAML reads it as true or false. | Write `auto`, or the number of the code (see [Two runs](analysis/amr.md#two-runs)). |
| `[FunFlux] parameters.amr.chroquetas_script_sha256: must be the SHA-256 checksum of ChroQueTas.sh, 64 hexadecimal characters (got: <value>)` (or `chroquetas_db_sha256: must be the SHA-256 checksum of the FungAMR database archive, ...`) | A shortened or mistyped checksum. | Copy the 64 characters from `config/config.yaml`, or leave the key out for the default. |
| `[FunFlux] links.chroquetas_link: required, but it is missing or empty in the config.` | The AMR screen is on and the link is not in the configuration. | Copy the line from `config/config.yaml`, or set `parameters.amr.chroquetas: off`. |

### Read-based modes

These messages concern the read and assembly settings of the read-based modes. `<mode>` is
`longread` or `hybrid`, and each of these blocks is checked only in its own mode. `spades_kmers`
is checked in `illumina` and `hybrid` mode.

| Message | Cause | Fix |
|---|---|---|
| `[FunFlux] parameters.spades_kmers: expected auto, or whole numbers separated by commas (got: <value>).` | A word other than `auto`, or a list with another separator. | Write `auto`, or for example `"21,33,55,77,99,127"`. |
| `[FunFlux] parameters.spades_kmers: every k must be odd and below 128, and the list ascending (got: <value>).` | An even k, a k below 1 or above 127, a k given twice, or a list out of order. | |
| `[FunFlux] parameters.<mode>.assembler: must be hifiasm | flye (got: <value>)` | Another word. | |
| `[FunFlux] parameters.<mode>.flye_input_mode: 'pacbio-hifi' is not a value any more. ...` | PacBio HiFi reads are recognised from the file name. | Name the reads `{sample}_pacbio.fastq.gz` and set `flye_input_mode: auto`. |
| `[FunFlux] parameters.<mode>.flye_input_mode: must be auto | nano-raw | nano-hq (got: <value>)` | Another word. | |
| `[FunFlux] parameters.<mode>.filtlong.target_bases: must be empty or a size such as 1500000000 or 1.5g (got: <value>)` | | |
| `[FunFlux] parameters.<mode>.hifiasm.genome_size: must be empty or a size such as 40000000 or 40m (got: <value>)` | | |
| `[FunFlux] parameters.<mode>.hifiasm.telomere_motif: must be empty or a DNA/IUPAC motif (got: <value>)` | | |

### Notes that announce a later failure

These NOTE lines do not stop the run. Each one warns of a problem that shows later, in a step or in
the result.

| Note | What follows |
|---|---|
| `[FunFlux] NOTE: directories.eggnog_db '<path>' does not exist on this machine; ...` | `eggnog_annotation` fails. |
| `[FunFlux] NOTE: files.iprscan '<path>' does not exist on this machine; ...` | `iprscan_annotation` fails. |
| `[FunFlux] NOTE: GeneMark-ES is enabled and needs a valid, unexpired licence key at <home>/.gm_key. ...` | An expired key stops the run at `genemark_check`. |
| `[FunFlux] NOTE: sample name '<name>' is not a valid NCBI locus-tag prefix ...` | The run finishes, but the gene IDs cannot be submitted to GenBank as they are. A valid prefix is 3 to 12 letters or digits, starting with a letter. |

!!! warning "`directories.eggnog_db` and `files.iprscan` are only tested for existence"

    An eggNOG folder without the data, or an InterProScan installation that does not start,
    passes these checks and fails when its rule runs.

## Failures while the workflow runs

These messages come from single steps after the jobs have started. Each entry gives the message
and its cause, and most also say what to do. The step's log holds the tool's own messages.

### GeneMark-ES licence

```text
[FunFlux] GeneMark-ES refused to run because of its licence (key file <home>/.gm_key). The key has probably expired: request a new one from the GeneMark site, or set parameters.genemark.mode: off. See <log>.
```

**Cause.** Rule `genemark_check` ran GeneMark-ES on a slice of the masked genome, and the output
mentioned the licence. The key in `~/.gm_key` has most likely expired (see
[The licence check before gene prediction](getting-started/external-tools.md#the-licence-check-before-gene-prediction)).

**Fix.** Download a new key from the GeneMark site and copy it to `~/.gm_key` of the user who runs
Snakemake. GeneMark's INSTALL file gives the two commands: `gunzip gm_key.gz`, then
`cp gm_key ~/.gm_key`. Then run the same Snakemake command again.

To run without GeneMark-ES, set `parameters.genemark.mode: off`. Gene prediction, and every step
after it, then runs again for the samples already predicted with GeneMark-ES.

```text
[FunFlux] GeneMark-ES could not start from <dir>: the Perl named in the first line of gmes_petap.pl (GeneMark ships it as /usr/bin/perl, the computer's own) is missing or lacks a Perl module GeneMark needs. funannotate starts GeneMark the same way, so its gene prediction would fail too. Fix it once, in that folder: perl change_path_in_perl_scripts.pl '/usr/bin/env perl' (GeneMark's scripts then use the Perl of FunFlux's funannotate environment, which has the modules), or install the missing modules into that Perl, or set parameters.genemark.mode: off. See <log>.
```

**Cause.** The test output contained `Can't locate`, which means a missing Perl module, or
`bad interpreter`, which means that the Perl named in the script's first line does not exist. The
check starts `gmes_petap.pl` the way funannotate predict does, so the first line of the script
chooses the Perl. GeneMark ships it as `#!/usr/bin/perl`, which often lacks GeneMark's modules. On
a cluster, the compute nodes can have another `/usr/bin/perl` than the login node.

**Fix.** Run `perl change_path_in_perl_scripts.pl "/usr/bin/env perl"` once in the GeneMark
folder. The scripts then use the Perl of the funannotate environment, which holds the modules.

```text
[FunFlux] GeneMark-ES could not start from <dir> (a program or a file is missing). See <log>.
```

**Cause.** The test output contained `command not found` or `No such file`.
`logs/genemark_check_{sample}.log` holds GeneMark's output and names what is missing.

Any other failure on the test slice, for example too little sequence to train on, is logged and
the run continues.

### SignalP 6

```text
[FunFlux] signalp6_annotation: SignalP 6 is not installed in the conda environment of this rule (the environment was probably rebuilt after SignalP 6 was installed). Delete the file <output_dir>/00.resources/signalp6/installation.json and run the workflow again, so that install_signalp6 installs SignalP 6 into the current environment.
```

**Cause.** SignalP 6 is installed into the signalp6 conda environment, and the record of that
installation, `installation.json`, lives in the output folder. A run with another
`--conda-prefix`, from another launch folder, or after `.snakemake/conda/` was deleted gets a new
environment without SignalP 6. The old record still keeps `install_signalp6` from running again.

**Fix.** Delete `00.resources/signalp6/installation.json` and run again. SignalP 6 and funannotate
annotate then run again for every sample of the run.

When the installation itself fails, `install_signalp6` stops with
`SignalP 6 installation failed: ...` in `logs/install_signalp6.log`. The rest of the message names
the cause:

| The message continues with | Cause and fix |
|---|---|
| `the supplied archive does not contain the '<mode>' model (...). Download the DTU archive for the selected mode.` | The archive holds the model of another mode. A fast archive such as `signalp-6.0i.fast.tar.gz` works only with `mode: fast`. Get the archive for the mode, or change the mode. |
| `configured archive_sha256 does not match the supplied archive (expected ..., observed ...)` | `files.signalp6_archive` names another file than the one the checksum was taken from. Check the path, or update `parameters.signalp6.archive_sha256` with what `sha256sum <archive>` prints. |
| `archive_sha256 must be a 64-character hexadecimal SHA-256 digest` | A typo in `parameters.signalp6.archive_sha256`. |
| `expected exactly one SignalP Python package ...`, `archive member has an unsafe path: ...`, `archive contains an unsupported link or special file: ...` | The file is not the archive as DTU distributes it. Use the downloaded file unchanged. |
| `the official archive did not install the signalp6 command` | The package installed without its command-line program. The pip output is above the message in the same log. |

### InterProScan

A failed InterProScan job is started once more before the run stops, and
`logs/iprscan_annotation_{sample}.log` keeps the last attempt only.

InterProScan is a Java program, and the memory it may use, its heap, is the `-Xmx` value on the
`java` line at the end of `interproscan.sh` in your InterProScan installation. FunFlux does not
change it, and `resources.ram_gb` does not reach it. If `logs/iprscan_annotation_{sample}.log`
shows `java.lang.OutOfMemoryError`, raise `-Xmx` there. Keep it below the memory the machine can
spare while eggNOG-mapper, antiSMASH and SignalP 6 run beside InterProScan.

After InterProScan ends, the rule checks that `iprscan/interproscan_out.xml` is complete and is an
InterProScan result, and stops otherwise:

```text
InterProScan XML validation failed for <path>: <reason>
```

**Cause.** InterProScan exited without an error but left a file that does not parse, or that is
not an InterProScan result. Its own messages are above this line in the same log. A full disk is
one cause to rule out, as the XML can be a few GB.

**Fix.** Run again. The rule deletes the old XML before it starts, so a stale file is never taken
for a new result.

`Could not find a *.proteins.fa file in <dir>/predict_results or update_results.` means that
funannotate predict wrote no protein file. Its messages are in
`logs/funannotate_prediction_{sample}.log`.

### Tiberius

Rule `tiberius_model` downloads the Tiberius model, and `tiberius_prediction` runs Tiberius on
each sample.

| Message | Log | Cause and fix |
|---|---|---|
| `[FunFlux] tiberius_model: the file downloaded from <link> does not match parameters.tiberius.model_sha256, so it was NOT unpacked ...` | `tiberius_model.log` | The download was cut short, or the authors replaced the model at the same address. Run again for a fresh download. If the checksum still differs, the model has changed, so check it and put its checksum into `parameters.tiberius.model_sha256`, or set `parameters.tiberius.mode: off` (see [Tiberius model](getting-started/databases.md#tiberius-model)). |
| `curl: (6) Could not resolve host` or another `curl` error | `tiberius_model.log` | No network access to `bioinf.uni-greifswald.de` (see [Network access](getting-started/installation.md#network)). |
| the job ends with `Killed`, or Snakemake reports that `tiberius_prediction` was terminated without a Python error | `tiberius_prediction_{sample}.log` | The job ran out of memory, as Tiberius needs about 16 GB. Free memory, run fewer jobs at once, or set `parameters.tiberius.mode: off`. |

### Telomere repeat database and tidk

Rule `telomere_db` downloads tidk's telomere repeat database, and `telomere_units` reads it when it
chooses the telomere unit of each sample.

| Message | Log | Cause and fix |
|---|---|---|
| `[FunFlux] telomere_db: the file downloaded from <link> does not match parameters.contig_summary.telomere_db_sha256, so it was NOT used ...` | `telomere_db.log` | The download was cut short, or the link and the checksum do not belong together. The template address names one commit of the database repository, so the file behind it does not change. The download is kept as `00.resources/telomeres/telo_database_v2.tsv.failed`. Run again for a fresh download. To accept a file you have looked at, put its checksum, which `sha256sum <file>` prints, into `parameters.contig_summary.telomere_db_sha256`. For another version of the database, set the link and its checksum together (see [Telomere repeat database](getting-started/databases.md#telomere-db)). |
| `curl: (6) Could not resolve host` or another `curl` error | `telomere_db.log` | No network access to `raw.githubusercontent.com` (see [Network access](getting-started/installation.md#network)). |
| `[FunFlux] telomere_units: error: ...` | `telomere_units_{sample}.log` | The script stopped on a missing input, a malformed database file or a tidk failure, and tidk's own message follows. A database file that does not parse means that the file at the link has changed. Delete `00.resources/telomeres/telo_database_v2.tsv` and run again, and the checksum then says whether it is still the expected file. |

### ChroQueTas (AMR screen)

With the AMR screen on, rule `chroquetas_db` downloads ChroQueTas and the FungAMR database, and
`amr_chroquetas` runs for each sample.

| Message | Log | Cause and fix |
|---|---|---|
| `[FunFlux] chroquetas_db: the file downloaded from <link> does not match parameters.amr.chroquetas_script_sha256, so it was NOT used ...` (or `... parameters.amr.chroquetas_db_sha256, so it was NOT unpacked ...`) | `chroquetas_db.log` | The download was cut short, or the link and the checksums do not belong together. The template address names one commit, so the files behind it do not change. The download is kept as `00.resources/chroquetas/ChroQueTas.sh.failed` or `db-v20250811.tgz.failed`. Run again for a fresh download. For another version of ChroQueTas, set the link and both checksums together (see [ChroQueTas and the FungAMR database](getting-started/databases.md#chroquetas)). |
| `[FunFlux] chroquetas_db: the database archive from <link> holds no species folder ...` | `chroquetas_db.log` | The link names a commit whose `FungAMR_db/db-v20250811.tgz` is not the database ChroQueTas expects. Use the template link. |
| `curl: (6) Could not resolve host` or another `curl` error | `chroquetas_db.log` | No network access to `raw.githubusercontent.com` (see [Network access](getting-started/installation.md#network)). |
| `[FunFlux] amr_chroquetas: error: ...` | `amr_chroquetas_{sample}.log` | The script stopped on a missing input file, a database folder without species folders, or a position list of the database (`<target>.txt`) that does not look like the one of ChroQueTas 1.0.1. It also stops on a path with a space, or when a ChroQueTas run left neither its summary file nor the message that no target could be extracted. The last 20 lines of ChroQueTas's log follow the message. A position list unlike that of ChroQueTas 1.0.1 means another version of the database, so use the template link and checksums. ChroQueTas cannot handle a path with a space, so move the output folder to a path without spaces. ChroQueTas also takes `bash`, `grep -P` and other shell commands from the machine, so check them as shown under [Commands outside the conda environments](getting-started/installation.md#host-commands). |

### funannotate database

Rule `funannotate_db` downloads the funannotate database and unpacks it into
`directories.funannotate_db`, and `normalize_funannotate_db` repairs the paths stored in it.

| Message | Log | Cause and fix |
|---|---|---|
| `[FunFlux] funannotate_db: <archive> does not match its checksum file <sha256>, so it was NOT unpacked ...` | `funannotate_db.log` | The download was cut short or damaged. Delete both files from `directories.funannotate_db` and run again. Both are then downloaded afresh. |
| `curl: (6) Could not resolve host` or another `curl` error | `funannotate_db.log` | No network access to `zenodo.org` (see [Network access](getting-started/installation.md#network)). Download the funannotate database by hand on a machine with network access and unpack it into `directories.funannotate_db`, as described under [Downloading it by hand](getting-started/databases.md#downloading-it-by-hand). Then run again. A folder that holds `funannotate-db-info.txt` is not downloaded again. |
| `ERROR: Funannotate trained-species directory is missing: ...`, `ERROR: Funannotate database index is missing: ...` | `normalize_funannotate_db.log` | `directories.funannotate_db` does not hold a complete, unpacked funannotate database. Point it at the folder that holds `trained_species/` and `funannotate-db-info.txt`. If the path is right but the database is incomplete, delete `funannotate-db-info.txt` and run again, as in the next row. |
| `ERROR: Trained-species asset is missing after extraction: ...`, `ERROR: Database asset is missing after extraction for archived path: ...` | `normalize_funannotate_db.log` | The unpacked database lacks files. Delete `funannotate-db-info.txt` from the database folder and run again. The archive kept there is then checked against its checksum and unpacked again, without a new download. |
| `ERROR: [Errno 13] Permission denied: ...` | `normalize_funannotate_db.log` | The database folder is read-only and its stored paths still need repair. Make it writable for one run. Once the paths are right, this step only reads it. |

!!! warning "Two runs repairing one funannotate database at once can damage it"

    Two runs repairing one database at the same time, or one database reached under two different
    paths, cause damage or repeated work. See
    [Reference databases](getting-started/databases.md#what-the-workflow-does).

### BUSCO downloads

Where BUSCO's lineage data come from depends on `parameters.busco.lineage`. With `auto`, every
BUSCO job needs network access, while a named lineage needs it once, on the computer that runs
Snakemake (see [BUSCO on compute nodes without internet](reference/running.md#cluster-busco)).

| Setting | Download | Without network access |
|---|---|---|
| `auto`, the default | every BUSCO job downloads the datasets for its sample into `02.assembly/{sample}/eval/busco/busco_tmp/`, which is deleted when the job ends, and samples do not share the download. A hypocrealean fungus fetches `eukaryota_odb12` (29.1 MB), six placement files (8.3 MB) and `hypocreales_odb12` (274.6 MB) | the BUSCO job fails and the run stops, but `--keep-going` lets the other jobs finish. BUSCO's messages are in `logs/completeness_{sample}.log` |
| a lineage, such as `hypocreales_odb12` | rule `busco_data` downloads it once per output folder into `00.resources/busco/`, on the computer that runs Snakemake, and every BUSCO job reads it offline | `busco_data` fails, and its messages are in `logs/busco_data.log`. The BUSCO jobs themselves need no network |

With a named lineage, a failed download stops rule `busco_data`:

| Message | Log | Cause and fix |
|---|---|---|
| `[FunFlux] busco_data: BUSCO did not download the lineage <lineage> (details in <log>). Check parameters.busco.lineage against the list that busco --list-datasets prints, or set it to auto.` | `busco_data.log` | A name BUSCO's server does not know, or no network access to `busco-data.ezlab.org`. |

When BUSCO publishes a newer dataset, a run with `auto` uses it, and its completeness figures can
differ from an earlier run of the same genome. A named lineage stays at the release in
`00.resources/busco/` for every later run into the same output folder.

### Medaka model

These messages apply to ONT samples in `longread` and `hybrid` mode. Rule `check_medaka_model`
checks the Medaka model before the assembler starts. When the check fails, it prints its log, with
the table of available models, to the console as well.

```text
ERROR: Medaka could not infer a model from the reads: their FASTQ headers carry no basecaller tag, or one this Medaka version does not recognise. Set parameters.<mode>.medaka_model to one of the consensus models below (pick by your flowcell / device / accuracy), or to FALSE to skip Medaka.
```

```text
ERROR: Configured medaka_model '<model>' is not available in this Medaka version.
```

**Fix.** Use one of these settings:

| To | Set |
|---|---|
| name the model | `parameters.<mode>.medaka_model`, for example `parameters.longread.medaka_model: r1041_e82_400bps_sup_v4.3.0` |
| skip Medaka | `parameters.<mode>.medaka_model: false` |
| fall back to the inferred model when a named one is not available | `parameters.<mode>.medaka_model_fallback_auto: true` |

```text
ERROR: Medaka model '<model>' exists, but its weights are not stored on this computer and could not be downloaded (Medaka downloads them from GitHub). ...
```

**Cause.** The check downloads the weights of a model Medaka does not ship, so that polishing
needs no network, and it runs on the computer that runs Snakemake.

**Fix.** Start Snakemake on a computer with internet, which on a cluster is the login node.
Otherwise name a model whose weights come with Medaka 2.2.2, `r1041_e82_400bps_hac_v6.0.0` or
`r1041_e82_400bps_sup_v5.2.0`, or set `parameters.<mode>.medaka_model: false`.

[The Medaka model](methods/medaka-model.md#what-can-fail) describes the other messages of the
check, what `auto` needs in the read headers and how to look at them.

### Other stops in the shared stages

The steps below run in more than one mode, and each can stop the run. Where a row mentions the
log, it is the log of the rule in the Rule column.

| Message | Rule | Cause and fix |
|---|---|---|
| `[FunFlux] filter_contigs: no contig of <sample> is left (kept 0 of M records). ...` | `filter_contigs` | In `contigs` mode, the FASTA is empty, or its headers are SPAdes-style and no contig reaches 500 bp and 2x coverage. The log says which filter ran. In `illumina` and `hybrid` mode, no SPAdes contig reached 500 bp and 2x coverage, so check `spades/spades.log` and the reads. |
| `[FunFlux] long_read_consensus: Medaka failed for sample <sample> with model '<model>' ...` | `long_read_consensus` | Medaka itself failed, for example for lack of memory or on a damaged read file, after the model had passed `check_medaka_model`. Medaka's messages are in the same log. `parameters.<mode>.medaka_model: false` runs without polishing. |
| the hifiasm rule stops after `hifiasm` finished, with no `[FunFlux]` line | `long_read_assembly` | hifiasm's primary graph holds no contig, so the step that turns it into `assembly.fasta` stopped the rule. hifiasm's own messages are in the log. |
| `No contigs were kept for sample '<sample>' with mode '<mode>'. ...` | `select_contigs` | The decontamination policy removed every contig, for example `include` with a genus that BLAST never assigned. Snakemake deletes the outputs of the failed job, the decision table included. `contaminants/bestscore.blob.blobDB.table.txt` still shows the genus of every contig, and `--keep-incomplete` keeps the failed job's files. See [Decontamination](analysis/decontamination.md). |
| `Mode 'include' requires at least one genus.` | `select_contigs` | `decontamination.mode: include` is set without a genus for the sample. The dry run stops on this first, with `[FunFlux] parameters.decontamination: mode include keeps only the genera you name, ...` (above), so the rule meets it only when a table was changed after the run started. |
| `ERROR: Mitochondrial split would leave an empty nuclear assembly.` | `separate_mitochondria` | Every contig was confirmed mitochondrial, for example in a FASTA that holds only a mitogenome. Run it without the screen, or with `mitochondria.mode: report`. |
| `[FunFlux] its_windows: ITSx's fungal profile file .../ITSx_db/HMMs/F.hmm is missing or empty. Is the job running in its itsx conda environment (--sdm conda), with the pinned ITSx 1.1.3?` | `its_windows` | The job did not run in the itsx environment, or that environment has no intact ITSx 1.1.3. Run with `--sdm conda`. If the run already does, delete that environment's folder, which `--list-conda-envs` shows, so that Snakemake builds it again from the pin file. |
| `[FunFlux] its_windows: no 5.8S profile found in .../ITSx_db/HMMs/F.hmm (no NAME line contains 58S). Does this ITSx version name its 5.8S profiles differently?` | `its_windows` | The ITSx in the itsx environment is not the pinned 1.1.3. Delete that environment's folder, which `--list-conda-envs` shows, so that Snakemake builds it again from the pin file. |
| `RepeatModeler did not produce a usable repeat library.` | `repeat_masking` | Neither `{sample}_rm-families.fa` nor `consensi.fa.classified` exists or holds sequence. RepeatModeler's messages are in the same log. `masking_method: tantan` does not use RepeatModeler. |
| `Not enough gene models <N> to train Augustus (200 required), exiting` (in `logs/funannotate_prediction_{sample}.log`) | `funannotate_prediction` | Column 4 of the sample is empty, so funannotate trains AUGUSTUS, and BUSCO found fewer than 200 usable training genes. The assembly may be fragmented, or the fungus far from the BUSCO set (`fungi` has 290 genes). Name a stored AUGUSTUS model in column 4 (see [The AUGUSTUS model](analysis/gene-prediction.md#augustus-model)), or use a BUSCO set closer to the fungus. |
| `Funannotate did not create annotations.secretome.txt from the supplied SignalP 6 results.` | `funannotate_annotation` | funannotate annotate did not use the SignalP 6 table. Its messages are in the same log. |
| any MultiQC error | `multiqc` | Unlike a failed QUAST report, a failed MultiQC report stops the run. See `logs/multiqc.log`. |

### Messages that do not stop the run

After these messages the run goes on, so a finished run can still lack a result. The Log column
says where each message is written.

| Message | Log | Meaning |
|---|---|---|
| `[FunFlux] genome_assembly_evaluation: QUAST exited with an error for sample <sample>. ...` | `genome_assembly_evaluation_{sample}.log` | QUAST is a report that no other step reads, so the run continues and the MultiQC page lacks the QUAST table for this sample. After fixing the cause, delete `02.assembly/{sample}/eval/quast/` and run again. |
| `[FunFlux] map_evaluation: Qualimap exited with an error for sample <sample>. ...` (in `hybrid` also `map_guide_evaluation`) | `map_evaluation_{sample}.log`, `map_guide_evaluation_{sample}.log` | Qualimap is a report, so the run continues and the MultiQC page lacks the mapping table for this sample. After fixing the cause, delete `02.assembly/{sample}/eval/qualimap/` (or `guide/eval/qualimap/`) and run again. |
| `[FunFlux] raw_long_read_qc: NanoPlot exited with an error for sample <sample>. ...` (or `filtered_long_read_qc`) | `raw_long_read_qc_{sample}.log`, `filtered_long_read_qc_{sample}.log` | NanoPlot is a report, so the run continues and the MultiQC page lacks the read-QC row for this sample. After fixing the cause, delete `01.reads/{sample}/long/raw_qc/` (or `filtered_qc/`) and run again. |
| `[FunFlux] its_extraction: no 5.8S window was found in the genome of <sample> ...` | `its_extraction_{sample}.log` | No 5.8S gene matched the 358 anchor profiles of ITSx. The ITS files are empty and `{sample}_ITS_copies.txt` holds only a comment line. A fungus whose 5.8S differs from all of them gives this result. |
| `No ITS1 sequence to classify for sample <sample> ...` (or ITS2) | `its_classification_{sample}.log` | ITSx found none, or every record failed the length guard, which accepts ITS1 and ITS2 of 50-600 bp. `{sample}_ITS_copies.txt` lists the rejected records. |
| `[FunFlux] telomere_units <sample>: no telomere repeat identified (best candidate <UNIT>: <n> ends in the telomere orientation, <m> opposite)` | `telomere_units_{sample}.log`, and the console | No unit reached 2 contig ends in the telomere orientation with more than twice as many as in the opposite one. `telomere_unit.txt` reads `none`, the telomere columns of the contig summary read `0`, and chromosomes come out as `nuclear_contig`. This is usual in a fragmented assembly or in a yeast with irregular telomeres. Set `parameters.contig_summary.telomere_motif` when the repeat is known (see [Telomeres](analysis/contig-summary.md#telomeres)). |
| `[FunFlux] NOTE telomere_units <sample>: ...` | `telomere_units_{sample}.log`, and the console | The chosen telomere unit is unknown to the literature and to the tidk database, or the unit set in `parameters.contig_summary.telomere_motif` is not the unit the data favour. `eval/telomere_candidates.tsv` shows every unit with its score. |
| `[FunFlux] NOTE amr_chroquetas <sample>: ...` | `amr_chroquetas_{sample}.log`, and the console | The species in column 2 of the annotation table matches no folder of the FungAMR database, so the sample is not screened. The NOTE names the folders of the same genus, if any. Correct the name when it is misspelt (see [The species](analysis/amr.md#species)). `00.resources/chroquetas/species.txt` lists the folders. Column 2 is also funannotate's species. A corrected name re-runs `augustus_model_check`, `funannotate_prediction` and everything after them, and it goes into the GenBank and GFF files and their names. Change it only when the new name is the right one for the annotation as well. |
| `WARNING: N genera each hold at least 5% of this assembly (...)` | `select_contigs_{sample}.log` | A mixed culture, or one genome split across related genera by BLAST. Read `contig_taxonomy_decisions.tsv` before trusting the assembly. |
| `WARNING: decontamination is discarding X% of the assembly (...)` | `select_contigs_{sample}.log` | More than 20 % of the bases were removed. If the kept assembly looks incomplete but not contaminated, the contig selector most likely removed genome. |
| `GeneMark-ES stopped on the test slice (exit status N) but not because of its licence or installation; ...` | `genemark_check_{sample}.log` | Gene prediction itself then decides whether GeneMark-ES suits the genome. |

## On a cluster { #on-a-cluster }

Most stops on a cluster come from a value larger than a node or a partition allows, or from a
compute node that cannot see a file or the internet. The step's log in `logs/` and its Slurm log
say which. How to run FunFlux on a cluster is described under
[On a cluster](reference/running.md#on-a-cluster).

| Message | Cause | Fix |
|---|---|---|
| `maximum number of parallel jobs/used nodes has to be specified for remote execution (use --jobs N with N being a number >= 1)` | `--executor slurm` without `--jobs`. | Add `--jobs N`, the number of jobs Snakemake may keep queued at once. |
| `SLURM job submission failed. The error message was sbatch: error: Batch job submission failed: Requested node configuration is not available` | A job asks for more CPUs (`resources.threads`) or more memory (`resources.ram_gb`) than any node of the partition has. | Lower the value to what `sinfo -o "%n %c"` shows for CPUs, or `sinfo -o "%n %m"` for memory in MB (see [How the settings reach Slurm](reference/running.md#cluster-resources)). |
| `sbatch: error: Batch job submission failed: Requested time limit is invalid (missing or exceeds some limit)`, or jobs that stay pending with the reason `PartitionTimeLimit` in `squeue` | `resources.runtime_hours` is longer than the time limit of the partition. | Lower it to what `sinfo -o "%P %l"` shows, or choose a partition with a longer limit through `slurm_partition` in `--default-resources`. |
| `SLURM-job '<id>' failed, SLURM status is: 'OUT_OF_MEMORY'` | The step used more memory than it asked for. | The `Error in rule <rule>:` block names the step. For a heavy step, raise `resources.ram_gb`. The heavy steps are listed under [How the settings reach Slurm](reference/running.md#cluster-resources). For a light step, raise `mem_mb` in `--default-resources`. |
| `SLURM-job '<id>' failed, SLURM status is: 'TIMEOUT'` | The step ran past its time limit. | For a heavy or long step, raise `resources.runtime_hours`. For a light step, raise `runtime` in `--default-resources`. |
| `No wall time information given. ...` and `No job memory information ('mem_mb' or 'mem_mb_per_cpu') is given - submitting without. ...` (warnings) | The light steps ask for no time and no memory, so Slurm's defaults apply. | Nothing, if those defaults are a few GB and an hour or more. Otherwise add `--default-resources mem_mb=8000 runtime=240`. |
| `MissingOutputException ... Missing files after 5 seconds. This might be due to filesystem latency.` (or after 60 seconds) | The shared storage shows a file on the login node a while after a compute node wrote it. | Add `--latency-wait 60`, as in the [Run A command](reference/running.md#cluster-run-a). If it is already there, raise it, for example to 120. |
| a dry run plans finished steps, or `download_all -n` plans downloads again, with the reason `Software environment definition has changed since last execution` | The run uses another `--conda-prefix` than the one the steps ran with. | Use the same `--conda-prefix` as before, from the same launch folder. |
| BUSCO (`completeness`) fails on the compute nodes and `logs/completeness_{sample}.log` shows a failed download | `parameters.busco.lineage: auto` needs internet in every BUSCO job. | Name a lineage, and run `download_all` on the login node (see [BUSCO on compute nodes without internet](reference/running.md#cluster-busco)). |
| `[FunFlux] GeneMark is enabled but there is no licence key at <home>/.gm_key` in a job's log | The compute node does not see the key in the home folder. | Copy `~/.gm_key` to the home folder the compute nodes see (see [The GeneMark-ES key](reference/running.md#cluster-genemark)). |
| `SSL certificate problem` (curl) or `CERTIFICATE_VERIFY_FAILED` (BUSCO, antiSMASH) in a download log | A network proxy re-signs HTTPS with an in-house certificate, and the conda environments trust only their own list. | Export `CURL_CA_BUNDLE` and `SSL_CERT_FILE` before starting Snakemake (see [Behind a proxy that re-signs HTTPS](getting-started/installation.md#proxy)). |
