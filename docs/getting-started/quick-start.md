# Quick start

This page goes step by step from cloning FunFlux to a finished run, then lists the main results
and what a run costs in time, memory and disk space. The example annotates one assembled genome in
`contigs` mode. For reads, step 4 sets another `mode` and input folder and step 5 needs a real
`directories.blast_db`, while the other steps are the same. If you are not sure which mode fits
your data, see [Choosing a mode](choosing-a-mode.md).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **conda environment** | a separate software installation for one tool, in which its rules run; with `--sdm conda`, Snakemake builds each one from `workflow/envs/` |
    | **rule** | one step of the workflow as Snakemake defines it, such as `filter_contigs`; the console and the log files in `logs/` use its name |
    | **the screen** | the contamination and mitochondrion screen: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **launch folder** | the folder Snakemake is started from, normally the FunFlux repository; relative paths are read from it, and `.snakemake/` is created there |
    | **funannotate database** | the frozen snapshot of funannotate's data, such as Swiss-Prot, Pfam, the BUSCO set and stored AUGUSTUS models, downloaded once from Zenodo into `directories.funannotate_db` |
    | **dry run** | a Snakemake run with `-n`: FunFlux checks the configuration and the inputs, and Snakemake lists the jobs it plans without running any; it takes seconds |
    | **job** | one run of a rule, for one sample or once for the whole run; the dry run lists the jobs it plans |
    | **AMR screen** | the optional resistance-mutation screen, off by default (`parameters.amr.chroquetas`): ChroQueTas looks in each genome for the changes FungAMR lists for the species; results in `05.amr/` |

## 1. Clone FunFlux and install Snakemake

Clone the repository and install Snakemake in an environment of its own:

```bash
git clone https://github.com/iLivius/FunFlux.git
cd FunFlux
conda create -c conda-forge -c bioconda -n snakemake "snakemake>=9.14.6,<10"
conda activate snakemake
```

conda must be 24.7.1 or newer (`conda --version`), because Snakemake 9.14.6 stops with an error on
an older one. On a cluster whose conda module is older, install Miniforge in your own folder and
use its conda. Set `conda config --set channel_priority strict` once, or Snakemake prints a
warning.

FunFlux v2 was developed and tested with Snakemake 9.14.6, and it stops at start-up on an older
version. Later 9.x releases should work, while 10 or later has not been tested. The analysis tools
are not installed in this environment: with `--sdm conda`, Snakemake builds the conda environment
of each tool from its pinned package list in `workflow/envs/` the first time a rule needs it.

Run every command below from the repository folder, where Snakemake finds `workflow/Snakefile`
without being told. See [Installation](installation.md) for the details.

## 2. Provide what FunFlux does not download

Apart from the resources below, FunFlux downloads what it needs itself (step 8). Download or
install these yourself before the first run that needs them, and give their paths in step 5.

| Config key | What | Size | Needed |
|---|---|---|---|
| `directories.eggnog_db` | eggNOG-mapper data of eggNOG 5.0.2 (`eggnog.db`, `eggnog_proteins.dmnd`) | 48 GB | always |
| `files.iprscan` | `interproscan.sh` of a local InterProScan 5.78-109.0 (not on bioconda) | about 50 GB installed | always |
| `directories.blast_db` | NCBI `core_nt`, with `taxdb.btd`, `taxdb.bti`, `nodes.dmp` and `names.dmp` in the same folder | 264 GB (build of 11 July 2026) | `illumina`, `longread`, `hybrid`; `contigs` with the screen |
| `directories.genemark_dir` | GeneMark-ES, with its licence key in `~/.gm_key` | | optional |
| `files.signalp6_archive` | the SignalP 6 archive from DTU | | optional |

`core_nt` is needed only when the contamination and mitochondrion screen runs. The two optional
tools, GeneMark-ES and SignalP 6, are licensed. How to get each resource is explained on
[Reference databases](databases.md) and [Licensed and external tools](external-tools.md).

## 3. Copy the configuration template

Make your own copy of the configuration template, for steps 4 and 5 to fill in:

```bash
cp config/config.yaml config/config_custom.yaml
```

Edit the copy. `config/config_custom.yaml` is listed in `.gitignore`, so it stays out of commits
and `git pull` leaves it alone.

## 4. Set the mode and the input folder

`mode` selects the input mode, and the input key of that mode names the folder of the input files
(`hybrid` has two such keys). The example reads assemblies:

```yaml
mode: contigs

input:
  contigs_dir: /path/to/assemblies      # SAMPLE1.fasta -> sample SAMPLE1
```

The four modes take their input from these keys and file names:

| `mode` | Input key | File names |
|---|---|---|
| `illumina` | `input.illumina_dir` | `{sample}_R1.<ext>` and `{sample}_R2.<ext>` |
| `longread` | `input.longread_dir` | `{sample}_ont.<ext>` or `{sample}_pacbio.<ext>` (PacBio HiFi) |
| `hybrid` | `input.illumina_dir` and `input.longread_dir` (may be the same folder) | both of the above |
| `contigs` | `input.contigs_dir` | `{sample}.fasta`, `{sample}.fa` or `{sample}.fna`, uncompressed |

`<ext>` is `fastq`, `fq`, `fastq.gz` or `fq.gz`, the same for all files of one kind in a batch.
The file name without its suffix is the sample name. It appears in every output path and becomes
the locus-tag prefix of every predicted gene. In `contigs` mode every file in the folder whose name
contains a dot is read as a sample, so keep only the FASTA files there. The naming rules and mixed
batches are described on [Choosing a mode](choosing-a-mode.md).

## 5. Set the output folder and the database paths

These keys tell FunFlux where to write the results, and where to find the databases, InterProScan
and the annotation parameter table of step 6:

```yaml
directories:
  output_dir:     /path/to/results/run1
  eggnog_db:      /path/to/eggnog_db
  funannotate_db: /path/to/funannotate_db
  blast_db:       /path/to/blast_db      # illumina, longread, hybrid; contigs only with the screen

files:
  annotation_params: /path/to/annotation_parameters.tsv
  iprscan:           /path/to/interproscan-5.78-109.0/interproscan.sh

resources:
  threads: 24
```

- `funannotate_db` may name a folder that does not exist yet, and the first run downloads the
  database into it (step 8).
- A relative path is read from the launch folder, the folder you start Snakemake in, not from the
  folder of the config file.
- `resources.threads` is the thread ceiling of the rules, and `--cores` should get the same number.
- The download links under `links:` and everything under `parameters:` can stay as they are.

!!! warning "Let one run repair the paths in `funannotate_db` before another run uses that folder"

    FunFlux repairs the absolute paths stored inside the `funannotate_db` folder without a lock.
    Let one run get past this repair, rule `normalize_funannotate_db`, before you start another
    one on the same folder, as
    [Reference databases](databases.md#what-the-workflow-does) explains.

The optional keys switch on the two licensed tools, point FunFlux at databases you already have,
or add the screen in `contigs` mode.

| To | Set |
|---|---|
| let funannotate use GeneMark-ES | `directories.genemark_dir: /path/to/gmes_linux_64_4` (absolute path; licence key in `~/.gm_key`) |
| predict signal peptides with SignalP 6 | `files.signalp6_archive: /path/to/signalp-6.0i.fast.tar.gz` (absolute path) and `parameters.signalp6.mode: fast` |
| use an antiSMASH database you already have | `directories.antismash_db: /path/to/antismash_db` (the folder that holds `clusterblast/`) |
| use a UNITE file you already have | `directories.unite_db: /path/to/unite` (the folder that holds an uncompressed `unite_its_sintax.fasta`) |
| find the mitochondrial genome and get a BLAST taxonomy per contig | `parameters.contigs.screen: true` and `directories.blast_db` |

Every key is listed on [Configuration](../reference/configuration.md).

## 6. Write the annotation parameter table

For each sample, the annotation parameter table gives the species, a protein FASTA of related
fungi and, optionally, a stored AUGUSTUS model. It has one line per sample, with three or four
columns separated by tabs, and lines starting with `#` are skipped. This command writes the table
of the example:

```bash
printf 'SAMPLE1\tFusarium graminearum\t/path/to/fusarium_proteins.faa\t-\n' \
    > /path/to/annotation_parameters.tsv
```

| Column | Content | In the example |
|---|---|---|
| 1 | the sample name from step 4 | `SAMPLE1` |
| 2 | the species, genus and species; the result files are named after it | `Fusarium graminearum`, giving `Fusarium_graminearum_SAMPLE1.gbk` |
| 3 | a protein FASTA of closely related fungi | `/path/to/fusarium_proteins.faa` |
| 4, optional | a stored AUGUSTUS model: `-` or no fourth column, the default, makes funannotate train AUGUSTUS on the genome's own conserved (BUSCO) genes; a name, one of the folder names in `/path/to/funannotate_db/trained_species/` such as `fusarium_graminearum`, makes AUGUSTUS use that stored model as it is instead | `-` |

Column 3 takes the predicted proteins of a species close to yours, for example the RefSeq
proteome of a species of the same genus. As evidence for its gene models, funannotate aligns them
to the genome, together with UniProtKB/Swiss-Prot, which it always adds. The file must exist and
must not be empty. To download a proteome from NCBI, put the accession of its genome in place of
`<accession>`, for example `GCF_000143535.2` for *Botrytis cinerea* B05.10:

```bash
curl -L -o proteins.zip \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/<accession>/download?include_annotation_type=PROT_FASTA"
unzip proteins.zip      # gives ncbi_dataset/data/<accession>/protein.faa
```

Whether to let AUGUSTUS train or to name a stored AUGUSTUS model is discussed on
[Gene prediction](../analysis/gene-prediction.md#augustus-model). What each column does, and what
the dry run checks, is described under
[Annotation parameter table](../reference/configuration.md#annotation-parameter-table).

## 7. Dry run

Before the real run, check the configuration and the plan with a dry run:

```bash
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml -n
```

`-n` builds the plan and stops. FunFlux checks the configuration and the input files before it
plans any job, so the dry run takes seconds and writes nothing to `output_dir`. It stops at the
first problem with a `[FunFlux]` line that names the key, file or sample to fix. For the example
it prints, after the banner:

```text
Mode: contigs — Annotation of fungal genome assemblies you already have.
[FunFlux] NOTE: contigs mode without parameters.contigs.screen: your FASTA gets no contamination or mitochondrion screen. If its headers are SPAdes-style (>NODE_<n>_length_<L>_cov_<C>), rule filter_contigs drops contigs under 500 bp or under 2x coverage; otherwise it only cuts each header at its first space or tab. A mitochondrial genome in your FASTA, if present, will be annotated as nuclear. Set parameters.contigs.screen: true to find and split it off (needs directories.blast_db).
Decontamination mode: off (discard_no_hit=false).
GeneMark-ES: off (funannotate runs its no-GeneMark branch).
Tiberius: on (its gene models join the others in EVidenceModeler, weight 1).
SignalP 6: off.
Contig summary: telomere motif auto (the repeat unit is chosen from the data), a telomere from 5 tandem copies within 25000 bp of a contig end, an rDNA array from 3 rDNA unit(s).
AMR screen: off.
Sample SAMPLE1 will be processed.
AUGUSTUS: SAMPLE1 is trained on its BUSCO genes (checked once the database is unpacked).
```

The `Decontamination mode` line is printed in every run, and without the screen it has no effect.

Snakemake then prints `Conda environment ... will be created.` for every environment not built
yet (15 here), followed by the jobs it plans:

```text
Job stats:
job                           count
--------------------------  -------
all                               1
antismash_annotation              1
antismash_db                      1
augustus_model_check              1
blast_adapters                    1
cleanup                           1
completeness                      1
contig_summary                    1
eggnog_annotation                 1
filter_contigs                    1
funannotate_annotation            1
funannotate_db                    1
funannotate_prediction            1
funannotate_preprocess            1
genome_assembly_evaluation        1
iprscan_annotation                1
its_classification                1
its_extraction                    1
its_windows                       1
multiqc                           1
normalize_funannotate_db          1
repeat_masking                    1
screen_adapters                   1
telomere_db                       1
telomere_units                    1
tiberius_gff3                     1
tiberius_model                    1
tiberius_prediction               1
unite_db                          1
unite_identified                  1
univec_illumina                   1
total                            31
```

One block per job follows, each with the reason it runs. Each further sample adds 21 jobs, and
other settings change the plan as follows:

| Setting | Jobs |
|---|---|
| funannotate database already in `funannotate_db` | no `funannotate_db` |
| `directories.antismash_db` or `directories.unite_db` set | `antismash_db_local` or `unite_db_local` instead of the download job |
| GeneMark-ES on | `genemark_check` per sample: GeneMark-ES on a 200 kb slice of the masked genome before gene prediction (about a minute), so that an expired licence key stops the run before `funannotate predict` starts instead of partway through it |
| SignalP 6 on | `install_signalp6` once, `signalp6_annotation` per sample |
| `parameters.tiberius.mode: off` | no `tiberius_model`, and no `tiberius_prediction` and `tiberius_gff3` per sample |
| `parameters.its.unite_identified_only: false` | no `unite_identified`: SINTAX classifies against the whole UNITE release |
| `parameters.contigs.screen: true` | 8 per sample: `map_contigs`, `blast_contigs`, `blob_json`, `blob_table`, `select_contigs`, `screen_mitochondria`, `finalize_assembly`, `separate_mitochondria` |
| `parameters.amr.chroquetas: on` | `chroquetas_db` once, `amr_chroquetas` per sample: the AMR screen for antifungal resistance mutations ([Antifungal resistance mutations](../analysis/amr.md)) |

With the funannotate database present, a local antiSMASH database, and GeneMark-ES and SignalP 6
on, the plan has 33 jobs. The read-based modes plan more jobs per sample, for the read processing,
the assembly, the polishing and the mapping QC:

| Mode | Jobs in the first run of one sample |
|---|---|
| `contigs`, the example | 31 |
| `illumina` | 47 |
| `longread` | 45 (PacBio HiFi) or 50 (ONT with Medaka) |
| `hybrid` | 62 (PacBio HiFi) or 67 (ONT with Medaka) |

How these counts are made up is explained on [Running FunFlux](../reference/running.md#dry-run).

!!! warning "A dry run that passes does not show that the eggNOG data, InterProScan and core_nt are usable"

    The dry run checks these paths only in part.

    | Path | What the dry run checks | What happens |
    |---|---|---|
    | `directories.eggnog_db`, `files.iprscan` | existence only | a missing `eggnog_db` folder or `iprscan` file prints a NOTE (`... does not exist on this machine; the eggNOG step will fail without it.`) and the dry run goes on; wrong content there goes unnoticed until eggNOG-mapper or InterProScan fails, after gene prediction, hours into the run |
    | `directories.antismash_db`, `directories.unite_db` | `clusterblast/` and `unite_its_sintax.fasta` | the dry run stops without them |
    | `directories.blast_db`, with the screen on | `nodes.dmp` and `names.dmp` | a folder without either stops the dry run with `MissingInputException in rule blob_json` (`blob_guide_json` in `hybrid` mode); a folder without core_nt passes the dry run and fails when `blast_contigs` starts |

## 8. Run

When the dry run shows the plan you expect, start the run with the same command without `-n`:

```bash
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml
```

`--cores` is the only CPU setting of this command: give it the value of `resources.threads`, and
do not add `--jobs` or `-j`. Each CPU-bound rule declares its threads, and Snakemake runs jobs side
by side while their threads fit into `--cores`. The caps of all rules are listed on
[CPU cores and threads](../reference/running.md#cpu-cores-and-threads).

On a cluster, either let Snakemake send each step to the scheduler (`--executor slurm --jobs N`),
or run this command inside one batch job on one node. Both ways, and the preparation for compute
nodes without internet, are described on [On a cluster](../reference/running.md#on-a-cluster).

Snakemake keeps its records and the conda environments in `.snakemake/` inside the launch folder.
Launch from the same folder every time, or add `--conda-prefix /path/to/conda_envs`, so that the
environments are built once.

After an interruption, run the same command: finished jobs are kept and the rest is redone.
Snakemake asks for `--unlock` after a killed run, and for `--rerun-incomplete` when a job left
incomplete output.

The first run downloads and builds what the table lists. The environments and the funannotate
database are built once and reused by every run that points at the same folder. The files under
`00.resources/` are fetched once per output folder, while with the automatic lineage the BUSCO
data are fetched again in every BUSCO job.

| What | Where | Size | When |
|---|---|---|---|
| conda environments (15 in the example) | `.snakemake/conda/`, or the `--conda-prefix` folder | about 20 GB, plus conda's package cache | first use of each |
| funannotate database snapshot (Zenodo record 23265189) and its `.sha256` file | `directories.funannotate_db` | 2.1 GB archive; 6.1 GB unpacked, 8.1 GB with the archive | when `funannotate-db-info.txt` is missing there |
| antiSMASH 8.0.4 databases | `00.resources/antismash/databases/` | about 9.4 GB | unless `directories.antismash_db` is set |
| UNITE SINTAX reference, release 19.02.2025, all eukaryotes | `00.resources/unite/unite_its_sintax.fasta` | 35.7 MB download, 207 MB unpacked | unless `directories.unite_db` is set |
| Tiberius `fungi` model | `00.resources/tiberius/model/` | 30 MB download, 31 MB unpacked | while `parameters.tiberius.mode` is `on`, the default |
| telomere repeat database of tidk | `00.resources/telomeres/telo_database_v2.tsv` | 0.6 MB | always |
| Illumina records of NCBI UniVec, build 10.1 | `00.resources/adapters/univec_illumina.fasta` | 1.7 MB download, about 31 kB kept (151 records) | always |
| BUSCO lineage datasets | `02.assembly/{sample}/eval/busco/busco_tmp/`, deleted after the job | eukaryota_odb12 29 MB, compressed, plus 8 MB of placement files for the automatic lineage step and the lineage dataset BUSCO places the genome in | every BUSCO job, for every sample, with `parameters.busco.lineage: auto`; a named lineage is fetched once per output folder into `00.resources/busco/` |

Four of these downloads are checked with `sha256sum -c` before they are used:

| Download | Checked against |
|---|---|
| funannotate archive | its `.sha256` file; a mismatch stops the run with the instruction to delete both files and start again |
| Tiberius model | `parameters.tiberius.model_sha256` |
| telomere repeat database | `parameters.contig_summary.telomere_db_sha256` |
| UniVec | `parameters.contig_summary.univec_sha256` |

For the UNITE file, `00.resources/unite/PROVENANCE.txt` records the URL, download date, md5 and
size. For the telomere repeat database, `00.resources/telomeres/PROVENANCE.txt` records the URL,
download date, checksum and row count. For UniVec, `00.resources/adapters/PROVENANCE.txt` records
the URL, download date, checksum and the number of Illumina records kept.

InterProScan computes every match itself and needs no internet: EBI's lookup service, which
returned matches already computed, did not answer InterProScan 5 in November 2026, so FunFlux does
not ask it
(see [Functional annotation](../analysis/functional-annotation.md#interproscan)). The read-based
modes also download the PhiX genome (`illumina`, `hybrid`) and the Medaka model for ONT samples.

## Results

The table lists the main result files. Paths are relative to `output_dir`,
`{sample}` is the sample name, and `<Genus_species>` is the `species` column of the annotation
table with an underscore for the space.

| What | File |
|---|---|
| Genome that QUAST, BUSCO, the ITS steps and funannotate read | `02.assembly/{sample}/contigs_final.fasta` |
| Assembly statistics (QUAST) | `02.assembly/{sample}/eval/quast/report.tsv` |
| Completeness (BUSCO) | `02.assembly/{sample}/eval/busco/short_summary.specific.<lineage>.{sample}.txt` |
| ITS1, ITS2 and full ITS, one record per rDNA unit | `03.taxonomy/{sample}/{sample}.ITS1.fasta`, `{sample}.ITS2.fasta`, `{sample}.full.fasta` |
| ITS copies found, and records rejected by the length check | `03.taxonomy/{sample}/{sample}_ITS_copies.txt` |
| ITS classification against UNITE | `03.taxonomy/{sample}/{sample}_ITS1_taxonomy.txt`, `{sample}_ITS2_taxonomy.txt` |
| Genome as annotated, contigs sorted by length and renamed `scaffold_1`, `scaffold_2`, ... | `04.annotation/{sample}/contigs_sort.fasta` |
| Contigs removed by `funannotate clean` before annotation | `04.annotation/{sample}/contigs_removed_by_clean.txt` |
| Gene models and functional annotation | `04.annotation/{sample}/funannotate/annotate_results/<Genus_species>_{sample}.gbk`, `.gff3`, `.tbl`, `.sqn`, `.proteins.fa`, `.annotations.txt` |
| InterProScan, eggNOG-mapper, antiSMASH | `04.annotation/{sample}/iprscan/interproscan_out.xml`, `eggnog/{sample}.emapper.annotations`, `antismash/index.html` |
| SignalP 6 | `04.annotation/{sample}/signalp6/prediction_results.txt` |
| Antifungal resistance mutations, with the AMR screen on | `05.amr/{sample}/amr_summary.tsv` |
| QUAST and BUSCO of all samples on one page | `06.report/multiqc_report.html` |
| Logs, one per rule and sample | `logs/` |

QUAST, BUSCO and the ITS files use the contig names of `contigs_final.fasta`, while the annotation
uses the `scaffold_N` names. The `annotated_as` column of
`02.assembly/{sample}/eval/contig_summary.tsv` links the two. Every file is described on
[Output files](../reference/output.md).

## Run time and disk

As the table below shows, plan for several hours per genome and for about 16 GB of memory for
Tiberius.

| What | Time, memory or space |
|---|---|
| One fungal genome in `contigs` mode | several hours |
| InterProScan | usually a few hours for a fungal proteome |
| Tiberius | about half an hour for a 40 Mb genome, with about 16 GB of memory |
| Output folder of one sample | in the order of 10 GB |
| A downloaded antiSMASH database | about 9.4 GB more in `00.resources/` |
| funannotate database, outside the output folder | 8.1 GB |
| conda environments, outside the output folder | about 20 GB for the 15 of the example, 22 GB with SignalP 6 |

How long the steps take is described on [Run time](../reference/running.md#run-time), and where
the space goes on [Disk use](../reference/output.md#disk-use).
