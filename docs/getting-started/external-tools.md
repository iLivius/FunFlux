# Licensed and external tools

Three tools are not installed by conda: InterProScan, which every run needs, and the optional
GeneMark-ES and SignalP 6. You obtain them yourself, and FunFlux calls them from its own
environments. Install InterProScan before the first run. For each tool, this page explains how to
obtain and install it, how FunFlux runs it and what happens without it.

| Tool | Needed | Licence | Config keys | Without it |
|---|---|---|---|---|
| [GeneMark-ES](#genemark-es) | optional | free for academic, non-profit and US government users; the key expires one year after download | `directories.genemark_dir`, `parameters.genemark.mode` | funannotate predict runs without GeneMark-ES |
| [SignalP 6](#signalp-6) | optional, off by default | DTU: download page for academic users, others contact DTU | `files.signalp6_archive`, `parameters.signalp6.*` | no secretome annotation |
| [InterProScan](#interproscan) 5.78-109.0 | required | Apache-2.0 | `files.iprscan` | the run stops: InterProScan cannot be switched off |

The licences of GeneMark-ES and SignalP 6 limit who may use them freely, and so do two of the
databases InterProScan searches in every run: commercial users of PROSITE pay a yearly licence fee,
and SMART is free for academic users only. The licence terms in full are on
[Licensing](../about/licensing.md).

The optional AMR screen, also called the resistance-mutation screen, runs ChroQueTas. ChroQueTas
is not installed by conda either, and its licence allows non-commercial use only. You do not
download it yourself. With `parameters.amr.chroquetas: on`, the workflow downloads it at run time,
as described in [ChroQueTas and the FungAMR database](databases.md#chroquetas).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **self-training gene finder** | a gene finder that learns its parameters from the genome it annotates, with no stored model; GeneMark-ES works this way |
    | **licence key** | the file `~/.gm_key` that GeneMark-ES needs; it expires one year after download, and an expired key stops the run at `genemark_check`, before gene prediction |
    | **secretome** | the predicted secreted proteins, which funannotate annotate lists from the SignalP 6 result in `annotations.secretome.txt`; none without SignalP 6 |
    | **model weights** | the trained numbers that make up a model, such as a SignalP 6 model |
    | **Sec/SPI** | signal peptides of the Sec pathway, cleaved by signal peptidase I; with `--organism eukarya`, SignalP 6 predicts only these |
    | **pin file** | the file `<name>.linux-64.pin.txt` beside each environment file: the exact package builds FunFlux v2 was tested with, installed without solving dependencies |
    | **`--conda-prefix`** | the Snakemake option naming the folder the conda environments are built in, instead of `.snakemake/conda/` in the launch folder; several clones or users can share it |
    | **lookup service** | the EBI online service that returns InterProScan matches already computed for known proteins; it receives checksums, not sequences, and saves run time |

## GeneMark-ES { #genemark-es }

GeneMark-ES is a self-training gene finder that funannotate predict adds to AUGUSTUS, SNAP and
GlimmerHMM. It is optional, and FunFlux uses it only when you name an installation.

### Licence and installation

GeneMark-ES is free ("royalty free") for academic and non-profit institutions and US government
agencies. Other users must contact the GeneMark team for a licence. Download the program and the
licence key from <https://genemark.bme.gatech.edu/license_download.cgi>. In October 2026, the
download for 64-bit Linux was GeneMark-ES/ET/EP+ 4.73.

The licence is valid for one year from the date of download, and the date of the key file does not
tell when it expires. A key file dated January 2026 worked in July and had expired in October.

Unpack the program, for example into `/path/to/gmes_linux_64_4`, and follow its `INSTALL` file:

```bash
gunzip gm_key.gz
cp gm_key ~/.gm_key
cd /path/to/gmes_linux_64_4
perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
```

The last command makes GeneMark's Perl scripts start the first `perl` on the `PATH` instead of
`/usr/bin/perl`. This matters because funannotate predict starts `gmes_petap.pl` as a program, so
the first line of the script chooses the Perl. Once the command has run, FunFlux runs GeneMark-ES
with the Perl 5.32.1 of its `funannotate` environment. That environment contains the Perl modules
the `INSTALL` file lists, so nothing has to be installed from CPAN. Without the command,
GeneMark-ES runs with the computer's `/usr/bin/perl`, which often lacks those modules.

The key must be `~/.gm_key` in the home folder of the account that runs Snakemake and its jobs. On
a cluster, it must be in the home folder that the compute nodes see, as described in
[The GeneMark-ES key](../reference/running.md#cluster-genemark). Do not commit or share the key
file.

```yaml
directories:
  genemark_dir: /path/to/gmes_linux_64_4   # absolute path, holds gmes_petap.pl
parameters:
  genemark:
    mode: auto     # auto: use GeneMark-ES when genemark_dir is set; off: never
```

With GeneMark-ES enabled, the run stops before any job starts when `genemark_dir` is not an
absolute path, when `gmes_petap.pl` is missing from it or not executable, or when `~/.gm_key` does
not exist. Whether the key has expired is checked later, before gene prediction starts, as the
next section describes.

### The licence check before gene prediction

GeneMark-ES itself checks only that the key file exists. An expired key is found inside GeneMark's
programs, partway through funannotate predict, which then fails with
`Your license period has ended`. Rule `genemark_check` makes this failure show within a minute
instead. It runs once per sample, after repeat masking. It starts GeneMark-ES the way funannotate
does, with the options funannotate uses, on the first 200,000 bytes of the masked genome
(`04.annotation/{sample}/contigs_mask.fasta`). funannotate predict waits for its result.

What FunFlux does next depends on what GeneMark-ES prints on this slice:

| GeneMark-ES output on the slice | FunFlux |
|---|---|
| mentions "license" | stops: `[FunFlux] GeneMark-ES refused to run because of its licence (key file <home>/.gm_key). The key has probably expired: request a new one from the GeneMark site, or set parameters.genemark.mode: off. See <log>.` |
| contains "can't locate" or "bad interpreter" | stops: `[FunFlux] GeneMark-ES could not start from <genemark_dir>: the Perl named in the first line of gmes_petap.pl (GeneMark ships it as /usr/bin/perl, the computer's own) is missing or lacks a Perl module GeneMark needs. ...`; the message names `change_path_in_perl_scripts.pl` as the fix |
| contains "command not found" or "no such file" | stops: `[FunFlux] GeneMark-ES could not start from <genemark_dir> (a program or a file is missing). See <log>.` |
| another error | continues, because the slice may be too short to train on; `genemark_check.txt` records the exit status |
| no error | continues; `genemark_check.txt` reads `GeneMark-ES ran on the test slice.` |

The check therefore stops the run only for an expired key or for a GeneMark-ES that cannot start.
`genemark_check.txt` is written to `04.annotation/{sample}/`, and the GeneMark-ES output to
`logs/genemark_check_{sample}.log`. After a stop, put the new key in `~/.gm_key` and start the same
command again. The check then runs again, because a failed check leaves no result.

### Without GeneMark-ES

With `genemark_dir` empty or `mode: off`, there is no `genemark_check` job, and funannotate predict
runs without GeneMark-ES. Rule `funannotate_prediction` clears `GENEMARK_PATH` and limits `PATH` to
the funannotate environment plus `/usr/bin` and `/bin`, so a GeneMark-ES installed elsewhere is not
picked up. FunFlux v2 was tested with GeneMark-ES on, and with `mode: off` on one public genome, as
reported in [Validation on public genomes](../about/validation.md#gene-finders).

## SignalP 6 { #signalp-6 }

SignalP 6 finds signal peptides in the predicted proteins, and funannotate annotate turns its
result into the secretome annotation,
`04.annotation/{sample}/funannotate/annotate_misc/annotations.secretome.txt`. SignalP 6 is
optional and off by default.

### Licence and download

The SignalP 6 page of DTU, <https://services.healthtech.dtu.dk/services/SignalP-6.0/>, states:
"There is a download page for academic users; other users are requested to contact DTU Health
Technology Software Package Manager". DTU offers version 6.0i as two packages, `fast` and
`slow_sequential`. FunFlux never downloads SignalP and contains no DTU address. It installs the
archive you name.

`parameters.signalp6.mode` chooses the model that SignalP 6 runs, and so the model that the
archive must hold:

| `parameters.signalp6.mode` | Model the archive must hold | Notes; speed and memory from the SignalP README |
|---|---|---|
| `off` (default) | none | |
| `fast` | `distilled_model_signalp6.pt` | a smaller model that approximates the full one; the DTU package `fast` holds only this model; the mode FunFlux v2 was tested with |
| `slow` | `ensemble_model_signalp6.pt` | the full model, run in parallel; needs more than 14 GB of RAM |
| `slow-sequential` | `sequential_models_signalp6` | the full model, run sequentially: the memory of `fast`, 6 times slower |

FunFlux looks for the model in `models/` and in `signalp/model_weights/` of the archive.

```yaml
files:
  signalp6_archive: /path/to/signalp-6.0i.fast.tar.gz   # absolute path
parameters:
  signalp6:
    mode: fast
    archive_sha256:          # optional: the 64-character checksum, the first field sha256sum <archive> prints
    torch_num_threads: 8     # threads of each signalp6_annotation job, at most resources.threads
    batch_size: 10           # proteins per batch (--bsize)
```

With a mode other than `off`, the run stops before any job starts when `signalp6_archive` is empty,
not an absolute path, or not a readable file.

### How FunFlux installs it

FunFlux installs SignalP 6 from your archive once per output folder, and then runs it once per
sample. The two steps appear in the log under these names:

1. `install_signalp6` runs once per output folder, before the first SignalP job. It checks the
   archive against `archive_sha256` when that is set, and refuses an archive with absolute paths,
   `..` or links in it. It then checks that the model of the chosen mode is there. It installs
   SignalP 6 into the `signalp6` conda environment (Python 3.10, PyTorch 1.13.1 for CPU), with
   nothing fetched from PyPI. It keeps the unpacked package, model weights included, in
   `00.resources/signalp6/package/` (1.6 GB with the fast model), and records what it installed
   in `00.resources/signalp6/installation.json`.
2. `signalp6_annotation` runs once per sample. It runs SignalP 6 on the predicted proteins with
   `--organism eukarya`, the chosen `--mode`, `torch_num_threads` and `batch_size` (`--bsize`),
   and writes the result to `04.annotation/{sample}/signalp6/prediction_results.txt`.
   `--organism eukarya` limits the prediction to the Sec/SPI signal peptides of eukaryotes.

The installation log is `logs/install_signalp6.log`, and its errors start with
`SignalP 6 installation failed:`. The package inside the 6.0i `fast` archive reports version
`6.0h`, and that is what `installation.json` records.

### When the conda environment is rebuilt

SignalP 6 is installed inside the `signalp6` conda environment, and a new environment does not
have it. Snakemake builds a new one when a FunFlux update changes `signalp6.yaml` or its pin file,
when `--conda-prefix` points elsewhere, or when `.snakemake/conda` is deleted. The next SignalP job
then stops with a message that names `00.resources/signalp6/installation.json`. Delete that file
and start the same command again, with the environment folder writable, and SignalP 6 is installed
again. See [Troubleshooting](../troubleshooting.md#signalp-6).

!!! warning "The SignalP 6 archive and its model weights must not be passed on"

    Both are licensed. Do not commit, upload or pass on the archive or the model weights.
    `00.resources/signalp6/package/` holds the weights, so remove it before you share an output
    folder. pip also keeps a built copy of the SignalP package in the pip cache (`~/.cache/pip`)
    of the account that ran the installation. `pip cache remove signalp6` deletes that copy. With
    the `fast` archive, whose weights sit in `models/`, that copy holds no weights. With an archive
    whose weights sit in `signalp/model_weights/`, pip puts the weights into the pip cache and
    into the `signalp6` conda environment as well, and so into a shared `--conda-prefix`.

### Without SignalP 6

With `mode: off`, no SignalP job runs and funannotate annotate gets no SignalP result. The log of
funannotate annotate then reads `Skipping secretome: neither SignalP nor Phobius searches were run`.
Rule `funannotate_annotation` limits `PATH` to the funannotate environment plus `/usr/bin` and
`/bin`, so a SignalP installed elsewhere is not picked up.

## InterProScan { #interproscan }

InterProScan assigns protein families, domains and GO terms to every predicted protein, and
funannotate annotate builds the domain and GO annotation from its XML. Every run needs it.
InterProScan 5 is not on bioconda, so FunFlux runs the `interproscan.sh` of your installation
inside its `interproscan` environment. That environment supplies Java (OpenJDK 11.0.30) and the
Python 3 (3.11.15) that InterProScan runs for seven of its analyses: Gene3D, FunFam, HAMAP,
MobiDB-lite, PANTHER, PIRSR and the PROSITE profiles.

!!! note "FunFlux v2 runs InterProScan 5, not InterProScan 6"

    FunFlux v2 was tested with InterProScan 5.78-109.0, the final release of InterProScan 5.
    Later InterPro releases need InterProScan 6, a different program that FunFlux v2 does not use.

### Installation

The InterProScan documentation gives these steps for 5.78-109.0, a 7.0 GB download that takes
50 GB once installed:

```bash
mkdir -p /path/to/software && cd /path/to/software
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.78-109.0-64-bit.tar.gz.md5
tar -pxvzf interproscan-5.78-109.0-64-bit.tar.gz
cd interproscan-5.78-109.0
python3 setup.py -f interproscan.properties
```

`setup.py` indexes the HMM databases, and `-f` re-indexes all of them. The InterProScan
documentation also requires Perl 5 and Python 3.8 or newer on the machine. In a FunFlux run, both
come from the `interproscan` environment (Python 3.11, Perl 5.32 with thread support), together
with Java 11.

To test the installation by hand, use the environment FunFlux runs InterProScan in. It is built
from the pinned package list in the FunFlux repository and holds the same Java, Python and Perl,
so the test runs as it would in FunFlux on any Linux system. On a system without Python 3.8 or
newer, run `setup.py` in this environment too:

```bash
conda create -n interproscan_test --file /path/to/FunFlux/workflow/envs/interproscan.linux-64.pin.txt
conda activate interproscan_test
```

The InterProScan documentation gives two test commands, run in the installation folder:

```bash
./interproscan.sh -i test_all_appl.fasta -f tsv -dp   # every analysis computed locally
./interproscan.sh -i test_all_appl.fasta -f tsv       # also contacts EBI's lookup service
```

Then point `files.iprscan` at the `interproscan.sh` of the installation:

```yaml
files:
  iprscan: /path/to/software/interproscan-5.78-109.0/interproscan.sh
```

### What FunFlux runs

Rule `iprscan_annotation` runs InterProScan once per sample, on the protein FASTA of funannotate
predict:

```text
interproscan.sh -i <proteins> -f XML -goterms -pa -cpu <threads>
                -T 04.annotation/{sample}/iprscan/tmp
                -o 04.annotation/{sample}/iprscan/interproscan_out.xml
```

- **Threads and run time.** InterProScan is capped at 10 threads (`-cpu`), because it gains little
  from more. Snakemake reserves the 10 threads for the whole job, and with `--cores 24`
  eggNOG-mapper (12 threads) runs beside it. InterProScan takes 1 to 5 hours for a fungal
  proteome, depending on how many of its proteins EBI's lookup service already knows. The XML can
  be a few GB.
- **The XML check.** FunFlux calls InterProScan directly, not through `funannotate iprscan`. Before
  funannotate annotate reads the XML, FunFlux checks that it parses and has the root element
  `protein-matches`. If it does not, the rule stops with
  `InterProScan XML validation failed for ...`.
- **The lookup service.** InterProScan asks EBI's lookup service for matches already computed for
  known proteins. The address of the service is set by `precalculated.match.lookup.service.url`
  in `interproscan.properties`. Without the service, InterProScan computes every match itself, which
  takes longer. See [Functional annotation](../analysis/functional-annotation.md#interproscan).
- **Licensed analyses.** Phobius, SignalP and TMHMM are switched off in the InterProScan download
  because they contain licensed components. FunFlux does not switch them on, and SignalP 6 is its
  own optional step.

### Missing InterProScan

InterProScan cannot be switched off, because funannotate annotate needs its result. What happens
without it depends on `files.iprscan`:

| `files.iprscan` | What happens |
|---|---|
| empty | the run stops before any job starts with `[FunFlux] files.iprscan: required, but it is missing or empty in the config.` |
| a path that does not exist | a NOTE before any job starts, `[FunFlux] NOTE: files.iprscan '<path>' does not exist on this machine; the InterProScan step will fail without it.`; then `iprscan_annotation` fails when it starts, after gene prediction |

!!! warning "A `files.iprscan` path that does not exist only prints a NOTE, and the run fails after gene prediction"

    A NOTE does not stop the run, so a missing InterProScan shows only when `iprscan_annotation`
    starts. Read the NOTEs of the [dry run](../reference/running.md#dry-run).
