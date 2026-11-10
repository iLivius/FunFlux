# Helper scripts: how to find the rule a script belongs to

Every script that a rule runs lives in a folder named after **the rule file that
runs it**. The folder name and the `.smk` file name are the same, so you can find
the rule file by eye, without searching:

| script folder | run by (path inside `workflow/`) | what its script does |
|---|---|---|
| `05_ont_adapters/` | `rules/shared/05_ont_adapters.smk` (longread and hybrid modes, ONT samples, unless `ont_adapter` is off) | `find_ont_adapter.py` reads every ONT read, finds the sequencing adapter and the barcode block that follows it, and writes the sequence fastplong removes (or that there is nothing to remove); its `check` subcommand counts the adapters left in the trimmed reads (`rules/shared/00_common.smk` also imports it, to check a sequence given in the configuration) |
| `10_decontam/` | `rules/shared/10_decontam.smk` | `select_contigs_by_taxonomy.py` reads the BlobTools genus table, decides which contigs to keep, and writes a table with the reason for every contig, kept or dropped |
| `12_medaka_check/` | `rules/shared/12_medaka_check.smk` (longread and hybrid modes, ONT samples with Medaka) | `medaka_model_check.py` checks the Medaka polishing model, or works it out from the reads, before the assembly starts, and makes Medaka store the model's weights there and then, so that the polishing needs no internet |
| `15_mitochondria/` | `rules/shared/15_mitochondria.smk`; also `rules/longread/20_assembly.smk` and `rules/hybrid/40_long_assembly.smk` (the `evidence-*` commands, in rule `assembly_evidence`) and `rules/hybrid/25_guide.smk` (`screen` on the guide assembly) | `mitochondria.py` finds the mitochondrial contig(s) and splits them off from the nuclear genome; in the long-read modes it first turns the assembler's own per-contig table (Flye's `assembly_info.txt`, hifiasm's primary graph) into the evidence table the screen reads |
| `25_contig_summary/` | `rules/shared/25_contig_summary.smk` | `telomere_units.py` chooses the telomere repeat unit of the assembly from the data (tidk lists the tandem repeats at the contig ends and calls the arrays of every candidate unit; the unit whose arrays point the way telomeres point wins, checked against the literature and tidk's telomere database) and writes the per-end telomere calls; `adapter_screen.py` turns the BLAST hits of adapter and barcode sequences on the draft contigs into adapter loci at contig ends and inside contigs, and writes the sample's row of the Adapters table of the report; `contig_summary.py` gathers the evidence about every contig of the draft assembly (those telomere calls, read depth, rDNA units, BUSCO markers, genes, the screens' verdicts, the adapter loci) into one table per sample and gives each contig a verdict with its reasons |
| `30_taxonomy/` | `rules/shared/30_taxonomy.smk` | `unite_identified.py` keeps the UNITE records whose taxonomy names a genus, so that a record identified only to the kingdom cannot hold the SINTAX result there; `its_windows.py` cuts one window around every rDNA copy, so ITSx sees one unit at a time; `its_distinct.py` collapses the identical ITS copies ITSx returns, counts them, and keeps records of implausible length away from the classifier |
| `40_annotation/` | `rules/shared/40_annotation.smk` | `normalize_funannotate_db_paths.py` rewrites the absolute paths stored inside a downloaded funannotate database, so that they point into the database folder on this computer; `check_augustus_model.py` works out, as funannotate does, whether AUGUSTUS uses a stored model or is trained on the genome, stops on a model name the database does not hold, and records the choice (`rules/shared/00_common.smk` also imports it, to make the same check at parse time) |
| `45_signalp6/` | `rules/shared/45_signalp6.smk` | `install_signalp6.py` installs SignalP 6 from the licensed archive the user supplies |
| `50_amr/` | `rules/shared/50_amr.smk` (only when `parameters.amr.chroquetas` is on) | `chroquetas_amr.py` runs ChroQueTas, which rule chroquetas_db downloads, to screen the genome for the antifungal resistance mutations of the FungAMR database: once on the nuclear genome, and once for cytochrome b on the mitochondrial contigs with the mitochondrial genetic code; it matches the species of the annotation table to a database folder, turns ChroQueTas's files into one row per target protein and one row per position checked, and gives every change a verdict from FungAMR's evidence for and against resistance |

**A script loose in this folder, not in a numbered one, is not part of a run.**
It is a maintenance tool you run by hand (`clean_funflux_output.sh`).

## Tests sit next to the script they test

`15_mitochondria/mitochondria.py` is tested by `15_mitochondria/test_mitochondria.py`,
and so on (`install_signalp6.py` has no tests). Because a test file and its script
are in the same folder, pytest finds the script on its own. Nine tests
(`test_find_ont_adapter.py`, `test_mitochondria.py`, `test_telomere_units.py`, `test_adapter_screen.py`, `test_contig_summary.py`, `test_its_windows.py`, `test_its_distinct.py`, `test_unite_identified.py`, `test_chroquetas_amr.py`) also put
their own folder first on Python's search path, so the import works whichever
folder pytest is started from.

Run the whole suite from the repository root. pytest is not part of Snakemake:
install it into the Snakemake environment first (`docs/requirements.txt` says how).

```bash
pytest workflow/scripts miscellaneous -q
```

Most tests of `50_amr/test_chroquetas_amr.py` run the script with a stand-in for
ChroQueTas that hands out copies of real ChroQueTas output, so they need nothing
else. Six tests run the real ChroQueTas, five of them on
`50_amr/test_data/b0510_beta_tubulin_region.fasta` (3.8 kb of the public *Botrytis
cinerea* B05.10 reference genome around the beta-tubulin gene), two of those with a
codon changed at test time. They are skipped unless the programs of the chroquetas
environment (`envs/chroquetas.yaml`) are on PATH and the environment variable
`FUNFLUX_CHROQUETAS_TEST_DIR` names a folder that holds `ChroQueTas.sh` and
`FungAMR_db/`, as rule chroquetas_db lays them out.

## Paths are declared once

No rule spells out a script path. Each one is a constant in
`rules/shared/00_common.smk` (`ONT_ADAPTER_SCRIPT`, `SELECT_TAXONOMY_SCRIPT`, `MITOCHONDRIA_SCRIPT`,
`MEDAKA_CHECK_SCRIPT`, `TELOMERE_UNITS_SCRIPT`, `ADAPTER_SCREEN_SCRIPT`, `CONTIG_SUMMARY_SCRIPT`, `UNITE_IDENTIFIED_SCRIPT`, `ITS_WINDOWS_SCRIPT`, `ITS_DISTINCT_SCRIPT`, `NORMALIZE_DB_SCRIPT`,
`AUGUSTUS_CHECK_SCRIPT`, `SIGNALP6_INSTALL_SCRIPT`, `CHROQUETAS_AMR_SCRIPT`), so moving a script means editing one line, not hunting
through the rule files. If you move something, run

```bash
python miscellaneous/check_comment_references.py
```

It fails when a comment names a repository path, a rule (written as "rule" followed
by the rule's name) or a dotted config key that cannot be found. It does not see
everything: it reads the `#` comment lines (not docstrings) of the Snakefiles, the
rule files, the Python and shell scripts and the two config templates, and the
Markdown pages; in a git checkout, only the files git tracks. It checks only rule
names that contain an underscore, and only paths that start with `workflow/`,
`config/`, `docs/` or `miscellaneous/`, or with `rules/`, `scripts/` or `envs/`,
which it looks up under `workflow/` (so the paths in the table at the top of this
page are checked too). In a git checkout, a path that git ignores, such as your own
`config/config_custom.yaml`, counts as present, and any other file must be tracked
by git. Lines that describe earlier versions or planned files are skipped (the top of
the script has the details). Names in capitals that no file assigns at the start of a
line (`NAME = ...`; a shell `export NAME=...` is not recognised), and names that look
like a misspelt rule, never make it fail: `--advisories` lists them.
`--config` names another config file to check the keys against; the default is
`config/config.yaml`.
