# Decontamination

The contamination and mitochondrion screen, called the screen on this page, examines the draft
assembly. This page covers its contamination part, which gives every contig of the draft a genus
and keeps or removes contigs by genus. The screen runs before QC and annotation. A contaminant
contig left in the assembly changes the QUAST, BUSCO and ITS results, and the fungal gene
predictor calls false genes on it.

The genus comes from a megablast search against NCBI `core_nt`. BlobTools turns the hits into one
genus per contig, and the contig selector applies `parameters.decontamination`, which decides what
is kept. The default, `mode: off`, removes nothing and still records a genus and a decision for
every contig, so you can read them before you choose a filter. The
[mitochondrion screen](mitochondria.md) reads the same BLAST results.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector |
    | **draft assembly (the draft)** | the assembly the screen examines, before any contig is removed: `contigs_filt.fasta` in `illumina` mode and in `contigs` mode with the screen, `assembly.fasta` in `longread` and `hybrid` |
    | **read-based modes** | the modes `illumina`, `longread` and `hybrid`, which start from reads, as against `contigs` |
    | **coverage track** | the reads aligned back to the draft, `contaminants/{sample}_map.bam`, from which BlobTools takes each contig's coverage; in `contigs` mode, the contigs aligned to themselves |
    | **contig selector** | rule `select_contigs`: it keeps or removes each contig by genus, as `parameters.decontamination` says; the default, `off`, keeps every contig |
    | **genus filter** | the setting `parameters.decontamination.mode: include` or `exclude`, under which the contig selector keeps or removes contigs by the genus BlobTools gave them |
    | **confirmed mitochondrial contig** | a contig the mitochondrion screen classes `confirmed_mitochondrial`; with the default `separate_confirmed`, it is kept whatever its genus and moved to `contigs_mitochondrial.fasta` |
    | **`--forcerun`** | a Snakemake option that runs the named rule again for every sample, with every step after it, even when Snakemake sees no change |

## When the screen runs

The screen always runs in the read-based modes, and in `contigs` mode only on request. It works the
same in every mode, except for the coverage track ([below](#coverage-track)).

| Mode | Screen | Screened draft, in `02.assembly/{sample}/` | Kept contigs then go through |
|---|---|---|---|
| `contigs` | with `parameters.contigs.screen: true` | `contigs_filt.fasta`: the input FASTA after the contig filter | nothing more |
| `illumina` | always | `contigs_filt.fasta`: the SPAdes contigs after the length and coverage filter | nothing more |
| `longread` | always | `assembly.fasta`: the hifiasm or Flye assembly | Medaka (ONT samples) |
| `hybrid` | always | `assembly.fasta`: the hifiasm or Flye assembly | Medaka (ONT samples), then Polypolish |

- **After the screen**, the kept contigs, polished where the table says so, become
  `contigs_all.fasta`. In every mode that file is then split into the nuclear genome
  `contigs_final.fasta` and `contigs_mitochondrial.fasta`, as described in
  [Mitochondria](mitochondria.md#the-split).
- **In `hybrid`**, the Illumina guide assembly gets a screen of its own, with the same
  `parameters.decontamination` and `parameters.mitochondria` settings. Its steps carry `guide` in
  their names, from `index_guide_contigs` to `select_guide_contigs`, and its files are under
  `02.assembly/{sample}/guide/`. The read pairs that map to its kept contigs, the selected pairs,
  go to Polypolish and, for ONT samples with `short_read_guidance: true` (off by default), to
  Filtlong ([Hybrid](../modes/hybrid.md)). A contig removed from the guide therefore also loses
  its Illumina reads for those steps.
- **With the screen on**, `directories.blast_db` must be set, and an empty value stops the run
  before any job starts. In `contigs` mode without the screen, the dry run prints a NOTE that the
  FASTA gets no contamination or mitochondrion screen
  ([Pre-assembled contigs](../modes/contigs.md)).

## Steps and files

Two steps gather the evidence, `map_contigs` and `blast_contigs`. Two BlobTools steps join it into
one row per contig, and the contig selector, `select_contigs`, decides.

```text
draft contigs
 ├── map_contigs ─────┐                                  coverage track (BAM)
 │                    ├── blob_json ── blob_table ──┐    one taxonomy row per contig
 ├── blast_contigs ───┘                             │
 │        └── screen_mitochondria ── confirmed IDs ─┤    see Mitochondria
 └──────────────────────────────────────────────────┴── select_contigs ── contigs_decontam.fasta
                                                                          + audit files
```

| Step | Tool | Does | Writes, in `02.assembly/{sample}/contaminants/` |
|---|---|---|---|
| `map_contigs` | `contigs`: minimap2 2.30, samtools 1.21; `illumina`: bowtie2 2.5.4 (index built by `index_contigs`), samtools 1.21; `longread`, `hybrid`: minimap2 2.30 with the sample's preset, samtools 1.21 | aligns the reads to the draft: the coverage track | `{sample}_map.bam`, `{sample}_map.bam.bai` |
| `blast_contigs` | BLAST+ 2.16.0 | searches every contig in `core_nt` | `blastout` |
| `blob_json` | BlobTools 1.1.1 | joins the BLAST hits and the coverage | `blob.blobDB.json`, `blob.{sample}_map.bam.cov` |
| `blob_table` | BlobTools 1.1.1 | gives each contig one taxon per rank | `bestscore.blob.blobDB.table.txt` |
| `select_contigs` | `select_contigs_by_taxonomy.py` (Python standard library), with the pinned Python 3.12 of the `python` conda environment | keeps or removes each contig by its genus | `contig_taxonomy_decisions.tsv`, `{sample}_composition.txt`, `contigs.list`, `contigs_decontam.fasta` |

The step names are those Snakemake prints on the console, and each step writes
`logs/<step>_{sample}.log`. BLAST uses at most 24 threads and the mapping `resources.threads`,
while every other step uses one. [Output files](#output-files) lists which files stay after the
run.

## Choosing a policy

Start with `off`, read the decisions, then filter. `off` is the default and removes nothing, while
`exclude` with a genus that core_nt gives to your own contigs removes genome. The table suggests a
setting for each situation.

| Situation | Setting |
|---|---|
| first run of any isolate | `mode: off`; read `contig_taxonomy_decisions.tsv` and the composition file ([Audit files](#audit-files)) |
| a contaminant genus (a bacterium or yeast grown with the culture) has contigs of its own, and core_nt does not give that genus to your fungus | `mode: exclude` with that genus |
| `contig_taxonomy_decisions.tsv` shows which genus labels your isolate's contigs carry, and every other label belongs to a contaminant | `mode: include` with the labels of your isolate's contigs, which need not include its true genus |
| one sample of a batch needs another policy | a `sample_overrides` row for that sample |
| only the contigs without any BLAST hit should go | `mode: exclude` with no genus listed, and `discard_no_hit: true` |

The screen can run alone, so that you see the decisions of a policy before QC and annotation
start:

```bash
snakemake --sdm conda --cores 24 --configfile config/config_custom.yaml --until select_contigs
```

## The keep/drop decision

The contig selector takes the genus of each contig from the BlobTools table
([BlobTools](#blobtools-one-genus-per-contig)) and applies `parameters.decontamination`:

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

The dry run prints the policy in use, for example
`Decontamination mode: off (discard_no_hit=false).`

### Modes

`parameters.decontamination.mode`, written `mode` on this page, chooses one of three policies.

| `mode` | Keeps | Removes |
|---|---|---|
| `off` (default) | every contig | nothing |
| `include` | the contigs of the listed genera | every other contig |
| `exclude` | every contig whose genus is not listed | the contigs of the listed genera |

For each contig the tests run in this order, and the first that applies decides.

1. A confirmed mitochondrial contig is kept
   ([Protected mitochondrial contigs](#protected-mitochondrial-contigs)).
2. With `mode: off`, the contig is kept.
3. With `discard_no_hit: true`, a contig whose genus is `no-hit` is removed. This applies in
   `include` and `exclude` alike.
4. Otherwise the genus test of the mode decides.

`include` must list the genera that BLAST gave the isolate's own contigs, which for a poorly
sequenced fungus need not be its true genus
([Why there is no auto mode](#why-there-is-no-auto-mode)). When the screen runs, a sample left in
`include` mode without any genus stops the run before any job starts:

```text
[FunFlux] parameters.decontamination: mode include keeps only the genera you name, but no genus is given for sample(s) S. Name the genus in include_genera, in include_genera_by_sample (one row per sample) or in the include_genera column of sample_overrides.
```

`discard_no_hit` accepts true/false, yes/no, on/off or 1/0. It changes which contigs are kept only
in `exclude` mode.

| `mode` | Effect of `discard_no_hit: true` |
|---|---|
| `off` | none |
| `include` | only the reason code of the `no-hit` contigs changes (`discarded_no_hit`); they go anyway because they are not on the list |
| `exclude` | the `no-hit` contigs are removed |

### Genus lists and aliases

`include_genera` and `exclude_genera` take one genus, a YAML list, or one string with the names
separated by `;` or `,`:

```yaml
    include_genera: [Botrytis, Alternaria]
    exclude_genera: "Pseudomonas;Bacillus"
```

Genus names are compared without regard to case, and `exclude` uses no aliases. `include` also
accepts two kinds of alias, which only ever keep more contigs.

- **A table of bacterial genus synonyms**, under which *Peribacillus* and *Priestia* match
  *Bacillus*. It contains no fungal genus.
- **The prefix rule**, which strips the prefixes brady, meso, neo, sino, aeri, caldi and geo when
  at least four letters remain. *Neocosmospora* and *Cosmospora* then match each other, and so do
  *Geosmithia* and *Smithia*.

For a fungus, only the prefix rule can apply.

### Per-sample settings

Three settings name files. The two genus files, `include_genera_by_sample` and
`exclude_genera_file`, add to the run-wide lists, and a `sample_overrides` row has the last word.
When the screen runs, the dry run checks the files before
any job starts, as described below, and `select_contigs` reads them again when it runs. A relative
path is resolved from the launch folder, the folder Snakemake is started from, so write absolute
paths.

| Key | File | Combines with the run-wide settings |
|---|---|---|
| `include_genera_by_sample` | TSV with the header columns `sample` and `genus`; several rows per sample allowed; a cell may hold `A;B` | adds the genera to `include_genera` for that sample |
| `exclude_genera_file` | text file, one or more genera per line, separated by `;` or `,` (not by spaces); lines that start with `#` are skipped, and a `#` after a genus is read as part of the name | adds the genera to `exclude_genera` for every sample |
| `sample_overrides` | TSV with the header columns `sample` and `mode` (required) and `include_genera`, `exclude_genera`, `discard_no_hit` (optional) | replaces, for that sample, each setting the row fills in |

The run-wide lists come first, the genera of the two genus files are added to them, and then the
`sample_overrides` row is applied. An override's `include_genera`
therefore replaces both the run-wide list and the `include_genera_by_sample` additions of that
sample.

In the example override table below, which is tab-separated, `isolate1` keeps every contig and
`isolate2` loses its *Pseudomonas* contigs and its contigs without hits. Every other sample follows
the run-wide policy.

```text
sample	mode	include_genera	exclude_genera	discard_no_hit
isolate1	off			
isolate2	exclude		Pseudomonas	true
```

An empty cell keeps the run-wide value. Every row needs all its tabs, including those of empty
cells, because a missing tab shifts the columns.

When the screen runs, a faulty file stops the run before any job starts. The dry run checks the
three files against the samples of the run and stops in these cases.

| Checked | Stops the run when |
|---|---|
| each file that is set | it is not an existing file, or is not UTF-8 text |
| the header line | `include_genera_by_sample` lacks the columns `sample` and `genus`, or `sample_overrides` lacks `sample` and `mode` |
| the mode of every `sample_overrides` row, whatever its sample | it is not `off`, `include`, `exclude` or blank (`auto` gets its own message, with the explanation quoted under [Why there is no auto mode](#why-there-is-no-auto-mode)) |
| the `discard_no_hit` cell of a sample of the run | it is not true/false, yes/no, on/off, 1/0 or blank |
| the genera of a sample of the run | the sample is left in `include` mode without a genus (the message above) |

The messages of the first four checks name the key and the file and, for a bad row, its sample,
for example
`[FunFlux] parameters.decontamination.sample_overrides: sample 'S' has mode 'M' in '<table>'; the mode column must be off, include, exclude or blank (blank keeps the run-wide mode).`
Without the screen, the files are not read and not checked. A table edited after the dry run is
checked again when `select_contigs` runs, and a faulty one then stops that step, with a Python
error in `logs/select_contigs_{sample}.log`.

A configuration that still has the FunFlux 1.x key `parameters.genus` and no `decontamination`
block is translated, with a NOTE, as described in
[Coming from v1](../getting-started/from-v1.md). A genus in the v1 key becomes `mode: include`
with that genus, and the NOTE says that `include` keeps only the contigs BLAST assigned to it. An
empty v1 key becomes `mode: off`, and the NOTE says that FunFlux 1.x then kept the genus with the
most contigs, and that v2 makes no automatic choice:

```text
[FunFlux] NOTE: parameters.genus (a v1 key) is empty. FunFlux 1.x kept the genus with the most contigs when parameters.genus was empty; FunFlux v2 has no automatic choice, so every contig is kept (decontamination mode off). To filter, read contig_taxonomy_decisions.tsv, then set parameters.decontamination: exclude the contaminant genera, or include the genera that BLAST assigned to the isolate's own contigs.
```

### Protected mitochondrial contigs

The contig selector also reads `mitochondria/confirmed_mitochondrial.ids`, the contigs that
`screen_mitochondria` confirmed from the same BLAST table. It keeps them before any other test,
`mode: off` and `discard_no_hit` included, with the reason `protected_confirmed_mitochondrial`.

The protection is needed because BLAST can give a fungal mitochondrion a genus other than the one
a genus filter keeps, or no hit at all. The nuclear contigs of a poorly sequenced fungus, for
example, can be labelled with a better-sequenced relative, while its mitochondrion, whose
reference genome is in core_nt, carries the true genus. Without the protection, `include` with the
relative alone would then remove the mitochondrion.

A protected contig that has no BlobTools row is kept with the genus `no-taxonomy`. A protected ID
that is not in the draft FASTA stops the step
(`Protected mitochondrial contig(s) were not found in FASTA: ...`), which means that the list and
the FASTA come from different assemblies.

!!! warning "Only the default `mitochondria.mode: separate_confirmed` protects the mitochondrion from a genus filter"

    With `mitochondria.mode: report` or `off`, the list of confirmed contigs is empty and nothing
    is protected. A contig that the contig selector removes cannot come back at the later split.

## Why there is no auto mode

FunFlux 1.x and BacFlux, the sister workflow for bacteria, have a fourth mode, `auto`, which keeps
one dominant genus and removes every other contig. FunFlux 1.x and BacFlux up to v2.0.1 chose the
genus that BlobTools gives to the most contigs. BacFlux v2.1.0 by default chooses the genus that
holds the most bases (`auto_vote: bp`). FunFlux v2 has no `auto` mode, because BLAST genus labels
are unreliable for poorly sequenced fungi.

When `core_nt` holds few sequences of the isolate's genus, BlobTools gives most of the isolate's
contigs the genus of better-sequenced relatives, often several of them, and many short contigs can
get no hit at all. In FunFlux 1.x, `auto` counted contigs, not bases, and the `no-hit` contigs took
part in the vote unless `discard_no_hit: true`. The winner could therefore be `no-hit` or a
mislabel, and every contig with another label, genuine genome included, was removed.

`auto` in `parameters.decontamination.mode` stops every run before any job starts:

```text
[FunFlux] parameters.decontamination.mode: 'auto' is not available in FunFlux. BLAST genus labels are unreliable for poorly sequenced fungi, so auto can remove most of the genome. Use off (the default), read contig_taxonomy_decisions.tsv, then exclude the contaminant genera or include the genera that BLAST assigned to the isolate's own contigs.
```

In the mode column of `sample_overrides`, `auto` stops the run with the same explanation, but only
when the screen runs, because without the screen that table is not read.

!!! warning "Removal by genus removes genome when the genus calls are wrong"

    `include` and `exclude` remain, because there you name the genera after reading
    `contig_taxonomy_decisions.tsv` and the composition file. The trap is the true genus. When
    `core_nt` holds few of its sequences, BLAST gives it to only a few of the isolate's contigs,
    and `include` with the true genus alone removes most of the genome. The warnings in
    `logs/select_contigs_{sample}.log` flag such a removal, but with `include` or `exclude` the
    run then goes on to annotate what is left.

    The calls also depend on the BLAST hit depth. Keeping more hits per contig changes the
    bit-score sums, and with them the genus of some contigs
    ([BLAST hit depth](../methods/blast-hit-depth.md)).

## Audit files

The contig selector records its decision on every contig in one table, the share of each genus in
a composition file, and its warnings in a log.

### contig_taxonomy_decisions.tsv

The decisions table has one row per contig of the BlobTools table, in table order, plus any
protected contig without a BlobTools row.

| Column | Content |
|---|---|
| `contig` | contig name, the first word of the FASTA header |
| `assigned_genus` | the BlobTools genus, `no-hit`, or `no-taxonomy` |
| `action` | `keep` or `remove` |
| `reason` | one of the codes below |

The `reason` column holds one of these codes.

| `reason` | Meaning |
|---|---|
| `protected_confirmed_mitochondrial` | confirmed mitochondrial contig, kept in every mode |
| `mode_off` | `mode: off` |
| `discarded_no_hit` | no BLAST hit, and `discard_no_hit: true` |
| `included_genus` | `include`; genus listed |
| `included_genus_alias` | `include`; an alias of a listed genus |
| `not_included_genus` | `include`; genus not listed |
| `excluded_genus` | `exclude`; genus listed |
| `not_excluded_genus` | `exclude`; genus not listed |

With the default settings, every row reads `keep` / `mode_off`, except those of the confirmed
mitochondrial contigs, which read `keep` / `protected_confirmed_mitochondrial`.

### Composition file

`{sample}_composition.txt` has one line per genus, over every contig of the BlobTools table, kept
or removed. Each line gives the genus share of the bases and of the contig count, as fractions with
two decimals, and the lines are sorted by bases:

```text
<genus>: bases <share of the bases>  contigs <share of the contigs>
```

The two shares can disagree, because a genus can hold many short contigs while another holds most
of the genome.

### Log messages

The contig selector writes these messages to `logs/select_contigs_{sample}.log`.

| Message | Printed when |
|---|---|
| `WARNING: N genera each hold at least 5% of this assembly (...). Either the sample is a mixed culture, or one genome is being split across related genera by the BLAST assignment. Check the kept/removed split in contig_taxonomy_decisions.tsv before trusting this assembly.` | two or more genera hold at least 5 % of the bases |
| `WARNING: decontamination is discarding X% of the assembly (R of T bp). If the kept assembly then looks incomplete but NOT contaminated, the filter has most likely removed genome rather than contamination.` | more than 20 % of the bases are removed |
| `WARNING: N contigs from BlobTools table were not found in FASTA.` | the BlobTools table names contigs the FASTA lacks |
| `Sample S decontamination mode: M` and `Kept K of N contig assignments.` | when the step succeeds, as the last two lines |

The warnings change no decision and do not stop the run. They are checked in every mode, so the
first one can appear with the default `mode: off` too.

The step fails in these cases.

| The step fails when | Message |
|---|---|
| no contig is kept | `No contigs were kept for sample 'S' with mode 'M'. Check taxonomy assignments and decontamination settings.` |
| a protected ID is missing from the FASTA | `Protected mitochondrial contig(s) were not found in FASTA: ...` |
| `include` has no genus (the dry run normally stops that case first) | `Mode 'include' requires at least one genus.` |
| the BlobTools table has no readable row | `No contig taxonomy records could be parsed from '<table>'.` |

Snakemake then deletes the outputs of the step, the decisions table included. Add
`--keep-incomplete` to the command line to keep them for inspection. Snakemake treats files kept
this way as finished, so delete them, or add `--forcerun select_contigs`, before the next run.

## Changing the policy after a run

Snakemake notices a changed decontamination setting, BLAST hit depth or mitochondrion threshold,
and re-runs the steps that use it and every step after them.

| Change | Re-runs |
|---|---|
| a new decontamination policy | `select_contigs`, the mitochondrion split, QC, ITS and the whole annotation of every sample: hours per sample |
| a new hit depth, `nt_version` or `blast_db` path | BLAST as well, and every step after it |
| a new mitochondrion threshold | the mitochondrion screen and every step after it |

!!! warning "A BLAST database or a file-setting table changed in place needs `--forcerun`"

    A changed setting re-runs everything after it by itself. A BLAST database updated in place
    re-runs nothing, so add `--forcerun blast_contigs`. A table of the three file settings edited
    in place re-runs nothing either, because Snakemake sees the path, not the content. Add
    `--forcerun select_contigs` to apply the edited table.

## Coverage track

`blobtools create` requires a BAM file of reads aligned to the draft. Only the read-based modes
give a real coverage track.

| Mode | Reads aligned to the draft | Aligner | BAM after the run |
|---|---|---|---|
| `contigs` | none: the contigs are aligned to themselves | `minimap2 -a`, no preset | deleted once BlobTools has read it |
| `illumina` | the trimmed read pairs | bowtie2 | kept; Qualimap reads it |
| `longread`, `hybrid` | the filtered long reads | minimap2, `map-ont` or `map-hifi` per sample | kept; Qualimap reads it |

In the read-based modes, bowtie2 runs in its default end-to-end mode and minimap2 with the
sample's preset. Qualimap's report on the kept BAM is described on
[Assembly QC](assembly-qc.md#qualimap-read-modes).

!!! note "In `contigs` mode the coverage values mean nothing"

    In `contigs` mode the BAM exists because BlobTools needs one. Each contig is one "read" that
    maps once, to itself, and BlobTools computes coverage as mapped reads × mean read length /
    contig length. The value is therefore the mean contig length divided by the contig's own
    length, below 1 for contigs longer than the mean and above 1 for shorter ones. The contig
    selector reads only the genus, in every mode.

## BLAST against core_nt

`blast_contigs` runs one megablast search per sample, here with the default settings:

```bash
BLASTDB=/path/to/blast_db \
blastn -task megablast \
  -query 02.assembly/{sample}/contigs_filt.fasta \
  -db /path/to/blast_db/core_nt \
  -outfmt '6 qseqid staxids bitscore pident evalue length qlen slen qcovs qcovhsp sskingdoms scomnames sscinames sblastnames stitle' \
  -num_threads 24 -evalue 1e-5 \
  -max_target_seqs 50 -max_hsps 5 \
  -out 02.assembly/{sample}/contaminants/blastout
```

BLAST's megablast task looks for highly similar sequences and is much faster than the more
sensitive `blastn` task. `directories.blast_db` must hold the files below.

| File | Read by | If missing |
|---|---|---|
| the BLAST database named by `parameters.nt_version` (default `core_nt`) | blastn | BLAST fails |
| `taxdb.btd`, `taxdb.bti` | blastn, for the organism-name columns | those columns read `N/A` |
| `nodes.dmp`, `names.dmp` (NCBI taxdump) | BlobTools | the dry run stops with `MissingInputException in rule blob_json` |

Without `taxdb` the run goes on, with `N/A` for the organism names. The downloads are described on
[Reference databases](../getting-started/databases.md).

`parameters.nt_version` (default `core_nt`) is the name of the database inside
`directories.blast_db`, either `core_nt` or `nt`. Two more settings decide how much BLAST keeps per
contig.

| Key | BLAST option | Keeps | Default | Allowed | The dry run prints a NOTE when |
|---|---|---|---|---|---|
| `parameters.blast_max_target_seqs` | `-max_target_seqs` | database sequences ("subjects") per contig | 50 | integer ≥ 1 | below 5 |
| `parameters.blast_max_hsps` | `-max_hsps` | alignments per subject | 5 | integer ≥ 1 | below 5, unless `mitochondria.mode` is `off` |

Keep the defaults, 50 subjects with 5 alignments each per contig, because smaller values lose
evidence in two ways.

- BLAST applies `-max_target_seqs` during the search, not as a final best-N filter, so a small
  value can lose better matches. BLAST recommends 5 or more.
- The `qcovs` column (query coverage per subject), which the mitochondrion screen reads, is
  computed from the kept alignments only. With `-max_hsps 1`, BLAST keeps only the single best
  local alignment per database sequence, so a complete mitochondrial genome that aligns in several
  blocks can look only half covered.

The evidence behind these defaults, and the effect of the hit depth on the genus calls, are
described on [BLAST hit depth](../methods/blast-hit-depth.md).

BLAST is the slow step of the screen. On 24 threads it took 40 to 51 minutes for the 42.6 Mb
*Botrytis cinerea* B05.10 assemblies of the [validation runs](../about/validation.md). In
`contigs` mode, every other step of the screen takes seconds.

## BlobTools: one genus per contig

BlobTools joins the BLAST hits and the coverage track into one table, with a taxon per contig and
rank. `blob_json` runs `blobtools create` on the draft, the BAM, `blastout`, `nodes.dmp` and
`names.dmp`, and `blob_table` flattens the result:

```bash
blobtools view --input blob.blobDB.json --out bestscore --taxrule bestsum --rank all --hits
```

`bestsum` adds up the bit scores of all hits per taxon, separately at each rank, and the taxon
with the highest sum wins. `bestscore.blob.blobDB.table.txt` has one row per contig, with these
columns.

| Columns | Content |
|---|---|
| 1 to 5 | name, length, GC, N count and coverage (`bam0`) |
| 6 to 29, four per rank: superkingdom, phylum, order, family, genus and species | the winning taxon, its summed score, the number of other taxa with hits and, added by `--hits`, every competing taxon with its own sum |

The contig selector reads only the genus, column 22, and a contig without a BLAST hit has the genus
`no-hit`. The competing taxa show how close a call was.

## Output files

The steps on this page write these files to `02.assembly/{sample}/contaminants/`. Everything the
decision rests on is kept. Only the BlobTools intermediates and, in `contigs` mode, the BAM are
deleted.

| File | Kept | Content |
|---|---|---|
| `blastout` | yes | megablast hits, one line per alignment, the 15 columns of `-outfmt` [above](#blast-against-core_nt) |
| `bestscore.blob.blobDB.table.txt` | yes | BlobTools table, one row per contig |
| `contig_taxonomy_decisions.tsv` | yes | genus, action and reason per contig |
| `{sample}_composition.txt` | yes | share of bases and contigs per genus |
| `contigs.list` | yes | names of the kept contigs |
| `contigs_decontam.fasta` | yes | kept contigs with their original headers, 80 bases per line |
| `{sample}_map.bam`, `{sample}_map.bam.bai` | `contigs` mode: no; read-based modes: yes | coverage track; in the read-based modes Qualimap reads it too |
| `blob.blobDB.json`, `blob.{sample}_map.bam.cov` | no | BlobTools intermediates |

The full output tree is on [Output files](../reference/output.md).
