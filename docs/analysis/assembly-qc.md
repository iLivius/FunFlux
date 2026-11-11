# Assembly QC

Three tools check the quality of each sample's assembly. QUAST measures its size and contiguity,
and BUSCO its completeness. In the read-based modes, Qualimap also shows how the sample's reads map
back to the draft. Read this page to interpret their results, to choose the BUSCO lineage of a
batch and to find their files.

QUAST and BUSCO run on `02.assembly/{sample}/contigs_final.fasta`, the nuclear genome that ITS
extraction and funannotate also read. When the contamination and mitochondrion screen runs with
the default `parameters.mitochondria.mode: separate_confirmed`, the mitochondrion has been split
off before this point, as described in [Mitochondria](mitochondria.md). If it were left in, its
genetic code, gene content and copy number would distort both the contig statistics and the BUSCO
counts of a genome that is annotated as nuclear. In `contigs` mode without the screen,
`contigs_final.fasta` is the input FASTA after the contig filter.

| Tool | Version | Measures | Rule | Output, in `02.assembly/{sample}/` | If the tool fails |
|---|---|---|---|---|---|
| QUAST | 5.3.0 | size and contiguity | `genome_assembly_evaluation` | `eval/quast/` | a line in the log; the run continues |
| BUSCO | 6.0.0 | completeness of conserved single-copy genes | `completeness` | `eval/busco/` | the job fails and the run stops (`--keep-going` lets the other jobs finish) |
| Qualimap | 2.3 | reads mapped back to the draft (read-based modes) | `map_evaluation` | `eval/qualimap/` | a line in the log; the run continues |

Only MultiQC and the contig summary read these results. The contig summary counts the complete
BUSCO markers on each contig and is written after the annotation
([Contig summary](contig-summary.md)).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **nuclear genome** | the delivered genome without the confirmed mitochondrial contigs, `contigs_final.fasta`; it holds every contig when nothing is split off (`parameters.mitochondria.mode` `off` or `report`, or no screen) |
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **draft assembly (the draft)** | the assembly the screen examines, before any contig is removed: `contigs_filt.fasta` in `illumina` mode and in `contigs` mode with the screen, `assembly.fasta` in `longread` and `hybrid` |
    | **guide assembly** | in `hybrid` mode, the SPAdes assembly of the Illumina reads; it is screened, and only the pairs that map to it go on |
    | **read-based modes (read modes)** | the modes `illumina`, `longread` and `hybrid`, which start from reads, as against `contigs` |
    | **coverage track** | the reads aligned back to the draft, `contaminants/{sample}_map.bam`; BlobTools takes each contig's coverage from it, the contig summary its read depth and Qualimap its report; in `contigs` mode, the contigs aligned to themselves |
    | **lineage dataset** | the BUSCO gene set of one taxonomic group, such as `hypocreales_odb12`, for the completeness check; with `lineage: auto`, BUSCO chooses one by placing the genome |
    | **contig summary** | the table `02.assembly/{sample}/eval/contig_summary.tsv`: one row per contig of the draft, with its evidence (telomeres, depth, rDNA units, genes) and a verdict |

## QUAST

QUAST measures the size and contiguity of the genome. Rule `genome_assembly_evaluation` runs it:

```bash
quast contigs_final.fasta -o eval/quast --labels {sample} --no-icarus -t <threads>
# hybrid: the screened guide assembly first, then the final genome, two labels
quast guide/contaminants/contigs_decontam.fasta contigs_final.fasta -o eval/quast \
  --labels {sample}_guide,{sample} --no-icarus -t <threads>
```

QUAST uses `min(resources.threads, 24)` threads. It runs without a reference genome, so it reports
size and contiguity only, with no misassemblies, no genome fraction and no NG50. Contigs shorter
than 500 bp, QUAST's default minimum, are left out of every statistic except the rows marked
`(>= 0 bp)`. `--labels` names the assembly after the sample, so QUAST and MultiQC show the sample
name. `--no-icarus` leaves out QUAST's interactive contig viewer.

When QUAST fails, its log ends with
`[FunFlux] genome_assembly_evaluation: QUAST exited with an error for sample S. ...`,
`eval/quast/` stays in place, and the MultiQC page has no QUAST row for that sample.

In `hybrid` mode, one QUAST call compares the screened Illumina guide assembly with the final
genome, labelled `{sample}_guide` and `{sample}` in that order. The guide is the short-read view
of the same isolate, so the two rows show what the long reads added.

### What the numbers mean

QUAST writes its report as `eval/quast/report.tsv`, also as `.txt`, `.html`, `.pdf` and `.tex`,
and transposed as `transposed_report.*`. These are some of its rows:

| Row | Meaning |
|---|---|
| `# contigs (>= 0 bp)` | all contigs |
| `# contigs (>= 1000 bp)` to `(>= 50000 bp)` | contigs of at least 1, 5, 10, 25 and 50 kb |
| `Total length (>= 0 bp)` | sum of all contig lengths (bp) |
| `Total length (>= 50000 bp)` | bases in contigs of 50 kb or more |
| `# contigs` | contigs of 500 bp or more |
| `Largest contig` | length of the longest contig (bp) |
| `Total length` | bases in contigs of 500 bp or more |
| `GC (%)` | GC content |
| `N50` | contigs of this length or longer hold half of the total length |
| `N90` | the same for 90 % of the total length |
| `auN` | sum of the squared contig lengths divided by the total length: the mean length of the contig that a random base lies in |
| `L50` | the number of contigs, longest first, that hold half of the total length |
| `L90` | the same for 90 % |
| `# N's per 100 kbp` | gap characters |

A higher N50 and auN and a lower L50 mean fewer and longer contigs.

## BUSCO

BUSCO searches the genome for genes that nearly every member of a lineage carries exactly
once (single-copy orthologs) and reports how many it finds whole, duplicated, in pieces or
not at all.

```bash
busco --mode genome -i contigs_final.fasta --out {sample} --out_path eval/busco \
  --download_path eval/busco/busco_tmp --auto-lineage-euk --opt-out-run-stats --cpu <threads>
```

Rule `completeness` uses `min(resources.threads, 24)` threads, as QUAST does. With a fixed
lineage, `-l <lineage>` replaces `--auto-lineage-euk`, and `--offline --download_path 00.resources/busco`
replaces the download folder of the sample. In a eukaryote genome, BUSCO 6.0.0 finds the genes
with miniprot unless another predictor is requested. FunFlux requests none, and the summary says
`BUSCO was run in mode: euk_genome_min`.

### Automatic or fixed lineage

The key `parameters.busco.lineage` names the dataset BUSCO uses, or leaves the choice to BUSCO:

```yaml
parameters:
  busco:
    lineage: auto        # or a dataset name such as hypocreales_odb12
```

| Setting | What BUSCO runs | Summaries in `eval/busco/` | MultiQC rows per sample |
|---|---|---|---|
| `auto` (default) | eukaryota_odb12, then the most specific lineage it reaches by placing the genome in the eukaryote tree | `short_summary.generic.eukaryota_odb12.{sample}.txt` and `short_summary.specific.<lineage>.{sample}.txt` | 2 |
| a dataset name | that dataset only | `short_summary.specific.<lineage>.{sample}.txt` | 1 |

!!! warning "For a batch, a fixed lineage keeps the BUSCO numbers comparable"

    Each sample is placed on its own, so `auto` can choose different lineages for the samples of
    one batch, and their numbers are then not comparable. A fixed lineage gives every sample the
    same dataset, skips the eukaryota run and the placement, and leaves a smaller result folder.

Each summary also exists as `.json`. The log names the lineage that the placement chose, in the
line `Lineage <lineage> is selected, supported by <n> markers out of <m>`. A dataset name must
match `[a-z0-9_]+_odb[0-9]+` after lower-casing, or the dry run stops. `busco --list-datasets`
prints the available names, each with its version after a dot, as in `hypocreales_odb12.2`.
Leave out the dot and the number after it.

### Reading the summary

Each summary has one result line of this form:

```text
C:<C>%[S:<S>%,D:<D>%],F:<F>%,M:<M>%,n:<n>,E:<E>%
```

C is the completeness of the genome, and D, F, M and E say what lies behind it:

| Letter | Meaning | What an unusual value can mean |
|---|---|---|
| C | complete: found over its full length, S + D | a lower C after a decontamination policy removed contigs points to removed genome rather than removed contamination ([Decontamination](decontamination.md#log-messages)) |
| S | complete, one copy | |
| D | complete, two or more copies | a high D can come from a recent whole-genome duplication, from haplotypes assembled separately, or from a mixed culture of two fungi of the same lineage, in which each gene is present twice |
| F | fragmented: found in part | the gene is only partly present, or no full-length gene model was built although the gene may be there |
| M | missing | the gene is absent, or the search found no significant match |
| n | genes in the dataset | |
| E | complete genes whose miniprot prediction contains an internal stop codon | caution before using BUSCO's predicted sequences downstream; those genes still count as complete |

The readings of D (genome duplication, haplotypes), F, M and E are those of the BUSCO user guide.

The header of each summary names the dataset and its creation date, for example
`hypocreales_odb12 (Creation date: 2026-05-22, number of genomes: 76, number of BUSCOs: 4349)`.
With `auto`, BUSCO downloads the current dataset at every run, and a newer dataset can change the
numbers. Compare samples only across runs with the same dataset date. A fixed lineage is
downloaded once per output folder, so every later run into that folder uses the same release.

### Downloads, disk and network

The lineage setting also decides where the datasets come from, and when BUSCO needs the network:

| Setting | Where the datasets come from | Network |
|---|---|---|
| `auto` | every BUSCO job downloads its datasets into `eval/busco/busco_tmp/`, which is deleted when the job ends. Each sample downloads its own copy, because two BUSCO jobs downloading into one folder can corrupt it. The job downloads eukaryota_odb12, the placement files and the chosen lineage | in every BUSCO job, on the computer that runs it |
| a fixed lineage | rule `busco_data` downloads it once per output folder into `00.resources/busco/lineages/<lineage>/`, on the computer that runs Snakemake (`logs/busco_data.log`). Every BUSCO job reads it with `--offline` and writes nothing there, so the samples share one copy | once; the BUSCO jobs need none |

!!! warning "With `auto`, the computer that runs BUSCO needs internet in every job"

    A fixed lineage needs internet once, on the computer that runs Snakemake, and the BUSCO jobs
    need none. This is how BUSCO runs on cluster nodes without internet
    ([Network access](../getting-started/installation.md#network)).

With a fixed lineage, the JSON summary records `offline: True` and the path of the shared lineage
folder, `00.resources/busco/lineages/<lineage>/`. The `.txt` summary and `full_table.tsv` are the
same as with a download per job.

The odb12 datasets of 2026-05-22 have these archive sizes:

| Dataset | Archive size (MB) |
|---|---|
| eukaryota | 29.1 |
| fungi | 125.8 |
| ascomycota | 310.5 |
| sordariomycetes | 311.7 |
| hypocreales | 274.6 |
| basidiomycota | 151.2 |
| agaricomycetes | 196.7 |
| mucoromycota | 132.1 |

With `auto`, a job therefore downloads tens to hundreds of MB, every time. This is what stays on
disk afterwards:

| In `eval/busco/` | After the run |
|---|---|
| the summaries | kept |
| `run_<lineage>/`: `full_table.tsv` (one line per gene with its status and location), `missing_busco_list.tsv`, `miniprot_output/`, BUSCO's own copy of the summary | kept |
| `auto_lineage/`, `logs/`, `tmp/` | kept |
| every `hmmer_output/` and `busco_sequences/` folder | deleted |

The folder keeps hundreds of MB per sample.

FunFlux passes `--opt-out-run-stats`, so BUSCO 6.0.0 sends no usage statistics to
busco-data.ezlab.org and makes no ipinfo.io lookup. Without the flag, BUSCO would send run details,
results, a checksum of the input FASTA and the country of the machine's IP address. The flag
changes no result.

## Qualimap (read modes)

Qualimap shows how the sample's reads map back to the draft. Rule `map_evaluation` runs it in
`illumina`, `longread` and `hybrid` mode:

```bash
unset DISPLAY
JAVA_TOOL_OPTIONS="-Djava.awt.headless=true" \
qualimap bamqc -bam contaminants/{sample}_map.bam --java-mem-size=<min(ram_gb, 64)>G \
  -nt <min(threads, 24)> -outdir eval/qualimap -outformat html
```

The BAM is the coverage track of the screen, the sample's reads aligned to the draft before
decontamination ([Decontamination](decontamination.md#coverage-track)). Every read-based mode gets
a report on the reads that built the draft, `hybrid` a second one for the Illumina guide, and
`contigs` mode none:

| Mode | Rule | Reads in the BAM | Report in `02.assembly/{sample}/` |
|---|---|---|---|
| `illumina` | `map_evaluation` | the trimmed Illumina pairs (bowtie2) | `eval/qualimap/` |
| `longread` | `map_evaluation` | the filtered long reads (minimap2) | `eval/qualimap/` |
| `hybrid` | `map_evaluation` | the filtered long reads (minimap2) | `eval/qualimap/` |
| `hybrid` | `map_guide_evaluation`, the same command | the Illumina pairs mapped to the guide assembly (`guide/contaminants/{sample}_map.bam`) | `guide/eval/qualimap/` |
| `contigs` | none | its only BAM is the contigs aligned to themselves | |

The Java heap is the smaller of `resources.ram_gb` and 64 GB, and a changed `ram_gb` does not
re-run the report ([CPU cores and threads](../reference/running.md#cpu-cores-and-threads)). The two
settings in front of the command let Qualimap draw its plots on a server without a display. The
Java runtime notes them with one `Picked up JAVA_TOOL_OPTIONS` line in the log, which is expected.

When Qualimap fails, the log ends with
`[FunFlux] map_evaluation: Qualimap exited with an error for sample S. ...`,
`eval/qualimap/` stays in place, empty, and the MultiQC page has no mapping table for the sample.

The report gives the share of mapped reads, insert sizes of read pairs, mapping quality, GC
content, error rate, the depth along the draft and, in `genome_results.txt`, the mean depth of
every contig. The mitochondrion stands out by its depth, because a cell holds many copies of it.
A contaminant often has a depth different from that of the isolate.

## In the MultiQC report

MultiQC reads `eval/quast/report.tsv`, the BUSCO summaries and, in the read-based modes,
`eval/qualimap/` of every sample into `06.report/multiqc_report.html`. In `hybrid` it also reads
`guide/eval/qualimap/`. Each tool has one row per sample, with a second QUAST and Qualimap row for
the guide in `hybrid` and a second BUSCO row with the automatic lineage:

| Mode | QUAST rows | BUSCO rows | Qualimap rows |
|---|---|---|---|
| `contigs` | 1 | 2 with the automatic lineage, 1 with a fixed one | none |
| `illumina`, `longread` | 1 | as above | 1: `mapping QC | <sample>` |
| `hybrid` | 2 | as above | 2: `mapping QC | <sample>` and `mapping Illumina QC | <sample>` |

The report also holds the contig summary table of every sample, with one row per contig of the
draft assembly ([Contig summary](contig-summary.md#multiqc)). [Report](report.md) describes the
row names, the General Statistics table and the BUSCO folders that MultiQC skips.

## Re-runs

Snakemake deletes `eval/quast/`, `eval/busco/` and `eval/qualimap/` before it runs the step
again, so files you add to these folders are lost. These changes re-run the following steps:

| Change | Re-runs |
|---|---|
| `contigs_final.fasta`, for example after a new decontamination policy or mitochondrion threshold | QUAST and BUSCO |
| the BAM, that is, a new draft or new reads | Qualimap |
| `busco.lineage` | BUSCO, the contig summary and MultiQC only, after `busco_data` has downloaded a newly named lineage |

A new BUSCO lineage therefore re-runs only BUSCO and the two steps that read it, not the assembly
or the annotation.

## Output files

Everything these steps write is kept, except BUSCO's download folder, `busco_tmp/`, and the
`hmmer_output/` and `busco_sequences/` folders inside `eval/busco/`.

| Path | Kept | Read by |
|---|---|---|
| `02.assembly/{sample}/eval/quast/` | yes | MultiQC (`report.tsv`) |
| `02.assembly/{sample}/eval/busco/` | yes | MultiQC (the `short_summary.*.txt` files) |
| `02.assembly/{sample}/eval/busco/busco_tmp/` (`auto` only) | no | BUSCO (downloaded datasets) |
| `00.resources/busco/lineages/<lineage>/`, `00.resources/busco/file_versions.tsv` (fixed lineage only) | yes | BUSCO, in every later run into the output folder |
| `02.assembly/{sample}/eval/qualimap/` | yes (read-based modes) | MultiQC |
| `02.assembly/{sample}/guide/eval/qualimap/` | yes (`hybrid`) | MultiQC |
| `02.assembly/{sample}/eval/telomere_candidates.tsv`, `telomere_unit.txt`, `telomere_ends.tsv` | yes | the contig summary; you. The telomere repeat unit of the assembly, from tidk ([Contig summary](contig-summary.md#telomeres)) |
| `02.assembly/{sample}/eval/contig_summary.tsv`, `contig_summary_mqc.tsv` | yes | you; MultiQC (`contig_summary_mqc.tsv`). Written after the annotation ([Contig summary](contig-summary.md)) |
| `06.report/multiqc_report.html`, `06.report/multiqc_data/`, `06.report/multiqc_config.yaml` | yes | you |
| `logs/genome_assembly_evaluation_{sample}.log`, `logs/completeness_{sample}.log`, `logs/busco_data.log` (fixed lineage), `logs/map_evaluation_{sample}.log`, `logs/map_guide_evaluation_{sample}.log`, `logs/telomere_units_{sample}.log`, `logs/contig_summary_{sample}.log`, `logs/multiqc.log` | yes | you |

The full output tree is on [Output files](../reference/output.md).
