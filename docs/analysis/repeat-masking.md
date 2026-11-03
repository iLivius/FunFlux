# Repeat masking

FunFlux prepares the nuclear genome for gene prediction in two steps, the same in every mode.
First, funannotate removes short and duplicated contigs and sorts and renames the rest. Then the
repeats are soft-masked. The files of both steps are in `04.annotation/{sample}/`. Read this page
to find out why a contig can be missing from the annotation, how the annotated contigs are named,
and which masking method to choose.

```text
contigs_final.fasta
 └─ funannotate clean ── contigs_clean.fasta    short and duplicated contigs removed
     └─ funannotate sort ── contigs_sort.fasta  longest first, renamed scaffold_N
         └─ masking ── contigs_mask.fasta       repeats in lowercase; on to gene prediction
```

`funannotate clean` and `funannotate sort` run in one rule, and the masking in another:

| Rule | Tools | Reads | Writes |
|---|---|---|---|
| `funannotate_preprocess` | funannotate 1.8.17: `clean`, then `sort` | `02.assembly/{sample}/contigs_final.fasta` | `contigs_clean.fasta`, `contigs_sort.fasta`, `contigs_removed_by_clean.txt` |
| `repeat_masking` | tantan 51, or RepeatModeler 2.0.8 and RepeatMasker 4.2.3 | `contigs_sort.fasta` | `contigs_mask.fasta`, `repeatmasking/` |

`contigs_mask.fasta` goes on to [gene prediction](gene-prediction.md), where funannotate predict,
[Tiberius](gene-prediction.md#tiberius) and, when GeneMark-ES is on, its licence check read it.
Tiberius gets the same file so that the scaffold names match. Its `fungi` model was trained
without repeat masking, so it reads lowercase and uppercase letters alike, and for Tiberius the
genome is in effect unmasked.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **soft-masking** | writing repeats in lower case, so that gene finders can treat them differently without losing the sequence |
    | **`scaffold_N`** | the name `funannotate sort` gives each contig, `scaffold_1` being the longest; nothing is joined, and the gene models refer to these names |
    | **haplotig** | a contig holding the second copy of a heterozygous region, assembled apart from the first; `funannotate clean` drops such duplicates when shorter than the N50 |
    | **scratch** | a working file or folder that the step which wrote it deletes when that step ends |
    | **EVidenceModeler** | the program inside `funannotate predict` that builds one consensus gene per locus from the gene finders and protein alignments, weighting each source's votes |
    | **repeat library** | a FASTA of a genome's repeat families; RepeatModeler builds it from the genome itself, and RepeatMasker then masks the copies it finds |
    | **random seed** | the number a program's random choices start from; fixing it makes a later run repeat the same choices |

## funannotate clean

`funannotate clean` drops contigs before anything is annotated:

```bash
funannotate clean -i contigs_final.fasta -o contigs_clean.fasta --cpus <threads>
```

FunFlux keeps funannotate's defaults. A contig is dropped when:

| Condition | Option (default) |
|---|---|
| it is shorter than 500 bp | `-m 500` |
| it contains a character that is not an IUPAC nucleotide code | none |
| it contains fewer than four different letters | none |
| it is shorter than the assembly N50, and one minimap2 2.31 alignment to a contig of at least N50 length, or to a shorter contig, covers more than 95 % of its length at more than 95 % identity | `-p 95`, `-c 95` |

Only a contig shorter than the N50 can be dropped as a duplicate. Contigs at least as long as the
N50 are not tested for duplication, because FunFlux does not use funannotate's `--exhaustive`
option. The duplication test removes, for example, haplotigs that an assembly from a sequencing
provider may still contain. Left in, the genes on them would be annotated twice.

The step declares 2 threads, or 1 when `resources.threads` is 1, and runs that many minimap2
searches at a time (`--cpus`). Each search starts minimap2's default 3 threads, so the job uses a
little more CPU than it declares. The cleaned genome does not depend on the thread count.

When `funannotate clean` is killed half-way, it leaves a scratch folder, `clean_<random>`, in
`04.annotation/{sample}/`. The [cleanup step](functional-annotation.md#cleanup) removes it.

`logs/funannotate_preprocess_{sample}.log` gives the number of input contigs, the number of those
longer than 500 bp, and the N50. It also says how many contigs shorter than the N50 were checked
for duplication, how many were duplicated, and how many were written to `contigs_clean.fasta`.

!!! note "QUAST, BUSCO and the ITS search can count contigs that the annotation leaves out"

    `contigs_removed_by_clean.txt` lists, one name per line and in input order, every contig of
    `contigs_final.fasta` that is missing from `contigs_clean.fasta`. An empty file means that
    nothing was dropped. A dropped contig stays in `contigs_final.fasta`, so QUAST, BUSCO and the
    ITS search still count it, but it is not annotated.

## funannotate sort

Next, `funannotate sort` gives the contigs the names that the annotation uses:

```bash
funannotate sort -i contigs_clean.fasta -o contigs_sort.fasta --minlen 0
```

`sort` orders the contigs from longest to shortest and renames them `scaffold_1`, `scaffold_2`,
and so on. According to funannotate's help, the contigs are renamed because AUGUSTUS can have
problems with complicated contig names. `--minlen 0` drops nothing, and the sequences are not
changed.

`contigs_sort.fasta` is the genome that the annotation describes, so the assembly results and the
annotation name the same contig differently. Which name a contig has depends on the file:

| Files | Contig names |
|---|---|
| GFF3, GenBank and antiSMASH files | `scaffold_N`: every coordinate refers to a contig of `contigs_sort.fasta` |
| SignalP 6, InterProScan, eggNOG-mapper | none: they report per protein, by transcript ID (`SAMPLE1_000001-T1`) |
| QUAST, BUSCO and the ITS files | the original contig names |
| the contig summary | both: its `annotated_as` column links each contig to its `scaffold_N` ([Contig summary](contig-summary.md)) |

Since only names and order change, a contig can also be matched to its scaffold by sequence, or by
length where lengths are unique.

## Masking

Soft-masking writes repeats in lowercase instead of replacing them by N, so no sequence is lost.
AUGUSTUS and GeneMark-ES treat the lowercase stretches with caution. SNAP and GlimmerHMM, which
funannotate also runs, and Tiberius read lowercase and uppercase letters alike. After
EVidenceModeler, funannotate's repeat filter removes the gene models that lie for at least 90 % in
a masked repeat or whose protein matches a transposable element, whatever gene finder they came
from. See [What funannotate predict does](gene-prediction.md#what-funannotate-predict-does).

`parameters.masking_method` chooses the tools: tantan, which is fast and masks simple repeats
only, or RepeatModeler and RepeatMasker, which also mask transposons and other interspersed
repeats and take far longer.

| `masking_method` | Commands | What is masked | Run time | Threads |
|---|---|---|---|---|
| `tantan` (default) | `tantan contigs_sort.fasta > contigs_mask.fasta` | tandem repeats and low-complexity sequence | seconds | 1 |
| `repeatmodeler_repeatmasker` | `BuildDatabase -name {sample}_rm`, `RepeatModeler -database {sample}_rm -threads {threads}`, then `RepeatMasker -pa {threads} -xsmall -lib {sample}_rm-families.fa` | the repeat families RepeatModeler finds in this genome (interspersed repeats such as transposons), plus simple repeats and low-complexity sequence | tens of minutes or more | at most 12 ([below](#repeatmodeler-settings)) |

RepeatMasker masks with the library RepeatModeler built from the genome itself, so no Dfam or
RepBase library is needed. The library is `{sample}_rm-families.fa`, or RepeatModeler's
`consensi.fa.classified` when that file is missing or empty. When neither holds a sequence, the
rule stops with `RepeatModeler did not produce a usable repeat library.`, so an unmasked genome
never passes as a masked one.

The share of the genome that was masked is in the line `Genome loaded: ... repeats masked` of
funannotate predict's log and, with RepeatMasker, in its summary table
`repeatmasking/contigs_sort.fasta.tbl`.

### RepeatModeler settings

These keys take effect only with `masking_method: repeatmodeler_repeatmasker`.

| Key | Default | Passed as | Effect, from the RepeatModeler 2.0.8 help |
|---|---|---|---|
| `repeatmodeler_quick` | `true` | `-quick` | uses the smaller genome samples of RepeatModeler versions before 2.0.4: faster, with the sensitivity of those versions |
| `repeatmodeler_ltrstruct` | `false` | `-LTRStruct` | adds the LTR structural pipeline (LTR_Harvest and LTR_retriever) |
| `repeatmodeler_seed` | empty | `-srand <seed>` | fixes the seed of the random choice of genome samples; empty: a seed taken from the clock |

!!! note "Without `repeatmodeler_seed`, two runs on one genome mask slightly differently"

    RepeatModeler samples the genome at random, so two runs without a fixed seed build slightly
    different libraries, and the gene models and everything computed from them differ slightly
    too. RepeatModeler writes the seed it drew near the top of `logs/repeat_masking_{sample}.log`,
    in the line `Random Number Seed: <number>`. Its help describes `-srand` as the way to
    reproduce the sample choice between runs. Whether a fixed seed also gives a byte-identical
    masked genome has not been tested.

`repeat_masking` declares at most 12 threads and passes them to RepeatModeler as `-threads` and to
RepeatMasker as `-pa`. RepeatMasker's `-pa` counts parallel batches, and each batch can run
several RMBlast threads, so RepeatMasker may use more CPUs than 12.

## Outputs

`contigs_mask.fasta` goes on to gene prediction, and the other files record what was dropped,
renamed and masked.

| File, in `04.annotation/{sample}/` | Content | Read by |
|---|---|---|
| `contigs_clean.fasta` | the genome after `funannotate clean`, original contig names | you |
| `contigs_removed_by_clean.txt` | the contigs `clean` dropped; empty if none | you |
| `contigs_sort.fasta` | sorted and renamed `scaffold_N`: the genome that is annotated | `repeat_masking`, `contig_summary` |
| `contigs_mask.fasta` | `contigs_sort.fasta` with repeats in lowercase | `funannotate_prediction`, `genemark_check`, `tiberius_prediction` |
| `repeatmasking/` | deleted and created anew at every start of the rule; empty with tantan. With RepeatModeler: the BLAST database `{sample}_rm.*`, the library `{sample}_rm-families.fa` with its alignments `.stk`, RepeatModeler's working folder `RM_*`, and RepeatMasker's `contigs_sort.fasta.out` (every masked repeat), `.tbl` (masked bases by repeat class) and `.masked` | you |

The logs of the two steps are `logs/funannotate_preprocess_{sample}.log` and
`logs/repeat_masking_{sample}.log`.
