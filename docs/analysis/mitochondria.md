# Mitochondria

A fungal mitochondrial genome, usually some tens of kb and present in many copies per cell,
assembles as one more contig, or sometimes as two or more. This page describes how FunFlux
recognises these contigs, sets them apart from the nuclear genome and reports each call. Read it
to interpret `mitochondria/mitochondrial_screen.tsv` or to find out why no mitochondrion was split
off.

With the default `parameters.mitochondria.mode: separate_confirmed`, FunFlux finds the
mitochondrial contig with BLAST, protects it from removal by the contig selector and writes it to
`contigs_mitochondrial.fasta`. QUAST, BUSCO, the ITS search and funannotate then read only the
nuclear genome, `contigs_final.fasta`.

The mitochondrion is kept apart because of the genetic code. The mitochondria of most fungi use
genetic code 4, in which UGA codes for tryptophan. The mitochondria of *Saccharomyces cerevisiae*
and related yeasts use code 3, which changes AUA and CUN too. The nuclear gene predictor reads UGA
as a stop codon and calls false genes on a mitochondrial contig.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against NCBI core_nt, BlobTools, the mitochondrion screen and the contig selector |
    | **mitochondrion screen** | rule `screen_mitochondria`, part of the screen: it confirms mitochondrial contigs from BLAST hits to mitochondrial genomes that are not an animal's or a plant's and, in the long-read modes, the assembler evidence |
    | **draft assembly (the draft)** | the assembly the screen examines, before any contig is removed: `contigs_filt.fasta` in `illumina` mode and in `contigs` mode with the screen, `assembly.fasta` in `longread` and `hybrid` |
    | **contig selector** | rule `select_contigs`, the last step of the screen: it keeps or removes each contig by genus, as `parameters.decontamination` says; the default, `off`, keeps every contig |
    | **genus filter** | the setting `parameters.decontamination.mode: include` or `exclude`, under which the contig selector keeps or removes contigs by the genus BlobTools gave them |
    | **direct hit** | a BLAST hit to a database sequence whose title says it is a mitochondrial genome; only such hits count as BLAST evidence for the mitochondrion |
    | **qcovs** | the BLAST query coverage per subject: the percentage of a contig covered by all kept alignments to one database sequence |
    | **coverage ratio** | a contig's depth as the assembler reported it, divided by the median depth of all contigs |

## Where it runs

The mitochondrion screen is one step of the contamination and mitochondrion screen, called the
screen on this page. It runs wherever the screen runs, as listed under
[When the screen runs](decontamination.md#when-the-screen-runs).

| Mode | Mitochondrion screen |
|---|---|
| `illumina`, `longread`, `hybrid` | always |
| `hybrid`, Illumina guide assembly | its own mitochondrion screen, rule `screen_guide_mitochondria`, which decides from BLAST alone ([Hybrid](../modes/hybrid.md)) |
| `contigs` with `parameters.contigs.screen: true` | yes |
| `contigs` without the screen | none: nothing is split off, a mitochondrial genome in the input FASTA is annotated as nuclear, and the dry run prints a NOTE saying so |

The mitochondrion screen decides on the draft, before the contig selector runs, so that a genus
filter cannot remove the mitochondrion before it is set apart. The split comes last, after
polishing. The steps, in order:

| Step | Rule | Reads | Writes, in `02.assembly/{sample}/` |
|---|---|---|---|
| 1 | `blast_contigs` | the draft contigs | `contaminants/blastout` |
| 2 | `screen_mitochondria` | `blastout`, the draft contigs, the NCBI taxonomy file `nodes.dmp` of the BLAST database folder; in `longread` and `hybrid` also `assembly_evidence.tsv`, from `assembly_evidence` ([below](#with-assembler-evidence-longread-hybrid)) | `mitochondria/mitochondrial_screen.tsv`, `mitochondria/confirmed_mitochondrial.ids` |
| 3 | `select_contigs` | the confirmed IDs, among others | `contaminants/contigs_decontam.fasta`, with the confirmed contigs kept whatever their genus |
| 4 | polishing (long-read modes, where it runs), then `finalize_assembly` | the kept contigs | `contigs_all.fasta` |
| 5 | `separate_mitochondria` | `contigs_all.fasta`, the confirmed IDs | `contigs_final.fasta`, `contigs_mitochondrial.fasta` |

In the long-read modes, `longread` and `hybrid`, Medaka (ONT samples, unless
`medaka_model: false`) and Polypolish (`hybrid`) polish the kept contigs before
`finalize_assembly`. PacBio HiFi samples in `longread` are not polished.

## Classification

The mitochondrion screen gives every contig of the draft a class and a reason. BLAST is the
evidence in every mode, and the long-read modes add the assembler's evidence. A last check, on the
organism of the best direct hit, keeps the mitochondrial DNA of other organisms from being
confirmed ([below](#foreign-mitochondrial)).

### BLAST evidence

The mitochondrion screen reads `contaminants/blastout`, the BLAST table of the screen described
under [BLAST against core_nt](decontamination.md#blast-against-core_nt), and uses only direct
organelle hits. In BLAST's terms, each database sequence that a contig matches is a subject. A
direct hit is a hit whose subject title says that the subject is a mitochondrial genome. The title
must contain one of these phrases, in upper or lower case:

- `mitochondrion, complete genome`
- `mitochondrion, complete sequence`
- `mitochondrion, partial genome`
- `mitochondrion, partial sequence`
- `mitochondrial DNA` or `mitochondrial genome`, but not followed by a hyphen
- `organelle: mitochondrion`, NCBI's wording for the mitochondrial sequence of a genome assembly

A title with the word `mRNA` or `protein` never counts, because it describes a nuclear gene, its
transcript or its protein, even when the gene's name contains one of the phrases. For example,
`mitochondrial DNA-directed RNA polymerase, putative, partial mRNA` and `Increased loss of
mitochondrial DNA protein 1` are nuclear genes and do not count, nor does `mitochondrial ribosomal
protein L5 mRNA`.

For each contig, the mitochondrion screen takes the direct hit with the highest `qcovs`, the
percentage of the contig covered by all kept alignments to that subject.
`parameters.blast_max_hsps` sets how many alignments per subject BLAST keeps. With 1, BLAST keeps
only the single best local alignment per database sequence. A complete mitochondrial genome that
aligns in several blocks can then look only half covered, and a mitogenome rearranged relative to
its closest reference can fall below the threshold. In a run with the screen and with
`parameters.mitochondria.mode` other than `off`, the dry run prints a NOTE when `blast_max_hsps`
is below 5. `parameters.blast_max_target_seqs` (default 10) sets how many subjects BLAST keeps per
contig. BLAST applies that limit during the search, so it can miss a mitochondrial genome even when
that genome scores higher than the subjects it kept.

### Without assembler evidence (contigs, illumina)

In `contigs` and `illumina` mode, and on the `hybrid` guide, BLAST is the only evidence, so a
direct hit covering at least 35 % of the contig confirms it.

| Best direct hit | Class | Reason |
|---|---|---|
| `qcovs` ≥ `min_direct_blast_qcov` (35) | `confirmed_mitochondrial` | `strong_direct_blast` |
| `qcovs` below 35 | `mitochondrial_candidate` | `weak_direct_blast` |
| none | `not_mitochondrial` | `no_mitochondrial_evidence` |

`very_strong_qcov` and `min_assembler_coverage_ratio` change nothing in these modes, and the
coverage that SPAdes writes into its contig names is not used.

### With assembler evidence (longread, hybrid)

In the long-read modes, rule `assembly_evidence`, a short step right after the assembly, writes
`assembly_evidence.tsv` beside the assembly. It holds the read depth and the circularity of each
contig as the assembler reported them.

| Assembler | Depth | Circularity | Also reports |
|---|---|---|---|
| hifiasm | `rd:i` tag of the primary contig graph | the contig name ends with `c` (`ptg000001c`) | |
| Flye | `assembly_info.txt` | `assembly_info.txt` | repeat status, alternative-haplotype group |

The coverage ratio of a contig is its depth divided by the median depth of all contigs. A clean
circle is circular, not a repeat and in no alternative group. In these modes a direct hit of at
least 35 % needs the support of a high depth or a clean circle, unless it covers 70 % or more. The
tests run in this order:

| Evidence | Class | Reason |
|---|---|---|
| `qcovs` ≥ 35 and coverage ratio ≥ `min_assembler_coverage_ratio` (3.0) | `confirmed_mitochondrial` | `direct_blast_plus_high_coverage` |
| `qcovs` ≥ 35 and a clean circle | `confirmed_mitochondrial` | `direct_blast_plus_clean_circularity` |
| `qcovs` ≥ `very_strong_qcov` (70) | `confirmed_mitochondrial` | `very_strong_direct_blast` |
| any other direct hit | `mitochondrial_candidate` | `direct_blast_without_support` |
| no direct hit, circular and coverage ratio ≥ 3.0 | `mitochondrial_candidate` | `circular_and_high_coverage_without_blast` |
| no direct hit, circular | `circular_unclassified` | `circularity_only` |
| no direct hit, coverage ratio ≥ 3.0 | `high_coverage_unclassified` | `coverage_only` |
| none of the above | `not_mitochondrial` | `no_mitochondrial_evidence` |

Circularity or depth alone never confirms a contig. Only `confirmed_mitochondrial` contigs are
protected and split off, and the other classes change nothing in the assembly.

### Mitochondrial DNA of another organism { #foreign-mitochondrial }

A direct hit shows that a contig matches a mitochondrial genome, but not whose. A few human reads
in an Illumina library, for example, can assemble into short contigs of human mitochondrial DNA,
which match human mitochondrial genomes over their whole length and pass the tests above. Wherever
the screen runs, the guide's screen of `hybrid` mode included, the mitochondrion screen therefore
checks the organism of the best direct hit. It takes the hit's taxonomy ID from `blastout` and
works out from `nodes.dmp`, the NCBI taxonomy file in the BLAST database folder that BlobTools also
reads, whether that organism is a fungus, an animal or a plant. A contig whose best direct hit is
an animal's or a plant's mitochondrial genome gets this class instead of `confirmed_mitochondrial`
or `mitochondrial_candidate`:

| Best direct hit | Class | Reason |
|---|---|---|
| an animal (Metazoa) or a plant (Viridiplantae) | `foreign_mitochondrial` | `best_direct_hit_animal_or_plant` |

Such a contig is not confirmed, not split off and not protected from the genus filter. It stays in
the assembly, where the contig selector treats it like any other contig, so with the default
`parameters.decontamination.mode: off` it is kept
([Protected mitochondrial contigs](decontamination.md#protected-mitochondrial-contigs)), and the
[contig summary](contig-summary.md) gives it the verdict `foreign_mitochondrion`. Any other answer
leaves the class as the tables above give it: a fungus, but also an organism that the taxonomy
places in none of the three kingdoms, such as an environmental record filed under Eukaryota only,
or an oomycete, which NCBI does not count among the fungi. So does a hit without a usable taxonomy
ID (`0`, when the BLAST database carries none). The column `best_mito_kingdom` of the screen table
records the answer.

FunFlux 2.0.0 had no such check. In the `illumina` run of the
[validation page](../about/validation.md#other-results), it confirmed three contigs of human
mitochondrial DNA, of 572 to 1,139 bp, and put them in `contigs_mitochondrial.fasta` beside the
fungal mitochondrion. The `hybrid` run assembles the same three contigs from the same Illumina
reads, in its guide assembly, and there the screen of FunFlux 2.1.0 calls them
`foreign_mitochondrial`. The check left the fungal mitochondria of the four runs with ONT reads
confirmed as before.

## Settings

The mitochondrion screen has four settings, under `parameters.mitochondria`. Two act wherever the
screen runs, and two only where the assembler reports depth and circularity.

```yaml
parameters:
  mitochondria:
    mode: separate_confirmed          # off | report | separate_confirmed
    min_direct_blast_qcov: 35
    min_assembler_coverage_ratio: 3.0
    very_strong_qcov: 70
```

| Key | Default | Allowed | Used in |
|---|---|---|---|
| `mode` | `separate_confirmed` | `off`, `report`, `separate_confirmed` | every mode with the screen |
| `min_direct_blast_qcov` | 35 | number > 0 and ≤ 100 | every mode with the screen |
| `min_assembler_coverage_ratio` | 3.0 | number > 0 | `longread`, `hybrid` |
| `very_strong_qcov` | 70 | number > 0 and ≤ 100 | `longread`, `hybrid` |

With the screen on, the dry run prints the settings in use:

```text
Mitochondrion screen: separate_confirmed (direct BLAST coverage >= 35%, assembler coverage ratio >= 3, near-complete hit >= 70%).
```

## Modes

`parameters.mitochondria.mode` decides what FunFlux does with the findings of the mitochondrion
screen.

| `mode` | Screen table | `confirmed_mitochondrial.ids` | Protected from the contig selector | `contigs_mitochondrial.fasta` |
|---|---|---|---|---|
| `off` | every row `not_screened`, reason `screening_disabled` | empty | no | empty file |
| `report` | classes and reasons | empty | no | empty file |
| `separate_confirmed` (default) | classes and reasons | the confirmed contigs | yes | the confirmed contigs |

!!! warning "Only `separate_confirmed` protects the mitochondrion and sets it apart"

    `report` shows the calls without changing the assembly. With `off` and `report`,
    `contigs_final.fasta` holds every contig, a mitochondrion included, and it goes through gene
    prediction. `report` also leaves the mitochondrion exposed to `include`, `exclude` and
    `discard_no_hit`.

`screen_mitochondria` and `separate_mitochondria` run with every setting of `mitochondria.mode`, so
a finished sample has the same set of files. `mitochondria.mode: off` does not switch BLAST off,
because the rest of the screen still needs it.

The log of `screen_mitochondria` counts confirmed and selected contigs:
`Mitochondrial screen: N confirmed call(s); M selected for separation (mode=...).` When the screen
finds mitochondrial DNA of another organism, a second line names those contigs:
`Mitochondrial screen: N contig(s) of another organism's mitochondrial DNA (foreign_mitochondrial), not set apart: ...`

## The screen table

`mitochondria/mitochondrial_screen.tsv` has one row per draft contig, in FASTA order.
`classification` and `reason` give the call, and the other columns are the evidence behind it.

| Column | Content |
|---|---|
| `contig` | contig name, the first word of the FASTA header |
| `length` | contig length (bp) |
| `classification` | the class, see the tables above |
| `reason` | the test that decided |
| `best_mito_qcovs` | `qcovs` of the best direct hit (%); empty without a direct hit |
| `direct_mito_hits` | number of BLAST alignments to mitochondrial-genome subjects |
| `best_mito_title` | subject title of the best direct hit |
| `assembler_source` | `hifiasm` or `flye` |
| `assembler_circular` | `Y` or `N` |
| `assembler_coverage` | read depth reported by the assembler |
| `assembler_coverage_ratio` | depth divided by the median depth of all contigs |
| `assembler_repeat` | `Y` or `N` (Flye) |
| `assembler_alt_group` | alternative-haplotype group (Flye); `*` for none |
| `best_mito_slen` | length of the best hit's subject, the reference mitochondrial genome (bp) |
| `length_ratio` | `length` / `best_mito_slen`, two decimals |
| `concatemer_flag` | `Y` when `length_ratio` is above 1.5 and the contig is not recorded as circular, otherwise `N`; information only |
| `best_mito_kingdom` | `Fungi`, `Metazoa` or `Viridiplantae`: the kingdom of the organism of the best direct hit ([Mitochondrial DNA of another organism](#foreign-mitochondrial)); empty without a direct hit, when the hit has no usable taxonomy ID (`0`, or an ID missing from `nodes.dmp`), or when the taxonomy places it in none of the three |

The six `assembler_*` columns are empty in `contigs` and `illumina` mode.

## The split

`separate_mitochondria` reads `contigs_all.fasta` and the confirmed IDs and writes each contig to
one of two files. Order and headers are kept, and sequences are written 80 bases per line.

| File, in `02.assembly/{sample}/` | Content | Read by |
|---|---|---|
| `contigs_all.fasta` | the decontaminated assembly (polished, where polishing runs), mitochondrion included | `separate_mitochondria` |
| `contigs_final.fasta` | the nuclear genome | QUAST, BUSCO, ITS extraction, funannotate; the AMR screen when it is on |
| `contigs_mitochondrial.fasta` | the confirmed mitochondrial contigs; an empty file when there are none | the AMR screen when it is on, for cytochrome b ([Two runs](amr.md#two-runs)); the contig summary, for the contig names |

The [contig summary](contig-summary.md) reads the draft and the annotated genome,
`04.annotation/{sample}/contigs_sort.fasta`, instead of `contigs_final.fasta`.

The IDs are contig names, the first word of the FASTA header, so every step between the
mitochondrion screen and the split must keep them. `separate_mitochondria` stops before writing
anything in two cases:

| The rule stops when | Message |
|---|---|
| a confirmed ID is missing from `contigs_all.fasta` | `ERROR: Confirmed mitochondrial contig(s) are absent from final assembly: ...` |
| no nuclear contig would remain | `ERROR: Mitochondrial split would leave an empty nuclear assembly.` |

Otherwise the log ends with `Mitochondrial split: X mitochondrial and Y nuclear contig(s).`

FunFlux does not annotate the mitochondrial genome.

## Concatemers

An assembler can join two copies of the mitochondrial genome end to end. Such a contig matches
the reference over its whole length and is confirmed. The screen table shows it:

| `length_ratio` | Reading | `concatemer_flag` |
|---|---|---|
| about 1 | one copy | `N` |
| about 2 | two copies joined end to end | `Y` (above 1.5) |
| above 1.5, and the assembler reported the contig as circular | a closed circle longer than the reference: a larger genome, or two copies that the assembler joined into one circle (hifiasm did so with the HiFi reads of the [validation page](../about/validation.md#hifi-assembly): 49,194 bp against 24,712 bp) | `N` |

`Y` marks a contig that likely holds more than one copy, but only on a confirmed row. The reason
is that the flag is computed for every contig with a direct hit, and a long nuclear contig with a
short mitochondrion-like stretch also has a high ratio. In `contigs` and `illumina` mode
circularity is unknown, and the flag depends on the ratio alone.

!!! warning "FunFlux does not trim a concatemer to one copy"

    The flag changes nothing: the contig is split off as it is.

## When no mitochondrion is confirmed

`contigs_mitochondrial.fasta` is then empty, and a mitochondrion that was assembled stays in
`contigs_final.fasta` and is annotated as nuclear. The screen table shows why:

| Finding | Meaning | What to do |
|---|---|---|
| a `mitochondrial_candidate` row with reason `weak_direct_blast` | a direct hit covers less than `min_direct_blast_qcov` of the contig | if the length (usually some tens of kb) and the title fit, set `min_direct_blast_qcov` to the row's `best_mito_qcovs` or lower and run again. Or set `blast_max_target_seqs: 50`: BLAST applies that limit during the search, so 10 subjects can miss mitochondrial genomes that score higher and cover more of the contig |
| a `foreign_mitochondrial` row | the best direct hit is the mitochondrial genome of an animal or a plant (`best_mito_kingdom` is `Metazoa` or `Viridiplantae`) | read `best_mito_title`: a contig of another organism's mitochondrial DNA, such as human DNA from the library, is rightly left out of `contigs_mitochondrial.fasta`; no setting changes this check |
| no row with a direct hit | core_nt holds no mitochondrial genome close enough; or the subjects BLAST kept for the contig (10 by default) are records of single mitochondrial genes, such as cox1, which do not count; or a mitochondrial genome's title uses none of the wordings above | read the titles of the contig's hits, the last column of `contaminants/blastout` (a table without a header line). If they are single-gene records, set `blast_max_target_seqs: 50` and run again. Otherwise nothing can be set in `contigs` and `illumina` mode; in the long-read modes depth and circularity can only give a candidate or `..._unclassified` row |
| `blast_max_hsps` below 5 | `qcovs` counts only the alignments BLAST kept per subject, with 1 the best one alone | set `blast_max_hsps: 5` |

A weak hit, a low `blast_max_hsps` or too few subjects can be fixed with a setting, but a missing
reference genome cannot. A changed threshold re-runs the mitochondrion screen and every step after
it. A changed `blast_max_hsps` or `blast_max_target_seqs` re-runs BLAST first (in `hybrid` also the
guide's BLAST), and with a genus filter it can change the genus calls too
([BLAST hit depth](../methods/blast-hit-depth.md#ten-subjects-what-the-genus-calls-need)). Both are
described under [Changing the policy after a run](decontamination.md#changing-the-policy-after-a-run).

## Output files

Each sample keeps these files after the run:

| File | Kept | Content |
|---|---|---|
| `02.assembly/{sample}/mitochondria/mitochondrial_screen.tsv` | yes | one row per draft contig |
| `02.assembly/{sample}/mitochondria/confirmed_mitochondrial.ids` | yes | confirmed contig names, one per line; empty unless `separate_confirmed` |
| `02.assembly/{sample}/contigs_all.fasta` | yes | assembly with the mitochondrion |
| `02.assembly/{sample}/contigs_final.fasta` | yes | nuclear genome |
| `02.assembly/{sample}/contigs_mitochondrial.fasta` | yes | mitochondrial contigs |
| `logs/screen_mitochondria_{sample}.log`, `logs/separate_mitochondria_{sample}.log` | yes | counts and errors |

Without the screen (`contigs` mode by default), only `contigs_final.fasta` exists of these. The
full output tree is on [Output files](../reference/output.md).
