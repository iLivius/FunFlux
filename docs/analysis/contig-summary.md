# Contig summary

A fungal assembly is a list of contigs, and the list alone does not say what each contig is.
For every sample, FunFlux writes one table that does:
`02.assembly/{sample}/eval/contig_summary.tsv`. It has one row per contig of the draft assembly,
with the evidence the run produced about that contig and a verdict with its reasons. The table
also appears in the MultiQC report. Read this page to interpret the verdicts, to see how the
telomere repeat unit is chosen, and to set the unit yourself when the data cannot choose it.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **draft assembly (the draft)** | the assembly before the screen removes any contig, such as `contigs_filt.fasta` or `assembly.fasta`; the contig summary describes it |
    | **the screen** | the contamination and mitochondrion screen of the draft: BLAST against core_nt, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **verdict** | the label the contig summary gives each draft contig, one of eleven such as `chromosome`, `chromosome_arm`, `rDNA_array` or `mitochondrion`, with its reasons |
    | **telomere repeat unit (telomere unit)** | the short sequence repeated in tandem at chromosome ends, `TTAGGG` in most filamentous fungi; FunFlux chooses it from each assembly unless `parameters.contig_summary.telomere_motif` names one |
    | **telomere array** | tandem copies of the telomere repeat unit at a contig end; by default a telomere needs 5 or more copies within 25,000 bp of the end |
    | **canonical form** | the spelling tidk gives a repeat unit: the alphabetically first of its rotations and those of its reverse complement, so `TTAGGG` is written `AACCCT` |
    | **rDNA unit** | one copy of the repeat SSU-ITS1-5.8S-ITS2-LSU; FunFlux counts one per 5.8S gene found, that is, one per ITS window |
    | **BUSCO marker** | a BUSCO gene found complete on a contig (`Complete` or `Duplicated`); the contig summary counts them per contig as evidence of nuclear sequence |

## What the contigs of a fungal assembly are

A long-read assembly of a filamentous fungus holds a few kinds of contig side by side, and a
short-read assembly adds two more kinds, by the hundred. Each kind has its own signature in the
evidence FunFlux already produces, as the table below shows. Telomeres, the first piece of
evidence, are tandem arrays of a short repeat unit that cap every chromosome end. The unit is
`TTAGGG` in most filamentous fungi. The left end of a chromosome reads `CCCTAA` on the forward
strand, and the right end reads `TTAGGG`.

| Kind of contig | Telomere arrays | Read depth | Genes and BUSCO markers | Other signs |
|---|---|---|---|---|
| **whole chromosome**, telomere to telomere | at both ends | the depth of the bulk of the genome | along its whole length | size says nothing: an accessory chromosome or a minichromosome of 200 to 300 kb is a chromosome like a 4-Mb one |
| **chromosome arm** | at one end only; the other end is where the assembler stopped, usually at a repeat it could not resolve: the rDNA array, a centromere | as a chromosome | as a chromosome | |
| **mitochondrion**, a small genome (82 kb in *Botrytis cinerea* B05.10, the [worked example](#worked-example)) present in many copies per cell | | several times the nuclear depth | | a long-read assembler may close it into a circle and flag it as circular, but it may also join copies of it end to end or leave it in pieces ([Concatemers](mitochondria.md#concatemers)) |
| **piece of the rDNA array** | can carry one at one end, when the array sits at a chromosome end | says little: the reads of all the copies are shared out over however many pieces the assembler made of the array | no BUSCO marker; few or no protein-coding genes | several 5.8S genes in a row (the ITS windows of the taxonomy step, one per rDNA unit) |
| **alternative copy of a region, or noise**: a second haplotype assembled on its own, a stretch assembled twice, or a leftover of the assembly graph | | well below the genome's, because the reads of the region were split between two copies, or were few | | |
| **contaminant**: a contig from another organism | | its own, not the genome's | | its own GC content; the screen labelled it with a genus from its BLAST hits and removed it |
| **another organism's mitochondrial DNA**, most often human DNA from the library | | its own, usually low | | its best mitochondrial BLAST hit is an animal's or a plant's; the mitochondrion screen calls it `foreign_mitochondrial` ([Mitochondria](mitochondria.md#foreign-mitochondrial)) |
| **collapsed repeat** (short-read assemblies) | | well above the genome's: the reads of all the copies pile up on it | no BUSCO marker, at most a gene call or two (predictors call open reading frames inside transposons) | several genome copies of a transposon, a satellite or another repeat family, assembled into one contig because short reads cannot tell the copies apart |
| **short piece of the nuclear genome** (short-read assemblies) | | the depth of the genome | none, because it is too short to carry a gene | hundreds to a few thousand bases; nuclear sequence all the same, left over where the assembly broke |

An empty cell means that the evidence is not part of that kind's signature. Taken together, the
columns sort the contigs. Telomeres tell a chromosome from an arm. Depth sets the mitochondrion,
collapsed repeats and alternative copies apart from the nuclear genome. Genes and BUSCO markers
tell nuclear sequence from rDNA and debris. A few points the table cannot hold:

- With the default settings, the mitochondrion screen confirms the mitochondrion from a BLAST hit
  to a mitochondrial genome that covers at least 35 % of the contig, unless that genome is an
  animal's or a plant's. In `contigs` and
  `illumina` mode, and on the Illumina guide assembly of `hybrid` mode, that hit is enough. On the
  long-read assembly of `longread` and `hybrid` mode, the hit must also be backed by the assembler's
  evidence, unless it covers 70 % or more. That support is an assembler depth of at least 3 times
  the median of all contigs, or a clean circle: a contig recorded as circular, not a repeat and in
  no alternative group. By default, the confirmed mitochondrion is then split off before the
  annotation. See [Mitochondria](mitochondria.md#classification).
- The ribosomal RNA genes sit in tens to hundreds of tandem copies of one unit,
  SSU-ITS1-5.8S-ITS2-LSU, that no assembler resolves whole.
- A long-read assembler resolves most repeats in place, so a long-read assembly has few collapsed
  repeats.

The contig summary gathers the evidence behind these signatures into one row per contig, and says
for each contig how it was weighed. The telomere arrays are found with tidk, in the sequence
itself. The other evidence is the read depth, the assembler's depth and circularity flag, the rDNA
units, the BUSCO markers, the predicted genes, and the calls of the contamination and mitochondrion
screen.

## The table

`02.assembly/{sample}/eval/contig_summary.tsv` has one row per contig of the draft assembly. Its
`verdict` and `reason` columns give the answer, `annotated_as` gives the contig's name in the
annotation, and the other columns hold the evidence:

| Column | Content | `NA` when |
|---|---|---|
| `contig` | the name in the draft assembly: the first word of the FASTA header, as the screen tables use it | never |
| `length` | bp | never |
| `gc` | GC fraction, 3 decimals | never |
| `depth` | mean read depth, measured with samtools coverage (`meandepth` of `contaminants/{sample}_depth.tsv`), 2 decimals | `contigs` mode; a contig without a row in the table |
| `depth_x_median` | `depth` divided by the median depth of the nuclear genome, 2 decimals ([the median](#median-depth)) | `depth` is `NA` |
| `assembler_depth` | the assembler's own depth, 1 decimal | `illumina` and `contigs` mode; an empty cell in `assembly_evidence.tsv` |
| `circular` | `Y` or `N`, the assembler's flag | `illumina` and `contigs` mode |
| `telomere_start`, `telomere_end` | copies of the telomere unit in the array at that end: the array length tidk measured, divided by the unit length and rounded; only an array that points the way a telomere points counts ([Telomeres](#telomeres)) | never (`0` when none) |
| `telomere_end_offset_start`, `telomere_end_offset_end` | distance (bp) from the contig end to the array | no telomere array at that end |
| `telomere_unit` | the telomere repeat unit of the assembly in tidk's canonical form (`AACCCT` for `TTAGGG`), the same on every row of a sample; `none` when no unit was found | never |
| `rdna_units` | ITS windows on the contig, that is, 5.8S genes found | never (`0` when none) |
| `busco_markers` | complete BUSCO markers on the contig (status `Complete` or `Duplicated`; a duplicated marker counts where each copy lies) | never (`0` when none) |
| `genes` | mRNA features of the annotation on the scaffold the contig became | the contig has no scaffold in `contigs_sort.fasta`: dropped by `funannotate clean`, split off as the mitochondrion, or removed by the screen |
| `genus` | `assigned_genus` of `contig_taxonomy_decisions.tsv` | no screen |
| `screen` | `kept (<reason>)`, `removed (<reason>)` or `mitochondrion (<reason>)`, with the reason of the contig selector or, for a mitochondrion, of the mitochondrion screen; a `mitochondrial_candidate`, `foreign_mitochondrial`, `circular_unclassified` or `high_coverage_unclassified` class of the mitochondrion screen is appended after `;` | no screen |
| `adapter` | what the contig adapter screen found: `none`, `at_end` (adapter or barcode within 100 bp of a contig end), `internal` (a whole adapter inside the contig) or `at_end+internal`; it changes no verdict ([The contig adapter screen](../methods/adapter-trimming.md#the-contig-adapter-screen)) | never |
| `annotated_as` | `scaffold_N`, the name of the contig in every annotation file; or `dropped by funannotate clean`, `not annotated (mitochondrion)`, `not annotated (removed by the screen)` | never |
| `verdict` | one of the eleven verdicts ([The verdict](#the-verdict)): `mitochondrion`, `contaminant`, `foreign_mitochondrion`, `rDNA_array`, `chromosome`, `chromosome_arm`, `low_depth`, `repeat`, `nuclear_contig`, `nuclear_fragment`, `unclassified` | never |
| `reason` | the evidence behind the verdict, in words, for example `telomere repeats at both ends (17 and 19 copies); depth 1.00x the median; 543 genes, 131 BUSCO markers` | never |

The table describes the draft assembly, which is `contigs_final.fasta` only in `contigs` mode
without the screen. Which file is the draft depends on the mode:

| Mode | The draft |
|---|---|
| `illumina`; `contigs` with the screen | `contigs_filt.fasta` |
| `longread`, `hybrid` | `assembly.fasta` |
| `contigs` without the screen | `contigs_final.fasta` |

!!! note "The contigs the screen removed and the mitochondrion have a row too"

    Because the table describes the draft, the contigs the screen removed and the mitochondrion
    get a row too, and nothing disappears silently. Evidence a mode cannot produce is `NA`, and
    the verdict rests on the rest.

## The verdict

Eleven verdict rules give each contig its verdict. They are tried in the order of the table below, and
the first that applies decides. Evidence a mode cannot produce is unknown (`NA` in the table) and
is not tested, so "depth within 0.5 to 1.5 times the median when known" passes when the depth is
unknown. A telomere is a contig end where tidk found an array of the telomere unit of at least
`min_telomere_copies` copies, pointing the way a telomere points.

| Rule | Verdict | The contig has | What it is |
|---|---|---|---|
| 1 | `mitochondrion` | the class `confirmed_mitochondrial` from the mitochondrion screen, or a place in `contigs_mitochondrial.fasta` | the mitochondrial genome |
| 2 | `contaminant` | the action `remove` from the screen; the reason names the genus | a contig of another organism |
| 3 | `foreign_mitochondrion` | the class `foreign_mitochondrial` from the mitochondrion screen, and a place in the assembly (the contig selector kept it; with the default `decontamination.mode: off` it keeps every contig) | mitochondrial DNA of an animal or a plant, most often human DNA from the library |
| 4 | `rDNA_array` | at least `parameters.contig_summary.min_rdna_units` rDNA units (3 by default), and fewer predicted genes than rDNA units, or no gene count | a piece of the rDNA array |
| 5 | `chromosome` | telomere repeats at both ends; depth within 0.5 to 1.5 times the median when known; at least one gene when known | a whole chromosome; size plays no part, so a small accessory chromosome qualifies |
| 6 | `chromosome_arm` | telomere repeats at one end; the rest as in rule 5 | a chromosome arm |
| 7 | `low_depth` | depth below 0.5 times the median | an alternative copy of a region (a second haplotype, a stretch assembled twice) or noise |
| 8 | `repeat` | no BUSCO marker; at most two predicted genes; depth above 1.5 times the median (when the depth is unknown this rule does not fire) | a collapsed repeat: several genome copies of a transposon, a satellite or another repeat family assembled into one contig, so the reads of all the copies pile up on it |
| 9 | `nuclear_contig` | at least one gene or BUSCO marker; depth within 0.5 to 1.5 times the median when known; no telomere repeat | a fragment of a chromosome |
| 10 | `nuclear_fragment` | no gene and no BUSCO marker; no telomere repeat; depth within 0.5 to 1.5 times the median, or unknown | a piece of the nuclear genome too short to carry a gene; a short-read assembly ends in many of them |
| 11 | `unclassified` | anything else | the reason says what was and was not seen |

The screen decides first, then the rDNA count, then the telomeres, and depth and genes sort out
the rest. Two verdict rules allow a gene call where none is expected. In verdict rule 4, a pure
rDNA array holds no protein-coding gene, but the gene predictor now and then calls a gene inside
an rDNA unit. A contig with more rDNA units than genes is still mostly rDNA. In verdict rule 8, a
gene call or two is allowed because gene predictors call open reading frames inside transposons.

"No BUSCO marker" (verdict rule 8) and "no gene and no BUSCO marker" (verdict rule 10) apply only
when at least one of the two counts is known, and then every known count they name must be `0`.
In a run both are always counted, because BUSCO and
funannotate run in every mode. Only a contig without a scaffold in `contigs_sort.fasta` lacks its
gene count.

The `reason` column spells out the evidence each verdict rule used, with the numbers, so a verdict
can be checked against the columns.

### The median depth { #median-depth }

The median depth is what `depth_x_median` and the verdict rules compare each contig with. It is
weighted by length, so a 4-Mb chromosome counts a hundred times as much as a 40-kb piece. It is
taken over the contigs of at least 100 kb, or over all contigs when none is that long, and it
leaves out the mitochondrion and the contigs the screen removed. Short debris, the deep
mitochondrion and a contaminant genome at its own depth therefore cannot pull the median away from
the depth of the bulk of the genome.

## How the table is made

### Where it runs

Three steps make the table, in every mode, and in the read-based modes a fourth measures the read
depth. The telomere unit is chosen as soon as the draft exists, and the table itself waits for the
annotation:

| Step | Runs | Does | Writes |
|---|---|---|---|
| `telomere_db` | once per output folder | downloads the tidk telomere repeat database and checks its SHA-256 ([Reference databases](../getting-started/databases.md#telomere-db)) | `00.resources/telomeres/telo_database_v2.tsv`, `PROVENANCE.txt` |
| `telomere_units` | once per sample, as soon as the draft exists | runs tidk on the draft and chooses the telomere repeat unit of the assembly ([Telomeres](#telomeres)); its own conda environment (tidk 0.3.0 with Python 3.12); seconds to a minute | `eval/telomere_candidates.tsv`, `telomere_unit.txt`, `telomere_ends.tsv` |
| `contig_depth` | once per sample, in the read-based modes, as soon as the reads are aligned to the draft | measures the mean read depth of every contig with samtools coverage, from the alignment of the sample's reads to the draft that the screen makes ([Coverage track](decontamination.md#coverage-track)) | `contaminants/{sample}_depth.tsv` |
| `contig_summary` | once per sample, after the annotation; one of the last jobs of a sample | gathers the evidence into the table | `eval/contig_summary.tsv`, `contig_summary_mqc.tsv` |

The MultiQC report reads the `_mqc` copy of the table, and no other step reads it.

### Where the evidence comes from

Telomeres, rDNA units, BUSCO markers and genes are known in every mode. The read depth is known
only in the read-based modes, the assembler's depth and circularity only in the long-read modes,
and the calls of the screen only where the screen runs. Each piece of evidence is read from these
files:

| Evidence | What it says | Taken from | Modes |
|---|---|---|---|
| telomere arrays at the contig ends | at both ends: a whole chromosome; at one end: a chromosome arm | `eval/telomere_unit.txt` and `eval/telomere_ends.tsv`: tidk on the draft sequence itself ([Telomeres](#telomeres)) | all |
| read depth | the nuclear contigs share one depth; the mitochondrion is several times deeper; an alternative copy of a region, or noise, is much shallower | `contaminants/{sample}_depth.tsv`, column `meandepth`: samtools coverage of the sample's reads on the draft ([Coverage track](decontamination.md#coverage-track)) | `illumina`, `longread`, `hybrid` |
| assembler depth and circularity | a circular contig at a high depth is the usual shape of a mitochondrion | `assembly_evidence.tsv`, columns `coverage` and `circular` | `longread`, `hybrid` |
| rDNA units | a contig made of rDNA units is a piece of the rDNA array, not a chromosome | `03.taxonomy/{sample}/its_windows.tsv`, one row per 5.8S gene found ([ITS taxonomy](its-taxonomy.md#step-1-one-window-per-rdna-copy)) | all |
| BUSCO markers | genuine nuclear sequence, as against rDNA and debris | `eval/busco/run_<lineage>/full_table.tsv` of the specific lineage, the markers with status `Complete` or `Duplicated` ([BUSCO](assembly-qc.md#busco)) | all |
| predicted genes | the same | the mRNA features of `04.annotation/{sample}/funannotate/annotate_results/<Genus_species>_{sample}.gff3`, joined to the draft contig through `04.annotation/{sample}/contigs_sort.fasta` ([the join](#genes)) | all |
| the screen | removed as a contaminant, with its genus; confirmed mitochondrial | `contaminants/contig_taxonomy_decisions.tsv`, `mitochondria/mitochondrial_screen.tsv`, `contigs_mitochondrial.fasta` | screen |

Paths without a stage prefix are under `02.assembly/{sample}/`. In `hybrid` mode the screen files
are those of the long-read draft, not of the Illumina guide.

In `contigs` mode the screen maps the contigs to themselves, so a depth would reflect the
self-alignment and nothing else ([Coverage track](decontamination.md#coverage-track)). The depth
is therefore not measured there, and the depth columns are `NA` in `contigs` mode, with or without
the screen.

!!! note "The read depth comes from samtools coverage, not from BlobTools"

    FunFlux 2.0.0 took the read depth from the coverage column of the BlobTools table (`bam0`).
    BlobTools estimates that column as the number of reads whose main (primary) alignment lies on
    the contig, times the mean length of all alignments in the BAM file, divided by the contig
    length. In that mean, secondary alignments, the other places where a read could also align,
    count as length 0 because they carry no sequence, and supplementary alignments, the further
    pieces of a read that aligns in parts, count only their aligned part. On long reads the
    estimate can therefore run low. On the *B. cinerea* B05.10 ONT reads of the
    [worked example](#worked-example), it came out almost five times too low on the untrimmed
    reads and 12 % too low on trimmed ones. A contig that reads reach only through supplementary
    alignments, such as a short piece of a mitochondrion, got a depth of 0. samtools coverage
    instead adds up the aligned bases at every position of the contig, with supplementary
    alignments counted and secondary ones left out. Run again with this change on their own
    files, the four 2.1.0 runs with ONT reads of the [validation page](../about/validation.md)
    kept every verdict. The BlobTools table is still made, and the contamination screen still
    takes its genus calls from it.

### Telomeres

The telomere repeat unit differs between fungal groups
([the table below](#which-fungi-carry-which-repeat)), so FunFlux does not assume one. It chooses
the unit from the assembly itself and says whether the unit it found is a known telomere repeat.
Rule `telomere_units` runs tidk (Brown et al. 2025) on the draft assembly, and tidk does the heavy
work. FunFlux uses two of its commands:

| tidk command | What it does | Used for |
|---|---|---|
| `tidk explore` | lists every tandem repeat unit of a sequence set with its copy number | collecting the candidate units (step 2 below) |
| `tidk ends` | takes one unit and reports, for both ends of every contig, whether an array of that unit sits there, how long it is, how far from the contig end, and which way it points | scoring each candidate (step 4) and the telomere columns of the table |

!!! note "tidk measures the arrays, and FunFlux chooses the unit"

    Neither command chooses the telomere unit of a genome. Most units `tidk explore` lists are
    microsatellites, and `tidk ends` needs the unit as input. The choice is FunFlux's, made from
    how each candidate unit behaves at the contig ends, and checked against the telomere repeat
    database that comes with tidk and against the literature.

tidk writes every unit in a canonical form, the alphabetically smallest of all rotations of the
unit and of its reverse complement, so `TTAGGG`, `GGGTTA` and `CCCTAA` are all `AACCCT`. The
`telomere_unit` column, the three telomere files and the log use that form, and the table of
database units under [Which fungi carry which repeat](#which-fungi-carry-which-repeat) gives the
published spelling next to it.

#### How the unit is chosen { #telomere-unit }

1. **The contig ends.** FunFlux writes a FASTA of terminal windows. Of every contig longer than
   two windows (`parameters.contig_summary.telomere_window`, 25,000 bp by default), it takes the
   first and the last window, as `<contig>__start` and `<contig>__end`. A contig of at most two
   windows goes in whole. The window is wide because the array is not always the last thing in a
   contig. When the assembler ran past the telomere, the array sits some kilobases inside the
   contig end, and a window that stops short of it misses the telomere altogether.

2. **The units found there.** `tidk explore` lists every unit of 5 to 30 bp that forms a tandem
   run of at least `min_telomere_copies` copies in the windows, with its copy number
   (`copies_at_ends`). It runs once more on the whole draft (`copies_genome_wide`). A telomere
   unit has nearly all its copies at the ends, while a satellite has most of them inside the
   contigs. Both numbers go into the candidates table, but neither enters the score, because a
   microsatellite can have nearly as many copies at the contig ends as the telomere unit. The copy
   number also does not say which way its arrays point.

3. **The candidates.** They are the units found at the ends, every unit the database lists with
   `strong` or `weak` support, every unit of the literature table and, when
   `parameters.contig_summary.telomere_motif` names one, that unit. All are written in canonical
   form and merged. Three filters drop what cannot be a telomere repeat, and the configured unit
   is exempt from them:

    | Filter | Drops, for example | Why |
    |---|---|---|
    | 5 to 30 bp long | `ACAT`, `AAAT` | no telomere repeat shorter than 5 bp is known, and 4-bp units are low-complexity sequence that forms arrays everywhere: `AAAT` alone gives 192 "telomere" ends on the Illumina (SPAdes) assembly of *Botrytis cinerea* B05.10 |
    | not a repeat of a shorter unit | `ACACAC` (`AC` three times), `TTAGGGTTAGGG` | the same array would be counted twice under two names |
    | the two strands differ: the unit holds more G than C, or more C than G | `AAAT`, `AGAGCC` | every telomerase-made repeat has a G-rich strand, and the orientation test of step 4 needs the two strands to be told apart |

4. **The score.** For every candidate, `tidk ends` runs on the whole draft. It searches the
   window at each contig end and accepts an array no further from the end than the window
   (`-w` and `--max-offset`, both `telomere_window`). An array must be at least
   `min_telomere_copies` copies long (`--min-length`, that is 30 bp for a 6-bp unit), and tidk's
   default mode tolerates sequencing errors inside the array. Each contig end gets a status:

    | Status | Meaning |
    |---|---|
    | `present` | an array pointing the way a telomere points: the C-rich strand at a contig start, the G-rich strand at a contig end |
    | `wrong_strand` | an array pointing the other way |
    | `not_terminal` | an array further from the end than the window |
    | `absent` | no array |
    | `other_end` | FunFlux's addition: on a contig shorter than two windows, one array is seen from both ends, once as `present` at its own end and once as `wrong_strand` from the other end; when the two rows describe the same stretch of the contig (their coordinates overlap), the row with the larger offset is relabelled and does not count |

    Per unit, `telomere_ends` is the number of `present` rows, `opposite_ends` the number of
    `wrong_strand` rows, `score` their difference, and `t2t_contigs` the contigs with `present`
    at both ends. The reasoning is that a telomere repeat points one way at every chromosome end,
    with its G-rich strand towards the end. A satellite that happens to sit at contig ends points
    either way with equal chance, so however often it occurs, its score stays near zero.

5. **The choice.** A unit is eligible with at least 2 telomere ends and more than twice as many
   telomere ends as opposite ends. The eligible units are ranked by score, then known units
   (literature or database) before unknown ones, then the shorter unit, then alphabetically. The
   first is the telomere unit of the assembly. When none is eligible, the unit is `none`, every
   telomere column reads `0`, and a chromosome comes out as `nuclear_contig`. The log then names
   the best candidate and its counts.

#### Which fungi carry which repeat

`TTAGGG`, the vertebrate telomere repeat, is the unit of most filamentous Ascomycota (the
Pezizomycotina), of most Basidiomycota and of the Mucoromycota. The departures are single species
or small groups, and the budding yeasts, the fission yeasts and the Microstromatales have no single
unit. The full references are on [Citation and references](../about/citation.md#references). The
table lists the known departures and the groups that no single unit describes:

| Fungi | Telomere repeat | Source |
|---|---|---|
| most Pezizomycotina (filamentous Ascomycota), most Basidiomycota, Mucoromycota; among them *Botrytis cinerea* and *Aspergillus fumigatus* | `TTAGGG` | Lue 2021; Sun et al. 2023; Levis et al. 1997; Bowyer et al. 2022 |
| the *Aspergillus flavus* group | `TTAGGGTCAACA` | Sun et al. 2023 |
| *Peltaster fructicola* | `TAGGG` | Sun et al. 2023 |
| *Drechmeria coniospora* | `CCGTTGCTGTTG` | Sun et al. 2023 |
| *Penicillium polonicum* | `TTTAGGG` | Sun et al. 2023 |
| a few Basidiomycota: *Malassezia*, *Hericium* | `TTAGGA` (*Malassezia furfur*), `TTAGTG` (*Malassezia restricta*, *M. globosa*), `TTGGA` (*Hericium erinaceus*) | Sun et al. 2023 |
| *Cryptococcus neoformans* | `TTA(G)4-6`: `TTA` followed by four to six `G`, a unit of varying length | Edman 1992 |
| budding yeasts (Saccharomycotina): *Candida albicans* (23 bp), *Kluyveromyces lactis* (25 bp) and most others | long and often irregular units of 8 to 25 bp; no single unit | McEachern and Blackburn 1994; Červenák et al. 2021; Lue 2021 |
| *Yarrowia lipolytica* | `TTAGTCAGGG`, a regular unit among the budding yeasts | Lue 2021; Červenák et al. 2021 |
| fission yeasts (*Schizosaccharomyces*) | degenerate `GGTTAC`-type repeats; no single unit | Sepsiova et al. 2016 |
| Microstromatales yeasts | no telomerase; a different tandem array at every chromosome end | Brejová et al. 2026 |

Every exact unit of this table is a candidate in every run. When the data choose one of them, the
candidates table marks it `in_literature`. The telomere repeat database of tidk
([Reference databases](../getting-started/databases.md#telomere-db)) adds the units reported with
`strong` or `weak` support in its 2,990 rows of all kingdoms, 256 distinct units. A unit the data
choose is marked `in_database`, with the number of species per kingdom that carry it. The 126
fungal rows of the database name these units with `strong` or `weak` support:

| Unit, canonical form | Fungal species in the database |
|---|---|
| `AACCCT` (`TTAGGG`) | 20 |
| `AAACCCT` (`TTTAGGG`) | 18 |
| `AACCT` (`TTAGG`) | 9 |
| `AACCTAACT` | 3 |
| `AAAAAACCCT` | 2 |
| `AAAACCCT`, `AACCCCT`, `AACACT`, `AACCC`, `AACCCTAT`, `AAAT` | 1 each |

`TTAGGG` and `TTTAGGG` account for most fungal species in the database. `AAAT` is 4 bp long and
has no G-rich strand, so it fails two of the three filters and is never scored. For a species in
none of these lists, TeloBase (Lyčka et al. 2024) and the Telomerase Database (Podlevsky et al.
2008) list the known telomere repeats.

#### Setting the unit { #setting-the-motif }

`parameters.contig_summary.telomere_motif` sets the telomere unit by hand. A configured unit is
always used, and the data's own choice is still reported next to it.

```yaml
parameters:
  contig_summary:
    telomere_motif: auto             # or a unit, such as TTAGGGTCAACA for the Aspergillus flavus group
```

| `telomere_motif` | Unit used | The three filters | Marked `chosen` in the candidates table |
|---|---|---|---|
| `auto` (default) | the unit chosen from the data, as described above | applied | the data's choice |
| a unit: letters A, C, G and T only, 4 to 30 long, upper or lower case | the unit as given, whatever its score: it is written to `telomere_unit.txt` and used for the telomere columns | skipped | the given unit |

Either strand and any rotation work, because `tidk ends` gives the same result for `TTAGGG`,
`CCCTAA` and `GGGTTA`. The ranking still runs, and a NOTE says when the data favour another unit or
give the configured one no support ([the NOTE](#telomere-note)). Set a unit when the data cannot
decide, for example in a fragmented assembly with few chromosome ends. A changed value of this key
re-runs `telomere_units`, `contig_summary` and MultiQC.

!!! warning "Setting one `telomere_motif` key does not set the other"

    hifiasm has a key of the same name, which is separate:

    | Key | Read by | The unit is written |
    |---|---|---|
    | `parameters.contig_summary.telomere_motif` | the contig summary | on either strand, in any rotation |
    | `parameters.longread.hifiasm.telomere_motif` or `parameters.hybrid.hifiasm.telomere_motif` (`--telo-m`) | the assembler, hifiasm | as it reads at the 5' end of a chromosome: `CCCTAA` for `TTAGGG` |

#### The three telomere files { #telomere-files }

Rule `telomere_units` writes three files to `02.assembly/{sample}/eval/`. They hold the numbers
behind the choice of the unit, and the arrays of that unit at every contig end.

`telomere_candidates.tsv` has one row per unit that was found at the contig ends or formed at
least one array when it was scored. A known unit that formed none is counted in the log only. The
rows are sorted by `score`, then `telomere_ends`, both descending, then `unit`:

| Column | Content |
|---|---|
| `unit` | the unit in canonical form |
| `unit_length` | bp |
| `copies_at_ends` | copies `tidk explore` counted in the terminal windows; `0` for a known unit that was not found there |
| `copies_genome_wide` | copies `tidk explore` counted on the whole draft |
| `telomere_ends` | contig ends with an array pointing the way a telomere points (`present`) |
| `opposite_ends` | contig ends with an array pointing the other way (`wrong_strand`) |
| `score` | `telomere_ends` minus `opposite_ends` |
| `t2t_contigs` | contigs with a `present` array at both ends |
| `in_literature` | the fungal groups or species of the literature table that carry the unit, with the published spelling and the references; `-` when none |
| `in_database` | species per kingdom with `strong` or `weak` support in the tidk database, for example `Fungi 20; Metazoa 420; Viridiplantae 11; (no kingdom) 2` (the few rows without a kingdom are counted as `(no kingdom)`); `-` when none |
| `chosen` | `yes` for the telomere unit of the assembly, `no` otherwise |

`telomere_unit.txt` holds one line, the unit in canonical form or `none`.

`telomere_ends.tsv` is the per-end table of the chosen unit after the same-array rule, the
relabelling to `other_end` in [step 4](#telomere-unit). It has tidk's ten columns and values, and
holds the header line alone when the unit is `none`:

| Column | Content |
|---|---|
| `id` | the contig |
| `sequence_length` | its length in bp |
| `side` | `start` or `end` |
| `status` | `present`, `wrong_strand`, `not_terminal`, `absent` or `other_end` ([step 4](#telomere-unit)) |
| `telomere_start`, `telomere_end` | the array's coordinates on the contig, 0-based; `NA` when absent |
| `telomere_length` | the array's length in bp |
| `offset` | distance from the contig end to the array, bp |
| `strand` | `C-rich` or `G-rich`; a telomere is C-rich at a contig start and G-rich at a contig end |
| `window_limited` | tidk's flag: `true` when the array reaches the edge of the window tidk searched, so it may continue beyond it |

The contig summary reads the rows of this table by their status, and only `present` rows are
telomeres:

| `status` | In the contig summary |
|---|---|
| `present` | a telomere: the copies in `telomere_start` and `telomere_end` are `telomere_length` divided by the unit length, rounded, and the two offset columns are tidk's `offset`. An array more than 1,000 bp inside the contig is also named in the `reason` column, for example `(the start array is 2814 bp from the end)`; the verdict does not change |
| `wrong_strand` | not a telomere; the reason names it (`array in the opposite orientation at the end`) |
| `other_end`, `not_terminal` | nothing |

An array more than 1,000 bp inside is named because a telomere is the very end of a chromosome, so
sequence beyond the array is a misjoin or an interstitial telomeric repeat.

!!! warning "An array further inside than the window is not seen at all"

    A wider `telomere_window` finds it, and the [worked example](#worked-example) shows an array
    that a narrow window would miss.

#### The log line and the NOTE { #telomere-note }

The log `logs/telomere_units_{sample}.log` holds one line that sums up the choice, followed by the
NOTE when there is one, and the rule copies that line to the Snakemake console:

```text
[FunFlux] telomere_units <sample>: telomere repeat unit AACCCT (TTAGGG read on the other strand): 34 contig ends, 15 contigs with a telomere at both ends; known from the literature (filamentous Ascomycota, most Basidiomycota, Mucoromycota; Botrytis cinerea) and from 20 fungal species in the database
```

| Text | Appears when |
|---|---|
| the part in parentheses after the unit | the published spelling differs from the canonical form; for a configured unit it says `set by telomere_motif` |
| the same line, reading `... no telomere repeat identified (best candidate <UNIT>: <n> ends in the telomere orientation, <m> opposite)` | no candidate is eligible |
| a second line, `[FunFlux] NOTE telomere_units <sample>: ...`, also copied to the console | the chosen unit is unknown to both the literature table and the database, so it may be a telomere repeat nobody has reported or a satellite that happens to point one way; or a configured `parameters.contig_summary.telomere_motif` is not the data's first choice, because another unit scores higher or the configured one has no support |

A NOTE is a reason to open the candidates table, which shows the numbers behind either case.

### Genes: the join to the scaffold names { #genes }

`genes` is the one column that cannot be joined to the draft by name, because `funannotate sort`
renames the contigs `scaffold_1`, `scaffold_2`, ... longest first
([funannotate sort](repeat-masking.md#funannotate-sort)). The script therefore joins the GFF3 to
the draft by sequence:

| Draft contig | Matched to its scaffold by |
|---|---|
| unchanged since the draft | the identical sequence: a draft contig and a scaffold with the identical sequence are the same contig |
| polished after the draft (Medaka for ONT samples, Polypolish in `hybrid` mode), so it differs from its scaffold by a base here and there | five 300-bp probes taken from inside the contig, searched in the scaffolds whose length is within 5 % of its own |

The `annotated_as` column gives the result of this join. It links the two sets of contig names,
which none of the other result files does.

## A worked example: *Botrytis cinerea* B05.10 { #worked-example }

The genome of the [validation runs](../about/validation.md), the public *Botrytis cinerea* B05.10
(RefSeq assembly GCF_000143535.2), shows how the table reads on a real assembly. It has 18 chromosomes: 16 core chromosomes of 2.0 to
4.1 Mb and two minichromosomes of 247 and 209 kb. Its mitochondrial genome is 82 kb, and the
reference assembly carries telomere repeats at 27 of its 36 chromosome ends.

A `longread` run of public ONT reads of the strain (about 136x), assembled with hifiasm, gave 33
contigs: 32 in the nuclear genome and the mitochondrion.
This run was made with FunFlux 2.0.0. FunFlux 2.1.0, which trims the ONT adapter from the
reads first, gives the same verdict counts on these reads with one more `low_depth` contig. That
run has a mitochondrial contig of 139,240 bp that is not circular
([Validation on public genomes](../about/validation.md)). The absolute depths of this example,
such as the median of 24.5, are FunFlux 2.0.0's BlobTools estimate, which ran almost five times
too low on these untrimmed reads. The ratios to the median, which the verdicts use, change much
less: measured with samtools coverage on the same reads, as FunFlux 2.1.0 does, the median is
about 117, the mitochondrion is at 6.5 instead of 7.4 times the median and the two
minichromosomes at 0.82, and no verdict changes
([Where the evidence comes from](#where-the-evidence-comes-from)). The telomere unit
chosen from the data is `AACCCT`, that is `TTAGGG`, with arrays of 102 to 126 bp (17 to 21
copies) at 34 of the 66 contig ends. Every one of these arrays points the way a telomere points,
with the C-rich strand at every contig start and the G-rich strand at every contig end, so there
is no opposite end. Of the 33 contigs, 15 have an array at both ends. With the default settings
and a median depth of 24.5, the table calls the 33 contigs:

| Verdict | Contigs | What they are |
|---|---|---|
| `chromosome` | 15 | 13 of the 16 core chromosomes, and the two minichromosomes |
| `chromosome_arm` | 3 | the other 3 core chromosomes, each with one end without a telomere array |
| `rDNA_array` | 13 | pieces of the rDNA array |
| `mitochondrion` | 1 | the mitochondrial genome |
| `low_depth` | 1 | an alternative copy of a region, or noise |

The 18 chromosomes of the reference come out as 18 contigs, 15 of them with telomere arrays at both
ends. The other 15 contigs are rDNA pieces, the mitochondrion and one low-depth contig. A few rows,
with the columns that decided:

| `contig` | `length` | `depth_x_median` | `telomere_start` | `telomere_end` | `rdna_units` | `busco_markers` | `genes` | `annotated_as` | `verdict` |
|---|---|---|---|---|---|---|---|---|---|
| `ptg000005l` | 4,109,401 | 1.00 | 18 | 17 | 0 | 287 | 1,127 | `scaffold_1` | `chromosome` |
| `ptg000007l` | 3,044,302 | 1.00 | 18 | 0 | 0 | 193 | 869 | `scaffold_4` | `chromosome_arm` |
| `ptg000014l` | 2,512,653 | 0.99 | 18 | 0 | 6 | 184 | 701 | `scaffold_9` | `chromosome_arm` |
| `ptg000022l` | 208,977 | 0.79 | 18 | 17 | 0 | 0 | 21 | `scaffold_18` | `chromosome` |
| `ptg000004l` | 282,192 | 0.77 | 17 | 17 | 0 | 0 | 22 | `scaffold_17` | `chromosome` |
| `ptg000030l` | 175,720 | 1.15 | 18 | 0 | 20 | 0 | 1 | `scaffold_19` | `rDNA_array` |
| `ptg000011l` | 58,954 | 0.16 | 0 | 0 | 6 | 0 | 0 | `scaffold_25` | `rDNA_array` |
| `ptg000015l` | 61,330 | 0.27 | 0 | 0 | 6 | 0 | `NA` | `dropped by funannotate clean` | `rDNA_array` |
| `ptg000028c` | 82,351 | 7.37 | 0 | 0 | 0 | 0 | `NA` | `not annotated (mitochondrion)` | `mitochondrion` |
| `ptg000029l` | 42,998 | 0.07 | 0 | 0 | 0 | 1 | 5 | `scaffold_27` | `low_depth` |

What the rows show:

- **The core chromosomes.** The 16 contigs of 2.0 to 4.1 Mb are the 16 core chromosomes. Of
  these, 13 are `chromosome` and 3 are `chromosome_arm`, each with one end without a telomere
  array. `ptg000014l` carries 6 rDNA units and 701 genes, so verdict rule 4 does not apply and it
  stays a chromosome arm. The start array of `ptg000007l` begins 2,814 bp inside the contig, and
  its reason says so:
  `telomere repeats at the start only (18 copies) (the start array is 2814 bp from the end)`.
  Three chromosomes have one array more than 1,000 bp inside the contig end, named in their
  reasons the same way: 1,134 bp in `ptg000002l`, 3,660 bp in `ptg000009l` and 12,300 bp in
  `ptg000017l`.
- **The minichromosomes.** The two contigs of 200 to 300 kb are the two minichromosomes,
  `chromosome` whatever their size. `ptg000022l` (209 kb) has its arrays right at the contig ends
  and 21 genes. `ptg000004l` (282 kb) has arrays of 17 copies that lie 19,746 and 15,234 bp inside
  the contig ends. The contig therefore carries about 35 kb beyond its telomeres, and the
  reference chromosome is 247 kb. Its reason names both offsets:
  `telomere repeats at both ends (17 and 17 copies) (the start array is 19746 bp from the end; the end array is 15234 bp from the end); depth 0.77x the median; 22 genes, 0 BUSCO markers`.
  With a window of 3,000 bp, neither array would be seen and the contig would read
  `nuclear_contig`, which is why the default window is 25 kb. The depths of the two
  minichromosomes, 0.79 and 0.77 times the median, are within the normal range.
- **The rDNA pieces.** `ptg000030l` is a 176-kb rDNA array with a telomere at its start (18
  copies), 20 rDNA units and one predicted gene, a gene call inside an rDNA unit. With fewer genes
  than rDNA units, it is an `rDNA_array` by verdict rule 4, with the telomere named in the reason.
  The 13 `rDNA_array` contigs are 36 to 176 kb of rDNA units at 0.16 to 1.15 times the median:
  the array is assembled in pieces and partly collapsed. Two of them, `ptg000015l` and
  `ptg000024l`, were dropped by `funannotate clean` as duplicates, so their `genes` cell is `NA`
  and `annotated_as` says so. The `genus` cell of the rDNA pieces names another genus than
  *Botrytis*, because a BLAST genus label is only as good as the closest sequence in the database
  ([Decontamination](decontamination.md#why-there-is-no-auto-mode)). The rDNA unit holds lone
  copies of `TTAGGG` and `CCCTAA`, but a telomere call needs an array of at least 5 copies, so
  none of the other rDNA pieces gets a telomere.
- **The mitochondrion.** `ptg000028c` is confirmed by the mitochondrion screen. It is circular, at
  7.4 times the median depth, and not annotated.
- **The low-depth contig.** `ptg000029l`, 43 kb at 0.07 times the median with 5 genes, is
  `low_depth`, an alternative copy of a region or noise.

The candidates table shows how the unit was chosen. `telomere_candidates.tsv` has 20 rows, and
these are the three units with the most copies at the contig ends:

| Unit | Copies at the contig ends | Ends in the telomere orientation | Opposite ends | Score |
|---|---|---|---|---|
| `AACCCT` | 657 | 34 | 0 | 34 |
| `ACCTCC` | 164 | 0 | 3 | -3 |
| `AGATGG` | 152 | 0 | 4 | -4 |

Only the telomere unit has arrays that point the telomere way. The two microsatellites that follow
it have three and four arrays at contig ends, all pointing the wrong way for a telomere. `AACCCT`
has 657 copies at the contig ends and 661 in the whole assembly, so nearly every copy sits at an
end, while `AGATGG` has 405 in the whole assembly. Seven units have a single end in the telomere
orientation and are not eligible. The table marks `AACCCT` known from the literature and from 20
fungal species in the database. The log of `telomere_units` ends with the summary line quoted
[above](#telomere-note) and has no NOTE.

A short-read assembly of the same strain chooses the same unit from thinner evidence. It was made
with SPAdes from public Illumina reads and has 413 contigs, 4 of them mitochondrial:

| Unit | Copies at the contig ends | Ends in the telomere orientation | Opposite ends | Score | Outcome |
|---|---|---|---|---|---|
| `AACCCT` | 227 | 13 | 1 | 12 | chosen; no contig with an array at both ends |
| `AAGAG` | 33 | 13 | 5 | 8 | eligible, ranks second |
| `AGAGGG` | 115 | 12 | 6 | 6 | not eligible: 12 is not more than twice 6 |
| `AAAGAG` | 179 | 5 | 5 | 0 | arrays point either way |
| `AAGGAG` | 171 | 6 | 7 | -1 | arrays point either way |

Here the copy number cannot choose the unit, but the orientation still can. The microsatellites
`AAAGAG` and `AAGGAG` have nearly as many copies at the contig ends as `AACCCT`, but their arrays
point either way, so they score 0 and -1. `tidk ends` first reports 7 opposite ends for `AACCCT`.
Six of them are the same array seen from the other end of a short contig of 1.4 to 12 kb, which
the same-array rule relabels `other_end`. The seventh is a 30-bp run 12.4 kb inside a long contig.

The gene count is FunFlux's own prediction. Where a published annotation exists, as here, its gene
count is an independent check of the run as a whole. The table itself reads no reference.

## Settings { #settings }

The contig summary has five settings in `parameters.contig_summary`, and the address of the
telomere database in `links`:

```yaml
parameters:
  contig_summary:
    telomere_motif: auto
    telomere_db_sha256: 1d4f5c8c2cce6a76aaa8de01602317a8e0b8d0401711b96f3c7a826c385da699
    telomere_window: 25000
    min_telomere_copies: 5
    min_rdna_units: 3
links:
  telomere_db_link: https://raw.githubusercontent.com/tolkit/a-telomeric-repeat-database/73c635c7b80e5ab369b4c7ec67d23db32dcaa1b8/data/telo_database_v2.tsv
```

| Key | Default | Allowed values | What it does |
|---|---|---|---|
| `telomere_motif` | `auto` | `auto` (any case), or letters A, C, G and T, 4 to 30 long; lower case is upper-cased | `auto`: the unit is chosen from the data. A unit: used as given, the data's choice reported next to it ([Setting the unit](#setting-the-motif)). |
| `telomere_db_sha256` | the checksum of the database file at the template address | 64 hexadecimal characters | The SHA-256 the downloaded database must have. Change it together with `links.telomere_db_link`. |
| `telomere_window` | 25000 | an integer from 500 to 100000 | The terminal window: the bases at each contig end in which the candidate units are collected, and how far from the end `tidk ends` accepts an array. The array can sit some kb inside a contig end when the assembler ran past the telomere, so the window must cover that; the offset columns show where each array lies. |
| `min_telomere_copies` | 5 | an integer ≥ 2 | Copies of the unit an array needs: `tidk ends` requires this many unit lengths (30 bp for a 6-bp unit), and `tidk explore` lists a unit only when it forms a tandem run of this many copies. |
| `min_rdna_units` | 3 | an integer ≥ 1 | rDNA units (ITS windows) from which a contig can be an `rDNA_array`. |

The five values are checked before any job starts, like the other blocks
([Configuration](../reference/configuration.md#contig_summary)). A value outside the range stops
the run with a `[FunFlux] parameters.contig_summary.<key>: ...` message. Every dry run prints the
settings in one line:
`Contig summary: telomere motif auto (the repeat unit is chosen from the data), a telomere from 5 tandem copies within 25000 bp of a contig end, an rDNA array from 3 rDNA unit(s).`
A changed value re-runs the steps that use it ([Re-runs](#re-runs)).

## In the MultiQC report { #multiqc }

For the report, `02.assembly/{sample}/eval/contig_summary_mqc.tsv` holds the same table prepared
for MultiQC, with the rows named `{sample}: {contig}` and a description that names the telomere
unit of the sample. MultiQC shows the tables of all samples of the run as one table in the section
**Contig summary** of `06.report/multiqc_report.html`. It writes the numbers to
`06.report/multiqc_data/multiqc_contig_summary.txt`. The table has all 21 columns, `telomere_unit`
among them, so samples of different species keep their units apart. The column settings of the
table show or hide columns.

!!! note "Above 500 rows, MultiQC shows a violin plot instead of the table"

    MultiQC shows a table with more than 500 rows (its `max_table_rows` setting) as a violin plot
    instead, so for a batch of fragmented assemblies read `contig_summary.tsv` itself.

## Re-runs

Each change in the left column re-runs the contig summary steps on the right:

| Change | Re-runs |
|---|---|
| a new draft | `telomere_units`, `contig_depth` (read-based modes) and `contig_summary` |
| a new database file | `telomere_units` |
| `parameters.contig_summary.telomere_motif`, `telomere_window` or `min_telomere_copies` | `telomere_units` and `contig_summary` |
| `links.telomere_db_link` or `telomere_db_sha256` | `telomere_db`, then `telomere_units` of every sample and the steps after it |
| a new telomere unit or ends table | `contig_summary` |
| `min_rdna_units` | `contig_summary` |
| `links.univec_link` or `univec_sha256` | `univec_illumina`, then `blast_adapters`, `screen_adapters` and `contig_summary` of every sample |
| a new adapter screen result | `contig_summary` |
| a new screen result, a new `busco.lineage` (the BUSCO markers are counted per contig), new ITS windows or a new annotation | `contig_summary` |

A re-run of these steps costs seconds to a minute and never touches the annotation. MultiQC
follows each of them.

## Limits

Some genomes, assemblies and settings leave a piece of evidence missing or unreliable. The table
lists these situations and what the contig summary does in each:

| Situation | What the table does |
|---|---|
| telomeres that vary from copy to copy: *Saccharomyces*, *Schizosaccharomyces*, *Cryptococcus* | only exact units are found, and no tandem run of one unit exists: `tidk explore` lists no unit, `tidk ends` finds no array, and the unit comes out as `none` even when the irregular repeat is written into `parameters.contig_summary.telomere_motif` |
| telomere repeats that differ between chromosome ends: the Microstromatales yeasts | no single unit covers more than a few of them |
| a unit that no paper and no database names | it is still found when it forms an array of at least `min_telomere_copies` copies within the window at the contig ends, because every tandem unit of the contig ends is a candidate; the NOTE then says that the unit is unknown |
| a genuine telomere unit that fails one of the three filters | it would be missed; none is known, and `parameters.contig_summary.telomere_motif` bypasses the filters |
| an assembly with few chromosome ends, such as a fragmented short-read assembly | few votes: a unit needs at least 2 telomere ends and more than twice as many telomere ends as opposite ends; when it falls short no unit is chosen, so set `parameters.contig_summary.telomere_motif` when the repeat is known |
| an array further inside a contig than `telomere_window` | it is not seen; the [worked example](#worked-example) shows the trade-off of a wider window |
| an array of the chosen unit pointing the other way at a contig end (`wrong_strand`) | not a telomere call; the reason names it, and the verdict rests on the rest |
| a contig that ends at a centromere | centromeres are not detected: it is a `chromosome_arm` like any other |
| `contigs` mode, where depth is unknown (it is known only in the read-based modes) | verdict rules 7 and 8 never apply; a mitochondrion that the screen did not confirm, or a contig at twice the depth of the genome, is judged on the other evidence alone: a gene-less contig at twice the depth is then a `nuclear_fragment`, not a `repeat` |
| `parameters.mitochondria.mode: off`, read-based modes | the screen confirms nothing, and the mitochondrion is judged on the rest of its evidence: its depth, several times the median, keeps it out of every verdict but `repeat` (when no gene was predicted on it) or `unclassified` |
| `parameters.mitochondria.mode: off`, `contigs` mode | without a depth, the mitochondrion is a `nuclear_contig` when genes were predicted on it |
| `parameters.mitochondria.mode: report` | the confirmed contig stays in the nuclear genome and is annotated; the table still calls it `mitochondrion` (verdict rule 1, from the screen's class), and its `annotated_as` cell names the `scaffold_N` under which funannotate annotated it, with its `genes` cell filled ([Mitochondria](mitochondria.md#modes)) |
| a fungus whose 5.8S matches none of ITSx's profiles | `rdna_units` counts the 5.8S genes that nhmmer found with ITSx's profiles ([Why the genome is cut into windows](its-taxonomy.md#why-the-genome-is-cut-into-windows)), so it is `0` on every contig, and the rDNA pieces are judged on depth and genes alone: `low_depth`, `repeat`, `nuclear_fragment`, `nuclear_contig` or `unclassified` |
| the gene count | it is the run's own gene prediction, counted per scaffold: evidence that a contig is nuclear sequence, not a measure of annotation quality |
| a polished contig (Medaka in `longread` and `hybrid` mode, Polypolish in `hybrid` mode) | its join to its `scaffold_N` rests on five 300-bp probes and a length within 5 %; `dropped by funannotate clean` on a contig that is in `contigs_sort.fasta` means the join failed, and the `genes` cell of that contig is then `NA` |

!!! note "Where evidence is missing or unreliable, the verdict rests on the evidence that remains"

    For the fungi of the first two rows, the table still gives the depth, rDNA, BUSCO and gene
    evidence, and the verdicts `chromosome` and `chromosome_arm` do not occur. The `annotated_as`
    column shows whether each contig found its scaffold.

## Output files

Every file is kept, and `contig_summary.tsv` is the one to read:

| Path | Content |
|---|---|
| `00.resources/telomeres/telo_database_v2.tsv` | the telomere repeat database of tidk, downloaded once per output folder |
| `00.resources/telomeres/PROVENANCE.txt` | source URL, download date, checksum and row count of the database |
| `02.assembly/{sample}/eval/telomere_candidates.tsv` | every candidate telomere unit with its counts ([The three telomere files](#telomere-files)) |
| `02.assembly/{sample}/eval/telomere_unit.txt` | the telomere unit of the assembly, or `none` |
| `02.assembly/{sample}/eval/telomere_ends.tsv` | the arrays of that unit at every contig end |
| `02.assembly/{sample}/contaminants/{sample}_depth.tsv` | the read depth of every contig, from samtools coverage; read-based modes only |
| `02.assembly/{sample}/eval/contig_summary.tsv` | the table, one row per contig of the draft |
| `02.assembly/{sample}/eval/contig_summary_mqc.tsv` | the same table for MultiQC |
| `logs/telomere_db.log` | the database download |
| `logs/telomere_units_{sample}.log` | the count line of the script, the one-line summary and the NOTE when there is one; tidk's own messages only when tidk fails |
| `logs/contig_depth_{sample}.log` | samtools' own messages; empty when the step succeeds |
| `logs/contig_summary_{sample}.log` | the one-line summary |

Each step writes its files only when every input was read without error. The job stops on a
missing file, a malformed database, a tidk failure, a table without an expected column, a depth
that is not a number or a FASTA with two records of the same name. It reports the cause as one line
in the log, `[FunFlux] telomere_units: error: ...` or `contig_summary.py: error: ...`, which for a
tidk failure includes tidk's own message.

The full output tree is on [Output files](../reference/output.md).
