# BLAST hit depth

This page explains how many BLAST hits FunFlux keeps for each contig, which e-value the search
uses, and why. Read it before you change `parameters.blast_max_target_seqs` or
`parameters.blast_max_hsps`, when you update from FunFlux 2.0.0, or when you compare genus calls
with those of FunFlux 2.0.0 or 1.1.x.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **the screen** | the contamination and mitochondrion screen of the draft assembly: BLAST, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **draft assembly** | the assembly the screen examines, before the screen removes any contig: `assembly.fasta` in `longread` and `hybrid` mode, and `contigs_filt.fasta`, the contigs after FunFlux's contig filter, in `illumina` mode, for the `hybrid` guide and in `contigs` mode; not yet polished |
    | **guide assembly** | in `hybrid` mode, the SPAdes assembly of the Illumina reads; it gets a screen of its own and is not delivered |
    | **contig selector** | rule `select_contigs`, the last step of the screen: it keeps or removes each contig by genus, as `parameters.decontamination` says; the default, `off`, keeps every contig |
    | **core_nt** | the NCBI core nucleotide BLAST database, which the screen searches every draft contig against; you download it, with its taxonomy files, into `directories.blast_db` |
    | **subject** | the BLAST word for a database sequence a contig matched, the contig being the query; BLAST keeps up to 10 per contig (`blast_max_target_seqs`; BLAST's own option, `-max_target_seqs`, calls them target sequences) |
    | **HSP** | high-scoring segment pair, one local alignment between a contig and a subject; BLAST keeps up to `blast_max_hsps` (5) per subject |
    | **hit** | one line of the BLAST table: one HSP between a contig and a subject |
    | **bit score** | BLAST's score for one HSP: the longer and closer the match, the higher the score |
    | **e-value** | BLAST's estimate of how often a match at least this good would turn up by chance in a search of the whole database; the search drops matches above its cut-off, so 1e-10 is stricter than 1e-5 |
    | **bestsum** | the rule BlobTools uses to give each contig a taxon: it adds up the bit scores of all the contig's hits per taxon and takes the taxon with the highest sum |
    | **direct hit** | a hit to a subject whose title names a mitochondrial genome, such as "mitochondrion, complete genome"; the mitochondrion screen counts only these ([Mitochondria](../analysis/mitochondria.md) lists the wordings) |
    | **FunFluxL, FunFluxL+** | two of the v1 workflows that FunFlux v2 replaces: FunFluxL for long reads and the unpublished FunFluxL+ for hybrid data |
    | **clean circle** | a contig the assembler recorded as circular, not as a repeat and not as an alternative copy of another contig |
    | **genus filter** | the setting `parameters.decontamination.mode: include` or `exclude`, under which the contig selector keeps or removes contigs by the genus BlobTools gave them |
    | **MGISEQ** | a short-read sequencer from MGI Tech that produces paired reads like Illumina's |

The contamination and mitochondrion screen runs one BLAST search per draft assembly. Rule
`blast_contigs` runs it with megablast against core_nt and an e-value of 1e-10. In `hybrid` mode,
rule `blast_guide_contigs` does the same for the guide assembly, with the same
`blast_max_target_seqs` and `blast_max_hsps`. Two programs read the result table, BlobTools and
the mitochondrion screen:

```mermaid
---
config:
  themeVariables:
    fontSize: 14px
  flowchart:
    nodeSpacing: 35
    rankSpacing: 35
---
flowchart TD
    C([draft contigs]) --> B[megablast against core_nt<br/>blast_contigs]
    B --> T([BLAST table<br/>by default at most 10 subjects<br/>per contig, 5 HSPs per subject])
    T --> BT[BlobTools<br/>assigns each contig a genus]
    T --> MS[mitochondrion screen<br/>looks for mitochondrial contigs]
```

Two BLAST options decide how much of the search is kept. By default FunFlux sets them to 10 and 5
in every mode.

| Key | Default | blastn option | blastn help text (BLAST+ 2.16.0) |
|---|---|---|---|
| `parameters.blast_max_target_seqs` | 10 | `-max_target_seqs` | "Maximum number of aligned sequences to keep (value of 5 or more is recommended)" |
| `parameters.blast_max_hsps` | 5 | `-max_hsps` | "Set maximum number of HSPs per subject sequence to save for each query" |

| Workflow or source | `max_target_seqs` / `max_hsps` | e-value |
|---|---|---|
| FunFlux 1.1.x (Illumina) | 10 / 1 | 1e-5 |
| FunFluxL and FunFluxL+ | 50 / 5 | 1e-5 |
| FunFlux 2.0.0, every mode | 50 / 5 | 1e-5 |
| FunFlux 2.1.0, every mode | 10 / 5 | 1e-10 |
| BlobTools manual, for its "most accurate results" | 10 / 1 | 1e-25 |
| BlobToolKit tutorial "Adding hits" | 10 / 1 | 1e-25 |
| Sanger Tree of Life BlobToolKit pipeline | 10 / 1 | 1e-10 |

FunFlux keeps 10 subjects per contig, as the three BlobTools sources do, and searches with the
e-value of the Sanger Tree of Life pipeline, 1e-10. Where the sources keep one HSP per subject,
FunFlux keeps 5, because its mitochondrion screen reads the same table and needs more than one HSP
per subject ([below](#five-hsps-what-the-mitochondrion-screen-needs)). The sources do not explain
their values.

FunFlux 2.0.0 kept 50 subjects with an e-value of 1e-5, the values of FunFluxL. With 50 subjects, a
genus with many database entries of one gene, such as the same gene sequenced from many isolates,
can outvote the genus of a contig's closest hits
([Ten subjects](#ten-subjects-what-the-genus-calls-need)). The e-value made almost no difference
in the [comparison](#funflux-200-and-210-compared): a close match over hundreds of bases has an
e-value far below 1e-5, 1e-10 and 1e-25 alike, so the cut-off drops only weak or short matches.

## Five HSPs: what the mitochondrion screen needs

The mitochondrion screen, rule `screen_mitochondria`, reads the BLAST column `qcovs` ("Query
Coverage Per Subject"). `qcovs` is the share of the contig covered by all the HSPs that BLAST saved
for one subject. With `-max_hsps 1` only the best HSP is saved, and `qcovs` equals `qcovhsp` ("Query
Coverage Per HSP"). Rearrangements, indels or a different start position on the circle can make a
mitochondrial genome align to its closest reference in several blocks. With one HSP per subject,
such a genome then looks partly covered.

For each contig, the mitochondrion screen takes the highest `qcovs` among its direct hits and
confirms the contig as follows:

| Mode | Evidence from the assembler | A contig is confirmed when its highest `qcovs` |
|---|---|---|
| `contigs`, `illumina`, the `hybrid` guide | none | reaches `mitochondria.min_direct_blast_qcov` (default 35). A contig with a weaker direct hit becomes a `mitochondrial_candidate`; in `contigs` and `illumina` mode it stays in `contigs_final.fasta` and goes through the nuclear gene prediction |
| `longread`, the `hybrid` long-read assembly | the depth, or a clean circle | reaches `mitochondria.very_strong_qcov` (70) on its own, or `min_direct_blast_qcov` (35) when the contig also has a high depth (at least `min_assembler_coverage_ratio`, 3.0, times the median) or is a clean circle |

!!! warning "One HSP per subject can leave a complete mitochondrial genome unconfirmed"

    With one HSP per subject, only the best block of a multi-block match counts. A complete
    mitochondrial genome that aligns to its nearest core_nt reference in many pieces can then fall
    below `min_direct_blast_qcov`. The mitochondrion screen is described on
    [Mitochondria](../analysis/mitochondria.md).

The default of 5 HSPs has not been compared with other values on a public genome.

## Ten subjects: what the genus calls need

BLAST applies `-max_target_seqs` during the search, not as a filter on the finished, ranked list
of hits. A small value can therefore drop subjects that would have ranked higher (Shah et al.
2019). NCBI's reply confirmed the mechanism, fixed a related bug in BLAST+ 2.8.1 and made BLAST+
warn when fewer than five subjects are requested (Madden et al. 2019). BLAST+ 2.16.0 prints:

```
Warning: [blastn] Examining 5 or more matches is recommended
```

FunFlux's 10 meets that advice. More subjects are not better for the genus calls, because bestsum
adds up the bit scores of every hit per genus. A genus with many database entries of one gene, such
as the same gene sequenced from many isolates of one species, then gains the scores of every entry,
while the genus with the closest hits may have only a few entries. With 10 subjects, at most 10
entries of one gene can enter the sum. That weakens the effect; it does not remove it.

!!! example "One contig of *Cordyceps cateniannulata* GXU-8616"

    The HiFi assembly of the [validation runs](../about/validation.md#pacbio-hifi-reads-cordyceps-cateniannulata-gxu-8616)
    has a contig of 4.9 Mb. With 50 subjects, 41 of them were *Beauveria* (38 of them *B.
    bassiana*), all the same gene (carbamoyl-phosphate synthase) sequenced from different
    isolates, with bit scores of 7,062 to 7,306 each. Six were *Cordyceps*, including the best
    single hit (bit score 10,567). The sums were 296,140 for *Beauveria* and 53,026 for
    *Cordyceps*, and BlobTools called the contig *Beauveria*. With 10 subjects (five *Cordyceps*,
    three *Beauveria*, two *Akanthomyces*), the sums were 45,921 for *Cordyceps* and 21,879 for
    *Beauveria*, and the call was *Cordyceps*, the strain's genus.

## FunFlux 2.0.0 and 2.1.0 compared

The table compares the setting of FunFlux 2.0.0 (50 subjects, e-value 1e-5) with that of FunFlux
2.1.0 (10 subjects, e-value 1e-10). Both keep 5 HSPs per subject. The four draft assemblies were
made from the public reads of the [validation runs](../about/validation.md), and the screen ran on
them in `contigs` mode, so BLAST was the only evidence of the mitochondrion screen. *B. cinerea* is
well represented in core_nt. *C. cateniannulata* itself has few sequences there, so its contigs are
called from its relatives: other *Cordyceps* species, but also other genera such as *Beauveria*.
"Genus changed" is the share of the assembly's bases in contigs whose genus differs between the two
settings.

| Draft assembly | Contigs | Genus changed | The strain's genus | Bases called as that genus, 2.0.0 (50 subjects) | 2.1.0 (10 subjects) |
|---|---|---|---|---|---|
| *B. cinerea* B05.10, Illumina reads, SPAdes | 413 | 0.02 % | *Botrytis* | 99.7 % | 99.6 % |
| *B. cinerea* B05.10, ONT reads, hifiasm | 34 | 0.17 % | *Botrytis* | 97.6 % | 97.7 % |
| *C. cateniannulata* GXU-8616, MGISEQ reads, SPAdes | 294 | 16.9 % | *Cordyceps* | 56.4 % | 72.8 % |
| *C. cateniannulata* GXU-8616, HiFi reads, hifiasm | 19 | 27.9 % | *Cordyceps* | 46.4 % | 74.3 % |

- **On *B. cinerea*** the two settings agree: a few small contigs change genus.
- **On *C. cateniannulata*** 10 subjects change the genus of 17 % (MGISEQ) and 28 % (HiFi) of the
  assembly, almost all of it from *Beauveria*, *Fusarium* or *Colletotrichum* to *Cordyceps*. Some
  calls are still wrong. Two HiFi contigs (4.8 and 3.1 Mb) and twelve MGISEQ contigs (four of 0.8
  to 2.8 Mb and eight small ones) stay *Beauveria*, about a quarter of each assembly. Each of the
  six large ones has its best single hit in *Cordyceps*, one of them in *C. cateniannulata*
  GXU-8616 itself. But four to nine of its 10 subjects are the same gene sequenced from many
  *Beauveria* isolates.
- **The e-value alone**, 1e-10 instead of 1e-5 with 50 subjects, changed the genus of at most
  0.03 % of the bases. Five short contigs (0.6 to 4.2 kb) lost their only hits and became
  `no-hit`. With the default, `mode: off`, they stay. With a genus filter they are now removed:
  `include` removes every `no-hit` contig, and `exclude` does so with `discard_no_hit: true`.
  FunFlux 2.0.0 had given four of them the strain's genus.
- **The mitochondrion screen** confirmed the same mitochondrial contigs, with the same coverage,
  under both settings. Both strains' own mitochondrial genomes are in core_nt, and their
  mitochondrial contigs match them best. Three contigs of 32 kb to 2 Mb were
  `mitochondrial_candidate` with 50 subjects and not with 10, because the subject behind the call,
  which covered 1 % of the contig or less, was no longer among the 10 kept. One of these subjects
  was a mitochondrial genome; the other two were nuclear genes whose names contain "mitochondrial
  DNA", which the screen counts as direct hits. A candidate stays in the nuclear genome, so no
  genome file changed; only those three rows of `mitochondrial_screen.tsv` did.
- **Without the strain's own species in core_nt**, as for a new species, the three mitochondrial
  contigs were searched again with that species left out. With it, their best direct hit covered
  99 to 100 % of them. Without it, all three were still confirmed with either setting, but 10
  subjects left less margin. The best direct hit of the GXU-8616 mitochondrion covered 46 % with
  50 subjects and 35 % with 10 in the HiFi draft, where the contig holds two copies of the
  mitochondrial genome end to end, and 70 and 57 % in the MGISEQ draft. 35 % is the threshold when
  BLAST is the only evidence; in `longread` mode the HiFi contig would also have needed its circle,
  which hifiasm recorded. The 10-subject search missed the mitochondrial genomes that covered most,
  although they scored higher than the ones it kept, because BLAST applies the limit during the
  search. The B05.10 mitochondrion reached 89 % with either setting. If no mitochondrion is
  confirmed, see [Mitochondria](../analysis/mitochondria.md#when-no-mitochondrion-is-confirmed).
- **Run time** was the same with either setting: 44 to 63 minutes per draft on 24 threads, with
  both searches running at the same time on a shared server. The BLAST table is three to five
  times smaller with 10 subjects.

The genus calls matter only with a genus filter, `decontamination.mode` `include` or `exclude`.
With `off`, the default, every contig is kept whatever its genus, so the genus calls change no
delivered genome file. See [Decontamination](../analysis/decontamination.md).

These are two genomes, so the numbers show the direction of the effect, not its size on another
fungus.

## Updating from FunFlux 2.0.0

- **Your configuration file.** The FunFlux 2.0.0 template, `config/config.yaml`, sets
  `blast_max_target_seqs: 50`, and updating FunFlux does not change the copy you made from it. Set
  it to 10, or delete the line, to search as FunFlux 2.1.0 does.
- **An output folder made with 2.0.0.** Started from the same launch folder as before, FunFlux
  2.1.0 runs the BLAST search there again, even if you keep `blast_max_target_seqs: 50`, because
  the BLAST command itself changed (its e-value) and Snakemake compares it with the one recorded
  in `.snakemake/` of the launch folder. Every later step follows for each sample with the screen,
  including the whole annotation: hours per sample. With the default, `mode: off`, every contig is
  still kept, but the genes are predicted again, and two annotations of the same genome differ at
  about 1 % of the gene models ([Run-to-run differences](run-to-run-differences.md)). Started from
  a new clone or another launch folder, Snakemake has no record of the old command and keeps the
  2.0.0 BLAST table; add `--forcerun blast_contigs` (in `hybrid` also `blast_guide_contigs`) to
  search again. To keep finished 2.0.0 results as they are, leave that folder to FunFlux 2.0.0.
- **A genus filter chosen from a 2.0.0 run** may no longer fit, because some contigs change genus.
  Read the new `contig_taxonomy_decisions.tsv` before you filter again
  ([Decontamination](../analysis/decontamination.md)).

## Notices and re-runs

Both values must be whole numbers of at least 1. With the screen on, two notices can print when
the workflow starts:

```
[FunFlux] NOTE: parameters.blast_max_hsps is 1. Below 5, BLAST reports the coverage of a multi-block mitochondrial match as that of its best block only, which weakens the mitochondrion screen (a complete mitogenome can fall under min_direct_blast_qcov).
[FunFlux] NOTE: parameters.blast_max_target_seqs is 3; BLAST recommends 5 or more.
```

The first notice, on `blast_max_hsps`, prints when that value is below 5 and `mitochondria.mode`
is not `off`. The second, on `blast_max_target_seqs`, prints when that value is below 5.

Changing either value runs the BLAST search again (in `hybrid` both searches) and every step after
it, including the whole annotation: hours per sample.

## Reproducing earlier versions

Apart from the mitochondrion case described on
[Mitochondria](../analysis/mitochondria.md#when-no-mitochondrion-is-confirmed), set these values
only to compare with results of an earlier version. The e-value is fixed at 1e-10 and cannot be
set, while the earlier versions searched with 1e-5; in the comparison above, the e-value changed the
genus of at most 0.03 % of the bases.

To come close to FunFlux 2.0.0, FunFluxL or FunFluxL+:

```yaml
parameters:
  blast_max_target_seqs: 50
  blast_max_hsps: 5
```

To come close to FunFlux 1.1.x:

```yaml
parameters:
  blast_max_target_seqs: 10
  blast_max_hsps: 1
```

FunFlux 1.1.x took 10 subjects and 1 HSP, but not the e-value, from the BlobToolKit tutorial
"Adding hits", which prepares a BLAST table for the taxonomy only. These values print the
`blast_max_hsps` notice, because the mitochondrion screen is weaker with one HSP. FunFlux v2 has
not been tested with them.

## References

- Shah N., Nute M.G., Warnow T., Pop M. (2019). Misunderstood parameter of NCBI BLAST impacts the
  correctness of bioinformatics workflows. *Bioinformatics* 35(9):1613-1614.
  <https://doi.org/10.1093/bioinformatics/bty833>
- Madden T.L., Busby B., Ye J. (2019). Reply to the paper: Misunderstood parameters of NCBI BLAST
  impacts the correctness of bioinformatics workflows. *Bioinformatics* 35(15):2699-2700.
  <https://doi.org/10.1093/bioinformatics/bty1026>
- BlobToolKit, "Adding hits":
  <https://blobtoolkit.genomehubs.org/blobtools2/blobtools2-tutorials/adding-data-to-a-dataset/adding-hits/>
- BlobTools manual, "Taxonomy file": <https://blobtools.readme.io/docs/taxonomy-file>
- Sanger Tree of Life, BlobToolKit pipeline, BLAST settings in `conf/modules.config`:
  <https://github.com/sanger-tol/blobtoolkit>
