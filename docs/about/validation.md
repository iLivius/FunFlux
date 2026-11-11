# Validation on public genomes

This page shows what FunFlux produced on two public genomes, and how to repeat the runs. Read it
to see, on real data, what each read-based mode delivers and how much the gene-finder settings
change the annotation.

*Botrytis cinerea* strain B05.10 was run from ONT and Illumina reads in the three read-based modes,
`illumina`, `longread` and `hybrid`. Five comparison runs then annotated the genome of its FunFlux
2.0.0 `hybrid` run again in `contigs` mode. *Cordyceps cateniannulata* strain GXU-8616 was run from PacBio HiFi
reads, alone and with MGISEQ short reads
([PacBio HiFi reads](#pacbio-hifi-reads-cordyceps-cateniannulata-gxu-8616)). Flye has run in
`longread` mode on both genomes but not in `hybrid` mode, where every run used hifiasm, the
default assembler. The results are measurements on these two genomes, so read the limitations
below before the numbers.

The four runs with ONT reads were made with FunFlux 2.1.0, whose main change, ONT adapter
trimming, applies only to ONT reads
([Adapter trimming in the runs with ONT reads](#adapter-trimming)). In those runs it changed the
contig ends and the mitochondrion much more than the chromosomes or the genes. The other nine
runs, with Illumina, MGISEQ or HiFi reads or in `contigs` mode, were made with FunFlux 2.0.0 and
not repeated. The five `contigs`-mode runs therefore annotated the genome of the 2.0.0 `hybrid`
run, not that of the 2.1.0 `hybrid` run in the assembly tables. Each table, or the text just
before it, says which version its numbers come from.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **Tiberius** | a deep-learning gene finder with a `fungi` model its authors trained on many fungal genomes; on by default, it adds its gene models to funannotate's evidence |
    | **annotation parameter table** | the tab-separated file named by `files.annotation_params`, one row per sample: name, species, a protein FASTA of related fungi and, optionally, an AUGUSTUS model |
    | **protein evidence** | proteins aligned to the genome to support the gene models: the FASTA in column 3 of the annotation table, plus UniProt Swiss-Prot, always added |
    | **stored AUGUSTUS model** | the AUGUSTUS parameters trained earlier for one species, kept in `trained_species/` of the funannotate database; named in column 4, they are used as they are |
    | **telomere array** | tandem copies of the telomere repeat unit at a contig end; by default 5 or more copies within 25,000 bp of the end make a telomere |
    | **contig summary** | the table `02.assembly/{sample}/eval/contig_summary.tsv`: one row per contig of the draft, with its evidence (telomeres, depth, rDNA units, genes) and a verdict such as `chromosome` or `rDNA_array` |
    | **adapter trimming** | removing the ONT sequencing adapter, and the barcode after it, from the reads before assembly; new in FunFlux 2.1.0 ([Adapter trimming and the adapter screen](../methods/adapter-trimming.md)) |
    | **barcode block** | in a native-barcoded ONT library, the 40 bp that follow the adapter: a short flank, the 24-bp barcode that marks the sample and a second flank; the adapter and this block are the 68 bp trimmed from these reads |
    | **concatemer** | a contig that holds more than one copy of the mitochondrial genome, joined end to end; flagged (`concatemer_flag` Y) when its `length_ratio`, the contig length divided by that of the mitochondrial genome it best matches in BLAST, is above 1.5 and the assembler did not report it circular; never trimmed |
    | **MGISEQ** | a short-read sequencer from MGI Tech (DNBSEQ technology) that produces paired reads like Illumina's; FunFlux ran such reads as Illumina pairs |

!!! warning "The numbers record what FunFlux does on a well-studied genome, not its accuracy on a new one"

    - **One genome** for ONT and Illumina reads (n = 1). Nothing here shows how FunFlux does with
      ONT or Illumina reads of other fungi. The HiFi runs, on another genome, have their own
      limitations ([below](#pacbio-hifi-reads-cordyceps-cateniannulata-gxu-8616)).
    - **Two versions of FunFlux.** The nine runs made with FunFlux 2.0.0 were not repeated with
      2.1.0. Version 2.1.0 changes nothing in how Illumina, MGISEQ and HiFi reads are filtered and
      assembled. Its new BLAST settings give some contigs another genus in the contamination
      screen, but with the default settings they change no genome file
      ([FunFlux 2.0.0 and 2.1.0 compared](../methods/blast-hit-depth.md#funflux-200-and-210-compared)).
      Its mitochondrion screen would change the files of the `illumina` run: three contigs of
      human mitochondrial DNA would stay in `contigs_final.fasta` instead of going to
      `contigs_mitochondrial.fasta` ([Other results](#other-results)). The contig summaries of
      the nine runs took the read depth from BlobTools, as FunFlux 2.0.0 did.
      The four runs with ONT reads were made before the last changes of 2.1.0: InterProScan
      without EBI's lookup service, how the mitochondrion screen reads the titles of BLAST
      records and checks the organism of its best match, and how the contig summary measures
      read depth. InterProScan asked the lookup service, got no answer and computed every match
      itself, as it now does from the start. Run again with the other changes on the runs' own
      files, the screen gave the same result in all four, and the contig summary the same
      verdicts.
    - **Tiberius** has seen genes of this species. *B. cinerea* is one of the species the Tiberius
      `fungi` model was trained on ([its training list](citation.md#validation-references)). What
      Tiberius adds here may be more than it adds on a species it has not seen.
    - **A stored AUGUSTUS model** for *B. cinerea*, `botrytis_cinerea`, is in the funannotate
      database. Only one comparison run used it. Every other run trained AUGUSTUS on the genome's
      own BUSCO genes.
    - **The protein evidence** given to funannotate was the reference annotation's own 13,703
      proteins, beside UniProt Swiss-Prot, which FunFlux always adds and which holds 194 reviewed
      B05.10 proteins. It pulls the gene models towards the reference. The agreement with the
      reference is therefore an upper bound, not what a new species gets with proteins of related
      fungi. Every run had the same evidence, so the runs stay comparable with each other.
    - **The reference annotation** is not a perfect answer key. NCBI lists it as submitted by
      Syngenta Biotechnology, Inc. and carried into RefSeq on 2018-04-17 by its annotation
      propagation pipeline. The RefSeq records are marked provisional (not yet given final NCBI
      review). The genome paper, whose authors include Syngenta staff, predicted the genes, used
      RNA-seq of more than 20 samples to check and improve them, and curated chromosome 1 by hand
      ([van Kan et al. 2017](citation.md#validation-references)). A FunFlux gene that differs from
      the reference is not necessarily wrong.
    - **Run times** come from one shared server, with other runs going on at the same time, so they
      are indicative only.

## What this page is { #scope }

The page is a record of what FunFlux 2.1.0 and 2.0.0 produced in these runs, with their settings.
It checks that each mode gives a sensible genome and gene set, and it measures how much two
annotations of the same genome, with the same gene-finder settings, differ. It is not a
benchmark, because it has one genome per read type and one replicate. Nor is it a comparison with
other tools or pipelines, or a measure of accuracy on a species the gene finders have not seen.

B05.10 was chosen because it has a gapless genome, a gene set and public reads of the strain
itself. The assembly, the genes, the mitochondrion and the β-tubulin change can therefore be set
against published records:

- **A gapless reference genome**, which its authors call near-finished: 18 chromosomes in 18
  sequences without gaps, 42.6 Mb.
- **A reference annotation** of 11,698 protein-coding genes.
- **Public reads of the strain**: ONT reads from the reference's own sample record, and Illumina
  pairs of a lab stock.
- **Features to check**: telomeres at both ends of 10 of the 18 reference chromosomes
  ([van Kan et al. 2017](citation.md#validation-references)), and the β-tubulin change E198A, which
  the reference protein carries too ([AMR worked example](../analysis/amr.md#worked-example)).

## Data and settings { #data }

The B05.10 runs used the public reference of the strain and two public read sets. Depth is the
number of bases divided by the reference length.

| Data | Accession | Details | Used by |
|---|---|---|---|
| Reference genome and annotation | GCF_000143535.2 (assembly ASM14353v4) | 18 chromosomes, 42,630,066 bp, PacBio reads assembled in 2015; 11,698 protein-coding genes with 13,703 mRNAs and proteins (RefSeq, 2018) | its proteins: column 3 of the [annotation parameter table](../reference/configuration.md#annotation-parameter-table) in every run; its genome and genes: the comparisons on this page |
| ONT reads | SRR40741277 | PromethION, 315,306 reads, 5.79 Gb (about 136x); BioSample SAMN02953628, the reference's own | `longread` 136x (hifiasm and Flye), `hybrid` |
| ONT subsample | the first 95,000 reads of SRR40741277 | 1.78 Gb (about 42x) | `longread` 42x + Medaka |
| Illumina reads | SRR10692911 | HiSeq 2000, 11,963,877 pairs, reads of up to 150 bp, 3.35 Gb (about 79x); BioSample SAMN13563037, a stock named B05.10-1 | `illumina`, `hybrid` |

The ONT reads come from the same sample record as the reference. The Illumina reads come from
another lab's stock of the strain, so a few real differences between the two stocks are possible.

The ONT reads come from a native-barcoded library. Reads 1 to about 95,000 of the file carry
almost no adapter or barcode, most likely because they were trimmed at basecalling, while most of
the later reads carry both. The 42x subsample is those first reads, so FunFlux 2.1.0 counted them
as already trimmed and left them as they were
([Adapter trimming in the runs with ONT reads](#adapter-trimming)).

The public ONT file names no basecaller model, so Medaka cannot work its model out: with `auto`, it
reads the model from the read headers ([The Medaka model](../methods/medaka-model.md)). The two
136x runs and the `hybrid` run therefore skip Medaka. The 42x run names a model,
`r1041_e82_400bps_sup_v5.2.0`, which is a guess, because the flow cell and the basecaller are not
known.

In the next table, the five mode runs are those in `illumina`, `longread` and `hybrid` mode, and
the five comparison runs are the annotations in `contigs` mode. The table lists the settings that
differ from the configuration template, apart from paths: the input and output folders, the
databases, the annotation table, InterProScan and the SignalP 6 archive. Keys missing from the run
configurations took their defaults. Each run used the template of its own version. Apart from the
new UniVec download address and its checksum, the two templates differ in two settings, which the
table leaves out: `blast_max_target_seqs`, 50 in FunFlux 2.0.0 and 10 in 2.1.0, and `ont_adapter`, new in
2.1.0. The four runs with ONT reads kept its default, `detect`, so FunFlux worked out from the
reads what to trim, if anything.

| Setting | Template default | The five mode runs | The five comparison runs |
|---|---|---|---|
| `mode` | `illumina` | `longread` (three times), `illumina`, `hybrid` | `contigs`, on the `contigs_final.fasta` of the 2.0.0 `hybrid` run |
| `assembler` (under `longread` or `hybrid`) | `hifiasm` | `hifiasm`; `flye` in one 136x run | not used |
| `medaka_model` (under `longread` or `hybrid`) | `auto` | `false` (the two 136x runs and `hybrid`); `r1041_e82_400bps_sup_v5.2.0` (42x) | not used |
| `short_read_guidance` (under `hybrid`) | `false` | `true` in the `hybrid` run: Filtlong guided by the Illumina reads, with `filtlong_split: 1000` | not used |
| `resources.threads` | 24 | 16 to 32 (the run-time table below) | 20 |
| `resources.ram_gb` | 128 | 256 | 256 |
| `directories.genemark_dir` | empty: no GeneMark-ES | a GeneMark-ES installation, so GeneMark-ES on | the same, but `genemark.mode: off` in two runs |
| `parameters.tiberius.mode` | `on` | `on` | `on`, but `off` in two runs |
| `parameters.signalp6.mode` | `off` | `fast` | `fast` |
| `parameters.amr.chroquetas` | `off` | `on` | `off` |
| Annotation table, column 3 (proteins) | proteins of related fungi | the reference's own proteins | the same |
| Annotation table, column 4 (AUGUSTUS model) | `-`: trained on the genome | `-` | `-`, but `botrytis_cinerea` in one run |

One of these changes shapes every gene result of the B05.10 runs: column 3 of the annotation table
gave the reference's own proteins as evidence (see the limitations above). The other changes
concern the mode, the assembler, the Medaka model, the Filtlong guidance of the `hybrid` run, and
the CPU and memory. Three optional tools were also switched on: GeneMark-ES, SignalP 6 and the AMR
screen. The five comparison runs, in addition, each changed one gene-finder setting, or two (the
run with neither GeneMark-ES nor Tiberius), or none (the replicate).

## Adapter trimming in the runs with ONT reads { #adapter-trimming }

Before Filtlong, FunFlux 2.1.0 looks for the ONT adapter at the start of every read and works out
the sequence to trim: the adapter and, in a barcoded library, the barcode block after it. When
fewer than 1 % of the reads start with the adapter, and fewer than 1 % with a barcode block, the
reads count as already trimmed and are left as they are; otherwise fastplong removes that
sequence. After assembly, the adapter screen searches every draft assembly for adapters and
barcodes ([Adapter trimming and the adapter screen](../methods/adapter-trimming.md)). The three
runs on all the ONT reads, the two 136x `longread` runs and the `hybrid` run, start from the same
read file and chose the same sequence to trim.

| What was counted | `longread` 136x (hifiasm and Flye) and `hybrid` | `longread` 42x + Medaka |
|---|---|---|
| Reads that start with the adapter | 121,568 of 315,306 (38.6 %), 92.9 % of them followed by native barcode NB24; 76,346 more start with the barcode block alone | 13 of 95,000 (0.01 %); 41 more start with a barcode block alone |
| Sequence trimmed | the adapter and the barcode block, 68 bp | none: fewer than 1 % of the reads start with the adapter or a barcode block, so the reads count as already trimmed |
| Reads trimmed by fastplong | 218,788; 20.3 Mb removed, 0.35 % of the bases; 72 reads cut in two where the sequence sat inside the read | none |
| Reads given to the assembler that still start with the adapter | 0.11 % (`hybrid`: 0.02 %) | 0.01 % |
| Contig ends with adapter or barcode | none | one: 33 bp of native barcode NB21 at the end of the chromosome 13 contig |

fastplong trims more reads than the scan counts, because it also searches the end of each read
and allows a quarter of the bases to differ, while the scan counts only exact matches to the core
of the adapter at the read start.

Trimming changed the nuclear genome only a little. The two 136x hifiasm runs, `longread` and
`hybrid`, still have all 18 chromosomes in one contig, but their rDNA pieces and some contig ends
differ, and the `longread` run has one contig more: a 70.9-kb second copy of part of
chromosome 7, at low read depth, with 18 genes.
hifiasm gives the same assembly again from the same reads, so these are effects of the trimming,
not run-to-run variation. The contig ends changed, and the mitochondrion most of all.

- **Contig ends.** The adapter screen found no adapter or barcode in the contigs of the three
  trimmed runs. The Flye contigs reach 21 telomeric chromosome ends in both versions; beyond the
  end of the reference chromosome, these carried a median of 72 bp in 2.0.0 and 8 bp in 2.1.0.
  The tails were most likely adapter and barcode, which untrimmed reads leave at contig ends, most
  often just outside the telomere
  ([Adapter trimming and the adapter screen](../methods/adapter-trimming.md)). Contigs with a
  telomere array at both ends numbered 15 in both versions in the hifiasm 136x run, 6 in both in
  the Flye run, and 16 against 14 in the `hybrid` run.
- **The 42x run** assembled the same reads as in 2.0.0 and gave the same assembly: its contigs
  are identical over the 99.9 % of their bases that align to the reference. The adapter screen
  still found 33 bp of native barcode NB21 at one contig end, probably from the 2 reads given to
  hifiasm that carry it; the library's own barcode is NB24. FunFlux reports such a contig end and
  does not trim it.
- **The mitochondrion** changed in all three 136x runs ([Other results](#other-results)).

## Assembly per mode { #assembly }

The first table sets the genome of each mode run against the 18 chromosomes of the reference.
Every run with long reads used hifiasm, the default assembler, except the one marked Flye, which ran
Flye 2.9.6 with `--nano-hq`, as `flye_input_mode: auto` chose.

| Run | FunFlux | Contigs | Total length (Mb) | N50 (Mb) | Chromosomes in one contig (of 18) | Share of each chromosome covered (%) | Contigs with telomeres at both ends |
|---|---|---|---|---|---|---|---|
| Reference (GCF_000143535.2) | – | 18 | 42.63 | 2.62 | 18 | – | 10 ([van Kan et al. 2017](citation.md#validation-references)) |
| `longread`, ONT 136x | 2.1.0 | 33 | 44.15 | 2.58 | 18 | 99.8 to 100 | 15 |
| `longread`, ONT 136x, Flye | 2.1.0 | 21 | 42.95 | 2.58 | 16 | 99.8 to 100 | 6 |
| `longread`, ONT 42x + Medaka | 2.1.0 | 18 | 42.97 | 2.63 | 18 | 99.8 to 100 | 14 |
| `illumina`, Illumina 79x | 2.0.0 | 409 | 42.23 | 0.31 | 0 | 96.0 to 99.8 | 0 |
| `hybrid`, ONT 136x + Illumina 79x | 2.1.0 | 34 | 44.24 | 2.58 | 18 | 99.8 to 100 | 16 |

On B05.10, the three hifiasm runs with long reads put each of the 18 chromosomes into one contig
and cover almost all of each. Flye, on the same 136x reads, is as contiguous (N50 2.58 Mb), but
puts 16 chromosomes into one contig and ends fewer contigs in a telomere array. The Illumina-only
run gives 409 contigs, no chromosome in one piece, and leaves up to 4 % of a chromosome uncovered.

- **Contigs and length** are those of `contigs_final.fasta`, the nuclear genome. The hifiasm 136x
  and `hybrid` assemblies are 1.5 and 1.6 Mb longer than the reference. Most of the extra length
  is in 13 and 14 contigs that the contig summary calls `rDNA_array` (1.02 Mb in the hifiasm 136x
  run, 1.08 Mb in the `hybrid` run), which are pieces of the rDNA array. The 42x and Flye runs
  have no such contigs.
- **The Flye run** has chromosomes 1 and 13 in two pieces each (2.72 and 1.39 Mb; 1.95 and 0.31
  Mb). Its genome is 0.32 Mb longer than the reference. It has four contigs fewer than the 2.0.0
  Flye run. Repeat pieces of 7.6 to 17.5 kb that were separate contigs in 2.0.0 are now part of
  the chromosome 5 and 16 contigs, or were not assembled again.
- **Telomeres** come from FunFlux's [contig summary](../analysis/contig-summary.md): arrays of the
  repeat unit it chose, `AACCCT`, at both contig ends. `AACCCT` is tidk's canonical form of
  `TTAGGG`, the same repeat read on the other strand. The hifiasm long-read runs have 14 to 16
  such contigs, and in each, 2 to 4 more chromosomes lack an array at one end. In the Flye run, 23
  contig ends carry an array, against 32 to 35 in the hifiasm long-read runs. The reference's 10
  comes from its paper, not from FunFlux's search, so it is not the same measure.
- **Filtlong guidance** was on in the `hybrid` run. The unguided `longread` 136x run, on the same
  ONT reads, has 33 contigs against 34, and all 18 chromosomes in one contig in both. That is one
  reason the guidance is off by default ([Filtlong guidance](../methods/filtlong-guidance.md)).

The second table gives the differences of each assembly from the reference, and the run time.

| Run | FunFlux | Substitutions per Mb | Small indels per Mb | Run time (h) | Cores |
|---|---|---|---|---|---|
| `longread`, ONT 136x | 2.1.0 | 35.1 | 35.2 | 10.6 | 28 |
| `longread`, ONT 136x, Flye | 2.1.0 | 35.3 | 33.3 | 9.8 | 20 |
| `longread`, ONT 42x + Medaka | 2.1.0 | 37.9 | 30.9 | 9.2 | 16 |
| `illumina`, Illumina 79x | 2.0.0 | 14.8 | 12.2 | 7.9 | 24 |
| `hybrid`, ONT 136x + Illumina 79x | 2.1.0 | 38.2 | 23.9 | 11.1 | 32 |

Every assembly differs from the reference by 15 to 38 substitutions and 12 to 35 small indels per
Mb. In the long-read assemblies most substitutions sit in a few regions, and in every run most
small indels sit in runs of identical bases. Each run took 8 to 11 hours.

- **Substitutions** are not spread out. In each long-read assembly 95 % or more of them fall in
  ten windows of 10 kb of the reference, and 79 to 88 % in the first 50 kb of chromosome 2. Six
  or seven other chromosome contigs run past the ends of their reference chromosomes, and that
  extra sequence aligns there, to the most similar stretch of the reference. These substitutions
  therefore most likely mark chromosome-end sequence that the reference lacks, not errors along
  the whole genome.
- **Between runs**, the substitution counts move with those contig ends. The chromosome 10
  contig runs 13 kb past the end of its reference chromosome in the 2.0.0 `longread` run and
  45 kb in the 2.1.0 `hybrid` run, far enough to add 110 substitutions at the start of
  chromosome 2; in the other two runs it stops 7 to 8 kb past the end. The rDNA pieces of the
  2.1.0 `longread` and `hybrid` runs add 19 and 86 more, where they align to an rDNA unit at the
  start of chromosome 4 (2 in each 2.0.0 run). That is most of why the `hybrid` run now has
  more substitutions than the hifiasm 136x run (38.2 against 35.1 per Mb), where it had fewer in
  FunFlux 2.0.0 (32.0 against 38.9).
- **Small indels**, of 1 to 10 bp, are spread out. In every run three quarters or more of them
  (76 to 86 %) sit in a run of at least 3 identical bases of the reference. The `hybrid` assembly,
  corrected with the Illumina reads by Polypolish, has about a third fewer than the hifiasm 136x
  assembly without Medaka. The 42x run with Medaka also differs from the hifiasm 136x run in
  depth, so it does not measure Medaka alone, while [Medaka on the 42x run](#medaka) does. The
  Illumina contigs differ least, but they cover less of each chromosome.
- **Not an error rate.** Some of these differences may be errors of the reference itself or real
  differences between lab stocks. The counts compare the runs with each other.
- **Run time** is that of the main Snakemake run, from its log. InterProScan took the largest
  share: 4.0 to 6.6 h of each 2.1.0 run, as it computed every match itself. EBI's lookup
  service, which returns matches already computed, did not answer InterProScan 5 in November
  2026, and FunFlux 2.1.0 does not ask it
  ([Functional annotation](../analysis/functional-annotation.md#interproscan)).

To count the chromosomes and the differences, each `contigs_final.fasta` was aligned to the
reference with minimap2 2.31 (`-x asm5 -c --cs`). A contig is placed on the chromosome that most of
its bases align to, and a chromosome counts as one contig when the largest contig placed on it is
at least 97 % of its length. The share of each chromosome covered counts alignments with a mapping
quality of at least 5. Differences are counted in the alignments with a mapping quality of at
least 30, per Mb of aligned contig.

### Medaka on the 42x run { #medaka }

The 42x run keeps both the assembly that Medaka polished (`contaminants/contigs_decontam.fasta`)
and Medaka's result (`medaka/consensus.fasta`). Medaka can therefore be measured alone, on the same
reads and the same assembly, before and after. QUAST 5.3.0 compared both with the reference, and
BUSCO 6.0.0 scored both, and the reference, with `ascomycota_odb12`. These numbers were measured
on the 2.0.0 run. The 2.1.0 run gave the same assembly before and after Medaka
([Adapter trimming in the runs with ONT reads](#adapter-trimming)), so they hold for it too.

| | Before Medaka | After Medaka | Reference |
|---|---|---|---|
| Indels | 1,879 | 1,590 | – |
| of them up to 5 bp | 1,717 | 1,458 | – |
| of them longer | 162 | 132 | – |
| Mismatches | 1,684 | 1,698 | – |
| Contigs, misassemblies, genome fraction | 19, 30, 99.98 % | the same | – |
| BUSCO complete (of 2,826) | 99.2 % | 99.2 % | 99.3 % |
| BUSCO genes with internal stop codons | 483 | 479 | 482 |

Medaka removed 289 indels, 15 %, most of them 5 bp or shorter, and added 14 mismatches. The
number of contigs and of misassemblies, the genome fraction and the BUSCO completeness did not
change. The BUSCO genes with internal stop codons went from 483 to 479, but Medaka's gain is too
small to show in the genes. The reference itself has as many such genes (482), so that count
measures BUSCO's alignments, not the assembly. Medaka took about 15 minutes on 12 threads.

- **The model**, `r1041_e82_400bps_sup_v5.2.0`, is a guess, because the public reads carry no
  basecaller tag in their headers. The reads look like R10.4.1 reads called with a "sup" model
  (median identity to the reference 98.97 %, median read quality Q20.5, on 20,000 reads), so the
  guess is plausible.
- **These counts** are QUAST's, over all aligned bases. The table of differences under
  [Assembly per mode](#assembly) counts differently (minimap2, mapping quality of at least 30,
  indels of 1 to 10 bp), so the two tables cannot be compared number for number.

## Gene models per mode { #gene-models }

The gene models of each mode run were set against the reference annotation, to see how much the
assembly changes the genes.

| Run | FunFlux | Protein-coding genes | BUSCO complete (of 1,312) | Reference genes found | FunFlux genes that match | Matching genes |
|---|---|---|---|---|---|---|
| Reference (RefSeq) | – | 11,698 | – | – | – | – |
| `longread`, ONT 136x | 2.1.0 | 11,837 | 1,299 (99.0 %) | 90.6 % | 89.6 % | 10,606 |
| `longread`, ONT 136x, Flye | 2.1.0 | 11,798 | 1,302 (99.2 %) | 90.2 % | 89.5 % | 10,564 |
| `longread`, ONT 42x + Medaka | 2.1.0 | 11,789 | 1,300 (99.1 %) | 90.2 % | 89.6 % | 10,564 |
| `illumina`, Illumina 79x | 2.0.0 | 11,704 | 1,302 (99.2 %) | 90.4 % | 90.5 % | 10,590 |
| `hybrid`, ONT 136x + Illumina 79x | 2.1.0 | 11,838 | 1,302 (99.2 %) | 90.5 % | 89.6 % | 10,602 |

!!! note "On B05.10, the assembly mode and the assembler made little difference to the genes"

    The runs predict 11,704 to 11,838 protein-coding genes, and 10,564 to 10,606 of them match a
    reference gene: a spread of 42, under 1 % of the genes.

- **BUSCO complete** is the number of BUSCO `dikarya_odb9` genes found complete among the
  predicted proteins, as funannotate annotate reports it. The five runs differ by at most three
  BUSCO genes.
- **Between versions**, each 2.1.0 run has 8 to 38 more genes than the 2.0.0 run of the same
  mode, and between 21 fewer and 15 more matching genes. The 42x run annotated the same genome in
  both versions and still differs by 11 genes and 10 matching genes, so differences of this size
  can arise without any change in the genome ([Repeated runs](#gene-finders)). In the hifiasm 136x
  run, 18 of the 35 extra genes lie on the low-depth copy of part of chromosome 7 that only the
  trimmed reads produced ([Adapter trimming in the runs with ONT reads](#adapter-trimming)).
- **Reference genes found** is the share of the reference's gene loci with a matching FunFlux gene
  (gffcompare's locus sensitivity). **FunFlux genes that match** is the share of FunFlux genes with
  a matching reference gene (locus precision). Both are about 90 % in every run.
- **A match** means that the coding part of the FunFlux gene has the same introns as a reference
  transcript, while its start and stop may differ. A gene without introns matches when the two
  overlap almost completely.
- **tRNA genes**, 211 in every run, are left out of all counts.

To compare the genes, Liftoff 1.6.3 carried the reference genes onto the genome FunFlux annotated:
all 11,698 of them, but 11,697 onto the Illumina assembly. Then gffcompare 0.12.10 compared the
coding parts of the two gene sets
([Shumate & Salzberg 2021; Pertea & Pertea 2020](citation.md#validation-references)). These tools
were used for this page only and are not part of FunFlux.

## Repeated runs and the gene finders { #gene-finders }

To tell run-to-run noise from the effect of each gene finder, one genome was annotated six times,
all with FunFlux 2.0.0. It is the genome of the 2.0.0 `hybrid` run (33 contigs, 44.05 Mb), not
the 2.1.0 `hybrid` genome in the tables above. That run annotated it, and so did five
`contigs`-mode runs on its `contigs_final.fasta`. The five ran at the same time, each with 20 CPU
threads instead of the `hybrid` run's 32. All six annotated the same masked genome, byte for
byte. The reference run, the 2.0.0 `hybrid` run, used Tiberius, GeneMark-ES and AUGUSTUS trained
on the genome. GeneMark-ES is off by default, so the run without it is the one with the default
settings.

| Run | What changed, besides the mode and the threads | Protein-coding genes | Genes identical to the reference run | Matching genes | Matching genes, difference from the reference run |
|---|---|---|---|---|---|
| Reference (the 2.0.0 `hybrid` run) | – | 11,830 | – | 10,623 | – |
| Replicate | nothing | 11,825 | 11,718 | 10,619 | −4 |
| No GeneMark-ES (the default settings) | `genemark.mode: off` | 11,834 | 11,253 | 10,590 | −33 |
| No Tiberius | `tiberius.mode: off` | 11,729 | 10,494 | 10,257 | −366 |
| Neither | both off; AUGUSTUS, SNAP and GlimmerHMM remain | 11,503 | 9,528 | 9,630 | −993 |
| Stored AUGUSTUS model | annotation table, column 4: `botrytis_cinerea` | 11,934 | 11,277 | 10,694 | +71 |

In these runs, two annotations of the same genome with the same gene-finder settings differ at
about 1 % of the genes, but by only 4 matching genes. Differences of a few hundred genes are
therefore unlikely to be noise. Tiberius matters most: without it about 3 % of the matching genes
are lost, and about 9 % without both Tiberius and GeneMark-ES. GeneMark-ES on its own makes a small
difference (33 genes), which one replicate cannot separate from noise. The stored *B. cinerea*
AUGUSTUS model adds 71. Read these numbers with the head start described in the note below.

- **The replicate.** 112 of the 11,830 genes of the reference run differ in the replicate. The
  matches lost (38) and gained (34) almost cancel, hence the difference of 4.
- **Noise.** One replicate is one measurement of the run-to-run noise, not its range, and it holds
  any effect of the thread count as well.
- **The stored AUGUSTUS model** is AUGUSTUS's model for *B. cinerea* itself.
- **Identical** means the same coding intervals: every exon, the start and the stop.
  [Run-to-run differences](../methods/run-to-run-differences.md) explains which steps differ
  between runs with the same input and settings, and why.
- **Run time.** Each comparison run took 7.6 to 9.5 h.

!!! note "Tiberius and the stored AUGUSTUS model had a head start on this species"

    Tiberius was trained on *B. cinerea* among many other fungi, and the stored AUGUSTUS model is a
    *B. cinerea* model. Both may therefore add more here than on a species they have not seen.
    Every run also had the reference's own proteins as evidence, which brings all runs closer to
    the reference and may hide part of what the gene finders add. These runs are consistent with
    Tiberius being on by default and GeneMark-ES being optional
    ([Gene prediction](../analysis/gene-prediction.md#tiberius)), but they do not measure what
    Tiberius adds on a new species.

## Other results { #other-results }

The ITS taxonomy and the AMR screen gave the same result in all five mode runs, and the four runs
with ONT reads gave the same result in both versions. In these runs they depend neither on the
mode nor on the assembler:

- **ITS taxonomy.** The main ITS2 sequence gives the genus *Botrytis* (confidence 0.99) and the
  species *B. caroliniana* (0.83), not the strain's species. The main ITS1 sequence stops at the
  family Sclerotiniaceae (0.93). See [ITS taxonomy](../analysis/its-taxonomy.md).
- **AMR screen.** It finds β-tubulin E198A, with the verdict `resistance` (carbendazim,
  iprodione). Bos1 and OS1 are wild type, and cytochrome b is not evaluated. See the
  [worked example](../analysis/amr.md#worked-example) of the AMR page.

The mitochondrion was found in every run ([Mitochondria](../analysis/mitochondria.md)). The next
table gives it per run, with the verdicts of the contig summary
([worked example](../analysis/contig-summary.md#worked-example), from the FunFlux 2.0.0
`longread` 136x run). The `illumina` column comes from FunFlux 2.0.0, the others from 2.1.0.

| Per run | `longread` 136x | `longread` 136x, Flye | `longread` 42x + Medaka | `illumina` (2.0.0) | `hybrid` |
|---|---|---|---|---|---|
| Mitochondrial contigs, in `contigs_mitochondrial.fasta` (bp) | 139,240: a full copy, then 41 kb of it again and a 15.8-kb inverted copy of part of it; flagged as a concatemer | 79,494 and 2,842: together one complete copy | 116,252: a full copy with 34 kb repeated at its ends; not flagged | 82,466, plus three contigs of human mitochondrial DNA, 572 to 1,139 bp | 122,980: a full copy with 41 kb repeated at its ends, one and a half copies; not flagged |
| Contig summary: `chromosome` | 15 | 6 | 14 | 0 | 16 |
| `chromosome_arm` | 3 | 11 | 4 | 10 | 2 |
| `rDNA_array` | 13 | 0 | 0 | 0 | 14 |
| `nuclear_contig` | 0 | 4 | 0 | 260 | 0 |
| `nuclear_fragment` | 0 | 0 | 0 | 31 | 0 |
| `repeat` | 0 | 0 | 0 | 73 | 0 |
| `low_depth` | 2 | 0 | 0 | 34 | 2 |
| `unclassified` | 0 | 0 | 0 | 1 | 0 |
| `mitochondrion` | 1 | 2 | 1 | 4 | 1 |
| All contigs | 34 | 23 | 19 | 413 | 35 |

No 2.1.0 run gives the mitochondrion as one circular contig of about the published length,
82,212 bp (KC832409). In FunFlux 2.0.0 both 136x `longread` runs did: 82,351 bp with hifiasm and
82,338 bp with Flye. In the 2.1.0 `longread` 136x run, hifiasm gives the 139,240-bp contig of the
table. Flye gives two linear contigs that together make one complete copy, 82,336 bp: a single
read adds a side branch to its assembly graph at both places where the two contigs meet, so Flye
did not close the circle. The `hybrid` run went from about two copies (175,274 bp) to one and a
half.

In the `longread` 136x run, the trimming caused the change. Run again on the same trimmed reads,
hifiasm gave the same assembly, base for base. Run on the untrimmed reads, filtered by Filtlong as
in 2.0.0, it gave the 2.0.0 assembly back, with the circle of 82,351 bp. Trimming 218,788 reads
(69 %) was enough to send hifiasm along another path through the tangled mitochondrial part of its
assembly graph; the nuclear genome hardly changed. For Flye this was not tested. Whether trimming
helps or harms the mitochondrion in general cannot be said from these runs: in the `hybrid` run
the contig went from about two copies to one and a half. On these reads, small changes in the
input changed the mitochondrion far more than the nuclear genome.

- **Check a mitochondrial contig** before you use it: compare its length with a published
  mitochondrial genome of the species, if there is one, and read `length_ratio` and
  `concatemer_flag` in `mitochondrial_screen.tsv`
  ([Concatemers](../analysis/mitochondria.md#concatemers)). The flag needs a `length_ratio` above
  1.5, so it misses contigs of about one and a half copies: the `hybrid` contig (ratio 1.50) and
  the 42x contig (1.41) are not flagged.
- **The three small contigs of the `illumina` run** are human mitochondrial DNA: their best BLAST
  hits, over their whole length, are human mitochondrial genomes. They most likely come from a few
  human reads in the Illumina library. FunFlux 2.0.0 confirmed them, because they match
  mitochondrial genomes, and put them in `contigs_mitochondrial.fasta` beside the fungal
  mitochondrion. FunFlux 2.1.0 also checks the organism of the best match. Run again on the guide
  assembly of the `hybrid` run, which holds the same three contigs because it is built from the
  same Illumina reads, its screen calls them `foreign_mitochondrial`. In an `illumina` run they
  would therefore stay in `contigs_final.fasta`, where the contamination screen treats them like
  any other contig and the contig summary calls them `foreign_mitochondrion`
  ([Mitochondria](../analysis/mitochondria.md#foreign-mitochondrial)).

In the contig summary, every hifiasm run with long reads has 18 contigs that are `chromosome` or
`chromosome_arm`, one per chromosome. The Flye run has 17 such contigs, some of them small pieces,
and four `nuclear_contig`s of 2.2 to 3.0 Mb with no telomere array. In 2.1.0 its chromosome 16
contig is 32 kb longer, and it counts as a `chromosome_arm` because it has a telomere array at
that new end. The 32 kb are subtelomeric repeat that the reference does not have at that end, and
the outer 6 kb, with the array, have about half the usual read depth, so this telomere cannot be
confirmed. The Illumina contigs are mostly `nuclear_contig`, with no telomere at both ends.

The contig summary describes the draft assembly, before the mitochondrion is split off and before
polishing. Each run therefore counts its mitochondrial contigs (one, two in the Flye run, or four
in `illumina`) on top of the nuclear contigs of the assembly table. The contig summary also lists
the 42x and `hybrid` mitochondria at their draft lengths, 116,228 and 122,952 bp.

## PacBio HiFi reads: *Cordyceps cateniannulata* GXU-8616 { #pacbio-hifi-reads-cordyceps-cateniannulata-gxu-8616 }

Three runs used PacBio HiFi reads of one strain of an insect-pathogenic fungus (family
Cordycipitaceae). Two ran in `longread` mode, one with hifiasm, the default assembler, and one
with Flye. The third ran in `hybrid` mode, which adds MGISEQ short reads of the same strain. On
this genome the runs show that the HiFi path works from reads to annotation, and how close its
genome comes to a chromosome-scale reference. The three runs were made with FunFlux 2.0.0;
FunFlux 2.1.0 does not trim HiFi reads, so its main change does not concern them. Their contig
summaries took the read depth from BlobTools' estimate, as in FunFlux 2.0.0; the depth ratios
below were not measured again with the samtools coverage of FunFlux 2.1.0
([Contig summary](../analysis/contig-summary.md#where-the-evidence-comes-from)).

!!! warning "The HiFi runs do not measure the accuracy of the genes"

    - **One genome** (n = 1). Nothing here shows how FunFlux does on HiFi reads of other fungi.
    - **The reference** comes from the same sample. It was assembled from HiFi reads with hifiasm
      0.19.8, polished with short reads and joined into 7 chromosomes with Hi-C data
      ([Liu et al. 2026](citation.md#validation-references)). The sample's only HiFi run in the
      archive is the one used here. The reference is chromosome-scale with the help of Hi-C, while
      reads alone may not give every chromosome in one piece. Agreement with it, above all for the
      hifiasm runs, is expected and is not an independent check.
    - **No gene comparison.** The reference's gene annotation (10,743 protein-coding genes at NCBI)
      was not compared with FunFlux's genes. The gene table gives counts and BUSCO only.
    - **Tiberius** has seen no genes of this genus, only genes of the same family. No *Cordyceps*
      species is among the training species of the Tiberius `fungi` model, but two species of the
      same family, *Akanthomyces muscarius* and *Beauveria bassiana*, are
      ([its training list](citation.md#validation-references)).
    - **The protein evidence** comes from another species, *Cordyceps militaris* CM01, as for a
      new genome. The B05.10 runs had the reference's own proteins, which pull the genes towards
      the reference.
    - **Run times** come from one shared server, with other runs going on at the same time, so they
      are indicative only.

### Data and settings { #hifi-data }

The reads and the reference come from one sample record. Depth is again the number of bases
divided by the reference length.

| Data | Accession | Details | Used by |
|---|---|---|---|
| Reference genome | GCA_056505955.1 (assembly ASM5650595v1) | 7 chromosomes, 32,707,717 bp, no gaps; BioSample SAMN48746374 | the comparisons below |
| PacBio HiFi reads | SRR34289959 | PacBio Revio, 147,159 reads, 2.77 Gb (about 85x), mean length 18.8 kb; the same BioSample | all three runs |
| MGISEQ short reads | SRR34289958 | MGISEQ-2000RS, 16,427,628 pairs of 2 x 150 bp, 4.93 Gb (about 151x); the same BioSample | `hybrid` |
| Protein evidence | GCF_000225605.1 | *Cordyceps militaris* CM01, 9,651 proteins (RefSeq) | column 3 of the annotation table, all three runs |

The MGISEQ short reads are DNBSEQ pairs, and `hybrid` mode reads them as it reads Illumina pairs.
In every run Filtlong kept 124,832 HiFi reads, 2.49 Gb, which is 90 % of the bases
([Read filter](../modes/longread.md#read-filter)).

The next table gives the settings of the three HiFi runs.

| Setting | `longread`, hifiasm | `longread`, Flye | `hybrid` |
|---|---|---|---|
| `mode` | `longread` | `longread` | `hybrid` |
| `assembler` (under `longread` or `hybrid`) | `hifiasm` (the default) | `flye`, which runs with `--pacbio-hifi` | `hifiasm` (the default), then Polypolish with the short reads |
| `resources.threads` | 24 | 20 | 24 |
| Annotation table, column 3 (proteins) | the *C. militaris* proteins | the same | the same |

Everything else is set as in the B05.10 mode runs ([settings table](#data)): 256 GB of memory,
GeneMark-ES and SignalP 6 `fast` on, the AMR screen on, and AUGUSTUS trained on the genome. HiFi
reads get no Medaka.

### Assembly { #hifi-assembly }

The table sets each HiFi assembly against the 7 chromosomes of the reference, with the run time of
each run.

| Run | Contigs | Total length (Mb) | N50 (Mb) | Reference chromosomes in one contig (of 7) | Contigs with telomeres at both ends | Run time (h) | Cores |
|---|---|---|---|---|---|---|---|
| Reference (GCA_056505955.1) | 7 | 32.71 | 4.79 | 7 | 7 | – | – |
| `longread`, hifiasm | 18 | 34.58 | 4.79 | 6 | 7 | 6.6 | 24 |
| `longread`, Flye | 12 | 34.01 | 4.79 | 6 | 6 | 7.3, in two attempts | 20 |
| `hybrid`, hifiasm + Polypolish | 18 | 34.58 | 4.79 | 6 | 7 | 7.8 | 24 |

On this genome, both assemblers give six of the seven reference chromosomes as one contig each,
with a telomere array at both ends, and cover all of every reference chromosome. The seventh breaks
in two at its rDNA array. Both assemblies also hold about 1.3 Mb that the reference lacks.

- **Chromosome 7** (3.58 Mb) comes as a piece of 3.10 Mb and one of 0.52 Mb, or 3.07 and 0.51 Mb
  with Flye. Each piece ends in a telomere array. The break lies in the tandem copies of the rDNA
  unit, the repeat that holds the ITS.
- **The extra 1.3 Mb.** hifiasm gives it as one contig of 1,329,053 bp with a telomere array at
  both ends, the seventh `chromosome` in the next table. Flye gives it as four contigs, 1.30 Mb
  together. Only 4 % of it aligns to the reference, in short pieces. Its read depth is close to
  that of the chromosomes (1.14 times the median), and FunFlux predicts 275 genes on it. The
  contamination and mitochondrion screen places it in *Cordyceps* from its BLAST hits. Whether it
  is a small chromosome that the reference lacks was not checked further.
- **Other contigs.** hifiasm adds 7 small contigs of rDNA units (`rDNA_array`, 0.34 Mb together)
  and 2 contigs with almost no read depth (`low_depth`, 0.17 Mb) that do not align to the
  reference. Flye has neither, and its genome is 0.57 Mb smaller.
- **`hybrid`.** Its draft assembly is the hifiasm run's, byte for byte, because both runs filter
  and assemble the same HiFi reads. Polypolish then changed 38 positions with the short reads.
- **Run time** is that of the main Snakemake run, from its log. The Flye run stopped once in
  InterProScan, with an error that did not come back when the run was started again by hand. Its
  time is the sum of both attempts, 3.7 and 3.6 h. The run was made without the second attempt
  that FunFlux makes by itself after a failed InterProScan job
  ([Running FunFlux](../reference/running.md#run-time)).

The assemblies were checked as for B05.10, with minimap2 2.31 (`-x asm5`), and a chromosome counts
as one contig when the largest contig placed on it is at least 97 % of its length. The reference's
telomeres were counted with FunFlux's own telomere search, run on the reference, which found
`AACCCT` arrays at both ends of all 7 chromosomes.

The next table gives the contig summary verdicts of each run and its mitochondrial contig.

| Per run | `longread`, hifiasm | `longread`, Flye | `hybrid` |
|---|---|---|---|
| Contig summary: `chromosome` | 7 | 6 | 7 |
| `chromosome_arm` | 2 | 2 | 2 |
| `rDNA_array` | 7 | 0 | 7 |
| `nuclear_contig` | 0 | 3 | 0 |
| `low_depth` | 2 | 0 | 2 |
| `unclassified` | 0 | 1 | 0 |
| `mitochondrion` | 1 | 1 | 1 |
| All contigs | 19 | 13 | 19 |
| Mitochondrial contig, in `contigs_mitochondrial.fasta` (bp) | 49,194, circular: about two copies of the genome end to end | 24,597, circular: one copy | as the hifiasm run |

With hifiasm, all 14 ends of the reference chromosomes and both ends of the extra contig carry a
telomere array. With Flye, the 14 reference ends do too, but the extra 1.3 Mb has an array at one
end only and comes as three `nuclear_contig`s and one `chromosome_arm`. Flye's 0.51-Mb piece of
chromosome 7 is `unclassified`, because its read depth, 1.74 times the median, is above the 1.5
limit.

!!! warning "A circular mitochondrial contig of two copies is not flagged as a concatemer"

    The mitochondrion published for this strain is 24,712 bp
    ([Liu et al. 2026](citation.md#validation-references)). Flye gives one copy. hifiasm marks its
    49,194-bp contig as circular, so the concatemer flag stays `N` although the contig holds about
    two copies. FunFlux sets the flag to `Y` only when `length_ratio` is above 1.5 and the
    assembler did not record the contig as circular. `length_ratio` is the contig length divided by
    the length of the mitochondrial genome it best matches in BLAST
    ([Concatemers](../analysis/mitochondria.md#concatemers)).

### Genes { #hifi-genes }

The table gives the gene count of each HiFi run and two BUSCO scores, one of the genome and one of
the gene models.

| Run | Protein-coding genes | BUSCO of the genome: complete, `cordycipitaceae_odb12` (of 4,122) | BUSCO of the gene models: complete, `dikarya_odb9` (of 1,312) |
|---|---|---|---|
| Reference (GCA_056505955.1) | 10,743 | – | – |
| `longread`, hifiasm | 10,924 | 4,112 (99.8 %) | 1,299 (99.0 %) |
| `longread`, Flye | 10,919 | 4,112 (99.8 %) | 1,300 (99.1 %) |
| `hybrid` | 10,921 | 4,112 (99.8 %) | 1,298 (98.9 %) |

The three runs agree, differing by at most 5 genes and 2 BUSCO genes. In each run 273 to 275 genes
lie on the 1.3 Mb that the reference lacks, so the counts cannot be set against the reference's
10,743 one to one.

- **BUSCO of the genome** comes from the [assembly QC](../analysis/assembly-qc.md#busco) on
  `contigs_final.fasta`. BUSCO chose `cordycipitaceae_odb12` itself (`busco.lineage: auto`).
- **BUSCO of the gene models** comes from funannotate annotate, with `dikarya_odb9`, as in the
  B05.10 gene table. The two columns use different gene sets and cannot be compared with each
  other.
- **`hybrid` and hifiasm** annotated nearly the same genome, 38 positions apart, and differ by 3
  genes. That is about as much as the B05.10 reference run and its replicate, which differ by 5
  genes ([Repeated runs](#gene-finders)).
- **tRNA genes**, 105 to 107 per run, are left out.

### AMR screen and ITS { #hifi-other }

The AMR screen and the ITS classification gave the same result in all three runs.

!!! note "A species outside the AMR database is reported, not failed, and the run goes on"

    *Cordyceps cateniannulata* is not among the 57 species of the ChroQueTas database, so the three
    runs report it as not screened. `amr_summary.tsv` has one row with status `not screened`, and a
    NOTE in the log gives the reason. The database holds no *Cordyceps*, and its only species of
    the same family is *Beauveria bassiana* ([The species](../analysis/amr.md#species)).

The ITS taxonomy is the same in all three runs. The main ITS2 sequence (158 bp) gives the strain's
own species, *C. cateniannulata* (confidence 0.96, and 0.97 for the genus *Cordyceps*). The main
ITS1 sequence (182 bp) stops at the phylum Ascomycota (class 0.60, below the cutoff of 0.8). The
main ITS sequences themselves are the same in all three runs. The hifiasm assemblies hold 43
copies of the main ITS1 and 44 of the main ITS2, the Flye assembly 4 of each. See
[ITS taxonomy](../analysis/its-taxonomy.md).

### To repeat the HiFi runs { #hifi-reproduce }

The HiFi runs can be repeated from public data, with FunFlux 2.0.0, which made them. First download
the reads, the protein evidence and the reference genome:

```bash
# the reads (ENA); ENA names the HiFi file _subreads, but it holds HiFi reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR342/059/SRR34289959/SRR34289959_subreads.fastq.gz \
     ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR342/058/SRR34289958/SRR34289958_1.fastq.gz \
     ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR342/058/SRR34289958/SRR34289958_2.fastq.gz

# the protein evidence and the reference genome (NCBI Datasets)
curl -L -o GCF_000225605.1.zip \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000225605.1/download?include_annotation_type=PROT_FASTA"
curl -L -o GCA_056505955.1.zip \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_056505955.1/download?include_annotation_type=GENOME_FASTA"
# unpack both into ncbi_dataset/data/<accession>/ (protein.faa; the genome FASTA);
# -o because both zip files hold files of the same name, such as README.md
unzip -o GCF_000225605.1.zip
unzip -o GCA_056505955.1.zip
```

Name the reads `GXU8616_pacbio.fastq.gz`, `GXU8616_R1.fastq.gz` and `GXU8616_R2.fastq.gz`. The
sample name, `GXU8616`, is the strain name without the hyphen. The `longread` runs had a folder
with the `_pacbio` file only.
Write the annotation table with the *C. militaris* proteins in column 3:

```text
#sample	species	proteins	augustus_species
GXU8616	Cordyceps cateniannulata	/path/to/ncbi_dataset/data/GCF_000225605.1/protein.faa	-
```

Then use the settings above, and run as in [How to reproduce](#reproduce), step 5.

## How to reproduce { #reproduce }

The B05.10 runs can be repeated from public data in five steps. Use the FunFlux version that made
each run: 2.1.0 for the four runs with ONT reads, 2.0.0 for the others. Even then the gene models
differ slightly from run to run ([Repeated runs](#gene-finders)), and a Flye assembly does not
repeat exactly ([Run-to-run differences](../methods/run-to-run-differences.md)). Another version
can differ more, as the mitochondrion of the 136x runs shows ([Other results](#other-results)).

1. Download the data:

    ```bash
    # the reference: genome, annotation and proteins (NCBI Datasets)
    curl -L -o GCF_000143535.2.zip \
      "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000143535.2/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=PROT_FASTA"
    unzip GCF_000143535.2.zip   # ncbi_dataset/data/GCF_000143535.2/: the genome FASTA, genomic.gff, protein.faa

    # the ONT reads (sra-tools)
    fasterq-dump SRR40741277 && gzip SRR40741277.fastq

    # the Illumina pairs (ENA)
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/011/SRR10692911/SRR10692911_1.fastq.gz \
         ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/011/SRR10692911/SRR10692911_2.fastq.gz
    ```

2. Name the reads `B0510_ont.fastq.gz`, `B0510_R1.fastq.gz` and `B0510_R2.fastq.gz`, as the modes
   expect ([Choosing a mode](../getting-started/choosing-a-mode.md)). Symbolic links are fine. For
   the 42x run, take the first 95,000 ONT reads, 4 lines each, as sample `B0510m`:

    ```bash
    zcat SRR40741277.fastq.gz | head -n 380000 | gzip > B0510m_ont.fastq.gz
    ```

3. Write the [annotation table](../reference/configuration.md#annotation-parameter-table),
   tab-separated, with the reference proteins in column 3. Add a line for `B0510m` for the 42x run,
   and write `botrytis_cinerea` in column 4 for the run with the stored AUGUSTUS model:

    ```text
    #sample	species	proteins	augustus_species
    B0510	Botrytis cinerea	/path/to/ncbi_dataset/data/GCF_000143535.2/protein.faa	-
    ```

4. Copy the configuration template and set the mode, the input folder, the databases and the
   settings in the [settings table](#data). For the comparison runs, run the `hybrid`
   configuration with FunFlux 2.0.0 as well, and use `contigs` mode with its
   `02.assembly/B0510/contigs_final.fasta`, named `B0510.fasta`.
5. Run each configuration, with `--cores` equal to its `resources.threads`
   ([Running FunFlux](../reference/running.md)):

    ```bash
    snakemake --sdm conda --cores 32 --configfile config/config_custom.yaml
    ```

To compare a B05.10 run with the reference, Liftoff carries the reference genes onto the genome
FunFlux annotated, `04.annotation/{sample}/contigs_sort.fasta`. Then gffcompare compares the coding
parts of the lifted reference genes with FunFlux's protein-coding genes. The commands are in the
box below.

??? note "The comparison commands are not part of FunFlux"

    ```bash
    REF=ncbi_dataset/data/GCF_000143535.2
    # 1. carry the reference genes onto the annotated genome (its contigs are named scaffold_N)
    liftoff -g $REF/genomic.gff -o lifted_reference.gff3 -u unmapped_reference_genes.txt -p 16 \
      -dir liftoff_tmp -polish 04.annotation/B0510/contigs_sort.fasta \
      $REF/GCF_000143535.2_ASM14353v4_genomic.fna
    # 2. rebuild the lifted reference transcripts from their coding parts only (CDS lines become
    #    exons), because FunFlux's genes have no untranslated regions
    awk -F'\t' 'BEGIN{OFS="\t"} /^#/ {print; next} $3=="exon" {next} $3=="CDS" {$3="exon"; print; next} {print}' \
      lifted_reference.gff3_polished | gffread -T -o lifted_reference_cds.gtf -
    # 3. FunFlux's protein-coding genes (-C leaves out the tRNA genes, which have no CDS)
    gffread -T -C -o funflux_coding.gtf \
      04.annotation/B0510/funannotate/annotate_results/Botrytis_cinerea_B0510.gff3
    # 4. compare: gffcmp_cds.stats holds the sensitivity, the precision and the matching counts
    gffcompare -r lifted_reference_cds.gtf -o gffcmp_cds funflux_coding.gtf

    # the chromosome check: align the nuclear genome to the reference
    minimap2 -x asm5 -c --cs -t 8 $REF/GCF_000143535.2_ASM14353v4_genomic.fna \
      02.assembly/B0510/contigs_final.fasta > contigs_vs_reference.paf
    ```

    All of these tools come from bioconda: Liftoff 1.6.3 (its own alignments made with minimap2
    2.24), gffread 0.12.9, gffcompare 0.12.10, and minimap2 2.31 for the chromosome check.
