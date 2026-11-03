# Run-to-run differences

Two FunFlux runs on the same input, with the same settings, databases and conda environments, do
not always write the same files. This page records which steps repeat exactly, which do not,
why, and what FunFlux does about it. It is a running record, and the page says so wherever a
step has not been checked yet.

The differences matter when you compare two runs, for example before and after changing a
setting. A difference that also appears when the same setting is run twice is not caused by the
setting. Compare genome files byte for byte, but compare gene models, and everything computed
from them, with some tolerance. When the difference matters, run each setting at least twice.

A step was checked by running it twice on the same input and comparing its output files, byte
for byte or, for gene models, model by model. A gene model counts as identical when all its
coding intervals (exons, start and stop) agree.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **identical gene model** | a gene model with the same coding intervals in two runs: every exon, the start and the stop |
    | **random seed** | the number a program's random choices start from; fixing it makes a later run repeat the same choices |
    | **BUSCO training genes** | the BUSCO genes funannotate looks for to train SNAP, GlimmerHMM and, without a stored AUGUSTUS model, AUGUSTUS |
    | **Tiberius** | a deep-learning gene finder with a `fungi` model its authors trained on many fungal genomes; on by default, it adds its gene models to funannotate's evidence |
    | **Perl hash order** | the order in which Perl lists the entries of a lookup table (a hash); it changes from run to run unless `PERL_HASH_SEED=0` fixes it |
    | **self-training gene finder** | a gene finder that learns its parameters from the genome it annotates, with no stored model; GeneMark-ES works this way |
    | **EVidenceModeler** | the program inside `funannotate predict` that builds one consensus gene per locus from the gene finders and protein alignments, weighting each source's votes |
    | **matching gene** | a FunFlux gene whose coding part has the same introns as a reference transcript, start and stop aside; an intronless gene matches on near-complete overlap |

## Steps that repeat exactly

The steps below gave the same result in both runs. Where one of their files is not identical byte
for byte, it differs only in a path, a start time or the order of its lines, as the remarks say.

| Step | Rule | What was compared | Remark |
|---|---|---|---|
| Contig filtering | `filter_contigs` | the FASTA, byte for byte | |
| funannotate clean and sort | `funannotate_preprocess` | the FASTA files, byte for byte | the result does not depend on the number of CPUs |
| tantan masking | `repeat_masking` with `masking_method: tantan` | the masked genome, byte for byte | |
| QUAST | `genome_assembly_evaluation` | the report files, byte for byte | |
| BUSCO of the genome | `completeness` | the scores | the summary files print the path of the input, which changes with the output folder |
| ITS windows, ITSx and SINTAX | `its_windows`, `its_extraction`, `its_classification` | the sequences and the classification | SINTAX runs with a fixed random seed (`--randseed 1`) on one thread; ITSx's summary file also prints its start time |
| funannotate's BUSCO training genes and SNAP | inside `funannotate_prediction` | `busco.final.gff3` and the SNAP predictions, byte for byte | |
| Protein alignments and AUGUSTUS | inside `funannotate_prediction` | the coordinates | the files list the same lines in a different order |
| Tiberius | `tiberius_prediction` | the GTF, byte for byte | two runs on the same genome, with 25 and 32 threads; TensorFlow warns at start-up that its CPU maths may give slightly different numbers, but here they did not change the gene models |
| Tiberius's GTF converted to GFF3 | `tiberius_gff3` | the GFF3, byte for byte | FunFlux fixes Perl's hash order (`PERL_HASH_SEED=0`, see EVidenceModeler below) |
| MultiQC | `multiqc` | the numbers | two of its data files hold absolute paths |

## Steps that do not repeat exactly

In two FunFlux runs on the same input, the steps below do not give the same result. The table
gives what differs, the cause where it was traced, and what FunFlux does about it.

| Step | What differs | Why | What FunFlux does |
|---|---|---|---|
| RepeatModeler (`masking_method: repeatmodeler_repeatmasker`) | the repeat library, and with it the masking and the genes predicted after it | RepeatModeler takes a random sample of the genome and seeds it from the clock | setting `parameters.repeatmodeler_seed` fixes the seed (`-srand`). Whether the whole annotation then repeats has not been checked. The default, tantan, repeats exactly. |
| GeneMark-ES, inside `funannotate_prediction` | its gene predictions | its self-training on the genome does not repeat exactly; the cause inside GeneMark-ES was not traced | nothing, because GeneMark-ES has no seed option |
| GlimmerHMM, inside `funannotate_prediction` | its gene predictions | not traced | nothing |
| EVidenceModeler 1.1.1, inside `funannotate_prediction` | the consensus gene at some loci, even when its inputs hold the same predictions | where two gene structures score the same, EVidenceModeler takes the one that comes first in a Perl hash (a lookup table), and Perl lists a hash in an order that changes from run to run. funannotate also writes EVidenceModeler's input files in an order that changes from run to run. | nothing yet. With Perl's hash order fixed (`PERL_HASH_SEED=0 PERL_PERTURB_KEYS=0`), EVidenceModeler gave the same result three times on one test region. Whether that makes the whole gene prediction repeat has not been checked. |
| Everything computed from the predicted proteins: InterProScan, eggNOG-mapper, antiSMASH, SignalP 6, `funannotate annotate` | the details | they inherit the differences of the gene models | nothing |
| Flye (`assembler: flye`) | the assembly | Flye does not repeat its assembly byte for byte | hifiasm is the default assembler |

!!! note "The final gene models of two runs on the same genome differ at a small share of the loci"

    Everything computed from the proteins inherits these differences. On the public genome
    *Botrytis cinerea* B05.10, two runs of the gene prediction on the same masked genome differed
    at 112 of 11,830 gene models, about 1 %. One run had 32 threads and the other 20. The number
    of genes that match the reference annotation went from 10,623 to 10,619. See
    [Validation on public genomes](../about/validation.md#gene-finders).

## Not checked yet

- The read-based front ends of v2 (`illumina`, `longread`, `hybrid`): read filtering, SPAdes,
  hifiasm, Medaka and Polypolish.
- The gene prediction without GeneMark-ES, and the effect of the number of threads on its own
  (the B05.10 runs above differ in the threads as well).
- Tiberius on other genomes, and with a graphics card (FunFlux runs it on the CPU).
- Genomes of other sizes and other fungi.
