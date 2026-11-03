# ITS on long-read contigs

This page explains why ITSx gets short windows instead of whole contigs, which can make it fail on
long-read assemblies. The ITS steps themselves, their settings and their output files are
described on [ITS taxonomy](../analysis/its-taxonomy.md).

FunFlux v2 cuts the ITS region out of the genome in two steps. First, nhmmer finds every 5.8S
gene. ITSx then runs on a short ITS window around each one, turned so that it reads from SSU to
LSU. Both steps run in every mode, on `contigs_final.fasta`.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **SSU, 5.8S and LSU** | the small-subunit, 5.8S and large-subunit rRNA genes of one rDNA unit, which reads SSU-ITS1-5.8S-ITS2-LSU |
    | **rDNA unit** | one copy of the repeat SSU-ITS1-5.8S-ITS2-LSU; FunFlux counts one per 5.8S gene found, that is, one per ITS window |
    | **rDNA array** | the tandem copies of the rDNA unit |
    | **ITS window** | the stretch cut around one 5.8S gene, 2,000 bp on each side by default (`its.window_flank`), turned to read from SSU to LSU before ITSx runs |
    | **anchor** | one of the four stretches ITSx finds, at the SSU end, both 5.8S ends and the LSU start; ITS1 and ITS2 are cut between them |
    | **SINTAX** | the VSEARCH classifier (`vsearch --sintax`) that names each ITS1 and ITS2 rank by rank against UNITE, with a confidence for each rank |
    | **v1 workflows** | the four workflows FunFlux 2.0.0 replaces: FunFlux 1.x (Illumina), Funnotator (assemblies), FunFluxL (long reads) and the unpublished FunFluxL+ (hybrid) |

The table compares the v1 workflows with FunFlux v2.

| | v1 workflows | FunFlux v2 |
|---|---|---|
| What ITSx gets | whole contigs, which can fail on long-read assemblies | one window of 4,001 bp (with the default flank) around each 5.8S gene, turned to read from SSU to LSU |
| Run time of ITSx | minutes on a whole genome | seconds, after a search for the 5.8S genes that takes about a minute for a 40-Mb genome |
| Names of the classified records | named after their contig | `{sample}_ITS1_1`, `{sample}_ITS2_1` and so on |

Because the windows are short and all read in the same direction, ITSx avoids the two failures
described below.

## Cause

On whole long-read contigs, ITSx can fail in two ways, and it then writes either nothing or a
wrong ITS2. ITSx finds four 45-bp anchors with nhmmer, at the end of the SSU, at both ends of the
5.8S and at the start of the LSU. It cuts ITS1 and ITS2 between them.

- **A failed strand decision.** On a multi-megabase contig ITSx also has to decide which strand
  to read, and that decision fails when the rDNA array lies on the minus strand. ITSx then reports
  `Domains found in wrong order, sequence may be chimeric`, with `SSU seems to overlap 5.8S` and
  `LSU seems to overlap 5.8S`, and writes nothing.
- **A far too long "ITS2" between distant anchors.** ITSx takes an SSU-end anchor and an
  LSU-start anchor that lie far apart on the contig. Between them lies an anchor for the end of a
  5.8S, but none for its start. ITSx then reports
  `Start of 5.8S sequence not found; ITS region is suspiciously long (> 1500 bp)` and writes no
  ITS1. Everything from that 5.8S end to the LSU start is written as ITS2, so SINTAX classifies a
  long stretch of genome as if it were an ITS2.

ITSx's `--allow_reorder` option changes the outcome in neither case.

## What the windows change

The turn of each window fixes the failed strand decision, and the window size, together with the
length limits, keeps a false ITS2 out of the classification.

- **The turn.** Each window is turned to read from SSU to LSU before ITSx sees it.
- **The window size.** A window of 4,001 bp holds part of one rDNA unit only. The 2 kb on each
  side of the 5.8S gene hold the SSU end and the LSU start, around an ITS1-5.8S-ITS2 region of
  0.4 to 1.5 kb.
- **The length limits.** ITS1 or ITS2 outside 50 to 600 bp is listed in `{sample}_ITS_copies.txt`
  and not classified, and this is what keeps a false ITS2 out of SINTAX. A full ITS1-5.8S-ITS2
  region longer than `its.max_full_length` (1,500 bp by default) is marked `rejected_too_long` in
  the same table. The full region itself is never classified. That default is the length above
  which ITSx itself calls an ITS region "suspiciously long".

!!! note "Windows from the minus strand fail too unless they are turned"

    Without the turn, windows from the minus strand fail in the same way as whole minus-strand
    contigs. The cause of that failure is how ITSx handles the strand, not the tandem copies or
    the window size.
