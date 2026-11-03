# ITS taxonomy

Stage `03.taxonomy` cuts the fungal ITS barcode (ITS1, 5.8S, ITS2) out of the nuclear genome,
`02.assembly/{sample}/contigs_final.fasta`, and classifies ITS1 and ITS2 against UNITE with
SINTAX. It runs in every mode, once per sample.

The ITS calls are for you to read, as a check on the identity of each isolate. No later step uses
them, and the [contig summary](contig-summary.md) only counts the 5.8S genes found.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **rDNA unit** | one copy of the repeat SSU-ITS1-5.8S-ITS2-LSU; FunFlux counts one per 5.8S gene found, that is, one per ITS window |
    | **SSU, 5.8S and LSU** | the small-subunit, 5.8S and large-subunit rRNA genes of one rDNA unit, which reads SSU-ITS1-5.8S-ITS2-LSU |
    | **ITS window** | the stretch cut around one 5.8S gene, 2,000 bp on each side by default (`its.window_flank`), turned to read from SSU to LSU before ITSx runs |
    | **profile HMM** | a hidden Markov model, a statistical model, that describes one sequence family; ITSx searches with such profiles |
    | **SINTAX** | the VSEARCH classifier (`vsearch --sintax`) that names each ITS1 and ITS2 rank by rank against UNITE, with a confidence for each rank |
    | **confidence (SINTAX)** | the share of SINTAX's repeated classifications, each on a random subset of k-mers, that give the same name at a rank; not a sequence identity |
    | **UNITE** | the reference database of ITS sequences, each with its taxonomy; by default, FunFlux uses the release for all eukaryotes and only its records named to genus |

For each sample, three rules run one after the other:

```text
contigs_final.fasta
 └─ its_windows ── one window per 5.8S gene, turned to read from SSU to LSU
     └─ its_extraction ── ITSx cuts ITS1, 5.8S and ITS2 out of each window
         └─ its_classification ── identical copies merged; ITS1 and ITS2 classified
                                   with SINTAX against UNITE (records with a genus)
```

Two more rules prepare the UNITE reference, once per run. The table lists all five, with the
tools they use and the files they read and write.

| Rule | Tools | Reads | Writes (in `03.taxonomy/{sample}/` unless a path is given) |
|---|---|---|---|
| `its_windows` | nhmmer (HMMER 3.4) with the 5.8S profiles of ITSx 1.1.3, then `its_windows.py` | `contigs_final.fasta` | `its_windows.fasta`, `its_windows.tsv` |
| `its_extraction` | ITSx 1.1.3 | `its_windows.fasta` | `{sample}.ITS1.fasta`, `{sample}.ITS2.fasta`, `{sample}.full.fasta`, `{sample}.summary.txt` |
| `its_classification` | `its_distinct.py`, then VSEARCH 2.30.0 `--sintax` | the three ITSx FASTA files and the UNITE reference | `{sample}_ITS_copies.txt`, `{sample}_ITS1_taxonomy.txt`, `{sample}_ITS2_taxonomy.txt` |
| `unite_db` or `unite_db_local`, once per run | curl and gzip, or a symbolic link | `links.unite_its_link` or `directories.unite_db` | `00.resources/unite/unite_its_sintax.fasta` and `PROVENANCE.txt` |
| `unite_identified`, once per run | `unite_identified.py` | `unite_its_sintax.fasta` | `00.resources/unite/unite_its_sintax.identified.fasta` and `identified_filter.txt` |

`its_windows.py`, `unite_identified.py` and `its_distinct.py` use only the Python standard
library. Each runs with the pinned Python of its rule's conda environment, Python 3.12.13 in itsx
and Python 3.12.15 in vsearch, so these steps need no `python` command on the machine.

## Why the genome is cut into windows

ITSx uses profile HMMs to find the end of the SSU, both ends of the 5.8S and the start of the
LSU, and cuts ITS1 and ITS2 out between them. The rDNA of a fungus is an array of tandem rDNA
units, each SSU-ITS1-5.8S-ITS2-LSU. A long-read assembly can place the whole array on a contig of
several megabases, often on the minus strand. Given such a contig whole, ITSx reports the domains
in the wrong order ("sequence may be chimeric") and writes nothing.

FunFlux therefore gives ITSx one short window per rDNA unit, turned so that it reads from SSU to
LSU. The cause and the design are explained in
[ITS on long-read contigs](../methods/its-long-read-contigs.md).

The windows are cut from `contigs_final.fasta` in every mode. In `contigs` mode that is the
filtered input FASTA. When the contamination and mitochondrion screen runs, the contigs removed
by the decontamination policy are no longer in it. With the default
`parameters.mitochondria.mode: separate_confirmed`, the mitochondrion is not in it either.

## Step 1: one window per rDNA copy

Rule `its_windows` finds the 5.8S genes of `contigs_final.fasta` and cuts a window around each,
with two commands:

```bash
nhmmer --cpu {threads} --dna -E 0.001 --tblout <hits.tbl> <5.8S profiles> contigs_final.fasta
python its_windows.py --tblout <hits.tbl> --fasta contigs_final.fasta --flank 2000 \
  --out-fasta its_windows.fasta --out-table its_windows.tsv
```

1. From ITSx's fungal profile file, `ITSx_db/HMMs/F.hmm` in the itsx environment, the rule takes
   every profile with `58S` in its NAME line. In ITSx 1.1.3 these are 358 of its 538 profiles,
   anchored at the start or the end of the 5.8S gene. The count is written to
   `logs/its_windows_{sample}.log`. nhmmer then searches `contigs_final.fasta` with these
   profiles.
2. `its_windows.py` keeps the hits with an E-value of at most 1e-5. The start and end profiles
   hit each 5.8S gene many times, so the script merges the hits into sites, one site per gene.
   On each contig and strand, a hit joins the current site when it starts at most 1,000 bp after
   the end of that site. Around the middle of each site the script cuts
   `parameters.its.window_flank` bases on each side (default 2,000), clipped at the contig ends,
   so a full window is 4,001 bp. A window on the minus strand is reverse-complemented.

The merge distance and the E-value limit are fixed, and `its.window_flank`, a whole number of at
least 200, is the only setting. The default 2,000 bp on each side hold the SSU end and the LSU
start that ITSx needs, around an ITS1-5.8S-ITS2 region of 0.4 to 1.5 kb.

`its_windows.fasta` names each window `contig:start-end_strand`, with 1-based positions on the
contig as it is stored in `contigs_final.fasta`, for example `ptg000001l:10001-14001_-`.
`its_windows.tsv` has one row per window:

| Column | Content |
|---|---|
| `window_id` | the FASTA header |
| `contig`, `start`, `end`, `strand` | where the window lies on the contig |
| `anchor_hits` | how many nhmmer hits were merged into the site |
| `best_evalue` | the smallest E-value among them |
| `window_length` | bases; 4,001 unless clipped at a contig end |

The table thus has one row per 5.8S gene found, and the [contig summary](contig-summary.md)
counts these rows per contig as rDNA units. A genome without any 5.8S hit gives an empty FASTA and
a table with only the header line, and the rule succeeds. The rule stops if a contig named in the
nhmmer table is missing from the FASTA or has another length there.

Two problems with ITSx's profile file stop the rule before the search, each with a message in the
log and on the console:

- a missing or empty `F.hmm`:
  `[FunFlux] its_windows: ITSx's fungal profile file <env>/bin/ITSx_db/HMMs/F.hmm is missing or empty. Is the job running in its itsx conda environment (--sdm conda), with the pinned ITSx 1.1.3?`
- no profile with `58S` in its NAME line:
  `[FunFlux] its_windows: no 5.8S profile found in <env>/bin/ITSx_db/HMMs/F.hmm (no NAME line contains 58S). Does this ITSx version name its 5.8S profiles differently?`

## Step 2: ITSx on the windows

Rule `its_extraction` runs ITSx on the windows:

```bash
ITSx -i its_windows.fasta -o 03.taxonomy/{sample}/{sample} --temp <scratch folder> \
  --preserve T --only_full F --graphical F --complement T --positions F --not_found F \
  -E 0.001 --nhmmer T --cpu {threads}
```

The options have the following effects:

| Option | Effect in ITSx 1.1.3 |
|---|---|
| `-E 0.001` | domain E-value cutoff (ITSx default 1e-5) |
| `--nhmmer T` | search with nhmmer instead of hmmsearch |
| `--complement T` | search both strands (ITSx default); the windows are already oriented |
| `--preserve T` | keep the window names as record names |
| `--only_full F` | report partial regions as well |
| `--graphical F`, `--positions F`, `--not_found F` | no graphical output, no table of positions, no list of sequences without ITS |

ITSx writes four files:

| File | Content |
|---|---|
| `{sample}.ITS1.fasta` | one ITS1 record per rDNA unit |
| `{sample}.ITS2.fasta` | one ITS2 record per rDNA unit |
| `{sample}.full.fasta` | ITS1, 5.8S and ITS2 together, one record per rDNA unit |
| `{sample}.summary.txt` | ITSx's own report: windows read, ITS detected on the main or the complementary strand, sequences judged chimeric |

The scratch folder given with `--temp` is deleted when the rule ends. When `its_windows.fasta` is
empty because no 5.8S gene was found, ITSx is not run. The four files are then written empty, and
`logs/its_extraction_{sample}.log` says so.

## Step 3: copies, distinct sequences and length limits

Rule `its_classification` first runs `its_distinct.py`, which reads the three ITSx files and does
three things for each region:

- It groups identical sequences and counts the rDNA units behind each. Upper and lower case, and
  line wrapping, are ignored.
- It checks the length against the limits below. A sequence outside them is listed as
  `rejected_too_short` or `rejected_too_long` and not classified.
- It names the accepted sequences `{sample}_ITS1_<n>`, `{sample}_ITS2_<n>` and
  `{sample}_full_<n>`, numbered in order of first appearance.

| Region | Accepted length | Classified |
|---|---|---|
| ITS1 | 50 to 600 bp (fixed) | yes |
| ITS2 | 50 to 600 bp (fixed) | yes |
| full (ITS1, 5.8S and ITS2) | at most `parameters.its.max_full_length` bp (default 1,500, a whole number of at least 400); no lower limit | no: it appears in the copies table only |

!!! note "The length limits keep a false ITS2 out of the classification"

    ITSx can join an SSU end and an LSU start that lie far apart, with only the end of a 5.8S
    between them, and so report a far too long "ITS2", to which SINTAX would still give a name. ITSx itself calls an ITS region
    suspiciously long above 1,500 bp, the default upper limit of the full region.

## Step 4: SINTAX against UNITE

Rule `its_classification` then classifies the distinct ITS2 sequences against UNITE:

```bash
vsearch --sintax <distinct ITS2 sequences> --db 00.resources/unite/unite_its_sintax.identified.fasta \
  --tabbedout 03.taxonomy/{sample}/{sample}_ITS2_taxonomy.txt \
  --sintax_cutoff 0.8 --randseed 1 --threads 1
```

The same command classifies ITS1. `--db` is the UNITE reference after
[the UNITE reference filter](#the-unite-reference-filter) described below. With
`parameters.its.unite_identified_only: false`, it is `00.resources/unite/unite_its_sintax.fasta`,
the whole release. `--sintax_cutoff` is `parameters.its_taxonomy_cutoff` (default 0.8), a number
above 0 and at most 1.

When no ITS1 or no ITS2 sequence is left, vsearch is not run for that region. The taxonomy file is
then created empty, and `logs/its_classification_{sample}.log` says
`No ITS1 sequence to classify for sample ... Skipping classification.`

SINTAX computes the confidence of each rank by classifying the query many times, each time from a
random subset of its k-mers (short words of the sequence), and reports how often the same name
came out. Without `--randseed`, vsearch draws a new random seed in every run, so the same sequence
gets slightly different confidences each time. The vsearch manual states that a fixed seed gives
the same result only with `--threads 1`, so FunFlux classifies on one thread.

## Reading the results

### The copies table

`{sample}_ITS_copies.txt` is tab-separated, with one row per distinct sequence: ITS1 rows first,
then ITS2, then full.

| Column | Content |
|---|---|
| `region` | `ITS1`, `ITS2` or `full` |
| `sequence_id` | the name of the distinct sequence, for example `SAMPLE1_ITS2_1` |
| `length` | bases |
| `copies` | how many rDNA units (windows) carry exactly this sequence |
| `status` | `ok`, `rejected_too_short` or `rejected_too_long` |
| `windows` | the names of those windows, separated by commas |

The number of rows per region and their `copies` show how alike the rDNA units of the assembly
are. Two or more rows for one region mean that the units differ, and each distinct sequence is
classified on its own. Rows with different classifications can point to a second organism,
because windows are cut from every contig of `contigs_final.fasta`. That file keeps contaminant
contigs unless the screen removed them, as described on [Decontamination](decontamination.md).

!!! note "`copies` counts the rDNA units of the assembly, not those of the genome"

    The `copies` column counts the rDNA units in the assembly that carry exactly this sequence.
    A short-read assembly collapses the array into one unit, while a long-read assembly can
    resolve several.

A `rejected_*` row has `-` as `sequence_id` and is not classified. A line starting with `#` after
the header means that no ITS was found, either because there was no 5.8S hit or because ITSx
found nothing in the windows.

### The taxonomy files

`{sample}_ITS1_taxonomy.txt` and `{sample}_ITS2_taxonomy.txt` are tab-separated, without a
header, with one line per classified sequence:

| Column | Content |
|---|---|
| 1 | the `sequence_id` of the copies table |
| 2 | the predicted taxonomy from the top rank (`d:`) down to `p:` phylum, `c:` class, `o:` order, `f:` family, `g:` genus, `s:` species, each rank with its confidence in brackets |
| 3 | the strand of the match |
| 4 | the ranks of column 2 with a confidence of at least `its_taxonomy_cutoff`, without the confidences |

Column 4 is the call at the cutoff, and column 2 shows how far the lower ranks got. When the
phylum is below the cutoff, for example, column 4 stops at `d:Fungi`, although column 2 may go on
to a genus and a species. A different `its_taxonomy_cutoff` changes column 4 only.

ITS1 and ITS2 are classified separately and can reach different depths, so read both files.

!!! note "A SINTAX confidence is not a sequence identity"

    A confidence measures how consistently the k-mer subsets point to one name. A confidence of
    1.00 at species level does not mean that the ITS equals a reference sequence.

The species name comes from UNITE and does not reach the annotation. The GenBank files and the
funannotate file names carry the species of the
[annotation parameter table](../reference/configuration.md#annotation-parameter-table), which can
differ from the ITS call.

## The UNITE reference filter

By default, FunFlux classifies against the UNITE records that name a genus, because the other
records can stop the call at `d:Fungi`. Every UNITE record carries its taxonomy in the header,
rank by rank:

```text
>LC784324|SH0155860.10FU;tax=d:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales,f:Sclerotiniaceae,g:Botrytis,s:Botrytis_cinerea;
```

Many records are identified only down to some rank above the genus, and their taxonomy then stops
there: `tax=d:Fungi;` alone, for example, or a string that ends with the family. In the default
release, about a third of the entries have no genus. SINTAX counts how often each name comes out
over its repeated classifications ([Step 4](#step-4-sintax-against-unite)). A record identified
only to the kingdom that is identical to the query is therefore followed as often as the
identified records. The result then stops at `d:Fungi`, even for a genome whose ITS matches
well-identified species perfectly. Classified against the records that carry a genus, the same
ITS gets its genus and species with the usual confidence values.

The filter is on by default, and the unfiltered result is one setting away:

| `parameters.its.unite_identified_only` | SINTAX reference | An ITS that matches well-identified species |
|---|---|---|
| `true` (default) | `00.resources/unite/unite_its_sintax.identified.fasta`, the records with a genus | gets its genus and species with the usual confidence values |
| `false` | `00.resources/unite/unite_its_sintax.fasta`, the whole release; rule `unite_identified` does not run | can stop at `d:Fungi` |

Once per run, rule `unite_identified` runs `unite_identified.py` on
`00.resources/unite/unite_its_sintax.fasta`. It takes a second or two and does three things:

- It keeps the records whose `tax=` field holds a genus rank, `g:` followed by a name. It drops
  the others, whose taxonomy stops above the genus or whose header has no `tax=` field.
- It writes the kept records, unchanged, to `00.resources/unite/unite_its_sintax.identified.fasta`,
  and the two counts, `kept_with_genus: <n>` and `dropped_without_genus: <n>`, to
  `00.resources/unite/identified_filter.txt`. `logs/unite_identified.log` holds a one-line summary
  and the same counts.
- It stops the run if no record at all carries a genus, which means the file is not a UNITE
  SINTAX release.

## The UNITE reference

FunFlux either downloads the UNITE reference or links a copy you already hold. A local copy in
`directories.unite_db` takes precedence over the download link.

- **Download.** While `directories.unite_db` is empty, rule `unite_db` fetches the gzip file
  named by `links.unite_its_link`. `curl -L --fail` downloads it and `gzip -dc` unpacks it on
  the fly into `00.resources/unite/unite_its_sintax.fasta`, so the `.gz` file is not kept. No
  checksum is checked. `PROVENANCE.txt` records the URL, the UTC date, and the md5 and size of
  the FASTA.
- **Local copy.** When `directories.unite_db` names a folder, rule `unite_db_local` runs instead.
  The folder must hold an uncompressed `unite_its_sintax.fasta`, otherwise the run stops before
  any job starts. `00.resources/unite/unite_its_sintax.fasta` becomes a symbolic link to it, so
  nothing is copied, and the folder is only read. `PROVENANCE.txt` records the folder, the UTC
  date, md5 and size. A NOTE says that `links.unite_its_link` is not used.

Both rules log to `logs/unite_db.log`. The files in `00.resources/unite/` are made once per output
folder, and later runs on the same output folder reuse them. Whichever rule made the FASTA,
`unite_identified` then filters it as described above.

The table compares the two UNITE releases of version 19.02.2025 (UNITE 10.0).

| Release | DOI | Sequences | Size | In FunFlux |
|---|---|---|---|---|
| USEARCH/UTAX/SINTAX release for eukaryotes (licence CC BY 4.0) | 10.15156/BIO/3301246 | 266,589, from all eukaryotes, not only fungi | 35.7 MB to download, 207 MB unpacked | the default link; FunFlux v2 was tested with it |
| Fungi-only release | 10.15156/BIO/3301245 | 168,030 | | |

FunFlux v2 was tested with the eukaryote release, but another release in SINTAX format can be
given through either key. How to download a local copy is described on
[Reference databases](../getting-started/databases.md).

## Threads and run time

Each of the three steps takes seconds to minutes. Their thread caps are low, so that they can run
beside InterProScan, which holds 10 threads for hours:

| Rule | Threads (at most) |
|---|---|
| `its_windows` | 8 |
| `its_extraction` | 8 |
| `its_classification` | 1 |

nhmmer and ITSx run with `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`, so they stay within
their thread count. Without these two settings, HMMER's OpenBLAS library starts one thread per CPU
core whatever `--cpu` says, which makes hundreds of threads on a large server.

## Limits

- A fungus whose 5.8S matches none of the 358 profiles gives no window and empty ITS files, and
  the copies table then holds the `#` line.
- The window names use the contig names of `contigs_final.fasta`. The annotation renames the
  contigs `scaffold_N`, as described under [funannotate sort](repeat-masking.md#funannotate-sort),
  and the `annotated_as` column of the [contig summary](contig-summary.md) links the two names.
