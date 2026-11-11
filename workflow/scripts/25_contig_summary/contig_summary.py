#!/usr/bin/env python3
"""Say what every contig of the draft assembly is: one table per sample, with reasons.

BIOLOGY AND WHY THIS SCRIPT EXISTS
----------------------------------
A finished fungal assembly is a list of contigs, and the list alone does not say what
each contig is. A long-read assembly of a filamentous fungus typically holds, side by
side: whole chromosomes (one contig from telomere to telomere), chromosome arms (the
assembler stopped at a repeat, often the rDNA array or the centromere), the
mitochondrial genome (a small circle present in many copies per cell), pieces of the
ribosomal DNA (rDNA) array (tens to hundreds of tandem copies of the
SSU-ITS1-5.8S-ITS2-LSU unit, which no assembler resolves completely), alternative
copies of a region (a second haplotype, or a stretch assembled twice), low-depth
debris, and sometimes a contaminant that the contamination screen removed. A
short-read assembly adds two more kinds, by the hundred: collapsed repeats (several
genome copies of a transposon or a satellite folded into one contig, so that the
reads of all the copies pile up on it) and short pieces of the nuclear genome with
no gene on them. FunFlux
already produces the evidence for all of this, spread over half a dozen files. This
script gathers it into one row per contig and gives every contig a verdict with its
reasons, so that a reader sees at a glance how complete the assembly is and which
contigs deserve a closer look.

The evidence, and what each piece says:
  * telomere repeats at the contig ends, called by tidk with the repeat unit that the
    telomere step chose for the sample: a telomere at both ends means a whole
    chromosome, at one end a chromosome arm;
  * read depth (the mean depth that samtools coverage measures from the sample's reads
    mapped to the draft): the nuclear contigs share one depth; the mitochondrion is
    several times deeper; an alternative copy of a region, or noise, is much shallower;
  * the assembler's own depth and its circularity flag (long-read modes);
  * rDNA units: the ITS windows that the taxonomy step cut, one per 5.8S gene found;
  * complete BUSCO markers and predicted protein-coding genes, which tell genuine
    nuclear sequence from rDNA and debris;
  * the verdicts of the two screens, the mitochondrion screen and the genus-based
    contamination screen.

WHERE THIS SCRIPT FITS IN THE WORKFLOW
--------------------------------------
    1. The screens, the assembly QC, the ITS step, the telomere step (rule
       telomere_units, script telomere_units.py) and funannotate all run first and
       write the files listed under INPUTS (the constants in capitals are the names
       used in workflow/rules/shared/00_common.smk).
    2. THIS SCRIPT is run once per sample by rule contig_summary
       (workflow/rules/shared/25_contig_summary.smk), in every mode, and writes the
       two OUTPUTS into the sample's eval folder. The draft FASTA and the two files
       of the telomere step are required: every other input is optional, and a mode
       that cannot produce a piece of evidence (no reads in contigs mode, say) simply
       leaves its columns as NA.
    3. Rule multiqc reads the *_mqc.tsv copy, so the table appears in the report;
       the table itself is a record for the reader. No rule reads it.

WHAT THIS SCRIPT DOES
---------------------
    1. Reads the draft assembly (DRAFT_CONTIGS): the kept contigs, the removed ones
       and the mitochondrion alike, so that nothing disappears silently.
    2. Reads the telomere unit of the sample and tidk's telomere calls at every
       contig end, both written by telomere_units.py (see TELOMERES below).
    3. Reads whichever optional tables were given, and joins them to the contigs by
       contig name. The gene prediction is the one input that cannot be joined by
       name, because funannotate renames the contigs scaffold_1, scaffold_2, ... in
       order of length; it is joined by sequence instead (see map_contigs_to_scaffolds).
    4. Works out the median read depth and every contig's depth relative to it.
    5. Applies the verdict rules (see RULES) to every contig.
    6. Writes the table and its MultiQC copy.

TELOMERES: THE UNIT AND THE CALLS COME FROM telomere_units.py
-------------------------------------------------------------
Telomeres are tandem arrays of a short repeat that cap every chromosome end. The
repeat unit differs between fungal groups: TTAGGG in most filamentous fungi, other
units in some groups (the Aspergillus flavus group, Malassezia, Yarrowia, ...), and
long, irregular repeats in most budding yeasts, which no single unit describes (Lue
2021; Sun et al. 2023). This script does not look for the unit itself. Rule
telomere_units runs telomere_units.py first: with tidk (Brown et al. 2025) it lists
the tandem repeat units found at the contig ends, scores every candidate unit (the
units found, the units of tidk's telomere database, the units of the literature and
the configured unit, if any) by how it behaves at the contig ends, takes the
best-supported one as the telomere unit of the sample, and runs tidk ends with that
unit over every contig. telomere_units.py and docs/analysis/contig-summary.md hold
the literature table and the choice rules. This script reads two of its files:

1. telomere_unit.txt (--telomere-unit-file): one line, the unit in tidk's canonical
   form (the alphabetically smallest of all rotations of the unit and of its reverse
   complement, so that TTAGGG, GGGTTA and CCCTAA are all written AACCCT), or the
   word none when no unit was eligible. It fills the telomere_unit column, the same
   on every row of the sample (the tables of several samples merge into one in the
   MultiQC report, and two species may use two units), and it is named in the
   MultiQC description.
2. telomere_ends.tsv (--telomere-ends): tidk's table for that unit, one row per
   contig end, with tidk's columns id, sequence_length, side, status,
   telomere_start, telomere_end, telomere_length, offset, strand, window_limited (a
   header line alone when the unit is none). side is start or end; telomere_length
   is the length of the array in bp and offset its distance from the contig end in
   bp. status says what tidk found within the telomere window of that end:
     present       an array of at least min_telomere_copies copies, on the strand a
                   telomere has at that end (on the forward strand of a contig the
                   array of most filamentous fungi reads CCCTAA at the start and
                   TTAGGG at the end). This is a telomere call.
     wrong_strand  an array on the other strand. A telomere repeat points one way
                   at every chromosome end, so this is NOT a telomere call; the
                   reason column still names it ("array in the opposite orientation
                   at the end"), because the contig does carry the repeat.
     other_end     one array seen from both ends of a short contig; telomere_units.py
                   gives this status to the row with the larger offset, so that the
                   array is counted once. Nothing is said about it.
     not_terminal  an array further inside the contig than the window allows.
                   Nothing is said about it.
     absent        no array. A contig without any row in the table has no telomere
                   at either end.

From this, per contig end: telomere_start and telomere_end hold the copies of the
unit in a present array (its length in bp divided by the unit length, rounded to the
nearest whole number; 0 for every other status), and a telomere is called from
min_telomere_copies copies (5 by default; telomere_units.py asks tidk for arrays of
at least that many copies, so a present row always passes when both steps use the
same value, as they do in the workflow). The two offset columns,
telomere_end_offset_start and telomere_end_offset_end, carry the offset of a
present array (NA otherwise). A telomere is the very end of a chromosome, so an
array that sits far inside the contig means sequence beyond the telomere: a misjoin,
or an interstitial telomeric repeat. An array more than FAR_FROM_END_BP (1,000) bases
from the end is therefore also named in the reason column. The verdict does not
change; the reader is told.

RULES: THE VERDICT, IN THIS ORDER (THE FIRST RULE THAT APPLIES DECIDES)
----------------------------------------------------------------------
Evidence a mode cannot produce is "unknown" (NA in the table), and a rule simply
does not test it: "depth in 0.5-1.5 when known" passes when the depth is unknown.
   1. mitochondrion     the mitochondrion screen confirmed the contig (classification
                        confirmed_mitochondrial), or the contig is in the
                        mitochondrial FASTA that the split step wrote.
   2. contaminant       the contamination screen removed it (action remove); the
                        reason names the genus.
   3. foreign_mitochondrion
                        the mitochondrion screen found the contig to be mitochondrial
                        DNA of an animal or a plant, not of the fungus (classification
                        foreign_mitochondrial): human DNA from the library, for example.
                        The contamination screen kept it (with the default,
                        decontamination mode off, it keeps every contig), so it stays in
                        the assembly, and the verdict says what it is.
   4. rDNA_array        at least min_rdna_units rDNA units (3 by default), and fewer
                        protein-coding genes than rDNA units (or the gene count is
                        unknown). A pure rDNA array holds no protein-coding gene, but
                        the gene predictor now and then calls one inside an rDNA
                        unit, so a strict zero would misname a real rDNA array that
                        carries a stray gene call; a contig with more rDNA units
                        than genes is mostly rDNA.
   5. chromosome        telomere repeats at both ends, depth 0.5-1.5 times the
                        median when known, at least one gene when known. Size plays
                        no part: small accessory chromosomes qualify.
   6. chromosome_arm    telomere repeats at one end, the rest as in 5.
   7. low_depth         depth below 0.5 times the median: an alternative copy of a
                        region (a second haplotype, a stretch assembled twice) or
                        noise.
   8. repeat            no BUSCO marker, at most two predicted genes, depth above 1.5
                        times the median: a collapsed repeat, that is, several genome
                        copies of a transposon, a satellite or another repeat family
                        that the assembler folded into one contig, so that the reads
                        of all the copies pile up on it. A gene call or two is allowed
                        because predictors call open reading frames inside
                        transposons. When the depth is unknown this rule does not fire.
   9. nuclear_contig    at least one gene or BUSCO marker, depth 0.5-1.5 when known,
                        no telomere repeat: a chromosome fragment.
  10. nuclear_fragment  no gene and no BUSCO marker, no telomere repeat, depth
                        0.5-1.5 or unknown: a piece of the nuclear genome too short
                        to carry a gene. A short-read assembly ends in many of them.
  11. unclassified      anything else; the reason says what was and was not seen.
"No gene and no BUSCO marker" (rules 8 and 10) means that every known count is 0 and
at least one of the two counts is known. In the workflow both are always counted
(BUSCO and funannotate run in every mode), and only a contig without a scaffold in
the annotated genome lacks its gene count. Run by hand without either file, the two
rules do not fire: nothing then says that the contig lacks genes, and it stays
unclassified.
The median depth is weighted by length and taken over the contigs of at least 100 kb
(all contigs when none is that long), leaving out the mitochondrion and the removed
contigs, so that short debris, the deep mitochondrion and a contaminant genome
cannot pull it.

INPUTS
------
    --contigs            DRAFT_CONTIGS, the draft assembly (FASTA). Required; the
                         table has one row per record, in file order. Names are the
                         first word of the header, as every other table uses them.
    --telomere-unit-file TELOMERE_UNIT_FILE, telomere_unit.txt of telomere_units.py:
                         the telomere repeat unit of the sample, or none. Required.
    --telomere-ends      TELOMERE_ENDS, telomere_ends.tsv of telomere_units.py: tidk's
                         table of what it found at every contig end (see TELOMERES).
                         Required; a row for a contig that is not in the draft stops
                         the script, a contig without rows has no telomere.
    --mito-contigs       MITO_CONTIGS, the mitochondrial contig(s) the split step wrote
                         (FASTA, possibly empty). Only the names are read.
    --decisions          CONTIG_DECISIONS, the contamination screen's audit table:
                         columns contig, assigned_genus, action (keep/remove), reason.
    --mito-screen        MITO_SCREEN_TSV, the mitochondrion screen's report: columns
                         contig, classification, reason (and more, not read).
    --depth-table        CONTIG_DEPTH, the samtools coverage table of rule
                         contig_depth: a header line starting with "#rname", then one
                         line per contig; the read depth is the column meandepth.
    --assembly-evidence  ASSEMBLY_EVIDENCE (long-read modes): columns contig, source,
                         circular (Y/N), coverage, repeat, alt_group; only contig,
                         circular and coverage are read.
    --its-windows        ITS_WINDOW_TABLE, its_windows.tsv: one row per rDNA unit found,
                         with a contig column.
    --busco-dir          BUSCO_DIR: the folder BUSCO wrote. The script reads the
                         full_table.tsv of the run_* folder with the largest marker
                         set (the specific lineage, in auto mode) and counts, per
                         contig, the markers with status Complete or Duplicated.
    --gff3               the funannotate GFF3 (annotate_results/<species>.gff3): the
                         protein-coding genes are its mRNA features, counted per
                         sequence name (scaffold_N).
    --annotated-fasta    CONTIGS_SORT, the genome as funannotate annotated it, with the
                         scaffold_N names. Needed to join the GFF3 to the draft contigs.
    --adapter-screen     ADAPTER_SCREEN, adapter_screen.tsv of the contig adapter screen
                         (adapter_screen.py, same folder): one row per adapter locus, with
                         the columns contig and position (end or internal). It fills the
                         adapter column and does not change any verdict.
    --min-telomere-copies, --min-rdna-units
                         the two tuning options (see TELOMERES and RULES); the
                         workflow passes the values of the configuration block
                         parameters.contig_summary (keys min_telomere_copies and
                         min_rdna_units).
    --sample             the sample name, for the MultiQC row names and the messages.

OUTPUTS
-------
    --out          CONTIG_SUMMARY, contig_summary.tsv: tab-separated, one header line,
                   one row per draft contig, the columns of COLUMNS below.
    --out-mqc      CONTIG_SUMMARY_MQC, contig_summary_mqc.tsv: the same table behind
                   MultiQC's custom-content header lines, with a first column that
                   holds the row name "<sample>: <contig>", so that the tables of
                   several samples merge into one in the report.
    standard error one summary line (the counts per verdict, the telomere unit and the
                   ends with a telomere, the BUSCO lineage used). Rule contig_summary
                   writes the stream to its log.

Exit status 0 when the table was written. A missing or unreadable file, a table
without the columns named above, a depth or count that is not a number, a FASTA
with two records of the same name, a GFF3 line with fewer than nine columns, or a
telomere table that names a contig the draft does not hold stops the script with
one readable message and exit status 1, and no output is written.

REFERENCES
----------
The three papers behind the telomere paragraphs above. docs/about/citation.md lists
them too, with the papers behind the telomere units of the fungal groups, which
telomere_units.py carries.
Lue NF (2021). Duplex telomere-binding proteins in fungi with canonical telomere
    repeats: new lessons in the rapid evolution of telomere proteins. Frontiers in
    Genetics 12: 638790. https://doi.org/10.3389/fgene.2021.638790
Sun Q, Wang H, Tao S, Xi X (2023). Large-scale detection of telomeric motif sequences
    in genomic data using TelFinder. Microbiology Spectrum 11(2): e03928-22.
    https://doi.org/10.1128/spectrum.03928-22
Brown MR, Gonzalez de la Rosa PM, Blaxter M (2025). tidk: a toolkit to rapidly
    identify telomeric repeats from genomic datasets. Bioinformatics 41(2): btaf049.
    https://doi.org/10.1093/bioinformatics/btaf049
"""

# Only Python's standard library is used, so the Python 3.12 that the vsearch conda
# environment holds (the one the ITS scripts run with) is all this script needs.
import argparse
import csv
import glob
import os
import re
import sys
from collections import OrderedDict


# ======================================================================
# SETTINGS
# ======================================================================

# Defaults of the two tuning options. The workflow passes the values of its
# contig_summary configuration block; a run by hand without the options uses these.
#   min telomere copies: copies of the telomere unit from which a telomere is called.
#                        telomere_units.py asks tidk for arrays of at least this many
#                        copies, so the two steps are meant to use the same value
#                        (the workflow passes one setting to both).
#   min rDNA units:      ITS windows from which a contig can be an rDNA array.
DEFAULT_MIN_TELOMERE_COPIES = 5
DEFAULT_MIN_RDNA_UNITS = 3

# The limits of the two options. They are the same limits the workflow applies when
# it reads the configuration, so that a run by hand is held to the same rules.
MIN_TELOMERE_COPIES_MIN = 2
MIN_RDNA_UNITS_MIN = 1

# The telomere table (see TELOMERES at the top): the two sides of a contig, and the
# statuses a row may carry (tidk's four, plus other_end from telomere_units.py). Any
# other text in the status column means that the wrong file was given, and the
# script stops.
TELOMERE_SIDES = ("start", "end")
TELOMERE_STATUSES = {"present", "absent", "wrong_strand", "not_terminal", "other_end"}

# A telomere array that begins more than this many bases inside the contig is named
# in the reason column ("the start array is 2814 bp from the end"): a telomere is the
# very end of a chromosome, so sequence beyond it points at a misjoin or an
# interstitial telomeric repeat. The verdict does not change; the reader is told.
FAR_FROM_END_BP = 1000

# Depth. The median is taken over the contigs of at least this length, so that short
# debris cannot pull it, and a contig is "at normal depth" between these two
# multiples of the median: below is an alternative copy or noise, above a collapsed
# repeat or an organelle.
MEDIAN_MIN_CONTIG_LENGTH = 100000
NORMAL_DEPTH_LOW = 0.5
NORMAL_DEPTH_HIGH = 1.5
# A collapsed repeat may still carry a gene call or two: gene predictors call open reading
# frames inside transposons. Up to this many predicted genes do not keep a deep, marker-less
# contig from being called a repeat.
REPEAT_MAX_GENES = 2

# The join of the draft contigs to funannotate's scaffold_N names is by identical
# sequence. When a contig was polished after the draft (Medaka for Nanopore samples,
# Polypolish in hybrid mode), the two sequences differ by a base here and there and
# the identity test fails; the fallback then takes PROBE_COUNT pieces of PROBE_LENGTH
# bases from the inside of the contig and looks for them in the scaffolds whose
# length is within PROBE_LENGTH_TOLERANCE of the contig's (polishing changes a length
# by far less). See map_contigs_to_scaffolds().
PROBE_LENGTH = 300
PROBE_COUNT = 5
PROBE_LENGTH_TOLERANCE = 0.05

# What an unknown value reads in the table.
NA = "NA"

# The columns of contig_summary.tsv, in this order.
COLUMNS = [
    "contig", "length", "gc",
    "depth", "depth_x_median", "assembler_depth", "circular",
    "telomere_start", "telomere_end",
    "telomere_end_offset_start", "telomere_end_offset_end",
    "telomere_unit",
    "rdna_units", "busco_markers", "genes",
    "genus", "screen", "adapter", "annotated_as",
    "verdict", "reason",
]


# ======================================================================
# STEP 1: READ A FASTA FILE
# ======================================================================

def read_fasta(path):
    """Read a FASTA file into {contig name: sequence in upper case}, in file order.

    Takes in:  the path of a FASTA file: the draft assembly (--contigs) or the
               annotated genome (--annotated-fasta). Sequences may be wrapped over
               several lines.
    Does:      reads the file line by line. The contig name is the first word of the
               ">" header, which is the name every table of the workflow uses (BLAST,
               BlobTools, the screens, BUSCO, the ITS table). Sequence lines are
               joined and put in upper case, so that soft-masked (lower-case) bases
               do not upset the identity test of map_contigs_to_scaffolds(), which
               must compare like with like.
    Produces:  an ordered dictionary {name: sequence}.
    Used by:   run(), for the draft and for the annotated genome.
    Stops with a clear error when two records share a name (the join to the tables
    would be ambiguous), when a sequence line comes before any header (not a FASTA
    file), or when the file holds no record at all.
    """
    sequences = OrderedDict()
    current_name = None
    current_lines = []

    def store_current_record():
        # Join the sequence lines of the record just read and file it under its name.
        if current_name is None:
            return
        if current_name in sequences:
            raise ValueError(
                f"{path}: the contig name {current_name!r} appears more than once, so "
                "the tables cannot be joined to the contigs"
            )
        sequences[current_name] = "".join(current_lines).upper()

    with open(path, "r", encoding="utf-8", errors="replace") as fasta:
        for line_number, line in enumerate(fasta, start=1):
            text = line.strip()
            if text == "":
                continue
            if text.startswith(">"):
                store_current_record()
                words = text[1:].split()
                current_name = words[0] if len(words) > 0 else ""
                current_lines = []
            elif current_name is None:
                raise ValueError(
                    f"{path}, line {line_number}: found sequence before any FASTA header "
                    f"line (starting with >), so this is not a FASTA file: {text[:40]!r}"
                )
            else:
                current_lines.append(text)
    store_current_record()

    if len(sequences) == 0:
        raise ValueError(f"{path}: no FASTA record found")
    return sequences


def read_fasta_names(path):
    """Return the set of contig names in a FASTA file, without reading the sequences.

    Takes in:  --mito-contigs, the mitochondrial FASTA the split step wrote. It is
               empty when no contig was confirmed, which gives an empty set.
    Does:      collects the first word of every ">" header line.
    Produces:  a set of names.
    Used by:   run(); the names mark contigs as mitochondrion in decide().
    """
    names = set()
    with open(path, "r", encoding="utf-8", errors="replace") as fasta:
        for line in fasta:
            if line.startswith(">"):
                words = line[1:].split()
                if len(words) > 0:
                    names.add(words[0])
    return names


def gc_fraction(sequence):
    """Return the GC content of a sequence as a fraction of its A, C, G and T bases.

    Takes in:  one contig sequence in upper case.
    Does:      counts G and C and divides by the count of the four plain bases, so
               that runs of N (unknown bases) do not dilute the value.
    Produces:  a number between 0 and 1, or 0.0 for a sequence without any plain base.
    Used by:   run(), for the gc column.
    """
    g_and_c = sequence.count("G") + sequence.count("C")
    plain_bases = g_and_c + sequence.count("A") + sequence.count("T")
    if plain_bases == 0:
        return 0.0
    return g_and_c / plain_bases


# ======================================================================
# STEP 2: THE TELOMERE UNIT AND THE TELOMERE CALLS OF telomere_units.py
# ======================================================================

def read_telomere_unit(path):
    """Read telomere_unit.txt: the telomere repeat unit of the sample, or none.

    Takes in:  --telomere-unit-file, the one-line file telomere_units.py wrote: the
               unit in canonical form (AACCCT for TTAGGG), or the word none when no
               unit was eligible.
    Does:      reads the one line, puts it in upper case and checks that it is a unit
               made of the letters A, C, G and T, or none.
    Produces:  the unit as text, or None for none.
    Used by:   run(): the unit length turns tidk's array lengths into copies, and the
               unit fills the telomere_unit column and the MultiQC description.
    Stops with a clear error on an empty file, on more than one line, and on
    anything else than a unit or none (the wrong file was given).
    """
    with open(path, "r", encoding="utf-8", errors="replace") as unit_file:
        lines = [line.strip() for line in unit_file if line.strip() != ""]
    if len(lines) != 1:
        raise ValueError(
            f"{path}: expected one line holding the telomere repeat unit or none, "
            f"found {len(lines)} non-empty lines"
        )
    text = lines[0].upper()
    if text == "NONE":
        return None
    if not re.fullmatch(r"[ACGT]+", text):
        raise ValueError(
            f"{path}: the telomere repeat unit must be made of the letters A, C, G and T, "
            f"or be the word none (got: {lines[0]!r})"
        )
    return text


def read_telomere_ends(path, contig_names):
    """Read telomere_ends.tsv: what tidk found at every contig end, for the chosen unit.

    Takes in:  --telomere-ends, the table telomere_units.py wrote (tidk's columns id,
               sequence_length, side, status, telomere_start, telomere_end,
               telomere_length, offset, strand, window_limited; see TELOMERES at the
               top), and the names of the draft contigs.
    Does:      reads the columns id, side, status, telomere_length and offset of every
               row; the other columns are not needed. The array length and the offset
               are read as numbers for a present row; the other rows carry NA there.
    Produces:  {contig: {side: (status, telomere_length, offset)}}, side being start
               or end. A contig without a row is absent from the dictionary, which
               means no telomere at either end. telomere_length and offset are None
               unless the status is present.
    Used by:   run(), through telomere_columns().
    Stops with a clear error on a row for a contig that is not in the draft (the
    table belongs to another assembly), on a side or a status that tidk does not
    write, on two rows for the same contig end, and on a present row without its
    numbers.
    """
    ends_by_contig = {}
    rows = read_tsv_rows(path, {"id", "side", "status", "telomere_length", "offset"})
    for line_number, row in enumerate(rows, start=2):   # line 1 is the header
        contig = (row["id"] or "").strip()
        if contig == "":
            continue
        if contig not in contig_names:
            raise ValueError(
                f"{path}, line {line_number}: the contig {contig!r} is not in the draft "
                "assembly, so this telomere table belongs to another assembly"
            )
        side = (row["side"] or "").strip()
        if side not in TELOMERE_SIDES:
            raise ValueError(
                f"{path}, line {line_number}: the side must be start or end (got: {side!r})"
            )
        status = (row["status"] or "").strip()
        if status not in TELOMERE_STATUSES:
            raise ValueError(
                f"{path}, line {line_number}: the status {status!r} is not one of "
                f"{', '.join(sorted(TELOMERE_STATUSES))}"
            )
        contig_ends = ends_by_contig.setdefault(contig, {})
        if side in contig_ends:
            raise ValueError(
                f"{path}, line {line_number}: a second row for the {side} of {contig!r}"
            )

        # Only a present array has a length and an offset worth reading.
        telomere_length = None
        offset = None
        if status == "present":
            telomere_length = number_from_cell(row["telomere_length"], path, "telomere_length")
            offset = number_from_cell(row["offset"], path, "offset")
            if telomere_length is None or offset is None:
                raise ValueError(
                    f"{path}, line {line_number}: a present array needs its telomere_length "
                    f"and offset ({contig!r}, {side})"
                )
            telomere_length = int(telomere_length)
            offset = int(offset)
        contig_ends[side] = (status, telomere_length, offset)
    return ends_by_contig


def telomere_columns(contig_ends, unit):
    """Turn the telomere rows of one contig into the telomere columns of the table.

    Takes in:  the rows of one contig from read_telomere_ends(), {side: (status,
               telomere_length, offset)}, empty for a contig without rows; and the
               telomere unit (None when none was identified, in which case no contig
               has a row).
    Does:      for each side, a present array gives its copies of the unit (the array
               length in bp divided by the unit length, rounded to the nearest whole
               number) and its offset from the contig end; every other status, and a
               missing row, gives 0 copies and no offset. The status itself is kept
               for the reason column, which names an array in the opposite
               orientation (status wrong_strand).
    Produces:  a dictionary with, for the start of the contig, telomere_start (the
               copies), telomere_end_offset_start (bp, or None) and
               telomere_status_start (the status, or None without a row), and the
               same three with "end" for the end of the contig.
    Used by:   run(), once per contig; decide() and telomere_text() read the fields.
    """
    columns = {}
    for side in TELOMERE_SIDES:
        status, telomere_length, offset = contig_ends.get(side, (None, None, None))
        copies = 0
        if status == "present":
            copies = round(telomere_length / len(unit))
        columns["telomere_" + side] = copies
        columns["telomere_end_offset_" + side] = offset
        columns["telomere_status_" + side] = status
    return columns


# ======================================================================
# STEP 3: READ THE OPTIONAL TABLES
# ======================================================================

def read_tsv_rows(path, required_columns):
    """Read a tab-separated table with a header line, checking the columns we need.

    Takes in:  the path, and the set of column names the caller will read.
    Does:      reads the file with csv.DictReader (each row becomes a dictionary
               column name -> cell text) and stops when a required column is missing
               from the header: a table from another tool, or a cut-off file, must not
               silently give empty columns.
    Produces:  a list of row dictionaries.
    Used by:   the readers of the telomere ends, decisions, mitochondrion screen,
               assembly evidence, ITS window and adapter tables.
    """
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as table:
        reader = csv.DictReader(table, delimiter="\t")
        found = set(reader.fieldnames or [])
        missing = sorted(required_columns - found)
        if len(missing) > 0:
            raise ValueError(f"{path}: missing the column(s) {', '.join(missing)}")
        return list(reader)


def number_from_cell(text, path, column):
    """Turn the text of a table cell into a number, or stop with a message naming the column.

    Takes in:  the cell text, plus the file and column names for the message.
    Produces:  a float. An empty cell gives None (the value was not recorded).
    Used by:   read_telomere_ends(), read_depths() and read_assembly_evidence().
    """
    text = (text or "").strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{path}: the {column} value {text!r} is not a number") from None


def read_decisions(path):
    """Read the contamination screen's audit table.

    Takes in:  --decisions (contig_taxonomy_decisions.tsv): columns contig,
               assigned_genus, action, reason, one row per contig the screen judged.
    Produces:  {contig: (genus, action, reason)}; action is keep or remove.
    Used by:   run(); decide() reads the action, the genus and screen columns the rest.
    """
    decisions = {}
    for row in read_tsv_rows(path, {"contig", "assigned_genus", "action", "reason"}):
        contig = (row["contig"] or "").strip()
        if contig == "":
            continue
        decisions[contig] = (
            (row["assigned_genus"] or "").strip(),
            (row["action"] or "").strip(),
            (row["reason"] or "").strip(),
        )
    return decisions


def read_mito_screen(path):
    """Read the mitochondrion screen's report.

    Takes in:  --mito-screen (mitochondrial_screen.tsv): one row per draft contig
               with its classification (confirmed_mitochondrial,
               mitochondrial_candidate, foreign_mitochondrial (mitochondrial DNA of
               an animal or a plant, such as human), circular_unclassified,
               high_coverage_unclassified, not_mitochondrial or not_screened) and
               the reason code; the other columns are not read.
    Produces:  {contig: (classification, reason)}.
    Used by:   run(); decide() tests the classification, the screen column shows it.
    """
    screen = {}
    for row in read_tsv_rows(path, {"contig", "classification", "reason"}):
        contig = (row["contig"] or "").strip()
        if contig == "":
            continue
        screen[contig] = ((row["classification"] or "").strip(), (row["reason"] or "").strip())
    return screen


def read_depths(path):
    """Read the mean read depth of every contig from the samtools coverage table.

    Takes in:  --depth-table (CONTIG_DEPTH, {sample}_depth.tsv of rule contig_depth): the
               output of samtools coverage on the screen's alignment of the sample's reads
               to the draft. Its first line names the columns and starts with "#rname";
               every other line is one contig. The depth is the column meandepth: the
               read bases aligned to the contig, added up over all its positions and
               divided by its length (supplementary alignments counted, secondary ones
               left out).
    Does:      finds the columns rname and meandepth in the header line and reads one
               number per contig.
    Produces:  {contig: depth}.
    Used by:   run(), for the depth column and the median.
    Stops with a clear error when the header line lacks either column, when a depth is
    not a number, when a contig line is too short, or when the table holds no contig.
    """
    depths = {}
    name_column = None
    depth_column = None
    with open(path, "r", encoding="utf-8", errors="replace") as table:
        for line_number, line in enumerate(table, start=1):
            text = line.rstrip("\n")
            if text.strip() == "":
                continue
            if text.startswith("#"):
                # The header line: "#rname startpos endpos numreads ... meandepth ...".
                # The # is dropped so that the first column name reads rname.
                columns = text[1:].split("\t")
                if "rname" not in columns or "meandepth" not in columns:
                    raise ValueError(
                        f"{path}, line {line_number}: expected the samtools coverage "
                        "header with the columns rname and meandepth"
                    )
                name_column = columns.index("rname")
                depth_column = columns.index("meandepth")
                continue
            if depth_column is None:
                raise ValueError(f"{path}, line {line_number}: a contig line before the header line")
            fields = text.split("\t")
            # A contig line must reach at least as far as the later of the two columns.
            needed_columns = max(name_column, depth_column) + 1
            if len(fields) < needed_columns:
                raise ValueError(
                    f"{path}, line {line_number}: expected at least "
                    f"{needed_columns} columns but found {len(fields)}"
                )
            contig = fields[name_column].strip()
            depth = number_from_cell(fields[depth_column], path, "meandepth")
            if depth is None:
                raise ValueError(f"{path}, line {line_number}: the meandepth of {contig} is empty")
            depths[contig] = depth
    if len(depths) == 0:
        raise ValueError(f"{path}: no contig row found in the depth table")
    return depths


def read_assembly_evidence(path):
    """Read the assembler's own depth and circularity flag per contig.

    Takes in:  --assembly-evidence (assembly_evidence.tsv, long-read modes): the
               table mitochondria.py writes from Flye's assembly_info.txt or hifiasm's
               graph, with the columns contig, source, circular (Y/N), coverage,
               repeat, alt_group. Only contig, circular and coverage are read.
    Produces:  {contig: (coverage or None, circular True/False)}.
    Used by:   run(), for the assembler_depth and circular columns and the reasons.
    """
    evidence = {}
    for row in read_tsv_rows(path, {"contig", "circular", "coverage"}):
        contig = (row["contig"] or "").strip()
        if contig == "":
            continue
        coverage = number_from_cell(row["coverage"], path, "coverage")
        circular = (row["circular"] or "").strip().upper() in {"Y", "YES", "TRUE", "1"}
        evidence[contig] = (coverage, circular)
    return evidence


def read_its_window_counts(path):
    """Count the rDNA units found on every contig.

    Takes in:  --its-windows (its_windows.tsv): one row per window that the ITS step
               cut around a 5.8S gene, that is, one row per rDNA unit found; column
               contig names the contig the window lies on.
    Produces:  {contig: number of windows}.
    Used by:   run(), for the rdna_units column.
    """
    counts = {}
    for row in read_tsv_rows(path, {"contig"}):
        contig = (row["contig"] or "").strip()
        if contig == "":
            continue
        counts[contig] = counts.get(contig, 0) + 1
    return counts


def find_busco_table(busco_dir):
    """Find the full_table.tsv of the BUSCO run with the largest marker set.

    Takes in:  --busco-dir, the folder BUSCO wrote for the sample. The workflow runs
               BUSCO with --auto-lineage-euk by default: BUSCO then places the
               genome, writes the broad eukaryote run under auto_lineage/ and the
               specific run (hypocreales_odb12, say, with thousands of markers) as
               run_<lineage>/ at the top; with a fixed lineage there is one run
               folder.
    Does:      looks for run_*/full_table.tsv directly in the folder and, when there
               is none, one level deeper. Each table's header line states the size of
               its marker set ("number of BUSCOs: 2826"); the largest set is the most
               specific lineage and is the one used. A table without that line is
               sized by counting its marker rows.
    Produces:  (path, lineage name).
    Used by:   read_busco_markers().
    Stops with a clear error when the folder holds no full_table.tsv.
    """
    candidates = sorted(glob.glob(os.path.join(busco_dir, "run_*", "full_table.tsv")))
    if len(candidates) == 0:
        candidates = sorted(glob.glob(os.path.join(busco_dir, "*", "run_*", "full_table.tsv")))
    if len(candidates) == 0:
        raise ValueError(f"{busco_dir}: no run_*/full_table.tsv found in the BUSCO folder")

    best_path = None
    best_lineage = None
    best_size = -1
    for path in candidates:
        size = None
        lineage = os.path.basename(os.path.dirname(path))[len("run_"):]
        with open(path, "r", encoding="utf-8", errors="replace") as table:
            marker_ids = set()
            for line in table:
                if line.startswith("#"):
                    found = re.search(r"number of BUSCOs:\s*(\d+)", line)
                    if found:
                        size = int(found.group(1))
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 0 and fields[0].strip() != "":
                    marker_ids.add(fields[0].strip())
        if size is None:
            size = len(marker_ids)
        if size > best_size:
            best_path, best_lineage, best_size = path, lineage, size
    return best_path, best_lineage


def read_busco_markers(busco_dir, contig_names):
    """Count the complete BUSCO markers located on every contig.

    Takes in:  --busco-dir and the names of the draft contigs.
    Does:      reads the full_table.tsv chosen by find_busco_table(). Its columns are
               Busco id, Status, Sequence, Gene Start, Gene End, Strand, Score, Length,
               OrthoDB url, Description; a marker with status Complete or Duplicated
               was found whole on the sequence named in the third column (a
               Duplicated marker is a complete marker found twice, and both copies
               are counted where they lie). Some BUSCO versions write the sequence as
               name:start-end; when such a name is not a draft contig, the part
               before the last colon is used.
    Produces:  ({contig: markers}, lineage name).
    Used by:   run(), for the busco_markers column.
    """
    path, lineage = find_busco_table(busco_dir)
    counts = {}
    with open(path, "r", encoding="utf-8", errors="replace") as table:
        for line in table:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            status = fields[1].strip()
            if status not in ("Complete", "Duplicated"):
                continue
            sequence_name = fields[2].strip()
            if sequence_name not in contig_names and ":" in sequence_name:
                sequence_name = sequence_name.rsplit(":", 1)[0]
            counts[sequence_name] = counts.get(sequence_name, 0) + 1
    return counts, lineage


def read_gene_counts(path):
    """Count the protein-coding genes (mRNA features) per sequence in a GFF3 file.

    Takes in:  --gff3, funannotate's annotation. A GFF3 line has nine tab-separated
               columns: sequence name, source, feature type, start, end, score,
               strand, phase, attributes. Each protein-coding gene has one or more
               mRNA features (one per transcript); counting mRNA features gives the
               number of gene models, which is what funannotate reports as genes.
               tRNA features and the gene/exon/CDS lines are not counted.
    Does:      counts the lines whose third column is mRNA, per first column (the
               scaffold_N name). Lines starting with # are comments; the optional
               ##FASTA section at the end of a GFF3 file ends the reading.
    Produces:  {scaffold name: genes}.
    Used by:   run(), for the genes column (joined to the draft contigs through the
               scaffold mapping).
    Stops with a clear error on a data line with fewer than nine columns.
    """
    counts = {}
    with open(path, "r", encoding="utf-8", errors="replace") as gff3:
        for line_number, line in enumerate(gff3, start=1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(
                    f"{path}, line {line_number}: a GFF3 line needs 9 tab-separated "
                    f"columns but this one has {len(fields)}"
                )
            if fields[2].strip() == "mRNA":
                sequence_name = fields[0].strip()
                counts[sequence_name] = counts.get(sequence_name, 0) + 1
    return counts


# ======================================================================
# STEP 4: JOIN THE DRAFT CONTIGS TO FUNANNOTATE'S SCAFFOLD NAMES
# ======================================================================

def read_adapter_loci(path):
    """Read the adapter loci of the contig adapter screen.

    Takes in:  --adapter-screen (adapter_screen.tsv): one row per adapter locus, with the
               columns contig and position (end or internal) among others.
    Produces:  a dict contig -> the set of positions found on it ({"end"}, {"internal"} or
               both). A contig without a locus is not in the dict.
    Used by:   run(); adapter_column_text() turns it into the adapter column.
    """
    loci = {}
    for row in read_tsv_rows(path, {"contig", "position"}):
        position = (row["position"] or "").strip()
        if position not in ("end", "internal"):
            raise ValueError(f"{path}: position must be end or internal, got {position!r}")
        loci.setdefault(row["contig"], set()).add(position)
    return loci


def adapter_column_text(adapter_loci, contig):
    """Build the adapter column: none, at_end, internal or at_end+internal.

    Takes in:  the dict of read_adapter_loci() (None when the screen was not given) and the
               contig name.
    Produces:  NA without a screen; otherwise none, at_end (adapter within 100 bp of a contig
               end), internal (a whole adapter inside the contig) or both. The verdicts do not
               use it: adapter at a contig end changes neither what the contig is nor where it
               came from.
    Used by:   run().
    """
    if adapter_loci is None:
        return NA
    positions = adapter_loci.get(contig, set())
    if positions == {"end", "internal"}:
        return "at_end+internal"
    if positions == {"end"}:
        return "at_end"
    if positions == {"internal"}:
        return "internal"
    return "none"


def find_scaffold_by_probes(sequence, scaffolds, claimed):
    """Find the scaffold that holds a polished copy of this contig.

    Takes in:  one draft contig sequence, the scaffolds {name: sequence}, and the
               set of scaffold names already matched to another contig.
    Does:      the fallback of map_contigs_to_scaffolds() for a contig whose sequence
               is not identical to any scaffold. Polishing changes a base here and
               there, so most 300-bp pieces of the contig survive unchanged. PROBE_COUNT
               pieces of PROBE_LENGTH bases are taken from evenly spaced positions
               inside the contig (the ends are avoided, because the assembler trims
               and the polisher edits them most), and each unclaimed scaffold whose
               length is within PROBE_LENGTH_TOLERANCE of the contig's length is
               searched for them. The scaffold that holds most pieces is the match,
               provided it holds more than half of them and no other scaffold ties.
               The length test also rejects the one trap of the method: a contig that
               funannotate clean dropped as a near-duplicate of a longer contig would
               find its pieces in that longer scaffold.
    Produces:  the scaffold name, or None when no scaffold qualifies.
    Used by:   map_contigs_to_scaffolds().
    """
    length = len(sequence)
    if length < PROBE_COUNT * PROBE_LENGTH:
        return None

    probes = []
    for probe_number in range(1, PROBE_COUNT + 1):
        start = (length * probe_number) // (PROBE_COUNT + 1) - PROBE_LENGTH // 2
        probe = sequence[start:start + PROBE_LENGTH]
        # A piece with unknown bases (N) cannot be expected to match; leave it out.
        if "N" not in probe:
            probes.append(probe)
    if len(probes) == 0:
        return None

    hits = {}
    for scaffold, scaffold_sequence in scaffolds.items():
        if scaffold in claimed:
            continue
        if abs(len(scaffold_sequence) - length) > PROBE_LENGTH_TOLERANCE * length:
            continue
        found = 0
        for probe in probes:
            if probe in scaffold_sequence:
                found += 1
        if found > 0:
            hits[scaffold] = found
    if len(hits) == 0:
        return None

    # The scaffold holding most pieces; between equal counts the alphabetically first
    # name, so that the answer is the same on every run.
    def most_pieces_first(scaffold_and_found):
        scaffold, found = scaffold_and_found
        return (-found, scaffold)

    best_scaffold, best_found = sorted(hits.items(), key=most_pieces_first)[0]
    if best_found * 2 <= len(probes):
        return None   # not even a majority of the pieces: no confident match
    ties = [name for name, found in hits.items() if found == best_found]
    if len(ties) > 1:
        return None   # two scaffolds hold the same number of pieces: ambiguous
    return best_scaffold


def map_contigs_to_scaffolds(contigs, scaffolds):
    """Work out which scaffold_N of the annotated genome each draft contig became.

    Takes in:  the draft contigs and the annotated genome (--annotated-fasta), both
               as {name: sequence}.
    Does:      funannotate clean drops contigs that are near-duplicates of longer ones
               (and the mitochondrion and the removed contigs never reach it), and
               funannotate sort renames the rest scaffold_1, scaffold_2, ... in order
               of length, with no record of the old names. The only safe join is by
               sequence. First every contig whose sequence is identical to a
               scaffold's is matched (this is the case in the modes without polishing
               after the draft: PacBio HiFi samples, illumina mode, contigs mode). The
               contigs left over are tried with find_scaffold_by_probes(), longest
               first, because funannotate clean keeps the longer of two near-identical
               contigs and the longer one should get the first claim. Each scaffold
               is matched to one contig at most.
    Produces:  {contig: scaffold}; a contig without a match (dropped by funannotate
               clean, or not annotated at all) is absent.
    Used by:   run(), for the genes and annotated_as columns.
    """
    scaffold_by_sequence = {}
    for scaffold, sequence in scaffolds.items():
        scaffold_by_sequence[sequence] = scaffold

    mapping = {}
    claimed = set()
    for contig, sequence in contigs.items():
        scaffold = scaffold_by_sequence.get(sequence)
        if scaffold is not None and scaffold not in claimed:
            mapping[contig] = scaffold
            claimed.add(scaffold)

    def longest_first(contig):
        return (-len(contigs[contig]), contig)

    unmatched = [contig for contig in contigs if contig not in mapping]
    for contig in sorted(unmatched, key=longest_first):
        scaffold = find_scaffold_by_probes(contigs[contig], scaffolds, claimed)
        if scaffold is not None:
            mapping[contig] = scaffold
            claimed.add(scaffold)
    return mapping


# ======================================================================
# STEP 5: THE MEDIAN DEPTH
# ======================================================================

def weighted_median_depth(rows):
    """The typical read depth of the nuclear genome: a length-weighted median.

    Takes in:  the contig rows built in run(), each with length, depth (None when
               unknown) and the screen's verdicts.
    Does:      takes the contigs with a known depth that are neither the mitochondrion
               nor removed by the screen, keeps those of at least
               MEDIAN_MIN_CONTIG_LENGTH bases (all of them when none is that long),
               sorts them by depth and finds the depth at which the summed length
               passes half of the total. Weighting by length means a 4-Mb chromosome
               counts 100 times as much as a 40-kb piece, so the median is the depth
               of the bulk of the genome whatever the debris does.
    Produces:  the median depth, or None when no contig has a depth.
    Used by:   run(), for the depth_x_median column.
    """
    candidates = []
    for row in rows:
        if row["depth"] is None or row["is_mitochondrion"] or row["is_removed"]:
            continue
        candidates.append((row["depth"], row["length"]))
    if len(candidates) == 0:
        return None

    long_candidates = [pair for pair in candidates if pair[1] >= MEDIAN_MIN_CONTIG_LENGTH]
    if len(long_candidates) > 0:
        candidates = long_candidates

    candidates.sort()
    total_length = 0
    for _, length in candidates:
        total_length += length
    summed = 0
    for depth, length in candidates:
        summed += length
        if summed * 2 >= total_length:
            return depth
    return candidates[-1][0]


# ======================================================================
# STEP 6: THE VERDICT
# ======================================================================

def telomere_notes_text(row):
    """Name the telomere arrays of one contig that deserve a second look.

    Two cases, both read from tidk's table (see TELOMERES at the top):
      * a present array that begins more than FAR_FROM_END_BP inside the contig. A
        telomere is the very end of a chromosome, so an array that begins kilobases
        before the contig end says that sequence lies beyond the telomere: a
        misjoined piece, or an interstitial telomeric repeat. The offset columns
        carry the number; this puts it into the reason, where a reader scanning the
        verdicts sees it.
      * an array in the opposite orientation (status wrong_strand). The contig does
        carry the repeat, but it points the wrong way for a telomere, so it is no
        telomere call and the counts say 0; the reader should still know.
    Returns "" when neither applies at either end.
    """
    parts = []
    for side in TELOMERE_SIDES:
        offset = row["telomere_end_offset_" + side]
        if offset is not None and offset > FAR_FROM_END_BP:
            parts.append(f"the {side} array is {offset} bp from the end")
        if row["telomere_status_" + side] == "wrong_strand":
            parts.append(f"array in the opposite orientation at the {side}")
    return "; ".join(parts)


def telomere_text(row, min_copies):
    """Describe the telomere evidence of one contig in words, for the reason column."""
    start_copies = row["telomere_start"]
    end_copies = row["telomere_end"]
    at_start = start_copies >= min_copies
    at_end = end_copies >= min_copies
    if at_start and at_end:
        text = f"telomere repeats at both ends ({start_copies} and {end_copies} copies)"
    elif at_start:
        text = f"telomere repeats at the start only ({start_copies} copies)"
    elif at_end:
        text = f"telomere repeats at the end only ({end_copies} copies)"
    else:
        text = "no telomere repeat at either end"
    notes = telomere_notes_text(row)
    if notes:
        text += f" ({notes})"
    return text


def depth_text(row):
    """Describe the depth of one contig relative to the median, for the reason column."""
    ratio = row["depth_x_median"]
    if ratio is None:
        return "depth unknown"
    text = f"depth {ratio:.2f}x the median"
    if ratio < NORMAL_DEPTH_LOW:
        text += f" (below {NORMAL_DEPTH_LOW}x)"
    elif ratio > NORMAL_DEPTH_HIGH:
        text += f" (above {NORMAL_DEPTH_HIGH}x)"
    return text


def genes_text(row):
    """Describe the gene evidence (predicted genes, BUSCO markers), for the reason column."""
    if row["genes"] is None:
        parts = ["genes unknown"]
    else:
        parts = [f"{row['genes']} gene" + ("" if row["genes"] == 1 else "s")]
    if row["busco_markers"] is not None:
        parts.append(f"{row['busco_markers']} BUSCO marker" + ("" if row["busco_markers"] == 1 else "s"))
    return ", ".join(parts)


def decide(row, min_telomere_copies, min_rdna_units):
    """Give one contig its verdict and the reason, by the rules at the top of the file.

    Takes in:  one contig row built in run() (Python values, None for unknown) and the
               two thresholds (--min-telomere-copies, --min-rdna-units).
    Does:      tries the eleven rules in order; the first one that applies decides.
               Unknown evidence is skipped, never counted against the contig.
    Produces:  (verdict, reason).
    Used by:   run().
    """
    depth_ratio = row["depth_x_median"]
    depth_normal = depth_ratio is None or NORMAL_DEPTH_LOW <= depth_ratio <= NORMAL_DEPTH_HIGH
    genes = row["genes"]
    busco = row["busco_markers"]
    rdna = row["rdna_units"]
    at_start = row["telomere_start"] >= min_telomere_copies
    at_end = row["telomere_end"] >= min_telomere_copies

    # 1. The mitochondrion screen confirmed it, or the split step wrote it out.
    if row["is_mitochondrion"]:
        if row["mito_classification"] is not None:
            reason = (
                f"mitochondrion screen: {row['mito_classification']} "
                f"({row['mito_reason']})"
            )
        else:
            reason = "listed in the mitochondrial FASTA of the split step"
        reason += "; " + depth_text(row)
        if row["circular"] is True:
            reason += "; circular (assembler)"
        return "mitochondrion", reason

    # 2. The contamination screen removed it.
    if row["is_removed"]:
        return "contaminant", (
            f"removed by the screen: genus {row['genus']} ({row['screen_reason']})"
        )

    # 3. Mitochondrial DNA of an animal or a plant (human DNA from the library, for
    #    example), which the contamination screen kept: say so, rather than let the
    #    rules below call it nuclear sequence of the fungus.
    # (.get: a row built by hand, as some tests do, may leave the screen fields out.)
    if row.get("mito_classification") == "foreign_mitochondrial":
        return "foreign_mitochondrion", (
            f"mitochondrion screen: foreign_mitochondrial ({row['mito_reason']}); "
            f"genus {row['genus']}; {depth_text(row)}"
        )

    # 4. An rDNA array: enough rDNA units, and fewer genes than units (or unknown).
    if rdna is not None and rdna >= min_rdna_units and (genes is None or genes < rdna):
        reason = f"{rdna} rDNA units; {genes_text(row)}"
        if at_start or at_end:
            reason += "; " + telomere_text(row, min_telomere_copies)
        reason += "; " + depth_text(row)
        return "rDNA_array", reason

    # 5 and 6. Telomeres at both ends, or at one end, with normal depth and genes.
    genes_present_or_unknown = genes is None or genes > 0
    if (at_start or at_end) and depth_normal and genes_present_or_unknown:
        reason = f"{telomere_text(row, min_telomere_copies)}; {depth_text(row)}; {genes_text(row)}"
        if at_start and at_end:
            return "chromosome", reason
        return "chromosome_arm", reason

    # 7. Too shallow: an alternative copy of a region, or noise.
    if depth_ratio is not None and depth_ratio < NORMAL_DEPTH_LOW:
        return "low_depth", (
            f"{depth_text(row)}; {genes_text(row)}; {telomere_text(row, min_telomere_copies)}"
        )

    # Rules 8 to 10 rest on whether the contig carries genes. "Has genes": a known gene
    # count or BUSCO count above 0. "No genes": every known count is 0, and at least one
    # of the two is known. With neither count there is no gene evidence at all, and
    # rules 8 and 10 must not call a contig gene-less for want of looking. In the
    # workflow the BUSCO count is always known, so only a run by hand gets here without.
    # Rule 8 is looser on the gene count: a deep contig with no BUSCO marker and at most
    # REPEAT_MAX_GENES predicted genes still counts as a repeat.
    has_genes = (genes is not None and genes > 0) or (busco is not None and busco > 0)
    gene_evidence_known = genes is not None or busco is not None
    no_genes = gene_evidence_known and not has_genes
    # rDNA units below the rDNA_array bar still belong in the reason of a gene-less
    # contig: a deep piece with an rDNA unit on it is most likely collapsed rDNA.
    rdna_note = ""
    if rdna is not None and rdna > 0:
        rdna_note = f"; {rdna} rDNA unit" + ("" if rdna == 1 else "s")

    # 8. Too deep, no BUSCO marker, at most a gene call or two: a collapsed repeat.
    #    Several genome copies of a transposon, a satellite or another repeat family
    #    were assembled into one contig, so the reads of all the copies pile up on it.
    #    Needs a known depth.
    no_markers = busco is None or busco == 0
    few_genes = genes is None or genes <= REPEAT_MAX_GENES
    if (gene_evidence_known and no_markers and few_genes
            and depth_ratio is not None and depth_ratio > NORMAL_DEPTH_HIGH):
        reason = f"{depth_text(row)}; {genes_text(row)}"
        if genes:
            reason += " (a gene call inside a repeat is common)"
        return "repeat", reason + f"; {telomere_text(row, min_telomere_copies)}" + rdna_note

    # 9. Genuine nuclear sequence without a telomere: a chromosome fragment.
    if has_genes and depth_normal and not (at_start or at_end):
        return "nuclear_contig", (
            f"{genes_text(row)}; {depth_text(row)}; {telomere_text(row, min_telomere_copies)}"
        )

    # 10. Normal (or unknown) depth, no telomere, no gene: a piece of the nuclear genome
    #    too short to carry a gene.
    if no_genes and depth_normal and not (at_start or at_end):
        return "nuclear_fragment", (
            f"{genes_text(row)}; {depth_text(row)}; {telomere_text(row, min_telomere_copies)}"
            + rdna_note
        )

    # 11. Nothing above fits: say what was seen.
    reason = f"{telomere_text(row, min_telomere_copies)}; {depth_text(row)}; {genes_text(row)}"
    if rdna is not None:
        reason += f"; {rdna} rDNA unit" + ("" if rdna == 1 else "s")
    return "unclassified", reason


# ======================================================================
# STEP 7: WRITE THE OUTPUT FILES
# ======================================================================

def format_cell(value, decimals=None):
    """Turn a value into table text: NA for unknown, a fixed number of decimals for a float."""
    if value is None:
        return NA
    if value is True:
        return "Y"
    if value is False:
        return "N"
    if decimals is not None:
        return f"{value:.{decimals}f}"
    return str(value)


def formatted_row(row):
    """Return the output cells of one contig row, in the order of COLUMNS.

    Takes in:  one contig row from run(), with Python values.
    Produces:  a list of strings, one per column of COLUMNS.
    Used by:   write_table() and write_mqc_table().
    """
    return [
        row["contig"],
        format_cell(row["length"]),
        format_cell(row["gc"], 3),
        format_cell(row["depth"], 2),
        format_cell(row["depth_x_median"], 2),
        format_cell(row["assembler_depth"], 1),
        format_cell(row["circular"]),
        format_cell(row["telomere_start"]),
        format_cell(row["telomere_end"]),
        format_cell(row["telomere_end_offset_start"]),
        format_cell(row["telomere_end_offset_end"]),
        row["telomere_unit"],
        format_cell(row["rdna_units"]),
        format_cell(row["busco_markers"]),
        format_cell(row["genes"]),
        format_cell(row["genus"]),
        format_cell(row["screen"]),
        format_cell(row["adapter"]),
        format_cell(row["annotated_as"]),
        row["verdict"],
        row["reason"],
    ]


def write_table(path, rows):
    """Write contig_summary.tsv: the header line and one row per draft contig.

    Used by:   run(). Nothing in the workflow reads it back; it is the deliverable.
    """
    with open(path, "w", encoding="utf-8") as table:
        table.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            table.write("\t".join(formatted_row(row)) + "\n")


def mqc_header_lines(telomere_unit):
    """The header lines that make MultiQC pick the table up as custom content.

    Takes in:  the telomere unit of the sample (None when none was identified).
    Does:      builds the five lines: the id, the section title, a description shown
               under it, the plot type and the table's namespace. The description
               names the unit the telomere calls rest on, so that a reader of the
               report knows which repeat was looked for. MultiQC finds every
               *_mqc.tsv file below the output folder on its own.
    Produces:  a list of five lines, without line breaks.
    Used by:   write_mqc_table().
    """
    if telomere_unit is None:
        unit_text = "no telomere repeat unit could be identified for this sample"
    else:
        unit_text = f"telomere repeat unit {telomere_unit}"
    return [
        "# id: contig_summary",
        "# section_name: 'Contig summary'",
        "# description: 'What each contig of the draft assembly is (chromosome, chromosome "
        "arm, rDNA array, mitochondrion, mitochondrial DNA of another organism, "
        "contaminant, low-depth piece, collapsed repeat, "
        f"nuclear contig or fragment), judged from telomere repeats ({unit_text}), read "
        "depth, rDNA units, BUSCO markers and the gene prediction. The reason column names "
        "the evidence behind every verdict.'",
        "# plot_type: 'table'",
        "# pconfig: {namespace: 'Contig summary'}",
    ]


def write_mqc_table(path, rows, sample, telomere_unit):
    """Write the MultiQC copy of the table.

    Takes in:  the rows, the sample name and the telomere unit (None for none).
    Does:      writes the lines of mqc_header_lines(), then the header with a first
               column for the row name, then one row per contig whose row name is
               "<sample>: <contig>". MultiQC takes the first column as the name of
               the row, so the names must differ between samples for the tables of a
               batch to merge into one.
    Used by:   run(). Rule multiqc lists the file among its inputs.
    """
    with open(path, "w", encoding="utf-8") as table:
        for line in mqc_header_lines(telomere_unit):
            table.write(line + "\n")
        table.write("\t".join(["Sample"] + COLUMNS) + "\n")
        for row in rows:
            table.write("\t".join([f"{sample}: {row['contig']}"] + formatted_row(row)) + "\n")


# ======================================================================
# PUTTING THE STEPS TOGETHER
# ======================================================================

def run(args):
    """Do the whole job: the draft and the evidence in, the two tables out.

    Takes in:  args, the parsed command-line options (see parse_arguments()).
    Does:      runs the steps described at the top of the file. Every input is read
               and checked before the first output file is written, so bad input
               never leaves a half-written table behind.
    Produces:  the two output files and the line on standard error. Raises OSError
               or ValueError on bad input; main() turns that into a message.
    Used by:   main().
    """
    # Step 1: the draft contigs, the two files of the telomere step, and the optional
    # inputs that were given.
    contigs = read_fasta(args.contigs)
    contig_names = set(contigs)

    telomere_unit = read_telomere_unit(args.telomere_unit_file)
    telomere_ends = read_telomere_ends(args.telomere_ends, contig_names)
    # Without a unit, telomere_units.py writes the header line alone: rows next to a
    # unit of none mean that the two files do not come from the same run.
    if telomere_unit is None and len(telomere_ends) > 0:
        raise ValueError(
            f"{args.telomere_ends}: holds rows although {args.telomere_unit_file} says "
            "none, so the two telomere files do not belong together"
        )
    telomere_unit_text = telomere_unit if telomere_unit is not None else "none"

    mito_names = read_fasta_names(args.mito_contigs) if args.mito_contigs else set()
    decisions = read_decisions(args.decisions) if args.decisions else None
    mito_screen = read_mito_screen(args.mito_screen) if args.mito_screen else None
    depths = read_depths(args.depth_table) if args.depth_table else None
    evidence = read_assembly_evidence(args.assembly_evidence) if args.assembly_evidence else None
    rdna_counts = read_its_window_counts(args.its_windows) if args.its_windows else None
    busco_counts = None
    busco_lineage = None
    if args.busco_dir:
        busco_counts, busco_lineage = read_busco_markers(args.busco_dir, contig_names)
    gene_counts = read_gene_counts(args.gff3) if args.gff3 else None
    adapter_loci = read_adapter_loci(args.adapter_screen) if args.adapter_screen else None
    scaffold_mapping = None
    if args.annotated_fasta:
        scaffolds = read_fasta(args.annotated_fasta)
        scaffold_mapping = map_contigs_to_scaffolds(contigs, scaffolds)

    # Step 2: one row per contig with the telomere columns and the evidence joined by
    # name. The telomere unit is the same on every row of the sample.
    rows = []
    ends_with_telomere = 0
    for contig, sequence in contigs.items():
        row = {"contig": contig, "length": len(sequence), "gc": gc_fraction(sequence)}
        row.update(telomere_columns(telomere_ends.get(contig, {}), telomere_unit))
        row["telomere_unit"] = telomere_unit_text
        for end_name in ("telomere_start", "telomere_end"):
            if row[end_name] >= args.min_telomere_copies:
                ends_with_telomere += 1

        # The screens. A contig is the mitochondrion when the screen confirmed it or the
        # split step wrote it out; it is removed when the decisions table says so.
        classification, mito_reason = (None, None)
        if mito_screen is not None and contig in mito_screen:
            classification, mito_reason = mito_screen[contig]
        row["mito_classification"] = classification
        row["mito_reason"] = mito_reason
        row["is_mitochondrion"] = classification == "confirmed_mitochondrial" or contig in mito_names

        genus, action, screen_reason = (None, None, None)
        if decisions is not None and contig in decisions:
            genus, action, screen_reason = decisions[contig]
        row["genus"] = genus
        row["screen_reason"] = screen_reason
        row["is_removed"] = action == "remove"
        row["screen"] = screen_column_text(row, action)
        row["adapter"] = adapter_column_text(adapter_loci, contig)

        # Depth and the assembler's evidence.
        row["depth"] = depths.get(contig) if depths is not None else None
        assembler_depth, circular = (None, None)
        if evidence is not None and contig in evidence:
            assembler_depth, circular = evidence[contig]
        row["assembler_depth"] = assembler_depth
        row["circular"] = circular

        # rDNA units, BUSCO markers and genes. A contig that is in the ITS table or the
        # BUSCO table with no row has none of them, which is 0, not unknown.
        row["rdna_units"] = rdna_counts.get(contig, 0) if rdna_counts is not None else None
        row["busco_markers"] = busco_counts.get(contig, 0) if busco_counts is not None else None
        scaffold = scaffold_mapping.get(contig) if scaffold_mapping is not None else None
        row["genes"] = None
        if gene_counts is not None and scaffold is not None:
            row["genes"] = gene_counts.get(scaffold, 0)
        row["annotated_as"] = annotated_as_text(row, scaffold, scaffold_mapping is not None)

        rows.append(row)

    # Step 3: the median depth, every contig's depth relative to it, and the verdicts.
    median_depth = weighted_median_depth(rows)
    verdict_counts = {}
    for row in rows:
        row["depth_x_median"] = None
        if row["depth"] is not None and median_depth:
            row["depth_x_median"] = row["depth"] / median_depth
        row["verdict"], row["reason"] = decide(row, args.min_telomere_copies, args.min_rdna_units)
        verdict_counts[row["verdict"]] = verdict_counts.get(row["verdict"], 0) + 1

    # Step 4: the two output files.
    write_table(args.out, rows)
    write_mqc_table(args.out_mqc, rows, args.sample, telomere_unit)

    # One summary line for the log.
    counts_text = ", ".join(
        f"{count} {verdict}" for verdict, count in sorted(verdict_counts.items())
    )
    if telomere_unit is None:
        telomere_summary = (
            f"no telomere repeat unit identified, so no telomere at any of the "
            f"{2 * len(rows)} contig ends"
        )
    else:
        telomere_summary = (
            f"telomere repeat unit {telomere_unit} at {ends_with_telomere} of "
            f"{2 * len(rows)} contig ends"
        )
    lineage_text = f"; BUSCO lineage {busco_lineage}" if busco_lineage else ""
    median_text = f"; median depth {median_depth:.2f}" if median_depth is not None else ""
    print(
        f"contig_summary.py: sample {args.sample}: {len(rows)} contigs: {counts_text}; "
        f"{telomere_summary}{median_text}{lineage_text}.",
        file=sys.stderr,
    )


def screen_column_text(row, action):
    """Build the screen column: kept / removed (with the reason) / mitochondrion.

    Takes in:  the contig row (its mitochondrion fields are already filled) and the
               action of the decisions table (keep, remove, or None without the table).
    Produces:  "mitochondrion (<reason>)" for a confirmed mitochondrion, "kept (<reason>)"
               or "removed (<reason>)" from the decisions table, with the mitochondrion
               screen's classification added when it is a candidate, an unclassified
               circle or deep contig, or another organism's mitochondrial DNA
               (foreign_mitochondrial): something the reader should know about. None
               when neither screen ran.
    Used by:   run().
    """
    if row["is_mitochondrion"]:
        if row["mito_reason"]:
            return f"mitochondrion ({row['mito_reason']})"
        return "mitochondrion"
    if action is None:
        return None
    text = "kept" if action == "keep" else "removed"
    if row["screen_reason"]:
        text += f" ({row['screen_reason']})"
    classification = row["mito_classification"]
    if classification not in (None, "", "not_mitochondrial", "not_screened"):
        text += f"; {classification} ({row['mito_reason']})"
    return text


def annotated_as_text(row, scaffold, annotated_fasta_given):
    """Build the annotated_as column: under which name funannotate annotated the contig.

    Takes in:  the contig row, its scaffold from the mapping (None when unmatched), and
               whether --annotated-fasta was given at all.
    Produces:  "scaffold_N" whenever the contig has a scaffold in the annotated genome,
               whatever the screens said about it: with parameters.mitochondria.mode
               set to report, the mitochondrion screen confirms the mitochondrion but
               does not split it off, so the contig reaches funannotate and IS
               annotated under some scaffold_N, and the column must say so. Without a
               scaffold: "not annotated (mitochondrion)"; "not annotated (removed by
               the screen)"; "dropped by funannotate clean" for a kept nuclear contig;
               None when the annotated genome was not given.
    Used by:   run().
    """
    if scaffold is not None:
        return scaffold
    if row["is_mitochondrion"]:
        return "not annotated (mitochondrion)"
    if row["is_removed"]:
        return "not annotated (removed by the screen)"
    if not annotated_fasta_given:
        return None
    return "dropped by funannotate clean"


# ======================================================================
# COMMAND LINE
# ======================================================================

def parse_arguments(argv):
    """Read the options given on the command line.

    Takes in:  argv, the list of words after the script name (None means "use the real
               command line"; the tests pass their own list).
    Does:      defines the options with argparse (Python's built-in reader of
               command-line options). The draft, the two telomere files, the sample
               name and the two outputs are required; every other input is
               optional, and an option given with an empty value counts as not
               given, because Snakemake writes an empty parameter as nothing. The
               two tuning options are checked against the same limits the workflow
               applies to the configuration.
    Produces:  an object whose attributes are the options (dashes become underscores:
               args.out_mqc, args.telomere_unit_file, ...).
    Used by:   main().
    """
    parser = argparse.ArgumentParser(
        prog="contig_summary.py",
        description=(
            "Gather the evidence FunFlux produced about every contig of the draft "
            "assembly into one table and give each contig a verdict with its reasons."
        ),
    )
    parser.add_argument("--contigs", required=True, help="the draft assembly (FASTA)")
    parser.add_argument(
        "--telomere-unit-file", required=True,
        help="telomere_unit.txt of telomere_units.py: the telomere repeat unit of the sample, or none",
    )
    parser.add_argument(
        "--telomere-ends", required=True,
        help="telomere_ends.tsv of telomere_units.py: tidk's table of what it found at every contig end",
    )
    # The ten optional inputs may be given with nothing after them: nargs="?" (zero or
    # one value) allows that, and the flag then takes its const value, an empty text,
    # meaning "not given". Snakemake's {...:q} quoting writes an empty value as nothing
    # at all, so a rule can write every flag and leave the value empty in the modes
    # that lack the file, as rule select_contigs does with its policy flags.
    optional_inputs = [
        ("--mito-contigs", "the mitochondrial contigs FASTA written by the split step"),
        ("--decisions", "contig_taxonomy_decisions.tsv of the contamination screen"),
        ("--mito-screen", "mitochondrial_screen.tsv of the mitochondrion screen"),
        ("--depth-table", "the samtools coverage table of the reads mapped to the draft "
                          "(read depth in column meandepth); not in contigs mode, where no "
                          "reads are mapped"),
        ("--assembly-evidence", "assembly_evidence.tsv (long-read modes)"),
        ("--its-windows", "its_windows.tsv, one row per rDNA unit found"),
        ("--busco-dir", "the BUSCO output folder of the sample"),
        ("--gff3", "the funannotate GFF3 annotation"),
        ("--annotated-fasta", "the genome as funannotate annotated it (scaffold_N names)"),
        ("--adapter-screen", "adapter_screen.tsv of the contig adapter screen"),
    ]
    for flag, help_text in optional_inputs:
        parser.add_argument(flag, nargs="?", const="", default="", help=help_text + " (optional)")
    parser.add_argument(
        "--min-telomere-copies", type=int, default=DEFAULT_MIN_TELOMERE_COPIES,
        help="copies of the telomere unit from which a telomere is called (default: %(default)s)",
    )
    parser.add_argument(
        "--min-rdna-units", type=int, default=DEFAULT_MIN_RDNA_UNITS,
        help="rDNA units from which a contig can be an rDNA array (default: %(default)s)",
    )
    parser.add_argument("--sample", required=True, help="the sample name, for the MultiQC row names")
    parser.add_argument("--out", required=True, help="output: contig_summary.tsv")
    parser.add_argument("--out-mqc", required=True, help="output: the MultiQC copy of the table")

    args = parser.parse_args(argv)

    # The same checks as the workflow makes on the configuration block.
    if args.min_telomere_copies < MIN_TELOMERE_COPIES_MIN:
        parser.error(f"--min-telomere-copies must be at least {MIN_TELOMERE_COPIES_MIN}")
    if args.min_rdna_units < MIN_RDNA_UNITS_MIN:
        parser.error(f"--min-rdna-units must be at least {MIN_RDNA_UNITS_MIN}")
    return args


def main(argv=None):
    """Run the script from the command line.

    Takes in:  argv, the command-line words (None = the real command line).
    Does:      parses the options and calls run(). If run() finds a problem with an
               input, prints one readable line to stderr instead of a Python traceback
               and exits with code 1, so Snakemake sees that the step failed.
    Produces:  exit code 0 on success, 1 on a problem with the input files, 2 on a bad
               command line (argparse's own code).
    Used by:   the last lines of this file when rule contig_summary runs the script,
               and by the tests, which pass their own argv.
    """
    args = parse_arguments(argv)
    try:
        run(args)
    except (OSError, ValueError) as problem:
        # OSError: a file cannot be opened or written (for example it does not exist).
        # ValueError: a table lacks a column, a number cannot be read, a FASTA has two
        # records of one name, a GFF3 line is short, the telomere table names a contig
        # the draft lacks or carries a status tidk does not write, or a file is not
        # UTF-8 text.
        print(f"contig_summary.py: error: {problem}", file=sys.stderr)
        sys.exit(1)


# Run main() only when this file is run as a script, not when the tests import it.
if __name__ == "__main__":
    main()
