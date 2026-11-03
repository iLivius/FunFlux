#!/usr/bin/env python3
"""Choose the telomere repeat unit of a draft assembly from the data, and call its telomeres.

BIOLOGY AND WHY THIS SCRIPT EXISTS
----------------------------------
Telomeres are tandem arrays of a short repeat that cap every chromosome end. A contig
that carries the array at both ends is a whole chromosome, one that carries it at one
end a chromosome arm: this is the evidence the contig summary (contig_summary.py)
builds its verdicts on. The repeat unit, however, differs between fungal groups. It is
TTAGGG, the vertebrate unit, in most filamentous Ascomycota, most Basidiomycota and the
Mucoromycota, but TTAGGGTCAACA in the Aspergillus flavus group, TTTAGGG in Penicillium
polonicum, TTAGTG in some Malassezia, and a 23- or 25-bp unit in the budding yeasts
Candida albicans and Kluyveromyces lactis (the constant LITERATURE lists the known cases).
Counting one assumed unit would call no telomere at all in a genome that uses another
unit, and the user would not know why. So this script chooses the unit from the
assembly itself, by how the candidate units behave at the contig ends, and then says
whether the chosen unit is a known telomere repeat: named in the literature, or carried
by species of the telomere repeat database that comes with tidk.

WHAT TIDK DOES AND WHAT THIS SCRIPT ADDS
----------------------------------------
tidk (Brown et al. 2025) does the two heavy parts:
  * `tidk explore` lists every repeat unit that forms a tandem run of at least a given
    number of copies anywhere in a set of sequences, with its copy number. Each unit is
    reported in tidk's canonical form: the alphabetically smallest of all rotations of
    the unit and of its reverse complement, so that TTAGGG, GGGTTA and CCCTAA (the same
    array read from another position, or from the other strand) are all written AACCCT.
  * `tidk ends` takes one unit and calls a telomere of that unit at both ends of every
    sequence: an array of at least a minimum length, within a maximum distance of the
    end, on the strand a telomere has there. The G-rich strand of a telomere runs
    towards the chromosome end, so the array reads CCCTAA (C-rich) at the start of a
    contig and TTAGGG (G-rich) at its end. tidk's status per contig end is `present`
    (such an array), `wrong_strand` (an array on the other strand), `not_terminal` (an
    array further inside than allowed) or `absent`.
The choice between the units, the check against what is known and the three output
files are FunFlux's own (THE METHOD).

WHERE THIS SCRIPT FITS IN THE WORKFLOW
--------------------------------------
    1. The draft assembly exists (DRAFT_CONTIGS: every contig the assembler made, the
       ones the screens later remove included) and the telomere repeat database of
       tidk was downloaded once per run by rule telomere_db
       (workflow/rules/shared/25_contig_summary.smk). The names in capitals are the
       names used in workflow/rules/shared/00_common.smk.
    2. THIS SCRIPT is run once per sample by rule telomere_units (the same rule file),
       in every mode, in the tidk conda environment (workflow/envs/tidk.yaml), and
       writes the three OUTPUTS into the sample's eval folder.
    3. Rule contig_summary reads telomere_unit.txt and telomere_ends.tsv and turns the
       per-end calls into the telomere columns and the verdicts of the contig summary.
       The candidates table is a record for the reader; no rule reads it.

THE METHOD
----------
Three settings: the window W (--telomere-window, 25,000 bp by default: how far from a
contig end an array may sit, because the assembler now and then runs past the telomere
and the array then lies some kilobases inside the contig end), the minimum copies C
(--min-telomere-copies, 5 by default: tidk's minimum array length is C times the unit
length, 30 bp of TTAGGG, and a unit must form a run of at least C copies somewhere to
be listed at all), and --telomere-motif, `auto` or a unit (A CONFIGURED UNIT below).
    1. Contig ends. The first and the last W bases of every contig longer than 2W are
       written to a FASTA file as <name>__start and <name>__end; a contig of at most 2W
       bases goes in whole, under its own name.
    2. Units found at the ends. `tidk explore` on that file lists every unit with a
       run of at least C copies and its copies in the windows (copies_at_ends); the
       same run on the whole draft gives copies_genome_wide. A telomere unit has
       nearly all of its copies at the ends; a satellite has most of them inside the
       contigs. The table shows both numbers, but they do not enter the score: a
       microsatellite such as AAAGAG can have nearly as many copies at the ends as
       the telomere unit, and the copy number does not say which way its arrays
       point.
    3. Candidates: every unit found at the ends in step 2, every unit of the database
       with strong or weak support, every unit of THE LITERATURE, and the configured
       unit when there is one; all in canonical form, duplicates merged. Three filters
       apply to every candidate except the configured unit. The unit is 5 to 30 bp long
       (no telomere repeat shorter than 5 bp is known, and 4-bp units are
       low-complexity artefacts that tidk ends finds at almost every contig end). It is
       not a repeat of a shorter unit (ACACAC is AC; TTAGGGTTAGGG would score the very
       arrays of TTAGGG a second time). Its two strands differ, that is, it holds more
       G than C or more C than G: every telomerase-made repeat has a G-rich strand, and
       that strand is what lets tidk tell the orientation of an array; AAAT has no
       G-rich strand, AGAGCC has two identical strands, and for such a unit every
       array counts as a telomere whichever way it points.
    4. Score. For every candidate, `tidk ends` on the whole draft with window W,
       maximum offset W and minimum array length C x unit length, in its default
       error-tolerant mode (a sequencing error inside the array does not break it).
       Then one correction per contig, the same-array rule: when the start row and the
       end row of a contig both report an array and the two coordinate intervals
       overlap, they are one array seen from both ends of a short contig (tidk looks W
       bases in from each end, and a contig shorter than 2W is looked through
       twice), and the row with the larger offset gets the status `other_end`, a
       FunFlux status, not tidk's. After that: telomere_ends = rows with status
       present, opposite_ends = rows with status wrong_strand, score = telomere_ends -
       opposite_ends, t2t_contigs = contigs with present at both ends. The score works
       because a telomere repeat points one way at every chromosome end, whereas a
       satellite that happens to sit at contig ends points either way with equal
       chance, so its score stays near zero however often it occurs.
    5. Choice. A candidate is eligible with at least 2 telomere ends and more than
       twice as many telomere ends as opposite ends. The candidates are ranked by
       score, then known units (literature or database) before unknown ones, then the
       shorter unit, then alphabetically; the first eligible candidate is the telomere
       unit of the assembly, and `none` when no candidate is eligible.
On the ONT assembly of the public benchmark genome (Botrytis cinerea B05.10, 33
contigs) this gives AACCCT, that is TTAGGG, at 34 contig ends with 15 contigs capped at
both ends, while the most copied microsatellites score about 0; on the short-read
assembly of the same genome (413 contigs) AACCCT has 13 telomere ends and 1 opposite
end, after six duplicates of short contigs became other_end.

A CONFIGURED UNIT
-----------------
When --telomere-motif names a unit instead of `auto`, that unit is the one written to
telomere_unit.txt and telomere_ends.tsv and marked chosen, whatever its score. The
ranking still runs, and a NOTE line says when the data favour another unit, or give the
configured one no support. The configured unit is exempt from the three filters.

THE LITERATURE
--------------
The units named for fungi in the literature, with the published spelling, the fungi
and the references, are the constant LITERATURE below; every one of them is a
candidate in every run. The budding yeasts with irregular repeats, the fission yeasts,
Cryptococcus (TTA followed by four to six G) and the telomerase-negative
Microstromatales have no exact unit and are not listed: their telomeres are not found
by an exact unit, and the candidates table then shows what the contig ends do carry.

INPUTS
------
    --contigs              DRAFT_CONTIGS, the draft assembly (FASTA). Names are the
                           first word of the header, as every other table uses them.
    --telomere-db          TELOMERE_DB, the telomere repeat database of tidk
                           (telo_database_v2.tsv): tab-separated, one row per species
                           with the columns species, kingdom, telomere_unit (canonical
                           form) and support (strong, weak, unchecked or none), among
                           others. Only the rows with strong or weak support count.
    --telomere-motif       auto, or the telomere repeat unit of the species, 4 to 30
                           letters A, C, G and T, in any rotation or orientation.
    --telomere-window      W, 500 to 100000 (default 25000).
    --min-telomere-copies  C, at least 2 (default 5).
    --sample               the sample name, for the messages.
    --tidk                 the tidk program (default: tidk, found on PATH).
The workflow passes the three settings from the configuration block
parameters.contig_summary (keys telomere_motif, telomere_window, min_telomere_copies).

OUTPUTS
-------
    --out-candidates   TELOMERE_CANDIDATES, telomere_candidates.tsv: one row per
                       candidate that was found at the contig ends or formed at least
                       one array (a known unit that did neither is counted in the log
                       only; the configured unit always has a row), the columns of
                       CANDIDATE_COLUMNS: unit, unit_length, copies_at_ends,
                       copies_genome_wide, telomere_ends, opposite_ends, score,
                       t2t_contigs, in_literature (the published spelling, the fungi
                       and the references, or -), in_database (species with strong or
                       weak support per kingdom, for example "Fungi 20; Metazoa 420;
                       Viridiplantae 11", with the few rows without a kingdom counted
                       as "(no kingdom)", or -) and chosen (yes or no). Sorted by
                       score, then telomere_ends, both descending, then unit.
    --out-unit         TELOMERE_UNIT_FILE, telomere_unit.txt: one line, the unit in
                       canonical form, or none.
    --out-ends         TELOMERE_ENDS, telomere_ends.tsv: the per-end table of the chosen
                       unit after the same-array rule: tidk's ten columns (id,
                       sequence_length, side, status, telomere_start, telomere_end,
                       telomere_length, offset, strand, window_limited) and values,
                       the status also other_end; the header line alone when the
                       unit is none.
    standard error     one plain line with the counts and the run time, then the
                       summary line "[FunFlux] telomere_units <sample>: ..." (the unit,
                       how many contig ends and contigs with both ends, and whether
                       the unit is known) or "... no telomere repeat identified (best
                       candidate ...)", then "[FunFlux] NOTE telomere_units <sample>:
                       ..." when the chosen unit is unknown to the literature and the
                       database, or when a configured unit is not the data's first
                       choice. Rule telomere_units writes the stream to its log and
                       copies the [FunFlux] lines to the console.

Exit status 0 when the three files were written. A missing file, a database without
the four columns, a tidk failure (the message includes tidk's own), a motif that is
neither auto nor a unit of 4 to 30 letters A, C, G and T, a window outside 500 to
100000 or fewer than 2 copies stops the script with one readable message starting
"[FunFlux] telomere_units:" and exit status 1, and no output is written. tidk's own
files go to a temporary folder that is removed at the end.

REFERENCES
----------
Brown MR, Gonzalez de la Rosa PM, Blaxter M (2025). tidk: a toolkit to rapidly
    identify telomeric repeats from genomic datasets. Bioinformatics 41(2): btaf049.
    https://doi.org/10.1093/bioinformatics/btaf049
Lue NF (2021). Duplex telomere-binding proteins in fungi with canonical telomere
    repeats: new lessons in the rapid evolution of telomere proteins. Frontiers in
    Genetics 12: 638790. https://doi.org/10.3389/fgene.2021.638790
Sun Q, Wang H, Tao S, Xi X (2023). Large-scale detection of telomeric motif sequences
    in genomic data using TelFinder. Microbiology Spectrum 11(2): e03928-22.
    https://doi.org/10.1128/spectrum.03928-22
Levis C, Giraud T, Dutertre M, Fortini D, Brygoo Y (1997). Telomeric DNA of Botrytis
    cinerea: a useful tool for strain identification. FEMS Microbiology Letters
    157(2): 267-272. https://doi.org/10.1111/j.1574-6968.1997.tb12783.x
Červenák F, Sepšiová R, Nosek J, Tomáška Ľ (2021). Step-by-step evolution of
    telomeres: lessons from yeasts. Genome Biology and Evolution 13(2): evaa268.
    https://doi.org/10.1093/gbe/evaa268
McEachern MJ, Blackburn EH (1994). A conserved sequence motif within the
    exceptionally diverse telomeric sequences of budding yeasts. Proceedings of the
    National Academy of Sciences USA 91(8): 3453-3457.
    https://doi.org/10.1073/pnas.91.8.3453
"""

# Only Python's standard library is used; tidk is run as a separate program with
# subprocess. The tidk conda environment holds both tidk and a Python 3.12.
import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import OrderedDict


# ======================================================================
# SETTINGS
# ======================================================================

# Defaults of the three settings and of the tidk program. The workflow passes the
# values of its contig_summary configuration block; a run by hand without the options
# uses these.
DEFAULT_TELOMERE_MOTIF = "auto"
DEFAULT_TELOMERE_WINDOW = 25000
DEFAULT_MIN_TELOMERE_COPIES = 5
DEFAULT_TIDK = "tidk"

# The limits of the settings. They are the same limits the workflow applies when it
# reads the configuration, so that a run by hand is held to the same rules. A configured
# unit may be 4 letters long (the filters below do not apply to it); the window has a
# ceiling because tidk looks that far in from every contig end, and a window of
# chromosome size would no longer describe the contig ends.
MOTIF_MIN_LENGTH = 4
MOTIF_MAX_LENGTH = 30
TELOMERE_WINDOW_MIN = 500
TELOMERE_WINDOW_MAX = 100000
MIN_TELOMERE_COPIES_MIN = 2

# The unit lengths tidk explore searches, and the length filter of the candidates: the
# same range. No telomere repeat shorter than 5 bp is known, and almost all are shorter
# than 30 (Brown et al. 2025); the longest fungal unit in the literature is 25 bp.
UNIT_MIN_LENGTH = 5
UNIT_MAX_LENGTH = 30

# tidk explore counts, by default, only the repeats in the outer 1 % of every sequence
# (its --distance option, a proportion of the sequence length from each end). The
# contig-end windows ARE the ends, every base of them must count, and the genome-wide
# count must cover the whole contigs: 0.5 from each end is the whole sequence.
EXPLORE_DISTANCE = 0.5

# When a candidate is eligible: at least this many telomere ends, and more than this
# many times as many telomere ends as opposite ends. 2 and 2: one end could be an
# interstitial array or a misjoin; a unit that points the wrong way at a third of its
# ends is not a telomere.
MIN_TELOMERE_ENDS = 2
OPPOSITE_ENDS_FACTOR = 2

# tidk's statuses, and the one this script adds.
STATUS_PRESENT = "present"
STATUS_ABSENT = "absent"
STATUS_WRONG_STRAND = "wrong_strand"
STATUS_OTHER_END = "other_end"

# What the unit file says when no unit was found, and what a table cell says when a
# unit is not known from the literature or the database.
NO_UNIT = "none"
NOT_KNOWN = "-"

# The database: the columns read, the support values that count, and the label of a
# row whose kingdom cell is empty (a few rows of the database have none).
DATABASE_COLUMNS = {"species", "kingdom", "telomere_unit", "support"}
DATABASE_SUPPORT = {"strong", "weak"}
NO_KINGDOM = "(no kingdom)"

# The columns of the two tidk tables, as tidk 0.3.0 writes them. A different header
# means another tidk version, and the script stops rather than misread the table.
EXPLORE_COLUMNS = ["canonical_repeat_unit", "copies", "copies_as_unit", "copies_as_revcomp"]
ENDS_COLUMNS = [
    "id", "sequence_length", "side", "status", "telomere_start", "telomere_end",
    "telomere_length", "offset", "strand", "window_limited",
]

# The columns of telomere_candidates.tsv, in this order.
CANDIDATE_COLUMNS = [
    "unit", "unit_length", "copies_at_ends", "copies_genome_wide",
    "telomere_ends", "opposite_ends", "score", "t2t_contigs",
    "in_literature", "in_database", "chosen",
]

# The telomere repeat units named for fungi in the literature: the published spelling,
# the fungi, and the references by author and year (THE LITERATURE at the top). Every
# unit here is a candidate in every run, and the in_literature column of the candidates
# table carries this text for the unit. The spelling is the one the papers use; the
# script compares units in canonical form, so any rotation or orientation of it is the
# same unit.
LITERATURE = [
    ("TTAGGG",
     "filamentous Ascomycota, most Basidiomycota, Mucoromycota; Botrytis cinerea",
     "Lue 2021; Sun et al. 2023; Levis et al. 1997"),
    ("TTAGGGTCAACA", "Aspergillus flavus group", "Sun et al. 2023"),
    ("TAGGG", "Peltaster fructicola", "Sun et al. 2023"),
    ("CCGTTGCTGTTG", "Drechmeria coniospora", "Sun et al. 2023"),
    ("TTTAGGG", "Penicillium polonicum", "Sun et al. 2023"),
    ("TTAGGA", "Malassezia furfur", "Sun et al. 2023"),
    ("TTAGTG", "Malassezia restricta, Malassezia globosa", "Sun et al. 2023"),
    ("TTGGA", "Hericium erinaceus", "Sun et al. 2023"),
    ("TTAGTCAGGG", "Yarrowia lipolytica", "Lue 2021; Červenák et al. 2021"),
    ("ACGGATGTCTAACTTCTTGGTGT", "Candida albicans", "McEachern and Blackburn 1994"),
    ("ACGGATTTGATTAGGTATGTGGTGT", "Kluyveromyces lactis", "McEachern and Blackburn 1994"),
]


# ======================================================================
# THE DRAFT ASSEMBLY AND ITS CONTIG ENDS
# ======================================================================

def read_fasta(path):
    """Read a FASTA file into {contig name: sequence in upper case}, in file order.

    Takes in:  the path of the draft assembly (--contigs). Sequences may be wrapped
               over several lines.
    Does:      reads the file line by line. The contig name is the first word of the
               ">" header, which is the name every table of the workflow uses (and the
               name tidk writes in its tables). Sequence lines are joined and put in
               upper case, so that soft-masked (lower-case) bases are the same bases.
    Produces:  an ordered dictionary {name: sequence}.
    Used by:   run(), to write the contig-end windows.
    Stops with a clear error when two records share a name (the per-end table could
    not be joined to the contigs), when a sequence line comes before any header (not
    a FASTA file), or when the file holds no record at all.
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
                "the per-end table could not be joined to the contigs"
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


def write_contig_ends(contigs, window, path):
    """Write the terminal windows of every contig as one FASTA file, for tidk explore.

    Takes in:  the draft contigs {name: sequence}, the window length W
               (--telomere-window) and the path to write.
    Does:      for a contig longer than 2W, writes its first W bases as <name>__start
               and its last W bases as <name>__end. A contig of at most 2W bases is
               written whole, under its own name: its two windows would overlap, and an
               array in the overlap would be counted twice. One line per sequence;
               tidk reads that, and nothing else reads the file.
    Produces:  the file; returns the number of records written.
    Used by:   run(), step 1 of the method. tidk explore on this file gives the units
               found at the contig ends and their copies_at_ends.
    """
    records = 0
    with open(path, "w", encoding="utf-8") as fasta:
        for name, sequence in contigs.items():
            if len(sequence) > 2 * window:
                fasta.write(f">{name}__start\n{sequence[:window]}\n")
                fasta.write(f">{name}__end\n{sequence[-window:]}\n")
                records += 2
            else:
                fasta.write(f">{name}\n{sequence}\n")
                records += 1
    return records


# ======================================================================
# REPEAT UNITS: CANONICAL FORM AND THE THREE FILTERS
# ======================================================================

# The swap table of reverse_complement(): A<->T and C<->G; N stays N. Any other
# letter is left unchanged.
COMPLEMENT_TABLE = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence):
    """Return the sequence as it reads on the opposite strand.

    Takes in:  a piece of DNA in upper case.
    Does:      swaps every base for its partner (A<->T, C<->G) and reads the result
               backwards, because the opposite strand runs the other way.
    Produces:  the reverse-complemented text, same length.
    Used by:   canonical_unit() and spelling_note().
    """
    return sequence.translate(COMPLEMENT_TABLE)[::-1]


def rotations(unit):
    """Return every rotation of a unit: the same tandem repeat read from each position.

    Takes in:  a repeat unit, for example TTAGGG.
    Produces:  the list [TTAGGG, TAGGGT, AGGGTT, GGGTTA, GGTTAG, GTTAGG].
    Used by:   canonical_unit() and spelling_note().
    """
    forms = []
    for shift in range(len(unit)):
        forms.append(unit[shift:] + unit[:shift])
    return forms


def canonical_unit(unit):
    """Return the canonical spelling of a tandem repeat unit, as tidk writes it.

    Takes in:  a repeat unit, for example TTAGGG.
    Does:      a tandem array can be read from any of its positions and from either
               strand, so the same array yields TTAGGG, TAGGGT, GGGTTA, ... on one
               strand and CCCTAA, CCTAAC, ... on the other. All rotations of the unit
               and of its reverse complement are listed, and the alphabetically
               smallest one is the canonical form (here AACCCT). This is how tidk
               collapses equivalent repeats (Brown et al. 2025), so the units of tidk
               explore, of the database (already canonical) and of the literature can
               be compared.
    Produces:  the canonical form, a string of the same length.
    Used by:   read_telomere_database(), literature_by_unit() and run() (the
               configured unit).
    """
    forms = rotations(unit) + rotations(reverse_complement(unit))
    return min(forms)


def is_repeat_of_shorter_unit(unit):
    """Tell whether a unit is itself a tandem repeat of a shorter unit.

    Takes in:  a candidate unit, for example TTAGGGTTAGGG or ACACAC.
    Does:      for every length that divides the unit's length, checks whether the
               unit is that many copies of its first bases. A homopolymer (AAAAAA) is
               the case of length 1.
    Produces:  True for such a unit, False for a genuine unit.
    Used by:   passes_candidate_filters(): the shorter unit is a candidate on its own,
               and the longer one would score the very same arrays a second time.
    """
    unit_length = len(unit)
    for shorter_length in range(1, unit_length):
        if unit_length % shorter_length != 0:
            continue
        if unit == unit[:shorter_length] * (unit_length // shorter_length):
            return True
    return False


def strands_differ(unit):
    """Tell whether the two strands of a unit can be told apart: more G than C, or more C than G.

    Takes in:  a candidate unit.
    Does:      compares the count of G with the count of C. Every telomerase-made
               repeat has a G-rich strand, and that strand is what lets tidk say which
               way an array points. A unit with as many G as C (AAAT: none of either;
               AGAGCC: two of each) has no G-rich strand, so every array of it counts
               as a telomere whichever way it points, and the score cannot separate a
               telomere from a satellite.
    Produces:  True when the counts differ.
    Used by:   passes_candidate_filters().
    """
    return unit.count("G") != unit.count("C")


def passes_candidate_filters(unit):
    """Apply the three filters of step 3 of the method to one candidate unit.

    Takes in:  a unit in canonical form.
    Does:      keeps the unit when it is UNIT_MIN_LENGTH to UNIT_MAX_LENGTH bases long,
               is not a repeat of a shorter unit, and has two strands that differ.
    Produces:  True for a plausible telomere unit, False for one to drop.
    Used by:   collect_candidates(), for every candidate except the configured unit.
    """
    if not UNIT_MIN_LENGTH <= len(unit) <= UNIT_MAX_LENGTH:
        return False
    if is_repeat_of_shorter_unit(unit):
        return False
    if not strands_differ(unit):
        return False
    return True


# ======================================================================
# WHAT IS KNOWN: THE LITERATURE AND THE DATABASE
# ======================================================================

def literature_by_unit():
    """Index the LITERATURE table by canonical unit.

    Produces:  {canonical unit: list of (published spelling, fungi, references)}; a
               list, because two published spellings could be the same unit.
    Used by:   run(); the entries fill the in_literature column and the summary line.
    """
    entries_by_unit = {}
    for spelling, fungi, references in LITERATURE:
        unit = canonical_unit(spelling)
        entries_by_unit.setdefault(unit, []).append((spelling, fungi, references))
    return entries_by_unit


def in_literature_text(entries):
    """Write the in_literature cell: the published spelling, the fungi and the references.

    Takes in:  the entries of one unit from literature_by_unit().
    Produces:  for example "TTAGGG: filamentous Ascomycota, most Basidiomycota,
               Mucoromycota; Botrytis cinerea (Lue 2021; Sun et al. 2023; Levis et al.
               1997)"; several entries are separated by " | ".
    Used by:   collect_candidates().
    """
    parts = []
    for spelling, fungi, references in entries:
        parts.append(f"{spelling}: {fungi} ({references})")
    return " | ".join(parts)


def spelling_note(spelling, canonical):
    """Say how a published or configured spelling relates to the canonical form.

    Takes in:  a spelling (TTAGGG) and the canonical form of the same unit (AACCCT).
    Does:      nothing when the two are the same. Otherwise says whether the canonical
               form is the spelling read from another position of the repeat (a
               rotation on the same strand) or read on the other strand.
    Produces:  "", "TTAGGG read from another position" or "TTAGGG read on the other
               strand".
    Used by:   summary_line(), so that a reader who knows the unit as TTAGGG
               understands AACCCT.
    """
    if spelling == canonical:
        return ""
    if canonical in rotations(spelling):
        return f"{spelling} read from another position"
    return f"{spelling} read on the other strand"


def read_telomere_database(path):
    """Read the telomere repeat database of tidk: which species carry which unit.

    Takes in:  --telomere-db, telo_database_v2.tsv: tab-separated, a header line, one
               row per species with its taxonomy, its telomere_unit (in tidk's
               canonical form) and the support of that unit: strong (the unit was
               confirmed at the sequence ends), weak, unchecked or none (no unit).
    Does:      keeps the rows with strong or weak support whose unit is made of A, C,
               G and T, and counts the species per kingdom for every unit. A row
               without a kingdom is counted under NO_KINGDOM. The unit is put in
               canonical form once more, in case a later database version is not.
    Produces:  {canonical unit: {kingdom: set of species names}}.
    Used by:   run(); the counts fill the in_database column and the summary line.
    Stops with a clear error when one of the four columns is missing (a cut-off file,
    or another table altogether) or when no row has strong or weak support.
    """
    species_by_unit = {}
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as table:
        reader = csv.DictReader(table, delimiter="\t")
        found = set(reader.fieldnames or [])
        missing = sorted(DATABASE_COLUMNS - found)
        if len(missing) > 0:
            raise ValueError(
                f"{path}: not the telomere repeat database of tidk: missing the "
                f"column(s) {', '.join(missing)}"
            )
        for row in reader:
            support = (row["support"] or "").strip()
            if support not in DATABASE_SUPPORT:
                continue
            unit = (row["telomere_unit"] or "").strip().upper()
            if not re.fullmatch(r"[ACGT]+", unit):
                continue
            kingdom = (row["kingdom"] or "").strip() or NO_KINGDOM
            species = (row["species"] or "").strip()
            unit = canonical_unit(unit)
            species_by_kingdom = species_by_unit.setdefault(unit, {})
            species_by_kingdom.setdefault(kingdom, set()).add(species)
    if len(species_by_unit) == 0:
        raise ValueError(f"{path}: no row with strong or weak support in the telomere database")
    return species_by_unit


def in_database_text(species_by_kingdom):
    """Write the in_database cell: the species with strong or weak support, per kingdom.

    Takes in:  {kingdom: set of species} of one unit, from read_telomere_database().
    Produces:  for example "Fungi 20; Metazoa 420; Viridiplantae 11", kingdoms in
               alphabetical order, the rows without a kingdom last.
    Used by:   collect_candidates().
    """
    def named_kingdoms_first(kingdom_and_species):
        # Sort key: False sorts before True, so the named kingdoms come first, in
        # alphabetical order, and the no-kingdom label last.
        kingdom, _ = kingdom_and_species
        return (kingdom == NO_KINGDOM, kingdom)

    parts = []
    for kingdom, species in sorted(species_by_kingdom.items(), key=named_kingdoms_first):
        parts.append(f"{kingdom} {len(species)}")
    return "; ".join(parts)


# ======================================================================
# RUNNING TIDK AND READING ITS TABLES
# ======================================================================

def tidk_error_text(stderr):
    """Pick the message out of tidk's error stream, for the error line of this script.

    Takes in:  what tidk wrote to its error stream.
    Does:      drops the progress lines ("[+] ...") and the warning that no clade
               database was fetched (tidk prints it on every run; this script never
               uses that database), and joins the rest.
    Produces:  one line of text.
    Used by:   run_tidk().
    """
    lines = []
    for line in stderr.splitlines():
        text = line.strip()
        if text == "" or text.startswith("[+]") or text.startswith("Warning! No clades"):
            continue
        lines.append(text)
    if len(lines) == 0:
        return "(no message from tidk)"
    return " ".join(lines)


def run_tidk(tidk, arguments):
    """Run one tidk command and stop with a readable message when it fails.

    Takes in:  the tidk program (--tidk) and the list of words after it.
    Does:      runs the command, capturing what it prints. A program that cannot be
               found, or an exit status other than 0, becomes a ValueError whose text
               names the command and carries tidk's own message.
    Produces:  what tidk wrote to standard output (tidk explore prints its table
               there; tidk ends writes files and prints nothing of use).
    Used by:   tidk_explore() and tidk_ends().
    """
    command = [tidk] + arguments
    try:
        finished = subprocess.run(command, capture_output=True, text=True)
    except FileNotFoundError:
        raise ValueError(
            f"the tidk program {tidk!r} was not found: install tidk 0.3.0 or name the "
            "program with --tidk"
        ) from None
    if finished.returncode != 0:
        raise ValueError(
            f"tidk {arguments[0]} failed (exit status {finished.returncode}) for the "
            f"command '{' '.join(command)}': {tidk_error_text(finished.stderr)}"
        )
    return finished.stdout


def parse_explore_table(text):
    """Read the table tidk explore prints: one row per repeat unit with its copies.

    Takes in:  tidk explore's standard output: possibly a few message lines, then the
               header line canonical_repeat_unit, copies, copies_as_unit,
               copies_as_revcomp, then one row per unit, most copies first.
    Does:      skips everything before the header, checks the header, and reads the
               unit and its copies (the sum of both strands) from every row. A unit
               with a letter other than A, C, G and T (from an N in the sequence)
               cannot be searched and is left out.
    Produces:  an ordered dictionary {canonical unit: copies}, in tidk's order.
    Used by:   tidk_explore().
    Stops with a clear error when the header is not tidk 0.3.0's or a copy number is
    not a whole number.
    """
    copies_by_unit = OrderedDict()
    header_seen = False
    for line in text.splitlines():
        fields = line.rstrip("\n").split("\t")
        # The lines before the header are tidk's messages: skip them.
        if not header_seen:
            if fields[0].strip() == EXPLORE_COLUMNS[0]:
                if fields != EXPLORE_COLUMNS:
                    raise ValueError(
                        f"tidk explore wrote the columns {fields}, not the expected "
                        f"{EXPLORE_COLUMNS}: another tidk version?"
                    )
                header_seen = True
            continue
        if line.strip() == "":
            continue
        if len(fields) < 2:
            raise ValueError(f"tidk explore wrote a row without a copy number: {line!r}")
        unit = fields[0].strip().upper()
        if not re.fullmatch(r"[ACGT]+", unit):
            continue
        try:
            copies_by_unit[unit] = int(fields[1])
        except ValueError:
            raise ValueError(
                f"tidk explore wrote the copy number {fields[1]!r} for {unit}, not a whole number"
            ) from None
    if not header_seen:
        raise ValueError("tidk explore printed no table (no header line canonical_repeat_unit)")
    return copies_by_unit


def tidk_explore(tidk, fasta, min_copies):
    """List the tandem repeat units of a FASTA file with tidk explore.

    Takes in:  the tidk program, the FASTA to search (the contig-end windows, or the
               whole draft) and C (--min-telomere-copies).
    Does:      runs tidk explore for unit lengths UNIT_MIN_LENGTH to UNIT_MAX_LENGTH
               with threshold C (a unit is listed when it forms a tandem run of at
               least C copies somewhere) over the whole of every sequence
               (EXPLORE_DISTANCE), and reads the table it prints.
    Produces:  {canonical unit: copies}.
    Used by:   run(), step 2 of the method, twice.
    """
    output = run_tidk(tidk, [
        "explore",
        "--minimum", str(UNIT_MIN_LENGTH),
        "--maximum", str(UNIT_MAX_LENGTH),
        "--threshold", str(min_copies),
        "--distance", str(EXPLORE_DISTANCE),
        fasta,
    ])
    return parse_explore_table(output)


def ends_table_path(folder, unit):
    """The path of the per-end table tidk ends writes for one unit in the temporary folder."""
    return os.path.join(folder, unit + ".ends.tsv")


def read_ends_table(path):
    """Read the per-end table tidk ends wrote: two rows per sequence, start and end.

    Takes in:  the path of a .ends.tsv file: the header ENDS_COLUMNS, then one row per
               sequence end with the status, the 0-based coordinates of the array on
               the sequence (telomere_start, telomere_end), its length, its offset
               (distance from that end, in bp), the strand (C-rich or G-rich) and
               window_limited (true when the array reaches the inner edge of the
               window tidk searched, so it may be longer than reported); NA in every
               column after the status when no array was found.
    Does:      checks the header and reads every row into a dictionary column -> text.
               The values stay text, so that telomere_ends.tsv can carry them exactly
               as tidk wrote them.
    Produces:  a list of row dictionaries, in tidk's order.
    Used by:   tidk_ends(), and run() for the chosen unit's table.
    Stops with a clear error when the header is not tidk 0.3.0's or a row is short.
    """
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as table:
        reader = csv.reader(table, delimiter="\t")
        header = next(reader, None)
        if header != ENDS_COLUMNS:
            raise ValueError(
                f"{path}: tidk ends wrote the columns {header}, not the expected "
                f"{ENDS_COLUMNS}: another tidk version?"
            )
        for fields in reader:
            if len(fields) == 0:
                continue
            if len(fields) != len(ENDS_COLUMNS):
                raise ValueError(
                    f"{path}: a row has {len(fields)} columns instead of {len(ENDS_COLUMNS)}: {fields}"
                )
            # Pair each column name with its cell: {"id": ..., "sequence_length": ..., ...}.
            rows.append(OrderedDict(zip(ENDS_COLUMNS, fields)))
    return rows


def tidk_ends(tidk, fasta, unit, window, min_copies, folder):
    """Call the telomeres of one unit at both ends of every contig with tidk ends.

    Takes in:  the tidk program, the draft assembly, the unit (canonical form), W and
               C, and the temporary folder for tidk's files.
    Does:      runs tidk ends with window W, maximum offset W (an array anywhere in
               the window is terminal enough: the assembler may have run past the
               telomere) and minimum array length C x unit length, in the default
               error-tolerant mode. tidk writes <folder>/<unit>.ends.tsv, .ends.json
               and .bed; the table is read back.
    Produces:  the rows of read_ends_table().
    Used by:   run(), step 4 of the method, once per candidate.
    """
    run_tidk(tidk, [
        "ends",
        "--string", unit,
        "--window", str(window),
        "--min-length", str(min_copies * len(unit)),
        "--max-offset", str(window),
        "--output", unit,
        "--dir", folder,
        fasta,
    ])
    return read_ends_table(ends_table_path(folder, unit))


# ======================================================================
# THE SAME-ARRAY RULE AND THE COUNTS PER UNIT
# ======================================================================

def coordinate(row, column):
    """Return one coordinate of a row that reports an array as a whole number.

    Takes in:  a row of the per-end table whose status is not absent, and the column
               (telomere_start, telomere_end or offset).
    Produces:  the number.
    Used by:   intervals_overlap() and mark_arrays_seen_from_both_ends().
    Stops with a clear error when the cell is not a number (NA in a row with an array
    would mean a tidk version that writes the table differently).
    """
    text = row[column]
    try:
        return int(text)
    except ValueError:
        raise ValueError(
            f"tidk ends wrote {text!r} as the {column} of {row['id']} ({row['side']}), "
            f"status {row['status']}: not a number"
        ) from None


def intervals_overlap(row_a, row_b):
    """Tell whether the arrays of two rows share at least one base of the contig."""
    a_before_b_ends = coordinate(row_a, "telomere_start") < coordinate(row_b, "telomere_end")
    b_before_a_ends = coordinate(row_b, "telomere_start") < coordinate(row_a, "telomere_end")
    return a_before_b_ends and b_before_a_ends


def mark_arrays_seen_from_both_ends(rows):
    """The same-array rule: one array seen from both ends of a short contig counts once.

    Takes in:  the rows of one tidk ends run (two per contig).
    Does:      tidk looks W bases in from each end, so a contig shorter than 2W is
               looked through twice, and one array near its start is reported by the
               start row (present, offset 0, say) and again by the end row (the same
               coordinates, an offset of some kb, and the opposite status, because
               the strand that is right for a start is wrong for an end). When both
               rows of a contig report an array and the two coordinate intervals
               overlap, the row with the larger offset (the far view of the array)
               gets the status other_end; on equal offsets the end row does. The
               rows are changed in place.
    Produces:  the same list.
    Used by:   run(), after every tidk ends run and for the chosen unit's table.
    """
    rows_by_contig = OrderedDict()
    for row in rows:
        rows_by_contig.setdefault(row["id"], {})[row["side"]] = row

    for sides in rows_by_contig.values():
        start_row = sides.get("start")
        end_row = sides.get("end")
        if start_row is None or end_row is None:
            continue
        if start_row["status"] == STATUS_ABSENT or end_row["status"] == STATUS_ABSENT:
            continue
        if not intervals_overlap(start_row, end_row):
            continue
        if coordinate(end_row, "offset") >= coordinate(start_row, "offset"):
            end_row["status"] = STATUS_OTHER_END
        else:
            start_row["status"] = STATUS_OTHER_END
    return rows


def count_ends(rows):
    """Count what the per-end table says about one unit, after the same-array rule.

    Takes in:  the rows of one unit.
    Does:      counts the ends with status present (the array points the way a
               telomere points: G-rich strand towards the contig end), the ends with
               status wrong_strand, and the contigs with present at both ends.
    Produces:  {"telomere_ends", "opposite_ends", "score", "t2t_contigs"}, score being
               telomere_ends minus opposite_ends.
    Used by:   run(), once per candidate.
    """
    telomere_ends = 0
    opposite_ends = 0
    present_sides = {}
    for row in rows:
        if row["status"] == STATUS_PRESENT:
            telomere_ends += 1
            present_sides.setdefault(row["id"], set()).add(row["side"])
        elif row["status"] == STATUS_WRONG_STRAND:
            opposite_ends += 1

    t2t_contigs = 0
    for sides in present_sides.values():
        if "start" in sides and "end" in sides:
            t2t_contigs += 1

    return {
        "telomere_ends": telomere_ends,
        "opposite_ends": opposite_ends,
        "score": telomere_ends - opposite_ends,
        "t2t_contigs": t2t_contigs,
    }


def formed_array(rows):
    """Tell whether tidk reported an array of the unit at any contig end (any status but absent)."""
    for row in rows:
        if row["status"] != STATUS_ABSENT:
            return True
    return False


# ======================================================================
# THE CANDIDATES, THEIR RANKING AND THE CHOICE
# ======================================================================

def collect_candidates(units_at_ends, database, literature, configured_unit):
    """Gather the candidate units from the four sources and apply the filters.

    Takes in:  the units tidk explore found at the contig ends ({unit: copies}), the
               database ({unit: {kingdom: species}}), the literature ({unit: entries})
               and the configured unit in canonical form (None when the motif is auto).
    Does:      takes every unit of the four sources once, in alphabetical order (so
               that tidk is run in the same order on every run), notes where it comes
               from and what is known about it, and drops the units that fail
               passes_candidate_filters(), except the configured unit.
    Produces:  (candidates, dropped): an ordered dictionary {unit: candidate} and the
               number of units the filters dropped. A candidate is a dictionary with
               unit, found_at_ends, in_literature and in_database (the table cells, or
               NOT_KNOWN), literature_entries, database_species, known and
               configured; run() adds the copies and the counts.
    Used by:   run(), step 3 of the method.
    """
    # Every unit of the three sources, each once (| joins sets).
    all_units = set(units_at_ends) | set(database) | set(literature)
    if configured_unit is not None:
        all_units.add(configured_unit)

    candidates = OrderedDict()
    dropped = 0
    for unit in sorted(all_units):
        candidate = {
            "unit": unit,
            "found_at_ends": unit in units_at_ends,
            "configured": unit == configured_unit,
            "in_literature": NOT_KNOWN,
            "in_database": NOT_KNOWN,
            "literature_entries": literature.get(unit, []),
            "database_species": database.get(unit, {}),
        }
        if unit in literature:
            candidate["in_literature"] = in_literature_text(literature[unit])
        if unit in database:
            candidate["in_database"] = in_database_text(database[unit])
        candidate["known"] = unit in literature or unit in database

        if not candidate["configured"] and not passes_candidate_filters(unit):
            dropped += 1
            continue
        candidates[unit] = candidate
    return candidates, dropped


def is_eligible(candidate):
    """Tell whether a scored candidate can be the telomere unit of the assembly.

    Takes in:  a candidate with its counts.
    Does:      asks for at least MIN_TELOMERE_ENDS telomere ends, and for more than
               OPPOSITE_ENDS_FACTOR times as many telomere ends as opposite ends.
    Produces:  True or False.
    Used by:   first_eligible().
    """
    enough_ends = candidate["telomere_ends"] >= MIN_TELOMERE_ENDS
    one_direction = candidate["telomere_ends"] > OPPOSITE_ENDS_FACTOR * candidate["opposite_ends"]
    return enough_ends and one_direction


def rank_candidates(candidates):
    """Order the scored candidates, best first.

    Takes in:  the candidates (a list of dictionaries with their counts).
    Does:      sorts by score (highest first), then known units before unknown ones,
               then the shorter unit, then alphabetically, so that the order is the
               same on every run.
    Produces:  a new, sorted list.
    Used by:   run().
    """
    def best_first(candidate):
        # Sort key: a tuple sorts by its first item and uses the next ones on ties; the
        # minus sign puts the highest score first.
        known_rank = 0 if candidate["known"] else 1
        return (-candidate["score"], known_rank, len(candidate["unit"]), candidate["unit"])

    return sorted(candidates, key=best_first)


def first_eligible(ranked):
    """Return the first eligible candidate of the ranking, or None when there is none."""
    for candidate in ranked:
        if is_eligible(candidate):
            return candidate
    return None


# ======================================================================
# WRITE THE OUTPUT FILES
# ======================================================================

def write_candidates_table(path, candidates):
    """Write telomere_candidates.tsv.

    Takes in:  the path and the scored candidates (with their chosen flag).
    Does:      writes one row per candidate that was found at the contig ends, formed
               at least one array, or is the configured unit; the other known units
               are left out (the log counts them). Sorted by score, then
               telomere_ends, both descending, then unit.
    Produces:  the file; returns the number of rows written.
    Used by:   run(). A record for the reader: it shows at once how every unit behaves
               at the contig ends and which one was chosen.
    """
    def table_order(candidate):
        # Sort key: highest score first, then most telomere ends, then the unit.
        return (-candidate["score"], -candidate["telomere_ends"], candidate["unit"])

    rows = []
    for candidate in candidates:
        if candidate["found_at_ends"] or candidate["formed_array"] or candidate["configured"]:
            rows.append(candidate)

    with open(path, "w", encoding="utf-8") as table:
        table.write("\t".join(CANDIDATE_COLUMNS) + "\n")
        for candidate in sorted(rows, key=table_order):
            cells = [
                candidate["unit"],
                str(len(candidate["unit"])),
                str(candidate["copies_at_ends"]),
                str(candidate["copies_genome_wide"]),
                str(candidate["telomere_ends"]),
                str(candidate["opposite_ends"]),
                str(candidate["score"]),
                str(candidate["t2t_contigs"]),
                candidate["in_literature"],
                candidate["in_database"],
                "yes" if candidate["chosen"] else "no",
            ]
            table.write("\t".join(cells) + "\n")
    return len(rows)


def write_unit_file(path, unit):
    """Write telomere_unit.txt: one line, the unit in canonical form or none.

    Used by:   run(). Rule contig_summary reads it.
    """
    with open(path, "w", encoding="utf-8") as unit_file:
        unit_file.write(unit + "\n")


def write_ends_table(path, rows):
    """Write telomere_ends.tsv: the per-end table of the chosen unit, after the same-array rule.

    Takes in:  the path and the rows (empty when the unit is none: the header alone).
    Does:      writes tidk's ten columns and values, the status possibly other_end.
    Used by:   run(). Rule contig_summary reads it.
    """
    with open(path, "w", encoding="utf-8") as table:
        table.write("\t".join(ENDS_COLUMNS) + "\n")
        for row in rows:
            table.write("\t".join(row[column] for column in ENDS_COLUMNS) + "\n")


# ======================================================================
# THE LINES FOR THE LOG AND THE CONSOLE
# ======================================================================

def known_text(candidate):
    """Say where a unit is known from, for the summary line.

    Takes in:  a candidate.
    Produces:  "known from the literature (<fungi>) and from 20 fungal species in the
               database", either half alone, "known from the database (Metazoa 10),
               from no fungus" when only other kingdoms carry it, or "not known from
               the literature or the database".
    Used by:   summary_line().
    """
    parts = []
    if len(candidate["literature_entries"]) > 0:
        fungi = []
        for _, entry_fungi, _ in candidate["literature_entries"]:
            fungi.append(entry_fungi)
        parts.append(f"the literature ({'; '.join(fungi)})")
    if len(candidate["database_species"]) > 0:
        fungal_species = len(candidate["database_species"].get("Fungi", set()))
        if fungal_species > 0:
            parts.append(f"{fungal_species} fungal species in the database")
        else:
            parts.append(f"the database ({candidate['in_database']}), from no fungus")
    if len(parts) == 0:
        return "not known from the literature or the database"
    return "known from " + " and from ".join(parts)


def summary_line(sample, chosen, best, configured_spelling):
    """Write the one-line summary that the rule copies to the console.

    Takes in:  the sample name; the chosen candidate (None when the unit is none); the
               best-ranked candidate among those that formed an array (None when no
               candidate did), named when nothing was chosen; and the spelling the
               user configured (None when the motif is auto).
    Produces:  "[FunFlux] telomere_units <sample>: telomere repeat unit AACCCT (TTAGGG
               read on the other strand): 34 contig ends, 15 contigs with a telomere at
               both ends; known from ..." or "[FunFlux] telomere_units <sample>: no
               telomere repeat identified (best candidate X: n ends in the telomere
               orientation, m opposite)". The part in brackets after the unit is given
               when a literature or configured spelling differs from the canonical
               form, and says "set by telomere_motif" for a configured unit.
    Used by:   run().
    """
    prefix = f"[FunFlux] telomere_units {sample}: "
    if chosen is None:
        if best is None:
            return prefix + (
                "no telomere repeat identified (no candidate unit formed an array at a contig end)"
            )
        return prefix + (
            f"no telomere repeat identified (best candidate {best['unit']}: "
            f"{best['telomere_ends']} ends in the telomere orientation, "
            f"{best['opposite_ends']} opposite)"
        )

    unit = chosen["unit"]
    remarks = []
    spelling = configured_spelling
    if len(chosen["literature_entries"]) > 0:
        spelling = chosen["literature_entries"][0][0]
    if spelling is not None:
        note = spelling_note(spelling, unit)
        if note != "":
            remarks.append(note)
    if chosen["configured"]:
        remarks.append("set by telomere_motif")
    unit_text = unit
    if len(remarks) > 0:
        unit_text += f" ({'; '.join(remarks)})"

    return prefix + (
        f"telomere repeat unit {unit_text}: {chosen['telomere_ends']} contig ends, "
        f"{chosen['t2t_contigs']} contigs with a telomere at both ends; {known_text(chosen)}"
    )


def note_lines(sample, chosen, data_choice, configured_spelling):
    """Write the NOTE lines, when there is something the user should know.

    Takes in:  the sample name, the chosen candidate (None when the unit is none), the
               data's first choice (the first eligible candidate, None when there is
               none) and the configured spelling (None when the motif is auto).
    Does:      one NOTE when the chosen unit is unknown to the literature and the
               database; one NOTE when a configured unit is not the data's first
               choice, naming the unit the data favour or saying that the configured
               unit has no support.
    Produces:  a list of lines starting "[FunFlux] NOTE telomere_units <sample>: ",
               possibly empty.
    Used by:   run().
    """
    prefix = f"[FunFlux] NOTE telomere_units {sample}: "
    notes = []
    if chosen is None:
        return notes

    if not chosen["known"]:
        notes.append(prefix + (
            f"the telomere repeat unit {chosen['unit']} is not known from the literature "
            "or the telomere database: look at the arrays in telomere_ends.tsv before "
            "trusting the telomere calls"
        ))

    if chosen["configured"] and data_choice is not chosen:
        configured_text = f"the configured unit {configured_spelling} ({chosen['unit']})"
        own_counts = (
            f"{chosen['telomere_ends']} contig ends in the telomere orientation and "
            f"{chosen['opposite_ends']} opposite"
        )
        if data_choice is None:
            notes.append(prefix + (
                f"{configured_text} has {own_counts}, and no other unit is supported "
                "either; set telomere_motif to auto to let the data decide, or look at "
                "telomere_candidates.tsv"
            ))
        else:
            known = "known" if data_choice["known"] else "not known from the literature or the database"
            notes.append(prefix + (
                f"the data favour {data_choice['unit']} ({data_choice['telomere_ends']} contig "
                f"ends in the telomere orientation, {data_choice['opposite_ends']} opposite; "
                f"{known}) over {configured_text}, which has {own_counts}; set "
                "telomere_motif to auto or to that unit and re-run"
            ))
    return notes


# ======================================================================
# PUTTING THE STEPS TOGETHER
# ======================================================================

def run(args):
    """Do the whole job: the draft and the database in, the three files out.

    Takes in:  args, the parsed and checked command-line options.
    Does:      the steps of THE METHOD, with tidk's files in a temporary folder that
               is removed at the end. Every input is read before tidk runs, and the
               output files are written last, so bad input never leaves a
               half-written file behind.
    Produces:  the three output files and the lines on standard error. Raises OSError
               or ValueError on bad input or a tidk failure; main() turns that into a
               message.
    Used by:   main().
    """
    started = time.monotonic()

    # The inputs: the draft contigs, the database, the literature, the configured unit.
    contigs = read_fasta(args.contigs)
    database = read_telomere_database(args.telomere_db)
    literature = literature_by_unit()
    configured_unit = None
    configured_spelling = None
    if args.telomere_motif != "auto":
        configured_spelling = args.telomere_motif
        configured_unit = canonical_unit(args.telomere_motif)

    with tempfile.TemporaryDirectory(prefix="telomere_units_") as folder:
        # Step 1 and 2: the contig-end windows, and the units tidk explore finds in
        # them and in the whole draft.
        ends_fasta = os.path.join(folder, "contig_ends.fasta")
        write_contig_ends(contigs, args.telomere_window, ends_fasta)
        copies_at_ends = tidk_explore(args.tidk, ends_fasta, args.min_telomere_copies)
        copies_genome_wide = tidk_explore(args.tidk, args.contigs, args.min_telomere_copies)

        # Step 3: the candidates.
        candidates, dropped = collect_candidates(copies_at_ends, database, literature, configured_unit)

        # Step 4: tidk ends for every candidate, the same-array rule, the counts.
        for unit, candidate in candidates.items():
            rows = tidk_ends(
                args.tidk, args.contigs, unit, args.telomere_window, args.min_telomere_copies, folder
            )
            candidate["formed_array"] = formed_array(rows)
            mark_arrays_seen_from_both_ends(rows)
            candidate.update(count_ends(rows))
            candidate["copies_at_ends"] = copies_at_ends.get(unit, 0)
            candidate["copies_genome_wide"] = copies_genome_wide.get(unit, 0)

        # Step 5: the ranking and the choice. A configured unit is the chosen one
        # whatever its rank; the data's first choice is still worked out for the NOTE.
        ranked = rank_candidates(list(candidates.values()))
        data_choice = first_eligible(ranked)
        if configured_unit is not None:
            chosen = candidates[configured_unit]
        else:
            chosen = data_choice
        # "is" asks whether it is the very same candidate, not one with the same values.
        for candidate in ranked:
            candidate["chosen"] = candidate is chosen
        # The best-ranked candidate that formed an array at all: named in the summary
        # when nothing was chosen (a unit without any array has nothing to be named for).
        best = None
        for candidate in ranked:
            if candidate["formed_array"]:
                best = candidate
                break

        # The per-end table of the chosen unit: tidk's table for it, read back from the
        # temporary folder, with the same-array rule applied.
        chosen_rows = []
        if chosen is not None:
            chosen_rows = read_ends_table(ends_table_path(folder, chosen["unit"]))
            mark_arrays_seen_from_both_ends(chosen_rows)

        # The output files.
        table_rows = write_candidates_table(args.out_candidates, ranked)
        write_unit_file(args.out_unit, chosen["unit"] if chosen is not None else NO_UNIT)
        write_ends_table(args.out_ends, chosen_rows)

    # The lines for the log: the counts, the summary, the NOTEs.
    left_out = len(ranked) - table_rows
    seconds = time.monotonic() - started
    print(
        f"telomere_units.py: sample {args.sample}: {len(contigs)} contigs; tidk explore "
        f"listed {len(copies_at_ends)} units at the contig ends and {len(copies_genome_wide)} "
        f"genome-wide; {len(ranked)} candidate units scored with tidk ends ({dropped} dropped "
        f"by the filters before scoring); {table_rows} rows in the candidates table, "
        f"{left_out} known units formed no array and are left out of it; {seconds:.1f} s.",
        file=sys.stderr,
    )
    print(summary_line(args.sample, chosen, best, configured_spelling), file=sys.stderr)
    for line in note_lines(args.sample, chosen, data_choice, configured_spelling):
        print(line, file=sys.stderr)


# ======================================================================
# COMMAND LINE
# ======================================================================

def parse_arguments(argv):
    """Read the options given on the command line.

    Takes in:  argv, the list of words after the script name (None means "use the real
               command line"; the tests pass their own list).
    Does:      defines the options with argparse (Python's built-in reader of
               command-line options). The limits of the three settings are checked in
               check_options(), so that they give the same kind of message as a
               missing file.
    Produces:  an object whose attributes are the options (dashes become underscores:
               args.out_unit, args.telomere_motif, ...).
    Used by:   main().
    """
    parser = argparse.ArgumentParser(
        prog="telomere_units.py",
        description=(
            "Choose the telomere repeat unit of a draft assembly from the data with tidk, "
            "and call the telomeres at both ends of every contig."
        ),
    )
    parser.add_argument("--contigs", required=True, help="the draft assembly (FASTA)")
    parser.add_argument("--telomere-db", required=True, help="the telomere repeat database of tidk (TSV)")
    parser.add_argument(
        "--telomere-motif", default=DEFAULT_TELOMERE_MOTIF,
        help="auto, or the telomere repeat unit of the species, letters ACGT only (default: %(default)s)",
    )
    parser.add_argument(
        "--telomere-window", type=int, default=DEFAULT_TELOMERE_WINDOW,
        help="bases searched at each contig end (default: %(default)s)",
    )
    parser.add_argument(
        "--min-telomere-copies", type=int, default=DEFAULT_MIN_TELOMERE_COPIES,
        help="tandem copies of the unit from which a telomere is called (default: %(default)s)",
    )
    parser.add_argument("--sample", required=True, help="the sample name, for the messages")
    parser.add_argument("--out-candidates", required=True, help="output: telomere_candidates.tsv")
    parser.add_argument("--out-unit", required=True, help="output: telomere_unit.txt")
    parser.add_argument("--out-ends", required=True, help="output: telomere_ends.tsv")
    parser.add_argument(
        "--tidk", default=DEFAULT_TIDK, help="the tidk program (default: %(default)s on PATH)",
    )
    return parser.parse_args(argv)


def check_options(args):
    """Hold the three settings to the limits the workflow applies to the configuration.

    Takes in:  the parsed options; the motif is normalised in place (auto in any case
               becomes "auto", a unit is upper-cased).
    Does:      stops with a ValueError, which main() turns into the error line, when
               the motif is neither auto nor a unit of MOTIF_MIN_LENGTH to
               MOTIF_MAX_LENGTH letters A, C, G and T, when the window is outside
               TELOMERE_WINDOW_MIN to TELOMERE_WINDOW_MAX, or when the copies are below
               MIN_TELOMERE_COPIES_MIN.
    Used by:   main().
    """
    motif = args.telomere_motif.strip()
    if motif.lower() == "auto":
        args.telomere_motif = "auto"
    else:
        args.telomere_motif = motif.upper()
        is_dna = re.fullmatch(r"[ACGT]+", args.telomere_motif) is not None
        right_length = MOTIF_MIN_LENGTH <= len(args.telomere_motif) <= MOTIF_MAX_LENGTH
        if not (is_dna and right_length):
            raise ValueError(
                "--telomere-motif must be auto (the unit is chosen from the data) or the "
                f"telomere repeat unit of the species, {MOTIF_MIN_LENGTH} to {MOTIF_MAX_LENGTH} "
                f"letters A, C, G and T, such as TTAGGG (got: {motif!r})"
            )
    if not TELOMERE_WINDOW_MIN <= args.telomere_window <= TELOMERE_WINDOW_MAX:
        raise ValueError(
            f"--telomere-window must be a whole number from {TELOMERE_WINDOW_MIN} to "
            f"{TELOMERE_WINDOW_MAX} (got: {args.telomere_window})"
        )
    if args.min_telomere_copies < MIN_TELOMERE_COPIES_MIN:
        raise ValueError(
            f"--min-telomere-copies must be a whole number of at least {MIN_TELOMERE_COPIES_MIN} "
            f"(got: {args.min_telomere_copies})"
        )


def main(argv=None):
    """Run the script from the command line.

    Takes in:  argv, the command-line words (None = the real command line).
    Does:      parses and checks the options and calls run(). A problem with an input,
               a setting or tidk becomes one readable line on standard error instead
               of a Python traceback, and exit status 1, so Snakemake sees that the
               step failed.
    Produces:  exit status 0 on success, 1 on a problem, 2 on a command line argparse
               cannot read (argparse's own code).
    Used by:   the last lines of this file when rule telomere_units runs the script,
               and by the tests, which pass their own argv.
    """
    args = parse_arguments(argv)
    try:
        check_options(args)
        run(args)
    except (OSError, ValueError) as problem:
        # OSError: a file cannot be opened or written (for example it does not exist).
        # ValueError: a setting outside its limits, a database without the columns, a
        # FASTA with two records of one name, a tidk failure or a tidk table that does
        # not look like tidk 0.3.0's.
        print(f"[FunFlux] telomere_units: error: {problem}", file=sys.stderr)
        sys.exit(1)


# Run main() only when this file is run as a script, not when the tests import it.
if __name__ == "__main__":
    main()
