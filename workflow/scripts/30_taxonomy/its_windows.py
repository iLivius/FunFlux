#!/usr/bin/env python3
"""Cut a window around every 5.8S gene so that ITSx sees one rDNA unit at a time.

BIOLOGY AND WHY THIS SCRIPT EXISTS
----------------------------------
The fungal ITS barcode (ITS1, 5.8S, ITS2) sits inside the ribosomal DNA (rDNA)
repeat. One repeat unit reads, in this order:

    SSU - ITS1 - 5.8S - ITS2 - LSU

and the unit is repeated head to tail, tens to hundreds of times (a "tandem array";
14 to 1,442 copies, median 82, estimated in 91 fungal genomes by Lofgren et al. 2019,
Molecular Ecology 28: 721-730). The workflow extracts the ITS from the assembled
genome with the program ITSx and classifies it against UNITE (a curated reference
database of fungal ITS sequences).

On long-read assemblies the whole rDNA array can sit on one multi-megabase contig,
often on the MINUS strand: the units then read along the reverse complement of the
contig as it is stored in the FASTA file. Given such a whole contig, ITSx reports
"domains found in wrong order, may be chimeric" and writes nothing, whereas the same
contig reverse-complemented works. So instead of whole contigs we hand ITSx short
windows, each holding ONE rDNA unit in the standard orientation (SSU first, LSU
last); v1 gave ITSx the whole genome. docs/methods/its-long-read-contigs.md
describes the design.

WHERE THIS SCRIPT FITS IN THE WORKFLOW
--------------------------------------
    1. Rule its_windows (workflow/rules/shared/30_taxonomy.smk) runs nhmmer (the DNA
       search program of the HMMER package) against the genome, contigs_final.fasta,
       with ITSx's own fungal 5.8S profiles. A profile is a statistical model of a
       sequence family. ITSx supplies short ones for the start and for the end of
       the 5.8S gene, so one gene is hit many times, by different profiles. nhmmer
       writes a hit table (nhmmer --tblout) with an E-value for each hit: the
       number of hits this good that would be expected by chance in a search of
       this size, so smaller is better. The rule runs nhmmer with -E 0.001, which
       already leaves out every hit with an E-value above 0.001. One 5.8S gene
       still gets many hits with an E-value of at most 1e-5, from different
       profiles.
    2. THIS SCRIPT, run next by the same rule, reads that table and the genome FASTA
       and cuts one window for every 5.8S gene.
    3. Rule its_extraction runs ITSx on the windows. Rule its_classification then
       collapses the identical ITS copies (its_distinct.py) and classifies them
       against UNITE with SINTAX.

WHAT THIS SCRIPT DOES
---------------------
    1. Reads the hit table and keeps the hits with an E-value of at most
       --max-evalue.
    2. Merges hits that lie close together on the same contig and strand into
       one "5.8S site": a hit joins the site when it starts at most
       --merge-distance bases after the end of the site built so far. This is
       needed because the start-anchor and end-anchor profiles hit the same gene
       many times.
    3. Cuts a window centred on the middle of each site, with --flank bases on
       each side, clipped at the contig ends.
    4. Reverse-complements the window when the site is on the minus strand, so
       that every window reads SSU - ITS1 - 5.8S - ITS2 - LSU.
    5. Writes the windows as a FASTA file and as a table with one row per window.

INPUTS
------
    --tblout  The nhmmer table. Columns are separated by whitespace and comment
              lines start with #. Columns used (counting from 0): 0 contig name,
              8 envelope start, 9 envelope end (the stretch of the contig that
              nhmmer considers to hold the hit), 10 contig length ("sq len"),
              11 strand (+ or -), 12 E-value.
    --fasta   The genome (the contigs) that the table was made from.

OUTPUTS
-------
    --out-fasta  The windows. Header format: >contig:start-end_strand, for example
                 >ptg000001l:10001-14001_-
                 The coordinates are 1-based and inclusive (the first base of a
                 contig is base 1) and are clipped to the contig ends. The strand
                 (+ or -) follows an underscore and is NOT put in parentheses,
                 because ITSx mishandles parentheses in headers. Sequence lines
                 are wrapped at 80 characters. The workflow keeps this file as
                 its_windows.fasta, and rule its_extraction gives it to ITSx.
    --out-table  Tab-separated, one row per window, with the columns:
                 window_id (the header without the ">"), contig, start, end,
                 strand, anchor_hits (how many nhmmer hits were merged into the
                 site), best_evalue (the smallest E-value among them) and
                 window_length (in bases). The workflow keeps it as
                 its_windows.tsv, for people who want to see where the rDNA
                 copies lie. No rule reads it.

If the genome has no recognisable 5.8S, the script still succeeds: it writes an
empty FASTA and a table that holds only the header line, and prints one line to
stderr saying so. The workflow must not fail on such a genome: rule its_extraction
sees the empty FASTA and writes empty ITS files instead of running ITSx.
"""

# Only Python's standard library is used, so the Python 3.12 that the itsx conda
# environment of rule its_windows holds is all this script needs (see the note at
# the top of 30_taxonomy.smk).
import argparse
import sys


# ======================================================================
# SETTINGS
# ======================================================================

# Default values of the three tuning options. Each can be changed on the command line.
#   flank:          bases kept on each side of the middle of the 5.8S site. 2000 bases
#                   on each side comfortably hold the end of the SSU, all of ITS1 and
#                   ITS2 and the start of the LSU, which ITSx needs to recognise the
#                   unit (a whole fungal ITS1-5.8S-ITS2 region is 0.4 to 1.5 kb long).
#   max E-value:    hits with a larger (worse) E-value are ignored as chance matches.
#   merge distance: a hit that starts at most this many bases after the end of the
#                   site built so far, on the same contig and strand, joins that site.
# Rule its_windows passes only --flank (parameters.its.window_flank in the config,
# 2000 by default, at least 200); the other two keep the defaults below. The rule's
# nhmmer search (-E 0.001) already drops weaker hits, so the 1e-5 below is the
# stricter of the two limits.
DEFAULT_FLANK = 2000
DEFAULT_MAX_EVALUE = 1e-5
DEFAULT_MERGE_DISTANCE = 1000

# Sequence lines in the output FASTA file are wrapped at this many characters.
FASTA_LINE_WIDTH = 80

# A data line of the nhmmer table must have at least 13 columns, because the E-value,
# the last column read here, is column 12 (counting from 0). Real lines have 16 or
# more: the score, the bias and the description of the contig ("-" when it has none)
# follow the E-value, and a description with spaces in it adds further columns. None
# of these is read.
MIN_TBLOUT_COLUMNS = 13


# ======================================================================
# REVERSE COMPLEMENT
# ======================================================================

# The DNA letters and their partners on the opposite strand.
# Plain bases pair A<->T and C<->G. The IUPAC ambiguity codes (single letters that
# stand for several possible bases) pair up in the same way:
#   R (A or G)    <-> Y (C or T)
#   K (G or T)    <-> M (A or C)
#   B (C, G or T) <-> V (A, C or G)
#   D (A, G or T) <-> H (A, C or T)
# S (C or G), W (A or T) and N (any base) are their own partners.
BASES = "ACGTRYKMBVDHSWN"
PARTNERS = "TGCAYRMKVBHDSWN"

# The swap table that reverse_complement() applies to a window with str.translate(),
# which swaps every letter of a text in one go. str.maketrans() builds the table from
# two strings of equal length: each letter of the first string is replaced by the
# letter at the same position in the second string. Lowercase letters (soft-masked
# bases, i.e. bases an earlier tool flagged as repeats) get lowercase partners, so the
# masking pattern survives. A character that is not in the table (a gap "-", a star,
# any odd letter) is left unchanged.
COMPLEMENT_TABLE = str.maketrans(
    BASES + BASES.lower(),
    PARTNERS + PARTNERS.lower(),
)


def reverse_complement(sequence):
    """Return the DNA sequence as it reads on the opposite strand.

    Takes in:  a piece of DNA as text, in practice one window cut from a contig by
               cut_window().
    Does:      replaces every base by its partner (see COMPLEMENT_TABLE), then reads
               the result backwards, because the opposite strand runs in the
               opposite direction.
    Produces:  the reverse-complemented text: same length, same upper/lower case
               pattern, and any character that is not a DNA letter is left as it was.
    Used by:   cut_window(), for every window whose 5.8S site is on the minus strand.
    """
    complemented = sequence.translate(COMPLEMENT_TABLE)

    # [::-1] means "the same text read from the last character to the first".
    reversed_and_complemented = complemented[::-1]
    return reversed_and_complemented


# ======================================================================
# STEP 1: READ THE NHMMER TABLE
# ======================================================================

def parse_hit_line(line, line_number, tblout_path):
    """Read one data line of the nhmmer table and return the values we need from it.

    Takes in:  one line of the table (not a comment line), plus its line number and
               the file name, which are only used to write a helpful error message.
    Does:      splits the line into columns and picks out the contig name (column
               0), the envelope start and end (columns 8 and 9), the contig length
               (column 10), the strand (column 11) and the E-value (column 12).
               On the minus strand nhmmer lists the envelope backwards (start larger
               than end), so the two numbers are put in ascending order here.
    Produces:  contig, strand, start, end, evalue, contig_length (start <= end).
    Used by:   read_hits().
    Stops with a clear error if the line is too short, a number cannot be read, or
    the strand is neither + nor -.
    """
    columns = line.split()

    # The error messages below are f-strings (an f just before the opening quote):
    # the values between curly braces are filled in. !r prints the offending line
    # between quotes, so that stray spaces or odd characters in it become visible.
    if len(columns) < MIN_TBLOUT_COLUMNS:
        raise ValueError(
            f"{tblout_path}, line {line_number}: expected an nhmmer --tblout line with at "
            f"least {MIN_TBLOUT_COLUMNS} columns but found {len(columns)}: {line.strip()!r}"
        )

    contig = columns[0]
    strand = columns[11]

    # int() and float() turn the text of a column into a number and fail on text that
    # is not a number. That failure becomes a message that says where the problem is.
    # ("from None" hides Python's own, less helpful, message about the failed conversion.)
    try:
        envelope_from = int(columns[8])
        envelope_to = int(columns[9])
        contig_length = int(columns[10])
        evalue = float(columns[12])
    except ValueError:
        raise ValueError(
            f"{tblout_path}, line {line_number}: could not read envelope start, envelope "
            "end, contig length and E-value (columns 9, 10, 11 and 13, counting from 1) "
            f"as numbers: {line.strip()!r}"
        ) from None

    if strand not in ("+", "-"):
        raise ValueError(
            f"{tblout_path}, line {line_number}: the strand column (column 12, counting "
            f"from 1) must be + or - but is {strand!r}"
        )

    # Smaller coordinate first, whatever the strand.
    start = min(envelope_from, envelope_to)
    end = max(envelope_from, envelope_to)

    return contig, strand, start, end, evalue, contig_length


def read_hits(tblout_path, max_evalue):
    """Read the nhmmer table and keep only the 5.8S hits that are good enough.

    Takes in:  the table written by `nhmmer --tblout` in rule its_windows, and the
               worst E-value we still accept (--max-evalue).
    Does:      reads every data line (skipping comment lines starting with # and
               blank lines). A hit whose E-value is above max_evalue is counted and
               thrown away. The remaining hits are filed under their (contig, strand)
               pair, because only hits on the same contig AND the same strand can
               belong to the same 5.8S gene.
    Produces:  three things, in this order:
                 hits_by_contig_and_strand  a dictionary (a lookup table) that maps
                     (contig, strand) to a list of hits, each hit being the three
                     numbers (start, end, evalue);
                 contig_lengths_in_table  a dictionary mapping each contig that has
                     a good hit to the length the table states for it;
                 n_weak_hits  how many hits were thrown away for being too weak.
    Used by:   run(), which passes the hits on to merge_hits_into_sites() and the
               contig lengths on to check_contigs_match_table().
    """
    hits_by_contig_and_strand = {}
    contig_lengths_in_table = {}
    n_weak_hits = 0

    with open(tblout_path, "r", encoding="utf-8") as tblout:
        for line_number, line in enumerate(tblout, start=1):
            # Comment lines (the column titles at the top and the run summary at the
            # bottom) start with "#". Blank lines carry nothing either.
            if line.startswith("#") or line.strip() == "":
                continue

            contig, strand, start, end, evalue, contig_length = parse_hit_line(
                line, line_number, tblout_path
            )

            # This filter comes BEFORE the merging on purpose: a weak, chance hit that
            # happens to lie near a real gene must not stretch that gene's site.
            if evalue > max_evalue:
                n_weak_hits += 1
                continue

            key = (contig, strand)
            if key not in hits_by_contig_and_strand:
                hits_by_contig_and_strand[key] = []
            hits_by_contig_and_strand[key].append((start, end, evalue))
            contig_lengths_in_table[contig] = contig_length

    return hits_by_contig_and_strand, contig_lengths_in_table, n_weak_hits


# ======================================================================
# STEP 2: MERGE HITS THAT BELONG TO THE SAME 5.8S GENE
# ======================================================================

def merge_hits_into_sites(hits_by_contig_and_strand, merge_distance):
    """Merge the hits that belong to the same 5.8S gene into one "site".

    Takes in:  the good hits from read_hits(), and --merge-distance.
    Does:      the start-anchor and the end-anchor profiles of ITSx hit the same 5.8S
               gene many times over (the real tables hold dozens of hits per gene),
               but we want ONE window per gene. So, for each (contig, strand), the
               hits are put in order from left to right along the contig and each
               hit is compared with the site built so far: if the hit starts no
               more than merge_distance bases after the end of the site, it joins
               the site and stretches it if it reaches further; otherwise it begins
               a new site. Neighbouring rDNA units are several thousand bases apart
               (the length of one rDNA unit), far more than the default 1000, so
               they stay separate. Hits on different strands are never merged with
               each other.
    Produces:  a list of sites. Each site is a dictionary with the keys contig,
               strand, start, end, anchor_hits (how many hits were merged) and
               best_evalue (the smallest E-value among them).
    Used by:   run(), which then reads the contigs and cuts one window per site.
    """
    sites = []

    # The hits are grouped by strand BEFORE they are sorted and merged. Sorting both
    # strands together and comparing every hit only with the previous site would go
    # wrong: a plus hit, a minus hit and another plus hit close together would not all
    # be merged. With one group per (contig, strand), the two plus hits form one site
    # and the minus hit forms its own.
    for (contig, strand), hits in hits_by_contig_and_strand.items():
        # Each hit is the triple (start, end, evalue), and sorting a list of triples
        # orders them by start first. This puts the hits in order along the contig.
        # sorted() makes a sorted copy, so the caller's lists keep their original order.
        sorted_hits = sorted(hits)

        current_site = None
        for start, end, evalue in sorted_hits:
            # A hit joins the current site if there is one and the hit starts no more
            # than merge_distance bases after the end of the site. Overlapping hits
            # give a negative distance, so they always join.
            joins_current_site = (
                current_site is not None
                and start - current_site["end"] <= merge_distance
            )

            if joins_current_site:
                # Stretch the site if this hit reaches further, count the hit, and
                # keep the best (smallest) E-value seen so far.
                current_site["end"] = max(current_site["end"], end)
                current_site["anchor_hits"] += 1
                current_site["best_evalue"] = min(current_site["best_evalue"], evalue)
            else:
                # Too far from the site so far (or the first hit): a new site starts.
                current_site = {
                    "contig": contig,
                    "strand": strand,
                    "start": start,
                    "end": end,
                    "anchor_hits": 1,
                    "best_evalue": evalue,
                }
                sites.append(current_site)

    return sites


# ======================================================================
# STEP 3: READ THE CONTIGS THAT CARRY A SITE AND CHECK THEY MATCH THE TABLE
# ======================================================================

def read_contigs(fasta_path, wanted_contigs):
    """Read the genome FASTA, keeping only the contigs we need.

    Takes in:  the genome FASTA that nhmmer searched, and the set (a collection
               without repeats) of contig names that have a 5.8S site.
    Does:      goes through the file line by line. The name of a contig is the first
               word of its ">" header line, which is also the name nhmmer writes in
               its table. The sequence lines of contigs we do not need are skipped,
               so a genome of 100 Mb costs no more memory than the few contigs that
               carry the rDNA.
    Produces:  a dictionary mapping contig name to its full sequence, with the
               contigs in the same order as in the FASTA file. Contigs that are in
               wanted_contigs but not in the file are simply absent from it;
               check_contigs_match_table() reports them.
    Used by:   run(), which passes it to check_contigs_match_table() and then to
               cut_window().
    Stops with a clear error if a wanted contig name appears twice in the FASTA,
    because the coordinates in the table would then be ambiguous.
    """
    sequence_lines_by_contig = {}
    current_contig = None   # the wanted contig whose lines are being read, else None

    with open(fasta_path, "r", encoding="utf-8") as fasta:
        for line in fasta:
            if line.startswith(">"):
                # A new contig starts here. Its name is the first word of the header.
                words = line[1:].split()
                if len(words) > 0:
                    name = words[0]
                else:
                    name = ""   # a header without a name cannot match any hit

                if name in wanted_contigs:
                    if name in sequence_lines_by_contig:
                        raise ValueError(
                            f"{fasta_path}: the contig name {name!r} appears more than "
                            "once, so the coordinates in the nhmmer table are ambiguous"
                        )
                    sequence_lines_by_contig[name] = []
                    current_contig = name
                else:
                    current_contig = None
            elif current_contig is not None:
                # A sequence line of a contig we want. strip() removes the end of line.
                sequence_lines_by_contig[current_contig].append(line.strip())

    # Join the lines of each contig into one long piece of text.
    sequences = {}
    for contig, sequence_lines in sequence_lines_by_contig.items():
        sequences[contig] = "".join(sequence_lines)

    return sequences


def check_contigs_match_table(sequences, wanted_contigs, contig_lengths_in_table,
                              fasta_path, tblout_path):
    """Stop with a clear error if the nhmmer table and the genome FASTA do not fit.

    Takes in:  the contigs read by read_contigs(), the names of the contigs that
               have a 5.8S site, the contig lengths stated in the nhmmer table, and
               the two file names (used only in the error messages).
    Does:      checks that every contig with a hit is present in the FASTA, and that
               its length in the FASTA equals the length nhmmer wrote in the table.
               If the table came from a different genome file, the hit coordinates
               would point at the wrong bases, and we would cut windows from the
               wrong places without noticing. Better to stop here.
    Produces:  nothing if all is well; otherwise raises ValueError, which main()
               turns into a one-line error message and a non-zero exit code.
    Used by:   run(), just before the windows are cut.
    """
    # First check: every contig with a site is in the FASTA. All missing names are
    # reported at once, in alphabetical order.
    missing_contigs = []
    for contig in sorted(wanted_contigs):
        if contig not in sequences:
            missing_contigs.append(contig)

    if len(missing_contigs) > 0:
        raise ValueError(
            f"these contigs have 5.8S hits in {tblout_path} but are not in {fasta_path}: "
            f"{', '.join(missing_contigs)}. The table must be made from this same FASTA file."
        )

    # Second check: each of them has the length the table states. The first mismatch,
    # in alphabetical order of the contig names, stops the script.
    for contig in sorted(wanted_contigs):
        length_in_fasta = len(sequences[contig])
        length_in_table = contig_lengths_in_table[contig]
        if length_in_fasta != length_in_table:
            raise ValueError(
                f"contig {contig} is {length_in_table} bases long in {tblout_path} but "
                f"{length_in_fasta} bases long in {fasta_path}. The table must be made "
                "from this same FASTA file, otherwise the window coordinates are wrong."
            )


# ======================================================================
# STEP 4: CUT ONE WINDOW AROUND EVERY SITE
# ======================================================================

def site_sort_key(site):
    """Sorting helper: tell list.sort() how to order the sites of one contig.

    Takes in:  one site (a dictionary, see merge_hits_into_sites()).
    Produces:  the values to sort by: start first, then end, then strand. Sorting by
               start puts the sites from left to right along the contig; end and
               strand only break ties, so the order is the same on every run.
    Used by:   order_sites().
    """
    return (site["start"], site["end"], site["strand"])


def order_sites(sites, sequences):
    """Put the sites in a predictable order, so that the output never changes between runs.

    Takes in:  the sites from merge_hits_into_sites(), and the contigs from
               read_contigs() (only the order of their names is used here).
    Does:      merge_hits_into_sites() returns the sites one (contig, strand) group
               at a time, in the order in which each group has its first good hit in
               the nhmmer table. nhmmer lists its hits profile by profile, not along
               the genome, so that order says nothing about where the sites lie. We
               sort the sites by contig, in the order the contigs appear in the
               genome FASTA, and then from left to right along each contig, with the
               plus- and minus-strand sites mixed together.
    Produces:  the same sites as a new list in that order.
    Used by:   run(), so that the windows are written in this order.
    """
    ordered_sites = []

    # Looping over a dictionary gives its keys, here the contig names in FASTA order.
    for contig in sequences:
        sites_on_contig = []
        for site in sites:
            if site["contig"] == contig:
                sites_on_contig.append(site)

        sites_on_contig.sort(key=site_sort_key)
        ordered_sites.extend(sites_on_contig)

    return ordered_sites


def cut_window(site, contig_sequence, flank):
    """Cut the window around one 5.8S site and turn it to the standard orientation.

    Takes in:  one site from merge_hits_into_sites(), the full sequence of the contig
               it sits on (from read_contigs()), and --flank.
    Does:      centres the window on the middle of the site and takes `flank` bases
               on each side, so a full-size window is 2 * flank + 1 bases long. Near
               the start or the end of a contig the window is clipped to what exists.
               If the site is on the minus strand, the window is reverse-complemented,
               so that it reads SSU - ITS1 - 5.8S - ITS2 - LSU, which is the
               orientation ITSx expects.
    Produces:  a dictionary with the keys window_id, contig, start, end, strand,
               anchor_hits, best_evalue, window_length (the columns of the output
               table) and sequence (the window itself). The window_id is
               contig:start-end_strand, with start and end being the 1-based,
               inclusive, clipped positions of the window on the contig as stored.
    Used by:   run(), which hands the windows to write_windows_fasta() and
               write_windows_table().
    """
    contig_length = len(contig_sequence)

    # The middle base of the site. "//" divides and rounds down to a whole number.
    middle = (site["start"] + site["end"]) // 2

    # Window limits, counted from 1, never running off either end of the contig.
    start = max(1, middle - flank)
    end = min(contig_length, middle + flank)

    # Python counts positions from 0 and a slice stops just before its end value.
    # So the bases from position `start` to position `end` (counted from 1, both
    # included) are the slice [start - 1 : end].
    window_sequence = contig_sequence[start - 1:end]

    if site["strand"] == "-":
        window_sequence = reverse_complement(window_sequence)

    # The strand follows an underscore and not parentheses, because ITSx mishandles
    # parentheses in sequence headers.
    window_id = f"{site['contig']}:{start}-{end}_{site['strand']}"

    window = {
        "window_id": window_id,
        "contig": site["contig"],
        "start": start,
        "end": end,
        "strand": site["strand"],
        "anchor_hits": site["anchor_hits"],
        "best_evalue": site["best_evalue"],
        "window_length": len(window_sequence),
        "sequence": window_sequence,
    }
    return window


# ======================================================================
# STEP 5: WRITE THE OUTPUT FILES
# ======================================================================

def write_windows_fasta(path, windows):
    """Write the windows as a FASTA file (this is the input of ITSx).

    Takes in:  the list of windows from cut_window(), possibly empty.
    Does:      writes ">window_id" as the header and the sequence below it, wrapped
               at FASTA_LINE_WIDTH characters per line.
    Produces:  the file --out-fasta. With an empty list the file is empty (0 bytes).
    Used by:   run(). The file is read next by rule its_extraction, which runs ITSx
               on it, or, when the file is empty, writes empty ITS files instead.
    """
    with open(path, "w", encoding="utf-8") as fasta:
        for window in windows:
            fasta.write(f">{window['window_id']}\n")

            sequence = window["sequence"]
            for line_start in range(0, len(sequence), FASTA_LINE_WIDTH):
                line_end = line_start + FASTA_LINE_WIDTH
                fasta.write(sequence[line_start:line_end] + "\n")


def write_windows_table(path, windows):
    """Write one row per window, so every window can be traced back to its evidence.

    Takes in:  the list of windows from cut_window(), possibly empty.
    Does:      writes a header line and then one tab-separated row per window. The
               best E-value is written in scientific notation with three
               significant digits (for example 3.10e-12).
    Produces:  the file --out-table. With an empty list it holds only the header line.
    Used by:   run(). Nothing in the workflow reads the table back; it is there so
               that a reader can trace every window to its contig, strand and 5.8S
               evidence.
    """
    with open(path, "w", encoding="utf-8") as table:
        table.write(
            "window_id\tcontig\tstart\tend\tstrand\tanchor_hits\tbest_evalue\twindow_length\n"
        )
        for window in windows:
            # One value per column, in the order of the header line above. ".2e" gives
            # the E-value form 3.10e-12 (two digits after the point).
            fields = [
                window["window_id"],
                window["contig"],
                str(window["start"]),
                str(window["end"]),
                window["strand"],
                str(window["anchor_hits"]),
                f"{window['best_evalue']:.2e}",
                str(window["window_length"]),
            ]
            table.write("\t".join(fields) + "\n")


# ======================================================================
# PUTTING THE STEPS TOGETHER
# ======================================================================

def run(args):
    """Do the whole job: nhmmer table and genome FASTA in, windows FASTA and table out.

    Takes in:  args, the parsed command-line options (see parse_arguments()).
    Does:      runs steps 1 to 5 in order. All checks of the input happen before the
               first output file is written, so bad input never leaves half-written
               outputs behind.
    Produces:  the two output files and a one-line summary on stderr (rule
               its_windows adds it to its log, logs/its_windows_{sample}.log).
               Raises ValueError or OSError on bad input; main() turns that into a
               message.
    Used by:   main().
    """
    # Step 1: read the hit table and keep the good hits.
    hits_by_contig_and_strand, contig_lengths_in_table, n_weak_hits = read_hits(
        args.tblout, args.max_evalue
    )

    # No good hit at all means that this genome has no recognisable 5.8S. That is a
    # legitimate result for some genomes and not an error, so write empty outputs
    # (an empty FASTA and a header-only table) and finish normally. The genome FASTA
    # is not needed in this case, so it is not even opened.
    # (":g" in the messages writes the E-value limit in a compact form, e.g. 1e-05.)
    if len(hits_by_contig_and_strand) == 0:
        write_windows_fasta(args.out_fasta, [])
        write_windows_table(args.out_table, [])
        print(
            "its_windows.py: no 5.8S found (no nhmmer hit with E-value <= "
            f"{args.max_evalue:g}; {n_weak_hits} weaker hit(s) ignored). "
            "Wrote an empty FASTA and a header-only table.",
            file=sys.stderr,
        )
        return

    # Step 2: merge the hits that belong to the same 5.8S gene into sites.
    sites = merge_hits_into_sites(hits_by_contig_and_strand, args.merge_distance)

    # Step 3: read only the contigs that carry a site, and check they fit the table.
    # The same loop also counts the good hits, for the summary line at the end.
    contigs_with_sites = set()
    n_good_hits = 0
    for site in sites:
        contigs_with_sites.add(site["contig"])
        n_good_hits += site["anchor_hits"]

    sequences = read_contigs(args.fasta, contigs_with_sites)
    check_contigs_match_table(
        sequences, contigs_with_sites, contig_lengths_in_table, args.fasta, args.tblout
    )

    # Step 4: cut one window around every site, in a predictable order.
    ordered_sites = order_sites(sites, sequences)
    windows = []
    for site in ordered_sites:
        window = cut_window(site, sequences[site["contig"]], args.flank)
        windows.append(window)

    # Step 5: write the two output files.
    write_windows_fasta(args.out_fasta, windows)
    write_windows_table(args.out_table, windows)

    # One summary line: windows written, contigs they lie on, hits used and ignored.
    print(
        f"its_windows.py: wrote {len(windows)} window(s) on {len(sequences)} contig(s), "
        f"merged from {n_good_hits} nhmmer hit(s) with E-value <= {args.max_evalue:g} "
        f"({n_weak_hits} weaker hit(s) ignored).",
        file=sys.stderr,
    )


# ======================================================================
# COMMAND LINE
# ======================================================================

def parse_arguments(argv):
    """Read the options given on the command line.

    Takes in:  argv, the list of words after the script name (None means "use the real
               command line"; the tests pass their own list).
    Does:      defines the options and checks them, using argparse (Python's built-in
               reader of command-line options). Stops with a usage message if a
               required option is missing or --flank is smaller than 1.
    Produces:  an object whose attributes are the options: tblout, fasta, out_fasta,
               out_table, flank, max_evalue and merge_distance (argparse turns the
               dash of --out-fasta, --out-table, --max-evalue and --merge-distance
               into an underscore).
    Used by:   main().
    """
    parser = argparse.ArgumentParser(
        prog="its_windows.py",
        description=(
            "Cut a window around every 5.8S gene found by nhmmer, so that ITSx sees one "
            "rDNA unit at a time, always in the SSU - ITS1 - 5.8S - ITS2 - LSU orientation."
        ),
    )
    parser.add_argument(
        "--tblout", required=True,
        help="nhmmer table (nhmmer --tblout) of the ITSx fungal 5.8S profiles against the genome",
    )
    parser.add_argument(
        "--fasta", required=True,
        help="genome FASTA (contigs) that the nhmmer table was made from",
    )
    parser.add_argument(
        "--out-fasta", required=True,
        help="output: the windows, as FASTA with headers like contig:start-end_strand",
    )
    parser.add_argument(
        "--out-table", required=True,
        help="output: one row per window (window_id, contig, start, end, strand, "
             "anchor_hits, best_evalue, window_length), tab-separated",
    )
    parser.add_argument(
        "--flank", type=int, default=DEFAULT_FLANK,
        help="bases kept on each side of the middle of the 5.8S site (default: %(default)s)",
    )
    parser.add_argument(
        "--max-evalue", type=float, default=DEFAULT_MAX_EVALUE,
        help="ignore nhmmer hits with an E-value above this (default: %(default)s)",
    )
    parser.add_argument(
        "--merge-distance", type=int, default=DEFAULT_MERGE_DISTANCE,
        help="hits on the same contig and strand that start at most this many bases after "
             "the end of the site built so far belong to the same 5.8S site "
             "(default: %(default)s)",
    )

    args = parser.parse_args(argv)

    # A flank below 1 would cut empty or nonsensical windows. parser.error() prints the
    # usage and the message and stops with exit code 2, like a missing option does.
    if args.flank < 1:
        parser.error("--flank must be at least 1")

    return args


def main(argv=None):
    """Run the script from the command line.

    Takes in:  argv, the command-line words (None = the real command line).
    Does:      parses the options and calls run(). If run() finds a problem with the
               input, prints one readable line to stderr instead of a long Python
               traceback and exits with code 1, so Snakemake sees that the step failed.
    Produces:  exit code 0 on success (including the "no 5.8S found" case), 1 on a
               problem with the input files, 2 on a bad command line.
    Used by:   the last lines of this file when rule its_windows runs the script, and
               by the tests, which pass their own argv.
    """
    args = parse_arguments(argv)
    try:
        run(args)
    except (OSError, ValueError) as problem:
        # OSError: a file cannot be opened or written (for example it does not exist).
        # ValueError: the content of an input file is wrong (raised by parse_hit_line(),
        # read_contigs() and check_contigs_match_table(), or by Python when a file is
        # not UTF-8 text).
        print(f"its_windows.py: error: {problem}", file=sys.stderr)
        sys.exit(1)


# Run main() only when this file is run as a script, not when the tests import it.
if __name__ == "__main__":
    main()
