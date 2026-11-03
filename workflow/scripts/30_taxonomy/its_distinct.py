#!/usr/bin/env python3
"""Collapse the ITS records of all rDNA copies into distinct sequences, count the
copies, and keep implausible records away from the classifier.

BIOLOGY AND WHY THIS SCRIPT EXISTS
----------------------------------
The fungal ITS barcode (ITS1, 5.8S, ITS2) sits inside the ribosomal DNA (rDNA)
repeat, and the repeat is present many times in one genome: tens to hundreds of
copies in a row (a "tandem array"; 14 to 1,442, median 82, estimated in 91 fungal
genomes by Lofgren et al. 2019, Molecular Ecology 28: 721-730), usually identical
or nearly so. A short-read assembly folds the copies into one; a long-read assembly
can resolve several. The workflow cuts one window around every 5.8S gene
(its_windows.py) and runs ITSx on the windows, so ITSx returns ONE RECORD PER COPY:
an assembly with five copies gives five ITS1 records, five ITS2 records and five
"full" (ITS1 + 5.8S + ITS2) records. Two things follow.

  * Classifying every copy against UNITE (a curated reference database of fungal
    ITS sequences) would repeat the same work and print the same answer once per
    copy. Identical sequences are therefore collapsed into ONE record for the
    classifier.
  * The number of copies is itself a result (it describes the rDNA array as the
    assembly holds it), so it is counted and written to a table.

The second job of the script is a safety net. ITS1 and ITS2 of fungi are about
50 to 600 bases long. ITSx can occasionally stitch across a damaged locus and report
a stretch of genome far longer than any ITS as an "ITS2", and the classifier
(SINTAX) would give that stretch a meaningless taxonomy. Such a record
is a mis-extraction, not a barcode. Records outside the length limits are NOT passed
on to the classifier; they are listed in the table together with the reason, so
nothing disappears silently.

Both jobs are new in v2 (docs/methods/its-long-read-contigs.md describes the
design); v1 classified every record that ITSx returned.

WHERE THIS SCRIPT FITS IN THE WORKFLOW
--------------------------------------
    1. Rule its_windows cuts one window per rDNA copy (its_windows.py).
    2. Rule its_extraction runs ITSx on the windows and writes {sample}.ITS1.fasta,
       {sample}.ITS2.fasta and {sample}.full.fasta, one record per copy. The header
       of each record is the window it came from, contig:start-end_strand. For a
       genome with no recognisable 5.8S the three files are empty.
    3. THIS SCRIPT is run by rule its_classification on those three files. It writes
       the table straight to {sample}_ITS_copies.txt, and the two FASTA files of
       distinct sequences into a temporary folder.
    4. The same rule then runs `vsearch --sintax` on each of those FASTA files that
       is not empty, and writes {sample}_ITS1_taxonomy.txt and
       {sample}_ITS2_taxonomy.txt (empty when there was nothing to classify). The
       temporary folder is deleted when the rule ends; the sequence names link each
       SINTAX result to its row in the table.

WHAT THIS SCRIPT DOES
---------------------
    1. Reads the three ITSx FASTA files (any of them may be empty).
    2. For each region (ITS1, ITS2, full) groups the records that have the same
       sequence. Capital and small letters are treated alike and line wrapping is
       ignored, so the sequence is compared in upper case with all line breaks removed.
    3. Judges the length of every distinct sequence:
         ITS1 and ITS2:  --min-part-length to --max-part-length bases are accepted.
         full:           at most --max-full-length bases is accepted (no minimum).
       A sequence that is too short or too long is "rejected", the others are "ok".
    4. Writes one FASTA record per ok sequence of ITS1 and of ITS2, named
       <sample>_<region>_<n> (n = 1, 2, ... in order of first appearance, counting
       only ok sequences of that region). The classifier output is then traceable to
       the table through this name. The "full" region is only reported in the table.
    5. Writes the table: one row per distinct sequence, ok or rejected.

INPUTS
------
    --sample           the sample name, used to name the output records.
    --its1 --its2 --full
                       the three ITSx FASTA files. A file that exists but holds no
                       record is fine; a file that does not exist is an error.
    --min-part-length  shortest ITS1 or ITS2 accepted (default 50)
    --max-part-length  longest ITS1 or ITS2 accepted (default 600)
    --max-full-length  longest "full" region accepted (default 1500, the length
                       above which ITSx 1.1.3 itself warns "ITS region is
                       suspiciously long"; see docs/methods/its-long-read-contigs.md)
                       The workflow always passes all three limits: the two part
                       limits are fixed in 00_common.smk, and the full-length limit
                       comes from parameters.its.max_full_length in the config (at
                       least 400). The defaults here are the workflow's defaults, for
                       running the script by hand.

OUTPUTS
-------
    --out-its1, --out-its2
        FASTA of the distinct ok sequences of that region, upper case, one line per
        sequence. Empty (0 bytes) when no sequence is ok.
    --out-table
        Tab-separated. The first line is always
            region<TAB>sequence_id<TAB>length<TAB>copies<TAB>status<TAB>windows
        followed by one row per distinct sequence, ordered ITS1, then ITS2, then
        full, and within a region in order of first appearance:
            region       ITS1, ITS2 or full
            sequence_id  <sample>_<region>_<n> for an ok sequence, "-" for a rejected one
            length       bases in the sequence
            copies       how many records (rDNA copies) carry this exact sequence
            status       ok, rejected_too_short or rejected_too_long
            windows      the window names of those records, separated by commas
        If all three input files hold no record at all, the table is the header line
        plus one comment line that starts with # and says that no ITS region was
        extracted, and both FASTA files are empty.

The script exits with status 0 whenever the input could be read (also when it finds
nothing, or when it rejects records), because a genome without a clean ITS is a result
and not a failure of the workflow. A file that does not exist, or a file that is not
FASTA (a sequence line before any ">" header, or a ">" with no name after it), stops
it with one readable message and exit status 1. Impossible limits on the command line
stop it with exit status 2, before any file is read.
"""

# Only Python's standard library is used, so the Python 3.12 that the vsearch conda
# environment of rule its_classification holds is all this script needs (see the
# note at the top of 30_taxonomy.smk).
import argparse
import sys


# ======================================================================
# SETTINGS
# ======================================================================

# The three regions that ITSx reports, in the order they are read and written to the
# table. "full" is the whole ITS1 + 5.8S + ITS2 stretch.
REGIONS = ("ITS1", "ITS2", "full")

# Default length limits in bases (see INPUTS above). Rule its_classification passes
# its own values on the command line: ITS_MIN_PART_LENGTH and ITS_MAX_PART_LENGTH,
# fixed in 00_common.smk, and ITS_MAX_FULL_LENGTH, read there from
# parameters.its.max_full_length in the config.
DEFAULT_MIN_PART_LENGTH = 50
DEFAULT_MAX_PART_LENGTH = 600
DEFAULT_MAX_FULL_LENGTH = 1500

# The three verdicts of the length guard, exactly as they are written to the table.
STATUS_OK = "ok"
STATUS_TOO_SHORT = "rejected_too_short"
STATUS_TOO_LONG = "rejected_too_long"

# What the sequence_id column holds for a sequence that was rejected.
NO_SEQUENCE_ID = "-"

# The columns of the table ({sample}_ITS_copies.txt), in order.
TABLE_COLUMNS = ("region", "sequence_id", "length", "copies", "status", "windows")


# ======================================================================
# STEP 1: READ AN ITSX FASTA FILE
# ======================================================================

def read_fasta_records(fasta_path):
    """Read one ITSx FASTA file and return its records as a list.

    Takes in:  the path of a FASTA file written by ITSx (rule its_extraction). It may
               be empty. The header of every record is the name of the window the
               sequence was cut from, for example ptg000001l:10001-14001_-
    Does:      reads the file line by line. The window name is the first word of the
               header line (the text after ">" up to the first space). The sequence
               lines that follow are joined into one piece of text, so line wrapping
               does not matter, and turned into capital letters, so that the same
               sequence written in small and in capital letters is the same sequence.
               Blank lines are skipped.
    Produces:  a list of (window_name, sequence) pairs in the order of the file; an
               empty list for an empty file.
    Used by:   run(), once for each of the three input files.
    Stops with a clear error if a sequence line comes before the first header (then
    the file is not FASTA) or if a header has no name.
    """
    records = []
    window_name = None       # the name of the record being read; None before the first header
    sequence_pieces = []     # the sequence lines of that record

    with open(fasta_path, "r", encoding="utf-8") as fasta:
        for line_number, line in enumerate(fasta, start=1):
            # strip() removes the end of line (also the Windows one) and blanks at both ends.
            text = line.strip()
            if text == "":
                continue

            if text.startswith(">"):
                # A new record starts here, so the one before it is complete: join its
                # sequence lines and turn them into capital letters.
                if window_name is not None:
                    sequence = "".join(sequence_pieces).upper()
                    records.append((window_name, sequence))

                words = text[1:].split()
                if len(words) == 0:
                    raise ValueError(
                        f"{fasta_path}, line {line_number}: a FASTA header line (starting "
                        "with >) has no name after the > sign"
                    )
                window_name = words[0]
                sequence_pieces = []
            else:
                # A sequence line. Before the first header it means the file is not FASTA.
                if window_name is None:
                    raise ValueError(
                        f"{fasta_path}, line {line_number}: found sequence before any FASTA "
                        f"header line (starting with >), so this is not a FASTA file: {text[:40]!r}"
                    )
                sequence_pieces.append(text)

    # The last record has no ">" after it to complete it.
    if window_name is not None:
        sequence = "".join(sequence_pieces).upper()
        records.append((window_name, sequence))

    return records


# ======================================================================
# STEP 2: COLLAPSE IDENTICAL SEQUENCES
# ======================================================================

def group_identical_sequences(records):
    """Put the records that carry the same sequence into one group.

    Takes in:  the list of (window_name, sequence) pairs of ONE region, from
               read_fasta_records().
    Does:      files each window name under its sequence. The sequences are already in
               capital letters without line breaks, so two copies of the rDNA unit
               that read the same are put together whatever the case or wrapping of
               the input file. Groups are kept in the order in which their sequence
               appears for the first time. (A Python dictionary remembers the order in
               which keys were added, which is what keeps the order here.)
    Produces:  a list of groups. Each group is a dictionary with the keys sequence and
               windows (the list of window names, in file order). The number of
               copies is the length of that list.
    Used by:   describe_region().
    """
    # A lookup table from each distinct sequence to the windows that carry it.
    window_names_by_sequence = {}
    for window_name, sequence in records:
        if sequence not in window_names_by_sequence:
            window_names_by_sequence[sequence] = []
        window_names_by_sequence[sequence].append(window_name)

    # The same content as a list of groups, in order of first appearance.
    groups = []
    for sequence, window_names in window_names_by_sequence.items():
        groups.append({"sequence": sequence, "windows": window_names})
    return groups


# ======================================================================
# STEP 3: THE LENGTH GUARD
# ======================================================================

def judge_length(region, length, min_part_length, max_part_length, max_full_length):
    """Decide whether a sequence of this length can be a real ITS region.

    Takes in:  the region (ITS1, ITS2 or full), the length of one distinct sequence in
               bases, and the three limits.
    Does:      ITS1 and ITS2 must be at least min_part_length and at most
               max_part_length bases long (both limits are accepted). The "full"
               region must be at most max_full_length bases long; it has no minimum,
               because its shortness is not what goes wrong when ITSx stitches across
               a damaged locus.
    Produces:  one of the three words ok, rejected_too_short, rejected_too_long.
    Used by:   describe_region().
    """
    # The full region: a maximum only.
    if region == "full":
        if length > max_full_length:
            return STATUS_TOO_LONG
        return STATUS_OK

    # ITS1 or ITS2: a minimum and a maximum.
    if length < min_part_length:
        return STATUS_TOO_SHORT
    if length > max_part_length:
        return STATUS_TOO_LONG
    return STATUS_OK


def describe_region(sample, region, records, min_part_length, max_part_length, max_full_length):
    """Turn the records of one region into one table row per distinct sequence.

    Takes in:  the sample name, the region, the (window_name, sequence) pairs of that
               region, and the three length limits.
    Does:      groups identical sequences (group_identical_sequences), judges the
               length of each group (judge_length) and gives the groups that pass a
               name <sample>_<region>_<n>. The counter n only moves on for sequences
               that pass, so the names of a region are always 1, 2, 3, ... with no
               gaps, whatever was rejected in between. A rejected sequence gets the
               name "-".
    Produces:  a list of rows in order of first appearance. Each row is a dictionary
               with the keys region, sequence_id, length, copies, status, windows (a
               list of window names) and sequence.
    Used by:   run().
    """
    rows = []
    n_ok = 0   # how many sequences of this region have passed the guard so far

    for group in group_identical_sequences(records):
        sequence = group["sequence"]
        length = len(sequence)
        status = judge_length(region, length, min_part_length, max_part_length, max_full_length)

        if status == STATUS_OK:
            n_ok += 1
            sequence_id = f"{sample}_{region}_{n_ok}"
        else:
            sequence_id = NO_SEQUENCE_ID

        rows.append({
            "region": region,
            "sequence_id": sequence_id,
            "length": length,
            "copies": len(group["windows"]),
            "status": status,
            "windows": group["windows"],
            "sequence": sequence,
        })

    return rows


# ======================================================================
# STEP 4: WRITE THE OUTPUT FILES
# ======================================================================

def write_table(table_path, sample, rows_by_region, no_record_at_all):
    """Write the table of all distinct sequences, ok and rejected.

    Takes in:  the path for --out-table, the sample name, the rows of every region
               (a dictionary from region to the list made by describe_region()), and
               a yes/no: did all three input files hold no record?
    Does:      writes the header line, then the rows ordered ITS1, ITS2, full. The
               window names of a row are joined with commas. If there was no record
               at all, one comment line (starting with #) follows the header and says
               so, because an empty table would otherwise look like a failure.
    Produces:  the file --out-table. Rule its_classification passes
               {sample}_ITS_copies.txt as that path, so this is the final table,
               read by people.
    Used by:   run().
    """
    with open(table_path, "w", encoding="utf-8") as table:
        table.write("\t".join(TABLE_COLUMNS) + "\n")

        if no_record_at_all:
            table.write(
                f"# Sample {sample}: no ITS region was extracted (no recognisable 5.8S in the "
                "genome, or ITSx found nothing).\n"
            )

        # One row per distinct sequence: the ITS1 rows first, then ITS2, then full. The
        # fields follow TABLE_COLUMNS; the numbers are turned into text so that they
        # can be joined with tabs.
        for region in REGIONS:
            for row in rows_by_region[region]:
                fields = (
                    row["region"],
                    row["sequence_id"],
                    str(row["length"]),
                    str(row["copies"]),
                    row["status"],
                    ",".join(row["windows"]),
                )
                table.write("\t".join(fields) + "\n")


def write_distinct_fasta(fasta_path, rows):
    """Write the ok sequences of one region as FASTA, for the classifier.

    Takes in:  the path for --out-its1 or --out-its2, and the rows of that region from
               describe_region().
    Does:      writes ">sequence_id" and the sequence on one line for every row whose
               status is ok. Rejected rows are not written. The sequence is in capital
               letters (see read_fasta_records()); ITS regions are short enough that
               wrapping the line would gain nothing.
    Produces:  the file. With no ok row the file is empty (0 bytes), and rule
               its_classification then skips vsearch for that region and writes an
               empty taxonomy file instead.
    Used by:   run(). The file is read by vsearch --sintax in rule its_classification.
    """
    with open(fasta_path, "w", encoding="utf-8") as fasta:
        for row in rows:
            if row["status"] != STATUS_OK:
                continue
            fasta.write(f">{row['sequence_id']}\n")
            fasta.write(row["sequence"] + "\n")


# ======================================================================
# PUTTING THE STEPS TOGETHER
# ======================================================================

def summarise_region(region, records, rows):
    """Return a short piece of text about one region for the one-line report on stderr.

    Takes in:  the region, its records, and its rows from describe_region().
    Produces:  for example "ITS2 5 record(s), 1 distinct, 1 ok, 0 rejected".
    Used by:   run().
    """
    n_ok = 0
    for row in rows:
        if row["status"] == STATUS_OK:
            n_ok += 1
    n_rejected = len(rows) - n_ok
    return (
        f"{region} {len(records)} record(s), {len(rows)} distinct, "
        f"{n_ok} ok, {n_rejected} rejected"
    )


def run(args):
    """Do the whole job: three ITSx FASTA files in, two FASTA files and a table out.

    Takes in:  args, the parsed command-line options (see parse_arguments()).
    Does:      runs steps 1 to 4 in order. All the input is read and checked before
               the first output file is opened, so bad input never leaves
               half-written outputs behind.
    Produces:  the three output files and a one-line report on stderr (rule
               its_classification writes it to its log,
               logs/its_classification_{sample}.log). Raises OSError or ValueError on
               bad input; main() turns that into a message.
    Used by:   main().
    """
    # Step 1: read the three ITSx files, in the order of REGIONS. A missing or broken
    # file stops the script here, before any output is written.
    input_paths = {"ITS1": args.its1, "ITS2": args.its2, "full": args.full}
    records_by_region = {}
    for region in REGIONS:
        records_by_region[region] = read_fasta_records(input_paths[region])

    # Steps 2 and 3: collapse identical sequences and judge their length. The records
    # of all three regions are also counted, to recognise the "nothing extracted" case.
    rows_by_region = {}
    n_records = 0
    for region in REGIONS:
        records = records_by_region[region]
        n_records += len(records)
        rows_by_region[region] = describe_region(
            args.sample, region, records,
            args.min_part_length, args.max_part_length, args.max_full_length,
        )

    # Step 4: write the outputs.
    no_record_at_all = (n_records == 0)
    write_table(args.out_table, args.sample, rows_by_region, no_record_at_all)
    write_distinct_fasta(args.out_its1, rows_by_region["ITS1"])
    write_distinct_fasta(args.out_its2, rows_by_region["ITS2"])

    # The report on stderr: one line saying that nothing was extracted, or else one
    # line with the counts of every region.
    if no_record_at_all:
        print(
            f"its_distinct.py: no ITS region was extracted for sample {args.sample} (no "
            "recognisable 5.8S in the genome, or ITSx found nothing). Wrote empty FASTA "
            "files and a table with only the header and a comment line.",
            file=sys.stderr,
        )
        return

    # For each region: the records read, the distinct sequences among them, and how
    # many of those passed or failed the length guard.
    region_summaries = []
    for region in REGIONS:
        region_summary = summarise_region(
            region, records_by_region[region], rows_by_region[region]
        )
        region_summaries.append(region_summary)
    report = "; ".join(region_summaries)
    print(f"its_distinct.py: sample {args.sample}: {report}.", file=sys.stderr)


# ======================================================================
# COMMAND LINE
# ======================================================================

def parse_arguments(argv):
    """Read the options given on the command line.

    Takes in:  argv, the list of words after the script name (None means "use the real
               command line"; the tests pass their own list).
    Does:      defines the options and checks them, using argparse (Python's built-in
               reader of command-line options). Stops with a usage message if a
               required option is missing, a limit is smaller than 1, or the shortest
               accepted ITS1/ITS2 is longer than the longest.
    Produces:  an object whose attributes are the options: sample, its1, its2, full,
               min_part_length, max_part_length, max_full_length, out_its1, out_its2
               and out_table (argparse turns the dashes into underscores).
    Used by:   main().
    """
    parser = argparse.ArgumentParser(
        prog="its_distinct.py",
        description=(
            "Collapse the identical ITS1, ITS2 and full-length records that ITSx returns for "
            "the copies of the rDNA array, count the copies, and keep records of implausible "
            "length away from the classifier."
        ),
    )
    parser.add_argument(
        "--sample", required=True,
        help="sample name; the distinct sequences are named <sample>_<region>_<n>",
    )
    parser.add_argument("--its1", required=True, help="ITSx ITS1 FASTA (may be empty)")
    parser.add_argument("--its2", required=True, help="ITSx ITS2 FASTA (may be empty)")
    parser.add_argument("--full", required=True, help="ITSx full-length (ITS1-5.8S-ITS2) FASTA (may be empty)")
    parser.add_argument(
        "--min-part-length", type=int, default=DEFAULT_MIN_PART_LENGTH,
        help="shortest ITS1 or ITS2 accepted, in bases (default: %(default)s)",
    )
    parser.add_argument(
        "--max-part-length", type=int, default=DEFAULT_MAX_PART_LENGTH,
        help="longest ITS1 or ITS2 accepted, in bases (default: %(default)s)",
    )
    parser.add_argument(
        "--max-full-length", type=int, default=DEFAULT_MAX_FULL_LENGTH,
        help="longest full-length region accepted, in bases; there is no minimum (default: %(default)s)",
    )
    parser.add_argument(
        "--out-its1", required=True,
        help="output: FASTA of the distinct ITS1 sequences that passed the length guard",
    )
    parser.add_argument(
        "--out-its2", required=True,
        help="output: FASTA of the distinct ITS2 sequences that passed the length guard",
    )
    parser.add_argument(
        "--out-table", required=True,
        help="output: one row per distinct sequence (region, sequence_id, length, copies, "
             "status, windows), tab-separated",
    )

    args = parser.parse_args(argv)

    # Refuse limits that make no sense. A minimum below 1 would accept empty ITS1 and
    # ITS2 records; a maximum below the minimum would reject every ITS1 and ITS2; a
    # full-length limit below 1 would reject every real full-length record.
    # parser.error() prints the usage and the message and stops with exit code 2, like
    # a missing option does.
    if args.min_part_length < 1:
        parser.error("--min-part-length must be at least 1")
    if args.max_part_length < args.min_part_length:
        parser.error("--max-part-length must not be smaller than --min-part-length")
    if args.max_full_length < 1:
        parser.error("--max-full-length must be at least 1")

    return args


def main(argv=None):
    """Run the script from the command line.

    Takes in:  argv, the command-line words (None = the real command line).
    Does:      parses the options and calls run(). If run() finds a problem with the
               input, prints one readable line to stderr instead of a long Python
               traceback and exits with code 1, so Snakemake sees that the step failed.
    Produces:  exit code 0 on success (including the "no ITS region" case and a run
               that rejects records), 1 on a problem with the input files, 2 on a bad
               command line.
    Used by:   the last lines of this file when rule its_classification runs the
               script, and by the tests, which pass their own argv.
    """
    args = parse_arguments(argv)
    try:
        run(args)
    except (OSError, ValueError) as problem:
        # OSError: a file cannot be opened or written (for example it does not exist).
        # ValueError: an input file is not FASTA (raised by read_fasta_records()) or
        # is not UTF-8 text (raised by Python while reading it).
        print(f"its_distinct.py: error: {problem}", file=sys.stderr)
        sys.exit(1)


# Run main() only when this file is run as a script, not when the tests import it.
if __name__ == "__main__":
    main()
