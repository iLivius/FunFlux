#!/usr/bin/env python3
"""Keep the UNITE reference records whose taxonomy names a genus, for SINTAX.

BIOLOGY AND WHY THIS SCRIPT EXISTS
----------------------------------
The workflow classifies the ITS1 and ITS2 of each genome with SINTAX (the classifier
that vsearch offers) against UNITE, the curated reference database of fungal ITS
sequences, in its "general FASTA release for SINTAX". Every record of that release
carries its taxonomy in the header, rank by rank:

    >ACC|SH0000001.10FU;tax=d:Fungi,p:Ascomycota,...,g:Botrytis,s:Botrytis_cinerea;

Many records are only identified down to some rank above the genus. Their taxonomy
then simply stops there: for example "tax=d:Fungi;" alone (kingdom only), or a string
that ends with the family. In the release the workflow downloads by default, about a
third of the entries name no genus.

SINTAX works out its confidence values by classifying the query many times, each time
from a random sample of its k-mers (short words of the sequence), and reporting how
often the same name came out. When a kingdom-only record is identical to the query,
SINTAX follows that record as often as the identified ones, and the result stops at
the kingdom even for a genome whose ITS matches well-identified species perfectly.
Classifying against the records that carry a genus gives the genus and species with
the usual confidence values. The records without a genus add nothing a reader can
use, so the workflow leaves them out before classifying.

WHERE THIS SCRIPT FITS IN THE WORKFLOW
--------------------------------------
    1. Rule unite_db (download) or unite_db_local (link to a copy the user holds)
       makes 00.resources/unite/unite_its_sintax.fasta, the whole release, once per
       run (workflow/rules/shared/30_taxonomy.smk).
    2. THIS SCRIPT is run once per run by rule unite_identified on that file. It
       writes 00.resources/unite/unite_its_sintax.identified.fasta and prints the two
       counts, which the rule saves as 00.resources/unite/identified_filter.txt.
    3. Rule its_classification hands the identified file to vsearch --sintax (--db)
       for every sample. The configuration key parameters.its.unite_identified_only
       (default true) switches this on; with false the rule above is not planned and
       the classification reads the whole release as before.

WHAT THIS SCRIPT DOES
---------------------
    Reads the FASTA record by record and copies a record to the output when its
    header has a "tax=" field whose taxonomy holds a genus rank: a field "g:" followed
    by a name, such as g:Botrytis. Everything else is dropped: a record whose
    taxonomy stops above the genus, a record with "g:" and nothing after it, and a
    header without a "tax=" field at all. Kept records are copied as they are, header
    and sequence lines alike, so a sequence spread over several lines stays intact.
    Sequences are never read into memory as a whole: the file is handled line by line.

INPUTS
------
    --fasta      the UNITE SINTAX FASTA (uncompressed), from rule unite_db or
                 unite_db_local.

OUTPUTS
-------
    --out-fasta  the kept records, in the order of the input. Read by vsearch --sintax
                 in rule its_classification.
    standard output: two lines, the counts, exactly as rule unite_identified saves
                 them in identified_filter.txt:
                     kept_with_genus: <n>
                     dropped_without_genus: <n>
    standard error: one summary line (the rule writes it to its log,
                 logs/unite_identified.log), which also says how many of the dropped
                 headers had no "tax=" field.

Exit status 0 when the file was read and at least one record was kept. A file that
cannot be opened, a file that is not FASTA (a sequence line before any ">" header), or a
file in which no record carries a genus (then it is not a UNITE SINTAX release, or an
empty download) stops the script with one readable message and exit status 1, and no
output FASTA is left behind.
"""

# Only Python's standard library is used, so the Python 3.12 that the vsearch conda
# environment holds is all this script needs (see the note at the top of 30_taxonomy.smk).
import argparse
import os
import sys


# ======================================================================
# STEP 1: DECIDE WHETHER ONE HEADER NAMES A GENUS
# ======================================================================

def header_names_a_genus(header_line):
    """Tell whether a UNITE SINTAX header carries a genus rank.

    Takes in:  one header line of the UNITE FASTA, starting with ">", for example
               >ACC|SH0000001.10FU;tax=d:Fungi,p:Ascomycota,g:Botrytis,s:Botrytis_cinerea;
    Does:      finds the "tax=" field. The taxonomy is the text after "tax=" up to the
               first ";" (the release ends the field with one; without it the rest of
               the line is taken). The taxonomy is split at the commas into its ranks,
               and the header names a genus when one rank starts with "g:" and has a
               name after the colon. Only the taxonomy is searched, so a "g:" anywhere
               else in the header (in the accession, say) does not count.
    Produces:  (names_a_genus, has_tax_field), two yes/no values. The second one only
               serves the count of headers without a "tax=" field in the summary line.
    Used by:   filter_records().
    """
    position = header_line.find("tax=")
    if position == -1:
        return False, False

    taxonomy = header_line[position + len("tax="):]
    end = taxonomy.find(";")
    if end != -1:
        taxonomy = taxonomy[:end]

    for rank in taxonomy.split(","):
        rank = rank.strip()
        if rank.startswith("g:") and len(rank) > len("g:"):
            return True, True
    return False, True


# ======================================================================
# STEP 2: COPY THE KEPT RECORDS, LINE BY LINE
# ======================================================================

def filter_records(fasta_path, out_fasta_path):
    """Copy the records that name a genus from the UNITE FASTA to the output FASTA.

    Takes in:  the path of the UNITE SINTAX FASTA and the path to write.
    Does:      reads the input line by line. A line that starts with ">" begins a
               record and decides, with header_names_a_genus(), whether that record is
               copied; every following line up to the next header belongs to the same
               record and is copied or skipped with it. Blank lines are skipped. The
               lines are written exactly as read, so the output is the input minus the
               dropped records.
    Produces:  the output FASTA, and three counts: records kept, records dropped,
               and, among the dropped, headers without a "tax=" field.
    Used by:   run().
    Stops with a clear error (ValueError) if a sequence line comes before the first
    header, because then the file is not FASTA.
    """
    n_kept = 0
    n_dropped = 0
    n_without_tax = 0
    # Whether the record being read is copied. None before the first header.
    keep_current_record = None

    with open(fasta_path, "r", encoding="utf-8") as fasta, \
            open(out_fasta_path, "w", encoding="utf-8") as out_fasta:
        for line_number, line in enumerate(fasta, start=1):
            # Keep the line's own end-of-line out of the test, but write the line as it is.
            text = line.strip()
            if text == "":
                continue

            if text.startswith(">"):
                names_a_genus, has_tax_field = header_names_a_genus(text)
                keep_current_record = names_a_genus
                if names_a_genus:
                    n_kept += 1
                else:
                    n_dropped += 1
                    if not has_tax_field:
                        n_without_tax += 1
            elif keep_current_record is None:
                # A sequence line before the first header: this is not a FASTA file.
                raise ValueError(
                    f"{fasta_path}, line {line_number}: found sequence before any FASTA "
                    f"header line (starting with >), so this is not a FASTA file: {text[:40]!r}"
                )

            if keep_current_record:
                out_fasta.write(line)

    return n_kept, n_dropped, n_without_tax


# ======================================================================
# PUTTING THE STEPS TOGETHER
# ======================================================================

def run(args):
    """Do the whole job: the UNITE FASTA in, the identified FASTA and the counts out.

    Takes in:  args, the parsed command-line options (see parse_arguments()).
    Does:      calls filter_records(), then checks that at least one record was kept:
               none at all means the input is not a UNITE SINTAX release (or an empty
               file), and classifying against an empty reference would silently give
               empty results. Prints the two count lines on standard output and the
               summary line on standard error.
    Produces:  the output FASTA and the printed lines. Raises OSError or ValueError on
               bad input; main() turns that into a message.
    Used by:   main().
    """
    n_kept, n_dropped, n_without_tax = filter_records(args.fasta, args.out_fasta)

    if n_kept == 0:
        raise ValueError(
            f"{args.fasta}: no record carries a genus rank (g:<name> in its tax= field). "
            "Is this the UNITE general FASTA release for SINTAX, and is the file complete?"
        )

    # The two count lines: rule unite_identified saves standard output as
    # identified_filter.txt, so nothing else may be printed there.
    print(f"kept_with_genus: {n_kept}")
    print(f"dropped_without_genus: {n_dropped}")

    n_records = n_kept + n_dropped
    print(
        f"unite_identified.py: {n_records} records read, {n_kept} kept (taxonomy names a "
        f"genus), {n_dropped} dropped (no genus; {n_without_tax} of them without a tax= "
        "field).",
        file=sys.stderr,
    )


# ======================================================================
# COMMAND LINE
# ======================================================================

def parse_arguments(argv):
    """Read the options given on the command line.

    Takes in:  argv, the list of words after the script name (None means "use the real
               command line"; the tests pass their own list).
    Does:      defines the two options with argparse (Python's built-in reader of
               command-line options) and stops with a usage message if one is missing.
    Produces:  an object with the attributes fasta and out_fasta.
    Used by:   main().
    """
    parser = argparse.ArgumentParser(
        prog="unite_identified.py",
        description=(
            "Keep the records of a UNITE SINTAX FASTA whose taxonomy names a genus "
            "(g:<name>), and print how many were kept and dropped."
        ),
    )
    parser.add_argument("--fasta", required=True, help="the UNITE SINTAX FASTA, uncompressed")
    parser.add_argument(
        "--out-fasta", required=True,
        help="output: the records whose taxonomy names a genus, in the input order",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Run the script from the command line.

    Takes in:  argv, the command-line words (None = the real command line).
    Does:      parses the options and calls run(). If run() finds a problem with the
               input, prints one readable line to stderr instead of a Python traceback,
               removes the output FASTA it may have started (so no half-written
               reference is left for vsearch to read), and exits with code 1, so
               Snakemake sees that the step failed.
    Produces:  exit code 0 on success, 1 on a problem with the input, 2 on a bad
               command line (argparse's own code).
    Used by:   the last lines of this file when rule unite_identified runs the script,
               and by the tests, which pass their own argv.
    """
    args = parse_arguments(argv)
    try:
        run(args)
    except (OSError, ValueError) as problem:
        # OSError: a file cannot be opened or written (for example it does not exist).
        # ValueError: the input is not FASTA or holds no record with a genus (raised in
        # filter_records() and run()), or is not UTF-8 text (raised by Python while
        # reading it).
        print(f"unite_identified.py: error: {problem}", file=sys.stderr)
        if os.path.exists(args.out_fasta):
            os.remove(args.out_fasta)
        sys.exit(1)


# Run main() only when this file is run as a script, not when the tests import it.
if __name__ == "__main__":
    main()
