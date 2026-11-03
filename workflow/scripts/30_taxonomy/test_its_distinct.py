"""Tests for its_distinct.py, the script that collapses the ITS records of the rDNA copies.

How to run (from the repository root):

    python -m pytest workflow/scripts/30_taxonomy -q

Only the Python standard library and pytest are needed. No real genome is needed
either: every input is made up while the tests run.

THE MADE-UP DATA
----------------
  * Records are (header, sequence) pairs. The header is a window name of the kind ITSx
    keeps from its input, contig:start-end_strand, for example ptg1:1-4001_-.
  * Random sequences of an exact length are made by dna(length, seed). A fixed seed
    gives the same sequence on every run, and two different seeds give two different
    sequences (the chance that two random 100-base sequences agree is negligible).
  * Three made-up sequences of realistic length are written into this file: an ITS1
    (165 bases), an ITS2 (170 bases) and the full-length region that holds them (490
    bases: ITS1, a 154-base stretch standing in for the 5.8S gene, ITS2 and one more
    base). Five identical copies of them are what a long-read assembly with five
    rDNA units gives, which is the case the script was written for.
  * The limits are the workflow's: ITS1 and ITS2 from 50 to 600 bases, full region at
    most 1500 bases.
"""

import os
import random
import subprocess
import sys
from types import SimpleNamespace

import pytest

# SimpleNamespace (used by make_inputs below) is a plain container whose fields are read
# with a dot, like made.its1.

# The script sits in the same folder as this test file. Putting that folder on the
# list of places where Python looks for modules makes "import its_distinct" work
# wherever pytest is started from.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TEST_DIR)

import its_distinct  # noqa: E402  (this import must come after the line above)

SCRIPT_PATH = os.path.join(TEST_DIR, "its_distinct.py")

SAMPLE = "S1"

# The first line of the table, exactly as the script must write it.
TABLE_HEADER_LINE = "region\tsequence_id\tlength\tcopies\tstatus\twindows"

# The made-up example sequences (see the top of this file). They were drawn at
# random once and pasted here, so they are the same on every run.
EXAMPLE_ITS1 = (
    "GTTAAGGTTACGTCCTGCTTAGCCAGCGCGCACAACGCTCATACACGTTAGCGCCGAGAC"
    "CGACTGTAAGACACTCATACGGCGAGGATCGTATTGAGGTAGCGAGACAACCCCTGCGAG"
    "TTTGGCGAGGGAGGCTGAAACGGGCTATAGGCCCATATCGTTAGA"
)
EXAMPLE_ITS2 = (
    "ATGAGGGTAGGCCCTGTTAACACGCGTCTTCATTCCGGTCCGGTCAGGGGCCCGCCAGGA"
    "TTTTACGGTCCGAACGGTATTAACAGCAGTCAGTTGTGCGTCTTGGGTAGTTATAAATTC"
    "AGGGCCATATGACCGCCAGAAAAAAACATTGCGCTCGGCGAGAGAGATTC"
)
EXAMPLE_FULL = (
    "GTTAAGGTTACGTCCTGCTTAGCCAGCGCGCACAACGCTCATACACGTTAGCGCCGAGAC"
    "CGACTGTAAGACACTCATACGGCGAGGATCGTATTGAGGTAGCGAGACAACCCCTGCGAG"
    "TTTGGCGAGGGAGGCTGAAACGGGCTATAGGCCCATATCGTTAGATTAACTCCGAGAGAT"
    "CGATTGAGCGCAGCCTTATACATTCGATAGTCGTAGTCTATGATATGTGGCTTTGTACTG"
    "CTGACCACGCTCGCGTTGCCAAAGGACTCGTACTCACTTCAACAAACGCAGACCCCATGC"
    "TAGAGCACTATGTTATGGAATGAGGGTAGGCCCTGTTAACACGCGTCTTCATTCCGGTCC"
    "GGTCAGGGGCCCGCCAGGATTTTACGGTCCGAACGGTATTAACAGCAGTCAGTTGTGCGT"
    "CTTGGGTAGTTATAAATTCAGGGCCATATGACCGCCAGAAAAAAACATTGCGCTCGGCGA"
    "GAGAGATTCG"
)


# ======================================================================
# HELPERS: MAKING INPUTS AND READING OUTPUTS
# ======================================================================

def dna(length, seed):
    """Return random A/C/G/T text of exactly `length` bases (the same text for the same seed)."""
    rng = random.Random(seed)
    return "".join(rng.choices("ACGT", k=length))


def window_names(contig, n):
    """Return n made-up window names on one contig, the way its_windows.py names them."""
    return [f"{contig}:{1 + 8000 * i}-{4001 + 8000 * i}_-" for i in range(n)]


def records_for(names, sequence):
    """Return one (window_name, sequence) record for every window name: identical copies."""
    return [(name, sequence) for name in names]


def write_fasta_file(path, records, line_width=70):
    """Write a FASTA file from a list of (header, sequence) pairs.

    The sequences are wrapped at 70 columns, so that the tests also show that
    multi-line records are read correctly. An empty list gives an empty file.
    """
    with open(path, "w") as fasta:
        for header, sequence in records:
            fasta.write(f">{header}\n")
            for line_start in range(0, len(sequence), line_width):
                fasta.write(sequence[line_start:line_start + line_width] + "\n")


def make_inputs(folder, its1=(), its2=(), full=()):
    """Write the three ITSx files into `folder` and return their paths and the output paths.

    Each of its1, its2 and full is either a list of (header, sequence) pairs (written as
    FASTA, wrapped at 70 columns) or a text that is written exactly as it is (for the
    tests that need to control the line breaks, blank lines or errors).
    """
    made = SimpleNamespace(
        its1=folder / "S1.ITS1.fasta",
        its2=folder / "S1.ITS2.fasta",
        full=folder / "S1.full.fasta",
        out_its1=folder / "distinct.ITS1.fasta",
        out_its2=folder / "distinct.ITS2.fasta",
        out_table=folder / "S1_ITS_copies.txt",
    )
    for path, content in ((made.its1, its1), (made.its2, its2), (made.full, full)):
        if isinstance(content, str):
            path.write_text(content)
        else:
            write_fasta_file(path, list(content))
    return made


def command_line(made, *extra_options):
    """Return the command-line words for `made`, followed by any extra options."""
    words = [
        "--sample", SAMPLE,
        "--its1", str(made.its1),
        "--its2", str(made.its2),
        "--full", str(made.full),
        "--out-its1", str(made.out_its1),
        "--out-its2", str(made.out_its2),
        "--out-table", str(made.out_table),
    ]
    words.extend(extra_options)
    return words


def run_script(made, *extra_options):
    """Run its_distinct.main() on the files in `made`, adding any extra options."""
    its_distinct.main(command_line(made, *extra_options))


def read_table(path):
    """Read the table back. Returns (all_lines, rows).

    rows holds one dictionary per data row (a line that does not start with #), from
    the column name to the text found in that column.
    """
    with open(path) as table:
        all_lines = table.read().splitlines()

    column_names = all_lines[0].split("\t")
    rows = []
    for line in all_lines[1:]:
        if line.startswith("#"):
            continue
        values = line.split("\t")
        row = {}
        for position, column_name in enumerate(column_names):
            row[column_name] = values[position]
        rows.append(row)
    return all_lines, rows


def read_fasta_back(path):
    """Read a FASTA output back. Returns a list of (header, sequence) pairs.

    The script writes each sequence on one line, so a record is two lines; the test
    for the format of the FASTA checks that separately.
    """
    with open(path) as fasta:
        lines = fasta.read().splitlines()

    records = []
    for line in lines:
        if line.startswith(">"):
            records.append([line[1:], ""])
        else:
            records[-1][1] += line
    return [(header, sequence) for header, sequence in records]


# ======================================================================
# IDENTICAL COPIES ARE COLLAPSED AND COUNTED
# ======================================================================

def test_identical_copies_collapse_into_one_row_with_the_copy_count(tmp_path):
    """Five identical ITS2 records give ONE row with copies 5 and the five windows listed."""
    names = window_names("ptg1", 5)
    made = make_inputs(tmp_path, its2=records_for(names, dna(170, 1)))

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert len(rows) == 1
    row = rows[0]
    assert row["region"] == "ITS2"
    assert row["sequence_id"] == "S1_ITS2_1"
    assert row["length"] == "170"
    assert row["copies"] == "5"
    assert row["status"] == "ok"
    # The windows are listed in the order of the file, separated by commas.
    assert row["windows"] == ",".join(names)


def test_distinct_fasta_holds_one_record_per_distinct_sequence(tmp_path):
    """The classifier is given one record, not five, and it carries the traceable name."""
    sequence = dna(170, 1)
    made = make_inputs(tmp_path, its2=records_for(window_names("ptg1", 5), sequence))

    run_script(made)

    assert read_fasta_back(made.out_its2) == [("S1_ITS2_1", sequence)]
    # ITS1 had no record, so its file exists but is empty.
    assert made.out_its1.read_text() == ""


def test_two_different_sequences_give_two_ids_in_order_of_first_appearance(tmp_path):
    """A B A B A: A appeared first, so A is _1 (3 copies) and B is _2 (2 copies)."""
    seq_a = dna(150, 10)
    seq_b = dna(160, 11)
    records = [
        ("w1", seq_a), ("w2", seq_b), ("w3", seq_a), ("w4", seq_b), ("w5", seq_a),
    ]
    made = make_inputs(tmp_path, its1=records)

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert [row["sequence_id"] for row in rows] == ["S1_ITS1_1", "S1_ITS1_2"]
    assert [row["length"] for row in rows] == ["150", "160"]
    assert [row["copies"] for row in rows] == ["3", "2"]
    assert [row["windows"] for row in rows] == ["w1,w3,w5", "w2,w4"]
    assert read_fasta_back(made.out_its1) == [("S1_ITS1_1", seq_a), ("S1_ITS1_2", seq_b)]


def test_the_same_sequence_in_two_regions_is_not_merged_across_regions(tmp_path):
    """Grouping is done per region: one sequence in ITS1 and in ITS2 gives a row in each."""
    sequence = dna(120, 5)
    made = make_inputs(
        tmp_path, its1=[("w1", sequence)], its2=[("w1", sequence), ("w2", sequence)]
    )

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert [(row["region"], row["sequence_id"], row["copies"]) for row in rows] == [
        ("ITS1", "S1_ITS1_1", "1"),
        ("ITS2", "S1_ITS2_1", "2"),
    ]


def test_case_and_line_wrapping_do_not_make_sequences_different(tmp_path):
    """Capital, small and mixed letters, wrapped in different places, are ONE sequence.

    The three records below spell the same 60 bases. The first is on one line in
    capitals, the second is on one line in small letters, the third is mixed case and
    wrapped at 17 columns with blanks at the ends of some lines.
    """
    sequence = dna(60, 3)
    mixed = sequence[:10].lower() + sequence[10:40] + sequence[40:].lower()
    wrapped_lines = [mixed[start:start + 17] for start in range(0, len(mixed), 17)]
    text = (
        f">w1\n{sequence}\n"
        f">w2\n{sequence.lower()}\n"
        f">w3\n" + "\n".join(wrapped_lines[:-1]) + " \n" + wrapped_lines[-1] + "\n"
    )
    made = make_inputs(tmp_path, its2=text)

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert len(rows) == 1
    assert rows[0]["copies"] == "3"
    assert rows[0]["length"] == "60"
    # The classifier is given the upper-case sequence on a single line.
    assert made.out_its2.read_text() == f">S1_ITS2_1\n{sequence}\n"


def test_windows_line_endings_are_understood(tmp_path):
    """A file written on Windows (lines end with \\r\\n) is read like any other."""
    sequence = dna(60, 4)
    made = make_inputs(tmp_path, its2=f">w1\r\n{sequence}\r\n>w2\r\n{sequence}\r\n")

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert rows[0]["copies"] == "2"
    assert rows[0]["windows"] == "w1,w2"


def test_window_name_is_the_first_word_of_the_header(tmp_path):
    """Anything after the first blank in a header line is a description and is not kept."""
    made = make_inputs(
        tmp_path, its2=[("ptg1:1-4001_+ some description words", dna(100, 2))]
    )

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert rows[0]["windows"] == "ptg1:1-4001_+"


# ======================================================================
# THE LENGTH GUARD
# ======================================================================

def test_too_long_its2_is_rejected_and_absent_from_the_classifier_fasta(tmp_path):
    """A mis-extracted 80,000-base "ITS2" is listed as rejected; the good ITS2 passes.

    Without the length guard such a record would go to the classifier, which would
    classify a long stretch of genome as if it were an ITS2. Here it must stay out of
    the FASTA that goes to the classifier but remain visible in the table.
    """
    good = dna(170, 20)
    giant = dna(80000, 21)
    made = make_inputs(
        tmp_path, its2=[("contig_1:1-80000_+", giant), ("contig_4:1-4001_-", good)]
    )

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert len(rows) == 2

    rejected_row = rows[0]
    assert rejected_row["sequence_id"] == "-"
    assert rejected_row["length"] == "80000"
    assert rejected_row["copies"] == "1"
    assert rejected_row["status"] == "rejected_too_long"
    assert rejected_row["windows"] == "contig_1:1-80000_+"

    good_row = rows[1]
    assert good_row["sequence_id"] == "S1_ITS2_1"      # the rejected one did not use up a number
    assert good_row["status"] == "ok"

    assert read_fasta_back(made.out_its2) == [("S1_ITS2_1", good)]


def test_too_short_part_is_rejected(tmp_path):
    """An ITS1 of 49 bases is shorter than the 50 bases of the guard."""
    made = make_inputs(tmp_path, its1=[("w1", dna(49, 30))])

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert rows[0]["status"] == "rejected_too_short"
    assert rows[0]["sequence_id"] == "-"
    assert rows[0]["length"] == "49"
    assert made.out_its1.read_text() == ""


@pytest.mark.parametrize("region", ["ITS1", "ITS2"])
@pytest.mark.parametrize(
    "length, expected_status",
    [
        (49, "rejected_too_short"),
        (50, "ok"),
        (600, "ok"),
        (601, "rejected_too_long"),
    ],
)
def test_part_length_limits_include_50_and_600(tmp_path, region, length, expected_status):
    """The limits are inclusive: 50 and 600 pass, 49 and 601 do not (ITS1 and ITS2 alike)."""
    records = [("w1", dna(length, length))]
    if region == "ITS1":
        made = make_inputs(tmp_path, its1=records)
    else:
        made = make_inputs(tmp_path, its2=records)

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert len(rows) == 1
    assert rows[0]["region"] == region
    assert rows[0]["status"] == expected_status


@pytest.mark.parametrize(
    "length, expected_status",
    [
        (1, "ok"),          # the full region has NO minimum
        (10, "ok"),
        (490, "ok"),
        (1500, "ok"),       # the limit itself is accepted
        (1501, "rejected_too_long"),
    ],
)
def test_full_region_has_a_maximum_but_no_minimum(tmp_path, length, expected_status):
    """The full region may be as short as it likes, but at most 1500 bases long."""
    made = make_inputs(tmp_path, full=[("w1", dna(length, 100 + length))])

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert len(rows) == 1
    assert rows[0]["region"] == "full"
    assert rows[0]["status"] == expected_status
    if expected_status == "ok":
        assert rows[0]["sequence_id"] == "S1_full_1"
    else:
        assert rows[0]["sequence_id"] == "-"


def test_full_region_is_only_reported_in_the_table(tmp_path):
    """The full-length sequence is not written to either classifier FASTA."""
    made = make_inputs(tmp_path, full=[("w1", dna(490, 40))])

    run_script(made)

    assert made.out_its1.read_text() == ""
    assert made.out_its2.read_text() == ""
    all_lines, rows = read_table(made.out_table)
    assert [row["region"] for row in rows] == ["full"]


def test_limits_can_be_changed_on_the_command_line(tmp_path):
    """With 100-200 for the parts and 300 for the full region, the verdicts move."""
    made = make_inputs(
        tmp_path,
        its1=[("w1", dna(150, 50))],      # inside 100-200: ok
        its2=[("w1", dna(250, 51))],      # above 200: too long
        full=[("w1", dna(350, 52))],      # above 300: too long
    )

    run_script(
        made, "--min-part-length", "100", "--max-part-length", "200", "--max-full-length", "300"
    )

    all_lines, rows = read_table(made.out_table)
    assert [(row["region"], row["status"]) for row in rows] == [
        ("ITS1", "ok"),
        ("ITS2", "rejected_too_long"),
        ("full", "rejected_too_long"),
    ]


def test_ids_count_only_sequences_that_pass_the_guard(tmp_path):
    """ok, too long, ok: the ids are _1, "-", _2 (rejected sequences leave no gap)."""
    seq_1 = dna(100, 60)
    seq_long = dna(700, 61)
    seq_2 = dna(110, 62)
    made = make_inputs(
        tmp_path, its1=[("w1", seq_1), ("w2", seq_long), ("w3", seq_2), ("w4", seq_long)]
    )

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert [row["sequence_id"] for row in rows] == ["S1_ITS1_1", "-", "S1_ITS1_2"]
    assert [row["copies"] for row in rows] == ["1", "2", "1"]
    assert [row["status"] for row in rows] == ["ok", "rejected_too_long", "ok"]
    assert [header for header, sequence in read_fasta_back(made.out_its1)] == [
        "S1_ITS1_1", "S1_ITS1_2"
    ]


def test_run_that_rejects_records_still_succeeds(tmp_path, capsys):
    """Rejecting every record is a result, not an error: normal return, empty FASTA, rows listed."""
    made = make_inputs(tmp_path, its2=[("w1", dna(700, 70))])

    run_script(made)   # must return normally: no SystemExit, so the exit code is 0

    assert made.out_its2.read_text() == ""
    all_lines, rows = read_table(made.out_table)
    assert [row["status"] for row in rows] == ["rejected_too_long"]
    # Some record existed, so this is NOT the "nothing extracted" case: no comment line.
    assert not any(line.startswith("#") for line in all_lines)
    assert "1 rejected" in capsys.readouterr().err


# ======================================================================
# THE FORMAT OF THE TABLE AND OF THE FASTA FILES
# ======================================================================

def test_table_format(tmp_path):
    """Header line, six tab-separated fields on every row, numbers that are numbers."""
    made = make_inputs(
        tmp_path,
        its1=records_for(window_names("c", 3), dna(100, 80)),
        its2=records_for(window_names("c", 3), dna(120, 81)) + [("x", dna(900, 82))],
        full=records_for(window_names("c", 3), dna(500, 83)),
    )

    run_script(made)

    text = made.out_table.read_text()
    assert text.endswith("\n")
    lines = text.splitlines()
    assert lines[0] == TABLE_HEADER_LINE

    allowed_statuses = {"ok", "rejected_too_short", "rejected_too_long"}
    for line in lines[1:]:
        fields = line.split("\t")
        assert len(fields) == 6
        region, sequence_id, length, copies, status, windows = fields
        assert region in {"ITS1", "ITS2", "full"}
        assert int(length) > 0
        assert int(copies) > 0
        assert status in allowed_statuses
        assert len(windows.split(",")) == int(copies)
        if status == "ok":
            assert sequence_id.startswith("S1_" + region + "_")
        else:
            assert sequence_id == "-"

    # No comment line when there are records.
    assert not any(line.startswith("#") for line in lines)


def test_rows_are_ordered_its1_then_its2_then_full_then_by_first_appearance(tmp_path):
    """The order is by region first, and inside a region by the first appearance."""
    made = make_inputs(
        tmp_path,
        its1=[("a", dna(100, 90)), ("b", dna(101, 91))],
        its2=[("a", dna(110, 92)), ("b", dna(111, 93)), ("c", dna(112, 94))],
        full=[("a", dna(300, 95)), ("b", dna(301, 96))],
    )

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert [(row["region"], row["length"]) for row in rows] == [
        ("ITS1", "100"), ("ITS1", "101"),
        ("ITS2", "110"), ("ITS2", "111"), ("ITS2", "112"),
        ("full", "300"), ("full", "301"),
    ]


def test_outputs_are_overwritten_and_not_appended(tmp_path):
    """Running twice on the same files leaves exactly the same outputs."""
    made = make_inputs(tmp_path, its2=records_for(window_names("c", 4), dna(130, 97)))

    run_script(made)
    first_table = made.out_table.read_text()
    first_fasta = made.out_its2.read_text()
    run_script(made)

    assert made.out_table.read_text() == first_table
    assert made.out_its2.read_text() == first_fasta


# ======================================================================
# NO RECORD AT ALL
# ======================================================================

def check_nothing_extracted(made, capsys):
    """Shared checks for a run in which all three input files hold no record.

    Expected: both FASTA files empty, a table of exactly two lines (the header and one
    comment line that starts with # and says that no ITS region was extracted), and
    one line on stderr.
    """
    assert made.out_its1.exists() and made.out_its1.read_text() == ""
    assert made.out_its2.exists() and made.out_its2.read_text() == ""

    lines = made.out_table.read_text().splitlines()
    assert len(lines) == 2
    assert lines[0] == TABLE_HEADER_LINE
    assert lines[1].startswith("#")
    assert "no ITS region was extracted" in lines[1]
    assert "5.8S" in lines[1]

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "no ITS region was extracted" in captured.err
    assert captured.err.count("\n") == 1


def test_three_empty_files_give_a_header_a_comment_line_and_empty_outputs(tmp_path, capsys):
    """A genome without a recognisable 5.8S must not make the workflow fail."""
    made = make_inputs(tmp_path)   # three empty files

    run_script(made)   # must return normally: no SystemExit, so the exit code is 0

    check_nothing_extracted(made, capsys)


def test_files_with_only_blank_lines_count_as_empty(tmp_path, capsys):
    """Blank lines carry no record."""
    made = make_inputs(tmp_path, its1="\n\n", its2="   \n", full="\n")

    run_script(made)

    check_nothing_extracted(made, capsys)


@pytest.mark.parametrize("region", ["its1", "its2", "full"])
def test_one_region_with_records_is_not_the_nothing_extracted_case(tmp_path, region):
    """It is enough that ONE of the three files holds a record: then there is no comment line."""
    made = make_inputs(tmp_path, **{region: [("w1", dna(100, 110))]})

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert not any(line.startswith("#") for line in all_lines)
    assert len(rows) == 1


# ======================================================================
# A REALISTIC CASE: FIVE rDNA COPIES (LONG READS) OR ONE (SHORT READS)
# ======================================================================

def test_the_example_sequences_have_the_expected_lengths():
    """Sanity check of the data written into this test file: ITS1 165, ITS2 170, full 490."""
    assert len(EXAMPLE_ITS1) == 165
    assert len(EXAMPLE_ITS2) == 170
    assert len(EXAMPLE_FULL) == 490


def test_five_identical_its2_records_give_one_row_with_five_copies(tmp_path):
    """Five identical 170-base ITS2 records are ONE sequence with copies 5."""
    names = window_names("ptg000001l", 5)
    made = make_inputs(tmp_path, its2=records_for(names, EXAMPLE_ITS2))

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert rows == [{
        "region": "ITS2",
        "sequence_id": "S1_ITS2_1",
        "length": "170",
        "copies": "5",
        "status": "ok",
        "windows": ",".join(names),
    }]
    assert read_fasta_back(made.out_its2) == [("S1_ITS2_1", EXAMPLE_ITS2)]


def test_long_read_result_is_five_copies_of_each_region(tmp_path):
    """A long-read assembly: five identical rDNA copies, 165 / 170 / 490 bases, all ok."""
    names = window_names("ptg000001l", 5)
    made = make_inputs(
        tmp_path,
        its1=records_for(names, EXAMPLE_ITS1),
        its2=records_for(names, EXAMPLE_ITS2),
        full=records_for(names, EXAMPLE_FULL),
    )

    run_script(made)

    joined = ",".join(names)
    expected_lines = [
        TABLE_HEADER_LINE,
        f"ITS1\tS1_ITS1_1\t165\t5\tok\t{joined}",
        f"ITS2\tS1_ITS2_1\t170\t5\tok\t{joined}",
        f"full\tS1_full_1\t490\t5\tok\t{joined}",
    ]
    assert made.out_table.read_text().splitlines() == expected_lines
    assert read_fasta_back(made.out_its1) == [("S1_ITS1_1", EXAMPLE_ITS1)]
    assert read_fasta_back(made.out_its2) == [("S1_ITS2_1", EXAMPLE_ITS2)]


def test_short_read_result_is_one_copy_of_each_region(tmp_path):
    """The short-read assembly folds the array into one unit: one copy, the same sequences."""
    name = "NODE_1_length_123456_cov_10.512345:10001-14001_+"
    made = make_inputs(
        tmp_path,
        its1=[(name, EXAMPLE_ITS1)], its2=[(name, EXAMPLE_ITS2)], full=[(name, EXAMPLE_FULL)],
    )

    run_script(made)

    all_lines, rows = read_table(made.out_table)
    assert [(row["region"], row["length"], row["copies"]) for row in rows] == [
        ("ITS1", "165", "1"), ("ITS2", "170", "1"), ("full", "490", "1"),
    ]
    assert read_fasta_back(made.out_its2) == [("S1_ITS2_1", EXAMPLE_ITS2)]


# ======================================================================
# BAD INPUT GIVES A CLEAR ERROR
# ======================================================================

def test_missing_input_file_is_a_clear_error(tmp_path, capsys):
    """An input file that does not exist gives a one-line message, not a traceback."""
    made = make_inputs(tmp_path, its1=[("w1", dna(100, 120))])
    made.its2.unlink()

    with pytest.raises(SystemExit) as stop:
        run_script(made)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("its_distinct.py: error:")
    assert "S1.ITS2.fasta" in message
    assert message.count("\n") == 1   # one line, no Python traceback

    # Nothing is written when the input is bad.
    assert not made.out_its1.exists()
    assert not made.out_its2.exists()
    assert not made.out_table.exists()


def test_sequence_before_any_header_is_a_clear_error(tmp_path, capsys):
    """A file that starts with sequence and not with ">" is not FASTA: stop and say where.

    The ITS1 and ITS2 files are fine; the full-length file starts with a sequence line.
    """
    made = make_inputs(
        tmp_path,
        its1=[("w1", dna(100, 121))],
        its2=[("w1", dna(100, 122))],
        full="ACGTACGTACGT\n>w1\nACGT\n",
    )

    with pytest.raises(SystemExit) as stop:
        run_script(made)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("its_distinct.py: error:")
    assert "S1.full.fasta" in message
    assert "line 1" in message
    assert "before any FASTA header" in message
    assert not made.out_table.exists()
    assert not made.out_its1.exists()


def test_header_without_a_name_is_a_clear_error(tmp_path, capsys):
    """A line with only ">" names no window, so the table could not say where a record came from."""
    made = make_inputs(tmp_path, its2=">w1\nACGT\n>\nACGT\n")

    with pytest.raises(SystemExit) as stop:
        run_script(made)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert "line 3" in message
    assert "no name" in message
    assert not made.out_table.exists()


@pytest.mark.parametrize(
    "bad_options, expected_word",
    [
        (["--min-part-length", "0"], "--min-part-length"),
        (["--min-part-length", "700"], "--max-part-length"),    # larger than the default maximum, 600
        (["--max-part-length", "10"], "--max-part-length"),     # smaller than the default minimum, 50
        (["--max-full-length", "0"], "--max-full-length"),
    ],
    ids=["min_is_zero", "min_above_max", "max_below_min", "full_limit_is_zero"],
)
def test_impossible_limits_are_refused_by_the_command_line(tmp_path, capsys, bad_options, expected_word):
    """Limits that would reject everything stop the script before anything is read or written."""
    made = make_inputs(tmp_path, its2=[("w1", dna(100, 122))])

    with pytest.raises(SystemExit) as stop:
        run_script(made, *bad_options)

    assert stop.value.code == 2   # argparse's code for a bad command line
    assert expected_word in capsys.readouterr().err
    assert not made.out_table.exists()


def test_missing_required_option_is_refused(tmp_path, capsys):
    """Leaving out --sample is a command-line error (exit code 2)."""
    made = make_inputs(tmp_path)
    words = command_line(made)
    position = words.index("--sample")
    del words[position:position + 2]

    with pytest.raises(SystemExit) as stop:
        its_distinct.main(words)

    assert stop.value.code == 2
    assert "--sample" in capsys.readouterr().err


# ======================================================================
# THE REPORT ON STDERR
# ======================================================================

def test_summary_line_goes_to_stderr_and_nothing_to_stdout(tmp_path, capsys):
    """One line on stderr says how many records, distinct sequences, ok and rejected per region."""
    made = make_inputs(
        tmp_path,
        its1=records_for(window_names("c", 5), dna(100, 130)),
        its2=records_for(window_names("c", 5), dna(120, 131)) + [("x", dna(800, 132))],
        full=records_for(window_names("c", 5), dna(500, 133)),
    )

    run_script(made)

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.count("\n") == 1
    assert "sample S1" in captured.err
    assert "ITS1 5 record(s), 1 distinct, 1 ok, 0 rejected" in captured.err
    assert "ITS2 6 record(s), 2 distinct, 1 ok, 1 rejected" in captured.err
    assert "full 5 record(s), 1 distinct, 1 ok, 0 rejected" in captured.err


# ======================================================================
# THE SCRIPT AS SNAKEMAKE WILL CALL IT
# ======================================================================

def run_from_command_line(made, *extra_options):
    """Run the script file in a separate Python process, as a Snakemake rule would."""
    command = [sys.executable, SCRIPT_PATH] + command_line(made, *extra_options)
    return subprocess.run(command, capture_output=True, text=True)


def test_command_line_run_succeeds(tmp_path):
    """Exit code 0, the table and the FASTA written, the summary on stderr."""
    made = make_inputs(tmp_path, its2=records_for(window_names("ptg1", 5), EXAMPLE_ITS2))

    finished = run_from_command_line(made)

    assert finished.returncode == 0
    assert finished.stdout == ""
    assert "ITS2 5 record(s), 1 distinct, 1 ok, 0 rejected" in finished.stderr
    assert read_fasta_back(made.out_its2) == [("S1_ITS2_1", EXAMPLE_ITS2)]


def test_command_line_run_without_any_record_still_exits_with_code_0(tmp_path):
    """The workflow must not fail on a genome without a recognisable 5.8S."""
    made = make_inputs(tmp_path)

    finished = run_from_command_line(made)

    assert finished.returncode == 0
    assert "no ITS region was extracted" in finished.stderr
    assert made.out_its1.read_text() == ""
    assert made.out_table.read_text().splitlines()[0] == TABLE_HEADER_LINE


def test_command_line_run_with_a_missing_file_exits_with_code_1(tmp_path):
    """A missing input fails the step (non-zero exit) with a readable message."""
    made = make_inputs(tmp_path)
    made.full.unlink()

    finished = run_from_command_line(made)

    assert finished.returncode == 1
    assert "its_distinct.py: error:" in finished.stderr
    assert "S1.full.fasta" in finished.stderr
    assert "Traceback" not in finished.stderr


def test_command_line_run_that_rejects_a_record_exits_with_code_0(tmp_path):
    """The length guard is a verdict listed in the table, not a failure of the step."""
    made = make_inputs(tmp_path, its2=[("w1", dna(900, 140))])

    finished = run_from_command_line(made)

    assert finished.returncode == 0
    assert "rejected_too_long" in made.out_table.read_text()
