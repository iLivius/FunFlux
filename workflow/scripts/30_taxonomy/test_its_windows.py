"""Tests for its_windows.py, the script that cuts one window around every 5.8S gene.

How to run (from the repository root):

    python -m pytest workflow/scripts/30_taxonomy -q

Only the Python standard library and pytest are needed. No real genome is needed
either: every input is made up while the tests run.

THE MADE-UP DATA
----------------
  * A synthetic contig, "contig_1", of 30,000 bases of random sequence. It carries
    two made-up rDNA units of 4,001 bases each (the size of one full window):
        - the "plus unit" is placed as it is, at bases 8079-12079;
        - the "minus unit" is the SAME unit reverse-complemented, at 18079-22079.
    Because both copies are the same unit, a correct run gives two windows that are
    identical to each other and to the original unit: the minus window must be turned
    back into the canonical orientation.
    The unit also holds a stretch of lowercase (soft-masked) bases and a run of N, to
    check that case and N survive the cutting and the reverse complement.
  * A hand-written nhmmer table (NHMMER_TABLE below) with several anchor hits for each
    unit, listed in no particular order, as nhmmer does. They have to be merged into
    one 5.8S site per unit. The table also holds one weak hit (E-value 4.0e-04, above
    the default limit of 1e-5) that must be ignored.
  * A few extra, tiny contigs and tables for the special cases (clipping at the
    contig ends, merging distance, strands, errors).
"""

import os
import random
import subprocess
import sys
from types import SimpleNamespace

import pytest

# SimpleNamespace (used by make_inputs below) is a plain container whose fields are read
# with a dot, like made.fasta.

# The script sits in the same folder as this test file. Putting that folder on the
# list of places where Python looks for modules makes "import its_windows" work
# wherever pytest is started from.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TEST_DIR)

import its_windows  # noqa: E402  (this import must come after the line above)

SCRIPT_PATH = os.path.join(TEST_DIR, "its_windows.py")


# ======================================================================
# THE MADE-UP DATA
# ======================================================================

CONTIG = "contig_1"
CONTIG_LENGTH = 30000

# With the default flank of 2000 bases a full window is 2 * 2000 + 1 = 4001 bases long.
WINDOW_LENGTH = 4001

# Where the two made-up rDNA units sit on contig_1 (counted from 1, both ends included).
PLUS_UNIT_START = 8079     # the unit as it is, bases 8079-12079
MINUS_UNIT_START = 18079   # the unit reverse-complemented, bases 18079-22079

# The columns of the output table, exactly as the script must write them.
TABLE_HEADER_LINE = (
    "window_id\tcontig\tstart\tend\tstrand\tanchor_hits\tbest_evalue\twindow_length"
)

# The two comment lines that nhmmer writes above the hits, and the lines below them.
# The script has to skip every line that starts with "#".
NHMMER_COMMENTS_TOP = (
    "# target name   accession  query name   accession  hmmfrom hmm to alifrom  ali to"
    " envfrom  env to  sq len strand   E-value  score  bias  description of target\n"
    "#  ------------ ---------- ------------ ---------- ------- ------- -------"
    " ------- ------- ------- ------- ------ --------- ------ ----- ------------------\n"
)
NHMMER_COMMENTS_BOTTOM = (
    "#\n"
    "# Program:         nhmmer\n"
    "# Pipeline mode:   SEARCH\n"
    "# [ok]\n"
)

# The hand-written table for contig_1 (30,000 bases long). The columns are, counting
# from 0: 0 contig, 1 accession, 2 profile, 3 accession, 4 hmmfrom, 5 hmmto, 6 alifrom,
# 7 alito, 8 envfrom, 9 envto, 10 contig length, 11 strand, 12 E-value, 13 score,
# 14 bias, 15 description. The description is the rest of the FASTA header line, here
# "length=30000 cov=55.2", so it can hold spaces; the script must not mind that.
#
# PLUS UNIT: six hits whose envelopes (columns 8 and 9) together cover 10001-10158,
#   so the middle of the site is (10001 + 10158) // 2 = 10079 and the window is
#   10079 +- 2000 = 8079-12079. The smallest E-value is 3.1e-12.
#   The alignment columns (6 and 7) are deliberately a little narrower than the
#   envelope (they start at 10002), so a script that used them by mistake would find
#   the middle 10080 and be caught.
# MINUS UNIT: four hits on the minus strand. nhmmer lists the coordinates of a minus
#   hit backwards (from > to). Their envelopes cover 20001-20158, so the middle is
#   20079 and the window is 18079-22079. The smallest E-value is 4.0e-13.
# WEAK HIT: E-value 4.0e-04 at 26001-26045 on the plus strand. It is above 1e-5, so it
#   must be ignored: no third window.
NHMMER_HIT_ROWS = """\
contig_1 - start_zero   -  2 45 20045 20002 20045 20001 30000 - 4.0e-13 54.1 0.5 length=30000 cov=55.2
contig_1 - end_seven    -  1 45 10114 10158 10114 10158 30000 + 2.0e-09 42.1 7.4 length=30000 cov=55.2
contig_1 - start_thirty -  1 45 10002 10045 10002 10045 30000 + 3.1e-12 51.2 0.4 length=30000 cov=55.2
contig_1 - end_ten      -  1 45 20158 20114 20158 20114 30000 - 9.8e-10 43.0 7.1 length=30000 cov=55.2
contig_1 - start_weak   -  2 45 26002 26045 26001 26045 30000 + 4.0e-04 24.7 0.6 length=30000 cov=55.2
contig_1 - start_cr     -  2 45 10002 10045 10001 10045 30000 + 1.2e-10 45.8 0.4 length=30000 cov=55.2
contig_1 - end_nine     -  1 44 20157 20114 20158 20114 30000 - 1.1e-07 35.2 8.0 length=30000 cov=55.2
contig_1 - end_fifteen  - 15 45 10128 10158 10115 10158 30000 + 5.5e-06 30.1 4.2 length=30000 cov=55.2
contig_1 - start_29     -  4 45 10004 10045 10003 10045 30000 + 4.5e-08 38.0 0.6 length=30000 cov=55.2
contig_1 - start_two    -  4 45 20043 20002 20044 20001 30000 - 2.9e-08 38.2 0.4 length=30000 cov=55.2
contig_1 - end_foura    -  1 44 10114 10157 10114 10158 30000 + 8.0e-07 33.0 8.0 length=30000 cov=55.2
"""

NHMMER_TABLE = NHMMER_COMMENTS_TOP + NHMMER_HIT_ROWS + NHMMER_COMMENTS_BOTTOM


def random_dna(rng, length):
    """Return random A/C/G/T text. rng is a seeded random.Random, so every run is identical."""
    letters = rng.choices("ACGT", k=length)
    return "".join(letters)


def reverse_complement_simple(sequence):
    """Reverse complement written the slow, obvious way, independent of the script.

    The tests use this to work out what the script SHOULD produce, so that the script
    is not checked against itself. It only knows A, C, G, T, N and their lowercase
    forms, which is all the made-up sequences contain.
    """
    partner = {
        "A": "T", "C": "G", "G": "C", "T": "A", "N": "N",
        "a": "t", "c": "g", "g": "c", "t": "a", "n": "n",
    }
    partner_bases = []
    for base in sequence:
        partner_bases.append(partner[base])
    partner_bases.reverse()   # the opposite strand is read in the opposite direction
    return "".join(partner_bases)


def build_synthetic_contig():
    """Build contig_1 and return (contig_sequence, canonical_unit).

    The canonical unit is the 4,001-base rDNA unit in its proper orientation. The
    contig holds it once as it is (plus strand) and once reverse-complemented
    (minus strand), with random sequence in between and around.
    """
    rng = random.Random(20260707)

    unit = random_dna(rng, WINDOW_LENGTH)
    # Make part of the unit lowercase and put a run of N in it, so that the tests
    # also check that case and N are kept. The total length stays 4,001.
    unit = unit[:300] + unit[300:360].lower() + "NNNNN" + unit[365:]

    before_plus_unit = random_dna(rng, PLUS_UNIT_START - 1)
    between_the_units = random_dna(rng, MINUS_UNIT_START - (PLUS_UNIT_START + WINDOW_LENGTH))
    after_minus_unit = random_dna(rng, CONTIG_LENGTH - (MINUS_UNIT_START + WINDOW_LENGTH - 1))

    contig = (
        before_plus_unit
        + unit
        + between_the_units
        + reverse_complement_simple(unit)
        + after_minus_unit
    )
    assert len(contig) == CONTIG_LENGTH
    return contig, unit


def hit_line(contig, contig_length, strand, low, high, evalue):
    """Return one nhmmer table line for a hit that covers bases low-high of a contig.

    For the minus strand nhmmer writes the coordinates backwards (from > to), so this
    helper does the same. Alignment and envelope coordinates are written identically.
    """
    if strand == "+":
        coordinate_from, coordinate_to = low, high
    else:
        coordinate_from, coordinate_to = high, low

    columns = [
        contig, "-", "some_anchor_profile", "-", "1", "45",
        str(coordinate_from), str(coordinate_to), str(coordinate_from), str(coordinate_to),
        str(contig_length), strand, f"{evalue:g}", "40.0", "1.0", "-",
    ]
    return " ".join(columns) + "\n"


def make_table(*hit_lines):
    """Return the text of an nhmmer table that holds the given hit lines."""
    return NHMMER_COMMENTS_TOP + "".join(hit_lines) + NHMMER_COMMENTS_BOTTOM


def write_fasta_file(path, records, line_width=70):
    """Write a FASTA file from a list of (header, sequence) pairs.

    The input is wrapped at 70 columns, not at the 80 that the script writes, so the
    tests also show that multi-line input is read correctly.
    """
    with open(path, "w") as fasta:
        for header, sequence in records:
            fasta.write(f">{header}\n")
            for line_start in range(0, len(sequence), line_width):
                fasta.write(sequence[line_start:line_start + line_width] + "\n")


def make_inputs(folder, records, table_text):
    """Write a genome FASTA and an nhmmer table into `folder` and return their paths.

    records is a list of (header, sequence) pairs. The returned object also holds the
    paths where the script is to write its two outputs.
    """
    fasta_path = folder / "genome.fasta"
    tblout_path = folder / "58S.tblout"
    write_fasta_file(fasta_path, records)
    tblout_path.write_text(table_text)
    return SimpleNamespace(
        fasta=fasta_path,
        tblout=tblout_path,
        out_fasta=folder / "windows.fasta",
        out_table=folder / "windows.tsv",
    )


@pytest.fixture
def inputs(tmp_path):
    """The standard made-up inputs: contig_1 with its two units, plus its nhmmer table.

    The genome also holds two unrelated contigs, one before and one after contig_1, to
    check that contigs without a 5.8S hit are skipped without trouble.
    """
    contig, unit = build_synthetic_contig()
    rng = random.Random(7)
    records = [
        ("other_contig_0 no rDNA here", random_dna(rng, 3000)),
        (f"{CONTIG} length={CONTIG_LENGTH} cov=55.2", contig),   # description after the name
        ("other_contig_2 no rDNA here either", random_dna(rng, 2500)),
    ]
    made = make_inputs(tmp_path, records, NHMMER_TABLE)
    made.contig = contig
    made.unit = unit
    return made


def run_script(made, *extra_options):
    """Run its_windows.main() on the files in `made`, adding any extra options."""
    command_line = [
        "--tblout", str(made.tblout),
        "--fasta", str(made.fasta),
        "--out-fasta", str(made.out_fasta),
        "--out-table", str(made.out_table),
    ]
    command_line.extend(extra_options)
    its_windows.main(command_line)


def read_output_fasta(path):
    """Read the script's FASTA output back. Returns a list with one dictionary per window.

    Each dictionary holds the header (without ">"), the list of sequence lines, and
    the full sequence.
    """
    windows = []
    with open(path) as fasta:
        for line in fasta:
            line = line.rstrip("\n")
            if line.startswith(">"):
                windows.append({"header": line[1:], "lines": []})
            else:
                windows[-1]["lines"].append(line)

    for window in windows:
        window["sequence"] = "".join(window["lines"])
    return windows


def read_output_table(path):
    """Read the script's table output back. Returns (list_of_column_names, list_of_rows).

    Each row is a dictionary from column name to the text found in that column.
    """
    with open(path) as table:
        lines = table.read().splitlines()

    column_names = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        values = line.split("\t")
        row = {}
        for position, column_name in enumerate(column_names):
            row[column_name] = values[position]
        rows.append(row)
    return column_names, rows


def headers_of(path):
    """Return the list of window headers in an output FASTA file."""
    headers = []
    for window in read_output_fasta(path):
        headers.append(window["header"])
    return headers


# ======================================================================
# REVERSE COMPLEMENT
# ======================================================================

def test_reverse_complement_handles_every_iupac_code_and_case():
    """Each base swaps with its partner; lowercase stays lowercase.

    The pairs are A/T, C/G, R/Y, K/M, B/V, D/H, while S, W and N are their own partner.
    """
    pairs = [
        ("A", "T"), ("C", "G"), ("G", "C"), ("T", "A"),
        ("R", "Y"), ("Y", "R"), ("K", "M"), ("M", "K"),
        ("B", "V"), ("V", "B"), ("D", "H"), ("H", "D"),
        ("S", "S"), ("W", "W"), ("N", "N"),
    ]
    for base, partner in pairs:
        assert its_windows.reverse_complement(base) == partner
        assert its_windows.reverse_complement(base.lower()) == partner.lower()

    # All 15 codes in one string: complement is TGCAYRMKVBHDSWN, then read backwards.
    assert its_windows.reverse_complement("ACGTRYKMBVDHSWN") == "NWSDHBVKMRYACGT"
    assert its_windows.reverse_complement("acgtrykmbvdhswn") == "nwsdhbvkmryacgt"

    # Mixed case: complement is TgCaNn, then read backwards.
    assert its_windows.reverse_complement("AcGtNn") == "nNaCgT"


def test_reverse_complement_twice_gives_back_the_original():
    """Reading the opposite strand of the opposite strand must return the first strand."""
    sequence = "ACGTRYKMBVDHSWNacgtrykmbvdhswnAACCGGTT"
    once = its_windows.reverse_complement(sequence)
    twice = its_windows.reverse_complement(once)
    assert once != sequence
    assert twice == sequence


def test_reverse_complement_leaves_other_characters_alone():
    """A character that is not a DNA letter (a gap, a star, an X) is kept unchanged."""
    # Complement is TX-A*, then read backwards.
    assert its_windows.reverse_complement("AX-T*") == "*A-XT"


# ======================================================================
# THE MAIN CASE: TWO RDNA UNITS, ONE ON EACH STRAND
# ======================================================================

def test_two_units_give_two_windows_and_the_weak_hit_is_ignored(inputs):
    """Ten good hits fall into two sites (one per unit); the weak hit makes no third window."""
    run_script(inputs)

    windows = read_output_fasta(inputs.out_fasta)
    assert len(windows) == 2


def test_window_coordinates_and_headers(inputs):
    """Plus window 8079-12079, minus window 18079-22079; headers end in _+ or _-."""
    run_script(inputs)

    headers = headers_of(inputs.out_fasta)
    assert headers == ["contig_1:8079-12079_+", "contig_1:18079-22079_-"]

    # The strand must follow an underscore. ITSx mishandles parentheses in headers.
    for header in headers:
        assert "(" not in header
        assert ")" not in header
        assert header.endswith("_+") or header.endswith("_-")


def test_plus_strand_window_is_the_contig_slice(inputs):
    """A plus-strand window is the contig text from base 8079 to base 12079, unchanged."""
    run_script(inputs)

    windows = read_output_fasta(inputs.out_fasta)

    # Python slices count from 0 and stop before the end, so bases 8079-12079 counted
    # from 1 are the slice [8078:12079].
    slice_of_contig = inputs.contig[8078:12079]
    assert windows[0]["sequence"] == slice_of_contig


def test_minus_strand_window_is_the_exact_reverse_complement_of_the_contig_slice(inputs):
    """A minus-strand window is the other strand of contig bases 18079-22079."""
    run_script(inputs)

    windows = read_output_fasta(inputs.out_fasta)

    slice_of_contig = inputs.contig[18078:22079]
    expected = reverse_complement_simple(slice_of_contig)
    assert windows[1]["sequence"] == expected
    assert windows[1]["sequence"] != slice_of_contig   # it really was flipped


def test_both_windows_hold_the_same_unit_in_the_canonical_orientation(inputs):
    """The made-up unit was placed once as it is and once reverse-complemented.

    After the script has turned the minus window around, both windows must be the
    original unit, including its lowercase stretch and its run of N.
    """
    run_script(inputs)

    windows = read_output_fasta(inputs.out_fasta)

    assert windows[0]["sequence"] == inputs.unit
    assert windows[1]["sequence"] == inputs.unit
    assert "NNNNN" in windows[1]["sequence"]
    assert windows[1]["sequence"].lower() != windows[1]["sequence"]   # lowercase kept


def test_fasta_sequence_lines_are_wrapped_at_80_columns(inputs):
    """4,001 bases are written as 50 lines of 80 bases and a last line of 1 base."""
    run_script(inputs)

    windows = read_output_fasta(inputs.out_fasta)
    for window in windows:
        line_lengths = []
        for line in window["lines"]:
            line_lengths.append(len(line))

        assert line_lengths == [80] * 50 + [1]


def test_table_has_the_expected_columns_and_values(inputs):
    """One row per window with the window id, coordinates, evidence and length."""
    run_script(inputs)

    column_names, rows = read_output_table(inputs.out_table)
    assert "\t".join(column_names) == TABLE_HEADER_LINE
    assert len(rows) == 2

    plus_row = rows[0]
    assert plus_row["window_id"] == "contig_1:8079-12079_+"
    assert plus_row["contig"] == "contig_1"
    assert plus_row["start"] == "8079"
    assert plus_row["end"] == "12079"
    assert plus_row["strand"] == "+"
    assert plus_row["anchor_hits"] == "6"           # six hits merged into the plus site
    assert plus_row["best_evalue"] == "3.10e-12"
    assert plus_row["window_length"] == "4001"

    minus_row = rows[1]
    assert minus_row["window_id"] == "contig_1:18079-22079_-"
    assert minus_row["contig"] == "contig_1"
    assert minus_row["start"] == "18079"
    assert minus_row["end"] == "22079"
    assert minus_row["strand"] == "-"
    assert minus_row["anchor_hits"] == "4"          # four hits merged into the minus site
    assert minus_row["best_evalue"] == "4.00e-13"
    assert minus_row["window_length"] == "4001"


def test_summary_line_goes_to_stderr(inputs, capsys):
    """The script reports on stderr how many windows it wrote, and nothing on stdout."""
    run_script(inputs)

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "wrote 2 window(s)" in captured.err
    assert captured.err.count("\n") == 1


# ======================================================================
# THE OPTIONS
# ======================================================================

def test_flank_option_changes_the_window_size(inputs):
    """With --flank 500 each window is the middle of the site +- 500, i.e. 1,001 bases."""
    run_script(inputs, "--flank", "500")

    # The middles are 10079 (plus site) and 20079 (minus site).
    headers = headers_of(inputs.out_fasta)
    assert headers == ["contig_1:9579-10579_+", "contig_1:19579-20579_-"]

    windows = read_output_fasta(inputs.out_fasta)
    assert windows[0]["sequence"] == inputs.contig[9578:10579]
    assert windows[1]["sequence"] == reverse_complement_simple(inputs.contig[19578:20579])


def test_max_evalue_option_lets_a_weaker_hit_in(inputs):
    """With --max-evalue 1e-3 the weak hit (4.0e-04) counts and forms a third window."""
    run_script(inputs, "--max-evalue", "1e-3")

    # The weak hit covers 26001-26045, so its middle is 26023 and its window is 24023-28023.
    headers = headers_of(inputs.out_fasta)
    assert headers == [
        "contig_1:8079-12079_+",
        "contig_1:18079-22079_-",
        "contig_1:24023-28023_+",
    ]

    column_names, rows = read_output_table(inputs.out_table)
    assert rows[2]["anchor_hits"] == "1"
    assert rows[2]["best_evalue"] == "4.00e-04"


def test_merge_distance_option_splits_sites_that_are_further_apart(inputs):
    """With --merge-distance 50, the start-anchor and end-anchor hits (69 bases apart) split.

    Each unit then gives two sites instead of one, so four windows in total.
    """
    run_script(inputs, "--merge-distance", "50")

    # Plus: start anchors cover 10001-10045 (middle 10023), end anchors 10114-10158 (10136).
    # Minus: start anchors cover 20001-20045 (middle 20023), end anchors 20114-20158 (20136).
    headers = headers_of(inputs.out_fasta)
    assert headers == [
        "contig_1:8023-12023_+",
        "contig_1:8136-12136_+",
        "contig_1:18023-22023_-",
        "contig_1:18136-22136_-",
    ]


def test_flank_below_one_is_refused(inputs):
    """A flank of 0 (or less) would cut nonsensical windows, so the command line refuses it."""
    with pytest.raises(SystemExit) as stop:
        run_script(inputs, "--flank", "0")

    assert stop.value.code == 2   # argparse's code for a bad command line
    assert not inputs.out_fasta.exists()


# ======================================================================
# MERGING AND FILTERING THE HITS
# ======================================================================

def test_hits_exactly_at_the_merge_distance_are_merged_and_one_base_further_are_not(tmp_path):
    """Default merge distance 1000: a gap of exactly 1000 merges, a gap of 1001 does not.

    The gap is the start of the next hit minus the end of the hit before it.
    """
    contig_name = "m"
    contig = random_dna(random.Random(1), 20000)
    records = [(contig_name, contig)]

    # Hit A covers 9001-9045. Hit B starts at 10045, which is 10045 - 9045 = 1000 later.
    # One site 9001-10089 with its middle at 9545: window 7545-11545, two hits merged.
    merged_table = make_table(
        hit_line(contig_name, 20000, "+", 9001, 9045, 1e-10),
        hit_line(contig_name, 20000, "+", 10045, 10089, 1e-10),
    )
    folder_a = tmp_path / "merged"
    folder_a.mkdir()
    made_a = make_inputs(folder_a, records, merged_table)
    run_script(made_a)

    assert headers_of(made_a.out_fasta) == ["m:7545-11545_+"]
    column_names, rows = read_output_table(made_a.out_table)
    assert rows[0]["anchor_hits"] == "2"

    # Hit B starts at 10046 this time, a gap of 1001, so there are two separate sites.
    # Site A: middle (9001 + 9045) // 2 = 9023, window 7023-11023.
    # Site B: middle (10046 + 10090) // 2 = 10068, window 8068-12068.
    split_table = make_table(
        hit_line(contig_name, 20000, "+", 9001, 9045, 1e-10),
        hit_line(contig_name, 20000, "+", 10046, 10090, 1e-10),
    )
    folder_b = tmp_path / "split"
    folder_b.mkdir()
    made_b = make_inputs(folder_b, records, split_table)
    run_script(made_b)

    assert headers_of(made_b.out_fasta) == ["m:7023-11023_+", "m:8068-12068_+"]


def test_hits_on_opposite_strands_are_never_merged_even_when_interleaved(tmp_path):
    """Hits are merged per strand: a plus hit, a minus hit and a plus hit, all close together.

    The two plus hits (9001-9045 and 9201-9245) form one site and the minus hit
    (9101-9145) forms another. Both sites have the middle 9123, so both windows are
    7123-11123, and they differ only by the strand. When the start and the end are
    equal, the plus site is listed first.
    """
    contig_name = "m"
    contig = random_dna(random.Random(2), 20000)
    table = make_table(
        hit_line(contig_name, 20000, "+", 9001, 9045, 1e-10),
        hit_line(contig_name, 20000, "-", 9101, 9145, 1e-10),
        hit_line(contig_name, 20000, "+", 9201, 9245, 1e-10),
    )
    made = make_inputs(tmp_path, [(contig_name, contig)], table)
    run_script(made)

    assert headers_of(made.out_fasta) == ["m:7123-11123_+", "m:7123-11123_-"]

    column_names, rows = read_output_table(made.out_table)
    assert rows[0]["anchor_hits"] == "2"
    assert rows[1]["anchor_hits"] == "1"


def test_merging_leaves_the_callers_hit_lists_unchanged():
    """merge_hits_into_sites() sorts a copy of each list of hits, never the list it is given.

    The two hits are handed over right to left. After the call the caller's list must
    still hold them in that order, and the merged site must still span both hits.
    """
    hits_right_to_left = [(9201, 9245, 1e-10), (9001, 9045, 1e-10)]
    hits_by_contig_and_strand = {("m", "+"): hits_right_to_left}
    order_before = list(hits_right_to_left)   # a copy, to compare with afterwards

    sites = its_windows.merge_hits_into_sites(hits_by_contig_and_strand, 1000)

    assert hits_by_contig_and_strand[("m", "+")] == order_before
    assert len(sites) == 1
    assert sites[0]["start"] == 9001
    assert sites[0]["end"] == 9245
    assert sites[0]["anchor_hits"] == 2


def test_hit_exactly_at_the_evalue_limit_is_kept(tmp_path):
    """The limit is inclusive: with the default 1e-5, a hit of 1e-5 counts but 1.1e-5 does not."""
    contig = random_dna(random.Random(9), 20000)

    # Hit at the limit: kept. It covers 9001-9045, so the middle is 9023 and the window
    # is 7023-11023.
    table_at_the_limit = make_table(hit_line("m", 20000, "+", 9001, 9045, 1e-5))
    made = make_inputs(tmp_path, [("m", contig)], table_at_the_limit)
    run_script(made)
    assert headers_of(made.out_fasta) == ["m:7023-11023_+"]

    # Hit just above the limit: ignored, so no window is left.
    table_above_the_limit = make_table(hit_line("m", 20000, "+", 9001, 9045, 1.1e-5))
    made.tblout.write_text(table_above_the_limit)
    run_script(made)
    assert headers_of(made.out_fasta) == []


def test_weak_hit_next_to_a_good_site_does_not_stretch_it(tmp_path):
    """Weak hits are dropped BEFORE merging, so they cannot move the middle of a site.

    The good hit covers 9001-9045 (middle 9023, window 7023-11023). The weak hit at
    9500-9544 lies within the merge distance; had it been merged, the site would end
    at 9544 and the window would be a different one.
    """
    contig_name = "m"
    contig = random_dna(random.Random(3), 20000)
    table = make_table(
        hit_line(contig_name, 20000, "+", 9001, 9045, 1e-10),
        hit_line(contig_name, 20000, "+", 9500, 9544, 4e-4),
    )
    made = make_inputs(tmp_path, [(contig_name, contig)], table)
    run_script(made)

    assert headers_of(made.out_fasta) == ["m:7023-11023_+"]
    column_names, rows = read_output_table(made.out_table)
    assert rows[0]["anchor_hits"] == "1"


# ======================================================================
# CLIPPING AND ORDER
# ======================================================================

def test_windows_are_clipped_at_the_contig_ends(tmp_path):
    """A window never runs off a contig; the header and the sequence show the clipped range.

    near_start (5,000 bases): plus hit 801-958, middle 879, so the window would start at
      879 - 2000 < 1; it is clipped to 1-2879.
    near_end (5,000 bases): minus hit 4401-4558, middle 4479, so the window would end at
      4479 + 2000 > 5000; it is clipped to 2479-5000 and then reverse-complemented.
    tiny (3,000 bases): plus hit 1401-1558, middle 1479; both ends are clipped, so the
      window is the whole contig, 1-3000.
    """
    rng = random.Random(4)
    near_start = random_dna(rng, 5000)
    near_end = random_dna(rng, 5000)
    tiny = random_dna(rng, 3000)
    records = [("near_start", near_start), ("near_end", near_end), ("tiny", tiny)]
    table = make_table(
        hit_line("tiny", 3000, "+", 1401, 1558, 1e-10),
        hit_line("near_end", 5000, "-", 4401, 4558, 1e-10),
        hit_line("near_start", 5000, "+", 801, 958, 1e-10),
    )
    made = make_inputs(tmp_path, records, table)
    run_script(made)

    assert headers_of(made.out_fasta) == [
        "near_start:1-2879_+",
        "near_end:2479-5000_-",
        "tiny:1-3000_+",
    ]

    windows = read_output_fasta(made.out_fasta)
    assert windows[0]["sequence"] == near_start[0:2879]
    assert windows[1]["sequence"] == reverse_complement_simple(near_end[2478:5000])
    assert windows[2]["sequence"] == tiny

    # window_length in the table is the clipped length, not 4001.
    column_names, rows = read_output_table(made.out_table)
    assert rows[0]["window_length"] == "2879"
    assert rows[1]["window_length"] == "2522"
    assert rows[2]["window_length"] == "3000"


def test_windows_follow_the_fasta_order_and_then_the_position(tmp_path):
    """Windows come out in FASTA contig order, then left to right along each contig.

    The order is neither alphabetical nor the order of the table. The FASTA lists zeta
    before alpha, but the table starts with alpha. On zeta the table lists the minus
    site first (middle 16079) and then the plus sites from right to left (middles 10079
    and 3079). The windows must still come out from left to right, whatever the strand:
    1079-5079, 8079-12079 and 14079-18079.
    """
    rng = random.Random(5)
    records = [("zeta", random_dna(rng, 20000)), ("alpha", random_dna(rng, 20000))]
    table = make_table(
        hit_line("alpha", 20000, "+", 10001, 10158, 1e-10),
        hit_line("zeta", 20000, "-", 16001, 16158, 1e-10),
        hit_line("zeta", 20000, "+", 10001, 10158, 1e-10),
        hit_line("zeta", 20000, "+", 3001, 3158, 1e-10),
    )
    made = make_inputs(tmp_path, records, table)
    run_script(made)

    assert headers_of(made.out_fasta) == [
        "zeta:1079-5079_+",
        "zeta:8079-12079_+",
        "zeta:14079-18079_-",
        "alpha:8079-12079_+",
    ]


# ======================================================================
# NO HITS
# ======================================================================

def check_empty_outputs(made, capsys):
    """Shared checks for a run that found no 5.8S.

    Expected: an empty FASTA, a table with only its header line, and one line on stderr.
    """
    assert made.out_fasta.exists()
    assert made.out_fasta.read_text() == ""
    assert made.out_table.read_text() == TABLE_HEADER_LINE + "\n"

    captured = capsys.readouterr()
    assert "no 5.8S found" in captured.err
    assert captured.err.count("\n") == 1   # a single line


def test_table_without_any_hit_gives_empty_outputs_and_no_error(tmp_path, capsys):
    """A genome without a recognisable 5.8S must not make the workflow fail."""
    rng = random.Random(6)
    made = make_inputs(tmp_path, [("m", random_dna(rng, 1000))], make_table())

    run_script(made)   # must return normally: no SystemExit, so the exit code is 0

    check_empty_outputs(made, capsys)


def test_table_with_only_weak_hits_gives_empty_outputs(inputs, capsys):
    """Every hit is above the E-value limit, so nothing is left: same empty outputs."""
    run_script(inputs, "--max-evalue", "1e-20")

    check_empty_outputs(inputs, capsys)


# ======================================================================
# BAD INPUT GIVES A CLEAR ERROR
# ======================================================================

def test_hit_on_a_contig_missing_from_the_fasta_is_a_clear_error(inputs, capsys):
    """The table names a contig that the FASTA does not contain: stop and say which one."""
    table = make_table(hit_line("contig_missing", 12000, "+", 5001, 5158, 1e-12))
    inputs.tblout.write_text(table)

    with pytest.raises(SystemExit) as stop:
        run_script(inputs)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("its_windows.py: error:")
    assert "contig_missing" in message
    assert "not in" in message
    assert str(inputs.fasta) in message
    assert message.count("\n") == 1   # one line, no Python traceback

    # Nothing is written when the input is inconsistent.
    assert not inputs.out_fasta.exists()
    assert not inputs.out_table.exists()


def test_contig_length_that_differs_from_the_table_is_a_clear_error(inputs, capsys):
    """The table says contig_1 has 29,999 bases but the FASTA has 30,000: wrong genome file."""
    table = make_table(hit_line(CONTIG, 29999, "+", 10001, 10158, 1e-12))
    inputs.tblout.write_text(table)

    with pytest.raises(SystemExit) as stop:
        run_script(inputs)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert "contig_1" in message
    assert "29999" in message
    assert "30000" in message
    assert not inputs.out_fasta.exists()


def test_duplicate_contig_name_in_the_fasta_is_a_clear_error(tmp_path, capsys):
    """Two contigs with the same name would make the table's coordinates ambiguous."""
    rng = random.Random(8)
    records = [("dup", random_dna(rng, 6000)), ("dup", random_dna(rng, 6000))]
    table = make_table(hit_line("dup", 6000, "+", 3001, 3158, 1e-10))
    made = make_inputs(tmp_path, records, table)

    with pytest.raises(SystemExit) as stop:
        run_script(made)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert "dup" in message
    assert "more than once" in message


@pytest.mark.parametrize(
    "bad_line, expected_words",
    [
        # Too few columns.
        ("contig_1 - some_profile - 1 45 10002 10045\n", ["line 3", "columns"]),
        # The envelope start is not a number.
        ("contig_1 - some_profile - 1 45 10002 10045 abc 10045 30000 + 1e-10 40.0 1.0 -\n",
         ["line 3", "numbers"]),
        # The strand is neither + nor -.
        ("contig_1 - some_profile - 1 45 10002 10045 10001 10045 30000 . 1e-10 40.0 1.0 -\n",
         ["line 3", "strand"]),
    ],
    ids=["too_few_columns", "text_instead_of_number", "strand_is_a_dot"],
)
def test_malformed_table_line_is_a_clear_error_with_its_line_number(
    inputs, capsys, bad_line, expected_words
):
    """A damaged table line stops the script with a message that names the line."""
    # The two comment lines of NHMMER_COMMENTS_TOP are lines 1 and 2, so the bad line is line 3.
    inputs.tblout.write_text(NHMMER_COMMENTS_TOP + bad_line + NHMMER_COMMENTS_BOTTOM)

    with pytest.raises(SystemExit) as stop:
        run_script(inputs)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    for word in expected_words:
        assert word in message


def test_missing_table_file_is_a_clear_error(inputs, capsys):
    """A table file that does not exist gives a one-line message, not a traceback."""
    inputs.tblout.unlink()

    with pytest.raises(SystemExit) as stop:
        run_script(inputs)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("its_windows.py: error:")
    assert "58S.tblout" in message


# ======================================================================
# THE SCRIPT AS SNAKEMAKE WILL CALL IT
# ======================================================================

def run_from_command_line(made):
    """Run the script file in a separate Python process, as a Snakemake rule would."""
    command = [
        sys.executable, SCRIPT_PATH,
        "--tblout", str(made.tblout),
        "--fasta", str(made.fasta),
        "--out-fasta", str(made.out_fasta),
        "--out-table", str(made.out_table),
    ]
    return subprocess.run(command, capture_output=True, text=True)


def test_command_line_run_succeeds(inputs):
    """Exit code 0, two windows written, summary on stderr."""
    finished = run_from_command_line(inputs)

    assert finished.returncode == 0
    assert finished.stdout == ""
    assert "wrote 2 window(s)" in finished.stderr
    assert headers_of(inputs.out_fasta) == ["contig_1:8079-12079_+", "contig_1:18079-22079_-"]


def test_command_line_run_without_any_5_8s_still_exits_with_code_0(tmp_path):
    """The workflow must not fail on a genome without a recognisable 5.8S."""
    made = make_inputs(tmp_path, [("m", "ACGT" * 100)], make_table())
    finished = run_from_command_line(made)

    assert finished.returncode == 0
    assert "no 5.8S found" in finished.stderr
    assert made.out_fasta.read_text() == ""


def test_command_line_run_with_a_missing_contig_exits_with_code_1(inputs):
    """An inconsistent input fails the step (non-zero exit) with a readable message."""
    table = make_table(hit_line("contig_missing", 12000, "+", 5001, 5158, 1e-12))
    inputs.tblout.write_text(table)

    finished = run_from_command_line(inputs)

    assert finished.returncode == 1
    assert "its_windows.py: error:" in finished.stderr
    assert "contig_missing" in finished.stderr
    assert "Traceback" not in finished.stderr
