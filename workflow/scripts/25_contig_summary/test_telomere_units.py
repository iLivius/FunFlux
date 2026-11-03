"""Tests for telomere_units.py, the script that chooses the telomere repeat unit from the data.

How to run (from the repository root):

    python -m pytest workflow/scripts/25_contig_summary -q

Only the Python standard library and pytest are needed for most tests. The tests that
run tidk itself are skipped when tidk is not on PATH (the tidk conda environment,
workflow/envs/tidk.yaml, holds it); to run them, put that environment's bin folder
first on PATH.

THE MADE-UP DATA
----------------
A small assembly of six contigs, built so that every part of the method shows:
  * chrom_a        12 kb of random sequence with a telomere array at each end: 15
                   copies of CCCTAA at the start (after 60 bases) and 12 copies of
                   TTAGGG at the end (before 30 bases): a whole chromosome.
  * chrom_b        8 kb with 10 copies of CCCTAA at the start and 10 of TTAGGG at the
                   end: a second whole chromosome.
  * short_dup      14 copies of CCCTAA, then 1.5 kb: a short contig whose one array is
                   seen from both ends (present at the start; from the end the same
                   array with the opposite status), the case of the same-array rule.
  * sat_one_way    8 copies of ATTCTT at the start, then 4 kb: a satellite that happens
                   to point the way a telomere points there.
  * sat_other_way  8 copies of AAGAAT at the start, then 4 kb: the same satellite the
                   other way round. Together the two give the unit AAGAAT one telomere
                   end and one opposite end: score 0, not eligible.
  * plain          4 kb, then 40 copies of ACAT: a 4-bp microsatellite at a contig end,
                   which the length filter drops.
The window is set to 3,000 bases, so that the two chromosomes are cut into start and
end windows while the short contigs go whole. The telomere database is a made-up table
with the real header: AACCCT with strong or weak support in two fungi and one animal,
AAACCCT in a plant (and in a fungus with support unchecked, which does not count),
AACCT in a species without a kingdom, and AAAT, which the filters drop.
"""

import csv
import os
import random
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

# SimpleNamespace (used by make_inputs below) is a plain container whose fields are
# read with a dot, like made.contigs.

# The script sits in the same folder as this test file. Putting that folder on the
# list of places where Python looks for modules makes "import telomere_units" work
# wherever pytest is started from.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TEST_DIR)

import telomere_units  # noqa: E402  (this import must come after the line above)

SCRIPT_PATH = os.path.join(TEST_DIR, "telomere_units.py")

# The tests that run tidk itself need the program on PATH.
needs_tidk = pytest.mark.skipif(shutil.which("tidk") is None, reason="tidk is not on PATH")

TELOMERE = "TTAGGG"
TELOMERE_RC = "CCCTAA"
SATELLITE = "AAGAAT"           # canonical form AAGAAT; G-rich strand
SATELLITE_RC = "ATTCTT"        # the C-rich strand, the one a start expects
SAMPLE = "SAMPLE1"
WINDOW = 3000


# ======================================================================
# THE MADE-UP DATA
# ======================================================================

def random_dna(length, seed):
    """Return `length` random bases; the seed makes the sequence the same on every run."""
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


# The six contigs (name -> sequence). Built once; every test reads them.
CONTIGS = {
    "chrom_a": random_dna(60, 1) + TELOMERE_RC * 15 + random_dna(12000, 2) + TELOMERE * 12 + random_dna(30, 3),
    "chrom_b": TELOMERE_RC * 10 + random_dna(8000, 5) + TELOMERE * 10,
    "short_dup": TELOMERE_RC * 14 + random_dna(1500, 6),
    "sat_one_way": SATELLITE_RC * 8 + random_dna(4000, 7),
    "sat_other_way": SATELLITE * 8 + random_dna(4000, 8),
    "plain": random_dna(4000, 9) + "ACAT" * 40,
}

# The real header of telo_database_v2.tsv, and the made-up rows: (species, kingdom,
# unit, support).
DATABASE_HEADER = [
    "taxon_id", "species", "genus", "family", "order", "class", "phylum", "kingdom",
    "superkingdom", "telomere_unit", "support", "ends_with_telomere", "ends_checked",
    "t2t_sequences", "source", "assembly", "telo_seq_2",
]
DATABASE_ROWS = [
    ("Rhizosphaera kalkhoffii", "Fungi", "AACCCT", "strong"),
    ("Ramularia vizellae", "Fungi", "AACCCT", "weak"),
    ("Osmia bicornis", "Metazoa", "AACCCT", "strong"),
    ("Osmia bicornis", "Metazoa", "AACCCT", "strong"),      # a second assembly of one species
    ("Arabidopsis thaliana", "Viridiplantae", "AAACCCT", "strong"),
    ("Some yeast", "Fungi", "AAACCCT", "unchecked"),         # does not count
    ("Unplaced species", "", "AACCT", "weak"),               # no kingdom
    ("Short unit species", "Metazoa", "AAAT", "strong"),     # dropped by the filters
    ("No unit species", "Metazoa", "", "none"),
]


def write_fasta(path, records):
    """Write (name, sequence) pairs as FASTA, wrapped at 80 bases."""
    with open(path, "w") as fasta:
        for name, sequence in records:
            fasta.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                fasta.write(sequence[start:start + 80] + "\n")


def write_database(path, rows=DATABASE_ROWS):
    """Write a telomere database in the real layout: the 17-column header, one row per species."""
    with open(path, "w") as table:
        table.write("\t".join(DATABASE_HEADER) + "\n")
        for number, (species, kingdom, unit, support) in enumerate(rows, start=1):
            cells = [str(number), species, "", "", "", "", "", kingdom, "Eukaryota", unit, support,
                     "", "", "", "tidk 0.3.0 assembly", "", ""]
            table.write("\t".join(cells) + "\n")


def make_inputs(folder):
    """Write the made-up assembly and database into `folder`; return every path."""
    paths = SimpleNamespace()
    paths.contigs = os.path.join(folder, "assembly.fasta")
    write_fasta(paths.contigs, list(CONTIGS.items()))
    paths.database = os.path.join(folder, "telo_database_v2.tsv")
    write_database(paths.database)
    paths.out_candidates = os.path.join(folder, "telomere_candidates.tsv")
    paths.out_unit = os.path.join(folder, "telomere_unit.txt")
    paths.out_ends = os.path.join(folder, "telomere_ends.tsv")
    return paths


def full_argv(paths, **changes):
    """The complete command line for the made-up assembly; `changes` replaces options."""
    options = {
        "--contigs": paths.contigs,
        "--telomere-db": paths.database,
        "--telomere-motif": "auto",
        "--telomere-window": WINDOW,
        "--min-telomere-copies": 5,
        "--sample": SAMPLE,
        "--out-candidates": paths.out_candidates,
        "--out-unit": paths.out_unit,
        "--out-ends": paths.out_ends,
    }
    options.update(changes)
    argv = []
    for option, value in options.items():
        if value is not None:
            argv.extend([option, str(value)])
    return argv


def read_candidates(path):
    """Read the candidates table back as {unit: row dictionary}, in file order."""
    with open(path, newline="") as table:
        rows = list(csv.DictReader(table, delimiter="\t"))
    return {row["unit"]: row for row in rows}


def read_ends(path):
    """Read a per-end table back as {(id, side): row dictionary}."""
    with open(path, newline="") as table:
        rows = list(csv.DictReader(table, delimiter="\t"))
    return {(row["id"], row["side"]): row for row in rows}


def read_unit(path):
    """Read telomere_unit.txt back: its one line."""
    with open(path) as unit_file:
        return unit_file.read().strip()


def candidate(unit, telomere_ends, opposite_ends, known=True, t2t=0, found=True, array=True, configured=False):
    """A scored candidate as run() builds it, for the tests of the ranking and the files."""
    return {
        "unit": unit,
        "found_at_ends": found,
        "configured": configured,
        "in_literature": "TTAGGG: some fungi (Someone 2020)" if known else "-",
        "in_database": "-",
        "literature_entries": [("TTAGGG", "some fungi", "Someone 2020")] if known else [],
        "database_species": {},
        "known": known,
        "formed_array": array,
        "telomere_ends": telomere_ends,
        "opposite_ends": opposite_ends,
        "score": telomere_ends - opposite_ends,
        "t2t_contigs": t2t,
        "copies_at_ends": 10,
        "copies_genome_wide": 20,
        "chosen": False,
    }


# The per-end table tidk ends 0.3.0 wrote for the unit TTAGGG on the ONT assembly of
# the public benchmark genome (window 25,000, minimum length 30, maximum offset
# 25,000): 33 contigs, 34 ends with a telomere, 15 contigs with one at both ends.
EXAMPLE_ENDS_ONT = """id	sequence_length	side	status	telomere_start	telomere_end	telomere_length	offset	strand	window_limited
ptg000001l	2019885	start	present	78	204	126	78	C-rich	false
ptg000001l	2019885	end	present	2019699	2019813	114	72	G-rich	false
ptg000002l	2354010	start	present	78	192	114	78	C-rich	false
ptg000002l	2354010	end	present	2352768	2352876	108	1134	G-rich	false
ptg000003l	49742	start	absent	NA	NA	NA	NA	NA	NA
ptg000003l	49742	end	absent	NA	NA	NA	NA	NA	NA
ptg000004l	282192	start	present	19746	19848	102	19746	C-rich	false
ptg000004l	282192	end	present	266856	266958	102	15234	G-rich	false
ptg000005l	4109401	start	present	78	186	108	78	C-rich	false
ptg000005l	4109401	end	present	4109215	4109317	102	84	G-rich	false
ptg000006l	2686745	start	present	72	192	120	72	C-rich	false
ptg000006l	2686745	end	present	2686571	2686679	108	66	G-rich	false
ptg000007l	3044302	start	present	2814	2922	108	2814	C-rich	false
ptg000007l	3044302	end	absent	NA	NA	NA	NA	NA	NA
ptg000008l	2183122	start	present	48	162	114	48	C-rich	false
ptg000008l	2183122	end	present	2182936	2183044	108	78	G-rich	false
ptg000009l	3230148	start	present	3660	3768	108	3660	C-rich	false
ptg000009l	3230148	end	present	3229956	3230070	114	78	G-rich	false
ptg000010l	2726025	start	present	108	222	114	108	C-rich	false
ptg000010l	2726025	end	present	2725833	2725947	114	78	G-rich	false
ptg000011l	58954	start	absent	NA	NA	NA	NA	NA	NA
ptg000011l	58954	end	absent	NA	NA	NA	NA	NA	NA
ptg000012l	158138	start	absent	NA	NA	NA	NA	NA	NA
ptg000012l	158138	end	absent	NA	NA	NA	NA	NA	NA
ptg000013l	3334582	start	present	78	198	120	78	C-rich	false
ptg000013l	3334582	end	present	3334402	3334510	108	72	G-rich	false
ptg000014l	2512653	start	present	84	192	108	84	C-rich	false
ptg000014l	2512653	end	absent	NA	NA	NA	NA	NA	NA
ptg000015l	61330	start	absent	NA	NA	NA	NA	NA	NA
ptg000015l	61330	end	absent	NA	NA	NA	NA	NA	NA
ptg000016l	2576789	start	present	84	186	102	84	C-rich	false
ptg000016l	2576789	end	present	2576603	2576717	114	72	G-rich	false
ptg000017l	2399821	start	present	84	198	114	84	C-rich	false
ptg000017l	2399821	end	present	2387413	2387521	108	12300	G-rich	false
ptg000018l	2027725	start	present	78	186	108	78	C-rich	false
ptg000018l	2027725	end	present	2027539	2027653	114	72	G-rich	false
ptg000019l	78179	start	absent	NA	NA	NA	NA	NA	NA
ptg000019l	78179	end	absent	NA	NA	NA	NA	NA	NA
ptg000020l	2630675	start	present	78	192	114	78	C-rich	false
ptg000020l	2630675	end	present	2630483	2630591	108	84	G-rich	false
ptg000021l	2257611	start	present	84	198	114	84	C-rich	false
ptg000021l	2257611	end	present	2257431	2257545	114	66	G-rich	false
ptg000022l	208977	start	present	78	186	108	78	C-rich	false
ptg000022l	208977	end	present	208797	208899	102	78	G-rich	false
ptg000023l	2432023	start	present	78	186	108	78	C-rich	false
ptg000023l	2432023	end	absent	NA	NA	NA	NA	NA	NA
ptg000024l	67439	start	absent	NA	NA	NA	NA	NA	NA
ptg000024l	67439	end	absent	NA	NA	NA	NA	NA	NA
ptg000025l	38975	start	absent	NA	NA	NA	NA	NA	NA
ptg000025l	38975	end	absent	NA	NA	NA	NA	NA	NA
ptg000026l	85908	start	absent	NA	NA	NA	NA	NA	NA
ptg000026l	85908	end	absent	NA	NA	NA	NA	NA	NA
ptg000027l	35521	start	absent	NA	NA	NA	NA	NA	NA
ptg000027l	35521	end	absent	NA	NA	NA	NA	NA	NA
ptg000028c	82351	start	absent	NA	NA	NA	NA	NA	NA
ptg000028c	82351	end	absent	NA	NA	NA	NA	NA	NA
ptg000029l	42998	start	absent	NA	NA	NA	NA	NA	NA
ptg000029l	42998	end	absent	NA	NA	NA	NA	NA	NA
ptg000030l	175720	start	present	78	186	108	78	C-rich	false
ptg000030l	175720	end	absent	NA	NA	NA	NA	NA	NA
ptg000031l	65582	start	absent	NA	NA	NA	NA	NA	NA
ptg000031l	65582	end	absent	NA	NA	NA	NA	NA	NA
ptg000032l	36585	start	absent	NA	NA	NA	NA	NA	NA
ptg000032l	36585	end	absent	NA	NA	NA	NA	NA	NA
ptg000033l	72367	start	absent	NA	NA	NA	NA	NA	NA
ptg000033l	72367	end	absent	NA	NA	NA	NA	NA	NA
"""


def ends_row(contig, side, status, start="NA", end="NA", offset="NA", strand="NA", length=5000):
    """One row of a per-end table, as read_ends_table() returns it (every value text)."""
    telomere_length = "NA" if start == "NA" else str(int(end) - int(start))
    return {
        "id": contig, "sequence_length": str(length), "side": side, "status": status,
        "telomere_start": str(start), "telomere_end": str(end), "telomere_length": telomere_length,
        "offset": str(offset), "strand": strand, "window_limited": "false",
    }


# ======================================================================
# THE BUILDING BLOCKS: CANONICAL FORM, THE FILTERS, THE SPELLING NOTE
# ======================================================================

@pytest.mark.parametrize("unit", ["TTAGGG", "GGGTTA", "CCCTAA", "TAACCC", "AGGGTT"])
def test_canonical_form_of_the_telomere_motif_is_aaccct(unit):
    """Every rotation of TTAGGG, on either strand, has the canonical form AACCCT."""
    assert telomere_units.canonical_unit(unit) == "AACCCT"


def test_canonical_form_is_the_same_for_rotations_and_the_other_strand():
    """A 12-bp unit, its rotations and its reverse complement share one canonical form."""
    unit = "ACGTTGCAGGTC"
    expected = telomere_units.canonical_unit(unit)
    for rotation in telomere_units.rotations(unit):
        assert telomere_units.canonical_unit(rotation) == expected
        assert telomere_units.canonical_unit(telomere_units.reverse_complement(rotation)) == expected
    all_forms = telomere_units.rotations(unit) + telomere_units.rotations(telomere_units.reverse_complement(unit))
    assert expected == min(all_forms)


@pytest.mark.parametrize(
    "unit, expected",
    [
        ("AACCCT", True),                      # TTAGGG
        ("AAACCCT", True),                     # TTTAGGG
        ("AACTTCTTGGTGTACGGATGTCT", True),     # the 23-bp Candida albicans unit
        ("ACAT", False),                       # 4 bp: too short
        ("A" * 31, False),                     # too long
        ("ACACAC", False),                     # a repeat of AC
        ("AACCCTAACCCT", False),               # a repeat of AACCCT
        ("AAAAAA", False),                     # a homopolymer
        ("AAAT", False),                       # too short, and no G-rich strand
        ("AAAAT", False),                      # no G and no C: the strands are the same
        ("AGAGCC", False),                     # two G and two C: the strands are the same
        ("AACGT", False),                      # one G, one C
    ],
)
def test_the_three_candidate_filters(unit, expected):
    """Length 5 to 30, not a repeat of a shorter unit, and two strands that differ."""
    assert telomere_units.passes_candidate_filters(unit) is expected


def test_strands_differ_compares_the_counts_of_g_and_c():
    """AACCCT has three C and no G; AGAGCC has two of each."""
    assert telomere_units.strands_differ("AACCCT") is True
    assert telomere_units.strands_differ("AGAGCC") is False
    assert telomere_units.is_repeat_of_shorter_unit("ACGACG") is True
    assert telomere_units.is_repeat_of_shorter_unit("ACGTTGCAGGTC") is False


def test_every_literature_unit_passes_the_filters_and_has_its_own_canonical_form():
    """The eleven published units are plausible candidates, and no two are the same unit."""
    by_unit = telomere_units.literature_by_unit()
    assert len(by_unit) == len(telomere_units.LITERATURE) == 11
    for spelling, fungi, references in telomere_units.LITERATURE:
        unit = telomere_units.canonical_unit(spelling)
        assert telomere_units.passes_candidate_filters(unit), spelling
        assert by_unit[unit] == [(spelling, fungi, references)]
    assert by_unit["AACCCT"][0][0] == "TTAGGG"
    assert "Botrytis cinerea" in by_unit["AACCCT"][0][1]
    assert telomere_units.in_literature_text(by_unit["AACCCT"]) == (
        "TTAGGG: filamentous Ascomycota, most Basidiomycota, Mucoromycota; Botrytis cinerea "
        "(Lue 2021; Sun et al. 2023; Levis et al. 1997)"
    )


def test_spelling_note_says_how_the_canonical_form_relates_to_the_published_one():
    """TTAGGG is AACCCT on the other strand; TTAGGGTCAACA is its canonical form read from another position."""
    assert telomere_units.spelling_note("TTAGGG", "AACCCT") == "TTAGGG read on the other strand"
    assert telomere_units.spelling_note("TTAGGGTCAACA", "AACATTAGGGTC") == "TTAGGGTCAACA read from another position"
    assert telomere_units.spelling_note("AACCCT", "AACCCT") == ""


# ======================================================================
# THE CONTIG ENDS, AND READING WHAT TIDK WRITES
# ======================================================================

def test_contig_ends_file_cuts_the_long_contigs_and_keeps_the_short_ones_whole(tmp_path):
    """A contig longer than 2W gives two windows; one of at most 2W goes in whole."""
    contigs = {"long": random_dna(7000, 30), "exact": random_dna(6000, 31), "short": random_dna(100, 32)}
    path = tmp_path / "ends.fasta"
    records = telomere_units.write_contig_ends(contigs, 3000, path)
    assert records == 4
    written = telomere_units.read_fasta(path)
    assert list(written) == ["long__start", "long__end", "exact", "short"]
    assert written["long__start"] == contigs["long"][:3000]
    assert written["long__end"] == contigs["long"][-3000:]
    assert written["exact"] == contigs["exact"]
    assert written["short"] == contigs["short"]


def test_explore_table_is_read_from_tidk_output_with_its_message_lines():
    """The rows after the header are read; messages before it and units with N are skipped."""
    text = (
        "Warning! No clades found in the database. Run 'tidk build' to fetch the latest data.\n"
        "[+]\tExploring genome for potential telomeric repeats between lengths 5 and 30.\n"
        "canonical_repeat_unit\tcopies\tcopies_as_unit\tcopies_as_revcomp\n"
        "ACAT\t1020\t447\t573\n"
        "AACCCT\t661\t374\t287\n"
        "AANCCT\t9\t9\t0\n"
        "\n"
    )
    copies = telomere_units.parse_explore_table(text)
    assert list(copies.items()) == [("ACAT", 1020), ("AACCCT", 661)]
    assert telomere_units.parse_explore_table("canonical_repeat_unit\tcopies\tcopies_as_unit\tcopies_as_revcomp\n") == {}

    with pytest.raises(ValueError, match="another tidk version"):
        telomere_units.parse_explore_table("canonical_repeat_unit\tcopies\n")
    with pytest.raises(ValueError, match="no table"):
        telomere_units.parse_explore_table("[+]\tnothing\n")
    with pytest.raises(ValueError, match="not a whole number"):
        telomere_units.parse_explore_table(
            "canonical_repeat_unit\tcopies\tcopies_as_unit\tcopies_as_revcomp\nAACCCT\tmany\t1\t1\n"
        )


def test_ends_table_of_the_benchmark_example_counts_34_ends_and_15_t2t_contigs(tmp_path):
    """The example table gives telomere_ends 34, opposite_ends 0, score 34, t2t 15, and the same-array rule changes nothing."""
    path = tmp_path / "example.ends.tsv"
    path.write_text(EXAMPLE_ENDS_ONT)
    rows = telomere_units.read_ends_table(path)
    assert len(rows) == 66
    assert list(rows[0]) == telomere_units.ENDS_COLUMNS
    assert rows[0]["id"] == "ptg000001l" and rows[0]["offset"] == "78"
    assert rows[4]["status"] == "absent" and rows[4]["offset"] == "NA"

    before = [row["status"] for row in rows]
    telomere_units.mark_arrays_seen_from_both_ends(rows)
    assert [row["status"] for row in rows] == before
    assert telomere_units.count_ends(rows) == {
        "telomere_ends": 34, "opposite_ends": 0, "score": 34, "t2t_contigs": 15,
    }
    assert telomere_units.formed_array(rows) is True


def test_ends_table_with_another_header_is_refused(tmp_path):
    """A table whose columns are not tidk 0.3.0's stops the script rather than being misread."""
    path = tmp_path / "odd.ends.tsv"
    path.write_text("id\tside\tstatus\nx\tstart\tabsent\n")
    with pytest.raises(ValueError, match="another tidk version"):
        telomere_units.read_ends_table(path)
    path.write_text(EXAMPLE_ENDS_ONT.splitlines()[0] + "\nptg000001l\t100\tstart\n")
    with pytest.raises(ValueError, match="3 columns instead of 10"):
        telomere_units.read_ends_table(path)


# ======================================================================
# THE SAME-ARRAY RULE AND THE COUNTS
# ======================================================================

def test_same_array_rule_marks_the_far_view_of_one_array_as_other_end():
    """One array seen from both ends of a short contig keeps the near row; the far row becomes other_end."""
    rows = [
        ends_row("short", "start", "present", 0, 84, 0, "C-rich", length=1584),
        ends_row("short", "end", "wrong_strand", 2, 86, 1498, "C-rich", length=1584),
        ends_row("chrom", "start", "present", 60, 150, 60, "C-rich", length=12000),
        ends_row("chrom", "end", "present", 11900, 11970, 30, "G-rich", length=12000),
        ends_row("empty", "start", "absent"),
        ends_row("empty", "end", "absent"),
        ends_row("one_side", "start", "wrong_strand", 100, 150, 100, "G-rich"),
        ends_row("one_side", "end", "absent"),
    ]
    telomere_units.mark_arrays_seen_from_both_ends(rows)
    statuses = {(row["id"], row["side"]): row["status"] for row in rows}
    assert statuses[("short", "start")] == "present"
    assert statuses[("short", "end")] == "other_end"       # the same bases, the larger offset
    assert statuses[("chrom", "start")] == "present"        # two different arrays: unchanged
    assert statuses[("chrom", "end")] == "present"
    assert statuses[("empty", "start")] == "absent"
    assert statuses[("one_side", "start")] == "wrong_strand"
    assert telomere_units.count_ends(rows) == {
        "telomere_ends": 3, "opposite_ends": 1, "score": 2, "t2t_contigs": 1,
    }


def test_same_array_rule_keeps_the_row_with_the_smaller_offset_whichever_side_it_is():
    """An array near the END of a short contig keeps the end row; equal offsets mark the end row."""
    rows = [
        ends_row("near_end", "start", "wrong_strand", 1400, 1500, 1400, "G-rich", length=1500),
        ends_row("near_end", "end", "present", 1398, 1500, 0, "G-rich", length=1500),
        ends_row("all_array", "start", "present", 0, 120, 0, "C-rich", length=120),
        ends_row("all_array", "end", "wrong_strand", 0, 120, 0, "C-rich", length=120),
    ]
    telomere_units.mark_arrays_seen_from_both_ends(rows)
    statuses = {(row["id"], row["side"]): row["status"] for row in rows}
    assert statuses[("near_end", "start")] == "other_end"
    assert statuses[("near_end", "end")] == "present"
    assert statuses[("all_array", "start")] == "present"
    assert statuses[("all_array", "end")] == "other_end"


def test_a_row_with_an_array_but_no_number_is_a_clear_error():
    """NA in a coordinate of a row that reports an array means another tidk: stop."""
    rows = [
        ends_row("odd", "start", "present", "NA", "NA", "NA", "C-rich"),
        ends_row("odd", "end", "present", 0, 50, 0, "G-rich"),
    ]
    with pytest.raises(ValueError, match="not a number"):
        telomere_units.mark_arrays_seen_from_both_ends(rows)


# ======================================================================
# ELIGIBILITY, RANKING AND THE CHOICE
# ======================================================================

@pytest.mark.parametrize(
    "telomere_ends, opposite_ends, expected",
    [(2, 0, True), (1, 0, False), (0, 0, False), (8, 5, False), (5, 2, True), (4, 2, False), (34, 0, True)],
)
def test_eligibility_needs_two_ends_and_more_than_twice_the_opposite_ends(telomere_ends, opposite_ends, expected):
    """At least 2 telomere ends, and telomere ends > 2 x opposite ends."""
    assert telomere_units.is_eligible(candidate("AACCCT", telomere_ends, opposite_ends)) is expected


def test_ranking_by_score_then_known_then_length_then_alphabet():
    """Score first; between equal scores a known unit, then a shorter unit, then the alphabet."""
    candidates = [
        candidate("AAGGAG", 3, 0, known=False),
        candidate("AACCCT", 5, 0, known=True),
        candidate("AAAGAG", 3, 0, known=True),
        candidate("AAACCCT", 3, 0, known=True),
        candidate("AAAAGAG", 3, 0, known=False),
        candidate("ACAGGG", 10, 7, known=False),      # score 3 too
    ]
    ranked = telomere_units.rank_candidates(candidates)
    assert [c["unit"] for c in ranked] == ["AACCCT", "AAAGAG", "AAACCCT", "AAGGAG", "ACAGGG", "AAAAGAG"]


def test_first_eligible_skips_a_better_scored_candidate_that_points_both_ways():
    """20 ends against 12 opposite scores 8 but is not eligible; 5 against 0 wins."""
    ranked = telomere_units.rank_candidates([
        candidate("AGATG", 20, 12), candidate("AACCCT", 5, 0), candidate("AAGAAT", 1, 1),
    ])
    assert ranked[0]["unit"] == "AGATG"
    assert telomere_units.first_eligible(ranked)["unit"] == "AACCCT"
    assert telomere_units.first_eligible([candidate("AAGAAT", 1, 1)]) is None


def test_collect_candidates_merges_the_sources_and_filters_all_but_the_configured_unit():
    """Units from the ends, the database and the literature are merged; ACAT and AAAT drop; a configured AAAT stays."""
    at_ends = {"ACAT": 100, "AACCCT": 50, "AAGAAT": 8}
    database = {"AACCCT": {"Fungi": {"a", "b"}}, "AAAT": {"Metazoa": {"c"}}, "AAACCCT": {"Fungi": {"d"}}}
    literature = telomere_units.literature_by_unit()

    candidates, dropped = telomere_units.collect_candidates(at_ends, database, literature, None)
    assert "ACAT" not in candidates and "AAAT" not in candidates
    assert dropped == 2
    assert list(candidates) == sorted(candidates)                       # alphabetical: the same order every run
    assert candidates["AACCCT"]["found_at_ends"] is True
    assert candidates["AACCCT"]["in_database"] == "Fungi 2"
    assert candidates["AACCCT"]["in_literature"].startswith("TTAGGG: ")
    assert candidates["AACCCT"]["known"] is True
    assert candidates["AAGAAT"]["known"] is False
    assert candidates["AAGAAT"]["in_database"] == "-" and candidates["AAGAAT"]["in_literature"] == "-"
    assert candidates["AAACCCT"]["found_at_ends"] is False
    # From the ends: AACCCT and AAGAAT (ACAT dropped). From the database: AAACCCT (AACCCT
    # is the same unit, AAAT dropped). From the literature: the nine units not yet listed.
    assert len(candidates) == 2 + 1 + 9

    candidates, dropped = telomere_units.collect_candidates(at_ends, database, literature, "AAAT")
    assert "AAAT" in candidates and candidates["AAAT"]["configured"] is True
    assert dropped == 1


def test_summary_line_names_the_unit_the_counts_and_what_is_known():
    """The summary gives the canonical unit, the published spelling, the counts and the sources."""
    chosen = candidate("AACCCT", 34, 0, t2t=15)
    chosen["database_species"] = {"Fungi": {"a", "b"}, "Metazoa": {"c"}}
    chosen["in_database"] = "Fungi 2; Metazoa 1"
    line = telomere_units.summary_line(SAMPLE, chosen, chosen, None)
    assert line == (
        "[FunFlux] telomere_units SAMPLE1: telomere repeat unit AACCCT (TTAGGG read on the "
        "other strand): 34 contig ends, 15 contigs with a telomere at both ends; known from "
        "the literature (some fungi) and from 2 fungal species in the database"
    )

    # Only other kingdoms in the database, and not in the literature.
    other = candidate("AAGGAG", 4, 0, known=False)
    other["database_species"] = {"Metazoa": {"c"}}
    other["in_database"] = "Metazoa 1"
    assert telomere_units.summary_line(SAMPLE, other, other, None).endswith(
        "4 contig ends, 0 contigs with a telomere at both ends; known from the database (Metazoa 1), from no fungus"
    )
    unknown = candidate("AAGGAG", 4, 0, known=False)
    assert telomere_units.summary_line(SAMPLE, unknown, unknown, None).endswith(
        "not known from the literature or the database"
    )

    # A configured unit: the user's spelling and "set by telomere_motif".
    configured = candidate("AAGGAG", 0, 0, known=False, configured=True)
    assert "AAGGAG (CTCCTT read on the other strand; set by telomere_motif): 0 contig ends" in (
        telomere_units.summary_line(SAMPLE, configured, None, "CTCCTT")
    )

    # Nothing chosen: the best candidate is named, or the absence of any array.
    best = candidate("AAGAAT", 1, 1, known=False)
    assert telomere_units.summary_line(SAMPLE, None, best, None) == (
        "[FunFlux] telomere_units SAMPLE1: no telomere repeat identified (best candidate "
        "AAGAAT: 1 ends in the telomere orientation, 1 opposite)"
    )
    assert telomere_units.summary_line(SAMPLE, None, None, None).endswith(
        "(no candidate unit formed an array at a contig end)"
    )


def test_note_lines_for_an_unknown_unit_and_for_a_configured_unit_the_data_do_not_support():
    """A NOTE names an unknown unit; a NOTE says when the data favour another unit than the configured one."""
    known_choice = candidate("AACCCT", 34, 0)
    assert telomere_units.note_lines(SAMPLE, known_choice, known_choice, None) == []
    assert telomere_units.note_lines(SAMPLE, None, None, None) == []

    unknown_choice = candidate("AAGGAG", 4, 0, known=False)
    notes = telomere_units.note_lines(SAMPLE, unknown_choice, unknown_choice, None)
    assert len(notes) == 1
    assert notes[0].startswith("[FunFlux] NOTE telomere_units SAMPLE1: the telomere repeat unit AAGGAG is not known")

    # The configured unit is the data's first choice: no NOTE.
    configured = candidate("AACCCT", 34, 0, configured=True)
    assert telomere_units.note_lines(SAMPLE, configured, configured, "TTAGGG") == []

    # The data favour another unit than the configured one.
    wrong = candidate("AACATTAGGGTC", 0, 0, configured=True)
    notes = telomere_units.note_lines(SAMPLE, wrong, known_choice, "TTAGGGTCAACA")
    assert len(notes) == 1
    assert notes[0] == (
        "[FunFlux] NOTE telomere_units SAMPLE1: the data favour AACCCT (34 contig ends in the "
        "telomere orientation, 0 opposite; known) over the configured unit TTAGGGTCAACA "
        "(AACATTAGGGTC), which has 0 contig ends in the telomere orientation and 0 opposite; "
        "set telomere_motif to auto or to that unit and re-run"
    )

    # No unit is supported, the configured one included.
    notes = telomere_units.note_lines(SAMPLE, wrong, None, "TTAGGGTCAACA")
    assert len(notes) == 1
    assert "has 0 contig ends in the telomere orientation and 0 opposite, and no other unit is supported either" in notes[0]

    # A configured unit that is unknown AND not the data's choice gets both NOTEs.
    odd = candidate("AAGGAG", 1, 0, known=False, configured=True)
    assert len(telomere_units.note_lines(SAMPLE, odd, known_choice, "AAGGAG")) == 2


# ======================================================================
# THE DATABASE AND THE OUTPUT FILES
# ======================================================================

def test_database_reading_counts_species_per_kingdom_with_strong_or_weak_support(tmp_path):
    """Strong and weak rows count, once per species; unchecked and none do not; no kingdom gets a label."""
    path = tmp_path / "telo_database_v2.tsv"
    write_database(path)
    database = telomere_units.read_telomere_database(path)
    assert set(database) == {"AACCCT", "AAACCCT", "AACCT", "AAAT"}
    assert database["AACCCT"] == {"Fungi": {"Rhizosphaera kalkhoffii", "Ramularia vizellae"}, "Metazoa": {"Osmia bicornis"}}
    assert database["AAACCCT"] == {"Viridiplantae": {"Arabidopsis thaliana"}}
    assert database["AACCT"] == {"(no kingdom)": {"Unplaced species"}}
    assert telomere_units.in_database_text(database["AACCCT"]) == "Fungi 2; Metazoa 1"
    assert telomere_units.in_database_text(database["AACCT"]) == "(no kingdom) 1"
    mixed = {"(no kingdom)": {"x", "y"}, "Viridiplantae": {"p"}, "Fungi": {"f"}}
    assert telomere_units.in_database_text(mixed) == "Fungi 1; Viridiplantae 1; (no kingdom) 2"


def test_database_without_the_needed_columns_is_refused(tmp_path):
    """A table that lacks the support column is not the database: a clear error, naming the column."""
    path = tmp_path / "broken.tsv"
    path.write_text("species\tkingdom\ttelomere_unit\nA b\tFungi\tAACCCT\n")
    with pytest.raises(ValueError, match="missing the column\\(s\\) support"):
        telomere_units.read_telomere_database(path)
    path.write_text("\t".join(DATABASE_HEADER) + "\n")
    with pytest.raises(ValueError, match="no row with strong or weak support"):
        telomere_units.read_telomere_database(path)


def test_candidates_table_has_the_columns_the_order_and_only_the_units_seen(tmp_path):
    """Header exactly as specified; score then telomere_ends descending then unit; known units without an array left out."""
    candidates = [
        candidate("AAGAAT", 1, 1, known=False),
        candidate("AACCCT", 5, 0, known=True, t2t=2),
        candidate("AAACCCT", 0, 0, known=True, found=False, array=False),       # known, never seen: left out
        candidate("ACCCT", 0, 0, known=True, found=False, array=True),          # formed an array: in
        candidate("AAGGAG", 2, 2, known=False),
        candidate("AATCC", 0, 0, known=True, found=False, array=False, configured=True),   # configured: in
    ]
    candidates[1]["chosen"] = True
    path = tmp_path / "telomere_candidates.tsv"
    written = telomere_units.write_candidates_table(path, candidates)
    assert written == 5

    with open(path) as table:
        lines = table.read().splitlines()
    assert lines[0] == "unit\tunit_length\tcopies_at_ends\tcopies_genome_wide\ttelomere_ends\topposite_ends\tscore\tt2t_contigs\tin_literature\tin_database\tchosen"
    rows = read_candidates(path)
    # Score 5 first; then the four at score 0 by telomere ends (2, 1, 0, 0), the last two alphabetically.
    assert list(rows) == ["AACCCT", "AAGGAG", "AAGAAT", "AATCC", "ACCCT"]
    assert rows["AACCCT"] == {
        "unit": "AACCCT", "unit_length": "6", "copies_at_ends": "10", "copies_genome_wide": "20",
        "telomere_ends": "5", "opposite_ends": "0", "score": "5", "t2t_contigs": "2",
        "in_literature": "TTAGGG: some fungi (Someone 2020)", "in_database": "-", "chosen": "yes",
    }
    assert rows["AAGAAT"]["chosen"] == "no" and rows["AAGAAT"]["in_literature"] == "-"
    assert "AAACCCT" not in rows


def test_candidates_table_sorts_equal_scores_by_telomere_ends(tmp_path):
    """Between equal scores the unit with more telomere ends comes first, then the alphabet."""
    candidates = [candidate("AAGAAT", 1, 1, known=False), candidate("AAGGAG", 2, 2, known=False), candidate("AAATC", 2, 2, known=False)]
    path = tmp_path / "t.tsv"
    telomere_units.write_candidates_table(path, candidates)
    assert list(read_candidates(path)) == ["AAATC", "AAGGAG", "AAGAAT"]


def test_unit_and_ends_files_are_written_as_specified(tmp_path):
    """telomere_unit.txt is one line; telomere_ends.tsv has tidk's header and the rows, or the header alone."""
    telomere_units.write_unit_file(tmp_path / "u.txt", "none")
    assert (tmp_path / "u.txt").read_text() == "none\n"
    rows = [ends_row("c", "start", "present", 0, 84, 0, "C-rich"), ends_row("c", "end", "other_end", 2, 86, 1498, "C-rich")]
    telomere_units.write_ends_table(tmp_path / "e.tsv", rows)
    lines = (tmp_path / "e.tsv").read_text().splitlines()
    assert lines[0] == "\t".join(telomere_units.ENDS_COLUMNS)
    assert lines[1] == "c\t5000\tstart\tpresent\t0\t84\t84\t0\tC-rich\tfalse"
    assert lines[2].split("\t")[3] == "other_end"
    telomere_units.write_ends_table(tmp_path / "empty.tsv", [])
    assert (tmp_path / "empty.tsv").read_text() == "\t".join(telomere_units.ENDS_COLUMNS) + "\n"


# ======================================================================
# BAD INPUT GIVES A CLEAR ERROR, AND NO OUTPUT FILE (NO TIDK NEEDED)
# ======================================================================

def assert_stopped_with_message(paths, argv, capsys, expected_text):
    """Run main(argv); it must exit with code 1, print one [FunFlux] error line and write nothing."""
    with pytest.raises(SystemExit) as stop:
        telomere_units.main(argv)
    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("[FunFlux] telomere_units: error:")
    assert expected_text in message
    assert message.count("\n") == 1
    assert not os.path.exists(paths.out_candidates)
    assert not os.path.exists(paths.out_unit)
    assert not os.path.exists(paths.out_ends)


def test_missing_contigs_file_is_a_clear_error(tmp_path, capsys):
    """An input that does not exist gives one readable line, not a traceback."""
    paths = make_inputs(tmp_path)
    missing = str(tmp_path / "nowhere.fasta")
    assert_stopped_with_message(paths, full_argv(paths, **{"--contigs": missing}), capsys, "nowhere.fasta")


def test_malformed_database_is_a_clear_error(tmp_path, capsys):
    """A database without the support column is refused, naming the column."""
    paths = make_inputs(tmp_path)
    with open(paths.database, "w") as table:
        table.write("species\tkingdom\ttelomere_unit\nA b\tFungi\tAACCCT\n")
    assert_stopped_with_message(paths, full_argv(paths), capsys, "missing the column(s) support")


@pytest.mark.parametrize(
    "option, value, expected_text",
    [
        ("--telomere-motif", "TTAGGX", "must be auto"),
        ("--telomere-motif", "TTA", "4 to 30 letters"),
        ("--telomere-motif", "", "must be auto"),
        ("--telomere-window", "100", "from 500 to 100000"),
        ("--telomere-window", "200000", "from 500 to 100000"),
        ("--min-telomere-copies", "1", "at least 2"),
    ],
)
def test_settings_outside_their_limits_are_refused(tmp_path, capsys, option, value, expected_text):
    """The three settings are held to the same limits as the configuration (exit code 1)."""
    paths = make_inputs(tmp_path)
    assert_stopped_with_message(paths, full_argv(paths, **{option: value}), capsys, expected_text)


def test_tidk_that_cannot_be_found_is_a_clear_error(tmp_path, capsys):
    """A wrong --tidk gives one readable line that names the program."""
    paths = make_inputs(tmp_path)
    missing_tidk = str(tmp_path / "no_such_folder" / "tidk")
    assert_stopped_with_message(paths, full_argv(paths, **{"--tidk": missing_tidk}), capsys, "was not found")


def test_duplicate_contig_name_stops_the_script(tmp_path, capsys):
    """Two FASTA records with the same name make the per-end table ambiguous: stop."""
    paths = make_inputs(tmp_path)
    write_fasta(paths.contigs, [("twice", "ACGT" * 300), ("twice", "TTTT" * 300)])
    assert_stopped_with_message(paths, full_argv(paths), capsys, "'twice' appears more than once")


# ======================================================================
# THE WHOLE RUN WITH TIDK ON THE MADE-UP ASSEMBLY
# ======================================================================

@needs_tidk
def test_the_data_choose_aaccct_on_the_made_up_assembly(tmp_path, capsys):
    """Five telomere ends, two T2T contigs, the short contig's far view other_end, the satellite at score 0."""
    paths = make_inputs(tmp_path)
    telomere_units.main(full_argv(paths))

    assert read_unit(paths.out_unit) == "AACCCT"

    rows = read_candidates(paths.out_candidates)
    assert list(rows)[0] == "AACCCT"
    assert rows["AACCCT"]["telomere_ends"] == "5"
    assert rows["AACCCT"]["opposite_ends"] == "0"
    assert rows["AACCCT"]["score"] == "5"
    assert rows["AACCCT"]["t2t_contigs"] == "2"
    assert rows["AACCCT"]["chosen"] == "yes"
    assert rows["AACCCT"]["in_literature"].startswith("TTAGGG: filamentous Ascomycota")
    assert rows["AACCCT"]["in_database"] == "Fungi 2; Metazoa 1"
    # 15 + 12 + 10 + 10 + 14 = 61 copies on both strands, all of them at the ends; tidk's
    # error-tolerant count can differ from the exact count by a copy or two, and the
    # count in the windows from the count in the whole contigs by as much.
    assert 55 <= int(rows["AACCCT"]["copies_at_ends"]) <= 65
    assert 55 <= int(rows["AACCCT"]["copies_genome_wide"]) <= 65
    assert rows["AAGAAT"]["telomere_ends"] == "1"
    assert rows["AAGAAT"]["opposite_ends"] == "1"
    assert rows["AAGAAT"]["score"] == "0"
    assert rows["AAGAAT"]["chosen"] == "no"
    assert rows["AAGAAT"]["in_literature"] == "-" and rows["AAGAAT"]["in_database"] == "-"
    assert "ACAT" not in rows                                       # 4 bp: dropped by the length filter
    assert "AAAT" not in rows
    for unit, row in rows.items():
        assert row["chosen"] == ("yes" if unit == "AACCCT" else "no")

    ends = read_ends(paths.out_ends)
    assert len(ends) == 2 * len(CONTIGS)
    # The arrays sit 60 bases from the start and 30 from the end; the error-tolerant
    # search may take a base or two of the flank into the array.
    assert ends[("chrom_a", "start")]["status"] == "present"
    assert 50 <= int(ends[("chrom_a", "start")]["offset"]) <= 60
    assert ends[("chrom_a", "end")]["status"] == "present"
    assert 20 <= int(ends[("chrom_a", "end")]["offset"]) <= 30
    assert ends[("chrom_b", "start")]["status"] == "present"
    assert ends[("chrom_b", "end")]["status"] == "present"
    assert ends[("short_dup", "start")]["status"] == "present"
    assert ends[("short_dup", "end")]["status"] == "other_end"
    assert ends[("sat_one_way", "start")]["status"] == "absent"
    assert ends[("plain", "end")]["status"] == "absent"
    assert ends[("plain", "end")]["offset"] == "NA"

    captured = capsys.readouterr()
    assert captured.out == ""
    lines = captured.err.splitlines()
    assert lines[0].startswith(f"telomere_units.py: sample {SAMPLE}: 6 contigs")
    assert lines[1] == (
        f"[FunFlux] telomere_units {SAMPLE}: telomere repeat unit AACCCT (TTAGGG read on the "
        "other strand): 5 contig ends, 2 contigs with a telomere at both ends; known from the "
        "literature (filamentous Ascomycota, most Basidiomycota, Mucoromycota; Botrytis cinerea) "
        "and from 2 fungal species in the database"
    )
    assert "[FunFlux] NOTE" not in captured.err


@needs_tidk
def test_a_wrong_configured_unit_is_written_and_the_note_names_the_data_choice(tmp_path, capsys):
    """With telomere_motif TTAGGGTCAACA the files carry that unit, and the NOTE says the data favour AACCCT."""
    paths = make_inputs(tmp_path)
    telomere_units.main(full_argv(paths, **{"--telomere-motif": "ttagggtcaaca"}))

    configured = telomere_units.canonical_unit("TTAGGGTCAACA")
    assert read_unit(paths.out_unit) == configured == "AACATTAGGGTC"
    rows = read_candidates(paths.out_candidates)
    assert rows[configured]["chosen"] == "yes"
    assert rows[configured]["telomere_ends"] == "0"
    assert rows["AACCCT"]["chosen"] == "no"
    assert rows["AACCCT"]["telomere_ends"] == "5"
    ends = read_ends(paths.out_ends)
    assert len(ends) == 2 * len(CONTIGS)
    assert all(row["status"] == "absent" for row in ends.values())

    captured = capsys.readouterr().err
    assert (
        "telomere repeat unit AACATTAGGGTC (TTAGGGTCAACA read from another position; set by "
        "telomere_motif): 0 contig ends, 0 contigs with a telomere at both ends; known from the "
        "literature (Aspergillus flavus group)"
    ) in captured
    note_lines = [line for line in captured.splitlines() if line.startswith("[FunFlux] NOTE")]
    assert len(note_lines) == 1
    assert note_lines[0] == (
        f"[FunFlux] NOTE telomere_units {SAMPLE}: the data favour AACCCT (5 contig ends in the "
        "telomere orientation, 0 opposite; known) over the configured unit TTAGGGTCAACA "
        "(AACATTAGGGTC), which has 0 contig ends in the telomere orientation and 0 opposite; "
        "set telomere_motif to auto or to that unit and re-run"
    )


@needs_tidk
def test_a_configured_unit_the_data_agree_with_gives_no_note(tmp_path, capsys):
    """telomere_motif CCCTAA (the other strand of TTAGGG) is AACCCT, the data's own choice: no NOTE."""
    paths = make_inputs(tmp_path)
    telomere_units.main(full_argv(paths, **{"--telomere-motif": "CCCTAA"}))
    assert read_unit(paths.out_unit) == "AACCCT"
    captured = capsys.readouterr().err
    assert "AACCCT (TTAGGG read on the other strand; set by telomere_motif): 5 contig ends" in captured
    assert "[FunFlux] NOTE" not in captured


@needs_tidk
def test_an_unknown_unit_is_chosen_from_the_data_and_gets_a_note(tmp_path, capsys):
    """A genome whose chromosome ends carry a made-up 8-bp unit: that unit wins, and the NOTE says it is unknown."""
    unit = "ACAGTTGG"                                   # three G, one C: a G-rich strand
    unit_c_rich = telomere_units.reverse_complement(unit)
    contigs = [
        ("odd_a", unit_c_rich * 10 + random_dna(5000, 60) + unit * 10),
        ("odd_b", unit_c_rich * 9 + random_dna(4000, 61) + unit * 8),
        ("odd_c", random_dna(3000, 62)),
    ]
    paths = make_inputs(tmp_path)
    write_fasta(paths.contigs, contigs)
    telomere_units.main(full_argv(paths))

    canonical = telomere_units.canonical_unit(unit)
    assert read_unit(paths.out_unit) == canonical
    rows = read_candidates(paths.out_candidates)
    assert rows[canonical]["telomere_ends"] == "4"
    assert rows[canonical]["t2t_contigs"] == "2"
    assert rows[canonical]["in_literature"] == "-" and rows[canonical]["in_database"] == "-"
    assert rows[canonical]["chosen"] == "yes"
    captured = capsys.readouterr().err
    assert f"telomere repeat unit {canonical}: 4 contig ends, 2 contigs with a telomere at both ends; not known from the literature or the database" in captured
    assert f"[FunFlux] NOTE telomere_units {SAMPLE}: the telomere repeat unit {canonical} is not known" in captured


@needs_tidk
def test_no_unit_when_no_candidate_is_eligible(tmp_path, capsys):
    """Random contigs without any array: the unit is none, the ends table is the header alone."""
    paths = make_inputs(tmp_path)
    write_fasta(paths.contigs, [("r1", random_dna(6000, 70)), ("r2", random_dna(4000, 71))])
    telomere_units.main(full_argv(paths))
    assert read_unit(paths.out_unit) == "none"
    with open(paths.out_ends) as table:
        assert table.read() == "\t".join(telomere_units.ENDS_COLUMNS) + "\n"
    with open(paths.out_candidates) as table:
        assert table.readline().startswith("unit\tunit_length")
    captured = capsys.readouterr().err
    assert f"[FunFlux] telomere_units {SAMPLE}: no telomere repeat identified" in captured
    assert "[FunFlux] NOTE" not in captured


@needs_tidk
def test_command_line_run_writes_the_three_files(tmp_path):
    """Exit code 0, the three files exist, the lines are on stderr and nothing on stdout."""
    paths = make_inputs(tmp_path)
    finished = subprocess.run(
        [sys.executable, SCRIPT_PATH] + full_argv(paths), capture_output=True, text=True,
    )
    assert finished.returncode == 0, finished.stderr
    assert finished.stdout == ""
    assert finished.stderr.startswith(f"telomere_units.py: sample {SAMPLE}: 6 contigs")
    assert f"[FunFlux] telomere_units {SAMPLE}: telomere repeat unit AACCCT" in finished.stderr
    assert os.path.exists(paths.out_candidates)
    assert read_unit(paths.out_unit) == "AACCCT"
    assert len(read_ends(paths.out_ends)) == 2 * len(CONTIGS)


@needs_tidk
def test_command_line_run_with_a_broken_database_exits_with_code_1(tmp_path):
    """The step must fail (non-zero exit) with a readable message and no traceback."""
    paths = make_inputs(tmp_path)
    with open(paths.database, "w") as table:
        table.write("species\tkingdom\n")
    finished = subprocess.run(
        [sys.executable, SCRIPT_PATH] + full_argv(paths), capture_output=True, text=True,
    )
    assert finished.returncode == 1
    assert "[FunFlux] telomere_units: error:" in finished.stderr
    assert "Traceback" not in finished.stderr
    assert not os.path.exists(paths.out_unit)
