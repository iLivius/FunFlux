"""Tests for unite_identified.py, the script that keeps the UNITE records with a genus.

How to run (from the repository root):

    python -m pytest workflow/scripts/30_taxonomy -q

Only the Python standard library and pytest are needed. No real UNITE file is needed:
every input is made up while the tests run.

THE MADE-UP DATA
----------------
  * A record is a (header, sequence) pair. The header follows the UNITE SINTAX
    release, accession|SH-code;tax=<ranks>; for example
    ACC1|SH0000001.10FU;tax=d:Fungi,p:Ascomycota,g:Botrytis,s:Botrytis_cinerea;
    The helper header(ranks) builds one from the ranks alone.
  * The sequences are short made-up DNA; the script never looks at them.
"""

import os
import subprocess
import sys

import pytest

# The script sits in the same folder as this test file. Putting that folder on the
# list of places where Python looks for modules makes "import unite_identified" work
# wherever pytest is started from.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TEST_DIR)

import unite_identified  # noqa: E402  (this import must come after the line above)

SCRIPT_PATH = os.path.join(TEST_DIR, "unite_identified.py")

# Four made-up taxonomies: a species (genus and species given), a genus without a
# species, a family (stops above the genus) and the kingdom alone.
RANKS_SPECIES = "d:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales,f:Sclerotiniaceae,g:Botrytis,s:Botrytis_cinerea"
RANKS_GENUS_ONLY = "d:Fungi,p:Ascomycota,c:Sordariomycetes,o:Hypocreales,f:Nectriaceae,g:Fusarium"
RANKS_FAMILY = "d:Fungi,p:Ascomycota,c:Leotiomycetes,o:Helotiales,f:Sclerotiniaceae"
RANKS_KINGDOM = "d:Fungi"


# ======================================================================
# HELPERS: MAKING INPUTS AND READING OUTPUTS
# ======================================================================

def header(ranks, accession="ACC1"):
    """Return a UNITE SINTAX header (without the ">") for these ranks."""
    return f"{accession}|SH0000001.10FU;tax={ranks};"


def write_fasta_file(path, records):
    """Write (header, sequence) pairs as FASTA, one sequence line per record."""
    with open(path, "w") as fasta:
        for record_header, sequence in records:
            fasta.write(f">{record_header}\n{sequence}\n")


def make_paths(folder):
    """Return the input and output paths inside `folder`."""
    return folder / "unite_its_sintax.fasta", folder / "unite_its_sintax.identified.fasta"


def run_script(fasta_path, out_path):
    """Run unite_identified.main() on the two paths."""
    unite_identified.main(["--fasta", str(fasta_path), "--out-fasta", str(out_path)])


def read_fasta_back(path):
    """Read a FASTA file back as a list of (header, joined sequence) pairs."""
    with open(path) as fasta:
        lines = fasta.read().splitlines()
    records = []
    for line in lines:
        if line.startswith(">"):
            records.append([line[1:], ""])
        else:
            records[-1][1] += line
    return [(record_header, sequence) for record_header, sequence in records]


# ======================================================================
# WHICH HEADERS NAME A GENUS
# ======================================================================

@pytest.mark.parametrize(
    "ranks, expected",
    [
        (RANKS_SPECIES, True),       # genus and species
        (RANKS_GENUS_ONLY, True),    # genus, no species
        (RANKS_FAMILY, False),       # stops at the family
        (RANKS_KINGDOM, False),      # kingdom alone: tax=d:Fungi;
        ("d:unidentified", False),
        ("d:Fungi,p:Ascomycota,g:", False),   # "g:" with nothing after it is no genus
        ("d:Fungi,p:Ascomycota,g:,s:Botrytis_cinerea", False),   # an empty genus before a species is no genus either
    ],
    ids=["species", "genus_only", "family", "kingdom", "unidentified", "empty_genus",
         "empty_genus_before_species"],
)
def test_header_names_a_genus(ranks, expected):
    """A header names a genus when a rank reads g:<name>; the tax= field is always there."""
    names_a_genus, has_tax_field = unite_identified.header_names_a_genus(">" + header(ranks))
    assert names_a_genus is expected
    assert has_tax_field is True


def test_header_without_a_tax_field():
    """A header without tax= names no genus, and the second value says the field is missing."""
    names_a_genus, has_tax_field = unite_identified.header_names_a_genus(">ACC1|SH0000001.10FU")
    assert names_a_genus is False
    assert has_tax_field is False


def test_g_colon_outside_the_taxonomy_does_not_count():
    """Only the ranks of the tax= field are searched: a g: in the accession is ignored."""
    line = ">g:fake|SH0000001.10FU;tax=d:Fungi;"
    names_a_genus, has_tax_field = unite_identified.header_names_a_genus(line)
    assert names_a_genus is False
    assert has_tax_field is True


def test_taxonomy_without_the_closing_semicolon_is_read_too():
    """A header that ends right after the species, without ";", still names its genus."""
    line = ">ACC1|SH0000001.10FU;tax=" + RANKS_SPECIES
    names_a_genus, has_tax_field = unite_identified.header_names_a_genus(line)
    assert names_a_genus is True


# ======================================================================
# THE FILTER ON A FILE
# ======================================================================

def test_records_with_a_genus_are_kept_and_the_others_dropped(tmp_path, capsys):
    """Species and genus-only records are kept; family, kingdom and no-tax records are dropped."""
    fasta_path, out_path = make_paths(tmp_path)
    records = [
        (header(RANKS_SPECIES, "ACC1"), "ACGTACGT"),
        (header(RANKS_FAMILY, "ACC2"), "ACGTACGA"),
        (header(RANKS_GENUS_ONLY, "ACC3"), "ACGTACGC"),
        (header(RANKS_KINGDOM, "ACC4"), "ACGTACGG"),
        ("ACC5|SH0000005.10FU", "ACGTACTT"),            # no tax= field at all
    ]
    write_fasta_file(fasta_path, records)

    run_script(fasta_path, out_path)

    assert read_fasta_back(out_path) == [records[0], records[2]]
    captured = capsys.readouterr()
    assert captured.out == "kept_with_genus: 2\ndropped_without_genus: 3\n"
    assert captured.err.count("\n") == 1
    assert "5 records read, 2 kept" in captured.err
    assert "3 dropped" in captured.err
    assert "1 of them without a tax= field" in captured.err


def test_kept_records_are_copied_as_they_are_including_multi_line_sequences(tmp_path):
    """A sequence spread over several lines is copied whole, line breaks and case included."""
    fasta_path, out_path = make_paths(tmp_path)
    text = (
        f">{header(RANKS_SPECIES, 'ACC1')}\n"
        "ACGTACGTAC\n"
        "gtacgtacgt\n"
        "ACG\n"
        f">{header(RANKS_KINGDOM, 'ACC2')}\n"
        "TTTT\n"
        "TTTT\n"
        f">{header(RANKS_GENUS_ONLY, 'ACC3')}\n"
        "GGGG\n"
    )
    fasta_path.write_text(text)

    run_script(fasta_path, out_path)

    expected = (
        f">{header(RANKS_SPECIES, 'ACC1')}\n"
        "ACGTACGTAC\n"
        "gtacgtacgt\n"
        "ACG\n"
        f">{header(RANKS_GENUS_ONLY, 'ACC3')}\n"
        "GGGG\n"
    )
    assert out_path.read_text() == expected


def test_blank_lines_are_skipped(tmp_path):
    """Blank lines inside the file are neither copied nor counted."""
    fasta_path, out_path = make_paths(tmp_path)
    fasta_path.write_text(f">{header(RANKS_SPECIES)}\n\nACGT\n\n")

    run_script(fasta_path, out_path)

    assert out_path.read_text() == f">{header(RANKS_SPECIES)}\nACGT\n"


def test_windows_line_endings_are_understood(tmp_path, capsys):
    """A file written on Windows (lines end with \\r\\n) is read like any other."""
    fasta_path, out_path = make_paths(tmp_path)
    fasta_path.write_bytes(
        f">{header(RANKS_SPECIES)}\r\nACGT\r\n>{header(RANKS_KINGDOM, 'ACC2')}\r\nACGT\r\n".encode()
    )

    run_script(fasta_path, out_path)

    assert capsys.readouterr().out == "kept_with_genus: 1\ndropped_without_genus: 1\n"
    assert read_fasta_back(out_path) == [(header(RANKS_SPECIES), "ACGT")]


def test_output_is_overwritten_and_not_appended(tmp_path):
    """Running twice on the same file leaves exactly the same output."""
    fasta_path, out_path = make_paths(tmp_path)
    write_fasta_file(fasta_path, [(header(RANKS_SPECIES), "ACGT")])

    run_script(fasta_path, out_path)
    first = out_path.read_text()
    run_script(fasta_path, out_path)

    assert out_path.read_text() == first


# ======================================================================
# BAD INPUT GIVES A CLEAR ERROR, AND NO OUTPUT FILE
# ======================================================================

def test_file_without_any_genus_is_an_error(tmp_path, capsys):
    """A file in which no record names a genus is not a usable reference: stop, say so."""
    fasta_path, out_path = make_paths(tmp_path)
    write_fasta_file(fasta_path, [(header(RANKS_KINGDOM), "ACGT"), (header(RANKS_FAMILY, "ACC2"), "ACGT")])

    with pytest.raises(SystemExit) as stop:
        run_script(fasta_path, out_path)

    assert stop.value.code == 1
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.startswith("unite_identified.py: error:")
    assert "no record carries a genus rank" in captured.err
    assert not out_path.exists()


def test_empty_file_is_an_error(tmp_path, capsys):
    """An empty file (a failed download, say) is refused the same way."""
    fasta_path, out_path = make_paths(tmp_path)
    fasta_path.write_text("")

    with pytest.raises(SystemExit) as stop:
        run_script(fasta_path, out_path)

    assert stop.value.code == 1
    assert "no record carries a genus rank" in capsys.readouterr().err
    assert not out_path.exists()


def test_sequence_before_any_header_is_an_error(tmp_path, capsys):
    """A file that starts with sequence and not with ">" is not FASTA: stop and say where."""
    fasta_path, out_path = make_paths(tmp_path)
    fasta_path.write_text(f"ACGTACGT\n>{header(RANKS_SPECIES)}\nACGT\n")

    with pytest.raises(SystemExit) as stop:
        run_script(fasta_path, out_path)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert "line 1" in message
    assert "before any FASTA header" in message
    assert not out_path.exists()


def test_missing_input_file_is_a_clear_error(tmp_path, capsys):
    """An input file that does not exist gives a one-line message, not a traceback."""
    fasta_path, out_path = make_paths(tmp_path)

    with pytest.raises(SystemExit) as stop:
        run_script(fasta_path, out_path)

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("unite_identified.py: error:")
    assert "unite_its_sintax.fasta" in message
    assert message.count("\n") == 1
    assert not out_path.exists()


def test_missing_option_is_refused(tmp_path, capsys):
    """Leaving out --out-fasta is a command-line error (exit code 2)."""
    fasta_path, out_path = make_paths(tmp_path)

    with pytest.raises(SystemExit) as stop:
        unite_identified.main(["--fasta", str(fasta_path)])

    assert stop.value.code == 2
    assert "--out-fasta" in capsys.readouterr().err


# ======================================================================
# THE SCRIPT AS SNAKEMAKE WILL CALL IT
# ======================================================================

def test_command_line_run_prints_the_two_count_lines(tmp_path):
    """Exit code 0; standard output holds exactly the two count lines the rule saves."""
    fasta_path, out_path = make_paths(tmp_path)
    write_fasta_file(fasta_path, [
        (header(RANKS_SPECIES, "ACC1"), "ACGT"),
        (header(RANKS_KINGDOM, "ACC2"), "ACGT"),
        (header(RANKS_GENUS_ONLY, "ACC3"), "ACGT"),
    ])

    finished = subprocess.run(
        [sys.executable, SCRIPT_PATH, "--fasta", str(fasta_path), "--out-fasta", str(out_path)],
        capture_output=True, text=True,
    )

    assert finished.returncode == 0
    assert finished.stdout == "kept_with_genus: 2\ndropped_without_genus: 1\n"
    assert "3 records read, 2 kept" in finished.stderr
    assert len(read_fasta_back(out_path)) == 2


def test_command_line_run_on_a_file_without_a_genus_exits_with_code_1(tmp_path):
    """The step must fail (non-zero exit) with a readable message and no output file."""
    fasta_path, out_path = make_paths(tmp_path)
    write_fasta_file(fasta_path, [(header(RANKS_KINGDOM), "ACGT")])

    finished = subprocess.run(
        [sys.executable, SCRIPT_PATH, "--fasta", str(fasta_path), "--out-fasta", str(out_path)],
        capture_output=True, text=True,
    )

    assert finished.returncode == 1
    assert "unite_identified.py: error:" in finished.stderr
    assert "Traceback" not in finished.stderr
    assert not out_path.exists()
