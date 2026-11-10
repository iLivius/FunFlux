"""Tests for contig_summary.py, the script that says what every contig of the assembly is.

How to run (from the repository root):

    python -m pytest workflow/scripts/25_contig_summary -q

Only the Python standard library and pytest are needed. No real assembly is needed:
every input is made up while the tests run.

THE MADE-UP DATA
----------------
A small assembly of nine contigs, each built to trigger one verdict rule:
  * chrom_both   12 kb of random sequence with a telomere array at each end: 15 copies
                 of CCCTAA at the start (after 60 bases of other sequence) and 12
                 copies of TTAGGG at the end (before 30 bases). Normal depth, 5 genes,
                 1 BUSCO marker, 1 rDNA unit (an rDNA locus on a chromosome must not
                 make it an rDNA array).
  * arm_one      8 kb with 8 copies of CCCTAA at the start only. Normal depth, 3 genes.
  * rdna         a made-up 1,000-base rDNA unit repeated 10 times in tandem, with 4
                 ITS windows on it, low depth, and no scaffold in the annotated
                 genome (funannotate clean dropped it).
  * mito         6 kb that the mitochondrion screen confirmed, circular, very deep.
  * removed      5 kb that the contamination screen removed as Escherichia.
  * plain        7 kb of random sequence, normal depth, 4 genes, 2 BUSCO markers.
  * lowdepth     5 kb at a tenth of the median depth.
  * repeat_deep  1.2 kb of random sequence at six times the median depth, with a
                 scaffold in the annotated genome but no gene and no BUSCO marker on
                 it: a collapsed repeat.
  * fragment     1.5 kb of random sequence at normal depth, with a scaffold but no
                 gene and no BUSCO marker: a piece too short to carry a gene.
The tables around them (decisions, mitochondrion screen, BlobTools table, assembler
evidence, ITS windows, BUSCO folder, GFF3, annotated genome) are written in the exact
layouts the workflow produces, and so are the two files of the telomere step: a
telomere_unit.txt holding AACCCT (the canonical form of TTAGGG) and a
telomere_ends.tsv in tidk's layout, with the arrays of chrom_both and arm_one as
present rows, absent rows for rdna and mito, and no row for the other contigs (a
contig without rows has no telomere at either end). The script never looks at the
sequence for telomeres: the arrays in the sequences only make the made-up assembly
look like a real one. make_inputs() builds all of it in a temporary folder.
"""

import csv
import os
import random
import subprocess
import sys
from types import SimpleNamespace

import pytest

# SimpleNamespace (used by make_inputs below) is a plain container whose fields are
# read with a dot, like made.contigs.

# The script sits in the same folder as this test file. Putting that folder on the
# list of places where Python looks for modules makes "import contig_summary" work
# wherever pytest is started from.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TEST_DIR)

import contig_summary  # noqa: E402  (this import must come after the line above)

SCRIPT_PATH = os.path.join(TEST_DIR, "contig_summary.py")

TELOMERE = "TTAGGG"
TELOMERE_RC = "CCCTAA"
TELOMERE_UNIT = "AACCCT"     # the canonical form of TTAGGG, as telomere_units.py writes it
SAMPLE = "SAMPLE1"

# The columns of telomere_ends.tsv, as tidk writes them.
TELOMERE_ENDS_COLUMNS = [
    "id", "sequence_length", "side", "status", "telomere_start", "telomere_end",
    "telomere_length", "offset", "strand", "window_limited",
]


# ======================================================================
# THE MADE-UP DATA
# ======================================================================

def random_dna(length, seed):
    """Return `length` random bases; the seed makes the sequence the same on every run."""
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


# The nine contigs (name -> sequence). Built once; every test reads them.
RDNA_UNIT = random_dna(1000, 6)
CONTIGS = {
    "chrom_both": random_dna(60, 1) + TELOMERE_RC * 15 + random_dna(12000, 2) + TELOMERE * 12 + random_dna(30, 3),
    "arm_one": random_dna(40, 4) + TELOMERE_RC * 8 + random_dna(8000, 5),
    "rdna": RDNA_UNIT * 10,
    "mito": random_dna(6000, 7),
    "removed": random_dna(5000, 8),
    "plain": random_dna(7000, 9),
    "lowdepth": random_dna(5000, 10),
    "repeat_deep": random_dna(1200, 70),
    "fragment": random_dna(1500, 71),
}

# Read depth per contig, as the BlobTools table reports it. The length-weighted
# median over the nuclear, kept contigs (all shorter than 100 kb, so all count) is
# 23.0: sorted by depth, lowdepth (5,000 bases), rdna (10,000) and fragment (1,500)
# add up to 16,500, and plain (7,000) carries the sum past half of the 45,040 bases.
DEPTHS = {
    "chrom_both": 25.0, "arm_one": 24.0, "rdna": 8.0, "mito": 180.0,
    "removed": 40.0, "plain": 23.0, "lowdepth": 2.0,
    "repeat_deep": 140.0, "fragment": 20.0,
}
MEDIAN_DEPTH = 23.0

# Which scaffold_N each contig became in the annotated genome (rdna was dropped by
# funannotate clean; mito and removed never reach funannotate). The two short contigs
# have a scaffold, as short contigs do in a real run, and no gene on it.
SCAFFOLDS = {
    "scaffold_1": "chrom_both", "scaffold_2": "arm_one", "scaffold_3": "plain",
    "scaffold_4": "lowdepth", "scaffold_5": "fragment", "scaffold_6": "repeat_deep",
}

# Genes (mRNA features) per scaffold in the GFF3. A scaffold not listed has none.
GENES = {"scaffold_1": 5, "scaffold_2": 3, "scaffold_3": 4}


def write_fasta(path, records):
    """Write (name, sequence) pairs as FASTA, wrapped at 80 bases."""
    with open(path, "w") as fasta:
        for name, sequence in records:
            fasta.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                fasta.write(sequence[start:start + 80] + "\n")
            if sequence == "":
                fasta.write("\n")


def write_tsv(path, header, rows):
    """Write a header line and rows (lists of values) as a tab-separated table."""
    with open(path, "w") as table:
        table.write("\t".join(header) + "\n")
        for row in rows:
            table.write("\t".join(str(value) for value in row) + "\n")


def write_blob_table(path, depths):
    """Write a BlobTools table in the real layout: ## lines, a # name header, one row per contig."""
    with open(path, "w") as table:
        table.write("## 1.1.1\n## assembly\t: assembly.fasta\n## taxrule\t: bestsum\n##\n")
        table.write("# name\tlength\tGC\tN\tbam0\tsuperkingdom.t.6%s\tsuperkingdom.s.7%s\n")
        for contig, depth in depths.items():
            table.write(f"{contig}\t{len(CONTIGS[contig])}\t0.4\t0\t{depth}\tEukaryota\t100.0\n")


def write_busco_dir(folder):
    """Write a BUSCO folder with two runs; the one with the larger marker set must be used."""
    specific = os.path.join(folder, "run_fungi_odb12")
    generic = os.path.join(folder, "run_eukaryota_odb12")
    os.makedirs(specific)
    os.makedirs(generic)
    busco_header = "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\tOrthoDB url\tDescription\n"
    with open(os.path.join(specific, "full_table.tsv"), "w") as table:
        table.write("# BUSCO version is: 6.0.0 \n")
        table.write("# The lineage dataset is: fungi_odb12 (Creation date: 2026-05-22, number of genomes: 100, number of BUSCOs: 4)\n")
        table.write(busco_header)
        table.write("1at4751\tComplete\tplain\t100\t900\t+\t500.0\t200\turl\tdesc\n")
        table.write("2at4751\tComplete\tplain\t2000\t2900\t-\t400.0\t180\turl\tdesc\n")
        table.write("3at4751\tDuplicated\tchrom_both\t3000\t3900\t+\t300.0\t150\turl\tdesc\n")
        table.write("3at4751\tDuplicated\tarm_one\t3000\t3900\t+\t300.0\t150\turl\tdesc\n")
        table.write("4at4751\tFragmented\tplain\t5000\t5300\t+\t50.0\t60\turl\tdesc\n")
        table.write("5at4751\tMissing\n")
    with open(os.path.join(generic, "full_table.tsv"), "w") as table:
        table.write("# BUSCO version is: 6.0.0 \n")
        table.write("# The lineage dataset is: eukaryota_odb12 (Creation date: 2026-05-22, number of genomes: 456, number of BUSCOs: 2)\n")
        table.write(busco_header)
        for marker in range(5):
            table.write(f"{marker}at2759\tComplete\tplain\t100\t900\t+\t500.0\t200\turl\tdesc\n")


def write_gff3(path, genes_per_scaffold):
    """Write a funannotate-like GFF3: gene, mRNA, exon and CDS lines per gene, plus a tRNA."""
    with open(path, "w") as gff3:
        gff3.write("##gff-version 3\n")
        for scaffold, n_genes in genes_per_scaffold.items():
            for number in range(1, n_genes + 1):
                gene_id = f"{scaffold}_g{number}"
                start, end = number * 1000, number * 1000 + 500
                gff3.write(f"{scaffold}\tfunannotate\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id};\n")
                gff3.write(f"{scaffold}\tfunannotate\tmRNA\t{start}\t{end}\t.\t+\t.\tID={gene_id}-T1;Parent={gene_id};\n")
                gff3.write(f"{scaffold}\tfunannotate\texon\t{start}\t{end}\t.\t+\t.\tID={gene_id}-T1.exon1;Parent={gene_id}-T1;\n")
                gff3.write(f"{scaffold}\tfunannotate\tCDS\t{start}\t{end}\t.\t+\t0\tID={gene_id}-T1.cds;Parent={gene_id}-T1;\n")
        gff3.write("scaffold_1\tfunannotate\tgene\t9000\t9072\t.\t-\t.\tID=trna1;\n")
        gff3.write("scaffold_1\tfunannotate\ttRNA\t9000\t9072\t.\t-\t.\tID=trna1-T1;Parent=trna1;\n")


def telomere_row(contig, contig_length, side, status, telomere_length=None, offset=None):
    """One row of telomere_ends.tsv in tidk's layout, for one end of one contig.

    A row with an array (present, wrong_strand, other_end, not_terminal) gets the
    array's coordinates worked out from its length and its distance from the contig
    end (the offset), and the strand tidk reports: C-rich at the start and G-rich at
    the end for a telomere, the other way round for wrong_strand. An absent row, or
    a row given without a length, holds NA.
    """
    if telomere_length is None:
        return [contig, contig_length, side, status, "NA", "NA", "NA", "NA", "NA", "NA"]
    if side == "start":
        array_start = offset
        strand = "C-rich"
    else:
        array_start = contig_length - offset - telomere_length
        strand = "G-rich"
    if status == "wrong_strand":
        strand = "G-rich" if strand == "C-rich" else "C-rich"
    array_end = array_start + telomere_length
    return [contig, contig_length, side, status, array_start, array_end, telomere_length, offset, strand, "false"]


def write_telomere_files(unit_path, ends_path, unit_text, rows):
    """Write telomere_unit.txt (one line) and telomere_ends.tsv (tidk's header and the rows)."""
    with open(unit_path, "w") as unit_file:
        unit_file.write(unit_text + "\n")
    write_tsv(ends_path, TELOMERE_ENDS_COLUMNS, rows)


# The telomere table of the made-up assembly: the arrays of chrom_both (15 and 12
# copies of a 6-bp unit, 60 and 30 bp from the ends) and arm_one (8 copies, 40 bp from
# the start), absent rows for rdna and mito, no row for the other contigs.
def made_up_telomere_rows():
    """The rows of telomere_ends.tsv for the made-up assembly."""
    return [
        telomere_row("chrom_both", len(CONTIGS["chrom_both"]), "start", "present", 15 * 6, 60),
        telomere_row("chrom_both", len(CONTIGS["chrom_both"]), "end", "present", 12 * 6, 30),
        telomere_row("arm_one", len(CONTIGS["arm_one"]), "start", "present", 8 * 6, 40),
        telomere_row("arm_one", len(CONTIGS["arm_one"]), "end", "absent"),
        telomere_row("rdna", len(CONTIGS["rdna"]), "start", "absent"),
        telomere_row("rdna", len(CONTIGS["rdna"]), "end", "absent"),
        telomere_row("mito", len(CONTIGS["mito"]), "start", "absent"),
        telomere_row("mito", len(CONTIGS["mito"]), "end", "absent"),
    ]


def make_inputs(folder):
    """Write every input file of the made-up assembly into `folder`; return their paths."""
    paths = SimpleNamespace()
    paths.contigs = os.path.join(folder, "assembly.fasta")
    write_fasta(paths.contigs, list(CONTIGS.items()))

    paths.telomere_unit = os.path.join(folder, "telomere_unit.txt")
    paths.telomere_ends = os.path.join(folder, "telomere_ends.tsv")
    write_telomere_files(paths.telomere_unit, paths.telomere_ends, TELOMERE_UNIT, made_up_telomere_rows())

    paths.mito_contigs = os.path.join(folder, "contigs_mitochondrial.fasta")
    write_fasta(paths.mito_contigs, [("mito", CONTIGS["mito"])])

    paths.decisions = os.path.join(folder, "contig_taxonomy_decisions.tsv")
    decision_rows = []
    for contig in CONTIGS:
        if contig == "removed":
            decision_rows.append([contig, "Escherichia", "remove", "excluded_genus"])
        elif contig == "mito":
            decision_rows.append([contig, "Botrytis", "keep", "protected_confirmed_mitochondrial"])
        else:
            decision_rows.append([contig, "Botrytis", "keep", "mode_off"])
    write_tsv(paths.decisions, ["contig", "assigned_genus", "action", "reason"], decision_rows)

    paths.mito_screen = os.path.join(folder, "mitochondrial_screen.tsv")
    screen_rows = []
    for contig in CONTIGS:
        if contig == "mito":
            screen_rows.append([contig, len(CONTIGS[contig]), "confirmed_mitochondrial", "direct_blast_plus_high_coverage", "99.00"])
        elif contig == "plain":
            screen_rows.append([contig, len(CONTIGS[contig]), "mitochondrial_candidate", "weak_direct_blast", "12.00"])
        else:
            screen_rows.append([contig, len(CONTIGS[contig]), "not_mitochondrial", "no_mitochondrial_evidence", ""])
    write_tsv(paths.mito_screen, ["contig", "length", "classification", "reason", "best_mito_qcovs"], screen_rows)

    paths.blob_table = os.path.join(folder, "bestscore.blob.blobDB.table.txt")
    write_blob_table(paths.blob_table, DEPTHS)

    paths.assembly_evidence = os.path.join(folder, "assembly_evidence.tsv")
    evidence_rows = []
    for contig in CONTIGS:
        circular = "Y" if contig == "mito" else "N"
        coverage = "420.000000" if contig == "mito" else "100.000000"
        evidence_rows.append([contig, "hifiasm", circular, coverage, "N", "*"])
    write_tsv(paths.assembly_evidence, ["contig", "source", "circular", "coverage", "repeat", "alt_group"], evidence_rows)

    paths.its_windows = os.path.join(folder, "its_windows.tsv")
    its_rows = []
    for number in range(4):
        start = number * 1000 + 1
        its_rows.append([f"rdna:{start}-{start + 999}_+", "rdna", start, start + 999, "+", 12, "1.00e-15", 1000])
    its_rows.append(["chrom_both:5001-9001_-", "chrom_both", 5001, 9001, "-", 30, "2.00e-15", 4001])
    write_tsv(paths.its_windows, ["window_id", "contig", "start", "end", "strand", "anchor_hits", "best_evalue", "window_length"], its_rows)

    paths.busco_dir = os.path.join(folder, "busco")
    write_busco_dir(paths.busco_dir)

    paths.gff3 = os.path.join(folder, "annotation.gff3")
    write_gff3(paths.gff3, GENES)

    paths.annotated_fasta = os.path.join(folder, "contigs_sort.fasta")
    write_fasta(paths.annotated_fasta, [(scaffold, CONTIGS[contig]) for scaffold, contig in SCAFFOLDS.items()])

    paths.out = os.path.join(folder, "contig_summary.tsv")
    paths.out_mqc = os.path.join(folder, "contig_summary_mqc.tsv")
    return paths


def full_argv(paths, **changes):
    """The complete command line for the made-up assembly; `changes` replaces options."""
    options = {
        "--contigs": paths.contigs,
        "--telomere-unit-file": paths.telomere_unit,
        "--telomere-ends": paths.telomere_ends,
        "--mito-contigs": paths.mito_contigs,
        "--decisions": paths.decisions,
        "--mito-screen": paths.mito_screen,
        "--blob-table": paths.blob_table,
        "--assembly-evidence": paths.assembly_evidence,
        "--its-windows": paths.its_windows,
        "--busco-dir": paths.busco_dir,
        "--gff3": paths.gff3,
        "--annotated-fasta": paths.annotated_fasta,
        "--sample": SAMPLE,
        "--out": paths.out,
        "--out-mqc": paths.out_mqc,
    }
    options.update(changes)
    argv = []
    for option, value in options.items():
        if value is not None:
            argv.extend([option, str(value)])
    return argv


def read_table(path):
    """Read the output table back as {contig: row dictionary}."""
    with open(path, newline="") as table:
        rows = list(csv.DictReader(table, delimiter="\t"))
    return {row["contig"]: row for row in rows}


# ======================================================================
# THE BUILDING BLOCKS: THE TWO TELOMERE FILES, THE MEDIAN DEPTH
# ======================================================================

def test_telomere_unit_file_gives_the_unit_or_none(tmp_path):
    """The one-line file reads as the unit in upper case, or as None for the word none."""
    unit_path = tmp_path / "telomere_unit.txt"
    unit_path.write_text("AACCCT\n")
    assert contig_summary.read_telomere_unit(unit_path) == "AACCCT"
    unit_path.write_text("  aaccct  \n\n")
    assert contig_summary.read_telomere_unit(unit_path) == "AACCCT"
    unit_path.write_text("none\n")
    assert contig_summary.read_telomere_unit(unit_path) is None


@pytest.mark.parametrize(
    "text, expected_text",
    [
        ("", "found 0 non-empty lines"),
        ("AACCCT\nAACCT\n", "found 2 non-empty lines"),
        ("TTAGGX\n", "letters A, C, G and T"),
    ],
)
def test_telomere_unit_file_that_is_not_a_unit_is_refused(tmp_path, text, expected_text):
    """An empty file, two lines, or a word that is no unit means the wrong file was given."""
    unit_path = tmp_path / "telomere_unit.txt"
    unit_path.write_text(text)
    with pytest.raises(ValueError) as problem:
        contig_summary.read_telomere_unit(unit_path)
    assert expected_text in str(problem.value)


def test_telomere_ends_table_is_read_per_contig_and_side(tmp_path):
    """Rows become {contig: {side: (status, length, offset)}}; only present rows carry numbers."""
    ends_path = tmp_path / "telomere_ends.tsv"
    write_tsv(ends_path, TELOMERE_ENDS_COLUMNS, [
        telomere_row("chrom_both", len(CONTIGS["chrom_both"]), "start", "present", 90, 60),
        telomere_row("chrom_both", len(CONTIGS["chrom_both"]), "end", "wrong_strand", 72, 30),
        telomere_row("arm_one", len(CONTIGS["arm_one"]), "start", "absent"),
    ])
    ends = contig_summary.read_telomere_ends(ends_path, set(CONTIGS))
    assert ends == {
        "chrom_both": {"start": ("present", 90, 60), "end": ("wrong_strand", None, None)},
        "arm_one": {"start": ("absent", None, None)},
    }
    # A header line alone (the unit is none) gives an empty dictionary.
    write_tsv(ends_path, TELOMERE_ENDS_COLUMNS, [])
    assert contig_summary.read_telomere_ends(ends_path, set(CONTIGS)) == {}


@pytest.mark.parametrize(
    "rows, expected_text",
    [
        ([telomere_row("ghost", 5000, "start", "present", 30, 0)], "'ghost' is not in the draft assembly"),
        ([telomere_row("plain", 7000, "middle", "absent")], "side must be start or end"),
        ([telomere_row("plain", 7000, "start", "found", 30, 0)], "'found' is not one of"),
        ([telomere_row("plain", 7000, "start", "present")], "telomere_length value 'NA' is not a number"),
        ([["plain", 7000, "start", "present", 0, 30, "", "", "C-rich", "false"]], "present array needs its telomere_length and offset"),
        ([telomere_row("plain", 7000, "end", "absent"), telomere_row("plain", 7000, "end", "absent")], "a second row for the end of 'plain'"),
    ],
)
def test_telomere_ends_table_that_does_not_fit_is_refused(tmp_path, rows, expected_text):
    """A row for an unknown contig, an unknown side or status, a present row with NA or empty numbers, a duplicate."""
    ends_path = tmp_path / "telomere_ends.tsv"
    write_tsv(ends_path, TELOMERE_ENDS_COLUMNS, rows)
    with pytest.raises(ValueError) as problem:
        contig_summary.read_telomere_ends(ends_path, set(CONTIGS))
    assert expected_text in str(problem.value)


def test_telomere_columns_turn_array_lengths_into_copies():
    """A present array gives its copies (length over unit length, rounded) and its offset; anything else 0 and None."""
    ends = {"start": ("present", 100, 60), "end": ("wrong_strand", None, None)}
    assert contig_summary.telomere_columns(ends, "AACCCT") == {
        "telomere_start": 17, "telomere_end_offset_start": 60, "telomere_status_start": "present",
        "telomere_end": 0, "telomere_end_offset_end": None, "telomere_status_end": "wrong_strand",
    }
    # 126 bp of a 6-bp unit are 21 copies.
    assert contig_summary.telomere_columns({"end": ("present", 126, 72)}, "AACCCT")["telomere_end"] == 21
    # A contig without rows, with or without a unit, has nothing at either end.
    empty = contig_summary.telomere_columns({}, None)
    assert empty["telomere_start"] == 0 and empty["telomere_end"] == 0
    assert empty["telomere_end_offset_start"] is None and empty["telomere_status_end"] is None


def test_weighted_median_uses_long_contigs_only_when_there_are_some():
    """The median is weighted by length over the contigs of at least 100 kb."""
    def row(length, depth, mito=False, removed=False):
        return {"length": length, "depth": depth, "is_mitochondrion": mito, "is_removed": removed}
    rows = [row(150000, 10.0), row(200000, 30.0), row(1000, 1000.0), row(120000, 500.0, mito=True)]
    assert contig_summary.weighted_median_depth(rows) == 30.0
    # Without any contig of 100 kb every contig counts; a removed contig never does.
    rows = [row(5000, 2.0), row(10000, 8.0), row(7000, 23.0), row(8000, 24.0), row(12000, 25.0), row(50000, 99.0, removed=True)]
    assert contig_summary.weighted_median_depth(rows) == 23.0
    assert contig_summary.weighted_median_depth([row(5000, None)]) is None


# ======================================================================
# THE WHOLE RUN ON THE MADE-UP ASSEMBLY
# ======================================================================

def test_every_contig_gets_the_expected_verdict(tmp_path, capsys):
    """The nine contigs trigger the nine rules they were built for."""
    paths = make_inputs(tmp_path)
    contig_summary.main(full_argv(paths))
    table = read_table(paths.out)

    assert list(table) == list(CONTIGS)   # one row per contig, in FASTA order
    assert table["chrom_both"]["verdict"] == "chromosome"
    assert table["arm_one"]["verdict"] == "chromosome_arm"
    assert table["rdna"]["verdict"] == "rDNA_array"
    assert table["mito"]["verdict"] == "mitochondrion"
    assert table["removed"]["verdict"] == "contaminant"
    assert table["plain"]["verdict"] == "nuclear_contig"
    assert table["lowdepth"]["verdict"] == "low_depth"
    assert table["repeat_deep"]["verdict"] == "repeat"
    assert table["fragment"]["verdict"] == "nuclear_fragment"

    # The reasons name the evidence.
    assert "telomere repeats at both ends (15 and 12 copies)" in table["chrom_both"]["reason"]
    assert "5 genes" in table["chrom_both"]["reason"]
    assert "telomere repeats at the start only (8 copies)" in table["arm_one"]["reason"]
    assert table["rdna"]["reason"].startswith("4 rDNA units; genes unknown")
    assert "confirmed_mitochondrial (direct_blast_plus_high_coverage)" in table["mito"]["reason"]
    assert "circular" in table["mito"]["reason"]
    assert table["removed"]["reason"] == "removed by the screen: genus Escherichia (excluded_genus)"
    assert table["plain"]["reason"].startswith("4 genes, 2 BUSCO markers")
    assert "below 0.5x" in table["lowdepth"]["reason"]
    assert table["repeat_deep"]["reason"] == (
        "depth 6.09x the median (above 1.5x); 0 genes, 0 BUSCO markers; "
        "no telomere repeat at either end"
    )
    assert table["fragment"]["reason"] == (
        "0 genes, 0 BUSCO markers; depth 0.87x the median; no telomere repeat at either end"
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "9 contigs" in captured.err
    assert (
        "1 chromosome, 1 chromosome_arm, 1 contaminant, 1 low_depth, 1 mitochondrion, "
        "1 nuclear_contig, 1 nuclear_fragment, 1 rDNA_array, 1 repeat"
    ) in captured.err
    assert "telomere repeat unit AACCCT at 3 of 18 contig ends" in captured.err
    assert "BUSCO lineage fungi_odb12" in captured.err


def test_the_evidence_columns_hold_the_right_values(tmp_path):
    """Length, GC, depth, depth ratio, assembler columns, telomere columns, counts, genus, screen, annotated_as."""
    paths = make_inputs(tmp_path)
    contig_summary.main(full_argv(paths))
    table = read_table(paths.out)

    assert list(table["plain"]) == contig_summary.COLUMNS
    assert table["plain"]["length"] == "7000"
    assert table["chrom_both"]["length"] == str(len(CONTIGS["chrom_both"]))
    assert float(table["plain"]["gc"]) == round(contig_summary.gc_fraction(CONTIGS["plain"]), 3)
    assert table["rdna"]["gc"] == f"{contig_summary.gc_fraction(RDNA_UNIT):.3f}"

    assert table["chrom_both"]["depth"] == "25.00"
    assert table["chrom_both"]["depth_x_median"] == f"{25.0 / MEDIAN_DEPTH:.2f}"
    assert table["mito"]["depth_x_median"] == f"{180.0 / MEDIAN_DEPTH:.2f}"
    assert table["mito"]["assembler_depth"] == "420.0"
    assert table["mito"]["circular"] == "Y"
    assert table["plain"]["circular"] == "N"

    # The telomere columns come from the two files of the telomere step: the copies
    # of a present array, its offset, and the unit, the same on every row.
    assert table["chrom_both"]["telomere_start"] == "15"
    assert table["chrom_both"]["telomere_end"] == "12"
    assert table["chrom_both"]["telomere_end_offset_start"] == "60"
    assert table["chrom_both"]["telomere_end_offset_end"] == "30"
    assert table["arm_one"]["telomere_end"] == "0"
    assert table["arm_one"]["telomere_end_offset_end"] == "NA"
    assert table["rdna"]["telomere_start"] == "0"              # an absent row
    assert table["plain"]["telomere_start"] == "0"             # no row at all
    assert table["plain"]["telomere_end_offset_start"] == "NA"
    for contig in CONTIGS:
        assert table[contig]["telomere_unit"] == "AACCCT", contig

    assert table["rdna"]["rdna_units"] == "4"
    assert table["chrom_both"]["rdna_units"] == "1"
    assert table["plain"]["rdna_units"] == "0"
    assert table["plain"]["busco_markers"] == "2"        # Complete x2; the Fragmented one does not count
    assert table["chrom_both"]["busco_markers"] == "1"   # a Duplicated marker is complete, counted where it lies
    assert table["chrom_both"]["genes"] == "5"
    assert table["plain"]["genes"] == "4"
    assert table["lowdepth"]["genes"] == "0"             # a scaffold without any mRNA
    assert table["rdna"]["genes"] == "NA"                # no scaffold: dropped by funannotate clean
    assert table["mito"]["genes"] == "NA"

    assert table["removed"]["genus"] == "Escherichia"
    assert table["removed"]["screen"] == "removed (excluded_genus)"
    assert table["chrom_both"]["screen"] == "kept (mode_off)"
    assert table["plain"]["screen"] == "kept (mode_off); mitochondrial_candidate (weak_direct_blast)"
    assert table["mito"]["screen"] == "mitochondrion (direct_blast_plus_high_coverage)"

    assert table["chrom_both"]["annotated_as"] == "scaffold_1"
    assert table["plain"]["annotated_as"] == "scaffold_3"
    assert table["rdna"]["annotated_as"] == "dropped by funannotate clean"
    assert table["mito"]["annotated_as"] == "not annotated (mitochondrion)"
    assert table["removed"]["annotated_as"] == "not annotated (removed by the screen)"


def test_mqc_copy_has_the_header_lines_and_sample_contig_row_names(tmp_path):
    """The MultiQC file starts with the custom-content header, names the unit, and names rows '<sample>: <contig>'."""
    paths = make_inputs(tmp_path)
    contig_summary.main(full_argv(paths))
    with open(paths.out_mqc) as mqc:
        lines = mqc.read().splitlines()

    assert lines[0] == "# id: contig_summary"
    assert lines[1] == "# section_name: 'Contig summary'"
    assert lines[2].startswith("# description: '")
    assert "telomere repeat unit AACCCT" in lines[2]
    assert lines[3] == "# plot_type: 'table'"
    assert lines[4] == "# pconfig: {namespace: 'Contig summary'}"
    assert lines[5] == "\t".join(["Sample"] + contig_summary.COLUMNS)
    assert lines[6].startswith(f"{SAMPLE}: chrom_both\tchrom_both\t")
    assert len(lines) == 6 + len(CONTIGS)

    # The data part is the same as the main table.
    with open(paths.out) as table:
        main_lines = table.read().splitlines()
    for mqc_line, main_line in zip(lines[6:], main_lines[1:]):
        assert mqc_line.split("\t", 1)[1] == main_line


def test_every_status_of_the_telomere_table_is_read_as_the_rules_say(tmp_path):
    """present, wrong_strand, other_end, not_terminal, absent, and a present array far from the end."""
    paths = make_inputs(tmp_path)
    rows = [
        # chrom_both: a telomere at both ends, the end array 2,500 bp inside the contig.
        telomere_row("chrom_both", len(CONTIGS["chrom_both"]), "start", "present", 90, 60),
        telomere_row("chrom_both", len(CONTIGS["chrom_both"]), "end", "present", 72, 2500),
        # arm_one: a telomere at the start, an array the wrong way round at the end.
        telomere_row("arm_one", len(CONTIGS["arm_one"]), "start", "present", 48, 40),
        telomere_row("arm_one", len(CONTIGS["arm_one"]), "end", "wrong_strand", 48, 1200),
        # plain: an array the wrong way round at the start, nothing at the end.
        telomere_row("plain", len(CONTIGS["plain"]), "start", "wrong_strand", 36, 0),
        telomere_row("plain", len(CONTIGS["plain"]), "end", "absent"),
        # lowdepth: one array seen from both ends of a short contig.
        telomere_row("lowdepth", len(CONTIGS["lowdepth"]), "start", "present", 30, 0),
        telomere_row("lowdepth", len(CONTIGS["lowdepth"]), "end", "other_end", 30, 4970),
        # fragment: an array too far inside the contig to count.
        telomere_row("fragment", len(CONTIGS["fragment"]), "start", "not_terminal", 30, 900),
        telomere_row("fragment", len(CONTIGS["fragment"]), "end", "absent"),
    ]
    write_tsv(paths.telomere_ends, TELOMERE_ENDS_COLUMNS, rows)
    contig_summary.main(full_argv(paths))
    table = read_table(paths.out)

    # present at both ends, one of them far from the end: a chromosome, and the reason says how far.
    assert table["chrom_both"]["verdict"] == "chromosome"
    assert table["chrom_both"]["telomere_end_offset_end"] == "2500"
    assert table["chrom_both"]["reason"].startswith(
        "telomere repeats at both ends (15 and 12 copies) (the end array is 2500 bp from the end)"
    )
    # wrong_strand is no telomere call: 0 copies, no offset, but the reason names the array.
    assert table["arm_one"]["verdict"] == "chromosome_arm"
    assert table["arm_one"]["telomere_end"] == "0"
    assert table["arm_one"]["telomere_end_offset_end"] == "NA"
    assert table["arm_one"]["reason"].startswith(
        "telomere repeats at the start only (8 copies) (array in the opposite orientation at the end)"
    )
    assert table["plain"]["verdict"] == "nuclear_contig"
    assert table["plain"]["telomere_start"] == "0"
    assert table["plain"]["reason"].endswith(
        "no telomere repeat at either end (array in the opposite orientation at the start)"
    )
    # other_end and not_terminal say nothing.
    assert table["lowdepth"]["verdict"] == "low_depth"
    assert table["lowdepth"]["telomere_start"] == "5"
    assert table["lowdepth"]["telomere_end"] == "0"
    assert table["lowdepth"]["reason"].endswith("telomere repeats at the start only (5 copies)")
    assert table["fragment"]["verdict"] == "nuclear_fragment"
    assert table["fragment"]["telomere_start"] == "0"
    assert table["fragment"]["reason"].endswith("no telomere repeat at either end")
    # A contig without any row has no telomere at either end.
    assert table["rdna"]["telomere_start"] == "0"
    assert table["rdna"]["telomere_end_offset_start"] == "NA"


def test_a_unit_of_none_leaves_every_telomere_column_empty(tmp_path, capsys):
    """No unit identified: counts 0, offsets NA, the unit column none, no contig a chromosome."""
    paths = make_inputs(tmp_path)
    write_telomere_files(paths.telomere_unit, paths.telomere_ends, "none", [])
    contig_summary.main(full_argv(paths))
    table = read_table(paths.out)

    for contig in CONTIGS:
        assert table[contig]["telomere_unit"] == "none", contig
        assert table[contig]["telomere_start"] == "0", contig
        assert table[contig]["telomere_end"] == "0", contig
        assert table[contig]["telomere_end_offset_start"] == "NA", contig
        assert table[contig]["telomere_end_offset_end"] == "NA", contig
    assert table["chrom_both"]["verdict"] == "nuclear_contig"
    assert table["arm_one"]["verdict"] == "nuclear_contig"
    assert "no telomere repeat unit identified, so no telomere at any of the 18 contig ends" in capsys.readouterr().err
    with open(paths.out_mqc) as mqc:
        description = mqc.read().splitlines()[2]
    assert "no telomere repeat unit could be identified for this sample" in description


def test_with_only_the_contigs_and_the_telomere_files_every_other_column_is_na(tmp_path, capsys):
    """Without the optional inputs the evidence columns read NA and the rules skip them."""
    paths = make_inputs(tmp_path)
    contig_summary.main([
        "--contigs", paths.contigs, "--telomere-unit-file", paths.telomere_unit,
        "--telomere-ends", paths.telomere_ends, "--sample", SAMPLE,
        "--out", paths.out, "--out-mqc", paths.out_mqc,
    ])
    table = read_table(paths.out)

    for column in ["depth", "depth_x_median", "assembler_depth", "circular", "rdna_units",
                   "busco_markers", "genes", "genus", "screen", "adapter", "annotated_as"]:
        assert table["plain"][column] == "NA", column
    assert table["chrom_both"]["verdict"] == "chromosome"
    assert "depth unknown" in table["chrom_both"]["reason"]
    assert "genes unknown" in table["chrom_both"]["reason"]
    assert table["arm_one"]["verdict"] == "chromosome_arm"
    # No depth, no genes, no rDNA count: nothing left to decide on.
    assert table["plain"]["verdict"] == "unclassified"
    assert table["plain"]["reason"] == "no telomere repeat at either end; depth unknown; genes unknown"
    assert table["mito"]["verdict"] == "unclassified"
    assert "median depth" not in capsys.readouterr().err


def test_a_bare_optional_flag_counts_as_not_given(tmp_path):
    """Snakemake writes an empty value as nothing, so '--blob-table' alone must mean 'no table'."""
    paths = make_inputs(tmp_path)
    argv = full_argv(paths, **{"--blob-table": None, "--assembly-evidence": None})
    argv += ["--blob-table", "--assembly-evidence"]      # both flags bare, at the very end
    contig_summary.main(argv)
    table = read_table(paths.out)
    assert table["plain"]["depth"] == "NA"
    assert table["mito"]["assembler_depth"] == "NA"
    assert table["mito"]["verdict"] == "mitochondrion"   # the screen report still names it
    # Without a depth nothing is "low" or "deep": the two gene-less contigs at an unknown
    # depth are nuclear fragments, which is what the gene evidence alone supports.
    assert table["lowdepth"]["verdict"] == "nuclear_fragment"
    assert table["repeat_deep"]["verdict"] == "nuclear_fragment"
    assert "depth unknown" in table["repeat_deep"]["reason"]


def test_the_adapter_column_shows_the_adapter_loci_and_changes_no_verdict(tmp_path):
    """The contig adapter screen fills the adapter column: none, at_end, internal or both."""
    paths = make_inputs(tmp_path)
    contig_summary.main(full_argv(paths))
    verdicts_without = {contig: row["verdict"] for contig, row in read_table(paths.out).items()}

    adapter_screen = os.path.join(tmp_path, "adapter_screen.tsv")
    write_tsv(adapter_screen, ["contig", "contig_length", "start", "end", "distance_to_end", "position",
                               "longest_hit", "identity", "adapters"], [
        ["chrom_both", 100, 1, 40, 0, "end", 40, "100.0", "ONT_native_barcode_NB24"],
        ["plain", 100, 50, 77, 23, "internal", 28, "100.0", "ONT_adapter_top"],
        ["arm_one", 100, 1, 30, 0, "end", 30, "100.0", "ONT_adapter_top"],
        ["arm_one", 100, 40, 67, 33, "internal", 28, "100.0", "ONT_adapter_top"],
    ])
    contig_summary.main(full_argv(paths, **{"--adapter-screen": adapter_screen}))
    table = read_table(paths.out)
    assert table["chrom_both"]["adapter"] == "at_end"
    assert table["plain"]["adapter"] == "internal"
    assert table["arm_one"]["adapter"] == "at_end+internal"
    assert table["mito"]["adapter"] == "none"
    assert {contig: row["verdict"] for contig, row in table.items()} == verdicts_without


def test_mito_fasta_alone_marks_the_mitochondrion(tmp_path):
    """Without the screen report, the mitochondrial FASTA still names the mitochondrion."""
    paths = make_inputs(tmp_path)
    contig_summary.main(full_argv(paths, **{"--mito-screen": None, "--decisions": None}))
    table = read_table(paths.out)
    assert table["mito"]["verdict"] == "mitochondrion"
    assert table["mito"]["reason"].startswith("listed in the mitochondrial FASTA")
    assert table["mito"]["screen"] == "mitochondrion"
    assert table["mito"]["annotated_as"] == "not annotated (mitochondrion)"
    assert table["removed"]["verdict"] != "contaminant"   # no decisions table: nothing is "removed"
    assert table["removed"]["screen"] == "NA"


def test_a_confirmed_mitochondrion_that_was_not_split_off_shows_its_scaffold(tmp_path):
    """Report mode of the mitochondrion screen: the contig is confirmed but not split off.

    It then reaches funannotate and is annotated like any other contig. The verdict stays
    mitochondrion (rule 1), but annotated_as must name its scaffold and genes must be
    filled, because the annotation files do hold it under that name.
    """
    paths = make_inputs(tmp_path)
    write_fasta(paths.mito_contigs, [])      # report mode: the split step wrote nothing
    annotated = [(scaffold, CONTIGS[contig]) for scaffold, contig in SCAFFOLDS.items()]
    write_fasta(paths.annotated_fasta, annotated + [("scaffold_7", CONTIGS["mito"])])
    genes = dict(GENES)
    genes["scaffold_7"] = 2
    write_gff3(paths.gff3, genes)
    contig_summary.main(full_argv(paths))
    table = read_table(paths.out)
    assert table["mito"]["verdict"] == "mitochondrion"
    assert table["mito"]["annotated_as"] == "scaffold_7"
    assert table["mito"]["genes"] == "2"
    assert table["mito"]["screen"] == "mitochondrion (direct_blast_plus_high_coverage)"


def test_an_rdna_contig_with_a_stray_gene_is_still_an_rdna_array(tmp_path):
    """Four rDNA units and one gene: an rDNA array; four units and four genes: not."""
    paths = make_inputs(tmp_path)
    # Give rdna a scaffold (scaffold_7) and one gene on it.
    write_fasta(paths.annotated_fasta, [(scaffold, CONTIGS[contig]) for scaffold, contig in SCAFFOLDS.items()] + [("scaffold_7", CONTIGS["rdna"])])
    genes = dict(GENES)
    genes["scaffold_7"] = 1
    write_gff3(paths.gff3, genes)
    contig_summary.main(full_argv(paths))
    table = read_table(paths.out)
    assert table["rdna"]["genes"] == "1"
    assert table["rdna"]["annotated_as"] == "scaffold_7"
    assert table["rdna"]["verdict"] == "rDNA_array"
    assert table["rdna"]["reason"].startswith("4 rDNA units; 1 gene")

    genes["scaffold_7"] = 4
    write_gff3(paths.gff3, genes)
    contig_summary.main(full_argv(paths))
    table = read_table(paths.out)
    assert table["rdna"]["verdict"] != "rDNA_array"
    assert table["rdna"]["verdict"] == "low_depth"   # 8.0 / 23.0 is below 0.5


def test_thresholds_change_the_calls(tmp_path):
    """A higher --min-telomere-copies turns arm_one into a plain nuclear contig; a higher --min-rdna-units frees rdna."""
    paths = make_inputs(tmp_path)
    contig_summary.main(full_argv(paths, **{"--min-telomere-copies": 10, "--min-rdna-units": 5}))
    table = read_table(paths.out)
    assert table["arm_one"]["verdict"] == "nuclear_contig"     # 8 copies are below 10
    assert table["chrom_both"]["verdict"] == "chromosome"      # 15 and 12 still pass
    assert table["rdna"]["verdict"] == "low_depth"             # 4 units are below 5


# ======================================================================
# THE JOIN TO FUNANNOTATE'S SCAFFOLD NAMES
# ======================================================================

def test_scaffold_mapping_by_identity_and_by_probes():
    """Identical sequences map directly; a polished copy maps by probes; a contained piece does not."""
    long_contig = random_dna(20000, 30)
    polished = list(long_contig)
    for position in (100, 7000, 9000, 15000, 19000):
        polished[position] = "A" if polished[position] != "A" else "C"
    polished = "".join(polished)
    assert polished != long_contig

    contigs = {
        "exact": random_dna(9000, 31),
        "polished": long_contig,
        "piece": long_contig[:3000],         # a near-duplicate funannotate clean would drop
        "gone": random_dna(4000, 32),         # no scaffold at all
    }
    scaffolds = {"scaffold_1": polished, "scaffold_2": contigs["exact"]}
    mapping = contig_summary.map_contigs_to_scaffolds(contigs, scaffolds)
    assert mapping == {"exact": "scaffold_2", "polished": "scaffold_1"}


def test_probe_search_needs_a_majority_and_a_matching_length():
    """A scaffold that holds only two of five probes, or has the wrong length, is not a match."""
    sequence = random_dna(20000, 33)
    # Keep only the middle 7,000 bases intact: probes 3 and part of the rest are destroyed.
    scrambled = random_dna(7000, 34) + sequence[7000:14000] + random_dna(6000, 35)
    assert contig_summary.find_scaffold_by_probes(sequence, {"s": scrambled}, set()) is None
    # The same sequence with a tail far longer than 5 percent of its length is rejected.
    too_long = sequence + random_dna(3000, 36)
    assert contig_summary.find_scaffold_by_probes(sequence, {"s": too_long}, set()) is None
    # A lightly polished copy of the right length is found; a claimed scaffold is skipped.
    polished = sequence[:5000] + "T" + sequence[5001:]
    assert contig_summary.find_scaffold_by_probes(sequence, {"s": polished}, set()) == "s"
    assert contig_summary.find_scaffold_by_probes(sequence, {"s": polished}, {"s"}) is None


# ======================================================================
# BAD INPUT GIVES A CLEAR ERROR, AND NO OUTPUT FILE
# ======================================================================

def test_depth_that_is_not_a_number_stops_the_script(tmp_path, capsys):
    """A BlobTools table with text in the depth column is refused with exit code 1."""
    paths = make_inputs(tmp_path)
    with open(paths.blob_table, "a") as table:
        table.write("broken\t1000\t0.4\t0\tabc\tEukaryota\t1.0\n")

    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths))

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("contig_summary.py: error:")
    assert "'abc' is not a number" in message
    assert not os.path.exists(paths.out)
    assert not os.path.exists(paths.out_mqc)


def test_table_without_the_needed_columns_stops_the_script(tmp_path, capsys):
    """A decisions table that lacks the action column is refused, naming the column."""
    paths = make_inputs(tmp_path)
    write_tsv(paths.decisions, ["contig", "assigned_genus", "reason"], [["plain", "Botrytis", "mode_off"]])

    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths))

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert "missing the column(s) action" in message
    assert not os.path.exists(paths.out)


def test_a_telomere_row_for_a_contig_the_draft_lacks_stops_the_script(tmp_path, capsys):
    """A telomere table that names a contig the draft does not hold belongs to another assembly."""
    paths = make_inputs(tmp_path)
    rows = made_up_telomere_rows() + [telomere_row("ghost", 5000, "start", "present", 30, 0)]
    write_tsv(paths.telomere_ends, TELOMERE_ENDS_COLUMNS, rows)

    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths))

    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("contig_summary.py: error:")
    assert "'ghost' is not in the draft assembly" in message
    assert not os.path.exists(paths.out)


def test_telomere_rows_next_to_a_unit_of_none_stop_the_script(tmp_path, capsys):
    """A unit of none comes with a header-only table; rows next to it mean mismatched files."""
    paths = make_inputs(tmp_path)
    write_telomere_files(paths.telomere_unit, paths.telomere_ends, "none", made_up_telomere_rows())

    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths))

    assert stop.value.code == 1
    assert "do not belong together" in capsys.readouterr().err
    assert not os.path.exists(paths.out)


def test_busco_folder_without_a_table_stops_the_script(tmp_path, capsys):
    """An empty BUSCO folder is a broken input, not 'no markers'."""
    paths = make_inputs(tmp_path)
    empty = tmp_path / "empty_busco"
    empty.mkdir()
    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths, **{"--busco-dir": str(empty)}))
    assert stop.value.code == 1
    assert "no run_*/full_table.tsv" in capsys.readouterr().err


def test_duplicate_contig_name_stops_the_script(tmp_path, capsys):
    """Two FASTA records with the same name make the join ambiguous: stop."""
    paths = make_inputs(tmp_path)
    write_fasta(paths.contigs, [("twice", "ACGT" * 300), ("twice", "TTTT" * 300)])
    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths))
    assert stop.value.code == 1
    assert "'twice' appears more than once" in capsys.readouterr().err


def test_missing_input_file_is_a_clear_error(tmp_path, capsys):
    """An input that does not exist gives one readable line, not a traceback."""
    paths = make_inputs(tmp_path)
    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths, **{"--its-windows": str(tmp_path / "nowhere.tsv")}))
    assert stop.value.code == 1
    message = capsys.readouterr().err
    assert message.startswith("contig_summary.py: error:")
    assert "nowhere.tsv" in message
    assert message.count("\n") == 1


@pytest.mark.parametrize(
    "option, value, expected_text",
    [
        ("--min-telomere-copies", "1", "at least 2"),
        ("--min-rdna-units", "0", "at least 1"),
    ],
)
def test_options_outside_their_limits_are_refused(tmp_path, capsys, option, value, expected_text):
    """The two tuning options are held to the same limits as the configuration (exit code 2)."""
    paths = make_inputs(tmp_path)
    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths, **{option: value}))
    assert stop.value.code == 2
    assert expected_text in capsys.readouterr().err


def test_the_telomere_files_are_required(tmp_path, capsys):
    """Without --telomere-unit-file or --telomere-ends argparse refuses the command line (exit code 2)."""
    paths = make_inputs(tmp_path)
    with pytest.raises(SystemExit) as stop:
        contig_summary.main(full_argv(paths, **{"--telomere-ends": None}))
    assert stop.value.code == 2
    assert "--telomere-ends" in capsys.readouterr().err


# ======================================================================
# THE SCRIPT AS SNAKEMAKE WILL CALL IT
# ======================================================================

def test_command_line_run_writes_the_two_files(tmp_path):
    """Exit code 0, the two files exist, the summary is on stderr and nothing on stdout."""
    paths = make_inputs(tmp_path)
    finished = subprocess.run(
        [sys.executable, SCRIPT_PATH] + full_argv(paths), capture_output=True, text=True,
    )
    assert finished.returncode == 0, finished.stderr
    assert finished.stdout == ""
    assert finished.stderr.startswith("contig_summary.py: sample SAMPLE1: 9 contigs")
    assert os.path.exists(paths.out)
    assert os.path.exists(paths.out_mqc)
    assert len(read_table(paths.out)) == len(CONTIGS)


def test_command_line_run_on_a_broken_table_exits_with_code_1(tmp_path):
    """The step must fail (non-zero exit) with a readable message and no traceback."""
    paths = make_inputs(tmp_path)
    write_tsv(paths.mito_screen, ["contig", "length"], [["mito", "6000"]])
    finished = subprocess.run(
        [sys.executable, SCRIPT_PATH] + full_argv(paths), capture_output=True, text=True,
    )
    assert finished.returncode == 1
    assert "contig_summary.py: error:" in finished.stderr
    assert "Traceback" not in finished.stderr
    assert not os.path.exists(paths.out)


def test_a_deep_contig_with_a_gene_call_or_two_is_a_repeat():
    """Rule 7 allows up to REPEAT_MAX_GENES predicted genes when there is no BUSCO marker."""
    def row(genes, busco, depth_x):
        # Only the keys the rules and the reason texts read.
        return {"is_mitochondrion": False, "is_removed": False, "rdna_units": 0,
                "telomere_start": 0, "telomere_end": 0,
                "telomere_end_offset_start": None, "telomere_end_offset_end": None,
                "telomere_status_start": None, "telomere_status_end": None,
                "depth_x_median": depth_x, "genes": genes, "busco_markers": busco}
    assert contig_summary.decide(row(2, 0, 3.4), 5, 3)[0] == "repeat"
    assert "gene call inside a repeat" in contig_summary.decide(row(1, 0, 3.4), 5, 3)[1]
    # Three genes, or a BUSCO marker, or a normal depth: not a repeat.
    assert contig_summary.decide(row(3, 0, 3.4), 5, 3)[0] == "unclassified"
    assert contig_summary.decide(row(1, 1, 3.4), 5, 3)[0] == "unclassified"
    assert contig_summary.decide(row(1, 0, 1.0), 5, 3)[0] == "nuclear_contig"
