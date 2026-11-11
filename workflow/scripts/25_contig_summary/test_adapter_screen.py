"""Tests for adapter_screen.py.

What the script does, in one sentence: it turns the BLAST hits of adapter sequences on the
draft contigs into adapter loci (at a contig end or inside a contig) and writes one row per
sample for the Adapters table of the MultiQC report.

The promises these tests pin down:
  * an end hit needs at least 20 bp and 90 % identity within 100 bp of a contig end; 19 bp,
    89 % identity, or 101 bp from the end do not count as an end hit;
  * an internal hit needs at least 28 bp (a whole ONT adapter); 27 bp does not count;
  * hits on the minus strand (BLAST writes their contig coordinates in reverse) are placed
    correctly;
  * overlapping hits of several adapter records form one locus, which names them all;
  * the MultiQC row reports the read-level numbers of a trimmed ONT sample and a verdict.

Every test writes a small FASTA and a made-up BLAST table to pytest's temporary folder.

Run: pytest workflow/scripts/25_contig_summary/test_adapter_screen.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import adapter_screen as screen

CONTIG_LENGTH = 100_000


def write_contigs(tmp_path, names=("ctg1", "ctg2")):
    # A FASTA of contigs of CONTIG_LENGTH bases each (the screen needs only their names and
    # lengths); returns its path, for --contigs.
    path = tmp_path / "contigs.fasta"
    with open(path, "w") as out:
        for name in names:
            out.write(f">{name}\n{'A' * CONTIG_LENGTH}\n")
    return str(path)


def blast_row(adapter, contig, identity, length, sstart, send, query_length=40):
    # One line of the BLAST table in the column order of rule blast_adapters.
    return "\t".join(str(value) for value in [
        adapter, contig, identity, length, 1, length, sstart, send, 1e-5, 40, query_length, CONTIG_LENGTH,
    ]) + "\n"


def run_screen(tmp_path, rows, choice=None, check=None):
    # Run the script as rule screen_adapters does: on the BLAST rows given, the two test
    # contigs and, when given, a choice and a check file of the read-level steps
    # ("name<TAB>value" lines). Returns the loci table as a list of dicts (one per locus)
    # and the MultiQC file as a list of lines.
    blast = tmp_path / "blast.tsv"
    blast.write_text("".join(rows))
    arguments = [
        "--blast", str(blast), "--contigs", write_contigs(tmp_path), "--sample", "S1",
        "--out", str(tmp_path / "loci.tsv"), "--out-mqc", str(tmp_path / "mqc.tsv"),
    ]
    if choice:
        choice_path = tmp_path / "choice.tsv"
        choice_path.write_text("".join(f"{k}\t{v}\n" for k, v in choice.items()))
        arguments += ["--read-choice", str(choice_path)]
    if check:
        check_path = tmp_path / "check.tsv"
        check_path.write_text("".join(f"{k}\t{v}\n" for k, v in check.items()))
        arguments += ["--read-check", str(check_path)]
    screen.main(arguments)
    lines = (tmp_path / "loci.tsv").read_text().strip().split("\n")
    header = lines[0].split("\t")
    loci = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    mqc = (tmp_path / "mqc.tsv").read_text().strip().split("\n")
    return loci, mqc


# ── End hits ──

def test_an_end_hit_of_20_bp_at_90_percent_counts(tmp_path):
    loci, _ = run_screen(tmp_path, [blast_row("NB24", "ctg1", 90.0, 20, 51, 70)])
    assert len(loci) == 1
    assert loci[0]["position"] == "end"
    assert loci[0]["distance_to_end"] == "50"


def test_shorter_or_less_similar_end_hits_do_not_count(tmp_path):
    # 19 bp at the very end of the contig; 30 bp, 50 bp from the end, but 89 % identical.
    rows = [
        blast_row("NB24", "ctg1", 95.0, 19, 1, 19),
        blast_row("NB24", "ctg1", 89.0, 30, 51, 80),
    ]
    loci, _ = run_screen(tmp_path, rows)
    assert loci == []


def test_a_hit_101_bp_from_the_end_is_judged_as_internal(tmp_path):
    # Starts at base 102: 101 bases lie before it, so it is an internal hit and 25 bp is too short.
    loci, _ = run_screen(tmp_path, [blast_row("NB24", "ctg1", 100.0, 25, 102, 126)])
    assert loci == []


def test_a_minus_strand_hit_at_the_right_end_is_placed_correctly(tmp_path):
    loci, _ = run_screen(tmp_path, [blast_row("ONT_adapter_top", "ctg2", 100.0, 28, CONTIG_LENGTH, CONTIG_LENGTH - 27)])
    assert loci[0]["start"] == str(CONTIG_LENGTH - 27)
    assert loci[0]["end"] == str(CONTIG_LENGTH)
    assert loci[0]["distance_to_end"] == "0"


# ── Internal hits ──

def test_an_internal_hit_needs_a_whole_adapter(tmp_path):
    rows = [
        blast_row("ONT_adapter_top", "ctg1", 100.0, 28, 50_001, 50_028),
        blast_row("ONT_adapter_top", "ctg2", 100.0, 27, 50_001, 50_027),
    ]
    loci, _ = run_screen(tmp_path, rows)
    assert [(locus["contig"], locus["position"]) for locus in loci] == [("ctg1", "internal")]


# ── Merging ──

def test_overlapping_hits_form_one_locus_naming_every_adapter(tmp_path):
    rows = [
        blast_row("ONT_adapter_top", "ctg1", 100.0, 28, 5, 32),
        blast_row("ONT_native_barcode_NB24", "ctg1", 100.0, 40, 33, 72),
        blast_row("S1_trim_1", "ctg1", 100.0, 68, 5, 72),
    ]
    loci, _ = run_screen(tmp_path, rows)
    assert len(loci) == 1
    assert loci[0]["start"] == "5" and loci[0]["end"] == "72"
    assert loci[0]["longest_hit"] == "68"
    assert loci[0]["adapters"] == "ONT_adapter_top,ONT_native_barcode_NB24,S1_trim_1"


# ── The MultiQC row ──

def test_the_mqc_row_reports_the_reads_and_a_clean_verdict(tmp_path):
    choice = {"setting": "detect", "percent_with_adapter": "38.56", "published_barcodes": "NB24 87.0 %",
              "summary": "adapter + native barcode NB24 (68 bp)"}
    check = {"percent_with_adapter": "0.12"}
    loci, mqc = run_screen(tmp_path, [], choice=choice, check=check)
    assert loci == []
    assert mqc[0] == "# id: adapter_screen"
    assert mqc[-1].split("\t") == [
        "S1", "38.56", "0.12", "NB24 87.0 %", "adapter + native barcode NB24 (68 bp)", "0", "0",
        "none of the listed adapters found",
    ]


def test_the_mqc_row_of_an_untrimmed_sample_counts_the_loci(tmp_path):
    rows = [
        blast_row("NB24", "ctg1", 100.0, 40, 1, 40),
        blast_row("NB24", "ctg2", 100.0, 40, CONTIG_LENGTH - 39, CONTIG_LENGTH),
        blast_row("ONT_adapter_top", "ctg2", 100.0, 28, 40_000, 40_027),
    ]
    _, mqc = run_screen(tmp_path, rows)
    assert mqc[-1].split("\t") == [
        "S1", "NA", "NA", "NA", "reads not trimmed by FunFlux", "2", "1",
        "adapter at 2 contig ends; 1 adapter match inside contigs to check",
    ]


def test_two_loci_at_one_contig_end_count_as_one_end(tmp_path):
    # Two separate loci near the start of ctg1 (more than 10 bp apart): one contig end with adapter.
    rows = [
        blast_row("NB24", "ctg1", 100.0, 40, 1, 40),
        blast_row("ONT_adapter_top", "ctg1", 100.0, 28, 61, 88),
    ]
    loci, mqc = run_screen(tmp_path, rows)
    assert len(loci) == 2
    assert mqc[-1].split("\t")[5] == "1"
