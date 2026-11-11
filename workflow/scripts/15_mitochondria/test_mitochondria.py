"""Tests for mitochondria.py, the mitochondrial screen of FunFlux v2.

Everything here runs on tiny hand-written files, with no assembler, no BLAST
database and no third-party package (only pytest). The tests are grouped like
the script:

  1. classify()          every branch of the decision function, one case each,
                         including the "very strong direct BLAST" cut-off with a
                         non-default value
  2. the evidence readers  Flye assembly_info.txt and hifiasm GFA, including the
                         hifiasm contig-name suffix (c = circular, l = linear)
  3. read_blast()        picks the best direct organelle hit and keeps both lengths
  4. the concatemer flag  Y for a contig that is two genome copies end to end,
                         N otherwise, and proof that it changes nothing else
  5. split               writes the right two FASTA files, fails clearly when a
                         confirmed ID is missing
  6. the command line    what the Snakemake rules actually run, including the
                         older flag name --min-flye-coverage-ratio

Why these tests matter: the script decides which contig is taken out of the
genome that goes into the nuclear gene predictor. A wrong call either sends the
mitochondrion through the wrong gene predictor or removes a nuclear contig from
the annotation, and neither is visible in the final files unless someone
looks. The numbers used below (a 30,000 bp reference, a 54,000 bp contig,
depths 20 and 150) are made up but of a realistic size.

Run:  python -m pytest workflow/scripts/15_mitochondria -q
"""

import csv
import os
import subprocess
import sys

import pytest

# Make "import mitochondria" find the script that sits next to this test file,
# whichever folder pytest is started from.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_PATH = os.path.join(SCRIPT_DIR, "mitochondria.py")
sys.path.insert(0, SCRIPT_DIR)

import mitochondria as mito  # noqa: E402  (has to come after the path line above)


# =============================================================================
# Helpers shared by the tests
# =============================================================================

# The report columns of version 0.2 of the screen, in their order.
# FunFlux v2 must keep them exactly where they were and add its own at the end.
OLD_REPORT_COLUMNS = [
    "contig", "length", "classification", "reason", "best_mito_qcovs",
    "direct_mito_hits", "best_mito_title", "assembler_source",
    "assembler_circular", "assembler_coverage", "assembler_coverage_ratio",
    "assembler_repeat", "assembler_alt_group",
]
NEW_REPORT_COLUMNS = ["best_mito_slen", "length_ratio", "concatemer_flag", "best_mito_kingdom"]

# Reference titles as they appear in the last BLAST column (stitle).
MITO_TITLE = "Example fungus mitochondrion, complete sequence"
NUCLEAR_TITLE = "Fusarium oxysporum chromosome 3, complete sequence"

# Typical read depth of a nuclear contig, and of the mitochondrial contig, in the
# classify() examples (150 / 20 = 7.5 times the median).
MEDIAN_DEPTH = 20.0
DEEP = 150.0

# A small taxonomy in the layout of NCBI's nodes.dmp (taxonomy ID, parent ID, rank,
# separated by a tab, a bar and a tab), enough for the species check: the root (1),
# Eukaryota (2759); Fungi (4751) with Ascomycota (4890) and a fungus (5062); Metazoa
# (33208) with Homo sapiens (9606); Viridiplantae (33090) with a plant (3702); an
# environmental record filed under Eukaryota only (100272, "uncultured eukaryote"); and
# an oomycete (4787, Phytophthora infestans) under the stramenopiles (33634), which NCBI
# does not place in Fungi.
TAXONOMY_LINES = [
    "1\t|\t1\t|\tno rank\t|",
    "2759\t|\t1\t|\tsuperkingdom\t|",
    "4751\t|\t2759\t|\tkingdom\t|",
    "4890\t|\t4751\t|\tphylum\t|",
    "5062\t|\t4890\t|\tspecies\t|",
    "33208\t|\t2759\t|\tkingdom\t|",
    "9606\t|\t33208\t|\tspecies\t|",
    "33090\t|\t2759\t|\tkingdom\t|",
    "3702\t|\t33090\t|\tspecies\t|",
    "100272\t|\t2759\t|\tspecies\t|",
    "33634\t|\t2759\t|\tclade\t|",
    "4787\t|\t33634\t|\tspecies\t|",
]
HUMAN_MITO_TITLE = "Homo sapiens isolate X1 mitochondrion, complete genome"


def write_text(path, text):
    """Write a small text file and give back its path as a string."""
    with open(str(path), "w", encoding="utf-8") as handle:
        handle.write(text)
    return str(path)


def blast_line(contig, title, qcovs, qlen=1000, slen=1000, taxids="5062"):
    """Return one line of the workflow's BLAST table (15 tab-separated columns).

    The columns are: qseqid staxids bitscore pident evalue length qlen slen
    qcovs qcovhsp sskingdoms scomnames sscinames sblastnames stitle.
    The script reads the contig name, the taxonomy IDs (for the species check),
    qlen, slen, qcovs and the title, so the other columns are filled with
    plausible dummy values. 5062 is a fungus in the small taxonomy of TAXONOMY_LINES.
    """
    fields = [
        contig, taxids, "30000", "99.9", "0.0", "28000", str(qlen), str(slen),
        str(qcovs), str(qcovs), "Eukaryota", "fungal sp.", "fungal sp.",
        "fungi", title,
    ]
    return "\t".join(fields)


def write_blast(path, lines):
    """Write a BLAST table made of the given lines."""
    return write_text(path, "\n".join(lines) + "\n")


def write_contigs(path, contig_lengths):
    """Write a FASTA with one contig per (name, length) pair; sequences are all A."""
    lines = []
    for name, length in contig_lengths:
        lines.append(">" + name)
        lines.append("A" * length)
    return write_text(path, "\n".join(lines) + "\n")


def read_table(path):
    """Read a tab-separated report into (column names, list of row dictionaries)."""
    with open(str(path), "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        column_names = reader.fieldnames
    return column_names, rows


def rows_by_contig(path):
    """Read a report into {contig name: row dictionary}."""
    column_names, rows = read_table(path)
    return {row["contig"]: row for row in rows}


def new_column_values(row):
    """The three concatemer values of one report row, as (best_mito_slen, length_ratio, flag)."""
    return (row["best_mito_slen"], row["length_ratio"], row["concatemer_flag"])


def run_cli(*arguments):
    """Run the script the way a Snakemake rule does and return the finished process.

    This is "python mitochondria.py <arguments>" in a separate process, so the
    exit status and the ERROR line on stderr can be checked as well.
    """
    command = [sys.executable, SCRIPT_PATH] + list(arguments)
    return subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )


def parse(*arguments):
    """Parse command-line arguments with the script's own parser (no process started)."""
    return mito.build_parser().parse_args(list(arguments))


def make_hit(qcov):
    """A stand-in for what read_blast() returns for one contig (classify reads best_qcov only)."""
    return {
        "best_qcov": qcov,
        "best_title": MITO_TITLE,
        "best_qlen": 29980,
        "best_slen": 30000,
        "direct_hits": 1,
    }


def make_assembler(coverage, circular=False, repeat=False, alt_group="*"):
    """A stand-in for one contig's assembler record (depth, circular, repeat, alt_group)."""
    return {
        "source": "flye",
        "coverage": coverage,
        "circular": circular,
        "repeat": repeat,
        "alt_group": alt_group,
    }


def run_classify(evidence, assembler_record, median_depth=MEDIAN_DEPTH,
                 minimum_qcov=35.0, minimum_ratio=3.0, mode="separate_confirmed",
                 very_strong_qcov=70.0):
    """Call classify() with the defaults of the screen subcommand unless told otherwise."""
    return mito.classify(
        evidence, assembler_record, median_depth, minimum_qcov, minimum_ratio,
        mode, very_strong_qcov,
    )


def test_tests_load_the_script_that_sits_next_to_them():
    # A guard: if this fails, the tests are checking some other copy of the script.
    assert os.path.samefile(mito.__file__, SCRIPT_PATH)


# =============================================================================
# 1. classify(): the decision function, every branch once
# =============================================================================
# classify() gives one contig a class and a reason code. Only the class
# "confirmed_mitochondrial" leads to the contig being split off the nuclear
# assembly; the reason code is what a person reads in the report to see why.
# The cases below follow the order of the branches in the function.

def test_classify_mode_off_does_not_screen():
    # With screening switched off nothing is examined, whatever the evidence says.
    result = run_classify(make_hit(100.0), make_assembler(DEEP, circular=True), mode="off")
    assert result == ("not_screened", "screening_disabled", None)


# --- No assembler evidence (short-read runs): BLAST is the only information ---

def test_classify_blast_only_strong_hit_is_confirmed():
    # A realistic short-read case: the best direct hit covers 52 percent of the
    # mitochondrial contig, against a threshold of 30.
    result = run_classify(make_hit(52.0), None, minimum_qcov=30.0)
    assert result == ("confirmed_mitochondrial", "strong_direct_blast", None)


def test_classify_blast_only_hit_exactly_at_the_threshold_counts():
    # The threshold is "at least", not "more than".
    result = run_classify(make_hit(35.0), None, minimum_qcov=35.0)
    assert result == ("confirmed_mitochondrial", "strong_direct_blast", None)


def test_classify_blast_only_weak_hit_is_only_a_candidate():
    # A direct organelle hit that covers only 20 percent of the contig is not enough.
    result = run_classify(make_hit(20.0), None, minimum_qcov=35.0)
    assert result == ("mitochondrial_candidate", "weak_direct_blast", None)


def test_classify_blast_only_without_a_hit_is_not_mitochondrial():
    result = run_classify(None, None)
    assert result == ("not_mitochondrial", "no_mitochondrial_evidence", None)


# --- With assembler evidence (long-read and hybrid runs) ---

def test_classify_direct_hit_plus_high_coverage_is_confirmed():
    # Hit covers 50 percent (above the threshold 35, below "very strong" 70) and
    # the contig is 150 / 20 = 7.5 times deeper than the median: confirmed
    # because BLAST and read depth agree.
    result = run_classify(make_hit(50.0), make_assembler(DEEP))
    assert result[:2] == ("confirmed_mitochondrial", "direct_blast_plus_high_coverage")
    assert result[2] == pytest.approx(DEEP / MEDIAN_DEPTH)


def test_classify_direct_hit_plus_clean_circle_is_confirmed():
    # Normal depth, but the assembler closed the contig into a circle (and it is
    # neither a graph repeat nor an alternative copy): BLAST and topology agree.
    result = run_classify(make_hit(50.0), make_assembler(20.0, circular=True))
    assert result[:2] == ("confirmed_mitochondrial", "direct_blast_plus_clean_circularity")
    assert result[2] == pytest.approx(1.0)


def test_classify_very_strong_hit_is_confirmed_without_any_support():
    # 80 percent of the contig matches an organelle genome; the assembler gives
    # no support (normal depth, linear) but a hit this strong is decisive.
    result = run_classify(make_hit(80.0), make_assembler(20.0))
    assert result[:2] == ("confirmed_mitochondrial", "very_strong_direct_blast")
    assert result[2] == pytest.approx(1.0)


def test_classify_very_strong_cutoff_can_be_changed():
    # The new --very-strong-qcov option. A contig whose hit covers 55 percent and
    # that has no assembler support is judged differently at different cut-offs.
    no_support = make_assembler(20.0)
    hit = make_hit(55.0)

    # Default cut-off 70: 55 is not very strong, so it stays a candidate.
    default_result = run_classify(hit, no_support)
    assert default_result[:2] == ("mitochondrial_candidate", "direct_blast_without_support")

    # Lowered to 50: 55 now counts as very strong.
    lowered_result = run_classify(hit, no_support, very_strong_qcov=50.0)
    assert lowered_result[:2] == ("confirmed_mitochondrial", "very_strong_direct_blast")

    # Raised to 60: 55 is below it again.
    raised_result = run_classify(hit, no_support, very_strong_qcov=60.0)
    assert raised_result[:2] == ("mitochondrial_candidate", "direct_blast_without_support")

    # A hit exactly at the cut-off counts ("at least").
    boundary_result = run_classify(make_hit(50.0), no_support, very_strong_qcov=50.0)
    assert boundary_result[:2] == ("confirmed_mitochondrial", "very_strong_direct_blast")


def test_classify_very_strong_cutoff_default_is_70():
    # With no explicit value the cut-off must still be the old fixed 70 percent.
    assert mito.DEFAULT_VERY_STRONG_QCOV == 70.0
    no_support = make_assembler(20.0)
    assert run_classify(make_hit(70.0), no_support)[1] == "very_strong_direct_blast"
    assert run_classify(make_hit(69.9), no_support)[1] == "direct_blast_without_support"
    # Calling classify() without the new argument at all gives the same answer as 70.
    old_style_call = mito.classify(make_hit(70.0), no_support, MEDIAN_DEPTH, 35.0, 3.0, "separate_confirmed")
    assert old_style_call[:2] == ("confirmed_mitochondrial", "very_strong_direct_blast")


def test_classify_very_strong_cutoff_can_be_raised_above_70():
    # Raising the cut-off also works: 75 percent is no longer "very strong" at 80.
    no_support = make_assembler(20.0)
    assert run_classify(make_hit(75.0), no_support)[1] == "very_strong_direct_blast"
    assert run_classify(make_hit(75.0), no_support, very_strong_qcov=80.0)[1] == "direct_blast_without_support"


def test_classify_direct_hit_without_support_is_a_candidate():
    # A direct hit is there, but it is not very strong and neither depth nor
    # circularity backs it up: worth a look, not enough to split the contig off.
    result = run_classify(make_hit(50.0), make_assembler(20.0))
    assert result[:2] == ("mitochondrial_candidate", "direct_blast_without_support")
    assert result[2] == pytest.approx(1.0)


def test_classify_circular_and_deep_without_blast_is_a_candidate():
    # Looks like a mitochondrion (closed circle, many times the normal depth) but
    # BLAST found no organelle reference: only a candidate.
    result = run_classify(None, make_assembler(DEEP, circular=True))
    assert result[:2] == ("mitochondrial_candidate", "circular_and_high_coverage_without_blast")
    assert result[2] == pytest.approx(DEEP / MEDIAN_DEPTH)


def test_classify_circularity_alone_never_confirms():
    result = run_classify(None, make_assembler(20.0, circular=True))
    assert result[:2] == ("circular_unclassified", "circularity_only")


def test_classify_high_coverage_alone_never_confirms():
    result = run_classify(None, make_assembler(DEEP))
    assert result[:2] == ("high_coverage_unclassified", "coverage_only")
    assert result[2] == pytest.approx(DEEP / MEDIAN_DEPTH)


def test_classify_unremarkable_contig_is_not_mitochondrial():
    result = run_classify(None, make_assembler(20.0))
    assert result[:2] == ("not_mitochondrial", "no_mitochondrial_evidence")
    assert result[2] == pytest.approx(1.0)


# --- Edge cases around the branches ---

def test_classify_missing_depth_gives_no_ratio():
    # A hifiasm GFA without rd:i tags has no depth. Then there is no ratio (None),
    # depth cannot support a hit, and only a very strong hit can still confirm.
    no_depth = make_assembler(None)
    middle_hit = run_classify(make_hit(50.0), no_depth)
    assert middle_hit == ("mitochondrial_candidate", "direct_blast_without_support", None)
    strong_hit = run_classify(make_hit(80.0), no_depth)
    assert strong_hit == ("confirmed_mitochondrial", "very_strong_direct_blast", None)
    # The same happens when the median depth of the whole assembly is unknown.
    unknown_median = run_classify(make_hit(50.0), make_assembler(DEEP), median_depth=None)
    assert unknown_median == ("mitochondrial_candidate", "direct_blast_without_support", None)


def test_classify_a_circle_that_is_a_repeat_or_an_alternative_is_not_clean():
    # "Clean" circularity needs: circular, not a graph repeat, not an alternative
    # copy of another contig. Otherwise a middle-strength hit gets no support.
    hit = make_hit(50.0)
    circular_repeat = make_assembler(20.0, circular=True, repeat=True)
    circular_alternative = make_assembler(20.0, circular=True, alt_group="1")
    for record in (circular_repeat, circular_alternative):
        result = run_classify(hit, record)
        assert result[:2] == ("mitochondrial_candidate", "direct_blast_without_support")


def test_classify_depth_ratio_exactly_at_the_threshold_counts_as_high():
    # 60 / 20 = 3.0, exactly the default ratio threshold ("at least").
    result = run_classify(make_hit(50.0), make_assembler(60.0), minimum_ratio=3.0)
    assert result[:2] == ("confirmed_mitochondrial", "direct_blast_plus_high_coverage")


def test_classify_report_mode_labels_contigs_like_separate_confirmed():
    # The two modes differ only in what screen() selects afterwards, not in labels.
    hit = make_hit(50.0)
    record = make_assembler(DEEP)
    assert run_classify(hit, record, mode="report") == run_classify(hit, record, mode="separate_confirmed")


def test_classify_recording_the_circle_does_not_change_a_long_read_call():
    # A long-read case: a hifiasm contig ending in c, hit covering 100 percent,
    # depth 150 against a median of 20. hifiasm used to be read as "linear"; now
    # the c in the name makes it circular. Depth already supports the hit and
    # gives the reason, so the call and its reason text must be the same either way.
    as_linear = run_classify(make_hit(100.0), make_assembler(DEEP, circular=False), minimum_qcov=30.0)
    as_circular = run_classify(make_hit(100.0), make_assembler(DEEP, circular=True), minimum_qcov=30.0)
    assert as_linear == as_circular
    assert as_circular[:2] == ("confirmed_mitochondrial", "direct_blast_plus_high_coverage")


# =============================================================================
# 2. The two evidence readers (Flye and hifiasm)
# =============================================================================
# Both turn an assembler's own output into the same small per-contig table, so
# that "screen" does not need to know which assembler made the contigs.

# A tiny Flye assembly_info.txt: four contigs with the real column layout.
#   contig_2 is a closed circle (the mitochondrion), contig_3 a closed circle
#   that Flye also marks as a graph repeat, contig_4 belongs to alternative group 1.
FLYE_INFO_TEXT = (
    "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
    "contig_1\t4400000\t20\tN\tN\t1\t*\t1\n"
    "contig_2\t29980\t150\tY\tN\t1\t*\t2\n"
    "contig_3\t7000\t800\tY\tY\t40\t*\t3\n"
    "contig_4\t120000\t11\tN\tN\t1\t1\t4\n"
)

# A tiny hifiasm primary GFA. The last letter of a contig name is hifiasm's own
# verdict: c = circular, l = linear. S lines are contigs; the A (read placement)
# and L (link) lines and the H (header) line must be ignored.
GFA_TEXT = (
    "H\tVN:Z:1.0\n"
    "S\tptg000001l\tACGTACGTAC\tLN:i:10\trd:i:20\n"
    "S\tptg000002c\tACGTACGTACGT\tLN:i:12\trd:i:150\n"
    "S\tptg000003c\tACGT\tLN:i:4\n"
    "S\th1tg000004l\tACGT\tLN:i:4\trd:i:25\n"
    "S\tctg000005l\tACGT\tLN:i:4\trd:i:30\n"
    "A\tptg000001l\t0\t-\tread1\t0\t100\tid:i:1\tHG:A:a\n"
    "L\tptg000001l\t+\tptg000002c\t+\t0M\n"
)


def test_read_flye_info_records_circularity_depth_repeat_and_group(tmp_path):
    info_path = write_text(tmp_path / "assembly_info.txt", FLYE_INFO_TEXT)
    records = mito.read_flye_info(info_path)

    assert list(records) == ["contig_1", "contig_2", "contig_3", "contig_4"]

    linear_nuclear = records["contig_1"]
    assert linear_nuclear["source"] == "flye"
    assert linear_nuclear["circular"] is False
    assert linear_nuclear["coverage"] == 20.0

    closed_circle = records["contig_2"]
    assert closed_circle["circular"] is True
    assert closed_circle["repeat"] is False
    assert closed_circle["alt_group"] == "*"
    assert closed_circle["coverage"] == 150.0

    circular_repeat = records["contig_3"]
    assert circular_repeat["circular"] is True
    assert circular_repeat["repeat"] is True

    alternative_copy = records["contig_4"]
    assert alternative_copy["circular"] is False
    assert alternative_copy["alt_group"] == "1"


def test_read_flye_info_without_contigs_fails_clearly(tmp_path):
    header_only = write_text(tmp_path / "assembly_info.txt", "#seq_name\tlength\tcov.\tcirc.\n")
    with pytest.raises(ValueError, match="No Flye records found"):
        mito.read_flye_info(header_only)


def test_evidence_flye_command_writes_the_common_table(tmp_path):
    info_path = write_text(tmp_path / "assembly_info.txt", FLYE_INFO_TEXT)
    output_path = tmp_path / "assembly_evidence.tsv"

    result = run_cli("evidence-flye", "--flye-info", info_path, "--output", str(output_path))
    assert result.returncode == 0, result.stderr
    assert "wrote 4 Flye contig record(s)" in result.stdout

    expected = (
        "contig\tsource\tcircular\tcoverage\trepeat\talt_group\n"
        "contig_1\tflye\tN\t20.000000\tN\t*\n"
        "contig_2\tflye\tY\t150.000000\tN\t*\n"
        "contig_3\tflye\tY\t800.000000\tY\t*\n"
        "contig_4\tflye\tN\t11.000000\tN\t1\n"
    )
    assert output_path.read_text() == expected


def test_read_hifiasm_gfa_c_suffix_means_circular(tmp_path):
    # CHANGED in FunFlux v2: hifiasm marks a circular contig by a final c.
    gfa_path = write_text(tmp_path / "hifiasm.p_ctg.gfa", GFA_TEXT)
    records = mito.read_hifiasm_gfa(gfa_path)

    assert records["ptg000002c"]["circular"] is True      # ends in c
    assert records["ptg000003c"]["circular"] is True      # ends in c (no depth tag)
    assert records["ptg000001l"]["circular"] is False     # ends in l
    assert records["h1tg000004l"]["circular"] is False    # haplotype-style name, ends in l


def test_read_hifiasm_gfa_only_the_last_letter_counts(tmp_path):
    # ctg000005l contains a c, but it is the last letter that carries the verdict.
    gfa_path = write_text(tmp_path / "hifiasm.p_ctg.gfa", GFA_TEXT)
    records = mito.read_hifiasm_gfa(gfa_path)
    assert records["ctg000005l"]["circular"] is False


def test_read_hifiasm_gfa_reads_depth_and_ignores_other_line_types(tmp_path):
    gfa_path = write_text(tmp_path / "hifiasm.p_ctg.gfa", GFA_TEXT)
    records = mito.read_hifiasm_gfa(gfa_path)

    # Only the five S lines became records; the H, A and L lines were skipped.
    assert list(records) == ["ptg000001l", "ptg000002c", "ptg000003c", "h1tg000004l", "ctg000005l"]
    # Depth comes from the rd:i tag; a segment without the tag has no depth.
    assert records["ptg000001l"]["coverage"] == 20.0
    assert records["ptg000002c"]["coverage"] == 150.0
    assert records["ptg000003c"]["coverage"] is None
    # hifiasm gives no repeat or alternative-group information.
    for record in records.values():
        assert record["source"] == "hifiasm"
        assert record["repeat"] is False
        assert record["alt_group"] == "*"


def test_read_hifiasm_gfa_without_segments_fails_clearly(tmp_path):
    no_segments = write_text(tmp_path / "hifiasm.p_ctg.gfa", "H\tVN:Z:1.0\nL\ta\t+\tb\t+\t0M\n")
    with pytest.raises(ValueError, match="No Hifiasm segment records found"):
        mito.read_hifiasm_gfa(no_segments)


def test_evidence_hifiasm_command_writes_the_common_table(tmp_path):
    gfa_path = write_text(tmp_path / "hifiasm.p_ctg.gfa", GFA_TEXT)
    output_path = tmp_path / "assembly_evidence.tsv"

    result = run_cli("evidence-hifiasm", "--gfa", gfa_path, "--output", str(output_path))
    assert result.returncode == 0, result.stderr
    assert "wrote 5 Hifiasm contig record(s); 4 include rd:i read-depth tags." in result.stdout

    expected = (
        "contig\tsource\tcircular\tcoverage\trepeat\talt_group\n"
        "ptg000001l\thifiasm\tN\t20.000000\tN\t*\n"
        "ptg000002c\thifiasm\tY\t150.000000\tN\t*\n"
        "ptg000003c\thifiasm\tY\t\tN\t*\n"
        "h1tg000004l\thifiasm\tN\t25.000000\tN\t*\n"
        "ctg000005l\thifiasm\tN\t30.000000\tN\t*\n"
    )
    assert output_path.read_text() == expected


def test_evidence_table_round_trip_gives_records_and_median_depth(tmp_path):
    # The table written by evidence-hifiasm is what "screen" reads back.
    gfa_path = write_text(tmp_path / "hifiasm.p_ctg.gfa", GFA_TEXT)
    table_path = str(tmp_path / "assembly_evidence.tsv")
    mito.write_assembly_evidence(table_path, mito.read_hifiasm_gfa(gfa_path))

    records, median_depth = mito.read_assembly_evidence(table_path)
    assert records["ptg000002c"]["circular"] is True
    assert records["ptg000001l"]["circular"] is False
    assert records["ptg000003c"]["coverage"] is None      # an empty cell stays "no depth"
    # Median of the four known depths 20, 150, 25, 30 is (25 + 30) / 2.
    assert median_depth == pytest.approx(27.5)


def test_evidence_table_with_a_missing_column_fails_clearly(tmp_path):
    table_path = write_text(tmp_path / "assembly_evidence.tsv", "contig\tsource\tcircular\nptg1\thifiasm\tY\n")
    with pytest.raises(ValueError, match="missing required column"):
        mito.read_assembly_evidence(table_path)


# =============================================================================
# 3. read_blast(): the best direct organelle hit per contig
# =============================================================================
# Only BLAST hits whose reference title says "this IS a mitochondrial genome"
# count. The best one (highest qcovs, the percentage of the contig covered by
# that reference) is kept, together with the contig length (qlen) and the
# reference length (slen) of that same hit.

def test_read_blast_picks_the_best_direct_hit_and_keeps_both_lengths(tmp_path):
    blast_path = write_blast(tmp_path / "blastout", [
        # contig_A has three direct hits; the second has the highest qcovs.
        blast_line("contig_A", "Fungus one mitochondrion, complete genome", 40, qlen=29980, slen=30000),
        blast_line("contig_A", "Fungus two mitochondrion, complete genome", 100, qlen=29980, slen=27900),
        blast_line("contig_A", "Fungus three mitochondrial DNA, complete genome", 100, qlen=29980, slen=29000),
        # contig_B only matches a nuclear chromosome.
        blast_line("contig_B", NUCLEAR_TITLE, 95, qlen=5000, slen=40000000),
        # contig_C matches a nuclear gene that merely has "mitochondrial" in its name.
        blast_line("contig_C", "Fusarium graminearum mitochondrial ribosomal protein L5 (MRPL5) mRNA", 90),
        # contig_D has one weak hit that uses the "mitochondrial genome" wording.
        blast_line("contig_D", "Neurospora crassa mitochondrial genome", 10, qlen=6000, slen=64840),
    ])
    evidence = mito.read_blast(blast_path)

    # Only the contigs with a direct organelle hit are present.
    assert sorted(evidence) == ["contig_A", "contig_D"]

    best = evidence["contig_A"]
    assert best["best_qcov"] == 100.0
    assert best["best_title"] == "Fungus two mitochondrion, complete genome"
    assert best["best_qlen"] == 29980
    assert best["best_slen"] == 27900          # the length of the BEST hit's reference
    assert best["direct_hits"] == 3            # all three direct hits are counted

    weak = evidence["contig_D"]
    assert weak["best_qcov"] == 10.0
    assert weak["best_qlen"] == 6000
    assert weak["best_slen"] == 64840
    assert weak["direct_hits"] == 1


def test_is_direct_mito_title_accepts_the_organelle_genome_wordings():
    # Every wording that says "this sequence IS mitochondrial DNA", in any case. The
    # last one is how NCBI titles the mitochondrial sequence of a genome assembly.
    titles = [
        "Fusarium graminearum mitochondrion, complete genome",
        "Example fungus MITOCHONDRION, COMPLETE SEQUENCE",
        "Example fungus mitochondrion, partial genome",
        "TPA_asm: Botrytis pseudocinerea mitochondrion, partial sequence",
        "Homo sapiens mitochondrial DNA, complete genome",
        "Neurospora crassa mitochondrial genome",
        "Cordyceps fumosorosea genome assembly, organelle: mitochondrion",
    ]
    for title in titles:
        assert mito.is_direct_mito_title(title), title


def test_is_direct_mito_title_rejects_nuclear_genes_named_after_mitochondrial_dna():
    # Records of nuclear genes, transcripts and proteins, some of them with
    # "mitochondrial DNA" or "mitochondrial genome" in their names. The first three
    # are NCBI titles that BLAST returned for fungal contigs.
    titles = [
        "Botrytis fragariae putative mitochondrial dna-directed rna polymerase protein (Bfra_005969), partial mRNA",
        "Cordyceps militaris CM01 mitochondrial DNA-directed RNA polymerase, putative (CCM_08799), partial mRNA",
        "Isaria fumosorosea ARSEF 2679 Increased loss of mitochondrial DNA protein 1 (ISF_05954), partial mRNA",
        "Example fungus mitochondrial genome maintenance protein MGM101",
        "Fusarium graminearum mitochondrial ribosomal protein L5 (MRPL5) mRNA",
        NUCLEAR_TITLE,
    ]
    for title in titles:
        assert not mito.is_direct_mito_title(title), title


def test_read_blast_counts_assembly_mitochondria_and_ignores_nuclear_gene_records(tmp_path):
    blast_path = write_blast(tmp_path / "blastout", [
        # contig_A: its best hit uses the genome-assembly wording; the second hit, a
        # nuclear gene named after mitochondrial DNA, covers more but must not count.
        blast_line("contig_A", "Cordyceps fumosorosea genome assembly, organelle: mitochondrion", 50,
                   qlen=49194, slen=27605),
        blast_line("contig_A", "Cordyceps militaris CM01 mitochondrial DNA-directed RNA polymerase, "
                   "putative (CCM_08799), partial mRNA", 90, qlen=49194, slen=4000),
        # contig_B only matches such a nuclear gene: no direct hit at all.
        blast_line("contig_B", "Isaria fumosorosea ARSEF 2679 Increased loss of mitochondrial DNA "
                   "protein 1 (ISF_05954), partial mRNA", 60, qlen=3000, slen=2500),
    ])
    evidence = mito.read_blast(blast_path)

    assert sorted(evidence) == ["contig_A"]
    best = evidence["contig_A"]
    assert best["best_qcov"] == 50.0
    assert best["best_title"] == "Cordyceps fumosorosea genome assembly, organelle: mitochondrion"
    assert best["best_slen"] == 27605
    assert best["direct_hits"] == 1


def test_read_blast_keeps_the_first_hit_when_two_are_equally_good(tmp_path):
    blast_path = write_blast(tmp_path / "blastout", [
        blast_line("contig_A", "First mitochondrion, complete sequence", 100, qlen=29980, slen=30000),
        blast_line("contig_A", "Second mitochondrion, complete sequence", 100, qlen=29980, slen=99999),
    ])
    best = mito.read_blast(blast_path)["contig_A"]
    assert best["best_title"] == "First mitochondrion, complete sequence"
    assert best["best_slen"] == 30000
    assert best["direct_hits"] == 2


def test_read_blast_skips_short_and_blank_lines(tmp_path):
    blast_path = write_text(
        tmp_path / "blastout",
        "\n"                                              # blank line
        "contig_A\t5062\t100\n"                           # only 3 columns: not a hit line
        + blast_line("contig_A", MITO_TITLE, 100, qlen=29980, slen=30000) + "\n",
    )
    evidence = mito.read_blast(blast_path)
    assert list(evidence) == ["contig_A"]
    assert evidence["contig_A"]["direct_hits"] == 1


def test_read_blast_fails_clearly_on_a_damaged_number(tmp_path):
    bad_qcovs = write_blast(tmp_path / "bad_qcovs", [blast_line("contig_A", MITO_TITLE, "abc")])
    with pytest.raises(ValueError, match="Invalid BLAST qcovs value"):
        mito.read_blast(bad_qcovs)

    bad_qlen = write_blast(tmp_path / "bad_qlen", [blast_line("contig_A", MITO_TITLE, 100, qlen="abc")])
    with pytest.raises(ValueError, match="Invalid BLAST qlen value"):
        mito.read_blast(bad_qlen)

    bad_slen = write_blast(tmp_path / "bad_slen", [blast_line("contig_A", MITO_TITLE, 100, slen="abc")])
    with pytest.raises(ValueError, match="Invalid BLAST slen value"):
        mito.read_blast(bad_slen)

    # A reference of length zero cannot be real and would break the length ratio.
    zero_slen = write_blast(tmp_path / "zero_slen", [blast_line("contig_A", MITO_TITLE, 100, slen=0)])
    with pytest.raises(ValueError, match="Invalid BLAST slen value"):
        mito.read_blast(zero_slen)


# =============================================================================
# 4. The concatemer flag (informational columns at the end of the report)
# =============================================================================
# A contig that is two copies of the mitochondrial genome end to end is a known
# assembler artefact. It still matches the organelle reference, so
# it is rightly confirmed as mitochondrial, but it is about twice as long as the
# reference. The report now says so in three extra columns. They are for the
# reader only: they must never change a classification, a reason, the list of
# confirmed IDs or the split.

def test_concatemer_columns_two_copies_end_to_end_is_flagged():
    # A 54,000 bp contig against a 30,000 bp reference: ratio 1.80.
    evidence = make_hit(100.0)
    assert mito.concatemer_columns(54000, evidence, None) == ("30000", "1.80", "Y")
    # A linear contig according to the assembler is flagged the same way.
    assert mito.concatemer_columns(54000, evidence, make_assembler(175.0, circular=False)) == ("30000", "1.80", "Y")


def test_concatemer_columns_one_copy_is_not_flagged():
    # One genome copy: a 29,980 bp contig against a 30,000 bp reference.
    assert mito.concatemer_columns(29980, make_hit(100.0), None) == ("30000", "1.00", "N")


def test_concatemer_columns_a_circular_contig_is_not_flagged():
    # The same long contig, but the assembler closed it into a circle: it is
    # taken to be a genuine larger genome, so no flag.
    result = mito.concatemer_columns(54000, make_hit(100.0), make_assembler(170.0, circular=True))
    assert result == ("30000", "1.80", "N")


def test_concatemer_columns_threshold_is_above_1_5():
    # A reference of 100 bases keeps the arithmetic easy to follow.
    reference_100 = {"best_slen": 100}
    # Exactly 1.5 times the reference is not "above 1.5".
    assert mito.concatemer_columns(150, reference_100, None) == ("100", "1.50", "N")
    # 1.51 times is.
    assert mito.concatemer_columns(151, reference_100, None) == ("100", "1.51", "Y")


def test_concatemer_flag_is_judged_on_the_rounded_ratio_that_is_printed():
    # 1504 / 1000 = 1.504 is shown as 1.50, so it is not flagged: the flag can
    # always be checked against the number a reader sees in the report.
    reference_1000 = {"best_slen": 1000}
    assert mito.concatemer_columns(1504, reference_1000, None) == ("1000", "1.50", "N")


def test_concatemer_columns_a_fragment_shorter_than_its_reference_is_not_flagged():
    assert mito.concatemer_columns(60, {"best_slen": 100}, None) == ("100", "0.60", "N")


def test_concatemer_columns_without_a_direct_hit_are_empty():
    assert mito.concatemer_columns(54000, None, None) == ("", "", "N")


def test_screen_report_ends_with_the_new_columns(tmp_path):
    # Three contigs: a nuclear contig without hit, a single-copy mitogenome, and a
    # double-copy mitogenome. BLAST-only (no assembler evidence), as in short-read runs.
    contigs_path = write_contigs(tmp_path / "contigs.fasta", [
        ("nuclear_1", 5000), ("mito_single", 29980), ("mito_double", 54000),
    ])
    blast_path = write_blast(tmp_path / "blastout", [
        blast_line("nuclear_1", NUCLEAR_TITLE, 90, qlen=5000, slen=40000000),
        blast_line("mito_single", MITO_TITLE, 100, qlen=29980, slen=30000),
        blast_line("mito_double", MITO_TITLE, 100, qlen=54000, slen=30000),
    ])
    report_path = tmp_path / "mitochondrial_screen.tsv"
    ids_path = tmp_path / "confirmed.ids"

    result = run_cli(
        "screen", "--blast", blast_path, "--contigs", contigs_path,
        "--report", str(report_path), "--confirmed-ids", str(ids_path),
    )
    assert result.returncode == 0, result.stderr

    column_names, rows = read_table(report_path)
    # The 13 older columns keep their places; the four new ones come last.
    assert column_names == OLD_REPORT_COLUMNS + NEW_REPORT_COLUMNS

    by_contig = rows_by_contig(report_path)
    # No direct organelle hit: nothing to compare with, so empty values and N.
    assert new_column_values(by_contig["nuclear_1"]) == ("", "", "N")
    # One copy of the genome: ratio 29,980 / 30,000 = 1.00, not flagged.
    assert new_column_values(by_contig["mito_single"]) == ("30000", "1.00", "N")
    # Two copies end to end: ratio 54,000 / 30,000 = 1.80, flagged.
    assert new_column_values(by_contig["mito_double"]) == ("30000", "1.80", "Y")


def test_concatemer_flag_changes_no_decision(tmp_path):
    # The same 54,000 bp contig with the same BLAST hit, once against a 30,000 bp
    # reference (ratio 1.80, flagged Y) and once against a 54,000 bp reference
    # (ratio 1.00, flagged N). Everything except the three new columns must be
    # identical: classification, reason, the confirmed-ID file and the split.
    contigs_path = write_contigs(tmp_path / "contigs.fasta", [("nuclear_1", 5000), ("mito_x", 54000)])
    outcomes = {}
    for label, reference_length in (("flagged", 30000), ("not_flagged", 54000)):
        blast_path = write_blast(tmp_path / ("blastout_" + label), [
            blast_line("mito_x", MITO_TITLE, 100, qlen=54000, slen=reference_length),
        ])
        report_path = tmp_path / ("report_" + label + ".tsv")
        ids_path = tmp_path / ("ids_" + label + ".txt")
        nuclear_path = tmp_path / ("nuclear_" + label + ".fasta")
        mito_path = tmp_path / ("mito_" + label + ".fasta")

        screen_result = run_cli(
            "screen", "--blast", blast_path, "--contigs", contigs_path,
            "--report", str(report_path), "--confirmed-ids", str(ids_path),
        )
        assert screen_result.returncode == 0, screen_result.stderr
        split_result = run_cli(
            "split", "--fasta", contigs_path, "--confirmed-ids", str(ids_path),
            "--nuclear-fasta", str(nuclear_path), "--mitochondrial-fasta", str(mito_path),
        )
        assert split_result.returncode == 0, split_result.stderr

        row = rows_by_contig(report_path)["mito_x"]
        outcomes[label] = {
            "old_columns": [row[name] for name in OLD_REPORT_COLUMNS],
            "flag": row["concatemer_flag"],
            "ids": ids_path.read_text(),
            "nuclear": nuclear_path.read_text(),
            "mito": mito_path.read_text(),
        }

    assert outcomes["flagged"]["flag"] == "Y"
    assert outcomes["not_flagged"]["flag"] == "N"
    # Still confirmed, still selected, still split off, with the same reason.
    assert outcomes["flagged"]["old_columns"][2:4] == ["confirmed_mitochondrial", "strong_direct_blast"]
    assert outcomes["flagged"]["old_columns"] == outcomes["not_flagged"]["old_columns"]
    assert outcomes["flagged"]["ids"] == outcomes["not_flagged"]["ids"] == "mito_x\n"
    assert outcomes["flagged"]["nuclear"] == outcomes["not_flagged"]["nuclear"]
    assert outcomes["flagged"]["mito"] == outcomes["not_flagged"]["mito"]
    assert outcomes["flagged"]["mito"].startswith(">mito_x\n")


def test_hifiasm_c_suffix_decides_the_concatemer_flag_end_to_end(tmp_path):
    # The whole long-read chain on hand-written files: hifiasm GFA ->
    # evidence-hifiasm -> screen. Three mitochondrial-looking contigs, all with a
    # full-length BLAST hit to a 30,000 bp reference:
    #   ptg000007c  29,980 bp   closed circle (c)     one copy        -> flag N
    #   ptg000008l  54,000 bp   linear (l)             two copies      -> flag Y
    #   ptg000009c  54,000 bp   closed circle (c)      long but a circle -> flag N
    # Six nuclear contigs at depth 20 come with them. They have to be the clear
    # majority, as in a real assembly, because the median depth of all contigs
    # is the "normal" depth that the three deep contigs are compared with.
    gfa_path = write_text(tmp_path / "hifiasm.p_ctg.gfa", (
        "S\tptg000001l\tACGT\tLN:i:4\trd:i:20\n"
        "S\tptg000002l\tACGT\tLN:i:4\trd:i:20\n"
        "S\tptg000003l\tACGT\tLN:i:4\trd:i:20\n"
        "S\tptg000004l\tACGT\tLN:i:4\trd:i:20\n"
        "S\tptg000005l\tACGT\tLN:i:4\trd:i:20\n"
        "S\tptg000006l\tACGT\tLN:i:4\trd:i:20\n"
        "S\tptg000007c\tACGT\tLN:i:4\trd:i:150\n"
        "S\tptg000008l\tACGT\tLN:i:4\trd:i:175\n"
        "S\tptg000009c\tACGT\tLN:i:4\trd:i:170\n"
    ))
    contigs_path = write_contigs(tmp_path / "contigs.fasta", [
        ("ptg000001l", 5000), ("ptg000002l", 6000), ("ptg000003l", 7000),
        ("ptg000004l", 8000), ("ptg000005l", 9000), ("ptg000006l", 10000),
        ("ptg000007c", 29980), ("ptg000008l", 54000), ("ptg000009c", 54000),
    ])
    blast_path = write_blast(tmp_path / "blastout", [
        blast_line("ptg000007c", MITO_TITLE, 100, qlen=29980, slen=30000),
        blast_line("ptg000008l", MITO_TITLE, 100, qlen=54000, slen=30000),
        blast_line("ptg000009c", MITO_TITLE, 100, qlen=54000, slen=30000),
    ])
    evidence_path = tmp_path / "assembly_evidence.tsv"
    report_path = tmp_path / "mitochondrial_screen.tsv"
    ids_path = tmp_path / "confirmed.ids"

    evidence_result = run_cli("evidence-hifiasm", "--gfa", gfa_path, "--output", str(evidence_path))
    assert evidence_result.returncode == 0, evidence_result.stderr
    screen_result = run_cli(
        "screen", "--blast", blast_path, "--contigs", contigs_path,
        "--assembly-evidence", str(evidence_path), "--report", str(report_path),
        "--confirmed-ids", str(ids_path), "--min-direct-blast-qcov", "30",
        "--min-assembler-coverage-ratio", "3.0",
    )
    assert screen_result.returncode == 0, screen_result.stderr

    rows = rows_by_contig(report_path)
    # Circularity now comes from the c in the name.
    assert rows["ptg000007c"]["assembler_circular"] == "Y"
    assert rows["ptg000008l"]["assembler_circular"] == "N"
    assert rows["ptg000009c"]["assembler_circular"] == "Y"
    # The concatemer flag follows from length ratio and circularity.
    assert rows["ptg000007c"]["length_ratio"] == "1.00"
    assert rows["ptg000007c"]["concatemer_flag"] == "N"
    assert rows["ptg000008l"]["length_ratio"] == "1.80"
    assert rows["ptg000008l"]["concatemer_flag"] == "Y"
    assert rows["ptg000009c"]["length_ratio"] == "1.80"
    assert rows["ptg000009c"]["concatemer_flag"] == "N"
    # All three are confirmed for the same reason: depth supports the BLAST hit.
    for contig in ("ptg000007c", "ptg000008l", "ptg000009c"):
        assert rows[contig]["classification"] == "confirmed_mitochondrial"
        assert rows[contig]["reason"] == "direct_blast_plus_high_coverage"
    # The nuclear contigs are untouched.
    assert rows["ptg000001l"]["classification"] == "not_mitochondrial"
    assert rows["ptg000001l"]["concatemer_flag"] == "N"
    assert ids_path.read_text() == "ptg000007c\nptg000008l\nptg000009c\n"


# =============================================================================
# 5. split: the nuclear and the mitochondrial FASTA
# =============================================================================
# split reads the final assembly and the list of confirmed IDs and writes two
# FASTA files: the nuclear contigs (which go on to assembly QC and gene
# prediction) and the mitochondrial contigs. It never changes a sequence.

SPLIT_INPUT_FASTA = (
    ">ctg1 first contig\n"
    "ACGTACGTAC\n"
    "GTACGTACGT\n"
    ">ctg2 the mitochondrion\n"
    "TTTTGGGGCC\n"
    ">ctg3\n"
    "AAAACCCCGG\n"
    ">ctg4 last contig\n"
    "GGGGTTTTAA\n"
)


def split_files(tmp_path, ids_text, fasta_text=SPLIT_INPUT_FASTA):
    """Write the split inputs and give back (arguments for the split command, output paths)."""
    fasta_path = write_text(tmp_path / "assembly.fasta", fasta_text)
    ids_path = write_text(tmp_path / "confirmed.ids", ids_text)
    nuclear_path = str(tmp_path / "out" / "nuclear.fasta")
    mito_path = str(tmp_path / "out" / "mito.fasta")
    arguments = [
        "split", "--fasta", fasta_path, "--confirmed-ids", ids_path,
        "--nuclear-fasta", nuclear_path, "--mitochondrial-fasta", mito_path,
    ]
    return arguments, nuclear_path, mito_path


def test_split_writes_the_right_two_fasta_files(tmp_path):
    arguments, nuclear_path, mito_path = split_files(tmp_path, "ctg2\n")
    args = parse(*arguments)
    args.func(args)

    # Nuclear: every contig that is not confirmed, in the original order, with
    # its full header; the sequence is re-written on one line (shorter than 80).
    expected_nuclear = (
        ">ctg1 first contig\n"
        "ACGTACGTACGTACGTACGT\n"
        ">ctg3\n"
        "AAAACCCCGG\n"
        ">ctg4 last contig\n"
        "GGGGTTTTAA\n"
    )
    expected_mito = (
        ">ctg2 the mitochondrion\n"
        "TTTTGGGGCC\n"
    )
    with open(nuclear_path, encoding="utf-8") as handle:
        assert handle.read() == expected_nuclear
    with open(mito_path, encoding="utf-8") as handle:
        assert handle.read() == expected_mito


def test_split_accepts_comments_and_blank_lines_in_the_id_file(tmp_path):
    arguments, nuclear_path, mito_path = split_files(tmp_path, "# confirmed by screen\n\nctg2\n\n")
    args = parse(*arguments)
    args.func(args)
    with open(mito_path, encoding="utf-8") as handle:
        assert handle.read().startswith(">ctg2 the mitochondrion\n")


def test_split_wraps_long_sequences_at_80_bases(tmp_path):
    long_fasta = ">ctg1\n" + "A" * 200 + "\n>ctg2\nCCCC\n"
    arguments, nuclear_path, mito_path = split_files(tmp_path, "ctg2\n", fasta_text=long_fasta)
    args = parse(*arguments)
    args.func(args)
    with open(nuclear_path, encoding="utf-8") as handle:
        nuclear_lines = handle.read().splitlines()
    assert nuclear_lines == [">ctg1", "A" * 80, "A" * 80, "A" * 40]


def test_split_with_no_confirmed_contig_gives_an_empty_mitochondrial_file(tmp_path):
    # "report" mode writes an empty ID list. The split rule still runs, and then
    # every contig stays nuclear and the mitochondrial file is empty.
    arguments, nuclear_path, mito_path = split_files(tmp_path, "")
    args = parse(*arguments)
    args.func(args)
    assert os.path.getsize(mito_path) == 0
    with open(nuclear_path, encoding="utf-8") as handle:
        nuclear_text = handle.read()
    assert nuclear_text.count(">") == 4


def test_split_fails_clearly_when_a_confirmed_id_is_not_in_the_fasta(tmp_path):
    # The split may run on a polished assembly rather than the FASTA that was
    # screened, so a confirmed contig that has disappeared must stop the run.
    arguments, nuclear_path, mito_path = split_files(tmp_path, "ctg2\nctg_missing\nanother_missing\n")
    args = parse(*arguments)
    with pytest.raises(ValueError) as error_info:
        args.func(args)
    message = str(error_info.value)
    assert "absent from final assembly" in message
    # Both missing IDs are named (sorted); the one that exists is not.
    assert "another_missing, ctg_missing" in message
    assert "ctg2" not in message
    # Nothing was written, so no half-finished output can be picked up later.
    assert not os.path.exists(nuclear_path)
    assert not os.path.exists(mito_path)


def test_split_fails_when_nothing_would_be_left_for_the_nuclear_assembly(tmp_path):
    arguments, nuclear_path, mito_path = split_files(tmp_path, "ctg1\nctg2\nctg3\nctg4\n")
    args = parse(*arguments)
    with pytest.raises(ValueError, match="empty nuclear assembly"):
        args.func(args)
    assert not os.path.exists(nuclear_path)


def test_read_fasta_rejects_duplicate_contig_names(tmp_path):
    duplicated = write_text(tmp_path / "dup.fasta", ">ctg1 one\nAAAA\n>ctg1 two\nCCCC\n")
    with pytest.raises(ValueError, match="Duplicate FASTA identifier: ctg1"):
        mito.read_fasta(duplicated)


# =============================================================================
# 6. The command line, as the Snakemake rules call it
# =============================================================================
# One small data set used by several tests: three nuclear contigs at depth 20
# and one mitochondrial candidate at depth 150 (7.5 times the median) whose
# BLAST hit covers 50 percent of it. 50 is above the default threshold of 35 but
# below the "very strong" 70, so the call depends on the depth ratio threshold.

CLI_FLYE_INFO = (
    "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
    "nuc_1\t50000\t20\tN\tN\t1\t*\t1\n"
    "nuc_2\t40000\t20\tN\tN\t1\t*\t2\n"
    "nuc_3\t30000\t20\tN\tN\t1\t*\t3\n"
    "mito_1\t28000\t150\tN\tN\t1\t*\t4\n"
)


def cli_dataset(tmp_path):
    """Write the small data set; give back the paths the screen command needs."""
    contigs_path = write_contigs(tmp_path / "contigs.fasta", [
        ("nuc_1", 5000), ("nuc_2", 4000), ("nuc_3", 3000), ("mito_1", 2800),
    ])
    blast_path = write_blast(tmp_path / "blastout", [
        blast_line("mito_1", MITO_TITLE, 50, qlen=2800, slen=2800),
    ])
    flye_info_path = write_text(tmp_path / "assembly_info.txt", CLI_FLYE_INFO)
    return contigs_path, blast_path, flye_info_path


def test_cli_parser_defaults_are_the_documented_ones():
    args = parse("screen", "--blast", "b", "--contigs", "c", "--report", "r", "--confirmed-ids", "i")
    assert args.mode == "separate_confirmed"
    assert args.min_direct_blast_qcov == 35.0
    assert args.min_assembler_coverage_ratio == 3.0
    assert args.very_strong_qcov == 70.0        # new option; default unchanged from v0.2
    assert args.assembly_evidence == ""
    assert args.flye_info == ""


def test_cli_old_flag_name_min_flye_coverage_ratio_still_works(tmp_path):
    # Backward compatibility: command lines written for version 0.1 of the script
    # call screen with --min-flye-coverage-ratio and --flye-info. The old name
    # must be accepted AND must act, not just be ignored. The same contig is
    # called differently at ratio threshold 3 and at 10 (its depth ratio is 7.5).
    contigs_path, blast_path, flye_info_path = cli_dataset(tmp_path)

    low_report = tmp_path / "low.tsv"
    low_ids = tmp_path / "low.ids"
    low_result = run_cli(
        "screen", "--blast", blast_path, "--contigs", contigs_path,
        "--flye-info", flye_info_path, "--report", str(low_report),
        "--confirmed-ids", str(low_ids), "--mode", "separate_confirmed",
        "--min-direct-blast-qcov", "35", "--min-flye-coverage-ratio", "3.0",
    )
    assert low_result.returncode == 0, low_result.stderr
    assert "1 confirmed call(s)" in low_result.stdout
    low_row = rows_by_contig(low_report)["mito_1"]
    assert low_row["classification"] == "confirmed_mitochondrial"
    assert low_row["reason"] == "direct_blast_plus_high_coverage"
    assert low_ids.read_text() == "mito_1\n"

    high_report = tmp_path / "high.tsv"
    high_ids = tmp_path / "high.ids"
    high_result = run_cli(
        "screen", "--blast", blast_path, "--contigs", contigs_path,
        "--flye-info", flye_info_path, "--report", str(high_report),
        "--confirmed-ids", str(high_ids), "--mode", "separate_confirmed",
        "--min-direct-blast-qcov", "35", "--min-flye-coverage-ratio", "10",
    )
    assert high_result.returncode == 0, high_result.stderr
    high_row = rows_by_contig(high_report)["mito_1"]
    assert high_row["classification"] == "mitochondrial_candidate"
    assert high_row["reason"] == "direct_blast_without_support"
    assert high_ids.read_text() == ""


def test_cli_old_and_new_flag_names_give_identical_reports(tmp_path):
    contigs_path, blast_path, flye_info_path = cli_dataset(tmp_path)
    old_report = tmp_path / "old.tsv"
    new_report = tmp_path / "new.tsv"
    common = ["screen", "--blast", blast_path, "--contigs", contigs_path, "--flye-info", flye_info_path]

    old_result = run_cli(*common, "--report", str(old_report), "--confirmed-ids", str(tmp_path / "o.ids"),
                         "--min-flye-coverage-ratio", "5")
    new_result = run_cli(*common, "--report", str(new_report), "--confirmed-ids", str(tmp_path / "n.ids"),
                         "--min-assembler-coverage-ratio", "5")
    assert old_result.returncode == 0, old_result.stderr
    assert new_result.returncode == 0, new_result.stderr
    assert old_report.read_text() == new_report.read_text()


def test_cli_flye_info_path_and_evidence_table_path_agree(tmp_path):
    # Two ways to give Flye information to screen: directly (--flye-info, the
    # compatibility path) or through the common table made by evidence-flye
    # (--assembly-evidence, what the v2 rules use). The report must be the same.
    contigs_path, blast_path, flye_info_path = cli_dataset(tmp_path)

    direct_report = tmp_path / "direct.tsv"
    direct_result = run_cli(
        "screen", "--blast", blast_path, "--contigs", contigs_path,
        "--flye-info", flye_info_path, "--report", str(direct_report),
        "--confirmed-ids", str(tmp_path / "direct.ids"),
    )
    assert direct_result.returncode == 0, direct_result.stderr

    evidence_path = tmp_path / "assembly_evidence.tsv"
    evidence_result = run_cli("evidence-flye", "--flye-info", flye_info_path, "--output", str(evidence_path))
    assert evidence_result.returncode == 0, evidence_result.stderr
    table_report = tmp_path / "table.tsv"
    table_result = run_cli(
        "screen", "--blast", blast_path, "--contigs", contigs_path,
        "--assembly-evidence", str(evidence_path), "--report", str(table_report),
        "--confirmed-ids", str(tmp_path / "table.ids"),
    )
    assert table_result.returncode == 0, table_result.stderr

    assert direct_report.read_text() == table_report.read_text()


def test_cli_very_strong_qcov_option_reaches_classify(tmp_path):
    # With a ratio threshold of 10 the depth does not support the 50 percent hit,
    # so the call depends only on the very-strong cut-off: 70 (default) leaves it a
    # candidate, 50 confirms it. This proves the option is passed all the way to
    # classify() and not silently replaced by the old fixed 70.
    contigs_path, blast_path, flye_info_path = cli_dataset(tmp_path)
    common = [
        "screen", "--blast", blast_path, "--contigs", contigs_path, "--flye-info", flye_info_path,
        "--min-assembler-coverage-ratio", "10",
    ]

    default_report = tmp_path / "default.tsv"
    default_result = run_cli(*common, "--report", str(default_report),
                             "--confirmed-ids", str(tmp_path / "default.ids"))
    assert default_result.returncode == 0, default_result.stderr
    assert rows_by_contig(default_report)["mito_1"]["reason"] == "direct_blast_without_support"

    lowered_report = tmp_path / "lowered.tsv"
    lowered_ids = tmp_path / "lowered.ids"
    lowered_result = run_cli(*common, "--very-strong-qcov", "50", "--report", str(lowered_report),
                             "--confirmed-ids", str(lowered_ids))
    assert lowered_result.returncode == 0, lowered_result.stderr
    lowered_row = rows_by_contig(lowered_report)["mito_1"]
    assert lowered_row["classification"] == "confirmed_mitochondrial"
    assert lowered_row["reason"] == "very_strong_direct_blast"
    assert lowered_ids.read_text() == "mito_1\n"


def test_cli_modes_decide_what_is_selected(tmp_path):
    # separate_confirmed: confirmed IDs are listed (the split uses them).
    # report: the same labels, but the list is empty so the assembly stays as it was.
    # off: nothing is screened at all.
    contigs_path, blast_path, flye_info_path = cli_dataset(tmp_path)

    expectations = {
        "separate_confirmed": ("confirmed_mitochondrial", "mito_1\n"),
        "report": ("confirmed_mitochondrial", ""),
        "off": ("not_screened", ""),
    }
    for mode, (expected_class, expected_ids) in expectations.items():
        report_path = tmp_path / (mode + ".tsv")
        ids_path = tmp_path / (mode + ".ids")
        result = run_cli(
            "screen", "--blast", blast_path, "--contigs", contigs_path,
            "--flye-info", flye_info_path, "--report", str(report_path),
            "--confirmed-ids", str(ids_path), "--mode", mode,
        )
        assert result.returncode == 0, result.stderr
        assert rows_by_contig(report_path)["mito_1"]["classification"] == expected_class
        assert ids_path.read_text() == expected_ids


def test_cli_screen_without_assembler_evidence_leaves_the_assembler_columns_empty(tmp_path):
    # Short-read use: BLAST only. The assembler columns stay empty and the call is
    # strong_direct_blast, as in the FunFlux short-read run.
    contigs_path, blast_path, flye_info_path = cli_dataset(tmp_path)
    report_path = tmp_path / "blast_only.tsv"
    result = run_cli(
        "screen", "--blast", blast_path, "--contigs", contigs_path,
        "--report", str(report_path), "--confirmed-ids", str(tmp_path / "blast_only.ids"),
        "--min-direct-blast-qcov", "30",
    )
    assert result.returncode == 0, result.stderr
    row = rows_by_contig(report_path)["mito_1"]
    assert (row["classification"], row["reason"]) == ("confirmed_mitochondrial", "strong_direct_blast")
    for column in ("assembler_source", "assembler_circular", "assembler_coverage",
                   "assembler_coverage_ratio", "assembler_repeat", "assembler_alt_group"):
        assert row[column] == ""


def test_cli_screen_with_a_missing_blast_file_ends_with_a_readable_error(tmp_path):
    contigs_path, blast_path, flye_info_path = cli_dataset(tmp_path)
    result = run_cli(
        "screen", "--blast", str(tmp_path / "does_not_exist"), "--contigs", contigs_path,
        "--report", str(tmp_path / "r.tsv"), "--confirmed-ids", str(tmp_path / "i.ids"),
    )
    assert result.returncode == 1
    assert result.stderr.startswith("ERROR:")
    assert "does_not_exist" in result.stderr
    assert "Traceback" not in result.stderr


def test_cli_split_with_a_missing_confirmed_id_ends_with_a_readable_error(tmp_path):
    arguments, nuclear_path, mito_path = split_files(tmp_path, "ctg_missing\n")
    result = run_cli(*arguments)
    assert result.returncode == 1
    assert result.stderr.startswith("ERROR: Confirmed mitochondrial contig(s) are absent from final assembly")
    assert "ctg_missing" in result.stderr
    assert "Traceback" not in result.stderr
    assert not os.path.exists(nuclear_path)
    assert not os.path.exists(mito_path)


def test_cli_split_accepts_the_short_option_names(tmp_path):
    # --nuclear and --mitochondrial are older spellings of --nuclear-fasta and
    # --mitochondrial-fasta and must keep working.
    fasta_path = write_text(tmp_path / "assembly.fasta", SPLIT_INPUT_FASTA)
    ids_path = write_text(tmp_path / "confirmed.ids", "ctg2\n")
    nuclear_path = tmp_path / "nuclear.fasta"
    mito_path = tmp_path / "mito.fasta"
    result = run_cli(
        "split", "--fasta", fasta_path, "--confirmed-ids", ids_path,
        "--nuclear", str(nuclear_path), "--mitochondrial", str(mito_path),
    )
    assert result.returncode == 0, result.stderr
    assert "1 mitochondrial and 3 nuclear contig(s)" in result.stdout
    assert mito_path.read_text().startswith(">ctg2 the mitochondrion\n")


# =============================================================================
# The species check (CHANGED in FunFlux 2.1.0): another organism's mitochondrion
# =============================================================================

def write_nodes(path):
    """Write the small taxonomy of TAXONOMY_LINES as a nodes.dmp file."""
    return write_text(path, "\n".join(TAXONOMY_LINES) + "\n")


def test_read_taxonomy_parents_reads_the_first_two_fields(tmp_path):
    parents = mito.read_taxonomy_parents(write_nodes(tmp_path / "nodes.dmp"))
    assert parents[5062] == 4890
    assert parents[9606] == 33208
    assert parents[1] == 1


def test_kingdom_of_walks_up_to_the_first_kingdom(tmp_path):
    parents = mito.read_taxonomy_parents(write_nodes(tmp_path / "nodes.dmp"))
    assert mito.kingdom_of("5062", parents) == "Fungi"
    assert mito.kingdom_of("9606", parents) == "Metazoa"
    assert mito.kingdom_of("3702", parents) == "Viridiplantae"
    # Several IDs: one fungal ID is enough.
    assert mito.kingdom_of("9606;5062", parents) == "Fungi"


def test_kingdom_of_cannot_tell_without_a_clear_kingdom(tmp_path):
    parents = mito.read_taxonomy_parents(write_nodes(tmp_path / "nodes.dmp"))
    # 0 (a BLAST database without taxonomy IDs), text, an empty cell, an unknown ID.
    assert mito.kingdom_of("0", parents) is None
    assert mito.kingdom_of("N/A", parents) is None
    assert mito.kingdom_of("", parents) is None
    assert mito.kingdom_of("123456", parents) is None
    # Filed under Eukaryota only, and an oomycete: neither is called foreign.
    assert mito.kingdom_of("100272", parents) is None
    assert mito.kingdom_of("4787", parents) is None
    # A damaged file with a loop: the walk stops, and the answer is "cannot tell".
    assert mito.kingdom_of("777", {777: 778, 778: 777}) is None


def screen_two_contigs(tmp_path, other_title, other_taxids, other_qcovs=100,
                       with_taxonomy=True, with_evidence=False):
    """Screen a fungal mitochondrial contig and one other; give back (rows, confirmed IDs, stdout).

    The other contig is 1,139 bp, like the human pieces of the B05.10 Illumina assembly.
    With with_evidence, an assembler table says both contigs are circular and deep, as a
    long-read assembler would.
    """
    contigs_path = write_contigs(tmp_path / "contigs.fasta", [
        ("mito_fungal", 29980), ("mito_other", 1139), ("nuclear_1", 50000),
    ])
    blast_path = write_blast(tmp_path / "blastout", [
        blast_line("mito_fungal", MITO_TITLE, 99, qlen=29980, slen=30000, taxids="5062"),
        blast_line("mito_other", other_title, other_qcovs, qlen=1139, slen=16569, taxids=other_taxids),
    ])
    report_path = tmp_path / "mitochondrial_screen.tsv"
    ids_path = tmp_path / "confirmed.ids"
    arguments = [
        "screen", "--blast", blast_path, "--contigs", contigs_path,
        "--report", str(report_path), "--confirmed-ids", str(ids_path),
    ]
    if with_taxonomy:
        arguments += ["--taxonomy-nodes", write_nodes(tmp_path / "nodes.dmp")]
    if with_evidence:
        evidence_path = write_text(tmp_path / "assembly_evidence.tsv",
            "contig\tsource\tcircular\tcoverage\trepeat\talt_group\n"
            "mito_fungal\thifiasm\tY\t400\tN\t*\n"
            "mito_other\thifiasm\tY\t300\tN\t*\n"
            "nuclear_1\thifiasm\tN\t50\tN\t*\n")
        arguments += ["--assembly-evidence", evidence_path]
    result = run_cli(*arguments)
    assert result.returncode == 0, result.stderr
    with open(str(ids_path), encoding="utf-8") as handle:
        confirmed = handle.read().split()
    return rows_by_contig(report_path), confirmed, result.stdout


def test_screen_reports_a_human_mitochondrion_as_foreign(tmp_path):
    rows, confirmed, stdout = screen_two_contigs(tmp_path, HUMAN_MITO_TITLE, "9606")
    assert rows["mito_fungal"]["classification"] == "confirmed_mitochondrial"
    assert rows["mito_fungal"]["best_mito_kingdom"] == "Fungi"
    assert rows["mito_other"]["classification"] == "foreign_mitochondrial"
    assert rows["mito_other"]["reason"] == "best_direct_hit_animal_or_plant"
    assert rows["mito_other"]["best_mito_kingdom"] == "Metazoa"
    # Only the fungus's own mitochondrion is set apart (and protected from the genus filter).
    assert confirmed == ["mito_fungal"]
    assert "foreign_mitochondrial" in stdout and "mito_other" in stdout


def test_screen_reports_a_plant_mitochondrion_as_foreign(tmp_path):
    plant_title = "Arabidopsis thaliana mitochondrion, complete genome"
    rows, confirmed, stdout = screen_two_contigs(tmp_path, plant_title, "3702")
    assert rows["mito_other"]["classification"] == "foreign_mitochondrial"
    assert rows["mito_other"]["best_mito_kingdom"] == "Viridiplantae"


def test_a_weak_human_hit_is_foreign_too(tmp_path):
    # 20 % query cover makes a candidate; with a human best hit it is foreign instead.
    rows, confirmed, stdout = screen_two_contigs(tmp_path, HUMAN_MITO_TITLE, "9606", other_qcovs=20)
    assert rows["mito_other"]["classification"] == "foreign_mitochondrial"


def test_a_deep_circular_human_contig_from_a_long_read_assembler_is_foreign(tmp_path):
    rows, confirmed, stdout = screen_two_contigs(tmp_path, HUMAN_MITO_TITLE, "9606", with_evidence=True)
    assert rows["mito_fungal"]["classification"] == "confirmed_mitochondrial"
    assert rows["mito_other"]["classification"] == "foreign_mitochondrial"
    assert rows["mito_other"]["assembler_circular"] == "Y"
    assert confirmed == ["mito_fungal"]


def test_screen_without_the_taxonomy_keeps_the_older_calls(tmp_path):
    rows, confirmed, stdout = screen_two_contigs(tmp_path, HUMAN_MITO_TITLE, "9606", with_taxonomy=False)
    assert rows["mito_other"]["classification"] == "confirmed_mitochondrial"
    assert rows["mito_other"]["best_mito_kingdom"] == ""
    assert confirmed == ["mito_fungal", "mito_other"]


def test_screen_leaves_the_call_alone_when_it_cannot_tell(tmp_path):
    # A BLAST database without taxonomy IDs writes 0 in staxids.
    rows, confirmed, stdout = screen_two_contigs(tmp_path, HUMAN_MITO_TITLE, "0")
    assert rows["mito_other"]["classification"] == "confirmed_mitochondrial"
    assert rows["mito_other"]["best_mito_kingdom"] == ""
    assert "mito_other" in confirmed


def test_an_oomycete_keeps_its_own_mitochondrion(tmp_path):
    # Phytophthora is not in NCBI's Fungi, but it is neither an animal nor a plant.
    oomycete_title = "Phytophthora infestans mitochondrion, complete genome"
    rows, confirmed, stdout = screen_two_contigs(tmp_path, oomycete_title, "4787")
    assert rows["mito_other"]["classification"] == "confirmed_mitochondrial"
    assert rows["mito_other"]["best_mito_kingdom"] == ""
