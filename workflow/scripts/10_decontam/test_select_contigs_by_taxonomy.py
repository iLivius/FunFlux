#!/usr/bin/env python3
"""Tests for select_contigs_by_taxonomy.py, the contig filter of the decontamination step.

What the script does. It reads a BlobTools table (one row per contig, with the
genus that BLAST assigned to it) and decides which contigs of the assembly to
keep. It writes four files: the kept contig IDs, the kept sequences (a FASTA), a
genus composition report, and contig_taxonomy_decisions.tsv, an audit table that
gives every contig an action (keep or remove) and a reason. There are three
modes: off (keep everything, the fungal default), include (keep only the listed
genera) and exclude (drop only the listed genera). BacFlux's fourth mode, auto
(keep the most common genus), is refused with an explanation in FunFlux.
--always-keep names contigs that must survive whatever the taxonomy says; in the
fungal workflow these are the confirmed mitochondrial contigs.

How these tests work. Each test writes a made-up assembly of six contigs into
pytest's temporary folder and runs the script on it the way the Snakemake rule
select_contigs does: as a command line, in its own Python process. That also
tests the real argument handling and the real exit status. Nothing here needs
BlobTools, BLAST or real data, and apart from pytest only the Python standard
library is used.

The six contigs (the genus is what BlobTools would report for them):

    contig_1   Fusarium      120 bp   the genus most contigs carry
    contig_2   Fusarium      100 bp
    contig_3   Fusarium       80 bp
    contig_4   Trichoderma    60 bp   a different fungus: contamination here
    contig_5   no-hit         40 bp   BLAST found nothing to call it
    contig_6   Botrytis       32 bp   stands in for the mitochondrion: BLAST
                                      gives it another genus than most contigs

contig_6 copies a case that can occur with a poorly sequenced fungus: its
confirmed mitochondrial contig is labelled with the true genus, while most of
the assembly's DNA is labelled with a better-sequenced genus. Here Botrytis
plays the true genus and Fusarium the majority label. Include Fusarium or
exclude Botrytis would remove the mitochondrion, which is why --always-keep
exists.

Run:
    python -m pytest workflow/scripts/10_decontam/test_select_contigs_by_taxonomy.py -q
"""

import subprocess
import sys
from pathlib import Path


# The script under test sits next to this file, wherever pytest is started from.
SCRIPT = Path(__file__).resolve().with_name("select_contigs_by_taxonomy.py")


# ── The toy assembly ─────────────────────────────────────────────────────────
# One row per contig: its ID, the genus BlobTools reports for it, the full FASTA
# header, and the sequence. Each sequence is made of one repeated pattern, so a
# contig is easy to recognise in an output file and its length is easy to check.
# contig_1 has a description after its ID (as many assemblers write); the script
# joins the BlobTools table and the FASTA on the first word only, but must keep
# the whole header in its output.
CONTIGS = [
    ("contig_1", "Fusarium",    "contig_1 length=120 cov=12.5", "A" * 120),
    ("contig_2", "Fusarium",    "contig_2",                     "C" * 100),
    ("contig_3", "Fusarium",    "contig_3",                     "G" * 80),
    ("contig_4", "Trichoderma", "contig_4",                     "T" * 60),
    ("contig_5", "no-hit",      "contig_5",                     "AC" * 20),
    ("contig_6", "Botrytis",    "contig_6",                     "GT" * 16),
]

# The IDs in table order, and the correct sequence of each, for the checks below.
ALL_IDS = [contig for contig, _, _, _ in CONTIGS]
SEQUENCES = {contig: sequence for contig, _, _, sequence in CONTIGS}

# The file names of the four outputs, modelled on the v1 fungal rule. The script
# writes to whatever paths it is given, so the names themselves are not tested.
LIST_NAME = "contigs.list"
FASTA_NAME = "assembly_decontam.fasta"
COMPOSITION_NAME = "toy_composition.txt"
DECISIONS_NAME = "contig_taxonomy_decisions.tsv"


# ── Building the two input files ─────────────────────────────────────────────

def make_blob_table(rows):
    """Return the text of a BlobTools table for a list of (contig, genus) rows.

    The real table (blobtools view --rank all) has 29 tab-separated columns: the
    contig name first, then length, GC and coverage, then four columns for each
    taxonomic rank. The script reads only column 1 (the contig) and column 22
    (the genus), so this fixture keeps that layout and fills every other column
    with a dash. Lines starting with # are the table's header, and the script
    must skip them.
    """
    lines = [
        "## 1.1.1",
        "## taxrule\t: bestsum",
        "# name\tlength\t(...)\tgenus.t.22%s\t(...)",
    ]
    for contig, genus in rows:
        fields = ["-"] * 29
        fields[0] = contig
        fields[1] = str(len(SEQUENCES[contig]))
        fields[21] = genus
        lines.append("\t".join(fields))
    return "\n".join(lines) + "\n"


def make_fasta():
    """Return the text of the toy assembly FASTA, wrapped at 50 columns.

    Wrapped on purpose: the script has to join the wrapped lines back into one
    sequence, and it writes its own output wrapped at 80 columns.
    """
    lines = []
    for _, _, header, sequence in CONTIGS:
        lines.append(">" + header)
        for start in range(0, len(sequence), 50):
            lines.append(sequence[start:start + 50])
    return "\n".join(lines) + "\n"


def write_always_keep(tmp_path, text):
    """Write a --always-keep file and return its path as text.

    In the workflow this file is confirmed_mitochondrial.ids, written by the
    mitochondrial screen: one contig ID per line.
    """
    path = tmp_path / "confirmed_mitochondrial.ids"
    path.write_text(text, encoding="utf-8")
    return str(path)


# ── Running the script and reading what it wrote ────────────────────────────

def run_selector(workdir, *flags, leave_out_of_table=()):
    """Write the inputs into workdir, run the script on them, return the process.

    workdir is a new folder for this run (it is created here, so every run in a
    test needs its own). The two inputs are written into it and the four outputs
    go to workdir/out/.

    flags are the policy options exactly as they would stand on the command line,
    for example "--mode", "exclude", "--exclude-genera", "Botrytis". The
    required options (input and output paths, sample name) are added here.

    leave_out_of_table gives a contig no row in the BlobTools table at all while
    it stays in the FASTA, for the test that needs a different toy assembly.

    The returned object has .returncode, .stdout and .stderr. The script
    fails by raising an error, which Python reports as exit status 1; a
    command-line error that argparse reports, such as --mode auto, gives 2.
    """
    workdir = Path(workdir)
    outputs = workdir / "out"
    outputs.mkdir(parents=True)

    rows = []
    for contig, genus, _, _ in CONTIGS:
        if contig in leave_out_of_table:
            continue
        rows.append((contig, genus))

    blob_path = workdir / "bestscore.blob.blobDB.table.txt"
    fasta_path = workdir / "assembly.fasta"
    blob_path.write_text(make_blob_table(rows), encoding="utf-8")
    fasta_path.write_text(make_fasta(), encoding="utf-8")

    command = [
        sys.executable, str(SCRIPT),
        "--bestscore", str(blob_path),
        "--contigs", str(fasta_path),
        "--sample", "toy",
        "--output-list", str(outputs / LIST_NAME),
        "--output-fasta", str(outputs / FASTA_NAME),
        "--composition", str(outputs / COMPOSITION_NAME),
        "--decisions", str(outputs / DECISIONS_NAME),
    ]
    command.extend(flags)
    # The output is decoded as UTF-8 explicitly, so the tests do not depend on the
    # language settings of the machine they run on.
    return subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
        errors="replace",
    )


def read_lines(path):
    """Return the lines of a text file, without their line endings."""
    return Path(path).read_text(encoding="utf-8").splitlines()


def read_fasta(path):
    """Read a FASTA file into a list of (header, sequence) pairs, in file order."""
    records = []
    header = None
    pieces = []
    for line in read_lines(path):
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(pieces)))
            header = line[1:]
            pieces = []
        else:
            pieces.append(line)
    if header is not None:
        records.append((header, "".join(pieces)))
    return records


def read_decisions(workdir):
    """Read the audit table into {contig: (assigned_genus, action, reason)}.

    Also checks that the table has exactly the four documented columns, since
    other steps and people read it by those names.
    """
    lines = read_lines(Path(workdir) / "out" / DECISIONS_NAME)
    assert lines[0] == "contig\tassigned_genus\taction\treason"
    decisions = {}
    for line in lines[1:]:
        contig, genus, action, reason = line.split("\t")
        decisions[contig] = (genus, action, reason)
    return decisions


def kept_ids(workdir):
    """Return the IDs of the kept contigs, after checking the outputs agree.

    The kept-ID list, the kept FASTA and the "keep" rows of the audit table are
    three views of one decision, and different steps read different ones, so they
    must never disagree. Every kept sequence must also be unchanged.
    """
    outputs = Path(workdir) / "out"
    from_list = read_lines(outputs / LIST_NAME)

    from_fasta = []
    for header, sequence in read_fasta(outputs / FASTA_NAME):
        contig_id = header.split()[0]
        from_fasta.append(contig_id)
        assert sequence == SEQUENCES[contig_id]

    from_decisions = []
    for contig, row in read_decisions(workdir).items():
        if row[1] == "keep":
            from_decisions.append(contig)

    assert from_fasta == from_list
    assert from_decisions == from_list
    return from_list


def without(decisions, contig):
    """Return a copy of a decisions dict with one contig left out.

    Used to show that protecting one contig changed nothing about the others.
    """
    others = dict(decisions)
    del others[contig]
    return others


# ── Mode off: keep everything ────────────────────────────────────────────────
# Off is the fungal default. Nothing is filtered, not even a contig BLAST could
# not place, but the audit files are still written so every run leaves the same
# trail whatever the mode.

def test_mode_off_keeps_every_contig(tmp_path):
    proc = run_selector(tmp_path, "--mode", "off")
    assert proc.returncode == 0, proc.stderr

    decisions = read_decisions(tmp_path)
    assert list(decisions) == ALL_IDS      # one row per contig, in table order
    for contig, genus, _, _ in CONTIGS:
        assert decisions[contig] == (genus, "keep", "mode_off")
    assert kept_ids(tmp_path) == ALL_IDS


def test_mode_off_ignores_discard_no_hit(tmp_path):
    # discard_no_hit is true in the fungal config. It must not matter when the
    # mode is off: the unplaced contig_5 is kept, because off bypasses all rules.
    proc = run_selector(tmp_path, "--mode", "off", "--discard-no-hit", "true")
    assert proc.returncode == 0, proc.stderr
    assert read_decisions(tmp_path)["contig_5"] == ("no-hit", "keep", "mode_off")
    assert kept_ids(tmp_path) == ALL_IDS


def test_without_the_mode_flag_the_mode_is_off(tmp_path):
    # Rule select_contigs always passes --mode, but a run by hand may leave it
    # out. The script then falls back to off, the FunFlux default (the BacFlux
    # copy falls back to auto), and keeps every contig.
    proc = run_selector(tmp_path)
    assert proc.returncode == 0, proc.stderr

    decisions = read_decisions(tmp_path)
    for contig, genus, _, _ in CONTIGS:
        assert decisions[contig] == (genus, "keep", "mode_off")
    assert kept_ids(tmp_path) == ALL_IDS


# ── The output files: FASTA layout and composition report ───────────────────

def test_kept_fasta_keeps_the_full_header_and_wraps_at_80_columns(tmp_path):
    proc = run_selector(tmp_path, "--mode", "off")
    assert proc.returncode == 0, proc.stderr

    lines = read_lines(tmp_path / "out" / FASTA_NAME)
    # The header keeps its description, not just the first word.
    assert lines[0] == ">contig_1 length=120 cov=12.5"
    # The input was wrapped at 50 columns; the output is wrapped at 80, so a
    # 120 bp contig becomes one line of 80 and one of 40.
    assert lines[1] == "A" * 80
    assert lines[2] == "A" * 40
    assert lines[3] == ">contig_2"


def test_composition_report_gives_bases_and_contigs_per_genus(tmp_path):
    # Each line gives a genus's share of the DNA ("bases") and of the contig
    # count ("contigs"), sorted by DNA. The two can disagree, which is exactly
    # what the report is there to show. Here Fusarium holds 300 of the 432 bp
    # and 3 of the 6 contigs. The report describes everything BlobTools saw,
    # whatever is kept, so the mode does not matter.
    proc = run_selector(tmp_path, "--mode", "off")
    assert proc.returncode == 0, proc.stderr
    assert read_lines(tmp_path / "out" / COMPOSITION_NAME) == [
        "Fusarium: bases 0.69  contigs 0.50",
        "Trichoderma: bases 0.14  contigs 0.17",
        "no-hit: bases 0.09  contigs 0.17",
        "Botrytis: bases 0.07  contigs 0.17",
    ]


# ── Mode auto: refused, with an explanation ──────────────────────────────────
# BacFlux has a fourth mode, auto, which keeps the genus with the most contigs.
# FunFlux does not: on a poorly sequenced fungus it can remove most of the
# genome (the script's docstring explains how). A request for auto, on the
# command line or in the override table, must stop the run before anything is
# written, with one message that says why and what to use instead (from the
# table, after the sample and the file name).

# The exact message (AUTO_MODE_REFUSAL in the script). It is copied here on
# purpose, so that a change of wording fails these tests and is made deliberately.
EXPECTED_AUTO_REFUSAL = (
    "decontamination mode 'auto' is not available in FunFlux: BLAST genus "
    "labels are unreliable for poorly sequenced fungi, so auto can remove most "
    "of the genome. Use off (the default), read "
    "contig_taxonomy_decisions.tsv, then exclude the contaminant genera or "
    "include the genera that BLAST assigned to the isolate's own contigs."
)


def test_mode_auto_on_the_command_line_is_refused_with_the_explanation(tmp_path):
    # argparse, which reads the command line, reports the refusal as an error on
    # the --mode argument: the last line of stderr is
    # "<script>: error: argument --mode: <message>". Without the explicit check
    # it would only say "invalid choice: 'auto'". "AUTO" in capitals is refused
    # the same way; each spelling runs in its own folder.
    for folder, spelling in [("lower_case", "auto"), ("capitals", "AUTO")]:
        proc = run_selector(tmp_path / folder, "--mode", spelling)
        assert proc.returncode != 0

        last_line = proc.stderr.strip().splitlines()[-1]
        assert last_line.endswith("error: argument --mode: " + EXPECTED_AUTO_REFUSAL)
        # The refusal comes before any input is read, so no output file exists.
        assert list((tmp_path / folder / "out").iterdir()) == []


def test_mode_auto_in_the_override_table_is_refused_with_the_same_message(tmp_path):
    # The run-wide mode is a valid one (off), but the override table asks for
    # auto. Two tables, each run in its own folder:
    #   own_row    the auto row is this sample's ("toy"), and the cell says
    #              "Auto": the script reads the mode column in lower case, so
    #              capitals must not slip past the refusal;
    #   other_row  the auto row belongs to another sample. It still stops this
    #              sample's job, because the script checks the mode of every
    #              row, so the message must say which row it means.
    cases = [("own_row", "toy", "Auto"), ("other_row", "OTHER", "auto")]
    for folder, row_sample, mode_cell in cases:
        overrides = tmp_path / f"{folder}_overrides.tsv"
        overrides.write_text(
            "sample\tmode\tinclude_genera\texclude_genera\tdiscard_no_hit\n"
            f"{row_sample}\t{mode_cell}\t\t\t\n",
            encoding="utf-8",
        )
        proc = run_selector(
            tmp_path / folder, "--mode", "off", "--sample-overrides", str(overrides)
        )
        assert proc.returncode != 0

        # The same message as on the command line, raised here as a ValueError,
        # with the sample of the row and the table's path in front of it.
        last_line = proc.stderr.strip().splitlines()[-1]
        assert last_line == (
            f"ValueError: Sample '{row_sample}' in override file '{overrides}': "
            + EXPECTED_AUTO_REFUSAL
        )
        # The override table is read before the BlobTools table, so again no
        # output file exists.
        assert list((tmp_path / folder / "out").iterdir()) == []


# ── Mode include: keep only the listed genera ────────────────────────────────

def test_include_mode_keeps_only_the_listed_genera(tmp_path):
    # Two genera are listed (separated by a semicolon, one of the accepted
    # separators); contigs of any other genus, and unplaced ones, go.
    proc = run_selector(
        tmp_path, "--mode", "include", "--include-genera", "Fusarium;Trichoderma"
    )
    assert proc.returncode == 0, proc.stderr

    decisions = read_decisions(tmp_path)
    for contig in ["contig_1", "contig_2", "contig_3"]:
        assert decisions[contig] == ("Fusarium", "keep", "included_genus")
    assert decisions["contig_4"] == ("Trichoderma", "keep", "included_genus")
    assert decisions["contig_5"] == ("no-hit", "remove", "not_included_genus")
    assert decisions["contig_6"] == ("Botrytis", "remove", "not_included_genus")
    assert kept_ids(tmp_path) == ["contig_1", "contig_2", "contig_3", "contig_4"]


# ── Mode exclude: drop only the listed genera ────────────────────────────────

def test_exclude_mode_removes_only_the_listed_genera(tmp_path):
    # The inverse of include: only Trichoderma goes. Everything else stays,
    # including the unplaced contig_5 (discard_no_hit is false here) and the
    # Botrytis contig_6.
    proc = run_selector(tmp_path, "--mode", "exclude", "--exclude-genera", "Trichoderma")
    assert proc.returncode == 0, proc.stderr

    decisions = read_decisions(tmp_path)
    assert decisions["contig_4"] == ("Trichoderma", "remove", "excluded_genus")
    for contig in ["contig_1", "contig_2", "contig_3", "contig_5", "contig_6"]:
        assert decisions[contig][1:] == ("keep", "not_excluded_genus")
    assert kept_ids(tmp_path) == [
        "contig_1", "contig_2", "contig_3", "contig_5", "contig_6"
    ]


# ── discard_no_hit: what happens to contigs BLAST could not place ───────────

def test_discard_no_hit_decides_the_fate_of_unplaced_contigs(tmp_path):
    # Exclude mode shows the difference most plainly, because the unplaced
    # contig_5 is not on the exclude list. With discard_no_hit false it stays.
    kept = run_selector(
        tmp_path / "false", "--mode", "exclude", "--exclude-genera", "Trichoderma",
        "--discard-no-hit", "false",
    )
    assert kept.returncode == 0, kept.stderr
    assert read_decisions(tmp_path / "false")["contig_5"] == (
        "no-hit", "keep", "not_excluded_genus"
    )
    assert "contig_5" in kept_ids(tmp_path / "false")

    # With discard_no_hit true it is removed, and the reason says why.
    removed = run_selector(
        tmp_path / "true", "--mode", "exclude", "--exclude-genera", "Trichoderma",
        "--discard-no-hit", "true",
    )
    assert removed.returncode == 0, removed.stderr
    decisions = read_decisions(tmp_path / "true")
    assert decisions["contig_5"] == ("no-hit", "remove", "discarded_no_hit")
    assert decisions["contig_4"] == ("Trichoderma", "remove", "excluded_genus")
    assert kept_ids(tmp_path / "true") == ["contig_1", "contig_2", "contig_3", "contig_6"]


def test_discard_no_hit_true_applies_in_include_mode_too(tmp_path):
    # The unplaced contig is removed in include mode either way (no-hit is not on
    # the include list), but the reason differs: with discard_no_hit true it is
    # "discarded_no_hit", checked before the genus rules of the mode are applied.
    include = run_selector(
        tmp_path / "include", "--mode", "include", "--include-genera", "Fusarium",
        "--discard-no-hit", "true",
    )
    assert include.returncode == 0, include.stderr
    assert read_decisions(tmp_path / "include")["contig_5"] == (
        "no-hit", "remove", "discarded_no_hit"
    )


# ── --always-keep: protect the confirmed mitochondrial contig ────────────────
# Each rescue test runs the same policy twice. The first run has no protection and
# shows that the filter really would remove the contig; the second lists the
# contig in a --always-keep file and shows that it now survives, with its own
# reason in the audit table and with every other decision unchanged.

def test_always_keep_rescues_a_contig_that_exclude_mode_would_remove(tmp_path):
    # The user excludes Botrytis, and contig_6 is labelled Botrytis.
    policy = ["--mode", "exclude", "--exclude-genera", "Botrytis"]

    plain = run_selector(tmp_path / "plain", *policy)
    assert plain.returncode == 0, plain.stderr
    plain_decisions = read_decisions(tmp_path / "plain")
    assert plain_decisions["contig_6"] == ("Botrytis", "remove", "excluded_genus")

    ids_file = write_always_keep(tmp_path, "contig_6\n")
    protected = run_selector(tmp_path / "protected", *policy, "--always-keep", ids_file)
    assert protected.returncode == 0, protected.stderr
    protected_decisions = read_decisions(tmp_path / "protected")
    assert protected_decisions["contig_6"] == (
        "Botrytis", "keep", "protected_confirmed_mitochondrial"
    )
    assert "contig_6" in kept_ids(tmp_path / "protected")
    assert without(protected_decisions, "contig_6") == without(plain_decisions, "contig_6")


def test_always_keep_rescues_a_contig_that_include_mode_would_remove(tmp_path):
    # Include keeps only Fusarium, so contig_6 would go as "not included".
    policy = ["--mode", "include", "--include-genera", "Fusarium"]

    plain = run_selector(tmp_path / "plain", *policy)
    assert plain.returncode == 0, plain.stderr
    plain_decisions = read_decisions(tmp_path / "plain")
    assert plain_decisions["contig_6"] == ("Botrytis", "remove", "not_included_genus")

    ids_file = write_always_keep(tmp_path, "contig_6\n")
    protected = run_selector(tmp_path / "protected", *policy, "--always-keep", ids_file)
    assert protected.returncode == 0, protected.stderr
    protected_decisions = read_decisions(tmp_path / "protected")
    assert protected_decisions["contig_6"] == (
        "Botrytis", "keep", "protected_confirmed_mitochondrial"
    )
    assert without(protected_decisions, "contig_6") == without(plain_decisions, "contig_6")


def test_always_keep_rescues_a_no_hit_contig_from_discard_no_hit(tmp_path):
    # A mitochondrion that BLAST cannot place at all is the other common case.
    # discard_no_hit true would remove it; protection has to win over that too.
    policy = [
        "--mode", "exclude", "--exclude-genera", "Trichoderma",
        "--discard-no-hit", "true",
    ]

    plain = run_selector(tmp_path / "plain", *policy)
    assert plain.returncode == 0, plain.stderr
    plain_decisions = read_decisions(tmp_path / "plain")
    assert plain_decisions["contig_5"] == ("no-hit", "remove", "discarded_no_hit")

    ids_file = write_always_keep(tmp_path, "contig_5\n")
    protected = run_selector(tmp_path / "protected", *policy, "--always-keep", ids_file)
    assert protected.returncode == 0, protected.stderr
    protected_decisions = read_decisions(tmp_path / "protected")
    assert protected_decisions["contig_5"] == (
        "no-hit", "keep", "protected_confirmed_mitochondrial"
    )
    assert "contig_5" in kept_ids(tmp_path / "protected")
    assert without(protected_decisions, "contig_5") == without(plain_decisions, "contig_5")


def test_always_keep_in_mode_off_records_the_protected_reason(tmp_path):
    # In mode off everything is kept anyway, so protection cannot change the
    # kept sequences (the FASTA is identical with and without it). It still
    # shows in the audit table: the protected contig says why it is special,
    # instead of "mode_off", and the others are untouched.
    plain = run_selector(tmp_path / "plain", "--mode", "off")
    assert plain.returncode == 0, plain.stderr

    ids_file = write_always_keep(tmp_path, "contig_6\n")
    protected = run_selector(
        tmp_path / "protected", "--mode", "off", "--always-keep", ids_file
    )
    assert protected.returncode == 0, protected.stderr

    decisions = read_decisions(tmp_path / "protected")
    assert decisions["contig_6"] == (
        "Botrytis", "keep", "protected_confirmed_mitochondrial"
    )
    for contig in ["contig_1", "contig_2", "contig_3", "contig_4", "contig_5"]:
        assert decisions[contig][1:] == ("keep", "mode_off")
    assert kept_ids(tmp_path / "protected") == ALL_IDS
    assert (tmp_path / "protected" / "out" / FASTA_NAME).read_bytes() == (
        tmp_path / "plain" / "out" / FASTA_NAME
    ).read_bytes()


def test_a_protected_contig_missing_from_the_blob_table_is_still_kept(tmp_path):
    # BlobTools has no row for contig_6 (say, it was never scored), but the
    # assembly FASTA still holds it. Without protection it is simply lost,
    # because the script only keeps contigs the table lists. The policy keeps
    # Fusarium only, so the kept list below shows the protected contig added
    # after the three Fusarium contigs.
    policy = ["--mode", "include", "--include-genera", "Fusarium"]

    plain = run_selector(
        tmp_path / "plain", *policy, leave_out_of_table=["contig_6"]
    )
    assert plain.returncode == 0, plain.stderr
    assert "contig_6" not in read_decisions(tmp_path / "plain")
    assert "contig_6" not in kept_ids(tmp_path / "plain")

    # With protection it is kept and written to the audit table with the
    # placeholder genus "no-taxonomy", after the contigs the table did list.
    ids_file = write_always_keep(tmp_path, "contig_6\n")
    protected = run_selector(
        tmp_path / "protected", *policy, "--always-keep", ids_file,
        leave_out_of_table=["contig_6"],
    )
    assert protected.returncode == 0, protected.stderr
    decisions = read_decisions(tmp_path / "protected")
    assert decisions["contig_6"] == (
        "no-taxonomy", "keep", "protected_confirmed_mitochondrial"
    )
    assert list(decisions)[-1] == "contig_6"
    assert kept_ids(tmp_path / "protected") == [
        "contig_1", "contig_2", "contig_3", "contig_6"
    ]


def test_a_protected_id_missing_from_the_fasta_fails_with_a_clear_message(tmp_path):
    # contig_94 to contig_99 exist nowhere, contig_6 is real. They are listed out
    # of order on purpose: the message names only the missing IDs, in sorted order
    # (the same on every run), so the mismatch between the mitochondrial list and
    # the assembly is easy to read. Six missing IDs, not two, because the script
    # holds them in a set, whose order changes from run to run; with six, an
    # unsorted message would almost never come out sorted by luck.
    ids_text = (
        "contig_97\ncontig_99\ncontig_6\ncontig_95\ncontig_98\ncontig_94\ncontig_96\n"
    )
    ids_file = write_always_keep(tmp_path, ids_text)
    proc = run_selector(tmp_path, "--mode", "off", "--always-keep", ids_file)

    assert proc.returncode != 0
    last_line = proc.stderr.strip().splitlines()[-1]
    assert last_line == (
        "ValueError: Protected mitochondrial contig(s) were not found in FASTA: "
        "contig_94, contig_95, contig_96, contig_97, contig_98, contig_99"
    )


def test_always_keep_file_may_hold_comments_blank_lines_and_stray_spaces(tmp_path):
    # The file format: one ID per line. Blank lines and lines starting with # are
    # skipped (an indented comment too), and spaces around an ID are trimmed.
    # If the comment text were mistaken for an ID, the run would fail because no
    # such contig exists in the FASTA, so a successful run proves it was skipped.
    ids_text = (
        "# confirmed mitochondrial contigs, one ID per line\n"
        "\n"
        "contig_6\n"
        "   contig_5   \n"
        "    # an indented comment is still a comment\n"
    )
    ids_file = write_always_keep(tmp_path, ids_text)
    proc = run_selector(
        tmp_path / "run", "--mode", "exclude", "--exclude-genera", "Botrytis",
        "--discard-no-hit", "true", "--always-keep", ids_file,
    )
    assert proc.returncode == 0, proc.stderr

    # contig_6 would be excluded and contig_5 discarded as no-hit; both survive.
    decisions = read_decisions(tmp_path / "run")
    assert decisions["contig_6"] == (
        "Botrytis", "keep", "protected_confirmed_mitochondrial"
    )
    assert decisions["contig_5"] == ("no-hit", "keep", "protected_confirmed_mitochondrial")
    assert kept_ids(tmp_path / "run") == ALL_IDS


def test_an_empty_or_bare_always_keep_flag_protects_nothing(tmp_path):
    # When there is nothing to protect, the Snakemake rule still writes the flag,
    # and the value reaches the script as an empty word. A bare flag with no value
    # must behave the same. Neither may protect or break anything: contig_6 is
    # still removed by the exclude list, and the run succeeds.
    policy = ["--mode", "exclude", "--exclude-genera", "Botrytis"]

    empty_word = run_selector(tmp_path / "empty_word", *policy, "--always-keep", "")
    assert empty_word.returncode == 0, empty_word.stderr
    assert read_decisions(tmp_path / "empty_word")["contig_6"] == (
        "Botrytis", "remove", "excluded_genus"
    )

    bare_flag = run_selector(tmp_path / "bare_flag", *policy, "--always-keep")
    assert bare_flag.returncode == 0, bare_flag.stderr
    assert read_decisions(tmp_path / "bare_flag")["contig_6"] == (
        "Botrytis", "remove", "excluded_genus"
    )


def test_the_snakemake_style_command_line_with_every_optional_flag_empty(tmp_path):
    # The rule select_contigs always passes all seven optional flags. Those not
    # set in the config reach the script as empty words, and the script must read
    # an empty word as "not configured" rather than stop with an argument error.
    proc = run_selector(
        tmp_path, "--mode", "off",
        "--include-genera", "",
        "--include-genera-by-sample", "",
        "--exclude-genera", "",
        "--exclude-genera-file", "",
        "--sample-overrides", "",
        "--discard-no-hit", "false",
        "--always-keep", "",
    )
    assert proc.returncode == 0, proc.stderr
    assert kept_ids(tmp_path) == ALL_IDS
