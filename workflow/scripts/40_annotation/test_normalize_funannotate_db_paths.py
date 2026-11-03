"""Tests for normalize_funannotate_db_paths.py.

What the script does, in one sentence: it rewrites the wrong absolute paths that
a downloaded funannotate database remembers (in trained_species/*/info.json and
in funannotate-db-info.txt) so that they point into the database folder on this
computer.

That is an edit INSIDE a shared database, so these tests pin down the promises
that matter most:
  * wrong paths are rewritten, and paths that are already right are left alone;
  * the counts in the summary are right;
  * a second run changes nothing (the script is "idempotent": repeating it is
    harmless);
  * a file that needs no change is not touched at all: same bytes, same
    modification time;
  * a problem (missing folder, missing file) stops the script with a clear
    message and writes no marker file;
  * the command-line wrapper writes the marker JSON and the log line exactly as
    the old Snakemake rule did.

Every test builds a tiny fake database in pytest's temporary folder (tmp_path),
so no real funannotate data and no network are needed.

Run: pytest workflow/scripts/40_annotation/test_normalize_funannotate_db_paths.py
"""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import normalize_funannotate_db_paths as nf


# ---------------------------------------------------------------------------
# Constants used to build the fake database
# ---------------------------------------------------------------------------

# Where the real snapshot was "built". None of these folders exists on the test
# computer, which is exactly the situation after downloading the real snapshot.
BUILDER_ROOT = "/home/builder/funannotate_db"
OTHER_ROOT = "/scratch/another_machine/db_copy"

# The gene finders listed in every info.json (same order as in the real snapshot).
GENE_FINDERS = ["augustus", "genemark", "codingquarry", "snap", "glimmerhmm"]

# Every file of the fake database is given this old modification time
# (1 January 2020, as seconds since 1970). A file that the script rewrites, even
# with identical content, would get today's time instead, so "was this file
# touched?" can be answered reliably.
OLD_TIME = 1577836800

# The script file itself, for the tests that start it like a real command.
SCRIPT = Path(__file__).with_name("normalize_funannotate_db_paths.py")


# ---------------------------------------------------------------------------
# Helper functions that build the fake database and look at it
# ---------------------------------------------------------------------------

def stamp_old_modification_times(folder):
    """Give every file under folder the old modification time (see OLD_TIME)."""
    for path in folder.rglob("*"):
        if path.is_file():
            os.utime(path, (OLD_TIME, OLD_TIME))


def snapshot(folder):
    """Record every file under folder as {relative path: (content, modification time)}.

    Comparing two snapshots shows which files were changed, added or touched.
    """
    record = {}
    for path in sorted(folder.rglob("*")):
        if path.is_file():
            relative_name = path.relative_to(folder).as_posix()
            record[relative_name] = (path.read_bytes(), path.stat().st_mtime_ns)
    return record


def read_json(path):
    """Read a JSON file into Python objects."""
    return json.loads(path.read_text(encoding="utf-8"))


def make_species(database, species, stored_paths, compact=False):
    """Create trained_species/<species>/ with an info.json and the folders it names.

    stored_paths  maps a gene finder (for example "augustus") to the path TEXT
                  written into info.json. This is where a wrong path, from
                  another computer, is planted. For every gene finder in it the
                  real folder is created under the database, because the script
                  checks that the repaired path exists.
    compact       True writes the JSON on one line instead of indented. The
                  script writes JSON with 2-space indentation, so a compact file
                  that comes through a run byte for byte unchanged proves that
                  it was not rewritten.
    """
    species_folder = database / "trained_species" / species
    species_folder.mkdir(parents=True)

    info = {}
    for gene_finder in GENE_FINDERS:
        if gene_finder in stored_paths:
            (species_folder / gene_finder).mkdir()
            info[gene_finder] = [
                {
                    "version": "fake 1.0",
                    "source": "test fixture",
                    "date": "2025-12-15",
                    "path": stored_paths[gene_finder],
                }
            ]
        else:
            # A gene finder that was not trained for this species has an empty entry.
            info[gene_finder] = [{}]

    info_path = species_folder / "info.json"
    if compact:
        info_path.write_text(json.dumps(info), encoding="utf-8")
    else:
        info_path.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    return info_path


def index_line(name, kind, path):
    """One line of funannotate-db-info.txt.

    The real columns are: name, type, PATH, version, date, number of records, md5.
    """
    columns = [
        name,
        kind,
        path,
        "1.0",
        "2025-12-15",
        "1000",
        "0123456789abcdef0123456789abcdef",
    ]
    return "\t".join(columns) + "\n"


def build_fake_database(tmp_path):
    """Build a small fake funannotate database with some wrong paths in it.

    Returns the database folder. What is inside:
      trained species
        Aspergillus_nidulans  1 wrong path (augustus)
        Fusarium_graminearum  2 wrong paths (augustus and snap)
        Neurospora_crassa     already correct, JSON on one line (compact)
      funannotate-db-info.txt
        merops, uniprot, pfam, busco_outgroups   wrong paths (4)
        go                                       already correct
        an empty line and a two-column line     not entries; must survive as they are
    So a first run must report 3 species files, 3 species paths and 4 index paths.
    """
    # resolve(): the script resolves the folder too, so compare like with like.
    database = tmp_path.resolve() / "funannotate_db"

    make_species(
        database,
        "Aspergillus_nidulans",
        {"augustus": f"{BUILDER_ROOT}/trained_species/Aspergillus_nidulans/augustus"},
    )
    make_species(
        database,
        "Fusarium_graminearum",
        {
            "augustus": f"{OTHER_ROOT}/trained_species/Fusarium_graminearum/augustus",
            "snap": f"{OTHER_ROOT}/trained_species/Fusarium_graminearum/snap",
        },
    )
    correct_path = database / "trained_species" / "Neurospora_crassa" / "augustus"
    make_species(
        database,
        "Neurospora_crassa",
        {"augustus": str(correct_path)},
        compact=True,
    )

    # Small stand-ins for the large database files that the index points to.
    for file_name in ["merops.dmnd", "uniprot.dmnd", "Pfam-A.hmm", "go.obo"]:
        (database / file_name).write_text("stand-in for a big file\n", encoding="utf-8")
    (database / "outgroups").mkdir()

    index_text = (
        index_line("merops", "diamond", f"{BUILDER_ROOT}/merops.dmnd")
        + index_line("uniprot", "diamond", f"{OTHER_ROOT}/uniprot.dmnd")
        + index_line("pfam", "hmmer3", f"{BUILDER_ROOT}/Pfam-A.hmm")
        + index_line("busco_outgroups", "outgroups", f"{BUILDER_ROOT}/outgroups")
        + index_line("go", "text", str(database / "go.obo"))
        + "\n"
        + "not_an_entry\tonly_two_columns\n"
    )
    (database / "funannotate-db-info.txt").write_text(index_text, encoding="utf-8")

    stamp_old_modification_times(database)
    return database


def expected_index_text(database):
    """The index file as it must look after a correct run on the fake database."""
    return (
        index_line("merops", "diamond", str(database / "merops.dmnd"))
        + index_line("uniprot", "diamond", str(database / "uniprot.dmnd"))
        + index_line("pfam", "hmmer3", str(database / "Pfam-A.hmm"))
        + index_line("busco_outgroups", "outgroups", str(database / "outgroups"))
        + index_line("go", "text", str(database / "go.obo"))
        + "\n"
        + "not_an_entry\tonly_two_columns\n"
    )


def make_minimal_database(tmp_path):
    """Make the smallest folder the script accepts.

    It holds an empty trained_species/ folder and an empty funannotate-db-info.txt.
    """
    database = tmp_path.resolve() / "funannotate_db"
    (database / "trained_species").mkdir(parents=True)
    (database / "funannotate-db-info.txt").write_text("", encoding="utf-8")
    return database


# ---------------------------------------------------------------------------
# One stored path: _rebase_funannotate_db_path
# ---------------------------------------------------------------------------
# This function turns ONE path from another computer into the matching path
# under the database folder, or refuses if the file is not there.

def test_rebase_moves_an_absolute_path_from_another_computer(tmp_path):
    database = tmp_path.resolve()
    (database / "uniprot.dmnd").write_text("x", encoding="utf-8")

    result = nf._rebase_funannotate_db_path(
        "/home/builder/funannotate_db/uniprot.dmnd", database
    )

    assert result == str(database / "uniprot.dmnd")


def test_rebase_accepts_a_relative_path(tmp_path):
    database = tmp_path.resolve()
    (database / "uniprot.dmnd").write_text("x", encoding="utf-8")

    result = nf._rebase_funannotate_db_path("uniprot.dmnd", database)

    assert result == str(database / "uniprot.dmnd")


def test_rebase_prefers_the_longest_trailing_part_that_exists(tmp_path):
    # Both <database>/sub/list.txt and <database>/list.txt exist. The old path
    # ended in sub/list.txt, so that longer match must win over the bare name.
    database = tmp_path.resolve()
    (database / "sub").mkdir()
    (database / "sub" / "list.txt").write_text("x", encoding="utf-8")
    (database / "list.txt").write_text("x", encoding="utf-8")

    result = nf._rebase_funannotate_db_path("/old/machine/sub/list.txt", database)

    assert result == str(database / "sub" / "list.txt")


def test_rebase_falls_back_to_a_shorter_part_when_the_longer_one_is_missing(tmp_path):
    # Only <database>/list.txt exists, so /old/machine/sub/list.txt maps to it.
    database = tmp_path.resolve()
    (database / "list.txt").write_text("x", encoding="utf-8")

    result = nf._rebase_funannotate_db_path("/old/machine/sub/list.txt", database)

    assert result == str(database / "list.txt")


def test_rebase_keeps_the_whole_trained_species_suffix(tmp_path):
    # Two species both have a folder called "augustus". The species name is part
    # of the suffix that must be kept, otherwise the wrong folder could be chosen.
    database = tmp_path.resolve()
    (database / "trained_species" / "Species_a" / "augustus").mkdir(parents=True)
    (database / "trained_species" / "Species_b" / "augustus").mkdir(parents=True)

    result = nf._rebase_funannotate_db_path(
        "/old/machine/db/trained_species/Species_b/augustus", database
    )

    assert result == str(database / "trained_species" / "Species_b" / "augustus")


def test_rebase_rejects_an_empty_path(tmp_path):
    # "/" has no pieces left once the root is dropped.
    with pytest.raises(ValueError, match="Invalid empty Funannotate database path"):
        nf._rebase_funannotate_db_path("/", tmp_path)


def test_rebase_reports_a_missing_trained_species_folder(tmp_path):
    with pytest.raises(ValueError, match="Trained-species asset is missing"):
        nf._rebase_funannotate_db_path(
            "/old/machine/db/trained_species/Nowhere/augustus", tmp_path
        )


def test_rebase_reports_a_missing_database_file(tmp_path):
    with pytest.raises(ValueError, match="Database asset is missing") as problem:
        nf._rebase_funannotate_db_path("/old/machine/db/missing.dmnd", tmp_path)

    # The message names the stored path, so the user can see which entry is the problem.
    assert "/old/machine/db/missing.dmnd" in str(problem.value)


# ---------------------------------------------------------------------------
# The whole database: normalize_funannotate_db_paths
# ---------------------------------------------------------------------------

def test_trained_species_paths_are_rewritten(tmp_path):
    database = build_fake_database(tmp_path)

    nf.normalize_funannotate_db_paths(database)

    species_folder = database / "trained_species"
    aspergillus = read_json(species_folder / "Aspergillus_nidulans" / "info.json")
    fusarium = read_json(species_folder / "Fusarium_graminearum" / "info.json")
    neurospora = read_json(species_folder / "Neurospora_crassa" / "info.json")

    assert aspergillus["augustus"][0]["path"] == str(
        species_folder / "Aspergillus_nidulans" / "augustus"
    )
    assert fusarium["augustus"][0]["path"] == str(
        species_folder / "Fusarium_graminearum" / "augustus"
    )
    assert fusarium["snap"][0]["path"] == str(
        species_folder / "Fusarium_graminearum" / "snap"
    )
    # Already correct before the run, so still the same.
    assert neurospora["augustus"][0]["path"] == str(
        species_folder / "Neurospora_crassa" / "augustus"
    )


def test_rewriting_changes_only_the_path_field(tmp_path):
    database = build_fake_database(tmp_path)
    info_path = database / "trained_species" / "Aspergillus_nidulans" / "info.json"
    before = read_json(info_path)

    nf.normalize_funannotate_db_paths(database)

    after = read_json(info_path)
    # Put the old path back into the new content: everything else (version,
    # source, date, the empty entries of the other gene finders) must be equal.
    after["augustus"][0]["path"] = before["augustus"][0]["path"]
    assert after == before


def test_database_index_paths_are_rewritten(tmp_path):
    database = build_fake_database(tmp_path)

    nf.normalize_funannotate_db_paths(database)

    # expected_index_text has the new paths in column 3. Every other column, the
    # already-correct "go" line, the empty line and the two-column line are
    # exactly as before.
    index_text = (database / "funannotate-db-info.txt").read_text(encoding="utf-8")
    assert index_text == expected_index_text(database)


def test_summary_counts_are_right(tmp_path):
    database = build_fake_database(tmp_path)

    summary = nf.normalize_funannotate_db_paths(database)

    assert summary == {
        # Three info.json files were examined (changed or not).
        "metadata_files": 3,
        # Aspergillus 1 + Fusarium 2 + Neurospora 0 (already correct).
        "metadata_paths_updated": 3,
        # merops, uniprot, pfam and busco_outgroups; "go" was already correct.
        "database_info_paths_updated": 4,
    }


def test_only_the_files_with_wrong_paths_are_changed(tmp_path):
    database = build_fake_database(tmp_path)
    before = snapshot(database)

    nf.normalize_funannotate_db_paths(database)

    after = snapshot(database)
    # No file was added (so no temporary file is left behind) or removed ...
    assert sorted(after) == sorted(before)
    # ... and exactly these three files were rewritten. Everything else - the
    # already-correct Neurospora file and every database file - has the same
    # content and the same old modification time.
    changed_files = sorted(name for name in before if after[name] != before[name])
    assert changed_files == [
        "funannotate-db-info.txt",
        "trained_species/Aspergillus_nidulans/info.json",
        "trained_species/Fusarium_graminearum/info.json",
    ]


def test_second_run_reports_zero_changes_and_touches_nothing(tmp_path):
    database = build_fake_database(tmp_path)
    nf.normalize_funannotate_db_paths(database)

    # Give every file the old time again, so that a rewrite during the second
    # run could not hide behind the time stamp of the first run.
    stamp_old_modification_times(database)
    after_first_run = snapshot(database)

    second_summary = nf.normalize_funannotate_db_paths(database)

    assert second_summary == {
        "metadata_files": 3,
        "metadata_paths_updated": 0,
        "database_info_paths_updated": 0,
    }
    assert snapshot(database) == after_first_run


def test_second_run_changes_nothing_when_the_database_is_reached_through_a_symlink(
    tmp_path,
):
    # The script follows symbolic links (it resolves the folder), so the paths
    # it writes are the real location. Running again through the same link must
    # still find nothing to change.
    database = build_fake_database(tmp_path)
    link = tmp_path / "database_link"
    link.symlink_to(database, target_is_directory=True)

    first_summary = nf.normalize_funannotate_db_paths(link)
    stamp_old_modification_times(database)
    after_first_run = snapshot(database)
    second_summary = nf.normalize_funannotate_db_paths(link)

    assert first_summary["metadata_paths_updated"] == 3
    assert second_summary["metadata_paths_updated"] == 0
    assert second_summary["database_info_paths_updated"] == 0
    assert snapshot(database) == after_first_run


def test_a_correct_index_is_not_rewritten_even_if_its_layout_is_unusual(tmp_path):
    database = build_fake_database(tmp_path)
    nf.normalize_funannotate_db_paths(database)

    # An index whose only line is correct but has NO final newline. If the script
    # rewrote it, it would add the newline, so unchanged bytes prove "not rewritten".
    index_path = database / "funannotate-db-info.txt"
    odd_but_correct = index_line("go", "text", str(database / "go.obo")).rstrip("\n")
    index_path.write_text(odd_but_correct, encoding="utf-8")
    stamp_old_modification_times(database)
    before = (index_path.read_bytes(), index_path.stat().st_mtime_ns)

    summary = nf.normalize_funannotate_db_paths(database)

    assert summary["database_info_paths_updated"] == 0
    assert (index_path.read_bytes(), index_path.stat().st_mtime_ns) == before


def test_entries_without_a_usable_path_are_ignored(tmp_path):
    database = make_minimal_database(tmp_path)
    species_folder = database / "trained_species" / "Odd_species"
    species_folder.mkdir()
    odd_info = {
        # A value that is not a list is skipped.
        "note": "text, not a list",
        # In a list, only a dictionary with a non-empty "path" counts. None of
        # these five is one: an empty dictionary, one without "path", one with
        # an empty path, a bare string and a number.
        "augustus": [{}, {"version": "no path here"}, {"path": ""}, "a string", 42],
    }
    info_path = species_folder / "info.json"
    info_path.write_text(json.dumps(odd_info), encoding="utf-8")
    before = info_path.read_bytes()

    summary = nf.normalize_funannotate_db_paths(database)

    assert summary == {
        "metadata_files": 1,
        "metadata_paths_updated": 0,
        "database_info_paths_updated": 0,
    }
    assert info_path.read_bytes() == before


def test_missing_database_directory_fails_clearly(tmp_path):
    missing = tmp_path.resolve() / "no_such_database"

    with pytest.raises(ValueError) as problem:
        nf.normalize_funannotate_db_paths(missing)

    # The message says what is missing and names the place that was searched.
    assert "trained-species directory is missing" in str(problem.value)
    assert str(missing) in str(problem.value)
    # Nothing was created as a side effect.
    assert not missing.exists()


def test_database_without_a_trained_species_folder_fails_clearly(tmp_path):
    # A folder that exists but is not a funannotate database (or an unpacking
    # that failed early) is reported the same way.
    empty_folder = tmp_path.resolve() / "not_a_database"
    empty_folder.mkdir()

    with pytest.raises(ValueError, match="trained-species directory is missing"):
        nf.normalize_funannotate_db_paths(empty_folder)


def test_missing_database_index_fails_clearly(tmp_path):
    database = tmp_path.resolve() / "funannotate_db"
    (database / "trained_species").mkdir(parents=True)

    with pytest.raises(ValueError, match="database index is missing") as problem:
        nf.normalize_funannotate_db_paths(database)

    assert "funannotate-db-info.txt" in str(problem.value)


def test_index_path_to_a_missing_file_fails_clearly(tmp_path):
    database = make_minimal_database(tmp_path)
    index_path = database / "funannotate-db-info.txt"
    index_path.write_text(
        index_line("uniprot", "diamond", f"{BUILDER_ROOT}/uniprot.dmnd"),
        encoding="utf-8",
    )
    before = index_path.read_bytes()

    with pytest.raises(ValueError, match="Database asset is missing"):
        nf.normalize_funannotate_db_paths(database)

    # The index was not changed by the failed run.
    assert index_path.read_bytes() == before


def test_rerun_after_an_error_finishes_the_job(tmp_path):
    database = build_fake_database(tmp_path)
    species_folder = database / "trained_species"
    index_path = database / "funannotate-db-info.txt"
    fusarium_info = species_folder / "Fusarium_graminearum" / "info.json"
    index_before = index_path.read_bytes()
    fusarium_before = fusarium_info.read_bytes()

    # Pretend the unpacking was incomplete: one trained-species folder is missing.
    missing_folder = species_folder / "Fusarium_graminearum" / "augustus"
    shutil.rmtree(missing_folder)

    with pytest.raises(ValueError, match="Trained-species asset is missing"):
        nf.normalize_funannotate_db_paths(database)

    # Species are handled in alphabetical order, so Aspergillus (before Fusarium)
    # was already repaired, and completely: the file is valid JSON with the new path.
    aspergillus = read_json(species_folder / "Aspergillus_nidulans" / "info.json")
    assert aspergillus["augustus"][0]["path"] == str(
        species_folder / "Aspergillus_nidulans" / "augustus"
    )
    # The file that caused the error and the index (done last) were not touched.
    assert fusarium_info.read_bytes() == fusarium_before
    assert index_path.read_bytes() == index_before

    # Repair the cause and run again: the rest is done, nothing is done twice.
    missing_folder.mkdir()
    summary = nf.normalize_funannotate_db_paths(database)

    assert summary == {
        "metadata_files": 3,
        "metadata_paths_updated": 2,        # only the two Fusarium paths were left
        "database_info_paths_updated": 4,
    }
    assert index_path.read_text(encoding="utf-8") == expected_index_text(database)


# ---------------------------------------------------------------------------
# The command-line wrapper: main()
# ---------------------------------------------------------------------------
# main() is what the Snakemake rule runs. It must write the marker JSON and the
# one-line log exactly as the old inline `run:` block did.

def run_main(database, output_folder):
    """Call main() with the three options; return (exit status, JSON path, log path)."""
    marker_name = "funannotate_db_paths.normalized.json"
    summary_path = output_folder / "00.resources" / marker_name
    log_path = output_folder / "logs" / "funannotate_db_paths.normalized.log"
    exit_status = nf.main(
        [
            "--database", str(database),
            "--summary-json", str(summary_path),
            "--log", str(log_path),
        ]
    )
    return exit_status, summary_path, log_path


def test_main_writes_the_marker_json_and_the_log_line(tmp_path):
    database = build_fake_database(tmp_path)

    exit_status, summary_path, log_path = run_main(database, tmp_path / "run_output")

    assert exit_status == 0
    expected_summary = {
        "metadata_files": 3,
        "metadata_paths_updated": 3,
        "database_info_paths_updated": 4,
        "database": str(database),
    }
    # json.dumps keeps the order of the keys it is given, so comparing the whole
    # text also checks the key order, the 2-space indent and the final newline.
    assert summary_path.read_text(encoding="utf-8") == (
        json.dumps(expected_summary, indent=2) + "\n"
    )
    assert log_path.read_text(encoding="utf-8") == (
        "[normalize_funannotate_db] normalized 3 trained-species path(s) "
        "across 3 metadata file(s), and 4 database-index path(s).\n"
    )
    # The temporary file used for the atomic write is gone.
    assert list((tmp_path / "run_output").rglob("*.tmp")) == []


def test_main_second_run_overwrites_the_marker_with_zero_counts(tmp_path):
    database = build_fake_database(tmp_path)
    output_folder = tmp_path / "run_output"
    run_main(database, output_folder)

    exit_status, summary_path, log_path = run_main(database, output_folder)

    assert exit_status == 0
    assert read_json(summary_path) == {
        "metadata_files": 3,
        "metadata_paths_updated": 0,
        "database_info_paths_updated": 0,
        "database": str(database),
    }
    assert log_path.read_text(encoding="utf-8") == (
        "[normalize_funannotate_db] normalized 0 trained-species path(s) "
        "across 3 metadata file(s), and 0 database-index path(s).\n"
    )


def test_main_accepts_relative_paths_and_creates_missing_folders(tmp_path, monkeypatch):
    build_fake_database(tmp_path)
    # Work from inside tmp_path so that all three options can be relative paths.
    monkeypatch.chdir(tmp_path)

    # marker.json is a bare file name (there is no folder to create for it);
    # logs/deeper/ does not exist yet and must be created.
    exit_status = nf.main(
        [
            "--database", "funannotate_db",
            "--summary-json", "marker.json",
            "--log", "logs/deeper/normalize.log",
        ]
    )

    assert exit_status == 0
    summary = read_json(tmp_path / "marker.json")
    # The summary records the ABSOLUTE path of the database, not the relative one.
    assert summary["database"] == str(tmp_path.resolve() / "funannotate_db")
    assert (tmp_path / "logs" / "deeper" / "normalize.log").is_file()


def test_main_fails_cleanly_when_the_database_is_missing(tmp_path, capsys):
    exit_status, summary_path, log_path = run_main(
        tmp_path / "no_such_database", tmp_path / "run_output"
    )

    assert exit_status == 1
    error_text = capsys.readouterr().err
    assert error_text.startswith("ERROR: ")
    assert "trained-species directory is missing" in error_text
    # Neither the marker nor the log was written, and no output folder was made,
    # so Snakemake sees a failed rule and starts nothing that waits for the marker.
    assert not summary_path.exists()
    assert not log_path.exists()
    assert not (tmp_path / "run_output").exists()


def test_main_requires_all_three_options():
    with pytest.raises(SystemExit) as problem:
        nf.main(["--database", "somewhere"])

    # 2 is argparse's exit status for a command-line mistake.
    assert problem.value.code == 2


# ---------------------------------------------------------------------------
# Starting the script like a real command
# ---------------------------------------------------------------------------

def run_script(*arguments):
    """Start the script in a separate Python process, as Snakemake does."""
    command = [sys.executable, str(SCRIPT)] + list(arguments)
    return subprocess.run(command, capture_output=True, text=True)


def test_script_runs_as_a_command_and_reports_failure_with_exit_status_1(tmp_path):
    database = build_fake_database(tmp_path)
    summary_path = tmp_path / "ok_run" / "marker.json"
    log_path = tmp_path / "ok_run" / "normalize.log"

    finished = run_script(
        "--database", str(database),
        "--summary-json", str(summary_path),
        "--log", str(log_path),
    )

    assert finished.returncode == 0, finished.stderr
    assert read_json(summary_path)["metadata_paths_updated"] == 3
    assert log_path.is_file()

    failed = run_script(
        "--database", str(tmp_path / "no_such_database"),
        "--summary-json", str(tmp_path / "bad_run" / "marker.json"),
        "--log", str(tmp_path / "bad_run" / "normalize.log"),
    )

    assert failed.returncode == 1
    assert "trained-species directory is missing" in failed.stderr
    assert not (tmp_path / "bad_run").exists()
