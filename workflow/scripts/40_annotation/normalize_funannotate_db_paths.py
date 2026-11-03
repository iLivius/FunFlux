#!/usr/bin/env python3
r"""Repair the file paths that a downloaded funannotate database remembers.

************************************************************************
* WARNING: THIS SCRIPT REWRITES FILES INSIDE THE FUNANNOTATE DATABASE. *
************************************************************************
The folder given with --database is normally a SHARED, long-lived resource: one
copy is used by many workflow runs, often by many users. Wherever a stored path
is wrong, this script overwrites the small metadata file that holds it, in
place, and keeps no copy of the old (wrong) content. Only the metadata files
described below are ever changed. The big data files (Pfam, UniProt, ...) are
never opened; they are only checked for existence. A rewritten file is a new
file: it belongs to the user who ran the script, and gets that user's default
permissions and a new modification time. There is no lock on the database, so
run one workflow at a time against a database that has not been repaired yet
(see "One run at a time" below).

Why this is needed
------------------
funannotate (the program that predicts the genes and annotates them) needs a
large set of reference data: protein and domain databases (UniProt, Pfam,
dbCAN, MEROPS, MIBiG, InterPro, ...) plus gene-finder settings trained for
about 166 species. funannotate can build this database itself (with
`funannotate setup`), but the workflow downloads a frozen copy of it instead:
one database snapshot published on Zenodo (links.funannotate_link in the
config). Two kinds of small text files inside the snapshot record where things
are, as ABSOLUTE paths: those of the computer that built the snapshot, or a
placeholder such as /path/to/funannotate_db (written by
miscellaneous/build_funannotate_db.sh):

  funannotate-db-info.txt
      The database index: one tab-separated line per database, with name,
      type, PATH, version, date, number of records and md5 checksum. Only the
      third column (the path) is rewritten.
  trained_species/<Species>/info.json
      One file per trained species. For each gene finder (augustus, genemark,
      codingquarry, snap, glimmerhmm) it holds a list of entries. The AUGUSTUS
      entry has a "path" to the folder with that species' trained AUGUSTUS
      parameters (trained_species/<Species>/augustus). A gene finder that was
      not trained for the species has an empty entry ({}), which is skipped.

Once the snapshot is unpacked on another computer, those paths point to folders
that do not exist there, and funannotate cannot find its own files.

What the script does
--------------------
For every stored path it works out where the same file is under --database and
writes that new path back. There are two cases (explained in front of the
functions below):
  * trained-species paths keep everything from "trained_species/" onwards:
        /home/builder/db/trained_species/Aspergillus_nidulans/augustus
     -> <database>/trained_species/Aspergillus_nidulans/augustus
  * database-index paths keep the longest end part of the old path that exists
    under --database:
        /home/builder/db/uniprot.dmnd  ->  <database>/uniprot.dmnd
The new path must exist. If it does not (the snapshot was only partly
downloaded or unpacked), the script stops with an error instead of writing a
path that cannot work.

Safe to run twice
-----------------
A second run changes nothing. A path that is already right is left alone, and a
file with nothing to correct is not rewritten at all: neither its content nor
its modification time changes. A file that does need correcting is first written
in full to "<name>.tmp", which is then renamed over the old file in one step, so
a file is never left half-written. (A ".tmp" file left behind by a killed run is
harmless: the next run overwrites it.) If the script stops part-way, for example
because a file is missing, the files already repaired stay repaired, and running
the script again finishes the job.

One run at a time
-----------------
Two copies of this script working on the same unrepaired database at the same
moment can write the same "<name>.tmp" file together and damage it. So start one
workflow, let it get past this step, and only then start another one that shares
the database. Once a database is repaired (and is still reached under the same
path, see below), this step only reads it, so later runs are harmless.

One path per database
---------------------
The new paths are built from the folder given with --database, after following
symbolic links. If the database is moved, or if people reach the same files
under different paths (for example through different mount points), the next run
rewrites the paths again to match ITS folder. That is correct, but two people
who use two different paths to one shared database would keep rewriting each
other's paths. Agree on one path.

Where it fits in the workflow
-----------------------------
Run by       rule normalize_funannotate_db (rules/shared/40_annotation.smk), once
             per run. That rule waits for rule funannotate_db, which downloads,
             checks and unpacks the snapshot when the folder does not hold one yet.
Takes in     --database      the unpacked funannotate database folder
                             (directories.funannotate_db in the config).
Writes       --summary-json  a small JSON file with the counts (how many
                             metadata files were examined, how many paths were
                             rewritten in them and in the database index) and the
                             absolute path of the database. It is the "database
                             is ready" marker (FUNANNOTATE_DB_MARKER,
                             00.resources/funannotate_db_paths.normalized.json):
                             rules funannotate_preprocess and
                             funannotate_prediction list it among their inputs,
                             so neither can start before this script has
                             finished successfully. Snakemake only needs the
                             file to exist; the numbers are for people.
             --log           the same counts as one line of text, for the rule's
                             log file (logs/normalize_funannotate_db.log).
Exit status  0 = finished.
             1 = something is wrong. A message starting with "ERROR:" goes to
             stderr; for an error the script does not expect, it is a Python
             traceback (Python's full error report) instead. The rule appends
             either one to the same log file. Neither the JSON nor the log line
             is written, so Snakemake sees the rule fail and starts nothing that
             needs the database.
             2 = a missing or unknown command-line option (argparse, Python's
             module for reading options, stops with a usage message).

Origin
------
The two functions that do the work (PART 1 below) were copied, with their logic
unchanged, from the Snakefile of FunFluxL+, the unpublished hybrid workflow,
where they were defined inline. New in FunFlux v2 are the command-line part
(PART 2: main() and the small functions it uses) and the comments.

Example
-------
    python normalize_funannotate_db_paths.py \
        --database /path/to/funannotate_db \
        --summary-json 00.resources/funannotate_db_paths.normalized.json \
        --log logs/normalize_funannotate_db.log

Only the Python standard library is used. The tests are in
test_normalize_funannotate_db_paths.py (same folder; run them with pytest).
"""

import argparse
import json
import os
import sys
from pathlib import Path

# CHANGED in FunFlux v2: "import json" and "from pathlib import Path" used to sit
# inside the two functions of PART 1 (in the Snakefile). They are now at the top
# of the file, together with the imports that the command-line part (PART 2)
# needs. Nothing else in the two functions was changed.


# =============================================================================
# PART 1: rewrite the stored paths (the two functions copied from the Snakefile)
# =============================================================================
# Apart from the two imports moved to the top of the file (see above), their
# code is unchanged. Kept word for word as well: the two docstrings and the two
# comments that sat inside _rebase_funannotate_db_path in FunFluxL+; each of those
# comments is followed by an "In plain words" version. Every other comment is new.
#
#   _rebase_funannotate_db_path(archived_path, database)
#       Works on ONE stored path. It answers the question "where is the file
#       that this old path was talking about, under my database folder?" and
#       returns the new path as text, or raises ValueError if the file is not
#       there.
#
#   normalize_funannotate_db_paths(database_root)
#       Works on the WHOLE database. It finds every stored path (in the species
#       info.json files and in the database index), asks the first function for
#       the right path, rewrites the files in which a path changed, and returns
#       the counts. main() in PART 2 calls it.


# ---- One stored path -> the right path on this computer ---------------------
# archived_path  one path exactly as written in a metadata file, usually an
#                absolute path on the computer that built the snapshot, e.g.
#                /home/builder/db/trained_species/Aspergillus_nidulans/augustus
#                It comes from an entry of a species' info.json, or from the
#                third column of funannotate-db-info.txt.
# database       the funannotate database folder on THIS computer, as a
#                pathlib.Path (normalize_funannotate_db_paths has already made
#                it absolute and followed symbolic links).
# returns        the matching path under `database`, as text. It is returned
#                only after checking that the path exists.
# raises         ValueError if the stored path is empty, or if the file or
#                folder it names cannot be found under the database.
# Called only by normalize_funannotate_db_paths, once for every stored path.
def _rebase_funannotate_db_path(archived_path, database):
    """Resolve an archived database path under the configured database root."""
    archived = Path(archived_path)

    # An absolute path starts with the root "/". Drop it, so that the remaining
    # pieces (home, builder, funannotate_db, trained_species, ...) can be
    # attached to a different folder. A relative path is used as it is.
    path_parts = archived.parts[1:] if archived.is_absolute() else archived.parts
    # Nothing left (the stored path was "/" or "."): there is no file to look for.
    if not path_parts:
        raise ValueError(f"Invalid empty Funannotate database path: {archived_path}")

    # Trained-species paths have a stable directory marker, so preserve the
    # complete suffix rather than relying on a potentially non-unique basename.
    # In plain words: every species has a folder called "augustus", so keeping
    # only that last name would be ambiguous. Instead everything from
    # "trained_species" onwards is kept, which includes the species name:
    # path_parts.index() finds the "trained_species" piece, and that piece and
    # all the pieces after it are attached to the database folder.
    if "trained_species" in path_parts:
        repaired_path = database.joinpath(
            *path_parts[path_parts.index("trained_species"):]
        )
        if repaired_path.exists():
            return str(repaired_path)
        # No second try with a shorter path: if this species folder is
        # missing, the snapshot is incomplete.
        raise ValueError(
            "Trained-species asset is missing after extraction: "
            f"{repaired_path}"
        )

    # Database-index entries normally point to files at the database root.
    # Find the longest trailing portion of the archived path that exists here,
    # without assuming the name of the machine or original database directory.
    # In plain words: for /home/builder/db/uniprot.dmnd the candidates are, in
    # this order, <database>/home/builder/db/uniprot.dmnd, then
    # <database>/builder/db/uniprot.dmnd, then <database>/db/uniprot.dmnd, and
    # last the bare file name <database>/uniprot.dmnd. The first candidate that
    # exists is used.
    for start in range(len(path_parts)):
        repaired_path = database.joinpath(*path_parts[start:])
        if repaired_path.exists():
            return str(repaired_path)

    # Not even the bare file name exists under the database folder.
    raise ValueError(
        "Database asset is missing after extraction for archived path: "
        f"{archived_path}"
    )


# ---- The whole database -------------------------------------------------------
# database_root  the funannotate database folder: the --database option of
#                main(), which the rule fills from directories.funannotate_db.
#                Rule funannotate_db has unpacked the snapshot into it, in this
#                run or an earlier one.
# does           A. goes through trained_species/*/info.json, then
#                B. goes through funannotate-db-info.txt.
#                In both it replaces every wrong path by the right one (from
#                _rebase_funannotate_db_path) and rewrites a file only if at
#                least one path in it changed.
# returns        a dictionary of counts, which main() writes to the summary
#                JSON and to the log line:
#                  metadata_files               info.json files examined
#                                               (changed or not)
#                  metadata_paths_updated       single paths rewritten in them
#                  database_info_paths_updated  paths rewritten in
#                                               funannotate-db-info.txt
# raises         ValueError if trained_species/ or funannotate-db-info.txt is
#                missing, or if a stored path cannot be matched to an existing
#                file (see _rebase_funannotate_db_path).
# Order of work  the info.json files one by one (alphabetical by species
#                folder), then the index. Every file is replaced on its own, so
#                an error part-way leaves the database partly repaired but never
#                corrupted. Fix the cause and run it again.
def normalize_funannotate_db_paths(database_root):
    """Rebase database-index and trained-species paths onto the configured root."""
    # Make the database folder an absolute path and follow symbolic links, so
    # that the new paths written into the files are the folder's real location.
    database = Path(database_root).resolve()

    # trained_species/ is part of every correctly unpacked snapshot. If it is
    # missing, the folder given is wrong or the unpacking failed: stop here.
    trained_species = database / "trained_species"
    if not trained_species.is_dir():
        raise ValueError(
            f"Funannotate trained-species directory is missing: {trained_species}"
        )

    # ---- Part A: the trained-species metadata (one info.json per species) ----
    metadata_files = 0   # info.json files examined, changed or not
    paths_updated = 0    # single paths rewritten inside those files
    # sorted(): species folders in alphabetical order, the same on every computer.
    for metadata_path in sorted(trained_species.glob("*/info.json")):
        metadata_files += 1
        with metadata_path.open() as handle:
            metadata = json.load(handle)

        # Remember whether this file needs rewriting at all. A file with
        # nothing to change is never touched.
        changed = False

        # info.json maps each gene finder (augustus, genemark, ...) to a list
        # of entries. A value that is not a list is not such an entry: skipped.
        for entries in metadata.values():
            if not isinstance(entries, list):
                continue
            for entry in entries:
                # Only dictionary entries with a non-empty "path" matter. The
                # empty placeholders ({}) of untrained gene finders are skipped.
                archived_path = entry.get("path") if isinstance(entry, dict) else None
                if not archived_path:
                    continue

                # Work out where this folder really is on this computer ...
                repaired_path = _rebase_funannotate_db_path(
                    archived_path, database
                )
                # ... and replace the stored path only if it differs. A path
                # that is already right is not counted; that is what makes a
                # second run report zero changes.
                if entry["path"] != repaired_path:
                    entry["path"] = repaired_path
                    paths_updated += 1
                    changed = True

        if changed:
            # Write the new content to "info.json.tmp" first, then rename it
            # over the old file in one step. An interrupted run therefore never
            # leaves a half-written info.json. The rewritten file has a 2-space
            # indent and a final newline, whatever its layout was before.
            temporary = metadata_path.with_name(metadata_path.name + ".tmp")
            with temporary.open("w") as handle:
                json.dump(metadata, handle, indent=2)
                handle.write("\n")
            temporary.replace(metadata_path)

    # ---- Part B: the database index (funannotate-db-info.txt) ----------------
    db_info_path = database / "funannotate-db-info.txt"
    if not db_info_path.is_file():
        raise ValueError(f"Funannotate database index is missing: {db_info_path}")

    # The whole index is read and checked before anything is written, so an
    # error in any line leaves the index file untouched.
    database_info_paths_updated = 0   # index lines whose path was rewritten
    normalized_lines = []             # the new file content, line by line
    with db_info_path.open() as handle:
        for line in handle:
            # One database per line, columns separated by tabs:
            # name, type, PATH, version, date, number of records, md5.
            fields = line.rstrip("\n").split("\t")

            # Lines with fewer than three columns, or an empty third column,
            # are not index entries and are kept exactly as they were.
            if len(fields) >= 3 and fields[2]:
                repaired_path = _rebase_funannotate_db_path(fields[2], database)
                if fields[2] != repaired_path:
                    fields[2] = repaired_path
                    database_info_paths_updated += 1
                # An entry line is put back together with a final newline, even
                # if it was the last line and had none. That matters only when
                # the file is rewritten below.
                line = "\t".join(fields) + "\n"
            normalized_lines.append(line)

    # Rewrite the index only if at least one path changed, again through a
    # temporary file that is renamed over the old one.
    if database_info_paths_updated:
        temporary = db_info_path.with_name(db_info_path.name + ".tmp")
        with temporary.open("w") as handle:
            handle.writelines(normalized_lines)
        temporary.replace(db_info_path)

    return {
        "metadata_files": metadata_files,
        "metadata_paths_updated": paths_updated,
        "database_info_paths_updated": database_info_paths_updated,
    }


# =============================================================================
# PART 2 (new in FunFlux v2): the command-line part
# =============================================================================
# In FunFluxL+ this work was done by the `run:` block (Python code written
# directly inside a rule) of its Snakemake rule normalize_funannotate_db. Now
# that rule runs this script instead, and the functions below do what the block
# did: run PART 1, then write the summary JSON and the one-line log, with the
# same layout and wording.


# Create the folder that will hold `file_path`, if it does not exist yet.
# Snakemake normally creates output folders itself, but this script may also be
# run by hand. A bare file name such as "marker.json" has no folder part, and
# os.makedirs("") would raise an error, so that case is skipped.
# Called by write_summary_json and write_log_line.
def make_parent_folder(file_path):
    """Make sure the folder for file_path exists."""
    parent_folder = os.path.dirname(file_path)
    if parent_folder:
        os.makedirs(parent_folder, exist_ok=True)


# Write the counts as a JSON file, through a temporary file.
# summary       the dictionary returned by normalize_funannotate_db_paths, plus
#               the "database" entry added by main().
# summary_path  where the "database is ready" marker goes
#               (00.resources/funannotate_db_paths.normalized.json in a run).
# The text goes to "<summary_path>.tmp" first, which is then renamed to the final
# name in one step, so Snakemake never sees a half-written marker. The layout
# (2-space indent, keys in the order they were added, final newline) is the same
# as the old rule produced. Called by main().
def write_summary_json(summary, summary_path):
    """Write the summary dictionary as JSON, via a temporary file."""
    make_parent_folder(summary_path)
    temporary = summary_path + ".tmp"
    with open(temporary, "w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    os.replace(temporary, summary_path)


# Write the one-line, human-readable result into the rule's log file
# (logs/normalize_funannotate_db.log in a run). The wording is the one
# FunFluxL+'s rule used. The file is overwritten, not added to. Called by main().
def write_log_line(summary, log_path):
    """Write the one-line log message."""
    make_parent_folder(log_path)
    with open(log_path, "w") as handle:
        handle.write(
            "[normalize_funannotate_db] normalized "
            f"{summary['metadata_paths_updated']} trained-species path(s) "
            f"across {summary['metadata_files']} metadata file(s), and "
            f"{summary['database_info_paths_updated']} database-index path(s).\n"
        )


# Define the three command-line options. All are required: the Snakemake rule
# passes them in. A missing one makes argparse stop with exit status 2.
# Called by main().
def build_parser():
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite the absolute paths stored in a downloaded funannotate "
            "database so that they point into the database folder on this "
            "computer. WARNING: this edits files inside the database folder, "
            "in place."
        )
    )
    parser.add_argument(
        "--database",
        required=True,
        metavar="DIR",
        help="the unpacked funannotate database folder (files in it are edited)",
    )
    parser.add_argument(
        "--summary-json",
        required=True,
        metavar="FILE",
        help="JSON file to write with the counts (the 'database is ready' marker)",
    )
    parser.add_argument(
        "--log",
        required=True,
        metavar="FILE",
        help="text file to write with a one-line summary",
    )
    return parser


# Entry point, used by the Snakemake rule normalize_funannotate_db and by the
# tests. argv is the list of command-line words (None means "read them from the
# real command line"). Returns the exit status: 0 = done, 1 = something is wrong.
def main(argv=None):
    """Repair the database, then write the summary JSON and the log line."""
    args = build_parser().parse_args(argv)

    try:
        # The only step that edits the database folder.
        summary = normalize_funannotate_db_paths(args.database)

        # Record which database this summary is about. abspath makes the path
        # absolute (relative to the current folder) but, unlike the resolve()
        # inside the function, does not follow symbolic links.
        summary["database"] = os.path.abspath(args.database)

        # Only after the database was repaired successfully are the marker and
        # the log written. If anything above fails they are not written (and
        # Snakemake removes the old output of a failed rule), so the rules that
        # wait for the marker do not start.
        write_summary_json(summary, args.summary_json)
        write_log_line(summary, args.log)
    except (ValueError, OSError) as error:
        # ValueError is what the copied functions raise for a problem with the
        # database: a missing folder, a missing file, a path that matches no
        # file. OSError covers files that cannot be read or written, for example
        # a database on a read-only disk. Both are shown as one clear line
        # instead of a traceback (Python's full error report). A damaged
        # info.json also arrives as a ValueError, because that is the JSON
        # reader's error type, but then the message does not name the file: the
        # copied function was left as it was. Any other error is not caught:
        # Python prints the traceback and also exits with status 1.
        sys.stderr.write(f"ERROR: {error}\n")
        return 1

    return 0


# Run main() only when the file is started as a program, not when the tests
# import it; its return value becomes the exit status.
if __name__ == "__main__":
    sys.exit(main())
