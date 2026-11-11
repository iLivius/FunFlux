#!/usr/bin/env python3
"""Check that FunFlux reads its configuration the way it should.

What this is for
    Almost all of FunFlux's input checking happens at "parse time", that is while
    Snakemake READS the workflow, before any job starts (hence this script's name):
    the mode, the input directories, the sample names, the annotation table and every
    tuning parameter are checked in workflow/rules/shared/00_common.smk. That makes a
    dry run (`snakemake -n -p`: plan the jobs, print the command each one would run,
    run none) the test. This script runs dry runs on purpose-built inputs and compares
    what Snakemake prints with what it should print. The cases come in seven groups:

      * every example block in config/config_test_examples.yaml must dry-run
        cleanly and print the mode, sample and decision lines a user expects to see;
      * a list of deliberately broken configurations must each be REFUSED, with the
        [FunFlux] message that tells the user what to fix;
      * in contigs mode, each optional feature (the contamination screen, SignalP 6,
        GeneMark-ES, Tiberius, local antiSMASH and UNITE copies, the UNITE reference filter,
        the AMR screen) must add or remove the expected rules in the planned jobs, and the
        number of jobs must be the expected one;
      * in the three read-based modes (illumina, longread, hybrid), the example blocks
        and a few variants of them (the AMR screen on among them) must plan the expected
        number of jobs, the two Medaka rules must be planned only for the ONT samples that
        use Medaka, and the commands Snakemake would run must carry each sample's settings
        (the Medaka model, the Filtlong guidance, the Flye read mode, the SPAdes memory
        limit, the assembler's file the evidence table is made from);
      * rule download_all, named on the command line, must plan exactly the download
        steps that each configuration defines, and a normal run must not plan it;
      * the heavy and the long steps must ask for the memory and the time limit that
        resources.ram_gb and resources.runtime_hours set (what a cluster scheduler is
        given), and the light steps for nothing;
      * a few unusual settings must be ACCEPTED and read the way the user meant them
        (some print a NOTE).

    It needs no real data and no databases: it builds stand-in files (empty, or a line
    or two) in a scratch folder, which is deleted at the end unless --keep is given.
    The dry runs are made without --sdm conda, so no conda environment is looked at or
    built and no conda program is needed: Snakemake alone is enough. Run it after every
    change to 00_common.smk, to a rule file or to the config files. It is a development
    tool: nothing in the workflow calls it.

How to run (from the repository root; use the Snakemake environment's python,
because PyYAML ships with Snakemake):
    /path/to/snakemake-env/bin/python miscellaneous/check_parse_time.py
Options: --snakefile (the Snakefile to test, relative to the repository root; default
workflow/Snakefile), --snakemake (the snakemake program to call; default: the one
on the PATH), --keep (do not delete the scratch folder, to look at the configurations
it wrote).

Output: a heading per group, then one PASS or FAIL line per case. A FAIL lists what
did not match and the last lines Snakemake printed. A dry run still running after
DRY_RUN_TIMEOUT seconds (300) is stopped and counted as a FAIL of its case. After the
last case comes one summary line; with --keep, one more line then names the kept
scratch folder.

Exit status: 0 if every case behaved as expected, 1 if any did not.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import yaml

# The repository root (the folder above miscellaneous/). Every dry run is started there.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The example configurations: one complete config per block (cut apart by extract_blocks).
EXAMPLES_FILE = os.path.join(REPO_ROOT, "config", "config_test_examples.yaml")

# The example blocks name stub (stand-in) files under /tmp/ff (see the header of the
# examples file). This script makes these stubs, and more, in its own scratch folder
# instead (build_stubs), and rewrites the paths on the way (extract_blocks).
EXAMPLE_ROOT = "/tmp/ff"

# Marker used in the edit lists below: "remove this key from the config". It is a
# unique object, not "" or False, because some cases SET a key to exactly those values.
DELETE = object()

# The longest a single dry run may take, in seconds. A dry run normally takes a few
# seconds; one that is still running after this long is stopped and counted as a FAIL
# of its case (see run_snakemake).
DRY_RUN_TIMEOUT = 300


def build_stubs(base):
    # Create the stand-in files that the example blocks and most cases share, inside the
    # scratch folder `base`. Some cases add their own input folder or annotation table
    # later (variant() and make_params()), and a few point at paths that are never
    # created, on purpose (nowhere, ref/none.tsv, ref/missing.faa, ref/no_overrides.tsv).
    # Called once by main(), before the example blocks are read; returns nothing.
    #
    # The input folders hold the same empty sample files that the header of
    # config_test_examples.yaml tells a user to create by hand. Sample discovery only
    # looks at file NAMES, so empty files are enough.
    layout = {
        "in/illumina": ["SAM1_R1.fastq.gz", "SAM1_R2.fastq.gz"],
        "in/longread": ["SAM1_ont.fastq.gz", "SAM2_pacbio.fastq.gz"],
        "in/hybrid_short": ["SAM1_R1.fastq.gz", "SAM1_R2.fastq.gz", "SAM2_R1.fastq.gz", "SAM2_R2.fastq.gz"],
        "in/contigs": ["SAM1.fasta"],
    }
    for directory, names in layout.items():
        os.makedirs(os.path.join(base, directory), exist_ok=True)
        for name in names:
            open(os.path.join(base, directory, name), "w").close()
    os.makedirs(os.path.join(base, "ref"), exist_ok=True)
    # Stand-ins for the optional features, so that a dry run can plan their jobs:
    #   * a SignalP 6 archive: 00_common.smk and the dry run only check that it exists;
    #   * a GeneMark folder with an executable gmes_petap.pl (00_common.smk checks that it
    #     exists and may be run), plus a fake home folder holding the licence key file
    #     .gm_key (see run_snakemake);
    #   * local copies of the antiSMASH and UNITE databases, each holding the one folder
    #     (clusterblast/) or file (unite_its_sintax.fasta) that 00_common.smk looks for to
    #     accept it as the real thing.
    open(os.path.join(base, "ref", "signalp6.tar.gz"), "w").close()
    os.makedirs(os.path.join(base, "ref", "genemark"), exist_ok=True)
    gmes_path = os.path.join(base, "ref", "genemark", "gmes_petap.pl")
    with open(gmes_path, "w") as handle:
        handle.write("#!/bin/sh\n")
    os.chmod(gmes_path, 0o755)
    os.makedirs(os.path.join(base, "home"), exist_ok=True)
    open(os.path.join(base, "home", ".gm_key"), "w").close()
    os.makedirs(os.path.join(base, "ref", "antismash_db", "clusterblast"), exist_ok=True)
    os.makedirs(os.path.join(base, "ref", "unite_db"), exist_ok=True)
    open(os.path.join(base, "ref", "unite_db", "unite_its_sintax.fasta"), "w").close()
    # The BLAST database stand-in. The contamination screen's BlobTools rule (blob_json)
    # lists nodes.dmp and names.dmp (the NCBI taxonomy dump that lives inside the BLAST
    # database folder) as INPUTS, on purpose, so that a missing dump is reported before any
    # job starts. A dry run that plans the screen therefore needs the two files to exist;
    # empty stand-ins are enough for planning the jobs.
    os.makedirs(os.path.join(base, "ref", "blastdb"), exist_ok=True)
    for dump_name in ("nodes.dmp", "names.dmp"):
        open(os.path.join(base, "ref", "blastdb", dump_name), "w").close()
    # Protein evidence for the annotation tables: one tiny protein (00_common.smk refuses an
    # empty protein file), and an empty file for the case that tests exactly that.
    with open(os.path.join(base, "ref", "proteins.faa"), "w") as handle:
        handle.write(">p1\nMKVLA\n")
    with open(os.path.join(base, "ref", "empty.faa"), "w"):
        pass
    # Annotation table rows with 2 and with 5 columns instead of 3 or 4, for the "wrong
    # number of columns" cases.
    with open(os.path.join(base, "ref", "badcols.tsv"), "w") as handle:
        handle.write("SAM1\tFusarium graminearum\n")
    with open(os.path.join(base, "ref", "fivecols.tsv"), "w") as handle:
        handle.write(f"SAM1\tFusarium graminearum\t{base}/ref/proteins.faa\tfusarium_graminearum\textra\n")
    # The annotation table the example blocks use (files.annotation_params), one row per
    # sample: sample, species, protein evidence, AUGUSTUS model. It needs a row for both
    # SAM1 and SAM2: the longread and hybrid examples have both samples, and 00_common.smk
    # refuses a sample that has no row.
    with open(os.path.join(base, "ref", "params.tsv"), "w") as handle:
        for sample in ("SAM1", "SAM2"):
            handle.write(f"{sample}\tFusarium graminearum\t{base}/ref/proteins.faa\tfusarium_graminearum\n")
    # Decontamination tables (parameters.decontamination) for the cases that test how
    # 00_common.smk checks them. Tab-separated, like the real ones:
    #   include_by_sample.tsv          : the per-sample include table; gives SAM1 its genus.
    #   include_without_header.tsv     : the same row without the header line (refused).
    #   overrides_auto.tsv             : an override row with the removed mode auto (refused).
    #   overrides_include_no_genus.tsv : an override row that sets include for SAM1 and
    #                                    leaves its include_genera cell empty (refused).
    #   overrides_SAM2_genus.tsv       : an override row that sets include for SAM2 and names
    #                                    two genera in one cell.
    #   overrides_bad_discard.tsv      : an override row for SAM1 whose discard_no_hit cell
    #                                    holds a word that is not yes or no (refused).
    #   include_latin1.tsv             : a per-sample include table saved in Latin-1, an old
    #                                    Windows encoding, instead of UTF-8 (refused).
    tables = {
        "include_by_sample.tsv": "sample\tgenus\nSAM1\tFusarium\n",
        "include_without_header.tsv": "SAM1\tFusarium\n",
        "overrides_auto.tsv": "sample\tmode\nSAM1\tauto\n",
        "overrides_include_no_genus.tsv": "sample\tmode\tinclude_genera\nSAM1\tinclude\t\n",
        "overrides_SAM2_genus.tsv": "sample\tmode\tinclude_genera\nSAM2\tinclude\tFusarium;Gibberella\n",
        "overrides_bad_discard.tsv": "sample\tmode\tdiscard_no_hit\nSAM1\texclude\tmaybe\n",
    }
    for name, text in tables.items():
        with open(os.path.join(base, "ref", name), "w") as handle:
            handle.write(text)
    # The genus "Fusárium" written in Latin-1, where the accented letter is one byte
    # (0xE1) that is not valid UTF-8.
    with open(os.path.join(base, "ref", "include_latin1.tsv"), "w", encoding="latin-1") as handle:
        handle.write("sample\tgenus\nSAM1\tFusárium\n")
    # A funannotate database stand-in that holds only its index file,
    # funannotate-db-info.txt: 00_common.smk then treats the database as present and
    # checks that the BUSCO set named in parameters.prediction.busco_db is a folder in it.
    os.makedirs(os.path.join(base, "ref", "fundb"), exist_ok=True)
    open(os.path.join(base, "ref", "fundb", "funannotate-db-info.txt"), "w").close()
    # A second stand-in for the AUGUSTUS model cases: the index file, the dikarya BUSCO set
    # (so the BUSCO check passes) and three trained species, as empty folders under
    # trained_species/. 00_common.smk checks a named model against these folder names.
    fundb_species = os.path.join(base, "ref", "fundb_species")
    os.makedirs(os.path.join(fundb_species, "dikarya"), exist_ok=True)
    open(os.path.join(fundb_species, "funannotate-db-info.txt"), "w").close()
    for species in ("fusarium_graminearum", "neurospora_crassa", "aspergillus_nidulans"):
        os.makedirs(os.path.join(fundb_species, "trained_species", species), exist_ok=True)
    # A third stand-in that holds a model stored under the name funannotate gives SAM1's own
    # training ("Fusarium graminearum" + "_" + "SAM1", in lower case): funannotate would use
    # it instead of training when column 4 is empty, and FunFlux says so in a NOTE.
    fundb_own = os.path.join(base, "ref", "fundb_own")
    os.makedirs(os.path.join(fundb_own, "dikarya"), exist_ok=True)
    open(os.path.join(fundb_own, "funannotate-db-info.txt"), "w").close()
    os.makedirs(os.path.join(fundb_own, "trained_species", "fusarium_graminearum_sam1"), exist_ok=True)
    # Annotation tables for two ways of writing rows: SAM1's fourth column left empty with a
    # trailing tab (the "left empty" form of the template), and a mixed batch in which SAM1
    # names a model and SAM2's fourth column is empty.
    with open(os.path.join(base, "ref", "params_trailing_tab.tsv"), "w") as handle:
        handle.write(f"SAM1\tFusarium graminearum\t{base}/ref/proteins.faa\t\n")
    with open(os.path.join(base, "ref", "params_mixed.tsv"), "w") as handle:
        handle.write(f"SAM1\tFusarium graminearum\t{base}/ref/proteins.faa\tfusarium_graminearum\n")
        handle.write(f"SAM2\tFusarium graminearum\t{base}/ref/proteins.faa\t\n")


def extract_blocks(base):
    # Cut config_test_examples.yaml into its standalone example configurations.
    # Takes `base`, the scratch folder build_stubs filled. Called once by main().
    # Returns {block name: config as a Python dict}, for example
    # {"mode = illumina": {...}, "mode = contigs": {...}, ...}.
    #
    # Read the file as text and point its paths at the scratch folder: first the stub
    # inputs (/tmp/ff/...), then the placeholder BLAST database path, which becomes the
    # stand-in folder that holds nodes.dmp and names.dmp (see build_stubs).
    with open(EXAMPLES_FILE) as handle:
        text = handle.read()
    text = text.replace(EXAMPLE_ROOT, base)
    text = text.replace("/abs/placeholder/db/NCBI_core_nt", base + "/ref/blastdb")

    # A block runs from a "BEGIN CONFIG: <title>" line to the next "END CONFIG" line (both
    # framed with ═ characters). The regular expression (a text pattern) captures the title
    # and the YAML text between the two lines; re.S lets "." match line breaks as well, so
    # one match can span the whole block.
    blocks = {}
    for match in re.finditer(r"# ═+ BEGIN CONFIG: (.+?) ═+\n(.*?)# ═+ END CONFIG", text, flags=re.S):
        title = match.group(1)
        block_text = match.group(2)
        # The block name is the title without its explanation in brackets:
        # "mode = longread (SAM1 is ONT, ...)" becomes "mode = longread".
        name = title.split(" (")[0].strip()
        blocks[name] = yaml.safe_load(block_text)
    return blocks


def set_value(config, dotted_key, value):
    # Apply one edit to a config dict: set the key named by a dotted path such as
    # "parameters.busco.lineage" (that is, config["parameters"]["busco"]["lineage"]),
    # creating any missing level on the way. The value DELETE removes the key instead
    # (a key that is not there is left alone). Changes `config` in place; returns nothing.
    # Called by variant() for each edit of a test case.
    parts = dotted_key.split(".")
    node = config
    # Walk down to the level that holds the last key. setdefault returns the existing
    # level, or adds an empty one when it is missing.
    for part in parts[:-1]:
        node = node.setdefault(part, {})
    if value is DELETE:
        node.pop(parts[-1], None)
    else:
        node[parts[-1]] = value


def run_snakemake(snakemake, snakefile, config_path, home=None, target=None):
    # Run one dry run (snakemake -n -p: plan the jobs, print their commands, run none)
    # from the repository root.
    #   snakemake   : the snakemake program to call (the --snakemake option)
    #   snakefile   : the Snakefile to test, relative to the repository root (--snakefile)
    #   config_path : the YAML file check() wrote, or None to run without any --configfile
    #   home        : a folder to use as HOME for this one run, or None to keep the real
    #                 HOME. GeneMark looks for its licence key in ~/.gm_key, and a test must
    #                 not depend on the real key of whoever runs it.
    #   target      : a rule to name on the command line (for example download_all), or
    #                 None to plan the default target, rule all, as a normal run does.
    # Returns (exit status, everything it printed: standard output, then standard error).
    # A dry run still running after DRY_RUN_TIMEOUT seconds is stopped (subprocess.run
    # then raises TimeoutExpired); it returns the exit status None and whatever the run
    # had printed by then, which may show where it hung. check() reports that as a FAIL
    # of the case, and the other cases still run. Called by check() in main().
    # -p prints the shell command of every planned job with its placeholders filled in,
    # so that a case can check what a rule hands a tool (the Medaka model name, the
    # Filtlong guidance options, the SPAdes memory limit). No --sdm conda: the
    # environments play no part in how a configuration is read, and the script must not
    # need the conda program.
    # The target comes right after the program name: an option that takes several values,
    # such as --configfile, would read a name written after it as one more of its values.
    command = [snakemake]
    if target:
        command.append(target)
    command += ["-n", "-p", "--cores", "1", "--snakefile", snakefile]
    if config_path:
        command += ["--configfile", config_path]
    if home:
        # A copy of this script's own environment variables, with HOME replaced.
        environment = dict(os.environ)
        environment["HOME"] = home
    else:
        # None: the dry run inherits this script's environment unchanged.
        environment = None
    try:
        result = subprocess.run(
            command, cwd=REPO_ROOT, capture_output=True, text=True, timeout=DRY_RUN_TIMEOUT, env=environment
        )
    except subprocess.TimeoutExpired as stopped:
        # What the run printed before it was stopped. TimeoutExpired holds each part as
        # raw bytes, even though text=True was asked for, or None when nothing came, so
        # each part is turned into text here.
        printed_before_stop = ""
        for part in (stopped.stdout, stopped.stderr):
            if part is None:
                continue
            if isinstance(part, bytes):
                part = part.decode("utf-8", errors="replace")
            printed_before_stop += part
        return None, printed_before_stop
    return result.returncode, result.stdout + result.stderr


def job_counts(output):
    # Read the "Job stats" table of a dry run into {rule name: number of jobs}, with the
    # key "total" for the whole run. Takes everything one dry run printed (from
    # run_snakemake); returns an empty dict when there is no table (a refused
    # configuration stops before any job is planned). Called by check() for the cases
    # that pass expect_jobs.
    #
    # The table looks like this: the rule name, spaces, the number of its jobs, and a
    # last line with the total. Snakemake prints it twice in a dry run (before the job
    # list and again after it, with the same numbers); the first one is read.
    #   Job stats:
    #   job                   count
    #   ------------------  -------
    #   all                       1
    #   filter_contigs            1
    #   ...
    #   total                    24
    counts = {}
    # From "Job stats:" to the first "total <number>" line after it. re.S lets "." match
    # line breaks, so the match spans the whole table; "?" after ".*" takes the shortest
    # match, so it stops at the first total line.
    table = re.search(r"Job stats:\n(.*?\ntotal\s+\d+)", output, flags=re.S)
    if table is None:
        return counts
    for line in table.group(1).splitlines():
        # A table row: a name, spaces, digits, nothing else. The header line ("job
        # count") and the dashes under it do not match.
        row = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\s+(\d+)\s*$", line)
        if row:
            counts[row.group(1)] = int(row.group(2))
    return counts


def job_resources(output):
    # Read what every planned job asks a cluster for into
    # {rule name: [{resource name: value as text}, one dict per job of that rule]}.
    # Takes everything one dry run printed (from run_snakemake); a rule with no job in the
    # plan is not in the result. Called by check() for the cases that pass
    # expect_resources.
    #
    # A dry run prints each job as a block that starts with "rule <name>:" (or
    # "localrule <name>:" for a job that runs on the launch machine) and holds one
    # resources line. Snakemake 9.14.6 writes the SPAdes job of the illumina example as
    #   rule illumina_assembly:
    #       ...
    #       resources: tmpdir=<TBD>, mem_gb=32, mem_mb=32768, mem_mib=31250, runtime=2880
    # and Snakemake 9.27.1 writes the same job as
    #       resources: tmpdir=/tmp, mem_gb=32, runtime=2880, mem_mb=32768, mem=32.77 GB, mem_mib=31250
    # that is, with the real temporary folder instead of <TBD>, in another order and with
    # one more entry (mem). Reading the line into names and values lets a case compare
    # the values FunFlux sets (mem_gb, mem_mb, java_mem, runtime) in any Snakemake version.
    jobs = {}
    current_rule = None
    for line in output.splitlines():
        stripped = line.strip()
        # The first line of a job block names its rule.
        start = re.match(r"^(?:local)?rule (\w+):$", stripped)
        if start:
            current_rule = start.group(1)
            continue
        # The resources line of that job: entries separated by ", ", each one name=value.
        if current_rule and stripped.startswith("resources: "):
            resources = {}
            for entry in stripped[len("resources: "):].split(", "):
                name, _, value = entry.partition("=")
                resources[name] = value
            jobs.setdefault(current_rule, []).append(resources)
            # This job is read; the next resources line belongs to the next job block.
            current_rule = None
    return jobs


def main():
    # Run every test case and report. Makes the scratch folder and fills it (build_stubs),
    # reads the example blocks (extract_blocks), then runs the cases below one by one;
    # check() prints PASS or FAIL for each. Ends with the summary line, and deletes the
    # scratch folder unless --keep was given. Returns the exit status, 0 when every case
    # behaved as expected and 1 otherwise, which the last lines of this file hand to the
    # shell.
    #
    # Command-line options (described in the module docstring at the top).
    parser = argparse.ArgumentParser(description="Dry-run FunFlux's example and deliberately broken configs.")
    parser.add_argument("--snakefile", default="workflow/Snakefile")
    parser.add_argument("--snakemake", default="snakemake")
    parser.add_argument("--keep", action="store_true", help="keep the scratch directory")
    args = parser.parse_args()

    # One scratch folder for the whole run. It holds the stand-in files, the configuration
    # each case writes, and the extra folders and tables some cases make. mkdtemp creates it
    # in the system's temporary folder, with a unique name that starts with funflux_parse_.
    base = tempfile.mkdtemp(prefix="funflux_parse_")
    build_stubs(base)
    blocks = extract_blocks(base)
    # The labels of the cases that did not behave as expected, for the summary at the end.
    failures = []

    # The three helpers below are defined inside main() so that they can use base, blocks,
    # failures and args directly.

    def check(label, config, expect_ok, expect_texts, forbid_texts=(), home=None, expect_jobs=None, target=None,
              expect_resources=None):
        # Run one test case and print PASS or FAIL for it.
        #   label        : the case name printed on the PASS/FAIL line
        #   config       : the configuration as a dict (an example block, or a variant() of
        #                  one), or None to run without any --configfile
        #   expect_ok    : True  -> the dry run must succeed (exit status 0);
        #                  False -> the dry run must be refused (any other exit status)
        #   expect_texts : texts that must all appear in what Snakemake prints
        #   forbid_texts : texts that must NOT appear
        #   home         : handed to run_snakemake (only the GeneMark cases use it)
        #   expect_jobs  : optional {rule name: number of jobs}, plus the key "total" for
        #                  the whole run, compared with the "Job stats" table (job_counts).
        #                  A rule that is not in the table counts as 0 jobs, so
        #                  {"long_read_consensus": 0} says "this rule must not be planned".
        #   target       : handed to run_snakemake (only the download_all cases use it)
        #   expect_resources : optional {rule name: {resource name: value}}, compared with
        #                  the resources line of every job of that rule (job_resources).
        #                  The value None says "this job must not set this resource" (a
        #                  light step asks the scheduler for nothing). Only the named
        #                  resources are compared, so the entries Snakemake adds by itself
        #                  (tmpdir, and mem and mem_mib, worked out from mem_mb) never
        #                  matter.
        # A failed case is added to `failures`.
        #
        # The texts are looked for in a copy of the output in which every run of spaces,
        # tabs and line breaks is one space. Snakemake prints a multi-line shell command
        # as one line with runs of spaces where the line continuations were, so a text
        # such as "-2 <file> --min_length 1000" would otherwise depend on the indentation
        # of the rule file. Every text a case names is written with single spaces.
        #
        # Write the configuration to a YAML file in the scratch folder, named after the case
        # label (every run of characters other than letters and digits becomes "_").
        if config is None:
            config_path = None
        else:
            config_path = os.path.join(base, re.sub(r"[^A-Za-z0-9]+", "_", label) + ".yaml")
            with open(config_path, "w") as handle:
                yaml.safe_dump(config, handle)
        status, output = run_snakemake(args.snakemake, args.snakefile, config_path, home, target)

        # Compare the outcome with the expectation, and collect every mismatch. A dry run
        # that run_snakemake had to stop (status None) has no exit status and no complete
        # output to compare: it is a FAIL, and the lines shown under it are what the run
        # had printed before it was stopped.
        problems = []
        if status is None:
            problems.append(f"the dry run did not finish within {DRY_RUN_TIMEOUT} seconds and was stopped")
        else:
            succeeded = (status == 0)
            if succeeded != expect_ok:
                if expect_ok:
                    expected_status = "0"
                else:
                    expected_status = "non-zero"
                problems.append(f"exit status {status}, expected {expected_status}")
            flat_output = re.sub(r"\s+", " ", output)
            for text in expect_texts:
                if text not in flat_output:
                    problems.append(f"missing text: {text!r}")
            for text in forbid_texts:
                if text in flat_output:
                    problems.append(f"unexpected text: {text!r}")
            if expect_jobs:
                counts = job_counts(output)
                if not counts:
                    problems.append("no Job stats table in the output")
                else:
                    for name, expected in expect_jobs.items():
                        planned = counts.get(name, 0)
                        if planned != expected:
                            problems.append(f"{name}: {planned} job(s) planned, expected {expected}")
            if expect_resources:
                planned_resources = job_resources(output)
                for rule_name, expected_resources in expect_resources.items():
                    jobs_of_rule = planned_resources.get(rule_name, [])
                    if not jobs_of_rule:
                        problems.append(f"{rule_name}: no job planned, so its resources cannot be checked")
                        continue
                    for resources in jobs_of_rule:
                        for name, expected in expected_resources.items():
                            # The dry run prints text, so a number in the case is compared
                            # as text; get() gives None when the job does not set it.
                            if expected is not None:
                                expected = str(expected)
                            found = resources.get(name)
                            if found != expected:
                                found_text = "not set" if found is None else found
                                expected_text = "not set" if expected is None else expected
                                problems.append(f"{rule_name}: {name} is {found_text}, expected {expected_text}")

        # Print the verdict. A FAIL also shows the last 8 lines Snakemake printed, which
        # usually hold the [FunFlux] message or the error.
        if problems:
            failures.append(label)
            print(f"FAIL  {label}")
            for problem in problems:
                print(f"        {problem}")
            print("        --- last lines of output ---")
            for line in output.strip().splitlines()[-8:]:
                print(f"        {line}")
        else:
            print(f"PASS  {label}")

    def make_params(rows, model="fusarium_graminearum"):
        # Write a small annotation table (files.annotation_params) into the scratch folder.
        #   rows  : [(sample, protein file), ...]; each row gets the same species as the
        #           examples' table
        #   model : the AUGUSTUS model of every row (column 4), the examples' one by default;
        #           None writes rows of three columns, with no model at all
        # Returns the path of the new table; the case then sets files.annotation_params to it.
        # mkstemp creates a new, empty file with a unique name (params_<random>.tsv) and
        # returns an open file number for it and its path. The file number is closed at
        # once, so that the file is not left open; the file is then written through its
        # path.
        file_number, path = tempfile.mkstemp(dir=base, prefix="params_", suffix=".tsv")
        os.close(file_number)
        with open(path, "w") as handle:
            for sample, protein_file in rows:
                if model is None:
                    handle.write(f"{sample}\tFusarium graminearum\t{protein_file}\n")
                else:
                    handle.write(f"{sample}\tFusarium graminearum\t{protein_file}\t{model}\n")
        return path

    def variant(block_name, edits, files=None):
        # Make a test configuration from an example block with some keys changed.
        #   block_name : the example block to start from, e.g. "mode = illumina"
        #   edits      : [(dotted key, new value or DELETE), ...], applied by set_value()
        #   files      : optional {input key: [file names]}. For each entry a new folder of
        #                empty files with exactly these names is made in the scratch folder,
        #                and config["input"][input key] points at it (to test sample discovery).
        # Returns the new config dict.
        #
        # Copy the block by writing it out as YAML text and reading it back, so that the
        # edits below change only this copy and never the shared example block.
        block_as_text = yaml.safe_dump(blocks[block_name])
        config = yaml.safe_load(block_as_text)
        for dotted_key, value in edits:
            set_value(config, dotted_key, value)
        if files:
            for input_key, names in files.items():
                directory = tempfile.mkdtemp(dir=base, prefix="in_")
                for name in names:
                    open(os.path.join(directory, name), "w").close()
                config["input"][input_key] = directory
        return config

    # ── The test cases ──
    # Each check(...) call below is one case. Its arguments, in order: the label, the
    # configuration, whether the dry run must succeed (True) or be refused (False), the
    # texts that must appear, and (optional) the texts that must not appear. Most texts are
    # pieces of the messages 00_common.smk prints, so a message reworded there makes its
    # case fail here: change both together.
    print("=== the example blocks must dry-run cleanly ===")
    # The start of the line Snakemake ends a successful dry run with.
    done = "This was a dry-run"
    # The illumina example has no long reads, so the PacBio NOTE must not appear.
    check("illumina example", blocks["mode = illumina"], True,
          ["Mode: illumina", "Sample SAM1 will be processed.", done], ["NOTE: samples with a _pacbio"])
    # The longread and hybrid examples mix an ONT sample (SAM1) and a PacBio HiFi sample
    # (SAM2). Each sample must get its own decision line: assembler flag, minimap2 preset,
    # Medaka or not and, in hybrid mode, whether Filtlong is guided by the Illumina reads.
    # The longread case also checks the NOTE about PacBio HiFi (CCS) reads.
    check("longread example (ONT + PacBio mixed)", blocks["mode = longread"], True,
          ["Mode: longread", "Sample SAM1 will be processed.", "Sample SAM2 will be processed.",
           "[FunFlux] SAM1: ont reads -> adapter check + fastplong, hifiasm --ont, minimap2 map-ont, Medaka.",
           "[FunFlux] SAM2: pacbio reads -> hifiasm, minimap2 map-hifi, no Medaka.",
           "HiFi (CCS) reads", done])
    check("hybrid example (ONT + PacBio mixed)", blocks["mode = hybrid"], True,
          ["Mode: hybrid",
           "[FunFlux] SAM1: ont reads -> adapter check + fastplong, Filtlong guided by the Illumina reads, hifiasm --ont, minimap2 map-ont, Medaka.",
           "[FunFlux] SAM2: pacbio reads -> Filtlong not guided, hifiasm, minimap2 map-hifi, no Medaka.",
           done])
    # Contigs mode without the screen prints a NOTE that says what rule filter_contigs does
    # to the FASTA and warns that a mitochondrial genome in it would be annotated as
    # nuclear (the old wording, "annotated as it is", was not true). With the screen it
    # prints the screen settings instead, and no such NOTE.
    check("contigs example (annotation only)", blocks["mode = contigs"], True,
          ["Mode: contigs", "contigs mode without parameters.contigs.screen",
           "rule filter_contigs drops contigs under 500 bp or under 2x coverage", done],
          ["annotated as it is"])
    check("contigs example with the screen", blocks["mode = contigs with the contamination and mitochondrion screen"], True,
          ["Mode: contigs", "Mitochondrion screen: separate_confirmed", done],
          ["contigs mode without parameters.contigs.screen"])

    print("=== deliberately broken configurations must be refused ===")
    # The configuration as a whole, and its mode. The first case runs with no --configfile.
    check("no configuration at all", None, False, ["No configuration supplied"])
    check("v1-style config without a mode", variant("mode = illumina", [("mode", DELETE)]), False,
          ["this configuration has no 'mode' key", "Coming from v1"])
    check("mode typo", variant("mode = illumina", [("mode", "ilumina")]), False,
          ["mode: must be one of illumina | longread | hybrid | contigs", "'ilumina'"])
    # Input folders and sample discovery. The cases that pass a {input key: [file names]}
    # dict point that input key at a new folder holding only those (empty) files.
    check("input directory missing", variant("mode = illumina", [("input.illumina_dir", DELETE)]), False,
          ["requires input.illumina_dir"])
    # "nowhere" is a path that is never created.
    check("input directory is not a directory", variant("mode = illumina", [("input.illumina_dir", base + "/nowhere")]), False,
          ["is not an existing directory"])
    check("R2 mate missing", variant("mode = illumina", [], {"illumina_dir": ["SAM1_R1.fastq.gz"]}), False,
          ["has an R1 file but no R2 mate"])
    # SAM1 ends in .fastq.gz and SAM2 in .fq: all files of one input folder must share one ending.
    check("two file extensions in one folder", variant("mode = illumina", [],
          {"illumina_dir": ["SAM1_R1.fastq.gz", "SAM1_R2.fastq.gz", "SAM2_R1.fq", "SAM2_R2.fq"]}), False,
          ["more than one file extension in this folder"])
    check("unrecognised extension", variant("mode = contigs", [], {"contigs_dir": ["SAM1.txt"]}), False,
          ["is not recognised"])
    check("sample under both _ont and _pacbio", variant("mode = longread", [],
          {"longread_dir": ["SAM1_ont.fastq.gz", "SAM1_pacbio.fastq.gz"]}), False,
          ["has both _ont and _pacbio"])
    # _hifi is not a recognised ending, so the folder holds no long-read file at all, and the
    # message must name the ending to use (_pacbio).
    check("longread directory without ONT or PacBio files", variant("mode = longread", [],
          {"longread_dir": ["SAM1_hifi.fastq.gz"]}), False,
          ["no input files found", "_pacbio"])
    # The hybrid example has Illumina reads for SAM1 and SAM2, and long reads for SAM1 (ONT)
    # and SAM2 (PacBio). Each of the next two cases replaces one side with SAM1 alone, so
    # SAM2 is missing from that side and the message must name it.
    check("hybrid sample without long reads", variant("mode = hybrid", [],
          {"longread_dir": ["SAM1_ont.fastq.gz"]}), False,
          ["needs Illumina and long reads for every sample", "No _ont or _pacbio file for: SAM2"])
    check("hybrid sample without Illumina reads", variant("mode = hybrid", [],
          {"illumina_dir": ["SAM1_R1.fastq.gz", "SAM1_R2.fastq.gz"]}), False,
          ["No Illumina _R1 file for: SAM2"])
    # A batch missing a different sample on each side: SAM3 has Illumina reads only, SAM4
    # long reads only. Both must be named, each under its own line.
    check("hybrid batch with a sample missing on each side", variant("mode = hybrid", [],
          {"illumina_dir": ["SAM1_R1.fastq.gz", "SAM1_R2.fastq.gz", "SAM2_R1.fastq.gz", "SAM2_R2.fastq.gz",
                            "SAM3_R1.fastq.gz", "SAM3_R2.fastq.gz"],
           "longread_dir": ["SAM1_ont.fastq.gz", "SAM2_pacbio.fastq.gz", "SAM4_ont.fastq.gz"]}), False,
          ["needs Illumina and long reads for every sample", "No _ont or _pacbio file for: SAM3",
           "No Illumina _R1 file for: SAM4"])
    # The space in "bad name" is one of the characters 00_common.smk rejects in a sample name.
    check("forbidden character in a sample name", variant("mode = contigs", [], {"contigs_dir": ["bad name.fasta"]}), False,
          ["contains characters that are not allowed"])
    # The annotation table (files.annotation_params) and the protein files it names. none.tsv
    # and missing.faa are never created.
    check("annotation table missing", variant("mode = illumina", [("files.annotation_params", base + "/ref/none.tsv")]), False,
          ["files.annotation_params", "not found"])
    # The longread example has samples SAM1 and SAM2; this table has a row for SAM1 only.
    check("sample without a row in the annotation table", variant("mode = longread",
          [("files.annotation_params", make_params([("SAM1", base + "/ref/proteins.faa")]))]), False,
          ["sample 'SAM2' has no row"])
    check("protein file named in the table is empty", variant("mode = illumina",
          [("files.annotation_params", make_params([("SAM1", base + "/ref/empty.faa")]))]), False,
          ["is empty"])
    check("protein file named in the table does not exist", variant("mode = illumina",
          [("files.annotation_params", make_params([("SAM1", base + "/ref/missing.faa")]))]), False,
          ["does not exist"])
    check("annotation table row with the wrong number of columns", variant("mode = illumina",
          [("files.annotation_params", base + "/ref/badcols.tsv")]), False,
          ["has 2 tab-separated columns; expected 3 or 4"])
    check("annotation table row with five columns", variant("mode = illumina",
          [("files.annotation_params", base + "/ref/fivecols.tsv")]), False,
          ["has 5 tab-separated columns; expected 3 or 4"])
    # The AUGUSTUS model of column 4. With the database stand-in on disk, a misspelt name
    # is refused and the closest stored name is suggested; funannotate itself would train
    # AUGUSTUS under the misspelt name without a word. Capitals count, as in funannotate.
    check("AUGUSTUS model the database does not hold", variant("mode = illumina",
          [("directories.funannotate_db", base + "/ref/fundb_species"),
           ("files.annotation_params", make_params([("SAM1", base + "/ref/proteins.faa")], "fusarium_graminaerum"))]), False,
          ["sample 'SAM1' names the AUGUSTUS model 'fusarium_graminaerum'", "Did you mean 'fusarium_graminearum'?"])
    check("AUGUSTUS model written with capitals", variant("mode = illumina",
          [("directories.funannotate_db", base + "/ref/fundb_species"),
           ("files.annotation_params", make_params([("SAM1", base + "/ref/proteins.faa")], "Fusarium_graminearum"))]), False,
          ["names the AUGUSTUS model 'Fusarium_graminearum'", "Did you mean 'fusarium_graminearum'?"])
    # Without a database on disk, only the characters of the name can be checked.
    check("AUGUSTUS model name that cannot be a folder name", variant("mode = illumina",
          [("files.annotation_params", make_params([("SAM1", base + "/ref/proteins.faa")], "fusarium/graminearum"))]), False,
          ["names the AUGUSTUS model 'fusarium/graminearum'", "is not a possible model name"])
    # Database paths, tuning parameters, licensed tools and download links. The screen always
    # runs in the read-based modes, so illumina mode needs directories.blast_db.
    check("required database path empty", variant("mode = illumina", [("directories.eggnog_db", "")]), False,
          ["directories.eggnog_db: required"])
    check("BLAST database missing where the screen runs", variant("mode = illumina", [("directories.blast_db", DELETE)]), False,
          ["directories.blast_db: required"])
    # PacBio HiFi reads are now recognised from the file name, so the message must say how to
    # name them ({sample}_pacbio.fastq.gz).
    check("flye_input_mode: pacbio-hifi is no longer a value", variant("mode = longread", [("parameters.longread.flye_input_mode", "pacbio-hifi")]), False,
          ["'pacbio-hifi' is not a value any more", "{sample}_pacbio.fastq.gz"])
    check("unknown assembler", variant("mode = longread", [("parameters.longread.assembler", "canu")]), False,
          ["parameters.longread.assembler: must be hifiasm | flye"])
    # "mod" instead of "mode": an unknown key inside a block must be named, not ignored.
    check("typo inside a parameter block", variant("mode = illumina", [("parameters.mitochondria.mod", "off")]), False,
          ["parameters.mitochondria: unsupported key(s): mod"])
    # medaka_model at the top of `parameters`, where the v1 FunFluxL workflow kept it; v2
    # reads it from the block named after the mode.
    check("a setting left at its v1 place", variant("mode = longread", [("parameters.medaka_model", "auto")]), False,
          ["parameters.medaka_model: this setting moved into the block named after your mode"])
    check("unknown top-level parameter", variant("mode = illumina", [("parameters.masking", "tantan")]), False,
          ["parameters.masking: unknown setting"])
    check("invalid masking method", variant("mode = illumina", [("parameters.masking_method", "repeatmasker")]), False,
          ["parameters.masking_method: must be tantan | repeatmodeler_repeatmasker"])
    check("invalid decontamination mode", variant("mode = illumina", [("parameters.decontamination", {"mode": "yes"})]), False,
          ["parameters.decontamination.mode: must be one of"])
    # Decontamination mode auto was removed from FunFlux, run-wide and in the
    # sample_overrides table; the message must say why (auto can remove most of the
    # genome) and what to do instead (start with off).
    check("decontamination mode auto is refused", variant("mode = illumina",
          [("parameters.decontamination", {"mode": "auto"})]), False,
          ["parameters.decontamination.mode: 'auto' is not available in FunFlux", "auto can remove most of the genome",
           "Use off (the default)"])
    check("mode auto in the sample_overrides table is refused", variant("mode = illumina",
          [("parameters.decontamination", {"sample_overrides": base + "/ref/overrides_auto.tsv"})]), False,
          ["parameters.decontamination.sample_overrides: sample 'SAM1' has mode 'auto'", "overrides_auto.tsv",
           "not available in FunFlux"])
    # The decontamination tables are checked before any job starts, but only where the
    # screen runs, so these cases use illumina mode. include mode needs a genus for every
    # sample it applies to: here the run-wide mode is include, and in the next case an
    # override row sets include for SAM1.
    check("include mode without any genus is refused", variant("mode = illumina",
          [("parameters.decontamination", {"mode": "include"})]), False,
          ["mode include keeps only the genera you name", "no genus is given for sample(s) SAM1"])
    check("include mode from an override row without a genus is refused", variant("mode = illumina",
          [("parameters.decontamination", {"sample_overrides": base + "/ref/overrides_include_no_genus.tsv"})]), False,
          ["no genus is given for sample(s) SAM1"])
    # no_overrides.tsv is never created.
    check("decontamination table that does not exist", variant("mode = illumina",
          [("parameters.decontamination", {"sample_overrides": base + "/ref/no_overrides.tsv"})]), False,
          ["parameters.decontamination.sample_overrides", "is not an existing file"])
    check("decontamination table without its header line", variant("mode = illumina",
          [("parameters.decontamination", {"include_genera_by_sample": base + "/ref/include_without_header.tsv"})]), False,
          ["parameters.decontamination.include_genera_by_sample", "header line that has the columns sample and genus"])
    # The contig selector reads the tables as UTF-8 text only; another encoding must stop
    # the dry run with a plain message, not a Python error.
    check("decontamination table that is not UTF-8 text", variant("mode = illumina",
          [("parameters.decontamination", {"include_genera_by_sample": base + "/ref/include_latin1.tsv"})]), False,
          ["parameters.decontamination.include_genera_by_sample", "is not UTF-8 text"], ["Traceback"])
    # SAM1's override row has discard_no_hit "maybe", which the selector cannot read as
    # yes or no.
    check("override row with a discard_no_hit that is not yes or no", variant("mode = illumina",
          [("parameters.decontamination", {"sample_overrides": base + "/ref/overrides_bad_discard.tsv"})]), False,
          ["sample 'SAM1' has discard_no_hit 'maybe'", "the discard_no_hit column must be true or false"])
    # parameters.contigs.screen acts only in contigs mode, but its value is checked in every
    # mode.
    check("contigs.screen with a wrong value in illumina mode", variant("mode = illumina",
          [("parameters.contigs", {"screen": "maybe"})]), False,
          ["parameters.contigs.screen: expected true or false"])
    # Python reads the word nan as a number ("not a number") that passes every range check;
    # it must be refused like any other value that is not a number.
    check("nan as a decimal setting", variant("mode = illumina", [("parameters.its_taxonomy_cutoff", "nan")]), False,
          ["parameters.its_taxonomy_cutoff: expected a number"])
    # The boolean True (what YAML makes of an unquoted true, yes or on) is not a mitochondria
    # mode, and the message must show the value it got.
    check("mitochondria mode written as a YAML true", variant("mode = illumina", [("parameters.mitochondria", {"mode": True})]), False,
          ["parameters.mitochondria.mode: must be one of off | report | separate_confirmed", "got: True"])
    # A BUSCO lineage must look like hypocreales_odb12; a bare order name is refused.
    check("invalid BUSCO lineage", variant("mode = illumina", [("parameters.busco.lineage", "Hypocreales")]), False,
          ["parameters.busco.lineage"])
    # busco --list-datasets prints each name with its version after a dot. A name copied with
    # that dot is refused, and the message says to leave the dot and the number out.
    check("BUSCO lineage copied with its version", variant("mode = illumina", [("parameters.busco.lineage", "hypocreales_odb12.2")]), False,
          ["parameters.busco.lineage: must be 'auto' or a BUSCO lineage name", "got: 'hypocreales_odb12.2'",
           "leave out the dot and the number after it"])
    # The contig summary (parameters.contig_summary). Its telomere motif is the word auto or a
    # DNA word of A, C, G and T, 4 to 30 letters long, which tidk then searches for: an N (any
    # base) or a three-letter word is refused. The checksum of the telomere database must look
    # like a SHA-256 (64 hexadecimal characters). The window has a floor and a ceiling (500 and
    # 100000), the two counts have floors (2 and 1), and a misspelt key stops the run like in
    # every other block.
    check("telomere motif with a letter that is not A, C, G or T", variant("mode = illumina",
          [("parameters.contig_summary", {"telomere_motif": "TTAGGN"})]), False,
          ["parameters.contig_summary.telomere_motif: must be auto (the unit is chosen from the data) "
           "or the telomere repeat unit of the species", "4 to 30 letters A, C, G and T", "'TTAGGN'"])
    check("telomere motif shorter than 4 letters", variant("mode = illumina",
          [("parameters.contig_summary", {"telomere_motif": "TTA"})]), False,
          ["parameters.contig_summary.telomere_motif: must be auto (the unit is chosen from the data) "
           "or the telomere repeat unit of the species", "'TTA'"])
    check("telomere database checksum that is not a SHA-256", variant("mode = illumina",
          [("parameters.contig_summary", {"telomere_db_sha256": "1d4f5c8c"})]), False,
          ["parameters.contig_summary.telomere_db_sha256: must be the SHA-256 checksum of the telomere database file",
           "64 hexadecimal characters", "'1d4f5c8c'"])
    check("telomere window below 500", variant("mode = illumina",
          [("parameters.contig_summary", {"telomere_window": 200})]), False,
          ["parameters.contig_summary.telomere_window: must be at least 500 (got 200)"])
    check("telomere window above 100000", variant("mode = illumina",
          [("parameters.contig_summary", {"telomere_window": 200000})]), False,
          ["parameters.contig_summary.telomere_window: must be at most 100000 (got 200000)"])
    check("min_telomere_copies below 2", variant("mode = illumina",
          [("parameters.contig_summary", {"min_telomere_copies": 1})]), False,
          ["parameters.contig_summary.min_telomere_copies: must be at least 2 (got 1)"])
    check("min_rdna_units below 1", variant("mode = illumina",
          [("parameters.contig_summary", {"min_rdna_units": 0})]), False,
          ["parameters.contig_summary.min_rdna_units: must be at least 1 (got 0)"])
    check("typo inside the contig_summary block", variant("mode = illumina",
          [("parameters.contig_summary", {"telomere_motive": "TTAGGG"})]), False,
          ["parameters.contig_summary: unsupported key(s): telomere_motive"])
    # The funannotate database stand-in (ref/fundb) has no fungi folder. The message must
    # give the setup command that worked when tested (-i busco, -w).
    check("BUSCO set missing from the funannotate database", variant("mode = illumina",
          [("directories.funannotate_db", base + "/ref/fundb"), ("parameters.prediction", {"busco_db": "fungi"})]), False,
          ["parameters.prediction.busco_db: 'fungi' is not in the funannotate database",
           "funannotate setup -i busco -b fungi -d " + base + "/ref/fundb -w"])
    check("SPAdes k-mer list with an even number", variant("mode = illumina", [("parameters.spades_kmers", "21,32,55")]), False,
          ["parameters.spades_kmers"])
    # auto, in any case, hands the k-mer choice to SPAdes: the command then reads -k auto.
    check("SPAdes k-mers set to auto", variant("mode = illumina", [("parameters.spades_kmers", "Auto")]), True,
          ["spades.py -k auto --isolate", done])
    check("SPAdes k-mers neither auto nor numbers", variant("mode = illumina", [("parameters.spades_kmers", "automatic")]), False,
          ["parameters.spades_kmers: expected auto, or whole numbers separated by commas"])
    # The time limit of the long steps on a cluster, in whole hours, at least 1.
    check("runtime_hours below 1", variant("mode = illumina", [("resources.runtime_hours", 0)]), False,
          ["resources.runtime_hours: must be at least 1 (got 0)"])
    check("runtime_hours that is not a whole number", variant("mode = illumina", [("resources.runtime_hours", 1.5)]), False,
          ["resources.runtime_hours: expected a whole number (got 1.5)"])
    # The scratch folder itself exists but holds no gmes_petap.pl.
    check("GeneMark directory without gmes_petap.pl", variant("mode = illumina", [("directories.genemark_dir", base)]), False,
          ["has no executable gmes_petap.pl"])
    check("SignalP 6 on without an archive", variant("mode = illumina", [("parameters.signalp6.mode", "fast")]), False,
          ["files.signalp6_archive: required"])
    # Tiberius is on by default. Its mode takes only on or off, its checksum must look like a
    # SHA-256 (64 characters 0-9 and a-f), and while it is on the model link is required.
    check("invalid Tiberius mode", variant("mode = illumina", [("parameters.tiberius", {"mode": "maybe"})]), False,
          ["parameters.tiberius.mode: must be on | off", "'maybe'"])
    check("Tiberius checksum that is not a SHA-256", variant("mode = illumina",
          [("parameters.tiberius", {"model_sha256": "4a309a63"})]), False,
          ["parameters.tiberius.model_sha256: expected the SHA-256 checksum"])
    check("Tiberius on without its model link", variant("mode = illumina", [("links.tiberius_model_link", DELETE)]), False,
          ["links.tiberius_model_link: required"])
    # Without links.unite_its_link and without directories.unite_db there is no UNITE
    # reference at all.
    check("no UNITE source at all", variant("mode = illumina", [("links.unite_its_link", DELETE)]), False,
          ["links.unite_its_link: required"])
    # The telomere database is downloaded in every mode (the contig summary runs in all of
    # them), so its link is required everywhere; contigs mode is the smallest case.
    check("telomere database link missing", variant("mode = contigs", [("links.telomere_db_link", DELETE)]), False,
          ["links.telomere_db_link: required"])
    check("funannotate link is not a Zenodo URL", variant("mode = illumina", [("links.funannotate_link", "https://example.org/db.tar.gz")]), False,
          ["expected a Zenodo URL"])
    check("Filtlong keep_percent out of range", variant("mode = longread", [("parameters.longread.filtlong", {"keep_percent": 150})]), False,
          ["parameters.longread.filtlong.keep_percent"])
    # split is not a key of the filtlong block (the hybrid block has its own filtlong_split).
    check("old Filtlong split key", variant("mode = longread", [("parameters.longread.filtlong", {"split": 1000})]), False,
          ["parameters.longread.filtlong: unsupported key(s): split"])
    # The ref/ stand-in folder exists but holds no clusterblast/ folder, so it is not an
    # antiSMASH database.
    check("optional antiSMASH database that is not there", variant("mode = illumina", [("directories.antismash_db", base + "/ref")]), False,
          ["directories.antismash_db", "has no 'clusterblast'"])
    # The AMR screen (parameters.amr). Its switch takes only on or off; the mitochondrial
    # genetic code is auto or the number of an NCBI genetic code (1-6, 9-16, 21-33: 7, 8
    # and 17 to 20 are not NCBI codes, and miniprot refuses them); the two checksums must look like a SHA-256 (64 characters 0-9 and a-f); a
    # misspelt key stops the run as in every other block; and while the step is on, the
    # address ChroQueTas and its database are downloaded from is required. The example
    # blocks name no such address, so switching the step on is enough for the last case.
    check("invalid AMR switch", variant("mode = illumina", [("parameters.amr", {"chroquetas": "maybe"})]), False,
          ["parameters.amr.chroquetas: must be on | off", "'maybe'"])
    check("mitochondrial genetic code 0", variant("mode = illumina", [("parameters.amr", {"mito_genetic_code": 0})]), False,
          ["parameters.amr.mito_genetic_code: must be auto or the number of an NCBI genetic code: "
           "1-6, 9-16 or 21-33 (got: 0)"])
    check("mitochondrial genetic code 34", variant("mode = illumina", [("parameters.amr", {"mito_genetic_code": 34})]), False,
          ["parameters.amr.mito_genetic_code: must be auto or the number of an NCBI genetic code", "(got: 34)"])
    check("mitochondrial genetic code 7, not an NCBI code", variant("mode = illumina",
          [("parameters.amr", {"mito_genetic_code": 7})]), False,
          ["parameters.amr.mito_genetic_code: must be auto or the number of an NCBI genetic code", "(got: 7)"])
    check("mitochondrial genetic code written as a word", variant("mode = illumina",
          [("parameters.amr", {"mito_genetic_code": "mold"})]), False,
          ["parameters.amr.mito_genetic_code: must be auto or the number of an NCBI genetic code", "'mold'"])
    check("mitochondrial genetic code with decimals", variant("mode = illumina",
          [("parameters.amr", {"mito_genetic_code": 4.5})]), False,
          ["parameters.amr.mito_genetic_code: must be auto or the number of an NCBI genetic code", "(got: 4.5)"])
    check("ChroQueTas script checksum that is not a SHA-256", variant("mode = illumina",
          [("parameters.amr", {"chroquetas_script_sha256": "d8961a0d"})]), False,
          ["parameters.amr.chroquetas_script_sha256: must be the SHA-256 checksum of ChroQueTas.sh",
           "64 hexadecimal characters", "'d8961a0d'"])
    check("FungAMR database checksum that is not a SHA-256", variant("mode = illumina",
          [("parameters.amr", {"chroquetas_db_sha256": "6f89ecb3"})]), False,
          ["parameters.amr.chroquetas_db_sha256: must be the SHA-256 checksum of the FungAMR database archive",
           "'6f89ecb3'"])
    check("typo inside the amr block", variant("mode = illumina", [("parameters.amr", {"chroquetaz": "on"})]), False,
          ["parameters.amr: unsupported key(s): chroquetaz"])
    check("AMR screen on without its link", variant("mode = contigs", [("parameters.amr", {"chroquetas": "on"})]), False,
          ["links.chroquetas_link: required"])

    print("=== the job graph each setting builds (contigs mode) ===")
    # The job graph is the set of jobs a dry run plans. Snakemake lists them in its "Job
    # stats" table, one line per rule with the number of its jobs, so a rule's name in the
    # output means the rule is in the graph, and expect_jobs (see check()) checks the
    # numbers. The optional features are tested in contigs mode, the smallest graph; the
    # read-based modes have their own group below. Names that are the beginning of
    # another name (antismash_db, antismash_db_local) are told apart by the job block of
    # the rule, which begins "rule <name>:".
    #
    # The one contigs sample of the example gives 31 jobs: 10 that run once per run (rule
    # all, the UNITE download and its reference filter, the antiSMASH download, the funannotate
    # database and its path normalisation, the Tiberius model, the telomere database, the
    # UniVec download of the adapter screen, MultiQC) and 21 for the sample, among them the
    # telomere repeat search (rule telomere_units), the adapter screen (rules blast_adapters
    # and screen_adapters) and, last, the contig summary (rule contig_summary; all of them in
    # every mode). The screen adds 8
    # jobs for the sample (its minimap2 self-mapping, BLAST, the two BlobTools steps, the
    # mitochondrion screen, the contig selection, finalize_assembly and the mitochondrion
    # split); Tiberius off takes 3 away (the model download and the two per-sample Tiberius
    # jobs); the UNITE reference filter off takes 1 away; the AMR screen, off by default, adds
    # 2 when it is on (the download of ChroQueTas and its database, and the sample's
    # ChroQueTas job); a named BUSCO lineage adds 1 (its download, rule busco_data).
    contigs_plain = blocks["mode = contigs"]
    # The output folder of the plain contigs example block is <scratch>/out/contigs, so a file
    # a command names can be written out in full. The variants below start from that block
    # (the screen case switches the screen on in it), so they write to the same folder.
    out_contigs = base + "/out/contigs"
    # The edits that switch the AMR screen on: the switch itself and the address of the pinned
    # ChroQueTas commit (the same address as in config/config.yaml), which the example
    # blocks do not carry. Used here and in the read-based modes below.
    chroquetas_link = "https://raw.githubusercontent.com/nmquijada/ChroQueTas/1aad051d6ad48f5e403636699d19095e655764f7"
    amr_on = [("parameters.amr", {"chroquetas": "on"}), ("links.chroquetas_link", chroquetas_link)]
    # The plain example plans the download rules for UNITE, antiSMASH and the telomere
    # database, the UNITE reference filter, the telomere repeat search, and none of the
    # optional rules. The telomere repeat search reads the draft and the downloaded database.
    # The contig summary gets the inputs every mode has (the ITS windows, the BUSCO folder,
    # the annotated genome with its GFF3, named after the species and the sample as
    # funannotate names it, the telomere unit and the telomere ends table) and none of the
    # five that depend on the mode: without the screen there are no screen files, and
    # without reads no depth table and no assembler table. Snakemake prints each command
    # on one line, so the five empty options stand next to each other.
    # BUSCO in its default auto mode downloads its lineage data inside the job, into the
    # sample's own busco_tmp folder, so there is no busco_data job. With SignalP 6 off,
    # funannotate annotate gets no --signalp option: --cpus follows --eggnog directly.
    check("annotation only: no screen, SignalP, GeneMark or local-database rules (31 jobs)", contigs_plain, True,
          [f"--download_path {out_contigs}/02.assembly/SAM1/eval/busco/busco_tmp --auto-lineage-euk",
           '--eggnog "$eggnog_annotations" --cpus',
           "filter_contigs", "its_windows", "its_classification", "funannotate_prediction", "multiqc",
           "rule unite_db:", "rule unite_identified:", "rule antismash_db:", "rule telomere_db:",
           "rule telomere_units:", "rule contig_summary:",
           f"--telomere-db {out_contigs}/00.resources/telomeres/telo_database_v2.tsv",
           f"--contigs {out_contigs}/02.assembly/SAM1/contigs_final.fasta",
           "--mito-contigs --decisions --mito-screen --depth-table --assembly-evidence --its-windows",
           f"--busco-dir {out_contigs}/02.assembly/SAM1/eval/busco",
           f"--gff3 {out_contigs}/04.annotation/SAM1/funannotate/annotate_results/Fusarium_graminearum_SAM1.gff3",
           f"--annotated-fasta {out_contigs}/04.annotation/SAM1/contigs_sort.fasta",
           f"--telomere-unit-file {out_contigs}/02.assembly/SAM1/eval/telomere_unit.txt",
           f"--telomere-ends {out_contigs}/02.assembly/SAM1/eval/telomere_ends.tsv", done],
          ["blast_contigs", "separate_mitochondria", "install_signalp6", "genemark_check",
           "unite_db_local", "antismash_db_local", "rule busco_data:", "--offline", "signalp_args"],
          expect_jobs={"total": 31, "filter_contigs": 1, "telomere_db": 1, "telomere_units": 1, "contig_summary": 1,
                       "multiqc": 1, "unite_db": 1, "unite_identified": 1, "antismash_db": 1, "funannotate_db": 1,
                       "normalize_funannotate_db": 1, "tiberius_model": 1})
    # The contigs example has no BLAST database path, so the screen case adds the stand-in
    # folder that holds nodes.dmp and names.dmp (see build_stubs). With the screen, the contig
    # summary reads the draft (contigs_filt.fasta) and gets the three screen files, but still
    # no depth table: in contigs mode the screen maps the contigs to themselves, and that
    # depth means nothing (the rule's comment in 00_common.smk, section 9h).
    check("the screen adds the screen rules (39 jobs)", variant("mode = contigs",
          [("parameters.contigs", {"screen": True}), ("directories.blast_db", base + "/ref/blastdb")]), True,
          ["map_contigs", "blast_contigs", "blob_json", "blob_table", "screen_mitochondria", "select_contigs",
           "finalize_assembly", "separate_mitochondria",
           f"--contigs {out_contigs}/02.assembly/SAM1/contigs_filt.fasta",
           f"--mito-contigs {out_contigs}/02.assembly/SAM1/contigs_mitochondrial.fasta",
           f"--decisions {out_contigs}/02.assembly/SAM1/contaminants/contig_taxonomy_decisions.tsv",
           f"--mito-screen {out_contigs}/02.assembly/SAM1/mitochondria/mitochondrial_screen.tsv "
           "--depth-table --assembly-evidence --its-windows", done],
          expect_jobs={"total": 39, "map_contigs": 1, "blast_contigs": 1, "blob_json": 1, "blob_table": 1,
                       "screen_mitochondria": 1, "select_contigs": 1, "finalize_assembly": 1,
                       "separate_mitochondria": 1, "telomere_units": 1, "contig_summary": 1})
    # SignalP 6 on: funannotate annotate gets --signalp and the sample's SignalP 6 table,
    # built in Snakemake (an empty bash list would stop bash older than 4.4).
    check("SignalP 6 adds its two rules", variant("mode = contigs",
          [("parameters.signalp6", {"mode": "fast"}), ("files.signalp6_archive", base + "/ref/signalp6.tar.gz")]), True,
          ["install_signalp6", "signalp6_annotation",
           f'--eggnog "$eggnog_annotations" --signalp {out_contigs}/04.annotation/SAM1/signalp6/prediction_results.txt --cpus',
           done])
    # The two GeneMark cases run with the fake home folder, so the stand-in licence key is
    # found there.
    # The check starts gmes_petap.pl as a program, the way funannotate predict does, so that
    # the Perl named in its first line answers, not the environment's Perl.
    check("GeneMark adds the licence check", variant("mode = contigs",
          [("directories.genemark_dir", base + "/ref/genemark")]), True,
          ["genemark_check", f'status=0 "{base}/ref/genemark/gmes_petap.pl" --ES --fungus', done],
          [f'perl "{base}/ref/genemark/gmes_petap.pl"'], home=base + "/home")
    check("GeneMark set to off adds nothing", variant("mode = contigs",
          [("directories.genemark_dir", base + "/ref/genemark"), ("parameters.genemark", {"mode": "off"})]), True,
          [done], ["genemark_check"], home=base + "/home")
    # Tiberius is on by default, so the plain example plans its three rules; off removes them.
    check("Tiberius, on by default, adds its three rules", contigs_plain, True,
          ["Tiberius: on", "rule tiberius_model:", "rule tiberius_prediction:", "rule tiberius_gff3:", done])
    check("Tiberius set to off adds nothing (28 jobs)", variant("mode = contigs", [("parameters.tiberius", {"mode": "off"})]), True,
          ["Tiberius: off.", done], ["tiberius_model", "tiberius_prediction", "tiberius_gff3"],
          expect_jobs={"total": 28})
    # A named BUSCO lineage is downloaded once per run by rule busco_data (busco --download,
    # into 00.resources/busco/), and the sample's BUSCO job reads it with --offline instead of
    # downloading into its own busco_tmp folder: 31 + 1 = 32 jobs.
    check("a named BUSCO lineage is downloaded once and read offline (32 jobs)", variant("mode = contigs",
          [("parameters.busco.lineage", "hypocreales_odb12")]), True,
          ["rule busco_data:",
           f"busco --download hypocreales_odb12 --download_path {out_contigs}/00.resources/busco",
           f"--offline --download_path {out_contigs}/00.resources/busco -l hypocreales_odb12", done],
          ["busco_tmp --", "--auto-lineage-euk"],
          expect_jobs={"total": 32, "busco_data": 1, "completeness": 1})
    # The contig adapter screen runs in contigs mode too. The example block has no
    # links.univec_link, so UniVec build 10.1 in NCBI's archive is downloaded and checked
    # against the default checksum; the contig summary gets the screen's loci.
    check("contig adapter screen: UniVec 10.1 by default, loci go to the contig summary", contigs_plain, True,
          ["https://ftp.ncbi.nlm.nih.gov/pub/UniVec/archive/10.1/UniVec",
           "687d8aa872d4ee215c9e0b3b004d0ee709fb243436d36d159a3f0ac5c25b5e36",
           f"--adapter-screen {out_contigs}/02.assembly/SAM1/eval/adapter_screen.tsv", done],
          expect_jobs={"univec_illumina": 1, "blast_adapters": 1, "screen_adapters": 1})
    check("univec_sha256 that is not a checksum is refused", variant("mode = contigs",
          [("parameters.contig_summary", {"univec_sha256": "not-a-checksum"})]), False,
          ["parameters.contig_summary.univec_sha256: must be the SHA-256 checksum of the UniVec file"])
    check("local antiSMASH and UNITE copies replace the downloads", variant("mode = contigs",
          [("directories.antismash_db", base + "/ref/antismash_db"), ("directories.unite_db", base + "/ref/unite_db")]), True,
          ["rule antismash_db_local:", "rule unite_db_local:", done], ["rule antismash_db:", "rule unite_db:"])
    # The UNITE reference filter (parameters.its.unite_identified_only) is on by default: rule
    # unite_identified is planned once, and the vsearch command of its_classification reads
    # the filtered file. The contigs example writes to <scratch>/out/contigs, so the --db
    # path can be named in full. With the key off the rule is not planned and vsearch reads
    # the whole release; a value that is not true or false stops the run (checked in every
    # mode, so the illumina example serves for that).
    check("UNITE reference filter, on by default: unite_identified planned, SINTAX reads the filtered file (31 jobs)",
          contigs_plain, True,
          ["rule unite_identified:", "unite_identified.py",
           f"--db {out_contigs}/00.resources/unite/unite_its_sintax.identified.fasta", done],
          [f"--db {out_contigs}/00.resources/unite/unite_its_sintax.fasta "],
          expect_jobs={"total": 31, "unite_identified": 1, "unite_db": 1})
    check("UNITE reference filter off: unite_identified not planned, SINTAX reads the whole release (30 jobs)",
          variant("mode = contigs", [("parameters.its.unite_identified_only", False)]), True,
          ["rule unite_db:", f"--db {out_contigs}/00.resources/unite/unite_its_sintax.fasta ", done],
          ["rule unite_identified:", "unite_identified.py", "unite_its_sintax.identified.fasta"],
          expect_jobs={"total": 30, "unite_identified": 0, "unite_db": 1})
    check("UNITE reference filter with a value that is not true or false",
          variant("mode = illumina", [("parameters.its.unite_identified_only", "maybe")]), False,
          ["parameters.its.unite_identified_only: expected true or false", "'maybe'"])
    # The AMR screen is off by default: the plain example plans neither of its rules and names
    # no 05.amr folder. MultiQC writes into 06.report, and no command names a 05.report
    # folder.
    check("AMR screen off by default: no AMR rules, no 05.amr, the report in 06.report (31 jobs)", contigs_plain, True,
          ["AMR screen: off.", f'--outdir "{out_contigs}/06.report"', done],
          ["rule chroquetas_db:", "rule amr_chroquetas:", "05.amr", "05.report"],
          expect_jobs={"total": 31, "chroquetas_db": 0, "amr_chroquetas": 0, "multiqc": 1})
    # On, it adds the download (rule chroquetas_db: the script, the database archive and the
    # licence from the pinned commit, each checked file against its default checksum) and
    # the sample's job (rule amr_chroquetas). The wrapper gets the downloaded script and the
    # database folder, the nuclear genome, the species of the annotation table (quoted,
    # because it holds a space), the mitochondrial code auto and its four output paths in
    # 05.amr/SAM1/. Without the screen nothing is split off, so --mito-contigs stands empty:
    # the mitochondrion, if any, is among the contigs of contigs_final.fasta. MultiQC reads
    # the AMR table after the contig summary.
    check("AMR screen on: the download and one ChroQueTas job per sample (33 jobs)", variant("mode = contigs", amr_on), True,
          ["AMR screen: ChroQueTas on (licence CC BY-NC-ND 4.0: non-commercial use).",
           "rule chroquetas_db:", "rule amr_chroquetas:",
           f'"{chroquetas_link}/bin/ChroQueTas.sh"', f'"{chroquetas_link}/FungAMR_db/db-v20250811.tgz"',
           f'"{chroquetas_link}/LICENSE"',
           f"d8961a0d311be57d4a9939bd9e4d90dbc796589b3298dd1835359ed776abfce3 {out_contigs}/00.resources/chroquetas/ChroQueTas.sh",
           f"6f89ecb3117110aad85b8b6ca7b8a3d3715cb5984b6cb75450a6cda5633d16e7 {out_contigs}/00.resources/chroquetas/db-v20250811.tgz",
           f"--strip-components=1 -C \"{out_contigs}/00.resources/chroquetas/FungAMR_db\"",
           "chroquetas_amr.py",
           f"--chroquetas-script {out_contigs}/00.resources/chroquetas/ChroQueTas.sh",
           f"--database {out_contigs}/00.resources/chroquetas/FungAMR_db",
           f"--genome {out_contigs}/02.assembly/SAM1/contigs_final.fasta --mito-contigs --species 'Fusarium graminearum' "
           "--mito-genetic-code auto",
           f"--sample SAM1 --work-dir {out_contigs}/05.amr/SAM1/chroquetas",
           f"--out-summary {out_contigs}/05.amr/SAM1/amr_summary.tsv",
           f"--out-positions {out_contigs}/05.amr/SAM1/amr_positions.tsv",
           f"--out-mqc {out_contigs}/05.amr/SAM1/amr_summary_mqc.tsv",
           "02.assembly/SAM1/eval/contig_summary_mqc.tsv 02.assembly/SAM1/eval/adapter_screen_mqc.tsv 05.amr/SAM1/amr_summary_mqc.tsv", done],
          expect_jobs={"total": 33, "chroquetas_db": 1, "amr_chroquetas": 1, "multiqc": 1})
    # With the screen (mitochondria mode separate_confirmed, the default) the mitochondrion is
    # split off into contigs_mitochondrial.fasta, which the job then waits for and hands to
    # the mitochondrial run: 39 + 2 = 41 jobs. In report mode the screen splits nothing off,
    # so --mito-contigs stands empty again and the mitochondrial run reads contigs_final.fasta.
    # (Rule contig_summary has an option of the same name, which still names the file in
    # report mode, so the text that must not appear includes the AMR script's next option.)
    contigs_screen_amr = [("parameters.contigs", {"screen": True}), ("directories.blast_db", base + "/ref/blastdb")] + amr_on
    check("AMR screen with the screen: the mitochondrial run reads the split-off mitochondrion (41 jobs)",
          variant("mode = contigs", contigs_screen_amr), True,
          [f"--genome {out_contigs}/02.assembly/SAM1/contigs_final.fasta "
           f"--mito-contigs {out_contigs}/02.assembly/SAM1/contigs_mitochondrial.fasta --species", done],
          expect_jobs={"total": 41, "chroquetas_db": 1, "amr_chroquetas": 1, "separate_mitochondria": 1})
    check("AMR screen with the mitochondrion screen in report mode: --mito-contigs stands empty (41 jobs)",
          variant("mode = contigs", contigs_screen_amr + [("parameters.mitochondria", {"mode": "report"})]), True,
          ["Mitochondrion screen: report", "--mito-contigs --species 'Fusarium graminearum'", done],
          [f"--mito-contigs {out_contigs}/02.assembly/SAM1/contigs_mitochondrial.fasta --species"],
          expect_jobs={"total": 41, "amr_chroquetas": 1})

    print("=== the job graph of the read-based modes ===")
    # The same kind of check for illumina, longread and hybrid mode, with the number of
    # jobs each time. How a count comes about: the jobs that run once per run (rule all,
    # the UNITE download and its reference filter, the antiSMASH download, the funannotate
    # database and its path normalisation, the Tiberius model, the telomere database, the
    # UniVec download, MultiQC: 10; the modes with Illumina reads add the PhiX download and its
    # index: 12), plus the per-sample jobs once per sample, plus,
    # for every ONT sample that uses Medaka, two more: the model check
    # (rule check_medaka_model) and the consensus (rule long_read_consensus), and for every
    # ONT sample whose reads are trimmed (parameters.<mode>.ont_adapter not off, the default),
    # three more: the adapter scan (rule find_ont_adapter), fastplong (rule
    # trim_ont_adapters) and the count of adapters left (rule check_ont_trimming). The AMR screen
    # is off in the examples; switched on, it adds one job per run (rule chroquetas_db) and
    # one per sample (rule amr_chroquetas), as the last case of each mode shows. The example
    # blocks have one illumina sample (SAM1) and two long-read samples, SAM1 with ONT and
    # SAM2 with PacBio HiFi reads, so the Medaka rules are planned once, for SAM1, and
    # SAM2's model file (SAM2_medaka_model.txt) is never asked for. Because the dry runs
    # print each job's command (-p), a case can also check the value a rule hands a tool.
    # The output folders of the illumina, longread and hybrid example blocks are
    # <scratch>/out/illumina, <scratch>/out/longread and <scratch>/out/hybrid, so the file a
    # command copies or reads can be named in full.
    out_illumina = base + "/out/illumina"
    out_longread = base + "/out/longread"
    out_hybrid = base + "/out/hybrid"

    # illumina: 12 once-per-run jobs + 35 for the one sample = 47. The chain is the PhiX
    # removal (download, index, bowtie2 with the alignments sent to /dev/null), fastp,
    # SPAdes with the memory limit of resources.ram_gb (32 in the example), the contig
    # filter, bowtie2 of the trimmed pairs on the draft (the coverage track of the
    # screen), Qualimap on that BAM with a Java heap of the same 32 GB, samtools coverage of
    # that BAM (rule contig_depth, the read depth for the contig summary), the screen, the
    # telomere repeat search on the draft, and no long-read rule at all. The contig summary
    # of this mode gets the three screen files and the depth table (real reads were
    # mapped), but no assembler table (SPAdes writes none): its --assembly-evidence stands
    # empty before --its-windows.
    check("illumina example plans the whole chain (one sample: 47 jobs)", blocks["mode = illumina"], True,
          ["rule download_phix:", "rule build_phix:", "rule map_phix:", "rule trim_adapters:",
           "rule illumina_assembly:", "rule filter_contigs:", "rule index_contigs:", "rule map_contigs:",
           "rule map_evaluation:", "rule screen_mitochondria:", "rule select_contigs:",
           "rule finalize_assembly:", "rule separate_mitochondria:", "rule telomere_db:",
           "rule telomere_units:", "rule contig_summary:", "rule contig_depth:",
           "-S /dev/null --local", "spades.py -k 21,33,55,77,99,127 --isolate", " -m 32 ",
           "--java-mem-size=32G",
           f"--contigs {out_illumina}/02.assembly/SAM1/contigs_filt.fasta",
           f"--depth-table {out_illumina}/02.assembly/SAM1/contaminants/SAM1_depth.tsv "
           "--assembly-evidence --its-windows", done],
          ["check_medaka_model", "long_read_consensus", "filter_long_reads", "minimap2"],
          expect_jobs={"total": 47, "download_phix": 1, "build_phix": 1, "map_phix": 1, "trim_adapters": 1,
                       "illumina_assembly": 1, "filter_contigs": 1, "index_contigs": 1, "map_contigs": 1,
                       "map_evaluation": 1, "blast_contigs": 1, "blob_json": 1, "blob_table": 1,
                       "screen_mitochondria": 1, "select_contigs": 1, "finalize_assembly": 1,
                       "separate_mitochondria": 1, "genome_assembly_evaluation": 1, "completeness": 1,
                       "telomere_db": 1, "telomere_units": 1, "contig_summary": 1, "contig_depth": 1,
                       "multiqc": 1, "chroquetas_db": 0, "amr_chroquetas": 0})
    # A second sample doubles the per-sample jobs and nothing else: 12 + 2 x 35 = 82.
    check("illumina with two samples (82 jobs)", variant("mode = illumina", [],
          {"illumina_dir": ["SAM1_R1.fastq.gz", "SAM1_R2.fastq.gz", "SAM2_R1.fastq.gz", "SAM2_R2.fastq.gz"]}), True,
          ["Sample SAM1 will be processed.", "Sample SAM2 will be processed.", done],
          expect_jobs={"total": 82, "download_phix": 1, "build_phix": 1, "map_phix": 2, "illumina_assembly": 2,
                       "telomere_db": 1, "telomere_units": 2, "contig_summary": 2, "contig_depth": 2,
                       "multiqc": 1})
    check("illumina with Tiberius off (44 jobs)", variant("mode = illumina", [("parameters.tiberius", {"mode": "off"})]), True,
          ["Tiberius: off.", done],
          expect_jobs={"total": 44, "tiberius_model": 0, "tiberius_prediction": 0, "tiberius_gff3": 0})

    # longread: 10 once-per-run jobs + 2 x 35 per-sample jobs + the two Medaka jobs and the
    # three adapter-trimming jobs of SAM1 = 85. The example names the Medaka model, so the check script gets it as it is
    # (--model "...") and may not fall back to working the model out from the reads
    # (--fallback-to-auto false, the default). Each sample's minimap2 preset and hifiasm
    # flag follow its reads: map-ont and --ont for SAM1, map-hifi and no flag for SAM2.
    # contigs_all.fasta of SAM1 is a copy of the Medaka consensus, SAM2's of the
    # decontaminated assembly (no Medaka for PacBio reads). The contig summary reads the
    # assembler's draft (assembly.fasta) and gets every optional file, the assembler's
    # evidence table included. medaka_consensus gets -x right after the draft (-d, SAM1's
    # decontaminated assembly): it rebuilds the draft's alignment index on every run. The
    # evidence table is made by a job of its own (rule assembly_evidence, one per sample),
    # which runs mitochondria.py evidence-hifiasm on the sample's primary graph; the
    # assembler's own command runs no Python (no "python" after its assembly.fasta). Both
    # NanoPlot reports of each sample are drawn without static PNG plots (--no_static), which
    # would need a Chrome browser on the computer.
    check("longread example: mixed batch with a named Medaka model (85 jobs)", blocks["mode = longread"], True,
          ["rule filter_long_reads:", "rule raw_long_read_qc:", "rule filtered_long_read_qc:",
           "rule check_medaka_model:", "rule long_read_assembly:", "rule assembly_evidence:",
           "rule long_read_consensus:",
           "rule map_contigs:", "rule map_evaluation:", "rule finalize_assembly:", "rule telomere_db:",
           "rule telomere_units:", "rule contig_summary:",
           "Medaka model: r1041_e82_400bps_sup_v4.3.0.",
           '--model "r1041_e82_400bps_sup_v4.3.0"', "--fallback-to-auto false", "SAM1_medaka_model.txt",
           "medaka_consensus -i", "-ax map-ont", "-ax map-hifi", "hifiasm -i --primary -f 0 --ont",
           "--loglength --no_static --prefix \"SAM1_\"", "--loglength --no_static --prefix \"SAM2_\"",
           f"-d {out_longread}/02.assembly/SAM1/contaminants/contigs_decontam.fasta -x -t",
           f"cp {out_longread}/02.assembly/SAM1/medaka/consensus.fasta",
           f"cp {out_longread}/02.assembly/SAM2/contaminants/contigs_decontam.fasta",
           f"--contigs {out_longread}/02.assembly/SAM1/assembly.fasta",
           f"--depth-table {out_longread}/02.assembly/SAM1/contaminants/SAM1_depth.tsv "
           f"--assembly-evidence {out_longread}/02.assembly/SAM1/assembly_evidence.tsv",
           f"mitochondria.py evidence-hifiasm --gfa {out_longread}/02.assembly/SAM1/hifiasm/hifiasm.p_ctg.gfa "
           f"--output {out_longread}/02.assembly/SAM1/assembly_evidence.tsv > ",
           f"mitochondria.py evidence-hifiasm --gfa {out_longread}/02.assembly/SAM2/hifiasm/hifiasm.p_ctg.gfa "
           f"--output {out_longread}/02.assembly/SAM2/assembly_evidence.tsv > ", done],
          ["SAM2_medaka_model.txt", "rule index_contigs:", "rule illumina_assembly:", "rule short_read_correction:",
           "download_phix", "--detect_adapter_for_pe", "spades.py", "assembly.fasta python"],
          expect_jobs={"total": 85, "filter_long_reads": 2, "raw_long_read_qc": 2, "filtered_long_read_qc": 2,
                       "long_read_assembly": 2, "assembly_evidence": 2, "check_medaka_model": 1,
                       "long_read_consensus": 1,
                       "map_contigs": 2, "map_evaluation": 2, "blast_contigs": 2, "blob_json": 2,
                       "screen_mitochondria": 2, "select_contigs": 2, "finalize_assembly": 2,
                       "separate_mitochondria": 2, "genome_assembly_evaluation": 2, "completeness": 2,
                       "telomere_db": 1, "telomere_units": 2, "contig_summary": 2, "contig_depth": 2,
                       "multiqc": 1,
                       "index_contigs": 0, "illumina_assembly": 0, "short_read_correction": 0,
                       "chroquetas_db": 0, "amr_chroquetas": 0})
    # Without a model name Medaka is still on: the check script works the model out from
    # the FASTQ headers (--model ""). The job count does not change.
    check("longread with the Medaka model worked out from the reads (85 jobs)", variant("mode = longread",
          [("parameters.longread.medaka_model", DELETE)]), True,
          ["Medaka model will be worked out from the FASTQ headers (checked before assembly).", '--model ""', done],
          ['--model "r1041'],
          expect_jobs={"total": 85, "check_medaka_model": 1, "long_read_consensus": 1})
    # medaka_model_fallback_auto reaches the check script as --fallback-to-auto true.
    check("longread with the fallback to auto switched on", variant("mode = longread",
          [("parameters.longread.medaka_model_fallback_auto", True)]), True,
          ['--model "r1041_e82_400bps_sup_v4.3.0"', "--fallback-to-auto true", done])
    # Medaka off: no model check, no consensus, and SAM1's contigs_all.fasta is a copy of
    # the decontaminated assembly like SAM2's: 85 - 2 = 83.
    check("longread with Medaka off (83 jobs)", variant("mode = longread", [("parameters.longread.medaka_model", False)]), True,
          ["Medaka disabled in the config. Skipping consensus polishing.",
           f"cp {out_longread}/02.assembly/SAM1/contaminants/contigs_decontam.fasta", done],
          ["check_medaka_model", "long_read_consensus", "medaka_consensus", "medaka_model_check.py"],
          expect_jobs={"total": 83, "check_medaka_model": 0, "long_read_consensus": 0})
    # A PacBio-only batch (SAM2 alone): the named model is ignored with a line that says
    # so, no Medaka rule and no adapter-trimming rule exists, and the one sample gives 10 + 35 = 45 jobs. The annotation
    # table still has its SAM1 row; a row for a sample that is not in the batch is fine.
    check("longread with a PacBio-only batch: no Medaka rules (45 jobs)", variant("mode = longread", [],
          {"longread_dir": ["SAM2_pacbio.fastq.gz"]}), True,
          ["Sample SAM2 will be processed.", "No ONT samples in this batch: any medaka_model setting is ignored.",
           "-ax map-hifi", done],
          ["Sample SAM1 will be processed.", "check_medaka_model", "long_read_consensus", "medaka_consensus",
           "-ax map-ont", "hifiasm -i --primary -f 0 --ont"],
          expect_jobs={"total": 45, "check_medaka_model": 0, "long_read_consensus": 0, "long_read_assembly": 1,
                       "assembly_evidence": 1})
    # Flye instead of hifiasm: the same jobs (one assembly rule either way), with Flye's
    # read mode per sample in the command: --nano-hq for ONT reads (the default model
    # means recent, high-accuracy reads), --pacbio-hifi for PacBio HiFi reads. The
    # evidence job runs evidence-flye on Flye's assembly_info.txt instead.
    check("longread with Flye (85 jobs)", variant("mode = longread",
          [("parameters.longread.assembler", "flye"), ("parameters.longread.medaka_model", DELETE)]), True,
          ["flye --nano-hq", "flye --pacbio-hifi", "--iterations 5",
           f"mitochondria.py evidence-flye --flye-info {out_longread}/02.assembly/SAM1/flye/assembly_info.txt "
           f"--output {out_longread}/02.assembly/SAM1/assembly_evidence.tsv > ", done],
          ["hifiasm -i", "evidence-hifiasm", "assembly.fasta python"],
          expect_jobs={"total": 85, "long_read_assembly": 2, "assembly_evidence": 2, "check_medaka_model": 1,
                       "long_read_consensus": 1})

    # ONT adapter trimming (parameters.<mode>.ont_adapter; rules/shared/05_ont_adapters.smk).
    # With the default, detect, the ONT sample SAM1 gets the three trimming jobs and the PacBio
    # sample SAM2 none; the scan is told to detect and reads ONT's published barcodes, fastplong
    # cuts exactly what it is given (--trimming_extension 0) and reads every sequence to trim
    # with -a, and Filtlong reads fastplong's output instead of the raw reads. The contig adapter
    # screen runs for both samples, and its Illumina records are downloaded once. SAM1's screen
    # also searches for the sequences trimmed from its reads, while SAM2's list of adapters ends
    # with the UniVec file (check() compares the output with every run of spaces made one).
    check("ONT adapter trimming by default: the ONT sample only, Filtlong reads the trimmed reads",
          blocks["mode = longread"], True,
          ["--setting detect", "--trimming_extension 0", "workflow/data/adapters/ont_chtd.fasta",
           f"-a {out_longread}/01.reads/SAM1/long/adapters/SAM1_adapter_sequence.fasta",
           f"{out_longread}/01.reads/SAM1/long/SAM1_trimmed.fastq.gz > {out_longread}/01.reads/SAM1/long/SAM1_filtered.fastq",
           f"univec_illumina.fasta {out_longread}/01.reads/SAM1/long/adapters/SAM1_adapter_sequence.fasta "
           f"> {out_longread}/02.assembly/SAM1/eval/adapter_query.fasta",
           f"univec_illumina.fasta > {out_longread}/02.assembly/SAM2/eval/adapter_query.fasta",
           done],
          ["SAM2_trimmed.fastq.gz", "SAM2_adapter_scan.tsv", "SAM2_adapter_sequence.fasta"],
          expect_jobs={"find_ont_adapter": 1, "trim_ont_adapters": 1, "check_ont_trimming": 1,
                       "univec_illumina": 1, "blast_adapters": 2, "screen_adapters": 2})
    # off: no trimming job, and Filtlong reads the raw reads, as in v2.0.0: 85 - 3 = 82.
    check("ONT adapter trimming off: no trimming jobs (82 jobs)", variant("mode = longread",
          [("parameters.longread.ont_adapter", "off")]), True,
          ["[FunFlux] SAM1: ont reads -> hifiasm --ont, minimap2 map-ont, Medaka.", done],
          ["rule find_ont_adapter:", "rule trim_ont_adapters:", "SAM1_trimmed.fastq.gz"],
          expect_jobs={"total": 82, "find_ont_adapter": 0, "trim_ont_adapters": 0, "check_ont_trimming": 0,
                       "blast_adapters": 2, "screen_adapters": 2})
    # A sequence is trimmed as given, in capitals; 68 bp is long enough for fastplong's read
    # splitting, so there is no NOTE.
    check("ONT adapter given as a sequence is used as it is", variant("mode = longread",
          [("parameters.longread.ont_adapter", "cctgtacttcgttcagttacgtattgctaaggttaagcatagttctgcatgatgggttagcagcacct")]),
          True, ["--setting CCTGTACTTCGTTCAGTTACGTATTGCTAAGGTTAAGCATAGTTCTGCATGATGGGTTAGCAGCACCT", done],
          ["NOTE: parameters.longread.ont_adapter"])
    check("ONT adapter shorter than 50 bp is accepted with a NOTE", variant("mode = longread",
          [("parameters.longread.ont_adapter", "CCTGTACTTCGTTCAGTTACGTATTGCT")]), True,
          ["NOTE: parameters.longread.ont_adapter is shorter than 50 bp", done])
    # A telomere repeat is never accepted as an adapter: trimming it would cut real genome.
    check("ONT adapter set to a telomere repeat is refused", variant("mode = longread",
          [("parameters.longread.ont_adapter", "TTAGGGTTAGGGTTAGGGTTAGGG")]), False,
          ["parameters.longread.ont_adapter: TTAGGGTTAGGGTTAGGGTTAGGG is a short unit repeated end to end"])
    check("ONT adapter with letters other than A, C, G, T is refused", variant("mode = longread",
          [("parameters.longread.ont_adapter", "CCTGTACTTCGTTCAGTTACNNNNNN")]), False,
          ["parameters.longread.ont_adapter: write detect, off, or the sequence to trim"])
    check("ONT adapter sequence in a PacBio-only batch is ignored with a NOTE", variant("mode = longread",
          [("parameters.longread.ont_adapter", "CCTGTACTTCGTTCAGTTACGTATTGCT")],
          {"longread_dir": ["SAM2_pacbio.fastq.gz"]}), True,
          ["NOTE: no ONT samples in this batch: parameters.longread.ont_adapter is ignored.", done],
          ["rule find_ont_adapter:"])
    # hybrid reads its own block; an unquoted off (YAML False) switches trimming off.
    check("hybrid with ONT adapter trimming off", variant("mode = hybrid",
          [("parameters.hybrid.ont_adapter", False)]), True,
          ["[FunFlux] SAM1: ont reads -> Filtlong guided by the Illumina reads, hifiasm --ont", done],
          ["rule trim_ont_adapters:"],
          expect_jobs={"find_ont_adapter": 0, "trim_ont_adapters": 0, "check_ont_trimming": 0})

    # hybrid: 12 once-per-run jobs + 2 x 50 per-sample jobs + the two Medaka jobs and the
    # three adapter-trimming jobs of SAM1 = 117. Each sample's Illumina reads are cleaned and assembled into a guide (SPAdes).
    # The guide is screened by its own rules (named *_guide_*), and the Illumina pairs
    # that map to the screened guide are selected (rule map_sel_contigs). The example switches
    # the Illumina guidance on (it is off by default), so Filtlong scores SAM1's ONT reads
    # against those pairs (-1 and -2, with --trim --split 1000 from filtlong_split); SAM2's
    # PacBio reads are never guided. The long-read assembly is
    # screened and polished, and Polypolish corrects it with the selected pairs (rule
    # short_read_correction; SAM1's starts from the Medaka consensus, SAM2's from the
    # decontaminated assembly). QUAST compares the screened guide with the final assembly
    # (two labels). The Medaka model is not named in the example, so the check works it
    # out (--model ""). The contig summary judges the long-read draft and gets ITS screen
    # files, never the guide's (which sit in the sample's guide/ folder). As in longread
    # mode, medaka_consensus gets -x right after the draft (-d): it rebuilds the draft's
    # alignment index on every run. As in longread mode, the evidence table is made by
    # rule assembly_evidence, a job of its own per sample.
    check("hybrid example plans the whole chain (117 jobs)", blocks["mode = hybrid"], True,
          ["rule download_phix:", "rule trim_adapters:", "rule illumina_assembly:", "rule filter_contigs:",
           "rule index_guide_contigs:", "rule map_guide_contigs:", "rule map_guide_evaluation:",
           "rule blast_guide_contigs:", "rule blob_guide_json:", "rule blob_guide_table:",
           "rule screen_guide_mitochondria:", "rule select_guide_contigs:", "rule index_selected_contigs:",
           "rule map_sel_contigs:", "rule filter_long_reads:", "rule check_medaka_model:",
           "rule long_read_assembly:", "rule assembly_evidence:", "rule long_read_consensus:",
           "rule short_read_correction:",
           "rule finalize_assembly:", "rule telomere_db:", "rule telomere_units:", "rule contig_summary:",
           "SAM1_selected_R2.fastq --trim --split 1000", '--model ""', "SAM1_medaka_model.txt",
           f"-d {out_hybrid}/02.assembly/SAM1/contaminants/contigs_decontam.fasta -x -t",
           f"cp {out_hybrid}/02.assembly/SAM1/medaka/consensus.fasta",
           f"cp {out_hybrid}/02.assembly/SAM2/contaminants/contigs_decontam.fasta",
           "polypolish polish", "--labels SAM1_guide,SAM1", "--labels SAM2_guide,SAM2",
           f"--contigs {out_hybrid}/02.assembly/SAM1/assembly.fasta",
           f"--decisions {out_hybrid}/02.assembly/SAM1/contaminants/contig_taxonomy_decisions.tsv",
           f"--depth-table {out_hybrid}/02.assembly/SAM1/contaminants/SAM1_depth.tsv "
           f"--assembly-evidence {out_hybrid}/02.assembly/SAM1/assembly_evidence.tsv",
           "-goterms -pa -dp -cpu",
           f"mitochondria.py evidence-hifiasm --gfa {out_hybrid}/02.assembly/SAM1/hifiasm/hifiasm.p_ctg.gfa "
           f"--output {out_hybrid}/02.assembly/SAM1/assembly_evidence.tsv > ", done],
          ["SAM2_selected_R2.fastq --trim", "SAM2_medaka_model.txt", "rule index_contigs:",
           f"--mito-screen {out_hybrid}/02.assembly/SAM1/guide/", f"--depth-table {out_hybrid}/02.assembly/SAM1/guide/",
           "assembly.fasta python"],
          expect_jobs={"total": 117, "download_phix": 1, "build_phix": 1, "map_phix": 2, "trim_adapters": 2,
                       "illumina_assembly": 2, "filter_contigs": 2, "index_guide_contigs": 2,
                       "map_guide_contigs": 2, "map_guide_evaluation": 2, "blast_guide_contigs": 2,
                       "blob_guide_json": 2, "blob_guide_table": 2, "screen_guide_mitochondria": 2,
                       "select_guide_contigs": 2, "index_selected_contigs": 2, "map_sel_contigs": 2,
                       "filter_long_reads": 2, "raw_long_read_qc": 2, "filtered_long_read_qc": 2,
                       "check_medaka_model": 1, "long_read_assembly": 2, "assembly_evidence": 2,
                       "map_contigs": 2, "map_evaluation": 2,
                       "blast_contigs": 2, "screen_mitochondria": 2, "select_contigs": 2,
                       "long_read_consensus": 1, "short_read_correction": 2, "finalize_assembly": 2,
                       "genome_assembly_evaluation": 2, "completeness": 2, "telomere_db": 1, "telomere_units": 2,
                       "contig_summary": 2, "contig_depth": 2, "multiqc": 1, "index_contigs": 0,
                       "chroquetas_db": 0, "amr_chroquetas": 0})
    # Guidance off: Filtlong gets no Illumina pairs for any sample, but the pairs are still
    # selected, because Polypolish needs them. The same 117 jobs.
    check("hybrid with guidance off still selects the Illumina pairs (117 jobs)", variant("mode = hybrid",
          [("parameters.hybrid.short_read_guidance", False)]), True,
          ["[FunFlux] SAM1: ont reads -> adapter check + fastplong, Filtlong not guided, hifiasm --ont", done],
          ["--trim ", "--split"],
          expect_jobs={"total": 117, "map_sel_contigs": 2, "short_read_correction": 2})
    # filtlong_split left empty (null): SAM1's Filtlong is guided by the pairs (-1, -2)
    # but neither trims nor splits the reads, so the next option after the pairs is
    # --min_length.
    check("hybrid with filtlong_split null: guided without --trim and --split (117 jobs)", variant("mode = hybrid",
          [("parameters.hybrid.filtlong_split", None)]), True,
          ["[FunFlux] SAM1: ont reads -> adapter check + fastplong, Filtlong guided by the Illumina reads, hifiasm --ont",
           "SAM1_selected_R2.fastq --min_length 1000", done],
          ["--trim ", "--split"],
          expect_jobs={"total": 117, "filter_long_reads": 2})
    # Medaka off in hybrid mode: 117 - 2 = 115, and Polypolish starts from the
    # decontaminated assembly for both samples.
    check("hybrid with Medaka off (115 jobs)", variant("mode = hybrid", [("parameters.hybrid.medaka_model", False)]), True,
          ["Medaka disabled in the config. Skipping consensus polishing.",
           f"cp {out_hybrid}/02.assembly/SAM1/contaminants/contigs_decontam.fasta", done],
          ["check_medaka_model", "long_read_consensus", "medaka_consensus"],
          expect_jobs={"total": 115, "check_medaka_model": 0, "long_read_consensus": 0, "short_read_correction": 2})
    check("hybrid with a named Medaka model (117 jobs)", variant("mode = hybrid",
          [("parameters.hybrid.medaka_model", "r1041_e82_400bps_sup_v4.3.0")]), True,
          ["Medaka model: r1041_e82_400bps_sup_v4.3.0.", '--model "r1041_e82_400bps_sup_v4.3.0"', done],
          ['--model ""'],
          expect_jobs={"total": 117, "check_medaka_model": 1, "long_read_consensus": 1})
    check("hybrid with Flye (117 jobs)", variant("mode = hybrid", [("parameters.hybrid.assembler", "flye")]), True,
          ["flye --nano-hq", "flye --pacbio-hifi",
           f"mitochondria.py evidence-flye --flye-info {out_hybrid}/02.assembly/SAM2/flye/assembly_info.txt "
           f"--output {out_hybrid}/02.assembly/SAM2/assembly_evidence.tsv > ", done],
          ["hifiasm -i", "evidence-hifiasm", "assembly.fasta python"],
          expect_jobs={"total": 117, "long_read_assembly": 2, "assembly_evidence": 2})

    # The AMR screen switched on in each read-based mode: the download once per run and one
    # ChroQueTas job per sample, so 47 + 2 = 49 (illumina, one sample), 85 + 3 = 88 (longread)
    # and 117 + 3 = 120 (hybrid, two samples each). The screen always runs in these modes and
    # splits the mitochondrion off (mitochondria mode separate_confirmed, the default), so
    # every job hands the mitochondrial run the sample's contigs_mitochondrial.fasta, and
    # MultiQC reads every sample's AMR table. The report folder is 06.report.
    check("illumina with the AMR screen on (49 jobs)", variant("mode = illumina", amr_on), True,
          ["rule chroquetas_db:", "rule amr_chroquetas:",
           f"--genome {out_illumina}/02.assembly/SAM1/contigs_final.fasta "
           f"--mito-contigs {out_illumina}/02.assembly/SAM1/contigs_mitochondrial.fasta",
           "05.amr/SAM1/amr_summary_mqc.tsv", f'--outdir "{out_illumina}/06.report"', done],
          ["05.report"],
          expect_jobs={"total": 49, "chroquetas_db": 1, "amr_chroquetas": 1, "multiqc": 1})
    check("longread with the AMR screen on (88 jobs)", variant("mode = longread", amr_on), True,
          [f"--mito-contigs {out_longread}/02.assembly/SAM1/contigs_mitochondrial.fasta",
           f"--mito-contigs {out_longread}/02.assembly/SAM2/contigs_mitochondrial.fasta",
           "05.amr/SAM1/amr_summary_mqc.tsv 05.amr/SAM2/amr_summary_mqc.tsv",
           f'--outdir "{out_longread}/06.report"', done],
          ["05.report"],
          expect_jobs={"total": 88, "chroquetas_db": 1, "amr_chroquetas": 2, "multiqc": 1})
    check("hybrid with the AMR screen on (120 jobs)", variant("mode = hybrid", amr_on), True,
          [f"--mito-contigs {out_hybrid}/02.assembly/SAM1/contigs_mitochondrial.fasta",
           f"--mito-contigs {out_hybrid}/02.assembly/SAM2/contigs_mitochondrial.fasta",
           "05.amr/SAM1/amr_summary_mqc.tsv 05.amr/SAM2/amr_summary_mqc.tsv",
           f'--outdir "{out_hybrid}/06.report"', done],
          ["05.report"],
          expect_jobs={"total": 120, "chroquetas_db": 1, "amr_chroquetas": 2, "multiqc": 1})

    print("=== download_all plans the download steps, and nothing else ===")
    # Rule download_all (workflow/Snakefile) asks for the files of every download step
    # the configuration defines (download_targets() in 00_common.smk), so that a run can be
    # prepared on a computer with internet. Named on the command line, it must plan those
    # steps, the path repair of the funannotate database (normalize_funannotate_db, which
    # follows the download) and itself, and no other job: no assembly, no annotation, no
    # rule all. Every case gives the total, so one extra job of any rule makes it fail. The
    # stand-in configurations name a funannotate folder that holds no database, so the
    # database download is planned in each of them.
    #
    # Every mode plans the funannotate database and its repair, the telomere database, the
    # UniVec download of the adapter screen and the Tiberius model, and, as the examples hold
    # no local copies, the UNITE and antiSMASH downloads: with download_all itself, 8 jobs. contigs mode adds nothing else.
    #
    # With parameters.busco.lineage: auto, as in the examples, the job of download_all
    # is shown by its message: the note that every BUSCO job will need internet, printed
    # in place of the usual "rule download_all:" description. With a named lineage there
    # is no note, and the usual description is printed.
    download_everywhere = {"download_all": 1, "funannotate_db": 1, "normalize_funannotate_db": 1,
                           "telomere_db": 1, "univec_illumina": 1, "tiberius_model": 1, "unite_db": 1,
                           "antismash_db": 1,
                           "all": 0}
    busco_auto_note = "[FunFlux] Note: parameters.busco.lineage is auto, so every BUSCO job downloads its own data"
    check("download_all in the contigs example plans the downloads and nothing else (8 jobs)",
          contigs_plain, True, [busco_auto_note, done],
          ["rule filter_contigs:", "rule unite_identified:", "rule completeness:"],
          expect_jobs=dict(download_everywhere, total=8, download_phix=0, chroquetas_db=0, busco_data=0,
                           check_medaka_model=0),
          target="download_all")
    # The AMR screen on adds the ChroQueTas download (8 + 1 = 9); off, as above, it is not
    # planned.
    check("download_all with the AMR screen on adds chroquetas_db (9 jobs)", variant("mode = contigs", amr_on), True,
          ["rule chroquetas_db:", done], ["rule amr_chroquetas:"],
          expect_jobs=dict(download_everywhere, total=9, chroquetas_db=1),
          target="download_all")
    # A named BUSCO lineage adds its download (8 + 1 = 9); BUSCO itself is not planned,
    # and there is no BUSCO note.
    check("download_all with a named BUSCO lineage adds busco_data (9 jobs)", variant("mode = contigs",
          [("parameters.busco.lineage", "hypocreales_odb12")]), True,
          [f"busco --download hypocreales_odb12 --download_path {out_contigs}/00.resources/busco",
           "rule download_all:", done],
          ["rule completeness:", busco_auto_note],
          expect_jobs=dict(download_everywhere, total=9, busco_data=1),
          target="download_all")
    # Local antiSMASH and UNITE copies: nothing to download for either, and their link
    # rules (antismash_db_local, unite_db_local) are not downloads, so 8 - 2 = 6.
    check("download_all with local antiSMASH and UNITE copies downloads neither (6 jobs)", variant("mode = contigs",
          [("directories.antismash_db", base + "/ref/antismash_db"), ("directories.unite_db", base + "/ref/unite_db")]),
          True, [done], ["rule antismash_db:", "rule unite_db:", "rule antismash_db_local:", "rule unite_db_local:"],
          expect_jobs={"total": 6, "download_all": 1, "funannotate_db": 1, "normalize_funannotate_db": 1,
                       "telomere_db": 1, "tiberius_model": 1},
          target="download_all")
    check("download_all with Tiberius off: no model download (7 jobs)", variant("mode = contigs",
          [("parameters.tiberius", {"mode": "off"})]), True,
          [done], ["rule tiberius_model:"],
          expect_jobs={"total": 7, "tiberius_model": 0, "download_all": 1, "funannotate_db": 1},
          target="download_all")
    # illumina mode adds the PhiX genome (8 + 1 = 9), but not its bowtie2 index (a step
    # that needs no internet).
    check("download_all in the illumina example adds download_phix (9 jobs)", blocks["mode = illumina"], True,
          ["rule download_phix:", done], ["rule build_phix:", "rule map_phix:"],
          expect_jobs=dict(download_everywhere, total=9, download_phix=1, build_phix=0, check_medaka_model=0),
          target="download_all")
    # longread mode adds the Medaka model check of the ONT sample SAM1 (8 + 1 = 9), and
    # none for the PacBio sample SAM2. The assembler that waits for that check is not planned.
    check("download_all in the longread example adds the Medaka check of the ONT sample (9 jobs)",
          blocks["mode = longread"], True,
          ["rule check_medaka_model:", "SAM1_medaka_model.txt", done],
          ["SAM2_medaka_model.txt", "rule long_read_assembly:", "rule download_phix:"],
          expect_jobs=dict(download_everywhere, total=9, check_medaka_model=1, download_phix=0),
          target="download_all")
    check("download_all in longread mode with Medaka off: no model check (8 jobs)", variant("mode = longread",
          [("parameters.longread.medaka_model", False)]), True,
          [done], ["rule check_medaka_model:"],
          expect_jobs=dict(download_everywhere, total=8, check_medaka_model=0),
          target="download_all")
    # hybrid mode has both: the PhiX genome and the Medaka check of SAM1 (8 + 2 = 10); with
    # the AMR screen on, ChroQueTas too (10 + 1 = 11).
    check("download_all in the hybrid example: PhiX and the Medaka check (10 jobs)", blocks["mode = hybrid"], True,
          ["rule download_phix:", "rule check_medaka_model:", "SAM1_medaka_model.txt", done],
          ["SAM2_medaka_model.txt", "rule illumina_assembly:"],
          expect_jobs=dict(download_everywhere, total=10, download_phix=1, check_medaka_model=1, chroquetas_db=0),
          target="download_all")
    check("download_all in hybrid mode with the AMR screen on (11 jobs)", variant("mode = hybrid", amr_on), True,
          ["rule chroquetas_db:", done], ["rule amr_chroquetas:"],
          expect_jobs=dict(download_everywhere, total=11, download_phix=1, check_medaka_model=1, chroquetas_db=1),
          target="download_all")
    # The other way round: a normal run (no rule named) plans rule all and never
    # download_all. This holds only while rule all is the first rule written in the
    # Snakefile, so the case guards that order.
    check("a normal run plans rule all, not download_all (31 jobs)", contigs_plain, True, [done],
          ["rule download_all:", busco_auto_note],
          expect_jobs={"total": 31, "all": 1, "download_all": 0})

    print("=== the memory and time the jobs ask a cluster for ===")
    # Under a cluster executor Snakemake hands each job's mem_mb (MB) and runtime (minutes)
    # resources to the scheduler. The heavy steps ask for resources.ram_gb (32 GB, so 32768
    # MB, in the illumina example) and resources.runtime_hours (48 h by default, so 2880
    # minutes); SPAdes (rule illumina_assembly) keeps its own mem_gb, which its -m option
    # reads. Qualimap (rule map_evaluation) asks for its Java heap plus 4 GB, but never more
    # than ram_gb: (32 + 4) x 1024 = 36864 MB would be more than the 32 GB set, so it asks
    # for 32768 MB. map_phix reads every read but needs little memory, so it asks only for
    # the time; fastp (rule trim_adapters) is light and asks for nothing. The cases compare
    # the values in each job's resources line (expect_resources), so they hold whatever
    # order and extra entries a Snakemake version prints there.
    check("heavy steps ask for ram_gb and 48 hours by default", blocks["mode = illumina"], True, [done],
          expect_resources={
              "illumina_assembly": {"mem_gb": 32, "mem_mb": 32768, "runtime": 2880},
              "map_evaluation": {"java_mem": 32, "mem_mb": 32768, "runtime": 2880},
              "map_phix": {"runtime": 2880, "mem_mb": None},
              "trim_adapters": {"mem_mb": None, "runtime": None},
          })
    # With ram_gb 100 the Qualimap heap stays at its cap of 64 GB, and the job asks for the
    # heap plus 4 GB, (64 + 4) x 1024 = 69632 MB, which is less than ram_gb, while the heavy
    # steps ask for 100 x 1024 = 102400 MB. No job may still ask for 36864 or 32768 MB.
    check("Qualimap asks for its heap plus 4 GB when ram_gb is larger", variant("mode = illumina",
          [("resources.ram_gb", 100)]), True, [done], ["mem_mb=36864", "mem_mb=32768"],
          expect_resources={
              "map_evaluation": {"java_mem": 64, "mem_mb": 69632},
              "illumina_assembly": {"mem_gb": 100, "mem_mb": 102400},
          })
    # A good value: 12 hours become 720 minutes, on every long step; no job may still ask
    # for 2880 minutes.
    check("runtime_hours 12 gives the long steps 720 minutes", variant("mode = illumina",
          [("resources.runtime_hours", 12)]), True, [done], ["runtime=2880"],
          expect_resources={
              "illumina_assembly": {"runtime": 720},
              "map_phix": {"runtime": 720},
          })

    print("=== accepted, with a NOTE ===")
    # Unusual settings that must be accepted. Some must print a NOTE; the others check
    # that the setting is read the way the user meant it.
    check("blast_max_hsps below 5 warns", variant("mode = illumina", [("parameters.blast_max_hsps", 1)]), True,
          ["NOTE: parameters.blast_max_hsps is 1"])
    # min_flye_coverage_ratio, the other name of the key, still works, and its value (4) is
    # the one used.
    check("other name of the coverage-ratio key still works", variant("mode = illumina",
          [("parameters.mitochondria", {"min_flye_coverage_ratio": 4})]), True,
          ["NOTE: parameters.mitochondria.min_flye_coverage_ratio is read as min_assembler_coverage_ratio",
           "assembler coverage ratio >= 4"])
    # A split length other than the default, written while the guidance is off. The split
    # cuts long reads where the Illumina reads stop supporting them, so without guidance it
    # is not used; a NOTE says so, because someone chose that number. The ONT sample's
    # Filtlong is not guided.
    check("hybrid split of 500 is ignored with a NOTE when guidance is off", variant("mode = hybrid",
          [("parameters.hybrid.short_read_guidance", False), ("parameters.hybrid.filtlong_split", 500)]), True,
          ["NOTE: parameters.hybrid.filtlong_split is ignored",
           "[FunFlux] SAM1: ont reads -> adapter check + fastplong, Filtlong not guided, hifiasm --ont"])
    # The configuration template writes the default split, 1000, with the guidance off. It is
    # not used either, but gets no NOTE: every default hybrid run would print one.
    check("hybrid split of 1000 with guidance off gives no NOTE", variant("mode = hybrid",
          [("parameters.hybrid.short_read_guidance", False)]), True,
          ["[FunFlux] SAM1: ont reads -> adapter check + fastplong, Filtlong not guided, hifiasm --ont"],
          ["filtlong_split is ignored"])
    # Without the short_read_guidance key the guidance is off, the default: no -1, -2, --trim
    # or --split for the ONT sample.
    check("hybrid without short_read_guidance: Filtlong not guided", variant("mode = hybrid",
          [("parameters.hybrid.short_read_guidance", DELETE)]), True,
          ["[FunFlux] SAM1: ont reads -> adapter check + fastplong, Filtlong not guided, hifiasm --ont"],
          ["--trim ", "--split", "SAM1_selected_R2.fastq --min_length"])
    # A sample name with an underscore is not a valid NCBI locus-tag prefix (3 to 12 letters
    # or digits, starting with a letter); such a sample still runs.
    check("sample name that is not a locus-tag prefix is accepted with a NOTE", variant("mode = illumina",
          [("files.annotation_params", make_params([("S_1", base + "/ref/proteins.faa")]))],
          {"illumina_dir": ["S_1_R1.fastq.gz", "S_1_R2.fastq.gz"]}), True,
          ["NOTE: sample name 'S_1' is not a valid NCBI locus-tag prefix", "Sample S_1 will be processed."])
    # YAML reads an unquoted `off` as the boolean False; the mode readers must treat it as the word "off".
    check("an unquoted mitochondria mode off (YAML False) is read as off", variant("mode = illumina",
          [("parameters.mitochondria", {"mode": False})]), True,
          ["Mitochondrion screen: off"])
    check("an unquoted decontamination mode off (YAML False) is read as off", variant("mode = illumina",
          [("parameters.decontamination", {"mode": False})]), True,
          ["Decontamination mode: off"])
    check("Medaka switched off", variant("mode = longread", [("parameters.longread.medaka_model", False)]), True,
          ["[FunFlux] SAM1: ont reads -> adapter check + fastplong, hifiasm --ont, minimap2 map-ont, no Medaka."])
    # A Medaka model with "fast" in its name means older, noisier ONT reads, so Flye gets
    # --nano-raw for the ONT sample; the PacBio sample always gets --pacbio-hifi.
    check("Flye with an explicit fast model gets --nano-raw", variant("mode = longread",
          [("parameters.longread.assembler", "flye"), ("parameters.longread.medaka_model", "r941_min_fast_g507")]), True,
          ["Flye --nano-raw", "Flye --pacbio-hifi"])
    # The hybrid example names no Medaka model, so the ONT sample gets --nano-hq.
    check("hybrid with Flye and the default model gets --nano-hq", variant("mode = hybrid",
          [("parameters.hybrid.assembler", "flye")]), True, ["Flye --nano-hq", "Flye --pacbio-hifi"])
    # A v1 configuration with an empty parameters.genus: FunFlux 1.x then kept the genus with
    # the most contigs; FunFlux v2 keeps every contig (mode off) and says so in a NOTE.
    check("an empty v1 genus becomes decontamination mode off, with a NOTE", variant("mode = illumina",
          [("parameters.genus", None)]), True,
          ["NOTE: parameters.genus (a v1 key) is empty", "Decontamination mode: off"])
    # include mode is accepted when every sample gets a genus from somewhere: SAM1 from the
    # include_genera_by_sample table, SAM2 from its row in the sample_overrides table.
    check("include mode with each sample's genus from a table is accepted", variant("mode = longread",
          [("parameters.decontamination", {"mode": "include",
                                           "include_genera_by_sample": base + "/ref/include_by_sample.tsv",
                                           "sample_overrides": base + "/ref/overrides_SAM2_genus.tsv"})]), True,
          ["Decontamination mode: include", done])
    # Without the screen there is no rule select_contigs, so the decontamination tables
    # are not checked: contigs mode without the screen accepts a table path that does
    # not exist (no_overrides.tsv is never created; the same path is refused above,
    # where the screen runs).
    check("decontamination tables are not checked without the screen", variant("mode = contigs",
          [("parameters.decontamination", {"sample_overrides": base + "/ref/no_overrides.tsv"})]), True,
          ["Decontamination mode: off", done], ["is not an existing file"])
    # An empty GeneMark mode (YAML null, or "") means auto, like the other mode keys: with
    # directories.genemark_dir set, GeneMark is switched on. These two cases run with the
    # fake home folder, so the stand-in licence key is found there.
    check("GeneMark mode left empty means auto", variant("mode = illumina",
          [("directories.genemark_dir", base + "/ref/genemark"), ("parameters.genemark", {"mode": None})]), True,
          ["NOTE: GeneMark-ES is enabled"], home=base + "/home")
    check("GeneMark mode \"\" means auto", variant("mode = illumina",
          [("directories.genemark_dir", base + "/ref/genemark"), ("parameters.genemark", {"mode": ""})]), True,
          ["NOTE: GeneMark-ES is enabled"], home=base + "/home")
    # The word null in quotes means off, like none (an unquoted null is an empty value,
    # which means: work the model out from the reads).
    check("Medaka model 'null' in quotes means off", variant("mode = longread",
          [("parameters.longread.medaka_model", "null")]), True,
          ["[FunFlux] SAM1: ont reads -> adapter check + fastplong, hifiasm --ont, minimap2 map-ont, no Medaka.",
           "Medaka disabled in the config"], ["Medaka model: null"])
    # Column 4 of the annotation table may be left out (three columns) or written "-": then
    # funannotate trains AUGUSTUS on the genome. A named model that the database holds is
    # used as it is. Each run says which, in one line.
    check("annotation table without a fourth column trains AUGUSTUS", variant("mode = illumina",
          [("files.annotation_params", make_params([("SAM1", base + "/ref/proteins.faa")], None))]), True,
          ["AUGUSTUS: SAM1 is trained on its BUSCO genes (checked once the database is unpacked).", done])
    check("annotation table with '-' as model trains AUGUSTUS", variant("mode = illumina",
          [("files.annotation_params", make_params([("SAM1", base + "/ref/proteins.faa")], "-"))]), True,
          ["AUGUSTUS: SAM1 is trained on its BUSCO genes (checked once the database is unpacked).", done])
    check("AUGUSTUS model the database holds is accepted", variant("mode = illumina",
          [("directories.funannotate_db", base + "/ref/fundb_species")]), True,
          ["AUGUSTUS: SAM1 uses the stored model fusarium_graminearum.", done])
    # Before the database is downloaded (the example blocks name a folder that does not
    # exist), a named model can only be announced, not confirmed.
    check("AUGUSTUS model named before the database is on disk", blocks["mode = illumina"], True,
          ["AUGUSTUS: SAM1 names the stored model fusarium_graminearum (checked once the database is unpacked).",
           done])
    check("annotation table row whose fourth column is an empty field", variant("mode = illumina",
          [("files.annotation_params", base + "/ref/params_trailing_tab.tsv")]), True,
          ["AUGUSTUS: SAM1 is trained on its BUSCO genes (checked once the database is unpacked).", done])
    check("mixed batch: one stored model, one trained", variant("mode = longread",
          [("directories.funannotate_db", base + "/ref/fundb_species"),
           ("files.annotation_params", base + "/ref/params_mixed.tsv")]), True,
          ["AUGUSTUS: SAM1 uses the stored model fusarium_graminearum; SAM2 is trained on its BUSCO genes.", done])
    check("a model stored under the sample's own name is used, with a NOTE", variant("mode = illumina",
          [("directories.funannotate_db", base + "/ref/fundb_own"),
           ("files.annotation_params", make_params([("SAM1", base + "/ref/proteins.faa")], None))]), True,
          ["AUGUSTUS: SAM1 uses the model fusarium_graminearum_sam1, stored earlier under its own name.",
           "NOTE: sample 'SAM1' leaves the AUGUSTUS model empty", done])
    check("AUGUSTUS optimisation for a sample without a model", variant("mode = illumina",
          [("parameters.prediction", {"augustus_optimize": True}),
           ("files.annotation_params", make_params([("SAM1", base + "/ref/proteins.faa")], None))]), True,
          ["AUGUSTUS: SAM1 is trained on its BUSCO genes, with the optimisation (checked once the database is unpacked).", done],
          ["NOTE: parameters.prediction.augustus_optimize"])
    check("AUGUSTUS optimisation with a named model gives a NOTE", variant("mode = illumina",
          [("parameters.prediction", {"augustus_optimize": True})]), True,
          ["NOTE: parameters.prediction.augustus_optimize is on, but it applies only to samples whose AUGUSTUS is trained",
           "SAM1 name(s) a stored model"])
    # The illumina example has resources.threads: 8, so 16 SignalP threads are lowered to 8.
    check("SignalP 6 threads above resources.threads are capped, with a NOTE", variant("mode = illumina",
          [("parameters.signalp6", {"mode": "fast", "torch_num_threads": 16}),
           ("files.signalp6_archive", base + "/ref/signalp6.tar.gz")]), True,
          ["NOTE: parameters.signalp6.torch_num_threads is 16, more than resources.threads (8); "
           "SignalP 6 will use 8 threads."])
    # YAML reads an unquoted `on` as the boolean True; Tiberius's mode reader must take it as on.
    # With Tiberius off the model link is not read, so it may be missing. A checksum in capitals
    # is accepted (it is lowered).
    check("Tiberius mode written as a YAML true is read as on", variant("mode = illumina",
          [("parameters.tiberius", {"mode": True})]), True, ["Tiberius: on", done])
    check("Tiberius off needs no model link", variant("mode = illumina",
          [("parameters.tiberius", {"mode": False}), ("links.tiberius_model_link", DELETE)]), True,
          ["Tiberius: off.", done])
    check("Tiberius checksum in capitals is accepted", variant("mode = illumina",
          [("parameters.tiberius", {"model_sha256": "4A309A630574F395D82D747DCF64716704DAB056DB534E427AE0CEC1B7410E20"})]), True,
          ["Tiberius: on", done])
    # The contig summary's settings (parameters.contig_summary). The defaults are printed in
    # one line; the motif, the window and the copies reach telomere_units.py as its options,
    # the copies and the rDNA units contig_summary.py. The motif auto is accepted in any case
    # and passed on as auto; a unit written in lower case is upper-cased; the floors of the
    # three numbers (window 500, copies 2, rDNA units 1) are accepted values; the database
    # checksum in capitals is accepted (it is lower-cased, and rule telomere_db checks the file
    # against the lower-cased form). Each command is printed on one line, so its options stand
    # next to each other.
    check("contig summary defaults are printed and passed on", contigs_plain, True,
          ["Contig summary: telomere motif auto (the repeat unit is chosen from the data), a telomere from 5 tandem "
           "copies within 25000 bp of a contig end, an rDNA array from 3 rDNA unit(s).",
           "--telomere-motif auto --telomere-window 25000 --min-telomere-copies 5",
           "--min-telomere-copies 5 --min-rdna-units 3", done])
    check("telomere motif auto written out is accepted", variant("mode = contigs",
          [("parameters.contig_summary", {"telomere_motif": "auto"})]), True,
          ["Contig summary: telomere motif auto (the repeat unit is chosen from the data),", "--telomere-motif auto", done])
    check("telomere motif Auto with a capital is read as auto", variant("mode = contigs",
          [("parameters.contig_summary", {"telomere_motif": "Auto"})]), True,
          ["Contig summary: telomere motif auto (the repeat unit is chosen from the data),", "--telomere-motif auto", done],
          ["--telomere-motif Auto"])
    check("lower-case telomere motif is upper-cased and used as given", variant("mode = contigs",
          [("parameters.contig_summary", {"telomere_motif": "ttaggg"})]), True,
          ["Contig summary: telomere motif TTAGGG (used as given; the data's own choice is reported next to it),",
           "--telomere-motif TTAGGG", done])
    check("contig summary settings at their floors (min_telomere_copies 2) are accepted", variant("mode = contigs",
          [("parameters.contig_summary", {"telomere_window": 500, "min_telomere_copies": 2, "min_rdna_units": 1})]), True,
          ["a telomere from 2 tandem copies within 500 bp of a contig end, an rDNA array from 1 rDNA unit(s).",
           "--telomere-window 500 --min-telomere-copies 2", "--min-telomere-copies 2 --min-rdna-units 1", done])
    # A longer unit is passed on as it is: the 12-letter unit of the Aspergillus flavus group.
    check("a 12-letter telomere motif is accepted", variant("mode = contigs",
          [("parameters.contig_summary", {"telomere_motif": "TTAGGGTCAACA"})]), True,
          ["Contig summary: telomere motif TTAGGGTCAACA (used as given", "--telomere-motif TTAGGGTCAACA", done])
    check("telomere database checksum in capitals is accepted", variant("mode = contigs",
          [("parameters.contig_summary",
            {"telomere_db_sha256": "1D4F5C8C2CCE6A76AAA8DE01602317A8E0B8D0401711B96F3C7A826C385DA699"})]), True,
          ["Contig summary: telomere motif auto",
           f"1d4f5c8c2cce6a76aaa8de01602317a8e0b8d0401711b96f3c7a826c385da699 {out_contigs}/00.resources/telomeres/telo_database_v2.tsv",
           done],
          ["1D4F5C8C"])
    # The AMR screen's settings (parameters.amr). An unquoted on (YAML True) switches it on, and
    # an unquoted off (YAML False) off, which needs no download address. The mitochondrial
    # code reaches the wrapper as auto or as the number: auto in any case and an empty value
    # both mean auto, a number is passed on as it is. The two checksums in capitals are
    # accepted and lower-cased (rule chroquetas_db checks the files against the lower-cased
    # form). A slash at the end of the download address is removed, so that the rule's three
    # addresses have no double slash.
    check("AMR switch written as a YAML true is read as on", variant("mode = contigs",
          [("parameters.amr", {"chroquetas": True}), ("links.chroquetas_link", chroquetas_link)]), True,
          ["AMR screen: ChroQueTas on (licence CC BY-NC-ND 4.0: non-commercial use).", "rule amr_chroquetas:", done])
    check("AMR switch written as a YAML false is read as off and needs no link", variant("mode = contigs",
          [("parameters.amr", {"chroquetas": False})]), True,
          ["AMR screen: off.", done], ["rule chroquetas_db:", "rule amr_chroquetas:"])
    check("mitochondrial genetic code 4 is passed on", variant("mode = contigs",
          [("parameters.amr", {"chroquetas": "on", "mito_genetic_code": 4}), ("links.chroquetas_link", chroquetas_link)]), True,
          ["--mito-genetic-code 4 ", done])
    check("mitochondrial genetic code Auto with a capital is read as auto", variant("mode = contigs",
          [("parameters.amr", {"chroquetas": "on", "mito_genetic_code": "Auto"}), ("links.chroquetas_link", chroquetas_link)]), True,
          ["--mito-genetic-code auto ", done], ["--mito-genetic-code Auto"])
    check("mitochondrial genetic code left empty means auto", variant("mode = contigs",
          [("parameters.amr", {"chroquetas": "on", "mito_genetic_code": None}), ("links.chroquetas_link", chroquetas_link)]), True,
          ["--mito-genetic-code auto ", done])
    check("AMR checksums in capitals are accepted", variant("mode = contigs",
          [("parameters.amr", {"chroquetas": "on",
                               "chroquetas_script_sha256": "D8961A0D311BE57D4A9939BD9E4D90DBC796589B3298DD1835359ED776ABFCE3",
                               "chroquetas_db_sha256": "6F89ECB3117110AAD85B8B6CA7B8A3D3715CB5984B6CB75450A6CDA5633D16E7"}),
           ("links.chroquetas_link", chroquetas_link)]), True,
          [f"d8961a0d311be57d4a9939bd9e4d90dbc796589b3298dd1835359ed776abfce3 {out_contigs}/00.resources/chroquetas/ChroQueTas.sh",
           f"6f89ecb3117110aad85b8b6ca7b8a3d3715cb5984b6cb75450a6cda5633d16e7 {out_contigs}/00.resources/chroquetas/db-v20250811.tgz",
           done],
          ["D8961A0D", "6F89ECB3"])
    check("a slash at the end of the ChroQueTas link is removed", variant("mode = contigs",
          [("parameters.amr", {"chroquetas": "on"}), ("links.chroquetas_link", chroquetas_link + "/")]), True,
          [f'"{chroquetas_link}/bin/ChroQueTas.sh"', done], ["//bin/ChroQueTas.sh"])

    # Summary line, then remove the scratch folder unless --keep was given.
    print()
    if failures:
        print(f"{len(failures)} case(s) did not behave as expected: " + "; ".join(failures))
    else:
        print("All cases behaved as expected.")
    if args.keep:
        print(f"Scratch directory kept: {base}")
    else:
        shutil.rmtree(base, ignore_errors=True)
    return 1 if failures else 0


# Run main() when the file is started as a script, and hand its result to the shell as
# the exit status.
if __name__ == "__main__":
    sys.exit(main())
