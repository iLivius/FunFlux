#!/usr/bin/env python3
"""Check the configured Medaka model, or work one out from the reads, before
the assembler starts. When the configured model is wrong, print the closest
available ones so that the fix is obvious.

Why check so early
------------------
Medaka polishes a long-read (ONT) assembly with the reads, using a model that
must match how the reads were basecalled. Polishing is one of the LAST steps of
a long-read run, because it needs the finished assembly. A model name this
Medaka version does not have (a typo, or a model dropped in a newer release)
therefore used to show up only after the assembler had already run for an hour
or more. Run before the assembler, this script stops such a run in seconds.

Where it fits in the workflow
-----------------------------
Run by rule check_medaka_model (rules/shared/12_medaka_check.smk), once per ONT
sample that uses Medaka, inside the Medaka conda environment, so `medaka` is on
PATH. That rule file, and the polishing rules that read the result, come with
the longread and hybrid modes and are not in the repository yet (the design is
in docs/methods/medaka-model.md). 00_common.smk already defines this
script's path (MEDAKA_CHECK_SCRIPT), its output file (MEDAKA_MODEL_RESOLVED) and
medaka_gate(), which the long-read assembler rule will list among its inputs so
that the assembler waits for that output.
CHANGED in FunFlux v2: BacFlux runs the check right after read filtering, on the
filtered reads. FunFlux runs it on the raw reads, so that it does not have to
wait for filtering (see --reads below).

Command-line options, all passed in by that rule
------------------------------------------------
  --model    the model name the user set in
             parameters.{longread|hybrid}.medaka_model. Empty means auto mode:
             work the model out from the reads.
  --reads    the RAW ONT read file the user gave: FASTQ, plain or
             gzip-compressed (.fastq or .fastq.gz). Auto mode reads the
             basecaller tag from the FASTQ headers of its first 100 reads.
             CHANGED in FunFlux v2: BacFlux passes FILT_LONG, the Filtlong
             output. FunFlux passes the raw file instead, because in hybrid mode
             the filtered file exists only hours into the run, and this check
             has to answer before the assembler starts. The file goes to Medaka
             as it is, .gz included (see resolve_from_reads()).
  --fallback-to-auto
             the opt-in middle way (config key medaka_model_fallback_auto,
             default off): when a NAMED model is not available, try auto mode
             instead of stopping, and use its answer with a loud warning.
  --out      the checked or worked-out model, written as one line. The
             polishing rule long_read_consensus (rules/longread/30_polish.smk,
             and its twin in rules/hybrid/50_polish.smk) reads it back with cat
             and passes it to `medaka_consensus -m`, so the model is decided
             once, not twice, and the two steps cannot disagree.

Exit status 0, with --out written, on success. On failure the script prints an
error and a table of suggested models to stderr, writes no --out, and exits
with status 1. A broken Medaka (no model list) also ends with status 1, with a
Python traceback instead of a table; a missing or unknown option ends with
status 2 (argparse, Python's module for reading options).

Auto mode
---------
Auto mode asks Medaka for the generic CONSENSUS model that matches the
basecaller named in the read headers (`--auto_model consensus` below): the
standard polishing model for that flowcell chemistry, accuracy level and
basecaller version.

CHANGED in FunFlux v2: BacFlux asks for `--auto_model consensus_bacteria`
instead, the research model Medaka offers for native BACTERIAL DNA (it corrects
the systematic errors that bacterial DNA methylation leaves in standard
basecalls). It is not meant for fungal DNA, so FunFlux asks for the generic
consensus model.

Auto mode also needs two things (CHANGED in FunFlux v2: this note is new):
  * Network access, once. Working out a model makes Medaka download its weights
    (from GitHub) the first time it meets a model that is not installed yet;
    later runs reuse the stored copy.
  * Exactly ONE basecaller model among the first 100 reads, read from the
    FASTQ header tag (basecall_model_version_id=...). Medaka stops when it
    finds none, or two different ones. ONT reads from public archives, or
    pooled from several basecalling runs, usually lack or mix tags: then the
    config must name the model (parameters.{longread|hybrid}.medaka_model). To
    see what your own reads carry: zcat reads.fastq.gz | head -n 1

Downloading a named model here
------------------------------
Medaka 2.2.2 ships the weights of four models only (all r1041_e82_400bps,
two of them for variant calling). For any other model, medaka_consensus
downloads the weights from GitHub the first time it is asked for them. In
FunFlux that would be rule long_read_consensus, the polishing step, which on a
cluster runs on a compute node, and compute nodes often have no internet. So
once a named model has passed the check, this script makes Medaka store its
weights right away (`medaka tools download_models --models <name>`), and rule
check_medaka_model runs on the computer that runs Snakemake itself, which has
internet (it is a local rule). Medaka keeps the weights inside its own conda
environment (or, when that folder cannot be written, in ~/.medaka/data), so
the polishing job finds them there without going to the internet. A model that
is already stored, the four shipped ones included, needs no network: Medaka
then only reports where the file is. When the weights cannot be fetched, the
script stops with a message, before the assembler starts, instead of letting
the polishing fail hours later. (CHANGED in FunFlux v2: new.) Auto mode needs
no extra step: working out the model already stores its weights (see above).

The suggestion tables
---------------------
The tables printed on failure list the consensus models of this Medaka version:
everything `medaka tools list_models` reports except the variant-calling models
(variant or snp in the name). That is a plain filter on the names, so the list
can also hold special-purpose models, for example the bacterial-methylation and
plant models; ignore those for fungi.

Testing
-------
Only list_available_models(), resolve_from_reads() and download_model() call
Medaka. Splitting a
model name into its parts, building the suggestion tables and main()'s choice
between a named model and auto mode do not. So the tests
(test_medaka_model_check.py, in the same folder) run with no Medaka installed:
they put stand-ins in place of the Medaka calls.
"""

import argparse
import os
import re
import subprocess
import sys


# ── Medaka calls (the only functions that touch the tool) ────────────────────

def list_available_models():
    """Return the model names this Medaka version knows, read from
    `medaka tools list_models`. That command prints one line
    'Available: m1, m2, ...' plus two 'Default ...:' lines, which are ignored.
    Asking the installed tool instead of carrying a hard-coded list is the whole
    point: the valid names depend on the Medaka version, which is exactly how a
    config model that worked last year stops working.

    Called once by main(), in every mode: the list is needed to check a named
    model and to print the suggestion tables. If the command prints no
    'Available:' line, Medaka itself is broken: the RuntimeError below is not
    caught, so Python prints a traceback (its full error report) and exits with
    status 1."""
    proc = subprocess.run(
        ["medaka", "tools", "list_models"],
        capture_output=True, text=True,
    )
    for line in proc.stdout.splitlines():
        if line.startswith("Available:"):
            # The rest of the line is a comma-separated list. Strip the spaces
            # around each name and drop empty pieces.
            names = line[len("Available:"):].split(",")
            return [n.strip() for n in names if n.strip()]
    # No "Available:" line means the tool itself failed (env broken, etc.).
    raise RuntimeError(
        "Could not read the Medaka model list ('medaka tools list_models' gave no "
        "'Available:' line). stderr:\n" + proc.stderr
    )


def resolve_from_reads(reads):
    """Ask Medaka to work out the consensus model from the basecaller tag in
    the reads (`medaka tools resolve_model --auto_model consensus <reads>`).

    Returns what Medaka prints: a PATH to the model file rather than a bare
    name, which is a valid `-m` argument all the same. Returns None when Medaka
    cannot work one out (no tag, two different tags, or a tag its lookup does
    not know): for this script that is an answer, not an error. Medaka's own
    messages are not shown, so a failed download of the model weights (no
    network) also ends as None. Called by main() in auto mode, and for a bad
    named model when --fallback-to-auto is on."""
    # CHANGED in FunFlux v2: the model type is "consensus" (Medaka's generic
    # consensus models), not "consensus_bacteria" as in BacFlux. Medaka 2.2.2
    # accepts "consensus", "variant" and "consensus_bacteria" here (class
    # AutoModel in medaka/medaka.py), so nothing else in the call changes.
    #
    # CHANGED in FunFlux v2: `reads` may be gzip-compressed, and is handed to
    # Medaka exactly as it is. In FunFlux it is the RAW file the user gave,
    # often reads.fastq.gz, and no plain-text copy is made because Medaka reads
    # gzip itself: models._model_from_fastq() opens the FASTQ with
    # pysam.FastxFile, which decompresses gzip on the fly, and it stops after
    # the first 100 records, so the size of the file does not matter. Checked
    # on 2026-10-07 against Medaka 2.2.2 (pysam 0.24.0):
    # `medaka tools resolve_model --auto_model consensus` returned the same
    # model for a plain FASTQ, a .fastq.gz, and two .fastq.gz files joined with
    # cat, and took about 3 seconds (mostly Medaka starting up) both for a
    # 5-read file and for a synthetic .fastq.gz holding 1 GB of FASTQ.
    proc = subprocess.run(
        ["medaka", "tools", "resolve_model", "--auto_model", "consensus", reads],
        capture_output=True, text=True,
    )
    # A non-zero exit status means Medaka could not work out a model.
    if proc.returncode != 0:
        return None
    # An empty answer counts as "could not work one out" as well.
    model = proc.stdout.strip()
    if not model:
        return None
    return model


def download_model(model):
    """Make Medaka store the weights of the NAMED model on this computer, so that
    the polishing step later finds them without the internet (see "Downloading a
    named model here" at the top of this file). Runs
    `medaka tools download_models --models <model>`, which downloads the weights
    from GitHub only when they are not stored yet, and otherwise just reports
    where they are.

    Returns None on success, or the last lines of Medaka's messages as text when
    the command failed (for example: no internet, or no folder Medaka may write
    to). Called by main() once a named model has passed the check."""
    proc = subprocess.run(
        ["medaka", "tools", "download_models", "--models", model],
        capture_output=True, text=True,
    )
    if proc.returncode == 0:
        return None
    # Medaka ends a failed download with a Python traceback whose last line says
    # what went wrong; the last three lines are enough for the message.
    last_lines = proc.stderr.strip().splitlines()[-3:]
    return "\n".join(last_lines)


# ── Reading a model name: plain text handling, tested without Medaka ─────────

# The parts of a model name a user picks a model by (called "axes" in the code
# below), in the words the Medaka names use. Everything else in a name (the
# chemistry tag e82, the speed 400bps, neural-network suffixes such as
# rl_lstm384_dwells) stays in the name but is not used for matching. The
# basecaller version is shown in the tables but not matched on either. It is
# either gNNN, the version of the Guppy basecaller (g507 = Guppy 5.0.7), or, in
# newer names, vX.Y.Z, the version of the basecaller model (v5.2.0).
# _DEVICES also gives the full device name for the tables (min -> MinION).
# "high" is the accuracy word of older, now deprecated models (r941_min_high_g303).
_DEVICES = {"min": "MinION", "prom": "PromethION"}
_ACCURACIES = {"fast", "hac", "sup", "high"}


def tokenize(model):
    """Split a Medaka model name into the parts a user picks by: flowcell,
    device, accuracy and version (each None when the name leaves it out), plus
    is_variant (True for the variant-calling models). Names are joined with
    underscores but irregular: some carry a device word and some do not, and
    versions are either gNNN or vX.Y.Z. So every piece is classified on its own,
    instead of being read from a fixed position. Example: r941_min_hac_g507 ->
    flowcell r941, device min, accuracy hac, version g507.

    Returns a dictionary with those five entries. Called by consensus_models(),
    _fmt_table() and suggest()."""
    parts = model.split("_")
    out = {"flowcell": None, "device": None, "accuracy": None,
           "version": None, "is_variant": False}
    for p in parts:
        # re.fullmatch tests the whole piece against a regular expression (a
        # text pattern). The flowcell is "r" plus digits only: r941, r103, r1041.
        if re.fullmatch(r"r\d+", p):
            out["flowcell"] = p
        elif p in _DEVICES:
            out["device"] = p
        elif p in _ACCURACIES:
            out["accuracy"] = p
        # A version is "g" plus digits (g507), or "v", a digit, and then more
        # letters, digits or dots (v5.2.0).
        elif re.fullmatch(r"g\d+", p) or re.fullmatch(r"v\d[\w.]*", p):
            # Only the first version-like piece is kept (no name in Medaka
            # 2.2.2 has two).
            if out["version"] is None:
                out["version"] = p
        elif p in {"variant", "snp"}:
            out["is_variant"] = True
    return out


def consensus_models(available):
    """The models usable for medaka_consensus polishing: all of `available`
    except the variant/snp models (they are for variant calling and would be
    the wrong choice here). Keeps the order of `available`. Called by suggest()
    and auto_failure_message()."""
    return [m for m in available if not tokenize(m)["is_variant"]]


def _fmt_table(models):
    """Return `models` as a text table, one row per model in alphabetical
    order, with the columns a user picks by: model, flowcell, device, accuracy,
    version. Device words are spelled out (min -> MinION), since that is what
    people recognise, and a part the name leaves out shows as "-". The columns
    are padded with spaces so that they line up. Called by suggest() and
    auto_failure_message()."""
    rows = []
    for m in sorted(models):
        t = tokenize(m)
        rows.append((
            m,
            t["flowcell"] or "-",
            _DEVICES.get(t["device"], "-"),
            t["accuracy"] or "-",
            t["version"] or "-",
        ))
    headers = ("model", "flowcell", "device", "accuracy", "version")

    # The width of each column is the length of its longest entry, the header
    # included. With no rows at all, the headers alone set the widths.
    widths = []
    for column in range(len(headers)):
        width = len(headers[column])
        for row in rows:
            width = max(width, len(row[column]))
        widths.append(width)

    # One line of text per row: every cell padded with spaces to the width of
    # its column, and two spaces between the columns.
    def format_line(cells):
        padded_cells = []
        for column, cell in enumerate(cells):
            padded_cells.append(cell.ljust(widths[column]))
        return "  ".join(padded_cells)

    # The header line first, then one line per model.
    lines = [format_line(headers)]
    for row in rows:
        lines.append(format_line(row))
    return "\n".join(lines)


def suggest(bad_model, available):
    """Build the suggestion text printed when a named model is not available.

    Narrows the consensus models to those that agree with the bad name on the
    parts it does specify (flowcell, device, accuracy). A typo such as
    'r941_min_hac_g508' narrows to the r941 / MinION / hac rows, where the real
    'r941_min_hac_g507' is the obvious pick. When the bad name has none of
    those parts, or no model agrees with it, the text shows the whole consensus
    table instead. Called by main() for a bad named model."""
    pool = consensus_models(available)
    want = tokenize(bad_model)

    # True when model m agrees with the bad name on every part the bad name
    # specified. A part the bad name leaves out (None) does not count.
    def matches(m):
        t = tokenize(m)
        for axis in ("flowcell", "device", "accuracy"):
            if want[axis] is not None and t[axis] != want[axis]:
                return False
        return True

    # True when the bad name specifies at least one of the three parts.
    bad_name_has_a_part = False
    for axis in ("flowcell", "device", "accuracy"):
        if want[axis] is not None:
            bad_name_has_a_part = True

    # Narrow the list only when there is something to narrow by.
    narrowed = []
    if bad_name_has_a_part:
        for m in pool:
            if matches(m):
                narrowed.append(m)

    if narrowed:
        return (
            f"Models matching your flowcell/device/accuracy "
            f"(flowcell={want['flowcell']}, device={_DEVICES.get(want['device'], want['device'])}, "
            f"accuracy={want['accuracy']}):\n\n" + _fmt_table(narrowed)
        )
    return (
        "Could not match your model name to a flowcell/device/accuracy. "
        "All consensus models available in this Medaka version:\n\n" + _fmt_table(pool)
    )


def auto_failure_message(available):
    """The text printed when AUTO mode fails, that is when resolve_from_reads()
    returned None. The docstring of resolve_from_reads() lists all the causes;
    this text names the two usual ones (no basecaller tag in the reads, or one
    Medaka's lookup does not know). It tells the user to set a named model and
    shows the whole consensus table to pick from. Called by main() in auto mode
    only (a failed --fallback-to-auto uses suggest() instead)."""
    return (
        "Medaka could not infer a model from the reads: their FASTQ headers carry "
        "no basecaller tag, or one this Medaka version does not recognise. Set "
        "parameters.<mode>.medaka_model to one of the consensus models below (pick "
        "by your flowcell / device / accuracy), or to FALSE to skip Medaka.\n\n"
        + _fmt_table(consensus_models(available))
    )


# ── Decide the model: auto, named, or named-then-auto ────────────────────────

# main() reads the command line, decides the model and writes it to --out.
# argv is the list of command-line words (the tests pass one; None means "read
# the real command line"). Returns the exit status: 0 = model written to --out,
# 1 = no usable model (message on stderr, no --out).
#
# Three ways in, tried in this order: an empty --model means auto mode, so the
# reads decide; a named model is accepted when it is a model file already on
# disk, or a name this Medaka version has; and only with --fallback-to-auto does
# a bad name turn to auto mode instead of stopping the run. Otherwise a named
# choice is never swapped for a guess: it is accepted or reported, which is why
# the option is off by default.
def main(argv=None):
    ap = argparse.ArgumentParser(description="Validate/resolve the Medaka model early.")
    ap.add_argument("--model", default="", help="explicit model, or empty for auto")
    # CHANGED in FunFlux v2: the raw long-read file (plain or .gz), not the
    # filtlong output — see the module docstring.
    ap.add_argument("--reads", required=True,
                    help="raw long-read FASTQ (plain or .gz) for auto-inference")
    ap.add_argument("--fallback-to-auto", default="false",
                    help="true: on an invalid explicit model, fall back to auto")
    ap.add_argument("--out", required=True, help="write the resolved model name here")
    args = ap.parse_args(argv)

    model = args.model.strip()
    # --fallback-to-auto arrives as text: true, 1, yes or on (upper or lower
    # case, spaces around it ignored) switch it on; anything else leaves it off.
    fallback = args.fallback_to_auto.strip().lower() in {"true", "1", "yes", "on"}

    # Success: print the optional warning (note) to stderr, write the model to
    # --out as one line, report it on stdout, and return exit status 0.
    def succeed(resolved, note=None):
        if note:
            sys.stderr.write(note + "\n")
        with open(args.out, "w") as fh:
            fh.write(resolved + "\n")
        print(f"Medaka model resolved to: {resolved}")
        return 0

    # Failure: print the message to stderr after "ERROR: " and return exit
    # status 1. Nothing is written to --out.
    def fail(message):
        sys.stderr.write("ERROR: " + message + "\n")
        return 1

    # The model names this Medaka version has: used to check a named model and
    # to build the suggestion tables.
    available = list_available_models()

    # ── AUTO mode ────────────────────────────────────────────────────────────
    # resolve_model returns a PATH to the model file (e.g. .../r941_min_hac_g507_
    # model_pt.tar.gz), NOT a bare name, and only on success: a non-zero exit
    # means it could not work one out. So any answer is already a valid, usable
    # `-m` argument; trust it exactly as v1 did. Do NOT check it against the list
    # of bare names (a path can never match, which would fail every auto run).
    if not model:
        resolved = resolve_from_reads(args.reads)
        if resolved:
            return succeed(resolved)
        return fail(auto_failure_message(available))

    # ── NAMED (explicit) model ───────────────────────────────────────────────
    # A path to a model file on disk is accepted as it is: Medaka's -m accepts a
    # file, and only NAMED models can be checked against this version's list.
    # (v1 did the same: its polishing rule, consensus_long, checked the name
    # against the list only when `[[ ! -e "$model" ]]`, that is when no such
    # file existed.)
    if os.path.exists(model):
        return succeed(model)
    if model in available:
        # The name is right; now make sure its weights are on this computer, so
        # that the polishing step needs no internet (CHANGED in FunFlux v2).
        download_problem = download_model(model)
        if download_problem is not None:
            return fail(
                f"Medaka model '{model}' exists, but its weights are not stored on "
                f"this computer and could not be downloaded (Medaka downloads them "
                f"from GitHub). Start Snakemake on a computer that can reach the "
                f"internet (on a cluster this check runs where Snakemake runs, not on "
                f"a compute node), or set parameters.<mode>.medaka_model to a model whose "
                f"weights come with the Medaka package (for Medaka 2.2.2: "
                f"r1041_e82_400bps_hac_v6.0.0 or r1041_e82_400bps_sup_v5.2.0), or "
                f"to false to skip Medaka. "
                f"Medaka said:\n{download_problem}"
            )
        return succeed(model)

    # A bad name. With --fallback-to-auto (the opt-in middle way), try auto mode
    # before failing. As in auto mode above, an answer is a usable path.
    if fallback:
        resolved = resolve_from_reads(args.reads)
        if resolved:
            return succeed(
                resolved,
                note=(
                    f"WARNING: configured medaka_model '{model}' is not available in "
                    f"this Medaka version; fell back to the auto-inferred model "
                    f"'{resolved}' (parameters.<mode>.medaka_model_fallback_auto is "
                    f"on). Set medaka_model to auto to silence this."
                ),
            )
        # Auto mode could not help either: stop, with the suggestion table.
        return fail(
            f"Configured medaka_model '{model}' is not available, and fallback "
            f"auto-inference could not infer a usable model from the reads.\n\n"
            + suggest(model, available)
        )

    # A bad name and the option is off: stop, with the suggestion table.
    return fail(
        f"Configured medaka_model '{model}' is not available in this Medaka version.\n\n"
        + suggest(model, available)
    )


# Run main() only when the file is started as a program, not when the tests
# import it; its return value becomes the exit status.
if __name__ == "__main__":
    sys.exit(main())
