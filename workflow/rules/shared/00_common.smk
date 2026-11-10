# FunFlux v2.0.0 — shared foundation module (rules/shared/00_common.smk)
#
# workflow/Snakefile includes this file FIRST, before any rule file. It reads the
# configuration (the YAML file given with --configfile), checks it, finds the
# samples, and defines the names the rule files share: the mode and its capability
# flags, the path of every file one stage passes to the next, the tool settings, and
# small helper functions. The rule files (rules/<mode>/ for the active mode's front
# end, rules/shared/ for the stages all modes share) read these names and never work
# them out again. A few rule files also name a file that only they use (the MultiQC
# configuration in 90_report.smk, for example). Apart from the Snakefile, which uses
# MODE to pick the front-end folder, this is the only place that looks at the mode
# name.
#
# No rule lives in this file. Everything here is plain Python that runs once, when
# Snakemake reads the workflow, before any job runs. (Snakemake calls that moment
# "parse time".) It also runs on every dry run (snakemake -n).
#
# That is why the configuration is CHECKED here. A wrong mode, a missing input
# directory, a bad sample name, a missing annotation row or a typo in a parameter
# name stops the run within seconds and says what to fix, instead of failing hours
# into an assembly. There are two kinds of message:
#   * a STOP — sys.exit("[FunFlux] <config key>: what is wrong and what to do").
#              The run does not start.
#   * a NOTE — print("[FunFlux] NOTE: ..."). Something worth knowing; the run goes
#              ahead.
# A few plain lines ("Mode: ...", "Sample SAMPLE1 will be processed.", ...) are printed
# too, so the top of the log says what the run is going to do.
# miscellaneous/check_parse_time.py dry-runs good and deliberately broken
# configurations and checks these messages: run it after every change to this file.
#
# Takes in: the configuration; the input directory of the active mode (listed, to
#           find the samples); the annotation table (files.annotation_params).
# Produces: no files. Python names (constants and functions) that every rule file
#           reads, and the messages above.
#
# Sections, in file order:
#
# 1.  Mode             : read and check config["mode"]
#                        (illumina | longread | hybrid | contigs).
# 1b. Config readers   : small functions that read a block, a text, a number or a
#                        yes/no from the config and stop with a plain message when
#                        it is wrong.
# 2.  Capability flags : HAS_SHORT_READS, HAS_LONG_READS, HAS_READS, IS_HYBRID,
#                        CONTIGS_SCREEN, SCREEN and DECONTAM_READS. The rule files
#                        test these, not the mode name, to decide which rules exist.
# 3.  Output root and reference databases : OUT, the database folders and the
#                        files the user supplies.
# 4.  Resources        : CPUS, RAM, RUNTIME_HOURS, capped_cpus(), and the memory and
#                        time a heavy step asks a cluster scheduler for
#                        (HEAVY_MEM_MB, LONG_RUNTIME_MINUTES).
# 5.  Stage layout     : the numbered output folders, and the path of every file
#                        one stage passes to the next.
# 6.  Scripts and shared constants : script paths, allowed sample-name characters
#                        and file endings, the shared awk programs.
# 7.  Helpers          : functions that the sections below need.
# 8.  Screening settings : decontamination policy, mitochondrion screen, BLAST depth.
#                        The check of the decontamination tables is defined here
#                        and called at the end of section 12.
# 9.  Annotation and QC settings : masking, ITS, BUSCO, funannotate, GeneMark,
#                        Tiberius, SignalP 6, SPAdes, and the contig summary and
#                        the AMR screen (ChroQueTas), each with the inputs it has
#                        in this mode.
# 10. Database download links.
# 11. Long-read options : assembler, Medaka, Flye read mode, Filtlong, ONT adapter
#                        trimming, hifiasm and the hybrid-only Illumina guidance
#                        (longread and hybrid only).
#     11b. A check of every top-level key under `parameters`: typos, and
#                        settings left at their old v1 place.
# 12. Sample discovery  : lists the input directory of the ACTIVE mode only, reads
#                        the annotation table, checks each sample's AUGUSTUS model and
#                        prints the "AUGUSTUS:" line, and, when the screen runs, checks
#                        the decontamination tables against the samples.
# 13. Per-sample technology decisions : which tool flags each long-read sample
#                        gets (an ONT sample and a PacBio sample need different ones),
#                        and which ONT samples have their adapters trimmed.
# 14. Wildcard constraint : {sample} may only stand for a real sample name.
# 15. Targets           : the files `rule all` asks for (the whole run), and the
#                        files `rule download_all` asks for (only the downloads).
#
# Where a comment names a page under docs/ (docs/methods/blast-hit-depth.md, for
# example), that page of the user documentation explains a choice made when the
# four fungal workflows (FunFlux, Funnotator, FunFluxL, FunFluxL+) were merged
# into this one.
#
# House rule for anyone editing this file: anything that depends on the INPUT DATA
# (listing an input directory, a message about samples, a stop on a bad sample
# name) must run only for the active mode: inside that mode's `if MODE == ...`
# branch (section 12), or after it, working only from what the branch found
# (SAMPLES, READ_TECH), as the annotation-table check, the decontamination-table
# check and section 13 do. So a run in one mode never looks at another mode's input
# directory or stops on its errors. Settings that do not depend on the input
# (database links, the decontamination policy) are read in every run.
#
# Every name defined here is visible in every rule file: Snakemake runs all included
# files in one shared space of names. A leading underscore (_note, _mito, ...) marks
# a name meant for this file only, and a new name must not reuse one from another file.

# Python's standard tools used below: csv reads the tab-separated decontamination
# tables; glob lists files that match a pattern; importlib.util loads
# check_augustus_model.py as a module (section 9d); math.isfinite tells a real number
# from nan and inf; os handles paths and tests files; re checks text against a
# regular expression (a pattern that text must match); sys.exit stops the run with
# a message; Path turns the output directory into an absolute path; unquote and
# urlparse take a download link apart.
import csv
import glob
import importlib.util
import math
import os
import re
import sys
from pathlib import Path
from urllib.parse import unquote, urlparse

# glob_wildcards is a Snakemake helper: given a pattern like "{sample}_R1.{extn}"
# it lists the files that match and returns what the wildcards matched. (A
# wildcard is a placeholder in braces, such as {sample}, that Snakemake fills in.)
# Imported explicitly so the reader can see where it comes from.
from snakemake.io import glob_wildcards


# ═══════════════════════════ 1. Mode ═══════════════════════════
# `config` is the dictionary Snakemake builds from the --configfile YAML file; it is
# available here without an import. The configuration is REQUIRED on the command
# line. The Snakefile deliberately declares no default `configfile:` (see the note
# there), so an empty config means the user forgot --configfile. Say so, rather than
# failing later with a Python error about a missing key.
if not config:
    sys.exit(
        "[FunFlux] No configuration supplied. Pass one explicitly, e.g.:\n"
        "  snakemake --sdm conda --cores N --configfile config/config.yaml"
    )

# The four accepted values of config["mode"].
MODES = ("illumina", "longread", "hybrid", "contigs")

# The mode is REQUIRED, and it is the one key that tells a v2 configuration from a
# v1 one: FunFlux, FunFluxL and FunFluxL+ configurations had no `mode:` key at all.
# So a missing mode gets a message that explains the change instead of just "not
# set". A missing, null or blank value counts as missing. strip() drops spaces
# around the word and lower() makes "Illumina" and "illumina" the same.
if not str(config.get("mode") or "").strip():
    sys.exit(
        "[FunFlux] mode: this configuration has no 'mode' key. FunFlux v2 needs one of "
        "illumina | longread | hybrid | contigs, and an 'input:' block that names the "
        "input directory for that mode. A configuration written for FunFlux, FunFluxL "
        "or FunFluxL+ v1 does not run unchanged: start from the config template shipped "
        "with FunFlux v2 (see 'Coming from v1' in the documentation)."
    )
MODE = str(config["mode"]).strip().lower()
if MODE not in MODES:
    sys.exit(
        "[FunFlux] mode: must be one of illumina | longread | hybrid | contigs "
        f"(got: {config['mode']!r})"
    )

# One line per mode, printed once so the top of the log says which pipeline is
# running. Purely cosmetic.
MODE_TAGLINE = {
    "illumina": "Genome assembly and annotation from Illumina paired-end reads",
    "longread": "Genome assembly and annotation from ONT or PacBio HiFi long reads",
    "hybrid":   "Genome assembly and annotation from long reads and Illumina paired-end reads together",
    "contigs":  "Annotation of fungal genome assemblies you already have",
}
print(f"Mode: {MODE} — {MODE_TAGLINE[MODE]}.")


# ═══════════════════════ 1b. Config readers ════════════════════
# Small functions that read one kind of value from the configuration and stop with
# a plain message when it is wrong. They exist so that every check below reads the
# same way and every message names the config key the user has to fix. The
# `dotted_key` that the yes/no and number readers take (for example
# "parameters.blast_max_hsps") is only used in that message.

def _note(message):
    # Print a NOTE: a message worth knowing that does not stop the run. Used all
    # through this file.
    print(f"[FunFlux] NOTE: {message}")


def _section(name):
    # Return the top-level block config[name] (for example "directories") as a
    # dict. `or {}` turns a block that is missing, or present but empty (YAML null),
    # into an empty dict. A block that is not a set of key: value lines at all (for
    # example `input: /some/dir`) is a mistake worth naming.
    # Called by _required_text and _optional_text, and for the `parameters` and
    # `resources` blocks.
    block = config.get(name, {}) or {}
    if not isinstance(block, dict):
        sys.exit(f"[FunFlux] {name}: must be a block of 'key: value' lines.")
    return block


def _required_text(section_name, key):
    # A text value that MUST be set (a path, a link), read from
    # config[section_name][key]. Returns it as text, with the spaces around it
    # removed. A key that is missing, null or empty stops the run: the shipped config
    # template carries empty placeholders on purpose, so an untouched template
    # fails here with the key's name rather than running against a wrong path.
    # Used in section 3 (output directory, databases, user files) and section 10
    # (download links).
    value = str(_section(section_name).get(key) or "").strip()
    if not value:
        sys.exit(f"[FunFlux] {section_name}.{key}: required, but it is missing or empty in the config.")
    return value


def _optional_text(section_name, key):
    # A text value that may be left empty, read the same way as _required_text.
    # Returns "" when it is not set; every caller treats "" as "not given". Used for
    # the optional folders and files (sections 3 and 9f), the UNITE link (section
    # 10) and the input directories (section 12).
    return str(_section(section_name).get(key) or "").strip()


def _config_bool(dotted_key, value, default=False):
    # Read a yes/no value predictably: true/yes/on/1 mean True, false/no/off/0 mean
    # False, a missing or null value (`key:` with nothing after it) means `default`.
    # Anything else, an empty quoted "" included, is a mistake and stops the run.
    # YAML 1.1 already turns a bare yes/no/on/off into True/False before we see
    # it, which is why a real True/False is accepted first.
    # Used for parameters.contigs.screen, parameters.decontamination.discard_no_hit,
    # the two RepeatModeler switches, parameters.its.unite_identified_only,
    # parameters.prediction.augustus_optimize, and the long-read yes/no settings
    # (section 11).
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "on"}:
        return True
    if text in {"false", "0", "no", "off"}:
        return False
    sys.exit(f"[FunFlux] {dotted_key}: expected true or false (got {value!r}).")


def _to_int(dotted_key, value, minimum=None, maximum=None):
    # Read a whole number and check its range (minimum and maximum are both
    # allowed values); stop with the key's name if it is not one. "12" and 12 are
    # both accepted; a decimal such as 2.5, an empty value and a True/False are
    # refused (Python would otherwise treat True as 1). Used for every whole-number
    # setting, for example the threads, the BLAST depths and the ITS lengths.
    if isinstance(value, bool):
        sys.exit(f"[FunFlux] {dotted_key}: expected a whole number (got {value!r}).")
    try:
        number = int(str(value).strip())
    except (TypeError, ValueError):
        sys.exit(f"[FunFlux] {dotted_key}: expected a whole number (got {value!r}).")
    if minimum is not None and number < minimum:
        sys.exit(f"[FunFlux] {dotted_key}: must be at least {minimum} (got {number}).")
    if maximum is not None and number > maximum:
        sys.exit(f"[FunFlux] {dotted_key}: must be at most {maximum} (got {number}).")
    return number


def _to_float(dotted_key, value, minimum=None, maximum=None, minimum_exclusive=False):
    # Same as _to_int for decimal numbers. minimum_exclusive=True means the value
    # must be strictly greater than `minimum` (used for "above zero"). Used for the
    # mitochondrion thresholds, the ITS cutoff and the Filtlong shares and weights.
    if isinstance(value, bool):
        sys.exit(f"[FunFlux] {dotted_key}: expected a number (got {value!r}).")
    try:
        number = float(str(value).strip())
    except (TypeError, ValueError):
        sys.exit(f"[FunFlux] {dotted_key}: expected a number (got {value!r}).")
    # float() also accepts the words nan ("not a number") and inf or infinity, with
    # or without a sign, and YAML reads .nan and .inf as such numbers. None of them
    # is a usable setting: nan passes every range test below (any comparison with
    # nan is False) and would show up as "nan%" in the mitochondrion summary line,
    # and inf is no threshold. So they get the same message as any other value that
    # is not a number.
    if not math.isfinite(number):
        sys.exit(f"[FunFlux] {dotted_key}: expected a number (got {value!r}).")
    if minimum is not None:
        # Two ways to be too small, with the matching wording for the message.
        if minimum_exclusive:
            too_small = number <= minimum
            relation = "above"
        else:
            too_small = number < minimum
            relation = "at least"
        if too_small:
            sys.exit(f"[FunFlux] {dotted_key}: must be {relation} {minimum} (got {number}).")
    if maximum is not None and number > maximum:
        sys.exit(f"[FunFlux] {dotted_key}: must be at most {maximum} (got {number}).")
    return number


# Every tool-specific setting lives under config["parameters"]. Read once here;
# sections 2, 8, 9, 11 and 11b take their settings from it.
PARAMETERS = _section("parameters")


def _check_block(dotted_name, block, allowed_keys):
    # Return `block` as a dict (empty if it is absent) and stop on any key that is not
    # in allowed_keys. An unknown key is almost always a typo, and a typo in a tuning
    # parameter silently falls back to the default, which is the worst kind of
    # mistake, so it is named together with the keys that do exist.
    # Called by _parameter_block below, and for the filtlong and hifiasm blocks
    # inside parameters.longread and parameters.hybrid (section 11).
    block = block or {}
    if not isinstance(block, dict):
        sys.exit(f"[FunFlux] {dotted_name}: must be a block of 'key: value' lines.")
    unknown = sorted(set(block) - set(allowed_keys))
    if unknown:
        sys.exit(
            f"[FunFlux] {dotted_name}: unsupported key(s): {', '.join(unknown)}. "
            f"Supported keys: {', '.join(sorted(allowed_keys))}."
        )
    return block


def _parameter_block(name, allowed_keys):
    # The block parameters.<name> as a checked dict: empty when the block is absent,
    # and a stop on any key it does not know. Used for the contigs, decontamination,
    # mitochondria, its, busco, prediction, genemark and signalp6 blocks, and for
    # parameters.longread or parameters.hybrid (section 11).
    return _check_block(f"parameters.{name}", PARAMETERS.get(name), allowed_keys)


# Words that mean "nothing / switched off" for a setting that is normally a
# number or a name. How YAML hands such a value over: an empty value or an
# unquoted null arrives as None; an unquoted off, no or false arrives as the
# boolean False; the word none (which YAML does not know) and anything in quotes
# ('null', '') arrive as text. All of them are treated the same.
EMPTY_WORDS = {"", "none", "null", "off", "false", "no"}


def _is_empty_value(value):
    # True when a setting means "nothing / switched off" (see EMPTY_WORDS just above).
    # A YAML True is NOT empty. Used for parameters.repeatmodeler_seed and, in
    # section 11, for filtlong.target_bases, hifiasm.genome_size,
    # hifiasm.telomere_motif and filtlong_split.
    if value is None:
        return True
    if isinstance(value, bool):
        return value is False
    return str(value).strip().lower() in EMPTY_WORDS


def _mode_word(value, default):
    # Read a setting whose value is a word such as off | auto | report. The trap: YAML
    # (version 1.1, which PyYAML and Snakemake read) turns an UNQUOTED off, no or false into
    # the boolean False, and on, yes or true into True, so `mode: off` in a configuration
    # arrives here as False and not as the text "off". Looking only at the text would reject
    # what the user clearly wrote. So False is read as "off", an empty value as the default,
    # and anything else as its lower-case text; a True is left as "true", which no mode
    # accepts, so the caller's own message says what the allowed words are.
    # Used for parameters.decontamination.mode (section 8a),
    # parameters.mitochondria.mode (section 8b) and parameters.genemark.mode
    # (section 9e).
    if value is None:
        return default
    if value is False:
        return "off"
    word = str(value).strip().lower()
    if word == "":
        return default
    return word


# ═══════════════════ 2. Capability flags ═══════════════════
# These True/False flags, not the mode name, decide which rules exist. A rule file
# wraps its rules in, for example, `if SCREEN:`, so the test says what the rules
# need instead of listing modes. A rule that needs something a mode does not have
# is never defined in that mode, so nothing waits for it. The comment after each
# flag names the steps that need it.
IS_HYBRID       = MODE == "hybrid"                 # the Illumina guide assembly, guided Filtlong, Polypolish
HAS_SHORT_READS = MODE in ("illumina", "hybrid")   # PhiX removal, fastp, SPAdes
HAS_LONG_READS  = MODE in ("longread", "hybrid")   # Filtlong, NanoPlot, hifiasm/Flye, Medaka
HAS_READS       = MODE != "contigs"                # anything that writes into 01.reads/ at all

# contigs mode (the old Funnotator) annotates a genome you already have. By default
# it does exactly that and nothing else. parameters.contigs.screen: true adds the
# same contamination and mitochondrion screen the read-based modes run, on your
# FASTA. Why it is optional: a provider's assembly usually still contains
# the mitochondrial genome, and annotating it with the nuclear gene predictor is
# wrong — but the screen needs the NCBI core_nt BLAST database, which contigs-mode
# users were told they did not need. The setting is read and checked in every mode,
# so a mistyped value stops the run even where it has no effect; it only changes
# something in contigs mode (the read-based modes always run the screen).
# docs/modes/contigs.md describes the mode with and without the screen.
_contigs_block = _parameter_block("contigs", {"screen"})
_contigs_screen_setting = _config_bool(
    "parameters.contigs.screen", _contigs_block.get("screen"), False
)
CONTIGS_SCREEN = MODE == "contigs" and _contigs_screen_setting

# SCREEN: True when this run includes the contamination + mitochondrion screen. Every
# read-based mode always does; contigs mode only when it was switched on above.
# Two shared modules (10_decontam, 15_mitochondria) define NO rules when this is
# False, and the BLAST database is only required when it is True (section 3).
SCREEN = HAS_READS or CONTIGS_SCREEN

# Which reads give BlobTools its coverage track in the contamination screen (BlobTools
# plots GC against coverage to separate genomes). 10_decontam.smk reads this to
# decide which version of rule map_contigs exists:
#   "illumina" — the Illumina reads, mapped to the SPAdes draft with bowtie2.
#   "long"     — the long reads, mapped to the long-read draft with minimap2.
#   "none"     — contigs mode has no reads, so the contigs are mapped to themselves.
#                The coverage BlobTools then reports reflects only contig length
#                (10_decontam.smk explains why) and says nothing about the
#                organism; only the BLAST taxonomy does real work there.
# The hybrid Illumina "guide" assembly is mapped by its own rules (index_guide_contigs and
# map_guide_contigs in rules/hybrid/25_guide.smk), always with bowtie2.
if MODE == "illumina":
    DECONTAM_READS = "illumina"
elif HAS_LONG_READS:
    DECONTAM_READS = "long"
else:
    DECONTAM_READS = "none"

# Without the screen, the only change made to the user's FASTA before QC and
# annotation is rule filter_contigs (contigs/20_assembly.smk), so the NOTE says
# what that rule does. It looks at the first header: SPAdes-style headers carry
# length and coverage, and contigs under 500 bp or under 2x coverage are dropped;
# for any other FASTA nothing is dropped, and each header is cut at its first space
# or tab.
if MODE == "contigs" and not CONTIGS_SCREEN:
    _note(
        "contigs mode without parameters.contigs.screen: your FASTA gets no contamination or "
        "mitochondrion screen. If its headers are SPAdes-style (>NODE_<n>_length_<L>_cov_<C>), "
        "rule filter_contigs drops contigs under 500 bp or under 2x coverage; otherwise it only "
        "cuts each header at its first space or tab. A mitochondrial genome in your FASTA, if "
        "present, will be annotated as nuclear. Set parameters.contigs.screen: true to find and "
        "split it off (needs directories.blast_db)."
    )


# ═══════════ 3. Output root and reference databases ═══════════
# OUT is the output directory as ONE absolute path, and every output path in section
# 5 is built from it. The v1 FunFlux workflows used Snakemake's `workdir:` directive
# instead, which changes the working directory while the workflow is still being
# READ, even for a dry run. Two problems followed: a placeholder output directory
# that did not exist yet crashed the dry run, and every input directory in the config
# had to be absolute, because a relative one was silently read from inside the output
# directory. With an absolute OUT the working directory never changes: input paths
# mean what the user expects (relative to where snakemake was started, or absolute),
# and output_dir only has to exist when a rule first writes there. resolve() makes
# the path absolute and also replaces any symbolic link in it by its target. OUT is
# kept as plain text, not as a Path object, because every path below is built from it
# with "+".
OUT = str(Path(_required_text("directories", "output_dir")).resolve())

# The three reference databases the user prepares before the first run.
#   blast_db       : NCBI nucleotide BLAST database (core_nt), used by the
#                    contamination and mitochondrion screen (rules blast_contigs and
#                    blob_json). Only needed when that screen runs (every read-based
#                    mode; contigs mode with screen), and only then required here.
#   eggnog_db      : eggNOG-mapper database (rule eggnog_annotation).
#   funannotate_db : the funannotate database snapshot (every funannotate rule). If
#                    the directory does not hold one yet the workflow downloads and
#                    checks it (rule funannotate_db, from links.funannotate_link), so
#                    this path may point at a place that does not exist yet.
# They are NOT opened here, so a dry run on a machine that does not hold the databases
# still works. A NOTE says when the eggNOG or the BLAST folder is not there (the
# funannotate folder may be missing on purpose: it is then downloaded).
#
# Each path is made ABSOLUTE right here (os.path.abspath, measured from the folder
# Snakemake was started in). A relative path would usually work, because the workflow
# itself never changes folder, but some rules change into another folder before they
# start their tool (funannotate_preprocess and funannotate_annotation do, and both use
# the funannotate database), and there a relative path points nowhere. One conversion
# here removes that whole class of problems.
EGGNOGDB       = os.path.abspath(_required_text("directories", "eggnog_db"))
FUNANNOTATEDB  = os.path.abspath(_required_text("directories", "funannotate_db"))
if SCREEN:
    BLASTDB = os.path.abspath(_required_text("directories", "blast_db"))
else:
    _blast_text = _optional_text("directories", "blast_db")
    BLASTDB = os.path.abspath(_blast_text) if _blast_text else ""

if not os.path.isdir(EGGNOGDB):
    _note(f"directories.eggnog_db '{EGGNOGDB}' does not exist on this machine; the eggNOG step will fail without it.")
if SCREEN and not os.path.isdir(BLASTDB):
    _note(f"directories.blast_db '{BLASTDB}' does not exist on this machine; the contamination screen will fail without it.")

# GeneMark-ES is an optional, separately licensed gene predictor (checked in
# section 9e). Its directory is read here, with the other paths. It is not made
# absolute: section 9e requires the user to give an absolute path.
GENEMARKDIR = _optional_text("directories", "genemark_dir")


def _resolve_optional_db_dir(config_key, probe_relpath, hint):
    # Read an OPTIONAL directories.<config_key> override and check it now, so a bad
    # path stops the run before any job starts instead of failing hours in.
    # probe_relpath is a file or folder that can only exist inside a real copy of
    # this particular database, so a directory that merely exists but holds the
    # wrong thing is caught too. hint is the last sentence of the stop message.
    # Returns the absolute path, or "" when the key is unset, which every caller
    # reads as "download it". Called twice, just below (antiSMASH and UNITE).
    path = _optional_text("directories", config_key)
    if not path:
        return ""
    if not os.path.isdir(path):
        sys.exit(f"[FunFlux] directories.{config_key}: '{path}' is not a directory. {hint}")
    if not os.path.exists(os.path.join(path, probe_relpath)):
        sys.exit(f"[FunFlux] directories.{config_key}: '{path}' has no '{probe_relpath}'. {hint}")
    # Absolute, because the rule that uses it builds symbolic links to files inside it, and a
    # link's relative target is read from the folder the link sits in, not from here.
    return os.path.abspath(path)


# Two databases the workflow can download itself, but that are large and shared by
# many projects. If the user already holds a copy, pointing at it means nothing is
# downloaded into every run's output directory. The workflow never writes to such a
# directory: the rule that uses it builds a local VIEW (symbolic links) inside
# 00.resources/ instead (rules antismash_db_local and unite_db_local). When the key is
# set it takes precedence over the download, and the download rule is not defined.
ANTISMASHDB = _resolve_optional_db_dir(
    "antismash_db", "clusterblast",
    "Point it at an antiSMASH databases directory (the one that holds clusterblast/, pfam/, ...)."
)
if ANTISMASHDB:
    print(f"Using the local antiSMASH database at '{ANTISMASHDB}'. Nothing will be downloaded.")

UNITEDB = _resolve_optional_db_dir(
    "unite_db", "unite_its_sintax.fasta",
    "Point it at a directory that holds the UNITE SINTAX reference as an uncompressed 'unite_its_sintax.fasta'."
)
if UNITEDB:
    print(f"Using the local UNITE database at '{UNITEDB}'. Nothing will be downloaded.")

# Files the user supplies (not directories).
#   iprscan           : the InterProScan launcher script, run by rule
#                       iprscan_annotation. InterProScan is not bioconda-installable
#                       and is several tens of GB, so it is installed by hand and
#                       called from here. A NOTE says when the file is not there.
#   annotation_params : one row per sample — species name, protein evidence and,
#                       optionally, the stored AUGUSTUS model to use (empty: AUGUSTUS
#                       is trained on the genome; section 9d). Read and checked in
#                       section 12.
IPRSCAN    = os.path.abspath(_required_text("files", "iprscan"))
ANNOTATION = os.path.abspath(_required_text("files", "annotation_params"))
if not os.path.isfile(IPRSCAN):
    _note(f"files.iprscan '{IPRSCAN}' does not exist on this machine; the InterProScan step will fail without it.")


# ════════════════════════ 4. Resources ═════════════════════════
# The single source of truth for compute limits. Every CPU-bound rule declares
# Snakemake's built-in `threads:` (capped_cpus(N) below), so `--cores` enforces
# the machine limit on its own; pass the SAME number to --cores that you set here.
# (The v1 workflows used a custom `resources: cpus`, which Snakemake does not count
# against --cores: a run started with --cores 8 still ran tools on 25 threads.)
# On a single computer, do not also pass --jobs/-j: there it is an alias for --cores.
# On a cluster (snakemake --executor slurm and the like) Snakemake needs --jobs: the
# number of jobs it may have queued at once. Each job then asks the scheduler for the
# CPUs of its rule's threads, and for the memory and time below.
#   CPUS          : resources.threads (default 24), the thread budget of the whole run.
#   RAM           : resources.ram_gb (default 128), the memory one heavy step may use,
#                   in GB. It is the cap SPAdes is given (-m, rule illumina_assembly)
#                   and, capped at 64 GB, the Java heap of Qualimap (rules
#                   map_evaluation and map_guide_evaluation). The rules take it as a
#                   Snakemake `resources:` entry, not as a param: a memory cap changes
#                   no result, so a changed ram_gb does not re-run them.
#   RUNTIME_HOURS : resources.runtime_hours (default 48), the longest a long step may
#                   run, in hours.
#
# What a cluster scheduler is told. Snakemake hands a job's `mem_mb` (memory in MB) and
# `runtime` (minutes) resources to the cluster scheduler, which reserves that much and
# stops the job when it uses more. A rule that declares neither asks for nothing
# (checked with Snakemake 9.14.6), so the scheduler's own defaults apply: often a GB or
# two of memory per CPU and the queue's default time limit. That is too little for the
# assemblers (hifiasm needed 17 GB on a 42 Mb genome, SPAdes about 15 GB and Flye 11 GB)
# and for the annotation tools, and InterProScan alone runs for hours. So the rules write
# one of two declarations:
#   * a heavy step (assemblers, read mapping with sorting, BLAST, BUSCO, Medaka,
#     Polypolish, repeat masking, the gene finders and the annotation tools) asks for
#     HEAVY_MEM_MB of memory and LONG_RUNTIME_MINUTES of time;
#   * a step that reads every read but holds little in memory (the PhiX filter, the
#     read selection of the hybrid mode, Filtlong in longread mode) asks only for
#     LONG_RUNTIME_MINUTES.
# Qualimap asks for its Java heap plus 4 GB for Java itself, but never for more than
# HEAVY_MEM_MB (rules map_evaluation and map_guide_evaluation). Every other step is
# light (seconds to minutes, a few GB at most) and declares nothing; on a cluster whose
# defaults are smaller than that, the user gives the light steps a default on the
# command line, for example --default-resources mem_mb=8000 runtime=240 (a rule's own
# declaration wins over it).
#
# On a single computer these two values change nothing: Snakemake's local runs count
# only threads against --cores, and ignore mem_mb and runtime unless the user limits
# them with --resources (for example --resources mem_mb=500000). A changed ram_gb or
# runtime_hours re-runs nothing either: Snakemake does not compare resources with
# those of the last run.
_resources = _section("resources")
CPUS = _to_int("resources.threads", _resources.get("threads", 24), minimum=1)
RAM  = _to_int("resources.ram_gb", _resources.get("ram_gb", 128), minimum=1)
RUNTIME_HOURS = _to_int("resources.runtime_hours", _resources.get("runtime_hours", 48), minimum=1)

# The two values the rules declare (see above): memory in MB (1 GB = 1024 MB) and
# time in minutes, the units Snakemake's mem_mb and runtime resources use.
HEAVY_MEM_MB = RAM * 1024
LONG_RUNTIME_MINUTES = RUNTIME_HOURS * 60


def capped_cpus(n):
    # Cap a rule's thread request at n: returns n, or CPUS if that is smaller. Some
    # tools stop scaling (or misbehave) past a certain thread count; rules write
    # `threads: capped_cpus(N)` (10_decontam, 20_qc, 30_taxonomy and 40_annotation
    # do) instead of spelling out min(CPUS, N) in dozens of places. Section 9f caps
    # the SignalP 6 thread count with it too.
    return min(CPUS, n)


# ═══════════ 5. Stage layout: output folders and the files passed between stages ═══════════
# The numbered output directories are the same in all four modes. A mode
# decides only how a set of contigs gets made; everything from taxonomy onwards
# lives in the same place whichever way the genome arrived. Every constant below
# is an ABSOLUTE path built from OUT (section 3). Rules use these names and append
# nothing but a file or folder name, so no rule hard-codes a stage number or the
# output root. docs/reference/output.md describes the layout.
#
# "{sample}" inside these strings is a Snakemake wildcard: a rule that lists the
# constant as an input or output gets one concrete path per sample, with the
# sample's name filled in. Python code outside a rule fills it in itself, with
# .format(sample=...) (sections 12 and 13) or expand() (section 15).
DIR_RESOURCES  = OUT + "/00.resources"   # run-level material, never per sample: downloaded databases, SignalP 6 install
DIR_READS      = OUT + "/01.reads"       # read QC and filtering; empty in contigs mode
DIR_ASSEMBLY   = OUT + "/02.assembly"    # assembly, polishing, decontamination, mitochondrion split, assembly QC
DIR_TAXONOMY   = OUT + "/03.taxonomy"    # ITS extraction and UNITE classification
DIR_ANNOTATION = OUT + "/04.annotation"  # repeat masking, gene prediction, functional annotation
DIR_AMR        = OUT + "/05.amr"         # antifungal resistance mutations (ChroQueTas); made only when the step is on
DIR_REPORT     = OUT + "/06.report"      # MultiQC
LOGS           = OUT + "/logs"           # one log per rule and sample (one per rule for the once-per-run rules)

# ── 00.resources: files that belong to the run, not to a sample ──
# Each is made once per run, and the jobs of every sample share it:
#   UNITE_FASTA, UNITE_PROVENANCE   : made by rule unite_db or unite_db_local; the
#                                     whole UNITE release. Rule unite_identified reads
#                                     the FASTA (and rule its_classification does, when
#                                     the filter below is off).
#   UNITE_FASTA_IDENTIFIED,
#   UNITE_IDENTIFIED_COUNTS         : made by rule unite_identified from UNITE_FASTA: the
#                                     UNITE records whose taxonomy names a genus, and a
#                                     two-line text with the number of records kept and
#                                     dropped. Rule its_classification reads the FASTA.
#                                     UNITE_REFERENCE, the file vsearch is given, is set
#                                     in section 9b to one of the two FASTA files, after
#                                     parameters.its.unite_identified_only is read.
#   ANTISMASH_DB_DIR                : made by rule antismash_db or antismash_db_local;
#                                     read by rule antismash_annotation.
#   FUNANNOTATE_DB_MARKER           : written by rule normalize_funannotate_db once the
#                                     funannotate database is ready; rules
#                                     funannotate_preprocess, augustus_model_check and
#                                     funannotate_prediction wait for it, and rule all
#                                     asks for it.
#   SIGNALP6_DIR, SIGNALP6_MANIFEST : rule install_signalp6 installs SignalP 6 into its
#                                     conda environment, moves the unpacked package with
#                                     the model weights here (signalp6/package) and
#                                     writes the manifest; rule signalp6_annotation
#                                     reads the manifest.
#   TIBERIUS_MODEL_DIR              : the Tiberius fungi model (the trained network),
#                                     downloaded as TIBERIUS_MODEL_ARCHIVE, checked and
#                                     unpacked by rule tiberius_model (42_tiberius.smk)
#                                     when Tiberius is on. TIBERIUS_MODEL_CONFIG and
#                                     TIBERIUS_MODEL_WEIGHTS are two of its three files,
#                                     the ones Snakemake checks for; rule
#                                     tiberius_prediction reads the whole folder.
#   TELOMERE_DB,
#   TELOMERE_DB_PROVENANCE          : the telomere repeat database of tidk, one row per
#                                     species whose chromosome ends have been examined
#                                     (its telomere repeat unit and how well it is
#                                     supported). Rule telomere_db (25_contig_summary.smk)
#                                     downloads it from links.telomere_db_link, checks it
#                                     against parameters.contig_summary.telomere_db_sha256
#                                     and writes the PROVENANCE file (the address, the
#                                     date, the checksum and the row count). Rule
#                                     telomere_units reads the database for every sample.
#                                     The database's repository carries no licence file,
#                                     so the file is fetched into the run, never shipped.
#   UNIVEC_ILLUMINA,
#   UNIVEC_PROVENANCE               : the Illumina adapter and primer records of NCBI's
#                                     UniVec (build 10.1), for the contig adapter screen.
#                                     Rule univec_illumina (25_contig_summary.smk) downloads
#                                     UniVec from links.univec_link, checks it against
#                                     parameters.contig_summary.univec_sha256, keeps the
#                                     records whose description names Illumina and writes
#                                     the PROVENANCE file. Illumina's oligonucleotide
#                                     sequences are copyrighted, so they are fetched into
#                                     the run, never shipped. Rule blast_adapters adds them
#                                     to the adapters shipped in workflow/data/adapters/.
#   CHROQUETAS_DIR and the six      : ChroQueTas and the FungAMR database of antifungal
#   CHROQUETAS_* paths below          resistance mutations, for the AMR screen when it is
#                                     on (section 9i). Rule chroquetas_db (50_amr.smk)
#                                     downloads them from links.chroquetas_link, checks the
#                                     script and the database archive against the two
#                                     checksums of parameters.amr, and unpacks the
#                                     database. Rule amr_chroquetas reads the script, the
#                                     species list and the database folder for every
#                                     sample. ChroQueTas and its database are licensed
#                                     CC BY-NC-ND 4.0, so they are fetched into the run,
#                                     never shipped.
#   PHIX_FASTA, PHIX_INDEX_PREFIX   : the PhiX174 genome, which Illumina libraries are
#                                     spiked with (rule download_phix fetches it from
#                                     links.phix_link), and its bowtie2 index (rule
#                                     build_phix). Rule map_phix maps every sample's
#                                     Illumina reads against it and keeps the reads
#                                     that do NOT match. Illumina and hybrid mode only.
UNITE_FASTA             = DIR_RESOURCES + "/unite/unite_its_sintax.fasta"
UNITE_PROVENANCE        = DIR_RESOURCES + "/unite/PROVENANCE.txt"   # where the UNITE file came from, when, and its md5
UNITE_FASTA_IDENTIFIED  = DIR_RESOURCES + "/unite/unite_its_sintax.identified.fasta"   # the records with a genus
UNITE_IDENTIFIED_COUNTS = DIR_RESOURCES + "/unite/identified_filter.txt"              # how many were kept and dropped
ANTISMASH_DB_DIR        = DIR_RESOURCES + "/antismash/databases"
FUNANNOTATE_DB_MARKER   = DIR_RESOURCES + "/funannotate_db_paths.normalized.json"
SIGNALP6_DIR            = DIR_RESOURCES + "/signalp6"
SIGNALP6_MANIFEST       = SIGNALP6_DIR + "/installation.json"
TIBERIUS_MODEL_ARCHIVE = DIR_RESOURCES + "/tiberius/model.tar.gz"     # deleted once it is unpacked
TIBERIUS_MODEL_DIR     = DIR_RESOURCES + "/tiberius/model"
TIBERIUS_MODEL_CONFIG  = TIBERIUS_MODEL_DIR + "/model_config.json"    # the network's settings
TIBERIUS_MODEL_WEIGHTS = TIBERIUS_MODEL_DIR + "/weights.h5"           # the trained weights
TELOMERE_DB            = DIR_RESOURCES + "/telomeres/telo_database_v2.tsv"   # tidk's telomere repeat database, one row per species
TELOMERE_DB_PROVENANCE = DIR_RESOURCES + "/telomeres/PROVENANCE.txt"         # where it came from, when, its checksum and row count
UNIVEC_ILLUMINA        = DIR_RESOURCES + "/adapters/univec_illumina.fasta"   # the Illumina records of UniVec 10.1
UNIVEC_PROVENANCE      = DIR_RESOURCES + "/adapters/PROVENANCE.txt"          # where they came from, when, the checksum and record count
CHROQUETAS_DIR          = DIR_RESOURCES + "/chroquetas"
CHROQUETAS_SCRIPT       = CHROQUETAS_DIR + "/ChroQueTas.sh"         # the program: one bash script, run as `bash ChroQueTas.sh ...`
CHROQUETAS_DB_ARCHIVE   = CHROQUETAS_DIR + "/db-v20250811.tgz"      # the database as downloaded; deleted once it is unpacked
CHROQUETAS_DB_DIR       = CHROQUETAS_DIR + "/FungAMR_db"            # the unpacked database: one folder per species
CHROQUETAS_SPECIES_LIST = CHROQUETAS_DIR + "/species.txt"           # the species folder names, one per line; written once the database is unpacked, so it marks it as ready
CHROQUETAS_PROVENANCE   = CHROQUETAS_DIR + "/PROVENANCE.txt"        # where the files came from, when, their checksums, the species count, the licence
CHROQUETAS_LICENSE      = CHROQUETAS_DIR + "/LICENSE"               # the licence text of ChroQueTas and its database (CC BY-NC-ND 4.0)
PHIX_DIR          = DIR_RESOURCES + "/phix"
PHIX_FASTA        = PHIX_DIR + "/phix.fna.gz"
PHIX_INDEX_PREFIX = PHIX_DIR + "/phix"       # bowtie2 writes phix.1.bt2 ... phix.rev.2.bt2

# ── 01.reads: what the read-QC steps pass to the assemblers ──
# One folder per sample and read type. The Illumina files are made by
# rules/illumina/10_reads.smk (illumina mode) and its copy rules/hybrid/10_reads.smk
# (hybrid mode); the long-read files by rules/longread/10_reads.smk and
# rules/hybrid/30_long_reads.smk.
READS_ILLUMINA_DIR = DIR_READS + "/{sample}/illumina"
READS_LONG_DIR     = DIR_READS + "/{sample}/long"
READS_SELECTED_DIR = DIR_READS + "/{sample}/selected"
# Illumina reads, in the order the steps run. Rule map_phix maps the raw pairs to the
# PhiX genome and keeps the pairs that do not match (NOPHIX_R1/R2, temporary). Rule
# trim_adapters (fastp) reads those, cuts adapters and low-quality ends, and writes the
# trimmed pairs plus fastp's two reports. The trimmed pairs go to SPAdes (rule
# illumina_assembly) and to the mapping that gives the contamination screen its coverage
# track: in illumina mode rule map_contigs (10_decontam.smk), in hybrid mode the guide
# rules in rules/hybrid/25_guide.smk.
NOPHIX_R1  = READS_ILLUMINA_DIR + "/{sample}_nophix_R1.fastq"
NOPHIX_R2  = READS_ILLUMINA_DIR + "/{sample}_nophix_R2.fastq"
TRIM_R1    = READS_ILLUMINA_DIR + "/{sample}_trim_R1.fastq"      # PhiX-free, adapter-trimmed (temporary)
TRIM_R2    = READS_ILLUMINA_DIR + "/{sample}_trim_R2.fastq"
FASTP_JSON = READS_ILLUMINA_DIR + "/{sample}_fastp.json"         # fastp's numbers, read by MultiQC
FASTP_HTML = READS_ILLUMINA_DIR + "/{sample}_fastp.html"         # the same as a page for people
# Long reads. Rule filter_long_reads (Filtlong) writes FILTERED_LONG, which the
# assembler (rule long_read_assembly), the minimap2 mapping of the draft (rule
# map_contigs) and Medaka (rule long_read_consensus) read. NanoPlot draws read-length
# and quality plots, with a NanoStats table of the numbers, once for the raw reads
# (rule raw_long_read_qc) and once for the filtered reads (rule filtered_long_read_qc);
# each run gets a folder of its own, and MultiQC reads both.
FILTERED_LONG         = READS_LONG_DIR + "/{sample}_filtered.fastq"   # Filtlong output (kept)
NANOPLOT_RAW_DIR      = READS_LONG_DIR + "/raw_qc"
NANOPLOT_FILTERED_DIR = READS_LONG_DIR + "/filtered_qc"
# ONT reads only, unless parameters.<mode>.ont_adapter is off (rules/shared/05_ont_adapters.smk).
# Rule find_ont_adapter reads the raw reads and decides which sequence (the ONT adapter, with the
# barcode block that follows it) fastplong removes; rule trim_ont_adapters writes the trimmed
# reads, which Filtlong reads instead of the raw ones (filtlong_input, section 13); rule
# check_ont_trimming counts the adapters left. The sequence file also goes to the contig adapter
# screen (rule blast_adapters), the choice and check files to its report (rule screen_adapters).
READS_ADAPTER_DIR     = READS_LONG_DIR + "/adapters"
ONT_ADAPTER_SCAN      = READS_ADAPTER_DIR + "/{sample}_adapter_scan.tsv"       # reads with the adapter, per 10,000 reads
ONT_ADAPTER_CHOICE    = READS_ADAPTER_DIR + "/{sample}_adapter_choice.tsv"     # the decision and the numbers behind it
ONT_FASTPLONG_OPTIONS = READS_ADAPTER_DIR + "/{sample}_fastplong_options.txt"  # one line: -s <sequence> ..., or -A
ONT_ADAPTER_SEQUENCE  = READS_ADAPTER_DIR + "/{sample}_adapter_sequence.fasta" # the trimmed sequence, or an empty file
TRIMMED_LONG          = READS_LONG_DIR + "/{sample}_trimmed.fastq.gz"         # fastplong output (temporary)
FASTPLONG_JSON        = READS_ADAPTER_DIR + "/{sample}_fastplong.json"         # fastplong's numbers
FASTPLONG_HTML        = READS_ADAPTER_DIR + "/{sample}_fastplong.html"         # the same as a page for people
ONT_ADAPTER_CHECK     = READS_ADAPTER_DIR + "/{sample}_adapter_check.tsv"      # reads still carrying the adapter
# hybrid only: the Illumina read pairs that map to the clean Illumina guide
# assembly (rule map_sel_contigs in rules/hybrid/25_guide.smk writes them). Filtlong
# uses them to score the long reads, Polypolish uses them to correct the final
# assembly.
SEL_R1 = READS_SELECTED_DIR + "/{sample}_selected_R1.fastq"
SEL_R2 = READS_SELECTED_DIR + "/{sample}_selected_R2.fastq"

# ── 02.assembly: the per-sample work area ──
ASSEMBLY_DIR = DIR_ASSEMBLY + "/{sample}"

# The three files a finished sample is judged by. All three sit directly in the
# sample's assembly directory.
#   CONTIGS_ALL   : the complete decontaminated, polished assembly, mitochondrion
#                   included. Written by rule finalize_assembly; read by rule
#                   separate_mitochondria.
#   FINAL_CONTIGS : the NUCLEAR genome (CONTIGS_ALL without the mitochondrion),
#                   written by rule separate_mitochondria. QUAST, BUSCO, ITS
#                   extraction and funannotate read this file and nothing else about
#                   the assembly. A fungal mitochondrion uses a different genetic
#                   code (code 4) and must not go through the nuclear gene predictor.
#   MITO_CONTIGS  : the confirmed mitochondrial contig(s); an empty file when none
#                   was confirmed. A delivered result, read by the contig summary (for
#                   the contig names; section 9h) and, when it is on, by the AMR screen
#                   (rule amr_chroquetas, for cytochrome b; section 9i).
# Without the screen (contigs mode, screen: false) only FINAL_CONTIGS exists, and
# the contig filter (rule filter_contigs) writes it directly.
CONTIGS_ALL   = ASSEMBLY_DIR + "/contigs_all.fasta"
FINAL_CONTIGS = ASSEMBLY_DIR + "/contigs_final.fasta"
MITO_CONTIGS  = ASSEMBLY_DIR + "/contigs_mitochondrial.fasta"

# DRAFT_CONTIGS is what goes INTO the contamination screen (10_decontam.smk and the
# mitochondrion screen read it). The mode's front end declares this exact path as
# the output of whichever rule makes it: rule filter_contigs in contigs and illumina
# mode, rule long_read_assembly in longread and hybrid mode.
#
# CAUTION: contigs_filt.fasta sits NEXT TO the spades/ directory, not inside it.
# The SPAdes rule (illumina_assembly) declares spades/ as a directory() output, and
# Snakemake deletes a directory() output before it re-runs that rule, so a file
# inside spades/ would be destroyed by any re-run of the assembler. (In v1 the whole
# 02.assembly/{sample}/ folder was that output.)
SPADES_DIR   = ASSEMBLY_DIR + "/spades"
HIFIASM_DIR  = ASSEMBLY_DIR + "/hifiasm"   # hifiasm's own files (graphs, .bin indexes), nothing else inside
FLYE_DIR     = ASSEMBLY_DIR + "/flye"      # Flye's own files
if MODE == "illumina":
    DRAFT_CONTIGS = ASSEMBLY_DIR + "/contigs_filt.fasta"    # SPAdes contigs, length/coverage filtered
elif HAS_LONG_READS:
    DRAFT_CONTIGS = ASSEMBLY_DIR + "/assembly.fasta"        # the long-read assembly (hifiasm or Flye)
elif CONTIGS_SCREEN:
    DRAFT_CONTIGS = ASSEMBLY_DIR + "/contigs_filt.fasta"    # the supplied FASTA after the contig filter
else:
    DRAFT_CONTIGS = FINAL_CONTIGS                           # contigs mode, no screen: the filter writes the final file directly

# Two small tables next to assembly.fasta (longread and hybrid only), both made
# right after the assembler has run:
#   ASSEMBLY_EVIDENCE   : per-contig read depth and circularity, which the
#                         mitochondrion screen (rule screen_mitochondria) reads. Rule
#                         assembly_evidence makes it by running mitochondria.py
#                         evidence-flye or evidence-hifiasm on the assembler's own
#                         output (Flye's assembly_info.txt, or hifiasm's primary graph).
#   ASSEMBLY_PROVENANCE : which tool, version and flags made the assembly, written
#                         by rule long_read_assembly itself.
ASSEMBLY_EVIDENCE   = ASSEMBLY_DIR + "/assembly_evidence.tsv"
ASSEMBLY_PROVENANCE = ASSEMBLY_DIR + "/assembly_provenance.tsv"

# The contamination screen (BLAST + BlobTools + the taxonomy selector, all in
# 10_decontam.smk) keeps all its working files in one directory. BLASTOUT, written
# by rule blast_contigs, is read by rule blob_json (BlobTools) and by the
# mitochondrion screen (rule screen_mitochondria). DECONTAM_CONTIGS is what the
# selector (rule select_contigs) writes: the kept contigs, which the next step of
# the mode's assembly chain reads (finalize_assembly, Medaka or Polypolish, by mode
# and sample; the table at the top of 10_decontam.smk lists each case). The other
# three files the selector writes are records for the reader.
DECONTAM_DIR     = ASSEMBLY_DIR + "/contaminants"
DECONTAM_CONTIGS = DECONTAM_DIR + "/contigs_decontam.fasta"
BLASTOUT         = DECONTAM_DIR + "/blastout"                   # the BLAST hits of every contig (tabular)
COMPOSITION      = DECONTAM_DIR + "/{sample}_composition.txt"   # genus share of the assembly
CONTIG_LIST      = DECONTAM_DIR + "/contigs.list"               # the IDs of the kept contigs
CONTIG_DECISIONS = DECONTAM_DIR + "/contig_taxonomy_decisions.tsv"  # the audit table: every contig, kept or dropped, and why

# The read alignment that gives BlobTools its coverage track (rule map_contigs makes
# it, rule blob_json reads it). In contigs mode map_contigs marks it temporary
# (10_decontam.smk says why). In the read-based modes it is kept: Qualimap (rule
# map_evaluation in 20_qc.smk) reads it too, and a re-run of BlobTools after a BLAST
# change should not have to map 40 million reads again.
# In illumina mode the reads are mapped with bowtie2, which needs an index of the
# draft first: rule index_contigs writes it next to the BAM under the prefix
# DECONTAM_INDEX_PREFIX (bowtie2 adds the six .bt2 endings itself), and rule
# map_contigs reads it.
# CAUTION: BlobTools names its coverage file after the BAM's file name, so BLOB_COV
# is DERIVED from DECONTAM_BAM instead of being typed out again. Renaming the BAM
# would otherwise silently break the declared output of rule blob_json.
# BLOB_PREFIX and BLOB_TABLE_PREFIX are the prefixes handed to BlobTools, which adds
# the endings itself: blob_json writes BLOB_JSON and BLOB_COV (both temporary), and
# blob_table writes BLOB_TABLE, the per-contig taxonomy that select_contigs reads.
DECONTAM_INDEX_PREFIX = DECONTAM_DIR + "/contigs"
DECONTAM_BAM      = DECONTAM_DIR + "/{sample}_map.bam"
BLOB_PREFIX       = DECONTAM_DIR + "/blob"
BLOB_JSON         = BLOB_PREFIX + ".blobDB.json"
BLOB_COV          = BLOB_PREFIX + "." + os.path.basename(DECONTAM_BAM) + ".cov"
BLOB_TABLE_PREFIX = DECONTAM_DIR + "/bestscore"
BLOB_TABLE        = BLOB_TABLE_PREFIX + ".blob.blobDB.table.txt"

# Mitochondrion screen (rule screen_mitochondria writes both files): which contigs
# are mitochondrial (the report, for people), and which IDs the taxonomy selector
# must protect from being discarded (select_contigs --always-keep). Rule
# separate_mitochondria reads the same ID list to split the finished genome.
MITO_DIR           = ASSEMBLY_DIR + "/mitochondria"
MITO_SCREEN_TSV    = MITO_DIR + "/mitochondrial_screen.tsv"
CONFIRMED_MITO_IDS = MITO_DIR + "/confirmed_mitochondrial.ids"

# Polishing, in the order it happens (longread and hybrid only). Medaka corrects ONT
# assemblies with the long reads (rule long_read_consensus, which writes its files
# into MEDAKA_DIR, and its alignment index of the draft next to the draft); Polypolish corrects hybrid assemblies with the Illumina reads
# (rule short_read_correction, into POLYPOLISH_DIR). Medaka's model is checked (or
# worked out from the read headers) BEFORE the assembler starts, by rule
# check_medaka_model in 12_medaka_check.smk, so a wrong model name fails in seconds;
# the result is kept in a small file (MEDAKA_MODEL_RESOLVED, see medaka_gate in
# section 13).
MEDAKA_MODEL_RESOLVED = ASSEMBLY_DIR + "/{sample}_medaka_model.txt"
MEDAKA_DIR            = ASSEMBLY_DIR + "/medaka"
MEDAKA_CONSENSUS      = MEDAKA_DIR + "/consensus.fasta"
POLYPOLISH_DIR        = ASSEMBLY_DIR + "/polypolish"
POLYPOLISH_CONTIGS    = POLYPOLISH_DIR + "/polished.fasta"

# hybrid only: the Illumina "guide" assembly (rules/hybrid/20_assembly.smk and
# rules/hybrid/25_guide.smk). SPAdes assembles the short reads into GUIDE_SPADES_DIR;
# the length-and-coverage filter writes GUIDE_DRAFT. The guide is then screened for
# contamination the same way the long-read draft is, by its own copies of the screen
# rules (index_guide_contigs, map_guide_contigs, blast_guide_contigs, blob_guide_json,
# blob_guide_table, screen_guide_mitochondria, select_guide_contigs), which write under
# GUIDE_DIR. The mitochondrion is protected here too, so that selection by genus can
# never throw its reads away. The read pairs that map to the screened guide are then
# selected for Filtlong and Polypolish (SEL_R1, SEL_R2 above), and QUAST compares
# GUIDE_DECONTAM with the final genome (rule genome_assembly_evaluation).
GUIDE_DIR             = ASSEMBLY_DIR + "/guide"
GUIDE_SPADES_DIR      = GUIDE_DIR + "/spades"
GUIDE_DRAFT           = GUIDE_DIR + "/contigs_filt.fasta"
GUIDE_DECONTAM_DIR    = GUIDE_DIR + "/contaminants"
GUIDE_BAM             = GUIDE_DECONTAM_DIR + "/{sample}_map.bam"        # the Illumina reads mapped to the guide
GUIDE_BLASTOUT        = GUIDE_DECONTAM_DIR + "/blastout"                # the BLAST hits of every guide contig
GUIDE_DECONTAM        = GUIDE_DECONTAM_DIR + "/contigs_decontam.fasta"  # the kept guide contigs
GUIDE_MITO_SCREEN_TSV = GUIDE_DIR + "/mitochondria/mitochondrial_screen.tsv"   # the guide's mitochondrion screen report
GUIDE_MITO_IDS        = GUIDE_DIR + "/mitochondria/confirmed_mitochondrial.ids"  # and its ID list
GUIDE_QUALIMAP_DIR    = GUIDE_DIR + "/eval/qualimap"                    # Qualimap on GUIDE_BAM (rule map_guide_evaluation)

# Assembly QC, per sample (20_qc.smk; MultiQC reads the results). QUAST
# (rule genome_assembly_evaluation) and BUSCO (rule completeness) each write a whole
# directory of files, and the rule declares that directory as its output, so a
# re-run starts clean, as in v1. QUALIMAP_DIR is the Qualimap report of the
# read-based modes (rule map_evaluation: the reads aligned to the draft, DECONTAM_BAM).
QUAST_DIR    = ASSEMBLY_DIR + "/eval/quast"
BUSCO_DIR    = ASSEMBLY_DIR + "/eval/busco"
QUALIMAP_DIR = ASSEMBLY_DIR + "/eval/qualimap"

# The contig summary, per sample (25_contig_summary.smk), in two steps. Rule telomere_units
# works out the telomere repeat unit of the assembly from the draft contigs with tidk and
# writes three files: every candidate unit with its score (TELOMERE_CANDIDATES), the chosen
# unit alone (TELOMERE_UNIT_FILE, one line: the unit, or none) and tidk's telomere call at
# both ends of every contig for that unit (TELOMERE_ENDS). Rule contig_summary reads the
# last two and joins them with the rest of the run's evidence into one table that says
# what each contig of the draft assembly is (a chromosome, a chromosome arm, the
# mitochondrion, a piece of the rDNA array, ...) with the evidence behind each verdict
# (CONTIG_SUMMARY), plus a copy in the form MultiQC reads (CONTIG_SUMMARY_MQC; rule
# multiqc lists it as an input). The files sit next to the QUAST and BUSCO folders,
# because the table is a judgement of the assembly, although the jobs run after the
# annotation, whose gene counts the table uses (section 9h has the settings and the
# inputs the summary has in each mode). The candidates table is for people; no rule
# reads it.
TELOMERE_CANDIDATES = ASSEMBLY_DIR + "/eval/telomere_candidates.tsv"   # every candidate unit, scored; the chosen one marked
TELOMERE_UNIT_FILE  = ASSEMBLY_DIR + "/eval/telomere_unit.txt"         # the chosen unit (canonical form), or none
TELOMERE_ENDS       = ASSEMBLY_DIR + "/eval/telomere_ends.tsv"         # tidk's call at both ends of every contig, for that unit
CONTIG_SUMMARY      = ASSEMBLY_DIR + "/eval/contig_summary.tsv"
CONTIG_SUMMARY_MQC  = ASSEMBLY_DIR + "/eval/contig_summary_mqc.tsv"    # the same table, for MultiQC

# The contig adapter screen, per sample and in every mode (25_contig_summary.smk). Rule
# blast_adapters searches the draft contigs for the adapters shipped in workflow/data/adapters/,
# the Illumina records of UniVec (UNIVEC_ILLUMINA) and, for a trimmed ONT sample, the sequence
# that was trimmed from its reads (ONT_ADAPTER_SEQUENCE); it writes the list it searched for
# (ADAPTER_QUERY, temporary) and BLAST's hits (ADAPTER_BLAST). Rule screen_adapters turns the
# hits into adapter loci near contig ends and inside contigs (ADAPTER_SCREEN; rule
# contig_summary adds them to its table as the adapter column) and writes one row per sample
# for the Adapters table of the MultiQC report (ADAPTER_SCREEN_MQC).
ADAPTER_QUERY       = ASSEMBLY_DIR + "/eval/adapter_query.fasta"       # the adapters searched for (temporary)
ADAPTER_BLAST       = ASSEMBLY_DIR + "/eval/adapter_blast.tsv"         # BLAST's hits on the draft contigs
ADAPTER_SCREEN      = ASSEMBLY_DIR + "/eval/adapter_screen.tsv"        # one row per adapter locus
ADAPTER_SCREEN_MQC  = ASSEMBLY_DIR + "/eval/adapter_screen_mqc.tsv"    # one row per sample, for MultiQC

# The names QUAST prints for the genome(s) it is given, one per input FASTA, in the
# same order. Every report that reads QUAST (MultiQC) finds the sample through these
# labels. In hybrid mode one QUAST call compares the Illumina "guide" assembly with
# the final genome, so there are two. The strings still contain the {sample}
# wildcard; quast_labels() in 20_qc.smk fills it in for each job.
if IS_HYBRID:
    QC_LABELS = ["{sample}_guide", "{sample}"]
else:
    QC_LABELS = ["{sample}"]

# ── 03.taxonomy: ITS extraction and classification (30_taxonomy.smk) ──
# Rule its_windows writes the two window files, rule its_extraction (ITSx) the four
# ITS files, and rule its_classification the copies table and the two taxonomy files.
TAXONOMY_DIR = DIR_TAXONOMY + "/{sample}"
ITS_WINDOWS  = TAXONOMY_DIR + "/its_windows.fasta"             # one window per rDNA copy (docs/methods/its-long-read-contigs.md)
ITS_WINDOW_TABLE = TAXONOMY_DIR + "/its_windows.tsv"           # the same windows as a table: contig, start, end, strand, anchor hits
ITS_COPIES   = TAXONOMY_DIR + "/{sample}_ITS_copies.txt"       # copies found, and any record the length guard rejected

# What ITSx writes from the windows (its -o prefix is TAXONOMY_DIR/{sample}): the whole
# ITS1-5.8S-ITS2 region and its two spacers, one record per rDNA copy, plus ITSx's own
# summary. SINTAX then classifies the distinct ITS1 and ITS2 sequences.
ITS_FULL      = TAXONOMY_DIR + "/{sample}.full.fasta"
ITS_ITS1      = TAXONOMY_DIR + "/{sample}.ITS1.fasta"
ITS_ITS2      = TAXONOMY_DIR + "/{sample}.ITS2.fasta"
ITS_SUMMARY   = TAXONOMY_DIR + "/{sample}.summary.txt"
ITS1_TAXONOMY = TAXONOMY_DIR + "/{sample}_ITS1_taxonomy.txt"
ITS2_TAXONOMY = TAXONOMY_DIR + "/{sample}_ITS2_taxonomy.txt"

# ── 04.annotation: the per-sample annotation work area (40_annotation.smk) ──
# The chain, in order: funannotate clean + sort (CONTIGS_CLEAN, CONTIGS_SORT) -> repeat
# masking (CONTIGS_MASK) -> gene prediction (FUNANNOTATE_DIR) -> InterProScan, eggNOG,
# antiSMASH and, if switched on, SignalP 6, each reading the gene prediction (the
# predicted proteins; antiSMASH reads the GenBank file) -> funannotate annotate,
# which merges all of them into the final annotation (ANNOTATION_COMPLETE) ->
# cleanup of the scratch files (CLEANUP_COMPLETE). The two *_COMPLETE files are
# marker files: empty files whose only job is to show that a step finished.
ANNOTATION_DIR   = DIR_ANNOTATION + "/{sample}"
CONTIGS_CLEAN    = ANNOTATION_DIR + "/contigs_clean.fasta"
CONTIGS_SORT     = ANNOTATION_DIR + "/contigs_sort.fasta"           # the annotated genome: contigs renamed scaffold_N
CONTIGS_MASK     = ANNOTATION_DIR + "/contigs_mask.fasta"           # CONTIGS_SORT with repeats in lowercase
CONTIGS_REMOVED  = ANNOTATION_DIR + "/contigs_removed_by_clean.txt" # what funannotate clean dropped, so nothing disappears silently
REPEATMASK_DIR   = ANNOTATION_DIR + "/repeatmasking"
FUNANNOTATE_DIR  = ANNOTATION_DIR + "/funannotate"
IPRSCAN_DIR      = ANNOTATION_DIR + "/iprscan"
IPRSCAN_XML      = IPRSCAN_DIR + "/interproscan_out.xml"
EGGNOG_DIR       = ANNOTATION_DIR + "/eggnog"
ANTISMASH_DIR    = ANNOTATION_DIR + "/antismash"
SIGNALP6_RESULTS = ANNOTATION_DIR + "/signalp6/prediction_results.txt"
# Three small files whose only job is to show that a step finished:
#   ANNOTATION_COMPLETE : the marker of rule funannotate_annotation. funannotate
#                         annotate writes its results into FUNANNOTATE_DIR, which
#                         is already the directory() output of rule
#                         funannotate_prediction, and Snakemake refuses an output
#                         inside another rule's directory() output; so this marker
#                         is the one output the rule can declare (the header of
#                         funannotate_annotation in 40_annotation.smk has the
#                         details). Rule cleanup waits for it.
#   GENEMARK_CHECK      : written by rule genemark_check (41_genemark.smk), a
#                         GeneMark-ES test run on a 200 kb slice, so a dead licence
#                         key shows up in a minute, not well into the gene
#                         prediction (see genemark_gate in section 12).
#   CLEANUP_COMPLETE    : the marker of rule cleanup, the last step of a sample;
#                         rule all asks for it.
ANNOTATION_COMPLETE = ANNOTATION_DIR + "/.annotation_complete.txt"
GENEMARK_CHECK      = ANNOTATION_DIR + "/genemark_check.txt"
CLEANUP_COMPLETE    = ANNOTATION_DIR + "/.cleanup_complete.txt"
# What AUGUSTUS does for this sample (a stored model, or training on the genome), as
# checked by rule augustus_model_check in 40_annotation.smk; funannotate_prediction
# waits for it and copies it into its log (section 9d explains the choice).
AUGUSTUS_MODEL_RECORD = ANNOTATION_DIR + "/augustus_model.txt"
# Tiberius's gene models of this sample (42_tiberius.smk; only made when Tiberius is
# on): TIBERIUS_GTF as Tiberius writes them (rule tiberius_prediction), and
# TIBERIUS_GFF3, the same models converted into the GFF3 that funannotate predict
# hands to EVidenceModeler (rule tiberius_gff3; funannotate_prediction reads it).
TIBERIUS_GTF  = ANNOTATION_DIR + "/tiberius/tiberius.gtf"
TIBERIUS_GFF3 = ANNOTATION_DIR + "/tiberius/tiberius.evm.gff3"

# ── 05.amr: the antifungal resistance screen (50_amr.smk; only when it is on) ──
# Rule amr_chroquetas runs ChroQueTas on the sample's nuclear genome and, for a species
# whose database folder holds the mitochondrial target cytochrome b, once more on its
# mitochondrial genome (section 9i says why). Each run has its own folder inside
# AMR_CHROQUETAS_DIR (nuclear/ and mitochondrial/, ChroQueTas's own files, untouched; the
# rule declares the whole folder as its output). From those files it writes three tables:
#   AMR_POSITIONS   : every position of a target protein that was checked, with the amino
#                     acid of the reference and of this genome, the call (FungAMR mutation,
#                     new mutation, no mutation), a verdict on resistance, the drugs with
#                     evidence for and against resistance, and the changes known to cause
#                     resistance at that position.
#   AMR_SUMMARY     : one row per target protein of the species: whether it was evaluated,
#                     its verdict, the mutations found, and the drugs with evidence for and
#                     against resistance. Rule all asks for it.
#   AMR_SUMMARY_MQC : the same table in the form MultiQC reads; rule multiqc lists it as an
#                     input.
AMR_DIR            = DIR_AMR + "/{sample}"
AMR_CHROQUETAS_DIR = AMR_DIR + "/chroquetas"
AMR_POSITIONS      = AMR_DIR + "/amr_positions.tsv"
AMR_SUMMARY        = AMR_DIR + "/amr_summary.tsv"
AMR_SUMMARY_MQC    = AMR_DIR + "/amr_summary_mqc.tsv"

# ── 06.report ──
MULTIQC_REPORT = DIR_REPORT + "/multiqc_report.html"   # the run's MultiQC page (rule multiqc in 90_report.smk)


# ═════════════ 6. Scripts and shared constants ═════════════
# Every helper script lives in a directory named after the rule file that runs it
# (workflow/scripts/README.md has the table), and its path is declared here ONCE.
# Moving a script means editing one line, not hunting through the rule files.
# workflow.basedir is the folder that holds the Snakefile (workflow/), as an
# absolute path, so these paths work wherever snakemake is started.
# The rule that runs each script:
#   ONT_ADAPTER_SCRIPT      : rules find_ont_adapter and check_ont_trimming
#                             (05_ont_adapters.smk), and the start-up check of a sequence
#                             given in parameters.<mode>.ont_adapter (section 11)
#   SELECT_TAXONOMY_SCRIPT  : rule select_contigs
#   MEDAKA_CHECK_SCRIPT     : rule check_medaka_model (12_medaka_check.smk)
#   MITOCHONDRIA_SCRIPT     : rules screen_mitochondria and separate_mitochondria
#   TELOMERE_UNITS_SCRIPT   : rule telomere_units (25_contig_summary.smk)
#   CONTIG_SUMMARY_SCRIPT   : rule contig_summary (25_contig_summary.smk)
#   ADAPTER_SCREEN_SCRIPT   : rule screen_adapters (25_contig_summary.smk)
#   UNITE_IDENTIFIED_SCRIPT : rule unite_identified
#   ITS_WINDOWS_SCRIPT      : rule its_windows
#   ITS_DISTINCT_SCRIPT     : rule its_classification
#   NORMALIZE_DB_SCRIPT     : rule normalize_funannotate_db
#   AUGUSTUS_CHECK_SCRIPT   : rule augustus_model_check, and the start-up check of the
#                             AUGUSTUS models (sections 9d and 12)
#   SIGNALP6_INSTALL_SCRIPT : rule install_signalp6
#   CHROQUETAS_AMR_SCRIPT   : rule amr_chroquetas (50_amr.smk)
WORKFLOW_DIR = workflow.basedir
SCRIPTS_DIR  = os.path.join(WORKFLOW_DIR, "scripts")
ONT_ADAPTER_SCRIPT      = os.path.join(SCRIPTS_DIR, "05_ont_adapters", "find_ont_adapter.py")
SELECT_TAXONOMY_SCRIPT  = os.path.join(SCRIPTS_DIR, "10_decontam", "select_contigs_by_taxonomy.py")
MEDAKA_CHECK_SCRIPT     = os.path.join(SCRIPTS_DIR, "12_medaka_check", "medaka_model_check.py")
MITOCHONDRIA_SCRIPT     = os.path.join(SCRIPTS_DIR, "15_mitochondria", "mitochondria.py")
TELOMERE_UNITS_SCRIPT   = os.path.join(SCRIPTS_DIR, "25_contig_summary", "telomere_units.py")
CONTIG_SUMMARY_SCRIPT   = os.path.join(SCRIPTS_DIR, "25_contig_summary", "contig_summary.py")
ADAPTER_SCREEN_SCRIPT   = os.path.join(SCRIPTS_DIR, "25_contig_summary", "adapter_screen.py")
# The adapter sequences shipped with FunFlux for the contig adapter screen (rule blast_adapters;
# workflow/data/adapters/README.md says where each comes from): ONT's adapters and barcodes from
# its Chemistry Technical Document, and the PacBio and MGI adapters of NCBI UniVec.
# ONT_BARCODE_LIST, the first of them, also goes to rule find_ont_adapter (05_ont_adapters.smk),
# which names the barcodes it finds in the reads from it and trims every barcode of a file that
# holds several.
ONT_BARCODE_LIST = os.path.join(WORKFLOW_DIR, "data", "adapters", "ont_chtd.fasta")
ADAPTER_LIST_FILES = [
    ONT_BARCODE_LIST,
    os.path.join(WORKFLOW_DIR, "data", "adapters", "univec_pacbio_mgi.fasta"),
]
UNITE_IDENTIFIED_SCRIPT = os.path.join(SCRIPTS_DIR, "30_taxonomy", "unite_identified.py")
ITS_WINDOWS_SCRIPT      = os.path.join(SCRIPTS_DIR, "30_taxonomy", "its_windows.py")
ITS_DISTINCT_SCRIPT     = os.path.join(SCRIPTS_DIR, "30_taxonomy", "its_distinct.py")
NORMALIZE_DB_SCRIPT     = os.path.join(SCRIPTS_DIR, "40_annotation", "normalize_funannotate_db_paths.py")
AUGUSTUS_CHECK_SCRIPT   = os.path.join(SCRIPTS_DIR, "40_annotation", "check_augustus_model.py")
SIGNALP6_INSTALL_SCRIPT = os.path.join(SCRIPTS_DIR, "45_signalp6", "install_signalp6.py")
CHROQUETAS_AMR_SCRIPT   = os.path.join(SCRIPTS_DIR, "50_amr", "chroquetas_amr.py")

# Sample names become file names, wildcards, report labels and the locus-tag prefix
# funannotate puts on every gene. These characters break paths or make wildcard
# matching ambiguous, so a sample whose name contains any of them is rejected
# (_check_sample_names, section 12). Underscore is ALLOWED: it is common in real
# sample names (the v1 fungal workflows rejected it; BacFlux v2 does not).
BAD_CHARS = set("*#@%^/! ?&:;|<>")

# Accepted input file extensions, checked by validate_extensions (section 7). Read
# files may be gzip-compressed; FASTA files may not (a .fasta.gz is refused).
FASTQ_SUFFIXES = ("fastq", "fq", "fastq.gz", "fq.gz")
FASTA_SUFFIXES = ("fasta", "fa", "fna")

# A valid NCBI locus-tag prefix is 3 to 12 letters or digits and starts with a letter.
# funannotate uses the sample name as that prefix (funannotate predict --name), so a
# sample name that does not fit gets a NOTE (it still runs, but the genome cannot be
# submitted to GenBank without renaming the genes). Checked in _check_sample_names.
LOCUS_TAG_PATTERN = r"[A-Za-z][A-Za-z0-9]{2,11}"

# Small awk programs shared by several rules, written once so the fragile
# one-liners are not copied around. Used as `awk {NAME:q}` (the :q lets Snakemake
# quote them for the shell). Today rule filter_contigs (contigs front end) uses all
# three. The r"""...""" form keeps every backslash and quote exactly as awk needs it.
#   FASTA_LIN_CMD  : turn a wrapped multi-line FASTA into one line per record.
#   FASTA_SEL_CMD  : keep SPAdes contigs with coverage >= 2.0 and length >= 500
#                    (fields 6 and 4 of the "_"-split SPAdes header
#                    NODE_<n>_length_<L>_cov_<C>); run with awk -F'_' on the
#                    one-line-per-record form.
#   FASTA_HEAD_CMD : cut each FASTA header at its first space or tab (with
#                    FS=" ", awk splits a line at runs of spaces and tabs).
#                    Polypolish appends " polypolish" to headers, and a header
#                    kept whole breaks BLAST, BlobTools and the mitochondrion
#                    screen, which identify a contig by the first word of its
#                    header.
FASTA_LIN_CMD  = r"""{if(NR==1) {printf "%s\n", $0} else {if(/^>/) {printf "\n%s\n", $0} else {printf $0}}}"""
FASTA_SEL_CMD  = r"""{if(/^>/ && $6>=2.0 && $4>=500) {printf "%s\n", $0; getline; print}}"""
FASTA_HEAD_CMD = r"""BEGIN{FS=" "} /^>/{print $1; next} {print}"""

# How contigs mode recognises an assembly made by SPAdes: its first header looks
# like ">NODE_12_length_45678_cov_23.4". Only then are length and coverage read out of
# the header and filtered; for any other assembler those numbers do not exist, and
# inventing them would silently delete real sequence. Strict on purpose: a header such
# as ">scaffold_1 length_12345" must not be mistaken for a SPAdes one.
SPADES_HEADER_REGEX = r"^>NODE_[0-9]+_length_[0-9]+_cov_"


# ═══════════════════════ 7. Helpers ═══════════════════════
# Functions that the sections below need.

def validate_extensions(where, expected, extensions, allowed_suffixes):
    # Check one class of discovered input files, from one input folder. All files of the
    # class must share ONE file extension, otherwise the {sample}.{extn} paths become
    # ambiguous.
    #   where            : the config key to name in a message, such as "input.illumina_dir"
    #   expected         : an example of a correct file name, for the "none found" message
    #   extensions       : the extensions glob_wildcards found
    #   allowed_suffixes : the acceptable endings (FASTQ_SUFFIXES or FASTA_SUFFIXES)
    # Returns the single validated extension, or stops with a message. Called by the
    # sample discovery of every mode (section 12).
    if not extensions:
        sys.exit(f"[FunFlux] {where}: no input files found. Expected names like {expected}.")
    distinct = sorted(set(extensions))
    if len(distinct) != 1:
        sys.exit(
            f"[FunFlux] {where}: more than one file extension in this folder "
            f"({', '.join(distinct)}). All files in it must use the same one."
        )
    extension = distinct[0]
    if not extension.endswith(allowed_suffixes):
        sys.exit(
            f"[FunFlux] {where}: file extension '{extension}' is not recognised "
            f"(expected one of: {', '.join(allowed_suffixes)})."
        )
    return extension


def medaka_disabled(value):
    # True only when the user explicitly turned Medaka OFF (value is the raw
    # medaka_model setting). An empty value does NOT count as disabled: a missing,
    # empty or "auto" medaka_model (an unquoted null too, which YAML hands over as
    # None) means "work the model out from the read headers", not "skip polishing".
    # The words none and null typed as text (none, or null in quotes) do mean off,
    # as they do for the other settings (EMPTY_WORDS in section 1b).
    # Used in section 11 to set USE_MEDAKA.
    if isinstance(value, bool):
        return value is False
    if value is None:
        return False
    return str(value).strip().lower() in {"false", "0", "no", "off", "none", "null"}


def has_explicit_medaka_model(value):
    # True only when value is a concrete Medaka model name. The words that mean
    # "work it out" (auto, true, yes, on, empty) and the words that mean "skip"
    # (false, no, off, none, null) are both False here. Used in section 11 (the
    # Flye read mode) and section 13 (the Medaka message).
    if isinstance(value, bool) or value is None:
        return False
    return str(value).strip().lower() not in {
        "", "true", "1", "yes", "on", "auto", "false", "0", "no", "off", "none", "null"
    }


def id_from_url(url):
    # A stable database ID from an archive URL: the file name without .tar.gz, .tgz
    # or its other extension. For the default link it turns
    #   https://zenodo.org/records/23265189/files/funannotate_db_1.8.17_2026-10-08.tar.gz?download=1
    # into funannotate_db_1.8.17_2026-10-08 (urlparse drops the "?download=1" part,
    # unquote turns a %20 back into a space). Used in section 10 to name the
    # funannotate database archive and its checksum file (rule funannotate_db).
    path = unquote(urlparse(url).path)
    filename = os.path.basename(path)
    if not filename:
        sys.exit(f"[FunFlux] links.funannotate_link: cannot read a file name from the URL: {url}")
    if filename.endswith(".tar.gz"):
        return filename[:-len(".tar.gz")]
    if filename.endswith(".tgz"):
        return filename[:-len(".tgz")]
    return os.path.splitext(filename)[0]


def zenodo_record_from_url(url):
    # The Zenodo record number from a standard https://zenodo.org/records/<id>/files/...
    # URL. The funannotate database snapshot is published on Zenodo together with a
    # .sha256 checksum file, and section 10 builds the checksum URL from this number.
    # For the default link the steps give:
    #   url_path : "/records/23265189/files/funannotate_db_1.8.17_2026-10-08.tar.gz"
    #   parts    : ["records", "23265189", "files", "funannotate_db_1.8.17_2026-10-08.tar.gz"]
    #   returned : "23265189"
    url_path = urlparse(url).path
    parts = url_path.strip("/").split("/")
    if len(parts) >= 3 and parts[0] == "records":
        return parts[1]
    sys.exit(
        "[FunFlux] links.funannotate_link: expected a Zenodo URL of the form "
        f"https://zenodo.org/records/<id>/files/<archive>, got: {url}"
    )


# ═════════════════ 8. Screening settings ═════════════════
# Everything that decides what the contamination and mitochondrion screen does
# (10_decontam.smk and 15_mitochondria.smk). These settings do not depend on the
# input, so they are read and checked in every run, even when the screen is off.
# The one exception is the check of the decontamination tables, which needs the
# samples and runs only with the screen (see "8a (continued)").

# ── 8a. Decontamination policy ──
# The taxonomy selector (select_contigs_by_taxonomy.py, run by rule select_contigs in
# 10_decontam.smk) takes plain text on its command line, so every setting is turned
# into text here, once.
def _as_decontam_text(value):
    # YAML may hand us None, a single value or a list; the selector wants one text:
    #   None           -> ""        (nothing set)
    #   list/tuple/set -> "A;B;C"   (semicolon-joined, blanks dropped)
    #   anything else  -> its stripped text
    # There is deliberately no branch for a mapping. Three of the fields
    # (include_genera_by_sample, exclude_genera_file, sample_overrides) are FILE
    # PATHS that the selector opens on disk, so they must arrive as plain paths; a
    # YAML mapping would turn into the text "{}" and be opened as a file of that name.
    # Called by _decontam_settings for each of the five text fields.
    if value is None:
        return ""
    if isinstance(value, (list, tuple, set)):
        names = []
        for item in value:
            name = str(item).strip()
            if name:
                names.append(name)
        return ";".join(names)
    return str(value).strip()


# Why FunFlux has no decontamination mode auto (removed on 2026-10-08; see
# docs/analysis/decontamination.md). auto kept the genus that BlobTools gave to
# the most contigs, and BLAST genus labels are unreliable for fungi with few
# relatives in NCBI nt. Many short contigs can get no hit at all, and these
# unplaced "no-hit" contigs can be the largest group although they hold little of
# the DNA: auto would then keep them and remove most of the genome. With no-hit
# left out of the vote, the winner can be a mislabel (a better-sequenced relative
# of the isolate), and auto would still remove the isolate's own contigs that
# BLAST labelled with other genera. The same labels make include with the
# isolate's true genus just as harmful: when BLAST gives that genus to only a few
# contigs, include with it alone removes most of the genome.
# So the advice at the end of the message is: start with off, read the genus
# BLAST gave each contig, then exclude the contaminant genera or include the
# genera BLAST gave the isolate's own contigs. Keep this advice the same as in
# the contig selector's own refusal of auto (select_contigs_by_taxonomy.py).
# A configuration that asks for auto, in parameters.decontamination.mode or in the
# mode column of the sample_overrides table, is stopped with this explanation:
_AUTO_MODE_REMOVED_TEXT = (
    "BLAST genus labels are unreliable for poorly sequenced fungi, so auto can remove most of "
    "the genome. Use off (the default), read contig_taxonomy_decisions.tsv, then exclude the "
    "contaminant genera or include the genera that BLAST assigned to the isolate's own contigs."
)


def _decontam_settings(default_mode):
    # Collect every decontamination choice into one dict that is passed unchanged to
    # the selector script, so all samples follow the same policy.
    #   Takes in : config["parameters"]["decontamination"], or, for a configuration
    #              written for v1, the old top-level parameters["genus"] key;
    #              default_mode is the mode used when nothing is set.
    #   Produces : a dict with a checked mode (off, include or exclude), five text
    #              fields and discard_no_hit as the text "true"/"false" (it is
    #              passed as a command-line argument, so it has to be text). Rule
    #              select_contigs hands each field to the selector as one
    #              command-line option.
    allowed = {
        "mode", "discard_no_hit", "include_genera", "include_genera_by_sample",
        "exclude_genera", "exclude_genera_file", "sample_overrides",
    }
    # A copy of the block, so that the changes below never touch config itself.
    settings = dict(_parameter_block("decontamination", allowed))

    # A v1 configuration named the target organism in one key, `genus`. Translate it
    # so such a configuration keeps working, and say so. This happens only when the
    # decontamination block is absent or empty:
    #   a genus filled in   -> mode include with that genus;
    #   a genus left empty  -> mode off. FunFlux 1.x then kept the genus with the
    #                          most contigs, as the removed mode auto did (see
    #                          above). FunFlux v2 makes no automatic choice, so every
    #                          contig is kept, and the NOTE says how to filter;
    #   no genus key at all -> default_mode.
    if not settings:
        legacy_genus = PARAMETERS.get("genus", None)
        if legacy_genus is not None and str(legacy_genus).strip():
            settings = {"mode": "include", "include_genera": legacy_genus}
            _note("parameters.genus is a v1 key; it was translated to parameters.decontamination "
                  "mode 'include'. Use parameters.decontamination in new configurations. include keeps "
                  "only the contigs BLAST assigned to that genus; for a fungus with few relatives in NCBI "
                  "nt that can be a small part of the genome, so read contig_taxonomy_decisions.tsv.")
        elif "genus" in PARAMETERS:
            settings = {"mode": "off"}
            _note("parameters.genus (a v1 key) is empty. FunFlux 1.x kept the genus with the most "
                  "contigs when parameters.genus was empty; FunFlux v2 has no automatic choice, so every "
                  "contig is kept (decontamination mode off). To filter, read contig_taxonomy_decisions.tsv, "
                  "then set parameters.decontamination: exclude the contaminant genera, or include the "
                  "genera that BLAST assigned to the isolate's own contigs.")
        else:
            settings = {"mode": default_mode}

    # Check the mode, then turn every other setting into the text the selector reads.
    # auto gets a message of its own, because users of FunFlux 1.x and of BacFlux
    # (which keeps auto) know it.
    mode = _mode_word(settings.get("mode"), default_mode)
    if mode == "auto":
        sys.exit(f"[FunFlux] parameters.decontamination.mode: 'auto' is not available in FunFlux. {_AUTO_MODE_REMOVED_TEXT}")
    if mode not in {"off", "include", "exclude"}:
        sys.exit(
            "[FunFlux] parameters.decontamination.mode: must be one of off | include | exclude "
            f"(got: {settings.get('mode')!r})"
        )
    settings["mode"] = mode
    settings["include_genera"]           = _as_decontam_text(settings.get("include_genera", ""))
    settings["include_genera_by_sample"] = _as_decontam_text(settings.get("include_genera_by_sample", ""))
    settings["exclude_genera"]           = _as_decontam_text(settings.get("exclude_genera", ""))
    settings["exclude_genera_file"]      = _as_decontam_text(settings.get("exclude_genera_file", ""))
    settings["sample_overrides"]         = _as_decontam_text(settings.get("sample_overrides", ""))
    settings["discard_no_hit"] = str(
        _config_bool("parameters.decontamination.discard_no_hit", settings.get("discard_no_hit"), False)
    ).lower()
    return settings


# The three modes, applied to the genus BlobTools gives each contig from its BLAST
# hits:
#   off     : keep every contig (the default);
#   include : keep only the contigs of the genera named. These must be the genera
#             BLAST gave the isolate's own contigs, which for a poorly sequenced
#             fungus are often not its true genus (see above);
#   exclude : remove only the contigs of the genera named (known contaminants).
# Whatever the mode, the confirmed mitochondrial contigs are kept. In exclude mode,
# discard_no_hit: true also removes the contigs BLAST could not place ("no-hit");
# include mode removes them anyway, because no-hit is not a genus on its list.
# select_contigs_by_taxonomy.py has the full rules. The fungal default is OFF:
# every contig is kept, and every contig still gets an audited verdict in
# contig_taxonomy_decisions.tsv. Fungal assemblies often need a human to read the
# BLAST/BlobTools result before anything is removed.
# DECONTAMINATION is read by the params of rule select_contigs, and by
# _check_decontam_tables below.
DECONTAMINATION = _decontam_settings("off")
print(
    "Decontamination mode: "
    f"{DECONTAMINATION['mode']} (discard_no_hit={DECONTAMINATION['discard_no_hit']})."
)


# ── 8a (continued). The decontamination tables, checked before any job starts ──
# parameters.decontamination may name three files:
#   include_genera_by_sample : tab-separated, columns sample and genus; each
#                              sample's own genera to keep, added to include_genera.
#   exclude_genera_file      : genera to remove, one per line, no header line.
#   sample_overrides         : tab-separated, columns sample and mode
#                              (include_genera, exclude_genera and discard_no_hit
#                              may follow); one row changes the policy of one sample.
# The selector opens these files only when rule select_contigs runs, after BLAST
# (often an hour or more into a run in contigs mode; in the read-based modes after
# the assembly as well), so a missing file or a table it cannot use would stop the
# run only then. _check_decontam_tables applies the selector's rules to them before
# any job starts. It needs the samples of this run, which are known only after
# sample discovery, so it is defined here and CALLED at the end of section 12. It
# runs only when the screen runs (SCREEN): without the screen there is no rule
# select_contigs, and the files are never read. A relative path is read from the
# folder Snakemake was started in, as the selector reads it (the workflow never
# changes folder; section 3).
#
# CAUTION: the four functions below repeat rules of select_contigs_by_taxonomy.py:
# read_lines_file, read_include_by_sample, read_overrides, split_genera, parse_bool,
# and the part of main() that applies a sample's override row. The selector applies
# the same rules again at run time, so the two must stay in step: change one,
# change the other.

def _split_genus_cell(value):
    # The selector's split_genera(): one cell, or one config value, may hold several
    # genera separated by commas, semicolons or tabs. Returns the names in the order
    # given, spaces trimmed, empty pieces dropped ([] for an empty cell, and for a
    # missing one, which arrives as None). Called by _check_decontam_tables.
    text = str(value or "").strip()
    text = text.replace(",", "\t")
    text = text.replace(";", "\t")
    genera = []
    for piece in text.split("\t"):
        name = piece.strip()
        if name:
            genera.append(name)
    return genera


def _stop_if_not_utf8(key, path):
    # The selector reads all three files as UTF-8 text and stops with a Python
    # error on any other encoding, for example a table that a spreadsheet program
    # saved as UTF-16 or in an old Windows encoding. Read the whole file once in
    # the same way, so that such a file stops the run now, with a plain message.
    #   key  : the config key, for the message
    #   path : a file that exists (the callers check that first)
    # Returns nothing when the file is UTF-8 text. Called by _read_decontam_table
    # and, for exclude_genera_file, by _check_decontam_tables.
    try:
        with open(path, "r", encoding="utf-8") as handle:
            handle.read()
    except UnicodeDecodeError:
        sys.exit(
            f"[FunFlux] parameters.decontamination.{key}: '{path}' is not UTF-8 text. "
            "Save it again as plain text with UTF-8 encoding."
        )


def _read_decontam_table(key, path, required_columns):
    # Read one of the two tables that have a header line (include_genera_by_sample,
    # sample_overrides) as the selector reads it: UTF-8 text, tab-separated, the
    # first line naming the columns.
    #   key              : the config key, for the messages
    #   path             : the table, as written in the config
    #   required_columns : the column names the header line must contain
    # Stops when the file does not exist, is not UTF-8 text, or lacks a required
    # column. Returns the rows, each a dict of column name -> cell text (None for a
    # cell the row does not have). Called twice by _check_decontam_tables.
    if not os.path.isfile(path):
        sys.exit(f"[FunFlux] parameters.decontamination.{key}: '{path}' is not an existing file.")
    _stop_if_not_utf8(key, path)
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        # fieldnames is the header line cut at its tabs (None for an empty file).
        columns = reader.fieldnames or []
        missing_columns = []
        for column in required_columns:
            if column not in columns:
                missing_columns.append(column)
        if missing_columns:
            sys.exit(
                f"[FunFlux] parameters.decontamination.{key}: '{path}' must start with a tab-separated "
                f"header line that has the columns {' and '.join(required_columns)} "
                f"(missing: {', '.join(missing_columns)})."
            )
        rows = list(reader)
    return rows


def _check_decontam_tables(samples):
    # Check the decontamination tables against the samples of this run, with the
    # selector's rules (see the CAUTION above).
    #   Takes in : DECONTAMINATION (the checked policy, above) and samples (SAMPLES,
    #              from section 12).
    #   Stops    : when a file that is set does not exist or is not UTF-8 text;
    #              when a table lacks its header columns; when a row of
    #              sample_overrides has a mode other than off, include, exclude or
    #              blank (auto gets the message that explains its removal); when the
    #              discard_no_hit cell of a sample of this run is not a yes/no word;
    #              and when a sample of this run would be filtered in include mode
    #              without a single genus, which the selector refuses ("Mode
    #              'include' requires at least one genus").
    #   Produces : nothing, and prints nothing when all is well.
    # Called once, at the end of section 12, and only when SCREEN is True.
    include_by_sample_path = DECONTAMINATION["include_genera_by_sample"]
    exclude_file_path = DECONTAMINATION["exclude_genera_file"]
    overrides_path = DECONTAMINATION["sample_overrides"]

    # exclude_genera_file: genera to remove, one per line. It only has to exist and be
    # UTF-8 text: the selector skips blank lines and lines starting with #, and reads
    # any other line as a list of genera.
    if exclude_file_path:
        if not os.path.isfile(exclude_file_path):
            sys.exit(f"[FunFlux] parameters.decontamination.exclude_genera_file: '{exclude_file_path}' is not an existing file.")
        _stop_if_not_utf8("exclude_genera_file", exclude_file_path)

    # include_genera_by_sample: collect the genera of each sample. A sample may have
    # several rows, and one genus cell may hold several genera.
    table_genera = {}   # sample name -> the genera its rows name, in file order
    if include_by_sample_path:
        rows = _read_decontam_table("include_genera_by_sample", include_by_sample_path, ["sample", "genus"])
        for row in rows:
            sample = str(row.get("sample") or "").strip()
            if sample not in table_genera:
                table_genera[sample] = []
            table_genera[sample].extend(_split_genus_cell(row.get("genus")))

    # sample_overrides: rows whose sample cell is empty or starts with # are skipped;
    # if a sample has two rows, the last one counts. The mode of every other row is
    # checked, whether or not its sample is in this run, because the selector reads
    # the whole table for every sample and stops on any bad mode. The discard_no_hit
    # cell is only kept here; it is checked below, for the samples of this run.
    override_mode = {}      # sample name -> its mode cell, lower case ("" keeps the run-wide mode)
    override_genera = {}    # sample name -> the genera of its include_genera cell
    override_discard = {}   # sample name -> its discard_no_hit cell, lower case ("" keeps the run-wide setting)
    if overrides_path:
        rows = _read_decontam_table("sample_overrides", overrides_path, ["sample", "mode"])
        for row in rows:
            sample = str(row.get("sample") or "").strip()
            if not sample or sample.startswith("#"):
                continue
            mode = str(row.get("mode") or "").strip().lower()
            if mode == "auto":
                sys.exit(
                    f"[FunFlux] parameters.decontamination.sample_overrides: sample '{sample}' has mode 'auto' "
                    f"in '{overrides_path}'. Mode auto is not available in FunFlux. {_AUTO_MODE_REMOVED_TEXT}"
                )
            if mode not in {"", "off", "include", "exclude"}:
                sys.exit(
                    f"[FunFlux] parameters.decontamination.sample_overrides: sample '{sample}' has mode '{mode}' "
                    f"in '{overrides_path}'; the mode column must be off, include, exclude or blank "
                    "(blank keeps the run-wide mode)."
                )
            override_mode[sample] = mode
            override_genera[sample] = _split_genus_cell(row.get("include_genera"))
            override_discard[sample] = str(row.get("discard_no_hit") or "").strip().lower()

    # For each sample of this run, the mode and the include list the selector will
    # use, and a note of the samples that would be filtered in include mode with
    # nothing to include:
    #   mode         : the run-wide mode, unless the sample's override row has a mode;
    #   include list : include_genera plus the sample's rows of include_genera_by_sample,
    #                  unless the sample's override row names genera, which then
    #                  replace that list (a blank cell keeps it).
    # The discard_no_hit cell of the sample's override row is checked here too: the
    # selector reads that cell for its own sample only, and stops on a word that is
    # not one of yes_no_words (its parse_bool; the same words as any yes/no setting of
    # the configuration, _config_bool in section 1b).
    yes_no_words = {"true", "1", "yes", "on", "false", "0", "no", "off"}
    run_wide_genera = _split_genus_cell(DECONTAMINATION["include_genera"])
    samples_without_genus = []
    for sample in samples:
        mode = DECONTAMINATION["mode"]
        include_list = run_wide_genera + table_genera.get(sample, [])
        if sample in override_mode:
            if override_mode[sample]:
                mode = override_mode[sample]
            if override_genera[sample]:
                include_list = override_genera[sample]
            discard_cell = override_discard[sample]
            if discard_cell and discard_cell not in yes_no_words:
                sys.exit(
                    f"[FunFlux] parameters.decontamination.sample_overrides: sample '{sample}' has discard_no_hit "
                    f"'{discard_cell}' in '{overrides_path}'; the discard_no_hit column must be true or false "
                    "(also yes/no, on/off, 1/0) or blank (blank keeps the run-wide setting)."
                )
        if mode == "include" and not include_list:
            samples_without_genus.append(sample)
    if samples_without_genus:
        sys.exit(
            "[FunFlux] parameters.decontamination: mode include keeps only the genera you name, but no genus "
            f"is given for sample(s) {', '.join(samples_without_genus)}. Name the genus in include_genera, "
            "in include_genera_by_sample (one row per sample) or in the include_genera column of sample_overrides."
        )


# ── 8b. Mitochondrion screen ──
# A fungal mitochondrial genome is small (usually some tens of kb) and arrives from
# every assembler as an ordinary contig. A
# contig is called mitochondrial when a BLAST hit to an organelle genome in NCBI
# core_nt covers enough of it (min_direct_blast_qcov percent). In the long-read
# modes, where the assembler also reports coverage and circularity, such a hit is
# confirmed only when that evidence agrees (coverage well above the
# nuclear contigs, or a circular contig), unless the hit is near-complete
# (very_strong_qcov). In contigs and illumina mode there is no such evidence, and
# the BLAST hit alone decides (15_mitochondria.smk has the details). Confirmed
# contigs are protected from the taxonomy selector and are written to their own
# file, so the nuclear gene predictor never sees them.
#   off                : no contig is examined and nothing is split off;
#                        contigs_final.fasta holds everything.
#   report             : write the screen table, split nothing off.
#   separate_confirmed : write the table and split confirmed contigs into their own
#                        file (the default).
# The four values go to rule screen_mitochondria through its params.
_mito = _parameter_block(
    "mitochondria",
    {"mode", "min_direct_blast_qcov", "min_assembler_coverage_ratio", "min_flye_coverage_ratio", "very_strong_qcov"},
)
MITO_MODE = _mode_word(_mito.get("mode"), "separate_confirmed")
if MITO_MODE not in {"off", "report", "separate_confirmed"}:
    sys.exit(
        "[FunFlux] parameters.mitochondria.mode: must be one of off | report | separate_confirmed "
        f"(got: {_mito.get('mode')!r})"
    )
MITO_MIN_QCOV = _to_float(
    "parameters.mitochondria.min_direct_blast_qcov", _mito.get("min_direct_blast_qcov", 35),
    minimum=0, minimum_exclusive=True, maximum=100,
)
# min_flye_coverage_ratio is accepted as another name of this key, with a NOTE; when
# both are given, min_assembler_coverage_ratio wins. The key does not name Flye because
# both long-read assemblers supply the number. Default 3.0 (contig depth at least 3
# times the median depth).
if "min_assembler_coverage_ratio" in _mito:
    _ratio_value = _mito["min_assembler_coverage_ratio"]
elif "min_flye_coverage_ratio" in _mito:
    _ratio_value = _mito["min_flye_coverage_ratio"]
    _note("parameters.mitochondria.min_flye_coverage_ratio is read as "
          "min_assembler_coverage_ratio; the value is used as it is. Write "
          "min_assembler_coverage_ratio in new configurations.")
else:
    _ratio_value = 3.0
MITO_MIN_RATIO = _to_float(
    "parameters.mitochondria.min_assembler_coverage_ratio", _ratio_value, minimum=0, minimum_exclusive=True
)
MITO_VERY_STRONG_QCOV = _to_float(
    "parameters.mitochondria.very_strong_qcov", _mito.get("very_strong_qcov", 70),
    minimum=0, minimum_exclusive=True, maximum=100,
)
# Say what the screen will do, but only in a run that has the screen. (:g prints
# 35.0 as 35.)
if SCREEN:
    print(
        f"Mitochondrion screen: {MITO_MODE} (direct BLAST coverage >= {MITO_MIN_QCOV:g}%, "
        f"assembler coverage ratio >= {MITO_MIN_RATIO:g}, near-complete hit >= {MITO_VERY_STRONG_QCOV:g}%)."
    )

# ── 8c. BLAST hit depth (docs/methods/blast-hit-depth.md) ──
# One BLAST run per contig feeds BOTH BlobTools (which sums hit scores per genus) and
# the mitochondrion screen (which reads how much of the contig an organelle genome
# covers, the BLAST column "qcovs"). qcovs is worked out from the alignments BLAST
# SAVED for each subject. With -max_hsps 1 only the single best local alignment
# survives, so a complete mitochondrial genome that aligns to its closest reference
# in several blocks is reported as only partly covered. That is why max_hsps
# defaults to 5. max_target_seqs is applied DURING the
# search, not as a final "best N" filter, so a small value can drop better subjects;
# BLAST itself recommends at least 5.
# All three values go to rule blast_contigs (10_decontam.smk). NT_VERSION is the
# name of the BLAST database inside directories.blast_db (core_nt unless changed).
# The two NOTEs below only appear in a run that has the screen.
NT_VERSION = str(PARAMETERS.get("nt_version") or "core_nt").strip()
BLAST_MAX_TARGET_SEQS = _to_int(
    "parameters.blast_max_target_seqs", PARAMETERS.get("blast_max_target_seqs", 50), minimum=1
)
BLAST_MAX_HSPS = _to_int(
    "parameters.blast_max_hsps", PARAMETERS.get("blast_max_hsps", 5), minimum=1
)
if SCREEN and BLAST_MAX_HSPS < 5 and MITO_MODE != "off":
    _note(f"parameters.blast_max_hsps is {BLAST_MAX_HSPS}. Below 5, BLAST reports the coverage of a "
          "multi-block mitochondrial match as that of its best block only, which weakens the "
          "mitochondrion screen (a complete mitogenome can fall under min_direct_blast_qcov).")
if SCREEN and BLAST_MAX_TARGET_SEQS < 5:
    _note(f"parameters.blast_max_target_seqs is {BLAST_MAX_TARGET_SEQS}; BLAST recommends 5 or more.")


# ═════════════ 9. Annotation and QC settings ═════════════
# The settings of the stages every mode shares after the assembly: repeat masking
# and gene prediction (40_annotation.smk, 41_genemark.smk), the SignalP 6 search for
# secreted proteins (45_signalp6.smk), ITS taxonomy (30_taxonomy.smk) and BUSCO
# (20_qc.smk), the contig summary (25_contig_summary.smk; section 9h, which also
# works out the inputs that rule has in this mode) and the screen for antifungal
# resistance mutations (50_amr.smk; section 9i, which also works out the inputs that
# rule has in this mode); plus the SPAdes k-mer list of the Illumina modes (illumina and
# hybrid; section 9g).

# ── 9a. Repeat masking ──
# Genes are predicted on a genome whose repeats are soft-masked (lowercase), so that
# repeat sequence does not produce false genes. Rule repeat_masking reads all
# four settings below.
#   tantan                     : fast; masks low-complexity and tandem repeats.
#   repeatmodeler_repeatmasker : builds a repeat library from THIS genome with
#                                RepeatModeler, then masks with RepeatMasker. Slow
#                                (hours) but finds interspersed repeats too.
# The two RepeatModeler switches add -quick (default on) and -LTRStruct (default
# off) to the RepeatModeler command; tantan ignores them.
MASKING_METHOD = str(PARAMETERS.get("masking_method") or "tantan").strip().lower()
if MASKING_METHOD not in {"tantan", "repeatmodeler_repeatmasker"}:
    sys.exit(
        "[FunFlux] parameters.masking_method: must be tantan | repeatmodeler_repeatmasker "
        f"(got: {PARAMETERS.get('masking_method')!r})"
    )
REPEATMODELER_QUICK    = _config_bool("parameters.repeatmodeler_quick", PARAMETERS.get("repeatmodeler_quick", True), True)
REPEATMODELER_LTRSTRUCT = _config_bool("parameters.repeatmodeler_ltrstruct", PARAMETERS.get("repeatmodeler_ltrstruct", False), False)
# RepeatModeler starts from a random seed taken from the clock, so two runs on the
# same genome build slightly different libraries. An integer here fixes the seed
# (RepeatModeler -srand) of its random choice of genome samples; whether the masking then
# repeats exactly has not been tested. Empty keeps the v1 behaviour.
# REPEATMODELER_SEED is text: "" (no seed) or the number written as text.
_seed_value = PARAMETERS.get("repeatmodeler_seed")
if _is_empty_value(_seed_value):
    REPEATMODELER_SEED = ""
else:
    REPEATMODELER_SEED = str(_to_int("parameters.repeatmodeler_seed", _seed_value, minimum=0))

# ── 9b. ITS extraction and classification ──
# The fungal ITS barcode is cut out of the assembled genome and classified against
# the UNITE reference with SINTAX. In every mode a window is first cut around each
# rDNA copy (rule its_windows), because ITSx fails on a multi-megabase contig when
# the array sits on the minus strand, as it can in a long-read assembly
# (docs/methods/its-long-read-contigs.md; the windows are made by
# workflow/scripts/30_taxonomy/its_windows.py).
#   ITS_TAXONOMY_CUTOFF : the SINTAX confidence a rank needs to be reported (above 0,
#                         at most 1, default 0.8; rule its_classification).
#   ITS_WINDOW_FLANK    : bases kept on each side of every 5.8S gene (default 2000,
#                         at least 200; rule its_windows).
#   ITS_MAX_FULL_LENGTH : the longest whole ITS1-5.8S-ITS2 record that is still
#                         accepted (default 1500, at least 400; rule its_classification).
#   UNITE_IDENTIFIED_ONLY : yes/no, default yes: classify against the UNITE records
#                         whose taxonomy names a genus, not against the whole release
#                         (rule unite_identified makes that file). About a third of
#                         the records of the default release stop above the genus
#                         ("tax=d:Fungi;" alone, for example); when such a record is
#                         identical to the query, SINTAX follows it and the result
#                         stops at the kingdom even for a well-identified species.
#   UNITE_REFERENCE     : the FASTA that rule its_classification hands to vsearch
#                         --db: UNITE_FASTA_IDENTIFIED when the filter is on,
#                         UNITE_FASTA (the whole release) when it is off. Both paths
#                         are set in section 5.
ITS_TAXONOMY_CUTOFF = _to_float(
    "parameters.its_taxonomy_cutoff", PARAMETERS.get("its_taxonomy_cutoff", 0.8),
    minimum=0, minimum_exclusive=True, maximum=1,
)
_its = _parameter_block("its", {"window_flank", "max_full_length", "unite_identified_only"})
ITS_WINDOW_FLANK    = _to_int("parameters.its.window_flank", _its.get("window_flank", 2000), minimum=200)
ITS_MAX_FULL_LENGTH = _to_int("parameters.its.max_full_length", _its.get("max_full_length", 1500), minimum=400)
UNITE_IDENTIFIED_ONLY = _config_bool(
    "parameters.its.unite_identified_only", _its.get("unite_identified_only"), True
)
if UNITE_IDENTIFIED_ONLY:
    UNITE_REFERENCE = UNITE_FASTA_IDENTIFIED
else:
    UNITE_REFERENCE = UNITE_FASTA
# Length guard applied before SINTAX (rule its_classification): ITS1 and ITS2 are
# 50-600 bp in fungi; a record outside that range is a mis-extraction (v1 could
# classify a huge "ITS2" as a non-fungus) and is listed in the ITS copies table
# with the reason instead. Fixed numbers, not settings.
ITS_MIN_PART_LENGTH = 50
ITS_MAX_PART_LENGTH = 600

# ── 9c. Genome completeness (BUSCO) ──
# auto lets BUSCO pick the lineage itself (--auto-lineage-euk), which is what the v1
# workflows did; a lineage name such as hypocreales_odb12 pins it. The name is
# lower-cased and must look like <name>_odb<number>. BUSCO_LINEAGE_ARGS is the
# matching option text, which rule completeness (20_qc.smk) puts on the BUSCO command
# line.
#
# Where the lineage data (BUSCO's gene sets) come from depends on that choice:
#   auto        : downloaded inside each BUSCO job, into that sample's own folder. The
#                 specific lineage is known only once BUSCO has placed the genome, and
#                 BUSCO has no download lock, so two samples fetching into one shared
#                 folder at the same time could corrupt it. The computer that runs the
#                 job needs internet.
#   a lineage   : downloaded ONCE per output folder, by rule busco_data (20_qc.smk),
#                 into BUSCO_DATA_DIR (00.resources/busco/), on the computer that runs
#                 Snakemake itself (a local rule). Every BUSCO job then reads it with
#                 --offline: it never downloads and never writes there, so samples can
#                 share it. A cluster whose compute nodes have no internet can run
#                 BUSCO this way.
# Three names carry that difference to the rules in 20_qc.smk:
#   BUSCO_DATA_READY    the file rule busco_data makes: dataset.cfg, the description file
#                       every BUSCO lineage folder holds.
#   BUSCO_DATA_INPUT    what rule completeness waits for: nothing in auto mode, that file
#                       when a lineage is named.
#   BUSCO_DOWNLOAD_ARGS the BUSCO options that say where the data are: in auto mode
#                       --download_path <the sample's busco_tmp folder> ({sample} in the
#                       path is filled in by Snakemake for each job); with a named lineage
#                       --offline --download_path <BUSCO_DATA_DIR>.
_busco = _parameter_block("busco", {"lineage"})
BUSCO_LINEAGE = str(_busco.get("lineage") or "auto").strip().lower()
# The message also explains the most likely mistake: `busco --list-datasets` prints each
# name with its dataset version after a dot (hypocreales_odb12.2), and a name copied from
# that list keeps the dot. BUSCO's -l option takes the name without it and fetches the
# newest version itself.
if BUSCO_LINEAGE != "auto" and not re.fullmatch(r"[a-z0-9_]+_odb[0-9]+", BUSCO_LINEAGE):
    sys.exit(
        "[FunFlux] parameters.busco.lineage: must be 'auto' or a BUSCO lineage name such as "
        f"hypocreales_odb12 (got: {_busco.get('lineage')!r}). busco --list-datasets prints "
        "each name with its version after a dot, as in hypocreales_odb12.2: leave out the "
        "dot and the number after it."
    )
BUSCO_LINEAGE_ARGS = "--auto-lineage-euk" if BUSCO_LINEAGE == "auto" else f"-l {BUSCO_LINEAGE}"
BUSCO_DATA_DIR = DIR_RESOURCES + "/busco"
BUSCO_DATA_READY = BUSCO_DATA_DIR + "/lineages/" + BUSCO_LINEAGE + "/dataset.cfg"
if BUSCO_LINEAGE == "auto":
    BUSCO_DATA_INPUT = []
    BUSCO_DOWNLOAD_ARGS = "--download_path " + BUSCO_DIR + "/busco_tmp"
else:
    BUSCO_DATA_INPUT = [BUSCO_DATA_READY]
    BUSCO_DOWNLOAD_ARGS = "--offline --download_path " + BUSCO_DATA_DIR

# ── 9d. funannotate gene prediction ──
# funannotate predict trains the ab-initio gene finders (programs that recognise genes
# from the sequence alone) SNAP and GlimmerHMM, and AUGUSTUS for the samples without a
# model in column 4 of the annotation table, on a set of conserved genes taken from a
# BUSCO dataset bundled in the funannotate database.
# Only 'dikarya' ships in the database snapshot; a fungus outside the Dikarya needs
# another set (fungi, for example) added to the database first, with the command
# the stop message below gives. PREDICTION_BUSCO_DB goes to rule
# funannotate_prediction (--busco_db). The name may only hold letters, digits and
# underscores.
_prediction = _parameter_block("prediction", {"busco_db", "augustus_optimize"})
PREDICTION_BUSCO_DB = str(_prediction.get("busco_db") or "dikarya").strip()
if not re.fullmatch(r"[A-Za-z0-9_]+", PREDICTION_BUSCO_DB):
    sys.exit(f"[FunFlux] parameters.prediction.busco_db: expected a plain dataset name such as dikarya (got: {PREDICTION_BUSCO_DB!r})")
# Check the dataset only when the database is already on disk; otherwise the
# database is about to be downloaded and there is nothing to look at yet.
# The stop message gives the command that adds a set:
#   funannotate setup -i busco -b <set> -d <database folder> -w
#   -i busco : download BUSCO sets only. Setup also copies into trained_species/
#              the AUGUSTUS species files the folder lacks, and rewrites
#              funannotate-db-info.txt from the lines it can read.
#   -w       : download with wget; funannotate's own downloader has failed on some
#              of its links (an HTTP 308 redirect, HTTP 403 from geneontology.org).
# Run it with FunFlux's funannotate environment activated, not only put on the PATH:
# setup needs AUGUSTUS_CONFIG_PATH, which activation sets. The message points to the
# documentation for finding that environment, because snakemake --list-conda-envs
# reads this configuration too and stops here until busco_db names a set the
# database holds. Tested on 2026-10-08 with the small microsporidia set, not with
# fungi (docs/getting-started/databases.md).
if os.path.isfile(os.path.join(FUNANNOTATEDB, "funannotate-db-info.txt")):
    if not os.path.isdir(os.path.join(FUNANNOTATEDB, PREDICTION_BUSCO_DB)):
        sys.exit(
            f"[FunFlux] parameters.prediction.busco_db: '{PREDICTION_BUSCO_DB}' is not in the funannotate "
            f"database at '{FUNANNOTATEDB}'. Add it with FunFlux's funannotate conda environment activated "
            "('Other BUSCO sets for gene prediction' in the documentation shows how to find it): "
            f"funannotate setup -i busco -b {PREDICTION_BUSCO_DB} -d {FUNANNOTATEDB} -w"
        )

# ── 9d (continued). The AUGUSTUS model: named, or trained on this genome ──
# AUGUSTUS weighs most in funannotate's final gene set: its models with strong evidence
# count double in EVidenceModeler. Column 4 of the annotation table (section 12) decides,
# per sample, which parameters it runs with:
#   * a name, such as fusarium_graminearum: funannotate copies that species' stored
#     AUGUSTUS parameters from the funannotate database (trained_species/<name>/) and uses
#     them AS THEY ARE. AUGUSTUS is not trained on this genome.
#   * empty or "-": funannotate trains AUGUSTUS on this genome's own conserved genes, the
#     BUSCO genes that also train SNAP and GlimmerHMM (parameters.prediction.busco_db).
#     It needs at least 200 of them (funannotate's --min_training_models) and stops
#     otherwise, late in the prediction.
# funannotate's log says which of the two it did ("augustus pretrained" or "busco"), but it
# gives no warning when a name is not in the database (a typo, or the wrong capitals): it
# then quietly trains AUGUSTUS under that name. So FunFlux checks every name against the
# database and stops on one that is not there, with the closest names it does hold. The
# check is made by workflow/scripts/40_annotation/check_augustus_model.py: in section 12,
# while Snakemake reads the workflow, when the database is already on disk, and otherwise
# by rule augustus_model_check, right after the database is unpacked (40_annotation.smk).
# The two choices can change a large share of the final gene models, so the choice is
# recorded:
# augustus_model_check writes it into AUGUSTUS_MODEL_RECORD, and the gene prediction
# copies that record into its log.
#
# parameters.prediction.augustus_optimize (default false) adds funannotate's
# --optimize_augustus, its "long" training, for the samples whose column 4 is empty. After
# the quick training (etraining, AUGUSTUS's training program, a few seconds),
# optimize_augustus.pl tunes the model's meta-parameters, the settings the quick training
# keeps fixed (such as the width of the windows around splice sites in which AUGUSTUS
# learns sequence patterns), by 8-fold cross-validation: the test genes (the tenth of the
# BUSCO models funannotate sets aside for testing) are cut into 8 parts, and
# each part is predicted by a model trained without it. Then AUGUSTUS is trained
# again. That adds about a quarter of an hour for a 40 Mb genome (funannotate lets it use at
# most 8 CPUs).
PREDICTION_AUGUSTUS_OPTIMIZE = _config_bool(
    "parameters.prediction.augustus_optimize", _prediction.get("augustus_optimize"), False
)

# check_augustus_model.py is the one place that knows how funannotate decides: section 12
# below and rule augustus_model_check both use it. importlib loads a Python file by its
# path as a module, so its functions can be called here (it is not an installed package).
# The three steps: find the file (spec), make an empty module from it, run the file's code
# into that module. AUGUSTUS_CHECK.model_problem and AUGUSTUS_CHECK.decide are then ready.
_augustus_check_spec = importlib.util.spec_from_file_location("check_augustus_model", AUGUSTUS_CHECK_SCRIPT)
AUGUSTUS_CHECK = importlib.util.module_from_spec(_augustus_check_spec)
_augustus_check_spec.loader.exec_module(AUGUSTUS_CHECK)


# ── 9e. GeneMark-ES (optional, separately licensed) ──
# GeneMark-ES is a self-training gene finder funannotate can use. It is licensed per
# user, so it is used only when the user names an installation (directories.genemark_dir);
# otherwise funannotate is shielded from any GeneMark that happens to be on the PATH and
# follows its supported no-GeneMark branch.
#   auto : use GeneMark when directories.genemark_dir is set. A mode key left
#          empty (`mode:`, `mode: ""` or an unquoted null) also means auto, as
#          for the other mode keys.
#   off  : never use it. false, no and 0 mean off too, and so do the words none
#          and 'null' (in quotes).
# The licence is a key file, ~/.gm_key, that EXPIRES. An expired key does not show
# up at start-up: gmes_petap.pl only checks that the file exists, and the failure
# arrives well into funannotate predict as "Your license period has ended".
# Rule genemark_check (41_genemark.smk) therefore test-runs GeneMark on a small
# slice of each genome first. When GeneMark is on, this section checks the folder,
# the gmes_petap.pl program and the key file now, before any job starts.
# GENEMARK_ENABLED decides whether 41_genemark.smk defines its rule and whether
# funannotate_prediction hands GeneMark to funannotate; GENEMARK_KEY is the key path
# named in its messages.
_genemark = _parameter_block("genemark", {"mode"})
# _mode_word (section 1b) turns a missing, null or blank value into "auto" and an
# unquoted off/no/false (a YAML False) into "off"; any other value arrives as its
# lower-case text, so the off words typed as text are listed here.
_genemark_mode = _mode_word(_genemark.get("mode"), "auto")
if _genemark_mode in {"off", "false", "0", "no", "none", "null"}:
    GENEMARK_ENABLED = False
elif _genemark_mode == "auto":
    GENEMARK_ENABLED = bool(GENEMARKDIR)
else:
    sys.exit(f"[FunFlux] parameters.genemark.mode: must be auto | off (got: {_genemark.get('mode')!r})")
GENEMARK_KEY = os.path.expanduser("~/.gm_key")
if GENEMARK_ENABLED:
    if not os.path.isabs(GENEMARKDIR):
        sys.exit("[FunFlux] directories.genemark_dir: must be an absolute path when GeneMark is enabled.")
    _gmes = os.path.join(GENEMARKDIR, "gmes_petap.pl")
    if not (os.path.isfile(_gmes) and os.access(_gmes, os.X_OK)):
        sys.exit(f"[FunFlux] directories.genemark_dir: '{GENEMARKDIR}' has no executable gmes_petap.pl.")
    if not os.path.isfile(GENEMARK_KEY):
        sys.exit(
            f"[FunFlux] GeneMark is enabled but there is no licence key at {GENEMARK_KEY}. Download your key "
            "from the GeneMark site, or set parameters.genemark.mode: off to run without GeneMark."
        )
    _note(f"GeneMark-ES is enabled and needs a valid, unexpired licence key at {GENEMARK_KEY}. An expired key "
          "stops funannotate predict late in the run; renew it before a long run.")
else:
    print("GeneMark-ES: off (funannotate runs its no-GeneMark branch).")

# ── 9e (continued). Tiberius, a second gene finder beside GeneMark-ES (on by default) ──
# Tiberius finds genes with a deep neural network joined to a hidden Markov model (the
# kind of model AUGUSTUS and GeneMark use for the order of exons, introns and the DNA
# between genes). It is not trained on the genome being annotated: its authors trained
# one "fungi" model on many fungal genomes, and FunFlux uses that model as it is. The
# three rules of 42_tiberius.smk download the model, run Tiberius on each sample's
# masked genome and convert its output; funannotate predict then hands Tiberius's gene
# models to EVidenceModeler, the step that combines all gene finders and protein
# alignments into one gene set, as one more source of evidence. GeneMark-ES
# is independent of this: with both on, both go into that combination. The section
# "Tiberius" of docs/analysis/gene-prediction.md describes it for users.
#   mode         : on (the default; an empty value, true, yes and 1 mean on too) | off
#                  (false, no, 0, none and 'null' mean off too). YAML turns an unquoted
#                  on into the boolean True, which _mode_word hands over as "true".
#   model_sha256 : the SHA-256 checksum (the fingerprint described at rule
#                  funannotate_db) that the downloaded model archive must have: 64
#                  characters 0-9 and a-f (capitals are accepted and lowered). The
#                  default is the fungi model of 8 February 2026, the one FunFlux v2 was
#                  tested with. The download address (links.tiberius_model_link, read in
#                  section 10) has no version number, so the authors can replace the file
#                  at the same address; the checksum turns that into a stop with a
#                  message instead of a silent change of the gene models.
# Two values are fixed, not settings:
#   TIBERIUS_SEQ_LEN    : Tiberius cuts the genome into pieces of this many bases for the
#                         network. 100080 is the fungi model's own value (it must be a
#                         multiple of 18). Tiberius 1.1.8 does not read it from the model
#                         and would use 500004, which stops it with the TensorFlow 2.17 of
#                         its environment (anything above 259992 does).
#   TIBERIUS_EVM_WEIGHT : the weight EVidenceModeler gives Tiberius's gene models, 1, the
#                         same as GeneMark-ES, SNAP, GlimmerHMM and the protein alignments
#                         (AUGUSTUS gets 1, and 2 for its models that the evidence
#                         supports well). Fixed until a public reference genome shows that
#                         another value gives better genes.
# TIBERIUS_ENABLED decides whether 42_tiberius.smk defines its rules and whether
# funannotate_prediction waits for and passes on their result (tiberius_gate and
# tiberius_arguments, section 12).
TIBERIUS_DEFAULT_MODEL_SHA256 = "4a309a630574f395d82d747dcf64716704dab056db534e427ae0cec1b7410e20"
TIBERIUS_SEQ_LEN = 100080
TIBERIUS_EVM_WEIGHT = 1
_tiberius = _parameter_block("tiberius", {"mode", "model_sha256"})
_tiberius_mode = _mode_word(_tiberius.get("mode"), "on")
if _tiberius_mode in {"on", "true", "yes", "1"}:
    TIBERIUS_ENABLED = True
elif _tiberius_mode in {"off", "false", "no", "0", "none", "null"}:
    TIBERIUS_ENABLED = False
else:
    sys.exit(f"[FunFlux] parameters.tiberius.mode: must be on | off (got: {_tiberius.get('mode')!r})")
TIBERIUS_MODEL_SHA256 = str(_tiberius.get("model_sha256") or TIBERIUS_DEFAULT_MODEL_SHA256).strip().lower()
if not re.fullmatch(r"[0-9a-f]{64}", TIBERIUS_MODEL_SHA256):
    sys.exit(
        "[FunFlux] parameters.tiberius.model_sha256: expected the SHA-256 checksum of the model archive, "
        f"64 characters 0-9 and a-f (got: {_tiberius.get('model_sha256')!r})"
    )
if TIBERIUS_ENABLED:
    print("Tiberius: on (its gene models join the others in EVidenceModeler, weight 1).")
else:
    print("Tiberius: off.")

# ── 9f. SignalP 6 (optional, licensed) ──
# SignalP 6 predicts which proteins are secreted. It is proprietary DTU software: a
# licensed user downloads the official archive by hand, and the workflow installs that
# exact local archive into its own conda environment. It never downloads SignalP and
# never contains a DTU URL.
#   mode              : off (the default; an empty value, false, no, none or null
#                       also mean off) | fast | slow | slow-sequential.
#   archive_sha256    : optional checksum of the archive, checked by the installer.
#   torch_num_threads : threads for SignalP's neural network, and the thread count of
#                       rule signalp6_annotation (default: CPUS or 8, whichever is
#                       smaller). A larger number may be set, up to CPUS; a number
#                       above CPUS is lowered to CPUS, with a NOTE.
#   batch_size        : proteins per batch (default 10).
# When SignalP 6 is on, the archive (files.signalp6_archive) must be an absolute path
# to an existing file; that is checked now. SIGNALP6_ENABLED decides whether
# 45_signalp6.smk defines its two rules (install_signalp6, signalp6_annotation),
# whether funannotate_annotation waits for their result, and whether rule all asks
# for it; the other values go to those two rules.
SIGNALP6_ARCHIVE = _optional_text("files", "signalp6_archive")
_signalp6 = _parameter_block("signalp6", {"mode", "archive_sha256", "torch_num_threads", "batch_size"})
_signalp6_mode = str(_signalp6.get("mode", "off")).strip().lower()
SIGNALP6_ARCHIVE_SHA256 = str(_signalp6.get("archive_sha256") or "").strip().lower()
if _signalp6_mode in {"", "off", "false", "0", "no", "none", "null"}:
    SIGNALP6_ENABLED = False
    SIGNALP6_MODE = ""
elif _signalp6_mode in {"fast", "slow", "slow-sequential"}:
    SIGNALP6_ENABLED = True
    SIGNALP6_MODE = _signalp6_mode
else:
    sys.exit(
        "[FunFlux] parameters.signalp6.mode: must be off | fast | slow | slow-sequential "
        f"(got: {_signalp6.get('mode')!r})"
    )
if SIGNALP6_ENABLED:
    if not SIGNALP6_ARCHIVE:
        sys.exit("[FunFlux] files.signalp6_archive: required when parameters.signalp6.mode is not off.")
    if not os.path.isabs(SIGNALP6_ARCHIVE):
        sys.exit("[FunFlux] files.signalp6_archive: must be an absolute path to the official DTU archive.")
    if not os.path.isfile(SIGNALP6_ARCHIVE):
        sys.exit(f"[FunFlux] files.signalp6_archive: '{SIGNALP6_ARCHIVE}' is not a readable file.")
# The thread count is capped at CPUS (resources.threads, the budget of the whole
# run) with capped_cpus (section 4), so that SignalP stays within that budget like
# every other rule. Snakemake itself only lowers a rule's threads to --cores, and
# --cores may be set higher than resources.threads, so without the cap SignalP could
# take more threads than the run was given. The NOTE appears only when SignalP 6 is
# on, since otherwise the setting is not used.
_torch_threads_requested = _to_int(
    "parameters.signalp6.torch_num_threads", _signalp6.get("torch_num_threads", min(CPUS, 8)), minimum=1
)
SIGNALP6_TORCH_NUM_THREADS = capped_cpus(_torch_threads_requested)
if SIGNALP6_ENABLED and _torch_threads_requested > CPUS:
    _note(f"parameters.signalp6.torch_num_threads is {_torch_threads_requested}, more than resources.threads "
          f"({CPUS}); SignalP 6 will use {CPUS} threads.")
SIGNALP6_BATCH_SIZE = _to_int("parameters.signalp6.batch_size", _signalp6.get("batch_size", 10), minimum=1)
print(f"SignalP 6: {SIGNALP6_MODE if SIGNALP6_ENABLED else 'off'}.")

# ── 9g. SPAdes k-mer sizes (Illumina assembly) ──
# The k-mer list SPAdes builds its graphs with. The default is the list the v1
# workflows hard-coded, 21 to 127: on a 40 Mb fungus sequenced at about 80x with
# 150 bp reads it gave fewer and longer contigs than SPAdes' own choice (21 to 77),
# at the same completeness. Each k must be odd and below 128, and the list ascending.
# `auto` (any case) hands the choice to SPAdes, which picks the list from the longest
# read (21,33,55 below 150 bp, 21,33,55,77 from 150 bp, 21 to 127 from 250 bp); with
# thin coverage (some 20-30x) the large k of the fixed list may hold too few reads.
# Read only in illumina and hybrid mode; rule illumina_assembly passes SPADES_KMERS
# to SPAdes (-k; SPAdes itself accepts `-k auto`). The setting may be text
# ("21,33,55") or a YAML list ([21, 33, 55]); either way it ends up as text with
# commas and no spaces.
SPADES_KMERS = "21,33,55,77,99,127"
if HAS_SHORT_READS:
    _kmers_value = PARAMETERS.get("spades_kmers", SPADES_KMERS)
    if isinstance(_kmers_value, (list, tuple)):
        _kmers_value = ",".join(str(k) for k in _kmers_value)
    _kmers_text = str(_kmers_value).replace(" ", "")
    _kmers_is_auto = _kmers_text.lower() == "auto"
if HAS_SHORT_READS and _kmers_is_auto:
    SPADES_KMERS = "auto"
elif HAS_SHORT_READS:
    # Every item must be a whole number ("21,,33" or "21;33" is refused here).
    try:
        _kmers = [int(k) for k in _kmers_text.split(",")]
    except ValueError:
        sys.exit(
            "[FunFlux] parameters.spades_kmers: expected auto, or whole numbers separated by commas "
            f"(got: {_kmers_value!r})."
        )
    # Every k must be odd, at least 1 and below 128 ...
    _kmer_not_allowed = any(k % 2 == 0 or k >= 128 or k < 1 for k in _kmers)
    # ... and the list strictly ascending: sorted, with no k given twice.
    _kmers_ascending = _kmers == sorted(set(_kmers))
    if _kmer_not_allowed or not _kmers_ascending:
        sys.exit(f"[FunFlux] parameters.spades_kmers: every k must be odd and below 128, and the list ascending (got: {_kmers_text}).")
    SPADES_KMERS = _kmers_text

# ── 9h. The contig summary (25_contig_summary.smk) ──
# After the annotation, two rules say what every contig of the draft assembly is. Rule
# telomere_units works out the telomere repeat unit of the assembly with tidk and calls
# the telomere at both ends of every contig; rule contig_summary joins that call with what
# the run learnt about every contig from files other rules wrote, in one table per sample
# (CONTIG_SUMMARY, section 5), and gives each contig a verdict: chromosome, chromosome arm,
# mitochondrion, rDNA array, low-depth copy, collapsed repeat, nuclear contig, nuclear
# fragment, contaminant or unclassified. The telomere, a tandem array of a short repeat
# unit that caps every chromosome end, is the one piece of evidence worked out from the
# sequence itself: a contig with the array at both ends is a whole chromosome, and at one
# end a chromosome arm. The unit differs between fungal groups. The settings:
#   telomere_motif      : auto (the default), or the repeat unit of the species, letters A,
#                         C, G and T only, 4 to 30 long, upper or lower case (upper-cased
#                         here, as the script does). auto: the unit is chosen from the data.
#                         Every tandem unit tidk finds at the contig ends, every unit of
#                         the tidk database and every unit the literature names is scored
#                         by how it behaves at the contig ends, and the best-supported one
#                         is the unit of the assembly (the top of 25_contig_summary.smk has
#                         the method). A unit: used as given, whatever the data say; the
#                         ranking still runs, and a NOTE in the log says when the data
#                         favour another unit or give the configured one no support.
#   telomere_db_sha256  : the SHA-256 checksum of the telomere database file at
#                         links.telomere_db_link (64 hexadecimal characters, upper or lower
#                         case; lower-cased here). Rule telomere_db checks the downloaded file
#                         against it before anything reads the file. The address is pinned
#                         to one revision of the database, so the two change together.
#   telomere_window     : the bases tidk looks at from each contig end: the window it
#                         scans, the farthest an array may lie from the end, and the length
#                         of the terminal windows in which the candidate units are
#                         collected (default 25000: the array can sit some kb inside a
#                         contig end when the assembler ran past the telomere, and the
#                         window must cover that; at least 500 and at most 100000: a wider
#                         window collects more unrelated tandem repeats as candidates, and
#                         every candidate costs one tidk run).
#   min_telomere_copies : the tandem copies of a unit that make a telomere (default 5, that
#                         is 30 bp of TTAGGG; at least 2). tidk accepts an array only from
#                         that length (the copies x the unit length) and lists a unit as a
#                         candidate only from that many tandem copies, and the summary
#                         counts a contig end as a telomere from that many copies.
#   min_rdna_units      : the rDNA units (the ITS windows of rule its_windows) from which a
#                         contig can be a piece of the rDNA array (default 3, at least 1).
#   univec_sha256       : the SHA-256 checksum of NCBI's UniVec file at links.univec_link
#                         (build 10.1 in NCBI's archive, which does not change). Rule
#                         univec_illumina checks the downloaded file against it before it
#                         takes the Illumina records out for the contig adapter screen.
# The limits are the ones the scripts enforce as well (telomere_units.py and
# contig_summary.py in workflow/scripts/25_contig_summary/), so a run of a script by hand
# is held to the same rules. The motif, the window and the copies go to rule telomere_units
# through its params, the copies and the rDNA units to rule contig_summary, so a changed
# value re-runs those steps, and MultiQC after them, and nothing else; the checksum goes to
# rule telomere_db. The line printed at the end says what the run will use.
_contig_summary = _parameter_block(
    "contig_summary",
    {"telomere_motif", "telomere_db_sha256", "telomere_window", "min_telomere_copies", "min_rdna_units",
     "univec_sha256"},
)
# The motif: the word auto in any case is stored as "auto"; anything else must be a unit.
CONTIG_SUMMARY_TELOMERE_MOTIF = str(_contig_summary.get("telomere_motif") or "auto").strip()
if CONTIG_SUMMARY_TELOMERE_MOTIF.lower() == "auto":
    CONTIG_SUMMARY_TELOMERE_MOTIF = "auto"
else:
    CONTIG_SUMMARY_TELOMERE_MOTIF = CONTIG_SUMMARY_TELOMERE_MOTIF.upper()
    if not re.fullmatch(r"[ACGT]{4,30}", CONTIG_SUMMARY_TELOMERE_MOTIF):
        sys.exit(
            "[FunFlux] parameters.contig_summary.telomere_motif: must be auto (the unit is chosen from "
            "the data) or the telomere repeat unit of the species, 4 to 30 letters A, C, G and T, such "
            f"as TTAGGG (got: {_contig_summary.get('telomere_motif')!r})"
        )
# The checksum of the database file at the address the config template pins; a copy of the
# template that leaves the key out gets it, like the other defaults of this block.
TELOMERE_DB_DEFAULT_SHA256 = "1d4f5c8c2cce6a76aaa8de01602317a8e0b8d0401711b96f3c7a826c385da699"
TELOMERE_DB_SHA256 = str(_contig_summary.get("telomere_db_sha256") or TELOMERE_DB_DEFAULT_SHA256).strip().lower()
if not re.fullmatch(r"[0-9a-f]{64}", TELOMERE_DB_SHA256):
    sys.exit(
        "[FunFlux] parameters.contig_summary.telomere_db_sha256: must be the SHA-256 checksum of the "
        f"telomere database file, 64 hexadecimal characters (got: {_contig_summary.get('telomere_db_sha256')!r})"
    )
# The checksum of UniVec build 10.1 (https://ftp.ncbi.nlm.nih.gov/pub/UniVec/archive/10.1/UniVec);
# a configuration without the key gets it, so a v2.0.0 configuration still runs.
UNIVEC_DEFAULT_SHA256 = "687d8aa872d4ee215c9e0b3b004d0ee709fb243436d36d159a3f0ac5c25b5e36"
UNIVEC_SHA256 = str(_contig_summary.get("univec_sha256") or UNIVEC_DEFAULT_SHA256).strip().lower()
if not re.fullmatch(r"[0-9a-f]{64}", UNIVEC_SHA256):
    sys.exit(
        "[FunFlux] parameters.contig_summary.univec_sha256: must be the SHA-256 checksum of the "
        f"UniVec file, 64 hexadecimal characters (got: {_contig_summary.get('univec_sha256')!r})"
    )
CONTIG_SUMMARY_TELOMERE_WINDOW = _to_int(
    "parameters.contig_summary.telomere_window", _contig_summary.get("telomere_window", 25000),
    minimum=500, maximum=100000,
)
CONTIG_SUMMARY_MIN_TELOMERE_COPIES = _to_int(
    "parameters.contig_summary.min_telomere_copies", _contig_summary.get("min_telomere_copies", 5), minimum=2
)
CONTIG_SUMMARY_MIN_RDNA_UNITS = _to_int(
    "parameters.contig_summary.min_rdna_units", _contig_summary.get("min_rdna_units", 3), minimum=1
)
if CONTIG_SUMMARY_TELOMERE_MOTIF == "auto":
    _telomere_motif_text = "auto (the repeat unit is chosen from the data)"
else:
    _telomere_motif_text = f"{CONTIG_SUMMARY_TELOMERE_MOTIF} (used as given; the data's own choice is reported next to it)"
print(
    f"Contig summary: telomere motif {_telomere_motif_text}, a telomere from "
    f"{CONTIG_SUMMARY_MIN_TELOMERE_COPIES} tandem copies within {CONTIG_SUMMARY_TELOMERE_WINDOW} bp "
    f"of a contig end, an rDNA array from {CONTIG_SUMMARY_MIN_RDNA_UNITS} rDNA unit(s)."
)

# ── 9h (continued). The inputs rule contig_summary has in this mode ──
# Three of the rule's inputs exist in every mode and are always given to the script: the ITS
# windows table, the BUSCO folder and the annotated genome with its GFF3. Five more depend
# on the mode, and they are worked out here, once, so that the rule file does not have to:
#   the three files of the screen   MITO_CONTIGS (the mitochondrial contigs the split wrote),
#                                   CONTIG_DECISIONS (every contig kept or removed, with its
#                                   genus) and MITO_SCREEN_TSV (the mitochondrion screen's
#                                   report): every run with the screen (SCREEN).
#   the BlobTools table             BLOB_TABLE, the read depth of every contig: only where
#                                   real reads were mapped to the draft (HAS_READS). In
#                                   contigs mode the screen maps the contigs to THEMSELVES
#                                   and the depth it reports reflects contig length and
#                                   nothing else (the contigs-mode rule map_contigs in
#                                   10_decontam.smk says why), so there the table is left
#                                   out even when the screen runs: given to the script, it
#                                   would make the longest contigs look shallow and call
#                                   them low-depth debris.
#   the assembler's evidence table  ASSEMBLY_EVIDENCE, the assembler's own depth and
#                                   circularity flag: the long-read modes (HAS_LONG_READS).
# Each of the five names below holds the file's path (still with the {sample} wildcard)
# where the mode makes it, and "" where it does not. The rule hands each to the script's
# matching option; an empty one reaches the script as an option without a value, which it
# reads as "not given", and the columns that need that file read NA. The list at the end
# holds the same files without the empty ones: the rule lists it as an input, so the job
# waits for them.
if SCREEN:
    CONTIG_SUMMARY_MITO_CONTIGS = MITO_CONTIGS
    CONTIG_SUMMARY_DECISIONS    = CONTIG_DECISIONS
    CONTIG_SUMMARY_MITO_SCREEN  = MITO_SCREEN_TSV
else:
    CONTIG_SUMMARY_MITO_CONTIGS = ""
    CONTIG_SUMMARY_DECISIONS    = ""
    CONTIG_SUMMARY_MITO_SCREEN  = ""
if HAS_READS:
    CONTIG_SUMMARY_BLOB_TABLE = BLOB_TABLE
else:
    CONTIG_SUMMARY_BLOB_TABLE = ""
if HAS_LONG_READS:
    CONTIG_SUMMARY_ASSEMBLY_EVIDENCE = ASSEMBLY_EVIDENCE
else:
    CONTIG_SUMMARY_ASSEMBLY_EVIDENCE = ""
CONTIG_SUMMARY_OPTIONAL_INPUTS = []
for _optional_path in (
    CONTIG_SUMMARY_MITO_CONTIGS, CONTIG_SUMMARY_DECISIONS, CONTIG_SUMMARY_MITO_SCREEN,
    CONTIG_SUMMARY_BLOB_TABLE, CONTIG_SUMMARY_ASSEMBLY_EVIDENCE,
):
    if _optional_path:
        CONTIG_SUMMARY_OPTIONAL_INPUTS.append(_optional_path)

# ── 9i. Antifungal resistance mutations with ChroQueTas (50_amr.smk; off by default) ──
# A fungus becomes resistant to an antifungal drug or a fungicide mostly through a point
# mutation in the protein the drug acts on (beta-tubulin for the benzimidazoles, cytochrome
# b for the QoI fungicides, Cyp51 for the azoles, succinate dehydrogenase for the SDHIs,
# ...), not by acquiring a resistance gene. ChroQueTas (Bédard et al. 2025) finds each
# target protein of the sample's species in the assembled genome and compares it, position
# by position, with the resistance mutations of the FungAMR database. The species is
# column 2 of the annotation table (species_name, section 12); a species the database does
# not hold is reported as "not screened", and the run goes on. Rule chroquetas_db downloads
# the program and the database once per run, and rule amr_chroquetas runs it once per
# sample; docs/analysis/amr.md describes the step for users.
#
# ChroQueTas and its database are licensed CC BY-NC-ND 4.0: non-commercial use only, and
# no modified copies may be passed on. FunFlux therefore never ships them: it downloads the
# unchanged files at run time from one pinned commit of the authors' repository and checks
# their checksums. The step is off by default, and switching it on is the user's
# acceptance of that licence, as with SignalP 6 (section 9f).
#
# The settings (parameters.amr):
#   chroquetas               : off (the default; an empty value, false, no, 0, none and
#                              'null' mean off too) | on (true, yes and 1 mean on too).
#                              YAML turns an unquoted on into the boolean True, which
#                              _mode_word hands over as "true".
#   mito_genetic_code        : the genetic code of the mitochondrial run. Cytochrome b, the
#                              target of the QoI fungicides, is encoded in the mitochondrial
#                              genome, and fungal mitochondria read some codons differently
#                              from the nucleus (in code 4, TGA codes for tryptophan, not
#                              for a stop). ChroQueTas takes one code for a whole run, so
#                              rule amr_chroquetas runs it on the nuclear genome with the
#                              code ChroQueTas picks for the species (12 for the yeasts of
#                              the CTG clade, 1 for the others), and, for a species whose
#                              database entry includes cytochrome b (15 of the 57), a second
#                              time on the mitochondrial genome with this one. For the other
#                              species the setting is not used.
#                              auto (the default; an empty value too): code 3 (yeast
#                              mitochondrial) for the Saccharomycetaceae species of the
#                              database, code 4 (mold mitochondrial) for every other species;
#                              the wrapper script chooses it per sample. Or the number of an
#                              NCBI genetic code, used for every sample: 1-6, 9-16 or 21-33
#                              (NCBI_GENETIC_CODES below; 7, 8 and 17 to 20 are not NCBI
#                              codes, and miniprot, which ChroQueTas uses to translate,
#                              refuses them).
#   chroquetas_script_sha256,
#   chroquetas_db_sha256     : the SHA-256 checksums (the fingerprint described at rule
#                              funannotate_db) that the downloaded ChroQueTas.sh and database
#                              archive must have: 64 characters 0-9 and a-f (capitals are
#                              accepted and lowered). The defaults are the checksums of the
#                              two files at the commit links.chroquetas_link pins in the
#                              config template (ChroQueTas v1.0.1), so a copy of the template
#                              that leaves the keys out gets them.
# Rule chroquetas_db reads the two checksums. Rule amr_chroquetas reads the mitochondrial
# code through its params, so a changed code re-runs that step and MultiQC and nothing
# else. AMR_ENABLED decides whether 50_amr.smk defines its two rules, whether
# links.chroquetas_link is required (section 10), whether rule multiqc reads the AMR table
# (90_report.smk) and whether rule all asks for it (section 15).
CHROQUETAS_DEFAULT_SCRIPT_SHA256 = "d8961a0d311be57d4a9939bd9e4d90dbc796589b3298dd1835359ed776abfce3"
CHROQUETAS_DEFAULT_DB_SHA256     = "6f89ecb3117110aad85b8b6ca7b8a3d3715cb5984b6cb75450a6cda5633d16e7"
_amr = _parameter_block(
    "amr", {"chroquetas", "mito_genetic_code", "chroquetas_script_sha256", "chroquetas_db_sha256"}
)
_chroquetas_mode = _mode_word(_amr.get("chroquetas"), "off")
if _chroquetas_mode in {"on", "true", "yes", "1"}:
    AMR_ENABLED = True
elif _chroquetas_mode in {"off", "false", "no", "0", "none", "null"}:
    AMR_ENABLED = False
else:
    sys.exit(f"[FunFlux] parameters.amr.chroquetas: must be on | off (got: {_amr.get('chroquetas')!r})")

# The mitochondrial code, kept as text: "auto", or the number written as text. The word auto
# is accepted in any case. Anything else must be the number of an NCBI genetic code; a YAML
# true/false (an unquoted yes, no, on, off) arrives as the text True/False and is refused
# with the rest.
NCBI_GENETIC_CODES = {
    1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
}
_mito_code_value = _amr.get("mito_genetic_code")
if _mito_code_value is None:
    _mito_code_text = ""
else:
    _mito_code_text = str(_mito_code_value).strip()
_mito_code_message = (
    "[FunFlux] parameters.amr.mito_genetic_code: must be auto or the number of an NCBI genetic code: "
    f"1-6, 9-16 or 21-33 (got: {_mito_code_value!r})"
)
if _mito_code_text == "" or _mito_code_text.lower() == "auto":
    AMR_MITO_GENETIC_CODE = "auto"
else:
    if not re.fullmatch(r"[0-9]+", _mito_code_text):
        sys.exit(_mito_code_message)
    _mito_code_number = int(_mito_code_text)
    if _mito_code_number not in NCBI_GENETIC_CODES:
        sys.exit(_mito_code_message)
    AMR_MITO_GENETIC_CODE = str(_mito_code_number)

CHROQUETAS_SCRIPT_SHA256 = str(
    _amr.get("chroquetas_script_sha256") or CHROQUETAS_DEFAULT_SCRIPT_SHA256
).strip().lower()
if not re.fullmatch(r"[0-9a-f]{64}", CHROQUETAS_SCRIPT_SHA256):
    sys.exit(
        "[FunFlux] parameters.amr.chroquetas_script_sha256: must be the SHA-256 checksum of ChroQueTas.sh, "
        f"64 hexadecimal characters (got: {_amr.get('chroquetas_script_sha256')!r})"
    )
CHROQUETAS_DB_SHA256 = str(
    _amr.get("chroquetas_db_sha256") or CHROQUETAS_DEFAULT_DB_SHA256
).strip().lower()
if not re.fullmatch(r"[0-9a-f]{64}", CHROQUETAS_DB_SHA256):
    sys.exit(
        "[FunFlux] parameters.amr.chroquetas_db_sha256: must be the SHA-256 checksum of the FungAMR database "
        f"archive, 64 hexadecimal characters (got: {_amr.get('chroquetas_db_sha256')!r})"
    )
if AMR_ENABLED:
    print("AMR screen: ChroQueTas on (licence CC BY-NC-ND 4.0: non-commercial use).")
else:
    print("AMR screen: off.")

# ── 9i (continued). The inputs rule amr_chroquetas has in this mode ──
# The nuclear run reads FINAL_CONTIGS in every mode: the nuclear genome after the screen,
# or, in contigs mode without the screen, every contig. Where the mitochondrial genome is
# depends on the run:
#   * with the screen and parameters.mitochondria.mode separate_confirmed (the default), the
#     confirmed mitochondrial contigs are in MITO_CONTIGS and no longer in FINAL_CONTIGS.
#     The wrapper reads MITO_CONTIGS for the mitochondrial run; when that file is empty
#     (no contig was confirmed) it does not run ChroQueTas on it and reports cytochrome b
#     as not evaluated.
#   * without the screen (contigs mode, screen: false), and with the screen in
#     mitochondria mode off or report, nothing is split off: MITO_CONTIGS is absent or
#     empty, and the mitochondrion, if the assembly holds it, is among the contigs of
#     FINAL_CONTIGS. The wrapper then reads FINAL_CONTIGS for the mitochondrial run too.
# AMR_MITO_CONTIGS holds the path of MITO_CONTIGS (still with the {sample} wildcard) in the
# first case and "" in the second. The rule hands it to the wrapper's --mito-contigs; an
# empty one reaches the script as an option without a value, which it reads as "not
# given". AMR_OPTIONAL_INPUTS holds the same file, or nothing: the rule lists it as an
# input, so the job waits for the file where there is one.
if SCREEN and MITO_MODE == "separate_confirmed":
    AMR_MITO_CONTIGS = MITO_CONTIGS
else:
    AMR_MITO_CONTIGS = ""
AMR_OPTIONAL_INPUTS = []
if AMR_MITO_CONTIGS:
    AMR_OPTIONAL_INPUTS.append(AMR_MITO_CONTIGS)


# ═══════════════ 10. Database download links ═══════════════
# Links the workflow downloads from when a database is not already on disk.
#   phix_link         : the PhiX174 genome, which Illumina libraries are spiked with;
#                       read pairs that align to it concordantly are removed.
#                       Short-read modes only (required there; rule download_phix
#                       fetches it as PHIX_FASTA).
#   funannotate_link  : the funannotate database snapshot (a Zenodo record). The
#                       checksum file next to it is found from the same record.
#                       Rule funannotate_db downloads both, as FUNANNOTATEID.tar.gz
#                       and FUNANNOTATEID.sha256, from FUNANNOTATEURL and
#                       FUNANNOTATESHAURL.
#   unite_its_link    : the UNITE SINTAX reference used for ITS classification
#                       (rule unite_db); not needed when directories.unite_db names
#                       a local copy, and ignored, with a NOTE, when both are given.
#   tiberius_model_link : the Tiberius fungi model (rule tiberius_model), checked
#                       against parameters.tiberius.model_sha256 (section 9e).
#                       Required when Tiberius is on, not read when it is off.
#   telomere_db_link  : the telomere repeat database of tidk (rule telomere_db
#                       fetches it as TELOMERE_DB), checked against
#                       parameters.contig_summary.telomere_db_sha256 (section 9h).
#                       Required in every mode: the contig summary runs in all of
#                       them, and the database is never shipped with FunFlux.
#   univec_link       : NCBI's UniVec file (rule univec_illumina keeps its Illumina records
#                       for the contig adapter screen), checked against
#                       parameters.contig_summary.univec_sha256 (section 9h). Optional: a
#                       configuration without the key gets the address of build 10.1 in
#                       NCBI's archive, so a v2.0.0 configuration still runs.
#   chroquetas_link   : the address of one pinned commit of the ChroQueTas repository on
#                       GitHub (raw files), from which rule chroquetas_db downloads three
#                       files: bin/ChroQueTas.sh, FungAMR_db/db-v20250811.tgz and LICENSE.
#                       Required when the AMR screen is on (section 9i), not read when it is
#                       off. A slash at the end is removed, so that the rule can add the
#                       three file paths to it.
if HAS_SHORT_READS:
    PHIXURL = _required_text("links", "phix_link")
else:
    PHIXURL = ""

FUNANNOTATEURL = _required_text("links", "funannotate_link")
FUNANNOTATEID = id_from_url(FUNANNOTATEURL)
ZENODORECORD = zenodo_record_from_url(FUNANNOTATEURL)
FUNANNOTATESHAURL = f"https://zenodo.org/records/{ZENODORECORD}/files/{FUNANNOTATEID}.sha256?download=1"

UNITEURL = _optional_text("links", "unite_its_link")
if UNITEDB and UNITEURL:
    _note("directories.unite_db is set, so links.unite_its_link is not used.")
if not UNITEDB and not UNITEURL:
    sys.exit("[FunFlux] links.unite_its_link: required, unless directories.unite_db names a local UNITE copy.")

if TIBERIUS_ENABLED:
    TIBERIUS_MODEL_URL = _required_text("links", "tiberius_model_link")
else:
    TIBERIUS_MODEL_URL = ""

TELOMERE_DB_URL = _required_text("links", "telomere_db_link")

UNIVEC_DEFAULT_URL = "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/archive/10.1/UniVec"
UNIVEC_URL = _optional_text("links", "univec_link") or UNIVEC_DEFAULT_URL

if AMR_ENABLED:
    CHROQUETAS_URL = _required_text("links", "chroquetas_link").rstrip("/")
else:
    CHROQUETAS_URL = ""


# ═══════════════ 11. Long-read options ═══════════════
# Read only in longread and hybrid mode, from the block NAMED AFTER THE MODE:
# parameters.longread or parameters.hybrid. Both blocks repeat the same keys (the
# hybrid one adds two of its own), so there is one rule for the user: fill in the
# block that has your mode's name. Illumina and contigs mode skip this section; the
# names it defines are set to harmless empty values first, so later sections can
# use them in every mode.
# The settings are checked here, section 13 turns them into one value per sample and
# prints what each sample gets, and the rules of the longread and hybrid front ends,
# of 05_ont_adapters.smk and of 12_medaka_check.smk read the results.
LONG_READ_ASSEMBLER = ""
MEDAKA_MODEL = None
MEDAKA_FALLBACK_AUTO = False
USE_MEDAKA = False
FLYE_ONT_FLAG = ""
HIFIASM_COMMON_ARGS = ""
FILTLONG_MIN_LENGTH = 1000
FILTLONG_KEEP_PERCENT = 90.0
FILTLONG_LENGTH_WEIGHT = 10.0
FILTLONG_TARGET_BASES = ""
FILTLONG_TARGET_ARG = ""
HYBRID_SHORT_READ_GUIDANCE = False
HYBRID_FILTLONG_SPLIT = None
ONT_ADAPTER_SETTING = "off"

if HAS_LONG_READS:
    # The keys the block may hold; hybrid adds its two guidance keys. _lr is the
    # checked block, and _lr_name ("parameters.longread" or "parameters.hybrid") is
    # how the messages below name it.
    _long_read_keys = {"assembler", "flye_input_mode", "medaka_model", "medaka_model_fallback_auto", "filtlong", "hifiasm", "ont_adapter"}
    if IS_HYBRID:
        _long_read_keys = _long_read_keys | {"short_read_guidance", "filtlong_split"}
    _lr = _parameter_block(MODE, _long_read_keys)
    _lr_name = f"parameters.{MODE}"

    # ── The assembler (docs/modes/longread.md) ──
    # hifiasm is the default; Flye is available as an alternative. One assembler
    # is used for the whole run. Both write the same standardised files
    # (assembly.fasta plus the evidence and provenance tables).
    LONG_READ_ASSEMBLER = str(_lr.get("assembler") or "hifiasm").strip().lower()
    if LONG_READ_ASSEMBLER not in {"hifiasm", "flye"}:
        sys.exit(f"[FunFlux] {_lr_name}.assembler: must be hifiasm | flye (got: {_lr.get('assembler')!r})")

    # ── Medaka (ONT only) ──
    # Medaka polishes an ONT assembly with the reads. Its model must match how the
    # reads were basecalled. auto (or empty) works the model out from the FASTQ
    # headers; a model name uses that model; false/off skips polishing. PacBio HiFi
    # reads are accurate enough and never go through Medaka, whatever is set here.
    # MEDAKA_MODEL keeps the value as YAML read it (text, a number, True/False or
    # None; "auto" when the key is missing).
    # medaka_model_fallback_auto lets the Medaka check fall back to auto when a named
    # model is not valid, instead of stopping; the check is explained in
    # workflow/scripts/12_medaka_check/medaka_model_check.py.
    MEDAKA_MODEL = _lr.get("medaka_model", "auto")
    MEDAKA_FALLBACK_AUTO = _config_bool(f"{_lr_name}.medaka_model_fallback_auto", _lr.get("medaka_model_fallback_auto"), False)
    USE_MEDAKA = not medaka_disabled(MEDAKA_MODEL)

    # ── Flye read mode (ONT samples; used only when the assembler is Flye) ──
    # PacBio samples always get Flye --pacbio-hifi. For ONT, auto picks --nano-hq
    # (reads from recent chemistry) unless an explicit Medaka model with "fast" in
    # its name says the reads are older and noisier, which means --nano-raw.
    # This couples two settings: naming a "fast" Medaka model changes the
    # ASSEMBLER's read mode, not only the polisher. Set flye_input_mode explicitly
    # to pin the Flye mode and still choose the Medaka model freely.
    # FLYE_MODE_MAP turns the config word into the Flye option; `{"auto",
    # *FLYE_MODE_MAP}` is the set of allowed words: auto plus the two keys of the map.
    # FLYE_ONT_FLAG is the option every ONT sample gets (section 13).
    FLYE_MODE_MAP = {"nano-raw": "--nano-raw", "nano-hq": "--nano-hq"}
    _flye_mode = str(_lr.get("flye_input_mode") or "auto").strip().lower()
    if _flye_mode == "pacbio-hifi":
        sys.exit(
            f"[FunFlux] {_lr_name}.flye_input_mode: 'pacbio-hifi' is not a value any more. FunFlux decides "
            "per sample from the file name: name PacBio HiFi reads {sample}_pacbio.fastq.gz and they get "
            "Flye --pacbio-hifi automatically."
        )
    if _flye_mode not in {"auto", *FLYE_MODE_MAP}:
        sys.exit(f"[FunFlux] {_lr_name}.flye_input_mode: must be auto | nano-raw | nano-hq (got: {_lr.get('flye_input_mode')!r})")
    if _flye_mode == "auto":
        if has_explicit_medaka_model(MEDAKA_MODEL) and "fast" in str(MEDAKA_MODEL).lower():
            FLYE_ONT_FLAG = "--nano-raw"
        else:
            FLYE_ONT_FLAG = "--nano-hq"
    else:
        FLYE_ONT_FLAG = FLYE_MODE_MAP[_flye_mode]
    if LONG_READ_ASSEMBLER == "hifiasm" and _flye_mode != "auto":
        _note(f"{_lr_name}.flye_input_mode is ignored because the assembler is hifiasm.")

    # ── Filtlong: which long reads to keep before assembly ──
    # Filtlong ranks reads by quality and length and keeps the best share. All of its
    # choices live here so a sample can be tuned without editing workflow code. The
    # defaults are sized for ONT reads; they are provisional for PacBio HiFi, which is
    # already high quality (docs/modes/longread.md).
    #   min_length    : reads shorter than this (bp) are dropped (default 1000).
    #   keep_percent  : the best this-many percent of the bases are kept (default 90).
    #   length_weight : how much read length counts, relative to read quality, when
    #                   reads are ranked (default 10; Filtlong's own default is 1).
    #   target_bases  : an absolute cap on the bases kept (see below).
    _filtlong = _check_block(
        f"{_lr_name}.filtlong", _lr.get("filtlong"),
        {"min_length", "keep_percent", "length_weight", "target_bases"},
    )
    FILTLONG_MIN_LENGTH    = _to_int(f"{_lr_name}.filtlong.min_length", _filtlong.get("min_length", 1000), minimum=1)
    FILTLONG_KEEP_PERCENT  = _to_float(f"{_lr_name}.filtlong.keep_percent", _filtlong.get("keep_percent", 90), minimum=0, minimum_exclusive=True, maximum=100)
    FILTLONG_LENGTH_WEIGHT = _to_float(f"{_lr_name}.filtlong.length_weight", _filtlong.get("length_weight", 10), minimum=0, minimum_exclusive=True)
    # target_bases is an absolute cap on the bases kept (for example 2.4g for a very
    # deep run). Empty (or 0) means no cap beyond keep_percent. Otherwise it must be a
    # number with an optional k, m or g (and an optional b) after it.
    _target_bases = _filtlong.get("target_bases")
    if _is_empty_value(_target_bases) or str(_target_bases).strip() == "0":
        FILTLONG_TARGET_BASES = ""
    else:
        FILTLONG_TARGET_BASES = str(_target_bases).strip()
        if not re.fullmatch(r"\d+(?:\.\d+)?(?:[kKmMgG](?:[bB])?)?", FILTLONG_TARGET_BASES):
            sys.exit(
                f"[FunFlux] {_lr_name}.filtlong.target_bases: must be empty or a size such as "
                f"1500000000 or 1.5g (got: {_target_bases!r})"
            )
    # The cap as Filtlong's command line takes it: "--target_bases 1.5g", or empty when
    # there is no cap. Rule filter_long_reads (longread and hybrid) puts
    # FILTLONG_TARGET_ARG into its command unquoted, so an empty value adds nothing and
    # a set value arrives as the two words Filtlong expects.
    if FILTLONG_TARGET_BASES != "":
        FILTLONG_TARGET_ARG = "--target_bases " + FILTLONG_TARGET_BASES

    # ── ONT adapter trimming (rules/shared/05_ont_adapters.smk; docs/methods/adapter-trimming.md) ──
    # Before Filtlong, fastplong removes the ONT sequencing adapter, and the barcode block that
    # follows it, from the ends of the ONT reads. PacBio HiFi reads are never trimmed.
    #   detect (default; also a missing or empty key, an unquoted null, or true) : rule
    #                 find_ont_adapter works the sequence out from the reads, or decides that
    #                 there is nothing to trim.
    #   off (also false, no, none, and null in quotes) : no trimming; Filtlong reads the raw
    #                 reads, as in v2.0.0.
    #                 These words are tested first: EMPTY_WORDS (section 1b) would count them as
    #                 "not set", which here would mean the default, detect.
    #   a sequence  : 20 to 200 bases of A, C, G and T, trimmed as given (fastplong searches only
    #                 the first and last 200 bp of a read). A short tandem repeat is refused, so
    #                 that a telomere repeat can never be set as the adapter by mistake.
    # ONT_ADAPTER_SETTING keeps "detect", "off" or the sequence in capitals; section 13 turns it
    # into the list of samples that are trimmed (ONT_TRIM_SAMPLES). The checks on a sequence use
    # find_ont_adapter.py itself, loaded as a module the way section 9d loads the AUGUSTUS check
    # (importlib: find the file, make an empty module from it, run the file's code into it).
    _ont_adapter_spec = importlib.util.spec_from_file_location("find_ont_adapter", ONT_ADAPTER_SCRIPT)
    ONT_ADAPTER_MODULE = importlib.util.module_from_spec(_ont_adapter_spec)
    _ont_adapter_spec.loader.exec_module(ONT_ADAPTER_MODULE)

    # A missing key or an unquoted empty value arrives as None (the default, detect); YAML turns
    # an unquoted true/false (also yes/no, on/off) into True/False.
    _ont_adapter = _lr.get("ont_adapter")
    if _ont_adapter is None:
        _ont_word = ""
    else:
        _ont_word = str(_ont_adapter).strip().lower()
    if _ont_word in {"off", "false", "no", "none", "null"}:
        ONT_ADAPTER_SETTING = "off"
    elif _ont_word in {"", "detect", "true", "yes", "on"}:
        ONT_ADAPTER_SETTING = "detect"
    else:
        _ont_sequence = ONT_ADAPTER_MODULE.clean_sequence(_ont_adapter)
        if _ont_sequence is None or not 20 <= len(_ont_sequence) <= 200:
            sys.exit(
                f"[FunFlux] {_lr_name}.ont_adapter: write detect, off, or the sequence to trim "
                f"(20 to 200 bases of A, C, G and T; got: {_ont_adapter!r})"
            )
        if ONT_ADAPTER_MODULE.is_tandem_repeat(_ont_sequence):
            sys.exit(
                f"[FunFlux] {_lr_name}.ont_adapter: {_ont_sequence} is a short unit repeated end to end "
                "(a telomere repeat, for example), not an adapter. Trimming it would cut real genome "
                "from the read ends."
            )
        if ONT_ADAPTER_MODULE.ADAPTER_CORE not in _ont_sequence:
            _note(
                f"{_lr_name}.ont_adapter does not contain the core of the ONT adapter "
                f"({ONT_ADAPTER_MODULE.ADAPTER_CORE}); it is trimmed as given."
            )
        if len(_ont_sequence) < ONT_ADAPTER_MODULE.MIN_SPLIT_LENGTH:
            _note(
                f"{_lr_name}.ont_adapter is shorter than {ONT_ADAPTER_MODULE.MIN_SPLIT_LENGTH} bp, so fastplong "
                "trims it at the read ends only and does not split reads where it occurs in the middle "
                "(a short sequence matches genome by chance)."
            )
        ONT_ADAPTER_SETTING = _ont_sequence

    # ── hifiasm: only the settings that have a clear biological or resource meaning ──
    # Graph-cleaning and overlap switches keep the defaults of the pinned version.
    # The settings are turned, one by one, into the list of hifiasm options
    # _hifiasm_args. --primary makes hifiasm write a primary and an alternate
    # assembly; the other options are explained where they are added.
    _hifiasm = _check_block(
        f"{_lr_name}.hifiasm", _lr.get("hifiasm"),
        {"genome_size", "bloom_filter_bits", "purge_dups", "homozygous_coverage", "telomere_motif"},
    )
    # bloom_filter_bits: hifiasm -f, the size of hifiasm's Bloom filter, a memory
    # table it uses while counting k-mers to set aside the k-mers seen only once. At
    # hifiasm's own default, 37, the table takes 16 GB at the start, and hifiasm's
    # manual advises 0 (no filter) for small genomes. 0 is the FunFlux default, which
    # config/config.yaml gives as right for small fungal genomes. Accepted here:
    # 0 to 63.
    _bloom_bits = _to_int(f"{_lr_name}.hifiasm.bloom_filter_bits", _hifiasm.get("bloom_filter_bits", 0), minimum=0, maximum=63)
    _hifiasm_args = ["--primary", "-f", str(_bloom_bits)]
    # genome_size: a haploid estimate such as 40m (--hg-size). Empty lets hifiasm infer the coverage.
    _genome_size = _hifiasm.get("genome_size")
    if not _is_empty_value(_genome_size) and str(_genome_size).strip() != "0":
        _genome_size = str(_genome_size).strip()
        if not re.fullmatch(r"\d+(?:\.\d+)?(?:[kKmMgG])?", _genome_size):
            sys.exit(f"[FunFlux] {_lr_name}.hifiasm.genome_size: must be empty or a size such as 40000000 or 40m (got: {_genome_size!r})")
        _hifiasm_args += ["--hg-size", _genome_size]
    # purge_dups: how aggressively duplicated haplotypes are purged (-l: 0 = off .. 3).
    # auto (or an empty value, default, none, null) keeps hifiasm's own choice: no -l.
    _purge = _hifiasm.get("purge_dups", "auto")
    _purge_is_auto = _purge is None or str(_purge).strip().lower() in {"", "auto", "default", "none", "null"}
    if not _purge_is_auto:
        _purge_level = _to_int(f"{_lr_name}.hifiasm.purge_dups", _purge, minimum=0, maximum=3)
        _hifiasm_args += ["-l", str(_purge_level)]
    # homozygous_coverage: the read depth of the homozygous peak (--hom-cov). auto (or
    # an empty value, default, none, null) lets hifiasm find it.
    _hom_cov = _hifiasm.get("homozygous_coverage", "auto")
    _hom_cov_is_auto = _hom_cov is None or str(_hom_cov).strip().lower() in {"", "auto", "default", "none", "null"}
    if not _hom_cov_is_auto:
        _hom_cov_depth = _to_int(f"{_lr_name}.hifiasm.homozygous_coverage", _hom_cov, minimum=1)
        _hifiasm_args += ["--hom-cov", str(_hom_cov_depth)]
    # telomere_motif: the telomere repeat as it reads at the 5' end of a chromosome
    # (--telo-m; CCCTAA for the common TTAGGG repeat), a DNA/IUPAC string,
    # upper-cased; never guessed.
    _telomere = _hifiasm.get("telomere_motif")
    if not _is_empty_value(_telomere):
        _telomere = str(_telomere).strip().upper()
        if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", _telomere):
            sys.exit(f"[FunFlux] {_lr_name}.hifiasm.telomere_motif: must be empty or a DNA/IUPAC motif (got: {_telomere!r})")
        _hifiasm_args += ["--telo-m", _telomere]
    # The part of the hifiasm command line that is the same for every sample. The one
    # per-sample difference (--ont for ONT reads) is added in section 13.
    HIFIASM_COMMON_ARGS = " ".join(_hifiasm_args)

    # ── hybrid only: Illumina guidance for Filtlong (docs/methods/filtlong-guidance.md) ──
    # With guidance on, Filtlong scores each ONT read against the Illumina read pairs
    # that map to the screened guide assembly: a long-read stretch that no short read
    # supports scores low. filtlong_split additionally CUTS reads where the short
    # reads stop supporting them (Filtlong --trim --split N, N = the shortest gap
    # that triggers a cut); null scores but never cuts. PacBio HiFi reads are never
    # guided or cut: they are already accurate, and cutting them at Illumina coverage
    # gaps would throw away what HiFi is for.
    #   short_read_guidance : default false. Guidance is off unless asked for, because
    #                         the Illumina reads cover repeats worst, and the long reads
    #                         that span repeats are the ones that hold an assembly
    #                         together: scoring them low, or cutting them there, can break
    #                         the assembly apart. Off, an ONT sample gets the plain length
    #                         and quality filter of longread mode.
    #   filtlong_split      : key absent -> 1000; empty or an "off" word (null, none,
    #                         off, false, no) -> None (no cut); otherwise a whole
    #                         number of at least 1.
    # A split number other than the default 1000 written while guidance is off is
    # ignored, with a NOTE: someone chose that number and expects it to act. The default
    # itself gets no note, because the configuration template writes it out, and every
    # run with guidance off would otherwise print it. Section 13 turns both settings into
    # the Filtlong options of each ONT sample.
    if IS_HYBRID:
        HYBRID_SHORT_READ_GUIDANCE = _config_bool(f"{_lr_name}.short_read_guidance", _lr.get("short_read_guidance"), False)
        if "filtlong_split" not in _lr:
            HYBRID_FILTLONG_SPLIT = 1000
        elif _is_empty_value(_lr["filtlong_split"]):
            HYBRID_FILTLONG_SPLIT = None
        else:
            HYBRID_FILTLONG_SPLIT = _to_int(f"{_lr_name}.filtlong_split", _lr["filtlong_split"], minimum=1)
        if not HYBRID_SHORT_READ_GUIDANCE and HYBRID_FILTLONG_SPLIT not in (None, 1000) and "filtlong_split" in _lr:
            _note(f"{_lr_name}.filtlong_split is ignored because short_read_guidance is false.")

    print(f"Long-read assembler: {LONG_READ_ASSEMBLER}.")


# ── 11b. parameters: keys in the wrong place, and typos ──
# Settings that the v1 FunFluxL workflow kept at the top of `parameters` now live in
# the block named after the mode. A key left at the old place would be silently
# ignored (the setting would fall back to its default), so name where it moved.
# MOVED_LONG_READ_KEYS maps each old key to its new name inside the mode's block.
# Any other top-level key that is not in KNOWN_PARAMETER_KEYS is a typo and stops
# the run too. This check runs in every mode, so a misspelt key is caught even
# when its section is not read. The keys are checked in alphabetical order, so the
# message names the first problem in that order. `genus` is accepted (the v1 key,
# section 8a) but left out of the list the message shows, so that new
# configurations do not start using it.
MOVED_LONG_READ_KEYS = {
    "long_read_assembler": "assembler",
    "flye_input_mode": "flye_input_mode",
    "medaka_model": "medaka_model",
    "filtlong": "filtlong",
    "hifiasm": "hifiasm",
}
KNOWN_PARAMETER_KEYS = {
    "nt_version", "blast_max_target_seqs", "blast_max_hsps", "spades_kmers", "decontamination",
    "mitochondria", "busco", "its", "its_taxonomy_cutoff", "masking_method", "repeatmodeler_quick",
    "repeatmodeler_ltrstruct", "repeatmodeler_seed", "prediction", "genemark", "tiberius", "signalp6",
    "contig_summary", "amr", "longread", "hybrid", "contigs", "genus",
}
for _key in sorted(PARAMETERS):
    if _key in MOVED_LONG_READ_KEYS:
        sys.exit(
            f"[FunFlux] parameters.{_key}: this setting moved into the block named after your mode: "
            f"parameters.longread.{MOVED_LONG_READ_KEYS[_key]} or parameters.hybrid.{MOVED_LONG_READ_KEYS[_key]}."
        )
    if _key not in KNOWN_PARAMETER_KEYS:
        sys.exit(
            f"[FunFlux] parameters.{_key}: unknown setting (a typo?). Known settings: "
            f"{', '.join(sorted(KNOWN_PARAMETER_KEYS - {'genus'}))}."
        )


# ═══════ 12. Sample discovery (per mode), the annotation table and the decontamination tables ═══════
# Only the ACTIVE mode's branch runs: only it touches the file system (it lists its
# own input directory) and only it can stop on "no input" or "bad sample name". The
# other modes' discovery never runs, which is what lets one Snakefile serve four
# kinds of input without an illumina run tripping over an empty longread_dir.
# The sample names come from the input file names. Every branch sets SAMPLES, the
# sorted list of sample names: `rule all` asks for the results of each (section
# 15), the {sample} wildcard may only stand for one of them (section 14), and
# 90_report.smk lists the QC folders of each. The end of this section checks the
# annotation table and the decontamination tables against them.
def _require_input_dir(config_key):
    # Read one input directory from config["input"][config_key] and stop clearly
    # if it is unset or does not exist. Called once per input directory the active
    # mode needs (hybrid needs two, which may be the same directory). The path is
    # returned as written: a relative path is fine, because the working directory
    # never changes (section 3).
    value = _optional_text("input", config_key)
    if not value:
        sys.exit(f"[FunFlux] mode={MODE} requires input.{config_key} to be set in the config.")
    if not os.path.isdir(value):
        sys.exit(f"[FunFlux] input.{config_key}: '{value}' is not an existing directory.")
    return value


def _check_sample_names(samples):
    # Reject sample names that would break paths or wildcards (BAD_CHARS), warn about
    # names the locus-tag convention does not fit (LOCUS_TAG_PATTERN), and echo each
    # accepted sample so the user can confirm the batch from the top of the log.
    # Samples are handled in alphabetical order. In the stop message each bad
    # character is shown in quotes (repr), so that a space is visible as ' '.
    # Called by every branch below.
    for sample in sorted(samples):
        bad = sorted(set(char for char in sample if char in BAD_CHARS))
        if bad:
            sys.exit(
                f"[FunFlux] sample name '{sample}' contains characters that are not allowed "
                f"({' '.join(repr(char) for char in bad)}). Rename the input files."
            )
        if not re.fullmatch(LOCUS_TAG_PATTERN, sample):
            _note(
                f"sample name '{sample}' is not a valid NCBI locus-tag prefix (3 to 12 letters or digits, "
                "starting with a letter). It runs, but the gene IDs funannotate writes cannot be submitted "
                "to GenBank as they are."
            )
        print(f"Sample {sample} will be processed.")


def _require_r2_mates(samples, illumina_dir, short_extn):
    # Illumina input is paired. Discovery starts from the R1 files, so check here that
    # each sample also has its R2 mate, with the same extension, instead of failing
    # deep inside a rule. Called by the illumina and hybrid branches.
    for sample in samples:
        r2_path = os.path.join(illumina_dir, f"{sample}_R2.{short_extn}")
        if not os.path.exists(r2_path):
            sys.exit(f"[FunFlux] input.illumina_dir: sample '{sample}' has an R1 file but no R2 mate: {r2_path}")


def _discover_long_reads(directory):
    # The file-name SUFFIX decides the technology: {sample}_ont.<ext> is Oxford
    # Nanopore, {sample}_pacbio.<ext> is PacBio HiFi (CCS) reads. A sample
    # may exist under one suffix only: both would make the assembler and polishing
    # choices undefined. The ONT files are listed first, so such a sample stops the
    # run when its PacBio file is found. See docs/modes/longread.md.
    #   Takes in: the long-read input directory.
    #   Returns : ({sample: (technology, extension)}, [every extension found]); the
    #             caller checks the extensions with validate_extensions.
    # Called by the longread and hybrid branches.
    found = {}
    extensions = []
    for technology in ("ont", "pacbio"):
        samples, exts = glob_wildcards(os.path.join(directory, "{sample}_" + technology + ".{extn}"))
        for sample, extension in zip(samples, exts):
            if sample in found:
                sys.exit(
                    f"[FunFlux] sample '{sample}' has both _{found[sample][0]} and _{technology} files "
                    f"in {directory}. Provide one long-read technology per sample."
                )
            found[sample] = (technology, extension)
            extensions.append(extension)
    return found, extensions


# Empty starting values, so every name exists in every mode; the branch of the
# active mode below fills in its own. CONTIGS_DIR and CONTIGS_EXTN are read by the
# contigs front end (workflow/rules/contigs/20_assembly.smk) to find each sample's
# FASTA; ILLUMINA_DIR and SHORT_EXTN by rule map_phix (workflow/rules/illumina/10_reads.smk
# and its hybrid copy) to name each sample's raw read pair; LONGREAD_DIR and LONG_EXTN
# only below, to build LONG_FASTQ, the per-sample raw long-read file that
# raw_long_fastq() hands to the long-read rules. READ_TECH is the basis of every
# per-sample decision in section 13.
ILLUMINA_DIR = ""
LONGREAD_DIR = ""
CONTIGS_DIR  = ""
SHORT_EXTN   = ""
LONG_EXTN    = ""
CONTIGS_EXTN = ""
LONG_FASTQ   = {}   # sample -> path of its long-read file (longread, hybrid)
READ_TECH    = {}   # sample -> "ont" or "pacbio"          (longread, hybrid)

if MODE == "illumina":
    # Illumina paired-end: samples are discovered from the R1 files.
    ILLUMINA_DIR = _require_input_dir("illumina_dir")
    _samples, _exts = glob_wildcards(os.path.join(ILLUMINA_DIR, "{sample}_R1.{extn}"))
    SHORT_EXTN = validate_extensions(
        "input.illumina_dir", "{sample}_R1.fastq.gz and {sample}_R2.fastq.gz", _exts, FASTQ_SUFFIXES
    )
    SAMPLES = sorted(set(_samples))
    _check_sample_names(SAMPLES)
    _require_r2_mates(SAMPLES, ILLUMINA_DIR, SHORT_EXTN)

elif MODE == "longread":
    # Long reads only: ONT or PacBio HiFi per sample, from the file-name suffix.
    LONGREAD_DIR = _require_input_dir("longread_dir")
    _found, _exts = _discover_long_reads(LONGREAD_DIR)
    LONG_EXTN = validate_extensions(
        "input.longread_dir", "{sample}_ont.fastq.gz or {sample}_pacbio.fastq.gz", _exts, FASTQ_SUFFIXES
    )
    SAMPLES = sorted(_found)
    _check_sample_names(SAMPLES)
    # _found[sample] is (technology, extension); [0] is the technology.
    for _sample in SAMPLES:
        READ_TECH[_sample] = _found[_sample][0]
        LONG_FASTQ[_sample] = os.path.join(LONGREAD_DIR, f"{_sample}_{READ_TECH[_sample]}.{LONG_EXTN}")

elif MODE == "hybrid":
    # Illumina AND long reads for the same isolates. The two are discovered
    # independently (they may sit in the same directory or in two) and then
    # intersected: a sample needs BOTH. A sample present on one side only stops the
    # run, with the missing files listed, instead of being dropped silently. Both
    # extensions are checked before the two sample lists are compared.
    ILLUMINA_DIR = _require_input_dir("illumina_dir")
    LONGREAD_DIR = _require_input_dir("longread_dir")
    _short_samples, _short_exts = glob_wildcards(os.path.join(ILLUMINA_DIR, "{sample}_R1.{extn}"))
    _found, _long_exts = _discover_long_reads(LONGREAD_DIR)
    SHORT_EXTN = validate_extensions(
        "input.illumina_dir", "{sample}_R1.fastq.gz and {sample}_R2.fastq.gz", _short_exts, FASTQ_SUFFIXES
    )
    LONG_EXTN = validate_extensions(
        "input.longread_dir", "{sample}_ont.fastq.gz or {sample}_pacbio.fastq.gz", _long_exts, FASTQ_SUFFIXES
    )
    _missing_long = sorted(set(_short_samples) - set(_found))
    _missing_short = sorted(set(_found) - set(_short_samples))
    if _missing_long or _missing_short:
        _lines = ["[FunFlux] hybrid mode needs Illumina and long reads for every sample."]
        if _missing_long:
            _lines.append("  No _ont or _pacbio file for: " + ", ".join(_missing_long))
        if _missing_short:
            _lines.append("  No Illumina _R1 file for: " + ", ".join(_missing_short))
        sys.exit("\n".join(_lines))
    SAMPLES = sorted(set(_short_samples))
    _check_sample_names(SAMPLES)
    _require_r2_mates(SAMPLES, ILLUMINA_DIR, SHORT_EXTN)
    # As in longread mode: [0] is the technology of the sample's long-read file.
    for _sample in SAMPLES:
        READ_TECH[_sample] = _found[_sample][0]
        LONG_FASTQ[_sample] = os.path.join(LONGREAD_DIR, f"{_sample}_{READ_TECH[_sample]}.{LONG_EXTN}")

else:
    # contigs: a genome you already have, one FASTA per sample. glob_wildcards is
    # GREEDY on {sample}, so the LAST dot splits off the extension: "my.genome.fasta"
    # is read as sample "my.genome" (the dot stays in the sample name).
    CONTIGS_DIR = _require_input_dir("contigs_dir")
    _samples, _exts = glob_wildcards(os.path.join(CONTIGS_DIR, "{sample}.{extn}"))
    CONTIGS_EXTN = validate_extensions("input.contigs_dir", "{sample}.fasta", _exts, FASTA_SUFFIXES)
    SAMPLES = sorted(set(_samples))
    _check_sample_names(SAMPLES)


# ── The annotation table: one row per sample ──
# files.annotation_params is a tab-separated table with three or four columns:
#   sample  species  proteins  [augustus_species]
# species is written into the GenBank/GFF files; proteins is a FASTA of proteins from
# a closely related fungus (funannotate aligns them to the genome as gene evidence);
# augustus_species, if given, names a trained AUGUSTUS model of the funannotate
# database, which AUGUSTUS then uses as it is. Left out, left empty or written "-",
# it means: train AUGUSTUS on this genome's BUSCO genes (section 9d explains both).
# Lines starting with # and blank lines are ignored.
def parse_annotation_file(annotation_file):
    # Read the table into {sample: {"species_name", "protein_file", "augustus_model"}};
    # augustus_model is "" for a sample whose AUGUSTUS is to be trained.
    # Stops when the file is missing, when a row does not have three or four
    # tab-separated columns, or when a sample has two rows. It does not look at the
    # protein files or the model names; the checks after the lookup functions below
    # do that, for the samples of this run only. Called once, just below (ANNOTATIONS).
    if not os.path.isfile(annotation_file):
        sys.exit(f"[FunFlux] files.annotation_params: '{annotation_file}' not found.")
    annotations = {}
    with open(annotation_file, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
            # Skip blank lines and comment lines; strip() also removes spaces and tabs
            # at the ends of the line.
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Split the row at its tabs and remove the spaces around each value.
            fields = [field.strip() for field in line.split("\t")]
            # strip() above has removed a trailing tab, so a row whose fourth column is
            # empty arrives here with three columns.
            if len(fields) not in (3, 4):
                sys.exit(
                    f"[FunFlux] files.annotation_params: line {line_number} has {len(fields)} tab-separated "
                    "columns; expected 3 or 4 (sample, species, proteins and, optionally, augustus_species)."
                )
            sample, species, protein_file = fields[:3]
            augustus_model = fields[3] if len(fields) == 4 else ""
            # "-" is the visible way of leaving the model out.
            if augustus_model == "-":
                augustus_model = ""
            # Absolute, like every other user-supplied path (see the databases in section 3).
            protein_file = os.path.abspath(protein_file)
            if sample in annotations:
                sys.exit(f"[FunFlux] files.annotation_params: sample '{sample}' appears more than once (line {line_number}).")
            annotations[sample] = {
                "species_name": species,
                "protein_file": protein_file,
                "augustus_model": augustus_model,
            }
    return annotations


ANNOTATIONS = parse_annotation_file(ANNOTATION)

# The annotation table is read through these small named functions, one value each, so a
# rule can say "the species of this job's sample" without a lambda (a nameless one-line
# function). Snakemake calls such a function once per job and hands it that job's
# wildcards (wildcards.sample). species_name and protein_evidence give rule
# funannotate_prediction its -s and --protein_evidence values; species_name,
# augustus_model and augustus_optimize_for_sample give rule augustus_model_check what to
# check, and augustus_arguments gives the prediction its AUGUSTUS options.
def species_name(wildcards):
    return ANNOTATIONS[wildcards.sample]["species_name"]


def protein_evidence(wildcards):
    return ANNOTATIONS[wildcards.sample]["protein_file"]


def augustus_model(wildcards):
    # The model name of column 4, or "" when AUGUSTUS is to be trained on this genome.
    return ANNOTATIONS[wildcards.sample]["augustus_model"]


def augustus_optimize_for_sample(wildcards):
    # "yes" when parameters.prediction.augustus_optimize is on AND this sample's column 4 is
    # empty (the only samples it applies to), "no" otherwise. Rule augustus_model_check
    # takes it as a param. It is worked out per sample on purpose: Snakemake re-runs a job
    # when one of its params changes, so switching the setting re-runs the check, and with
    # it the gene prediction, only for the samples whose AUGUSTUS is trained.
    if PREDICTION_AUGUSTUS_OPTIMIZE and not ANNOTATIONS[wildcards.sample]["augustus_model"]:
        return "yes"
    return "no"


def augustus_arguments(wildcards):
    # The AUGUSTUS options of this job's funannotate predict command (section 9d):
    #   a model named : --augustus_species <name>   (the stored model, no training)
    #   no model      : nothing, so funannotate trains AUGUSTUS on the BUSCO genes,
    #                   plus --optimize_augustus when parameters.prediction.
    #                   augustus_optimize is true
    model = ANNOTATIONS[wildcards.sample]["augustus_model"]
    if model:
        return f"--augustus_species {model}"
    if PREDICTION_AUGUSTUS_OPTIMIZE:
        return "--optimize_augustus"
    return ""


def funannotate_results_stem(wildcards):
    # funannotate names every file of its gene prediction <Genus_species>_<strain>.<ext>,
    # where the strain is the sample name (the rules pass --strain {sample}) and the
    # species has its spaces turned into underscores. antiSMASH, eggNOG and SignalP read
    # two of these files and the contig summary a third, so the stem is worked out once
    # here: for the example row of config/annotation_parameters.tsv (sample ARSEF3097,
    # Beauveria bassiana) it is "Beauveria_bassiana_ARSEF3097". Called by
    # funannotate_proteins, funannotate_gbk and funannotate_gff3.
    species = ANNOTATIONS[wildcards.sample]["species_name"].replace(" ", "_")
    return f"{species}_{wildcards.sample}"


def funannotate_proteins(wildcards):
    # The predicted proteins (FASTA) of this job's sample, inside the folder rule
    # funannotate_prediction writes: the input of eggNOG and SignalP 6 (rules
    # eggnog_annotation and signalp6_annotation read it through their params).
    prediction_dir = FUNANNOTATE_DIR.format(sample=wildcards.sample)
    file_stem = funannotate_results_stem(wildcards)
    return prediction_dir + "/predict_results/" + file_stem + ".proteins.fa"


def genemark_gate(wildcards):
    # The file that proves GeneMark-ES started on a test slice of this sample's genome.
    # funannotate_prediction lists it as an input when GeneMark is switched on, so an
    # expired licence key stops the run before gene prediction begins instead of well
    # into it. With GeneMark off it returns an empty list, which Snakemake
    # reads as "no input here", so nothing waits for it.
    if GENEMARK_ENABLED:
        return GENEMARK_CHECK.format(sample=wildcards.sample)
    return []


def tiberius_gate(wildcards):
    # Tiberius's gene models of this sample, converted for EVidenceModeler (rule
    # tiberius_gff3). funannotate_prediction lists them as an input when Tiberius is on,
    # so the prediction waits for them. With Tiberius off it returns an empty list
    # ("no input here"), as genemark_gate does.
    if TIBERIUS_ENABLED:
        return TIBERIUS_GFF3.format(sample=wildcards.sample)
    return []


def tiberius_arguments(wildcards):
    # The option of this job's funannotate predict command that hands it Tiberius's gene
    # models: --other_gff <file>:<weight>. funannotate passes the file on to
    # EVidenceModeler as "other_pred1", a source of the kind EVidenceModeler calls
    # OTHER_PREDICTION, with that weight (TIBERIUS_EVM_WEIGHT, section 9e). Empty text
    # when Tiberius is off, so the command is then the one without Tiberius.
    if TIBERIUS_ENABLED:
        return f"--other_gff {TIBERIUS_GFF3.format(sample=wildcards.sample)}:{TIBERIUS_EVM_WEIGHT}"
    return ""


def funannotate_gbk(wildcards):
    # The predicted genes of this job's sample in GenBank format, from the same
    # folder: the input of antiSMASH (rule antismash_annotation, through its params).
    prediction_dir = FUNANNOTATE_DIR.format(sample=wildcards.sample)
    file_stem = funannotate_results_stem(wildcards)
    return prediction_dir + "/predict_results/" + file_stem + ".gbk"


def funannotate_gff3(wildcards):
    # The final annotation of this job's sample in GFF3 format, in the annotate_results
    # folder that funannotate annotate writes inside the same prediction folder (rule
    # funannotate_annotation): the input of rule contig_summary (25_contig_summary.smk),
    # which counts the protein-coding genes (the mRNA features) on each scaffold. Reached
    # through that rule's params like the two files above, because the folder is the
    # directory() output of funannotate_prediction; the rule lists ANNOTATION_COMPLETE, the
    # marker funannotate_annotation writes last, as its input instead, so that the file is
    # there when the job starts.
    prediction_dir = FUNANNOTATE_DIR.format(sample=wildcards.sample)
    file_stem = funannotate_results_stem(wildcards)
    return prediction_dir + "/annotate_results/" + file_stem + ".gff3"


# Every sample of this batch needs a row, and the protein evidence it names must
# exist and not be empty. Only the samples being run are checked, so a leftover row
# for a sample that is not in this batch cannot stop the run.
for _sample in SAMPLES:
    if _sample not in ANNOTATIONS:
        sys.exit(f"[FunFlux] files.annotation_params: sample '{_sample}' has no row in '{ANNOTATION}'.")
    _protein_file = ANNOTATIONS[_sample]["protein_file"]
    if not os.path.isfile(_protein_file):
        sys.exit(f"[FunFlux] files.annotation_params: protein file '{_protein_file}' of sample '{_sample}' does not exist.")
    if os.path.getsize(_protein_file) == 0:
        sys.exit(f"[FunFlux] files.annotation_params: protein file '{_protein_file}' of sample '{_sample}' is empty.")

# The AUGUSTUS model of each sample of this batch (section 9d). When the database is
# already on disk (it has its index file), check_augustus_model.py decides now what
# funannotate will do, exactly as rule augustus_model_check will, and a named model that
# is not there stops the run. When the database is not on disk yet, only the characters
# of a name can be checked now; the rule checks the rest right after the download. One
# line then says what AUGUSTUS will do for each sample.
_funannotate_db_on_disk = os.path.isfile(os.path.join(FUNANNOTATEDB, "funannotate-db-info.txt"))
_augustus_plan = []
_samples_with_model = []
for _sample in SAMPLES:
    _model = ANNOTATIONS[_sample]["augustus_model"]
    if _model:
        _samples_with_model.append(_sample)
    if not _funannotate_db_on_disk:
        if _model:
            _problem = AUGUSTUS_CHECK.model_problem(_sample, _model, None)
            if _problem:
                sys.exit(_problem)
            _augustus_plan.append(f"{_sample} names the stored model {_model} (checked once the database is unpacked)")
        # Even an empty column 4 is decided only then: funannotate would use a model the
        # database holds under the sample's own name instead of training.
        elif PREDICTION_AUGUSTUS_OPTIMIZE:
            _augustus_plan.append(f"{_sample} is trained on its BUSCO genes, with the optimisation (checked once the database is unpacked)")
        else:
            _augustus_plan.append(f"{_sample} is trained on its BUSCO genes (checked once the database is unpacked)")
        continue
    _problem, _mode, _model_used = AUGUSTUS_CHECK.decide(
        _sample, ANNOTATIONS[_sample]["species_name"], _model, FUNANNOTATEDB, PREDICTION_AUGUSTUS_OPTIMIZE
    )
    if _problem:
        sys.exit(_problem)
    if _mode == AUGUSTUS_CHECK.MODE_STORED:
        _augustus_plan.append(f"{_sample} uses the stored model {_model_used}")
    elif _mode == AUGUSTUS_CHECK.MODE_OWN_STORED:
        _augustus_plan.append(f"{_sample} uses the model {_model_used}, stored earlier under its own name")
        # funannotate looks for this name before it trains; someone added an earlier
        # result to the database with `funannotate species`.
        _note(
            f"sample '{_sample}' leaves the AUGUSTUS model empty, but the funannotate database holds a model "
            f"under the name funannotate gives this sample's own training, '{_model_used}', and funannotate "
            f"will use it as it is instead of training. Delete '{FUNANNOTATEDB}/trained_species/{_model_used}' "
            "to train afresh, or name that model in column 4 to say so."
        )
    elif _mode == AUGUSTUS_CHECK.MODE_TRAINED_OPTIMIZED:
        _augustus_plan.append(f"{_sample} is trained on its BUSCO genes, with the optimisation")
    else:
        _augustus_plan.append(f"{_sample} is trained on its BUSCO genes")
if SAMPLES:
    print("AUGUSTUS: " + "; ".join(_augustus_plan) + ".")
if PREDICTION_AUGUSTUS_OPTIMIZE and _samples_with_model:
    _note(
        "parameters.prediction.augustus_optimize is on, but it applies only to samples whose AUGUSTUS is "
        f"trained; {', '.join(_samples_with_model)} name(s) a stored model in files.annotation_params, which "
        "AUGUSTUS uses as it is."
    )


# ── The decontamination tables (the check is defined in section 8a) ──
# Called here, not in section 8a, because it checks include mode sample by sample,
# and SAMPLES is known only now. Only when the screen runs: without it there is no
# rule select_contigs, and the tables are never read.
if SCREEN:
    _check_decontam_tables(SAMPLES)


# ═══════════ 13. Per-sample technology decisions ═══════════
# A longread or hybrid batch can hold S1_ont.fastq.gz next to S2_pacbio.fastq.gz. Six
# things then differ per sample, and two of them (adapter trimming and Medaka) change the
# job graph (the set of jobs Snakemake plans: a PacBio sample has no trimming and no
# Medaka jobs at all):
#   * the Flye read-mode flag                 (--nano-hq / --nano-raw  vs  --pacbio-hifi)
#   * the hifiasm technology flag             (--ont  vs  nothing: native HiFi mode)
#   * the minimap2 preset                     (map-ont  vs  map-hifi)
#   * whether the adapters are trimmed        (ONT only, unless ont_adapter is off)
#   * whether Filtlong is Illumina-guided     (ONT in a hybrid run with guidance on only)
#   * whether Medaka runs                     (ONT only)
# READ_TECH is the one source of truth; this section turns it into one value per
# sample, with one assignment per line, and prints the result so the log shows what
# each sample will get. The rules read the values through the small named functions
# below the loop.
# This serves the longread and hybrid modes. In illumina and contigs mode READ_TECH
# is empty and nothing here does anything.
#   ONT_SAMPLES, PACBIO_SAMPLES : the samples of each technology, in SAMPLES order.
#   MEDAKA_SAMPLES              : the ONT samples when Medaka is on, otherwise none.
#   ONT_TRIM_SAMPLES            : the ONT samples whose reads fastplong trims before Filtlong
#                                 (all of them unless parameters.<mode>.ont_adapter is off).
ONT_SAMPLES    = [sample for sample in SAMPLES if READ_TECH.get(sample) == "ont"]
PACBIO_SAMPLES = [sample for sample in SAMPLES if READ_TECH.get(sample) == "pacbio"]
MEDAKA_SAMPLES = ONT_SAMPLES if USE_MEDAKA else []
if ONT_ADAPTER_SETTING == "off":
    ONT_TRIM_SAMPLES = []
else:
    ONT_TRIM_SAMPLES = ONT_SAMPLES

# One dictionary per decision: sample name -> value.
FLYE_READ_FLAG = {}
HIFIASM_TECH_FLAG = {}
MINIMAP_PRESET = {}
FILTLONG_GUIDANCE_FLAGS = {}
POLISH_SOURCE = {}
FINALIZE_SOURCE = {}

# Only the long-read modes have anything to decide here; illumina and contigs mode
# loop over an empty list.
LONG_READ_SAMPLES = SAMPLES if HAS_LONG_READS else []

for _sample in LONG_READ_SAMPLES:
    # The tool flags of the sample's technology. PacBio HiFi: Flye --pacbio-hifi,
    # hifiasm in its native HiFi mode, minimap2 preset map-hifi, Filtlong never
    # guided. ONT: the Flye mode chosen in section 11, hifiasm --ont, map-ont, and in
    # a hybrid run with guidance on, the Filtlong options that score the reads
    # against this sample's selected Illumina pairs (-1 and -2 name the two read
    # files; --trim --split N is added when filtlong_split is set).
    if READ_TECH[_sample] == "pacbio":
        FLYE_READ_FLAG[_sample] = "--pacbio-hifi"
        HIFIASM_TECH_FLAG[_sample] = ""
        MINIMAP_PRESET[_sample] = "map-hifi"
        FILTLONG_GUIDANCE_FLAGS[_sample] = ""
    else:
        FLYE_READ_FLAG[_sample] = FLYE_ONT_FLAG
        HIFIASM_TECH_FLAG[_sample] = "--ont"
        MINIMAP_PRESET[_sample] = "map-ont"
        FILTLONG_GUIDANCE_FLAGS[_sample] = ""
        if IS_HYBRID and HYBRID_SHORT_READ_GUIDANCE:
            _guide = "-1 " + SEL_R1.format(sample=_sample) + " -2 " + SEL_R2.format(sample=_sample)
            if HYBRID_FILTLONG_SPLIT is not None:
                _guide += f" --trim --split {HYBRID_FILTLONG_SPLIT}"
            FILTLONG_GUIDANCE_FLAGS[_sample] = _guide

    # What the long-read assembly looks like once the long reads have polished it:
    # the Medaka consensus for ONT samples that use Medaka, otherwise the
    # decontaminated assembly as it is. POLISH_SOURCE and FINALIZE_SOURCE hold paths
    # that still contain {sample}; polish_source() and finalize_source() below fill
    # it in.
    if _sample in MEDAKA_SAMPLES:
        POLISH_SOURCE[_sample] = MEDAKA_CONSENSUS
    else:
        POLISH_SOURCE[_sample] = DECONTAM_CONTIGS

    # What becomes contigs_all.fasta: in longread mode the long-read-polished
    # assembly; in hybrid mode that assembly after Polypolish has corrected it with the
    # Illumina reads.
    if IS_HYBRID:
        FINALIZE_SOURCE[_sample] = POLYPOLISH_CONTIGS
    else:
        FINALIZE_SOURCE[_sample] = POLISH_SOURCE[_sample]

    # One line per sample for the log, for example
    #   [FunFlux] SAM1: ont reads -> hifiasm --ont, minimap2 map-ont, Medaka.
    # First the assembler with its technology flag, and whether Medaka runs.
    if LONG_READ_ASSEMBLER == "hifiasm":
        _assembler_text = "hifiasm" if READ_TECH[_sample] == "pacbio" else "hifiasm --ont"
    else:
        _assembler_text = "Flye " + FLYE_READ_FLAG[_sample]
    _medaka_text = "Medaka" if _sample in MEDAKA_SAMPLES else "no Medaka"

    # In a hybrid run, also say whether Filtlong is guided by the Illumina reads. That
    # is the one hybrid decision that differs between samples (ONT yes, PacBio never)
    # and that a user cannot otherwise see without reading the config. Longread mode
    # has no such choice, so its line stays shorter.
    if IS_HYBRID and FILTLONG_GUIDANCE_FLAGS[_sample]:
        _filtlong_text = "Filtlong guided by the Illumina reads, "
    elif IS_HYBRID:
        _filtlong_text = "Filtlong not guided, "
    else:
        _filtlong_text = ""

    # ONT samples whose reads are trimmed first say so at the front of the line.
    if _sample in ONT_TRIM_SAMPLES:
        _trim_text = "adapter check + fastplong, "
    else:
        _trim_text = ""
    print(f"[FunFlux] {_sample}: {READ_TECH[_sample]} reads -> {_trim_text}{_filtlong_text}{_assembler_text}, minimap2 {MINIMAP_PRESET[_sample]}, {_medaka_text}.")

# A batch with PacBio samples gets one NOTE: _pacbio must mean HiFi (CCS) reads
# (docs/modes/longread.md).
if PACBIO_SAMPLES:
    _note("samples with a _pacbio file are assumed to hold PacBio HiFi (CCS) reads: " + ", ".join(PACBIO_SAMPLES)
          + ". Continuous long reads (CLR) are not supported by this workflow.")

# A sequence given in parameters.<mode>.ont_adapter has nothing to act on in a batch without ONT
# reads; say so, as for the Medaka model below.
if HAS_LONG_READS and not ONT_SAMPLES and ONT_ADAPTER_SETTING not in ("detect", "off"):
    _note(f"no ONT samples in this batch: parameters.{MODE}.ont_adapter is ignored.")

# One line about Medaka for the whole run. If any sample uses Medaka: which model, or
# that it will be worked out from the reads. Otherwise, in a long-read run: why
# Medaka is not used (switched off, or no ONT sample in the batch). Illumina and
# contigs runs print nothing here.
if MEDAKA_SAMPLES:
    if has_explicit_medaka_model(MEDAKA_MODEL):
        print(f"Medaka model: {MEDAKA_MODEL}.")
    else:
        print("Medaka model will be worked out from the FASTQ headers (checked before assembly).")
elif HAS_LONG_READS:
    if ONT_SAMPLES and not USE_MEDAKA:
        print("Medaka disabled in the config. Skipping consensus polishing.")
    elif not ONT_SAMPLES:
        print("No ONT samples in this batch: any medaka_model setting is ignored.")


# The rules ask for these per-sample values through one-line functions. Snakemake
# calls such a function once per job with the job's wildcards. They are for the
# long-read rules: the longread and hybrid front ends, rule check_medaka_model in
# 12_medaka_check.smk, and the long-read version of rule map_contigs in
# 10_decontam.smk:
#   raw_long_fastq    : the sample's long-read input file
#   filtlong_input    : the reads Filtlong starts from (trimmed for ONT samples, otherwise raw)
#   read_tech         : ont or pacbio, the sample's technology (for the provenance table)
#   flye_read_flag    : --pacbio-hifi, --nano-hq or --nano-raw
#   hifiasm_tech_flag : --ont, or nothing for HiFi reads
#   minimap_preset    : map-ont or map-hifi
#   filtlong_guidance : the Filtlong options for Illumina guidance, or nothing
#   filtlong_guide_reads : the selected Illumina pairs those options name, or nothing
#   polish_source     : the assembly after long-read polishing (Medaka or not)
#   finalize_source   : the file that becomes contigs_all.fasta
#   medaka_gate       : the Medaka model check the assembler waits for, if any
#   trimmed_adapter_sequence, trimmed_adapter_inputs, trimmed_adapter_choice,
#   trimmed_adapter_check : the read-level files of a trimmed ONT sample that the contig
#                       adapter screen uses (rules blast_adapters and screen_adapters in
#                       25_contig_summary.smk), or nothing for any other sample
def raw_long_fastq(wildcards):
    return LONG_FASTQ[wildcards.sample]


def filtlong_input(wildcards):
    # The reads rule filter_long_reads (longread and hybrid) starts from: for an ONT sample
    # that is trimmed, fastplong's output (TRIMMED_LONG, rule trim_ont_adapters in
    # 05_ont_adapters.smk); for a PacBio sample, or with parameters.<mode>.ont_adapter off,
    # the raw file. The raw NanoPlot report and the Medaka model check keep reading the raw
    # file: they describe the reads as delivered.
    if wildcards.sample in ONT_TRIM_SAMPLES:
        return TRIMMED_LONG.format(sample=wildcards.sample)
    return LONG_FASTQ[wildcards.sample]


def read_tech(wildcards):
    return READ_TECH[wildcards.sample]


def flye_read_flag(wildcards):
    return FLYE_READ_FLAG[wildcards.sample]


def hifiasm_tech_flag(wildcards):
    return HIFIASM_TECH_FLAG[wildcards.sample]


def minimap_preset(wildcards):
    return MINIMAP_PRESET[wildcards.sample]


def filtlong_guidance(wildcards):
    return FILTLONG_GUIDANCE_FLAGS[wildcards.sample]


def filtlong_guide_reads(wildcards):
    # The files behind the -1 and -2 of a guided Filtlong command: the sample's own
    # selected Illumina pairs (SEL_R1, SEL_R2; rule map_sel_contigs in
    # rules/hybrid/25_guide.smk). Rule filter_long_reads of the hybrid front end lists
    # them as an input, so Filtlong waits for them. A sample whose command carries no
    # guidance (PacBio HiFi, or guidance off) gets an empty list and waits for nothing:
    # its Filtlong and assembly jobs do not wait for the guide screen, whose result they
    # do not use.
    if FILTLONG_GUIDANCE_FLAGS[wildcards.sample]:
        return [SEL_R1.format(sample=wildcards.sample), SEL_R2.format(sample=wildcards.sample)]
    return []


def polish_source(wildcards):
    return POLISH_SOURCE[wildcards.sample].format(sample=wildcards.sample)


def finalize_source(wildcards):
    return FINALIZE_SOURCE[wildcards.sample].format(sample=wildcards.sample)


def trimmed_adapter_sequence(wildcards):
    # The sequences fastplong trimmed from the reads of an ONT sample (ONT_ADAPTER_SEQUENCE, rule
    # find_ont_adapter in 05_ont_adapters.smk), as a list of one file, which rule blast_adapters
    # adds to the adapters it searches the contigs for. Any other sample (PacBio, Illumina,
    # contigs mode, or trimming off) gets an empty list, which the rule's command writes as
    # nothing at all.
    if wildcards.sample in ONT_TRIM_SAMPLES:
        return [ONT_ADAPTER_SEQUENCE.format(sample=wildcards.sample)]
    return []


def trimmed_adapter_inputs(wildcards):
    # The choice and check files of rules find_ont_adapter and check_ont_trimming that rule
    # screen_adapters reports in the Adapters table, for an ONT sample whose reads were trimmed.
    # Listed as inputs so that the screen waits for them; any other sample gets an empty list.
    if wildcards.sample in ONT_TRIM_SAMPLES:
        return [
            ONT_ADAPTER_CHOICE.format(sample=wildcards.sample),
            ONT_ADAPTER_CHECK.format(sample=wildcards.sample),
        ]
    return []


# The same two files one by one, for the params of rule screen_adapters: the path, or an empty
# text. Snakemake writes an empty text as nothing at all, so the script then sees the option
# without a value and reads it as "not given".


def trimmed_adapter_choice(wildcards):
    if wildcards.sample in ONT_TRIM_SAMPLES:
        return ONT_ADAPTER_CHOICE.format(sample=wildcards.sample)
    return ""


def trimmed_adapter_check(wildcards):
    if wildcards.sample in ONT_TRIM_SAMPLES:
        return ONT_ADAPTER_CHECK.format(sample=wildcards.sample)
    return ""


def medaka_gate(wildcards):
    # The file that proves the Medaka model was checked (MEDAKA_MODEL_RESOLVED, written
    # by rule check_medaka_model). Rule long_read_assembly lists it as an input for ONT
    # samples that use Medaka, so a wrong model name stops the run BEFORE a
    # multi-hour assembly. Any other sample gets an empty list and never waits.
    if wildcards.sample in MEDAKA_SAMPLES:
        return MEDAKA_MODEL_RESOLVED.format(sample=wildcards.sample)
    return []


# ═════════════ 14. Global {sample} wildcard constraint ═════════════
# Pin the {sample} wildcard to the exact set of discovered sample names. Without
# this, {sample} matches anything (even "/"), which lets one rule's output pattern
# match a file that belongs to another rule. For example, in hybrid mode the guide
# assembly's 02.assembly/S1/guide/contaminants/contigs_decontam.fasta also fits the
# output pattern of rule select_contigs, with sample = "S1/guide". Restricting
# {sample} to real names removes that whole class of accidental matches in one place.
# The constraint is a regular expression: the sample names joined by | ("or"), each
# passed through re.escape() so that a dot in a name only matches a dot.
# SAMPLES is never empty here (discovery stops when it finds nothing); the `.+`
# (any name) used when it is empty only keeps the directive valid in any case.
if SAMPLES:
    _SAMPLE_CONSTRAINT = "|".join(re.escape(s) for s in SAMPLES)
else:
    _SAMPLE_CONSTRAINT = r".+"


# A Snakemake directive, declared outside any rule, so it applies to every rule of
# the workflow.
wildcard_constraints:
    sample=_SAMPLE_CONSTRAINT


# ═══════════════════ 15. Targets ═══════════════════
# all_targets() is what `rule all` asks for (at the end of workflow/Snakefile). It
# is the one place 00_common has to know the names of the final outputs: the files a
# finished run delivers. Snakemake works backwards from them and runs every rule
# needed to make them. download_targets(), at the end of this section, is what
# `rule download_all` asks for: only the files the download steps make.
#
# EMPTY RULE FOLDER: asking Snakemake to build a file that no rule produces makes even a
# dry run fail. Every mode has its rule files in rules/<MODE>/, so a run plans the
# mode's whole job list. Should that directory ever hold no .smk file (a checkout with
# a front end removed, or a new mode whose rules are still being written),
# all_targets() returns an empty list instead, so `snakemake -n` plans one empty
# `all` job (Snakemake 9.14.6 lists it under "neither input nor output") rather than
# failing, and the configuration can still be checked.
def _rule_modules_present():
    # True when rules/<MODE>/ holds at least one .smk file. Called by all_targets()
    # below.
    rule_file_pattern = os.path.join(WORKFLOW_DIR, "rules", MODE, "*.smk")
    rule_files = glob.glob(rule_file_pattern)
    return len(rule_files) > 0


def _frontend_targets():
    # The files each mode's front end promises, one block per mode: the genome files
    # the mode's assembly chain ends with, the kept read and assembly files a user looks
    # at (the filtered long reads, the selected Illumina pairs, the polished assemblies,
    # the assembler's provenance table) and the QC results MultiQC reads (the fastp
    # report, the NanoPlot folders, the Qualimap folders). Rule multiqc asks for the QC
    # results too; they are listed here as well, as in v1's rule all, so that a run that
    # leaves the report out (--omit-from multiqc) still makes them. The intermediate
    # files (the trimmed reads, the draft, the BAM, ...) are not listed; Snakemake
    # builds them because the listed files need them. expand() fills {sample} with each
    # name in SAMPLES and returns one path per sample. Called by all_targets() below.
    targets = []
    if MODE == "contigs":
        # With the screen off, filter_contigs writes the final genome directly. With it
        # on, the genome is split into nuclear and mitochondrial files, and the screen's
        # table is delivered next to them.
        targets += expand(FINAL_CONTIGS, sample=SAMPLES)
        if CONTIGS_SCREEN:
            targets += expand(MITO_CONTIGS, sample=SAMPLES)
            targets += expand(MITO_SCREEN_TSV, sample=SAMPLES)
    elif MODE == "illumina":
        # The screen always runs, so the three genome files always exist. The fastp
        # report (rule trim_adapters) and Qualimap on the trimmed reads mapped to the
        # draft (rule map_evaluation) are the QC results of this mode.
        targets += expand(FINAL_CONTIGS, sample=SAMPLES)
        targets += expand(MITO_CONTIGS, sample=SAMPLES)
        targets += expand(MITO_SCREEN_TSV, sample=SAMPLES)
        targets += expand(FASTP_JSON, sample=SAMPLES)
        targets += expand(QUALIMAP_DIR, sample=SAMPLES)
    elif MODE == "longread":
        # The three genome files; the filtered long reads, which are kept (rule
        # filter_long_reads); NanoPlot on the raw and on the filtered reads; the
        # assembler's provenance table (rule long_read_assembly); Qualimap on the
        # filtered reads mapped to the draft; and the Medaka consensus of the ONT samples
        # that use Medaka (MEDAKA_SAMPLES; the list is empty when Medaka is off or the
        # batch holds no ONT sample, and expand() then adds nothing).
        targets += expand(FINAL_CONTIGS, sample=SAMPLES)
        targets += expand(MITO_CONTIGS, sample=SAMPLES)
        targets += expand(MITO_SCREEN_TSV, sample=SAMPLES)
        targets += expand(FILTERED_LONG, sample=SAMPLES)
        targets += expand(NANOPLOT_RAW_DIR, sample=SAMPLES)
        targets += expand(NANOPLOT_FILTERED_DIR, sample=SAMPLES)
        targets += expand(ASSEMBLY_PROVENANCE, sample=SAMPLES)
        targets += expand(QUALIMAP_DIR, sample=SAMPLES)
        targets += expand(MEDAKA_CONSENSUS, sample=MEDAKA_SAMPLES)
    else:
        # hybrid: everything longread mode delivers, plus the fastp report of the
        # Illumina reads, the Illumina pairs selected on the screened guide assembly
        # (SEL_R1/R2: Polypolish needs them whether or not Filtlong is guided), the
        # Polypolish result (rule short_read_correction) and Qualimap on the Illumina
        # reads mapped to the guide (rule map_guide_evaluation).
        targets += expand(FINAL_CONTIGS, sample=SAMPLES)
        targets += expand(MITO_CONTIGS, sample=SAMPLES)
        targets += expand(MITO_SCREEN_TSV, sample=SAMPLES)
        targets += expand(FILTERED_LONG, sample=SAMPLES)
        targets += expand(NANOPLOT_RAW_DIR, sample=SAMPLES)
        targets += expand(NANOPLOT_FILTERED_DIR, sample=SAMPLES)
        targets += expand(ASSEMBLY_PROVENANCE, sample=SAMPLES)
        targets += expand(QUALIMAP_DIR, sample=SAMPLES)
        targets += expand(MEDAKA_CONSENSUS, sample=MEDAKA_SAMPLES)
        targets += expand(FASTP_JSON, sample=SAMPLES)
        targets += expand(SEL_R1, sample=SAMPLES)
        targets += expand(SEL_R2, sample=SAMPLES)
        targets += expand(POLYPOLISH_CONTIGS, sample=SAMPLES)
        targets += expand(GUIDE_QUALIMAP_DIR, sample=SAMPLES)
    return targets


def _shared_targets():
    # The deliverables every mode ends with, in the order the stages run: assembly QC
    # (with the contig summary, made last but delivered next to it), the ITS barcode and
    # its classification, the annotation, the AMR screen, and the MultiQC page.
    # Intermediate files are not listed; Snakemake builds them because these need them.
    # FUNANNOTATE_DB_MARKER and MULTIQC_REPORT are one file per run, not per sample.
    # The SignalP 6 result is only asked for when SignalP 6 is switched on, and the AMR
    # summary only when the AMR screen is.
    # Called by all_targets() below.
    targets = []
    targets += expand(QUAST_DIR, sample=SAMPLES)
    targets += expand(BUSCO_DIR, sample=SAMPLES)
    # The contig summary (rule contig_summary, 25_contig_summary.smk): what each contig of
    # the draft is, one table per sample. It needs the ITS windows, BUSCO, the finished
    # annotation and the telomere call of rule telomere_units, so it is one of the last
    # jobs of a sample. The same job writes the table's MultiQC copy; the unit file and
    # the ends table are its inputs, and the candidates table comes from the same job
    # as those two, so none of them needs a listing of its own.
    targets += expand(CONTIG_SUMMARY, sample=SAMPLES)
    # The contig adapter screen (rule screen_adapters, 25_contig_summary.smk), in every mode;
    # the contig summary reads it too. For the ONT samples whose reads fastplong trimmed, the
    # count of adapters left in the trimmed reads (rule check_ont_trimming, 05_ont_adapters.smk).
    targets += expand(ADAPTER_SCREEN, sample=SAMPLES)
    if HAS_LONG_READS and ONT_TRIM_SAMPLES:
        targets += expand(ONT_ADAPTER_CHECK, sample=ONT_TRIM_SAMPLES)

    targets += expand(ITS_FULL, sample=SAMPLES)
    targets += expand(ITS_ITS1, sample=SAMPLES)
    targets += expand(ITS_ITS2, sample=SAMPLES)
    targets += expand(ITS_SUMMARY, sample=SAMPLES)
    targets += expand(ITS_COPIES, sample=SAMPLES)
    targets += expand(ITS1_TAXONOMY, sample=SAMPLES)
    targets += expand(ITS2_TAXONOMY, sample=SAMPLES)

    targets += [FUNANNOTATE_DB_MARKER]
    targets += expand(CONTIGS_CLEAN, sample=SAMPLES)
    targets += expand(CONTIGS_SORT, sample=SAMPLES)
    targets += expand(CONTIGS_MASK, sample=SAMPLES)
    targets += expand(FUNANNOTATE_DIR, sample=SAMPLES)
    targets += expand(IPRSCAN_XML, sample=SAMPLES)
    targets += expand(EGGNOG_DIR, sample=SAMPLES)
    targets += expand(ANTISMASH_DIR, sample=SAMPLES)
    if SIGNALP6_ENABLED:
        targets += expand(SIGNALP6_RESULTS, sample=SAMPLES)
    targets += expand(CLEANUP_COMPLETE, sample=SAMPLES)

    # The AMR summary of every sample (rule amr_chroquetas, 50_amr.smk). The same job writes
    # the positions table, the MultiQC copy and ChroQueTas's own folder, so the summary
    # alone stands for all four.
    if AMR_ENABLED:
        targets += expand(AMR_SUMMARY, sample=SAMPLES)

    targets += [MULTIQC_REPORT]
    return targets


def all_targets():
    # Called once, when Snakemake reads `rule all` at the very end of the Snakefile, so
    # SAMPLES, the mode flags and every constant above are final by then. Returns the
    # front-end files followed by the shared ones, or an empty list (see EMPTY RULE
    # FOLDER above).
    if not _rule_modules_present():
        return []
    return _frontend_targets() + _shared_targets()


# ── The download targets: what rule download_all asks for ──
# download_targets() is what `rule download_all` asks for (workflow/Snakefile, after
# rule all). Naming that rule on the command line runs the download steps of this
# configuration and nothing else, so that a run can be prepared on a computer with
# internet before the real run goes to compute nodes that have none (the Snakefile says
# how to use it).
#
# The list holds, for each download rule that THIS configuration defines, the file
# Snakemake checks to know that the rule has run. Each `if` below is the same switch
# that decides whether the rule exists at all, so the list never asks for a file that
# no rule of this run can make (Snakemake would stop with an error). Rule by rule:
#   download_phix       PHIX_FASTA. Illumina and hybrid mode only (HAS_SHORT_READS):
#                       the rule sits in rules/illumina/ and rules/hybrid/.
#   funannotate_db      FUNANNOTATE_DB_MARKER, in every mode. This is the marker that
#                       rule normalize_funannotate_db writes, not the index file of the
#                       database itself: the marker exists only once the database has
#                       been downloaded and unpacked (if its folder held none yet) AND
#                       the paths stored inside it have been repaired for this computer.
#                       The repair needs no internet; asking for it here means that the
#                       shared database folder is finished before any job of the real
#                       run reads it.
#   unite_db            UNITE_FASTA, unless directories.unite_db names a copy the user
#                       already holds (UNITEDB is then not empty, and nothing is
#                       downloaded).
#   telomere_db         TELOMERE_DB, in every mode (the contig summary runs in all).
#   univec_illumina     UNIVEC_ILLUMINA, in every mode (the contig adapter screen runs in
#                       all).
#   tiberius_model      TIBERIUS_MODEL_WEIGHTS, when Tiberius is on (TIBERIUS_ENABLED).
#   chroquetas_db       CHROQUETAS_SPECIES_LIST, when the AMR screen is on
#                       (AMR_ENABLED). The rule writes it last, once the database is
#                       unpacked.
#   antismash_db        ANTISMASH_DB_DIR, unless directories.antismash_db names a copy
#                       the user already holds (ANTISMASHDB is then not empty).
#   busco_data          BUSCO_DATA_READY, only when parameters.busco.lineage names a
#                       lineage. With auto (the default) BUSCO picks the lineage only
#                       after it has looked at each genome, so there is nothing to fetch
#                       in advance: every BUSCO job downloads its own data.
#   check_medaka_model  MEDAKA_MODEL_RESOLVED, once for every ONT sample polished with
#                       Medaka (MEDAKA_SAMPLES; longread and hybrid mode only). The check
#                       stores the weights of the sample's model, so that polishing needs
#                       no internet. It reads the sample's raw reads, so these must be
#                       in place.
# Not in the list: what the user installs and names by path (the BLAST and eggNOG
# databases, InterProScan, GeneMark-ES, the SignalP 6 archive), and the conda
# environments.
#
# Returns a list of paths, with {sample} filled in for the Medaka files. Called once,
# when Snakemake reads rule download_all, so every name above is final by then.
def download_targets():
    targets = []

    # The PhiX genome (Illumina and hybrid mode). Its rule sits in the mode's own folder,
    # so it is asked for only when that folder holds rule files (EMPTY RULE FOLDER, at
    # the top of this section).
    if HAS_SHORT_READS and _rule_modules_present():
        targets.append(PHIX_FASTA)

    # The funannotate database, downloaded if needed and then repaired.
    targets.append(FUNANNOTATE_DB_MARKER)

    # UNITE, unless the user holds a copy.
    if not UNITEDB:
        targets.append(UNITE_FASTA)

    # The telomere repeat database.
    targets.append(TELOMERE_DB)

    # The Illumina adapter records of UniVec, for the contig adapter screen.
    targets.append(UNIVEC_ILLUMINA)

    # The Tiberius fungi model.
    if TIBERIUS_ENABLED:
        targets.append(TIBERIUS_MODEL_WEIGHTS)

    # ChroQueTas and its FungAMR database.
    if AMR_ENABLED:
        targets.append(CHROQUETAS_SPECIES_LIST)

    # The antiSMASH database, unless the user holds a copy.
    if not ANTISMASHDB:
        targets.append(ANTISMASH_DB_DIR)

    # The BUSCO lineage, when one is named.
    if BUSCO_LINEAGE != "auto":
        targets.append(BUSCO_DATA_READY)

    # The checked Medaka model of each ONT sample polished with Medaka: the same test as
    # the `if` around rule check_medaka_model (12_medaka_check.smk).
    if HAS_LONG_READS and MEDAKA_SAMPLES:
        targets += expand(MEDAKA_MODEL_RESOLVED, sample=MEDAKA_SAMPLES)

    return targets
