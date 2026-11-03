"""Every rule that downloads must run where Snakemake runs, with pinned download tools.

What this test protects
-----------------------
Two things go wrong with downloads when FunFlux leaves the computer it was written on.

1. Cluster compute nodes often have no internet. On a cluster, Snakemake sends every
   job to a compute node, except the rules marked as local rules (`localrule: True`
   in the rule, or a `localrules:` line naming them), which it runs itself, on the
   computer it was started on. So every rule that fetches something from the
   internet must be a local rule.
2. curl, tar and gzip are not on every Linux computer (a minimal Debian or Ubuntu
   system has no curl, a minimal Red Hat system no tar). A rule that runs them must
   get them from a pinned conda environment: its pin list (<name>.linux-64.pin.txt,
   the exact package list Snakemake builds the environment from) must hold each of
   the three programs the shell runs (tar also needs gzip, which it calls to unpack
   a .tar.gz archive).

So this test reads every rule of every .smk file under workflow/rules/, with the
rule reader of test_python_environments.py (the same folder), and requires:
  * for each rule whose shell runs curl, tar or gzip: a conda environment whose pin
    list holds that program (and gzip for tar), and that the rule is local when it
    runs curl;
  * for the rules that download by other means (download-antismash-databases,
    busco --download, and the Medaka model check, which makes Medaka download the
    weights of a named model): that they are local.
A rule that breaks this is reported with its name, its file and the reason.

It also checks rule download_all, which a user names on the command line to fetch
everything a run downloads before the run goes to compute nodes without internet. Its
input list is download_targets() in workflow/rules/shared/00_common.smk. The test
requires:
  * that the list names an output of every download rule found above (for
    funannotate_db: the marker of normalize_funannotate_db, the step that follows it),
    so that a new download rule cannot be left out of it (download_all would then miss
    it, and the real run would try to download on a compute node);
  * that download_all is a local rule whose input is download_targets().
Whether each name sits under the right switch (PhiX only with Illumina reads, and so
on) is checked by dry runs in miscellaneous/check_parse_time.py, which this test,
reading the files as text, cannot do.

Run with:  python -m pytest miscellaneous -q
"""

import os
import re

# The repository root, the rule folder, the rule reader and the conda-line reader of the
# Python-environment test, which reads the rule files the same way (as plain text,
# without Snakemake).
from test_python_environments import REPO_ROOT, RULES_DIR, conda_path, read_rules


# The rules that fetch from the internet without running curl themselves, by name:
#   antismash_db       download-antismash-databases, a command of the antiSMASH package
#   busco_data         busco --download, for a named BUSCO lineage
#   check_medaka_model medaka_model_check.py, which makes Medaka store the weights of a
#                      named model (and, in auto mode, of the model it works out)
OTHER_DOWNLOAD_RULES = {"antismash_db", "busco_data", "check_medaka_model"}

# The three programs, each used as a command word: not inside a longer name such as
# libcurl, curl_args or tar.gz, and not as part of a path (envs/download.yaml).
#   (?<![\w./-])  the character in front is not a letter, digit, _, ., / or -
#   (?![\w./-])   nor is the character after it
PROGRAM_PATTERNS = {
    name: re.compile(r"(?<![\w./-])" + name + r"(?![\w./-])")
    for name in ("curl", "tar", "gzip")
}

# The shell line of a rule, and a `localrules:` line of a file (the older way of naming
# local rules, which 40_annotation.smk uses for augustus_model_check).
SHELL_LINE = re.compile(r"^\s*shell:")
LOCALRULES_LINE = re.compile(r"^\s*localrules:\s*(.+)$")
LOCALRULE_FLAG = re.compile(r"^\s*localrule:\s*True\s*$")


def all_rule_files():
    """Every .smk file under workflow/rules/, sorted, as full paths."""
    found = []
    for folder, _subfolders, names in os.walk(RULES_DIR):
        for name in names:
            if name.endswith(".smk"):
                found.append(os.path.join(folder, name))
    return sorted(found)


def shell_lines(rule_lines):
    """The lines of one rule from its `shell:` line to its end, comment lines left out."""
    lines = []
    in_shell = False
    for line in rule_lines:
        if SHELL_LINE.match(line):
            in_shell = True
        if in_shell and not line.strip().startswith("#"):
            lines.append(line)
    return lines


def local_names_of_file(smk_path):
    """The rule names a `localrules:` line of this file names (an empty set if none)."""
    names = set()
    with open(smk_path, encoding="utf-8") as handle:
        for line in handle:
            match = LOCALRULES_LINE.match(line)
            if match:
                for name in match.group(1).split(","):
                    names.add(name.strip())
    return names


def is_local(rule_name, rule_lines, file_local_names):
    """True when the rule carries `localrule: True` or its file's localrules line names it."""
    if rule_name in file_local_names:
        return True
    for line in rule_lines:
        if LOCALRULE_FLAG.match(line):
            return True
    return False


def pin_list_text(smk_path, environment):
    """The text of the pin list next to the rule's environment file, or None if absent.

    The environment path is written relative to the rule file, for example
    ../../envs/download.yaml; its pin list is ../../envs/download.linux-64.pin.txt."""
    yaml_path = os.path.normpath(os.path.join(os.path.dirname(smk_path), environment))
    pin_path = yaml_path[: -len(".yaml")] + ".linux-64.pin.txt"
    if not os.path.isfile(pin_path):
        return None
    with open(pin_path, encoding="utf-8") as handle:
        return handle.read()


def pin_holds(pin_text, program):
    """True when the pin list holds the conda package of that name (for example
    .../linux-64/curl-8.22.0-ha042cf0_0.conda for curl)."""
    return re.search(r"/" + program + r"-\d[^/\s]*\.(conda|tar\.bz2)\s*$", pin_text, flags=re.M) is not None


def problems_in(smk_path):
    """Every rule of one file that breaks the rules above, as readable sentences."""
    problems = []
    relative = os.path.relpath(smk_path, os.path.dirname(RULES_DIR))
    file_local_names = local_names_of_file(smk_path)
    for rule_name, rule_lines in read_rules(smk_path):
        shell_text = "\n".join(shell_lines(rule_lines))
        programs = [name for name, pattern in PROGRAM_PATTERNS.items() if pattern.search(shell_text)]
        # tar unpacks a .tar.gz archive by calling gzip, so it needs gzip too.
        if "tar" in programs and "gzip" not in programs:
            programs.append("gzip")
        local = is_local(rule_name, rule_lines, file_local_names)

        if programs:
            environment = conda_path(rule_lines)
            if environment is None:
                problems.append(
                    f"rule {rule_name} ({relative}): its shell runs {', '.join(programs)} but it has no "
                    "conda environment, so these come from the computer, which may lack them"
                )
            else:
                pin_text = pin_list_text(smk_path, environment)
                if pin_text is None:
                    problems.append(f"rule {rule_name} ({relative}): {environment} has no pin list next to it")
                else:
                    for program in programs:
                        if not pin_holds(pin_text, program):
                            problems.append(
                                f"rule {rule_name} ({relative}): its shell runs {program}, but the pin list "
                                f"of {environment} holds no {program} package"
                            )
        if ("curl" in programs or rule_name in OTHER_DOWNLOAD_RULES) and not local:
            problems.append(
                f"rule {rule_name} ({relative}): it downloads from the internet but is not a local rule "
                "(localrule: True), so on a cluster it would run on a compute node, which may have no internet"
            )
    return problems


def test_download_rules_are_local_and_pinned():
    problems = []
    for smk_path in all_rule_files():
        problems.extend(problems_in(smk_path))
    assert problems == [], "\n".join(problems)


def test_the_known_download_rules_are_found():
    # A guard on the guard: if the reader stopped finding the rules (a renamed rule, a
    # changed layout), the test above would pass on nothing. These rules download in
    # every FunFlux release and must be seen as downloading and local.
    expected = {"download_phix", "funannotate_db", "unite_db", "telomere_db", "tiberius_model",
                "chroquetas_db", "antismash_db", "busco_data", "check_medaka_model"}
    seen_local_downloaders = set()
    for smk_path in all_rule_files():
        file_local_names = local_names_of_file(smk_path)
        for rule_name, rule_lines in read_rules(smk_path):
            shell_text = "\n".join(shell_lines(rule_lines))
            downloads = PROGRAM_PATTERNS["curl"].search(shell_text) or rule_name in OTHER_DOWNLOAD_RULES
            if downloads and is_local(rule_name, rule_lines, file_local_names):
                seen_local_downloaders.add(rule_name)
    assert expected <= seen_local_downloaders, sorted(expected - seen_local_downloaders)


# ── rule download_all ──

# The file that defines download_targets(), and the entry point (workflow/Snakefile)
# that defines rule download_all.
COMMON_FILE = os.path.join(RULES_DIR, "shared", "00_common.smk")
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")

# A name of a constant as the rule files write it: capital letters, digits and _, with
# an optional leading _ (PHIX_FASTA, BUSCO_DATA_READY, _FUNANNOTATE_DB_INFO).
CONSTANT_NAME = re.compile(r"\b_?[A-Z][A-Z0-9_]*\b")

# The output line of a rule; an output may follow on the same line.
OUTPUT_LINE = re.compile(r"^(\s*)output:\s*(.*)$")

# The one download rule whose output download_all does not name directly. download_all
# asks for the funannotate database through FUNANNOTATE_DB_MARKER, the file that rule
# normalize_funannotate_db writes once it has repaired the paths inside the database, and
# that rule takes the output of funannotate_db as its input. (That output,
# _FUNANNOTATE_DB_INFO, is a name of 40_annotation.smk alone, so 00_common.smk could not
# name it anyway.)
REACHED_THROUGH = {"funannotate_db": "FUNANNOTATE_DB_MARKER"}


def download_targets_names():
    """The constant names used in the body of download_targets() in 00_common.smk.

    The body runs from the `def download_targets():` line to the next line of code at the
    left margin (or the end of the file); comment lines are left out."""
    with open(COMMON_FILE, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    names = set()
    inside = False
    for line in lines:
        stripped = line.strip()
        if line.startswith("def download_targets():"):
            inside = True
            continue
        if not inside or stripped == "" or stripped.startswith("#"):
            continue
        # A line of code at the left margin: the function has ended.
        if not line[0].isspace():
            break
        names.update(CONSTANT_NAME.findall(line))
    return names


def output_names(rule_lines):
    """The constant names on the output lines of one rule (from read_rules).

    The outputs run from the `output:` line to the next line indented no deeper than it
    (the next directive, such as params: or conda:)."""
    names = set()
    output_indent = None
    for line in rule_lines:
        stripped = line.strip()
        # Before the output line: look for it, and take any output written on it.
        if output_indent is None:
            match = OUTPUT_LINE.match(line)
            if match:
                output_indent = len(match.group(1))
                names.update(CONSTANT_NAME.findall(match.group(2)))
            continue
        if stripped == "" or stripped.startswith("#"):
            continue
        # The next directive (params:, conda:, ...) ends the outputs.
        indent = len(line) - len(line.lstrip())
        if indent <= output_indent:
            break
        names.update(CONSTANT_NAME.findall(line))
    return names


def local_download_rules():
    """{rule name: its output names} for every rule that downloads and is local (the same
    test as test_the_known_download_rules_are_found)."""
    found = {}
    for smk_path in all_rule_files():
        file_local_names = local_names_of_file(smk_path)
        for rule_name, rule_lines in read_rules(smk_path):
            shell_text = "\n".join(shell_lines(rule_lines))
            downloads = PROGRAM_PATTERNS["curl"].search(shell_text) or rule_name in OTHER_DOWNLOAD_RULES
            if downloads and is_local(rule_name, rule_lines, file_local_names):
                found[rule_name] = output_names(rule_lines)
    return found


def test_download_all_asks_for_every_download_rule():
    asked_for = download_targets_names()
    assert asked_for, "download_targets() was not found in 00_common.smk, or names nothing"
    problems = []
    # One output of a rule is enough: asking for it makes Snakemake run the rule.
    for rule_name, outputs in sorted(local_download_rules().items()):
        if rule_name in REACHED_THROUGH:
            wanted = {REACHED_THROUGH[rule_name]}
        else:
            wanted = outputs
        # wanted & asked_for: the names that are in both sets (none: the rule is missed).
        if not wanted & asked_for:
            problems.append(
                f"rule {rule_name}: download_targets() in 00_common.smk names none of {sorted(wanted)}, "
                "so snakemake download_all would not run it"
            )
    assert problems == [], "\n".join(problems)


def test_download_all_is_a_local_rule_asking_for_download_targets():
    rules = dict(read_rules(SNAKEFILE))
    assert "download_all" in rules, f"{SNAKEFILE} defines no rule download_all"
    rule_lines = rules["download_all"]
    assert any(LOCALRULE_FLAG.match(line) for line in rule_lines), "rule download_all is not a local rule"
    assert any(line.strip() == "download_targets()" for line in rule_lines), (
        "rule download_all does not take download_targets() as its input"
    )
