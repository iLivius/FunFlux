"""Every rule that starts Python must get it from a pinned conda environment.

What this test protects
-----------------------
Many FunFlux rules run a small Python script with a shell command such as
`python mitochondria.py screen ...`. The shell takes `python` from the PATH of
the job (the list of folders it searches for a program). When the rule has a
conda environment, Snakemake activates it first, and `python` is that
environment's own, pinned Python. When the rule has none, `python` is whatever
the computer offers, and that can be missing (Ubuntu 22.04 and later, Debian 12
and RHEL 9 install only `python3`), an old Python 2.7 or 3.6 on an old cluster,
or a different Python on each cluster node. A conda environment without Python
in it is no better: activating it takes the Snakemake environment off the PATH,
so `python` again comes from the computer.

So this test reads every rule of every .smk file under workflow/rules/ and, for
each rule whose shell command starts Python (as `python`, `python3`, a versioned
name such as `python3.12`, or a full path such as `/usr/bin/python3`), requires:
  * a `conda:` line, and
  * a pin list next to that environment file (<name>.linux-64.pin.txt, the exact
    package list Snakemake builds the environment from) that holds a Python 3
    package.
A rule that breaks this is reported with its name, its file and the reason.

How it reads the rule files
---------------------------
As plain text, without Snakemake, so it needs no configuration and no input
files. A rule starts at a line `rule <name>:` and ends at the next line of code
that is indented no deeper than that line. Lines inside a triple-quoted string
(a shell command, which may hold a here-document whose lines start at the left
margin) never end a rule. Only the lines from `shell:` to the end of the rule
are searched for `python`, and lines that start with # (comments) are skipped.
A rule that starts Python in another way (a `script:` line, a `run:` block) is
not looked at.

Run with:  python -m pytest miscellaneous -q
"""

import os
import re


# The repository root (the folder above miscellaneous/) and the rule files under it.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RULES_DIR = os.path.join(REPO_ROOT, "workflow", "rules")

# The first line of a rule, for example "    rule select_contigs:". The first group is
# the indentation in front of it, the second the rule's name.
RULE_START = re.compile(r"^(\s*)rule\s+(\w+)\s*:")

# The conda line of a rule, with or without the path on the same line.
CONDA_LINE = re.compile(r"^\s*conda:\s*(.*)$")

# The shell line of a rule; a short command may follow on the same line.
SHELL_LINE = re.compile(r"^\s*shell:\s*(.*)$")

# Python used as a program: `python`, `python3` or a versioned name such as
# `python3.12`, on its own or at the end of a full path such as /usr/bin/python3
# (both start the computer's own Python, which is what this test must catch).
# Not counted: the word inside a longer name or a folder name, so {params.python},
# envs/python.yaml, python_version, my-python and lib/python3.12/ do not count.
# How the pattern says this:
#   (?<![\w.-])       the character in front is not a letter, digit, _, . or -
#                     (a space, "/", "$(", "!" or ";" are fine, as is the line start)
#   python(?:3(?:\.\d+)?)?   python, python3 or python3.<digits>
#   (?![\w./-])       the character after is not a letter, digit, _, ., / or -
PYTHON_COMMAND = re.compile(r"(?<![\w.-])python(?:3(?:\.\d+)?)?(?![\w./-])")

# A Python 3 package in a pin list, for example
# https://conda.anaconda.org/conda-forge/linux-64/python-3.12.15-h5f976f7_1_cpython.conda
PYTHON3_PACKAGE = re.compile(r"/python-3\.\d+[^/]*\.(conda|tar\.bz2)$")


def read_rules(smk_path):
    """Cut one rule file into its rules.

    Input:  the path of a .smk file.
    Output: a list with one (rule name, list of the rule's lines) pair per rule,
            in file order. The first line of each list is the `rule <name>:` line.
    """
    with open(smk_path, encoding="utf-8") as handle:
        lines = handle.read().splitlines()

    rules = []
    current_name = None      # the rule being read, or None between rules
    current_indent = 0       # the indentation of its `rule <name>:` line
    current_lines = []
    inside_string = False    # True while inside a triple-quoted string

    for line in lines:
        stripped = line.strip()
        is_comment_or_blank = stripped == "" or stripped.startswith("#")

        if not inside_string and not is_comment_or_blank:
            rule_match = RULE_START.match(line)
            indent = len(line) - len(line.lstrip())
            # A line of code no deeper than the rule line ends the rule being read.
            if current_name is not None and indent <= current_indent:
                rules.append((current_name, current_lines))
                current_name = None
                current_lines = []
            # A new rule starts here.
            if rule_match:
                current_name = rule_match.group(2)
                current_indent = len(rule_match.group(1))

        if current_name is not None:
            current_lines.append(line)

        # An odd number of triple quotes on a line opens or closes a string. Comment
        # lines outside a string are skipped: some comments mention the quotes.
        if inside_string or not is_comment_or_blank:
            if line.count('"""') % 2 == 1:
                inside_string = not inside_string

    if current_name is not None:
        rules.append((current_name, current_lines))
    return rules


def conda_path(rule_lines):
    """Return the environment path written in the rule's conda line, or None.

    Input:  the lines of one rule (from read_rules). The path may stand on the
            conda line itself or on the next line that is not blank or a comment.
    Output: the path as written, without its quotes (for example
            "../../envs/python.yaml"), or None when the rule has no conda line.
    """
    for position, line in enumerate(rule_lines):
        match = CONDA_LINE.match(line)
        if not match:
            continue
        value = match.group(1).strip()
        if value == "":
            for next_line in rule_lines[position + 1:]:
                next_stripped = next_line.strip()
                if next_stripped != "" and not next_stripped.startswith("#"):
                    value = next_stripped
                    break
        return value.strip().strip(",").strip('"').strip("'")
    return None


def shell_runs_python(rule_lines):
    """Return True when the rule's shell command starts Python (see PYTHON_COMMAND).

    Input:  the lines of one rule (from read_rules).
    Looks at the lines from `shell:` to the end of the rule (the rest of the shell
    line included) and skips comment lines.
    """
    in_shell = False
    for line in rule_lines:
        if not in_shell:
            match = SHELL_LINE.match(line)
            if match:
                in_shell = True
                line = match.group(1)
            else:
                continue
        if line.strip().startswith("#"):
            continue
        if PYTHON_COMMAND.search(line):
            return True
    return False


def pin_list_has_python3(pin_path):
    """Return True when the pin list holds a Python 3 package (see PYTHON3_PACKAGE)."""
    with open(pin_path, encoding="utf-8") as handle:
        for line in handle:
            if PYTHON3_PACKAGE.search(line.strip()):
                return True
    return False


def python_rules(rules_dir):
    """List every rule under rules_dir whose shell command starts Python.

    Input:  a folder of rule files (normally workflow/rules/); every .smk file in
            it and its subfolders is read.
    Output: a list of (path of the .smk file, rule name, conda path or None).
    """
    found = []
    for folder, _subfolders, file_names in sorted(os.walk(rules_dir)):
        for file_name in sorted(file_names):
            if not file_name.endswith(".smk"):
                continue
            smk_path = os.path.join(folder, file_name)
            for rule_name, rule_lines in read_rules(smk_path):
                if shell_runs_python(rule_lines):
                    found.append((smk_path, rule_name, conda_path(rule_lines)))
    return found


def find_problems(rules_dir):
    """Check every rule that starts Python, and describe each one that breaks the rule.

    Input:  a folder of rule files (see python_rules).
    Output: one line of text per problem: the rule, its file and the reason. An
            empty list means every such rule gets a pinned Python 3.
    """
    problems = []
    for smk_path, rule_name, environment in python_rules(rules_dir):
        # The file is named from the rules folder on, for example
        # rules/shared/15_mitochondria.smk.
        short_path = os.path.relpath(smk_path, os.path.dirname(rules_dir))
        where = f"rule {rule_name} ({short_path})"
        if environment is None:
            problems.append(f"{where}: its shell runs python but it has no conda environment, "
                            f"so python would be whatever the computer has")
            continue
        if not environment.endswith(".yaml"):
            problems.append(f"{where}: the conda line '{environment}' is not a path to a .yaml file, "
                            f"so this test cannot find its pin list")
            continue
        # Snakemake reads a conda path relative to the rule file, and finds the pin list
        # next to the .yaml file under the name <name>.linux-64.pin.txt.
        yaml_path = os.path.normpath(os.path.join(os.path.dirname(smk_path), environment))
        pin_path = yaml_path[:-len(".yaml")] + ".linux-64.pin.txt"
        if not os.path.isfile(yaml_path):
            problems.append(f"{where}: its environment file {yaml_path} does not exist")
        elif not os.path.isfile(pin_path):
            problems.append(f"{where}: its environment has no pin list {pin_path}")
        elif not pin_list_has_python3(pin_path):
            problems.append(f"{where}: its shell runs python, but the pin list {pin_path} holds no "
                            f"Python 3, so python would be whatever the computer has")
    return problems


# =============================================================================
# The real rule files
# =============================================================================

def test_every_rule_that_runs_python_has_a_pinned_python():
    # The check itself: no rule of the workflow may start a Python it does not pin.
    problems = find_problems(RULES_DIR)
    assert problems == [], "\n".join(problems)


def test_the_reader_finds_the_rules_that_run_python():
    # A reader that found no rule at all would pass the test above without checking
    # anything. These rules are known to run a FunFlux script with `python`, so each
    # must be found, with the environment it is meant to use.
    found = {}
    for _smk_path, rule_name, environment in python_rules(RULES_DIR):
        found[rule_name] = environment
    expected = {
        "select_contigs": "../../envs/python.yaml",
        "screen_mitochondria": "../../envs/python.yaml",
        "separate_mitochondria": "../../envs/python.yaml",
        "normalize_funannotate_db": "../../envs/python.yaml",
        "augustus_model_check": "../../envs/python.yaml",
        "assembly_evidence": "../../envs/python.yaml",
        "screen_guide_mitochondria": "../../envs/python.yaml",
        "select_guide_contigs": "../../envs/python.yaml",
        "contig_summary": "../../envs/vsearch.yaml",
        "check_medaka_model": "../../envs/medaka.yaml",
    }
    for rule_name, environment in expected.items():
        assert rule_name in found, f"rule {rule_name} was not recognised as a rule that runs python"
        assert found[rule_name] == environment, (
            f"rule {rule_name}: expected the environment {environment}, found {found[rule_name]}")
    # The assembler runs no Python of its own: its evidence table is made by rule
    # assembly_evidence, so a failed table can never cost the assembly.
    assert "long_read_assembly" not in found


def test_the_pattern_tells_python_commands_from_other_words():
    # Shell lines that start a Python: each must be recognised, or a rule written
    # this way would slip past the check above.
    starts_python = [
        "python script.py",
        "python3 script.py",
        "python3.12 script.py",
        "/usr/bin/python3 script.py",
        "x=$(python -c 'print(1)')",
        "python - <<'PY'",
    ]
    for line in starts_python:
        assert PYTHON_COMMAND.search(line), f"not recognised as starting Python: {line}"
    # Lines that only mention the word: none may count, or the check would complain
    # about rules that start no Python at all.
    mentions_only = [
        "{params.python} script.py",
        "cp envs/python.yaml backup/",
        "echo python_version",
        "my-python script.py",
        "ls $CONDA_PREFIX/lib/python3.12/site-packages",
    ]
    for line in mentions_only:
        assert not PYTHON_COMMAND.search(line), f"wrongly counted as starting Python: {line}"


# =============================================================================
# A small made-up rule folder: the test must catch each kind of mistake
# =============================================================================

def write_file(path, text):
    """Write text into path, creating the folders above it as needed."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


# Four made-up rules: one that is right, one with no environment, one whose
# environment has no Python, and one that mentions python only in a comment. The
# here-document in the first rule starts at the left margin, as the real one in
# 40_annotation.smk does; it must not end the rule early.
FAKE_RULES = '''
rule good:
    output: "a.txt"
    conda:
        "../../envs/withpython.yaml"
    shell:
        """
        python - <<'PY'
import sys
PY
        python script.py > {output}
        """

if True:
    rule no_environment:
        output: "b.txt"
        shell:
            """
            python script.py > {output}
            """

rule environment_without_python:
    output: "c.txt"
    conda:
        "../../envs/nopython.yaml"
    shell:
        """
        x=$(python3 -c 'print(1)')
        """

rule python_only_in_a_comment:
    output: "d.txt"
    shell:
        """
        # python is not used here
        touch {output}
        """
'''


def make_fake_rules(tmp_path):
    """Build workflow/rules/shared/fake.smk and two environments under tmp_path.

    Output: the made-up rules folder, ready for find_problems().
    """
    workflow = tmp_path / "workflow"
    write_file(str(workflow / "rules" / "shared" / "fake.smk"), FAKE_RULES)
    write_file(str(workflow / "envs" / "withpython.yaml"), "dependencies:\n  - python=3.12\n")
    write_file(str(workflow / "envs" / "withpython.linux-64.pin.txt"),
               "@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/python-3.12.15-h5f976f7_1_cpython.conda\n")
    write_file(str(workflow / "envs" / "nopython.yaml"), "dependencies:\n  - hifiasm=0.25.0\n")
    write_file(str(workflow / "envs" / "nopython.linux-64.pin.txt"),
               "@EXPLICIT\nhttps://conda.anaconda.org/bioconda/linux-64/hifiasm-0.25.0-h5ca1c30_0.tar.bz2\n")
    return str(workflow / "rules")


def test_made_up_rules_are_read_correctly(tmp_path):
    # Three of the four rules start Python; the comment does not count. Each keeps its
    # own conda line (or none).
    rules_dir = make_fake_rules(tmp_path)
    found = {}
    for _smk_path, rule_name, environment in python_rules(rules_dir):
        found[rule_name] = environment
    assert found == {
        "good": "../../envs/withpython.yaml",
        "no_environment": None,
        "environment_without_python": "../../envs/nopython.yaml",
    }


def test_made_up_mistakes_are_reported(tmp_path):
    # The rule without an environment and the rule whose environment has no Python
    # are reported, each once and by name; the correct rule is not.
    rules_dir = make_fake_rules(tmp_path)
    problems = find_problems(rules_dir)
    assert len(problems) == 2
    assert problems[0].startswith("rule no_environment (")
    assert "has no conda environment" in problems[0]
    assert problems[1].startswith("rule environment_without_python (")
    assert "holds no Python 3" in problems[1]
