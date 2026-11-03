#!/usr/bin/env python3
"""Check that the things FunFlux's comments name actually exist.

    (Adopted from BacFlux, where it was written. Deliberate differences from BacFlux's
    copy:
      * --config chooses the config file that dotted keys are checked against. The
        default, config/config.yaml, is the one file BacFlux's copy reads.
      * Added on 2026-10-08, to close these blind spots:
          - in a git checkout only tracked files are read, and a path that git
            ignores counts as present;
          - the comments of the two config templates are read;
          - in Markdown the history skip ends at a blank line, and "1.x" counts as a
            word about history;
          - a line that says "planned", "moved from" or "moved into" is skipped
            (that line only);
          - the comment lines inside shell: commands are counted correctly;
          - paths written from inside workflow/ (rules/..., scripts/..., envs/...)
            are checked;
          - path templates and web addresses are left alone;
          - the summary says what each count holds.)

WHY THIS EXISTS
    Comments go out of date without anyone noticing. A rule gets renamed, a stage
    directory is renumbered, a constant moves to another file, and every comment
    pointing at the old name keeps sitting there looking authoritative. Nothing
    fails, no test goes red, and the reader who trusts it loses an afternoon.

    A large comment rewrite of BacFlux (2026-08-14), the sibling workflow this script
    comes from, turned up dozens of exactly that: a cited `rule mobilome_name_transposons`
    that was never built, v1 stage paths (`05.annotation` where v2 has `04.annotation`),
    a count that said four where the code had five. Every one of them could have been
    found by a script. This script is that check, so the next batch gets caught in
    seconds instead of by someone reading the file two years later.

WHAT IT READS
    In a git checkout (the normal case) only the files git tracks are read. A file is
    tracked from the moment it is added with git add, whether or not it has been
    committed since. A personal config, a scratch folder or any other file never added
    with git add does not change the result. When the repository folder is not the top
    folder of a git checkout (a copied or unpacked folder), or when git is not
    installed, every file is found by walking the folders, except in folders named
    .git, .snakemake, __pycache__ or ref.

    Of those files it reads
      * the comment lines (lines whose first character other than a space is #) of
        every .smk, .py and .sh file, of the file named Snakefile, and of the two
        config templates, config/config.yaml and config/config_test_examples.yaml;
      * every line of every .md file: a Markdown page is prose from top to bottom.
    Not read: the conda environment files in workflow/envs/. They hold version pins
    with short notes, and their text, comments included, is part of the fingerprint
    of each environment (the hash Snakemake names the environment folder after), so
    rewording a note would rebuild that environment. They are left alone, and a rule,
    key or path named in their notes is not checked. Not read either: the comments of
    this script itself.

    A line that speaks about history ("v1", "FunFlux 1.x", "renamed", "no longer", ...)
    is skipped, together with the lines right after it: a former name is correct there
    precisely because it does not exist any more. In code the skip ends at the first
    line that is not a comment; in Markdown it ends at the first blank line, the end
    of the paragraph. A line that names a file not written yet or the old place of
    something that has moved ("planned", "moved from", "moved into") is skipped too,
    but only that line.

    A web address (http://, https://, ftp://) is removed from a line before anything
    in it is checked: in https://geneontology.org/docs/go-citation-policy/ the part
    docs/go-citation-policy is a page of that web site, not a folder of this
    repository.

WHAT IT CHECKS
    Failures, each one a yes-or-no question (the thing is there or it is not, with no
    judgement involved). Any failure makes the exit status 1.

      rules    every "rule <name>" or "rules <name>" with a snake_case name (one that
               contains an underscore) must be a rule defined in a .smk file or in a
               file named Snakefile
      paths    every repository path must exist: one that starts with workflow/,
               config/, docs/ or miscellaneous/, and one written from inside
               workflow/ (rules/..., scripts/..., envs/...), which is looked up under
               workflow/. In a git checkout a path exists when git tracks it (a file,
               or a folder that holds tracked files) or when git ignores it: the pages
               name config/config_custom.yaml, the user's own copy of the template,
               which .gitignore keeps out of the repository. A file that is on this
               computer but neither tracked nor ignored is reported, because a fresh
               clone would not have it. A tracked file deleted with plain rm still
               counts as present until the deletion is staged (git rm, or git add
               after the rm), because until then git still lists it. Outside a git
               checkout a path must simply exist on disk (so there
               config/config_custom.yaml is reported unless the file is present).
               A path template, where a { or a * follows the
               name (config/config_{mode}.yaml, workflow/scripts/its_*.py), is skipped;
               in workflow/rules/{mode}/ the folder workflow/rules is still checked.
      config   every dotted config key (parameters.busco.lineage) must exist in the
               config file named by --config (default config/config.yaml)

    Advisories, counted in the summary and listed with --advisories, never a failure:
      * a snake_case name in backticks that is neither a rule nor defined anywhere
        else, but is close to a rule name: a probable typo;
      * an ALL_CAPS name with an underscore that no file read assigns at the start
        of a line (NAME = ... in Python or Snakemake, NAME=... in a shell script).
        An assignment further along a line, such as the shell's export NAME=..., is
        not recognised, so the list also holds names that rule files set that way.

    For information it also counts the comment lines INSIDE the command of a shell:
    directive, between its opening and its closing triple quotes. They are part of the
    command text Snakemake compares, so rewording one makes Snakemake re-run the rule:
    harmless for a database download, expensive for gene prediction or InterProScan.
    Comment lines written between shell: and the opening quotes are not part of the
    command: they are not counted, and rewording them re-runs nothing.

WHAT IT DOES NOT CHECK
    Whether an explanation is TRUE. A comment saying "we filter here because X" is
    beyond any script; only a reader who knows the biology can judge it. This tool
    clears away the mechanical errors so that human attention goes where it is the
    only thing that works.

RUN
    python miscellaneous/check_comment_references.py           # from the repo root
    python miscellaneous/check_comment_references.py --quiet   # only the failures

    It is a development tool: nothing in the workflow calls it. Run it by hand after
    editing comments, renaming a rule or moving a file (workflow/scripts/README.md asks
    for it after moving a script), and before every release.

    It prints one "path:line: message" line per failure, then a short summary of the
    counts unless --quiet is given. With --advisories it also lists the advisories, in
    the same form.

    Exit code 0 = every reference resolves. 1 = at least one does not. 2 = the file
    named by --config does not exist.
"""

import argparse
import difflib
import os
import re
import subprocess
import sys


# The repository root: the folder above miscellaneous/. Every path is taken relative to it.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Files whose comments we read: names that end in one of SOURCE_SUFFIXES or that are one
# of SOURCE_NAMES, and the YAML files listed in CONFIG_TEMPLATES.
SOURCE_SUFFIXES = (".smk", ".py", ".sh", ".md")
# workflow/Snakefile is the entry point, the file Snakemake reads first.
SOURCE_NAMES = ("Snakefile",)
# The only YAML files read: the two config templates, whose comments explain every key
# to the user. No other YAML file is read. The conda environment files in workflow/envs/
# hold version pins with short notes; their text, comments included, is part of the
# fingerprint of each environment, so rewording a note would rebuild the environment.
# They are left alone, and a name in their notes goes unchecked (their pin files end in
# .txt and are not read either).
CONFIG_TEMPLATES = ("config/config.yaml", "config/config_test_examples.yaml")
# Folders that are never entered when the files are found by walking the folders (outside
# a git checkout, see source_files): git's own data (.git), Snakemake's own logs, job
# records and conda environments (.snakemake), Python's cache of compiled files
# (__pycache__), and ref (a name kept from BacFlux; FunFlux has no such folder). In a git
# checkout none of them holds a tracked file, so there the list is not needed.
SKIP_DIRS = {".git", ".snakemake", "__pycache__", "ref"}

# ── What a claim looks like ──────────────────────────────────────────────────
# Each pattern below is a regular expression (a text pattern) that finds one kind of name
# in a comment line.
#
# A rule: "rule foo_bar", "rules foo_bar" (only the name right after the word is
# captured). The UNDERSCORE is required, and that is a deliberate trade. Without it the
# pattern matches ordinary English ("the rule always fires", "the rule asks for") and
# produces so much noise that a real finding is invisible, which is how checkers end up
# ignored. The cost is that a reference to a single-word rule (multiqc, cleanup) is not
# checked. Most FunFlux rule names are snake_case, and every stale reference the
# 2026-08-14 BacFlux comment pass turned up was of that shape, so this catches the class
# that actually occurs.
RULE_REFERENCE = re.compile(r"\brules?\s+([a-z][a-z0-9]*(?:_[a-z0-9]+)+)\b")

# A repository path. Anchored on the four real top-level directories so that
# "AMR/mapping" or "identity/coverage" in prose is not mistaken for a file. A path must
# end in a letter, digit or underscore, so a full stop after it is not taken as part of
# it. It stops at a { or a *: in "config/config_{mode}.yaml" the match is
# "config/config_", which check() then recognises as the start of a template (see
# is_path_template).
PATH_REFERENCE = re.compile(
    r"\b((?:workflow|config|docs|miscellaneous)/[A-Za-z0-9_./-]*[A-Za-z0-9_])"
)

# A path written from inside workflow/, the way the rule files and
# workflow/scripts/README.md write them: rules/shared/10_decontam.smk,
# scripts/30_taxonomy/its_distinct.py, envs/itsx.yaml. check() looks it up under
# workflow/. "(?<![...])" means "not right after one of these characters": the
# "rules/..." inside "workflow/rules/..." belongs to a path PATH_REFERENCE already checks,
# and the "envs/..." in "../../envs/x.yaml" is relative to the folder of a rule file, so
# neither is read a second time here.
WORKFLOW_PATH_REFERENCE = re.compile(
    r"(?<![A-Za-z0-9_./-])((?:rules|scripts|envs)/[A-Za-z0-9_./-]*[A-Za-z0-9_])"
)

# A web address: http://, https:// or ftp:// and everything up to the next space. It is
# removed from a line before anything in it is checked (see the module docstring).
WEB_ADDRESS = re.compile(r"(?:https?|ftp)://\S+")

# An ALL_CAPS name. It must contain an underscore (so it is at least three characters
# long), so that shouted prose (NOT, ONLY, THIS) is not treated as a name.
GLOBAL_REFERENCE = re.compile(r"\b([A-Z][A-Z0-9]*(?:_[A-Z0-9]+)+)\b")

# How close a name must be to a real rule before it is called a probable typo: the
# similarity score of Python's difflib, from 0 (nothing in common) to 1 (identical).
# 0.70 is empirical: 'check_medaka_check' scores 0.78 against 'check_medaka_model', the
# actual mistake this check exists for, so a cutoff of 0.78 or more misses it; lowering
# it much below 0.70 starts pairing unrelated names that merely share a prefix.
NEAR_MISS_CUTOFF = 0.70

# A snake_case name written in backticks (`blob_json`) or as the first line of a node
# label in a Mermaid diagram ([blob_json<br/>... or [blob_json]): the two places
# BacFlux's documentation, where this script was written, names a rule. A label in
# quotes (["blob_json<br/>..."]) does not match. These are the names the near-miss
# advisory looks at.
SNAKE_NAME = re.compile(r"`([a-z][a-z0-9]*(?:_[a-z0-9]+)+)`"
                        r"|\[([a-z][a-z0-9]*(?:_[a-z0-9]+)+)(?:<br/>|\])")

# A dotted config key, all lowercase: parameters.busco.lineage, files.annotation_params.
CONFIG_REFERENCE = re.compile(r"\b([a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)+)\b")

# A comment that is deliberately describing what something USED to be called is
# correct precisely because the name no longer resolves. The workflows do this
# constantly ("(v1 message: ...)", "v1 rule was called map_qc", "these used to be
# flat"), and without this the checker flags the repo's own change log as broken.
# "1.x" covers "FunFlux 1.x", the pages' name for the v1 workflows. A line carrying one
# of these markers is skipped, and so are the lines right after it (see check()).
HISTORICAL_MARKERS = re.compile(
    r"\b(v1|1\.x|used to|was called|were called|renamed|formerly|old name|previously|"
    r"superseded|no longer|used to be|predates|pre-dates|stale|"
    r"that path|never built|does not exist|was never)\b", re.IGNORECASE
)

# Words that mark only their OWN line as naming something that is not in the repository
# as it is now: a file that is planned but not written yet ("rules/shared/12_medaka_check.smk
# (planned: ...)"), or the old place of something that has moved ("input: # MOVED from
# directories.input_dir", "parameters.filtlong: this setting moved into the block named
# after your mode"). Unlike HISTORICAL_MARKERS they do not carry over to the lines after
# them: they usually sit in one row of a table or one line of a YAML example, and the next
# row is about something else (in workflow/scripts/README.md the rows after the planned
# 12_medaka_check.smk name rule files that exist and must be checked). On 2026-10-08 every
# "moved from" and "moved into" in the files this tool reads described the change from v1
# to v2. Plain "moved" is not a marker, because it is also used for the present ("the
# output folder can be moved"), and "removed from" does not match (\b, a word boundary,
# must come right before "moved").
ONE_LINE_MARKERS = re.compile(r"\b(planned|moved from|moved into)\b", re.IGNORECASE)

# Snakemake's own names: the parts of a rule (input, output, params, wildcards, log,
# resources, threads), config (the loaded configuration) and rules (as in
# rules.filter_contigs.output, the output of another rule). In a comment, "input.contigs"
# nearly always means a rule's own input: entry, not config["input"]["contigs"], although
# the config really does have an input: section (and a resources: section). So a dotted
# name that begins with one of these words is never checked against the config. The
# cost: a typo in an input or resources key named in a comment is not caught.
SNAKEMAKE_DIRECTIVES = {"input", "output", "params", "wildcards", "log",
                        "resources", "threads", "config", "rules"}

# Endings that mark a dotted name as a file name, not a config key
# (contigs_final.fasta, its_windows.py).
CONFIG_SKIP_SUFFIXES = (
    ".py", ".sh", ".smk", ".tsv", ".txt", ".fasta", ".fa", ".fna", ".faa", ".gff",
    ".gff3", ".gbff", ".json", ".yaml", ".yml", ".log", ".html", ".zip", ".gz",
    ".bz2", ".md", ".fastq", ".fq", ".bam", ".sam", ".ms", ".gbk", ".fas", ".xz",
)


def git_tracked_files():
    """The files git tracks in the repository, as a list of paths relative to REPO_ROOT.

    Returns None when REPO_ROOT is not the top folder of a git checkout: git is not
    installed, the folder is not under git at all (a copied or unpacked repository),
    or it sits inside some OTHER git checkout, whose file list would not describe this
    repository. Called by source_files() and check().
    """
    # Ask git which folder is the top of the checkout REPO_ROOT is in; it fails when
    # there is none.
    try:
        top_level = subprocess.run(["git", "-C", REPO_ROOT, "rev-parse", "--show-toplevel"],
                                   capture_output=True, text=True)
    except FileNotFoundError:
        return None                                     # git is not installed
    if top_level.returncode != 0:
        return None                                     # not a git checkout
    # realpath resolves symbolic links, so one folder reached by two routes compares equal.
    if os.path.realpath(top_level.stdout.strip()) != os.path.realpath(REPO_ROOT):
        return None                                     # a folder inside another checkout

    # The tracked files. -z separates the names with a zero byte instead of a new line,
    # so that git prints every name as it is (otherwise it puts a name with spaces or
    # accented letters in quotes).
    listing = subprocess.run(["git", "-C", REPO_ROOT, "ls-files", "-z"],
                             capture_output=True, text=True, check=True)
    tracked = []
    for name in listing.stdout.split("\0"):
        if name != "":
            tracked.append(name)
    return tracked


def is_source_file(relative_path):
    """True when this tool reads the comments of the file at relative_path.

    relative_path is relative to the repository root, with / between folders. A file is
    read when its name ends in one of SOURCE_SUFFIXES or is one of SOURCE_NAMES, or when
    it is one of CONFIG_TEMPLATES. Called by source_files().
    """
    if relative_path in CONFIG_TEMPLATES:
        return True
    name = os.path.basename(relative_path)
    if name in SOURCE_NAMES:
        return True
    return name.endswith(SOURCE_SUFFIXES)


def source_files():
    """Every file in the repo whose comments this tool should read, as a sorted list.

    In a git checkout the candidates are the files git tracks (git_tracked_files), so a
    file never added with git add (a personal config, scratch output) does not change
    the result; a tracked file that has been deleted from the disk but whose deletion
    git has not been told about (no git rm yet) is left out, as there is nothing to
    read. Outside a git checkout the candidates are
    found by walking the folders, without entering those named in SKIP_DIRS. A
    candidate is kept when is_source_file() says so. Returns full paths. Called by
    collect_defined_rules(), collect_defined_globals() and check().
    """
    candidates = []                     # paths relative to the repository root
    tracked = git_tracked_files()
    if tracked is not None:
        for relative in tracked:
            if os.path.isfile(os.path.join(REPO_ROOT, relative)):
                candidates.append(relative)
    else:
        for root, dirs, files in os.walk(REPO_ROOT):
            # Leave out the folders to skip. Changing the list `dirs` in place (dirs[:] = ...)
            # is how os.walk is told not to go into them.
            dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
            for name in files:
                candidates.append(os.path.relpath(os.path.join(root, name), REPO_ROOT))

    found = []
    for relative in candidates:
        if is_source_file(relative):
            found.append(os.path.join(REPO_ROOT, relative))
    return sorted(found)


def is_snakemake_file(path):
    """True for a file of workflow code: a .smk rule file, or an entry point named in
    SOURCE_NAMES (Snakefile). Only these define rules and shell: commands.
    Called by collect_defined_rules() and comment_lines().
    """
    return path.endswith(".smk") or os.path.basename(path) in SOURCE_NAMES


def opens_shell_command(stripped):
    """True when a line opens a shell: command that goes on over the next lines.

    Takes a line without the spaces around it. Such a line starts with three double
    quotes, or with r and three double quotes (a raw string, as rule multiqc writes
    it), and holds no second set of quotes that would close the command on the same
    line. Called by comment_lines().
    """
    starts_with_quotes = stripped.startswith('"""') or stripped.startswith('r"""')
    return starts_with_quotes and stripped.count('"""') == 1


def comment_lines(path):
    """Yield (line_number, text, inside_shell_command) for every comment in a file.

    Takes the full path of one file from source_files(); check() calls it for each
    file. A comment is a line whose first character other than a space is #, and
    `text` is that line without the leading # characters and spaces. A Markdown (.md)
    file is prose from top to bottom: every line that is not blank is yielded as it
    is. Blank lines are not yielded, so the end of a paragraph shows as a jump in the
    line numbers, the same way the end of a block of comment lines does in code
    (check() relies on this to end its history skip).

    `inside_shell_command` is True for a comment line inside the command of a
    Snakemake shell: directive, between its opening and its closing triple quotes:
    such a line is part of the command text Snakemake compares (see the module
    docstring). It is tracked only in Snakemake files (is_snakemake_file), where
    shell: is a directive; in every other file it is always False. After a shell:
    line the first line that is neither blank nor a comment is the command: either
    the opening quotes (opens_shell_command) of a command of several lines, which the
    next line holding three double quotes closes, or a one-line command with nothing
    inside to count. Comment lines before the opening quotes are notes about the
    command, outside it. Known limit: a command in triple single quotes (''') is not
    recognised and its comment lines count as outside; FunFlux has no such command.
    """
    is_markdown = path.endswith(".md")
    track_shell = is_snakemake_file(path)

    # Where the current line stands relative to a shell: directive:
    #   "outside"  not in a shell: directive
    #   "waiting"  after a shell: line, before the command itself
    #   "inside"   between the opening and the closing quotes of the command
    shell_state = "outside"

    # errors="replace": a byte that is not valid UTF-8 becomes a placeholder character
    # instead of stopping the run.
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for number, raw in enumerate(handle, start=1):
            stripped = raw.strip()

            if is_markdown:
                # Every line that is not blank, as it stands: nothing to strip and no
                # shell: command to be inside.
                if stripped != "":
                    yield number, raw.rstrip("\n"), False
                continue

            if track_shell:
                if shell_state == "outside" and stripped.startswith("shell:"):
                    # The command may begin on the shell: line itself (shell: """).
                    after_shell = stripped[len("shell:"):].strip()
                    if after_shell == "" or after_shell.startswith("#"):
                        shell_state = "waiting"         # the command is on a later line
                    elif opens_shell_command(after_shell):
                        shell_state = "inside"
                    # Anything else is a one-line command, already complete.
                    continue
                if shell_state == "waiting" and stripped != "" and not stripped.startswith("#"):
                    # The first line of the command: the opening quotes, or a one-line
                    # command that is already complete.
                    if opens_shell_command(stripped):
                        shell_state = "inside"
                    else:
                        shell_state = "outside"
                    continue
                if shell_state == "inside" and '"""' in stripped:
                    shell_state = "outside"             # the closing quotes
                    continue

            if stripped.startswith("#"):
                yield number, stripped.lstrip("#").strip(), shell_state == "inside"


def collect_defined_rules():
    """Every rule name defined in the repo, as a set.

    Reads the Snakemake files among source_files() (the .smk files and the files named
    in SOURCE_NAMES, see is_snakemake_file) and keeps the name from every line of the
    form "rule <name>:". Called by check().
    """
    names = set()
    pattern = re.compile(r"^\s*rule\s+([a-z_][a-z0-9_]*)\s*:")
    for path in source_files():
        if not is_snakemake_file(path):
            continue
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                match = pattern.match(line)
                if match:
                    names.add(match.group(1))
    return names


def collect_defined_globals():
    """Every ALL_CAPS name assigned at the start of a line in the files read, as a set.

    Reads every file from source_files(), this one included, and keeps the name
    from every line that starts with an ALL_CAPS name followed by "=" (spaces
    allowed before both): NAME = ... in Python or Snakemake, NAME=... in a shell
    script. An assignment further along a line, such as the shell's
    export NAME=value, is not recognised. Called by check() for the ALL_CAPS
    advisory.
    """
    names = set()
    pattern = re.compile(r"^\s*([A-Z][A-Z0-9_]*)\s*=")
    for path in source_files():
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                match = pattern.match(line)
                if match:
                    names.add(match.group(1))
    return names


def collect_config_keys(config_path):
    """Every dotted key path present in the config file `config_path`, as a set.

    A file with
        parameters:
          busco:
            lineage: auto
    gives "parameters", "parameters.busco" and "parameters.busco.lineage". Parsed
    with PyYAML when available, and by indentation otherwise, so the check still
    works in a bare interpreter with no third-party packages. Called by check() with
    the file named by --config; main() has already stopped with exit status 2 if that
    file does not exist.
    """
    keys = set()
    try:
        # PyYAML is imported here and not at the top, so that a Python without it
        # still reaches the indentation reader below.
        import yaml
        with open(config_path, encoding="utf-8") as handle:
            data = yaml.safe_load(handle)

        def walk(node, prefix):
            # Record the full dotted path of every key in one block, then do the same
            # for each block nested inside it (walk calls itself, one level deeper).
            if not isinstance(node, dict):
                return
            for key, value in node.items():
                dotted = f"{prefix}.{key}" if prefix else str(key)
                keys.add(dotted)
                walk(value, dotted)

        walk(data, "")
    except Exception:                                           # noqa: BLE001
        # No PyYAML, or a file it cannot parse: read the keys from the indentation
        # instead, two spaces per level. `stack` holds the keys of the enclosing
        # blocks; a key at level d replaces everything from level d on. (The noqa
        # mark tells a code checker that catching every kind of error is deliberate.)
        stack = []
        with open(config_path, encoding="utf-8") as handle:
            for line in handle:
                match = re.match(r"^(\s*)([a-z][a-z0-9_]*)\s*:", line)
                if not match:
                    continue
                depth = len(match.group(1)) // 2
                stack = stack[:depth] + [match.group(2)]
                keys.add(".".join(stack))
    return keys


def tracked_paths(tracked_files):
    """Every tracked file and every folder that holds one, as a set of relative paths.

    git tracks files, not folders, so a folder counts as present when at least one
    tracked file lies inside it: "workflow/rules/shared/00_common.smk" gives itself,
    "workflow/rules/shared", "workflow/rules" and "workflow". Takes the list from
    git_tracked_files(). Called by check(); path_problem() looks paths up in the set.
    """
    paths = set()
    for relative in tracked_files:
        paths.add(relative)
        folder = os.path.dirname(relative)
        while folder != "":
            paths.add(folder)
            folder = os.path.dirname(folder)
    return paths


def ignored_by_git(relative_path):
    """True when git's ignore rules cover relative_path (the file need not exist).

    `git check-ignore` applies the repository's .gitignore files and also the user's
    own ignore rules (.git/info/exclude, and the file named by git's core.excludesFile
    setting, by default ~/.config/git/ignore). -q prints nothing: the exit status gives
    the answer, 0 for ignored. Called by path_problem() in a git checkout.
    """
    result = subprocess.run(["git", "-C", REPO_ROOT, "check-ignore", "-q", relative_path],
                            capture_output=True)
    return result.returncode == 0


def path_problem(relative_path, tracked):
    """Why a path named in a comment is not valid, or None when it is.

    relative_path is relative to the repository root; tracked is the set from
    tracked_paths() in a git checkout, and None outside one. Outside a git checkout a
    path is valid when it exists on disk. In a git checkout it is valid when git tracks
    it (a file, or a folder that holds tracked files) or ignores it (see the module
    docstring); a file that is on this computer but neither tracked nor ignored is
    reported, because a fresh clone would not have it. A tracked path is not looked
    for on disk: a file deleted with plain rm stays tracked, and so valid here, until
    the deletion is staged. Called by check().
    """
    on_disk = os.path.exists(os.path.join(REPO_ROOT, relative_path))
    if tracked is None:
        if on_disk:
            return None
        return "path does not exist"
    # normpath writes the path the way git lists it ("workflow/rules/../envs" becomes
    # "workflow/envs").
    if os.path.normpath(relative_path) in tracked:
        return None
    if ignored_by_git(relative_path):
        return None
    if on_disk:
        return "path is not tracked by git"
    return "path does not exist"


def is_path_template(text, match):
    """True when a path found in a line is only the start of a path template.

    A template names a family of files with a placeholder: config/config_{mode}.yaml,
    workflow/scripts/its_*.py. The path patterns stop at the { or the *, so the match
    is the fixed start ("config/config_"), which is not a real path. This looks at the
    one character right after the match. (A lookahead in the pattern itself, "not
    followed by { or *", would not do: the pattern would give up one more character and
    report "config/config" instead.) A placeholder after a / is not a template here: in
    workflow/rules/{mode}/ the match is the folder workflow/rules, which must exist.
    Takes the line and a match of PATH_REFERENCE or WORKFLOW_PATH_REFERENCE in it.
    Called by check().
    """
    position_after = match.end()
    if position_after >= len(text):
        return False                        # the path is the last thing on the line
    return text[position_after] in "{*"


def check(config_path):
    # Read every comment in the repo and test each name it mentions.
    # Takes the full path of the config file named by --config. Returns four lists:
    #   problems        : "path:line: message" for every reference that does not resolve
    #                     (any one makes the run fail)
    #   rule_typos      : "path:line: message" for every snake_case name in backticks that
    #                     is close to a rule name: a probable typo (information only)
    #   unknown_globals : "path:line: message" for every ALL_CAPS name that no file read
    #                     assigns at the start of a line (see collect_defined_globals;
    #                     information only)
    #   shell_comments  : "path:line" of every comment line inside a shell: command
    # Called by main(), which prints the problems, the two advisory lists (only with
    # --advisories) and, in its summary, the counts.
    rules = collect_defined_rules()
    globals_defined = collect_defined_globals()
    config_keys = collect_config_keys(config_path)

    # In a git checkout: the tracked files and folders, which a path named in a comment
    # may be (see path_problem). None outside a git checkout.
    tracked_files = git_tracked_files()
    if tracked_files is None:
        tracked = None
    else:
        tracked = tracked_paths(tracked_files)

    # Every snake_case name that legitimately is NOT a rule: config keys, output
    # column names, filenames, constants. Without this the near-miss check would
    # flag things like `amr_gene_family` (a BacFlux example) for merely resembling a
    # rule. Two sources:
    #   * every lowercase name at the start of a line followed by ":" or "=" in any file
    #     read, this file included: a rule's input, output or param name, a key in a
    #     YAML example inside a Markdown file or in a config template, a Python or shell
    #     variable, a keyword argument on its own line;
    #   * the last part of every config key (lineage, from parameters.busco.lineage).
    defined_elsewhere = set()
    for source_path in source_files():
        try:
            with open(source_path, encoding="utf-8", errors="replace") as handle:
                source_text = handle.read()
        except Exception:                                   # noqa: BLE001
            continue                        # a file that cannot be read is skipped
        # re.M lets ^ match at the start of every line, not only at the start of the file.
        assigned_names = re.findall(r"^\s*([a-z][a-z0-9_]+)\s*[:=]", source_text, re.M)
        defined_elsewhere.update(assigned_names)
    for config_key in config_keys:
        defined_elsewhere.add(config_key.split(".")[-1])

    # The top-level sections of the config file (mode, input, directories, parameters, ...).
    # A dotted name that starts with anything else (os.path, re.sub) is not a config key.
    top_level_sections = set()
    for config_key in config_keys:
        top_level_sections.add(config_key.split(".")[0])

    problems = []
    rule_typos = []
    unknown_globals = []
    shell_comments = []

    for path in source_files():
        relative = os.path.relpath(path, REPO_ROOT)
        if os.path.abspath(path) == os.path.abspath(__file__):
            continue                        # its own comments quote wrong names on purpose
        # The history skip (see HISTORICAL_MARKERS): a line with a marker starts a block,
        # and the block goes on over the lines that follow it, up to the first gap in the
        # line numbers. comment_lines() leaves a gap at every line of code in a code file
        # and at every blank line in a Markdown file, so the block ends there.
        # previous_number starts at -10 so that the first line of a file never continues
        # a block.
        in_historical_block = False
        previous_number = -10
        for number, line_text, in_shell_command in comment_lines(path):
            where = f"{relative}:{number}"
            # Counted before the history skip below: a comment inside a shell: command is
            # part of the command text whatever it says.
            if in_shell_command:
                shell_comments.append(where)

            # Web addresses out first (see WEB_ADDRESS): nothing inside one is a rule, a
            # repository path, a config key or a word about history.
            text = WEB_ADDRESS.sub(" ", line_text)

            if number != previous_number + 1:
                in_historical_block = False     # a gap ends the block
            previous_number = number
            if HISTORICAL_MARKERS.search(text):
                in_historical_block = True
            if in_historical_block:
                continue                    # deliberately describing a former name
            if ONE_LINE_MARKERS.search(text):
                continue                    # names a planned file or an old place, this line only

            # Rules: every "rule <name>" must name a defined rule. RULE_REFERENCE only
            # captures names with an underscore, so the English after "rule" in "the rule
            # always fires" never reaches this test.
            for name in RULE_REFERENCE.findall(text):
                if name not in rules:
                    problems.append(f"{where}: no rule named '{name}'")

            # Paths from the repository root: every one must exist (see path_problem).
            for match in PATH_REFERENCE.finditer(text):
                if is_path_template(text, match):
                    continue
                candidate = match.group(1)
                if candidate.endswith("envs/x.yaml"):
                    continue                                    # the made-up example a comment
                                                                # uses to explain how the
                                                                # ../../envs/ hop works
                reason = path_problem(candidate, tracked)
                if reason is not None:
                    problems.append(f"{where}: {reason}: {candidate}")

            # Paths written from inside workflow/ (rules/..., scripts/..., envs/...): the same
            # test, on the path with workflow/ put in front.
            for match in WORKFLOW_PATH_REFERENCE.finditer(text):
                if is_path_template(text, match):
                    continue
                candidate = match.group(1)
                path_from_root = "workflow/" + candidate
                reason = path_problem(path_from_root, tracked)
                if reason is not None:
                    problems.append(f"{where}: {reason}: {candidate} (read as {path_from_root})")

            # Probable rule-name typos (advisory): a snake_case name in backticks that is
            # neither a rule nor defined elsewhere, but close to a rule name.
            for match in SNAKE_NAME.findall(text):
                # SNAKE_NAME has two alternatives (backticks, diagram label), so findall
                # gives a pair per match: the name in one half, "" in the other.
                bare = match[0] or match[1]
                if bare in rules or bare in defined_elsewhere:
                    continue
                # The one rule name most similar to it, if its score reaches the cutoff.
                near = difflib.get_close_matches(bare, sorted(rules), n=1, cutoff=NEAR_MISS_CUTOFF)
                if near:
                    # ADVISORY, not a failure. Measured on BacFlux, where this was written,
                    # the ratio was about nine false alarms to one real find: output column
                    # names (plasmid_score, same_replicon) and config keys are legitimately
                    # snake_case and resemble rule names closely. As a failure it would be
                    # switched off within a week; as a list you read after editing the
                    # documentation it is worth having, because it does catch the real
                    # thing: 'check_medaka_check' for 'check_medaka_model', which no
                    # exact-match test can see.
                    rule_typos.append(
                        f"{where}: '{bare}' is not a rule; did you mean '{near[0]}'?")

            # ALL_CAPS names (advisory): each should be assigned somewhere. Only an
            # assignment at the start of a line is recognised (see collect_defined_globals).
            for name in GLOBAL_REFERENCE.findall(text):
                if name in globals_defined:
                    continue
                unknown_globals.append(
                    f"{where}: '{name}' is not assigned at the start of a line in any file read")

            # Config keys: a dotted name whose first part is a section of the config file
            # must be a key of that file.
            for key in CONFIG_REFERENCE.findall(text):
                if key.endswith(CONFIG_SKIP_SUFFIXES):
                    continue                                    # a file name, not a key
                first_part = key.split(".")[0]
                if first_part in SNAKEMAKE_DIRECTIVES:
                    continue                                    # a rule's own input:/params:
                if first_part not in top_level_sections:
                    continue                                    # not a config path at all
                if key not in config_keys:
                    problems.append(f"{where}: no config key '{key}'")

    return problems, rule_typos, unknown_globals, shell_comments


def main(arguments=None):
    # Read the options, run check(), print the results, and return the exit status:
    # 0 every reference resolves, 1 at least one does not, 2 no such config file.
    # `arguments` is the list of command-line words; None (the normal case) means the
    # words the script was started with. The tests pass their own list.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true",
                        help="print only the failures, not the summary counts")
    parser.add_argument("--advisories", action="store_true",
                        help="also list the advisories: snake_case names close to a rule "
                             "name, and ALL_CAPS names that no file read assigns at the "
                             "start of a line (noisy by nature)")
    parser.add_argument("--config", default=os.path.join("config", "config.yaml"),
                        help="the config file that dotted keys in comments are checked "
                             "against, relative to the repository root (default: config/config.yaml)")
    args = parser.parse_args(arguments)
    # --config is taken relative to the repository root, not to the current folder.
    config_path = os.path.join(REPO_ROOT, args.config)
    if not os.path.exists(config_path):
        print(f"error: config file not found: {config_path}")
        return 2

    problems, rule_typos, unknown_globals, shell_comments = check(config_path)

    # The failures are always printed, one per line; the advisories only with --advisories.
    for problem in problems:
        print(problem)
    if args.advisories:
        for advisory in rule_typos:
            print(advisory)
        for advisory in unknown_globals:
            print(advisory)

    # The summary, unless --quiet: one count per advisory kind, the shell: comment count,
    # and the verdict.
    if not args.quiet:
        print()
        print(f"{len(rule_typos)} snake_case name(s) in backticks or diagram labels are close")
        print("  to a rule name but are not one: a misspelt rule, or a planned rule, a value")
        print("  or a file name that only looks like one — listed with --advisories, not failed.")
        print()
        print(f"{len(unknown_globals)} ALL_CAPS name(s) in comments are not assigned at the")
        print("  start of a line (NAME = ...) in any file read; a shell 'export NAME=...'")
        print("  is not recognised. Most are environment or shell variables, gene, contig")
        print("  or read-file names, or accessions, which is fine — listed with")
        print("  --advisories, not failed.")
        print()
        print(f"{len(shell_comments)} comment line(s) sit inside a shell: command.")
        print("  Those are part of the command text Snakemake compares: rewording one")
        print("  makes the rule re-run. Cheap for a download, expensive for gene")
        print("  prediction or InterProScan.")
        print()
        if problems:
            print(f"FAIL: {len(problems)} reference(s) in comments do not resolve.")
        else:
            print("OK: every rule, path and config key named in a comment exists.")

    return 1 if problems else 0


# Run main() when the file is started as a script, and hand its result to the shell as
# the exit status.
if __name__ == "__main__":
    sys.exit(main())
