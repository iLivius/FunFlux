"""Tests for check_comment_references.py.

What these tests protect
------------------------
check_comment_references.py reads the comments and the documentation pages of
the repository and fails when one names a rule, a path or a config key that
does not exist. A checker with a blind spot gives false confidence: it says OK
while broken references sit in the part it does not read, or it cries wolf on
correct text until nobody runs it. The reviews of 2026-10-08 found such blind
spots (they are listed at the top of the checker). Each section below pins one
fix, so that the blind spot cannot come back unnoticed, and the last section
checks that a real broken reference is still reported.

How the tests work
------------------
Every test builds a tiny fake repository in a throw-away temporary folder
(pytest's tmp_path): a rule file, a config file, a page or two. The checker
takes every path relative to its REPO_ROOT, so each test points REPO_ROOT at
the fake repository (pytest's monkeypatch puts the real value back after the
test). Most tests then call check() and look at the lists it returns; the tests
of the printed text call main() exactly as the command line would. The tests of
the git behaviour turn the fake repository into a git checkout first (git init,
git add). git is given a stand-in home folder for these tests, so the tester's
own git settings, such as a personal list of ignored files, cannot change the
outcome.

Run with:  python -m pytest miscellaneous -q
"""

import os
import subprocess

import check_comment_references as ccr


# =============================================================================
# Helpers shared by the tests
# =============================================================================

# The config file of every fake repository: three top-level sections, so that a
# dotted name starting with mode, directories or parameters is checked against it.
CONFIG_TEXT = (
    "mode: contigs\n"
    "directories:\n"
    "  antismash_db: \"\"\n"
    "parameters:\n"
    "  busco:\n"
    "    lineage: auto\n"
)


def write_file(root, relative_path, text):
    """Write text into root/relative_path, creating folders as needed.

    Input:  a root folder (the fake repository), a path relative to it, and the
            file content. Returns the full path of the file.
    """
    full_path = os.path.join(str(root), relative_path)
    folder = os.path.dirname(full_path)
    if folder != "":
        os.makedirs(folder, exist_ok=True)
    with open(full_path, "w") as handle:
        handle.write(text)
    return full_path


def make_repository(tmp_path, monkeypatch, files):
    """Build a fake repository in tmp_path/repo and point the checker at it.

    Input:  files maps a path inside the repository to the text of that file.
            config/config.yaml (CONFIG_TEXT) is always written first, because
            main() stops when the file named by --config is missing; `files`
            may replace it.
    Output: the repository folder.
    """
    repo = tmp_path / "repo"
    repo.mkdir()
    write_file(repo, "config/config.yaml", CONFIG_TEXT)
    for relative_path, text in files.items():
        write_file(repo, relative_path, text)
    monkeypatch.setattr(ccr, "REPO_ROOT", str(repo))
    return repo


def isolate_git(tmp_path, monkeypatch):
    """Give git a stand-in home folder and no system-wide or inherited settings.

    monkeypatch.setenv and delenv change the environment for this test only, and
    the checker's own git calls see the same environment.
    """
    home = tmp_path / "home"
    home.mkdir(exist_ok=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    for name in ("XDG_CONFIG_HOME", "GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE"):
        monkeypatch.delenv(name, raising=False)


def run_git(folder, *arguments):
    """Run one git command in `folder`; the test stops with an error if git fails."""
    subprocess.run(["git", "-C", str(folder)] + list(arguments),
                   check=True, capture_output=True, text=True)


def make_git_checkout(repo, tmp_path, monkeypatch):
    """Turn the fake repository into a git checkout and track every file in it.

    Nothing is committed: git ls-files lists a file as soon as it is added.
    Files written after this call are untracked.
    """
    isolate_git(tmp_path, monkeypatch)
    run_git(repo, "init", "-q")
    run_git(repo, "add", "-A")


def run_check(repo):
    """Run check() on the fake repository, against its config/config.yaml.

    Output: the four lists check() returns: problems, rule_typos,
            unknown_globals, shell_comments.
    """
    return ccr.check(str(repo / "config" / "config.yaml"))


def run_main(capsys, *options):
    """Run main() with its default --config (config/config.yaml) and any extra options.

    Output: (exit_status, output_lines): what main() returned and the lines it
            printed.
    """
    exit_status = ccr.main(list(options))
    # main() must hand back a real integer, so that the shell gets an exit status.
    assert isinstance(exit_status, int)
    output_lines = capsys.readouterr().out.splitlines()
    return exit_status, output_lines


# =============================================================================
# 1a. Comment lines inside a shell: command are counted, and only those
# =============================================================================


def test_comments_between_the_quotes_are_counted_not_the_notes_before_them(tmp_path, monkeypatch):
    """The first quotes after shell: open the command; the next ones close it.

    Before the fix the opening quotes closed the block at once: the comments
    inside the command were missed, and the note written between shell: and the
    quotes, which is not part of the command, was counted instead.
    """
    rule_text = (
        "rule demo_step:\n"                                  # line 1
        "    shell:\n"                                       # line 2
        "        # a note about the command, outside it\n"   # line 3
        "        \"\"\"\n"                                   # line 4: opens
        "        # step 1, inside the command\n"             # line 5
        "        echo one\n"                                 # line 6
        "        # step 2, inside the command\n"             # line 7
        "        \"\"\"\n"                                   # line 8: closes
        "# a comment after the rule\n"                       # line 9
    )
    repo = make_repository(tmp_path, monkeypatch, {"workflow/rules/demo.smk": rule_text})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert shell_comments == ["workflow/rules/demo.smk:5", "workflow/rules/demo.smk:7"]


def test_raw_string_command_and_comment_about_v1_inside_it(tmp_path, monkeypatch):
    """r\"\"\" opens a command too (rule multiqc writes it so).

    A comment inside the command that mentions v1 is still counted: the history
    skip decides what is checked, not what belongs to the command text.
    """
    rule_text = (
        "rule demo_report:\n"                   # line 1
        "    shell:\n"                          # line 2
        "        # (1) a long note\n"            # line 3
        "        #\n"                            # line 4
        "        r\"\"\"\n"                      # line 5: opens
        "        printf \"%s\\n\" x\n"           # line 6
        "        # kept as in v1\n"              # line 7
        "        \"\"\"\n"                       # line 8: closes
    )
    repo = make_repository(tmp_path, monkeypatch, {"workflow/rules/demo.smk": rule_text})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert shell_comments == ["workflow/rules/demo.smk:7"]


def test_one_line_command_and_command_on_the_shell_line(tmp_path, monkeypatch):
    """A one-line command holds no comment; shell: \"\"\" opens a command at once."""
    rule_text = (
        "rule demo_touch:\n"                                     # line 1
        "    shell:\n"                                           # line 2
        "        \"touch {output}\"\n"                           # line 3
        "# after the rule, outside any command\n"                # line 4
        "rule demo_echo:\n"                                      # line 5
        "    shell: \"\"\"\n"                                    # line 6: opens
        "        # inside a command that starts on the shell: line\n"   # line 7
        "        echo two\n"                                     # line 8
        "        \"\"\"\n"                                       # line 9: closes
        "# after the second rule\n"                              # line 10
    )
    repo = make_repository(tmp_path, monkeypatch, {"workflow/rules/demo.smk": rule_text})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert shell_comments == ["workflow/rules/demo.smk:7"]


def test_one_line_command_in_triple_quotes_and_note_on_the_shell_line(tmp_path, monkeypatch):
    """shell: \"\"\"touch x\"\"\" is a complete one-line command: the opening and the
    closing quotes are on the same line, so the comment after it is outside any
    command. A note on the shell: line itself (shell:  # ...) is not a command: the
    command starts on a later line, and the comments inside it count."""
    rule_text = (
        "rule demo_touch:\n"                                         # line 1
        "    shell: \"\"\"touch {output}\"\"\"\n"                    # line 2: one line
        "# after the rule, outside any command\n"                    # line 3
        "rule demo_note:\n"                                          # line 4
        "    shell:  # what the command below does\n"                # line 5: a note
        "        \"\"\"\n"                                           # line 6: opens
        "        # inside the command\n"                             # line 7
        "        \"\"\"\n"                                           # line 8: closes
    )
    repo = make_repository(tmp_path, monkeypatch, {"workflow/rules/demo.smk": rule_text})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert shell_comments == ["workflow/rules/demo.smk:7"]


def test_summary_prints_the_number_of_comments_inside_commands(tmp_path, monkeypatch, capsys):
    """The printed count is the count of comment lines inside the commands."""
    rule_text = (
        "rule demo_step:\n"
        "    shell:\n"
        "        # outside\n"
        "        \"\"\"\n"
        "        # inside 1\n"
        "        # inside 2\n"
        "        \"\"\"\n"
    )
    make_repository(tmp_path, monkeypatch, {"workflow/rules/demo.smk": rule_text})
    exit_status, output_lines = run_main(capsys)
    assert exit_status == 0
    assert "2 comment line(s) sit inside a shell: command." in output_lines


# =============================================================================
# 1b. Which files are read: the Snakefile and the two config templates yes,
#     the conda environment files and copies of a template kept elsewhere no
# =============================================================================


def test_snakefile_and_the_config_templates_are_read(tmp_path, monkeypatch):
    """A broken path in any of the three files is reported, and a rule defined in
    the Snakefile counts as a rule."""
    files = {
        "workflow/Snakefile": "# see docs/missing_one.md\nrule entry_point:\n    input: []\n",
        "config/config.yaml": CONFIG_TEXT + "# see docs/missing_two.md\n",
        "config/config_test_examples.yaml": "# see docs/missing_three.md\nmode: contigs\n",
        "docs/index.md": "The last rule is rule entry_point.\n",
    }
    repo = make_repository(tmp_path, monkeypatch, files)
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == [
        "config/config.yaml:7: path does not exist: docs/missing_two.md",
        "config/config_test_examples.yaml:1: path does not exist: docs/missing_three.md",
        "workflow/Snakefile:1: path does not exist: docs/missing_one.md",
    ]


def test_environment_files_and_template_copies_are_not_read(tmp_path, monkeypatch):
    """Comments in workflow/envs/ are never read, nor those of a copy of a config
    template kept elsewhere: only the two template paths count."""
    files = {
        "workflow/envs/tool.yaml": "# see docs/missing_env.md\nname: tool\n",
        "workflow/envs/tool.linux-64.pin.txt": "# see docs/missing_pin.md\n",
        "configs/config.yaml": "# see docs/missing_copy.md\nmode: x\n",
    }
    repo = make_repository(tmp_path, monkeypatch, files)
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == []


# =============================================================================
# 1c. In a git checkout only tracked files are read, and ignored paths are valid
# =============================================================================


def test_untracked_file_is_not_read_in_a_git_checkout(tmp_path, monkeypatch):
    """A local file never added with git add does not change the result; once it
    is added (no commit needed), the same file is read."""
    repo = make_repository(tmp_path, monkeypatch, {"docs/index.md": "Nothing to check.\n"})
    make_git_checkout(repo, tmp_path, monkeypatch)
    write_file(repo, "docs/scratch_notes.md", "See docs/missing.md.\n")

    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == []

    run_git(repo, "add", "docs/scratch_notes.md")
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == ["docs/scratch_notes.md:1: path does not exist: docs/missing.md"]


def test_tracked_files_and_folders_count_as_present(tmp_path, monkeypatch):
    """In a git checkout a path is valid when git tracks it: a file, a folder that
    holds tracked files (git lists only files), a path written from inside
    workflow/, and a path with ../ in it that leads to a tracked file (the way a
    rule file reaches the environment files). The broken path on the last line
    shows that the page was read."""
    page = (
        "Start at docs/index.md.\n"                                          # line 1: a file
        "The rule files sit in workflow/rules and workflow/rules/shared.\n"  # line 2: folders
        "One of them is rules/shared/real.smk.\n"                            # line 3: rules/...
        "It reaches workflow/rules/shared/../../envs/tool.yaml.\n"           # line 4: with ../
        "See docs/missing.md.\n"                                             # line 5: broken
    )
    files = {
        "docs/index.md": page,
        "workflow/rules/shared/real.smk": "rule real_step:\n    input: []\n",
        "workflow/envs/tool.yaml": "name: tool\n",
    }
    repo = make_repository(tmp_path, monkeypatch, files)
    make_git_checkout(repo, tmp_path, monkeypatch)
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == ["docs/index.md:5: path does not exist: docs/missing.md"]


def test_tracked_file_deleted_with_rm_is_not_read(tmp_path, monkeypatch):
    """A tracked file removed with plain rm (the deletion not yet staged) is still
    in git's list, but there is nothing left to read: it is skipped, and the run
    does not stop on it."""
    files = {
        "docs/index.md": "Nothing to check.\n",
        "docs/old_notes.md": "See docs/missing.md.\n",
    }
    repo = make_repository(tmp_path, monkeypatch, files)
    make_git_checkout(repo, tmp_path, monkeypatch)
    os.remove(repo / "docs" / "old_notes.md")
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == []


def test_path_ignored_by_git_counts_as_present(tmp_path, monkeypatch):
    """config/config_custom.yaml, listed in .gitignore, may be named; a path that is
    neither tracked nor ignored may not."""
    files = {
        ".gitignore": "config/config_custom.yaml\n",
        "docs/index.md": "Copy the template to config/config_custom.yaml.\n"
                         "\n"
                         "Then read config/config_other.yaml.\n",
    }
    repo = make_repository(tmp_path, monkeypatch, files)
    make_git_checkout(repo, tmp_path, monkeypatch)
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == ["docs/index.md:3: path does not exist: config/config_other.yaml"]


def test_file_on_disk_but_not_tracked_is_reported_as_such(tmp_path, monkeypatch):
    """A fresh clone would not have it, so naming it is an error, with its own message."""
    repo = make_repository(tmp_path, monkeypatch, {"docs/index.md": "See docs/local_notes.md.\n"})
    make_git_checkout(repo, tmp_path, monkeypatch)
    write_file(repo, "docs/local_notes.md", "my notes\n")
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == ["docs/index.md:1: path is not tracked by git: docs/local_notes.md"]


def test_outside_a_git_checkout_every_file_is_read(tmp_path, monkeypatch):
    """Without git, the folders are walked: every file counts, and an ignored path
    must exist on disk like any other."""
    files = {
        ".gitignore": "config/config_custom.yaml\n",
        "docs/scratch_notes.md": "See docs/missing.md.\n"
                                 "\n"
                                 "Copy the template to config/config_custom.yaml.\n",
    }
    repo = make_repository(tmp_path, monkeypatch, files)
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == [
        "docs/scratch_notes.md:1: path does not exist: docs/missing.md",
        "docs/scratch_notes.md:3: path does not exist: config/config_custom.yaml",
    ]


def test_folder_inside_another_git_checkout_is_walked(tmp_path, monkeypatch):
    """When the repository folder sits inside some other checkout (which does not
    track it), that checkout's empty file list must not be used: the folders are
    walked instead, so the page is read and its broken path reported."""
    isolate_git(tmp_path, monkeypatch)
    run_git(tmp_path, "init", "-q")
    repo = make_repository(tmp_path, monkeypatch, {"docs/index.md": "See docs/missing.md.\n"})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == ["docs/index.md:1: path does not exist: docs/missing.md"]


# =============================================================================
# 1d. Markdown: the history skip ends at a blank line; "1.x" is history
# =============================================================================


def test_markdown_history_skip_ends_at_a_blank_line(tmp_path, monkeypatch):
    """Before the fix one 'v1' switched checking off for the rest of the page."""
    page = (
        "In v1 the page was docs/old_page.md.\n"        # line 1: history
        "It also linked docs/old_index.md.\n"           # line 2: same paragraph, skipped
        "\n"                                            # line 3: the paragraph ends
        "The current page is docs/missing_page.md.\n"   # line 4: checked
    )
    repo = make_repository(tmp_path, monkeypatch, {"docs/page.md": page})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == ["docs/page.md:4: path does not exist: docs/missing_page.md"]


def test_funflux_1x_counts_as_history(tmp_path, monkeypatch):
    """'FunFlux 1.x' and '1.x' mark a line about the v1 workflows, like 'v1'."""
    page = (
        "FunFlux 1.x wrote docs/old_layout.md.\n"
        "\n"
        "Version 1.x kept docs/older_layout.md.\n"
    )
    repo = make_repository(tmp_path, monkeypatch, {"docs/page.md": page})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == []


def test_planned_and_moved_mark_only_their_own_line(tmp_path, monkeypatch):
    """'planned', 'moved from' and 'moved into' skip their line, not the next one:
    in a table the next row is about something else."""
    page = (
        "| `rules/shared/later.smk` (planned) |\n"                            # line 1
        "| `rules/shared/missing.smk` |\n"                                    # line 2
        "input:   # MOVED from directories.input_dir\n"                       # line 3
        "parameters.genus: this setting moved into parameters.busco.lineage\n"   # line 4
        "parameters.busco.lineages is a typo\n"                               # line 5
    )
    repo = make_repository(tmp_path, monkeypatch, {"docs/page.md": page})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == [
        "docs/page.md:2: path does not exist: rules/shared/missing.smk "
        "(read as workflow/rules/shared/missing.smk)",
        "docs/page.md:5: no config key 'parameters.busco.lineages'",
    ]


# =============================================================================
# 1e. Paths written from inside workflow/ are looked up under workflow/
# =============================================================================


def test_paths_written_from_inside_workflow_are_checked(tmp_path, monkeypatch):
    """rules/..., scripts/... and envs/... are read as workflow/rules/... and so on.

    The 'rules/...' inside 'workflow/rules/...' is not checked a second time, and
    the 'envs/...' in '../../envs/...' (relative to a rule file) is not read.
    """
    readme = (
        "| `rules/shared/real.smk` | exists |\n"                       # line 1
        "| `rules/shared/missing.smk` | missing |\n"                   # line 2
        "| `scripts/none.py` | missing |\n"                            # line 3
        "| `envs/none.yaml` | missing |\n"                             # line 4
        "The full path workflow/rules/shared/real.smk exists.\n"       # line 5
        "A rule file reaches ../../envs/other.yaml from its folder.\n"  # line 6
    )
    files = {
        "workflow/rules/shared/real.smk": "rule real_step:\n    input: []\n",
        "workflow/scripts/README.md": readme,
    }
    repo = make_repository(tmp_path, monkeypatch, files)
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == [
        "workflow/scripts/README.md:2: path does not exist: rules/shared/missing.smk "
        "(read as workflow/rules/shared/missing.smk)",
        "workflow/scripts/README.md:3: path does not exist: scripts/none.py "
        "(read as workflow/scripts/none.py)",
        "workflow/scripts/README.md:4: path does not exist: envs/none.yaml "
        "(read as workflow/envs/none.yaml)",
    ]


# =============================================================================
# 1f. Path templates are skipped
# =============================================================================


def test_path_templates_are_skipped(tmp_path, monkeypatch):
    """A { or * right after the match marks a template. Before the fix
    'config/config_{mode}.yaml' failed as 'config/config_'. A placeholder after a
    / leaves the folder before it checked."""
    page = (
        "Copy config/config_{mode}.yaml, run workflow/scripts/its_*.py.\n"   # line 1
        "\n"
        "The rule files sit in workflow/rulez/{mode}/.\n"                     # line 3
    )
    repo = make_repository(tmp_path, monkeypatch, {"docs/page.md": page})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == ["docs/page.md:3: path does not exist: workflow/rulez"]


# =============================================================================
# 1g. Text inside a web address is not checked
# =============================================================================


def test_paths_inside_web_addresses_are_ignored(tmp_path, monkeypatch):
    """https://geneontology.org/docs/go-citation-policy/ is not the repository folder
    docs/go-citation-policy; the rest of the line is still checked. A 'v1' inside a
    web address does not make the line a line about history."""
    page = (
        "Policy: https://geneontology.org/docs/go-citation-policy/, "
        "http://example.org/workflow/x.smk, ftp://ftp.example.org/config/y.yaml "
        "and docs/missing.md.\n"                                                    # line 1
        "\n"                                                                        # line 2
        "The notes at https://example.org/v1/notes.html are in docs/notes.md.\n"    # line 3
    )
    repo = make_repository(tmp_path, monkeypatch, {"docs/page.md": page})
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert problems == [
        "docs/page.md:1: path does not exist: docs/missing.md",
        "docs/page.md:3: path does not exist: docs/notes.md",
    ]


# =============================================================================
# 1h. The shell: tracker runs only on Snakemake files, never on Markdown
# =============================================================================


def test_shell_tracker_does_not_run_on_markdown_or_python(tmp_path, monkeypatch):
    """In a page, every line that is not blank is read, the shell: example included;
    a Python file with a 'shell:' line opens no command either."""
    page = (
        "```\n"                              # line 1
        "shell:\n"                           # line 2
        "    \"\"\"\n"                       # line 3
        "    # a comment in an example\n"    # line 4
        "    \"\"\"\n"                       # line 5
        "```\n"                              # line 6
    )
    # In a Python file, "shell:" followed by three quotes is only text in a docstring:
    # the comment after it must not be counted as inside a command.
    script = (
        "def example():\n"
        "    \"\"\"Show a rule.\n"
        "shell:\n"
        "    \"\"\"\n"
        "    # a plain comment\n"
        "    return 1\n"
    )
    repo = make_repository(tmp_path, monkeypatch, {"docs/page.md": page,
                                                   "workflow/scripts/demo.py": script})
    page_lines = list(ccr.comment_lines(str(repo / "docs" / "page.md")))
    line_numbers = [number for number, text, inside in page_lines]
    assert line_numbers == [1, 2, 3, 4, 5, 6]
    problems, rule_typos, unknown_globals, shell_comments = run_check(repo)
    assert shell_comments == []


# =============================================================================
# 1i. The filters that were removed could never have fired
# =============================================================================


def test_removed_filters_could_never_trigger(tmp_path, monkeypatch, capsys):
    """Proof, kept as a test, that each removed filter was dead code.

    * RULE_REFERENCE captures only names with an underscore, so none of the
      English words the rule check used to filter out can ever be captured.
    * The path patterns cannot capture a { or a * (their character class has
      neither), so a test for them inside the match never fires.
    * CONFIG_REFERENCE cannot capture a /.
    * main() stops with exit status 2 before check() runs when the config file is
      missing, so collect_config_keys() never sees a missing file.
    """
    english_words = ["all", "above", "below", "that", "the", "and", "in", "for",
                     "is", "it", "this", "which", "with"]
    for word in english_words:
        assert ccr.RULE_REFERENCE.findall("rule " + word) == []
        assert ccr.RULE_REFERENCE.findall("rules " + word) == []

    template_text = "config/config_{mode}.yaml workflow/scripts/its_*.py rules/{mode}/x.smk"
    for pattern in (ccr.PATH_REFERENCE, ccr.WORKFLOW_PATH_REFERENCE):
        for match in pattern.finditer(template_text):
            assert "{" not in match.group(1)
            assert "*" not in match.group(1)

    for key in ccr.CONFIG_REFERENCE.findall("parameters.busco/lineage docs/a.b.c x.y/z"):
        assert "/" not in key

    def check_must_not_run(config_path):
        raise AssertionError("check() ran although the config file is missing")

    make_repository(tmp_path, monkeypatch, {})
    monkeypatch.setattr(ccr, "check", check_must_not_run)
    exit_status = ccr.main(["--config", "config/no_such_file.yaml"])
    assert exit_status == 2
    assert "error: config file not found" in capsys.readouterr().out


# =============================================================================
# 1j. The printed summary says what each count holds
# =============================================================================


def test_summary_counts_each_advisory_kind_on_its_own_line(tmp_path, monkeypatch, capsys):
    """One probable rule typo and two ALL_CAPS names give the counts 1 and 2; no
    BacFlux tool is named, and the OK line does not claim the ALL_CAPS names are
    checked.

    The texts say that only an assignment at the start of a line is recognised:
    the name a Python file assigns (ASSIGNED_SETTING = 1) is not listed, while
    the name a shell script sets with export is, like the one set nowhere.
    """
    files = {
        "workflow/rules/demo.smk": "rule check_medaka_model:\n    input: []\n",
        "workflow/scripts/demo.py": "ASSIGNED_SETTING = 1\n",
        "workflow/scripts/demo.sh": "export EXPORTED_SETTING=1\n",
        "docs/page.md": "Run `check_medaka_check` with SOME_SETTING, ASSIGNED_SETTING "
                        "and EXPORTED_SETTING set.\n",
    }
    make_repository(tmp_path, monkeypatch, files)
    exit_status, output_lines = run_main(capsys)
    assert exit_status == 0
    output_text = "\n".join(output_lines)
    assert "1 snake_case name(s) in backticks or diagram labels are close" in output_text
    assert "2 ALL_CAPS name(s) in comments are not assigned at the" in output_text
    assert "a shell 'export NAME=...'" in output_text
    assert "GTDB-Tk" not in output_text
    assert "expensive for gene" in output_text
    assert output_lines[-1] == "OK: every rule, path and config key named in a comment exists."

    # With --advisories the advisories are listed, the typo first.
    exit_status, output_lines = run_main(capsys, "--advisories", "--quiet")
    assert output_lines == [
        "docs/page.md:1: 'check_medaka_check' is not a rule; did you mean 'check_medaka_model'?",
        "docs/page.md:1: 'SOME_SETTING' is not assigned at the start of a line "
        "in any file read",
        "docs/page.md:1: 'EXPORTED_SETTING' is not assigned at the start of a line "
        "in any file read",
    ]


# =============================================================================
# 2. A real broken reference is still reported
# =============================================================================


def test_real_broken_references_are_reported(tmp_path, monkeypatch, capsys):
    """A wrong rule name, a missing file and a wrong config key, in a git checkout:
    three failures, exit status 1."""
    rule_text = (
        "# rule filter_contig writes docs/missing.md, set by parameters.busco.lineages\n"
        "rule filter_contigs:\n"
        "    input: []\n"
    )
    repo = make_repository(tmp_path, monkeypatch, {"workflow/rules/demo.smk": rule_text})
    make_git_checkout(repo, tmp_path, monkeypatch)
    exit_status, output_lines = run_main(capsys)
    assert exit_status == 1
    assert output_lines[:3] == [
        "workflow/rules/demo.smk:1: no rule named 'filter_contig'",
        "workflow/rules/demo.smk:1: path does not exist: docs/missing.md",
        "workflow/rules/demo.smk:1: no config key 'parameters.busco.lineages'",
    ]
    assert output_lines[-1] == "FAIL: 3 reference(s) in comments do not resolve."
