"""Install a user-supplied official SignalP 6 archive for one workflow run.

SignalP 6 finds signal peptides: the short stretch at the start of a protein
that sends it into the secretory pathway, the route that secreted proteins take
out of the cell. funannotate annotate uses its result to list the secreted
proteins. SignalP 6 is licensed software of the Technical University of Denmark
(DTU), handed out after the licence has been accepted, so this script never
downloads, redistributes, or hard-codes a DTU URL. It installs exactly the
archive the user names in files.signalp6_archive.

Where it fits in the workflow
-----------------------------
Run by rule install_signalp6 (rules/shared/45_signalp6.smk), once per run and
only when parameters.signalp6.mode is not off. It is a ``script:`` rule:
Snakemake runs this file itself, inside the rule's conda environment
(envs/signalp6.yaml: PyTorch, NumPy and the other packages SignalP needs, but
not SignalP itself), and hands it an object called ``snakemake`` that holds the
rule's input, output, params and log.

What it does, in order
----------------------
  1. checks the settings, computes the archive's SHA-256 checksum and, if the
     user gave one (parameters.signalp6.archive_sha256), compares the two;
  2. unpacks the archive into a temporary folder, refusing any member that
     could write outside it;
  3. finds the SignalP Python package in it, and the model for the chosen mode;
  4. installs only that package into the rule's conda environment (pip, with
     nothing downloaded), and starts ``signalp6 --help`` once to show that the
     new command runs;
  5. moves the unpacked package, which holds the large model file, to the run's
     ``00.resources/signalp6/package`` folder;
  6. writes the manifest, ``00.resources/signalp6/installation.json``, which
     rule signalp6_annotation reads to find the model.

The manifest records where the model is as a path relative to the manifest
itself, so the run's output folder can be moved to another computer.

The manifest also records ``conda_prefix``, the folder of the conda
environment the package was installed into. That one is an absolute
path on the computer that ran the installer and is kept for diagnosis only;
nothing uses it to find files.

Everything this script prints, and everything the programs started through
run() print (pip, signalp6 --help), is also written to the rule's log file,
because Snakemake does not do that for a ``script:`` rule.

Needs Python 3.9 or newer (the environment has 3.10): the type hint list[str] in
run() is evaluated as soon as the script starts, and older Pythons stop there
with a TypeError.
"""

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path, PurePosixPath


# ── Helper functions (used by the installation steps at the bottom) ──────────


class Tee:
    """Write to two places at once: the original stream (the screen) and the rule's log file."""

    # stream    the original sys.stdout or sys.stderr
    # log_file  the rule's log file, already open for writing
    # Created twice at the start of the installation (below), to replace
    # sys.stdout and sys.stderr.
    def __init__(self, stream, log_file):
        self.stream = stream
        self.log_file = log_file

    # Called by print() and by Python's error report. Returns the number of
    # characters written, as a real file does.
    def write(self, text):
        self.stream.write(text)
        self.log_file.write(text)
        return len(text)

    # Push anything still held in memory out to the screen and to the log file.
    def flush(self):
        self.stream.flush()
        self.log_file.flush()


# Stop the installation with a message that starts "SignalP 6 installation
# failed:". It raises an error and never returns. Nothing in this script catches
# that error, so Python prints it (through the Tee, so into the log as well) and
# exits with status 1, and Snakemake marks the rule as failed. Called by every
# check in this file.
def fail(message: str) -> None:
    raise RuntimeError(f"SignalP 6 installation failed: {message}")


# Return the SHA-256 checksum of a file, as 64 hexadecimal characters. A checksum
# is a fingerprint of a file: if even one byte differs, the checksum differs.
# The file is read 1 MB at a time, so a large archive never has to fit in
# memory: iter(function, b"") calls the reading function again and again until
# it returns b"", the empty read at the end of the file. Called once, on the
# archive.
def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_extract(archive: Path, destination: Path) -> None:
    """Unpack the .tar.gz archive into `destination`, folders and regular files only.

    EVERY member is checked before anything is unpacked. A member whose path
    is absolute, or contains ".." (one folder up), could be written outside
    the destination folder; a link could point outside it; a special file (a
    device, a pipe) has no place in a software archive. Any of these stops the
    installation. Called once, with the temporary folder as destination.
    """
    with tarfile.open(archive, "r:gz") as tar:
        members = tar.getmembers()
        for member in members:
            member_path = PurePosixPath(member.name)
            if member_path.is_absolute() or ".." in member_path.parts:
                fail(f"archive member has an unsafe path: {member.name!r}")
            if not (member.isdir() or member.isfile()):
                fail(
                    "archive contains an unsupported link or special file: "
                    f"{member.name!r}"
                )
        tar.extractall(destination, members=members)


# Return the folder of SignalP's Python package inside the unpacked archive: the
# one folder, at any depth, that holds setup.py, requirements.txt and a signalp/
# subfolder (the DTU archive wraps it in further folders, for example
# signalp6_fast/signalp-6-package/). No such folder, or more than one, stops the
# installation. Called once, on the unpacked archive.
def find_package_root(extracted: Path) -> Path:
    candidates = sorted(
        path.parent
        for path in extracted.rglob("setup.py")
        if (path.parent / "requirements.txt").is_file()
        and (path.parent / "signalp").is_dir()
    )
    if len(candidates) != 1:
        fail(
            "expected exactly one SignalP Python package (setup.py, "
            f"requirements.txt and signalp/); found {len(candidates)}"
        )
    return candidates[0]


def find_model_dir(package_root: Path, mode: str) -> Path:
    """Return the folder inside the package that holds the model for `mode`.

    What is looked for depends on the mode: distilled_model_signalp6.pt
    (fast), ensemble_model_signalp6.pt (slow) or sequential_models_signalp6
    (slow-sequential). DTU documents two places for it, and both are tried:
    models/ in the main package, and signalp/model_weights/ in a model
    downloaded separately for another mode (the README inside the package
    describes both). Both are inside the unpacked archive; no path on this
    computer is assumed. If neither holds the model, the installation stops
    (fail() never returns, so the function always returns a folder). Called
    once, after find_package_root.
    """
    required = {
        "fast": "distilled_model_signalp6.pt",
        "slow": "ensemble_model_signalp6.pt",
        "slow-sequential": "sequential_models_signalp6",
    }[mode]
    candidates = [
        package_root / "models",
        package_root / "signalp" / "model_weights",
    ]
    for candidate in candidates:
        if (candidate / required).exists():
            return candidate
    searched = ", ".join(str(candidate.relative_to(package_root)) for candidate in candidates)
    fail(
        f"the supplied archive does not contain the {mode!r} model ({required}) "
        f"in {searched}. Download the DTU archive for the selected mode."
    )


# Run one program. `command` is the program and its arguments as a list of
# words. The command is printed first, after "+ " (as bash does with set -x).
# Its error messages are merged into its normal output, which is passed on line
# by line. A non-zero exit status raises CalledProcessError, which stops the
# installation. Called for pip and for signalp6 --help.
def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    # The program's output is read line by line and printed, so that it reaches the log file as
    # well as the screen (subprocess.run on its own would write it to the screen only).
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    for line in process.stdout:
        print(line, end="", flush=True)
    process.wait()
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, command)


# ── The installation, run from top to bottom when Snakemake starts the script ─
# `snakemake` is not imported: Snakemake puts it into the script before running
# it. Its input, output, params and log are those of rule install_signalp6.

# Copy everything printed from here on into the rule's log file (see the module docstring).
# buffering=1 writes every finished line to the log file at once, so the log is
# complete up to the last line even if the job is killed. The rule always
# declares a log; without one, the output goes to the screen only.
if snakemake.log:
    log_file = open(str(snakemake.log[0]), "w", buffering=1, encoding="utf-8")
    sys.stdout = Tee(sys.stdout, log_file)
    sys.stderr = Tee(sys.stderr, log_file)

# The archive (files.signalp6_archive) and the manifest to write
# (00.resources/signalp6/installation.json), both as absolute paths; then the
# mode (fast, slow or slow-sequential) and the checksum the user expects (empty
# when none was given), in lower case.
archive = Path(str(snakemake.input.archive)).resolve()
manifest = Path(str(snakemake.output.manifest)).resolve()
params = snakemake.params
mode = str(params.mode).strip().lower()
expected_sha256 = str(params.archive_sha256 or "").strip().lower()

# Step 1: check the settings. 00_common.smk has already checked the archive and
# the mode when Snakemake read the workflow; the form of the checksum (64
# characters, 0-9 and a-f) is checked only here.
if not archive.is_file():
    fail(f"configured archive does not exist or is not a regular file: {archive}")
if mode not in {"fast", "slow", "slow-sequential"}:
    fail(f"unsupported mode {mode!r}")
if expected_sha256 and (
    len(expected_sha256) != 64
    or any(character not in "0123456789abcdef" for character in expected_sha256)
):
    fail("archive_sha256 must be a 64-character hexadecimal SHA-256 digest")

# Compute the archive's checksum before installing anything. It is compared with
# the checksum the user gave (if any), and it is written into the manifest, so
# the run records exactly which archive was installed.
archive_sha256 = sha256sum(archive)
if expected_sha256 and archive_sha256 != expected_sha256:
    fail(
        "configured archive_sha256 does not match the supplied archive "
        f"(expected {expected_sha256}, observed {archive_sha256})"
    )

# 00.resources/signalp6, the manifest's folder, receives the files this script
# leaves in the output folder: the manifest, package/ and .matplotlib/.
resource_root = manifest.parent
resource_root.mkdir(parents=True, exist_ok=True)
# Matplotlib (the plotting library in this environment) keeps a cache folder, by
# default in the home folder, which may be read-only. Point it at
# 00.resources/signalp6/.matplotlib instead; programs started from here inherit
# the setting.
os.environ["MPLCONFIGDIR"] = str(resource_root / ".matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(exist_ok=True)
# A temporary working folder, .signalp6-install-<random letters> inside
# 00.resources/signalp6. The `finally:` at the end removes it, whether the
# installation worked or not.
staging_root = Path(tempfile.mkdtemp(prefix=".signalp6-install-", dir=resource_root))

try:
    # Steps 2 and 3: unpack the archive safely, then find the package and the
    # model for the chosen mode.
    extracted = staging_root / "extracted"
    extracted.mkdir()
    safe_extract(archive, extracted)
    package_root = find_package_root(extracted)
    model_dir = find_model_dir(package_root, mode)

    # Step 4: install SignalP's own package from the unpacked archive, and
    # nothing else. The NumPy and PyTorch versions it works with are already in
    # the environment (envs/signalp6.yaml); nothing comes from PyPI, the public
    # Python package index. sys.executable is the Python running this script,
    # the one of the rule's conda environment, so pip installs into that
    # environment. The pip options: --no-deps installs no other package;
    # --force-reinstall replaces a SignalP left from an earlier installation;
    # --no-build-isolation builds the package with the tools already in the
    # environment instead of downloading build tools; --no-input and
    # --disable-pip-version-check keep pip from asking questions or looking
    # online for a newer pip.
    python = sys.executable
    run(
        [
            python,
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            "--no-input",
            "--no-deps",
            "--force-reinstall",
            "--no-build-isolation",
            str(package_root),
        ]
    )
    # Ask the newly installed package for its version, for the manifest. Its
    # answer is kept, not printed.
    signalp_version = subprocess.check_output(
        [python, "-c", "import signalp; print(signalp.__version__)"], text=True
    ).strip()
    # pip put the signalp6 command into the environment's bin/ folder, which is
    # on PATH while the environment is active. Start it once with --help to show
    # that it runs; the help text goes to the log.
    signalp_executable = shutil.which("signalp6")
    if not signalp_executable:
        fail("the official archive did not install the signalp6 command")
    run([signalp_executable, "--help"])

    # Step 5: pip has already copied the Python package into the conda
    # environment. The unpacked package folder is kept as well, for its model
    # file only: it is moved to 00.resources/signalp6/package, replacing the copy
    # from an earlier installation, and the manifest stores the model folder
    # relative to the manifest (package/models, or package/signalp/model_weights),
    # so the output folder can be moved. That relative path is worked out before
    # the move, while model_dir still points into the temporary folder.
    model_relative = Path("package") / model_dir.relative_to(package_root)
    installed_package = resource_root / "package"
    if installed_package.exists():
        shutil.rmtree(installed_package)
    shutil.move(str(package_root), str(installed_package))

    # The manifest also records "conda_prefix", the folder of the conda
    # environment this installer ran in.
    #
    # Why: the pip command above installed the signalp package inside that conda
    # environment, but this manifest (installation.json) is written in the run's
    # output folder and outlives the environment. Snakemake builds one conda
    # environment per recipe and names its folder after a fingerprint (hash) of
    # the recipe, so the environment is rebuilt from scratch after the recipe
    # changes, after --conda-prefix is pointed somewhere else, or after
    # .snakemake/conda is deleted. The new environment has no signalp package,
    # but the old manifest is still there, so without a check the
    # signalp6_annotation rule would fail hours into a run. That rule therefore
    # checks that SignalP 6 is still installed before it starts; if it is not, it
    # stops at once and tells the user to delete this manifest, so that this
    # installer runs again. conda_prefix says which environment the package was
    # installed into, so that such a failure can be traced.
    #
    # CONDA_PREFIX is set by Conda while the rule's environment is active. If this
    # script is ever started outside Conda it is not set, and sys.prefix (the
    # folder of the Python that is running this script) is recorded instead.
    #
    # Unlike model_dir, which is a path relative to this manifest, conda_prefix
    # is an absolute path on the computer that ran the installer. It is for
    # diagnosis only: nothing uses it to locate files, and once the output
    # folder has been moved to another computer it is just an out-of-date hint.
    conda_prefix = os.environ.get("CONDA_PREFIX", sys.prefix)

    # Step 6: the content of the manifest. schema_version numbers the layout of
    # this file (1). The archive's name and checksum record which archive was
    # installed; model_dir is the model folder, relative to the manifest.
    payload = {
        "schema_version": 1,
        "archive": {"filename": archive.name, "sha256": archive_sha256},
        "mode": mode,
        "signalp_version": signalp_version,
        "model_dir": model_relative.as_posix(),
        "conda_prefix": conda_prefix,
    }
    # The JSON (keys in alphabetical order, 2-space indent) goes to the hidden
    # file .installation.json.tmp next to the manifest, which is then renamed to
    # installation.json in one step. A failed job therefore leaves a complete
    # manifest or none at all, never a half-written one.
    temporary_manifest = manifest.with_name(f".{manifest.name}.tmp")
    with temporary_manifest.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary_manifest, manifest)
finally:
    # Remove the temporary working folder in every case. ignore_errors: a
    # problem while cleaning up must not hide the error that stopped the
    # installation.
    shutil.rmtree(staging_root, ignore_errors=True)
