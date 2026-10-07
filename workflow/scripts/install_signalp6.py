"""Install a user-supplied official SignalP 6 archive for one workflow run.

SignalP 6 is distributed by DTU after licence acceptance, so this script never
downloads, redistributes, or hard-codes a DTU URL.  It safely unpacks the
archive selected in the config, installs only its Python package into this
rule's dedicated Conda environment, and keeps the large model file under the
run's ``00.resources/signalp6`` directory.  The manifest contains relative
paths only; it is therefore not tied to the machine that created the run.
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


def fail(message: str) -> None:
    raise RuntimeError(f"SignalP 6 installation failed: {message}")


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_extract(archive: Path, destination: Path) -> None:
    """Extract ordinary files only, rejecting archive paths that escape dest."""
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
    """Support both documented archive layouts without assuming a local path."""
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


def run(command: list[str]) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


archive = Path(str(snakemake.input.archive)).resolve()
manifest = Path(str(snakemake.output.manifest)).resolve()
params = snakemake.params
mode = str(params.mode).strip().lower()
expected_sha256 = str(params.archive_sha256 or "").strip().lower()

if not archive.is_file():
    fail(f"configured archive does not exist or is not a regular file: {archive}")
if mode not in {"fast", "slow", "slow-sequential"}:
    fail(f"unsupported mode {mode!r}")
if expected_sha256 and (
    len(expected_sha256) != 64
    or any(character not in "0123456789abcdef" for character in expected_sha256)
):
    fail("archive_sha256 must be a 64-character hexadecimal SHA-256 digest")

# Hash the exact local archive before installation. This both verifies an optional
# user-supplied checksum and records enough provenance to reproduce the run.
archive_sha256 = sha256sum(archive)
if expected_sha256 and archive_sha256 != expected_sha256:
    fail(
        "configured archive_sha256 does not match the supplied archive "
        f"(expected {expected_sha256}, observed {archive_sha256})"
    )

resource_root = manifest.parent
resource_root.mkdir(parents=True, exist_ok=True)
# Avoid writing Matplotlib's cache into an arbitrary (possibly read-only) HOME.
os.environ["MPLCONFIGDIR"] = str(resource_root / ".matplotlib")
Path(os.environ["MPLCONFIGDIR"]).mkdir(exist_ok=True)
staging_root = Path(tempfile.mkdtemp(prefix=".signalp6-install-", dir=resource_root))

try:
    extracted = staging_root / "extracted"
    extracted.mkdir()
    safe_extract(archive, extracted)
    package_root = find_package_root(extracted)
    model_dir = find_model_dir(package_root, mode)

    # The compatible NumPy/PyTorch stack is declared in envs/signalp6.yaml.
    # Install only the official local archive package here, never from PyPI.
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
    signalp_version = subprocess.check_output(
        [python, "-c", "import signalp; print(signalp.__version__)"], text=True
    ).strip()
    signalp_executable = shutil.which("signalp6")
    if not signalp_executable:
        fail("the official archive did not install the signalp6 command")
    run([signalp_executable, "--help"])

    # pip has already copied the importable package into the dedicated Conda
    # environment. Keep this separate archive-derived copy only as run-local
    # model data, so the manifest can use a portable relative path.
    model_relative = Path("package") / model_dir.relative_to(package_root)
    installed_package = resource_root / "package"
    if installed_package.exists():
        shutil.rmtree(installed_package)
    shutil.move(str(package_root), str(installed_package))

    payload = {
        "schema_version": 1,
        "archive": {"filename": archive.name, "sha256": archive_sha256},
        "mode": mode,
        "signalp_version": signalp_version,
        "model_dir": model_relative.as_posix(),
    }
    # Write beside the final manifest and rename atomically. A failed job then
    # leaves either a complete manifest or none at all, never misleading JSON.
    temporary_manifest = manifest.with_name(f".{manifest.name}.tmp")
    with temporary_manifest.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary_manifest, manifest)
finally:
    shutil.rmtree(staging_root, ignore_errors=True)
