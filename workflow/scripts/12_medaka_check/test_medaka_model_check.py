"""Unit tests for medaka_model_check.py — the pure name-handling logic, plus
main()'s branching with the two Medaka helpers monkeypatched.

Only list_available_models(), resolve_from_reads() and download_model() shell
out to the tool, and every test that reaches one of them swaps in a stub, so
this file runs with no Medaka install, no basecalled reads and no model
download.

Why the branching is worth pinning down: rule check_medaka_model runs the script
on the raw reads, before the assembler starts, and the assembler waits on its
output, so what main() decides here is whether a long-read run stops in seconds
or spends an hour in Flye first and only then discovers the model is unusable.

CHANGED in FunFlux v2: BacFlux runs the script after read filtering, FunFlux on
the raw reads. The last section of this file is new: it runs
resolve_from_reads() itself, with only subprocess.run faked, to pin down the
`consensus` model type and the untouched hand-over of a gzip-compressed file.

Run: pytest workflow/scripts/12_medaka_check/test_medaka_model_check.py
"""

# CHANGED in FunFlux v2: gzip and subprocess are only used by the new tests at
# the end of this file.
import gzip
import subprocess

import medaka_model_check as mc


# ── The model menu a real Medaka install prints ─────────────────────────────
# A small but representative slice of a real `medaka tools list_models`, covering
# the naming shapes the tokenizer must handle: old gNNN names with/without a
# device token, new vX.Y.Z names with pore/speed, and variant/snp models.
AVAILABLE = [
    "r103_fast_g507", "r103_hac_g507", "r103_sup_g507",
    "r103_hac_variant_g507",
    "r941_min_fast_g507", "r941_min_hac_g507", "r941_min_sup_g507",
    "r941_min_hac_variant_g507", "r941_min_hac_snp_g507",
    "r941_prom_fast_g507", "r941_prom_hac_g507", "r941_prom_sup_g507",
    "r1041_e82_400bps_fast_g615", "r1041_e82_400bps_hac_v5.2.0",
    "r1041_e82_400bps_sup_v5.2.0", "r1041_e82_400bps_sup_variant_v5.0.0",
]


# ── Splitting a model name into flowcell, device, accuracy and version ──────
# Medaka model names are underscore-joined but not fixed-width: some carry a
# min/prom device token and some do not, and the version tag is either gNNN or
# vX.Y.Z. tokenize() therefore classifies token by token, and these tests cover
# each shape the real menu above contains. Everything the tokenizer does not
# recognise (pore e82, speed 400bps, network suffixes) stays in the name but is
# not matched on.

def test_tokenize_old_style_with_device():
    t = mc.tokenize("r941_min_hac_g507")
    assert t["flowcell"] == "r941"
    assert t["device"] == "min"
    assert t["accuracy"] == "hac"
    assert t["version"] == "g507"
    assert t["is_variant"] is False


def test_tokenize_no_device_token():
    # r103 names carry no min/prom token; device must come back None, not crash.
    t = mc.tokenize("r103_fast_g507")
    assert t["flowcell"] == "r103"
    assert t["device"] is None
    assert t["accuracy"] == "fast"
    assert t["version"] == "g507"


def test_tokenize_new_style_version_and_pore():
    t = mc.tokenize("r1041_e82_400bps_sup_v5.2.0")
    assert t["flowcell"] == "r1041"
    assert t["accuracy"] == "sup"
    # First version-like token (v5.2.0) wins; pore/speed are ignored for matching.
    assert t["version"] == "v5.2.0"


def test_tokenize_flags_variant_and_snp():
    assert mc.tokenize("r941_min_hac_variant_g507")["is_variant"] is True
    assert mc.tokenize("r941_min_hac_snp_g507")["is_variant"] is True


# ── Narrowing the menu when the configured model is wrong ───────────────────
# The whole value of failing early is lost if the error just says "not available"
# against a list of 80 names. suggest() keeps only the models sharing whatever
# axes the bad name did get right, so the fix is usually the single line above or
# below. variant/snp models are stripped first: they are for variant CALLING, and
# offering one as a polishing model would be a wrong answer, not a near miss.

def test_consensus_models_drops_variant_and_snp():
    kept = mc.consensus_models(AVAILABLE)
    assert "r941_min_hac_g507" in kept
    assert "r941_min_hac_variant_g507" not in kept
    assert "r941_min_hac_snp_g507" not in kept
    assert "r1041_e82_400bps_sup_variant_v5.0.0" not in kept


def test_suggest_narrows_on_a_version_typo():
    # A plausible typo: right flowcell/device/accuracy, wrong version tag.
    text = mc.suggest("r941_min_hac_g999", AVAILABLE)
    assert "r941_min_hac_g507" in text          # the obvious correct pick is shown
    assert "r941_prom_hac_g507" not in text     # different device — excluded
    assert "r103_fast_g507" not in text         # different flowcell — excluded
    assert "variant" not in text                # variant models never suggested


def test_suggest_falls_back_to_full_table_when_unrecognisable():
    text = mc.suggest("totally-bogus-name", AVAILABLE)
    # No recognisable axis — show the whole consensus menu rather than nothing.
    assert "r941_min_hac_g507" in text
    assert "r1041_e82_400bps_sup_v5.2.0" in text


def test_suggest_narrows_on_flowcell_only():
    # Only the flowcell is recognisable — narrow by flowcell alone.
    text = mc.suggest("r1041_somethingwrong", AVAILABLE)
    assert "r1041_e82_400bps_hac_v5.2.0" in text
    assert "r941_min_hac_g507" not in text


# ── main(): an explicit model out of the config ─────────────────────────────
# The user set parameters.<mode>.medaka_model to a name. main() writes the model
# to --out on success and the Medaka polishing step reads it back, so the model is
# resolved once for the run rather than a second time inside the polishing step.
# An invalid name must exit 1 BEFORE the assembler starts, and must write no
# --out file at all: that file IS the model name Medaka will be handed, so there
# is no such thing as a half-valid one.

def test_main_explicit_valid(tmp_path, monkeypatch):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "download_model", lambda model: None)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g507", "--reads", "x.fastq",
                  "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == "r941_min_hac_g507"


def test_main_explicit_local_file_accepted_as_is(tmp_path, monkeypatch):
    # A path to a real file on disk is accepted without list validation
    # (Medaka -m accepts a file). Preserves v1's escape hatch.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    model_file = tmp_path / "custom_model.tar.gz"
    model_file.write_text("x")
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", str(model_file), "--reads", "x.fastq", "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == str(model_file)


def test_main_explicit_invalid_fails_with_suggestion(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g999", "--reads", "x.fastq",
                  "--out", str(out)])
    assert rc == 1
    assert not out.exists()
    assert "r941_min_hac_g507" in capsys.readouterr().err


# ── main(): auto mode, where Medaka answers with a PATH ─────────────────────
# Auto mode means the user left medaka_model empty and Medaka infers the model
# from the basecaller tag in the headers of the raw long-read FASTQ (CHANGED in
# FunFlux v2: BacFlux read the filtlong output here).
#
# The realistic thing `medaka tools resolve_model` prints: a PATH to the model
# file, NOT a bare name — so it is NEVER a member of the (bare-name) AVAILABLE
# list. Auto mode must accept it anyway (it is a valid `-m` argument). This
# guards against a regression: an earlier version required the
# resolved value to be in AVAILABLE, which failed every real auto run.
RESOLVED_PATH = "/opt/medaka/data/r941_min_hac_g507_model_pt.tar.gz"


def test_main_auto_success_accepts_a_path_not_in_the_name_list(tmp_path, monkeypatch):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: RESOLVED_PATH)
    assert RESOLVED_PATH not in AVAILABLE          # guard: it is a path, not a name
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "", "--reads", "x.fastq", "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == RESOLVED_PATH


def test_main_auto_failure_shows_menu(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: None)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "", "--reads", "x.fastq", "--out", str(out)])
    assert rc == 1
    assert not out.exists()
    err = capsys.readouterr().err
    assert "could not infer" in err
    assert "r941_min_hac_g507" in err


# ── main(): the opt-in fallback from a bad explicit model to auto ───────────
# parameters.<mode>.medaka_model_fallback_auto, default off. It only ever loosens
# the failure: an invalid explicit model may be replaced by an auto-inferred one,
# loudly, and if auto cannot help either the run still stops. The warning is part
# of the contract — silently polishing with a model the user did not choose is
# the thing this whole script exists to prevent.

def test_main_middle_option_falls_back_to_auto(tmp_path, monkeypatch, capsys):
    # Invalid explicit model, fallback on, auto succeeds — use auto and warn.
    # Auto returns a PATH (as real medaka does), which must be accepted.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: RESOLVED_PATH)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g999", "--reads", "x.fastq",
                  "--fallback-to-auto", "true", "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == RESOLVED_PATH
    assert "WARNING" in capsys.readouterr().err


def test_main_middle_option_but_auto_also_fails(tmp_path, monkeypatch):
    # Invalid explicit, fallback on, but auto can't help either — still fail.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "resolve_from_reads", lambda reads: None)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g999", "--reads", "x.fastq",
                  "--fallback-to-auto", "true", "--out", str(out)])
    assert rc == 1
    assert not out.exists()


# ── CHANGED in FunFlux v2: what main() asks Medaka, and what it hands over ──
# Every test above replaces the two Medaka helpers, so none of them looks at the
# command line itself. FunFlux made two deliberate changes that live exactly
# there, so the tests below run resolve_from_reads() for real and replace only
# subprocess.run, the one call that would start the actual tool. The stand-in
# records the command it is asked to run and answers the way Medaka does when it
# succeeds: a model path on stdout and exit code 0. Still no Medaka install, no
# real reads and no model download.
#
#   1. The model type is `consensus` (Medaka's generic consensus models), not
#      `consensus_bacteria` as in BacFlux.
#   2. The reads file is the RAW file the user gave and may be gzip-compressed.
#      Medaka reads gzip FASTQ itself (it opens the file with pysam), so the
#      script passes the path on untouched instead of unpacking a copy first.
#
# The last test of this section backs a promise in the module docstring: an
# explicit model is never resolved from the reads; only auto mode does that.

# The basecaller tag exactly as MinKNOW / Dorado write it into every FASTQ header.
BASECALLER_TAG = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"


def write_tiny_gzip_fastq(path):
    """Write a three-read gzip-compressed FASTQ whose headers carry the
    basecaller tag, like the raw reads a FunFlux user hands in. Returns the path
    as a string, which is how main() receives --reads."""
    records = ""
    for number in (1, 2, 3):
        header = (f"@read{number} runid=abc123 read={number} ch=7 "
                  f"basecall_model_version_id={BASECALLER_TAG}")
        records += header + "\nACGTACGTAC\n+\nIIIIIIIIII\n"
    with gzip.open(path, "wt") as handle:
        handle.write(records)
    return str(path)


def test_resolve_from_reads_asks_for_the_generic_consensus_model(monkeypatch):
    commands_run = []

    def fake_run(command, **kwargs):
        commands_run.append(command)
        return subprocess.CompletedProcess(
            command, 0, stdout=RESOLVED_PATH + "\n", stderr="")

    monkeypatch.setattr(mc.subprocess, "run", fake_run)
    resolved = mc.resolve_from_reads("reads.fastq")

    assert resolved == RESOLVED_PATH
    # One call, and the model type is `consensus`.
    assert commands_run == [["medaka", "tools", "resolve_model",
                             "--auto_model", "consensus", "reads.fastq"]]
    # Spelled out because this is the line that differs from BacFlux.
    assert "consensus_bacteria" not in commands_run[0]


def test_resolve_from_reads_passes_a_gzip_file_on_untouched(tmp_path, monkeypatch):
    gzip_reads = write_tiny_gzip_fastq(tmp_path / "reads.fastq.gz")
    # The file we wrote really is a gzip FASTQ carrying the basecaller tag.
    with gzip.open(gzip_reads, "rt") as handle:
        assert BASECALLER_TAG in handle.read()

    handed_over = {}

    def fake_run(command, **kwargs):
        # The last word of the command is the file Medaka would open. Look at it
        # now, while "Medaka" is running.
        path_given_to_medaka = command[-1]
        with open(path_given_to_medaka, "rb") as handle:
            first_two_bytes = handle.read(2)
        handed_over["path"] = path_given_to_medaka
        handed_over["first_two_bytes"] = first_two_bytes
        return subprocess.CompletedProcess(
            command, 0, stdout=RESOLVED_PATH + "\n", stderr="")

    monkeypatch.setattr(mc.subprocess, "run", fake_run)
    resolved = mc.resolve_from_reads(gzip_reads)

    assert resolved == RESOLVED_PATH
    # The same file, not an unpacked temporary copy ...
    assert handed_over["path"] == gzip_reads
    # ... and it is still gzip (every gzip file starts with the bytes 1f 8b).
    assert handed_over["first_two_bytes"] == b"\x1f\x8b"


def test_main_auto_mode_hands_the_raw_gzip_reads_to_medaka(tmp_path, monkeypatch):
    # Whole auto-mode path with a gzip --reads file: main() -> resolve_from_reads()
    # -> the Medaka command. Only the model list and subprocess.run are faked.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    gzip_reads = write_tiny_gzip_fastq(tmp_path / "reads.fastq.gz")
    commands_run = []

    def fake_run(command, **kwargs):
        commands_run.append(command)
        return subprocess.CompletedProcess(
            command, 0, stdout=RESOLVED_PATH + "\n", stderr="")

    monkeypatch.setattr(mc.subprocess, "run", fake_run)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "", "--reads", gzip_reads, "--out", str(out)])

    assert rc == 0
    assert out.read_text().strip() == RESOLVED_PATH
    assert commands_run == [["medaka", "tools", "resolve_model",
                             "--auto_model", "consensus", gzip_reads]]


def test_main_explicit_model_never_asks_medaka_to_resolve_a_model(tmp_path, monkeypatch):
    # An explicit model is compared with the list of names and its weights are
    # stored (download_model, stubbed here); Medaka is never asked to work a
    # model out from the reads.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "download_model", lambda model: None)

    def must_not_be_called(reads):
        raise AssertionError("resolve_from_reads() was called for an explicit model")

    monkeypatch.setattr(mc, "resolve_from_reads", must_not_be_called)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g507", "--reads", "x.fastq",
                  "--out", str(out)])
    assert rc == 0
    assert out.read_text().strip() == "r941_min_hac_g507"


# ── CHANGED in FunFlux v2: a named model's weights are stored by the check ──
# Medaka downloads the weights of a model it does not ship the first time it is
# asked for them. The check runs on the computer that runs Snakemake (a local
# rule), so it fetches them there, and the polishing job on a compute node,
# which may have no internet, finds them stored. The tests below pin down when
# the download is asked for, what is asked, and that a failed download stops the
# check before the assembler starts.

def test_main_explicit_model_asks_medaka_to_store_its_weights(tmp_path, monkeypatch):
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    asked_for = []

    def fake_download(model):
        asked_for.append(model)
        return None

    monkeypatch.setattr(mc, "download_model", fake_download)
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g507", "--reads", "x.fastq", "--out", str(out)])
    assert rc == 0
    assert asked_for == ["r941_min_hac_g507"]
    assert out.read_text().strip() == "r941_min_hac_g507"


def test_main_explicit_model_whose_weights_cannot_be_fetched_stops(tmp_path, monkeypatch, capsys):
    # No internet on the computer that runs the check: stop now, write no --out,
    # and say what to do, with Medaka's own last words.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)
    monkeypatch.setattr(mc, "download_model",
                        lambda model: "RuntimeError: ... could not be downloaded. Check you are connected to the internet")
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", "r941_min_hac_g507", "--reads", "x.fastq", "--out", str(out)])
    assert rc == 1
    assert not out.exists()
    err = capsys.readouterr().err
    assert "could not be downloaded" in err
    assert "r1041_e82_400bps_sup_v5.2.0" in err


def test_main_model_file_on_disk_is_not_downloaded(tmp_path, monkeypatch):
    # A path to a model file needs nothing from the internet.
    monkeypatch.setattr(mc, "list_available_models", lambda: AVAILABLE)

    def must_not_be_called(model):
        raise AssertionError("download_model() was called for a model file on disk")

    monkeypatch.setattr(mc, "download_model", must_not_be_called)
    model_file = tmp_path / "custom_model.tar.gz"
    model_file.write_text("x")
    out = tmp_path / "model.txt"
    rc = mc.main(["--model", str(model_file), "--reads", "x.fastq", "--out", str(out)])
    assert rc == 0


def test_download_model_runs_medaka_download_models(monkeypatch):
    # download_model() itself, with only subprocess.run faked: the command, and
    # None for success.
    commands_run = []

    def fake_run(command, **kwargs):
        commands_run.append(command)
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="Model ... is stored at: ...")

    monkeypatch.setattr(mc.subprocess, "run", fake_run)
    assert mc.download_model("r1041_e82_400bps_sup_v4.3.0") is None
    assert commands_run == [["medaka", "tools", "download_models", "--models", "r1041_e82_400bps_sup_v4.3.0"]]


def test_download_model_returns_the_last_lines_of_a_failure(monkeypatch):
    def fake_run(command, **kwargs):
        stderr = "Traceback (most recent call last):\n  File x\n  line 2\nRuntimeError: could not be downloaded\n"
        return subprocess.CompletedProcess(command, 1, stdout="", stderr=stderr)

    monkeypatch.setattr(mc.subprocess, "run", fake_run)
    problem = mc.download_model("r1041_e82_400bps_sup_v4.3.0")
    assert problem.splitlines()[-1] == "RuntimeError: could not be downloaded"
    assert len(problem.splitlines()) == 3
