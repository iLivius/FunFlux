"""Tests for check_augustus_model.py.

What the script does, in one sentence: it works out, the way funannotate predict does,
whether AUGUSTUS will use a stored model or be trained on the genome for one sample,
stops on a model name the database does not hold, and writes that decision down.

The promises these tests pin down:
  * a stored model is accepted only under its exact folder name (capitals count), and a
    hidden folder (name starting with ".") does not count, as in funannotate;
  * a misspelt name stops with the closest stored names, also when only the capitals
    are wrong;
  * a name with characters no folder name has is refused even before the database exists;
  * an empty column 4 means training (optimised or not), unless a model is stored under
    the sample's own name, which funannotate would then use instead;
  * the command line writes the record and exits 0, or prints the message to stderr,
    writes nothing and exits 1.

Every test builds a tiny fake database in pytest's temporary folder (tmp_path): only the
folder names under trained_species/ matter here.

Run: pytest workflow/scripts/40_annotation/test_check_augustus_model.py
"""

import check_augustus_model as cam


def make_database(tmp_path, species_names):
    # A fake funannotate database: one empty folder per stored model under trained_species/.
    database = tmp_path / "funannotate_db"
    for name in species_names:
        (database / "trained_species" / name).mkdir(parents=True)
    return str(database)


# ── model_problem: is a named model usable? ──

def test_a_stored_model_is_accepted(tmp_path):
    database = make_database(tmp_path, ["fusarium_graminearum", "neurospora_crassa"])
    assert cam.model_problem("SAMPLE1", "fusarium_graminearum", database) is None


def test_a_misspelt_model_is_refused_with_the_closest_name(tmp_path):
    database = make_database(tmp_path, ["fusarium_graminearum", "neurospora_crassa"])
    message = cam.model_problem("SAMPLE1", "fusarium_graminaerum", database)
    assert "sample 'SAMPLE1' names the AUGUSTUS model 'fusarium_graminaerum'" in message
    assert "Did you mean 'fusarium_graminearum'?" in message


def test_capitals_count_but_the_suggestion_ignores_them(tmp_path):
    # funannotate compares folder names exactly, so 'Fusarium_graminearum' is not found;
    # the suggestion still finds the stored name.
    database = make_database(tmp_path, ["fusarium_graminearum"])
    message = cam.model_problem("SAMPLE1", "Fusarium_graminearum", database)
    assert "Did you mean 'fusarium_graminearum'?" in message


def test_a_hidden_folder_is_not_a_model(tmp_path):
    database = make_database(tmp_path, [".hidden_model"])
    message = cam.model_problem("SAMPLE1", ".hidden_model", database)
    assert "is not a trained species" in message


def test_a_database_without_trained_species_names_the_folder(tmp_path):
    database = tmp_path / "empty_db"
    database.mkdir()
    message = cam.model_problem("SAMPLE1", "fusarium_graminearum", str(database))
    assert "The names are the folder names in" in message


def test_without_a_database_only_the_characters_are_checked():
    # Before the database is downloaded, a well-formed name passes (it is checked later)
    # and a path does not.
    assert cam.model_problem("SAMPLE1", "fusarium_graminaerum", None) is None
    message = cam.model_problem("SAMPLE1", "fusarium/graminearum", None)
    assert "is not a possible model name" in message


# ── decide: what funannotate will do ──

def test_an_empty_model_means_training(tmp_path):
    # The stored fusarium_graminearum model is not used although the species matches:
    # with column 4 empty, funannotate looks only for "<species>_<strain>".
    database = make_database(tmp_path, ["fusarium_graminearum"])
    assert cam.decide("SAMPLE1", "Fusarium graminearum", "", database, False) == (None, cam.MODE_TRAINED, "-")
    assert cam.decide("SAMPLE1", "Fusarium graminearum", "", database, True) == (None, cam.MODE_TRAINED_OPTIMIZED, "-")


def test_a_model_stored_under_the_samples_own_name_is_used(tmp_path):
    # funannotate names a model it trains "<species>_<strain>" in lower case and uses a
    # stored folder of that name instead of training (also with the optimisation on).
    database = make_database(tmp_path, ["fusarium_graminearum_sample1"])
    assert cam.own_model_name("Fusarium graminearum", "SAMPLE1") == "fusarium_graminearum_sample1"
    assert cam.decide("SAMPLE1", "Fusarium graminearum", "", database, True) == (
        None, cam.MODE_OWN_STORED, "fusarium_graminearum_sample1")


def test_a_named_model_is_reported_as_stored(tmp_path):
    database = make_database(tmp_path, ["fusarium_graminearum"])
    assert cam.decide("SAMPLE1", "Fusarium graminearum", "fusarium_graminearum", database, True) == (
        None, cam.MODE_STORED, "fusarium_graminearum")


# ── main: the command line used by rule augustus_model_check ──

def test_the_command_line_writes_the_record(tmp_path, capsys):
    database = make_database(tmp_path, ["fusarium_graminearum"])
    record = tmp_path / "augustus_model.txt"
    status = cam.main(["--database", database, "--sample", "SAMPLE1", "--species", "Fusarium graminearum",
                       "--model", "fusarium_graminearum", "--optimize", "no", "--record", str(record)])
    assert status == 0
    assert record.read_text() == (
        f"sample\tSAMPLE1\naugustus_mode\t{cam.MODE_STORED}\naugustus_model\tfusarium_graminearum\n")
    assert "augustus_model\tfusarium_graminearum" in capsys.readouterr().out


def test_the_command_line_stops_on_a_missing_model(tmp_path, capsys):
    database = make_database(tmp_path, ["fusarium_graminearum"])
    record = tmp_path / "augustus_model.txt"
    status = cam.main(["--database", database, "--sample", "SAMPLE1", "--species", "Fusarium graminearum",
                       "--model", "fusarium_graminaerum", "--optimize", "no", "--record", str(record)])
    assert status == 1
    assert not record.exists()
    assert "Did you mean 'fusarium_graminearum'?" in capsys.readouterr().err
