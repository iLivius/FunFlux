#!/usr/bin/env python3
r"""Work out which AUGUSTUS parameters funannotate predict will use for one sample.

Why this is needed
------------------
AUGUSTUS is the gene finder that weighs most in funannotate's final gene set: its models
with strong evidence count double when EVidenceModeler combines the gene finders. Column 4
of FunFlux's annotation table (files.annotation_params) decides, per sample, which
parameters it runs with:

  * a model name, such as fusarium_graminearum: funannotate copies that species' stored
    AUGUSTUS parameters from the funannotate database (trained_species/<name>/) and uses
    them AS THEY ARE. AUGUSTUS is not trained on this genome;
  * empty or "-": funannotate trains AUGUSTUS on this genome's own conserved genes, the
    BUSCO genes that also train SNAP and GlimmerHMM.

funannotate gives no warning about two cases, and both change the result quietly:

  * a name that is not in the database (a typo, or capitals the folder name does not
    have) makes it train AUGUSTUS under that name instead of stopping;
  * with no name, it first looks for a stored model under the name it gives a newly
    trained one, "<species>_<sample>" in lower case (beauveria_bassiana_arsef3097). Such a
    folder exists only if someone added an earlier result to the database with
    `funannotate species`, as funannotate's own end-of-run message suggests; funannotate
    then uses that stored model as it is and does not train.

This script makes the same decisions funannotate makes (predict.py and
CheckFunannotateSpecies in its library, version 1.8.17), stops on a name the database does
not hold, with the closest names it does hold, and writes down what will happen.

Usage (rule augustus_model_check, workflow/rules/shared/40_annotation.smk):
    python check_augustus_model.py --database DIR --sample NAME --species "Genus species" \
        --model MODEL_OR_EMPTY --optimize yes|no --record FILE

  --database  the funannotate database folder (directories.funannotate_db), already
              unpacked and repaired (rule normalize_funannotate_db)
  --sample    the sample name; funannotate predict gets it as --strain
  --species   column 2 of the annotation table; funannotate predict gets it as -s
  --model     column 4 of the annotation table, "" when it is empty or "-"
  --optimize  yes when parameters.prediction.augustus_optimize is on for this sample
  --record    the file to write (AUGUSTUS_MODEL_RECORD)

Output: the record, three tab-separated lines, which rule funannotate_prediction copies to
the top of its log:
    sample          ARSEF3097
    augustus_mode   pretrained (the stored model, used as it is; AUGUSTUS is not trained on this genome)
    augustus_model  fusarium_graminearum
The same lines are printed (the rule sends them to its log). Exit status 0; 1 when the named
model is not in the database, with the [FunFlux] message on stderr and no record written.

00_common.smk imports this file too, to make the same check while Snakemake reads the
workflow when the database is already on disk; the functions below are written for both.
Uses only the Python standard library. Tests: test_check_augustus_model.py (same folder).
"""

import argparse
import difflib
import os
import re
import sys

# The characters a model name may hold. The 166 trained species of the funannotate
# database use letters, digits, "_" and "-" (for example saccharomyces_cerevisiae_rm11-1a_1);
# "." is allowed too, although no stored name uses it. A "/" or a space can never name a
# folder there.
MODEL_NAME_PATTERN = r"[A-Za-z0-9_.-]+"

# The mode texts of the record, one per case.
MODE_STORED = "pretrained (the stored model, used as it is; AUGUSTUS is not trained on this genome)"
MODE_OWN_STORED = ("pretrained (a model stored earlier under this sample's own name, used as it is; "
                   "AUGUSTUS is not trained again)")
MODE_TRAINED = "trained on this genome's BUSCO genes"
MODE_TRAINED_OPTIMIZED = "trained on this genome's BUSCO genes, with --optimize_augustus"


def known_models(database):
    """The stored AUGUSTUS models of a funannotate database, sorted.

    They are the folder names in <database>/trained_species/, the list funannotate itself
    reads (CheckFunannotateSpecies); like funannotate, names that start with "." (hidden
    folders) are left out. Returns None when the database has no trained_species folder.
    """
    species_folder = os.path.join(database, "trained_species")
    if not os.path.isdir(species_folder):
        return None
    models = []
    for entry in sorted(os.listdir(species_folder)):
        if not entry.startswith(".") and os.path.isdir(os.path.join(species_folder, entry)):
            models.append(entry)
    return models


def own_model_name(species, sample):
    """The name funannotate gives a model it trains for this sample.

    funannotate predict builds it from -s and --strain: "<species>_<strain>", spaces turned
    into "_", in lower case (predict.py, aug_species). FunFlux passes the sample name as
    --strain, so "Beauveria bassiana" and ARSEF3097 give beauveria_bassiana_arsef3097.
    """
    return (species + "_" + sample).replace(" ", "_").lower()


def model_problem(sample, model, database):
    """The [FunFlux] message why the named model cannot be used, or None when it can.

    sample   : the sample name, for the message
    model    : the name from column 4 (not empty)
    database : the funannotate database folder, or None when it is not on disk yet; then
               only the characters of the name are checked
    """
    if not re.fullmatch(MODEL_NAME_PATTERN, model):
        return (
            f"[FunFlux] files.annotation_params: sample '{sample}' names the AUGUSTUS model '{model}', "
            "which is not a possible model name (letters, digits, '_', '-' and '.' only). Leave column 4 "
            "empty, or write '-', to train AUGUSTUS on this genome instead."
        )
    if database is None:
        return None
    models = known_models(database)
    if models is not None and model in models:
        return None
    species_folder = os.path.join(database, "trained_species")
    # difflib.get_close_matches ranks the stored names by how many letters they share with
    # the typed one; up to three that share at least 60 % are suggested. The comparison
    # ignores capitals, so 'Fusarium_graminearum' still finds fusarium_graminearum.
    suggestions = []
    if models:
        lowercase_names = {}
        for name in models:
            lowercase_names.setdefault(name.lower(), name)
        for close in difflib.get_close_matches(model.lower(), list(lowercase_names), n=3, cutoff=0.6):
            suggestions.append(lowercase_names[close])
    if suggestions:
        hint = "Did you mean " + " or ".join(f"'{name}'" for name in suggestions) + "?"
    else:
        hint = f"The names are the folder names in '{species_folder}'."
    return (
        f"[FunFlux] files.annotation_params: sample '{sample}' names the AUGUSTUS model '{model}', "
        f"which is not a trained species of the funannotate database at '{database}'. {hint} "
        "Leave column 4 empty, or write '-', to train AUGUSTUS on this genome instead. (funannotate "
        "itself would not stop: it would silently train AUGUSTUS under the name you typed.)"
    )


def decide(sample, species, model, database, optimize):
    """What funannotate predict will do with AUGUSTUS for this sample.

    Takes the same values as the command line (database must be on disk; optimize is True
    or False). Returns (problem, mode, model_used): problem is the [FunFlux] message when
    the named model cannot be used (mode and model_used are then None); otherwise problem
    is None, mode is one of the MODE_* texts and model_used the name of the stored model,
    or "-" when AUGUSTUS is trained.
    """
    if model:
        problem = model_problem(sample, model, database)
        if problem:
            return problem, None, None
        return None, MODE_STORED, model
    own_name = own_model_name(species, sample)
    models = known_models(database) or []
    if own_name in models:
        return None, MODE_OWN_STORED, own_name
    if optimize:
        return None, MODE_TRAINED_OPTIMIZED, "-"
    return None, MODE_TRAINED, "-"


def build_parser():
    """The command-line options (see the usage at the top of this file)."""
    parser = argparse.ArgumentParser(
        description="Check and record which AUGUSTUS parameters funannotate predict will use for one sample."
    )
    parser.add_argument("--database", required=True, metavar="DIR", help="the funannotate database folder")
    parser.add_argument("--sample", required=True, help="the sample name (funannotate's --strain)")
    parser.add_argument("--species", required=True, help="column 2 of the annotation table (funannotate's -s)")
    parser.add_argument("--model", required=True, help='column 4 of the annotation table; "" when empty')
    parser.add_argument("--optimize", required=True, choices=["yes", "no"],
                        help="yes when the long AUGUSTUS training is on for this sample")
    parser.add_argument("--record", required=True, metavar="FILE", help="the record to write")
    return parser


def main(argv=None):
    """Decide, then write and print the record. Returns the exit status (0 or 1)."""
    args = build_parser().parse_args(argv)
    problem, mode, model_used = decide(args.sample, args.species, args.model.strip(), args.database,
                                       args.optimize == "yes")
    if problem:
        print(problem, file=sys.stderr)
        return 1
    lines = [f"sample\t{args.sample}", f"augustus_mode\t{mode}", f"augustus_model\t{model_used}"]
    with open(args.record, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
