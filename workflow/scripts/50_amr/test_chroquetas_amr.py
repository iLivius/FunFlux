"""Tests for chroquetas_amr.py, the script that screens a genome for antifungal resistance mutations.

How to run (from the repository root):

    python -m pytest workflow/scripts/50_amr -q

Only the Python standard library, pytest and bash are needed for most tests. They read
copies of real ChroQueTas 1.0.1 output, pasted below as text, and run the whole script
with a stand-in for ChroQueTas.sh that writes those copies (see THE STAND-IN).

The tests that run the real ChroQueTas are skipped unless two things are in place:
  * the programs blastp, miniprot, mafft, seqkit and seqret on PATH (the chroquetas
    conda environment, workflow/envs/chroquetas.yaml: put its bin folder on PATH,
    after the one that holds pytest);
  * the environment variable FUNFLUX_CHROQUETAS_TEST_DIR naming a folder that holds
    ChroQueTas.sh and FungAMR_db/ with the species folders, the layout rule
    chroquetas_db makes (workflow/rules/shared/50_amr.smk).
Most of them run on test_data/b0510_beta_tubulin_region.fasta, 3,765 bp of the public
reference genome of Botrytis cinerea B05.10 (RefSeq GCF_000143535.2, chromosome 1,
bases 2,811,131-2,814,895: the beta-tubulin gene BctubA with 1 kb on each side, on the
forward strand of this piece), one on random sequence.
B05.10 carries alanine at position 198 of beta-tubulin, the change E198A that FungAMR
links to carbendazim resistance. One test changes that codon to valine at test time
(E198V, which FungAMR links to resistance to three drugs and to no resistance to
iprodione), to see the verdicts come through the real tool.

THE DATABASE CELLS USED HERE
----------------------------
The verdict tests use drug cells copied from the database of ChroQueTas 1.0.1: the
Botrytis cinerea beta-tubulin changes E198A and E198V; Candida albicans Cyp51 Q21L
and I55A (evidence of no resistance only), A61V (some drugs with both kinds of evidence,
such as Itraconazole(3/-2)) and K128Y (no drug with resistance evidence alone);
Candida albicans Upc2 I142S (one drug, both kinds); and the Cryptococcus deuterogattii
Uxs1 row, which names no drug.

THE REAL OUTPUTS COPIED HERE
---------------------------
The files and logs below are what ChroQueTas 1.0.1 wrote for the hybrid assembly of
B05.10 made by FunFlux: a nuclear run on the nuclear contigs (Beta-tubulin with E198A,
Bos1 and OS1 without mutation, Cytochrome_b not extracted) and a mitochondrial run
with genetic code 4 on the mitochondrial contig (Cytochrome_b extracted, but with
insertions against the reference), and the log of a run on random sequence, where no
target was found.

THE STAND-IN
------------
make_stand_in() writes a small bash script that takes ChroQueTas's options, copies the
prepared files of the run (chosen by the input's name, genome.fasta or mito.fasta)
into the -o folder, prints the prepared log, records its own command line, and exits
with the prepared exit status. The script under test cannot tell it from ChroQueTas,
so every step except the alignments themselves is tested without the tools.
"""

import csv
import os
import random
import shutil
import stat
import subprocess
import sys

import pytest

# The script sits in the same folder as this test file. Putting that folder on the
# list of places where Python looks for modules makes "import chroquetas_amr" work
# wherever pytest is started from.
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TEST_DIR)

import chroquetas_amr  # noqa: E402  (this import must come after the line above)

SCRIPT_PATH = os.path.join(TEST_DIR, "chroquetas_amr.py")
REGION_FASTA = os.path.join(TEST_DIR, "test_data", "b0510_beta_tubulin_region.fasta")
SAMPLE = "SAMPLE1"

# The real ChroQueTas needs its five programs on PATH, and the folder with the script
# and the database named by FUNFLUX_CHROQUETAS_TEST_DIR.
CHROQUETAS_TEST_DIR = os.environ.get("FUNFLUX_CHROQUETAS_TEST_DIR", "")
CHROQUETAS_PROGRAMS = ["blastp", "miniprot", "mafft", "seqkit", "seqret"]


def real_chroquetas_available():
    """True when the five programs are on PATH and the test folder holds ChroQueTas and its database."""
    for program in CHROQUETAS_PROGRAMS:
        if shutil.which(program) is None:
            return False
    if CHROQUETAS_TEST_DIR == "":
        return False
    script = os.path.join(CHROQUETAS_TEST_DIR, "ChroQueTas.sh")
    botrytis = os.path.join(CHROQUETAS_TEST_DIR, "FungAMR_db", "Botrytis_cinerea")
    return os.path.isfile(script) and os.path.isdir(botrytis)


needs_chroquetas = pytest.mark.skipif(
    not real_chroquetas_available(),
    reason="ChroQueTas's programs are not on PATH or FUNFLUX_CHROQUETAS_TEST_DIR is not set",
)

# The 57 species folders of the database of ChroQueTas 1.0.1 (archive db-v20250811.tgz).
DATABASE_FOLDERS = [
    "Arthroderma_vanbreuseghemii", "Aspergillus_awamori", "Aspergillus_fumigatus",
    "Aspergillus_niger", "Aspergillus_terreus", "Aspergillus_tubingensis",
    "Beauveria_bassiana", "Blumeria_graminis", "Blumeria_graminis_hordei",
    "Blumeria_graminis_tritici", "Botrytis_cinerea", "Candida_albicans",
    "Candida_dubliniensis", "Candida_metapsilosis", "Candida_orthopsilosis",
    "Candida_parapsilosis", "Candida_tropicalis", "Candidozyma_auris",
    "Cercospora_beticola", "Clavispora_lusitaniae", "Colletotrichum_acutatum",
    "Colletotrichum_fructicola", "Colletotrichum_gloeosporioides",
    "Colletotrichum_graminicola", "Colletotrichum_siamense", "Corynespora_cassiicola",
    "Cryptococcus_deuterogattii", "Cryptococcus_neoformans", "Erysiphe_necator",
    "Fusarium_fujikuroi", "Fusarium_graminearum", "Fusarium_incarnatum",
    "Fusarium_pseudograminearum", "Histoplasma_capsulatum", "Kluyveromyces_marxianus",
    "Magnaporthe_grisea", "Monilinia_fructicola", "Mycosphaerella_fijiensis",
    "Nakaseomyces_glabratus", "Passalora_fulva", "Phaeosphaeria_nodorum",
    "Phakopsora_pachyrhizi", "Pichia_kudriavzevii", "Pneumocystis_jirovecii",
    "Puccinia_triticina", "Pyrenopeziza_brassicae", "Saccharomyces_cerevisiae",
    "Saccharomyces_paradoxus", "Tapesia_acuformis", "Tapesia_yallundae",
    "Trichophyton_indotineae", "Trichophyton_interdigitale", "Trichophyton_mentagrophytes",
    "Trichophyton_rubrum", "Ustilaginoidea_virens", "Venturia_inaequalis",
    "Zymoseptoria_tritici",
]
BOTRYTIS_TARGETS = ["Beta-tubulin", "Bos1", "Cytochrome_b", "OS1"]

# The position lists (<Target>.txt) of Botrytis cinerea in that database, copied as
# they are. The made-up database of the stand-in tests uses them, so that the
# known_resistance_changes column is filled as in a real run.
POSITION_LIST_HEADER = "Position\tReference\tMutation\tFungicide_resistance\n"
BOTRYTIS_POSITION_LISTS = {
    "Beta-tubulin": POSITION_LIST_HEADER + (
        "198\tE\tA\tCarbendazim(7/NA),Iprodione(8/NA)\n"
        "198\tE\tG\tCarbendazim(8/NA)\n"
        "198\tE\tK\tCarbendazim(7/NA),Diethofencarb(8/NA),Iprodione(8/NA)\n"
        "198\tE\tV\tBenzimidazole(8/NA),Carbendazim(7/NA),Diethofencarb(8/NA),Iprodione(NA/-8)\n"
        "200\tF\tY\tCarbendazim(7/NA),Diethofencarb(8/NA)\n"
        "351\tT\tI\tZoxamide(8/NA)\n"
    ),
    "Bos1": POSITION_LIST_HEADER + (
        "345\tE\tG\tDicarboximide_fungicides(8/NA)\n"
        "346\tV\tA\tDicarboximide_fungicides(8/NA)\n"
        "369\tQ\tP\tProcymidone(8/NA)\n"
        "373\tN\tS\tFludioxonil,_Iprodione(8/NA),Iprodione(8/NA)\n"
        "386\tL\tF\tDicarboximide_fungicides(8/NA)\n"
        "531\tS\tN\tDicarboximide_fungicides(8/NA)\n"
        "538\tG\tR\tDicarboximide_fungicides(8/NA)\n"
        "570\tL\tW\tDicarboximide_fungicides(8/NA)\n"
        "581\tT\tP\tPyrrolnitrin(8/NA)\n"
        "585\tA\tH\tDicarboximide_fungicides(8/NA)\n"
        "1158\tD\tN\tFludioxonil(1/NA)\n"
    ),
    "Cytochrome_b": POSITION_LIST_HEADER + (
        "143\tG\tA\tAzoxystrobin(8/NA),Boscalid(8/NA),Carbendazim(7/NA),Iprodione(8/NA),"
        "Pyraclostrobin(7/NA),Strobilurin(8/NA)\n"
    ),
    "OS1": POSITION_LIST_HEADER + (
        "365\tI\tN\tFludioxonil,_Iprodione(8/NA)\n"
        "365\tI\tS\tCarbendazim(8/NA),Fludioxonil,_Iprodione(8/NA),Iprodione(7/NA),Vinclozolin(7/NA)\n"
        "373\tN\tS\tFludioxonil,_Iprodione(8/NA),Iprodione(8/NA)\n"
    ),
}

# Drug cells of single database rows (THE DATABASE CELLS USED HERE at the top).
CELL_E198A = "Carbendazim(7/NA),Iprodione(8/NA)"
CELL_E198V = "Benzimidazole(8/NA),Carbendazim(7/NA),Diethofencarb(8/NA),Iprodione(NA/-8)"
CELL_CYP51_Q21L = "Fluconazole(NA/-1)"
CELL_CYP51_I55A = (
    "Clotrimazole(NA/-3),Fluconazole(NA/-3),Isavuconazole(NA/-3),Itraconazole(NA/-3),"
    "Posaconazole(NA/-3),Voriconazole(NA/-3)"
)
CELL_CYP51_A61V = (
    "Clotrimazole(3/NA),Fluconazole(NA/-2),Isavuconazole(3/NA),Itraconazole(3/-2),"
    "Posaconazole(2/NA),Voriconazole(3/-2)"
)
CELL_CYP51_K128Y = (
    "Clotrimazole(NA/-3),Fluconazole(3/-8),Isavuconazole(NA/-3),Itraconazole(NA/-3),"
    "Posaconazole(NA/-3),Voriconazole(NA/-3)"
)
CELL_UPC2_I142S = "Fluconazole(8/-7)"
CELL_UXS1_NO_DRUG = ""

# The words of the levels as the cells write them.
WORDS_2 = "the change made in the gene, expressed in another species"
WORDS_3 = "bulk competition assay"
WORDS_7 = "population association study"
WORDS_8 = "natural isolates, not validated further"


# ======================================================================
# THE REAL OUTPUTS, COPIED AS TEXT
# ======================================================================

# The banner ChroQueTas prints at the start of every run.
BANNER = r"""
Put some...

   ____ _                ___            _____
  / ___| |__  _ __ ___  / _ \ _   _  __|_   _|_ _ ___
 | |   | '_ \| '__/ _ \| | | | | | |/ _ \| |/ _` / __|
 | |___| | | | | | (_) | |_| | |_| |  __/| | (_| \__ \
  \____|_| |_|_|  \___/ \__\__\__,_|\___||_|\__,_|___/   v1.0.1


                                                            ...IN YOUR LIFE!

Four steps and you are there!

The species Botrytis_cinerea has 4 proteins associated with AMR
"""

ENDING = """
PLEASE NOTE: the mutations reported by ChroQueTas must be considered according their confidence scores, where different combinations might exist.

See the main repositories for further details:
 -https://github.com/nmquijada/ChroQueTas
 -https://card.mcmaster.ca/fungamrhome
"""

NUCLEAR_LOG = BANNER + """
PLEASE NOTE: as -c/--trans_code has not been specified the genetic code has been automatically set to standard genetic code (1)

Running protein prediction and extraction (step 1/4)
WARNING: We could not extract Cytochrome_b from the genome...
Protein will be discarded for further screening
Done! (step 1/4)
Calculating protein similarity with reference (step 2/4)
Done! (step 2/4)
Performing alignment with reference (step 3/4)
Done! (step 3/4)
Inspecting mutations potentially causing AMR (step 4/4)
Done! (step 4/4)

 -FungAMR mutations found: 1
 -New mutations found: 0

Thanks for using ChroQueTas!

""" + ENDING

MITO_LOG = BANNER + """
Running protein prediction and extraction (step 1/4)
WARNING: We could not extract Beta-tubulin from the genome...
Protein will be discarded for further screening
WARNING: We could not extract Bos1 from the genome...
Protein will be discarded for further screening
WARNING: We could not extract OS1 from the genome...
Protein will be discarded for further screening
Done! (step 1/4)
Calculating protein similarity with reference (step 2/4)
Done! (step 2/4)
Performing alignment with reference (step 3/4)
Done! (step 3/4)
Inspecting mutations potentially causing AMR (step 4/4)
Done! (step 4/4)

No mutations found

Thanks for using ChroQueTas!

""" + ENDING

# A run on 50 kb of random sequence: no target found, exit status 1, and only the tmp/
# folder (with miniprot's index genome.mpi) left in the output folder.
NOTHING_FOUND_LOG = BANNER + """
PLEASE NOTE: as -c/--trans_code has not been specified the genetic code has been automatically set to standard genetic code (1)

Running protein prediction and extraction (step 1/4)
WARNING: We could not extract Beta-tubulin from the genome...
Protein will be discarded for further screening
WARNING: We could not extract Bos1 from the genome...
Protein will be discarded for further screening
WARNING: We could not extract Cytochrome_b from the genome...
Protein will be discarded for further screening
WARNING: We could not extract OS1 from the genome...
Protein will be discarded for further screening

WARNING: ChroQuetas could not extract any query proteins from the genome...
Exiting...
Still... Thanks for using ChroQueTas!

"""

BLASTP_HEADER = (
    "Query\tReference\tperc_identity\tperc_coverage\tevalue\tquery_length\tsubject_length\t"
    "length_alignment\tgaps\tbitscore\tqstart\tqend\tsstart\tsend\n"
)

# The files of the nuclear run that the script looks at. The .faa files only have to
# exist (the script never reads them), so they hold a placeholder sequence.
NUCLEAR_FILES = {
    "genome.ChroQueTaS.AMR_summary.txt": (
        "Protein\tFragment\tPosition_reference\tAA_reference\tAA_query\tFungicide_resistance\n"
        "Beta-tubulin\t1\t198\tE\tA\tCarbendazim(7/NA),Iprodione(8/NA)\n"
    ),
    "genome.ChroQueTaS.AMR_stats.txt": (
        "Protein\tFungAMR_mutations\tNew_mutations\n"
        "Beta-tubulin\t1\t0\n"
        "Bos1\t0\t0\n"
        "OS1\t0\t0\n"
    ),
    "genome.ChroQueTaS.Beta-tubulin.1.tsv": (
        "Position\tReference\tQuery\tResult\tFungicides\n"
        "198\tE\tA\tFungAMR MUTATION\tCarbendazim(7/NA),Iprodione(8/NA)\n"
        "200\tF\tF\tNo mutation\tNA\n"
        "351\tT\tT\tNo mutation\tNA\n"
    ),
    "genome.ChroQueTaS.Bos1.1.tsv": (
        "Position\tReference\tQuery\tResult\tFungicides\n"
        "345\tE\tE\tNo mutation\tNA\n"
        "346\tV\tV\tNo mutation\tNA\n"
        "369\tQ\tQ\tNo mutation\tNA\n"
        "373\tN\tN\tNo mutation\tNA\n"
        "386\tL\tL\tNo mutation\tNA\n"
        "531\tS\tS\tNo mutation\tNA\n"
        "538\tG\tG\tNo mutation\tNA\n"
        "570\tL\tL\tNo mutation\tNA\n"
        "581\tT\tT\tNo mutation\tNA\n"
        "585\tA\tA\tNo mutation\tNA\n"
        "1158\tD\tD\tNo mutation\tNA\n"
    ),
    "genome.ChroQueTaS.OS1.1.tsv": (
        "Position\tReference\tQuery\tResult\tFungicides\n"
        "365\tI\tI\tNo mutation\tNA\n"
        "373\tN\tN\tNo mutation\tNA\n"
    ),
    "genome_Beta-tubulin.1.blastp.txt": BLASTP_HEADER + (
        "genome_Beta-tubulin.1\tBeta-tubulin\t99.776\t100\t0.0\t448\t447\t447\t0\t935\t1\t447\t1\t447\n"
    ),
    "genome_Bos1.1.blastp.txt": BLASTP_HEADER + (
        "genome_Bos1.1\tBos1\t99.924\t100\t0.0\t1316\t1315\t1315\t0\t2698\t1\t1315\t1\t1315\n"
    ),
    "genome_OS1.1.blastp.txt": BLASTP_HEADER + (
        "genome_OS1.1\tOS1\t99.924\t100\t0.0\t1316\t1315\t1315\t0\t2698\t1\t1315\t1\t1315\n"
    ),
    "genome_Beta-tubulin.1.faa": ">genome_Beta-tubulin.1\nMREIVHLQTGQCGNQIGAAFWQ\n",
    "genome_Bos1.1.faa": ">genome_Bos1.1\nMADTAP\n",
    "genome_OS1.1.faa": ">genome_OS1.1\nMADTAP\n",
}

# The files of the mitochondrial run (genetic code 4).
MITO_FILES = {
    "mito.ChroQueTaS.AMR_summary.txt": "No AMR mutations found in mito\n",
    "mito.ChroQueTaS.Cytochrome_b.1.tsv": (
        "WARNING: genome contains insertions for Cytochrome_b and thus AMR are not calculated\n"
    ),
    "mito_Cytochrome_b.1.blastp.txt": BLASTP_HEADER + (
        "mito_Cytochrome_b.1\tCytochrome_b\t96.231\t100\t0.0\t399\t391\t398\t7\t747\t1\t398\t1\t391\n"
    ),
    "mito_Cytochrome_b.1.faa": ">mito_Cytochrome_b.1\nMRIFKSHPLLKLVNSYMIDSPQ\n",
}

# The files left by a run that found nothing: miniprot's index only.
NOTHING_FOUND_FILES = {"tmp/genome.mpi": "index\n"}


# ======================================================================
# HELPERS: FILES, THE MADE-UP DATABASE AND THE STAND-IN
# ======================================================================

def write_files(folder, files):
    """Write {relative path: text} into `folder`, making sub-folders as needed."""
    for relative_path, text in files.items():
        path = os.path.join(folder, relative_path)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as handle:
            handle.write(text)


def write_fasta(path, records):
    """Write (name, sequence) pairs as FASTA, wrapped at 60 bases."""
    with open(path, "w") as fasta:
        for name, sequence in records:
            fasta.write(f">{name}\n")
            for start in range(0, len(sequence), 60):
                fasta.write(sequence[start:start + 60] + "\n")


def random_dna(length, seed):
    """Return `length` random bases; the seed makes the sequence the same on every run."""
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def make_database(folder, species_targets):
    """Write a made-up database: {species folder: [targets]}, one .faa and .txt per target.

    The .txt position lists of Botrytis_cinerea are the real ones (BOTRYTIS_POSITION_LISTS);
    every other species gets the header line alone.
    """
    for species, targets in species_targets.items():
        species_folder = os.path.join(folder, species)
        os.makedirs(species_folder, exist_ok=True)
        for target in targets:
            with open(os.path.join(species_folder, target + ".faa"), "w") as faa:
                faa.write(f">{target}\nMADTAP\n")
            position_list = POSITION_LIST_HEADER
            if species == "Botrytis_cinerea" and target in BOTRYTIS_POSITION_LISTS:
                position_list = BOTRYTIS_POSITION_LISTS[target]
            with open(os.path.join(species_folder, target + ".txt"), "w") as txt:
                txt.write(position_list)
    return folder


def make_stand_in(folder, runs):
    """Write the stand-in for ChroQueTas.sh (THE STAND-IN at the top) and the files it hands out.

    Takes in:  a folder, and {"genome" or "mito": {"files": {...}, "log": text,
               "exit": status}} for the runs the test expects.
    Produces:  the path of the stand-in script. After a run, <folder>/<name>.args holds
               the options the stand-in was given.
    """
    source = os.path.join(folder, "stand_in_runs")
    os.makedirs(source, exist_ok=True)
    for name, run in runs.items():
        write_files(os.path.join(source, name), run.get("files", {}))
        with open(os.path.join(source, name + ".log"), "w") as handle:
            handle.write(run.get("log", ""))
        with open(os.path.join(source, name + ".exit"), "w") as handle:
            handle.write(str(run.get("exit", 0)) + "\n")
    script = os.path.join(folder, "ChroQueTas.sh")
    with open(script, "w") as handle:
        handle.write(
            "#!/bin/bash\n"
            "# A stand-in for ChroQueTas.sh: hands out the prepared files of the run.\n"
            f"source_folder='{source}'\n"
            'all_options="$*"\n'
            'while [ $# -gt 0 ]; do\n'
            '  case "$1" in\n'
            '    -g) genome="$2"; shift 2 ;;\n'
            '    -o) output="$2"; shift 2 ;;\n'
            '    *) shift ;;\n'
            '  esac\n'
            'done\n'
            'name=$(basename "$genome" .fasta)\n'
            'echo "$all_options" > "$source_folder/$name.args"\n'
            'echo "TERM=$TERM" >> "$source_folder/$name.args"\n'
            'mkdir -p "$output/tmp"\n'
            'if [ -d "$source_folder/$name" ]; then cp -r "$source_folder/$name/." "$output/"; fi\n'
            'cat "$source_folder/$name.log"\n'
            'exit $(cat "$source_folder/$name.exit")\n'
        )
    os.chmod(script, os.stat(script).st_mode | stat.S_IXUSR)
    return script


def stand_in_options(folder, name):
    """Read back the options the stand-in was given for one run, or None when it did not run."""
    path = os.path.join(folder, "stand_in_runs", name + ".args")
    if not os.path.isfile(path):
        return None
    with open(path) as handle:
        return handle.read()


def b0510_runs():
    """The two B05.10 runs (THE REAL OUTPUTS COPIED HERE), for the stand-in."""
    return {
        "genome": {"files": NUCLEAR_FILES, "log": NUCLEAR_LOG, "exit": 0},
        "mito": {"files": MITO_FILES, "log": MITO_LOG, "exit": 0},
    }


def full_argv(folder, script, database, **changes):
    """The complete command line of a run in `folder`; `changes` replaces options (None drops one)."""
    options = {
        "--chroquetas-script": script,
        "--database": database,
        "--genome": os.path.join(folder, "contigs_final.fasta"),
        "--mito-contigs": os.path.join(folder, "contigs_mitochondrial.fasta"),
        "--species": "Botrytis cinerea",
        "--mito-genetic-code": "auto",
        "--threads": "2",
        "--sample": SAMPLE,
        "--work-dir": os.path.join(folder, "amr", "chroquetas"),
        "--out-summary": os.path.join(folder, "amr", "amr_summary.tsv"),
        "--out-positions": os.path.join(folder, "amr", "amr_positions.tsv"),
        "--out-mqc": os.path.join(folder, "amr", "amr_summary_mqc.tsv"),
    }
    options.update(changes)
    argv = []
    for option, value in options.items():
        if value is not None:
            argv.extend([option, str(value)])
    return argv


def make_case(folder, runs=None, species_targets=None):
    """Write the genome, the mitochondrial contigs, a database and the stand-in; return (script, database)."""
    write_fasta(os.path.join(folder, "contigs_final.fasta"), [("contig_1", random_dna(300, 1))])
    write_fasta(os.path.join(folder, "contigs_mitochondrial.fasta"), [("mito_1", random_dna(300, 2))])
    if species_targets is None:
        species_targets = {"Botrytis_cinerea": BOTRYTIS_TARGETS, "Fusarium_graminearum": ["Beta-tubulin"]}
    database = make_database(os.path.join(folder, "FungAMR_db"), species_targets)
    script = make_stand_in(folder, b0510_runs() if runs is None else runs)
    return script, database


def read_table(path):
    """Read a tab-separated table back as a list of row dictionaries."""
    with open(path, newline="") as table:
        return list(csv.DictReader(table, delimiter="\t"))


def summary_by_target(path):
    """Read amr_summary.tsv back as {target: row}."""
    rows = {}
    for row in read_table(path):
        rows[row["target"]] = row
    return rows


# ======================================================================
# THE SPECIES: FOLDER, FIRST TWO WORDS, SYNONYMS, NO MATCH
# ======================================================================

@pytest.mark.parametrize(
    "name, folder, how",
    [
        ("Botrytis cinerea", "Botrytis_cinerea", chroquetas_amr.MATCH_NAME),
        ("botrytis_CINEREA", "Botrytis_cinerea", chroquetas_amr.MATCH_NAME),
        ("  Botrytis   cinerea ", "Botrytis_cinerea", chroquetas_amr.MATCH_NAME),
        ("Blumeria graminis f. sp. hordei", "Blumeria_graminis_hordei", chroquetas_amr.MATCH_NAME),
        ("Blumeria graminis f.sp. tritici", "Blumeria_graminis_tritici", chroquetas_amr.MATCH_NAME),
        ("Blumeria graminis F. SP. hordei", "Blumeria_graminis_hordei", chroquetas_amr.MATCH_NAME),
        ("Blumeria graminis", "Blumeria_graminis", chroquetas_amr.MATCH_NAME),
        ("Blumeria graminis f. sp. secalis", "Blumeria_graminis", chroquetas_amr.MATCH_FIRST_TWO_WORDS),
        ("Botrytis cinerea B05.10", "Botrytis_cinerea", chroquetas_amr.MATCH_FIRST_TWO_WORDS),
        ("Candida auris", "Candidozyma_auris", chroquetas_amr.MATCH_SYNONYM),
        ("candida glabrata CBS138", "Nakaseomyces_glabratus", chroquetas_amr.MATCH_SYNONYM),
        ("Cladosporium fulvum", "Passalora_fulva", chroquetas_amr.MATCH_SYNONYM),
        ("Pichia kudriavzevii", "Pichia_kudriavzevii", chroquetas_amr.MATCH_NAME),
        ("Issatchenkia orientalis", "Pichia_kudriavzevii", chroquetas_amr.MATCH_SYNONYM),
    ],
)
def test_species_is_matched_by_name_first_two_words_or_synonym(name, folder, how):
    """Case, underscores, extra spaces and "f. sp." do not matter; a strain name after the species is ignored."""
    assert chroquetas_amr.resolve_species(name, DATABASE_FOLDERS) == (folder, how)


def test_every_synonym_points_at_a_folder_of_the_database():
    """Each synonym leads to a real folder, and none of its names is a folder itself."""
    for current_name, folder in chroquetas_amr.SYNONYMS.items():
        assert folder in DATABASE_FOLDERS, folder
        assert current_name.replace(" ", "_") not in [name.lower() for name in DATABASE_FOLDERS]
        assert chroquetas_amr.resolve_species(current_name, DATABASE_FOLDERS) == (folder, chroquetas_amr.MATCH_SYNONYM)


def test_a_synonym_counts_only_when_the_database_holds_its_folder():
    """A database without Candidozyma_auris does not screen Candida auris."""
    assert chroquetas_amr.resolve_species("Candida auris", ["Botrytis_cinerea"]) == (None, None)


@pytest.mark.parametrize(
    "name, relatives",
    [
        ("Fusarium oxysporum f. sp. lycopersici",
         ["Fusarium_fujikuroi", "Fusarium_graminearum", "Fusarium_incarnatum", "Fusarium_pseudograminearum"]),
        ("Pyricularia oryzae", []),       # not a synonym of Magnaporthe grisea: another species
        ("Botrytis", ["Botrytis_cinerea"]),   # a genus alone is not a species
        ("", []),
    ],
)
def test_no_match_names_the_folders_of_the_same_genus(name, relatives):
    """An unknown species is not matched; the folders of its genus are offered, if any."""
    assert chroquetas_amr.resolve_species(name, DATABASE_FOLDERS) == (None, None)
    assert chroquetas_amr.same_genus_folders(name, DATABASE_FOLDERS) == relatives


def test_not_screened_note_names_the_species_and_the_genus_folders():
    """The note says the species is not among the database's species and lists the relatives."""
    assert chroquetas_amr.not_screened_note("Pyricularia oryzae", DATABASE_FOLDERS) == (
        "the species 'Pyricularia oryzae' is not among the 57 species of the ChroQueTas database"
    )
    assert chroquetas_amr.not_screened_note("Botrytis pseudocinerea", DATABASE_FOLDERS).endswith(
        "; the database holds these of the same genus: Botrytis_cinerea"
    )


def test_match_text_explains_a_match_that_is_not_exact():
    """Nothing for the name as given; the first two words or the synonym are named."""
    assert chroquetas_amr.match_text("Botrytis cinerea", chroquetas_amr.MATCH_NAME) == ""
    assert chroquetas_amr.match_text("Botrytis cinerea B05.10", chroquetas_amr.MATCH_FIRST_TWO_WORDS) == (
        " (from the first two words of 'Botrytis cinerea B05.10')"
    )
    assert chroquetas_amr.match_text("Candida auris", chroquetas_amr.MATCH_SYNONYM) == (
        " (the database name of 'Candida auris', a synonym in NCBI Taxonomy)"
    )


# ======================================================================
# COMPARTMENTS AND GENETIC CODES
# ======================================================================

def test_only_cytochrome_b_is_mitochondrial():
    """Cytochrome_b goes to the mitochondrial run, every other target to the nuclear run."""
    assert chroquetas_amr.compartment_of("Cytochrome_b") == "mitochondrial"
    for target in ["Beta-tubulin", "Cyp51", "Sdh", "Fks", "Bos1", "Atp1"]:
        assert chroquetas_amr.compartment_of(target) == "nuclear"


@pytest.mark.parametrize(
    "folder, setting, code",
    [
        ("Botrytis_cinerea", "auto", 4),
        ("Candida_albicans", "auto", 4),
        ("Saccharomyces_cerevisiae", "auto", 3),
        ("Saccharomyces_paradoxus", "auto", 3),
        ("Nakaseomyces_glabratus", "auto", 3),
        ("Kluyveromyces_marxianus", "auto", 3),
        ("Botrytis_cinerea", "1", 1),
        ("Saccharomyces_cerevisiae", "4", 4),
    ],
)
def test_mitochondrial_code_is_3_for_the_saccharomycetaceae_and_4_otherwise(folder, setting, code):
    """auto gives the yeast mitochondrial code to the four Saccharomycetaceae; a number is used as given."""
    assert chroquetas_amr.mitochondrial_code(folder, setting) == code


def test_command_puts_the_database_before_the_species_and_adds_the_code_when_given():
    """-f must come before -s; -c only for the mitochondrial run."""
    command = chroquetas_amr.chroquetas_command("/x/ChroQueTas.sh", "/t/genome.fasta", "/db", "Botrytis_cinerea", "/o/nuclear", 4, None)
    assert command == ["bash", "/x/ChroQueTas.sh", "-g", "/t/genome.fasta", "-f", "/db", "-s", "Botrytis_cinerea",
                       "-o", "/o/nuclear", "-t", "4"]
    assert command.index("-f") < command.index("-s")
    command = chroquetas_amr.chroquetas_command("/x/ChroQueTas.sh", "/t/mito.fasta", "/db", "Botrytis_cinerea", "/o/mitochondrial", 4, 4)
    assert command[-2:] == ["-c", "4"]


# ======================================================================
# READING CHROQUETAS'S FILES
# ======================================================================

def test_summary_file_with_mutations_and_without(tmp_path):
    """The table of the nuclear run gives one row; the one-line file of the mitochondrial run none."""
    path = tmp_path / "summary.txt"
    path.write_text(NUCLEAR_FILES["genome.ChroQueTaS.AMR_summary.txt"])
    rows = chroquetas_amr.read_summary_file(path)
    assert rows == [{
        "Protein": "Beta-tubulin", "Fragment": "1", "Position_reference": "198", "AA_reference": "E",
        "AA_query": "A", "Fungicide_resistance": "Carbendazim(7/NA),Iprodione(8/NA)",
    }]
    path.write_text(MITO_FILES["mito.ChroQueTaS.AMR_summary.txt"])
    assert chroquetas_amr.read_summary_file(path) == []
    path.write_text("Protein\tPosition\n")
    with pytest.raises(ValueError, match="not the summary file of ChroQueTas 1.0.1"):
        chroquetas_amr.read_summary_file(path)


def test_stats_file_counts_per_target_and_a_missing_file_means_none(tmp_path):
    """Three rows from the nuclear run; no file gives an empty dictionary."""
    path = tmp_path / "stats.txt"
    path.write_text(NUCLEAR_FILES["genome.ChroQueTaS.AMR_stats.txt"])
    assert chroquetas_amr.read_stats_file(path) == {"Beta-tubulin": (1, 0), "Bos1": (0, 0), "OS1": (0, 0)}
    assert chroquetas_amr.read_stats_file(tmp_path / "absent.txt") == {}
    path.write_text("Protein\tFungAMR_mutations\tNew_mutations\nBeta-tubulin\tone\t0\n")
    with pytest.raises(ValueError, match="cannot read the row"):
        chroquetas_amr.read_stats_file(path)


def test_per_target_table_rows_and_calls(tmp_path):
    """The Beta-tubulin table gives three rows with their calls; ChroQueTas's words become the calls of the script."""
    path = tmp_path / "genome.ChroQueTaS.Beta-tubulin.1.tsv"
    path.write_text(NUCLEAR_FILES["genome.ChroQueTaS.Beta-tubulin.1.tsv"])
    kind, rows = chroquetas_amr.read_target_table(path)
    assert kind == "compared"
    assert rows[0] == {
        "position": "198", "reference": "E", "found": "A", "call": "FungAMR mutation",
        "fungicides": "Carbendazim(7/NA),Iprodione(8/NA)", "chroquetas_call": "FungAMR mutation",
    }
    assert [row["call"] for row in rows] == ["FungAMR mutation", "no mutation", "no mutation"]

    # The two rarer words ChroQueTas writes, and a new mutation (its drugs are "Unknown").
    # The last two rows are the ones ChroQueTas 1.0.1 writes for the database rows without
    # a position number: the insertion of Uxs1, and the translocation of Ncp1, for which it
    # read the first residue of the protein. Both are "position not found", without amino acids.
    path.write_text(
        "Position\tReference\tQuery\tResult\tFungicides\n"
        "200\tL\tF\tNew mutation\tUnknown\n"
        "351\tT\t-\tNo alignment\tNA\n"
        "1182ins(fs)\t\t\tPosition not found\tNA\n"
        "ranslocation\tM\tM\tNo mutation\tNA\n"
    )
    kind, rows = chroquetas_amr.read_target_table(path)
    assert [row["call"] for row in rows] == ["new mutation", "no alignment", "position not found", "position not found"]
    assert (rows[3]["position"], rows[3]["reference"], rows[3]["found"]) == ("ranslocation", "", "")
    assert rows[3]["chroquetas_call"] == "no mutation"


def test_per_target_table_corrects_an_unreadable_codon_and_the_reference_amino_acid(tmp_path):
    """X found is not compared; the Hap1 N832N row of the database gives no mutation, not resistance."""
    path = tmp_path / "genome.ChroQueTaS.Hap1.1.tsv"
    path.write_text(
        "Position\tReference\tQuery\tResult\tFungicides\n"
        # Saccharomyces cerevisiae Hap1: the database lists N832N, Clotrimazole(4/NA), so
        # ChroQueTas calls the wild-type asparagine a FungAMR mutation.
        "832\tN\tN\tFungAMR MUTATION\tClotrimazole(4/NA)\n"
        # A codon with N bases, translated X: where FungAMR lists a stop (also written X),
        # ChroQueTas calls it a FungAMR mutation; elsewhere a new mutation.
        "198\tE\tX\tFungAMR MUTATION\tCarbendazim(4/NA)\n"
        "200\tF\tX\tNew mutation\tUnknown\n"
    )
    kind, rows = chroquetas_amr.read_target_table(path)
    assert kind == "compared"
    assert [row["call"] for row in rows] == ["no mutation", "position not found", "position not found"]
    # ChroQueTas's own calls stay, for the comparison with its counts; the X stays too.
    assert [row["chroquetas_call"] for row in rows] == ["FungAMR mutation", "FungAMR mutation", "new mutation"]
    assert [row["found"] for row in rows] == ["N", "X", "X"]
    verdicts = []
    for row in rows:
        verdicts.append(chroquetas_amr.position_verdict(row["call"], row["fungicides"]))
    assert verdicts == ["wild type", "not compared", "not compared"]


def test_per_target_warning_line_and_bad_tables(tmp_path):
    """The insertions warning is recognised; another header, a short row or an unknown Result stop the script."""
    path = tmp_path / "mito.ChroQueTaS.Cytochrome_b.1.tsv"
    path.write_text(MITO_FILES["mito.ChroQueTaS.Cytochrome_b.1.tsv"])
    assert chroquetas_amr.read_target_table(path) == ("warning", [])
    path.write_text("Position\tReference\tQuery\n198\tE\tA\n")
    with pytest.raises(ValueError, match="unexpected first line"):
        chroquetas_amr.read_target_table(path)
    path.write_text("Position\tReference\tQuery\tResult\tFungicides\n198\tE\tA\n")
    with pytest.raises(ValueError, match="3 columns instead of 5"):
        chroquetas_amr.read_target_table(path)
    path.write_text("Position\tReference\tQuery\tResult\tFungicides\n198\tE\tA\tMaybe\tNA\n")
    with pytest.raises(ValueError, match="unknown Result 'Maybe'"):
        chroquetas_amr.read_target_table(path)


def test_targets_not_extracted_are_read_from_the_log():
    """The nuclear run misses Cytochrome_b; the mitochondrial run the three nuclear targets."""
    assert chroquetas_amr.targets_not_extracted(NUCLEAR_LOG) == {"Cytochrome_b"}
    assert chroquetas_amr.targets_not_extracted(MITO_LOG) == {"Beta-tubulin", "Bos1", "OS1"}
    assert chroquetas_amr.targets_not_extracted(NOTHING_FOUND_LOG) == set(BOTRYTIS_TARGETS)


# ======================================================================
# DRUGS, EVIDENCE LEVELS AND VERDICTS
# ======================================================================

@pytest.mark.parametrize(
    "text, drugs",
    [
        ("Carbendazim(7/NA),Iprodione(8/NA)", [("Carbendazim", 7, None), ("Iprodione", 8, None)]),
        ("Benzimidazole(8/NA),Iprodione(NA/-8)", [("Benzimidazole", 8, None), ("Iprodione", None, -8)]),
        ("Fluconazole(2/-3)", [("Fluconazole", 2, -3)]),
        # Names with commas and brackets, as FungAMR writes some of them.
        ("Fludioxonil,_Iprodione(8/NA),Iprodione(8/NA)", [("Fludioxonil,_Iprodione", 8, None), ("Iprodione", 8, None)]),
        ("Azole(s)(7/NA),Azoxystrobin(8/NA)", [("Azole(s)", 7, None), ("Azoxystrobin", 8, None)]),
        ("1,10-Phenanthroline(NA/-8),Amphotericin_B(8/NA)", [("1,10-Phenanthroline", None, -8), ("Amphotericin_B", 8, None)]),
        ("Potassium_cyanide_(KCN)(NA/-8)", [("Potassium_cyanide_(KCN)", None, -8)]),
        ("NA", []),
        ("Unknown", []),
        ("", []),
    ],
)
def test_drug_lists_are_split_at_the_scores(text, drugs):
    """Every drug ends at its (positive/negative) score; NA is None."""
    assert chroquetas_amr.parse_drugs(text) == drugs


def test_unreadable_drug_list_is_a_clear_error():
    """A text without the score in brackets is not FungAMR's: stop."""
    with pytest.raises(ValueError, match="cannot read the drugs"):
        chroquetas_amr.parse_drugs("Carbendazim")


def test_level_words_name_the_kind_of_study_for_both_signs():
    """Eight levels; 8 and -8 get the same words."""
    assert sorted(chroquetas_amr.EVIDENCE_LEVELS) == [1, 2, 3, 4, 5, 6, 7, 8]
    assert chroquetas_amr.EVIDENCE_LEVELS == {
        1: "the change made in the species and measured",
        2: WORDS_2,
        3: WORDS_3,
        4: "experimental evolution",
        5: "gene deletion",
        6: "gene overexpression or duplication",
        7: WORDS_7,
        8: WORDS_8,
    }
    assert chroquetas_amr.evidence_words(8) == chroquetas_amr.evidence_words(-8) == WORDS_8


def test_each_drug_goes_to_the_list_of_its_evidence():
    """Positive only -> resistance; negative only -> no resistance; both -> both lists, each with its level."""
    drugs = chroquetas_amr.parse_drugs(CELL_CYP51_A61V)
    resistance, no_resistance = chroquetas_amr.split_by_evidence(drugs)
    assert resistance == [
        ("Clotrimazole", 3), ("Isavuconazole", 3), ("Itraconazole", 3), ("Posaconazole", 2), ("Voriconazole", 3),
    ]
    assert no_resistance == [("Fluconazole", -2), ("Itraconazole", -2), ("Voriconazole", -2)]


@pytest.mark.parametrize(
    "fungicides, verdict, resistance_to, no_resistance_to",
    [
        # B. cinerea beta-tubulin E198A: resistance only.
        (CELL_E198A, "resistance",
         f"Carbendazim (7: {WORDS_7}); Iprodione (8: {WORDS_8})", "-"),
        # E198V: resistance to three drugs, evidence of no resistance to iprodione.
        (CELL_E198V, "resistance",
         f"Benzimidazole (8: {WORDS_8}); Carbendazim (7: {WORDS_7}); Diethofencarb (8: {WORDS_8})",
         f"Iprodione (-8: {WORDS_8})"),
        # C. albicans Cyp51 Q21L and I55A: evidence of no resistance only.
        (CELL_CYP51_Q21L, "no resistance",
         "-", "Fluconazole (-1: the change made in the species and measured)"),
        (CELL_CYP51_I55A, "no resistance",
         "-",
         f"Clotrimazole (-3: {WORDS_3}); Fluconazole (-3: {WORDS_3}); Isavuconazole (-3: {WORDS_3}); "
         f"Itraconazole (-3: {WORDS_3}); Posaconazole (-3: {WORDS_3}); Voriconazole (-3: {WORDS_3})"),
        # A61V: three drugs with resistance evidence alone make it resistance, although
        # itraconazole and voriconazole also carry evidence of no resistance.
        (CELL_CYP51_A61V, "resistance",
         f"Clotrimazole (3: {WORDS_3}); Isavuconazole (3: {WORDS_3}); Itraconazole (3: {WORDS_3}); "
         f"Posaconazole (2: {WORDS_2}); Voriconazole (3: {WORDS_3})",
         f"Fluconazole (-2: {WORDS_2}); Itraconazole (-2: {WORDS_2}); Voriconazole (-2: {WORDS_2})"),
        # K128Y: the one drug with resistance evidence has evidence against it as well.
        (CELL_CYP51_K128Y, "conflicting",
         f"Fluconazole (3: {WORDS_3})",
         f"Clotrimazole (-3: {WORDS_3}); Fluconazole (-8: {WORDS_8}); Isavuconazole (-3: {WORDS_3}); "
         f"Itraconazole (-3: {WORDS_3}); Posaconazole (-3: {WORDS_3}); Voriconazole (-3: {WORDS_3})"),
        # Upc2 I142S: one drug, both kinds.
        (CELL_UPC2_I142S, "conflicting", f"Fluconazole (8: {WORDS_8})", f"Fluconazole (-7: {WORDS_7})"),
        # Uxs1: no drug named.
        (CELL_UXS1_NO_DRUG, "unknown", "-", "-"),
        ("NA", "unknown", "-", "-"),
    ],
)
def test_verdict_and_drug_cells_of_a_fungamr_change(fungicides, verdict, resistance_to, no_resistance_to):
    """Real database cells: the verdict of the change and its two drug cells."""
    drugs = chroquetas_amr.parse_drugs(fungicides)
    assert chroquetas_amr.mutation_verdict(drugs) == verdict
    assert chroquetas_amr.position_verdict("FungAMR mutation", fungicides) == verdict
    resistance, no_resistance = chroquetas_amr.split_by_evidence(drugs)
    assert chroquetas_amr.evidence_cell(resistance) == resistance_to
    assert chroquetas_amr.evidence_cell(no_resistance) == no_resistance_to


def test_drug_names_with_commas_are_kept_whole_and_separated_by_semicolons():
    """The FungAMR drug "Fludioxonil,_Iprodione" stays whole; the cell separates drugs with "; "."""
    resistance, _ = chroquetas_amr.split_by_evidence(
        chroquetas_amr.parse_drugs("Fludioxonil,_Iprodione(8/NA),Iprodione(8/NA)")
    )
    assert chroquetas_amr.evidence_cell(resistance) == (
        f"Fludioxonil,_Iprodione (8: {WORDS_8}); Iprodione (8: {WORDS_8})"
    )


@pytest.mark.parametrize(
    "call, fungicides, verdict",
    [
        ("new mutation", "Unknown", "unknown"),
        ("no mutation", "NA", "wild type"),
        ("no alignment", "NA", "not compared"),
        ("position not found", "NA", "not compared"),
    ],
)
def test_verdict_of_the_other_calls(call, fungicides, verdict):
    """A change FungAMR does not list is unknown; the reference amino acid is wild type; the rest not compared."""
    assert chroquetas_amr.position_verdict(call, fungicides) == verdict


@pytest.mark.parametrize(
    "position_verdicts, verdict",
    [
        (["wild type", "resistance", "conflicting", "unknown", "no resistance"], "resistance"),
        (["wild type", "no resistance", "unknown", "conflicting"], "conflicting"),
        (["no resistance", "wild type", "unknown"], "unknown"),
        (["wild type", "no resistance", "not compared"], "no resistance"),
        (["wild type", "not compared", "wild type"], "wild type"),
        (["not compared"], "-"),
        ([], "-"),
    ],
)
def test_target_verdict_takes_the_first_that_applies(position_verdicts, verdict):
    """resistance > conflicting > unknown > no resistance > wild type; nothing compared gives -."""
    assert chroquetas_amr.target_verdict(position_verdicts) == verdict


def test_change_log_text_gives_the_verdict_and_the_levels():
    """The console names each FungAMR change with its verdict and the levels of its drugs."""
    assert chroquetas_amr.change_log_text("E198A", "resistance", CELL_E198A) == (
        "E198A (resistance: Carbendazim 7, Iprodione 8)"
    )
    assert chroquetas_amr.change_log_text("E198V", "resistance", CELL_E198V) == (
        "E198V (resistance: Benzimidazole 8, Carbendazim 7, Diethofencarb 8; no resistance: Iprodione -8)"
    )
    assert chroquetas_amr.change_log_text("Q21L", "no resistance", CELL_CYP51_Q21L) == (
        "Q21L (no resistance: Fluconazole -1)"
    )
    assert chroquetas_amr.change_log_text("I142S", "conflicting", CELL_UPC2_I142S) == (
        "I142S (conflicting; resistance: Fluconazole 8; no resistance: Fluconazole -7)"
    )
    assert chroquetas_amr.change_log_text("C1182X", "unknown", CELL_UXS1_NO_DRUG) == "C1182X (unknown)"


# ======================================================================
# THE KNOWN RESISTANCE CHANGES OF A POSITION
# ======================================================================

def test_known_resistance_changes_of_the_botrytis_beta_tubulin(tmp_path):
    """Position 198 -> E198A, E198G, E198K, E198V; 200 -> F200Y; 351 -> T351I; any other position -."""
    path = tmp_path / "Beta-tubulin.txt"
    path.write_text(BOTRYTIS_POSITION_LISTS["Beta-tubulin"])
    known = chroquetas_amr.read_known_resistance_changes(str(path))
    assert known == {"198": ["E198A", "E198G", "E198K", "E198V"], "200": ["F200Y"], "351": ["T351I"]}
    assert chroquetas_amr.known_changes_text(known, "198") == "E198A, E198G, E198K, E198V"
    assert chroquetas_amr.known_changes_text(known, "200") == "F200Y"
    assert chroquetas_amr.known_changes_text(known, "351") == "T351I"
    assert chroquetas_amr.known_changes_text(known, "167") == "-"


def test_known_changes_leave_out_no_resistance_rows_rows_without_a_position_and_no_change(tmp_path):
    """Q21L and I55A carry only evidence of no resistance; insertion, duplication and translocation rows have no position; N832N is no change."""
    path = tmp_path / "Cyp51.txt"
    path.write_text(
        POSITION_LIST_HEADER
        + f"21\tQ\tL\t{CELL_CYP51_Q21L}\n"
        + f"55\tI\tA\t{CELL_CYP51_I55A}\n"
        + f"61\tA\tV\t{CELL_CYP51_A61V}\n"
        + f"128\tK\tY\t{CELL_CYP51_K128Y}\n"
        # The three database rows without a position number (Uxs1, Pdr1, Ncp1), as written there.
        + "1182ins(fs)\tC\t\t\n"
        + "\tc.3308dupT\t\tEchinocandin(8/NA),Fluconazole(8/NA)\n"
        + "ranslocation\tT\t\t5-fluorocytosine(4/NA),Amphotericin_B(4/NA)\n"
        # The Saccharomyces cerevisiae Hap1 row whose new amino acid is the reference one.
        + "832\tN\tN\tClotrimazole(4/NA)\n"
    )
    known = chroquetas_amr.read_known_resistance_changes(str(path))
    # K128Y has a positive level (Fluconazole 3, with -8 as well): it is listed.
    assert known == {"61": ["A61V"], "128": ["K128Y"]}


def test_a_position_list_with_another_header_or_a_short_row_stops(tmp_path):
    """Another first line, a row with fewer columns or a missing file: a clear error."""
    path = tmp_path / "Beta-tubulin.txt"
    path.write_text("Position\tReference\tQuery\tResult\tFungicides\n198\tE\tA\tFungAMR MUTATION\tNA\n")
    with pytest.raises(ValueError, match="not a position list of the ChroQueTas 1.0.1 database"):
        chroquetas_amr.read_known_resistance_changes(str(path))
    path.write_text(POSITION_LIST_HEADER + "198\tE\tA\n")
    with pytest.raises(ValueError, match="3 columns instead of 4"):
        chroquetas_amr.read_known_resistance_changes(str(path))
    with pytest.raises(ValueError, match="does not exist"):
        chroquetas_amr.read_known_resistance_changes(str(tmp_path / "absent.txt"))


# ======================================================================
# WHAT A RUN LEFT: THE COPIES OF EACH TARGET
# ======================================================================

def test_run_results_tell_compared_insertions_below_thresholds_and_skipped(tmp_path):
    """Each kind of copy is recognised from the files ChroQueTas leaves."""
    files = dict(NUCLEAR_FILES)
    # Beta-tubulin found a second time, with insertions (the BLASTP file exists).
    files["genome_Beta-tubulin.2.faa"] = ">genome_Beta-tubulin.2\nMREIV\n"
    files["genome_Beta-tubulin.2.blastp.txt"] = BLASTP_HEADER + "genome_Beta-tubulin.2\tBeta-tubulin\t95\t100\t0.0\t450\t447\t450\t3\t900\t1\t450\t1\t447\n"
    files["genome.ChroQueTaS.Beta-tubulin.2.tsv"] = "WARNING: genome contains insertions for Beta-tubulin and thus AMR are not calculated\n"
    # Bos1 found a second time, below the identity filter (no BLASTP file, the same warning).
    files["genome_Bos1.2.faa"] = ">genome_Bos1.2\nMADTAP\n"
    files["genome.ChroQueTaS.Bos1.2.tsv"] = "WARNING: genome contains insertions for Bos1 and thus AMR are not calculated\n"
    # OS1 extracted but skipped by ChroQueTas: no table at all.
    del files["genome.ChroQueTaS.OS1.1.tsv"]
    # Left-over files of a target the log names as not extracted are not read.
    files["genome_Cytochrome_b.1.faa"] = ">genome_Cytochrome_b.1\n\n"
    write_files(tmp_path, files)

    results = chroquetas_amr.read_run_results(str(tmp_path), "genome", BOTRYTIS_TARGETS, {"Cytochrome_b"})
    kinds = {}
    for target, copies in results.items():
        kinds[target] = [(copy["fragment"], copy["kind"]) for copy in copies]
    assert kinds == {
        "Beta-tubulin": [(1, "compared"), (2, "insertions")],
        "Bos1": [(1, "compared"), (2, "below thresholds")],
        "Cytochrome_b": [],
        "OS1": [(1, "skipped")],
    }
    assert len(results["Bos1"][0]["rows"]) == 11


def test_a_gap_at_a_listed_deletion_counts_only_within_the_copy(tmp_path):
    """A copy that lacks the first 52 amino acids: its gap at 20 is the missing start, not the deletion listed there."""
    # The BLASTP row is the one of a real run on the B05.10 beta-tubulin gene with codons
    # 1 to 53 cut off (reference positions 53 to 447 aligned). The position list behind
    # the table is made up: deletions (written -) listed at 20, outside that stretch, and
    # at 105, inside it, where the copy really lacks the amino acid.
    files = {
        "genome_Beta-tubulin.1.faa": ">genome_Beta-tubulin.1\nMREIV\n",
        "genome_Beta-tubulin.1.blastp.txt": BLASTP_HEADER + (
            "genome_Beta-tubulin.1\tBeta-tubulin\t99.494\t88.3669\t0.0\t396\t447\t395\t0\t823\t1\t395\t53\t447\n"
        ),
        "genome.ChroQueTaS.Beta-tubulin.1.tsv": (
            "Position\tReference\tQuery\tResult\tFungicides\n"
            "20\tF\t-\tFungAMR MUTATION\tCarbendazim(4/NA)\n"
            "105\tH\t-\tFungAMR MUTATION\tCarbendazim(4/NA)\n"
            "198\tE\tA\tFungAMR MUTATION\tCarbendazim(7/NA)\n"
        ),
    }
    write_files(tmp_path, files)
    results = chroquetas_amr.read_run_results(str(tmp_path), "genome", ["Beta-tubulin"], set())
    rows = results["Beta-tubulin"][0]["rows"]
    assert [row["call"] for row in rows] == ["no alignment", "FungAMR mutation", "FungAMR mutation"]
    assert [row["chroquetas_call"] for row in rows] == ["FungAMR mutation"] * 3


def test_aligned_reference_stretch_reads_the_blastp_file_by_column_name(tmp_path):
    """The smallest sstart and the largest send over the rows; a missing file or column is a clear error."""
    path = tmp_path / "genome_Beta-tubulin.1.blastp.txt"
    path.write_text(
        BLASTP_HEADER
        + "genome_Beta-tubulin.1\tBeta-tubulin\t99\t80\t0.0\t400\t447\t360\t0\t700\t1\t360\t60\t420\n"
        + "genome_Beta-tubulin.1\tBeta-tubulin\t99\t78\t0.0\t400\t447\t350\t0\t690\t40\t390\t90\t440\n"
    )
    assert chroquetas_amr.aligned_reference_stretch(str(path)) == (60, 440)
    path.write_text("Query\tReference\tqstart\tqend\nx\ty\t1\t400\n")
    with pytest.raises(ValueError, match="not the BLASTP file of ChroQueTas 1.0.1"):
        chroquetas_amr.aligned_reference_stretch(str(path))
    with pytest.raises(ValueError, match="does not exist"):
        chroquetas_amr.aligned_reference_stretch(str(tmp_path / "absent.blastp.txt"))


def test_fragment_numbers_do_not_mix_targets_whose_names_overlap(tmp_path):
    """Vma1's copies are not Vma16's, and Fcy2's are not Fcy21's."""
    write_files(tmp_path, {
        "genome_Vma1.1.faa": "x", "genome_Vma16.1.faa": "x", "genome_Vma16.2.faa": "x",
        "genome_Fcy21.1.faa": "x",
    })
    assert chroquetas_amr.fragment_numbers(str(tmp_path), "genome", "Vma1") == [1]
    assert chroquetas_amr.fragment_numbers(str(tmp_path), "genome", "Vma16") == [1, 2]
    assert chroquetas_amr.fragment_numbers(str(tmp_path), "genome", "Fcy2") == []


def test_run_judged_by_its_files_not_its_exit_status(tmp_path):
    """The summary file means finished; nothing extracted needs the log and miniprot's index; anything else stops."""
    targets = BOTRYTIS_TARGETS
    finished = tmp_path / "finished"
    write_files(finished, NUCLEAR_FILES)
    assert chroquetas_amr.check_run_finished("nuclear", str(finished), "genome", NUCLEAR_LOG, 1, targets) == "finished"

    nothing = tmp_path / "nothing"
    write_files(nothing, NOTHING_FOUND_FILES)
    assert chroquetas_amr.check_run_finished("nuclear", str(nothing), "genome", NOTHING_FOUND_LOG, 1, targets) == (
        "nothing extracted"
    )

    # The same log without miniprot's index: miniprot did not read the genome (a genetic
    # code it does not know, for example), and the error says so.
    no_index = tmp_path / "no_index"
    no_index.mkdir()
    with pytest.raises(ValueError, match="did not finish") as problem:
        chroquetas_amr.check_run_finished("nuclear", str(no_index), "genome", NOTHING_FOUND_LOG, 1, targets)
    assert "miniprot built no index of the input" in str(problem.value)
    assert "Check the genetic code (-c) and the input FASTA" in str(problem.value)

    # A run that stopped early: the error carries the last lines of its output.
    crashed_log = "Running protein prediction and extraction (step 1/4)\nERROR: miniprot is required and not installed\n"
    with pytest.raises(ValueError) as problem:
        chroquetas_amr.check_run_finished("nuclear", str(no_index), "genome", crashed_log, 1, targets)
    assert "exit status 1" in str(problem.value)
    assert "    ERROR: miniprot is required and not installed" in str(problem.value)


def test_last_lines_keeps_the_end_of_the_output_indented():
    """The last 20 non-empty lines, indented by four spaces."""
    text = "\n".join(f"line {number}" for number in range(1, 31)) + "\n\n"
    block = chroquetas_amr.last_lines(text, 20)
    assert block.splitlines()[0] == "    line 11"
    assert block.splitlines()[-1] == "    line 30"


def test_counts_must_agree_between_the_tables_and_chroquetas_summary(tmp_path):
    """The B05.10 files agree; a summary row more, or a stats row that differs, stops the script."""
    write_files(tmp_path, NUCLEAR_FILES)
    results = chroquetas_amr.read_run_results(str(tmp_path), "genome", BOTRYTIS_TARGETS, {"Cytochrome_b"})
    summary_rows = chroquetas_amr.read_summary_file(tmp_path / "genome.ChroQueTaS.AMR_summary.txt")
    stats = chroquetas_amr.read_stats_file(tmp_path / "genome.ChroQueTaS.AMR_stats.txt")
    chroquetas_amr.check_counts_agree("nuclear", results, summary_rows, stats)

    extra_row = dict(summary_rows[0])
    extra_row["Fungicide_resistance"] = "New mutation"
    with pytest.raises(ValueError, match="per-target tables give"):
        chroquetas_amr.check_counts_agree("nuclear", results, summary_rows + [extra_row], stats)
    with pytest.raises(ValueError, match="stats file gives Bos1"):
        chroquetas_amr.check_counts_agree("nuclear", results, summary_rows, {"Bos1": (1, 0)})


# ======================================================================
# THE STATUS AND THE SUMMARY OF ONE TARGET
# ======================================================================

def compared_copy(fragment, rows):
    """A compared copy as read_run_results() returns it; rows are (position, reference, found, call, fungicides)."""
    row_dictionaries = []
    for position, reference, found, call, fungicides in rows:
        row_dictionaries.append({
            "position": position, "reference": reference, "found": found, "call": call, "fungicides": fungicides,
        })
    return {"fragment": fragment, "kind": "compared", "rows": row_dictionaries}


E198A = ("198", "E", "A", "FungAMR mutation", "Carbendazim(7/NA),Iprodione(8/NA)")
F200F = ("200", "F", "F", "no mutation", "NA")
T351T = ("351", "T", "T", "no mutation", "NA")


@pytest.mark.parametrize(
    "kinds, status, note",
    [
        ([], "not found", chroquetas_amr.NOTE_NOT_FOUND),
        (["below thresholds"], "not found", chroquetas_amr.NOTE_NOT_FOUND),
        (["insertions"], "not evaluated", chroquetas_amr.NOTE_INSERTIONS),
        (["below thresholds", "insertions"], "not evaluated", chroquetas_amr.NOTE_INSERTIONS),
        (["skipped"], "not evaluated", chroquetas_amr.NOTE_SKIPPED),
    ],
)
def test_status_of_a_target_without_a_compared_copy(kinds, status, note):
    """No copy, or copies below the thresholds: not found; insertions or a skipped copy: not evaluated."""
    copies = []
    for number, kind in enumerate(kinds, start=1):
        copies.append({"fragment": number, "kind": kind, "rows": []})
    summary = chroquetas_amr.summarise_target(copies, False)
    assert (summary["status"], summary["note"]) == (status, note)
    assert summary["positions_checked"] == 0
    assert summary["verdict"] == "-"
    assert summary["fungamr_mutations"] == summary["new_mutations"] == "-"
    assert summary["resistance_to"] == summary["no_resistance_to"] == "-"


def test_no_mitochondrial_contig_leaves_cytochrome_b_not_evaluated():
    """The screen found no mitochondrion: not evaluated, no verdict, with the note that says why."""
    summary = chroquetas_amr.summarise_target([], True)
    assert summary["status"] == "not evaluated"
    assert summary["verdict"] == "-"
    assert summary["note"] == chroquetas_amr.NOTE_NO_MITO_CONTIG
    assert summary["note"].startswith("no mitochondrial contig")


def test_evaluated_target_lists_its_verdict_mutations_drugs_and_positions():
    """The B05.10 beta-tubulin: E198A, verdict resistance, carbendazim 7 and iprodione 8, three positions, no note."""
    summary = chroquetas_amr.summarise_target([compared_copy(1, [E198A, F200F, T351T])], False)
    assert summary["status"] == "evaluated"
    assert summary["positions_checked"] == 3
    assert summary["verdict"] == "resistance"
    assert summary["fungamr_mutations"] == "E198A (resistance)"
    assert summary["new_mutations"] == "-"
    assert summary["resistance_to"] == "Carbendazim (E198A, 7); Iprodione (E198A, 8)"
    assert summary["no_resistance_to"] == "-"
    assert summary["note"] == "-"


def test_a_target_without_any_change_is_wild_type():
    """Every compared position carries the reference amino acid."""
    summary = chroquetas_amr.summarise_target([compared_copy(1, [F200F, T351T])], False)
    assert (summary["status"], summary["verdict"]) == ("evaluated", "wild type")
    assert summary["fungamr_mutations"] == summary["resistance_to"] == summary["no_resistance_to"] == "-"


def test_second_copy_mutations_drugs_and_notes():
    """Copy 2 is named; one change on two copies names its drugs once; uncompared parts are noted."""
    new_200 = ("200", "F", "L", "new mutation", "Unknown")
    deleted = ("351", "T", "-", "no alignment", "NA")
    copies = [
        compared_copy(1, [E198A, F200F, T351T]),
        compared_copy(2, [E198A, new_200, deleted]),
        {"fragment": 3, "kind": "insertions", "rows": []},
    ]
    summary = chroquetas_amr.summarise_target(copies, False)
    assert summary["verdict"] == "resistance"
    assert summary["fungamr_mutations"] == "E198A (resistance), E198A (copy 2) (resistance)"
    assert summary["new_mutations"] == "F200L (copy 2)"
    assert summary["resistance_to"] == "Carbendazim (E198A, 7); Iprodione (E198A, 8)"
    assert summary["positions_checked"] == 3
    assert summary["note"] == (
        "found 2 times (copies 1, 2): a duplicated gene, or a region the assembly holds more than once; "
        "copy 3 not compared: insertions against the reference; positions not compared: 351 (no alignment, copy 2)"
    )


def test_a_drug_named_by_two_changes_is_listed_once_per_change():
    """E198V on one copy and F200Y on another: carbendazim twice, with each change; iprodione -8 in no_resistance_to."""
    e198v = ("198", "E", "V", "FungAMR mutation", CELL_E198V)
    f200y = ("200", "F", "Y", "FungAMR mutation", "Carbendazim(7/NA),Diethofencarb(8/NA)")
    summary = chroquetas_amr.summarise_target(
        [compared_copy(1, [e198v, F200F, T351T]), compared_copy(2, [E198A, f200y, T351T])], False,
    )
    assert summary["fungamr_mutations"] == (
        "E198V (resistance), E198A (copy 2) (resistance), F200Y (copy 2) (resistance)"
    )
    assert summary["resistance_to"] == (
        "Benzimidazole (E198V, 8); Carbendazim (E198V, 7); Diethofencarb (E198V, 8); "
        "Carbendazim (E198A, 7); Iprodione (E198A, 8); Carbendazim (F200Y, 7); Diethofencarb (F200Y, 8)"
    )
    assert summary["no_resistance_to"] == "Iprodione (E198V, -8)"
    assert summary["fungamr_log"] == [
        "E198V (resistance: Benzimidazole 8, Carbendazim 7, Diethofencarb 8; no resistance: Iprodione -8)",
        "E198A (copy 2) (resistance: Carbendazim 7, Iprodione 8)",
        "F200Y (copy 2) (resistance: Carbendazim 7, Diethofencarb 8)",
    ]


@pytest.mark.parametrize(
    "rows, verdict, fungamr_mutations, resistance_to, no_resistance_to",
    [
        # A change with evidence of no resistance only, next to wild-type positions.
        ([("21", "Q", "L", "FungAMR mutation", CELL_CYP51_Q21L), F200F], "no resistance",
         "Q21L (no resistance)", "-", "Fluconazole (Q21L, -1)"),
        # A new mutation outranks a no-resistance change.
        ([("21", "Q", "L", "FungAMR mutation", CELL_CYP51_Q21L), ("200", "F", "L", "new mutation", "Unknown")],
         "unknown", "Q21L (no resistance)", "-", "Fluconazole (Q21L, -1)"),
        # A conflicting change outranks a new mutation; its drug sits in both lists.
        ([("142", "I", "S", "FungAMR mutation", CELL_UPC2_I142S), ("200", "F", "L", "new mutation", "Unknown")],
         "conflicting", "I142S (conflicting)", "Fluconazole (I142S, 8)", "Fluconazole (I142S, -7)"),
        # A FungAMR change without drugs is unknown.
        ([("1182", "C", "Y", "FungAMR mutation", CELL_UXS1_NO_DRUG), T351T], "unknown", "C1182Y (unknown)", "-", "-"),
    ],
)
def test_target_verdict_from_its_changes(rows, verdict, fungamr_mutations, resistance_to, no_resistance_to):
    """The target takes the first verdict that applies among its positions; each change shows its own."""
    summary = chroquetas_amr.summarise_target([compared_copy(1, rows)], False)
    assert summary["verdict"] == verdict
    assert summary["fungamr_mutations"] == fungamr_mutations
    assert summary["resistance_to"] == resistance_to
    assert summary["no_resistance_to"] == no_resistance_to


def test_a_copy_without_any_compared_position_is_not_evaluated():
    """Every position in a gap (a gene cut at a contig end, a premature stop): not evaluated, no verdict, positions named."""
    gaps = [("420", "S", "-", "no alignment", "NA"), ("440", "E", "-", "no alignment", "NA")]
    summary = chroquetas_amr.summarise_target([compared_copy(1, gaps)], False)
    assert summary["status"] == "not evaluated"
    assert summary["verdict"] == "-"
    assert summary["positions_checked"] == 0
    assert summary["note"] == (
        chroquetas_amr.NOTE_NOTHING_COMPARED + "; positions not compared: 420 (no alignment), 440 (no alignment)"
    )


# ======================================================================
# THE OUTPUT TABLES AND THE LOG LINE
# ======================================================================

def test_the_three_tables_have_exactly_these_columns():
    """The column names and their order are part of the output format."""
    assert chroquetas_amr.POSITIONS_COLUMNS == [
        "sample", "target", "compartment", "fragment", "position", "reference_aa", "found_aa", "call",
        "verdict", "resistance_to", "no_resistance_to", "known_resistance_changes", "fungamr_drugs",
    ]
    assert chroquetas_amr.SUMMARY_COLUMNS == [
        "sample", "species", "target", "compartment", "status", "positions_checked", "verdict",
        "fungamr_mutations", "new_mutations", "resistance_to", "no_resistance_to", "note",
    ]


def test_position_rows_carry_verdict_drug_lists_known_changes_and_drugs_unchanged():
    """One row per position: the verdict, the two drug lists for FungAMR mutations, the known changes, ChroQueTas's text."""
    known = {"198": ["E198A", "E198G", "E198K", "E198V"], "200": ["F200Y"], "351": ["T351I"]}
    e198v = ("198", "E", "V", "FungAMR mutation", CELL_E198V)
    new_200 = ("200", "F", "L", "new mutation", "Unknown")
    deleted = ("351", "T", "-", "no alignment", "NA")
    rows = chroquetas_amr.position_rows(SAMPLE, "Beta-tubulin", "nuclear", [
        compared_copy(1, [E198A, F200F, T351T]),
        compared_copy(2, [e198v, new_200, deleted]),
        {"fragment": 3, "kind": "insertions", "rows": []},
    ], known)
    assert rows == [
        [SAMPLE, "Beta-tubulin", "nuclear", "1", "198", "E", "A", "FungAMR mutation", "resistance",
         f"Carbendazim (7: {WORDS_7}); Iprodione (8: {WORDS_8})", "-", "E198A, E198G, E198K, E198V", CELL_E198A],
        [SAMPLE, "Beta-tubulin", "nuclear", "1", "200", "F", "F", "no mutation", "wild type", "-", "-", "F200Y", "-"],
        [SAMPLE, "Beta-tubulin", "nuclear", "1", "351", "T", "T", "no mutation", "wild type", "-", "-", "T351I", "-"],
        [SAMPLE, "Beta-tubulin", "nuclear", "2", "198", "E", "V", "FungAMR mutation", "resistance",
         f"Benzimidazole (8: {WORDS_8}); Carbendazim (7: {WORDS_7}); Diethofencarb (8: {WORDS_8})",
         f"Iprodione (-8: {WORDS_8})", "E198A, E198G, E198K, E198V", CELL_E198V],
        [SAMPLE, "Beta-tubulin", "nuclear", "2", "200", "F", "L", "new mutation", "unknown", "-", "-", "F200Y", "-"],
        # A position that was not compared gets no known changes.
        [SAMPLE, "Beta-tubulin", "nuclear", "2", "351", "T", "-", "no alignment", "not compared", "-", "-", "-", "-"],
    ]


def test_mqc_table_has_the_custom_content_header_and_named_rows(tmp_path):
    """Five header lines, then the columns with a first Sample column, rows named "<sample>: <target>"."""
    path = tmp_path / "amr_summary_mqc.tsv"
    row = [SAMPLE, "Botrytis_cinerea", "Bos1", "nuclear", "evaluated", "11", "wild type", "-", "-", "-", "-", "-"]
    chroquetas_amr.write_mqc_table(path, [row], SAMPLE)
    lines = path.read_text().splitlines()
    assert lines[0] == "# id: amr_chroquetas"
    assert lines[1] == "# section_name: 'AMR mutations (ChroQueTas)'"
    assert lines[2].startswith("# description: '") and "FungAMR" in lines[2] and lines[2].endswith("'")
    description = lines[2][len("# description: '"):-1]
    assert "'" not in description
    # The description says what each verdict means.
    for verdict in ("resistance (", "conflicting (", "unknown (", "no resistance (", "wild type ("):
        assert verdict in description
    assert lines[3] == "# plot_type: 'table'"
    assert lines[4].startswith("# pconfig:")
    assert lines[5] == "Sample\t" + "\t".join(chroquetas_amr.SUMMARY_COLUMNS)
    assert lines[6] == f"{SAMPLE}: Bos1\t" + "\t".join(row)


def test_summary_line_names_the_folder_the_statuses_the_verdicts_and_the_mutations():
    """The console line counts the targets per status and per verdict, and gives each FungAMR change with its verdict."""
    e198v = ("198", "E", "V", "FungAMR mutation", CELL_E198V)
    q21l = ("21", "Q", "L", "FungAMR mutation", CELL_CYP51_Q21L)
    summaries = {
        "Beta-tubulin": chroquetas_amr.summarise_target([compared_copy(1, [E198A, F200F, T351T])], False),
        "Bos1": chroquetas_amr.summarise_target([compared_copy(1, [F200F])], False),
        "Cytochrome_b": chroquetas_amr.summarise_target([{"fragment": 1, "kind": "insertions", "rows": []}], False),
        "OS1": chroquetas_amr.summarise_target([], False),
    }
    line = chroquetas_amr.summary_line(SAMPLE, "Botrytis_cinerea", "", BOTRYTIS_TARGETS, summaries)
    assert line == (
        f"[FunFlux] amr_chroquetas {SAMPLE}: species folder Botrytis_cinerea; 4 targets: 2 evaluated "
        "(Beta-tubulin, Bos1), 1 not evaluated (Cytochrome_b), 1 not found (OS1); verdicts: 1 resistance "
        "(Beta-tubulin), 1 wild type (Bos1); FungAMR mutations: Beta-tubulin E198A (resistance: Carbendazim 7, "
        "Iprodione 8); new mutations: none"
    )

    summaries = {
        "Beta-tubulin": chroquetas_amr.summarise_target([compared_copy(1, [e198v])], False),
        "Cyp51": chroquetas_amr.summarise_target([compared_copy(1, [q21l])], False),
    }
    line = chroquetas_amr.summary_line(SAMPLE, "Fakeus_species", "", ["Beta-tubulin", "Cyp51"], summaries)
    assert line.endswith(
        "verdicts: 1 resistance (Beta-tubulin), 1 no resistance (Cyp51); FungAMR mutations: Beta-tubulin E198V "
        "(resistance: Benzimidazole 8, Carbendazim 7, Diethofencarb 8; no resistance: Iprodione -8); "
        "Cyp51 Q21L (no resistance: Fluconazole -1); new mutations: none"
    )

    nothing = {"OS1": chroquetas_amr.summarise_target([], False)}
    line = chroquetas_amr.summary_line(SAMPLE, "Botrytis_cinerea", "", ["OS1"], nothing)
    assert "1 target: 0 evaluated, 0 not evaluated, 1 not found (OS1); verdicts: none;" in line


# ======================================================================
# THE WHOLE SCRIPT WITH THE STAND-IN (NO TOOLS NEEDED)
# ======================================================================

# The three tables of the B05.10 runs, as the script must write them.
EXPECTED_SUMMARY = [
    "sample\tspecies\ttarget\tcompartment\tstatus\tpositions_checked\tverdict\tfungamr_mutations\tnew_mutations\t"
    "resistance_to\tno_resistance_to\tnote",
    f"{SAMPLE}\tBotrytis_cinerea\tBeta-tubulin\tnuclear\tevaluated\t3\tresistance\tE198A (resistance)\t-\t"
    "Carbendazim (E198A, 7); Iprodione (E198A, 8)\t-\t-",
    f"{SAMPLE}\tBotrytis_cinerea\tBos1\tnuclear\tevaluated\t11\twild type\t-\t-\t-\t-\t-",
    f"{SAMPLE}\tBotrytis_cinerea\tCytochrome_b\tmitochondrial\tnot evaluated\t0\t-\t-\t-\t-\t-\t"
    + chroquetas_amr.NOTE_INSERTIONS,
    f"{SAMPLE}\tBotrytis_cinerea\tOS1\tnuclear\tevaluated\t2\twild type\t-\t-\t-\t-\t-",
]
EXPECTED_POSITIONS_START = [
    "sample\ttarget\tcompartment\tfragment\tposition\treference_aa\tfound_aa\tcall\tverdict\tresistance_to\t"
    "no_resistance_to\tknown_resistance_changes\tfungamr_drugs",
    f"{SAMPLE}\tBeta-tubulin\tnuclear\t1\t198\tE\tA\tFungAMR mutation\tresistance\t"
    f"Carbendazim (7: {WORDS_7}); Iprodione (8: {WORDS_8})\t-\tE198A, E198G, E198K, E198V\t{CELL_E198A}",
    f"{SAMPLE}\tBeta-tubulin\tnuclear\t1\t200\tF\tF\tno mutation\twild type\t-\t-\tF200Y\t-",
    f"{SAMPLE}\tBeta-tubulin\tnuclear\t1\t351\tT\tT\tno mutation\twild type\t-\t-\tT351I\t-",
    f"{SAMPLE}\tBos1\tnuclear\t1\t345\tE\tE\tno mutation\twild type\t-\t-\tE345G\t-",
]
# The summary line of the B05.10 runs.
EXPECTED_CONSOLE_LINE = (
    f"[FunFlux] amr_chroquetas {SAMPLE}: species folder Botrytis_cinerea; 4 targets: 3 evaluated "
    "(Beta-tubulin, Bos1, OS1), 1 not evaluated (Cytochrome_b), 0 not found; verdicts: 1 resistance "
    "(Beta-tubulin), 2 wild type (Bos1, OS1); FungAMR mutations: Beta-tubulin E198A (resistance: "
    "Carbendazim 7, Iprodione 8); new mutations: none"
)


def test_the_b0510_runs_give_the_expected_tables(tmp_path, capsys):
    """Nuclear and mitochondrial runs, the compartment rule, the three tables, the console line."""
    script, database = make_case(tmp_path)
    chroquetas_amr.main(full_argv(tmp_path, script, database))

    amr = tmp_path / "amr"
    assert (amr / "amr_summary.tsv").read_text().splitlines() == EXPECTED_SUMMARY
    positions = (amr / "amr_positions.tsv").read_text().splitlines()
    assert positions[:5] == EXPECTED_POSITIONS_START
    assert len(positions) == 1 + 3 + 11 + 2
    mqc = (amr / "amr_summary_mqc.tsv").read_text().splitlines()
    assert mqc[6] == f"{SAMPLE}: Beta-tubulin\t" + EXPECTED_SUMMARY[1]

    # ChroQueTas's own folders, one per run; the mitochondrial run used -c 4, the
    # nuclear run no -c; both ran with TERM=dumb and the database before the species.
    assert sorted(os.listdir(amr / "chroquetas")) == ["mitochondrial", "nuclear"]
    assert (amr / "chroquetas" / "nuclear" / "genome.ChroQueTaS.Beta-tubulin.1.tsv").exists()
    nuclear_options = stand_in_options(tmp_path, "genome")
    mito_options = stand_in_options(tmp_path, "mito")
    assert "-c" not in nuclear_options.split()
    assert mito_options.split("\n")[0].endswith("-c 4")
    assert "-f " + str(database) + " -s Botrytis_cinerea" in nuclear_options
    assert "TERM=dumb" in nuclear_options

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err.count("[FunFlux]") == 1
    assert EXPECTED_CONSOLE_LINE in captured.err.splitlines()
    assert "----- ChroQueTas output, nuclear run -----" in captured.err


def test_a_target_found_by_the_other_run_is_ignored(tmp_path):
    """Cytochrome_b in the nuclear run (a nuclear copy of mitochondrial DNA) does not count; the mitochondrial run decides."""
    nuclear_files = dict(NUCLEAR_FILES)
    nuclear_files["genome_Cytochrome_b.1.faa"] = ">genome_Cytochrome_b.1\nMRIF\n"
    nuclear_files["genome_Cytochrome_b.1.blastp.txt"] = BLASTP_HEADER + "x\tCytochrome_b\t99\t100\t0.0\t391\t391\t391\t0\t700\t1\t391\t1\t391\n"
    nuclear_files["genome.ChroQueTaS.Cytochrome_b.1.tsv"] = "Position\tReference\tQuery\tResult\tFungicides\n143\tG\tG\tNo mutation\tNA\n"
    nuclear_log = NUCLEAR_LOG.replace("WARNING: We could not extract Cytochrome_b from the genome...\n", "")
    runs = b0510_runs()
    runs["genome"] = {"files": nuclear_files, "log": nuclear_log, "exit": 0}
    script, database = make_case(tmp_path, runs)
    chroquetas_amr.main(full_argv(tmp_path, script, database))
    rows = summary_by_target(tmp_path / "amr" / "amr_summary.tsv")
    assert rows["Cytochrome_b"]["status"] == "not evaluated"
    assert rows["Cytochrome_b"]["note"] == chroquetas_amr.NOTE_INSERTIONS
    for row in read_table(tmp_path / "amr" / "amr_positions.tsv"):
        assert row["target"] != "Cytochrome_b"


def test_without_the_mitochondrion_screen_the_mitochondrial_run_reads_the_genome(tmp_path):
    """--mito-contigs given empty (or with no value): both runs read the genome, as genome.fasta and mito.fasta."""
    script, database = make_case(tmp_path)
    argv = full_argv(tmp_path, script, database, **{"--mito-contigs": None}) + ["--mito-contigs"]
    chroquetas_amr.main(argv)
    assert stand_in_options(tmp_path, "mito") is not None
    assert summary_by_target(tmp_path / "amr" / "amr_summary.tsv")["Beta-tubulin"]["status"] == "evaluated"


def test_an_empty_mitochondrial_file_means_no_mitochondrial_run(tmp_path):
    """The screen found no mitochondrion: one run only, Cytochrome_b not evaluated."""
    script, database = make_case(tmp_path)
    open(tmp_path / "contigs_mitochondrial.fasta", "w").close()
    chroquetas_amr.main(full_argv(tmp_path, script, database))
    assert stand_in_options(tmp_path, "mito") is None
    assert os.listdir(tmp_path / "amr" / "chroquetas") == ["nuclear"]
    row = summary_by_target(tmp_path / "amr" / "amr_summary.tsv")["Cytochrome_b"]
    assert (row["status"], row["note"]) == ("not evaluated", chroquetas_amr.NOTE_NO_MITO_CONTIG)


def test_a_mitochondrial_code_given_by_hand_is_passed_on(tmp_path):
    """--mito-genetic-code 1 reaches ChroQueTas as -c 1."""
    script, database = make_case(tmp_path)
    chroquetas_amr.main(full_argv(tmp_path, script, database, **{"--mito-genetic-code": "1"}))
    assert stand_in_options(tmp_path, "mito").split("\n")[0].endswith("-c 1")


def test_a_species_with_only_cytochrome_b_has_no_nuclear_run(tmp_path):
    """Colletotrichum graminicola has only Cytochrome_b in the database: the nuclear run is not needed."""
    script, database = make_case(tmp_path, species_targets={"Colletotrichum_graminicola": ["Cytochrome_b"]})
    chroquetas_amr.main(full_argv(tmp_path, script, database, **{"--species": "Colletotrichum graminicola"}))
    assert stand_in_options(tmp_path, "genome") is None
    assert os.listdir(tmp_path / "amr" / "chroquetas") == ["mitochondrial"]
    assert list(summary_by_target(tmp_path / "amr" / "amr_summary.tsv")) == ["Cytochrome_b"]


def test_nothing_found_in_either_run_gives_four_targets_not_found(tmp_path):
    """Both runs end with "could not extract any query proteins": every target not found, no failure."""
    nothing = {"files": NOTHING_FOUND_FILES, "log": NOTHING_FOUND_LOG, "exit": 1}
    mito_nothing = {"files": {"tmp/mito.mpi": "index\n"}, "log": NOTHING_FOUND_LOG, "exit": 1}
    script, database = make_case(tmp_path, {"genome": nothing, "mito": mito_nothing})
    chroquetas_amr.main(full_argv(tmp_path, script, database))
    rows = summary_by_target(tmp_path / "amr" / "amr_summary.tsv")
    for target in BOTRYTIS_TARGETS:
        assert (rows[target]["status"], rows[target]["note"]) == ("not found", chroquetas_amr.NOTE_NOT_FOUND)
    assert len(read_table(tmp_path / "amr" / "amr_positions.tsv")) == 0


def test_a_species_not_in_the_database_is_not_screened(tmp_path, capsys):
    """One summary row with target -, an empty positions table, one NOTE line, ChroQueTas not run."""
    script, database = make_case(tmp_path)
    chroquetas_amr.main(full_argv(tmp_path, script, database, **{"--species": "Fusarium oxysporum"}))
    lines = (tmp_path / "amr" / "amr_summary.tsv").read_text().splitlines()
    assert lines[0] == "\t".join(chroquetas_amr.SUMMARY_COLUMNS)
    assert lines[1] == (
        f"{SAMPLE}\tFusarium oxysporum\t-\t-\tnot screened\t0\t-\t-\t-\t-\t-\tthe species 'Fusarium oxysporum' is "
        "not among the 2 species of the ChroQueTas database; the database holds these of the same genus: "
        "Fusarium_graminearum"
    )
    assert (tmp_path / "amr" / "amr_positions.tsv").read_text() == "\t".join(chroquetas_amr.POSITIONS_COLUMNS) + "\n"
    mqc = (tmp_path / "amr" / "amr_summary_mqc.tsv").read_text().splitlines()
    assert mqc[6].startswith(f"{SAMPLE}: -\t")
    assert os.listdir(tmp_path / "amr" / "chroquetas") == []
    assert stand_in_options(tmp_path, "genome") is None
    err = capsys.readouterr().err
    assert err.startswith(f"[FunFlux] NOTE amr_chroquetas {SAMPLE}: not screened: the species 'Fusarium oxysporum'")
    assert err.count("\n") == 1


# ======================================================================
# BAD INPUT AND FAILED RUNS: A CLEAR ERROR AND NO TABLE
# ======================================================================

def assert_stopped_with_message(tmp_path, argv, capsys, expected_text):
    """Run main(argv); it must exit with code 1, print the [FunFlux] error line and write no table."""
    with pytest.raises(SystemExit) as stop:
        chroquetas_amr.main(argv)
    assert stop.value.code == 1
    err = capsys.readouterr().err
    error_lines = [line for line in err.splitlines() if line.startswith("[FunFlux] amr_chroquetas: error:")]
    assert len(error_lines) == 1
    assert expected_text in err
    for name in ("amr_summary.tsv", "amr_positions.tsv", "amr_summary_mqc.tsv"):
        assert not (tmp_path / "amr" / name).exists()


def test_a_run_that_crashes_stops_the_script_with_its_last_lines(tmp_path, capsys):
    """A run without the summary file and without the "nothing extracted" message is a failure."""
    crashed = {"files": {}, "log": BANNER + "\nERROR: miniprot is required and not installed\n", "exit": 1}
    script, database = make_case(tmp_path, {"genome": crashed, "mito": b0510_runs()["mito"]})
    assert_stopped_with_message(tmp_path, full_argv(tmp_path, script, database), capsys,
                                "    ERROR: miniprot is required and not installed")


def test_a_run_with_exit_status_0_but_no_summary_file_is_a_failure(tmp_path, capsys):
    """The exit status is not trusted: no summary file, no success."""
    silent = {"files": {}, "log": "Done! (step 4/4)\n", "exit": 0}
    script, database = make_case(tmp_path, {"genome": silent, "mito": b0510_runs()["mito"]})
    assert_stopped_with_message(tmp_path, full_argv(tmp_path, script, database), capsys,
                                "the nuclear run of ChroQueTas did not finish (exit status 0)")


def test_an_existing_run_folder_is_refused(tmp_path, capsys):
    """An old nuclear/ folder must not be read as a new result."""
    script, database = make_case(tmp_path)
    os.makedirs(tmp_path / "amr" / "chroquetas" / "nuclear")
    assert_stopped_with_message(tmp_path, full_argv(tmp_path, script, database), capsys, "already exists")


@pytest.mark.parametrize(
    "option, value, expected_text",
    [
        ("--genome", "nowhere.fasta", "does not exist"),
        ("--mito-contigs", "nowhere_mito.fasta", "the mitochondrial contigs file"),
        ("--chroquetas-script", "nowhere.sh", "the ChroQueTas script"),
        ("--mito-genetic-code", "34", "NCBI genetic code: 1-6, 9-16 or 21-33 (got: '34')"),
        ("--mito-genetic-code", "7", "NCBI genetic code: 1-6, 9-16 or 21-33 (got: '7')"),
        ("--mito-genetic-code", "four", "NCBI genetic code: 1-6, 9-16 or 21-33 (got: 'four')"),
        ("--threads", "0", "at least 1"),
    ],
)
def test_bad_options_are_refused(tmp_path, capsys, option, value, expected_text):
    """Missing files and settings out of range give one readable error."""
    script, database = make_case(tmp_path)
    if value.endswith((".fasta", ".sh")):
        value = str(tmp_path / value)
    assert_stopped_with_message(tmp_path, full_argv(tmp_path, script, database, **{option: value}), capsys, expected_text)


def test_a_database_one_level_too_high_is_refused(tmp_path, capsys):
    """The folder above the species folders holds none itself."""
    script, database = make_case(tmp_path)
    assert_stopped_with_message(tmp_path, full_argv(tmp_path, script, database, **{"--database": str(tmp_path)}), capsys,
                                "holds no species folder")


def test_a_position_list_that_cannot_be_read_stops_before_chroquetas_runs(tmp_path, capsys):
    """Another header in a target's .txt file of the database: a clear error, and no run."""
    script, database = make_case(tmp_path)
    with open(os.path.join(database, "Botrytis_cinerea", "Bos1.txt"), "w") as txt:
        txt.write("Pos\tRef\tAlt\tDrugs\n")
    assert_stopped_with_message(tmp_path, full_argv(tmp_path, script, database), capsys,
                                "not a position list of the ChroQueTas 1.0.1 database")
    assert stand_in_options(tmp_path, "genome") is None


def test_an_empty_genome_is_refused(tmp_path, capsys):
    """A genome file without a sequence cannot be screened."""
    script, database = make_case(tmp_path)
    open(tmp_path / "contigs_final.fasta", "w").close()
    assert_stopped_with_message(tmp_path, full_argv(tmp_path, script, database), capsys, "holds no sequence")


def test_a_path_with_a_space_is_refused(tmp_path, capsys):
    """ChroQueTas does not quote its paths: an output folder with a space stops the script before it runs."""
    script, database = make_case(tmp_path)
    work_dir = str(tmp_path / "with space" / "chroquetas")
    assert_stopped_with_message(tmp_path, full_argv(tmp_path, script, database, **{"--work-dir": work_dir}), capsys,
                                "contains a space")
    assert stand_in_options(tmp_path, "genome") is None


def test_command_line_run_exits_0_and_a_failure_exits_1(tmp_path):
    """Run as a program: exit status 0 and nothing on standard output; a crashed run: exit status 1, no traceback."""
    script, database = make_case(tmp_path)
    finished = subprocess.run([sys.executable, SCRIPT_PATH] + full_argv(tmp_path, script, database),
                              capture_output=True, text=True)
    assert finished.returncode == 0, finished.stderr
    assert finished.stdout == ""
    assert f"[FunFlux] amr_chroquetas {SAMPLE}: species folder Botrytis_cinerea" in finished.stderr

    crashed_case = tmp_path / "crashed"
    crashed_case.mkdir()
    crashed = {"files": {}, "log": "ERROR: something\n", "exit": 1}
    script, database = make_case(crashed_case, {"genome": crashed, "mito": crashed})
    finished = subprocess.run([sys.executable, SCRIPT_PATH] + full_argv(crashed_case, script, database),
                              capture_output=True, text=True)
    assert finished.returncode == 1
    assert "[FunFlux] amr_chroquetas: error:" in finished.stderr
    assert "Traceback" not in finished.stderr


# ======================================================================
# THE REAL CHROQUETAS ON A PIECE OF THE B05.10 GENOME
# ======================================================================

def read_region():
    """Return the sequence of the B05.10 beta-tubulin region (test_data)."""
    lines = []
    with open(REGION_FASTA) as fasta:
        for line in fasta:
            if not line.startswith(">"):
                lines.append(line.strip())
    return "".join(lines)


def real_paths():
    """The real ChroQueTas script and database of FUNFLUX_CHROQUETAS_TEST_DIR."""
    return (os.path.join(CHROQUETAS_TEST_DIR, "ChroQueTas.sh"), os.path.join(CHROQUETAS_TEST_DIR, "FungAMR_db"))


def reference_protein(target):
    """Return the header and the sequence of a Botrytis cinerea reference protein of the real database."""
    _, database = real_paths()
    with open(os.path.join(database, "Botrytis_cinerea", target + ".faa")) as faa:
        lines = faa.read().splitlines()
    return lines[0], "".join(lines[1:])


def add_species(database, folder, proteins, positions_text):
    """Write a made-up species folder: {target: protein sequence}, each with the same position list."""
    os.makedirs(os.path.join(database, folder))
    for target, sequence in proteins.items():
        with open(os.path.join(database, folder, target + ".faa"), "w") as faa:
            faa.write(f">{target}\n{sequence}\n")
        with open(os.path.join(database, folder, target + ".txt"), "w") as txt:
            txt.write(positions_text)


def run_real(tmp_path, genome_records, species, database=None, mito_contigs=""):
    """Run the script with the real ChroQueTas on a genome built from (name, sequence) pairs; return the summary rows."""
    script, real_database = real_paths()
    os.makedirs(tmp_path, exist_ok=True)
    genome = tmp_path / "genome_input.fasta"
    write_fasta(genome, genome_records)
    argv = full_argv(tmp_path, script, database or real_database, **{
        "--genome": str(genome), "--mito-contigs": None, "--species": species,
    })
    if mito_contigs == "":
        argv.append("--mito-contigs")
    else:
        argv.extend(["--mito-contigs", mito_contigs])
    chroquetas_amr.main(argv)
    return summary_by_target(tmp_path / "amr" / "amr_summary.tsv")


@needs_chroquetas
def test_real_run_finds_e198a_in_the_b0510_beta_tubulin(tmp_path, capsys):
    """Beta-tubulin evaluated with E198A; the other targets are not in the 3.8-kb piece."""
    rows = run_real(tmp_path, [("tub_region", read_region())], "Botrytis cinerea")
    assert rows["Beta-tubulin"]["status"] == "evaluated"
    assert rows["Beta-tubulin"]["verdict"] == "resistance"
    assert rows["Beta-tubulin"]["fungamr_mutations"] == "E198A (resistance)"
    assert rows["Beta-tubulin"]["resistance_to"] == "Carbendazim (E198A, 7); Iprodione (E198A, 8)"
    assert rows["Beta-tubulin"]["no_resistance_to"] == "-"
    assert rows["Beta-tubulin"]["positions_checked"] == "3"
    for target in ["Bos1", "OS1", "Cytochrome_b"]:
        assert (rows[target]["status"], rows[target]["verdict"]) == ("not found", "-")
    # Both runs were made: without the screen the mitochondrial run reads the genome too,
    # and finds nothing (which must not be a failure).
    assert sorted(os.listdir(tmp_path / "amr" / "chroquetas")) == ["mitochondrial", "nuclear"]
    positions = read_table(tmp_path / "amr" / "amr_positions.tsv")
    assert [(row["position"], row["call"], row["verdict"], row["known_resistance_changes"]) for row in positions] == [
        ("198", "FungAMR mutation", "resistance", "E198A, E198G, E198K, E198V"),
        ("200", "no mutation", "wild type", "F200Y"),
        ("351", "no mutation", "wild type", "T351I"),
    ]
    assert positions[0]["resistance_to"] == f"Carbendazim (7: {WORDS_7}); Iprodione (8: {WORDS_8})"
    assert positions[0]["fungamr_drugs"] == CELL_E198A
    assert "FungAMR mutations: Beta-tubulin E198A (resistance: Carbendazim 7, Iprodione 8)" in capsys.readouterr().err


@needs_chroquetas
def test_real_run_on_e198v_gives_resistance_and_no_resistance_to_iprodione(tmp_path, capsys):
    """Codon 198 changed to valine: FungAMR mutation E198V, verdict resistance, iprodione in no_resistance_to."""
    # In the region file the gene lies on the forward strand. Its sixth exon starts at
    # base 1,525 with codon 54, so codon 198 is at bases 1,957-1,959 (index 1956 counted
    # from 0): GCG, alanine, in B05.10. GTG makes it valine. The mutated genome is built
    # here, so that no second data file is needed.
    region = read_region()
    assert region[1956:1959] == "GCG"
    with_valine = region[:1956] + "GTG" + region[1959:]
    rows = run_real(tmp_path, [("tub_region", with_valine)], "Botrytis cinerea")

    positions = read_table(tmp_path / "amr" / "amr_positions.tsv")
    assert (positions[0]["position"], positions[0]["reference_aa"], positions[0]["found_aa"]) == ("198", "E", "V")
    assert positions[0]["call"] == "FungAMR mutation"
    assert positions[0]["verdict"] == "resistance"
    assert positions[0]["resistance_to"] == (
        f"Benzimidazole (8: {WORDS_8}); Carbendazim (7: {WORDS_7}); Diethofencarb (8: {WORDS_8})"
    )
    assert positions[0]["no_resistance_to"] == f"Iprodione (-8: {WORDS_8})"
    assert positions[0]["known_resistance_changes"] == "E198A, E198G, E198K, E198V"
    assert positions[0]["fungamr_drugs"] == CELL_E198V

    tubulin = rows["Beta-tubulin"]
    assert tubulin["verdict"] == "resistance"
    assert tubulin["fungamr_mutations"] == "E198V (resistance)"
    assert tubulin["resistance_to"] == "Benzimidazole (E198V, 8); Carbendazim (E198V, 7); Diethofencarb (E198V, 8)"
    assert tubulin["no_resistance_to"] == "Iprodione (E198V, -8)"
    assert (
        "FungAMR mutations: Beta-tubulin E198V (resistance: Benzimidazole 8, Carbendazim 7, Diethofencarb 8; "
        "no resistance: Iprodione -8)"
    ) in capsys.readouterr().err


@needs_chroquetas
def test_real_run_on_two_copies_of_the_region_reports_copy_2(tmp_path):
    """The region on two contigs: fragments 1 and 2, both with E198A."""
    region = read_region()
    rows = run_real(tmp_path, [("copy_a", region), ("copy_b", region)], "Botrytis cinerea")
    assert rows["Beta-tubulin"]["fungamr_mutations"] == "E198A (resistance), E198A (copy 2) (resistance)"
    assert rows["Beta-tubulin"]["resistance_to"] == "Carbendazim (E198A, 7); Iprodione (E198A, 8)"
    assert rows["Beta-tubulin"]["note"].startswith("found 2 times (copies 1, 2)")


@needs_chroquetas
def test_real_run_on_random_sequence_finds_nothing(tmp_path):
    """50 kb of random sequence: every target not found, and no failure although ChroQueTas exits with 1."""
    rows = run_real(tmp_path, [("random_contig", random_dna(50000, 5))], "Botrytis cinerea")
    for target in BOTRYTIS_TARGETS:
        assert rows[target]["status"] == "not found"


@needs_chroquetas
def test_real_run_new_mutation_insertions_below_thresholds_and_skipped(tmp_path):
    """Made-up species built from the Botrytis beta-tubulin show the four other outcomes."""
    _, tubulin = reference_protein("Beta-tubulin")
    _, cytochrome_b = reference_protein("Cytochrome_b")
    database = str(tmp_path / "db")
    os.makedirs(database)
    positions = (
        "Position\tReference\tMutation\tFungicide_resistance\n"
        "198\tE\tA\tCarbendazim(7/NA),Iprodione(8/NA)\n"
        "200\tL\tY\tCarbendazim(8/NA)\n"
    )
    # A reference with L at 200 where B05.10 has F: L200F is a change FungAMR does not list.
    add_species(database, "Fakeus_newmut", {"Beta-tubulin": tubulin[:199] + "L" + tubulin[200:]}, positions)
    # Three residues fewer in the reference: the B05.10 protein has an insertion against it.
    add_species(database, "Fakeus_insertion", {"Beta-tubulin": tubulin[:99] + tubulin[102:]}, positions)
    # Every seventh residue changed: about 86 % identity, below the 90 % filter.
    distant = list(tubulin)
    for index in range(3, len(distant), 7):
        distant[index] = "W" if distant[index] != "W" else "C"
    add_species(database, "Fakeus_distant", {"Beta-tubulin": "".join(distant)}, positions)
    # Vma1 is not in the region; ChroQueTas then cuts "Vma1" out of "Vma16" and skips it.
    add_species(database, "Fakeus_overlap", {"Vma1": cytochrome_b, "Vma16": tubulin}, positions)

    region = [("tub_region", read_region())]
    rows = run_real(tmp_path / "newmut", region, "Fakeus newmut", database)
    assert rows["Beta-tubulin"]["fungamr_mutations"] == "E198A (resistance)"
    assert rows["Beta-tubulin"]["new_mutations"] == "L200F"
    assert rows["Beta-tubulin"]["verdict"] == "resistance"
    positions = read_table(tmp_path / "newmut" / "amr" / "amr_positions.tsv")
    assert [(row["position"], row["verdict"], row["known_resistance_changes"]) for row in positions] == [
        ("198", "resistance", "E198A"), ("200", "unknown", "L200Y"),
    ]
    rows = run_real(tmp_path / "insertion", region, "Fakeus insertion", database)
    assert (rows["Beta-tubulin"]["status"], rows["Beta-tubulin"]["note"]) == ("not evaluated", chroquetas_amr.NOTE_INSERTIONS)
    rows = run_real(tmp_path / "distant", region, "Fakeus distant", database)
    assert (rows["Beta-tubulin"]["status"], rows["Beta-tubulin"]["note"]) == ("not found", chroquetas_amr.NOTE_NOT_FOUND)
    rows = run_real(tmp_path / "overlap", region, "Fakeus overlap", database)
    assert rows["Vma1"]["status"] == "not found"
    assert (rows["Vma16"]["status"], rows["Vma16"]["note"]) == ("not evaluated", chroquetas_amr.NOTE_SKIPPED)



@needs_chroquetas
def test_real_run_premature_stop_gives_no_alignment_and_not_evaluated(tmp_path):
    """A stop codon at a FungAMR position: ChroQueTas sees a gap there, so the target is not evaluated."""
    _, tubulin = reference_protein("Beta-tubulin")
    database = str(tmp_path / "db")
    os.makedirs(database)
    # The only listed position is 351, as a nonsense mutation (X), the way FungAMR writes one.
    add_species(database, "Fakeus_stop", {"Beta-tubulin": tubulin},
                "Position\tReference\tMutation\tFungicide_resistance\n351\tT\tX\tCarbendazim(4/NA)\n")
    # Codon 351 of the gene (threonine, ACC) is at bases 2472-2474 of the region (index
    # 2471 counted from 0); TAA makes it a stop.
    region = read_region()
    assert region[2471:2474] == "ACC"
    with_stop = region[:2471] + "TAA" + region[2474:]
    rows = run_real(tmp_path / "stop", [("tub_region", with_stop)], "Fakeus stop", database)
    assert rows["Beta-tubulin"]["status"] == "not evaluated"
    assert rows["Beta-tubulin"]["note"] == chroquetas_amr.NOTE_NOTHING_COMPARED + "; positions not compared: 351 (no alignment)"
    positions = read_table(tmp_path / "stop" / "amr" / "amr_positions.tsv")
    assert [(row["position"], row["found_aa"], row["call"]) for row in positions] == [("351", "-", "no alignment")]


@needs_chroquetas
def test_real_run_gene_cut_at_a_contig_end_unreadable_codon_and_real_deletion(tmp_path):
    """A gap beyond the found copy and an X are not compared; a deletion inside the copy stays a FungAMR mutation."""
    _, tubulin = reference_protein("Beta-tubulin")
    database = str(tmp_path / "db")
    os.makedirs(database)
    # A made-up position list: deletions (written -) at 20 and 105, E198A and a stop at
    # 198 (FungAMR writes a stop as X), and F200Y.
    add_species(database, "Fakeus_gaps", {"Beta-tubulin": tubulin},
                "Position\tReference\tMutation\tFungicide_resistance\n"
                "20\tF\t-\tCarbendazim(4/NA)\n"
                "105\tH\t-\tCarbendazim(4/NA)\n"
                "198\tE\tA\tCarbendazim(7/NA)\n"
                "198\tE\tX\tCarbendazim(4/NA)\n"
                "200\tF\tY\tCarbendazim(7/NA)\n")
    region = read_region()
    # Codon 105 (histidine, CAT) is at bases 1678-1680 of the region (index 1677
    # counted from 0); codon 198 at index 1956 (see the E198V test).
    assert region[1677:1680] == "CAT"
    assert region[1956:1959] == "GCG"

    # The gene cut at a contig end: the contig starts at index 1520, 4 bases before
    # codon 54, so the copy lacks amino acids 1 to 53. Codon 198 is NNN, an unreadable
    # codon (index 1956 - 1520 = 436 in the cut contig).
    cut = region[1520:]
    cut = cut[:436] + "NNN" + cut[439:]
    rows = run_real(tmp_path / "cut", [("cut_region", cut)], "Fakeus gaps", database)
    positions = read_table(tmp_path / "cut" / "amr" / "amr_positions.tsv")
    assert [(row["position"], row["found_aa"], row["call"], row["verdict"]) for row in positions] == [
        ("20", "-", "no alignment", "not compared"),
        ("105", "H", "no mutation", "wild type"),
        ("198", "X", "position not found", "not compared"),
        ("200", "F", "no mutation", "wild type"),
    ]
    assert (rows["Beta-tubulin"]["status"], rows["Beta-tubulin"]["verdict"]) == ("evaluated", "wild type")
    assert rows["Beta-tubulin"]["note"] == "positions not compared: 20 (no alignment), 198 (position not found)"

    # A real deletion: codon 105 removed from the whole gene. The gap lies inside the
    # copy, so it is the deletion FungAMR lists.
    with_deletion = region[:1677] + region[1680:]
    rows = run_real(tmp_path / "deletion", [("tub_region", with_deletion)], "Fakeus gaps", database)
    positions = read_table(tmp_path / "deletion" / "amr" / "amr_positions.tsv")
    assert [(row["position"], row["found_aa"], row["call"], row["verdict"]) for row in positions] == [
        ("20", "F", "no mutation", "wild type"),
        ("105", "-", "FungAMR mutation", "resistance"),
        ("198", "A", "FungAMR mutation", "resistance"),
        ("200", "F", "no mutation", "wild type"),
    ]
    assert rows["Beta-tubulin"]["fungamr_mutations"] == "H105- (resistance), E198A (resistance)"
