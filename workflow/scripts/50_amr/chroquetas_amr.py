#!/usr/bin/env python3
"""Screen a fungal genome for the antifungal resistance mutations of FungAMR with ChroQueTas.

BIOLOGY AND WHY THIS SCRIPT EXISTS
----------------------------------
Fungi become resistant to antifungal drugs and agricultural fungicides mostly through
point mutations in the protein the drug binds: beta-tubulin for the benzimidazoles
(carbendazim), cytochrome b for the QoI fungicides (strobilurins), Cyp51 (Erg11) for
the azoles, succinate dehydrogenase (Sdh) for the SDHIs, the glucan synthase Fks for
the echinocandins, and so on. A changed amino acid at one known position of that
protein is enough. Resistance genes acquired from elsewhere, the rule in bacteria,
are rare in fungi, so a catalogue of resistance genes finds nothing here: the
question is, for each target protein of the species, which amino acid the sample
carries at each position where a change is known to matter.

FungAMR (Bédard et al. 2025) collects those mutations from the literature, each with
the drugs it was studied with and, per drug, a confidence score for the evidence that
it causes resistance or that it does not (EVIDENCE LEVELS AND VERDICTS below).
ChroQueTas, described in the same paper, screens a genome for them.

WHAT CHROQUETAS DOES
--------------------
ChroQueTas 1.0.1 and its FungAMR database cover 57 species and 8,285 mutations: only
single mutations, tested on their own rather than in combination, and only species
with a trusted reference protein for the target. An entry can carry evidence for
resistance, against it, or both: 4,655 of the 8,285 carry only evidence against, most
of them from deep mutational scanning of Candida albicans Cyp51 and Saccharomyces
cerevisiae Fcy1.

For every target protein of the species (one <Target>.faa file, with a <Target>.txt
list of positions, in the species' folder of the database), ChroQueTas
  1. finds the gene in the genome with miniprot (Li 2023), which aligns the
     reference protein to the DNA across introns, and translates every copy found
     (fragment 1, 2, ...; more than one copy means a duplicated gene, or a region the
     assembly holds twice);
  2. keeps a copy when BLASTP gives it more than 90 % identity with the reference
     over more than 75 % of the reference's length;
  3. aligns the kept copy with the reference with MAFFT; when the alignment puts a
     gap into the reference (the found protein has extra residues: an "insertion"),
     the positions cannot be numbered reliably and are not compared;
  4. otherwise reads, at every FungAMR position, the amino acid of the reference and
     of the sample, and calls it "FungAMR MUTATION" (the change is in FungAMR, with
     its drugs, whatever its evidence says), "New mutation" (another change at a
     known position) or "No mutation"; "No alignment" or "Position not found" when
     the position could not be read (see CALLS below).
Not detected: changes in promoters, overexpression, gene duplication or loss as
such, aneuploidy, and combinations of mutations.

TWO RUNS: THE NUCLEAR GENOME AND THE MITOCHONDRION
--------------------------------------------------
A ChroQueTas run translates the whole genome with one genetic code. For the nuclear
genome ChroQueTas chooses it itself from the species: code 12 for the eight species
of the CTG clade (Candida albicans and relatives, which read CTG as serine), code 1
otherwise. One target is encoded in the mitochondrion: cytochrome b (Cytochrome_b,
in 15 species of the database). Fungal mitochondria use another code: table 4 (mold
mitochondrial code, where TGA codes for tryptophan instead of stop), or table 3 in
the Saccharomycetaceae (yeast mitochondrial code). Translated with the nuclear code,
every TGA of the cytochrome b gene, a tryptophan in the mitochondrion, would read as a
stop. So this script runs ChroQueTas twice, each time into its own folder:
  * nuclear: on the genome, without -c, so ChroQueTas picks the nuclear code; it
    gives the results of every target except Cytochrome_b;
  * mitochondrial: only when the species has Cytochrome_b, with -c set to the
    mitochondrial code (--mito-genetic-code; auto gives 3 for the four
    Saccharomycetaceae of the database and 4 for every other species), on the
    mitochondrial contigs when the run's screen splits the mitochondrion off,
    otherwise on the whole genome; it gives the result of Cytochrome_b only.
A target found by the other run (cytochrome b in a nuclear copy of mitochondrial
DNA, say) is ignored. When the screen ran and found no mitochondrial contig, there is
no mitochondrial run, and Cytochrome_b is "not evaluated".

THE SPECIES
-----------
The species comes from the annotation table (column 2), and must be matched to a
folder of the database. Tried in this order, ignoring upper and lower case:
  1. the whole name, spaces as underscores and "f. sp." left out, so that "Blumeria
     graminis f. sp. hordei" finds Blumeria_graminis_hordei;
  2. the first two words (a strain name after the species is common);
  3. SYNONYMS below: other names of the same species, older or newer, whose database
     folder carries a different name (Candida auris -> Candidozyma_auris), each
     checked against NCBI Taxonomy.
No match: the sample is "not screened", which is not a failure. A NOTE line says so,
and names the database folders of the same genus, if there are any.

WHAT EACH TARGET GETS: THE STATUSES
-----------------------------------
    evaluated        at least one FungAMR position of at least one copy of the
                     protein was compared; the positions are in amr_positions.tsv.
    not evaluated    no position was compared: the protein was found with
                     insertions against the reference (step 3 above; a common cause
                     is an intron the protein aligner could not place, such as an
                     intron of a mitochondrial cytochrome b gene); or it was found,
                     but every FungAMR position lies in a gap of the found protein, at
                     a codon that could not be read, or outside the protein; or
                     ChroQueTas skipped the target (see
                     read_run_results()); or, for Cytochrome_b, there was no
                     mitochondrial contig to search.
    not found        no copy passed the identity and coverage filter of step 2.
    not screened     the species is not in the database (one row, target "-").

EVIDENCE LEVELS AND VERDICTS: DOES A CHANGE CAUSE RESISTANCE?
-------------------------------------------------------------
ChroQueTas calls every change FungAMR lists a "FungAMR mutation", but FungAMR also lists
changes that were tested and found NOT to cause resistance: more than half of the rows
of the database (see WHAT CHROQUETAS DOES above). This script therefore says, for every
change, which way the evidence points. Evidence of no resistance is the harder kind to
judge: a susceptible result can come from an assay that was not sensitive enough, and
such assays rarely test a resistant control alongside. The FungAMR authors kept these
reports in the database but left them out of their own analyses for that reason. So
"no resistance" means that the reports found no resistance to the drugs named, nothing
more.

FungAMR gives each drug of a change two scores, written Drug(positive/negative)
(ChroQueTas wiki, "Confidence score for antimicrobial resistance"): a positive level,
1 (strongest) to 8 (weakest), for evidence that the change causes resistance to that
drug; a negative level, -1 to -8, for evidence of the same kind that it does not; NA
when there is no report of that kind. The level names the kind of study, so one set of
words (EVIDENCE_LEVELS) reads right after either sign: after a positive level the
study found resistance, after a negative one it found none. The tables carry the words
without the explanations in brackets below.
    1  the change made in the species and measured (in a susceptible strain)
    2  the change made in the gene, expressed in another species (and measured)
    3  bulk competition assay (such as deep mutational scanning)
    4  experimental evolution (the change arose in independent replicates)
    5  gene deletion
    6  gene overexpression or duplication
    7  population association study (GWAS, QTL or classical genetics)
    8  natural isolates, not validated further (a clinical or field strain)
Levels 7 and 8 rest on association, not on a test of the change itself. Field isolates
are often resistant to several classes of fungicides at once, so a drug at level 7 or 8
that acts on another target can record a resistance found together with the change
rather than caused by it.

Each drug of a FungAMR change goes, with its level, into one or both lists:
    a positive level only   -> resistance_to
    a negative level only   -> no_resistance_to
    both                    -> both lists, each with its own level (a conflicting drug)
The verdict of one position (amr_positions.tsv):
    FungAMR mutation        resistance      at least one drug has a positive level only;
                            conflicting     otherwise, at least one drug has both;
                            no resistance   otherwise, every drug named has a negative
                                            level only;
                            unknown         no drug is named
    new mutation            unknown         a change FungAMR does not list, at a position
                                            where it lists others
    no mutation             wild type       the reference amino acid
    no alignment,           not compared
    position not found
The verdict of a target (amr_summary.tsv) is the first of these found among the
verdicts of its compared positions, all copies together: resistance, conflicting,
unknown, no resistance, wild type; "-" for a target that was not evaluated, not found or
not screened (its status says why).
Example: Botrytis cinerea beta-tubulin E198V, Benzimidazole(8/NA),Carbendazim(7/NA),
Diethofencarb(8/NA),Iprodione(NA/-8): evidence of resistance to three drugs and of no
resistance to iprodione; the verdict is resistance.

THE TABLES
----------
"-" marks an empty cell. Drugs are separated by "; ", because a FungAMR drug name can
hold a comma ("Fludioxonil,_Iprodione"); FungAMR's names are kept as they are.
amr_positions.tsv, one row per FungAMR position of every compared copy of a target
(POSITIONS_COLUMNS):
    sample, target, compartment (nuclear or mitochondrial), fragment (the copy),
    position, reference_aa, found_aa, call (CALLS below), verdict;
    resistance_to, no_resistance_to   for a FungAMR mutation, each drug with its level
        and the words of that level: "Carbendazim (7: population association study)",
        "Iprodione (-8: natural isolates, not validated further)";
    known_resistance_changes   for every compared row, the changes the database lists
        at this position with a positive level for at least one drug, in the order of
        the database file: "E198A, E198G, E198K, E198V". On a wild-type row they are the
        changes that would matter there; on an unknown row, the changes known there;
    fungamr_drugs   ChroQueTas's drug text, unchanged; "-" where it names no drug (NA,
        or Unknown for a new mutation).
amr_summary.tsv, one row per target of the species (SUMMARY_COLUMNS):
    sample, species (the database folder), target, compartment, status,
    positions_checked, verdict;
    fungamr_mutations   each FungAMR change with its verdict: "E198A (resistance)",
        "E198A (copy 2) (resistance)" for a second copy;
    new_mutations   "F200L", "F200L (copy 2)";
    resistance_to, no_resistance_to   the drugs of the target's FungAMR changes, with
        the change and the level: "Carbendazim (E198A, 7); Iprodione (E198A, 8)"; a
        drug named by two changes appears once for each;
    note.
amr_summary_mqc.tsv: the summary with MultiQC's header lines.

WHERE THIS SCRIPT FITS IN THE WORKFLOW
--------------------------------------
    1. Rule chroquetas_db (workflow/rules/shared/50_amr.smk) downloads the unmodified
       ChroQueTas script and FungAMR database once per run (CHROQUETAS_SCRIPT,
       CHROQUETAS_DB_DIR; the names in capitals are the names used in
       workflow/rules/shared/00_common.smk). The step is off unless
       parameters.amr.chroquetas is on, because of the licence (CC BY-NC-ND 4.0).
    2. THIS SCRIPT is run once per sample by rule amr_chroquetas (the same rule file),
       in the chroquetas conda environment (workflow/envs/chroquetas.yaml), which
       holds the programs ChroQueTas calls. It runs ChroQueTas as
       "bash <script> ..." and writes the three tables below.
    3. Rule multiqc (workflow/rules/shared/90_report.smk) reads the _mqc copy, so the
       table appears in the report; the other two tables are the deliverables.

INPUTS
------
    --chroquetas-script  CHROQUETAS_SCRIPT, ChroQueTas.sh as downloaded.
    --database           CHROQUETAS_DB_DIR, the folder that directly holds the species
                         folders of the FungAMR database. Besides handing it to
                         ChroQueTas, this script reads the position list of each target
                         (<species folder>/<Target>.txt) for known_resistance_changes.
    --genome             FINAL_CONTIGS, the nuclear genome after the screens (in contigs
                         mode without the screens, all contigs).
    --mito-contigs       MITO_CONTIGS, the mitochondrial contig(s) the screen set aside
                         (possibly an empty file), or an empty value when the run does
                         not split the mitochondrion off (no screen, or
                         parameters.mitochondria.mode off or report).
    --species            the species name of the annotation table.
    --mito-genetic-code  auto, or the number of an NCBI genetic code (1-6, 9-16 or
                         21-33) for every sample (parameters.amr.mito_genetic_code).
    --threads            threads for miniprot and MAFFT (default 4).
    --sample             the sample name, for the tables and the messages.

OUTPUTS
-------
    --work-dir       AMR_CHROQUETAS_DIR, the chroquetas/ folder this script fills:
                     nuclear/ and mitochondrial/, ChroQueTas's own output of each run,
                     untouched (a run that was not needed has no folder).
    --out-positions  AMR_POSITIONS, amr_positions.tsv: one row per FungAMR position of
                     every compared copy (THE TABLES above).
    --out-summary    AMR_SUMMARY, amr_summary.tsv: one row per target of the species
                     (THE TABLES above).
    --out-mqc        AMR_SUMMARY_MQC, amr_summary_mqc.tsv: the summary behind
                     MultiQC's custom-content header lines, rows named
                     "<sample>: <target>".
    standard error   per run, the command and ChroQueTas's own messages; then the
                     summary line "[FunFlux] amr_chroquetas <sample>: ..." (the species
                     folder, the targets per status, a tally of the target verdicts,
                     each FungAMR mutation with its verdict and drug levels, such as
                     "Beta-tubulin E198A (resistance: Carbendazim 7, Iprodione 8)"),
                     or, for a species not in the database, one
                     "[FunFlux] NOTE amr_chroquetas <sample>: ..." line. Rule
                     amr_chroquetas writes the stream to its log and copies the
                     [FunFlux] lines to the console.

Exit status 0 when the three tables were written. A missing input, a database folder
without species folders, a position list of the database that does not look like the
one of ChroQueTas 1.0.1, a setting out of range, an output folder of a run that
already exists, a path with a space (ChroQueTas cannot handle one), or a ChroQueTas
run that did not finish stops the script with a message starting
"[FunFlux] amr_chroquetas: error:" and exit status 1, and no table is written. For a
run that did not finish, the message carries the last lines of ChroQueTas's output.

REFERENCES
----------
Bédard C, et al. (2025). FungAMR: a comprehensive database for investigating fungal
    mutations associated with antimicrobial resistance. Nature Microbiology 10(9):
    2338-2352. https://doi.org/10.1038/s41564-025-02084-7 (FungAMR and ChroQueTas)
Li H (2023). Protein-to-genome alignment with miniprot. Bioinformatics 39(1):
    btad014. https://doi.org/10.1093/bioinformatics/btad014
Camacho C, et al. (2009). BLAST+: architecture and applications. BMC Bioinformatics
    10: 421. https://doi.org/10.1186/1471-2105-10-421
Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7:
    improvements in performance and usability. Molecular Biology and Evolution 30(4):
    772-780. https://doi.org/10.1093/molbev/mst010
"""

# Only Python's standard library is used; ChroQueTas is run as a separate program with
# subprocess. The chroquetas conda environment holds a Python 3.12 and the programs
# ChroQueTas calls.
import argparse
import os
import re
import subprocess
import sys
import tempfile
import time


# ======================================================================
# SETTINGS
# ======================================================================

# The prefix of the messages that rule amr_chroquetas copies to the console.
MESSAGE_NAME = "amr_chroquetas"

# The two runs: the name of each run's output folder inside --work-dir, and the name the
# input FASTA gets for ChroQueTas. ChroQueTas names its files after the input file, up to
# its last dot (genome.fasta -> genome.ChroQueTaS.AMR_summary.txt), so a genome file
# called B0510.final.fasta would give the prefix "B0510.final"; the input is therefore
# handed over as a link with a fixed name, and the prefix is always "genome" or "mito".
RUN_NUCLEAR = "nuclear"
RUN_MITOCHONDRIAL = "mitochondrial"
RUN_INPUT_NAMES = {RUN_NUCLEAR: "genome.fasta", RUN_MITOCHONDRIAL: "mito.fasta"}

# The targets of the database encoded in the mitochondrion. In the database of
# ChroQueTas 1.0.1 (the archive db-v20250811.tgz) only cytochrome b; every other target
# is nuclear.
MITOCHONDRIAL_TARGETS = {"Cytochrome_b"}

# The mitochondrial genetic code that --mito-genetic-code auto gives. NCBI Taxonomy
# (checked 2026-10-11) lists code 3 (yeast mitochondrial) for these four species of the
# database, all Saccharomycetaceae, and code 4 (mold mitochondrial) for every other
# species of the database.
YEAST_MITOCHONDRIAL_CODE = 3
MOLD_MITOCHONDRIAL_CODE = 4
YEAST_MITOCHONDRIAL_CODE_FOLDERS = {
    "Saccharomyces_cerevisiae",
    "Saccharomyces_paradoxus",
    "Nakaseomyces_glabratus",
    "Kluyveromyces_marxianus",
}
# The genetic codes a number given by hand may name: the NCBI genetic codes. The numbers
# 7, 8 and 17 to 20 are not NCBI codes, and miniprot, which ChroQueTas uses to translate,
# refuses them.
NCBI_GENETIC_CODES = {
    1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33,
}
NCBI_GENETIC_CODES_TEXT = "1-6, 9-16 or 21-33"
DEFAULT_THREADS = 4

# Other names of species whose database folder carries a different name of the same
# taxon. Some keys are the older name (Candida glabrata, now Nakaseomyces glabratus),
# some the newer one (Parastagonospora nodorum, whose folder carries the older
# Phaeosphaeria nodorum). Each entry was checked against NCBI Taxonomy on 2026-10-11:
# the record of the taxon (its id after the name) holds both names, one as its
# scientific name and the other as a synonym. The keys are the first two words of the
# name, in lower case.
SYNONYMS = {
    "candida glabrata": "Nakaseomyces_glabratus",          # taxid 5478
    "candida auris": "Candidozyma_auris",                  # taxid 498019
    "parastagonospora nodorum": "Phaeosphaeria_nodorum",   # taxid 13684
    "pseudocercospora fijiensis": "Mycosphaerella_fijiensis",  # taxid 1873960
    "fulvia fulva": "Passalora_fulva",                     # taxid 5499
    "cladosporium fulvum": "Passalora_fulva",              # taxid 5499
    "oculimacula yallundae": "Tapesia_yallundae",          # taxid 86028
    "oculimacula acuformis": "Tapesia_acuformis",          # taxid 110162
    "issatchenkia orientalis": "Pichia_kudriavzevii",      # taxid 4909
    "candida krusei": "Pichia_kudriavzevii",               # taxid 4909
    "candida lusitaniae": "Clavispora_lusitaniae",         # taxid 36911
    "lodderomyces metapsilosis": "Candida_metapsilosis",   # taxid 273372
    "lodderomyces orthopsilosis": "Candida_orthopsilosis", # taxid 273371
    "pyricularia grisea": "Magnaporthe_grisea",            # taxid 148305
    "blumeria hordei": "Blumeria_graminis_hordei",         # taxid 2867405
}
# How a species was matched to its folder (THE SPECIES at the top, steps 1 to 3).
MATCH_NAME = "name"
MATCH_FIRST_TWO_WORDS = "first two words"
MATCH_SYNONYM = "synonym"

# The words of the eight evidence levels (EVIDENCE LEVELS AND VERDICTS at the top), for
# the resistance_to and no_resistance_to cells of amr_positions.tsv; 8 and -8 get the
# same words.
EVIDENCE_LEVELS = {
    1: "the change made in the species and measured",
    2: "the change made in the gene, expressed in another species",
    3: "bulk competition assay",
    4: "experimental evolution",
    5: "gene deletion",
    6: "gene overexpression or duplication",
    7: "population association study",
    8: "natural isolates, not validated further",
}

# The verdicts (EVIDENCE LEVELS AND VERDICTS at the top). The first five can be the
# verdict of a target; VERDICT_ORDER gives their precedence: a target takes the first of
# them found among its compared positions.
VERDICT_RESISTANCE = "resistance"
VERDICT_CONFLICTING = "conflicting"
VERDICT_UNKNOWN = "unknown"
VERDICT_NO_RESISTANCE = "no resistance"
VERDICT_WILD_TYPE = "wild type"
VERDICT_NOT_COMPARED = "not compared"
VERDICT_ORDER = [
    VERDICT_RESISTANCE, VERDICT_CONFLICTING, VERDICT_UNKNOWN, VERDICT_NO_RESISTANCE,
    VERDICT_WILD_TYPE,
]

# How the drugs of one cell are separated. Not a comma: some FungAMR drug names hold one
# ("Fludioxonil,_Iprodione").
DRUG_SEPARATOR = "; "

# The words ChroQueTas writes in the Result column of a per-target table, and the call
# this script writes for each. The first three are the usual ones. ChroQueTas calls a
# position a FungAMR mutation when the amino acid it found is one of the new amino
# acids the database lists there.
# "No alignment": the found protein has a gap at the position. That is a deletion, a
# part of the gene the assembly lacks (a gene cut at a contig end), or a premature stop
# codon: miniprot writes a stop as *, which MAFFT turns into a gap. So a nonsense
# mutation that FungAMR lists (written X; 242 rows in the database of ChroQueTas 1.0.1,
# most of them in Saccharomyces cerevisiae Fcy1 and Fur1) comes out as "no alignment",
# not as a FungAMR mutation; the copy's protein (<prefix>_<Target>.<n>.faa) shows the
# stop as *. Where FungAMR lists a deletion (written -, 22 rows), ChroQueTas matches the
# gap to it and calls a FungAMR mutation instead, also where the gap only means that
# the found protein ends early; correct_deletions_beyond_the_copy() handles that case.
# "Position not found": the position lies outside the alignment, or is not a number.
# Three database rows have no position number. ChroQueTas reports the insertion
# (Cryptococcus deuterogattii Uxs1, "1182ins(fs)") this way, and skips the duplication
# (Nakaseomyces glabratus Pdr1, an empty position) without a word. It reads the
# translocation (Candidozyma auris Ncp1, "ranslocation") as position 1 and calls it like
# any other position; read_target_table() turns that call into "position not found".
# read_target_table() corrects two more of ChroQueTas's calls:
#   * X found: miniprot writes X for a codon it cannot read, one with N bases (a gap in
#     a scaffold, say). That is no amino acid at all: neither the stop FungAMR writes as
#     X (miniprot writes a stop as *) nor a new mutation. ChroQueTas calls it a FungAMR
#     mutation where FungAMR lists a stop, otherwise a new mutation; this script calls it
#     "position not found", and found_aa keeps the X.
#   * the reference amino acid found, called a FungAMR mutation: the database lists one
#     "change" to the reference amino acid itself, Saccharomyces cerevisiae Hap1 N832N,
#     so ChroQueTas calls the wild type there a FungAMR mutation. This script calls it
#     "no mutation".
# The call as ChroQueTas made it is kept as well (chroquetas_call), because its own
# summary and stats files count those calls (check_counts_agree()).
CALL_FUNGAMR = "FungAMR mutation"
CALL_NEW = "new mutation"
CALL_NONE = "no mutation"
CALL_NO_ALIGNMENT = "no alignment"
CALL_NOT_FOUND = "position not found"
CALLS = {
    "FungAMR MUTATION": CALL_FUNGAMR,
    "New mutation": CALL_NEW,
    "No mutation": CALL_NONE,
    "No alignment": CALL_NO_ALIGNMENT,
    "Position not found": CALL_NOT_FOUND,
}
# The calls for which the two amino acids were really compared.
COMPARED_CALLS = {CALL_FUNGAMR, CALL_NEW, CALL_NONE}
# Two letters of the found amino acid that are not amino acids (CALLS above): the gap
# MAFFT writes where the found protein has nothing, and the X miniprot writes for a
# codon it cannot read.
GAP_LETTER = "-"
UNREADABLE_CODON_LETTER = "X"

# The statuses of a target, and the notes that go with the last three.
STATUS_EVALUATED = "evaluated"
STATUS_NOT_EVALUATED = "not evaluated"
STATUS_NOT_FOUND = "not found"
STATUS_NOT_SCREENED = "not screened"
NOTE_INSERTIONS = (
    "the extracted protein has insertions against the reference, so the positions were "
    "not compared; a common cause is an intron the protein aligner could not place, such "
    "as an intron of a mitochondrial cytochrome b gene"
)
# The note of a target whose copies were read position by position, but whose every
# FungAMR position lies in a gap of the found protein or outside it; the note then
# lists those positions.
NOTE_NOTHING_COMPARED = "found, but none of its FungAMR positions could be compared"
# ChroQueTas keeps a copy only ABOVE both thresholds (more than 90 %, more than 75 %),
# so a copy at exactly 90 % identity counts as not found.
NOTE_NOT_FOUND = "not found in the genome at more than 90 % identity and more than 75 % coverage"
NOTE_NO_MITO_CONTIG = (
    "no mitochondrial contig: the mitochondrion screen set none aside, so this "
    "mitochondrially encoded target was not searched"
)
NOTE_SKIPPED = (
    "found, but ChroQueTas did not compare its positions; ChroQueTas 1.0.1 drops a target "
    "whose name contains the name of a target it did not find"
)

# The text that marks an empty cell in the output tables.
EMPTY = "-"

# The columns of the three output tables, in this order (THE TABLES at the top).
POSITIONS_COLUMNS = [
    "sample", "target", "compartment", "fragment", "position", "reference_aa",
    "found_aa", "call", "verdict", "resistance_to", "no_resistance_to",
    "known_resistance_changes", "fungamr_drugs",
]
SUMMARY_COLUMNS = [
    "sample", "species", "target", "compartment", "status", "positions_checked",
    "verdict", "fungamr_mutations", "new_mutations", "resistance_to", "no_resistance_to",
    "note",
]

# The files of ChroQueTas 1.0.1 as this script reads them. A header that differs means
# another ChroQueTas, and the script stops rather than misread the file.
TARGET_TABLE_COLUMNS = ["Position", "Reference", "Query", "Result", "Fungicides"]
# The position list of a target in the database (<species folder>/<Target>.txt): one
# row per change FungAMR lists, with its drugs and scores.
DATABASE_TABLE_COLUMNS = ["Position", "Reference", "Mutation", "Fungicide_resistance"]
SUMMARY_FILE_COLUMNS = [
    "Protein", "Fragment", "Position_reference", "AA_reference", "AA_query", "Fungicide_resistance",
]
STATS_FILE_COLUMNS = ["Protein", "FungAMR_mutations", "New_mutations"]
# The two columns of ChroQueTas's BLASTP table of a copy (<prefix>_<Target>.<n>.blastp.txt)
# that give the first and the last position of the reference protein the copy was
# aligned to. The file is read by column name, from its header line.
BLASTP_REFERENCE_START = "sstart"
BLASTP_REFERENCE_END = "send"
# The one line of a summary file without any mutation ("No AMR mutations found in genome"),
# and the one line of a per-target table whose positions were not compared.
NO_MUTATIONS_TEXT = "No AMR mutations found in "
INSERTIONS_WARNING = "WARNING: genome contains insertions for "
# How the summary file marks a new mutation: in place of the drugs.
SUMMARY_NEW_MUTATION = "New mutation"
# The two log lines that say a target, or every target, could not be extracted.
NOT_EXTRACTED_LINE = re.compile(r"WARNING: We could not extract (\S+) from the genome")
NOTHING_EXTRACTED_TEXT = "could not extract any query proteins from the genome"

# How many of ChroQueTas's last output lines an error message carries.
LOG_TAIL_LINES = 20

# One drug of a Fungicides cell: the name, then (positive/negative). The name may itself
# hold commas and brackets ("1,10-Phenanthroline", "Azole(s)",
# "Fludioxonil,_Iprodione"), so an item ends only where a score in brackets is followed
# by a comma or the end of the text.
DRUG_ITEM = re.compile(r"(.*?)\((NA|[0-9]+)/(NA|-[0-9]+)\)(?:,|$)")


# ======================================================================
# THE DATABASE AND THE SPECIES
# ======================================================================

def list_species_folders(database):
    """List the species folders of the FungAMR database.

    Takes in:  --database, the folder that should directly hold the species folders.
    Does:      keeps the sub-folders that hold at least one target protein (.faa
               file); hidden folders are skipped.
    Produces:  the folder names, sorted (Aspergillus_fumigatus, Botrytis_cinerea, ...).
    Used by:   run(), to match the species.
    Stops with a clear error when there is no such folder: the path is wrong, or it
    points one level too high (at the folder that holds FungAMR_db, say).
    """
    if not os.path.isdir(database):
        raise ValueError(f"the database folder {database} does not exist")
    folders = []
    for name in sorted(os.listdir(database)):
        path = os.path.join(database, name)
        if name.startswith(".") or not os.path.isdir(path):
            continue
        holds_protein = False
        for file_name in os.listdir(path):
            if file_name.endswith(".faa"):
                holds_protein = True
                break
        if holds_protein:
            folders.append(name)
    if len(folders) == 0:
        raise ValueError(
            f"{database} holds no species folder (a folder with .faa files): give the folder "
            "that directly holds the species folders of the FungAMR database"
        )
    return folders


def list_targets(database, folder):
    """List the target proteins of one species: the names of its .faa files, sorted."""
    targets = []
    for file_name in sorted(os.listdir(os.path.join(database, folder))):
        if file_name.endswith(".faa"):
            targets.append(file_name[:-len(".faa")])
    return targets


def read_known_resistance_changes(path):
    """Read, from the database, the changes known to cause resistance at each position of a target.

    Takes in:  the position list of one target, <database>/<species folder>/<Target>.txt,
               the file ChroQueTas itself reads: the header DATABASE_TABLE_COLUMNS, then
               one row per change, for example "198  E  V  Benzimidazole(8/NA),...".
    Does:      keeps the changes with a positive level for at least one drug (evidence
               that they cause resistance), written reference amino acid, position,
               new amino acid ("E198V"), in the order of the file, each once. Changes
               with negative levels only are left out. Rows without a position number
               (an insertion, a duplication, a translocation; three in the database of
               ChroQueTas 1.0.1) are left out too: no compared position can match them.
    Produces:  {position as text: [changes]}, for example {"198": ["E198A", "E198G",
               "E198K", "E198V"], "200": ["F200Y"], "351": ["T351I"]} for the beta-tubulin
               of Botrytis cinerea. The one row whose new amino acid is the reference
               one (Saccharomyces cerevisiae Hap1 N832N, CALLS below) is left out as
               well: it is no change.
    Used by:   run(), for the known_resistance_changes column (position_rows()).
    Stops with a clear error when the file is missing, its first line is not
    DATABASE_TABLE_COLUMNS, or a row has another number of columns.
    """
    if not os.path.isfile(path):
        raise ValueError(f"the database file {path} does not exist")
    lines = read_lines(path)
    expected_header = "\t".join(DATABASE_TABLE_COLUMNS)
    if len(lines) == 0 or lines[0].split("\t") != DATABASE_TABLE_COLUMNS:
        raise ValueError(
            f"{path}: not a position list of the ChroQueTas 1.0.1 database (the first line "
            f"should be {expected_header!r})"
        )
    changes_by_position = {}
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(DATABASE_TABLE_COLUMNS):
            raise ValueError(
                f"{path}: a row has {len(cells)} columns instead of {len(DATABASE_TABLE_COLUMNS)}: {line!r}"
            )
        position, reference, mutation, fungicides = cells
        if not position.isdigit() or reference == "" or mutation == "":
            continue
        if mutation == reference:
            continue
        has_positive_level = False
        for _, positive, _ in parse_drugs(fungicides):
            if positive is not None:
                has_positive_level = True
        if not has_positive_level:
            continue
        change = f"{reference}{position}{mutation}"
        if position not in changes_by_position:
            changes_by_position[position] = []
        if change not in changes_by_position[position]:
            changes_by_position[position].append(change)
    return changes_by_position


def species_words(name):
    """Split a species name into words, leaving out "f. sp.".

    Takes in:  the species as the annotation table gives it, for example "Blumeria
               graminis f. sp. hordei" or "Botrytis_cinerea".
    Does:      reads underscores as spaces, removes "f. sp." and "f.sp." (forma
               specialis, a host-specific form; the database folders leave it out) and
               splits at spaces.
    Produces:  a list of words: ["Blumeria", "graminis", "hordei"].
    Used by:   resolve_species() and same_genus_folders().
    """
    text = name.replace("_", " ")
    text = re.sub(r"\bf\.\s*sp\.", " ", text, flags=re.IGNORECASE)
    return text.split()


def resolve_species(name, folders):
    """Find the database folder of a species, by name, first two words or synonym.

    Takes in:  the species name (--species) and the species folders of the database.
    Does:      tries the three steps of THE SPECIES at the top, ignoring case: the
               whole name joined by underscores, its first two words, and the
               SYNONYMS table (looked up with the first two words). A synonym counts
               only when the database holds its folder.
    Produces:  (folder, how) with how = MATCH_NAME, MATCH_FIRST_TWO_WORDS or
               MATCH_SYNONYM; (None, None) when nothing matches.
    Used by:   run().
    """
    folder_by_lower_name = {}
    for folder in folders:
        folder_by_lower_name[folder.lower()] = folder

    words = species_words(name)
    if len(words) == 0:
        return None, None

    whole_name = "_".join(words).lower()
    if whole_name in folder_by_lower_name:
        return folder_by_lower_name[whole_name], MATCH_NAME

    first_two = "_".join(words[:2]).lower()
    if first_two in folder_by_lower_name:
        return folder_by_lower_name[first_two], MATCH_FIRST_TWO_WORDS

    synonym_key = " ".join(words[:2]).lower()
    synonym_folder = SYNONYMS.get(synonym_key)
    if synonym_folder is not None and synonym_folder.lower() in folder_by_lower_name:
        return folder_by_lower_name[synonym_folder.lower()], MATCH_SYNONYM

    return None, None


def match_text(species, how):
    """Say how the species was matched, for the log lines; nothing for an exact match.

    Produces:  "", " (from the first two words of '<species>')" or " (the database name
               of '<species>', a synonym in NCBI Taxonomy)".
    Used by:   run() and summary_line().
    """
    if how == MATCH_FIRST_TWO_WORDS:
        return f" (from the first two words of '{species}')"
    if how == MATCH_SYNONYM:
        return f" (the database name of '{species}', a synonym in NCBI Taxonomy)"
    return ""


def same_genus_folders(name, folders):
    """List the database folders of the same genus as a species that was not matched.

    Takes in:  the species name and the species folders.
    Produces:  the folders whose first part is the genus of the name (case ignored),
               for example the four Fusarium folders for "Fusarium oxysporum"; an
               empty list when there is none.
    Used by:   not_screened_note().
    """
    words = species_words(name)
    if len(words) == 0:
        return []
    genus = words[0].lower()
    matches = []
    for folder in folders:
        if folder.split("_")[0].lower() == genus:
            matches.append(folder)
    return matches


def compartment_of(target):
    """Say where a target is encoded: "mitochondrial" for cytochrome b, else "nuclear"."""
    if target in MITOCHONDRIAL_TARGETS:
        return RUN_MITOCHONDRIAL
    return RUN_NUCLEAR


def mitochondrial_code(folder, setting):
    """Work out the genetic code of the mitochondrial run.

    Takes in:  the species folder and --mito-genetic-code ("auto" or a number as text).
    Produces:  3 for the Saccharomycetaceae of YEAST_MITOCHONDRIAL_CODE_FOLDERS and 4
               for every other species when the setting is auto; the number otherwise.
    Used by:   run().
    """
    if setting == "auto":
        if folder in YEAST_MITOCHONDRIAL_CODE_FOLDERS:
            return YEAST_MITOCHONDRIAL_CODE
        return MOLD_MITOCHONDRIAL_CODE
    return int(setting)


# ======================================================================
# THE INPUT FASTA FILES
# ======================================================================

def fasta_has_sequence(path):
    """Tell whether a FASTA file holds at least one sequence line.

    Takes in:  a FASTA path (the genome, or the mitochondrial contigs, which the
               screen writes as an empty file when it found none).
    Does:      reads until the first line that is neither empty nor a header.
    Produces:  True or False.
    Used by:   run(): an empty genome is an error; empty mitochondrial contigs mean
               there is nothing to search for cytochrome b.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fasta:
        for line in fasta:
            text = line.strip()
            if text != "" and not text.startswith(">"):
                return True
    return False


# ======================================================================
# RUNNING CHROQUETAS
# ======================================================================

def chroquetas_command(script, genome, database, folder, out_dir, threads, genetic_code):
    """Build the ChroQueTas command line of one run.

    Takes in:  the script, the input FASTA (the link with its fixed name), the database
               folder, the species folder, the run's output folder, the threads, and
               the genetic code (None: let ChroQueTas choose the nuclear code).
    Does:      puts -f before -s: ChroQueTas checks the species while it reads its
               options, against the database it knows at that moment (by default one
               next to the script, which does not exist here).
    Produces:  the list of words, starting with "bash" and the script.
    Used by:   run().
    """
    command = [
        "bash", script,
        "-g", genome,
        "-f", database,
        "-s", folder,
        "-o", out_dir,
        "-t", str(threads),
    ]
    if genetic_code is not None:
        command.extend(["-c", str(genetic_code)])
    return command


def run_chroquetas(command, folder):
    """Run one ChroQueTas command and capture everything it prints.

    Takes in:  the command and a folder to run it in (a temporary folder, so that a
               stray file could not land among the outputs).
    Does:      runs the command with TERM=dumb (ChroQueTas colours its messages with
               tput, which prints errors when no terminal type is set), and collects
               standard output and standard error together, in order.
    Produces:  (exit status, the printed text, seconds taken). The exit status is
               recorded but not trusted: see check_run_finished().
    Used by:   run().
    """
    environment = dict(os.environ)
    environment["TERM"] = "dumb"
    started = time.monotonic()
    try:
        finished = subprocess.run(
            command, cwd=folder, env=environment, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace",
        )
    except FileNotFoundError:
        raise ValueError("the program bash was not found, so ChroQueTas cannot be run") from None
    seconds = time.monotonic() - started
    return finished.returncode, finished.stdout, seconds


def last_lines(text, count):
    """Return the last `count` non-empty lines of a text, each indented, as one block."""
    lines = []
    for line in text.splitlines():
        if line.strip() != "":
            lines.append("    " + line.rstrip())
    return "\n".join(lines[-count:])


def targets_not_extracted(log_text):
    """Return the targets ChroQueTas says it could not extract, from its printed text.

    Takes in:  what a run printed. For a target miniprot did not find, ChroQueTas
               prints "WARNING: We could not extract <Target> from the genome...".
    Produces:  the set of those target names.
    Used by:   check_run_finished() and read_run_results().
    """
    names = set()
    for line in log_text.splitlines():
        match = NOT_EXTRACTED_LINE.search(line)
        if match is not None:
            names.add(match.group(1))
    return names


def summary_file_path(out_dir, prefix):
    """The path of ChroQueTas's summary file of one run."""
    return os.path.join(out_dir, f"{prefix}.ChroQueTaS.AMR_summary.txt")


def check_run_finished(run_name, out_dir, prefix, log_text, exit_status, targets):
    """Judge whether a ChroQueTas run finished, from its files rather than its exit status.

    Takes in:  the run's name, output folder and file prefix, what it printed, its exit
               status and the targets of the species.
    Does:      ChroQueTas's exit status cannot be trusted: it carries on after a
               failed step, and it exits with 1 when it found no target at all. So:
               the run finished when its folder holds the summary file. It also
               counts as finished, with every target not found, when the log names
               every target as not extracted and miniprot's index of the genome
               (tmp/<prefix>.mpi, which ChroQueTas deletes only at the end of a run
               that went on) exists: without that file miniprot did not even read the
               genome, and "not found" would be a false result.
    Produces:  "finished" or "nothing extracted".
    Used by:   run().
    Stops with an error naming the run, its exit status and its last
    LOG_TAIL_LINES lines in any other case. When the log says nothing was extracted
    but the index is missing, the error says so: miniprot refused the genome or the
    genetic code, and ChroQueTas hides miniprot's own messages.
    """
    if os.path.isfile(summary_file_path(out_dir, prefix)):
        return "finished"

    missing = targets_not_extracted(log_text)
    every_target_missing = set(targets).issubset(missing)
    index_file = os.path.join(out_dir, "tmp", prefix + ".mpi")
    if every_target_missing and NOTHING_EXTRACTED_TEXT in log_text and os.path.isfile(index_file):
        return "nothing extracted"

    # The log says nothing was extracted, but without the index miniprot never read the
    # genome (a genetic code it does not know, or an input it cannot read): say so.
    if every_target_missing and NOTHING_EXTRACTED_TEXT in log_text:
        raise ValueError(
            f"the {run_name} run of ChroQueTas did not finish (exit status {exit_status}): it "
            f"extracted no target, but miniprot built no index of the input ({index_file} is "
            "missing), so the genome was not searched; ChroQueTas hides miniprot's own messages. "
            "Check the genetic code (-c) and the input FASTA. The last "
            f"{LOG_TAIL_LINES} lines of its output:\n" + last_lines(log_text, LOG_TAIL_LINES)
        )

    raise ValueError(
        f"the {run_name} run of ChroQueTas did not finish (exit status {exit_status}): "
        f"{out_dir} holds no {prefix}.ChroQueTaS.AMR_summary.txt, and its output does not say "
        f"that no target was found. The last {LOG_TAIL_LINES} lines of its output:\n"
        + last_lines(log_text, LOG_TAIL_LINES)
    )


# ======================================================================
# READING WHAT CHROQUETAS WROTE
# ======================================================================

def read_lines(path):
    """Read a text file into a list of lines without their line ends; empty lines are dropped."""
    lines = []
    with open(path, "r", encoding="utf-8", errors="replace") as text_file:
        for line in text_file:
            text = line.rstrip("\r\n")
            if text.strip() != "":
                lines.append(text)
    return lines


def read_summary_file(path):
    """Read ChroQueTas's summary of one run: every FungAMR and new mutation it found.

    Takes in:  <prefix>.ChroQueTaS.AMR_summary.txt: either the single line "No AMR
               mutations found in <prefix>", or the header SUMMARY_FILE_COLUMNS and one
               row per mutation (the last cell holds the drugs, or "New mutation").
    Produces:  a list of row dictionaries (column -> text); empty for no mutation.
    Used by:   run(), to check the counts against the per-target tables.
    Stops with a clear error on any other content.
    """
    lines = read_lines(path)
    if len(lines) == 1 and lines[0].startswith(NO_MUTATIONS_TEXT):
        return []
    if len(lines) == 0 or lines[0].split("\t") != SUMMARY_FILE_COLUMNS:
        raise ValueError(f"{path}: not the summary file of ChroQueTas 1.0.1 (unexpected first line)")
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(SUMMARY_FILE_COLUMNS):
            raise ValueError(
                f"{path}: a row has {len(cells)} columns instead of {len(SUMMARY_FILE_COLUMNS)}: {line!r}"
            )
        # Pair each column name with its cell: {"Protein": ..., "Fragment": ..., ...}.
        rows.append(dict(zip(SUMMARY_FILE_COLUMNS, cells)))
    return rows


def read_stats_file(path):
    """Read ChroQueTas's count of FungAMR and new mutations per target, when it wrote one.

    Takes in:  <prefix>.ChroQueTaS.AMR_stats.txt: the header STATS_FILE_COLUMNS and one
               row per target. ChroQueTas writes it only when it found a mutation, and
               then only for the targets it handled from the first mutation on, so a
               missing file or a missing row means no mutation.
    Produces:  {target: (FungAMR mutations, new mutations)}; empty without the file.
    Used by:   run(), to check the counts against the per-target tables.
    Stops with a clear error on another header or a count that is not a number.
    """
    if not os.path.isfile(path):
        return {}
    lines = read_lines(path)
    if len(lines) == 0 or lines[0].split("\t") != STATS_FILE_COLUMNS:
        raise ValueError(f"{path}: not the stats file of ChroQueTas 1.0.1 (unexpected first line)")
    counts = {}
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(STATS_FILE_COLUMNS):
            raise ValueError(f"{path}: cannot read the row {line!r}")
        target, fungamr_count, new_count = cells
        if not fungamr_count.isdigit() or not new_count.isdigit():
            raise ValueError(f"{path}: cannot read the row {line!r}")
        counts[target] = (int(fungamr_count), int(new_count))
    return counts


def read_target_table(path):
    """Read the per-target table of one copy of one target.

    Takes in:  <prefix>.ChroQueTaS.<Target>.<fragment>.tsv: either the single line
               "WARNING: genome contains insertions for <name> and thus AMR are not
               calculated", or the header TARGET_TABLE_COLUMNS and one row per FungAMR
               position.
    Does:      turns ChroQueTas's Result into the call of this script (CALLS) and
               keeps every cell as text, with three corrections (CALLS explains them):
                 a position that is not a whole number points at no residue of the
                 protein, so its row becomes "position not found", without amino
                 acids, whatever ChroQueTas called it;
                 an X found (an unreadable codon) that ChroQueTas called a FungAMR or
                 a new mutation becomes "position not found";
                 the reference amino acid found, but called a FungAMR mutation,
                 becomes "no mutation".
    Produces:  ("warning", []) for the single warning line, or ("compared", rows) with
               one dictionary per row: position, reference, found, call, fungicides,
               and chroquetas_call, the call as ChroQueTas made it (for
               check_counts_agree(), which compares with ChroQueTas's own counts).
    Used by:   read_run_results().
    Stops with a clear error on another header, a short row or a Result it does not
    know.
    """
    lines = read_lines(path)
    if len(lines) == 1 and lines[0].startswith(INSERTIONS_WARNING):
        return "warning", []
    if len(lines) == 0 or lines[0].split("\t") != TARGET_TABLE_COLUMNS:
        raise ValueError(f"{path}: not a per-target table of ChroQueTas 1.0.1 (unexpected first line)")
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(TARGET_TABLE_COLUMNS):
            raise ValueError(
                f"{path}: a row has {len(cells)} columns instead of {len(TARGET_TABLE_COLUMNS)}: {line!r}"
            )
        position, reference, found, result, fungicides = cells
        if result not in CALLS:
            raise ValueError(f"{path}: unknown Result {result!r} in the row {line!r}")
        chroquetas_call = CALLS[result]
        call = chroquetas_call
        if not position.isdigit():
            # The translocation row of the database: ChroQueTas read the first residue
            # of the protein for it. That is not a comparison, so it is not counted as one.
            call = CALL_NOT_FOUND
            reference = ""
            found = ""
        elif found == UNREADABLE_CODON_LETTER and chroquetas_call in (CALL_FUNGAMR, CALL_NEW):
            # A codon with N bases: the amino acid is unknown, so the position was not
            # compared. found_aa keeps the X, so the reader sees why.
            call = CALL_NOT_FOUND
        elif chroquetas_call == CALL_FUNGAMR and found == reference:
            # Hap1 N832N: the reference amino acid, which ChroQueTas matched to a
            # database row that lists it as the new one. No change, so no mutation.
            call = CALL_NONE
        rows.append({
            "position": position,
            "reference": reference,
            "found": found,
            "call": call,
            "fungicides": fungicides,
            "chroquetas_call": chroquetas_call,
        })
    return "compared", rows


def fragment_numbers(out_dir, prefix, target):
    """List the copies (fragments) of a target that miniprot extracted, by their .faa files.

    Takes in:  the run's output folder and prefix, and the target name.
    Does:      looks for <prefix>_<Target>.<n>.faa; n is the fragment number.
    Produces:  the numbers, sorted ([1], or [1, 2] for a target found twice).
    Used by:   read_run_results().
    """
    pattern = re.compile(re.escape(f"{prefix}_{target}.") + r"([0-9]+)\.faa")
    numbers = []
    for file_name in os.listdir(out_dir):
        match = pattern.fullmatch(file_name)
        if match is not None:
            numbers.append(int(match.group(1)))
    return sorted(numbers)


def aligned_reference_stretch(path):
    """Read which stretch of the reference protein BLASTP aligned a found copy to.

    Takes in:  <prefix>_<Target>.<n>.blastp.txt, ChroQueTas's table of the BLASTP
               comparison of one copy with the reference protein: a header line, then
               one row per aligned stretch (usually one). The columns sstart and send
               give the first and the last position of the reference in that stretch.
    Produces:  (first, last): the smallest sstart and the largest send. For example
               (53, 447) for a copy of a 447-amino-acid protein that lacks the first 52
               amino acids, because its gene was cut at a contig end.
    Used by:   correct_deletions_beyond_the_copy().
    Stops with a clear error when the file is missing, lacks the two columns, has no
    row, or a row cannot be read.
    """
    if not os.path.isfile(path):
        raise ValueError(f"the BLASTP file {path} does not exist")
    lines = read_lines(path)
    header = []
    if len(lines) > 0:
        header = lines[0].split("\t")
    if BLASTP_REFERENCE_START not in header or BLASTP_REFERENCE_END not in header or len(lines) < 2:
        raise ValueError(
            f"{path}: not the BLASTP file of ChroQueTas 1.0.1 (a header with the columns "
            f"{BLASTP_REFERENCE_START} and {BLASTP_REFERENCE_END}, then at least one row)"
        )
    start_column = header.index(BLASTP_REFERENCE_START)
    end_column = header.index(BLASTP_REFERENCE_END)

    starts = []
    ends = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(header):
            raise ValueError(f"{path}: cannot read the row {line!r}")
        start = cells[start_column]
        end = cells[end_column]
        if not start.isdigit() or not end.isdigit():
            raise ValueError(f"{path}: cannot read the row {line!r}")
        starts.append(int(start))
        ends.append(int(end))
    return min(starts), max(ends)


def correct_deletions_beyond_the_copy(rows, blastp_path):
    """Undo a FungAMR deletion that only marks where a found copy ends.

    Takes in:  the rows of one compared copy (read_target_table()) and the copy's BLASTP
               file.
    Does:      FungAMR lists 22 deletions of one amino acid, written "-" (CALLS above),
               and ChroQueTas calls any gap at such a position a FungAMR mutation. Inside
               the found protein a gap is a real deletion. Beyond its ends, a gap only
               means that the copy stops early: a gene cut at a contig end. So a gap
               called a FungAMR mutation keeps that call only within the stretch of the
               reference that BLASTP aligned the copy to (aligned_reference_stretch());
               outside it, the call becomes "no alignment". chroquetas_call stays as
               ChroQueTas made it, for check_counts_agree().
    Produces:  nothing; changes the call of those rows in place.
    Used by:   read_run_results().
    """
    first, last = aligned_reference_stretch(blastp_path)
    for row in rows:
        gap_called_deletion = row["call"] == CALL_FUNGAMR and row["found"] == GAP_LETTER
        if not gap_called_deletion or not row["position"].isdigit():
            continue
        position = int(row["position"])
        if position < first or position > last:
            row["call"] = CALL_NO_ALIGNMENT


def read_run_results(out_dir, prefix, targets, not_extracted):
    """Collect, per target, what one finished run found for every copy.

    Takes in:  the run's output folder and file prefix, the targets of the species, and
               the targets the run's output names as not extracted
               (targets_not_extracted()).
    Does:      a target named as not extracted has no copy, whatever files are left
               (ChroQueTas gives up on a target whose first copy is empty, and drops it
               from every later step). For every copy of the other targets
               (fragment_numbers()), decides what ChroQueTas did with it, from the
               files it left:
                 compared          the per-target table holds the positions (a gap
                                   at a deletion FungAMR lists, beyond the stretch
                                   the copy covers, becomes "no alignment":
                                   correct_deletions_beyond_the_copy());
                 insertions        the table holds the warning line and the BLASTP
                                   file <prefix>_<Target>.<n>.blastp.txt exists: the
                                   copy passed the identity and coverage filter, and
                                   the alignment put gaps into the reference;
                 below thresholds  the warning line without the BLASTP file:
                                   ChroQueTas writes the same warning for a copy the
                                   identity and coverage filter dropped, so the copy
                                   was in fact not close enough to the reference;
                 skipped           no per-target table at all. ChroQueTas 1.0.1
                                   removes a target it did not find from its list of
                                   targets with a text replacement that also cuts
                                   that name out of longer target names (Vma1 out of
                                   Vma16, Pdr1 out of Pdr15), and the longer target
                                   is then silently skipped.
    Produces:  {target: list of {"fragment", "kind", "rows"}}; a target without any
               copy has an empty list.
    Used by:   run().
    """
    results = {}
    for target in targets:
        copies = []
        if target in not_extracted:
            results[target] = copies
            continue
        for number in fragment_numbers(out_dir, prefix, target):
            table = os.path.join(out_dir, f"{prefix}.ChroQueTaS.{target}.{number}.tsv")
            blastp = os.path.join(out_dir, f"{prefix}_{target}.{number}.blastp.txt")
            if not os.path.isfile(table):
                copies.append({"fragment": number, "kind": "skipped", "rows": []})
                continue
            kind, rows = read_target_table(table)
            if kind == "compared":
                correct_deletions_beyond_the_copy(rows, blastp)
            if kind == "warning":
                if os.path.isfile(blastp):
                    kind = "insertions"
                else:
                    kind = "below thresholds"
            copies.append({"fragment": number, "kind": kind, "rows": rows})
        results[target] = copies
    return results


def check_counts_agree(run_name, results, summary_rows, stats):
    """Check that the per-target tables, the summary file and the stats file tell the same story.

    Takes in:  the run's name, the results of read_run_results() for every target of
               the species, the rows of the summary file and the stats counts.
    Does:      counts the FungAMR and new mutations per target in the per-target
               tables, and again in the summary file, and compares; then compares the
               stats file with the summary for every target it lists. The tables of
               this script are built from the per-target tables, so a disagreement
               means they were misread, or ChroQueTas misbehaved. The counts use the
               calls as ChroQueTas made them (chroquetas_call), because its summary
               file counts those.
    Used by:   run().
    Stops with a clear error naming the run and the counts when they differ.
    """
    from_tables = {}
    for target, copies in results.items():
        fungamr = 0
        new = 0
        for copy in copies:
            for row in copy["rows"]:
                if row["chroquetas_call"] == CALL_FUNGAMR:
                    fungamr += 1
                elif row["chroquetas_call"] == CALL_NEW:
                    new += 1
        if fungamr > 0 or new > 0:
            from_tables[target] = (fungamr, new)

    from_summary = {}
    for row in summary_rows:
        fungamr, new = from_summary.get(row["Protein"], (0, 0))
        if row["Fungicide_resistance"] == SUMMARY_NEW_MUTATION:
            new += 1
        else:
            fungamr += 1
        from_summary[row["Protein"]] = (fungamr, new)

    if from_tables != from_summary:
        raise ValueError(
            f"the {run_name} run: the per-target tables give the (FungAMR, new) mutations "
            f"{from_tables}, ChroQueTas's summary file {from_summary}"
        )
    for target, counts in stats.items():
        if counts != from_summary.get(target, (0, 0)):
            raise ValueError(
                f"the {run_name} run: ChroQueTas's stats file gives {target} the (FungAMR, new) "
                f"mutations {counts}, its summary file {from_summary.get(target, (0, 0))}"
            )


# ======================================================================
# DRUGS, EVIDENCE LEVELS AND VERDICTS
# ======================================================================

def parse_drugs(text):
    """Split a Fungicides cell into its drugs and their two evidence scores.

    Takes in:  ChroQueTas's text, for example "Carbendazim(7/NA),Iprodione(8/NA)" or
               "Benzimidazole(8/NA),Iprodione(NA/-8)"; "NA", "Unknown" (a new
               mutation) or an empty cell name no drug. The database's position lists
               (<Target>.txt) write the drugs the same way.
    Does:      reads one DRUG_ITEM after the other, from the start of the text.
    Produces:  a list of (drug, positive, negative): positive 1-8 or None, negative
               -1 to -8 or None; an empty list when no drug is named.
    Used by:   the verdict and table functions below, and read_known_resistance_changes().
    Stops with a clear error when the text cannot be read to its end.
    """
    if text in ("", "NA", "Unknown"):
        return []
    drugs = []
    position = 0
    while position < len(text):
        match = DRUG_ITEM.match(text, position)
        if match is None or match.end() == position:
            raise ValueError(f"cannot read the drugs and scores in {text!r}")
        name = match.group(1)
        positive = None if match.group(2) == "NA" else int(match.group(2))
        negative = None if match.group(3) == "NA" else int(match.group(3))
        drugs.append((name, positive, negative))
        position = match.end()
    return drugs


def split_by_evidence(drugs):
    """Sort the drugs of one FungAMR change by the kind of evidence each carries.

    Takes in:  the drugs of one change, as parse_drugs() returns them: (drug, positive,
               negative).
    Does:      the rule of EVIDENCE LEVELS AND VERDICTS at the top: a drug with a
               positive level goes to the resistance list, a drug with a negative level
               to the no-resistance list, and a drug with both (a conflicting drug) to
               both lists, each time with its own level.
    Produces:  (resistance, no_resistance): two lists of (drug, level), in FungAMR's
               order; the levels of the second list keep their minus sign.
    Used by:   position_rows(), summarise_target() and change_log_text().
    """
    resistance = []
    no_resistance = []
    for name, positive, negative in drugs:
        if positive is not None:
            resistance.append((name, positive))
        if negative is not None:
            no_resistance.append((name, negative))
    return resistance, no_resistance


def evidence_words(level):
    """Return the words of an evidence level, the same for 8 and -8 (EVIDENCE_LEVELS)."""
    return EVIDENCE_LEVELS.get(abs(level), "level not described")


def evidence_cell(drug_levels):
    """Write a resistance_to or no_resistance_to cell of amr_positions.tsv.

    Takes in:  one of the two lists of split_by_evidence(): (drug, level) pairs.
    Produces:  "Carbendazim (7: population association study); Iprodione (8: natural
               isolates, not validated further)", or, for a negative level, "Iprodione
               (-8: natural isolates, not validated further)"; EMPTY for an empty list.
    Used by:   position_rows().
    """
    parts = []
    for name, level in drug_levels:
        parts.append(f"{name} ({level}: {evidence_words(level)})")
    if len(parts) == 0:
        return EMPTY
    return DRUG_SEPARATOR.join(parts)


def mutation_verdict(drugs):
    """Give the verdict of a FungAMR change from the evidence of its drugs.

    Takes in:  the drugs of the change, as parse_drugs() returns them.
    Does:      takes the first that applies:
                 a drug with a positive level only       -> resistance
                 a drug with both kinds of level         -> conflicting
                 the drugs have negative levels only     -> no resistance
                 no drug named                           -> unknown
               A drug written NA/NA would carry no evidence and is not counted (the
               database of ChroQueTas 1.0.1 has none).
    Produces:  one of the VERDICT_ texts.
    Used by:   position_verdict().
    """
    has_resistance_only = False
    has_both = False
    has_no_resistance_only = False
    for _, positive, negative in drugs:
        if positive is not None and negative is not None:
            has_both = True
        elif positive is not None:
            has_resistance_only = True
        elif negative is not None:
            has_no_resistance_only = True

    if has_resistance_only:
        return VERDICT_RESISTANCE
    if has_both:
        return VERDICT_CONFLICTING
    if has_no_resistance_only:
        return VERDICT_NO_RESISTANCE
    return VERDICT_UNKNOWN


def position_verdict(call, fungicides):
    """Give the verdict of one position from its call and ChroQueTas's drug text.

    Takes in:  the call of the position (CALLS) and its Fungicides cell.
    Produces:  FungAMR mutation -> the verdict of mutation_verdict();
               new mutation     -> unknown (a change FungAMR does not list);
               no mutation      -> wild type;
               no alignment, position not found -> not compared.
    Used by:   position_rows() and summarise_target().
    """
    if call == CALL_FUNGAMR:
        return mutation_verdict(parse_drugs(fungicides))
    if call == CALL_NEW:
        return VERDICT_UNKNOWN
    if call == CALL_NONE:
        return VERDICT_WILD_TYPE
    return VERDICT_NOT_COMPARED


def target_verdict(position_verdicts):
    """Give the verdict of a target from the verdicts of its positions.

    Takes in:  the verdicts of every position of every compared copy of the target.
    Does:      returns the first of VERDICT_ORDER (resistance, conflicting, unknown, no
               resistance, wild type) found among them; positions that were not
               compared do not count.
    Produces:  that verdict, or EMPTY when no position was compared.
    Used by:   summarise_target().
    """
    for verdict in VERDICT_ORDER:
        if verdict in position_verdicts:
            return verdict
    return EMPTY


def change_log_text(change, verdict, fungicides):
    """Write one FungAMR change for the console line, with its verdict and the levels of its drugs.

    Takes in:  the change ("E198A", or "E198A (copy 2)"), its verdict and its Fungicides
               cell.
    Produces:  "E198A (resistance: Carbendazim 7, Iprodione 8)";
               "E198V (resistance: Benzimidazole 8, Carbendazim 7, Diethofencarb 8;
               no resistance: Iprodione -8)"; "Q21L (no resistance: Fluconazole -1)";
               "I142S (conflicting; resistance: Fluconazole 8; no resistance:
               Fluconazole -7)"; "(unknown)" after the change when no drug is named.
    Used by:   summarise_target(), for summary_line().
    """
    resistance, no_resistance = split_by_evidence(parse_drugs(fungicides))
    parts = []
    if len(resistance) > 0:
        items = []
        for name, level in resistance:
            items.append(f"{name} {level}")
        parts.append("resistance: " + ", ".join(items))
    if len(no_resistance) > 0:
        items = []
        for name, level in no_resistance:
            items.append(f"{name} {level}")
        parts.append("no resistance: " + ", ".join(items))
    details = "; ".join(parts)

    if details == "":
        # No drug named: the verdict alone, "(unknown)".
        inside = verdict
    elif verdict in (VERDICT_RESISTANCE, VERDICT_NO_RESISTANCE):
        # The details already start with the verdict: a resistance change always has a
        # drug with a positive level, a no-resistance change only negative levels.
        inside = details
    else:
        # A conflicting change: the verdict, then both lists.
        inside = f"{verdict}; {details}"
    return f"{change} ({inside})"


# ======================================================================
# THE OUTPUT ROWS
# ======================================================================

def cell(text):
    """Return a table cell: the text, or EMPTY when it is empty."""
    if text is None or text == "":
        return EMPTY
    return text


def mutation_label(row, fragment):
    """Name a change the usual way: "E198A", with " (copy 2)" when it lies on a second copy."""
    label = f"{row['reference']}{row['position']}{row['found']}"
    if fragment > 1:
        label += f" (copy {fragment})"
    return label


def known_changes_text(known_changes, position):
    """Write the known_resistance_changes cell of one position: "E198A, E198G, E198K, E198V", or EMPTY.

    Takes in:  the known changes of the target (read_known_resistance_changes()) and the
               position, as text.
    Used by:   position_rows().
    """
    changes = known_changes.get(position, [])
    if len(changes) == 0:
        return EMPTY
    return ", ".join(changes)


def position_rows(sample, target, compartment, copies, known_changes):
    """Build the amr_positions.tsv rows of one target: every position of every compared copy.

    Takes in:  the sample, the target, its compartment, its copies (read_run_results())
               and the changes the database lists as causing resistance at each of its
               positions (read_known_resistance_changes()).
    Does:      one row per row of each compared copy's table, in ChroQueTas's order,
               with the cells of THE TABLES at the top: the verdict (position_verdict());
               for a FungAMR mutation, resistance_to and no_resistance_to
               (split_by_evidence(), evidence_cell()); for a compared position, the known
               resistance changes there; and ChroQueTas's drug text unchanged, or EMPTY
               when it names no drug (NA; Unknown for a new mutation).
    Produces:  a list of lists of cells, in the order of POSITIONS_COLUMNS.
    Used by:   run().
    """
    rows = []
    for copy in copies:
        if copy["kind"] != "compared":
            continue
        for row in copy["rows"]:
            verdict = position_verdict(row["call"], row["fungicides"])
            drugs = parse_drugs(row["fungicides"])

            # The two drug lists belong to the changes FungAMR lists; every other call
            # leaves them empty.
            resistance_to = EMPTY
            no_resistance_to = EMPTY
            if row["call"] == CALL_FUNGAMR:
                resistance, no_resistance = split_by_evidence(drugs)
                resistance_to = evidence_cell(resistance)
                no_resistance_to = evidence_cell(no_resistance)

            # The known resistance changes, for the positions that were compared.
            known = EMPTY
            if verdict != VERDICT_NOT_COMPARED:
                known = known_changes_text(known_changes, row["position"])

            fungamr_drugs = EMPTY
            if len(drugs) > 0:
                fungamr_drugs = row["fungicides"]

            rows.append([
                sample, target, compartment, str(copy["fragment"]), cell(row["position"]),
                cell(row["reference"]), cell(row["found"]), row["call"], verdict,
                resistance_to, no_resistance_to, known, fungamr_drugs,
            ])
    return rows


def copy_reason(kind):
    """Say in a few words why a copy was not compared, for the note of an evaluated target."""
    if kind == "insertions":
        return "insertions against the reference"
    if kind == "below thresholds":
        return "not above 90 % identity and 75 % coverage"
    return "skipped by ChroQueTas"


def summary_drug_cell(entries):
    """Write a resistance_to or no_resistance_to cell of amr_summary.tsv.

    Takes in:  a list of (drug, change, level), each once, in order of appearance.
    Produces:  "Carbendazim (E198A, 7); Iprodione (E198A, 8)", or "Iprodione (E198V, -8)";
               EMPTY for an empty list.
    Used by:   summarise_target().
    """
    parts = []
    for name, change, level in entries:
        parts.append(f"{name} ({change}, {level})")
    if len(parts) == 0:
        return EMPTY
    return DRUG_SEPARATOR.join(parts)


def summarise_target(copies, no_mito_contig):
    """Give one target its status, verdict, mutations, drugs and note for amr_summary.tsv.

    Takes in:  the target's copies from the run of its compartment (an empty list
               when it had none, or when that run was not done), and no_mito_contig:
               True when the run was not done because the screen set no
               mitochondrial contig aside.
    Does:      applies the statuses of WHAT EACH TARGET GETS at the top:
                 a compared copy               -> evaluated, or not evaluated
                                                  (NOTE_NOTHING_COMPARED) when none
                                                  of its positions could be compared
                 else a copy with insertions   -> not evaluated (NOTE_INSERTIONS)
                 else a copy ChroQueTas skipped -> not evaluated (NOTE_SKIPPED)
                 else (no copy, or every copy below the thresholds)
                                               -> not found (NOTE_NOT_FOUND)
               For a target with a compared copy it gives the verdict of the target
               (target_verdict()); lists the FungAMR changes with their verdicts and the
               new mutations (fragment > 1 noted as "(copy n)"); lists the drugs of the
               FungAMR changes with the change and the level, in resistance_to and
               no_resistance_to (a change found on two copies names its drugs once; a
               drug named by two changes appears once for each); counts the positions
               compared; and notes the copies or positions that were not compared.
    Produces:  a dictionary with status, positions_checked, verdict, fungamr_mutations,
               new_mutations, resistance_to, no_resistance_to, note, and two lists for
               the console line: fungamr_log (change_log_text() of each FungAMR change)
               and new.
    Used by:   run().
    """
    summary = {
        "positions_checked": 0, "verdict": EMPTY, "fungamr_log": [], "new": [],
        "fungamr_mutations": EMPTY, "new_mutations": EMPTY,
        "resistance_to": EMPTY, "no_resistance_to": EMPTY,
    }
    if no_mito_contig:
        summary["status"] = STATUS_NOT_EVALUATED
        summary["note"] = NOTE_NO_MITO_CONTIG
        return summary

    kinds = []
    for copy in copies:
        kinds.append(copy["kind"])
    if "compared" not in kinds:
        if "insertions" in kinds:
            summary["status"] = STATUS_NOT_EVALUATED
            summary["note"] = NOTE_INSERTIONS
        elif "skipped" in kinds:
            summary["status"] = STATUS_NOT_EVALUATED
            summary["note"] = NOTE_SKIPPED
        else:
            summary["status"] = STATUS_NOT_FOUND
            summary["note"] = NOTE_NOT_FOUND
        return summary

    summary["status"] = STATUS_EVALUATED
    compared_positions = set()
    not_compared_positions = []
    position_verdicts = []
    fungamr_changes = []          # "E198A (resistance)", one per change and copy
    resistance_entries = []       # (drug, change, level), each once
    no_resistance_entries = []    # (drug, change, negative level), each once
    notes = []
    compared_copies = []
    for copy in copies:
        if copy["kind"] != "compared":
            notes.append(f"copy {copy['fragment']} not compared: {copy_reason(copy['kind'])}")
            continue
        compared_copies.append(str(copy["fragment"]))
        for row in copy["rows"]:
            verdict = position_verdict(row["call"], row["fungicides"])
            position_verdicts.append(verdict)
            if row["call"] in COMPARED_CALLS:
                compared_positions.add(row["position"])
            else:
                # For example "351 (no alignment)", or "351 (no alignment, copy 2)".
                where = row["call"]
                if copy["fragment"] > 1:
                    where += f", copy {copy['fragment']}"
                not_compared_positions.append(f"{cell(row['position'])} ({where})")

            if row["call"] == CALL_FUNGAMR:
                label = mutation_label(row, copy["fragment"])
                fungamr_changes.append(f"{label} ({verdict})")
                summary["fungamr_log"].append(change_log_text(label, verdict, row["fungicides"]))
                # The drugs, named with the change as such ("E198A", without the copy):
                # the same change on a second copy adds nothing new.
                change = mutation_label(row, 1)
                resistance, no_resistance = split_by_evidence(parse_drugs(row["fungicides"]))
                for name, level in resistance:
                    entry = (name, change, level)
                    if entry not in resistance_entries:
                        resistance_entries.append(entry)
                for name, level in no_resistance:
                    entry = (name, change, level)
                    if entry not in no_resistance_entries:
                        no_resistance_entries.append(entry)
            elif row["call"] == CALL_NEW:
                summary["new"].append(mutation_label(row, copy["fragment"]))

    summary["positions_checked"] = len(compared_positions)
    summary["verdict"] = target_verdict(position_verdicts)
    if len(fungamr_changes) > 0:
        summary["fungamr_mutations"] = ", ".join(fungamr_changes)
    if len(summary["new"]) > 0:
        summary["new_mutations"] = ", ".join(summary["new"])
    summary["resistance_to"] = summary_drug_cell(resistance_entries)
    summary["no_resistance_to"] = summary_drug_cell(no_resistance_entries)

    if len(compared_copies) > 1:
        notes.insert(0, (
            f"found {len(compared_copies)} times (copies {', '.join(compared_copies)}): a "
            "duplicated gene, or a region the assembly holds more than once"
        ))
    if len(not_compared_positions) > 0:
        notes.append("positions not compared: " + ", ".join(not_compared_positions))

    # Copies were read position by position, but no FungAMR position fell in a part of
    # the found protein that could be compared (every one lies in a gap, as at a gene cut
    # at a contig end or a premature stop codon, at a codon that could not be read, or
    # outside the alignment). Nothing was
    # checked, so the target is not evaluated (and its verdict is EMPTY); the notes list
    # the positions.
    if len(compared_positions) == 0:
        summary["status"] = STATUS_NOT_EVALUATED
        notes.insert(0, NOTE_NOTHING_COMPARED)

    summary["note"] = "; ".join(notes) if len(notes) > 0 else EMPTY
    return summary


def summary_cells(sample, species, target, compartment, summary):
    """Put one target's summary into the cells of amr_summary.tsv (order of SUMMARY_COLUMNS)."""
    return [
        sample, species, target, compartment, summary["status"], str(summary["positions_checked"]),
        summary["verdict"], summary["fungamr_mutations"], summary["new_mutations"],
        summary["resistance_to"], summary["no_resistance_to"], summary["note"],
    ]


def not_screened_note(species, folders):
    """Write the note of a sample whose species is not in the database.

    Takes in:  the species name and the species folders of the database.
    Produces:  "the species 'X' is not among the 57 species of the ChroQueTas
               database", followed by the database folders of the same genus when
               there are any.
    Used by:   run(), for the summary row and the NOTE line.
    """
    note = (
        f"the species '{species}' is not among the {len(folders)} species of the ChroQueTas database"
    )
    relatives = same_genus_folders(species, folders)
    if len(relatives) > 0:
        note += f"; the database holds these of the same genus: {', '.join(relatives)}"
    return note


# ======================================================================
# WRITE THE OUTPUT FILES
# ======================================================================

def write_table(path, columns, rows):
    """Write a tab-separated table: the header line, then one line per row of cells."""
    with open(path, "w", encoding="utf-8") as table:
        table.write("\t".join(columns) + "\n")
        for row in rows:
            table.write("\t".join(row) + "\n")


def mqc_header_lines():
    """The header lines that make MultiQC pick the summary up as custom content.

    Does:      builds the five lines: the id, the section title, a description shown
               under it (what was screened, what the verdicts mean and how to read the
               evidence levels), the plot type and the table's namespace. The
               description sits between single quotes, so it must not hold one itself.
               MultiQC finds every *_mqc.tsv file below the output folder on its own.
    Produces:  a list of five lines, without line breaks.
    Used by:   write_mqc_table().
    """
    return [
        "# id: amr_chroquetas",
        "# section_name: 'AMR mutations (ChroQueTas)'",
        "# description: 'Mutations of the FungAMR database in the drug target proteins of "
        "the species, screened with ChroQueTas (Bédard et al. 2025); a target is evaluated, "
        "not evaluated, not found, or the sample is not screened. The verdict of a target is "
        "the first that applies: resistance (a change with, for at least one drug, evidence "
        "of resistance and none against it), conflicting (a change with evidence both for and "
        "against resistance to the same drug), unknown (a change FungAMR does not list, or "
        "lists without drugs), no resistance (only reports that found no resistance, to the "
        "drugs named), wild type (the reference amino acid at every compared position). Each "
        "drug carries its evidence level, "
        "from 1 (strongest: the change made in the species and measured) to 8 (weakest: "
        "natural isolates, not validated further); a minus sign marks evidence of no "
        "resistance.'",
        "# plot_type: 'table'",
        "# pconfig: {namespace: 'AMR mutations'}",
    ]


def write_mqc_table(path, rows, sample):
    """Write the MultiQC copy of the summary.

    Takes in:  the summary rows (lists of cells) and the sample name.
    Does:      writes the lines of mqc_header_lines(), then the header with a first
               column for the row name, then every row named "<sample>: <target>".
               MultiQC takes the first column as the name of the row, so the names
               must differ between samples for the tables of a batch to merge into one.
    Used by:   run(). Rule multiqc lists the file among its inputs.
    """
    with open(path, "w", encoding="utf-8") as table:
        for line in mqc_header_lines():
            table.write(line + "\n")
        table.write("\t".join(["Sample"] + SUMMARY_COLUMNS) + "\n")
        for row in rows:
            target = row[SUMMARY_COLUMNS.index("target")]
            table.write("\t".join([f"{sample}: {target}"] + row) + "\n")


# ======================================================================
# THE LINES FOR THE LOG AND THE CONSOLE
# ======================================================================

def summary_line(sample, folder, how_text, targets, summaries):
    """Write the one-line summary that the rule copies to the console.

    Takes in:  the sample, the species folder, how it was matched (match_text()), the
               targets in order and their summaries.
    Does:      counts and names the targets per status, then per verdict (only the
               verdicts that occur, in the order of VERDICT_ORDER), then gives every
               FungAMR change with its verdict and the levels of its drugs
               (change_log_text()), and the new mutations.
    Produces:  "[FunFlux] amr_chroquetas <sample>: species folder Botrytis_cinerea; 4
               targets: 3 evaluated (Beta-tubulin, Bos1, OS1), 1 not evaluated
               (Cytochrome_b), 0 not found; verdicts: 1 resistance (Beta-tubulin),
               2 wild type (Bos1, OS1); FungAMR mutations: Beta-tubulin E198A
               (resistance: Carbendazim 7, Iprodione 8); new mutations: none".
    Used by:   run().
    """
    by_status = {STATUS_EVALUATED: [], STATUS_NOT_EVALUATED: [], STATUS_NOT_FOUND: []}
    by_verdict = {}
    for verdict in VERDICT_ORDER:
        by_verdict[verdict] = []
    fungamr_parts = []
    new_parts = []
    for target in targets:
        summary = summaries[target]
        by_status[summary["status"]].append(target)
        if summary["verdict"] in by_verdict:
            by_verdict[summary["verdict"]].append(target)
        for change_text in summary["fungamr_log"]:
            fungamr_parts.append(f"{target} {change_text}")
        if len(summary["new"]) > 0:
            new_parts.append(f"{target} {', '.join(summary['new'])}")

    status_parts = []
    for status, names in by_status.items():
        text = f"{len(names)} {status}"
        if len(names) > 0:
            text += f" ({', '.join(names)})"
        status_parts.append(text)

    # The tally of the verdicts: only those that occur; "none" when no target was
    # evaluated.
    verdict_parts = []
    for verdict, names in by_verdict.items():
        if len(names) > 0:
            verdict_parts.append(f"{len(names)} {verdict} ({', '.join(names)})")
    verdict_text = ", ".join(verdict_parts) if len(verdict_parts) > 0 else "none"

    fungamr_text = "; ".join(fungamr_parts) if len(fungamr_parts) > 0 else "none"
    new_text = "; ".join(new_parts) if len(new_parts) > 0 else "none"
    targets_text = "1 target" if len(targets) == 1 else f"{len(targets)} targets"
    return (
        f"[FunFlux] {MESSAGE_NAME} {sample}: species folder {folder}{how_text}; {targets_text}: "
        f"{', '.join(status_parts)}; verdicts: {verdict_text}; FungAMR mutations: {fungamr_text}; "
        f"new mutations: {new_text}"
    )


# ======================================================================
# PUTTING THE STEPS TOGETHER
# ======================================================================

def check_no_space(label, path):
    """Stop on a path with a space or a tab: ChroQueTas does not quote its paths."""
    if re.search(r"\s", os.path.abspath(path)):
        raise ValueError(
            f"the {label} {os.path.abspath(path)!r} contains a space; ChroQueTas 1.0.1 cannot "
            "work with such a path: use a folder whose path has no space"
        )


def write_not_screened(args, folders):
    """Write the three tables of a sample whose species is not in the database, and the NOTE line.

    Takes in:  the options and the species folders of the database.
    Does:      one summary row with target "-", status not screened, verdict "-" and
               the note of not_screened_note(); a positions table with the header
               alone; the MultiQC copy. ChroQueTas is not run, and --work-dir stays empty.
    Used by:   run().
    """
    note = not_screened_note(args.species, folders)
    summary = {
        "status": STATUS_NOT_SCREENED, "positions_checked": 0, "verdict": EMPTY,
        "fungamr_mutations": EMPTY, "new_mutations": EMPTY, "resistance_to": EMPTY,
        "no_resistance_to": EMPTY, "note": note,
    }
    summary_rows = [summary_cells(args.sample, args.species, EMPTY, EMPTY, summary)]
    write_table(args.out_positions, POSITIONS_COLUMNS, [])
    write_table(args.out_summary, SUMMARY_COLUMNS, summary_rows)
    write_mqc_table(args.out_mqc, summary_rows, args.sample)
    print(f"[FunFlux] NOTE {MESSAGE_NAME} {args.sample}: not screened: {note}", file=sys.stderr)


def plan_runs(targets, folder, args):
    """Decide which ChroQueTas runs a species needs, on which FASTA, with which genetic code.

    Takes in:  the targets of the species, its folder, and the options (--genome,
               --mito-contigs, --mito-genetic-code).
    Does:      TWO RUNS at the top: a nuclear run when the species has a nuclear
               target; a mitochondrial run when it has Cytochrome_b, on the
               mitochondrial contigs, or on the genome when the run does not split
               the mitochondrion off (--mito-contigs empty). When the screen ran and set
               no mitochondrial contig aside there is nothing to search, and the
               mitochondrial run is left out.
    Produces:  (runs, mito_skipped): a list of (run name, input FASTA, genetic code,
               None for ChroQueTas's own nuclear choice), and True when the
               mitochondrial run was left out for want of a mitochondrial contig.
    Used by:   run().
    """
    has_nuclear_target = False
    has_mitochondrial_target = False
    for target in targets:
        if compartment_of(target) == RUN_MITOCHONDRIAL:
            has_mitochondrial_target = True
        else:
            has_nuclear_target = True

    runs = []
    mito_skipped = False
    if has_nuclear_target:
        runs.append((RUN_NUCLEAR, args.genome, None))
    if has_mitochondrial_target:
        code = mitochondrial_code(folder, args.mito_genetic_code)
        if args.mito_contigs == "":
            # Nothing split off in this run: the mitochondrion is among the contigs.
            runs.append((RUN_MITOCHONDRIAL, args.genome, code))
        elif fasta_has_sequence(args.mito_contigs):
            runs.append((RUN_MITOCHONDRIAL, args.mito_contigs, code))
        else:
            mito_skipped = True
    return runs, mito_skipped


def do_run(run_name, input_fasta, code, folder, targets, temporary, args):
    """Run ChroQueTas once and read what it found for every target of the species.

    Takes in:  the run's name, its input FASTA and genetic code (from plan_runs()), the
               species folder and its targets, the temporary folder, and the options.
    Does:      links the input into the temporary folder under its fixed name (see
               RUN_INPUT_NAMES), runs ChroQueTas into <--work-dir>/<run name>, copies
               the command and ChroQueTas's output to standard error (the rule's log),
               judges the run by its files (check_run_finished()), and reads the copies
               of every target, checking the counts against ChroQueTas's own summary
               and stats files. After a run that found nothing, every target has no copy.
    Produces:  {target: copies}, for every target of the species; run() keeps those of
               this run's compartment.
    Used by:   run().
    """
    link = os.path.join(temporary, RUN_INPUT_NAMES[run_name])
    os.symlink(os.path.abspath(input_fasta), link)
    prefix = RUN_INPUT_NAMES[run_name].split(".")[0]
    out_dir = os.path.join(os.path.abspath(args.work_dir), run_name)
    command = chroquetas_command(
        os.path.abspath(args.chroquetas_script), link, os.path.abspath(args.database),
        folder, out_dir, args.threads, code,
    )
    print(
        f"chroquetas_amr.py: {run_name} run on {input_fasta} (as {RUN_INPUT_NAMES[run_name]}): "
        f"{' '.join(command)}",
        file=sys.stderr,
    )
    exit_status, log_text, seconds = run_chroquetas(command, temporary)
    print(f"----- ChroQueTas output, {run_name} run -----", file=sys.stderr)
    print(log_text.rstrip("\n"), file=sys.stderr)
    print(
        f"----- end of the {run_name} run: exit status {exit_status}, {seconds:.1f} s -----",
        file=sys.stderr,
    )

    outcome = check_run_finished(run_name, out_dir, prefix, log_text, exit_status, targets)
    if outcome == "nothing extracted":
        results = {}
        for target in targets:
            results[target] = []
        return results

    results = read_run_results(out_dir, prefix, targets, targets_not_extracted(log_text))
    chroquetas_summary = read_summary_file(summary_file_path(out_dir, prefix))
    stats = read_stats_file(os.path.join(out_dir, f"{prefix}.ChroQueTaS.AMR_stats.txt"))
    check_counts_agree(run_name, results, chroquetas_summary, stats)
    return results


def run(args):
    """Do the whole job: the genome in, ChroQueTas's folders and the three tables out.

    Takes in:  args, the parsed and checked command-line options.
    Does:      checks the inputs; matches the species (a species not in the database
               ends here, write_not_screened()); reads the known resistance changes of
               every target from the database (read_known_resistance_changes()); plans
               the runs (plan_runs()); runs ChroQueTas for each (do_run()) and keeps,
               for every target, the result of the run of its compartment; works out
               the verdicts and writes the three tables last, so that a failure never
               leaves a half-written table behind.
    Produces:  the folders, the tables and the lines on standard error. Raises OSError
               or ValueError on bad input or a failed run; main() turns that into a
               message.
    Used by:   main().
    """
    started = time.monotonic()

    # The inputs: ChroQueTas, the genome, the mitochondrial contigs, the database.
    for label, path in (("ChroQueTas script", args.chroquetas_script), ("genome", args.genome)):
        if not os.path.isfile(path):
            raise ValueError(f"the {label} {path} does not exist")
    if not fasta_has_sequence(args.genome):
        raise ValueError(f"the genome {args.genome} holds no sequence")
    if args.mito_contigs != "" and not os.path.isfile(args.mito_contigs):
        raise ValueError(f"the mitochondrial contigs file {args.mito_contigs} does not exist")
    folders = list_species_folders(args.database)

    # The folder this script fills. ChroQueTas refuses an output folder that already
    # exists, and an old result must never be read as a new one.
    os.makedirs(args.work_dir, exist_ok=True)
    for run_name in (RUN_NUCLEAR, RUN_MITOCHONDRIAL):
        old_folder = os.path.join(args.work_dir, run_name)
        if os.path.exists(old_folder):
            raise ValueError(f"{old_folder} already exists: remove it, or give another --work-dir")

    # The species, its targets and the runs they need.
    folder, how = resolve_species(args.species, folders)
    if folder is None:
        write_not_screened(args, folders)
        return
    targets = list_targets(args.database, folder)

    # The changes the database lists as causing resistance at each position of each
    # target, for the known_resistance_changes column. Read before ChroQueTas runs, so
    # that a position list this script cannot read stops the job at once.
    known_by_target = {}
    for target in targets:
        position_list = os.path.join(args.database, folder, target + ".txt")
        known_by_target[target] = read_known_resistance_changes(position_list)

    planned_runs, mito_skipped = plan_runs(targets, folder, args)
    how_text = match_text(args.species, how)
    print(
        f"chroquetas_amr.py: sample {args.sample}: species {args.species!r} -> database folder "
        f"{folder}{how_text}; targets: {', '.join(targets)}",
        file=sys.stderr,
    )
    check_no_space("output folder", args.work_dir)
    check_no_space("database folder", args.database)

    # Run ChroQueTas for each planned run. Every target keeps the result of the run of
    # its compartment; a target found by the other run is ignored.
    copies_by_target = {}
    with tempfile.TemporaryDirectory(prefix="chroquetas_amr_") as temporary:
        check_no_space("temporary folder", temporary)
        for run_name, input_fasta, code in planned_runs:
            results = do_run(run_name, input_fasta, code, folder, targets, temporary, args)
            for target in targets:
                if compartment_of(target) == run_name:
                    copies_by_target[target] = results[target]

    # The rows of the three tables, target by target in alphabetical order.
    positions = []
    summary_rows = []
    summaries = {}
    for target in targets:
        compartment = compartment_of(target)
        no_mito_contig = mito_skipped and compartment == RUN_MITOCHONDRIAL
        copies = copies_by_target.get(target, [])
        summary = summarise_target(copies, no_mito_contig)
        summaries[target] = summary
        positions.extend(position_rows(args.sample, target, compartment, copies, known_by_target[target]))
        summary_rows.append(summary_cells(args.sample, folder, target, compartment, summary))

    write_table(args.out_positions, POSITIONS_COLUMNS, positions)
    write_table(args.out_summary, SUMMARY_COLUMNS, summary_rows)
    write_mqc_table(args.out_mqc, summary_rows, args.sample)

    seconds = time.monotonic() - started
    print(
        f"chroquetas_amr.py: sample {args.sample}: {len(planned_runs)} ChroQueTas run(s), "
        f"{len(positions)} positions in the positions table; {seconds:.1f} s.",
        file=sys.stderr,
    )
    console_line = summary_line(args.sample, folder, how_text, targets, summaries)
    print(console_line, file=sys.stderr)


# ======================================================================
# COMMAND LINE
# ======================================================================

def parse_arguments(argv):
    """Read the options given on the command line.

    Takes in:  argv, the list of words after the script name (None means "use the real
               command line"; the tests pass their own list).
    Does:      defines the options with argparse (Python's built-in reader of
               command-line options). --mito-contigs may be given without a value:
               Snakemake's {...:q} quoting writes an empty value as nothing at all, and
               an empty value means "this run does not split the mitochondrion off".
    Produces:  an object whose attributes are the options (dashes become underscores:
               args.mito_contigs, args.out_summary, ...).
    Used by:   main().
    """
    parser = argparse.ArgumentParser(
        prog="chroquetas_amr.py",
        description=(
            "Screen a fungal genome for the antifungal resistance mutations of FungAMR with "
            "ChroQueTas: a nuclear run and, for cytochrome b, a mitochondrial run, then one "
            "summary row per target and one row per position checked."
        ),
    )
    parser.add_argument("--chroquetas-script", required=True, help="ChroQueTas.sh, as downloaded")
    parser.add_argument(
        "--database", required=True,
        help="the folder that directly holds the species folders of FungAMR",
    )
    parser.add_argument("--genome", required=True, help="the nuclear genome (FASTA)")
    parser.add_argument(
        "--mito-contigs", nargs="?", const="", default="",
        help=(
            "the mitochondrial contigs (FASTA, possibly empty); empty when the run does not "
            "split the mitochondrion off"
        ),
    )
    parser.add_argument("--species", required=True, help="the species name of the annotation table")
    parser.add_argument(
        "--mito-genetic-code", default="auto",
        help=(
            f"auto, or the NCBI genetic code ({NCBI_GENETIC_CODES_TEXT}) of the mitochondrial "
            "run (default: %(default)s)"
        ),
    )
    parser.add_argument(
        "--threads", type=int, default=DEFAULT_THREADS,
        help="threads for miniprot and MAFFT (default: %(default)s)",
    )
    parser.add_argument("--sample", required=True, help="the sample name, for the tables and messages")
    parser.add_argument(
        "--work-dir", required=True,
        help="the chroquetas/ folder for ChroQueTas's own output",
    )
    parser.add_argument("--out-summary", required=True, help="output: amr_summary.tsv")
    parser.add_argument("--out-positions", required=True, help="output: amr_positions.tsv")
    parser.add_argument("--out-mqc", required=True, help="output: amr_summary_mqc.tsv")
    return parser.parse_args(argv)


def check_options(args):
    """Hold the settings to the limits the workflow applies to the configuration.

    Takes in:  the parsed options; the genetic code is normalised in place ("Auto"
               becomes "auto", " 4 " becomes "4").
    Does:      stops with a ValueError, which main() turns into the error line, when
               the genetic code is neither auto nor one of NCBI_GENETIC_CODES, or when
               the threads are below 1.
    Used by:   main().
    """
    code = args.mito_genetic_code.strip()
    if code.lower() == "auto":
        args.mito_genetic_code = "auto"
    elif code.isdigit() and int(code) in NCBI_GENETIC_CODES:
        args.mito_genetic_code = str(int(code))
    else:
        raise ValueError(
            "--mito-genetic-code must be auto or the number of an NCBI genetic code: "
            f"{NCBI_GENETIC_CODES_TEXT} (got: {args.mito_genetic_code!r})"
        )
    if args.threads < 1:
        raise ValueError(f"--threads must be at least 1 (got: {args.threads})")


def main(argv=None):
    """Run the script from the command line.

    Takes in:  argv, the command-line words (None = the real command line).
    Does:      parses and checks the options and calls run(). A problem with an input,
               a setting or a ChroQueTas run becomes a readable message on standard
               error instead of a Python traceback, and exit status 1, so Snakemake
               sees that the step failed.
    Produces:  exit status 0 on success, 1 on a problem, 2 on a command line argparse
               cannot read (argparse's own code).
    Used by:   the last lines of this file when rule amr_chroquetas runs the script,
               and by the tests, which pass their own argv.
    """
    args = parse_arguments(argv)
    try:
        check_options(args)
        run(args)
    except (OSError, ValueError) as problem:
        # OSError: a file cannot be opened or written (for example it does not exist).
        # ValueError: a setting outside its limits, a database without species folders,
        # a run that did not finish, or a ChroQueTas file that does not look like the
        # files of ChroQueTas 1.0.1.
        print(f"[FunFlux] {MESSAGE_NAME}: error: {problem}", file=sys.stderr)
        sys.exit(1)


# Run main() only when this file is run as a script, not when the tests import it.
if __name__ == "__main__":
    main()
