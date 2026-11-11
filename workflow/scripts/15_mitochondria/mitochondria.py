#!/usr/bin/env python3
"""Find the mitochondrial genome among assembled contigs and set it apart.

Version 0.2 of FunFlux's mitochondrial screen, also used by the unpublished
FunFluxL+, with three deliberate differences marked "CHANGED in FunFlux v2" in
the code, and two more, the stricter title check and the species check, marked
"CHANGED in FunFlux 2.1.0". It still accepts every command line of version 0.1
of the script.

Four subcommands; the first word after the script name chooses the job.
screen: decides which contigs are mitochondrial, from the workflow's BLAST
table (hits against the NCBI nucleotide database, core_nt by default) and, in
the long-read modes, the assembler's evidence table. It writes a report and the
list of confirmed contig names.
split: writes the nuclear contigs and the mitochondrial contigs to two separate
FASTA files; the input assembly is not changed.
evidence-flye and evidence-hifiasm: turn the assembler's own output into that
per-contig evidence table (read depth, circular or linear).
A contig is never called mitochondrial from circularity or read depth alone, and
a contig whose best mitochondrial match is an animal or a plant (human
mitochondrial DNA from the library, for example) is reported as foreign, not set
apart as the sample's mitochondrion.
"""

# =============================================================================
# WHY THIS SCRIPT EXISTS, AND HOW THE DATA FLOWS THROUGH IT
# =============================================================================
# A fungal cell carries two genomes: the nuclear genome (the chromosomes, tens of
# megabases) and a small mitochondrial genome (usually some tens of kb) that
# exists in many copies per cell. An assembler
# cannot tell them apart: it delivers the mitochondrial genome as one more
# ordinary contig. If that contig stayed in the assembly, funannotate would run
# the nuclear gene predictor on it, which is wrong because fungal mitochondria
# translate with a different genetic code (code 4 for the filamentous fungi
# studied here: UGA means tryptophan, not "stop"). So this script finds the
# mitochondrial contig, and the workflow sets it aside before assembly QC and the
# nuclear annotation.
#
# How the pieces fit together. The rules named here are in
# workflow/rules/shared/15_mitochondria.smk and 10_decontam.smk.
#
#   Flye assembly_info.txt ---> evidence-flye ------\
#   hifiasm primary GFA ------> evidence-hifiasm ----+--> assembly_evidence.tsv
#     (long-read modes only; see the note below)                  |
#                                                                 v
#   blastout, the BLAST table (rule blast_contigs) -----------> screen
#   draft contigs, the FASTA that BLAST ran on ---------------> (rule screen_mitochondria)
#   nodes.dmp, the NCBI taxonomy (BLAST database folder) -----> (the species check)
#                                                                 |
#        mitochondrial_screen.tsv     <---------------------------+  the report, for people
#        confirmed_mitochondrial.ids  <---------------------------+  one contig name per line
#             |
#             +--> rule select_contigs keeps these contigs, whatever the
#             |    genus filter of the decontamination says
#             v
#   contigs_all.fasta (the finished assembly) + the list -----> split
#                                                              (rule separate_mitochondria)
#        contigs_final.fasta          the nuclear genome -> QUAST, BUSCO, the ITS
#                                     search and funannotate
#        contigs_mitochondrial.fasta  the mitochondrial contig(s), a deliverable
#                                     (empty when none); only the optional AMR
#                                     screen reads it, for cytochrome b
#
# The two evidence subcommands serve the long-read modes: rule assembly_evidence
# (rules/longread/20_assembly.smk, and its copy rules/hybrid/40_long_assembly.smk)
# runs evidence-hifiasm or evidence-flye once the assembler has finished and writes
# assembly_evidence.tsv next to assembly.fasta; rule screen_mitochondria
# (rules/shared/15_mitochondria.smk) then passes that file as --assembly-evidence.
# In contigs and illumina mode there is no assembler evidence and the rule passes
# an empty value: the decision rests on BLAST alone.
#
# Words used all through this file:
#   hit                   one line of the BLAST table: one local alignment
#                         between a contig and a reference sequence. A reference
#                         that aligns to the contig in several pieces gives
#                         several lines.
#   qcovs                 "query coverage per subject": the percentage of the
#                         CONTIG that BLAST aligned to one particular reference
#                         sequence: the contig positions covered by any saved
#                         alignment to that reference, each position counted
#                         once, so overlapping alignments are not added up.
#                         Every line of that reference shows the same number.
#                         100 means the whole contig matches that reference.
#   direct organelle hit  a hit whose reference title says the reference IS a
#                         mitochondrial genome (see is_direct_mito_title()).
#
# This script uses only the Python standard library. Every rule that runs it
# does so in workflow/envs/python.yaml, a conda environment that holds a pinned
# Python 3.12 and no other package (no Biopython, no pandas). It is written for
# Python 3.8 or newer, so it can also be run by hand with any Python from 3.8 on.

# Standard-library modules only.
import argparse                      # reads the command line (the subcommands and their options)
import csv                           # reads and writes the tab-separated tables
import os                            # makes output folders
import re                            # the pattern that recognises organelle titles
import statistics                    # the median read depth
import sys                           # error messages and exit status
from collections import OrderedDict  # FASTA records that keep the order of the input file


# ── Fixed settings: the organelle title check, the modes, two cut-offs ───────

# Recognise BLAST reference titles that describe an organelle genome itself.
# The title is the last BLAST column (stitle), free text typed by whoever
# deposited the sequence in NCBI. The decision is made in is_direct_mito_title()
# below, with two regular expressions (text patterns), each searched for anywhere
# in the title, upper or lower case alike. The match is on the wording of the
# title, not on the species.
#
# CHANGED in FunFlux 2.1.0: the check is stricter. FunFlux 2.0.0 and FunFluxL+
# used one pattern, with the first two wordings below only. The hyphen rule, the
# "organelle: mitochondrion" wording and NUCLEAR_GENE_RECORD are new: a nuclear
# gene named after mitochondrial DNA no longer counts as a mitochondrial genome,
# and the mitochondrial sequence of a genome assembly now does.
#
# MITO_GENOME_WORDING: the wordings that say "this sequence IS mitochondrial DNA":
#   - "mitochondrion," followed by "complete" or "partial" and then "genome" or
#     "sequence": "Fusarium graminearum mitochondrion, complete genome";
#   - "mitochondrial DNA" or "mitochondrial genome": "Homo sapiens mitochondrial
#     DNA, complete genome". Not when a hyphen follows, as in "mitochondrial
#     DNA-directed RNA polymerase", the name of a nuclear gene;
#   - "organelle: mitochondrion", the wording NCBI uses for the mitochondrial
#     sequence of a genome assembly: "Cordyceps fumosorosea genome assembly,
#     organelle: mitochondrion".
# Do not match the much broader word "mitochondrial": nuclear-encoded proteins
# routinely use it in their descriptions (for example "mitochondrial ribosomal
# protein L5 mRNA"), and a contig that matches such a gene is nuclear and must
# stay in the nuclear genome.
MITO_GENOME_WORDING = re.compile(
    r"\bmitochondrion,\s*(?:complete|partial)\s+(?:genome|sequence)\b"
    r"|\bmitochondrial\s+(?:DNA|genome)\b(?!-)"
    r"|\borganelle:\s*mitochondrion\b",
    re.IGNORECASE,
)

# NUCLEAR_GENE_RECORD: the words "mRNA" and "protein", which mark the record of a
# nuclear gene, its transcript or its protein. Some such genes carry one of the
# wordings above in their names, for example "Increased loss of mitochondrial DNA
# protein 1 (ISF_05954), partial mRNA" or "mitochondrial genome maintenance protein";
# a title with either word is never an organelle genome. The organelle-encoded
# genes themselves (cox1, cytb and so on) are deposited as genomic DNA, with
# "gene" and "cds" in the title, not "mRNA".
NUCLEAR_GENE_RECORD = re.compile(r"\bmRNA\b|\bprotein\b", re.IGNORECASE)


def is_direct_mito_title(title):
    """Say whether a BLAST reference title describes a mitochondrial genome.

    Takes in:  one title (the stitle column of the BLAST table).
    Does:      True when the title has one of the MITO_GENOME_WORDING wordings and
               neither word of NUCLEAR_GENE_RECORD.
    Used by:   read_blast(), which counts only the hits whose title passes this test.
    """
    return bool(MITO_GENOME_WORDING.search(title)) and not NUCLEAR_GENE_RECORD.search(title)


# The three values allowed for --mode (set by parameters.mitochondria.mode in
# the workflow configuration):
#   off                 do not screen: every contig is labelled "not_screened".
#   report              screen and write the report, but select nothing, so the
#                       assembly goes on exactly as it was.
#   separate_confirmed  screen, and pass the confirmed mitochondrial contigs on
#                       to be split off (the default).
VALID_MODES = {"off", "report", "separate_confirmed"}

# CHANGED in FunFlux v2: the 70 percent "very strong direct BLAST" cut-off was a
# number typed directly into classify(). It now lives here, is the default of
# the new --very-strong-qcov option, and is passed on to classify(), so it can
# be changed from the command line (the workflow passes
# parameters.mitochondria.very_strong_qcov). With the default, behaviour is
# unchanged.
DEFAULT_VERY_STRONG_QCOV = 70.0

# CHANGED in FunFlux v2: a contig that is more than this many times longer than
# the reference mitochondrial genome it matches is flagged as a possible
# concatemer (see concatemer_columns). One copy of the genome gives a ratio near
# 1.0, two copies end to end give about 2.0; 1.5 sits safely between the two.
CONCATEMER_RATIO = 1.5

# CHANGED in FunFlux 2.1.0: the three kingdoms of the species check, by their NCBI
# taxonomy IDs. A contig whose best direct organelle hit comes from an animal (Metazoa)
# or a plant (Viridiplantae) is taken to be the mitochondrial DNA of another organism,
# such as human mitochondrial DNA that came into the library with the sample, and
# screen() reports it as foreign_mitochondrial instead of confirming it (see
# kingdom_of()). Only these two kingdoms count as foreign: a best hit placed nowhere
# more precise than "Eukaryota" (an environmental or unclassified record) means
# "cannot tell", and so does a fungus-like organism outside Fungi, such as an oomycete,
# whose own mitochondrion must not be called foreign.
KINGDOMS = {
    4751: "Fungi",
    33208: "Metazoa",
    33090: "Viridiplantae",
}
FOREIGN_KINGDOMS = {"Metazoa", "Viridiplantae"}


# ── Small helpers: FASTA and ID files, output folders, table cells ───────────

def fasta_id(header):
    """Return the contig name: the first word of a FASTA header line.

    A header such as "NODE_1_length_123456_cov_10.512345 some comment" becomes
    "NODE_1_length_123456_cov_10.512345". BLAST (its qseqid column) and the
    assemblers use this first word as the contig name, so it is what links a
    BLAST hit to a FASTA record. Called by read_fasta().
    """
    return header.split()[0]


def ensure_parent(path):
    """Create the folder that an output file will go into, if it is missing.

    Snakemake normally creates output folders itself; this makes the script
    equally safe to run by hand. Called by write_fasta(),
    write_assembly_evidence() and screen() before they open an output file.
    """
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def read_fasta(path):
    """Read a FASTA file in order, and refuse two records with the same contig name.

    Takes in a FASTA file: the draft assembly screened by "screen", or the
    finished assembly (contigs_all.fasta) cut up by "split". Gives back an
    ordered dictionary {contig name: (full header text without ">", whole
    sequence on one line)}. The order of the file is kept, so the split files
    list contigs in the same order as the assembly. Two records with the same
    contig name would make the BLAST hits ambiguous, so that stops the run; so
    does a file without any record. Lines before the first header are ignored,
    and bytes that are not valid UTF-8 text are replaced instead of stopping
    the run.
    """
    records = OrderedDict()
    header = None
    sequence = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                # A new header: the record collected so far is complete, store it.
                if header is not None:
                    identifier = fasta_id(header)
                    if identifier in records:
                        raise ValueError(f"Duplicate FASTA identifier: {identifier}")
                    records[identifier] = (header, "".join(sequence))
                header = line[1:].strip()
                sequence = []
            elif header is not None:
                # A sequence line of the current record; the pieces are joined
                # into one string when the record is stored.
                sequence.append(line.strip())
    # The last record has no header after it, so it is stored here.
    if header is not None:
        identifier = fasta_id(header)
        if identifier in records:
            raise ValueError(f"Duplicate FASTA identifier: {identifier}")
        records[identifier] = (header, "".join(sequence))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def write_fasta(path, records):
    """Write FASTA records, 80 bases per line.

    Takes in a list of (header text, sequence) pairs: the nuclear or the
    mitochondrial contigs picked by "split". 80 bases per line is the usual
    FASTA layout, and QUAST, BUSCO, the ITS search and funannotate all read it.
    The sequence itself is never changed; only the line width may differ from
    the input file. An empty list gives an empty file (the mitochondrial file
    when no contig was confirmed).
    """
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def read_id_file(path):
    """Read one contig name per line, skipping blank lines and comments.

    Takes in confirmed_mitochondrial.ids, the list "screen" wrote. Lines
    starting with "#" and empty lines are ignored. Gives back a set of names
    that split() uses to decide which contigs are mitochondrial. An empty file
    (the case in "report" and "off" mode) gives an empty set, which is valid.
    """
    identifiers = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                identifiers.add(value)
    return identifiers


def is_yes(value):
    """Read a Y/N table cell: Y, YES, TRUE or 1 (any capitalisation) mean yes.

    Anything else, including an empty cell, means no. Used by read_flye_info()
    and read_assembly_evidence() for the circular and repeat columns of the
    assembler tables.
    """
    return str(value or "").strip().upper() in {"Y", "YES", "TRUE", "1"}


def as_float(value, label):
    """Turn table text into a number, or stop with a message naming the column.

    "label" says which column the text came from, so a damaged input file gives
    "Invalid BLAST qcovs value: 'abc'" instead of a bare Python traceback. Used
    by every reader of a table in this file.
    """
    try:
        return float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {label} value: {value!r}") from error


# ── BLAST evidence: the best direct organelle hit of each contig ─────────────

def read_blast(path):
    """Find, for each contig, its best hit to a mitochondrial genome.

    Takes in the workflow's BLAST table ("blastout", from rule blast_contigs),
    the same file that BlobTools uses for decontamination. It is tab-separated
    with 15 columns:
      1 qseqid    2 staxids    3 bitscore   4 pident     5 evalue
      6 length    7 qlen       8 slen       9 qcovs     10 qcovhsp
     11 sskingdoms  12 scomnames  13 sscinames  14 sblastnames  15 stitle
    (qlen is the length of the contig, slen the length of the reference
    sequence; the columns are counted from one.) Only columns 1, 2, 7, 8, 9 and
    15 are read. Column 2 (staxids, the taxonomy ID of the reference's organism)
    is kept for the best hit only, for the species check in screen() (CHANGED in
    FunFlux 2.1.0).

    Only direct organelle hits count (see is_direct_mito_title()). For each contig
    the best one is kept: the hit with the highest qcovs; on a tie the earlier
    line wins. This is the cautious choice: taking the single best reference,
    rather than adding up several, can only underestimate how much of the
    contig is organelle DNA. qcovs is only as complete as the alignments BLAST
    saved for each reference, which is why the workflow keeps up to 5 per
    reference by default (-max_hsps 5, parameters.blast_max_hsps); with 1, a
    mitochondrial contig that aligns to its reference in several blocks looks
    only partly covered (docs/methods/blast-hit-depth.md).

    Gives back {contig name: {"best_qcov", "best_title", "best_qlen",
    "best_slen", "best_taxids", "direct_hits"}}. direct_hits counts every direct organelle line
    of the contig, so a reference saved in five pieces counts five times.
    Contigs without a direct organelle hit are simply absent. screen() looks
    each contig up here.

    Only direct organelle lines are checked: there, a qcovs, qlen or slen that
    is not a number (text such as "abc", or an empty cell) stops the run with
    a message naming the column. The numbers on all other lines are never read.
    """
    evidence = {}
    # Bytes that are not valid UTF-8 (possible in the free-text titles) are
    # replaced instead of stopping the run.
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            # A line with fewer than 15 columns (a blank line, for example) is
            # not a hit line; skip it.
            if len(row) < 15:
                continue
            contig = row[0].strip()
            title = row[14].strip()
            # Only hits whose reference IS an organelle genome count here.
            if not contig or not is_direct_mito_title(title):
                continue
            qcov = as_float(row[8], "BLAST qcovs")
            # CHANGED in FunFlux v2: also read the contig length (qlen, column 7)
            # and the length of the reference organelle genome (slen, column 8).
            # concatemer_columns() divides the contig length by slen to notice a
            # contig that is two genome copies end to end. It is given the contig
            # length counted from the FASTA, so best_qlen is stored but not used
            # by this script. A reference length of zero or less cannot be real
            # and would break that division, so it stops the run with a clear
            # message, just as text that is not a number does.
            query_length = as_float(row[6], "BLAST qlen")
            subject_length = as_float(row[7], "BLAST slen")
            if subject_length <= 0:
                raise ValueError(f"Invalid BLAST slen value: {row[7]!r}")

            # Count every direct organelle line of this contig, the best one or not.
            previous = evidence.get(contig)
            if previous is None:
                hit_count = 1
            else:
                hit_count = previous["direct_hits"] + 1

            # The contig's first line, or a line that covers more of the contig
            # than the best so far, becomes its best hit.
            if previous is None or qcov > previous["best_qcov"]:
                evidence[contig] = {
                    "best_qcov": qcov,
                    "best_title": title,
                    # CHANGED in FunFlux v2: lengths of the contig and of the
                    # reference for the best hit, kept as whole numbers of bases.
                    "best_qlen": int(query_length),
                    "best_slen": int(subject_length),
                    # CHANGED in FunFlux 2.1.0: the NCBI taxonomy ID(s) of the
                    # organism of the best hit, as BLAST wrote them: one ID ("5180",
                    # Sclerotinia sclerotiorum), several joined by ";", or "0" when the
                    # BLAST database carries no taxonomy ID for that sequence.
                    "best_taxids": row[1].strip(),
                    "direct_hits": hit_count,
                }
            else:
                # Not better than the hit we already hold (a tie included); only count it.
                previous["direct_hits"] = hit_count
    return evidence


# ── The organism of the best hit: fungus, animal or plant? (CHANGED in FunFlux 2.1.0) ──
# The title check above says whether a reference IS a mitochondrial genome, not whose.
# A contig of human mitochondrial DNA (a few human reads in an Illumina library can
# assemble into one) matches human mitochondrial genomes over its whole length and
# would pass every check of classify(). The two functions below look the organism of
# the best hit up in the NCBI taxonomy, which the BLAST database folder carries as
# nodes.dmp (the same file BlobTools reads), and screen() reports such contigs as
# foreign_mitochondrial instead of confirming them: they stay in the assembly.

def read_taxonomy_parents(path):
    """Read the parent of every taxon from the NCBI taxonomy file nodes.dmp.

    Takes in:  nodes.dmp of the NCBI taxonomy dump (--taxonomy-nodes): one taxon
               per line, the fields separated by a tab, a bar and a tab; field 1 is
               the taxonomy ID and field 2 the ID of its parent (the root, 1, is its
               own parent).
    Does:      keeps those two numbers of every line. A line that does not start
               with two numbers is skipped.
    Gives back {taxonomy ID: parent ID}, used by kingdom_of(). The full NCBI file
    has about 2.9 million taxa (July 2026) and is read in a few seconds.
    """
    parents = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            # Only the first two fields are needed, so the line is cut at its first two
            # separators only (the 2); the rest of the line stays in one piece, unread.
            fields = line.split("\t|\t", 2)
            if len(fields) < 2:
                continue
            try:
                parents[int(fields[0])] = int(fields[1])
            except ValueError:
                continue
    return parents


def kingdom_of(taxids_text, parents):
    """Say whether the organism of a BLAST hit is a fungus, an animal or a plant.

    Takes in:  the staxids cell of the hit (one taxonomy ID, several joined by ";",
               or "0" when the BLAST database carries no taxonomy ID for that
               sequence) and the parent table of read_taxonomy_parents().
    Does:      walks from each ID up through its parents towards the root, and stops
               at the first of the KINGDOMS it meets (Fungi, Metazoa, Viridiplantae).
    Gives back "Fungi" when at least one ID lies inside Fungi; otherwise "Metazoa" or
    "Viridiplantae" when an ID lies in one of those; otherwise None, which means
    "cannot tell": an empty cell, "0", text, an ID missing from the taxonomy file, an
    organism placed in none of the three kingdoms (an environmental record filed
    under Eukaryota, an oomycete), or a damaged file whose walk never ends. screen()
    then leaves the call as classify() made it.
    """
    found = None
    for text in taxids_text.split(";"):
        text = text.strip()
        if not text.isdigit():
            continue
        taxid = int(text)
        # Walk up to the root. The step limit only guards against a damaged file with
        # a loop in it: the deepest taxon of the NCBI tree is about 40 levels down.
        for _ in range(200):
            if taxid in KINGDOMS:
                if KINGDOMS[taxid] == "Fungi":
                    return "Fungi"
                found = KINGDOMS[taxid]
                break
            parent = parents.get(taxid)
            if parent is None or parent == taxid:
                break
            taxid = parent
    return found


# ── Assembler evidence: read depth and circularity from Flye or hifiasm ──────

# The columns of the common assembler-evidence table (assembly_evidence.tsv,
# ASSEMBLY_EVIDENCE in 00_common.smk). Both evidence-flye and evidence-hifiasm
# write exactly these columns, so "screen" does not need to know which assembler
# made the contigs.
#   contig    contig name (matches the FASTA header and the BLAST qseqid)
#   source    "flye" or "hifiasm"
#   circular  Y if the assembler reports the contig as a closed circle
#   coverage  assembler-reported read depth (mean for Flye, rd:i tag for hifiasm);
#             empty when the assembler gave none
#   repeat    Y if Flye marks the contig as a repeat in its assembly graph
#   alt_group label of an alternative-haplotype group, "*" when there is none
EVIDENCE_FIELDS = [
    "contig", "source", "circular", "coverage", "repeat", "alt_group",
]


def write_assembly_evidence(path, records):
    """Write the common assembler-evidence table that "screen" reads.

    Takes in the per-contig records made by read_flye_info() or
    read_hifiasm_gfa(), and writes them as a tab-separated file with the columns
    of EVIDENCE_FIELDS: yes/no values become Y/N, read depth gets six decimals,
    and an unknown depth gives an empty cell. Called by evidence_flye() and
    evidence_hifiasm(). "screen" reads the file back with
    read_assembly_evidence() when it is given --assembly-evidence.
    """
    ensure_parent(path)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EVIDENCE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for contig, record in records.items():
            writer.writerow({
                "contig": contig,
                "source": record["source"],
                "circular": "Y" if record["circular"] else "N",
                "coverage": "" if record["coverage"] is None else f"{record['coverage']:.6f}",
                "repeat": "Y" if record["repeat"] else "N",
                "alt_group": record["alt_group"],
            })


def median_coverage(records):
    """Return the median read depth over the contigs that have one, or None.

    Most contigs are nuclear, so this median (one value per contig, whatever
    its length) is the "normal" nuclear depth. A mitochondrial genome is present
    in many copies per cell, so its depth is several times higher. classify()
    compares each contig with this median. Takes in {contig name: record} from
    read_assembly_evidence() or, on the --flye-info route of screen(), from
    read_flye_info(). Returns None when no contig has a depth (for example a
    hifiasm graph whose contigs carry no rd:i tag).
    """
    depths = []
    for record in records.values():
        if record["coverage"] is not None:
            depths.append(record["coverage"])
    if not depths:
        return None
    return statistics.median(depths)


def read_flye_info(path):
    """Read Flye circularity, graph flags, and long-read coverage by contig.

    Takes in Flye's assembly_info.txt, written next to assembly.fasta. It is
    tab-separated with one line per contig and the columns:
      #seq_name  length  cov.  circ.  repeat  mult.  alt_group  graph_path
    cov. is the mean read depth, circ. is Y when Flye closed the contig into a
    circle, repeat is Y when the contig is a repeat in the assembly graph, and
    alt_group names a group of alternative copies ("*" when there is none).
    Gives back {contig name: record} with the keys source, coverage, circular,
    repeat and alt_group. Called by evidence_flye(), which saves the records
    with write_assembly_evidence(), and by screen() when it is given
    --flye-info. Stops with a message if a depth is not a number or if no
    contig is found.
    """
    records = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            # Flye calls the name column "#seq_name"; "seq_name" is accepted as
            # well. A line without a name is skipped.
            contig = (row.get("#seq_name") or row.get("seq_name") or "").strip()
            if not contig:
                continue
            # A table without a "cov." column is read from a "coverage" column
            # instead (and one without "circ." from "circular", below). A line
            # with no depth at all stops the run.
            coverage = as_float(row.get("cov.", row.get("coverage")), "Flye coverage")
            records[contig] = {
                "source": "flye",
                "coverage": coverage,
                "circular": is_yes(row.get("circ.", row.get("circular"))),
                "repeat": is_yes(row.get("repeat")),
                "alt_group": (row.get("alt_group") or "*").strip(),
            }
    if not records:
        raise ValueError(f"No Flye records found in {path}")
    return records


def read_hifiasm_gfa(path):
    """Read Hifiasm primary GFA segments and their optional rd:i read depth.

    Takes in hifiasm's primary-contig GFA (hifiasm.p_ctg.gfa). A GFA
    (Graphical Fragment Assembly) file has one line per item, and the first
    letter says what it is. Only the "S" lines are used: they are the contigs
    ("segments"), laid out as
      S  <contig name>  <sequence>  <tags such as LN:i:30000 and rd:i:150>
    where the rd:i tag is the read depth (some hifiasm runs omit it). The other
    line types (links, read placements) are ignored.

    hifiasm gives no repeat or alternative-group flags, so those are filled
    with "no" and "*". Gives back {contig name: record}, in the same layout as
    read_flye_info(). Called by evidence_hifiasm(). Stops with a message if no
    S line is found.
    """
    records = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            # Keep only S lines that have at least a name and a sequence field.
            if len(fields) < 3 or fields[0] != "S":
                continue
            contig = fields[1].strip()
            if not contig:
                continue
            # The tags follow the sequence. The depth is the number after the
            # last colon of the rd:i tag (rd:i:150 gives 150); no tag, no depth.
            coverage = None
            for tag in fields[3:]:
                if tag.startswith("rd:i:"):
                    coverage = as_float(tag.rsplit(":", 1)[-1], "Hifiasm rd:i coverage")
                    break
            # CHANGED in FunFlux v2: hifiasm puts its own verdict on circularity
            # in the last letter of every contig name, c for circular and l for
            # linear (for example ptg000002c). v0.2 ignored it and recorded False
            # for every contig, so a hifiasm-closed mitochondrial circle was
            # reported as linear. Now a name ending in c is recorded as circular.
            is_circular = contig.endswith("c")
            records[contig] = {
                "source": "hifiasm",
                "coverage": coverage,
                "circular": is_circular,
                "repeat": False,
                "alt_group": "*",
            }
    if not records:
        raise ValueError(f"No Hifiasm segment records found in {path}")
    return records


def read_assembly_evidence(path):
    """Read the common assembler-evidence table and work out its median depth.

    Takes in the assembly_evidence.tsv written by evidence-flye or
    evidence-hifiasm (columns: see EVIDENCE_FIELDS). Stops with a message if a
    required column is missing, a depth is not a number, or no contig is found.
    Gives back two things: {contig name: record} and the median read depth over
    all contigs (None if none has a depth). Called by screen() when it is given
    --assembly-evidence; screen() hands both to classify().
    """
    records = {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        # Check the header first: all six columns must be there.
        required_columns = {"contig", "source", "circular", "coverage", "repeat", "alt_group"}
        found_columns = set(reader.fieldnames or [])
        missing = required_columns - found_columns
        if missing:
            raise ValueError(
                f"Assembly evidence is missing required column(s): {', '.join(sorted(missing))}"
            )
        for row in reader:
            contig = (row.get("contig") or "").strip()
            if not contig:
                continue
            # An empty cell means "not recorded": the source becomes "unknown",
            # the depth None and the alternative group "*" (none).
            source = (row.get("source") or "").strip()
            if not source:
                source = "unknown"
            coverage_text = (row.get("coverage") or "").strip()
            if coverage_text:
                coverage = as_float(coverage_text, "assembler coverage")
            else:
                coverage = None
            alt_group = (row.get("alt_group") or "").strip()
            if not alt_group:
                alt_group = "*"
            records[contig] = {
                "source": source,
                "coverage": coverage,
                "circular": is_yes(row.get("circular")),
                "repeat": is_yes(row.get("repeat")),
                "alt_group": alt_group,
            }
    if not records:
        raise ValueError(f"No assembly evidence records found in {path}")
    return records, median_coverage(records)


# ── The evidence-flye and evidence-hifiasm subcommands ───────────────────────

def evidence_flye(args):
    """Turn Flye's assembly_info.txt into the common assembler-evidence table.

    The "evidence-flye" subcommand: assembly_info.txt in (--flye-info),
    assembly_evidence.tsv out (--output), then one summary line. Rule
    assembly_evidence runs it when the run uses Flye (see the note at the top
    of the file); "screen" reads its table through --assembly-evidence.
    """
    records = read_flye_info(args.flye_info)
    write_assembly_evidence(args.output, records)
    print(f"Assembly evidence: wrote {len(records)} Flye contig record(s).")


def evidence_hifiasm(args):
    """Turn hifiasm's primary GFA into the common assembler-evidence table.

    The "evidence-hifiasm" subcommand: the primary GFA in (--gfa),
    assembly_evidence.tsv out (--output). Rule assembly_evidence runs it when
    the run uses hifiasm (see the note at the top of the file); "screen" reads
    its table through --assembly-evidence. The line it prints says how many
    contigs carried a read-depth tag.
    """
    records = read_hifiasm_gfa(args.gfa)
    write_assembly_evidence(args.output, records)
    # Count the contigs that had an rd:i tag, for the summary line.
    coverage_count = 0
    for record in records.values():
        if record["coverage"] is not None:
            coverage_count += 1
    print(
        f"Assembly evidence: wrote {len(records)} Hifiasm contig record(s); "
        f"{coverage_count} include rd:i read-depth tags."
    )


# ── Deciding one contig: classify() and the concatemer columns ───────────────

def classify(evidence, assembler_record, median_depth, minimum_qcov, minimum_ratio, mode,
             very_strong_qcov=DEFAULT_VERY_STRONG_QCOV):
    """Decide, for ONE contig, whether it is mitochondrial, and say why.

    Called by screen() once per contig. Takes in:
      evidence          the best direct organelle BLAST hit for this contig as
                        made by read_blast(), or None if it had none
      assembler_record  what the assembler recorded for this contig (depth,
                        circular, repeat, alt_group), or None when there is
                        none (no assembler evidence in contigs and illumina
                        modes, or a contig the assembler table does not list)
      median_depth      the median read depth over all contigs, the "normal"
                        nuclear depth (None when unknown)
      minimum_qcov      percent of the contig a direct organelle hit must cover
                        to count as confirming (--min-direct-blast-qcov)
      minimum_ratio     how many times the median depth a contig must reach to
                        count as "high coverage" (--min-assembler-coverage-ratio)
      mode              off, report or separate_confirmed
      very_strong_qcov  CHANGED in FunFlux v2 (was the literal 70.0): percent
                        coverage from which a direct hit alone confirms a contig
                        even without supporting assembler evidence
                        (--very-strong-qcov)
    Gives back three values (class, reason, coverage_ratio):
      class           confirmed_mitochondrial, mitochondrial_candidate,
                      circular_unclassified, high_coverage_unclassified,
                      not_mitochondrial or not_screened
      reason          a short code saying which rule decided; it is written to the
                      "reason" column of the report, so every call can be checked
      coverage_ratio  contig depth divided by the median depth, or None
    Only "confirmed_mitochondrial" leads to a contig being split off. Circularity
    or high depth alone only ever give a "...unclassified" or "candidate" label.
    CHANGED in FunFlux 2.1.0: screen() turns a confirmed or candidate call into
    foreign_mitochondrial when the best direct hit is an animal or a plant (see
    kingdom_of()).
    The checks below run in this order, and the first one that applies decides.
    """
    # Screening switched off: nothing is examined, nothing is confirmed.
    if mode == "off":
        return "not_screened", "screening_disabled", None

    # How much of the contig its best direct organelle hit covers; 0 without one.
    best_qcov = evidence["best_qcov"] if evidence else 0.0

    # No assembler record for this contig: BLAST is the only information. This
    # is always so in contigs and illumina modes, and for a contig the
    # assembler table does not list.
    if assembler_record is None:
        # The hit covers enough of the contig: confirmed.
        if best_qcov >= minimum_qcov:
            return "confirmed_mitochondrial", "strong_direct_blast", None
        # There is a direct hit, but it covers too little of the contig: only a candidate.
        if evidence:
            return "mitochondrial_candidate", "weak_direct_blast", None
        # No organelle hit at all: leave the contig in the nuclear assembly.
        return "not_mitochondrial", "no_mitochondrial_evidence", None

    # From here on the assembler's evidence is available as well.
    # Depth relative to the typical nuclear contig; stays None if either value is missing.
    coverage_ratio = None
    if assembler_record["coverage"] is not None and median_depth:
        coverage_ratio = assembler_record["coverage"] / median_depth
    # A "clean" circle: closed by the assembler, not a graph repeat, and not an
    # alternative copy of another contig.
    clean_circular = (
        assembler_record["circular"]
        and not assembler_record["repeat"]
        and assembler_record["alt_group"] == "*"
    )
    high_coverage = coverage_ratio is not None and coverage_ratio >= minimum_ratio

    # BLAST says organelle AND the assembler agrees (high depth or a clean
    # circle): confirmed. The reason names the support; a contig that is both
    # deep and a clean circle gets direct_blast_plus_high_coverage.
    if best_qcov >= minimum_qcov and (high_coverage or clean_circular):
        support = "high_coverage" if high_coverage else "clean_circularity"
        return "confirmed_mitochondrial", f"direct_blast_plus_{support}", coverage_ratio
    # A near-complete direct organelle match is decisive even when the assembler
    # supplies no supporting evidence (for example, a hifiasm GFA without rd:i).
    # CHANGED in FunFlux v2: the cut-off is the very_strong_qcov argument
    # (default 70.0, set with --very-strong-qcov) instead of a fixed 70.0.
    if best_qcov >= very_strong_qcov:
        return "confirmed_mitochondrial", "very_strong_direct_blast", coverage_ratio
    # A direct hit exists but is not enough: it covers less than the minimum, or
    # the assembler does not back it up (no high depth, no clean circle), and it
    # is not very strong either. Worth a look, not enough to split the contig off.
    if evidence:
        return "mitochondrial_candidate", "direct_blast_without_support", coverage_ratio
    # No organelle hit from here on. A contig that is circular AND deep looks like
    # a mitochondrion, but BLAST did not confirm it, so it stays a candidate.
    if assembler_record["circular"] and high_coverage:
        return "mitochondrial_candidate", "circular_and_high_coverage_without_blast", coverage_ratio
    # Circular only, or deep only: a hint for the reader, never a call.
    if assembler_record["circular"]:
        return "circular_unclassified", "circularity_only", coverage_ratio
    if high_coverage:
        return "high_coverage_unclassified", "coverage_only", coverage_ratio
    # Nothing points to an organelle.
    return "not_mitochondrial", "no_mitochondrial_evidence", coverage_ratio


def concatemer_columns(contig_length, evidence, assembler_record):
    """Work out the three concatemer columns of the screen report for one contig.

    CHANGED in FunFlux v2: this whole function is new.

    Background: an assembler can occasionally join two copies of the
    mitochondrial genome end to end into one contig (a tandem duplicate, also
    called a concatemer). Such a contig still matches the organelle
    reference (so it is rightly confirmed as mitochondrial) but it is about
    twice as long as the reference, which is how it can be noticed.

    Called by screen() once per contig. Takes in:
      contig_length     length of the contig in bases, counted from the FASTA
      evidence          the dictionary read_blast() made for this contig
                        (None when there was no direct organelle hit)
      assembler_record  the assembler's record for this contig (None when there
                        is none, as for every contig in contigs and illumina
                        modes)
    Gives back three pieces of text for the report, in this order:
      best_mito_slen   length of the reference genome of the best direct hit
      length_ratio     contig length divided by that length, rounded to 2 decimals
      concatemer_flag  Y when the ratio is above CONCATEMER_RATIO and the contig
                       is not recorded as circular (a contig that the assembler
                       closed into a circle is most likely a genuine, larger
                       mitochondrial genome), otherwise N
    Without a direct organelle hit there is no reference to compare with: the
    first two are empty and the flag is N. The flag is judged on the rounded
    ratio, so it always agrees with the number written in the report.

    These columns are information for the person reading the report. screen()
    never uses them to classify, to confirm or to split.
    """
    if not evidence:
        return "", "", "N"

    reference_length = evidence["best_slen"]
    length_ratio = round(contig_length / reference_length, 2)

    # No assembler record for this contig (always so in contigs and illumina
    # modes) means circularity is unknown, which counts as "not recorded as
    # circular".
    recorded_circular = bool(assembler_record and assembler_record["circular"])
    if length_ratio > CONCATEMER_RATIO and not recorded_circular:
        concatemer_flag = "Y"
    else:
        concatemer_flag = "N"

    return str(reference_length), f"{length_ratio:.2f}", concatemer_flag


# ── The screen and split subcommands ─────────────────────────────────────────

def screen(args):
    """Classify every contig and write the report and the confirmed-ID list.

    The "screen" subcommand, run by rule screen_mitochondria. Takes in:
      --contigs            the draft assembly FASTA that BLAST was run on
                           (contigs_filt.fasta in contigs mode)
      --blast              the BLAST table of that assembly (see read_blast)
      --assembly-evidence  optional: the table from evidence-flye or
                           evidence-hifiasm (long-read modes; empty in contigs mode)
      --flye-info          optional, older way: Flye's assembly_info.txt directly
      --taxonomy-nodes     optional, CHANGED in FunFlux 2.1.0: nodes.dmp of the NCBI
                           taxonomy, for the species check (see kingdom_of); without
                           it the check is skipped
      --mode and the thresholds (see classify)
    Writes:
      --report         mitochondrial_screen.tsv, one line per contig, in the
                       order of the FASTA. Columns:
                       contig, length, classification, reason; best_mito_qcovs
                       (percent of the contig covered by the best direct organelle
                       hit), direct_mito_hits (how many direct organelle lines the
                       BLAST table has for the contig; see read_blast),
                       best_mito_title (the reference named by the best hit);
                       assembler_source, assembler_circular, assembler_coverage,
                       assembler_coverage_ratio, assembler_repeat,
                       assembler_alt_group (empty without assembler evidence);
                       and, CHANGED in FunFlux v2, three informational columns at
                       the end: best_mito_slen, length_ratio, concatemer_flag
                       (see concatemer_columns); and, CHANGED in FunFlux 2.1.0,
                       best_mito_kingdom: Fungi, Metazoa or Viridiplantae for the
                       organism of the best direct hit, empty when that cannot be
                       told or the check did not run.
      --confirmed-ids  confirmed_mitochondrial.ids, one confirmed contig per line
                       (empty unless --mode is separate_confirmed). Read by rule
                       select_contigs, which keeps these contigs whatever the
                       genus filter says, and by rule separate_mitochondria,
                       which splits them off.
    Prints one summary line (in the rule's log): how many contigs were
    confirmed, and how many were selected for separation; and, when there are
    any, a second line that names the foreign_mitochondrial contigs.
    """
    # 1. Read the inputs: the contigs, the BLAST evidence and, when it was given,
    #    the assembler evidence (with the median read depth over all contigs).
    contigs = read_fasta(args.contigs)
    evidence_by_contig = read_blast(args.blast)
    assembler_by_contig = None
    median_depth = None
    if args.assembly_evidence:
        assembler_by_contig, median_depth = read_assembly_evidence(args.assembly_evidence)
    elif args.flye_info:
        # The older way, kept from script v0.1: read Flye's assembly_info.txt
        # directly. No FunFlux v2 rule uses it; --assembly-evidence is the
        # current way.
        assembler_by_contig = read_flye_info(args.flye_info)
        median_depth = median_coverage(assembler_by_contig)
    # CHANGED in FunFlux 2.1.0: the NCBI taxonomy, for the species check below.
    taxonomy_parents = read_taxonomy_parents(args.taxonomy_nodes) if args.taxonomy_nodes else None

    # 2. Classify every contig and write the report, one line per contig.
    ensure_parent(args.report)
    ensure_parent(args.confirmed_ids)
    # The report gets every call, in every mode. Which calls are acted on is
    # decided only in step 3, so "report" mode shows the calls without changing
    # the nuclear assembly.
    confirmed = []
    foreign = []
    with open(args.report, "w", encoding="utf-8", newline="") as report:
        writer = csv.writer(report, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "contig", "length", "classification", "reason", "best_mito_qcovs",
            "direct_mito_hits", "best_mito_title", "assembler_source",
            "assembler_circular", "assembler_coverage", "assembler_coverage_ratio",
            "assembler_repeat", "assembler_alt_group",
            # CHANGED in FunFlux v2: three informational columns added at the END,
            # so the position of every older column stays the same.
            "best_mito_slen", "length_ratio", "concatemer_flag",
            # CHANGED in FunFlux 2.1.0: one more column at the end.
            "best_mito_kingdom",
        ])
        # One report line per contig, in the order of the FASTA file.
        for contig, (_, sequence) in contigs.items():
            # The contig's best direct organelle hit (None without one) and its
            # assembler record (None without assembler evidence). A contig that
            # the assembler table does not list also gets None, and is then
            # judged on BLAST alone.
            evidence = evidence_by_contig.get(contig)
            assembler_record = assembler_by_contig.get(contig) if assembler_by_contig else None
            classification, reason, coverage_ratio = classify(
                evidence,
                assembler_record,
                median_depth,
                args.min_direct_blast_qcov,
                args.min_assembler_coverage_ratio,
                args.mode,
                # CHANGED in FunFlux v2: hand the --very-strong-qcov value on
                # to classify().
                args.very_strong_qcov,
            )
            # CHANGED in FunFlux 2.1.0: the species check. A contig whose best direct
            # organelle hit comes from an animal or a plant is another organism's
            # mitochondrial DNA. A confirmed or candidate call then becomes
            # foreign_mitochondrial: the contig is not split off, and it is not protected
            # from the genus filter of the contamination screen, which treats it like any
            # other contig. kingdom_of() answers None when it cannot tell, and then the
            # call stays as it is and the column stays empty.
            best_mito_kingdom = ""
            if taxonomy_parents is not None and evidence:
                kingdom = kingdom_of(evidence["best_taxids"], taxonomy_parents)
                if kingdom is not None:
                    best_mito_kingdom = kingdom
                if kingdom in FOREIGN_KINGDOMS and classification in (
                        "confirmed_mitochondrial", "mitochondrial_candidate"):
                    foreign.append(contig)
                    classification = "foreign_mitochondrial"
                    reason = "best_direct_hit_animal_or_plant"
            if classification == "confirmed_mitochondrial":
                confirmed.append(contig)
            # CHANGED in FunFlux v2: compute the three informational concatemer
            # columns. They are written to the report only; nothing below uses them.
            best_mito_slen, length_ratio, concatemer_flag = concatemer_columns(
                len(sequence), evidence, assembler_record
            )

            # The three BLAST columns: empty (and 0 hits) without a direct organelle hit.
            if evidence:
                best_mito_qcovs = f"{evidence['best_qcov']:.2f}"
                direct_mito_hits = evidence["direct_hits"]
                best_mito_title = evidence["best_title"]
            else:
                best_mito_qcovs = ""
                direct_mito_hits = 0
                best_mito_title = ""

            # The assembler columns: Y/N and depths when the contig has an
            # assembler record, all empty when it has none.
            if assembler_record:
                assembler_source = assembler_record["source"]
                assembler_circular = "Y" if assembler_record["circular"] else "N"
                if assembler_record["coverage"] is not None:
                    assembler_coverage = f"{assembler_record['coverage']:.2f}"
                else:
                    assembler_coverage = ""
                assembler_repeat = "Y" if assembler_record["repeat"] else "N"
                assembler_alt_group = assembler_record["alt_group"]
            else:
                assembler_source = ""
                assembler_circular = ""
                assembler_coverage = ""
                assembler_repeat = ""
                assembler_alt_group = ""
            # The depth ratio worked out by classify(); empty when there is none.
            if coverage_ratio is not None:
                assembler_coverage_ratio = f"{coverage_ratio:.2f}"
            else:
                assembler_coverage_ratio = ""

            # The values in the same order as the header line above.
            writer.writerow([
                contig,
                len(sequence),
                classification,
                reason,
                best_mito_qcovs,
                direct_mito_hits,
                best_mito_title,
                assembler_source,
                assembler_circular,
                assembler_coverage,
                assembler_coverage_ratio,
                assembler_repeat,
                assembler_alt_group,
                best_mito_slen,
                length_ratio,
                concatemer_flag,
                best_mito_kingdom,
            ])

    # 3. Write the list of confirmed contig names that the later rules use.
    # Only separate_confirmed acts on the calls. In report mode (and off) the
    # list is written empty, so no contig is protected from the genus filter
    # and none is split off: the assembly goes on unchanged.
    selected = confirmed if args.mode == "separate_confirmed" else []
    with open(args.confirmed_ids, "w", encoding="utf-8") as handle:
        for contig in selected:
            handle.write(f"{contig}\n")

    print(
        f"Mitochondrial screen: {len(confirmed)} confirmed call(s); "
        f"{len(selected)} selected for separation (mode={args.mode})."
    )
    # CHANGED in FunFlux 2.1.0: name the contigs of another organism, if any.
    if foreign:
        print(
            f"Mitochondrial screen: {len(foreign)} contig(s) of another organism's "
            f"mitochondrial DNA (foreign_mitochondrial), not set apart: {', '.join(foreign)}"
        )


def split(args):
    """Write the nuclear and the mitochondrial contigs to two FASTA files.

    The "split" subcommand, run by rule separate_mitochondria. Takes in:
      --fasta          the finished assembly (contigs_all.fasta: decontaminated
                       and, in the long-read modes, polished)
      --confirmed-ids  confirmed_mitochondrial.ids from "screen"
    Writes:
      --nuclear-fasta       every contig that is not on the list: the nuclear
                            genome (contigs_final.fasta) that goes on to QUAST,
                            BUSCO, the ITS search and gene prediction
      --mitochondrial-fasta the listed contigs (contigs_mitochondrial.fasta; an
                            empty file when the list is empty)
    The contigs keep their order and their headers. Stops with a clear message,
    before writing anything, if a listed name is missing from the FASTA or if
    nothing would be left for the nuclear assembly. Prints one summary line.
    """
    contigs = read_fasta(args.fasta)
    confirmed = read_id_file(args.confirmed_ids)
    # The list was made on the draft, but the split runs on the finished
    # assembly, which the long-read modes also polish in between. A listed name
    # that is not found here means the contig was renamed or lost on the way, or
    # the list belongs to another assembly: stop, rather than risk leaving a
    # renamed mitochondrion in the nuclear genome.
    missing = sorted(confirmed - set(contigs))
    if missing:
        raise ValueError(
            "Confirmed mitochondrial contig(s) are absent from final assembly: "
            + ", ".join(missing)
        )

    # Sort every contig into the nuclear or the mitochondrial list, keeping the
    # original order.
    nuclear = []
    mitochondrial = []
    for contig, record in contigs.items():
        if contig in confirmed:
            mitochondrial.append(record)
        else:
            nuclear.append(record)
    if not nuclear:
        raise ValueError("Mitochondrial split would leave an empty nuclear assembly.")

    # Both files are written only now, after every check has passed.
    write_fasta(args.nuclear, nuclear)
    write_fasta(args.mitochondrial, mitochondrial)
    print(
        f"Mitochondrial split: {len(mitochondrial)} mitochondrial and "
        f"{len(nuclear)} nuclear contig(s)."
    )


# ── The command line, and what runs when the script starts ───────────────────

def build_parser():
    """Define the command line: one subcommand per job, each with its options.

    Gives back the parser (argparse is Python's standard command-line reader).
    The workflow rules call the script as
    "python mitochondria.py <subcommand> --option value ...". Each subcommand
    records the function that does its job (set_defaults(func=...)), and main()
    calls that function. The tests use this parser directly too.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    # evidence-flye: Flye's assembly_info.txt -> common evidence table.
    flye_evidence_parser = commands.add_parser(
        "evidence-flye", help="normalize Flye assembly_info.txt for mitochondrial screening"
    )
    flye_evidence_parser.add_argument("--flye-info", required=True)
    flye_evidence_parser.add_argument("--output", required=True)
    flye_evidence_parser.set_defaults(func=evidence_flye)

    # evidence-hifiasm: hifiasm primary GFA -> common evidence table.
    hifiasm_evidence_parser = commands.add_parser(
        "evidence-hifiasm", help="normalize Hifiasm primary GFA read depth for mitochondrial screening"
    )
    hifiasm_evidence_parser.add_argument("--gfa", required=True)
    hifiasm_evidence_parser.add_argument("--output", required=True)
    hifiasm_evidence_parser.set_defaults(func=evidence_hifiasm)

    # screen: BLAST table + contigs (+ assembler evidence) -> report + confirmed IDs.
    screen_parser = commands.add_parser("screen", help="classify mitochondrial candidates")
    screen_parser.add_argument("--blast", required=True)
    screen_parser.add_argument("--contigs", required=True)
    screen_parser.add_argument("--report", required=True)
    screen_parser.add_argument("--confirmed-ids", required=True)
    # Assembler evidence for the long-read modes: the common table made by
    # evidence-flye or evidence-hifiasm. Left empty in contigs and illumina
    # modes (BLAST only).
    screen_parser.add_argument("--assembly-evidence", default="")
    # The older way to give Flye's table, kept from script v0.1 so that a command
    # line written for it still works. No FunFlux v2 rule uses it.
    screen_parser.add_argument("--flye-info", default="")
    # CHANGED in FunFlux 2.1.0: nodes.dmp of the NCBI taxonomy (the BLAST database
    # folder holds it), for the species check; empty skips the check.
    screen_parser.add_argument("--taxonomy-nodes", default="")
    # off, report or separate_confirmed (see VALID_MODES at the top of the file).
    screen_parser.add_argument("--mode", choices=sorted(VALID_MODES), default="separate_confirmed")
    # Percent of the contig a direct organelle hit must cover to count as confirming.
    screen_parser.add_argument("--min-direct-blast-qcov", type=float, default=35.0)
    # CHANGED in FunFlux v2: the cut-off for "very strong direct BLAST" is now an
    # option instead of a fixed 70.0 inside classify(). Default behaviour is unchanged.
    screen_parser.add_argument(
        "--very-strong-qcov", type=float, default=DEFAULT_VERY_STRONG_QCOV,
        help="percent of the contig a direct organelle hit must cover to confirm "
             "it on its own, without supporting assembler evidence (default: %(default)s)",
    )
    # How many times the median depth a contig needs to count as "high coverage".
    # The older name --min-flye-coverage-ratio (the name in script v0.1) still works
    # and means the same.
    screen_parser.add_argument(
        "--min-assembler-coverage-ratio", "--min-flye-coverage-ratio",
        dest="min_assembler_coverage_ratio", type=float, default=3.0,
    )
    screen_parser.set_defaults(func=screen)

    # split: final contigs FASTA + confirmed IDs -> nuclear FASTA + mitochondrial FASTA.
    # --nuclear and --mitochondrial are shorter spellings of the two output
    # options; both spellings work.
    split_parser = commands.add_parser("split", help="separate confirmed mitochondrial contigs")
    split_parser.add_argument("--fasta", required=True)
    split_parser.add_argument("--confirmed-ids", required=True)
    split_parser.add_argument("--nuclear-fasta", "--nuclear", dest="nuclear", required=True)
    split_parser.add_argument(
        "--mitochondrial-fasta", "--mitochondrial", dest="mitochondrial", required=True
    )
    split_parser.set_defaults(func=split)
    return parser


def main():
    """Read the command line and run the chosen subcommand."""
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    # A missing file or a damaged value becomes one readable "ERROR:" line (in the
    # rule's log, or in the Snakemake output) and exit status 1, which makes
    # Snakemake mark the rule as failed, instead of a long Python traceback.
    # Any other kind of error still shows the full traceback.
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
