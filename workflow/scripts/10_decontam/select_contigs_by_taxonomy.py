#!/usr/bin/env python3
"""Keep the assembly contigs that belong to the isolate, drop the rest, and record
a reason for every contig either way.

Where it sits in the workflow
-----------------------------
Rule select_contigs (shared/10_decontam.smk) runs this script as the last step
of the contamination screen (every read-based mode; contigs mode only with
parameters.contigs.screen: true). Before it, rule blast_contigs compares every
draft contig with NCBI nt, and BlobTools (rules blob_json and blob_table) turns
the BLAST hits into one taxon per contig and rank. The genus it reports is the
"bestsum" call: the genus whose hits add up to the highest total bit score
(BLAST's measure of how good a match is). The same table holds each contig's
coverage, which this script never reads. After it, the kept contigs
(DECONTAM_CONTIGS) go to rule finalize_assembly in contigs and illumina mode.
In longread and hybrid mode the polishing steps come first: Medaka for ONT
samples (PacBio HiFi samples skip it) and, in hybrid mode, Polypolish. See the
table at the top of shared/10_decontam.smk; only contigs mode is built so far.

The screen runs before annotation because a contaminant contig left in the
assembly would get false genes from the fungal gene predictor and would distort
QUAST, BUSCO and the ITS search. Even so, FunFlux removes nothing by default:
mode "off" keeps every contig and still writes a verdict for each one.

Shared with BacFlux
-------------------
BacFlux, the sister workflow for bacterial genomes, runs the same script: "which
contigs belong to the isolate?" is the same question for both. This copy is the
BacFlux v2.0.1 script with one fungal feature added, --always-keep, which
protects the mitochondrial contig (see below), and one BacFlux mode taken out,
auto (see "Why FunFlux has no auto mode"). Every place where the code differs
from the BacFlux copy carries a comment starting with "FUNFLUX CHANGE", so a
later fix can be carried from one workflow to the other by searching for it.
The changes marked "readability only" were rewritten in plainer steps (loops,
named sort rules) in the October 2026 readability pass and compute exactly
what the BacFlux code computes. The comments were rewritten for FunFlux and
are not marked, except one note marked "comment only" on what the genus tables
do to fungal names. The bacterial cases the comments cite (Bacillus, CheckM
...) are the BacFlux evidence the code was built and tuned on. The default
mode differs too: BacFlux's is "auto", FunFlux's is "off".

What it reads
-------------
  --bestscore    BLOB_TABLE, from rule blob_table: one row per contig with its
                 best-scoring taxon at every rank. Only column 1 (the contig
                 ID) and column 22 (the genus) are read.
  --contigs      DRAFT_CONTIGS, the draft assembly BlobTools judged (FASTA,
                 wrapped or one line per sequence).
  --always-keep  confirmed_mitochondrial.ids (CONFIRMED_MITO_IDS), from rule
                 screen_mitochondria: IDs of contigs to keep whatever their
                 genus, one per line; blank lines and lines starting with # are
                 skipped. An empty file or value protects nothing. BacFlux has
                 no such flag.
  --sample       the sample name: picks this sample's rows in the two
                 per-sample tables, and appears in the messages.
  the policy     seven flags carrying parameters.decontamination, turned into
                 text once by _decontam_settings() in shared/00_common.smk so a
                 whole batch follows one policy: --mode (off | include |
                 exclude), --discard-no-hit, two inline genus lists
                 (--include-genera, --exclude-genera) and three table paths
                 (--include-genera-by-sample, --exclude-genera-file,
                 --sample-overrides).

What it writes
--------------
  --output-fasta  DECONTAM_CONTIGS: the kept contigs with their full original
                  headers, sequence wrapped at 80 characters per line.
  --output-list   CONTIG_LIST: the kept contig IDs, one per line. A record
                  only; no rule reads it.
  --composition   COMPOSITION: each genus's share of the assembly's DNA and of
                  its contig count. No FunFlux rule reads it; in BacFlux, rule
                  annotation parses it for Bakta's --genus hint, so read
                  write_composition() before changing its line format.
  --decisions     CONTIG_DECISIONS, contig_taxonomy_decisions.tsv: the audit
                  table, one row per contig, four columns:
                    contig          the contig ID, from column 1 of the
                                    BlobTools table (the first word of the
                                    FASTA header)
                    assigned_genus  the genus from column 22 of the BlobTools
                                    table, or "no-hit"; "no-taxonomy" for a
                                    protected contig the table has no row for
                    action          keep or remove (from decide())
                    reason          one word saying why; decide() lists them all
  Warnings and a two-line summary go to stdout and stderr (the normal output
  and the error stream), which rule select_contigs writes to
  logs/select_contigs_{sample}.log. The script stops with an error, so the
  rule fails, when no contig is kept, when a protected ID is missing from the
  FASTA, or when an input or a setting cannot be used (the functions below say
  when).

Protecting the fungal mitochondrion (--always-keep)
---------------------------------------------------
A fungal mitochondrial genome usually assembles as one small contig, at much
higher coverage than the nuclear contigs, and BLAST can give it a different
genus from the one most nuclear contigs get, or no hit at all. For a poorly
sequenced fungus, for example, nt may hold the mitochondrial genome of its
genus but few of its nuclear sequences: the mitochondrial contig is then
labelled with the right genus, while most of the assembly's DNA is labelled
with a better-sequenced relative.
Include and exclude could both remove such a contig, and so could
discard_no_hit, before the step that sets the mitochondrion apart has seen it.

--always-keep closes that gap. Rule screen_mitochondria writes the IDs of the
contigs it confirmed as mitochondrial to confirmed_mitochondrial.ids (a contig
is confirmed when its best BLAST hit to a mitochondrial genome covers enough
of it; shared/15_mitochondria.smk gives the exact rules). This script keeps
every contig on that list whatever the mode, the genus or discard_no_hit say,
with the reason protected_confirmed_mitochondrial. The kept FASTA is what the
later split into nuclear and mitochondrial contigs (rule
separate_mitochondria) is made from, so a contig lost here cannot come back
there. Two unusual cases are handled on purpose:
  - a protected contig with no row in the BlobTools table is still kept, and
    appears in the decisions file with the genus "no-taxonomy";
  - a protected ID that is not in the --contigs FASTA stops the script with an
    error. It means the list and the FASTA come from different assemblies (or
    the headers were renamed in between), and carrying on would quietly hand
    the next step an assembly without its mitochondrion.

Why FunFlux has no auto mode
----------------------------
FUNFLUX CHANGE. BacFlux has a fourth mode, "auto", its default: keep one
dominant genus and remove every other contig. Up to BacFlux v2.0.1, auto kept
the genus that holds the most CONTIGS; BacFlux v2.1.0 by default keeps the
genus that holds the most bases. FunFlux leaves auto out, because BLAST genus
labels are unreliable for poorly sequenced fungi. When nt holds few sequences
of the isolate's genus, two things go wrong with a vote by contigs:
  - "no-hit" takes part in the vote unless discard_no_hit is true. Many short
    contigs can get no hit at all, so no-hit can be the largest group by
    contig count while holding only a small share of the DNA. Auto then keeps
    those unplaced contigs and removes most of the genome.
  - With discard_no_hit true, no-hit leaves the vote, and the winner can be a
    mislabel: the better-sequenced relative that BLAST gives most of the
    isolate's contigs. Auto then removes the no-hit contigs and every contig
    labelled with another genus, genuine genome included. Being part of one
    genome, those contigs are as deeply covered as the contigs of the winning
    genus (compare the coverage SPAdes writes into the contig names, cov_, a
    k-mer coverage, lower than read depth).
Include and exclude remain, because there a person names the genera after
reading contig_taxonomy_decisions.tsv and the composition report. Include must
list the genera BLAST gave the isolate's contigs, not the isolate's true
genus: when BLAST gives the true genus to only a few contigs, include with
that genus alone removes most of the genome as well. The user-facing account
is in docs/analysis/decontamination.md ("Why there is no auto mode").
Asked for auto, by --mode or in the override table, the script stops with
AUTO_MODE_REFUSAL (from the table, after the sample and the file name). In a
normal run shared/00_common.smk has already refused it when it read the
configuration.

A dropped contig in hybrid mode
-------------------------------
In hybrid mode a dropped contig can take its ONT reads with it. With
parameters.hybrid.short_read_guidance on (off by default; ONT samples only), only
the Illumina reads that map to the kept contigs of the Illumina guide assembly
go on to guide Filtlong, so long reads from a dropped contig score badly and
can be filtered out. The sequence then disappears from the assembly, not merely
from the taxonomy table. BacFlux's config/config.yaml (decontamination block)
explains this chain and the fixes. Its measured case, plasmid pMYS of
K. pneumoniae ATCC BAA-2146, shows how a contig gets dropped in the first
place: in a contigs-mode run, BacFlux's auto mode removed it as Escherichia,
because E. coli fills so much of nt that its hits outscored Klebsiella in the
bestsum.

Uses only the Python standard library. Rules select_contigs and
select_guide_contigs run it in workflow/envs/python.yaml, a conda environment
that holds a pinned Python 3.12 and nothing else, like mitochondria.py.
"""

import argparse
import csv
import sys
from collections import Counter, OrderedDict


# ── Fixed settings: modes, yes/no spellings, genus aliases ───────────────────

# The three decontamination modes. --mode must be one of them; the mode column
# of the sample override table must be one of them or blank (blank keeps the
# run-wide mode):
#   off      keep every contig; the audit files are still written (FunFlux default)
#   include  keep only the genera the user listed
#   exclude  remove only the genera the user listed
# FUNFLUX CHANGE: no auto mode (BacFlux also has "auto", its default, which
# keeps the genus with the most contigs; see "Why FunFlux has no auto mode" at
# the top of this file).
VALID_MODES = {"include", "exclude", "off"}

# FUNFLUX CHANGE: no auto mode. The one message for a request for auto, written
# once here and used in both places that can receive it: the --mode value
# (refuse_auto_mode()) and the mode column of the override table
# (read_overrides(), which puts the sample and the file name in front of it).
# Without it, auto would stop the run with a bare "invalid choice" or "Invalid
# mode" that does not say why. The message gives no figure: the reasons
# behind it are described at the top of this file.
AUTO_MODE_REFUSAL = (
    "decontamination mode 'auto' is not available in FunFlux: BLAST genus "
    "labels are unreliable for poorly sequenced fungi, so auto can remove most "
    "of the genome. Use off (the default), read "
    "contig_taxonomy_decisions.tsv, then exclude the contaminant genera or "
    "include the genera that BLAST assigned to the isolate's own contigs."
)

# The accepted spellings of yes and no for discard_no_hit, compared in lower
# case. The value arrives as text: "true" or "false" from shared/00_common.smk,
# or whatever was typed into the discard_no_hit column of the override table.
# parse_bool() stops the run on any other word.
FALSE_VALUES = {"false", "0", "no", "off"}
TRUE_VALUES = {"true", "1", "yes", "on"}

# FUNFLUX CHANGE (comment only): what the two genus tables below mean for fungi.
# Both were written for bacteria, and both are harmless for fungal names:
#   - The curated aliases are all bacterial genera. An alias only matters when
#     the target genus (an include entry) is one of them, so with a fungal
#     target it does nothing.
#   - The prefix list looks only at the first letters of a name, so it can
#     shorten a fungal genus that starts with them: Neocosmospora gives
#     cosmospora, Geosmithia gives smithia. The effect is small and only ever
#     keeps more: include mode then also keeps a contig labelled Cosmospora
#     when the target is Neocosmospora (or the other way round). No
#     contig is removed because of an alias: exclude mode ignores aliases, and
#     mode off, the FunFlux default, never looks at them.
#
# Genus names that BLAST and BlobTools often split ONE organism across. Each
# genus on the left was split off from the genus on the right (its species
# were once placed in that genus). Two genera match when they share a name in
# this table, so Bacillus, Paenibacillus, Peribacillus and Priestia all match
# one another. This guards against throwing out a genuine contig of the
# isolate as a contaminant; it is not an attempt to reconcile taxonomies.
#
# Peribacillus and Priestia were split off Bacillus sensu lato in 2020, so hits
# for one genome come back under either name. Without the table, a genome
# whose contigs are split between Bacillus and Peribacillus can lose its
# chromosome: a tie on contig count is broken alphabetically, Bacillus wins, and
# the chromosome labelled Peribacillus leaves as contamination. Both names are
# listed: fixing only one would hit the same bug from the other side.
#
# Pseudoarthrobacter (with an 'o') is left out on purpose: it is not a validly
# published name. LPSN (the List of Prokaryotic names with Standing in
# Nomenclature) lists only Pseudarthrobacter (Busse 2016): when Latin joins
# two word parts, the final 'o' of pseudo- is dropped before a vowel. Do not
# add the spelling with 'o' back.
GENUS_EQUIVALENCE_ALIASES = {
    "paenibacillus": ("bacillus",),
    "peribacillus": ("bacillus",),
    "priestia": ("bacillus",),
    "pseudarthrobacter": ("arthrobacter",),
    "paenarthrobacter": ("arthrobacter",),
    "paraburkholderia": ("burkholderia",),
}

# The same idea, blunter: a known prefix is cut off to reach the parent genus.
# Bradyrhizobium, Mesorhizobium, Neorhizobium and Sinorhizobium all lead to
# rhizobium; Aeribacillus, Caldibacillus and Geobacillus to bacillus. Carried
# over from the rule of thumb the original BacFlux selector used (a regular
# expression, that is, a text pattern), and kept because it covers splits the
# curated table has no entry for.
GENUS_EQUIVALENCE_PREFIXES = (
    "brady",
    "meso",
    "neo",
    "sino",
    "aeri",
    "caldi",
    "geo",
)


# ── Normalise genus names and fold known genus splits into one ───────────────

def normalize(value):
    """Return value as text with spaces trimmed at both ends.

    None (a TSV cell that is missing altogether) and an empty cell both become
    "". Used on the settings read from the command line and the TSV tables.
    """
    return str(value or "").strip()


def normalize_genus(value):
    """Lower-case a genus name, so that "Fusarium" as BlobTools writes it and
    "fusarium" as someone typed it into the config compare equal.

    Returns the name with spaces trimmed, in lower case ("" for None or an
    empty cell). Called by genus_aliases(), is_no_hit() and decide().
    """
    return normalize(value).lower()


def genus_aliases(value):
    """Return the set of names a genus may match: the genus itself and its aliases.

    Always the genus itself, in lower case, plus the name the curated table
    (GENUS_EQUIVALENCE_ALIASES) gives it, plus what is left after cutting off a
    known prefix (GENUS_EQUIVALENCE_PREFIXES):

    - Peribacillus -> peribacillus, bacillus
    - Paraburkholderia -> paraburkholderia, burkholderia
    - Bradyrhizobium -> bradyrhizobium, rhizobium
    - Fusarium -> fusarium (no alias, the usual result for a fungal genus)

    A prefix is cut off only when at least four letters remain, so a short name
    is never cut down to a stub that would match many unrelated names.

    Called by genus_matches(), which include mode uses. Exclude mode matches
    exact names only, because a too-broad alias there would remove contigs
    instead of rescuing them.
    """
    genus = normalize_genus(value)

    # FUNFLUX CHANGE (readability only): the steps below are spelled out; BacFlux
    # has the same logic in shorter expressions, with the same result.
    # The genus itself (an empty genus gets no names, so it matches nothing).
    aliases = set()
    if genus:
        aliases.add(genus)

    # The names the curated table gives it; none for a fungal genus.
    curated_names = GENUS_EQUIVALENCE_ALIASES.get(genus, ())
    aliases.update(curated_names)

    # The name left after a known prefix, when at least four letters remain.
    for prefix in GENUS_EQUIVALENCE_PREFIXES:
        if not genus.startswith(prefix):
            continue
        rest_of_name = genus[len(prefix):]
        if len(rest_of_name) >= 4:
            aliases.add(rest_of_name)
    return aliases


def genus_matches(genus, targets):
    """Return True when a genus shares a name with any of the target genera,
    directly or through an alias (see genus_aliases()).

    Takes in a contig's genus and a list of targets, the include list. Called
    by decide() in include mode, after the exact comparison has failed, to tell
    an alias match from a miss.
    """
    # FUNFLUX CHANGE (readability only): a plain loop; BacFlux writes the same
    # test in one line, with the same result.
    genus_names = genus_aliases(genus)
    for target in targets:
        target_names = genus_aliases(target)
        shared_names = genus_names.intersection(target_names)
        if shared_names:
            return True
    return False


# ── Read the keep/drop policy from config, genus files and overrides ─────────

# FUNFLUX CHANGE (comment only): shared/00_common.smk checks the three policy
# files before any job starts (_check_decontam_tables and the functions next to
# it), with the rules of split_genera(), parse_bool(), read_lines_file(),
# read_overrides(), read_include_by_sample() and the override part of main().
# Change a rule here and change it there too, or the dry run and this script
# will disagree about what is allowed.

def refuse_auto_mode(mode):
    """Stop on --mode auto with AUTO_MODE_REFUSAL; return any other mode unchanged.

    FUNFLUX CHANGE: no auto mode (BacFlux has no such check). main() hands
    this function to argparse (the standard Python module that reads the
    command line) as type= of --mode. argparse runs it on the text given
    after --mode BEFORE it compares that text with VALID_MODES, so --mode auto
    stops here with the explanation instead of argparse's bare "invalid
    choice: 'auto'". argparse prints the message after the usage line and
    ends the script with exit status 2. Capitals are ignored ("AUTO" is
    refused the same way), as in the mode column of the override table. Any
    other text passes through, and argparse checks it against VALID_MODES
    exactly as before.
    """
    if mode.lower() == "auto":
        raise argparse.ArgumentTypeError(AUTO_MODE_REFUSAL)
    return mode


def parse_bool(value, default=False):
    """Turn a yes/no setting into True or False.

    Takes in the text of --discard-no-hit, or of a discard_no_hit cell in the
    override table. Accepts true/false, yes/no, on/off and 1/0 in any mix of
    capitals (TRUE_VALUES, FALSE_VALUES). An empty value gives default. Any
    other word stops the run with an error, because reading a typo as False
    would silently change whether unplaced contigs are discarded. Called by
    main().
    """
    text = normalize(value).lower()
    if not text:
        return default
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise ValueError(f"Cannot parse boolean value '{value}'.")


def split_genera(value):
    """Split one genus list into single genus names.

    The list can arrive in any shape the config and the tables allow: one genus,
    "GenusA;GenusB" (how shared/00_common.smk passes a YAML list),
    "GenusA,GenusB", or names separated by tabs. Commas and semicolons are
    turned into tabs, the text is cut at the tabs, and empty pieces are dropped.
    Returns the names in the order given, spaces trimmed, capitals unchanged
    ([] for an empty value). Called by main() for --include-genera and
    --exclude-genera, and by the three table readers below.
    """
    text = normalize(value)
    if not text:
        return []
    for separator in [",", ";"]:
        text = text.replace(separator, "\t")
    return [item.strip() for item in text.split("\t") if item.strip()]


def read_lines_file(path):
    """Read the --exclude-genera-file: genera to remove, one per line.

    Blank lines and lines starting with # are skipped, so the file can carry
    comments; a line may also hold several genera in any shape split_genera()
    accepts. An empty path means the option was not set, and gives []. Called
    by main(), which adds these genera to --exclude-genera.
    """
    if not path:
        return []
    genera = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genera.extend(split_genera(line))
    return genera


def read_overrides(path):
    """Read the --sample-overrides table: per-sample changes to the policy.

    Optional: an empty path gives {}. A tab-separated table with a header line.
    Columns sample and mode are required (a table without them stops the run);
    include_genera, exclude_genera and discard_no_hit may follow. A sample's row
    replaces the run-wide value of each field it fills in, and leaves every
    other sample alone: the usual case is one run-wide policy plus one awkward
    isolate. Rows whose sample cell is empty or starts with # are skipped. If a
    sample has two rows, the last one counts. A mode outside VALID_MODES, in
    any row, stops the run, even for a sample that has no row of its own
    (auto with AUTO_MODE_REFUSAL, after the sample and the file name); a blank
    mode keeps the run-wide one.

    The trap: an EMPTY include or exclude cell means "keep the run-wide list",
    NOT "match nothing". A row must therefore still carry every tab where fields
    are blank, or the values slide into the wrong columns.

    Returns {sample: {"mode", "include", "exclude", "discard_no_hit"}}. Called
    by main().
    """
    if not path:
        return {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
        # Each row is read as a dictionary: column name (from the header line)
        # to cell text.
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "mode"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Override file '{path}' must contain at least columns: sample, mode."
            )
        overrides = {}
        for row in reader:
            sample = normalize(row.get("sample"))
            if not sample or sample.startswith("#"):
                continue
            mode = normalize(row.get("mode")).lower()
            # FUNFLUX CHANGE: no auto mode (BacFlux accepts it here). An auto
            # cell, in any capitals, gets the same explanation as --mode auto
            # rather than the bare "Invalid mode" below. Like that error, it
            # names the sample and the file: any row stops every sample's job,
            # so the log of one sample may be reporting another sample's row.
            if mode == "auto":
                raise ValueError(
                    f"Sample '{sample}' in override file '{path}': {AUTO_MODE_REFUSAL}"
                )
            if mode and mode not in VALID_MODES:
                raise ValueError(
                    f"Invalid mode '{mode}' for sample '{sample}' in '{path}'."
                )
            # The two genus columns become lists here, but discard_no_hit stays
            # raw TEXT. Turning it into True/False now would need a default, and
            # the obvious False would make a blank cell look like an explicit
            # "false" and quietly override the run-wide setting. Keeping the raw
            # text lets main() tell "not given" from "given as false".
            overrides[sample] = {
                "mode": mode,
                "include": split_genera(row.get("include_genera")),
                "exclude": split_genera(row.get("exclude_genera")),
                "discard_no_hit": row.get("discard_no_hit"),
            }
    return overrides


def read_include_by_sample(path, sample):
    """Read this sample's genera from the --include-genera-by-sample table.

    Two tab-separated columns with a header line, sample and genus (a table
    without them stops the run). A sample may have several rows, and a genus
    cell may hold several genera. The compact alternative to the override table
    when many samples each need their own include list but nothing else differs.
    Returns this sample's genera in file order ([] for an empty path, or for a
    sample with no rows). Called by main(), which adds them to --include-genera.
    """
    if not path:
        return []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "genus"}.issubset(reader.fieldnames):
            raise ValueError(
                f"Include file '{path}' must contain columns: sample, genus."
            )
        genera = []
        for row in reader:
            if normalize(row.get("sample")) == sample:
                genera.extend(split_genera(row.get("genus")))
        return genera


def read_contig_ids(path):
    """Read the --always-keep file into the set of contig IDs that must be kept.

    FUNFLUX CHANGE: added for the fungal workflow.

    Takes in: a text file with one contig ID per line; in FunFlux,
              confirmed_mitochondrial.ids from rule screen_mitochondria. The IDs
              are the first-word contig names that the FASTA and the BlobTools
              table also use. Blank lines and lines starting with # are skipped.
              Every other line is taken whole, spaces trimmed at both ends:
              "NODE_12 mitochondrion" would be read as ONE ID with a space in
              it, would match no contig, and would stop the run as an ID missing
              from the FASTA.
    Returns:  a set of IDs (a collection without repeats that answers "is this
              ID in it?" at once). An empty path means --always-keep was not
              given, and gives an empty set: nothing is protected.
    Next:     main() hands the set to decide(), which keeps every contig in it,
              and checks after writing the outputs that every ID is in the FASTA.
    """
    if not path:
        return set()
    contig_ids = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            contig_id = line.strip()
            # Skip blank lines and comment lines.
            if not contig_id or contig_id.startswith("#"):
                continue
            contig_ids.add(contig_id)
    return contig_ids


# ── Read the BlobTools genus table and the assembly FASTA ────────────────────

def read_bestscore(path):
    """Read the BlobTools table into (contig, genus) rows and contigs per genus.

    Takes in BLOB_TABLE (from rule blob_table). Header lines (#), blank lines,
    rows with fewer than 22 columns and rows without a contig ID are skipped.
    Returns:
      records  [(contig, genus), ...] in table order, so that the decisions file
               lists the contigs in the order BlobTools reported them;
      counts   {genus: number of contigs}, what the contig column of the
               composition report is built from.
    Called by main().

    If no row can be read, it stops HERE rather than later on "No contigs were
    kept": that message would send the reader to the decontamination settings,
    when the real fault is an unreadable or cut-off BlobTools table.
    """
    records = []
    counts = Counter()
    # errors="replace": a byte that is not valid UTF-8 becomes a placeholder
    # character instead of stopping the run.
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 22:
                continue
            contig = fields[0].strip()

            # The genus is column 22 of `blobtools view --rank all` (genus.t in
            # its header), fields[21] counting from zero; that is why a row with
            # fewer than 22 fields is skipped above rather than half-read. An
            # empty cell means BLAST placed nothing, and becomes "no-hit", so the
            # rest of the script has one spelling to test for.
            genus = fields[21].strip() or "no-hit"
            if not contig:
                continue

            records.append((contig, genus))
            counts[genus] += 1
    if not records:
        raise ValueError(f"No contig taxonomy records could be parsed from '{path}'.")
    return records, counts


def read_fasta(path):
    """Read the assembly FASTA into {full header: sequence}, in file order.

    Takes in DRAFT_CONTIGS. Sequence lines are collected and joined at the end
    of each record, so a wrapped FASTA (several lines per sequence) is read
    directly. v1 FunFluxL first wrote a one-line-per-sequence copy
    (assembly_linearized.fasta) for the selector, as BacFlux's v1 FastaFlux did
    (contigs_filt_lin.fasta); v2 dropped that step as unnecessary. A file with
    no FASTA record stops the run. Called by main().

    The file order kept here does not set the order of the outputs: main()
    writes the kept contigs in the order of the BlobTools table rows, then any
    protected contigs the table has no row for, sorted by ID. BlobTools
    normally lists the contigs in FASTA order, so the two orders usually agree.
    """
    records = OrderedDict()
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # A new header closes the previous record.
                if header is not None:
                    records[header] = "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
    if header is not None:
        # The last record has no following header to close it, so store it here.
        records[header] = "".join(seq_chunks)
    if not records:
        raise ValueError(f"No FASTA records could be parsed from '{path}'.")
    return records


def fasta_key(header):
    """Return the contig ID: the first word of a FASTA header.

    A header can carry text after the ID, as in "contig_1 length=1234", while the
    BlobTools table and confirmed_mitochondrial.ids name only "contig_1". That
    first word is what the FASTA and the table are joined on, here and in
    mitochondria.py. Called by main() to index the FASTA.
    """
    return header.split()[0]


# ── Decide keep or remove, one contig at a time ──────────────────────────────

def is_no_hit(genus):
    """Return True for a contig BLAST could not place.

    A test for the text "no-hit" anywhere in the label, not for equality:
    read_bestscore() turns an empty genus cell into the bare "no-hit", but
    BlobTools' own label is not guaranteed to be exactly that text, and any
    label containing it means the same here. Called by decide().
    """
    return "no-hit" in normalize_genus(genus)


# FUNFLUX CHANGE: no auto mode (BacFlux has a function choose_auto_genus() here,
# which picks the genus with the most contigs for auto mode; FunFlux leaves it
# out, see "Why FunFlux has no auto mode" at the top of this file).


def decide(contig, genus, mode, include, exclude, discard_no_hit, protected):
    """Decide whether one contig is kept or removed, and say why.

    Takes in: one row of the BlobTools table (contig ID and genus); the policy
    (mode, the include and exclude lists, discard_no_hit); and protected, the
    set of IDs read from --always-keep (empty when the flag was not given).
    Called by main() for every row of the table, and for protected contigs the
    table has no row for.

    Returns ("keep" or "remove", reason). The rules are tried in this order, and
    the first that applies decides: the protected list, mode off,
    discard_no_hit, then the genus test of the mode.

    The reason becomes the fourth column of contig_taxonomy_decisions.tsv and is
    the only record of why a contig went. All the reasons, so that the audit file
    can be read without opening this function:

      protected_confirmed_mitochondrial
                             listed in --always-keep; wins over every rule below
      mode_off               mode is off; nothing is filtered
      discarded_no_hit       not placed by BLAST, and discard_no_hit is true
      included_genus         named in the include list
      included_genus_alias   an alias of a genus in the include list
      not_included_genus     not in the include list
      excluded_genus         named in the exclude list (exact name, any capitals)
      not_excluded_genus     not in the exclude list, so kept
    """
    # FUNFLUX CHANGE (--always-keep): a contig
    # on the protected list is kept before any other rule is asked, ahead of mode
    # off, discard_no_hit and every genus test below. The protected contigs are
    # the confirmed mitochondrial ones, and BlobTools can give a mitochondrion a
    # genus other than the one the policy aims at, or no hit at all, so each of
    # those rules could otherwise say "remove". The test sits here, not in
    # main(), so that decide() stays the one place where keep or remove is
    # settled. Asking it even before mode off also makes the audit file give the
    # real reason for keeping the contig instead of "mode_off". With an empty
    # protected set this test never fires.
    if contig in protected:
        return "keep", "protected_confirmed_mitochondrial"

    # Lower-case the genus and both lists, so that Fusarium and fusarium match.
    genus_key = normalize_genus(genus)
    include_keys = {normalize_genus(item) for item in include}
    exclude_keys = {normalize_genus(item) for item in exclude}

    # Mode off keeps everything, no-hit contigs included.
    if mode == "off":
        return "keep", "mode_off"

    # discard_no_hit is applied BEFORE the genus tests, so it removes unplaced
    # contigs in include and exclude mode alike; only mode off (and the
    # protected list) escape it. An unplaced contig is often short and poorly
    # covered, but it can also be real sequence of the isolate that nt has no
    # close relative for (in BacFlux, for example a small plasmid), which is why
    # removing it is a setting and not fixed.
    if is_no_hit(genus) and discard_no_hit:
        return "remove", "discarded_no_hit"

    # FUNFLUX CHANGE: no auto mode (BacFlux tests mode auto here, between
    # discard_no_hit and include, with the reasons auto_genus,
    # auto_genus_alias, not_auto_genus and auto_no_genus).

    # Include mode is strict: only the listed genera (and their aliases) stay.
    # An empty list here means include was chosen without naming a genus, so
    # stop with an error rather than remove the whole assembly.
    if mode == "include":
        if not include_keys:
            raise ValueError("Mode 'include' requires at least one genus.")
        if genus_key in include_keys:
            return "keep", "included_genus"
        if genus_matches(genus, include):
            return "keep", "included_genus_alias"
        return "remove", "not_included_genus"

    # Exclude mode is the reverse: only the listed genera go, everything else
    # stays. Matching here is EXACT, without aliases: in include mode an alias
    # rescues a contig, but here it would remove contigs the user never named.
    if mode == "exclude":
        if genus_key in exclude_keys:
            return "remove", "excluded_genus"
        return "keep", "not_excluded_genus"

    # Only reached with a mode outside VALID_MODES, which the command line and
    # read_overrides() already refuse.
    raise ValueError(f"Unsupported decontamination mode '{mode}'.")


# ── Write the composition report and warn when the call looks wrong ──────────

def write_composition(path, counts, total, bases=None, total_bases=0):
    """Write each genus's share of the assembly to COMPOSITION, for people to read.

    Takes in counts (contigs per genus, from read_bestscore()) and total (the
    number of table rows), plus bases (bp per genus) and total_bases (bp of all
    table contigs found in the FASTA), both added up in main(). Writes one line
    per genus, over every contig in the BlobTools table, kept or removed:

        Fusarium: bases 0.92  contigs 0.61

      bases    the genus's share of the assembly's DNA
      contigs  the genus's share of the contig COUNT

    Nothing here changes which contigs are kept. The two figures can disagree,
    and that is the point: a genus can hold many short contigs while another
    holds most of the genome. BacFlux's auto mode votes on the contig count
    (FunFlux has no auto mode), and a genus can win that vote while a related
    genus holds more of the DNA; the wrong part of one genome is then kept, and
    a report of contig counts alone gives no hint. In FunFlux the report helps a person
    choose the genera for include or exclude mode.

    Lines are sorted by DNA, the fairer answer to "what is this sample mostly
    made of"; genera with equal amounts are sorted by name. When no table contig
    was found in the FASTA there are no lengths to divide by, and the v1
    one-figure form is written instead, sorted by contig count:
    "Fusarium: 0.61".

    No FunFlux rule reads this file. In BacFlux, rule annotation parses it to
    pick Bakta's --genus hint, and adding the second figure already changed
    what that parse returns: it now picks the alphabetically last genus, not
    the one with the most DNA (the comment on the genus= line of BacFlux's
    annotation rule explains why). Do not change the line format again
    without reading that comment. Called by main().
    """
    # FUNFLUX CHANGE (readability only): a named sort rule and named shares;
    # BacFlux has the same logic with lambdas as sort rules and the divisions
    # written inside the output text.
    # Sort rule for (genus, amount) pairs: the largest amount first (the minus
    # sign puts it at the top), and equal amounts by genus name.
    def largest_first(genus_and_amount):
        genus, amount = genus_and_amount
        return (-amount, genus)

    with open(path, "w", encoding="utf-8") as handle:
        if bases and total_bases:
            for genus, bp in sorted(bases.items(), key=largest_first):
                base_share = bp / total_bases
                contig_share = counts.get(genus, 0) / total
                handle.write(
                    f"{genus}: bases {base_share:.2f}  "
                    f"contigs {contig_share:.2f}\n"
                )
        else:
            # The v1 one-figure form, written when there are no lengths to
            # divide by. Same "Genus: 0.87" shape, ranked by contig count.
            for genus, count in sorted(counts.items(), key=largest_first):
                contig_share = count / total
                handle.write(f"{genus}: {contig_share:.2f}\n")


# Thresholds for the two warnings of warn_if_selection_looks_wrong(). They
# change no keep/remove decision; they only decide whether a warning is printed.
#
# MAJOR_GENUS_BASE_FRACTION: a genus holding at least this share of the DNA is a
# "major" genus. CONFUSED_GENUS_COUNT: this many major genera or more trigger
# the first warning. Two is deliberate: a clean isolate has ONE genus holding
# essentially all its DNA, so a second genus above 5% is already unusual.
#
# The thresholds were set in BacFlux on a batch of bacterial isolates. A clean
# isolate has one genus holding essentially all the DNA and loses only a few
# percent of the assembly. The samples that trip the warnings are of two
# opposite kinds, which is why the warning names both causes:
#   - a genome that BLAST scattered over several genus labels, so that auto mode
#     removed a large part of the assembly and completeness fell;
#   - a genuine two-organism culture, which CheckM reports as complete and
#     heavily contaminated: two genomes in one.
# It is one batch: treat the numbers as a well-supported starting point, not a
# universal constant. The calibration used bacteria only; a poorly
# sequenced fungus can trip the first warning even in mode off, because BLAST
# spreads its one genome over several related genera.
MAJOR_GENUS_BASE_FRACTION = 0.05
CONFUSED_GENUS_COUNT = 2

# LARGE_REMOVAL_BASE_FRACTION: removing more than this share of the assembly
# triggers the second warning. A truly contaminated culture can exceed it; so
# can a filter that has just thrown away the target genome. The warning cannot
# tell the two apart; it asks a person to look. This threshold rests on
# few samples, so treat it as a first cut, not a calibrated cutoff.
LARGE_REMOVAL_BASE_FRACTION = 0.20


def warn_if_selection_looks_wrong(bases, total_bases, removed_bases):
    """Print a warning when the genus picture, or the amount being removed,
    suggests that a person should check the keep/remove calls.

    Why. A genus filter assumes that BLAST labels the isolate's contigs with
    the genera the policy expects, and that every other genus is
    contamination. That breaks when BLAST spreads ONE genome across several
    related genera, which happens when the organism has few relatives in the
    nucleotide database: different contigs match different relatives. The
    script cannot tell that from real contamination, and without these warnings
    it would carry on silently either way. The thresholds, and how they were
    chosen, are described above MAJOR_GENUS_BASE_FRACTION.

    Takes in: bases (bp per genus), total_bases and removed_bases (bp), all
              added up in main() over the table contigs found in the FASTA.
    Prints:   up to two warnings to stdout, which rule select_contigs writes to
              logs/select_contigs_{sample}.log. Returns nothing and changes no
              decision. Called by main() in every mode, mode off included.
    """
    if not total_bases:
        return

    # FUNFLUX CHANGE (readability only): plain loops and a named sort rule;
    # BacFlux writes the same steps as one-line expressions and a lambda.
    # The "major" genera: those holding at least MAJOR_GENUS_BASE_FRACTION of
    # the DNA, in the order they were first met in the BlobTools table.
    major = []
    for genus, bp in bases.items():
        if bp / total_bases >= MAJOR_GENUS_BASE_FRACTION:
            major.append(genus)
    removed_fraction = removed_bases / total_bases

    if len(major) >= CONFUSED_GENUS_COUNT:
        # Name them largest first, each with its share of the DNA, as in
        # "Fusarium 81%, Aspergillus 9%". Genera with exactly the same
        # number of bases stay in table order.
        def most_dna_first(genus):
            return -bases[genus]

        listed_parts = []
        for genus in sorted(major, key=most_dna_first):
            share = bases[genus] / total_bases
            listed_parts.append(f"{genus} {share:.0%}")
        listed = ", ".join(listed_parts)
        print(
            f"WARNING: {len(major)} genera each hold at least "
            f"{MAJOR_GENUS_BASE_FRACTION:.0%} of this assembly ({listed}). "
            "Either the sample is a mixed culture, or one genome is being split "
            "across related genera by the BLAST assignment. Check the kept/removed "
            "split in contig_taxonomy_decisions.tsv before trusting this assembly."
        )

    if removed_fraction > LARGE_REMOVAL_BASE_FRACTION:
        print(
            f"WARNING: decontamination is discarding {removed_fraction:.0%} of the "
            f"assembly ({removed_bases:,} of {total_bases:,} bp). If the kept "
            "assembly then looks incomplete but NOT contaminated, the filter has "
            "most likely removed genome rather than contamination."
        )


# ── Apply the policy and write the four output files ─────────────────────────

def main():
    # Read the policy and the two inputs, decide every contig, write the four
    # output files, then stop with an error if the result cannot be right.
    #
    # Run by rule select_contigs (shared/10_decontam.smk), which fills in every
    # flag: the input and output paths, the sample name (the {sample} wildcard,
    # the part of each path that Snakemake fills with one sample's name) and the
    # policy from parameters.decontamination. The steps, in order:
    #   1. the run-wide policy, then this sample's override row;
    #   2. the BlobTools table, the FASTA and the protected IDs;
    #   3. one pass over the table: decide, write the audit row, add up bases;
    #   4. protected contigs that the table has no row for;
    #   5. the composition report, the warnings, the kept IDs, the kept FASTA;
    #   6. the checks that end the run with an error, then a short summary.

    # The command line. The seven optional flags (the six policy settings after
    # --mode, and --always-keep) may be given with nothing after them:
    # nargs="?" (zero or one value) allows that, and the flag then takes its
    # const value, an empty text ("false" for --discard-no-hit), meaning "not
    # set". This matters because Snakemake's {...:q} quoting writes an empty
    # value as nothing at all, so a setting left blank in the config leaves its
    # flag bare on the command line, which argparse would otherwise refuse with
    # "expected one argument". --mode falls back to off, the FunFlux default,
    # only when the flag is missing; rule select_contigs always passes it, with
    # the mode that _decontam_settings("off") in shared/00_common.smk read from
    # the config (off when the config gives none). config/config.yaml also
    # writes it out as "mode: off".
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bestscore", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--sample", required=True)
    # FUNFLUX CHANGE: no auto mode (BacFlux: default "auto" and no type=).
    # type=refuse_auto_mode makes --mode auto stop with AUTO_MODE_REFUSAL
    # instead of argparse's bare "invalid choice"; see refuse_auto_mode().
    parser.add_argument(
        "--mode", type=refuse_auto_mode, choices=sorted(VALID_MODES), default="off"
    )
    parser.add_argument("--include-genera", nargs="?", const="", default="")
    parser.add_argument("--include-genera-by-sample", nargs="?", const="", default="")
    parser.add_argument("--exclude-genera", nargs="?", const="", default="")
    parser.add_argument("--exclude-genera-file", nargs="?", const="", default="")
    parser.add_argument("--sample-overrides", nargs="?", const="", default="")
    parser.add_argument("--discard-no-hit", nargs="?", const="false", default="false")
    # FUNFLUX CHANGE: --always-keep is not in the BacFlux copy. It is FunFlux's
    # own flag, and it is the seventh optional flag (BacFlux has six). Unlike the others
    # it does not come from parameters.decontamination: rule select_contigs
    # passes CONFIRMED_MITO_IDS, the confirmed_mitochondrial.ids written by rule
    # screen_mitochondria. The contigs listed there are kept whatever BlobTools
    # called them, so that the later split into nuclear and mitochondrial
    # contigs still finds the mitochondrion. Left empty, it protects nothing.
    parser.add_argument("--always-keep", nargs="?", const="", default="")
    parser.add_argument("--output-list", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--composition", required=True)
    parser.add_argument("--decisions", required=True)
    args = parser.parse_args()

    # The run-wide policy. Each genus list is the inline config value PLUS what
    # its optional table adds (this sample's rows of --include-genera-by-sample,
    # every line of --exclude-genera-file), so a long shared exclude list can
    # live in a file while one or two genera stay visible in the config.
    mode = args.mode
    include = split_genera(args.include_genera)
    include.extend(read_include_by_sample(args.include_genera_by_sample, args.sample))
    exclude = split_genera(args.exclude_genera)
    exclude.extend(read_lines_file(args.exclude_genera_file))
    discard_no_hit = parse_bool(args.discard_no_hit, default=False)

    # This sample's row in the override table, if it has one, REPLACES each
    # value it fills in, and only those. A blank cell keeps the run-wide value,
    # which is how one awkward isolate gets its own policy without changing the
    # rest of the batch.
    overrides = read_overrides(args.sample_overrides)
    if args.sample in overrides:
        override = overrides[args.sample]
        # FUNFLUX CHANGE (readability only): an if, where BacFlux writes
        # mode = override["mode"] or mode; the result is the same.
        if override["mode"]:
            mode = override["mode"]
        if override["include"]:
            include = override["include"]
        if override["exclude"]:
            exclude = override["exclude"]
        # discard_no_hit is still raw text here (see read_overrides()): None when
        # the column or the cell is missing, blank when the cell is empty. Both
        # keep the run-wide value.
        if override["discard_no_hit"] is not None and normalize(override["discard_no_hit"]):
            discard_no_hit = parse_bool(override["discard_no_hit"], default=discard_no_hit)

    # The two real inputs. They are joined on the contig ID, so a renamed contig
    # anywhere in an earlier step shows up as the missing-contig warning at the
    # end of this function.
    records, counts = read_bestscore(args.bestscore)
    fasta_records = read_fasta(args.contigs)

    # Index the sequences by contig ID, the first word of the header, which is
    # what the BlobTools table uses. Each ID keeps its full header, so the kept
    # FASTA carries the input's headers unchanged. In FunFlux the two are
    # normally the same text: SPAdes writes length and coverage into the ID
    # itself (NODE_1_length_123456_cov_10.512345), and filter_contigs cuts any
    # other header at its first space or tab. If two headers share a first word, the
    # later record wins.
    # FUNFLUX CHANGE (readability only): a plain loop; BacFlux builds the same
    # dictionary in one line.
    fasta_by_id = {}
    for header, sequence in fasta_records.items():
        contig_id = fasta_key(header)
        fasta_by_id[contig_id] = (header, sequence)

    # FUNFLUX CHANGE: read the protected IDs once, here, and hand the set to
    # decide() for every contig. Without --always-keep the set is empty and the
    # run behaves exactly like the BacFlux script.
    protected = read_contig_ids(args.always_keep)

    # The number of table rows (what the contig shares are divided by).
    # FUNFLUX CHANGE: no auto mode (BacFlux also works out the auto winner
    # here, with choose_auto_genus(), and hands it to decide()).
    total = sum(counts.values())

    # One pass over the BlobTools rows: for each contig, decide, write its audit
    # row, add its length to the base totals, and collect it if it is kept.
    kept = []
    seen = set()
    # The base totals: how much DNA each genus holds, and how much is being
    # removed. They feed the composition report and the two warnings only; no
    # keep/remove decision uses them. The lengths come from the FASTA already in
    # memory, so nothing extra is read from disk.
    bases_by_genus = Counter()
    total_bases = 0
    removed_bases = 0
    with open(args.decisions, "w", encoding="utf-8") as decisions:
        decisions.write("contig\tassigned_genus\taction\treason\n")
        for contig, genus in records:
            # FUNFLUX CHANGE: decide() takes the protected set as its last
            # argument, so that a protected contig is kept (see decide()).
            action, reason = decide(
                contig, genus, mode, include, exclude, discard_no_hit, protected
            )
            decisions.write(f"{contig}\t{genus}\t{action}\t{reason}\n")

            # A contig with no FASTA record adds no bases; the missing-contig
            # warning further down reports it.
            if contig in fasta_by_id:
                length = len(fasta_by_id[contig][1])
                bases_by_genus[genus] += length
                total_bases += length
                if action != "keep":
                    removed_bases += length

            # Collect a kept contig once, and only if the FASTA has it. `seen`
            # guards against a BlobTools table that lists a contig on two rows:
            # writing it twice would put a duplicated sequence into the assembly.
            if action == "keep" and contig in fasta_by_id and contig not in seen:
                kept.append(contig)
                seen.add(contig)

        # FUNFLUX CHANGE: a protected contig
        # may have no row in the BlobTools table at all, so the loop above never
        # met it. It is kept anyway, so that saving the mitochondrion does not
        # depend on BlobTools having scored it. decide() still makes the call
        # (it says "keep" for any protected contig). The audit row gets the
        # placeholder genus "no-taxonomy", because there is no table row to read
        # one from. These contigs are not added to the composition or to the
        # base totals, which describe only what BlobTools scored. A protected ID
        # missing from the FASTA is skipped here and stops the run further down
        # (missing_protected). sorted() gives the same order on every run.
        for contig in sorted(protected - seen):
            if contig in fasta_by_id:
                placeholder_genus = "no-taxonomy"
                action, reason = decide(
                    contig, placeholder_genus, mode, include, exclude,
                    discard_no_hit, protected,
                )
                decisions.write(f"{contig}\t{placeholder_genus}\t{action}\t{reason}\n")
                # Always true today, because every contig in this loop is
                # protected; the test keeps main() obeying decide() should its
                # rules ever change.
                if action == "keep":
                    kept.append(contig)
                    seen.add(contig)

    # The composition report covers every contig in the BlobTools table, kept
    # and removed alike; that is what makes it useful for judging whether the
    # policy was right, even when almost nothing survived.
    write_composition(args.composition, counts, total, bases_by_genus, total_bases)

    # Advisory warnings (to stdout, so into the rule's log); no decision above
    # depends on them.
    warn_if_selection_looks_wrong(bases_by_genus, total_bases, removed_bases)

    # The kept IDs as plain text (CONTIG_LIST). A declared output of rule
    # select_contigs, but only a record: no rule reads it.
    with open(args.output_list, "w", encoding="utf-8") as handle:
        for contig in kept:
            handle.write(f"{contig}\n")

    # The kept contigs (DECONTAM_CONTIGS) with their full headers, sequence
    # wrapped at 80 characters per line. The table at the top of
    # shared/10_decontam.smk says which rule reads this file next in each mode:
    # finalize_assembly in contigs and illumina mode; in longread and hybrid
    # mode the polishing steps come first (Medaka for ONT samples, Polypolish
    # in hybrid mode).
    with open(args.output_fasta, "w", encoding="utf-8") as handle:
        for contig in kept:
            header, sequence = fasta_by_id[contig]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    # Contigs that the BlobTools table names but the FASTA does not have, each
    # counted once. The two are joined on the contig ID, so this catches a
    # mismatched pair of inputs, or an earlier step that rewrote the headers.
    # Only a warning, on stderr, which rule select_contigs also writes to
    # logs/select_contigs_{sample}.log.
    missing = sorted({contig for contig, _ in records if contig not in fasta_by_id})
    if missing:
        print(
            f"WARNING: {len(missing)} contigs from BlobTools table were not found in FASTA.",
            file=sys.stderr,
        )

    # FUNFLUX CHANGE: unlike the missing
    # BlobTools contigs above, which only earn a warning, a protected ID that is
    # not in the FASTA is an error. The list comes from the mitochondrial screen
    # of this same assembly, so a missing ID means the list and the FASTA come
    # from different files (or the headers were renamed in between), and
    # carrying on would hand the next step an assembly without its
    # mitochondrion. The check sits here, after the output files are written
    # and next to the other test for mismatched inputs, as in the fungal script.
    # The error ends the script with a non-zero exit status, so the rule fails.
    missing_protected = sorted(protected - set(fasta_by_id))
    if missing_protected:
        raise ValueError(
            "Protected mitochondrial contig(s) were not found in FASTA: "
            + ", ".join(missing_protected)
        )
    if not kept:
        # Nothing survived. Stop here, with a message that names the sample and
        # the mode, rather than leave an empty FASTA: Snakemake would count the
        # rule as done, and the run would fail only later (in contigs mode at
        # separate_mitochondria, whose script stops on an empty FASTA), with a
        # message that does not point at the decontamination settings. The
        # decisions file gives the reason each contig went, which tells a policy
        # naming the wrong genus from a sample the taxonomy screen could not
        # place at all. Snakemake deletes the outputs of a failed job, so to
        # read that file, run the script by hand or re-run Snakemake with
        # --keep-incomplete.
        raise ValueError(
            f"No contigs were kept for sample '{args.sample}' with mode '{mode}'. "
            "Check taxonomy assignments and decontamination settings."
        )
    # A two-line summary for the log.
    print(f"Sample {args.sample} decontamination mode: {mode}")
    print(f"Kept {len(kept)} of {len(records)} contig assignments.")


# Run main() when the file is started as a script.
if __name__ == "__main__":
    main()
