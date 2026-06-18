#!/usr/bin/env python3
"""Select assembly contigs using BlobTools genus assignments.

The script reads the best genus assignment for each contig, applies the
configured decontamination strategy, and writes a filtered contig FASTA plus
small audit files explaining what was kept or removed.

Expected inputs:
- a BlobTools bestscore table, where each row links one contig to its best
  taxonomic assignment;
- the assembly FASTA from which contigs should be retained or discarded;
- optional global and sample-specific include/exclude genus settings.

Main idea:
1. Convert the config/TSV settings into Python lists, booleans, and strings.
2. Read BlobTools genus assignments and count how common each genus is.
3. Decide whether each contig should be kept or removed.
4. Write both the filtered FASTA and audit files that make the decision
   transparent.
"""

import argparse
import csv
import sys
from collections import Counter, OrderedDict


# The workflow supports four decontamination modes:
# - auto: keep the most abundant genus inferred from BlobTools;
# - include: keep only user-provided genera;
# - exclude: remove user-provided genera;
# - off: keep everything, useful when decontamination is disabled.
VALID_MODES = {"auto", "include", "exclude", "off"}

# Boolean values may arrive from YAML, TSV, or the command line as strings.
# These sets make the accepted spelling explicit and prevent silent surprises.
FALSE_VALUES = {"false", "0", "no", "off"}
TRUE_VALUES = {"true", "1", "yes", "on"}

# Some bacterial genera are frequently split across related names in BLAST /
# BlobTools output, for example Bacillus/Paenibacillus,
# Arthrobacter/Pseudarthrobacter, or Burkholderia/Paraburkholderia. Curated
# aliases cover high-confidence cases observed in real runs. The prefix tuple
# below preserves a limited part of the old BacFlux empirical regex behavior.
# Together, these are heuristic safeguards against false contig removal, not a
# formal taxonomic reconciliation system.
GENUS_EQUIVALENCE_ALIASES = {
    "paenibacillus": ("bacillus",),
    "pseudarthrobacter": ("arthrobacter",),
    "pseudoarthrobacter": ("arthrobacter",),
    "paenarthrobacter": ("arthrobacter",),
    "paraburkholderia": ("burkholderia",),
}

GENUS_EQUIVALENCE_PREFIXES = (
    "brady",
    "meso",
    "neo",
    "sino",
    "aeri",
    "caldi",
    "geo",
)


# Small normalization helpers keep comparisons robust against empty values,
# extra spaces, and different capitalization in config files or TSV inputs.
def normalize(value):
    """Convert None or any value into a stripped string.

    This keeps downstream code simple: instead of repeatedly checking for None,
    every parser function can work with a clean string.
    """
    return str(value or "").strip()


def normalize_genus(value):
    """Normalize genus names for case-insensitive comparisons.

    BlobTools output and user config should normally agree on capitalization,
    but using lowercase internally avoids fragile matches such as
    "Pseudomonas" versus "pseudomonas".
    """
    return normalize(value).lower()


def genus_aliases(value):
    """Return comparable aliases for a genus name.

    The first alias is the normalized original genus. Additional aliases come
    from curated common reclassification cases plus a small set of legacy
    prefixes retained from the original BacFlux selector. For example:

    - Paenibacillus -> paenibacillus, bacillus
    - Pseudarthrobacter -> pseudarthrobacter, arthrobacter
    - Paenarthrobacter -> paenarthrobacter, arthrobacter
    - Paraburkholderia -> paraburkholderia, burkholderia

    These heuristic aliases are used for auto/include matching only. Exclude
    mode remains exact because over-broad contaminant removal is riskier there.
    """
    genus = normalize_genus(value)
    aliases = {genus} if genus else set()
    aliases.update(GENUS_EQUIVALENCE_ALIASES.get(genus, ()))
    for prefix in GENUS_EQUIVALENCE_PREFIXES:
        if genus.startswith(prefix) and len(genus) > len(prefix) + 3:
            aliases.add(genus[len(prefix):])
    return aliases


def genus_matches(genus, targets):
    """Return True when a genus matches any target directly or by alias."""
    genus_keys = genus_aliases(genus)
    return any(genus_keys & genus_aliases(target) for target in targets)


def parse_bool(value, default=False):
    """Parse config-style boolean values such as true/false, yes/no, 1/0.

    Empty values are allowed and fall back to the provided default. Invalid
    values raise an error because a misspelled boolean could change whether
    no-hit contigs are removed.
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
    """Split a genus list into a Python list of genus strings.

    The config may contain a single genus, a semicolon-separated list, a
    comma-separated list, or a tab-separated TSV cell. Internally all separators
    are converted to tabs first, then empty items are discarded.
    """
    text = normalize(value)
    if not text:
        return []
    for separator in [",", ";"]:
        text = text.replace(separator, "\t")
    return [item.strip() for item in text.split("\t") if item.strip()]


def read_lines_file(path):
    """Read an optional plain-text genus list.

    This is used for external include/exclude files. Each line may contain one
    genus or a small list accepted by split_genera(). Blank lines and lines
    starting with # are ignored so the file can be commented.
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
    """Read sample-specific decontamination settings from a TSV override file.

    The override file is optional. When present, each row can modify the mode,
    include list, exclude list, or discard_no_hit behavior for one sample.

    Important TSV detail:
    empty include/exclude cells are meaningful. They mean "do not replace the
    global list". Therefore the row must still contain the correct number of tab
    separators, even when some fields are blank.
    """
    if not path:
        return {}
    with open(path, "r", encoding="utf-8", newline="") as handle:
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
            if mode and mode not in VALID_MODES:
                raise ValueError(
                    f"Invalid mode '{mode}' for sample '{sample}' in '{path}'."
                )
            # Store parsed list columns as lists, but keep discard_no_hit as the
            # raw string for now so an empty cell can mean "do not override".
            overrides[sample] = {
                "mode": mode,
                "include": split_genera(row.get("include_genera")),
                "exclude": split_genera(row.get("exclude_genera")),
                "discard_no_hit": row.get("discard_no_hit"),
            }
    return overrides


def read_include_by_sample(path, sample):
    """Read optional per-sample include genera from a two-column TSV file.

    This helper supports a compact file with columns sample and genus. It is
    useful when many samples need different include lists but a full override
    table would be unnecessarily verbose.
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


def read_bestscore(path):
    """Parse BlobTools bestscore output and count genus assignments.

    The returned records list contains tuples: (contig_id, assigned_genus).
    A tuple is a small fixed-size grouping of values; here each tuple keeps the
    contig name and its genus together.

    The Counter named counts records how many contigs were assigned to each
    genus. This is later used for the composition summary and for auto mode.
    """
    records = []
    counts = Counter()
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 22:
                continue
            contig = fields[0].strip()

            # BlobTools genus is expected in column 22 in this table format.
            # Python uses zero-based indexing, so column 22 is fields[21].
            # Empty genus values are treated as no-hit/unclassified.
            genus = fields[21].strip() or "no-hit"
            if not contig:
                continue

            # Append a (contig, genus) tuple to preserve row order for the
            # decisions output, and update the per-genus count at the same time.
            records.append((contig, genus))
            counts[genus] += 1
    if not records:
        raise ValueError(f"No contig taxonomy records could be parsed from '{path}'.")
    return records, counts


def read_fasta(path):
    """Read the input contig FASTA while preserving the original record order.

    OrderedDict behaves like a dictionary, but remembers insertion order. The
    keys are FASTA headers and the values are full sequence strings. Preserving
    order keeps the filtered FASTA close to the original assembly layout.
    """
    records = OrderedDict()
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # A new header means the previous record is complete. Join all
                # accumulated sequence chunks before starting the next contig.
                if header is not None:
                    records[header] = "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                # FASTA sequences may span many lines, so collect them in a
                # list and join once. This is faster and cleaner than repeated
                # string concatenation.
                seq_chunks.append(line.strip())
    if header is not None:
        # Store the final record after the loop; there is no next header to
        # trigger the normal save step.
        records[header] = "".join(seq_chunks)
    if not records:
        raise ValueError(f"No FASTA records could be parsed from '{path}'.")
    return records


def fasta_key(header):
    """Match FASTA records by the first header token, as most tools do.

    If a FASTA header is "contig_1 length=1234", BlobTools usually refers only
    to "contig_1". Splitting on whitespace keeps these IDs compatible.
    """
    return header.split()[0]


def is_no_hit(genus):
    """Treat BlobTools no-hit labels as unclassified contigs.

    The string test is intentionally broad because BlobTools labels can vary
    slightly, but they usually contain no-hit when no confident taxonomy was
    assigned.
    """
    return "no-hit" in normalize_genus(genus)


def choose_auto_genus(counts, discard_no_hit):
    """Choose the most abundant assigned genus for auto decontamination mode.

    Auto mode assumes that the dominant assigned genus is the intended organism.
    If discard_no_hit is true, no-hit contigs are excluded from this choice so
    "no-hit" cannot accidentally become the selected target.
    """
    candidates = [
        (genus, count)
        for genus, count in counts.items()
        if not discard_no_hit or not is_no_hit(genus)
    ]
    if not candidates:
        return None

    # Sort primarily by descending count. If two genera have the same count,
    # sort alphabetically to make the result deterministic across runs.
    candidates.sort(key=lambda item: (-item[1], normalize_genus(item[0])))
    return candidates[0][0]


def decide(contig, genus, mode, include, exclude, discard_no_hit, auto_genus):
    """Return the keep/remove decision and a short reason for one contig.

    The result is a two-item tuple: ("keep" or "remove", reason). Writing the
    reason later into contig_taxonomy_decisions.tsv makes the filtering auditable
    instead of hiding decisions inside the code.
    """
    genus_key = normalize_genus(genus)
    include_keys = {normalize_genus(item) for item in include}
    exclude_keys = {normalize_genus(item) for item in exclude}

    # Mode off is an explicit bypass: keep everything, including no-hit contigs.
    if mode == "off":
        return "keep", "mode_off"

    # no-hit removal is checked before mode-specific genus logic. This means
    # discard_no_hit=true removes unclassified contigs in auto/include/exclude
    # modes unless mode is off.
    if is_no_hit(genus) and discard_no_hit:
        return "remove", "discarded_no_hit"

    # Auto mode keeps only contigs assigned to the most abundant genus chosen
    # above. This is a convenience mode for relatively clean single-organism
    # assemblies.
    if mode == "auto":
        if auto_genus is None:
            return "remove", "auto_no_genus"
        if genus_key == normalize_genus(auto_genus):
            return "keep", f"auto_genus:{auto_genus}"
        if genus_matches(genus, [auto_genus]):
            return "keep", f"auto_genus_alias:{auto_genus}"
        return "remove", f"not_auto_genus:{auto_genus}"

    # Include mode is strict: a user-provided target genus list is mandatory,
    # and every contig outside that list is removed.
    if mode == "include":
        if not include_keys:
            raise ValueError("Mode 'include' requires at least one genus.")
        if genus_key in include_keys:
            return "keep", "included_genus"
        if genus_matches(genus, include):
            return "keep", "included_genus_alias"
        return "remove", "not_included_genus"

    # Exclude mode is the inverse: only the listed contaminant genera are
    # removed, while all other assigned genera are retained.
    if mode == "exclude":
        if genus_key in exclude_keys:
            return "remove", "excluded_genus"
        return "keep", "not_excluded_genus"
    raise ValueError(f"Unsupported decontamination mode '{mode}'.")


def write_composition(path, counts, total):
    """Write genus composition as relative frequencies for quick inspection.

    The output is not used for filtering. It is a human-readable summary that
    helps decide whether auto/include/exclude settings make biological sense.
    """
    with open(path, "w", encoding="utf-8") as handle:
        # Sort by decreasing abundance so the dominant assignments appear first.
        for genus, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            handle.write(f"{genus}: {count / total:.2f}\n")


def main():
    # Command-line arguments are supplied by Snakemake from the workflow config.
    # Optional text arguments use nargs="?" and const="" so a command like
    # --include-genera with an empty YAML value is interpreted as an empty string
    # instead of causing argparse to fail with "expected one argument".
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bestscore", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--mode", choices=sorted(VALID_MODES), default="auto")
    parser.add_argument("--include-genera", nargs="?", const="", default="")
    parser.add_argument("--include-genera-by-sample", nargs="?", const="", default="")
    parser.add_argument("--exclude-genera", nargs="?", const="", default="")
    parser.add_argument("--exclude-genera-file", nargs="?", const="", default="")
    parser.add_argument("--sample-overrides", nargs="?", const="", default="")
    parser.add_argument("--discard-no-hit", nargs="?", const="false", default="false")
    parser.add_argument("--output-list", required=True)
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--composition", required=True)
    parser.add_argument("--decisions", required=True)
    args = parser.parse_args()

    # Start from global config values: mode, include/exclude lists, optional
    # genus list files, and whether unclassified/no-hit contigs are discarded.
    # Lists can be built from direct config values plus optional external files.
    mode = args.mode
    include = split_genera(args.include_genera)
    include.extend(read_include_by_sample(args.include_genera_by_sample, args.sample))
    exclude = split_genera(args.exclude_genera)
    exclude.extend(read_lines_file(args.exclude_genera_file))
    discard_no_hit = parse_bool(args.discard_no_hit, default=False)

    # Sample overrides replace the relevant global settings only for that sample.
    # Empty include/exclude override columns leave the global lists untouched.
    # This allows a run-wide exclude list plus special behavior for one sample.
    overrides = read_overrides(args.sample_overrides)
    if args.sample in overrides:
        override = overrides[args.sample]
        mode = override["mode"] or mode
        if override["include"]:
            include = override["include"]
        if override["exclude"]:
            exclude = override["exclude"]
        if override["discard_no_hit"] is not None and normalize(override["discard_no_hit"]):
            discard_no_hit = parse_bool(override["discard_no_hit"], default=discard_no_hit)

    # Load taxonomy assignments and FASTA records, then index FASTA by contig ID
    # so BlobTools rows can be connected back to the original sequences.
    records, counts = read_bestscore(args.bestscore)
    fasta_records = read_fasta(args.contigs)

    # fasta_by_id is a dictionary. The key is the short contig ID used by
    # BlobTools, and the value is a tuple containing the original full FASTA
    # header and sequence. Keeping the full header avoids losing metadata.
    fasta_by_id = {fasta_key(header): (header, sequence) for header, sequence in fasta_records.items()}
    total = sum(counts.values())
    auto_genus = choose_auto_genus(counts, discard_no_hit) if mode == "auto" else None

    # Apply the decision logic contig by contig. The decisions TSV is an audit
    # trail that explains every keep/remove call made by the workflow.
    kept = []
    seen = set()
    with open(args.decisions, "w", encoding="utf-8") as decisions:
        decisions.write("contig\tassigned_genus\taction\treason\n")
        for contig, genus in records:
            action, reason = decide(contig, genus, mode, include, exclude, discard_no_hit, auto_genus)
            decisions.write(f"{contig}\t{genus}\t{action}\t{reason}\n")

            # A contig is added to the output only if it was marked keep, exists
            # in the FASTA, and has not already been added. The seen set prevents
            # duplicate FASTA records if the taxonomy table contains duplicates.
            if action == "keep" and contig in fasta_by_id and contig not in seen:
                kept.append(contig)
                seen.add(contig)

    # The composition report is written after decisions so it is generated even
    # when the selected contig list is small; it summarizes the original input
    # taxonomy, not only the kept contigs.
    write_composition(args.composition, counts, total)

    # Write the list of kept contig IDs for tools that prefer a text selection.
    with open(args.output_list, "w", encoding="utf-8") as handle:
        for contig in kept:
            handle.write(f"{contig}\n")

    # Write the filtered assembly FASTA, wrapping sequences at 80 characters.
    with open(args.output_fasta, "w", encoding="utf-8") as handle:
        for contig in kept:
            header, sequence = fasta_by_id[contig]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")

    # Warn if BlobTools reported contigs that are not present in the input FASTA;
    # this usually indicates mismatched inputs or changed contig headers.
    missing = sorted({contig for contig, _ in records if contig not in fasta_by_id})
    if missing:
        print(
            f"WARNING: {len(missing)} contigs from BlobTools table were not found in FASTA.",
            file=sys.stderr,
        )
    if not kept:
        # An empty selected assembly is almost always a configuration problem,
        # so fail loudly instead of letting downstream rules consume nothing.
        raise ValueError(
            f"No contigs were kept for sample '{args.sample}' with mode '{mode}'. "
            "Check taxonomy assignments and decontamination settings."
        )
    print(f"Sample {args.sample} decontamination mode: {mode}")
    print(f"Kept {len(kept)} of {len(records)} contig assignments.")


if __name__ == "__main__":
    main()
