#!/usr/bin/env python3
r"""Turn BLAST hits of adapter sequences on the draft contigs into adapter loci, and report them.

Why this is needed
------------------
Sequencing adapters and barcodes that stay in the reads can end up in the assembly, most often
at the contig ends, just outside the telomere. NCBI screens every genome submitted to GenBank
for adapter sequence and asks for it to be removed. FunFlux trims ONT reads before assembly
(rules/shared/05_ont_adapters.smk); this screen checks the result in every mode, including
assemblies made elsewhere (contigs mode) and assemblies of PacBio HiFi or Illumina reads.

What goes in
------------
Rule blast_adapters searched the draft contigs with blastn -task blastn-short (made for short
queries) for: the ONT adapters and barcodes and the PacBio and MGI adapters shipped in
workflow/data/adapters/, the Illumina records of NCBI UniVec, and, for a trimmed ONT sample,
the very sequence that was trimmed from its reads. That last one catches a barcode that no
published list holds.

The rules (fixed; no setting changes them)
------------------------------------------
A hit counts when at least 90 % of its aligned bases are identical and:
  * at a contig end: it lies within 100 bp of the end (the distance NCBI's FCS-adaptor uses to
    call a hit terminal) and is at least 20 bp long. There is no rule on how much of the
    adapter record the hit covers: an adapter is often cut short at the very end of a contig,
    and a barcode left behind after its adapter was trimmed must still count;
  * inside a contig: it is at least 28 bp long, the length of a whole ONT adapter. Shorter
    internal matches turn up by chance in a 40 Mb genome.
Hits on the same contig that overlap (or lie within 10 bp of each other) are merged into one
locus, as several adapter records often match the same stretch.

What comes out
--------------
  --out      one row per locus: contig, contig length, start and end (1-based, inclusive),
             distance to the nearer contig end, position (end or internal), the longest hit
             and its identity, and the adapters that matched. A header line alone when nothing
             was found. Rule contig_summary reads it for its adapter column.
  --out-mqc  one row per sample for the Adapters table of the MultiQC report: for a trimmed ONT
             sample, the share of reads that started with the adapter before and after trimming,
             the published ONT barcodes found in the reads, and what was trimmed (from
             --read-choice and --read-check); for every sample, the contig ends with adapter,
             the loci inside contigs and a verdict.
The verdict never says that a genome is ready for submission: the lists cover the common
adapters only. NCBI's own screen, FCS-adaptor, is the check before a submission.

Usage (rule screen_adapters, workflow/rules/shared/25_contig_summary.smk):
    python adapter_screen.py --blast TSV --contigs FASTA --sample NAME \
        [--read-choice TSV] [--read-check TSV] --out TSV --out-mqc TSV

  --blast        the BLAST table of rule blast_adapters, columns: qseqid sseqid pident length
                 qstart qend sstart send evalue bitscore qlen slen
  --contigs      the draft assembly the hits are on (for the contig names and lengths)
  --read-choice  optional: the choice file of rule find_ont_adapter (ONT samples)
  --read-check   optional: the check file of rule check_ont_trimming (ONT samples)

One [FunFlux] summary line is printed for the console (the rule copies it to the log).
Exit status 0; 1 with a message on stderr when an input file cannot be read or is not in the
expected format (a BLAST line with too few columns, for example).

Uses only the Python standard library. Tests: test_adapter_screen.py (same folder).
"""

import argparse
import sys

END_DISTANCE = 100        # bp: a hit this close to a contig end is an end hit
MIN_END_LENGTH = 20       # bp aligned, for an end hit
MIN_INTERNAL_LENGTH = 28  # bp aligned, for a hit inside a contig
MIN_IDENTITY = 90.0       # percent identical bases
MERGE_GAP = 10            # bp: hits closer than this on one contig form one locus
MAX_NAMES = 5             # adapter names listed per locus; more are counted

LOCUS_COLUMNS = [
    "contig", "contig_length", "start", "end", "distance_to_end", "position",
    "longest_hit", "identity", "adapters",
]


def read_contig_lengths(path):
    """Contig name -> length, in file order, from a FASTA file."""
    lengths = {}
    name = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line)
    return lengths


def read_hits(path, contig_lengths):
    """Read the BLAST table and keep the hits that pass the rules.

    Takes in:  the BLAST table and the contig lengths.
    Does:      for every hit, works out its start and end on the contig (BLAST writes them in
               reverse order for a hit on the minus strand), its distance to the nearer end,
               and whether it is an end hit or an internal one; drops it when it is too short
               or too different for its position.
    Produces:  a list of dicts: contig, start, end, distance, position, length, identity, name.
    Used by:   main().
    """
    hits = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            name, contig, identity, length = fields[0], fields[1], float(fields[2]), int(fields[3])
            start = min(int(fields[6]), int(fields[7]))
            end = max(int(fields[6]), int(fields[7]))
            contig_length = contig_lengths.get(contig, int(fields[11]))
            distance = min(start - 1, contig_length - end)
            if distance <= END_DISTANCE:
                position = "end"
                long_enough = length >= MIN_END_LENGTH
            else:
                position = "internal"
                long_enough = length >= MIN_INTERNAL_LENGTH
            if identity < MIN_IDENTITY or not long_enough:
                continue
            hits.append({
                "contig": contig, "start": start, "end": end, "distance": distance,
                "position": position, "length": length, "identity": identity, "name": name,
            })
    return hits


def merge_loci(hits, contig_lengths):
    """Merge overlapping hits on one contig into loci.

    Takes in:  the kept hits and the contig lengths.
    Does:      sorts the hits by contig (in assembly order) and start, and joins a hit to the
               locus before it when it starts at most MERGE_GAP bp after that locus ends. A
               locus is an end locus when any of its hits is an end hit.
    Produces:  a list of loci (dicts with the LOCUS_COLUMNS).
    Used by:   main().
    """
    order = {name: index for index, name in enumerate(contig_lengths)}
    hits = sorted(hits, key=lambda hit: (order.get(hit["contig"], len(order)), hit["start"]))
    loci = []
    for hit in hits:
        if loci and loci[-1]["contig"] == hit["contig"] and hit["start"] <= loci[-1]["end"] + MERGE_GAP:
            locus = loci[-1]
            locus["end"] = max(locus["end"], hit["end"])
            locus["distance_to_end"] = min(locus["distance_to_end"], hit["distance"])
            if hit["position"] == "end":
                locus["position"] = "end"
            if hit["length"] > locus["longest_hit"]:
                locus["longest_hit"] = hit["length"]
                locus["identity"] = hit["identity"]
            locus["names"].add(hit["name"])
        else:
            loci.append({
                "contig": hit["contig"],
                "contig_length": contig_lengths.get(hit["contig"], "NA"),
                "start": hit["start"],
                "end": hit["end"],
                "distance_to_end": hit["distance"],
                "position": hit["position"],
                "longest_hit": hit["length"],
                "identity": hit["identity"],
                "names": {hit["name"]},
            })
    for locus in loci:
        names = sorted(locus.pop("names"))
        listed = ",".join(names[:MAX_NAMES])
        if len(names) > MAX_NAMES:
            listed += f" and {len(names) - MAX_NAMES} more"
        locus["adapters"] = listed
        locus["identity"] = f"{locus['identity']:.1f}"
    return loci


def read_key_values(path):
    """Read "name<TAB>value" lines (the read-level files of 05_ont_adapters.smk) as a dict."""
    values = {}
    if not path:
        return values
    with open(path) as handle:
        for line in handle:
            name, _, value = line.rstrip("\n").partition("\t")
            values[name] = value
    return values


def trimmed_text(choice):
    """What was trimmed from the reads, for the MultiQC table.

    Takes in:  the choice file of rule find_ont_adapter as a dict (empty for a sample whose reads
               FunFlux did not trim: PacBio, Illumina, contigs mode, or ont_adapter off).
    Produces:  the scan's own one-line summary, for example "adapter + native barcode NB24
               (68 bp)" or "nothing (no adapter found in the reads)".
    """
    if not choice:
        return "reads not trimmed by FunFlux"
    return choice.get("summary", "see the adapter choice file")


def contig_ends_with_adapter(loci):
    """Count the contig ends that carry adapter, as (contig, start or end) pairs.

    Takes in:  the loci of merge_loci().
    Does:      gives every end locus the side it is nearer to (the start or the end of its contig)
               and counts each contig end once: two separate loci at one end, or a short contig with
               a locus near each end, are then counted correctly.
    Used by:   write_mqc() and main().
    """
    ends = set()
    for locus in loci:
        if locus["position"] != "end":
            continue
        from_start = locus["start"] - 1
        from_end = locus["contig_length"] - locus["end"] if locus["contig_length"] != "NA" else from_start + 1
        side = "start" if from_start <= from_end else "end"
        ends.add((locus["contig"], side))
    return len(ends)


def verdict_text(end_loci, internal_loci):
    """The verdict column of the MultiQC table."""
    parts = []
    if end_loci:
        parts.append(f"adapter at {end_loci} contig end{'s' if end_loci != 1 else ''}")
    if internal_loci:
        parts.append(f"{internal_loci} adapter match{'es' if internal_loci != 1 else ''} inside contigs to check")
    if not parts:
        return "none of the listed adapters found"
    return "; ".join(parts)


def write_loci(path, loci):
    """Write the loci table (--out): a header line, then one line per locus.

    The columns are LOCUS_COLUMNS; a header line alone means that no adapter was found.
    Rule contig_summary reads the file for its adapter column.
    Used by:   main().
    """
    with open(path, "w") as out:
        out.write("\t".join(LOCUS_COLUMNS) + "\n")
        for locus in loci:
            out.write("\t".join(str(locus[column]) for column in LOCUS_COLUMNS) + "\n")


def write_mqc(path, sample, loci, choice, check):
    """Write the sample's row of the Adapters table, behind MultiQC's custom-content header lines.

    Rule multiqc lists the file among its inputs; the rows of all samples merge into one table
    because they share the id below and their row names (the sample names) differ.
    """
    end_loci = contig_ends_with_adapter(loci)
    internal_loci = sum(1 for locus in loci if locus["position"] == "internal")
    before = choice.get("percent_with_adapter", "NA") if choice else "NA"
    after = check.get("percent_with_adapter", "NA") if check else "NA"
    barcodes = choice.get("published_barcodes", "NA") if choice else "NA"
    with open(path, "w") as out:
        out.write("# id: adapter_screen\n")
        out.write("# section_name: 'Adapters'\n")
        out.write(
            "# description: 'ONT adapter trimming of the reads (ONT samples only) and the search of "
            "the draft contigs for ONT, PacBio, MGI and Illumina adapters and barcodes. An end hit lies "
            "within 100 bp of a contig end; adapter left there usually sits just outside the telomere. "
            "The lists cover the common adapters only: before a GenBank submission, also run NCBI FCS-adaptor.'\n"
        )
        out.write("# plot_type: 'table'\n")
        out.write("# pconfig: {namespace: 'Adapters'}\n")
        out.write("\t".join([
            "Sample", "Reads with adapter (%)", "Reads with adapter after trimming (%)",
            "ONT barcodes in the reads", "Trimmed from the reads", "Contig ends with adapter",
            "Adapter matches inside contigs", "Verdict",
        ]) + "\n")
        out.write("\t".join([
            sample, before, after, barcodes, trimmed_text(choice), str(end_loci), str(internal_loci),
            verdict_text(end_loci, internal_loci),
        ]) + "\n")


def build_parser():
    """Define the command line; the text at the top of this file says what each option holds.

    --read-choice and --read-check may come without a value: for a sample whose reads were
    not trimmed, Snakemake writes the empty path as nothing, and the option reads as "not given".
    Used by:   main().
    """
    parser = argparse.ArgumentParser(description="Report adapter loci on the draft contigs.")
    parser.add_argument("--blast", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--read-choice", nargs="?", const="", default="")
    parser.add_argument("--read-check", nargs="?", const="", default="")
    parser.add_argument("--out", required=True)
    parser.add_argument("--out-mqc", required=True)
    return parser


def main(argv=None):
    """Read the inputs, merge the hits into loci, write both tables and print one summary line.

    Takes in:  argv, the command line without the script name (None means the real command
               line; the tests pass a list).
    Does:      read_contig_lengths(), read_hits(), read_key_values(), merge_loci(),
               write_loci() and write_mqc(). An input that cannot be read ends the run with
               exit status 1 and one [FunFlux] line on stderr, instead of a Python traceback.
    Used by:   rule screen_adapters (25_contig_summary.smk).
    """
    args = build_parser().parse_args(argv)
    try:
        contig_lengths = read_contig_lengths(args.contigs)
        hits = read_hits(args.blast, contig_lengths)
        choice = read_key_values(args.read_choice)
        check = read_key_values(args.read_check)
    except (OSError, ValueError, IndexError) as error:
        sys.exit(f"[FunFlux] screen_adapters {args.sample}: cannot read an input file: {error}")
    loci = merge_loci(hits, contig_lengths)
    write_loci(args.out, loci)
    write_mqc(args.out_mqc, args.sample, loci, choice, check)
    end_loci = contig_ends_with_adapter(loci)
    internal_loci = sum(1 for locus in loci if locus["position"] == "internal")
    print(
        f"[FunFlux] screen_adapters {args.sample}: {verdict_text(end_loci, internal_loci)} "
        f"({len(contig_lengths)} contigs searched)."
    )


if __name__ == "__main__":
    main()
