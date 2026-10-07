#!/usr/bin/env python3
"""Screen assembled contigs for mitochondrial DNA and split confirmed calls.

Short-read implementation: mitochondrial screen v0.1.

FunFlux normally uses direct mitochondrial BLAST evidence only: the short-read
assembly has no long-read assembler metadata to corroborate a call. The optional
Flye parser is retained solely for compatible standalone use. ``split`` preserves
the input assembly while writing nuclear and mitochondrial FASTA files. It never
confirms a contig from circularity or coverage alone.
"""

import argparse
import csv
import os
import re
import statistics
import sys
from collections import OrderedDict


# Match titles that describe an organelle sequence itself.  Do not match the
# much broader word "mitochondrial": nuclear-encoded proteins routinely use it
# in their descriptions and must remain in the nuclear genome.
DIRECT_MITO_TITLE = re.compile(
    r"\bmitochondrion,\s*(?:complete|partial)\s+(?:genome|sequence)\b"
    r"|\bmitochondrial\s+(?:DNA|genome)\b",
    re.IGNORECASE,
)

VALID_MODES = {"off", "report", "separate_confirmed"}


def fasta_id(header):
    """Return the first FASTA-header token used by BLAST and Flye."""
    return header.split()[0]


def ensure_parent(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def read_fasta(path):
    """Read FASTA records in order and reject duplicate short identifiers."""
    records = OrderedDict()
    header = None
    sequence = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    identifier = fasta_id(header)
                    if identifier in records:
                        raise ValueError(f"Duplicate FASTA identifier: {identifier}")
                    records[identifier] = (header, "".join(sequence))
                header = line[1:].strip()
                sequence = []
            elif header is not None:
                sequence.append(line.strip())
    if header is not None:
        identifier = fasta_id(header)
        if identifier in records:
            raise ValueError(f"Duplicate FASTA identifier: {identifier}")
        records[identifier] = (header, "".join(sequence))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def write_fasta(path, records):
    """Write ordered FASTA records with conventional 80-base wrapping."""
    ensure_parent(path)
    with open(path, "w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def read_id_file(path):
    """Read one contig ID per line, accepting blank lines and comments."""
    identifiers = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if value and not value.startswith("#"):
                identifiers.add(value)
    return identifiers


def is_yes(value):
    return str(value or "").strip().upper() in {"Y", "YES", "TRUE", "1"}


def as_float(value, label):
    try:
        return float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {label} value: {value!r}") from error


def read_blast(path):
    """Collect direct mitochondrial-reference evidence from existing blastout.

    The workflow's BLAST format stores qcovs in column 9 and stitle in column
    15 (one-based).  qcovs is already aggregated by BLAST for each subject, so
    the largest direct organelle hit is a practical conservative signal.
    """
    evidence = {}
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 15:
                continue
            contig = row[0].strip()
            title = row[14].strip()
            if not contig or not DIRECT_MITO_TITLE.search(title):
                continue
            qcov = as_float(row[8], "BLAST qcovs")
            previous = evidence.get(contig)
            if previous is None or qcov > previous["best_qcov"]:
                evidence[contig] = {
                    "best_qcov": qcov,
                    "best_title": title,
                    "direct_hits": 1 if previous is None else previous["direct_hits"] + 1,
                }
            else:
                previous["direct_hits"] += 1
    return evidence


def read_flye_info(path):
    """Read Flye circularity, graph flags, and long-read coverage by contig."""
    records = {}
    coverage_values = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            contig = (row.get("#seq_name") or row.get("seq_name") or "").strip()
            if not contig:
                continue
            coverage = as_float(row.get("cov.", row.get("coverage")), "Flye coverage")
            records[contig] = {
                "coverage": coverage,
                "circular": is_yes(row.get("circ.", row.get("circular"))),
                "repeat": is_yes(row.get("repeat")),
                "alt_group": (row.get("alt_group") or "*").strip(),
            }
            coverage_values.append(coverage)
    if not records:
        raise ValueError(f"No Flye records found in {path}")
    return records, statistics.median(coverage_values)


def classify(evidence, flye_record, median_coverage, minimum_qcov, minimum_ratio, mode):
    """Return a conservative class and an audit-friendly explanation."""
    if mode == "off":
        return "not_screened", "screening_disabled", None

    best_qcov = evidence["best_qcov"] if evidence else 0.0
    if flye_record is None:
        if best_qcov >= minimum_qcov:
            return "confirmed_mitochondrial", "strong_direct_blast", None
        if evidence:
            return "mitochondrial_candidate", "weak_direct_blast", None
        return "not_mitochondrial", "no_mitochondrial_evidence", None

    coverage_ratio = flye_record["coverage"] / median_coverage
    clean_circular = (
        flye_record["circular"]
        and not flye_record["repeat"]
        and flye_record["alt_group"] == "*"
    )
    high_coverage = coverage_ratio >= minimum_ratio

    if best_qcov >= minimum_qcov and (high_coverage or clean_circular):
        support = "high_coverage" if high_coverage else "clean_circularity"
        return "confirmed_mitochondrial", f"direct_blast_plus_{support}", coverage_ratio
    # A near-complete direct organellar match remains decisive when no Flye
    # metadata was supplied, as is normally the case for FunFlux short reads.
    if best_qcov >= 70.0:
        return "confirmed_mitochondrial", "very_strong_direct_blast", coverage_ratio
    if evidence:
        return "mitochondrial_candidate", "direct_blast_without_support", coverage_ratio
    if flye_record["circular"] and high_coverage:
        return "mitochondrial_candidate", "circular_and_high_coverage_without_blast", coverage_ratio
    if flye_record["circular"]:
        return "circular_unclassified", "circularity_only", coverage_ratio
    if high_coverage:
        return "high_coverage_unclassified", "coverage_only", coverage_ratio
    return "not_mitochondrial", "no_mitochondrial_evidence", coverage_ratio


def screen(args):
    contigs = read_fasta(args.contigs)
    evidence_by_contig = read_blast(args.blast)
    flye_by_contig = None
    median_coverage = None
    if args.flye_info:
        flye_by_contig, median_coverage = read_flye_info(args.flye_info)

    ensure_parent(args.report)
    ensure_parent(args.confirmed_ids)
    # Keep biological classifications separate from downstream selection: report
    # mode must be informative without silently changing the nuclear assembly.
    confirmed = []
    with open(args.report, "w", encoding="utf-8", newline="") as report:
        writer = csv.writer(report, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "contig", "length", "classification", "reason", "best_mito_qcovs",
            "direct_mito_hits", "best_mito_title", "flye_circular", "flye_coverage",
            "flye_coverage_ratio", "flye_repeat", "flye_alt_group",
        ])
        for contig, (_, sequence) in contigs.items():
            evidence = evidence_by_contig.get(contig)
            flye_record = flye_by_contig.get(contig) if flye_by_contig else None
            classification, reason, coverage_ratio = classify(
                evidence,
                flye_record,
                median_coverage,
                args.min_direct_blast_qcov,
                args.min_flye_coverage_ratio,
                args.mode,
            )
            if classification == "confirmed_mitochondrial":
                confirmed.append(contig)
            writer.writerow([
                contig,
                len(sequence),
                classification,
                reason,
                f"{evidence['best_qcov']:.2f}" if evidence else "",
                evidence["direct_hits"] if evidence else 0,
                evidence["best_title"] if evidence else "",
                "Y" if flye_record and flye_record["circular"] else "N" if flye_record else "",
                f"{flye_record['coverage']:.2f}" if flye_record else "",
                f"{coverage_ratio:.2f}" if coverage_ratio is not None else "",
                "Y" if flye_record and flye_record["repeat"] else "N" if flye_record else "",
                flye_record["alt_group"] if flye_record else "",
            ])

    # Report mode is non-destructive: it records calls but deliberately leaves
    # the downstream assembly unchanged by writing an empty selection list.
    selected = confirmed if args.mode == "separate_confirmed" else []
    with open(args.confirmed_ids, "w", encoding="utf-8") as handle:
        for contig in selected:
            handle.write(f"{contig}\n")

    print(
        f"Mitochondrial screen: {len(confirmed)} confirmed call(s); "
        f"{len(selected)} selected for separation (mode={args.mode})."
    )


def split(args):
    contigs = read_fasta(args.fasta)
    confirmed = read_id_file(args.confirmed_ids)
    # The split can receive a selected or polished assembly, so fail explicitly
    # if a confirmed ID is no longer present rather than mislabelling output.
    missing = sorted(confirmed - set(contigs))
    if missing:
        raise ValueError(
            "Confirmed mitochondrial contig(s) are absent from final assembly: "
            + ", ".join(missing)
        )

    nuclear = []
    mitochondrial = []
    for contig, record in contigs.items():
        if contig in confirmed:
            mitochondrial.append(record)
        else:
            nuclear.append(record)
    if not nuclear:
        raise ValueError("Mitochondrial split would leave an empty nuclear assembly.")

    write_fasta(args.nuclear, nuclear)
    write_fasta(args.mitochondrial, mitochondrial)
    print(
        f"Mitochondrial split: {len(mitochondrial)} mitochondrial and "
        f"{len(nuclear)} nuclear contig(s)."
    )


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    screen_parser = commands.add_parser("screen", help="classify mitochondrial candidates")
    screen_parser.add_argument("--blast", required=True)
    screen_parser.add_argument("--contigs", required=True)
    screen_parser.add_argument("--report", required=True)
    screen_parser.add_argument("--confirmed-ids", required=True)
    screen_parser.add_argument("--flye-info", default="")
    screen_parser.add_argument("--mode", choices=sorted(VALID_MODES), default="separate_confirmed")
    screen_parser.add_argument("--min-direct-blast-qcov", type=float, default=35.0)
    screen_parser.add_argument("--min-flye-coverage-ratio", type=float, default=3.0)
    screen_parser.set_defaults(func=screen)

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
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
