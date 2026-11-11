#!/usr/bin/env python3
r"""Find the Oxford Nanopore adapter, and the barcode after it, at the start of ONT reads; check trimmed reads.

Why this is needed
------------------
An ONT read starts with the sequencing adapter that brought the DNA molecule to the pore. In a
barcoded library the adapter is followed by the barcode and its two short flanking sequences.
Dorado removes both when it is run by hand (`dorado basecaller`, kit14 kits), but reads basecalled
live in MinKNOW, or with `dorado duplex`, usually keep them, and nothing in the FASTQ file says
which. Left in the reads, they end up at the contig ends, just outside the telomere: in the
hifiasm assembly of the public Botrytis cinerea B05.10 reads (SRR40741277), every contig end with a
telomere carried about 80 bp of adapter and native barcode beyond it.

fastplong removes such a sequence well when it is told what the sequence is. Left to guess, it
looks only at the first reads of the file and takes the most common sequence at the read ends,
which in fungal reads can be the telomere repeat (TTAGGG): it then cuts the telomeres away. This
script therefore works the sequence out itself, from all the reads, and gives fastplong the result.

How the sequence is found (subcommand scan)
-------------------------------------------
1. The anchor. The first 150 bp of every read are searched for the 20 bp core of the ONT adapter
   top strand, TTCGTTCAGTTACGTATTGC. The whole adapter is CCTGTACTTCGTTCAGTTACGTATTGCT; the
   ligation, native and rapid adapters share it (ONT Chemistry Technical Document, section 15). A
   read holding the core exactly is "anchored". Reads whose adapter carries a basecalling error in
   the core are missed, so the share of anchored reads is a lower bound.
2. What follows the adapter. In a barcoded library every anchored read continues with a barcode
   block (native barcoding: AAGGTTAA, the 24 bp barcode, CAGCACCT: 40 bp). In a library without
   barcodes it continues straight into the genome, which differs from read to read. For block
   lengths of 8, 12, ... 120 bp, the script finds the most common sequence that follows the
   adapter, and keeps the longest length whose most common sequence is shared by at least 60 % of
   the anchored reads. On B05.10 the 40 bp native block was shared by 87 % of the anchored reads,
   and the most common 44 bp sequence by 1 %.
3. ONT's published barcodes. The sequence after the adapter is also compared with the native
   barcodes NB01-96 and the rapid barcodes RB01-96 that ONT publishes in its Chemistry Technical
   Document (section 14; the rule passes FunFlux's copy, workflow/data/adapters/ont_chtd.fasta).
   A read counts for a barcode when the 24 bases where the barcode belongs differ from it at 3
   positions or fewer. This names the barcode of the sample ("native barcode NB24" on B05.10) and
   finds every barcode of a file that holds several. The list only sorts the reads by barcode:
   the block trimmed for each barcode is again the one found in its own reads (step 2 on those
   reads alone), and the published block is used only when its reads share none. ONT publishes
   the rapid barcodes only in part, so the block in the reads is the safer source. A kit or
   barcode that is not in the list is trimmed all the same, as found in the reads.
4. The decision. fastplong gets the most common sequence with -s, on one line of options, and
   every sequence to trim, in both orientations, in a FASTA file it reads with -a:
   * fewer than 1 % of reads anchored, and fewer than 1 % starting with a barcode block either:
     the reads are already trimmed (or are not ONT reads). fastplong runs with -A, which copies
     the reads unchanged.
   * at least 1 % anchored, with published barcodes in the reads (each in at least 1 % of the
     anchored reads, and at least 10 reads): the adapter followed by each barcode's block. With
     two or more barcodes above 10 % a NOTE says so, because such a file may hold several
     samples; a barcode below 10 % is usually a few reads of a neighbouring sample assigned to
     this one, and it is trimmed without a note. fastplong also searches the read END for the
     reverse complement of every sequence, which is how the barcode and adapter sit at the 3' end
     of a native-barcoded read.
   * at least 1 % anchored, no published barcode: the adapter followed by the block of step 2,
     when there is one, or the adapter alone (a library without barcodes).
   * fewer than 1 % anchored, but at least 1 % of reads start with a native or rapid barcode
     block with no adapter in front (the adapter was trimmed earlier, the barcode left): each such
     block found in at least 1 % of the reads.
   * a sequence written in the configuration (parameters.<mode>.ont_adapter): that sequence.
   fastplong also splits a read where it finds a sequence in the middle (two molecules joined
   into one read). It allows the same share of mismatches there as at the ends, so a short
   sequence matches genome by chance: with the 28 bp adapter alone it cut 1,041 B05.10 reads,
   1,001 of them in error. Splitting is therefore switched off (--disable_chimera_splitting) when
   any sequence to trim is shorter than 50 bp. With the 68 bp adapter + barcode it cut 73 reads,
   all of them real joins.
5. Safety checks. A block that is a tandem repeat of a short unit (10 bp or less, a telomere for
   example) is never added. When the anchored reads continue with two or more different
   sequences that are NOT published barcodes, only the part they share is trimmed, with a NOTE:
   those barcodes then stay in the reads, and the contig adapter screen reports them. A block
   with the flanks of a native or rapid barcode around a barcode that is not in the list is
   trimmed as found, with a NOTE.

The check (subcommand check)
----------------------------
The reads Filtlong kept from the trimmed reads, the reads the assembler gets, are searched again:
for the adapter core in their first 150 bp and, when a block was trimmed, for the last 20 bp of
the most common trimmed sequence. More than 1 % of reads still starting with the adapter prints a
NOTE; it does not stop the run.

Usage (rules find_ont_adapter and check_ont_trimming, workflow/rules/shared/05_ont_adapters.smk):
    python find_ont_adapter.py scan --reads FASTQ --setting detect|SEQUENCE --sample NAME \
        [--barcode-list FASTA] --out-scan TSV --out-choice TSV --out-options TXT \
        --out-sequence FASTA
    python find_ont_adapter.py check --reads TRIMMED_FASTQ --choice TSV --sample NAME --out TSV

  --reads          the sample's raw ONT reads (scan), or the reads Filtlong kept after fastplong
                   (check); FASTQ, plain or gzip-compressed (.gz)
  --setting        detect, or the sequence written in parameters.<mode>.ont_adapter
  --barcode-list   ONT's published barcodes (workflow/data/adapters/ont_chtd.fasta); without it,
                   no barcode is named and step 4 never trims several barcodes
  --out-scan       one row per 10,000 reads: how many start with the adapter, and how many with
                   a barcode block and no adapter; a batch of trimmed reads in an otherwise
                   untrimmed file shows up here
  --out-choice     the decision and the numbers behind it, one "name<TAB>value" line each
  --out-options    the fastplong options, one line: "-s <most common sequence>", with
                   "--disable_chimera_splitting" added when a sequence is shorter than 50 bp, or
                   "-A" when nothing is trimmed
  --out-sequence   every sequence to trim, each in both orientations (fastplong reads it with -a,
                   and -a does not add the reverse complement itself), or an empty file; the
                   contig adapter screen adds it to its list of adapters
  --choice         (check) the --out-choice file of the scan
  --out            (check) "name<TAB>value" lines: reads, and how many still carry the adapter

Lines starting with [FunFlux] are printed for the console (the rules copy them from the log).
Exit status 0; 1 with a message on stderr when a file cannot be read, or when --setting is
neither detect nor a sequence of A, C, G and T only.

00_common.smk imports this file too, for is_tandem_repeat(), clean_sequence(), ADAPTER_CORE and
MIN_SPLIT_LENGTH, to check a sequence given in the configuration while Snakemake reads the workflow.
Uses only the Python standard library. Tests: test_find_ont_adapter.py (same folder).
"""

import argparse
import collections
import gzip
import random
import re
import sys

# The ONT adapter top strand and its 20 bp core (ONT Chemistry Technical Document, section 15).
ADAPTER = "CCTGTACTTCGTTCAGTTACGTATTGCT"
ADAPTER_CORE = "TTCGTTCAGTTACGTATTGC"
# Where the adapter ends, counted from the start of the core: the core is followed by one more
# adapter base (T) before the block or the genome begins.
ADAPTER_TAIL_AFTER_CORE = len(ADAPTER) - (ADAPTER.index(ADAPTER_CORE) + len(ADAPTER_CORE))

# The flanks around the barcodes (ONT Chemistry Technical Document, section 14). Native:
# AAGGTTAA, 24 bp barcode, CAGCACCT (40 bp). Rapid: GCTTGGGTGTTTAACC, 24 bp barcode, then a 50 bp
# constant sequence (90 bp). The two patterns find a barcode block with no adapter in front; for
# the rapid block only the first 16 bp of the 50 bp sequence are required, and the whole 50 bp
# are added to the sequence to trim.
NATIVE_BLOCK = re.compile(r"AAGGTTAA[ACGT]{24}CAGCACCT")
RAPID_TAIL = "GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
RAPID_BLOCK = re.compile(r"GCTTGGGTGTTTAACC[ACGT]{24}" + RAPID_TAIL[:16])
# Where the 24 bp barcode starts inside each kind of block: after the 8 bp native flank, after
# the 16 bp rapid flank.
BARCODE_LENGTH = 24
BARCODE_START = {"native": 8, "rapid": 16}

START_WINDOW = 150            # bp at the read start searched for the adapter or a block
BLOCK_LENGTHS = range(8, 121, 4)
MIN_READ_SHARE = 0.01         # below this share of reads, the adapter (or block) counts as absent
MIN_BLOCK_SHARE = 0.60        # share of anchored reads one block must reach to be trimmed
SECOND_BLOCK_SHARE = 0.10     # a second block this common means several barcodes in one file
MIN_BARCODE_SHARE = 0.01      # share of anchored reads a published barcode needs to be trimmed
MIN_BARCODE_READS = 10        # ... and at least this many reads
MAX_BARCODE_MISMATCHES = 3    # positions in which a read's barcode may differ from the published one
MIN_SPLIT_LENGTH = 50         # shorter sequences are trimmed at the ends only
MAX_TANDEM_UNIT = 10          # a block repeating a unit this short or shorter is refused
STORED_FOLLOWERS = 100_000    # how many "what follows the adapter" sequences are kept for steps 2-3
READS_PER_ROW = 10_000        # reads per row of the --out-scan table
CHECK_TAIL_LENGTH = 20        # bp at the end of the trimmed sequence searched for by the check
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    """Return the sequence of the other DNA strand, read 5' to 3' (A, C, G, T, upper case).

    Used by:   write_sequences(), which writes every sequence to trim in both orientations,
               because fastplong's -a does not add the reverse complement itself.
    """
    return sequence.translate(COMPLEMENT)[::-1]


def clean_sequence(text):
    """Upper-case a sequence and remove spaces; return None if anything but A, C, G, T is left.

    Takes in:  text, as written in the configuration.
    Used by:   scan() for --setting, and 00_common.smk (parameters.<mode>.ont_adapter).
    """
    sequence = re.sub(r"\s+", "", str(text)).upper()
    if not sequence or not re.fullmatch(r"[ACGT]+", sequence):
        return None
    return sequence


def is_tandem_repeat(sequence, max_unit=MAX_TANDEM_UNIT):
    """True if the sequence is a short unit (max_unit bp or less) repeated end to end.

    TTAGGGTTAGGGTTAGGG is TTAGGG three times and returns True; an adapter or a barcode does not.
    A partial copy at the end counts too (TTAGGGTTAGGGTT), but the unit must fit at least twice,
    so a sequence is never a repeat of itself. Sequences of 8 bp or less are never called a
    repeat, so that a short block is judged by the other rules.
    Used by:   choose_block() and 00_common.smk.
    """
    if len(sequence) <= 8:
        return False
    # The unit must fit at least twice. Otherwise a 9 or 10 bp sequence would count as one
    # copy of a 9 or 10 bp "unit", and a short block that the scan finds after the adapter
    # would be refused as a repeat.
    longest_unit = min(max_unit, len(sequence) // 2)
    for unit_length in range(1, longest_unit + 1):
        unit = sequence[:unit_length]
        repeated = unit * (len(sequence) // unit_length + 1)
        if repeated[:len(sequence)] == sequence:
            return True
    return False


def read_sequences(path):
    """Yield the sequence line of every FASTQ record (the second of its four lines).

    Takes in:  path, a FASTQ file; gzip-compressed when the name ends in .gz.
    Does:      reads the file line by line, so a file of any size fits in memory.
    Used by:   scan() and check().
    """
    if path.endswith(".gz"):
        handle = gzip.open(path, "rt")
    else:
        handle = open(path)
    with handle:
        for line_number, line in enumerate(handle):
            if line_number % 4 == 1:
                yield line.strip().upper()


# ======================================================================
# ONT'S PUBLISHED BARCODES (step 3)
# ======================================================================

def read_barcode_list(path):
    """Read the native and rapid barcode blocks of ONT's published list.

    Takes in:  path, workflow/data/adapters/ont_chtd.fasta (or "" for no list). Its records are
               named ONT_native_barcode_NB01 ... NB96 (AAGGTTAA, barcode, CAGCACCT) and
               ONT_rapid_barcode_RB01 ... RB96 (GCTTGGGTGTTTAACC, barcode, 50 bp sequence); the
               adapter records in the same file are skipped.
    Produces:  a list of dicts, one per barcode: name (NB24), kind (native or rapid), block (the
               barcode with its flanks, the sequence that follows the adapter in a read) and
               barcode (the 24 bp alone).
    Used by:   scan().
    """
    barcodes = []
    if not path:
        return barcodes
    name = None
    parts = []
    records = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts).upper()))

    for record_name, block in records:
        match = re.fullmatch(r"ONT_(native|rapid)_barcode_([A-Z]+\d+)", record_name)
        if not match:
            continue
        kind, barcode_name = match.group(1), match.group(2)
        start = BARCODE_START[kind]
        barcodes.append({
            "name": barcode_name,
            "kind": kind,
            "block": block,
            "barcode": block[start:start + BARCODE_LENGTH],
        })
    return barcodes


def mismatches(first, second):
    """The number of positions at which two sequences of the same length differ."""
    return sum(1 for a, b in zip(first, second) if a != b)


def count_barcodes(followers, barcodes):
    """Sort the followers by the published barcode their block holds.

    Takes in:  followers, the sequences that follow the adapter in a sample of the anchored reads;
               barcodes, the list of read_barcode_list().
    Does:      looks at the 24 bases where a native barcode would sit (after 8 bases of flank) and
               where a rapid barcode would sit (after 16), and finds the published barcode that
               differs from them at MAX_BARCODE_MISMATCHES positions or fewer (the closest one).
               To avoid comparing every read with all 192 barcodes, it first cuts the 24 bases
               into four pieces of 6 and looks each piece up in a table of the barcodes' own
               pieces at the same place. A read with up to 3 differences keeps at least one of
               its four pieces unchanged, so every barcode within the tolerance is found this way.
    Produces:  a dict: barcode name -> the followers that carry it (scan() counts them and finds
               each barcode's block in them). A follower that matches no barcode is left out.
    Used by:   scan().
    """
    piece_length = BARCODE_LENGTH // 4
    by_piece = collections.defaultdict(list)
    for barcode in barcodes:
        sequence = barcode["barcode"]
        for offset in range(0, BARCODE_LENGTH, piece_length):
            by_piece[(barcode["kind"], offset, sequence[offset:offset + piece_length])].append(barcode)

    followers_of = collections.defaultdict(list)
    for follower in followers:
        best = None
        best_mismatches = MAX_BARCODE_MISMATCHES + 1
        for kind, start in BARCODE_START.items():
            piece = follower[start:start + BARCODE_LENGTH]
            if len(piece) < BARCODE_LENGTH:
                continue
            for offset in range(0, BARCODE_LENGTH, piece_length):
                for barcode in by_piece.get((kind, offset, piece[offset:offset + piece_length]), []):
                    differences = mismatches(piece, barcode["barcode"])
                    if differences < best_mismatches:
                        best = barcode
                        best_mismatches = differences
        if best is not None:
            followers_of[best["name"]].append(follower)
    return followers_of


def name_block(block, barcodes):
    """Name the published barcode a block holds, for example "native barcode NB24", or return "".

    Takes in:  a barcode block as found in the reads (native: AAGGTTAA, barcode, CAGCACCT; rapid:
               GCTTGGGTGTTTAACC, barcode, 50 bp sequence) and the list of read_barcode_list().
    Does:      tells the kind of block from its flanks first, and then compares the 24 bases where
               the barcode of that kind sits with the barcodes of that kind only (up to
               MAX_BARCODE_MISMATCHES differences). The kind must come first: ONT uses the same 24 bp
               sequences for many native and rapid barcodes (RB13-RB96 equal NB13-NB96), and only
               the flanks tell them apart.
    Used by:   scan(), for the barcode blocks found with no adapter in front.
    """
    if NATIVE_BLOCK.match(block):
        kind = "native"
    elif RAPID_BLOCK.match(block):
        kind = "rapid"
    else:
        return ""
    start = BARCODE_START[kind]
    piece = block[start:start + BARCODE_LENGTH]
    for barcode in barcodes:
        if barcode["kind"] == kind and mismatches(piece, barcode["barcode"]) <= MAX_BARCODE_MISMATCHES:
            return f"{kind} barcode {barcode['name']}"
    return ""


def looks_like_barcode_block(block):
    """True if the block carries the flanks of a native or rapid barcode around 24 bases."""
    return bool(NATIVE_BLOCK.search(block) or RAPID_BLOCK.search(block))


# ======================================================================
# THE BLOCK THAT FOLLOWS THE ADAPTER (step 2)
# ======================================================================

def choose_block(followers):
    """Find the block that follows the adapter in most anchored reads (step 2 of the method).

    Takes in:  followers, the sequences that follow the adapter (up to 120 bp each) in a sample of
               the anchored reads.
    Does:      for each length in BLOCK_LENGTHS, counts the most common sequence of that length
               and its share of all followers (a follower shorter than the length cannot match).
               The share can only fall as the length grows, because a longer sequence is an
               extension of a shorter one, so the search stops at the first length below
               MIN_BLOCK_SHARE. That next length also tells whether the reads continue with two
               or more different sequences of at least SECOND_BLOCK_SHARE each (several barcodes).
    Produces:  (block, share, several_blocks): the block ("" when none reaches the threshold),
               the share of followers that carry it, and True when several different blocks
               continue after it. A block that is a tandem repeat is returned as "".
    Used by:   scan().
    """
    block = ""
    block_share = 0.0
    several_blocks = False
    if not followers:
        return block, block_share, several_blocks
    for length in BLOCK_LENGTHS:
        counts = collections.Counter(f[:length] for f in followers if len(f) >= length)
        if not counts:
            break
        ranked = counts.most_common(2)
        top_share = ranked[0][1] / len(followers)
        if top_share >= MIN_BLOCK_SHARE:
            block = ranked[0][0]
            block_share = top_share
            continue
        # The first length below the threshold: do two different sequences continue here?
        if len(ranked) == 2 and ranked[1][1] / len(followers) >= SECOND_BLOCK_SHARE:
            several_blocks = True
        break
    # The lengths go in steps of 4 bp, so the block can end up to 3 bp short of its true end (a
    # 90 bp rapid block would stop at 88). Try the next three lengths one base at a time.
    if block:
        for length in range(len(block) + 1, len(block) + 4):
            counts = collections.Counter(f[:length] for f in followers if len(f) >= length)
            if not counts:
                break
            longer, count = counts.most_common(1)[0]
            if count / len(followers) < MIN_BLOCK_SHARE:
                break
            block, block_share = longer, count / len(followers)
    if block and is_tandem_repeat(block):
        return "", 0.0, several_blocks
    return block, block_share, several_blocks


# ======================================================================
# FILES
# ======================================================================

def write_key_values(path, pairs):
    """Write "name<TAB>value" lines, in the given order."""
    with open(path, "w") as out:
        for name, value in pairs:
            out.write(f"{name}\t{value}\n")


def read_key_values(path):
    """Read the "name<TAB>value" lines that write_key_values() wrote, as a dict."""
    values = {}
    with open(path) as handle:
        for line in handle:
            name, _, value = line.rstrip("\n").partition("\t")
            values[name] = value
    return values


def write_sequences(path, sample, sequences):
    """Write the sequences to trim, each in both orientations, as FASTA (an empty file for none).

    Takes in:  sequences, a list of (sequence, label) pairs, the most common first.
    Does:      writes <sample>_trim_<n> and <sample>_trim_<n>_reverse_complement for each one; the
               label (for example "adapter + native barcode NB24") follows the name in the header.
    Used by:   scan().
    """
    with open(path, "w") as out:
        for number, (sequence, label) in enumerate(sequences, start=1):
            out.write(f">{sample}_trim_{number} {label}\n{sequence}\n")
            out.write(f">{sample}_trim_{number}_reverse_complement {label}\n{reverse_complement(sequence)}\n")


# ======================================================================
# THE TWO SUBCOMMANDS
# ======================================================================

def read_the_reads(path):
    """Read every read once and collect what the decision needs (steps 1 and 2).

    Takes in:  path, the raw ONT reads.
    Does:      for every read, looks for the adapter core in the first START_WINDOW bp (an
               anchored read) and, if it is not there, for a native or rapid barcode block. The
               sequence after the adapter of an anchored read goes into a sample of at most
               STORED_FOLLOWERS sequences drawn evenly from the whole file (reservoir sampling
               with a fixed seed, so that every run gives the same result); the barcode blocks
               without adapter are counted.
    Produces:  a dict: reads, anchored, block_only (counts), followers (the sample), block_counts
               (Counter of barcode blocks without adapter) and rows (the --out-scan rows).
    Used by:   scan().
    """
    sampler = random.Random(20261109)
    found = {"reads": 0, "anchored": 0, "block_only": 0, "followers": [],
             "block_counts": collections.Counter(), "rows": []}
    row_reads = row_anchored = row_block_only = 0

    for sequence in read_sequences(path):
        found["reads"] += 1
        row_reads += 1
        start = sequence[:START_WINDOW]
        core_position = start.find(ADAPTER_CORE)
        if core_position != -1:
            found["anchored"] += 1
            row_anchored += 1
            adapter_end = core_position + len(ADAPTER_CORE) + ADAPTER_TAIL_AFTER_CORE
            follower = sequence[adapter_end:adapter_end + max(BLOCK_LENGTHS)]
            # Reservoir sampling: keep the first STORED_FOLLOWERS, then replace a kept one with
            # falling probability, so that every anchored read had the same chance to be kept.
            if len(found["followers"]) < STORED_FOLLOWERS:
                found["followers"].append(follower)
            else:
                slot = sampler.randrange(found["anchored"])
                if slot < STORED_FOLLOWERS:
                    found["followers"][slot] = follower
        else:
            native = NATIVE_BLOCK.search(start)
            rapid = RAPID_BLOCK.search(start)
            if native:
                found["block_only"] += 1
                row_block_only += 1
                found["block_counts"][native.group(0)] += 1
            elif rapid:
                found["block_only"] += 1
                row_block_only += 1
                found["block_counts"][rapid.group(0)[:-16] + RAPID_TAIL] += 1
        if row_reads == READS_PER_ROW:
            found["rows"].append((found["reads"] - row_reads + 1, found["reads"], row_reads, row_anchored, row_block_only))
            row_reads = row_anchored = row_block_only = 0
    if row_reads:
        found["rows"].append((found["reads"] - row_reads + 1, found["reads"], row_reads, row_anchored, row_block_only))
    return found


def scan(args):
    """Subcommand scan: read every read, decide what fastplong removes, write the four files.

    Takes in:  args.reads (the raw ONT reads), args.setting (detect, or a given sequence),
               args.barcode_list (ONT's published barcodes, or ""), args.sample, and the four
               output paths.
    Does:      read_the_reads(), then the decision of step 4 of the method.
    Produces:  the scan table, the choice file, the options line and the sequence FASTA; one
               [FunFlux] summary line and, if needed, [FunFlux] NOTE lines on standard output.
    Used by:   main().
    """
    given = None
    if args.setting.strip().lower() != "detect":
        given = clean_sequence(args.setting)
        if given is None:
            sys.exit(f"[FunFlux] find_ont_adapter {args.sample}: the given sequence is not A, C, G, T only: {args.setting!r}")

    barcodes = read_barcode_list(args.barcode_list)
    found = read_the_reads(args.reads)
    reads = found["reads"]
    followers = found["followers"]
    anchored_share = found["anchored"] / reads if reads else 0.0
    block_only_share = found["block_only"] / reads if reads else 0.0
    block, block_share, several_blocks = choose_block(followers)

    # Step 3: the published barcodes among the anchored reads, most common first, and those
    # common enough to be trimmed.
    barcode_by_name = {barcode["name"]: barcode for barcode in barcodes}
    followers_of = count_barcodes(followers, barcodes)
    present = []
    for name in sorted(followers_of, key=lambda name: len(followers_of[name]), reverse=True):
        count = len(followers_of[name])
        share = count / len(followers)
        if count >= MIN_BARCODE_READS and share >= MIN_BARCODE_SHARE:
            present.append((name, share))
    barcodes_text = ", ".join(f"{name} {100 * share:.1f} %" for name, share in present) or "none found"

    # Step 4: the sequences to trim, as (sequence, label) pairs, the most common first, and the
    # NOTEs for the console.
    to_trim = []
    notes = []
    if given:
        to_trim = [(given, "the sequence given in the configuration")]
    elif anchored_share >= MIN_READ_SHARE and present:
        # Published barcodes in the reads: the adapter followed by each barcode's block, found in
        # the reads of that barcode alone (step 2 on them). When those reads share no block that
        # holds the barcode, the published block is used instead.
        for name, share in present:
            barcode = barcode_by_name[name]
            own_block = choose_block(followers_of[name])[0]
            # Keep the block of the reads when it holds this barcode where the barcode belongs,
            # allowing the same few differences as when the reads were counted.
            start = BARCODE_START[barcode["kind"]]
            own_barcode = own_block[start:start + BARCODE_LENGTH]
            if len(own_barcode) < BARCODE_LENGTH or mismatches(own_barcode, barcode["barcode"]) > MAX_BARCODE_MISMATCHES:
                own_block = barcode["block"]
            to_trim.append((ADAPTER + own_block, f"adapter + {barcode['kind']} barcode {name}"))
        common = [name for name, share in present if share >= SECOND_BLOCK_SHARE]
        if len(common) >= 2:
            notes.append(
                f"the reads carry {len(common)} barcodes in at least {SECOND_BLOCK_SHARE:.0%} of the reads "
                f"with the adapter ({barcodes_text}); all of them are trimmed. If they are different samples, "
                "demultiplex the reads and run each sample on its own."
            )
    elif anchored_share >= MIN_READ_SHARE:
        # No published barcode: the adapter followed by the common block, if there is one.
        if block:
            to_trim = [(ADAPTER + block, f"adapter + the {len(block)} bp block after it")]
            # Barcode flanks around a barcode that the list does not hold (only said when the
            # list was given: without it, every barcode would look unknown).
            if barcodes and looks_like_barcode_block(block):
                notes.append(
                    "the barcode after the adapter is not in ONT's published list; it is trimmed as "
                    "found in the reads."
                )
        else:
            to_trim = [(ADAPTER, "adapter (no common block after it)")]
        if several_blocks:
            notes.append(
                "the reads continue with several different sequences after the adapter that are not "
                "published ONT barcodes; only the part they share is trimmed, so those sequences stay in "
                "the reads and the contig adapter screen will report them. Write the whole sequence into "
                "parameters.<mode>.ont_adapter to trim it."
            )
    elif block_only_share >= MIN_READ_SHARE:
        # Barcode blocks with no adapter in front: each block found in at least 1 % of the reads.
        # The barcodes reported are then these, counted over all reads (no read has the adapter).
        named = []
        for found_block, count in found["block_counts"].most_common():
            if count >= MIN_BARCODE_READS and count / reads >= MIN_READ_SHARE:
                label = name_block(found_block, barcodes)
                if label:
                    named.append(f"{label.split()[-1]} {100 * count / reads:.1f} % of all reads, no adapter")
                else:
                    label = "a barcode not in ONT's published list"
                to_trim.append((found_block, f"{label}, with no adapter in front"))
        if named:
            barcodes_text = ", ".join(named)

    # The fastplong options: -s with the most common sequence (fastplong searches the read end
    # for its reverse complement), and --disable_chimera_splitting when any sequence to trim is
    # short; -A when there is nothing to trim. The rule always adds -a with the sequence FASTA.
    if to_trim:
        main_sequence, main_label = to_trim[0]
        shortest = min(len(sequence) for sequence, label in to_trim)
        if shortest < MIN_SPLIT_LENGTH:
            options = f"-s {main_sequence} --disable_chimera_splitting"
            splitting = "off"
        else:
            options = f"-s {main_sequence}"
            splitting = "on"
        if len(to_trim) == 1:
            summary = f"{main_label} ({len(main_sequence)} bp)"
        else:
            summary = f"{len(to_trim)} sequences: " + "; ".join(label for sequence, label in to_trim)
    else:
        main_sequence = ""
        options = "-A"
        splitting = "off"
        summary = "nothing (no adapter found in the reads)"

    with open(args.out_scan, "w") as out:
        out.write("first_read\tlast_read\treads\treads_with_adapter\treads_with_barcode_block_only\n")
        for row in found["rows"]:
            out.write("\t".join(str(value) for value in row) + "\n")
    write_key_values(args.out_choice, [
        ("sample", args.sample),
        ("setting", "given sequence" if given else "detect"),
        ("reads", reads),
        ("reads_with_adapter", found["anchored"]),
        ("percent_with_adapter", f"{100 * anchored_share:.2f}"),
        ("reads_with_barcode_block_only", found["block_only"]),
        ("block", block or "none"),
        ("published_barcodes", barcodes_text),
        ("trimmed_sequence", main_sequence or "none"),
        ("trimmed_length", len(main_sequence)),
        ("trimmed_sequences", len(to_trim)),
        ("chimera_splitting", splitting),
        ("summary", summary),
        ("note", " ".join(notes) or "none"),
    ])
    with open(args.out_options, "w") as out:
        out.write(options + "\n")
    write_sequences(args.out_sequence, args.sample, to_trim)

    print(
        f"[FunFlux] find_ont_adapter {args.sample}: {found['anchored']:,} of {reads:,} reads "
        f"({100 * anchored_share:.1f} %) start with the ONT adapter; fastplong removes: {summary}."
    )
    for note in notes:
        print(f"[FunFlux] NOTE find_ont_adapter {args.sample}: {note}")


def check(args):
    """Subcommand check: how many trimmed reads still start with the adapter (or the block)?

    Takes in:  args.reads (the reads Filtlong kept after fastplong), args.choice (the scan's
               choice file), args.sample, args.out.
    Does:      counts the reads with the adapter core in their first START_WINDOW bp and, when
               the most common trimmed sequence was longer than the adapter, the reads with its
               last CHECK_TAIL_LENGTH bp there (barcode left behind).
    Produces:  args.out, "name<TAB>value" lines; a [FunFlux] NOTE when more than MIN_READ_SHARE
               of the reads still start with the adapter.
    Used by:   main().
    """
    choice = read_key_values(args.choice)
    construct = choice.get("trimmed_sequence", "none")
    tail = ""
    if construct not in ("none", "") and len(construct) > len(ADAPTER):
        tail = construct[-CHECK_TAIL_LENGTH:]

    reads = 0
    with_adapter = 0
    with_tail = 0
    for sequence in read_sequences(args.reads):
        reads += 1
        start = sequence[:START_WINDOW]
        if ADAPTER_CORE in start:
            with_adapter += 1
        if tail and tail in start:
            with_tail += 1

    adapter_share = with_adapter / reads if reads else 0.0
    tail_share = with_tail / reads if reads else 0.0
    write_key_values(args.out, [
        ("sample", args.sample),
        ("reads", reads),
        ("reads_with_adapter", with_adapter),
        ("percent_with_adapter", f"{100 * adapter_share:.2f}"),
        ("reads_with_end_of_trimmed_sequence", with_tail if tail else "NA"),
        ("percent_with_end_of_trimmed_sequence", f"{100 * tail_share:.2f}" if tail else "NA"),
    ])
    print(
        f"[FunFlux] check_ont_trimming {args.sample}: after fastplong and Filtlong, {with_adapter:,} of "
        f"{reads:,} reads ({100 * adapter_share:.2f} %) still start with the ONT adapter."
    )
    if adapter_share > MIN_READ_SHARE:
        print(
            f"[FunFlux] NOTE check_ont_trimming {args.sample}: more than {MIN_READ_SHARE:.0%} of the reads "
            "kept after trimming still start with the ONT adapter. The adapter may differ from the one FunFlux searched "
            "for; the contig adapter screen shows whether it reached the assembly."
        )


def build_parser():
    """Define the command line: the subcommands scan and check, and their options.

    The text at the top of this file says what each option holds.
    Used by:   main().
    """
    parser = argparse.ArgumentParser(description="Find the ONT adapter in reads, or check trimmed reads.")
    subcommands = parser.add_subparsers(dest="command", required=True)

    scan_parser = subcommands.add_parser("scan", help="decide what fastplong removes")
    scan_parser.add_argument("--reads", required=True)
    scan_parser.add_argument("--setting", required=True, help="detect, or the sequence to trim")
    scan_parser.add_argument("--sample", required=True)
    scan_parser.add_argument("--barcode-list", default="", help="ONT's published barcodes (FASTA)")
    scan_parser.add_argument("--out-scan", required=True)
    scan_parser.add_argument("--out-choice", required=True)
    scan_parser.add_argument("--out-options", required=True)
    scan_parser.add_argument("--out-sequence", required=True)

    check_parser = subcommands.add_parser("check", help="count the adapters left after fastplong")
    check_parser.add_argument("--reads", required=True)
    check_parser.add_argument("--choice", required=True)
    check_parser.add_argument("--sample", required=True)
    check_parser.add_argument("--out", required=True)
    return parser


def main(argv=None):
    """Run the subcommand named on the command line.

    Takes in:  argv, the command line without the script name (None means the real command
               line; the tests pass a list).
    Does:      scan() or check(). A file that cannot be read ends the run with exit status 1
               and one [FunFlux] line on stderr, instead of a Python traceback.
    Used by:   rules find_ont_adapter and check_ont_trimming (05_ont_adapters.smk).
    """
    args = build_parser().parse_args(argv)
    try:
        if args.command == "scan":
            scan(args)
        else:
            check(args)
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        sys.exit(f"[FunFlux] {args.command} {args.sample}: cannot read an input file: {error}")


if __name__ == "__main__":
    main()
