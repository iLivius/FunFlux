"""Tests for find_ont_adapter.py.

What the script does, in one sentence: it finds the Oxford Nanopore adapter, and the barcode
block after it, at the start of ONT reads, decides which sequence fastplong removes, and later
counts how many trimmed reads still carry the adapter.

The promises these tests pin down:
  * native-barcoded reads give the adapter + the 40 bp barcode block (68 bp), with fastplong's
    read splitting on, and the barcode is named from ONT's published list (NB24); reads without
    barcodes give the adapter alone, with splitting off;
  * a rapid barcode is named too, and its 90 bp block is found to the base;
  * a read whose barcode carries a few basecalling errors still counts for that barcode;
  * several published barcodes in one file are all trimmed (fastplong gets every sequence in both
    orientations through -a), with a NOTE when two or more are common; a rare second barcode is
    trimmed without a NOTE;
  * barcodes that ONT does not publish are trimmed as found in the reads, with a NOTE; several of
    them give the shared part only, with a NOTE;
  * reads that are already trimmed give -A (fastplong copies them unchanged);
  * the adapter is found even when only the last part of the file carries it (fastplong's own
    guess looks only at the first reads);
  * a telomere repeat is never chosen: not at read starts without adapter, and not as the block
    after the adapter;
  * a barcode block with no adapter in front is trimmed on its own, and a rapid one is named as
    rapid although its 24 bp barcode equals a native one;
  * a sequence given in the configuration is used as it is;
  * the check counts the reads that still start with the adapter, with a NOTE above 1 %.

Every test writes made-up reads to pytest's temporary folder (tmp_path). The genome parts of
the reads are random sequence from a seeded generator, so every run sees the same reads. The
barcode list is the one FunFlux ships, workflow/data/adapters/ont_chtd.fasta.

Run: pytest workflow/scripts/05_ont_adapters/test_find_ont_adapter.py
"""

import gzip
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import find_ont_adapter as foa

BARCODE_LIST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "adapters", "ont_chtd.fasta")
ADAPTER = foa.ADAPTER
# Native barcode NB24 with its flanks, as in the public B05.10 reads (ONT Chemistry Technical
# Document, section 14).
NB24_BLOCK = "AAGGTTAA" + "GCATAGTTCTGCATGATGGGTTAG" + "CAGCACCT"
NB01_BLOCK = "AAGGTTAA" + "CACAAAGACACCGACAACTTTCTT" + "CAGCACCT"
# Rapid barcode RB05 with its flanks: GCTTGGGTGTTTAACC, the barcode, the 50 bp sequence after it.
RB05_BLOCK = "GCTTGGGTGTTTAACC" + "CTTGTCCAGGGTTTGTGTAACCTT" + foa.RAPID_TAIL
# Made-up barcodes that ONT does not publish, between the native flanks.
UNKNOWN_BLOCK_1 = "AAGGTTAA" + "ACGTACGTTTGGCCAAGGTTCCAA" + "CAGCACCT"
UNKNOWN_BLOCK_2 = "AAGGTTAA" + "TTTTGGGGCCCCAAAATTGCATGC" + "CAGCACCT"
TELOMERE = "CCCTAA" * 12


def genome(generator, length):
    # Random sequence standing in for genomic DNA, from the test's seeded generator
    # (random.Random(n)), so that every run of the tests builds the same reads.
    return "".join(generator.choice("ACGT") for _ in range(length))


def write_fastq(path, sequences, compress=False):
    # One FASTQ record per sequence, with a constant quality line of the same length.
    opener = gzip.open if compress else open
    with opener(path, "wt") as out:
        for number, sequence in enumerate(sequences):
            out.write(f"@read{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
    return str(path)


def run_scan(tmp_path, reads_path, setting="detect", barcode_list=BARCODE_LIST):
    # Run the scan subcommand and return the choice file, the options line and the sequence FASTA.
    outputs = {name: str(tmp_path / name) for name in ("scan.tsv", "choice.tsv", "options.txt", "construct.fasta")}
    foa.main([
        "scan", "--reads", reads_path, "--setting", setting, "--sample", "S1",
        "--barcode-list", barcode_list,
        "--out-scan", outputs["scan.tsv"], "--out-choice", outputs["choice.tsv"],
        "--out-options", outputs["options.txt"], "--out-sequence", outputs["construct.fasta"],
    ])
    choice = foa.read_key_values(outputs["choice.tsv"])
    with open(outputs["options.txt"]) as handle:
        options = handle.read().strip()
    with open(outputs["construct.fasta"]) as handle:
        construct_fasta = handle.read()
    return choice, options, construct_fasta


# ── The decision ──

def test_native_barcoded_reads_give_adapter_plus_barcode_block(tmp_path):
    generator = random.Random(1)
    reads = [genome(generator, 7) + ADAPTER + NB24_BLOCK + genome(generator, 400) for _ in range(300)]
    choice, options, construct_fasta = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert choice["trimmed_sequence"] == ADAPTER + NB24_BLOCK
    assert choice["trimmed_length"] == "68"
    assert options == "-s " + ADAPTER + NB24_BLOCK
    assert choice["chimera_splitting"] == "on"
    assert choice["published_barcodes"] == "NB24 100.0 %"
    assert choice["summary"] == "adapter + native barcode NB24 (68 bp)"
    # Both orientations, because fastplong's -a does not add the reverse complement itself.
    assert ADAPTER + NB24_BLOCK in construct_fasta
    assert foa.reverse_complement(ADAPTER + NB24_BLOCK) in construct_fasta


def test_reads_without_barcodes_give_the_adapter_alone_with_splitting_off(tmp_path):
    generator = random.Random(2)
    reads = [genome(generator, 5) + ADAPTER + genome(generator, 400) for _ in range(300)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert choice["block"] == "none"
    assert options == "-s " + ADAPTER + " --disable_chimera_splitting"


def test_trimmed_reads_are_passed_through(tmp_path):
    generator = random.Random(3)
    reads = [genome(generator, 400) for _ in range(300)]
    choice, options, construct_fasta = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert options == "-A"
    assert choice["trimmed_sequence"] == "none"
    assert construct_fasta == ""


def test_adapter_only_in_the_last_part_of_the_file_is_found(tmp_path):
    # 80 % trimmed reads first, then untrimmed ones: fastplong's own guess would see only the first.
    generator = random.Random(4)
    reads = [genome(generator, 400) for _ in range(400)]
    reads += [genome(generator, 6) + ADAPTER + NB24_BLOCK + genome(generator, 400) for _ in range(100)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq.gz", reads, compress=True))
    assert options == "-s " + ADAPTER + NB24_BLOCK
    assert choice["percent_with_adapter"] == "20.00"


def test_telomere_at_read_starts_is_never_chosen(tmp_path):
    # Trimmed reads from chromosome ends start with the telomere repeat: no adapter, nothing to trim.
    generator = random.Random(5)
    reads = [TELOMERE + genome(generator, 400) for _ in range(300)]
    _, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert options == "-A"


def test_a_telomere_after_the_adapter_is_not_added_as_a_block(tmp_path):
    generator = random.Random(6)
    reads = [genome(generator, 5) + ADAPTER + TELOMERE + genome(generator, 300) for _ in range(300)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert choice["block"] == "none"
    assert options == "-s " + ADAPTER + " --disable_chimera_splitting"


def test_several_published_barcodes_are_all_trimmed_with_a_note(tmp_path, capsys):
    generator = random.Random(7)
    reads = [genome(generator, 5) + ADAPTER + NB24_BLOCK + genome(generator, 300) for _ in range(160)]
    reads += [genome(generator, 5) + ADAPTER + NB01_BLOCK + genome(generator, 300) for _ in range(140)]
    choice, options, construct_fasta = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    # The more common barcode goes to -s; both go to the FASTA that fastplong reads with -a.
    assert options == "-s " + ADAPTER + NB24_BLOCK
    assert choice["trimmed_sequences"] == "2"
    for block in (NB24_BLOCK, NB01_BLOCK):
        assert ADAPTER + block in construct_fasta
        assert foa.reverse_complement(ADAPTER + block) in construct_fasta
    assert choice["published_barcodes"] == "NB24 53.3 %, NB01 46.7 %"
    assert "2 barcodes" in capsys.readouterr().out


def test_a_rare_second_barcode_is_trimmed_without_a_note(tmp_path, capsys):
    # 3 % of the reads carry a neighbouring sample's barcode, as after demultiplexing.
    generator = random.Random(13)
    reads = [genome(generator, 5) + ADAPTER + NB24_BLOCK + genome(generator, 300) for _ in range(388)]
    reads += [genome(generator, 5) + ADAPTER + NB01_BLOCK + genome(generator, 300) for _ in range(12)]
    choice, options, construct_fasta = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert options == "-s " + ADAPTER + NB24_BLOCK
    assert ADAPTER + NB01_BLOCK in construct_fasta
    assert "NOTE" not in capsys.readouterr().out


def test_a_barcode_with_basecalling_errors_still_counts(tmp_path):
    # Two of the 24 barcode bases changed in every read: still NB24, at most 3 may differ.
    generator = random.Random(14)
    damaged = list(NB24_BLOCK)
    for position in (10, 20):   # both inside the barcode, which spans positions 8 to 31
        damaged[position] = "G" if damaged[position] == "T" else "T"
    damaged = "".join(damaged)
    assert foa.mismatches(damaged, NB24_BLOCK) == 2
    reads = [genome(generator, 5) + ADAPTER + damaged + genome(generator, 300) for _ in range(300)]
    choice, _, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert choice["published_barcodes"] == "NB24 100.0 %"
    # The block trimmed is the one in the reads, damaged bases and all.
    assert choice["trimmed_sequence"] == ADAPTER + damaged


def test_a_rapid_barcode_is_named_and_its_block_found_to_the_base(tmp_path):
    generator = random.Random(15)
    reads = [genome(generator, 6) + ADAPTER + RB05_BLOCK + genome(generator, 300) for _ in range(300)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    # 90 bp: the lengths are tried in steps of 4 (88) and then one base at a time.
    assert choice["trimmed_sequence"] == ADAPTER + RB05_BLOCK
    assert choice["summary"] == "adapter + rapid barcode RB05 (118 bp)"


def test_an_unpublished_barcode_is_trimmed_as_found_with_a_note(tmp_path, capsys):
    generator = random.Random(16)
    reads = [genome(generator, 5) + ADAPTER + UNKNOWN_BLOCK_1 + genome(generator, 300) for _ in range(300)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert options == "-s " + ADAPTER + UNKNOWN_BLOCK_1
    assert choice["published_barcodes"] == "none found"
    assert "not in ONT's published list" in capsys.readouterr().out


def test_several_unpublished_barcodes_give_the_shared_part_and_a_note(tmp_path, capsys):
    generator = random.Random(17)
    reads = [genome(generator, 5) + ADAPTER + UNKNOWN_BLOCK_1 + genome(generator, 300) for _ in range(150)]
    reads += [genome(generator, 5) + ADAPTER + UNKNOWN_BLOCK_2 + genome(generator, 300) for _ in range(150)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert choice["trimmed_sequence"] == ADAPTER + "AAGGTTAA"
    assert "--disable_chimera_splitting" in options
    assert "only the part they share is trimmed" in capsys.readouterr().out


def test_without_the_list_the_block_is_still_trimmed_and_nothing_is_named(tmp_path, capsys):
    generator = random.Random(18)
    reads = [genome(generator, 7) + ADAPTER + NB24_BLOCK + genome(generator, 400) for _ in range(300)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads), barcode_list="")
    assert options == "-s " + ADAPTER + NB24_BLOCK
    assert choice["published_barcodes"] == "none found"
    assert "NOTE" not in capsys.readouterr().out


def test_a_rapid_block_without_adapter_is_named_as_rapid(tmp_path):
    # RB20's 24 bp barcode equals NB20's: the flanks must decide, and the summary say "rapid".
    generator = random.Random(19)
    barcodes = {b["name"]: b for b in foa.read_barcode_list(BARCODE_LIST)}
    assert barcodes["RB20"]["barcode"] == barcodes["NB20"]["barcode"]
    reads = [genome(generator, 4) + barcodes["RB20"]["block"] + genome(generator, 400) for _ in range(300)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert choice["summary"] == "rapid barcode RB20, with no adapter in front (90 bp)"
    assert choice["published_barcodes"] == "RB20 100.0 % of all reads, no adapter"


def test_the_shipped_list_holds_the_published_barcodes():
    barcodes = foa.read_barcode_list(BARCODE_LIST)
    names = {barcode["name"]: barcode for barcode in barcodes}
    assert len(barcodes) == 192
    assert names["NB24"]["block"] == NB24_BLOCK
    assert names["RB05"]["block"] == RB05_BLOCK


def test_a_barcode_block_without_adapter_is_trimmed_on_its_own(tmp_path):
    generator = random.Random(8)
    reads = [genome(generator, 4) + NB24_BLOCK + genome(generator, 400) for _ in range(300)]
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    assert choice["trimmed_sequence"] == NB24_BLOCK
    assert options == "-s " + NB24_BLOCK + " --disable_chimera_splitting"


def test_a_given_sequence_is_used_as_it_is(tmp_path):
    generator = random.Random(9)
    reads = [genome(generator, 400) for _ in range(50)]
    given = ADAPTER + NB24_BLOCK
    choice, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads), setting=given.lower())
    assert options == "-s " + given
    assert choice["setting"] == "given sequence"


def test_an_empty_file_gives_pass_through(tmp_path):
    _, options, _ = run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", []))
    assert options == "-A"


def test_the_scan_table_counts_reads_in_rows_of_ten_thousand(tmp_path):
    generator = random.Random(10)
    reads = [ADAPTER + genome(generator, 60) for _ in range(10_005)]
    run_scan(tmp_path, write_fastq(tmp_path / "reads.fastq", reads))
    with open(tmp_path / "scan.tsv") as handle:
        rows = [line.split("\t") for line in handle.read().strip().split("\n")]
    assert rows[1][:4] == ["1", "10000", "10000", "10000"]
    assert rows[2][:3] == ["10001", "10005", "5"]


# ── Small helpers ──

def test_tandem_repeats_are_recognised():
    assert foa.is_tandem_repeat("TTAGGGTTAGGGTTAGGG")
    assert foa.is_tandem_repeat("TTAGGGTTAGGGTT")
    assert not foa.is_tandem_repeat(ADAPTER)
    assert not foa.is_tandem_repeat(NB24_BLOCK)
    assert not foa.is_tandem_repeat("AAGGTTAA")
    # The unit must fit at least twice: a 9 to 12 bp block is not a repeat of itself,
    # while a short unit repeated over 9 bp still is.
    assert not foa.is_tandem_repeat("AAGGTTAAC")
    assert not foa.is_tandem_repeat("CCTGTACTTCCC")
    assert foa.is_tandem_repeat("ATATATATA")


def test_clean_sequence_upper_cases_and_refuses_other_letters():
    assert foa.clean_sequence(" cctg tac ") == "CCTGTAC"
    assert foa.clean_sequence("CCTGNAC") is None
    assert foa.clean_sequence("") is None


# ── The check ──

def test_the_check_counts_reads_that_still_carry_the_adapter(tmp_path, capsys):
    generator = random.Random(11)
    choice_path = tmp_path / "choice.tsv"
    foa.write_key_values(choice_path, [("trimmed_sequence", ADAPTER + NB24_BLOCK)])
    reads = [genome(generator, 300) for _ in range(95)]
    reads += [ADAPTER + NB24_BLOCK + genome(generator, 300) for _ in range(5)]
    reads_path = write_fastq(tmp_path / "trimmed.fastq", reads)
    out = tmp_path / "check.tsv"
    foa.main(["check", "--reads", reads_path, "--choice", str(choice_path), "--sample", "S1", "--out", str(out)])
    result = foa.read_key_values(out)
    assert result["reads"] == "100"
    assert result["reads_with_adapter"] == "5"
    assert result["reads_with_end_of_trimmed_sequence"] == "5"
    assert "NOTE" in capsys.readouterr().out


def test_the_check_is_quiet_when_the_reads_are_clean(tmp_path, capsys):
    generator = random.Random(12)
    choice_path = tmp_path / "choice.tsv"
    foa.write_key_values(choice_path, [("trimmed_sequence", "none")])
    reads_path = write_fastq(tmp_path / "trimmed.fastq", [genome(generator, 300) for _ in range(100)])
    out = tmp_path / "check.tsv"
    foa.main(["check", "--reads", reads_path, "--choice", str(choice_path), "--sample", "S1", "--out", str(out)])
    result = foa.read_key_values(out)
    assert result["reads_with_adapter"] == "0"
    assert result["reads_with_end_of_trimmed_sequence"] == "NA"
    assert "NOTE" not in capsys.readouterr().out
