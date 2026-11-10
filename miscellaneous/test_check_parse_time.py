"""job_resources() in check_parse_time.py must read the resources of a dry run in any Snakemake version.

What this test protects
-----------------------
check_parse_time.py checks the memory and time limits FunFlux's jobs ask a cluster for
(resources.ram_gb and resources.runtime_hours) in the resources line Snakemake prints
for every job of a dry run. Snakemake versions print that line differently: 9.14.6 writes
tmpdir=<TBD> and lists mem_mb before runtime, 9.27.1 writes the real temporary folder,
lists runtime first and adds an entry of its own (mem=32.77 GB). The automatic checks on
GitHub run the dry runs with 9.14.6 only, so this test feeds job_resources() the same
two jobs written both ways and requires the same values for the entries FunFlux sets.
It needs no Snakemake.

Run with: pytest miscellaneous
"""

import os
import sys

# check_parse_time.py sits next to this file; put this folder on Python's search path
# so that it can be imported by name.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from check_parse_time import job_resources


# Two jobs of the illumina example's dry run, shortened, as Snakemake 9.14.6 prints them:
# map_phix asks only for the time, the SPAdes job for memory and time.
OUTPUT_9_14_6 = """
rule map_phix:
    input: SAM1_R1.fastq.gz, SAM1_R2.fastq.gz
    jobid: 12
    wildcards: sample=SAM1
    resources: tmpdir=<TBD>, runtime=2880

Shell command: bowtie2 -x phix ...

rule illumina_assembly:
    jobid: 14
    wildcards: sample=SAM1
    resources: tmpdir=<TBD>, mem_gb=32, mem_mb=32768, mem_mib=31250, runtime=2880

Shell command: OMP_NUM_THREADS=1 spades.py ...
"""

# The same two jobs as Snakemake 9.27.1 prints them.
OUTPUT_9_27_1 = """
rule map_phix:
    input: SAM1_R1.fastq.gz, SAM1_R2.fastq.gz
    jobid: 12
    wildcards: sample=SAM1
    resources: tmpdir=/tmp, runtime=2880

Shell command: bowtie2 -x phix ...

rule illumina_assembly:
    jobid: 14
    wildcards: sample=SAM1
    resources: tmpdir=/tmp, mem_gb=32, runtime=2880, mem_mb=32768, mem=32.77 GB, mem_mib=31250

Shell command: OMP_NUM_THREADS=1 spades.py ...
"""


def test_both_versions_give_the_same_values():
    # The values FunFlux sets must come out the same, whatever the order and the extra
    # entries of each version.
    for output in (OUTPUT_9_14_6, OUTPUT_9_27_1):
        jobs = job_resources(output)
        spades = jobs["illumina_assembly"][0]
        assert spades["mem_gb"] == "32"
        assert spades["mem_mb"] == "32768"
        assert spades["runtime"] == "2880"
        phix = jobs["map_phix"][0]
        assert phix["runtime"] == "2880"
        # map_phix sets no memory: check_parse_time's cases read a missing entry as
        # "not set".
        assert "mem_mb" not in phix


def test_one_dict_per_job_and_localrule_blocks():
    # Two jobs of one rule (two samples) give two dicts, in the order printed; a job
    # that runs on the launch machine starts with "localrule" and is read the same way.
    output = (
        "rule trim_adapters:\n    wildcards: sample=SAM1\n    resources: tmpdir=<TBD>\n\n"
        "rule trim_adapters:\n    wildcards: sample=SAM2\n    resources: tmpdir=<TBD>\n\n"
        "localrule download_phix:\n    resources: tmpdir=/tmp\n"
    )
    jobs = job_resources(output)
    assert len(jobs["trim_adapters"]) == 2
    assert jobs["trim_adapters"][1] == {"tmpdir": "<TBD>"}
    assert jobs["download_phix"] == [{"tmpdir": "/tmp"}]


def test_no_jobs_gives_an_empty_result():
    # A refused configuration stops before any job is planned.
    assert job_resources("[FunFlux] mode: write illumina, longread, hybrid or contigs") == {}
