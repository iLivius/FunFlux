# Output files

This page lists the files and folders a FunFlux run writes, and whether each one is still there
when the run has finished. Every path on this page is relative to `directories.output_dir`, and a
relative value of that key is read from the launch folder, the folder Snakemake is started from.
The numbered folders are the same in every mode, and a mode with nothing to put into a folder
does not create it.

Most tables name the rule that writes each file, as the console and `logs/` name it. They also
give the file's status, which says whether it stays, and the modes or settings under which it is
written. The terms below define the labels of the Status and Modes columns.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **deliverable** | a result to read, keep or pass on; a later step may read it too |
    | **intermediate** | a working file that one step hands to the next; kept after the run |
    | **temporary** | deleted by Snakemake once the last step that reads it has finished; absent from a finished run |
    | **scratch** | deleted by the step that wrote it, when that step ends; a download whose checksum did not match is kept instead, for you to inspect |
    | **cleanup** | deleted by `cleanup`, the last step of each sample |
    | **screen** | the contamination and mitochondrion screen, always run in `illumina`, `longread` and `hybrid` mode, and in `contigs` mode only with `parameters.contigs.screen: true` |
    | **SignalP, GeneMark, Tiberius** | for SignalP, `parameters.signalp6.mode` is not `off`; for GeneMark, GeneMark-ES is enabled (`directories.genemark_dir` set, `parameters.genemark.mode: auto`); for Tiberius, `parameters.tiberius.mode` is `on`, the default |
    | **AMR** | antimicrobial resistance; as a label, the AMR screen is switched on: `parameters.amr.chroquetas` is `on` (off by default) |

## Files to read first

These are the files that hold the results of a run, with the question each one answers and the
modes that write it.

| File | What it answers | Modes |
|---|---|---|
| `02.assembly/{sample}/contigs_final.fasta` | The nuclear genome. QUAST, BUSCO, the ITS search and funannotate read this file and no other version of the assembly. | all |
| `02.assembly/{sample}/contigs_mitochondrial.fasta` | The confirmed mitochondrial contig(s); an empty file when none was confirmed. | screen |
| `02.assembly/{sample}/mitochondria/mitochondrial_screen.tsv` | Why each contig was or was not called mitochondrial. | screen |
| `02.assembly/{sample}/contaminants/contig_taxonomy_decisions.tsv` | Genus, keep or remove, and the reason, for every contig. | screen |
| `02.assembly/{sample}/eval/busco/short_summary.specific.<lineage>.{sample}.txt` | Completeness of the genome. | all |
| `02.assembly/{sample}/eval/quast/report.tsv` | Contig count, total length, N50, GC content. | all |
| `02.assembly/{sample}/eval/contig_summary.tsv` | What each contig of the draft assembly is (chromosome, chromosome arm, mitochondrion, rDNA array, low-depth piece, collapsed repeat, nuclear contig or fragment, contaminant, another organism's mitochondrial DNA, unclassified), with the evidence; it also links each contig to its `scaffold_N` name in the annotation. | all |
| `03.taxonomy/{sample}/{sample}_ITS1_taxonomy.txt`, `{sample}_ITS2_taxonomy.txt` | UNITE classification of the ITS barcode. | all |
| `03.taxonomy/{sample}/{sample}_ITS_copies.txt` | How many rDNA units were found, and which records the length guard rejected. | all |
| `04.annotation/{sample}/funannotate/annotate_results/` | The annotated genome: GenBank, GFF3, proteins, transcripts, and `.annotations.txt`, the functional annotation with one row per transcript. | all |
| `04.annotation/{sample}/contigs_sort.fasta` | The genome as annotated. | all |
| `05.amr/{sample}/amr_summary.tsv` | The FungAMR mutations found in the genome, one row per target protein of the species, with a verdict on resistance and the drugs with evidence for and against it. | AMR |
| `06.report/multiqc_report.html` | QUAST and BUSCO of every sample of the run on one page; in the read-based modes also the read QC (fastp, NanoPlot) and the mapping QC (Qualimap). | all |

!!! note "The annotation and the assembly use different contig names"

    `funannotate sort` renames the contigs `scaffold_1`, `scaffold_2`, ..., longest first, before
    repeat masking, so the masked genome and every annotation file use these names. The screen,
    QUAST, BUSCO and the ITS files use the names of `contigs_final.fasta`. The `annotated_as`
    column of `02.assembly/{sample}/eval/contig_summary.tsv` links the two, as described on
    [Contig summary](../analysis/contig-summary.md), and
    `04.annotation/{sample}/contigs_clean.fasta` holds the same sequences under their original
    names.

## Layout

The tree shows the folders of an output folder and the main files in them. A label in square
brackets names the modes that write a line, or the setting it needs. On the `eval/` line, the
label covers only `qualimap/`. A line without a label is written in every mode, although what
`00.resources/` holds depends on the mode and the settings, as its table shows.

```text
<output_dir>/
├── 00.resources/                        run-level material: databases, Tiberius model, SignalP 6 installation, ChroQueTas, a named BUSCO lineage
├── 01.reads/{sample}/                   [illumina, longread, hybrid]
├── 02.assembly/{sample}/
│   ├── contigs_final.fasta              nuclear genome
│   ├── contigs_all.fasta                screened genome, mitochondrion included   [screen]
│   ├── contigs_mitochondrial.fasta                                                [screen]
│   ├── contigs_filt.fasta               the draft the screen examines   [contigs + screen, illumina]
│   ├── contaminants/                    BLAST, BlobTools, keep/remove decisions   [screen]
│   ├── mitochondria/                    mitochondrion screen                      [screen]
│   ├── spades/                                                                    [illumina]
│   ├── hifiasm/ or flye/, assembly.fasta, assembly_evidence.tsv, assembly_provenance.tsv   [longread, hybrid]
│   ├── {sample}_medaka_model.txt, medaka/            [longread, hybrid: ONT samples with Medaka]
│   ├── polypolish/, guide/                                                        [hybrid]
│   └── eval/                            quast/, busco/, telomere_*, adapter_*, contig_summary*.tsv; qualimap/ [illumina, longread, hybrid]
├── 03.taxonomy/{sample}/                ITS windows, ITS sequences, UNITE classification
├── 04.annotation/{sample}/              masking, gene prediction, functional annotation
├── 05.amr/{sample}/                     antifungal resistance mutations (ChroQueTas)   [AMR]
├── 06.report/                           MultiQC
└── logs/                                one log per rule and sample
```

### Disk use

One sample in `contigs` mode takes in the order of 10 GB. The table shows where the space goes.

| Folder | What takes the space | Size |
|---|---|---|
| `01.reads/{sample}/`, `02.assembly/{sample}/` (read-based modes only) | the kept reads in `01.reads/` (the filtered long reads, uncompressed; in `hybrid` also the selected pairs of Illumina reads), the assembler's working folder (`spades/`, `hifiasm/` or `flye/`), the read-mapping BAM (two in `hybrid`) and the polishing output | 3.7 to 30 GB for the two folders together in the runs [below](#disk-use-per-run) |
| `04.annotation/{sample}/` | most of the rest: the InterProScan XML, there twice (`iprscan/interproscan_out.xml` and funannotate's copy `funannotate/annotate_misc/iprscan.xml`), and funannotate's working files in `funannotate/predict_misc/` and `annotate_misc/` | the XML can be a few GB |
| `02.assembly/{sample}/` | the BUSCO result tree | 0.65 to 1.8 GB with the automatic lineage in the runs [below](#disk-use-per-run) |
| `00.resources/` | the antiSMASH 8.0.4 database, unless `directories.antismash_db` is set | 9.4 GB |
| `00.resources/` | the SignalP 6 package, with the fast model | 1.6 GB |
| `00.resources/` | UNITE: the reference and, by default, its copy of the records that name a genus | 354 MB (207 + 147 MB) |
| `00.resources/` | the Tiberius model | 31 MB |
| `00.resources/` | the BUSCO lineage data, when `parameters.busco.lineage` names a lineage | a few hundred MB |
| `00.resources/` | ChroQueTas with its FungAMR database (AMR screen on) | 3 MB |
| `00.resources/` | the telomere repeat database | 0.6 MB |

In the long-read modes, the kept long reads are the largest single item. Apart from them, the
annotation and the assembly folder take most of a sample's space. The first keeps the
InterProScan XML twice, the second the read-mapping BAM and the assembler's working files. Of
the material shared by the run, the antiSMASH database is the largest item. With
`directories.antismash_db` set, `00.resources/` holds only links to your copy.

To free space once you have checked the results, run
`workflow/scripts/clean_funflux_output.sh <output folder>`. The script lists the bulky working
files that no result needs, and deletes them only with `--run`:

- assembler caches and read mappings;
- the working files and copies of Medaka, BUSCO and funannotate;
- the antiSMASH and UNITE databases and the BUSCO lineage data FunFlux downloaded;
- the licensed SignalP 6 package.

The reads kept in `01.reads/` go only with `--include-reads`. The results, the logs and tbl2asn's
validation reports stay.

!!! warning "Cleaning is for a folder that FunFlux will not use again"

    The script refuses a folder in which a sample has not finished, but it cannot tell whether
    another FunFlux run is using the folder. Do not give a cleaned folder to FunFlux again: the
    deleted files would be made again and most of each sample repeated.

The conda environments are not counted here. They are in `.snakemake/` of the launch folder, or
in the `--conda-prefix` folder, as listed under
[Outside the output folder](#outside-the-output-folder).

#### Disk use per run { #disk-use-per-run }

The sizes below were measured on two public genomes, with one sample per run, and are given as
`du -sh` prints them. The runs are described on
[Validation on public genomes](../about/validation.md). The four B05.10 runs with ONT reads were
made with FunFlux 2.1.0, the others with FunFlux 2.0.0.

| Genome | Run | `01.reads/` | `02.assembly/` | Whole output folder |
|---|---|---|---|---|
| *Botrytis cinerea* B05.10, 42.6 Mb | `illumina`, Illumina 79x | 0.6 MB | 3.7 GB | 13 GB |
| | `longread`, ONT 42x, hifiasm, Medaka | 3.1 GB | 7.1 GB | 20 GB |
| | `longread`, ONT 136x, hifiasm | 9.8 GB | 12 GB | 31 GB |
| | `longread`, ONT 136x, Flye | 9.8 GB | 6.6 GB | 26 GB |
| | `hybrid`, ONT 136x plus Illumina 79x, hifiasm | 16 GB | 14 GB | 38 GB |
| *Cordyceps cateniannulata* GXU-8616, 32.7 Mb | `longread`, PacBio HiFi 85x, hifiasm | 4.7 GB | 3.1 GB | 17 GB |
| | `longread`, PacBio HiFi 85x, Flye | 4.7 GB | 2.3 GB | 16 GB |
| | `hybrid`, PacBio HiFi 85x plus MGISEQ short reads 151x, hifiasm | 14 GB | 5.7 GB | 29 GB |

A whole sample took 13 to 38 GB in these runs, and the kept long reads were the largest single
file of every run with long reads. Their size follows the bases sequenced: 9.7 GB from the
5.8 Gb of ONT reads of B05.10, 4.7 GB from the 2.8 Gb of HiFi reads of GXU-8616.

- In `hybrid`, `01.reads/` also holds the selected pairs, uncompressed: 9.3 GB of the 14 GB for
  GXU-8616.
- `02.assembly/` was smaller with Flye than with hifiasm on both genomes, because hifiasm keeps
  larger working files: `hifiasm/` took 5.5 GB on B05.10 and 1.5 GB on GXU-8616, `flye/` 860 MB
  and 570 MB.
- The whole output folder of every run includes about 2 GB of `00.resources/`, most of it the
  SignalP 6 package. It does not include the antiSMASH database, because these runs set
  `directories.antismash_db`.

## 00.resources

`00.resources/` holds the material that belongs to the run, not to a sample. It is kept, so a
later run on the same output folder downloads or installs none of it again unless the settings
behind it change.

| Path | Content | Written by | Status | Modes |
|---|---|---|---|---|
| `unite/unite_its_sintax.fasta` | UNITE reference for SINTAX, unpacked from `links.unite_its_link` (207 MB, 266,589 sequences with the default link); with `directories.unite_db`, a symbolic link to your copy | `unite_db` or `unite_db_local` | intermediate | all |
| `unite/PROVENANCE.txt` | source URL (or source folder), UTC date, md5 and size in bytes of that FASTA | same | deliverable | all |
| `unite/unite_its_sintax.identified.fasta` | the records of that FASTA whose taxonomy names a genus (147 MB, 180,343 sequences with the default link; 86,246 are dropped); what SINTAX classifies against while `parameters.its.unite_identified_only` is `true` | `unite_identified` | intermediate | all |
| `unite/identified_filter.txt` | two lines: `kept_with_genus: <n>` and `dropped_without_genus: <n>`, the records the filter kept and dropped (180,343 and 86,246 with the default link) | `unite_identified` | deliverable | all |
| `antismash/databases/` | antiSMASH database from `download-antismash-databases`; with `directories.antismash_db`, one symbolic link per entry of your copy | `antismash_db` or `antismash_db_local` | intermediate | all |
| `funannotate_db_paths.normalized.json` | marker that the paths inside the funannotate database were repaired, with the number of files checked and paths rewritten | `normalize_funannotate_db` | intermediate | all |
| `signalp6/installation.json` | archive name and SHA-256, SignalP version, mode, model folder, and the conda environment SignalP 6 was installed into | `install_signalp6` | intermediate | SignalP |
| `signalp6/package/` | the SignalP 6 package with the model weights of the DTU archive (1.6 GB, fast model) | `install_signalp6` | intermediate | SignalP |
| `signalp6/.matplotlib/` | matplotlib cache | `install_signalp6` | intermediate | SignalP |
| `signalp6/.signalp6-install-*/` | the unpacked archive during installation | `install_signalp6` | scratch | SignalP |
| `tiberius/model/` | the Tiberius `fungi` model, unpacked from `links.tiberius_model_link` after its SHA-256 matched `parameters.tiberius.model_sha256`: `model_config.json`, `model_layers.json`, `weights.h5` (31 MB) | `tiberius_model` | intermediate | Tiberius |
| `tiberius/model.tar.gz` | the downloaded archive; deleted once unpacked, kept when its checksum does not match | `tiberius_model` | scratch | Tiberius |
| `telomeres/telo_database_v2.tsv` | the telomere repeat database of tidk, downloaded from `links.telomere_db_link` after its SHA-256 matched `parameters.contig_summary.telomere_db_sha256`: 4,051 rows, one per examined species, with the telomere repeat unit and its support (0.6 MB) | `telomere_db` | intermediate | all |
| `telomeres/PROVENANCE.txt` | source URL, download date, SHA-256 and row count of that file | `telomere_db` | deliverable | all |
| `telomeres/telo_database_v2.tsv.failed` | the downloaded file when its checksum did not match; absent otherwise | `telomere_db` | scratch | all |
| `adapters/univec_illumina.fasta` | the 151 Illumina adapter and primer records of NCBI UniVec build 10.1, downloaded from `links.univec_link` after its SHA-256 matched `parameters.contig_summary.univec_sha256`; the contig adapter screen searches for them | `univec_illumina` | intermediate | all |
| `adapters/PROVENANCE.txt` | source URL, download date, SHA-256 and the number of Illumina records kept | `univec_illumina` | deliverable | all |
| `adapters/univec_illumina.fasta.download.failed` | the downloaded file when its checksum did not match; absent otherwise | `univec_illumina` | scratch | all |
| `chroquetas/ChroQueTas.sh` | the ChroQueTas script, version 1.0.1, downloaded unchanged from `links.chroquetas_link` after its SHA-256 matched `parameters.amr.chroquetas_script_sha256` (29 kB) | `chroquetas_db` | intermediate | AMR |
| `chroquetas/FungAMR_db/` | the FungAMR database of ChroQueTas: one folder per species (57), each with a `<target>.faa` (the reference protein) and a `<target>.txt` (its resistance positions) per target; unpacked from `db-v20250811.tgz` after its SHA-256 matched `parameters.amr.chroquetas_db_sha256` | `chroquetas_db` | intermediate | AMR |
| `chroquetas/db-v20250811.tgz` | the downloaded database archive; deleted once unpacked | `chroquetas_db` | scratch | AMR |
| `chroquetas/species.txt` | the database's species folders, one per line, sorted | `chroquetas_db` | intermediate | AMR |
| `chroquetas/LICENSE` | the licence of ChroQueTas and its database, CC BY-NC-ND 4.0 | `chroquetas_db` | deliverable | AMR |
| `chroquetas/PROVENANCE.txt` | the three source addresses, the UTC date, the checksums, the species count and the licence | `chroquetas_db` | deliverable | AMR |
| `chroquetas/ChroQueTas.sh.failed`, `chroquetas/db-v20250811.tgz.failed` | a download whose checksum did not match; absent otherwise | `chroquetas_db` | scratch | AMR |
| `phix/` | `phix.fna.gz`, the PhiX genome from `links.phix_link`, still compressed, and its six bowtie2 index files `phix.*.bt2` | `download_phix`, `build_phix` | intermediate | illumina, hybrid |
| `busco/lineages/<lineage>/` | the BUSCO dataset named in `parameters.busco.lineage`, downloaded once and read by every BUSCO job with `--offline` (a few hundred MB); later runs into the folder use the same release | `busco_data` | intermediate | all, with a named lineage |
| `busco/file_versions.tsv` | BUSCO's list of the dataset versions it downloaded | `busco_data` | deliverable | all, with a named lineage |

The folders named in `directories.unite_db` and `directories.antismash_db` are only read.
`signalp6/package/` is licensed DTU material, so remove it from any copy of the output folder you
share, as explained on [Licensing](../about/licensing.md). `chroquetas/` is CC BY-NC-ND 4.0
material: it may be passed on unchanged, for non-commercial purposes and with credit.

## 01.reads

`01.reads/` holds the processed reads and their QC reports, in one folder per sample and read
type. It is not created in `contigs` mode.

| Path | Content | Written by | Status | Modes |
|---|---|---|---|---|
| `{sample}/illumina/{sample}_nophix_R1.fastq`, `_R2.fastq` | the read pairs that did not align to the PhiX genome, uncompressed | `map_phix` | temporary | illumina, hybrid |
| `{sample}/illumina/{sample}_trim_R1.fastq`, `_R2.fastq` | the PhiX-free pairs after fastp: adapters and low-quality ends cut, reads under 75 bp dropped | `trim_adapters` | temporary | illumina, hybrid |
| `{sample}/illumina/{sample}_fastp.html`, `{sample}_fastp.json` | fastp's read-QC report as a page and as numbers; MultiQC reads the JSON | `trim_adapters` | deliverable | illumina, hybrid |
| `{sample}/long/adapters/{sample}_adapter_scan.tsv` | ONT samples: the reads that start with the ONT adapter, and those with a barcode block and no adapter, in every block of 10,000 reads | `find_ont_adapter` | deliverable | longread, hybrid |
| `{sample}/long/adapters/{sample}_adapter_choice.tsv` | ONT samples: what fastplong removes and why, with the numbers behind the decision and ONT's published barcodes found in the reads | `find_ont_adapter` | deliverable | longread, hybrid |
| `{sample}/long/adapters/{sample}_fastplong_options.txt`, `{sample}_adapter_sequence.fasta` | ONT samples: the options handed to fastplong, and every trimmed sequence in both orientations as FASTA (empty when nothing was trimmed); fastplong reads the FASTA with `-a`, and the contig adapter screen searches for those sequences too | `find_ont_adapter` | intermediate | longread, hybrid |
| `{sample}/long/{sample}_trimmed.fastq.gz` | ONT samples: the reads after fastplong; Filtlong reads them | `trim_ont_adapters` | temporary | longread, hybrid |
| `{sample}/long/adapters/{sample}_fastplong.html`, `{sample}_fastplong.json` | ONT samples: fastplong's report, as a page and as numbers (reads trimmed, reads split) | `trim_ont_adapters` | deliverable | longread, hybrid |
| `{sample}/long/adapters/{sample}_adapter_check.tsv` | ONT samples: the trimmed reads that still start with the adapter | `check_ont_trimming` | deliverable | longread, hybrid |
| `{sample}/long/{sample}_filtered.fastq` | the long reads Filtlong kept, uncompressed; the assembler, the coverage mapping and Medaka read it | `filter_long_reads` | intermediate | longread, hybrid |
| `{sample}/long/raw_qc/`, `filtered_qc/` | NanoPlot on the raw and on the filtered reads: `{sample}_NanoStats.txt` (MultiQC reads it), an HTML report, the plots as HTML pages (no PNG, `--no_static`) | `raw_long_read_qc`, `filtered_long_read_qc` | deliverable | longread, hybrid |
| `{sample}/selected/{sample}_selected_R1.fastq`, `_R2.fastq` | the trimmed Illumina pairs properly paired on the screened guide assembly; Polypolish reads them, and so does Filtlong for ONT samples with `short_read_guidance: true` | `map_sel_contigs` | intermediate | hybrid |

## 02.assembly

`02.assembly/` holds one folder per sample with the genome files, the draft and the assembler's
files, the results of the screen and the assembly QC. In `contigs` mode without the screen, the
folder holds only `contigs_final.fasta` and `eval/`.

| Path | Content | Written by | Status | Modes |
|---|---|---|---|---|
| `contigs_final.fasta` | the nuclear genome | `filter_contigs` (no screen) or `separate_mitochondria` | deliverable | all |
| `contigs_all.fasta` | the screened genome, mitochondrion included | `finalize_assembly` | deliverable | screen |
| `contigs_mitochondrial.fasta` | the confirmed mitochondrial contigs, headers unchanged; an empty file when none was confirmed or `parameters.mitochondria.mode` is `off` or `report`. A contig whose best mitochondrial BLAST hit is an animal's or a plant's, such as a contig of human mitochondrial DNA, is classed `foreign_mitochondrial` and is not confirmed: it stays in the assembly, where the contamination screen keeps or removes it like any other contig | `separate_mitochondria` | deliverable | screen |
| `contigs_filt.fasta` | the draft the screen examines: your FASTA after the contig filter (`contigs`), or the SPAdes contigs of at least 500 bp and 2x coverage (`illumina`) | `filter_contigs` | intermediate | contigs + screen, illumina |
| `spades/` | SPAdes' working folder: the unfiltered `contigs.fasta`, `scaffolds.fasta`, the assembly graph, one `K<k>/` folder per k-mer size, `spades.log`; emptied whenever SPAdes runs again | `illumina_assembly` | intermediate | illumina |
| `hifiasm/` or `flye/` | the assembler's own files: hifiasm's graphs (`hifiasm.p_ctg.gfa`, the primary contigs; `hifiasm.a_ctg.gfa`) and `.bin` files, or Flye's `assembly.fasta`, `assembly_info.txt`, graph and working files; emptied whenever the assembler runs again | `long_read_assembly` | intermediate | longread, hybrid |
| `assembly.fasta` | the long-read draft the screen examines, one line per contig (hifiasm: the primary contigs; Flye: a copy of `flye/assembly.fasta`) | `long_read_assembly` | intermediate | longread, hybrid |
| `assembly_evidence.tsv` | per contig: `contig`, `source`, `circular`, `coverage`, `repeat`, `alt_group`, taken from `hifiasm/hifiasm.p_ctg.gfa` or `flye/assembly_info.txt`; the mitochondrion screen and the contig summary read it | `assembly_evidence` | deliverable | longread, hybrid |
| `assembly_provenance.tsv` | `field` and `value` columns: assembler, version, read type, input file, command line | `long_read_assembly` | deliverable | longread, hybrid |
| `{sample}_medaka_model.txt` | one line: the Medaka model name after the check, or the path of the model file `auto` found; the assembler waits for it, Medaka reads it | `check_medaka_model` | intermediate | longread, hybrid: ONT samples with Medaka |
| `medaka/consensus.fasta` | the screened draft polished by Medaka; `medaka/` also holds Medaka's working files and is emptied whenever Medaka runs again | `long_read_consensus` | intermediate | longread, hybrid: ONT samples with Medaka |
| `polypolish/polished.fasta` | the draft after correction with the selected pairs, headers cut at the first blank, one line per sequence | `short_read_correction` | intermediate | hybrid |
| `guide/` | the Illumina guide assembly and its screen ([below](#guide)) | | | hybrid |

What the contig filter, rule `filter_contigs`, does depends on the mode and, in `contigs` mode,
on the first FASTA header:

| Mode and draft | What `filter_contigs` does |
|---|---|
| `contigs`, SPAdes-style headers (`>NODE_12_length_45678_cov_23.4`) | the SPAdes filter: keeps contigs of at least 500 bp and 2x coverage |
| `contigs`, any other header | cuts each header at its first space or tab; removes nothing |
| `illumina`, SPAdes' `contigs.fasta` | the same SPAdes filter; an empty result stops the run |

An assembly from another assembler therefore loses no contig at this step, only the text after
the first space or tab of its headers. See [Pre-assembled contigs](../modes/contigs.md).

### contaminants/

`contaminants/` holds the contamination part of the screen: the BLAST hits, the BlobTools table,
the decisions of the contig selector (`select_contigs`) and the coverage track. In the read-based
modes it also holds the read depth of every contig, measured from the coverage track, for the
contig summary.

| Path | Content | Written by | Status |
|---|---|---|---|
| `blastout` | megablast hits of every contig, 15 tab-separated columns | `blast_contigs` | intermediate |
| `bestscore.blob.blobDB.table.txt` | BlobTools table: per contig, length, GC, coverage and the best-sum taxon at every rank | `blob_table` | deliverable |
| `contig_taxonomy_decisions.tsv` | `contig`, `assigned_genus`, `action` (`keep`, `remove`), `reason` | `select_contigs` | deliverable |
| `{sample}_composition.txt` | share of the DNA and of the contigs per genus, over every contig, kept or not | `select_contigs` | deliverable |
| `contigs.list` | IDs of the kept contigs | `select_contigs` | deliverable |
| `contigs_decontam.fasta` | the kept contigs | `select_contigs` | intermediate |
| `{sample}_map.bam`, `{sample}_map.bam.bai` | the coverage track: the trimmed Illumina pairs (`illumina`, bowtie2) or the filtered long reads (`longread`, `hybrid`, minimap2) aligned to the draft; in `contigs` mode the contigs mapped to themselves, only because BlobTools needs a BAM, and the coverage BlobTools derives from it means nothing | `map_contigs` | kept in the read-based modes (Qualimap and `contig_depth` read it too); temporary in `contigs` mode |
| `{sample}_depth.tsv` | the read depth of every contig, measured by `samtools coverage` from `{sample}_map.bam`: `#rname` (the contig), `startpos`, `endpos`, `numreads`, `covbases`, `coverage` (the percentage of the contig covered), `meandepth` (the mean read depth), `meanbaseq`, `meanmapq`. The depth in the contig summary is `meandepth`; FunFlux 2.0.0 took it from the coverage column of the BlobTools table, an estimate that can be several times too low with long reads. Not made in `contigs` mode, where no depth is measured | `contig_depth` | intermediate; `contig_summary` reads it |
| `contigs.*.bt2` (six files) | the bowtie2 index of the draft | `index_contigs` | temporary (`illumina`) |
| `blob.blobDB.json`, `blob.{sample}_map.bam.cov` | BlobTools database and coverage file | `blob_json` | temporary |

With the default `decontamination.mode: off`, every contig is kept and still gets a row in
`contig_taxonomy_decisions.tsv`. The decontamination policy and the meaning of each `reason` are
explained on [Decontamination](../analysis/decontamination.md).

### mitochondria/

`mitochondria/` holds the result of the mitochondrion screen.

| Path | Content | Written by | Status |
|---|---|---|---|
| `mitochondrial_screen.tsv` | one row per contig: class, reason, best organelle hit and its query cover, the assembler columns, `concatemer_flag`, and `best_mito_kingdom` (the kingdom of the organism of the best hit: `Fungi`, `Metazoa` or `Viridiplantae`; empty when that cannot be told) | `screen_mitochondria` | deliverable |
| `confirmed_mitochondrial.ids` | the confirmed contig IDs, one per line; empty unless `mitochondria.mode: separate_confirmed` | `screen_mitochondria` | intermediate |

The columns, classes and reasons are explained on [Mitochondria](../analysis/mitochondria.md).

### eval/

`eval/` holds the assembly QC: QUAST, BUSCO and, in the read-based modes, Qualimap, together with
the telomere search, the adapter screen and the contig summary.

| Path | Content | Written by | Status | Modes |
|---|---|---|---|---|
| `quast/` | `report.tsv`, `report.txt`, `report.html`, `report.pdf`, `report.tex`, `transposed_report.*`, `basic_stats/`, `quast.log` | `genome_assembly_evaluation` | deliverable | all |
| `busco/short_summary.specific.<lineage>.{sample}.txt`, `.json` | BUSCO result for the most specific lineage | `completeness` | deliverable | all |
| `busco/short_summary.generic.<lineage>.{sample}.txt`, `.json` | BUSCO result for the broad lineage; only with `busco.lineage: auto` | `completeness` | deliverable | all |
| `busco/run_<lineage>/` | `full_table.tsv` (every expected gene: complete, duplicated, fragmented or missing, and where), `missing_busco_list.tsv` | `completeness` | deliverable | all |
| `busco/auto_lineage/`, `busco/logs/`, `busco/tmp/` | BUSCO's lineage-placement runs, its logs, its `tmp/` folder | `completeness` | intermediate | all |
| `busco/busco_tmp/` | the lineage datasets downloaded for this sample; only with `busco.lineage: auto` | `completeness` | temporary | all |
| `hmmer_output/`, `busco_sequences/` in the BUSCO tree | per-gene search results and gene sequences | `completeness` | scratch | all |
| `qualimap/` | Qualimap on `contaminants/{sample}_map.bam`, the sample's reads on the draft: `qualimapReport.html`, `genome_results.txt` (mapping rate, mean depth per contig), `raw_data_qualimapReport/`; an empty folder when Qualimap failed | `map_evaluation` | deliverable | illumina, longread, hybrid |
| `telomere_candidates.tsv` | every telomere repeat unit weighed for the assembly, the units found at the contig ends and those of the tidk database and of the literature: copies, contig ends in the telomere orientation and in the opposite one, score, whether the unit is known, and which one was chosen | `telomere_units` | deliverable | all |
| `telomere_unit.txt` | one line: the telomere repeat unit of the assembly in canonical form, or `none` | `telomere_units` | intermediate; `contig_summary` reads it | all |
| `telomere_ends.tsv` | per contig end, tidk's call for the chosen unit: status, the array's coordinates and length, its distance from the end, its strand | `telomere_units` | deliverable; `contig_summary` reads it | all |
| `adapter_query.fasta` | the adapter and barcode sequences the screen searches for | `blast_adapters` | temporary | all |
| `adapter_blast.tsv` | BLAST's hits of the adapter and barcode sequences on the draft contigs | `blast_adapters` | intermediate | all |
| `adapter_screen.tsv` | one row per adapter locus on the draft contigs: contig, contig length, start and end of the locus, its distance to the nearer contig end, `end` or `internal`, the longest hit and its identity, and the adapters that matched; the header line alone when none was found | `screen_adapters` | deliverable; `contig_summary` reads it | all |
| `adapter_screen_mqc.tsv` | the sample's row of the Adapters table of the report | `screen_adapters` | deliverable; MultiQC reads it | all |
| `contig_summary.tsv` | one row per contig of the draft assembly: length, GC, read depth and depth relative to the median, the assembler's depth and circularity, the telomere array at each end (copies) and its distance from the end, the telomere unit, rDNA units, BUSCO markers, genes, genus, the screens' calls, the adapter loci, the `scaffold_N` name, the verdict and its reason; written after the annotation | `contig_summary` | deliverable | all |
| `contig_summary_mqc.tsv` | the same table behind MultiQC's header lines, rows named `{sample}: {contig}` | `contig_summary` | deliverable; MultiQC reads it | all |

The columns, the verdicts and how the telomere unit is chosen are explained on
[Contig summary](../analysis/contig-summary.md).

### guide/ (hybrid) { #guide }

In `hybrid` mode, `guide/` holds the Illumina guide assembly and its own screen, laid out like the
sample folder itself. The guide is not delivered. QUAST reports its screened version under the
label `{sample}_guide`, and the Illumina pairs that map to it are the selected pairs in
`01.reads/`.

| Path | Content | Written by | Status |
|---|---|---|---|
| `guide/spades/` | SPAdes' working folder for the guide, as `spades/` above | `illumina_assembly` | intermediate |
| `guide/contigs_filt.fasta` | the guide after the 500 bp / 2x filter | `filter_contigs` | intermediate |
| `guide/contaminants/{sample}_map.bam`, `.bai` | the trimmed Illumina pairs aligned to the guide; Qualimap reads it too | `map_guide_contigs` | intermediate |
| `guide/contaminants/filtered.*.bt2`, `selected.*.bt2` | the bowtie2 indexes of the filtered and of the screened guide | `index_guide_contigs`, `index_selected_contigs` | temporary |
| `guide/contaminants/blastout`, `bestscore.blob.blobDB.table.txt`, `contig_taxonomy_decisions.tsv`, `{sample}_composition.txt`, `contigs.list` | the same files as in `contaminants/` above, for the guide | `blast_guide_contigs`, `blob_guide_table`, `select_guide_contigs` | as in `contaminants/` |
| `guide/contaminants/contigs_decontam.fasta` | the screened guide: QUAST reports it as `{sample}_guide`, and the read selection maps to it | `select_guide_contigs` | intermediate |
| `guide/contaminants/blob.blobDB.json`, `blob.{sample}_map.bam.cov`, `selected_map.sam` | BlobTools' working files; the alignments of the read selection | `blob_guide_json`, `map_sel_contigs` | temporary |
| `guide/mitochondria/mitochondrial_screen.tsv`, `confirmed_mitochondrial.ids` | the mitochondrion screen of the guide, from BLAST alone | `screen_guide_mitochondria` | deliverable, intermediate |
| `guide/eval/qualimap/` | Qualimap on the Illumina pairs mapped to the guide (MultiQC row `mapping Illumina QC`) | `map_guide_evaluation` | deliverable |

## 03.taxonomy

`03.taxonomy/` holds the ITS taxonomy, in one folder per sample. The method is described on
[ITS taxonomy](../analysis/its-taxonomy.md).

| Path | Content | Written by | Status |
|---|---|---|---|
| `its_windows.fasta` | one window per rDNA unit, `its.window_flank` bp (2,000) on each side of the 5.8S gene, turned to read SSU-ITS1-5.8S-ITS2-LSU; headers `contig:start-end_strand` | `its_windows` | intermediate |
| `its_windows.tsv` | `window_id`, `contig`, `start`, `end`, `strand`, `anchor_hits`, `best_evalue`, `window_length` | `its_windows` | deliverable |
| `{sample}.full.fasta`, `{sample}.ITS1.fasta`, `{sample}.ITS2.fasta` | ITSx output, one record per rDNA unit, named after its window | `its_extraction` | deliverable |
| `{sample}.summary.txt` | ITSx summary | `its_extraction` | deliverable |
| `{sample}_ITS_copies.txt` | one row per distinct sequence: `region`, `sequence_id` (`-` when rejected), `length`, `copies`, `status` (`ok`, `rejected_too_short`, `rejected_too_long`), `windows` | `its_classification` | deliverable |
| `{sample}_ITS1_taxonomy.txt`, `{sample}_ITS2_taxonomy.txt` | SINTAX result for the distinct sequences (`{sample}_ITS1_1`, ...); no header; record name, taxonomy with a confidence per rank, strand, taxonomy kept at `its_taxonomy_cutoff` | `its_classification` | deliverable |
| `{sample}.ITSx_tmp/` | ITSx working folder | `its_extraction` | scratch |

The full ITS region is listed in the copies table but not classified. A genome without a
recognisable 5.8S gene gives an empty windows file, empty ITS and taxonomy files and a note in
the logs, and the run continues.

## 04.annotation

`04.annotation/` holds the annotation, in one folder per sample: the cleaned and renamed genome,
the repeat masking, the gene prediction and the functional annotation. The methods are described
on [Repeat masking](../analysis/repeat-masking.md),
[Gene prediction](../analysis/gene-prediction.md) and
[Functional annotation](../analysis/functional-annotation.md).

| Path | Content | Written by | Status | Modes |
|---|---|---|---|---|
| `contigs_clean.fasta` | `contigs_final.fasta` after `funannotate clean`, which removes contigs under 500 bp and duplicated contigs shorter than the N50 ([the rules](../analysis/repeat-masking.md#funannotate-clean)); original names | `funannotate_preprocess` | deliverable | all |
| `contigs_sort.fasta` | the same contigs sorted longest first and renamed `scaffold_N`: the genome that is annotated | `funannotate_preprocess` | deliverable | all |
| `contigs_removed_by_clean.txt` | names of the contigs `funannotate clean` removed; an empty file when none | `funannotate_preprocess` | deliverable | all |
| `contigs_mask.fasta` | `contigs_sort.fasta` with repeats in lowercase | `repeat_masking` | intermediate | all |
| `repeatmasking/` | RepeatModeler library and RepeatMasker reports; empty with `masking_method: tantan`; emptied whenever the rule starts | `repeat_masking` | intermediate | all |
| `genemark_check.txt` | date and outcome of the GeneMark-ES test on a 200 kb slice of `contigs_mask.fasta` | `genemark_check` | intermediate | GeneMark |
| `tiberius/tiberius.gtf` | Tiberius's gene models (GTF): one transcript per gene, with exons, CDS, start and stop codons | `tiberius_prediction` | intermediate | Tiberius |
| `tiberius/tiberius.evm.gff3` | the same models as GFF3 for EVidenceModeler, the file `funannotate predict` gets with `--other_gff` | `tiberius_gff3` | intermediate | Tiberius |
| `augustus_model.txt` | what AUGUSTUS does for the sample: the stored AUGUSTUS model it uses as it is, or training on the genome's BUSCO genes (optimised or not); copied to the top of `logs/funannotate_prediction_{sample}.log` | `augustus_model_check` | intermediate | all |
| `funannotate/predict_results/` | the gene prediction: `<Genus_species>_{sample}.gbk`, `.gff3`, `.proteins.fa`, `.cds-transcripts.fa`, `.mrna-transcripts.fa`, `.discrepency.report.txt`, `.error.summary.txt`, and `<model>.parameters.json` (the name of the stored AUGUSTUS model, or `<genus>_<species>_<sample>` in lower case for a model trained on the genome) | `funannotate_prediction` | intermediate | all |
| `funannotate/predict_misc/`, `funannotate/logfiles/` | funannotate's working files and logs | `funannotate_prediction` | intermediate | all |
| `funannotate/annotate_results/` | the final annotation: `<Genus_species>_{sample}.gbk`, `.gff3`, `.tbl`, `.sqn`, `.agp`, `.contigs.fsa`, `.scaffolds.fa`, `.proteins.fa`, `.cds-transcripts.fa`, `.mrna-transcripts.fa`, `.annotations.txt`, `.clusters.txt`, `.stats.json`, `.discrepency.report.txt`, and `Gene2Products.must-fix.txt`, `Gene2Products.need-curating.txt`, `Gene2Products.new-names-passed.txt` | `funannotate_annotation` | deliverable | all |
| `funannotate/annotate_misc/` | working files of funannotate annotate, among them `annotations.secretome.txt` (with SignalP), `genome.fixedproducts` and `iprscan.xml`, funannotate's copy of the InterProScan XML | `funannotate_annotation` | intermediate | all |
| `iprscan/interproscan_out.xml` | InterProScan matches with GO terms and pathways (can be a few GB) | `iprscan_annotation` | deliverable | all |
| `iprscan/tmp/` | InterProScan's working folder, left in place; empty or nearly so after a finished run | `iprscan_annotation` | intermediate | all |
| `eggnog/{sample}.emapper.annotations`, `.emapper.hits`, `.emapper.seed_orthologs` | eggNOG-mapper results | `eggnog_annotation` | deliverable | all |
| `eggnog/eggnog_tmp/` | eggNOG-mapper working folder | `eggnog_annotation` | temporary | all |
| `antismash/` | `index.html` (the cluster report), `antismash.gbk` (read by funannotate annotate), `antismash.json`, `antismash.zip`, `regions.js`, one `<scaffold>.regionNNN.gbk` per region | `antismash_annotation` | deliverable | all |
| `signalp6/` | `prediction_results.txt` (read by funannotate annotate), `output.gff3`, `output.json`, `region_output.gff3`, `processed_entries.fasta`, `proteins.for_signalp.faa`, `.matplotlib/` | `signalp6_annotation` | deliverable | SignalP |
| `.signalp6-{sample}.XXXXXX/` | SignalP working folder; it replaces `signalp6/` only after the result passed its checks | `signalp6_annotation` | scratch | SignalP |
| `.annotation_complete.txt`, `.cleanup_complete.txt` | empty markers: funannotate annotate finished; the cleanup finished | `funannotate_annotation`, `cleanup` | intermediate | all |
| `clean_*/`, `*.log`, `genome.fixedproducts` directly in this folder | funannotate scratch files, if any are left; listed in `logs/cleanup_{sample}.log` | | cleanup | all |

In the file names, `<Genus_species>` is the `species` column of the annotation parameter table
with spaces turned into underscores, for example `Fusarium_graminearum_SAMPLE1.gff3`. The sample
name is also the locus-tag prefix of every gene ID.

## 05.amr

`05.amr/` is written only when the AMR screen is on, and holds its results in one folder per
sample. The method, the target statuses and the columns are explained on
[Antifungal resistance mutations](../analysis/amr.md).

| Path | Content | Written by | Status |
|---|---|---|---|
| `amr_summary.tsv` | one row per target protein of the species: compartment, status (`evaluated`, `not evaluated`, `not found`, `not screened`), positions compared, the verdict (`resistance`, `conflicting`, `unknown`, `no resistance`, `wild type`), the FungAMR mutations with their verdicts and the new mutations found, the drugs with evidence of resistance and of no resistance with their levels, a note | `amr_chroquetas` | deliverable |
| `amr_positions.tsv` | one row per FungAMR position of every copy read position by position: target, compartment, copy, position, reference and found amino acid, call (`FungAMR mutation`, `new mutation`, `no mutation`; or `no alignment` or `position not found` where the position could not be compared), verdict, the drugs with evidence of resistance and of no resistance with their levels in words, the changes with evidence of resistance at that position, ChroQueTas's drug text | `amr_chroquetas` | deliverable |
| `amr_summary_mqc.tsv` | the summary behind MultiQC's header lines, rows named `{sample}: {target}` | `amr_chroquetas` | deliverable; MultiQC reads it |
| `chroquetas/nuclear/`, `chroquetas/mitochondrial/` | ChroQueTas's own output of the two runs, unchanged: its summary and statistics, one table per target copy, the extracted proteins, miniprot's alignments, the BLASTP comparisons, `tmp/`. A run that was not made has no folder: `mitochondrial/` exists only when the species has cytochrome b and its mitochondrial run was made | `amr_chroquetas` | deliverable |

## 06.report

`06.report/` holds the MultiQC report of the run and the files behind it.

| Path | Content | Status |
|---|---|---|
| `multiqc_report.html` | QUAST and BUSCO results, the contig summary table and the Adapters table of every sample of the run, and the AMR table with the AMR screen on; in the read-based modes also fastp (`illumina`, `hybrid`), NanoPlot (`longread`, `hybrid`) and Qualimap | deliverable |
| `multiqc_data/` | the tables behind the plots, for example `multiqc_general_stats.txt`, `multiqc_busco.txt`, `multiqc_quast.txt`, `multiqc_contig_summary.txt`, `multiqc_adapter_screen.txt` | deliverable |
| `multiqc_config.yaml` | the MultiQC configuration of the run, with the sample-name rules | intermediate |

Rows are named after the tool and the sample, for example `QUAST assembly QC | <sample>` and
`BUSCO completeness QC (<lineage>) | <sample>`. A sample has two BUSCO rows with the automatic
lineage and one with a named lineage. The rows of the contig summary table are named
`{sample}: {contig}`, and those of the AMR table `{sample}: {target}`. The rows of each mode are
listed under [Report](../analysis/report.md#sample-names). The report covers the samples of the
current run only and is rewritten when they change.

## logs

`logs/` holds a log for every per-sample rule and sample, and the logs of the steps that run once
for the whole run.

| File | Written by |
|---|---|
| `<rule>_{sample}.log` | every per-sample rule, for example `funannotate_prediction_SAMPLE1.log` |
| `funannotate_db.log`, `normalize_funannotate_db.log` | the funannotate database download and path repair |
| `unite_db.log` | `unite_db` or `unite_db_local` |
| `busco_data.log` | `busco_data`: the download of a named BUSCO lineage |
| `unite_identified.log` | `unite_identified`: the one-line summary of the filter and the two counts |
| `telomere_db.log` | `telomere_db`: the download of the telomere repeat database and its checksum |
| `univec_illumina.log` | `univec_illumina`: the download of UniVec, its checksum and the Illumina records kept |
| `chroquetas_db.log` | `chroquetas_db`: the download of ChroQueTas and its database, the two checksums, the unpacking (AMR screen on) |
| `tiberius_model.log` | `tiberius_model`: the download of the Tiberius `fungi` model, its checksum and the unpacking (Tiberius on) |
| `antismash_db.log` or `antismash_db_local.log` | the antiSMASH database rule of the run |
| `install_signalp6.log` | the SignalP 6 installation |
| `download_phix.log`, `build_phix.log` | the PhiX genome download and its bowtie2 index (`illumina`, `hybrid`) |
| `multiqc.log` | the report |

Snakemake deletes the log of a job before it runs that job again, so each log holds the latest
attempt only. Copy a log you want to keep before you run again. Some messages go only to the logs:

| Only in the log | Log |
|---|---|
| the warnings of the contig selector | `select_contigs_{sample}.log` |
| a QUAST, Qualimap or NanoPlot failure | `genome_assembly_evaluation_{sample}.log`, `map_evaluation_{sample}.log`, `map_guide_evaluation_{sample}.log`, `raw_long_read_qc_{sample}.log`, `filtered_long_read_qc_{sample}.log` |
| a genome without 5.8S | `its_extraction_{sample}.log` |
| the share of PhiX in an Illumina library | `map_phix_{sample}.log` |
| fastp's and Filtlong's read totals | `trim_adapters_{sample}.log`, `filter_long_reads_{sample}.log` |
| how many SPAdes contigs the filter kept | `filter_contigs_{sample}.log` |
| which file became `contigs_all.fasta` | `finalize_assembly_{sample}.log` |
| tidk's messages and the candidate units behind the telomere unit (the `[FunFlux] telomere_units` summary line and its NOTE reach the console too) | `telomere_units_{sample}.log` |
| the one-line summary of the contig summary: contigs per verdict, contig ends with a telomere, median depth, BUSCO lineage | `contig_summary_{sample}.log` |
| ChroQueTas's messages of both runs, among them the targets it could not extract (the `[FunFlux] amr_chroquetas` summary line and its NOTE reach the console too) | `amr_chroquetas_{sample}.log` |

Snakemake writes its own log of each run to `.snakemake/log/` in the launch folder. See
[Troubleshooting](../troubleshooting.md).

## Temporary files and folders a re-run replaces

The files and folders below are deleted during the run, either by Snakemake once no later step
needs them or by the step that made them. A finished run has none of them, so they take space
only while the run goes on.

| Path | Modes |
|---|---|
| `01.reads/{sample}/illumina/{sample}_nophix_R1.fastq`, `_R2.fastq` | illumina, hybrid |
| `01.reads/{sample}/illumina/{sample}_trim_R1.fastq`, `_R2.fastq` | illumina, hybrid |
| `01.reads/{sample}/long/{sample}_trimmed.fastq.gz` | longread, hybrid (ONT samples) |
| `02.assembly/{sample}/eval/adapter_query.fasta` | all |
| `02.assembly/{sample}/contaminants/{sample}_map.bam` and `.bai` | contigs + screen (kept in the read-based modes) |
| `02.assembly/{sample}/contaminants/contigs.*.bt2` | illumina |
| `02.assembly/{sample}/contaminants/blob.blobDB.json`, `blob.{sample}_map.bam.cov` | screen |
| `02.assembly/{sample}/guide/contaminants/filtered.*.bt2`, `selected.*.bt2`, `blob.blobDB.json`, `blob.{sample}_map.bam.cov`, `selected_map.sam` | hybrid |
| `02.assembly/{sample}/polypolish/draft.fasta`, its bwa index files, `align_1.sam`, `align_2.sam`, `filtered_align_1.sam`, `filtered_align_2.sam`, `polished_raw.fasta` | hybrid |
| `02.assembly/{sample}/eval/busco/busco_tmp/` (`busco.lineage: auto`), and BUSCO's `hmmer_output/` and `busco_sequences/` | all |
| `03.taxonomy/{sample}/{sample}.ITSx_tmp/` | all |
| `04.annotation/{sample}/eggnog/eggnog_tmp/` | all |
| `04.annotation/{sample}/.signalp6-{sample}.XXXXXX/`, the SignalP working folder, and `00.resources/signalp6/.signalp6-install-*/` | SignalP |
| a folder in the system's temporary directory, for `its_windows` and `its_classification` | all |
| a folder in the system's temporary directory, for `genemark_check` | GeneMark |

The folders below are outputs that Snakemake deletes before it runs their rule again.

| Folder | Rule |
|---|---|
| `01.reads/{sample}/long/raw_qc/`, `filtered_qc/` | `raw_long_read_qc`, `filtered_long_read_qc` |
| `02.assembly/{sample}/spades/`, `02.assembly/{sample}/guide/spades/` | `illumina_assembly` |
| `02.assembly/{sample}/hifiasm/` or `flye/` | `long_read_assembly` |
| `02.assembly/{sample}/medaka/` | `long_read_consensus` |
| `02.assembly/{sample}/eval/quast/` | `genome_assembly_evaluation` |
| `02.assembly/{sample}/eval/busco/` | `completeness` |
| `02.assembly/{sample}/eval/qualimap/`, `02.assembly/{sample}/guide/eval/qualimap/` | `map_evaluation`, `map_guide_evaluation` |
| `04.annotation/{sample}/funannotate/` | `funannotate_prediction` |
| `04.annotation/{sample}/antismash/` | `antismash_annotation` |
| `04.annotation/{sample}/eggnog/` | `eggnog_annotation` |
| `05.amr/{sample}/chroquetas/` | `amr_chroquetas` |
| `00.resources/antismash/databases/` | `antismash_db` or `antismash_db_local` |
| `06.report/multiqc_data/` | `multiqc` |

!!! warning "Files you put into these folders are lost on a re-run"

    The same holds for `repeatmasking/`, which `repeat_masking` empties itself, and for
    `signalp6/`, which `signalp6_annotation` replaces as a whole. Because funannotate annotate
    writes into `funannotate/`, a re-run of `funannotate_prediction` also deletes
    `annotate_results/`, and the annotation runs again. With `antismash_db_local`, a re-run
    deletes only the links, never your database.

`.snakemake_timestamp` inside such a folder is Snakemake's time stamp for it, so leave it in
place. What triggers a re-run is explained on
[Running FunFlux](running.md#what-makes-snakemake-re-run-a-rule).

## Outside the output folder

A run also writes outside `directories.output_dir`, in the places below.

| Where | What |
|---|---|
| `directories.funannotate_db` | the funannotate database: archive (`funannotate_db_1.8.17_2026-10-08.tar.gz` by default, 2.1 GB) and its `.sha256`, both kept after unpacking; the unpacked database, whose stored paths `normalize_funannotate_db` rewrites in place; 8.1 GB in total |
| `.snakemake/` in the launch folder | Snakemake's records of what was made and how, locks, markers of unfinished jobs, run logs, and the conda environments unless `--conda-prefix` names another folder |
| the BlobTools conda environment | `nodesDB.txt`, built by BlobTools from `nodes.dmp` and `names.dmp` on first use |
| the Medaka conda environment, or `~/.medaka/data` when it cannot be written | the weights of each Medaka model, downloaded by the model check the first time a model is used |
| the pip cache of the user | the wheel pip builds while installing SignalP 6 |

The databases are described on [Reference databases](../getting-started/databases.md), and the
environments and the launch folder on [Running FunFlux](running.md).
