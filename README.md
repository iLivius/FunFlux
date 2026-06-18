# FunFlux
Integrated workflow for fungal short-read genome assembly and annotation.

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.14.6-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13612159.svg)](https://doi.org/10.5281/zenodo.13612159)

---
```bash
__________             _______________
___  ____/___  ___________  ____/__  /___  _____  __
__  /_   _  / / /_  __ \_  /_   __  /_  / / /_  |/_/
_  __/   / /_/ /_  / / /  __/   _  / / /_/ /__>  <
/_/      \__,_/ /_/ /_//_/      /_/  \__,_/ /_/|_|

FunFlux v1.1.0

June 2026
```
---

## Synopsis
`FunFlux` is a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow designed for genome assembly and annotation of fungal short reads sequenced with Illumina technology. Pre-assembled fungal genomes can be analyzed with the bundled `Funnotator` flavour. The workflow includes read preprocessing, assembly, contig selection and decontamination, assembly quality control, genome completeness assessment, ITS extraction and taxonomic assignment, repeat masking, gene prediction, functional annotation, and summary reporting.

## Table of Contents
- [Rationale](#rationale)
- [Description](#description)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running FunFlux](#running-funflux)
- [Output](#output)
- [Acknowledgements](#acknowledgements)
- [Citation](#citation)
- [References](#references)

## Rationale
The analysis of fungal whole-genome sequencing (WGS) data involves a complex series of bioinformatic steps that can be challenging to execute manually. This process is time-consuming, prone to errors, and difficult to reproduce. `FunFlux` addresses these problems by providing an automated [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow for fungal genome assembly and annotation.

`FunFlux` is designed to streamline the annotation process with [funannotate](https://github.com/nextgenusfs/funannotate) in the absence of RNA sequencing evidence. It relies on both *ab initio* annotation and protein FASTA sequences from organisms of the same species or genus to enhance the accuracy of gene prediction and annotation.

## Description
Here's a breakdown of the `FunFlux` workflow:

01. **Preprocessing:**
    * Raw reads are checked for Illumina phiX contamination using [bowtie2](https://github.com/BenLangmead/bowtie2).
    * Adapters are removed and reads are filtered using [fastp](https://github.com/OpenGene/fastp).

02. **Assembly:**
    * Filtered reads are assembled into contigs with [SPAdes](https://github.com/ablab/spades).

03. **QC, Decontamination, Completeness Assessment, and ITS extraction:**
    * Contigs are filtered based on minimum length and coverage.
    * Filtered reads are mapped back to contigs using [bowtie2](https://github.com/BenLangmead/bowtie2) and [samtools](https://github.com/samtools/samtools). The resulting BAM file is analyzed with [QualiMap](http://qualimap.conesalab.org/).
    * Local alignments of contigs are performed against the [NCBI core nt](https://ftp.ncbi.nlm.nih.gov/blast/db/) database using [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).
    * Contaminant contigs are checked with [BlobTools](https://github.com/DRL/blobtools). Contig selection is handled by a dedicated taxonomy selector script. The selector can keep all contigs, keep the most abundant assigned genus, include specified genera, exclude specified genera, and optionally discard `no-hit` contigs.
    * Genome assembly quality is evaluated with [QUAST](https://github.com/ablab/quast).
    * Genome completeness is assessed with [BUSCO](https://busco.ezlab.org/).
    * ITS markers are detected and extracted with [ITSx](https://microbiology.se/software/itsx/).
    * ITS taxonomic assignment is performed with the [SINTAX](https://www.drive5.com/sintax/) classifier in [VSEARCH](https://github.com/torognes/vsearch) using the [UNITE](https://unite.ut.ee/repository.php) database.

04. **Gene Prediction:**

    `FunFlux` is optimized to leverage [funannotate](https://github.com/nextgenusfs/funannotate) when RNA sequencing data is not available. Instead, it uses external protein evidence and *ab initio* predictors to produce fungal gene models. The workflow splits the previous monolithic prediction step into:

    ```text
    funannotate_preprocess -> repeat_masking -> funannotate_prediction
    ```

    This split makes repeat masking replaceable without changing the downstream funannotate prediction and annotation logic.

    - Preprocessing the genome assembly
        - N50 calculation and contig duplication checking are performed by `funannotate clean`.
        - Contigs are sorted and headers are renamed with `funannotate sort`.

    - Repeat masking
        - The default strategy is direct [tantan](https://gitlab.com/mcfrith/tantan) softmasking.
        - The optional advanced strategy runs [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) and [RepeatMasker](https://www.repeatmasker.org/) with the de novo repeat library produced by the former.

    - Incorporating protein evidence
        - [DIAMOND](https://github.com/bbuchfink/diamond) is used by `funannotate` to search for homology between the genome and provided protein sequences from related taxa, as well as the `UniProt` database bundled in the configured `funannotate` database snapshot. Matches are refined by [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate).

    - *Ab initio* gene prediction
        - [GeneMark-ES](https://genemark.bme.gatech.edu/gmes_instructions.html) is made available to `funannotate` and contributes self-trained *ab initio* gene predictions from the genome sequence.
        - [BUSCO](https://busco.ezlab.org/) conserved genes are passed to [Augustus](https://github.com/Gaius-Augustus/Augustus) to improve training.
        - [SNAP](https://github.com/KorfLab/SNAP), [GlimmerHMM](https://ccb.jhu.edu/software/glimmerhmm/), and other `funannotate`-supported predictors contribute to consensus model building.

    - Combining and refining predictions
        - [EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler) combines evidence into final gene models.
        - tRNA genes are predicted with [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE).
        - NCBI-compatible annotation files are generated by `funannotate`.

05. **Gene Annotation:**

    Gene annotation integrates multiple tools and culminates in a final `funannotate annotate` step:

    - [InterProScan](https://github.com/ebi-pf-team/interproscan) is expected as an external local installation.
    - [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) is used for orthology and functional annotation.
    - [antiSMASH](https://github.com/antismash/antismash) detects secondary metabolite biosynthetic gene clusters. The `antiSMASH` database directory is a temporary Snakemake output and should be removed after antiSMASH jobs complete.
    - The configured funannotate database snapshot contains resources such as `UniProt`, `MEROPS`, `dbCAN`, `Pfam`, `GO`, `MIBiG`, `InterPro`, and `BUSCO` outgroups.

06. **Report:**
    * Results are parsed and aggregated with [MultiQC](https://github.com/MultiQC/MultiQC).

07. **Funnotator flavour:**
    * `Funnotator` is the annotation-only `FunFlux` flavour for pre-assembled fungal genome FASTA files. It skips Illumina preprocessing and `SPAdes` assembly, then applies the same fungal annotation logic where relevant.

[⬆ Back to Table of Contents](#table-of-contents)

## Installation
`FunFlux` automatically downloads most Conda-managed dependencies and several workflow databases. Some external databases and licensed tools still require manual installation before running the workflow.

1. **Download FunFlux:**

    ```bash
    git clone https://github.com/iLivius/FunFlux.git
    ```

2. **Install Snakemake:**

    `FunFlux` relies on [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to manage workflow execution.

    ```bash
    conda create -c conda-forge -c bioconda -n snakemake snakemake
    conda activate snakemake
    ```

3. **Databases and external software:**

    * `NCBI core nt` database:

        ```bash
        rsync --list-only rsync://ftp.ncbi.nlm.nih.gov/blast/db/core_nt.*.gz | grep '.tar.gz' | awk '{print "ftp.ncbi.nlm.nih.gov/blast/db/" $NF}' > nt_links.list
        cat nt_links.list | parallel -j4 'rsync -h --progress rsync://{} .'
        find . -name '*.gz' | parallel -j4 'echo {}; tar -zxf {}'

        wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        tar -zxvf taxdump.tar.gz

        wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'
        tar -zxvf taxdb.tar.gz

        wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        gunzip nucl_gb.accession2taxid.gz
        ```

        *NOTE: the complete NCBI core nt database and taxonomy-related files require more than 200 GB of disk space. It is not needed when using Funnotator on already assembled FASTA files.*

    * `UNITE` database:

        Manual download is not required for the standard workflow. The config contains:

        ```yaml
        links:
          unite_its_link: https://s3.hpc.ut.ee/plutof-public/original/338a1413-6039-4e00-b5cf-410346a1e366.gz
        ```

        The workflow downloads and decompresses this file automatically into:

        ```text
        03.post-processing/ITS_extraction/unite_its_sintax.fasta
        ```

    * `eggNOG diamond` database:

        ```bash
        conda create -n eggnog-mapper eggnog-mapper=2.1.13
        conda activate eggnog-mapper
        mkdir /data/eggnog_db
        download_eggnog_data.py --data_dir /data/eggnog_db -y
        ```

        *NOTE: the eggNOG database requires roughly 50 GB of disk space.*

    * Download and set up `GeneMark-ES/ET`:

        - Visit the [GeneMark](http://topaz.gatech.edu/GeneMark/license_download.cgi) download page.
        - Download `GeneMark-ES/ET`.
        - Change Perl script shebangs if necessary:

            ```bash
            cd /path/to/gmes_linux_64_4
            find . -type f -name "*.pl" -exec sed -i '1s|^#!/usr/bin/perl|#!/usr/bin/env perl|' {} +
            ./gmes_petap.pl
            ```

    * Download and set up `InterProScan`:

        The version tested was `v5.77-108.0`, which is not downloaded by the workflow. The official download instructions are available [here](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html). Download the tested archive and checksum file:

        ```bash
        wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.77-108.0/interproscan-5.77-108.0-64-bit.tar.gz
        wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.77-108.0/interproscan-5.77-108.0-64-bit.tar.gz.md5
        md5sum -c interproscan-5.77-108.0-64-bit.tar.gz.md5
        ```

        Then extract and initialize `InterProScan`:

        ```bash
        tar -pxvzf interproscan-5.77-108.0-64-bit.tar.gz
        cd interproscan-5.77-108.0
        python3 setup.py -f interproscan.properties
        ./interproscan.sh
        ```

    * Repeat masking tools:

        `tantan`, `RepeatModeler`, and `RepeatMasker` are installed through the workflow Conda environment. The advanced masking mode uses the custom RepeatModeler library with `RepeatMasker -lib`, so no separate curated RepeatMasker database is required for that mode.

[⬆ Back to Table of Contents](#table-of-contents)

## Configuration
Before running `FunFlux`, edit `config/config.yaml`. The file is organized into `links`, `directories`, `files`, `resources`, and `parameters`.

- `links`

    - [phix_link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015): PhiX genome reference used as Illumina sequencing control.
    - [funannotate_link](https://zenodo.org/records/18271295/files/funannotate_db_1.8.15_2025-12-15.tar.gz?download=1): URL to a frozen funannotate database snapshot. This is a database archive, not the funannotate executable version.
    - [unite_its_link](https://s3.hpc.ut.ee/plutof-public/original/338a1413-6039-4e00-b5cf-410346a1e366.gz): URL to the compressed UNITE SINTAX FASTA used for ITS classification.

- `directories`

    - **input_dir**: Directory containing paired-end FASTQ reads. Requirements:
        1. Files can only have `fastq`, `fq`, `fastq.gz`, or `fq.gz` extensions.
        2. All files in a run must use the same extension.
        3. Sample names must not contain underscores or these characters: `_*#@%^/! ?&:;|<>`.
        4. Use `_R1` and `_R2` to define paired-end reads, for example `strain-42_R1.fastq.gz` and `strain-42_R2.fastq.gz`.

        For `Funnotator`, provide assembled FASTA files instead. FASTA files can use `fasta`, `fa`, or `fna` extensions and sample names must not contain underscores.

    - **output_dir**: Directory where output files, `.snakemake` metadata, and Conda environments are stored. Reusing the same output directory avoids reinstalling environments.
    - **blast_db**: Path to the NCBI core nt database and taxonomy files.
    - **eggnog_db**: Path to the eggNOG-mapper database.
    - **genemark_dir**: Path to the `gmes_linux_64_4` directory.
    - **funannotate_db**: Path where the funannotate database snapshot is installed or already available.

- `files`

    - **annotation_params**: Path to a tab-delimited annotation parameter file. An example is provided in `config/annotation_parameters.tsv`.

        | #Sample        | Species                    | Proteins              | Model                |
        |----------------|----------------------------|-----------------------|----------------------|
        | ARSEF3097      | Beauveria bassiana         | /path/to/proteins.faa | fusarium_graminearum |
        | strain-42      | Lecanicillium fungicola    | /path/to/proteins.faa | fusarium_graminearum |

    - **iprscan**: Path to the InterProScan shell script.

- `resources`

    - **threads**: Maximum CPUs passed to individual tools inside rules. Some tools are capped internally where higher values are not useful or can be unstable.
    - **ram_gb**: Maximum RAM value used by memory-aware tools such as SPAdes and QualiMap.

    `--cores` is Snakemake's scheduler limit. `resources: threads` controls tool-level thread arguments.

- `parameters`

    **Decontamination**

    ```yaml
    decontamination:
      mode: off
      discard_no_hit: true
      include_genera:
      include_genera_by_sample:
      exclude_genera:
      exclude_genera_file:
      sample_overrides:
    ```

    Available modes:

    - `off`: keep all contigs. `discard_no_hit` is ignored.
    - `auto`: keep the most abundant assigned genus.
    - `include`: keep only listed genera.
    - `exclude`: remove listed genera.

    `discard_no_hit: true` removes BLAST `no-hit` contigs only when the mode is `auto`, `include`, or `exclude`.
    In `auto` and `include` modes, the selector can treat selected genus aliases and retained legacy prefixes as equivalent. This is mainly a safeguard against false contig removal when BLAST/BlobTools assigns related or recently reclassified genera inconsistently. Although FunFlux targets fungal genomes, bacterial genera may appear here because bacterial contamination can occur in fungal WGS assemblies. Alias-based decisions are recorded in `contig_taxonomy_decisions.tsv` with reasons such as `auto_genus_alias` or `included_genus_alias`. `exclude` mode remains exact. These aliases are heuristic safeguards, not a formal taxonomic reconciliation system.

    Genera can be supplied directly:

    ```yaml
    exclude_genera: Acidovorax;Pseudomonas;Sphingomonas
    ```

    or through a one-genus-per-line file:

    ```yaml
    exclude_genera_file: /path/to/exclude_genera.txt
    ```

    Optional sample overrides use a tab-separated file:

    ```text
    sample<TAB>mode<TAB>include_genera<TAB>exclude_genera<TAB>discard_no_hit
    strain-42<TAB>exclude<TAB><TAB><TAB>true
    ```

    Use real tab characters, not the literal string `<TAB>`. A sample-specific include or exclude list replaces the global list for that sample.

    **ITS taxonomy**

    ```yaml
    its_taxonomy_cutoff: 0.8
    ```

    **Repeat masking**

    ```yaml
    masking_method: tantan
    repeatmodeler_quick: true
    repeatmodeler_ltrstruct: false
    ```

    To use the advanced RepeatModeler + RepeatMasker strategy, change only `masking_method`:

    ```yaml
    masking_method: repeatmodeler_repeatmasker
    ```

    Available masking methods:

    - `tantan`: default lightweight softmasking.
    - `repeatmodeler_repeatmasker`: runs `BuildDatabase -> RepeatModeler -> RepeatMasker -lib <RepeatModeler library> -xsmall`.

    `repeatmodeler_quick: true` adds `RepeatModeler -quick`. `repeatmodeler_ltrstruct: true` adds `RepeatModeler -LTRStruct`.

[⬆ Back to Table of Contents](#table-of-contents)

## Running FunFlux
`FunFlux` can be executed as a Snakemake workflow.

```bash
conda activate snakemake
snakemake --configfile config/config.yaml --sdm conda --cores 12 --jobs 2
```

If you resume an interrupted run, keep using the same configuration file. When Snakemake reports incomplete output after a stopped job, rerun with `--rerun-incomplete`. If you intentionally updated workflow code but want to continue based only on file timestamps, add `--rerun-triggers mtime`.

```bash
snakemake --configfile config/config.yaml --unlock
snakemake --configfile config/config.yaml --sdm conda --cores 12 --jobs 2 --rerun-triggers mtime --rerun-incomplete
```

To analyze pre-assembled fungal genomes with `Funnotator`:

```bash
conda activate snakemake
snakemake --snakefile workflow/Funnotator --configfile config/config.yaml --sdm conda --cores 24 --jobs 2
```

After a successful run, optional cleanup of bulky intermediate files can be inspected with:

```bash
workflow/scripts/clean_funflux_output.sh --target /path/to/output_dir
```

To actually remove the listed files:

```bash
workflow/scripts/clean_funflux_output.sh --run --target /path/to/output_dir
```

The cleanup script is dry-run by default and refuses targets that do not look like `FunFlux` or `Funnotator` output directories.

[⬆ Back to Table of Contents](#table-of-contents)

## Output
Here's a breakdown of the sub-directories created by `FunFlux` within the main output folder. `Funnotator` produces a similar but simplified annotation-only output.

```text
├── 01.pre-processing
├── 02.assembly
├── 03.post-processing
├── 04.annotation
├── logs
└── report
```

- `01.pre-processing`: QC and statistics of raw and trimmed reads, produced by [fastp](https://github.com/OpenGene/fastp) v1.0.1.

- `02.assembly`: Output from [SPAdes](https://github.com/ablab/spades) v4.2.0. This directory contains raw contigs, filtered contigs, and selected/decontaminated contigs.

- `03.post-processing`: Contains:
    - **mapping_evaluation**: [QualiMap](http://qualimap.conesalab.org/) v2.3 output.
    - **contaminants**: BLAST+ v2.16.0 and BlobTools v1.1.1 decontamination output, including genus composition and `contig_taxonomy_decisions.tsv`.
    - **assembly_evaluation**: [QUAST](https://github.com/ablab/quast) v5.3.0 output.
    - **completeness_evaluation**: [BUSCO](https://busco.ezlab.org/) v6.0.0 output from `--auto-lineage-euk`.
    - **ITS_extraction**: [ITSx](https://microbiology.se/software/itsx/) v1.1.3 output and [VSEARCH](https://github.com/torognes/vsearch) v2.30.0 SINTAX classification against the automatically downloaded UNITE reference.

- `04.annotation`: Contains:
    - **repeatmasking**: [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) v2.0.8 and [RepeatMasker](https://www.repeatmasker.org/) v4.2.3 output, present inside each sample when `repeatmodeler_repeatmasker` is selected.
    - **iprscan**: [InterProScan](https://github.com/ebi-pf-team/interproscan) v5.77-108.0 XML output.
    - **eggnog**: [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) v2.1.13 annotation output.
    - **antismash**: [antiSMASH](https://github.com/antismash/antismash) v8.0.4 secondary metabolite output.
    - **funannotate**: Prediction and annotation directories from [funannotate](https://github.com/nextgenusfs/funannotate) v1.8.17.

        ```text
        ├── annotate_misc
        ├── annotate_results
        ├── logfiles
        ├── predict_misc
        └── predict_results
        ```

- `report`: [MultiQC](https://github.com/MultiQC/MultiQC) v1.33 report aggregating fastp, QualiMap, QUAST, BUSCO, and other supported outputs.

[⬆ Back to Table of Contents](#table-of-contents)

## Acknowledgements
This work was originally supported by the [BeXyl project](https://cordis.europa.eu/project/id/101060593) (Beyond Xylella, Integrated Management Strategies for Mitigating *Xylella fastidiosa* impact in Europe), funded under the HORIZON-CL6-2021-FARM2FORK-01-04 programme (grant agreement No. 101060593).

## Citation
If you use `FunFlux`, please cite:

Antonielli, L., Brader, G., & Compant, S. (2024). FunFlux: Integrated workflow for fungal genome assembly and annotation. Zenodo. https://doi.org/10.5281/zenodo.13612159

## References
01. Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology, 19(5), 455-477. https://doi.org/10.1089/cmb.2012.0021

02. Bengtsson-Palme, J., Ryberg, M., Hartmann, M., Branco, S., Wang, Z., Godhe, A., De Wit, P., Sánchez-García, M., Ebersberger, I., de Sousa, F., Amend, A., Jumpponen, A., Unterseher, M., Kristiansson, E., Abarenkov, K., Bertrand, Y. J. K., Sanli, K., Eriksson, K. M., Vik, U., ... Nilsson, R. H. (2013). Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods in Ecology and Evolution, 4(10), 914-919. https://doi.org/10.1111/2041-210X.12073

03. Blin, K., et al. (2025). antiSMASH 8.0: extended gene cluster detection capabilities and analyses of chemistry, enzymology and regulation. Nucleic Acids Research, 53(W1), W32-W38. https://doi.org/10.1093/nar/gkaf334

04. Blum, M., Chang, H.-Y., Chuguransky, S., Grego, T., Kandasaamy, S., Mitchell, A., Nuka, G., Paysan-Lafosse, T., Qureshi, M., Raj, S., Richardson, L., Salazar, G. A., Williams, L., Bork, P., Bridge, A., Gough, J., Haft, D. H., Letunic, I., Marchler-Bauer, A., ... Finn, R. D. (2021). The InterPro protein families and domains database: 20 years on. Nucleic Acids Research, 49(D1), D344-D354. https://doi.org/10.1093/nar/gkaa977

05. Borodovsky, M., & Lomsadze, A. (2011). Eukaryotic Gene Prediction Using GeneMark.hmm-E and GeneMark-ES. Current Protocols in Bioinformatics, Unit 4.6. https://doi.org/10.1002/0471250953.bi0406s35

06. Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12(1), 59-60. https://doi.org/10.1038/nmeth.3176

07. Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10, 421. https://doi.org/10.1186/1471-2105-10-421

08. Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Molecular Biology and Evolution, 38(12), 5825-5829. https://doi.org/10.1093/molbev/msab293

09. Challis, R., Richards, E., Rajan, J., Cochrane, G., & Blaxter, M. (2020). BlobToolKit - Interactive Quality Assessment of Genome Assemblies. G3 Genes|Genomes|Genetics, 10(4), 1361-1374. https://doi.org/10.1534/g3.119.400908

10. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560

11. Edgar, R. C. (2016). SINTAX: A simple non-Bayesian taxonomy classifier for 16S and ITS sequences. bioRxiv. https://doi.org/10.1101/074161

12. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. https://doi.org/10.1093/bioinformatics/btw354

13. Flynn, J. M., Hubley, R., Goubert, C., Rosen, J., Clark, A. G., Feschotte, C., & Smit, A. F. (2020). RepeatModeler2 for automated genomic discovery of transposable element families. Proceedings of the National Academy of Sciences, 117(17), 9451-9457. https://doi.org/10.1073/pnas.1921046117

14. Frith, M. C. (2011). A new repeat-masking method enables specific detection of homologous sequences. Nucleic Acids Research, 39(4), e23. https://doi.org/10.1093/nar/gkq1212

15. Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072-1075. https://doi.org/10.1093/bioinformatics/btt086

16. Haas, B. J., Salzberg, S. L., Zhu, W., Pertea, M., Allen, J. E., Orvis, J., White, O., Buell, C. R., & Wortman, J. R. (2008). Automated eukaryotic gene structure annotation using EVidenceModeler and the Program to Assemble Spliced Alignments. Genome Biology, 9(1), R7. https://doi.org/10.1186/gb-2008-9-1-r7

17. Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S. K., Cook, H., Mende, D. R., Letunic, I., Rattei, T., Jensen, L. J., von Mering, C., & Bork, P. (2019). eggNOG 5.0. Nucleic Acids Research, 47(D1), D309-D314. https://doi.org/10.1093/nar/gky1085

18. Jonathan M. Palmer, & Jason Stajich. (2020). Funannotate v1.8.1: Eukaryotic genome annotation [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.4054262

19. Jones, P., Binns, D., Chang, H.-Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G., Pesseat, S., Quinn, A. F., Sangrador-Vegas, A., Scheremetjew, M., Yong, S.-Y., Lopez, R., & Hunter, S. (2014). InterProScan 5: Genome-scale protein function classification. Bioinformatics, 30(9), 1236-1240. https://doi.org/10.1093/bioinformatics/btu031

20. Köster, J., & Rahmann, S. (2012). Snakemake - A scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520-2522. https://doi.org/10.1093/bioinformatics/bts480

21. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357-359. https://doi.org/10.1038/nmeth.1923

22. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079. https://doi.org/10.1093/bioinformatics/btp352

23. Nilsson, R. H., Larsson, K.-H., Taylor, A. F. S., Bengtsson-Palme, J., Jeppesen, T. S., Schigel, D., Kennedy, P., Picard, K., Glöckner, F. O., Tedersoo, L., Saar, I., Kõljalg, U., & Abarenkov, K. (2019). The UNITE database for molecular identification of fungi: Handling dark taxa and parallel taxonomic classifications. Nucleic Acids Research, 47(D1), D259-D264. https://doi.org/10.1093/nar/gky1022

24. Abarenkov, K., et al. (2024). The UNITE database for molecular identification and taxonomic communication of fungi and other eukaryotes: sequences, taxa and classifications reconsidered. Nucleic Acids Research, 52(D1), D791-D797. https://doi.org/10.1093/nar/gkad1039

25. Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: Advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics, 32(2), 292-294. https://doi.org/10.1093/bioinformatics/btv566

26. Rawlings, N. D., Waller, M., Barrett, A. J., & Bateman, A. (2014). MEROPS. Nucleic Acids Research, 42(D1), D503-D509. https://doi.org/10.1093/nar/gkt953

27. Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: A versatile open source tool for metagenomics. PeerJ, 4, e2584. https://doi.org/10.7717/peerj.2584

28. Smit, A. F. A., Hubley, R., & Green, P. RepeatMasker Open-4.0. http://www.repeatmasker.org

29. Stanke, M., Keller, O., Gunduz, I., Hayes, A., Waack, S., & Morgenstern, B. (2006). AUGUSTUS. Nucleic Acids Research, 34(Web Server issue), W435-W439. https://doi.org/10.1093/nar/gkl200

30. Tegenfeldt, F., Kuznetsov, D., Manni, M., Berkeley, M., Zdobnov, E. M., & Kriventseva, E. V. (2025). OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes. Nucleic Acids Research, 53(D1), D516-D522. https://doi.org/10.1093/nar/gkae987

31. The UniProt Consortium. (2023). UniProt: The Universal Protein Knowledgebase in 2023. Nucleic Acids Research, 51(D1), D523-D531. https://doi.org/10.1093/nar/gkac1052

32. Zheng, J., Ge, Q., Yan, Y., Zhang, X., Huang, L., & Yin, Y. (2023). dbCAN3: Automated carbohydrate-active enzyme and substrate annotation. Nucleic Acids Research, 51(W1), W115-W121. https://doi.org/10.1093/nar/gkad328
