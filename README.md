# FunFlux
<<<<<<< HEAD
FunFlux: A dedicated workflow for fungal genome assembly from short reads, decontamination, completeness validation, and comprehensive gene annotation.
=======
Integrated workflow for fungal genome assembly and annotation.

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.16.0-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)

```bash
__________             _______________              
___  ____/___  ___________  ____/__  /___  _____  __
__  /_   _  / / /_  __ \_  /_   __  /_  / / /_  |/_/
_  __/   / /_/ /_  / / /  __/   _  / / /_/ /__>  <  
/_/      \__,_/ /_/ /_//_/      /_/  \__,_/ /_/|_|  
                                                                                                                                      
FunFlux v1.0.3

August 2024
```

## Authors and Contributors
[AIT Austrian Institute of Technology, Center for Health & Bioresources](https://www.ait.ac.at/en/research-topics/bioresources)

- Livio Antonielli
- Günter Brader
- Stéphane Compant

## Synopsis
`FunFlux` is a [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow designed for the genome assembly and annotation of fungal short reads sequenced with Illumina technology.  It also supports the analysis of pre-assembled contigs. The workflow includes features such as contig selection and decontamination, genome completeness assessment, ITS extraction with taxonomic assignment, and precise gene prediction and annotation.

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
The analysis of fungal whole-genome sequencing (WGS) data involves a complex series of bioinformatic steps that can be challenging to execute manually. This process is often time-consuming, prone to errors, and difficult to reproduce. `FunFlux` addresses these challenges by offering a comprehensive and automated [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow specifically designed for fungal genomic data analysis.

`FunFlux` is designed to streamline the annotation process with [funannotate](https://github.com/nextgenusfs/funannotate), even in the absence of RNA sequencing evidence. It relies on *ab initio* annotation and incorporates protein FASTA sequences from organisms of the same species or genus to enhance the accuracy of gene prediction and annotation. 

## Description
Here's a breakdown of the `FunFlux` workflow:

01. **Preprocessing:**
    * Raw reads are checked for Illumina phiX contamination using [bowtie2](https://github.com/BenLangmead/bowtie2).

    * Adapters are removed and reads are filtered using [fastp](https://github.com/OpenGene/fastp).

02. **Assembly:**
    * Filtered reads are assembled into contigs with [SPAdes](https://github.com/ablab/spades).

03. **QC, Decontamination, Completeness Assessment, and ITS extraction:**
    * Contigs are filtered based on a minimum length of 500 bp and a coverage of 2x.
    * Filtered reads are mapped back to contigs using [bowtie2](https://github.com/BenLangmead/bowtie2) and [samtools](https://github.com/samtools/samtools). The resulting BAM file is analyzed with [QualiMap](http://qualimap.conesalab.org/).
    * Local alignments of contigs are performed against the [NCBI core nt](https://ftp.ncbi.nlm.nih.gov/blast/db/) database using [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/).
    * Contaminant contigs are checked with [BlobTools](https://github.com/DRL/blobtools). Unless otherwise specified (see [configuration](#configuration) section for more details), the output of this step will be parsed automatically to discard contaminants based on the relative taxonomic composition of the contigs.   
    * Genome assembly quality is evaluated with [Quast](https://github.com/ablab/quast).
    * Genome completeness is assessed with [BUSCO](https://busco.ezlab.org/) using taxon-specific markers.
    * ITS markers are detected and extracted with [ITSx](https://microbiology.se/software/itsx/).
    * ITS2 taxonomic assignment is performed with [SINTAX](https://www.drive5.com/sintax/) re-implemented in [VSEARCH](https://github.com/torognes/vsearch) using the [UNITE](https://unite.ut.ee/repository.php) database as reference. 

04. **Gene Prediction:**

    `FunFlux` is optimized to leverage the [funannotate](https://github.com/nextgenusfs/funannotate) pipeline in cases where RNA sequencing data is not available. Instead, it utilizes external protein evidence along with robust *ab initio* prediction methods to produce accurate gene models for fungal genomes. Below is a step-by-step breakdown of the workflow:

    - Preprocessing the genome assembly 
        - N50 calculation and contig duplication checking: As part of the cleaning process, the N50 value is calculated, and contigs shorter than this value are checked for duplication. Only unique, non-redundant contigs are retained, ensuring that the assembly is as clean and representative as possible.

        - Sorting and renaming FASTA headers: The assembled contigs are sorted by length and headers are renamed to ensure compatibility with follow-up tools. 
       
        - Repeat masking: Before gene prediction, the genome assembly is softmasked using the [tantan](https://gitlab.com/mcfrith/tantan) software to obscure repetitive elements, which helps in preventing spurious gene predictions in these regions.

    - Incorporating protein evidence 
        - Protein alignment: [DIAMOND](https://github.com/bbuchfink/diamond) is used to quickly search for homologies between the genome and provided protein sequences of closely related taxa, as well as the [UniProt](https://www.uniprot.org/) database. These matches are then refined with [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate), which aligns the protein sequences to the genome with high precision, providing evidence for gene structures.

    - *Ab initio* gene prediction
        - [GeneMark-ES](https://genemark.bme.gatech.edu/gmes_instructions.html): This tool performs self-training on the genome sequence to predict genes without the need for external training data, making it especially useful for identifying genes in regions lacking homology-based evidence.

    - Ortholog detection and model training
        - [BUSCO](https://busco.ezlab.org/): Based on conserved orthologous genes, it provides high-quality evidence for training gene prediction tools. Conserved genes are passed to [Augustus](https://github.com/Gaius-Augustus/Augustus) to improve its predictive accuracy.

        - [Augustus](https://github.com/Gaius-Augustus/Augustus) training: It works with the closest taxon model available, as well as the evidence from [BUSCO](https://busco.ezlab.org/), [DIAMOND](https://github.com/bbuchfink/diamond)/[Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate), and the outputs from other *ab initio* predictors like [SNAP](https://github.com/KorfLab/SNAP) and [GlimmerHMM](https://ccb.jhu.edu/software/glimmerhmm/). This comprehensive training enables Augustus to generate highly accurate gene predictions.

    - Combining predictions with [EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler)
        - [EVidenceModeler (EVM)](https://github.com/EVidenceModeler/EVidenceModeler): The predictions from various *ab initio* tools, such as [Augustus](https://github.com/Gaius-Augustus/Augustus), [SNAP](https://github.com/KorfLab/SNAP), [GlimmerHMM](https://ccb.jhu.edu/software/glimmerhmm/) are combined to generate consensus gene models.

    - Refining steps
        - Gene model filtering: The gene models generated by [EVM]((https://github.com/EVidenceModeler/EVidenceModeler)) are subjected to further filtering to remove short, low-confidence predictions, models spanning gaps, and potential transposable elements.

        - tRNA prediction: tRNA genes are predicted using [tRNAscan-SE](https://github.com/UCSC-LoweLab/tRNAscan-SE), ensuring comprehensive annotation of both protein-coding and non-coding genes.

        - NCBI submission preparation: Generation of an NCBI-compatible annotation table (.tbl format) and conversion to GenBank format using tbl2asn. The workflow also includes a validation step to parse NCBI error reports and alert users to any gene models that need manual correction.

05. **Gene Annotation:**

    A comprehensive gene annotation process assigns functional information to the identified genes. This process integrates multiple annotation tools and culminates in a final annotation round performed by [funannotate](https://github.com/nextgenusfs/funannotate). Below is an overview of the workflow:

    - [InterProScan](https://github.com/ebi-pf-team/interproscan) (v5.65-97.0): This tool is employed to assign protein domains and predict functional sites within the gene models. It integrates data from multiple databases such as `Pfam`, `SMART`, `PANTHER` and `PROSITE`, providing a rich set of functional annotations.

    - [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) (v2.1.12): This software is used to predict orthology and functional annotations based on the `EggNOG` database (v5.0). It helps in assigning Gene Ontology (GO) terms, enzyme codes, and pathway annotations to the gene models, offering insights into the biological roles of the proteins.
        
    - [antiSMASH](https://github.com/antismash/antismash) (v7.1): For fungal genomes, secondary metabolite gene clusters related to antibiotics or toxins are of particolar interest.

    - HMMer for [PFAM](http://pfam.xfam.org/) database (v36.0)

    - [DIAMOND](https://github.com/bbuchfink/diamond) for [UniProt](https://www.uniprot.org/) (v2024_01)

    - [DIAMOND](https://github.com/bbuchfink/diamond) for [MEROPS](https://www.ebi.ac.uk/merops/) database (v12.0)

    - CAZyme annotation with [dbCAN](https://bcb.unl.edu/dbCAN2/) (v12.0).

06. **Report:**
    * Results are parsed and aggregated to generate a report using [MultiQC](https://github.com/MultiQC/MultiQC).

## Installation
`FunFlux` automatically downloads all dependencies and several databases. However, some external databases require manual download before running the workflow.

1. **Download FunFlux:**

    Download via command line as:
    ```bash
    # Clone the directory
    git clone https://github.com/iLivius/FunFlux.git
    ```

2. **Install Snakemake:**

    `FunFlux` relies on [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to manage the workflow execution. Find the official and complete set of instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). To install Snakemake as a Conda environment:
    ```bash
    # Install Snakemake in a new Conda environment
    mamba create -c conda-forge -c bioconda -n snakemake snakemake
    ```

3. **Databases:**

    While `FunFlux` automates the installation of all software dependencies, some external databases need to be downloaded manually. If you have already installed these databases, you can skip this paragraph and proceed to the  [configuration](#configuration) section.
    
    Here are the required databases and software that need manual installation.

    * `NCBI core nt` database, adapted from [here](https://gist.github.com/ppflrs/336e49f8ae3843dc06cc3925940f3024):
        ```bash
        # Create a list of all core nt links in the directory designated to host the database (recommended)
        rsync --list-only rsync://ftp.ncbi.nlm.nih.gov/blast/db/core_nt.*.gz | grep '.tar.gz' | awk '{print "ftp.ncbi.nlm.nih.gov/blast/db/" $NF}' > nt_links.list
        
        # Alternatively, create a list of nt links for bacteria only 
        rsync --list-only rsync://ftp.ncbi.nlm.nih.gov/blast/db/nt_prok.*.gz | grep '.tar.gz' | awk '{print "ftp.ncbi.nlm.nih.gov/blast/db/" $NF}' > nt_prok_links.list
       
        # Download in parallel, without overdoing it
        cat nt*.list | parallel -j4 'rsync -h --progress rsync://{} .'

        # Decompress with multiple CPUs
        find . -name '*.gz' | parallel -j4 'echo {}; tar -zxf {}'

        # Get NCBI taxdump
        wget -c 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        tar -zxvf taxdump.tar.gz

        # Get NCBI BLAST taxonomy
        wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz'
        tar -zxvf taxdb.tar.gz

        # Get NCBI accession2taxid file
        wget -c 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'
        gunzip nucl_gb.accession2taxid.gz
        ```
        *NOTE: **Skip the download if you provide previously assembled contigs as input.** The complete NCBI core nt database and taxonomy-related files should take around 223 GB of hard drive space.*

    * `UNITE` database:

        - Visit the [UNITE USEARCH/UTAX release for eukaryotes](https://doi.plutof.ut.ee/doi/10.15156/BIO/2959341).
        - Download the `utax_reference_dataset_all_04.04.2024.fasta.gz` file and decompress it.
        - Add the PATH to the FASTA file in the `config.yaml` file. See the [configuration](#configuration) section.
   
    * `eggNOG diamond` database:
        ```bash
       # Create a Conda environment with eggnog-mapper, first
       conda create -n eggnog-mapper eggnog-mapper=2.1.12

       # Activate the environment
       conda activate eggnog-mapper

       # Create a directory where you want to install the diamond database for eggnog-mapper (example) 
       mkdir /data/eggnog_db

       # Finally, download the diamond db in the newly created directory 
       download_eggnog_data.py --data_dir /data/eggnog_db -y
        ```
        *NOTE: the eggNOG database requires ~50 GB of space.*

    * Download and set up `Genemark-ES/ET`:

        - Visit the `GeneMark` download page [here](http://topaz.gatech.edu/GeneMark/license_download.cgi).

        - Follow the instructions to download `GeneMark-ES/ET`.

        - Change the shebang line in Perl scripts, as follows:
            ```bash
            # After downloading, navigate to the GeneMark directory (example):
            cd /gmes_linux_64_4

            # Change the shebang line in all Perl scripts to use /usr/bin/env perl:
            find . -type f -name "*.pl" -exec sed -i '1s|^#!/usr/bin/perl|#!/usr/bin/env perl|' {} +

            # Test the software:
            ./gmes_petap.pl
            ```
    * Download and set up `InterProScan`:
        ```bash
        # Download this version although probably also more recent ones should work: 
        wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.65-97.0/interproscan-5.65-97.0-64-bit.tar.gz

        # Exctract the tarball:
        tar -pxvzf interproscan-5.65-97.0-64-bit.tar.gz

        # From inside the iprscan dir, index the hmm models:
        python3 setup.py -f interproscan.properties

        # Check the shell script, inside the iprscan dir:
        ./interproscan.sh
        ```

## Configuration
Before running `FunFlux`, you must edit the `config.yaml` file with a text editor. The file is organized in different sections: `links`, `directories`, `files`, `resources` and `parameters`, respectively.

If you have the FASTA files of previously assembled genomes as input, you must edit the `config_funnotator.yaml` file, instead.
 
- `links`

    This section should work fine as it is, therefore it is recommanded to change the `links` only if necessary or to update the database versions:

    - [phix_link](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015): Path to the PhiX genome reference used by Illumina for sequencing control. It is not needed in `config/config_funnotator.yaml`.

- `directories`

    Update paths based on your file system: 

    - **fastq_dir**: Directory containing the paired-end reads of your sequenced strains, in FASTQ format. You can provide as many as you like but at the following conditions:
        1. Files can only have the following extensions: `fastq`, `fq`, `fastq.gz`, or `fq.gz`.
        2. You can provide multiple samples but the extension should be the same for all files. So, don't mix files with different extensions.
        3. `No underscores` are allowed in sample names.
        4. Use `_R1` and `_R2` to define PE reads of each sample. i.e. `sample_R1.fastq.gz`, `sample_R2.fastq.gz`.

    - **input_dir** : Alternatively, if you want to analyze already available contigs, provide the directory path to the FASTA files, in the `config/config_funnotator.yaml`. Also in this case, you can provide as many genomes as you like but at the following conditions:
        1. Files can only have the following extensions: `fasta`, `fa`, or `fna`.
        2. You can provide multiple samples but the extension should be the same for all files. So, don't mix files with different extensions.
        3. `No underscores` are allowed in sample names. See example: `sample-1.fasta`, `sample2.fasta`.
                     
    - **out_dir**: This directory will store all output files generated by `FunFlux`. Additionally, by default, `FunFlux` will install required software and databases here, within Conda environments. Reusing this output directory for subsequent runs avoids reinstalling everything from scratch.

    - **blast_db**: Path to the whole `NCBI core nt`. See [installation](#installation). Not needed for `config/config_funnotator.yaml`, if you work with FASTA files of previouslsy assembled contigs. 

    - **eggnog_db**: Path to the diamond database for `eggNOG`. Download details in [installation](#installation).

    - **genemark_dir**: Path to `gmes_linux_64_4` directory. To get and configure `GeneMark-ES/ET`, see the [installation](#installation), above.

    - **funannotate_db**: Provide a path for [funannotate](https://github.com/nextgenusfs/funannotate) to automatically install the following databases:

        ```bash
        $ funannotate database

        Funannotate Databases currently installed:

        Database          Type        Version      Date         Num_Records   Md5checksum                     
        merops            diamond     12.0         2017-10-04          5009   a6dd76907896708f3ca5335f58560356
        uniprot           diamond     2024_01      2024-01-24        570830   c7507ea16b3c4807971c663994cad329
        dbCAN             hmmer3      12.0         2023-08-02           699   fb112af319a5001fbf547eac29e7c3b5
        pfam              hmmer3      36.0         2023-07            20795   0725495ccf049a4f198fcc0a92f7f38c
        repeats           diamond     1.0          2024-03-12         11950   4e8cafc3eea47ec7ba505bb1e3465d21
        go                text        2024-01-17   2024-01-17         47729   7e6b9974184dda306e6e07631f1783af
        mibig             diamond     1.4          2024-03-12         31023   118f2c11edde36c81bdea030a0228492
        interpro          xml         98.0         2024-01-25         40768   502ea05009761b893dedb56d5ea89c48
        busco_outgroups   outgroups   1.0          2024-03-12             8   6795b1d4545850a4226829c7ae8ef058
        gene2product      text        1.92         2023-10-02         34459   32a4a80987720e0872377de3207dc0f5
        ``` 

- `files`

    - **its_db**: Path to the [UNITE USEARCH/UTAX release for eukaryotes](https://doi.plutof.ut.ee/doi/10.15156/BIO/2959341) decompressed FASTA. Find the download instructions in the [installation](#installation) paragraph.

    - **annotation_params**: Path to a tab-delimited annotation parameter file, as displayed below. 
    
        An example is provided in the `config` directory as `annotation_parameters.tsv`:

        | #Sample        | Species                    | Proteins              | Model                |
        |:--------------:|:--------------------------:|:---------------------:|:--------------------:|
        | ARSEF3097      | Beauveria bassiana         | /path/to/proteins.faa | fusarium_graminearum |
        | 150-1          | Lecanicillium fungicola    | /path/to/proteins.faa | fusarium_graminearum |
        | FJII-L10-SW-P1 | Parengyodontium torokii    | /path/to/proteins.faa | fusarium_graminearum |
        | HWLR35         | Lecanicillium psalliotae   | /path/to/proteins.faa | fusarium_graminearum |
        | JC-1038        | Gamszarea kalimantanensis  | /path/to/proteins.faa | fusarium_graminearum |
        | MBC-099        | Lecanicillium aphanocladii | /path/to/proteins.faa | fusarium_graminearum |
        | MBC-350        | Akanthomyces uredinophilus | /path/to/proteins.faa | fusarium_graminearum |
        | MBC-401        | Cordyceps farinosa         | /path/to/proteins.faa | fusarium_graminearum |
        | MBC-695        | Akanthomyces uredinophilus | /path/to/proteins.faa | fusarium_graminearum |
        | MBC-701        | Akanthomyces dipterigenus  | /path/to/proteins.faa | fusarium_graminearum |

    - **iprscan**: Path to the [InterProScan](https://github.com/ebi-pf-team/interproscan) shell script. See the [installation](#installation) section for more details.
  
- `resources`
    
    In this section you can specify the hardware resources available to the workflow: 

  - threads: max number of CPUs used by each rule
  - ram_gb: max amount of RAM used by [SPAdes](https://github.com/ablab/spades). Not necessary in `config/config_funnotator.yaml`.

- `parameters`

    **Genus filtering**: `FunFlux` includes an optional parameter to specify the bacterial `genus` of contigs you wish to retain in the final assembly. If left blank, `FunFlux` will automatically keep contigs associated with the most abundant taxon, based on relative composition determined through `BLAST` analysis. While this approach generally works well, it has limitations, such as reduced resolution at the species level due to reliance on the cumulative best scores of `BLAST` hits. Additionally, this method may be problematic if the contaminant organism belongs to the same genus as your target organism, or if you are working with co-cultured closely related species or strains. If the `genus` parameter introduces more issues than benefits, simply remove the `genus` option from the `config.yaml` file. This feature is not available in `config/config_funnotator.yaml`.

    - **Using** the `genus` parameter: if a contaminant is ascertained to be more abundant than your target organism, you can re-run the workflow after reviewing the assembly [output](#output). Specify the `genus` of the desired bacterial taxon you want to keep in during the re-run. 
    
    - **Disabling** the `genus` filtering: if either the automatic inference of contaminant contigs or the manual selection of the desired taxon are still not working for you, simply delete the `genus` option from the `parameters`. In this case, only contigs tagged as "no-hit" after `BLAST` search will be filtered out.

## Running FunFlux
`FunFlux` can be executed as simply as a `Snakefile`. Please refer to the official [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html) for more details.
```bash
# First, activate the Snakemake Conda environment.
conda activate snakemake

# Navigate inside the FunFlux downloaded directory.

# Customize the "config.yaml" configuration file in the "config" sub-directory

# Launch the workflow
snakemake --sdm conda --cores 50 --jobs 2
```
*NOTE:  Starting from Snakemake version 8.4.7, the --use-conda option has been deprecated. Instead, you should now use --software-deployment-method conda or --sdm conda.*

**IMPORTANT**: If you need to analyze previously assembled fungal genomes, provided as FASTA files, use `Funnotator`, instead:
```bash
# Activate the Snakemake Conda environment
conda activate snakemake

# Navigate inside the FunFlux directory

# Customize the "config_funnotator.yaml" configuration file in the "config" sub-directory 

# Launch the workflow specifying the "Funannotator" Snakefile
snakemake --snakefile workflow/Funannotator --sdm conda --cores 50 --jobs 2
```

## Output
Here's a breakdown of the sub-directories created by `FunFlux` within the main output folder, along with explanations of their contents. Please notice that `Funnotator` will produce a similar, simplified output.
```
├── 01.pre-processing
├── 02.assembly
├── 03.post-processing
├── 04.annotation
├── logs
└── report
```
- `01.pre-processing`: QC and statistics of raw reads and trimmed reads, produced by [fastp](https://github.com/OpenGene/fastp) (v0.23.4).

- `02.assembly`: Content output by [SPAdes](https://github.com/ablab/spades) (v4.0.0). In addition to the raw contigs, you will also find the filtered contigs (>500bp and at least 2x) and the selected contigs, which are the contigs selected after BLAST search and decontamination (see `parameters` in the [configuration](#configuration) section above). The follow-up applications used during the worflow will either use selected contigs (i.e. for annotation purposes) or raw, filtered and selected contigs (i.e. to evaluate the genome completenness and contamination).

- `03.post-processing`: Contains the following sub-directories:
    - **mapping_evaluation**: [QualiMap](http://qualimap.conesalab.org/) (v2.3) output based on filtered contigs.
    - **contaminants**: Contig selection based on [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/) (v2.15.0) search and [BlobTools](https://github.com/DRL/blobtools) (1.1.1) analysis. Check the `composition` text file for a quick overview of the relative composition of your assembly.
    - **assembly_evaluation**: [Quast](https://github.com/ablab/quast) (v5.2.0) output based on selected contigs.
    - **completenness_evaluation**: [BUSCO](https://busco.ezlab.org/) (v5.5.0) output based on selected contigs.
     - **ITS_extraction**: [ITSx](https://microbiology.se/software/itsx/) (v1.1.3) output based on raw contigs and classified using the [SINTAX](https://www.drive5.com/sintax/) algorithm re-implemented in [VSEARCH](https://github.com/torognes/vsearch) (v2.28.1). It is recommendable to use the latest [UNITE](https://unite.ut.ee/repository.php) database, as reference.
 
- `04.annotation`: Contains the following sub-directories:
    - **iprscan**: Annotation output by [InterProScan](https://github.com/ebi-pf-team/interproscan) (v5.65-97.0), in XML format. 
    - **eggnog**: Functional annotation produced by [eggNOG](https://github.com/eggnogdb) mapper (v2.1.12).
    - **antismash**: Secondary metabolites inferred by [antiSMASH](https://github.com/antismash/antismash) (v7.1.0).
    - **funannotate**: Prediction and annotation directories output by [funannotate](https://github.com/nextgenusfs/funannotate) (v1.8.15).
        ```
        ├── annotate_misc
        ├── annotate_results
        ├── logfiles
        ├── predict_misc
        └── predict_results
        ```
   
- `report`: [MultiQC](https://github.com/MultiQC/MultiQC) (v1.23) is used to parse and aggregate the results of the following tools:
    1. [fastp](https://github.com/OpenGene/fastp) (v0.23.4)
    2. [QualiMap](http://qualimap.conesalab.org/) (v2.3)
    3. [Quast](https://github.com/ablab/quast) (v5.2.0)
    4. [BUSCO](https://busco.ezlab.org/) (v5.5.0)

## Acknowledgements
This work was supported by [BeXyl](https://cordis.europa.eu/project/id/101060593) (Beyond Xylella, Integrated Management Strategies for Mitigating *Xylella fastidiosa* impact in Europe) - HORIZON-CL6-2021-FARM2FORK-01-04, grant ID 101060593.

## Citation
In preparation.

## References
01. Bankevich, A., Nurk, S., Antipov, D., Gurevich, A. A., Dvorkin, M., Kulikov, A. S., Lesin, V. M., Nikolenko, S. I., Pham, S., Prjibelski, A. D., Pyshkin, A. V., Sirotkin, A. V., Vyahhi, N., Tesler, G., Alekseyev, M. A., & Pevzner, P. A. (2012). SPAdes: A New Genome Assembly Algorithm and Its Applications to Single-Cell Sequencing. Journal of Computational Biology, 19(5), 455–477. https://doi.org/10.1089/cmb.2012.0021

02. Bengtsson-Palme, J., Ryberg, M., Hartmann, M., Branco, S., Wang, Z., Godhe, A., De Wit, P., Sánchez-García, M., Ebersberger, I., de Sousa, F., Amend, A., Jumpponen, A., Unterseher, M., Kristiansson, E., Abarenkov, K., Bertrand, Y. J. K., Sanli, K., Eriksson, K. M., Vik, U., … Nilsson, R. H. (2013). Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. Methods in Ecology and Evolution, 4(10), 914–919. https://doi.org/10.1111/2041-210X.12073

03. Blin, K., Shaw, S., Augustijn, H. E., Reitz, Z. L., Biermann, F., Alanjary, M., Fetter, A., Terlouw, B. R., Metcalf, W. W., Helfrich, E. J. N., van Wezel, G. P., Medema, M. H., & Weber, T. (2023). antiSMASH 7.0: New and improved predictions for detection, regulation, chemical structures and visualisation. Nucleic Acids Research, 51(W1), W46–W50. https://doi.org/10.1093/nar/gkad344

04. Blum, M., Chang, H.-Y., Chuguransky, S., Grego, T., Kandasaamy, S., Mitchell, A., Nuka, G., Paysan-Lafosse, T., Qureshi, M., Raj, S., Richardson, L., Salazar, G. A., Williams, L., Bork, P., Bridge, A., Gough, J., Haft, D. H., Letunic, I., Marchler-Bauer, A., … Finn, R. D. (2021). The InterPro protein families and domains database: 20 years on. Nucleic Acids Research, 49(D1), D344–D354. https://doi.org/10.1093/nar/gkaa977

05. Borodovsky, M., & Lomsadze, A. (2011). Eukaryotic Gene Prediction Using GeneMark.hmm-E and GeneMark-ES. Current Protocols in Bioinformatics / Editoral Board, Andreas D. Baxevanis ... [et Al.], CHAPTER, Unit-4.610. https://doi.org/10.1002/0471250953.bi0406s35

06. Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature Methods, 12(1), 59–60. https://doi.org/10.1038/nmeth.3176

07. Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10, 421. https://doi.org/10.1186/1471-2105-10-421

08. Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Molecular Biology and Evolution, 38(12), 5825–5829. https://doi.org/10.1093/molbev/msab293

09. Challis, R., Richards, E., Rajan, J., Cochrane, G., & Blaxter, M. (2020). BlobToolKit – Interactive Quality Assessment of Genome Assemblies. G3 Genes|Genomes|Genetics, 10(4), 1361–1374. https://doi.org/10.1534/g3.119.400908

10. Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

11. Edgar, R. C. (2016). SINTAX: A simple non-Bayesian taxonomy classifier for 16S and ITS sequences (p. 074161). bioRxiv. https://doi.org/10.1101/074161

12. Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: Summarize analysis results for multiple tools and samples in a single report. Bioinformatics (Oxford, England), 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354

13. Finn, R. D., Bateman, A., Clements, J., Coggill, P., Eberhardt, R. Y., Eddy, S. R., Heger, A., Hetherington, K., Holm, L., Mistry, J., Sonnhammer, E. L. L., Tate, J., & Punta, M. (2014). Pfam: The protein families database. Nucleic Acids Research, 42(Database issue), D222–D230. https://doi.org/10.1093/nar/gkt1223

14. Frith, M. C. (2011). A new repeat-masking method enables specific detection of homologous sequences. Nucleic Acids Research, 39(4), e23. https://doi.org/10.1093/nar/gkq1212

15. Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: Quality assessment tool for genome assemblies. Bioinformatics (Oxford, England), 29(8), 1072–1075. https://doi.org/10.1093/bioinformatics/btt086

16. Haas, B. J., Salzberg, S. L., Zhu, W., Pertea, M., Allen, J. E., Orvis, J., White, O., Buell, C. R., & Wortman, J. R. (2008). Automated eukaryotic gene structure annotation using EVidenceModeler and the Program to Assemble Spliced Alignments. Genome Biology, 9(1), R7. https://doi.org/10.1186/gb-2008-9-1-r7

17. Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S. K., Cook, H., Mende, D. R., Letunic, I., Rattei, T., Jensen, L. J., von Mering, C., & Bork, P. (2019). eggNOG 5.0: A hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Research, 47(D1), D309–D314. https://doi.org/10.1093/nar/gky1085

18. Jonathan M. Palmer, & Jason Stajich. (2020). Funannotate v1.8.1: Eukaryotic genome annotation [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.4054262

19. Jones, P., Binns, D., Chang, H.-Y., Fraser, M., Li, W., McAnulla, C., McWilliam, H., Maslen, J., Mitchell, A., Nuka, G., Pesseat, S., Quinn, A. F., Sangrador-Vegas, A., Scheremetjew, M., Yong, S.-Y., Lopez, R., & Hunter, S. (2014). InterProScan 5: Genome-scale protein function classification. Bioinformatics, 30(9), 1236–1240. https://doi.org/10.1093/bioinformatics/btu031

20. KorfLab/SNAP. (2024). [C]. The Korf Lab. https://github.com/KorfLab/SNAP (Original work published 2017)

21. Köster, J., & Rahmann, S. (2012). Snakemake—A scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520–2522. https://doi.org/10.1093/bioinformatics/bts480

22. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), Article 4. https://doi.org/10.1038/nmeth.1923

23. Letunic, I., Khedkar, S., & Bork, P. (2021). SMART: Recent updates, new developments and status in 2020. Nucleic Acids Research, 49(D1), D458–D460. https://doi.org/10.1093/nar/gkaa937

24. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

25. Lowe, T. M., & Eddy, S. R. (1997). tRNAscan-SE: A program for improved detection of transfer RNA genes in genomic sequence. Nucleic Acids Research, 25(5), 955–964.

26. Majoros, W. H., Pertea, M., & Salzberg, S. L. (2004). TigrScan and GlimmerHMM: Two open source ab initio eukaryotic gene-finders. Bioinformatics (Oxford, England), 20(16), 2878–2879. https://doi.org/10.1093/bioinformatics/bth315

27. Nilsson, R. H., Larsson, K.-H., Taylor, A. F. S., Bengtsson-Palme, J., Jeppesen, T. S., Schigel, D., Kennedy, P., Picard, K., Glöckner, F. O., Tedersoo, L., Saar, I., Kõljalg, U., & Abarenkov, K. (2019). The UNITE database for molecular identification of fungi: Handling dark taxa and parallel taxonomic classifications. Nucleic Acids Research, 47(D1), D259–D264. https://doi.org/10.1093/nar/gky1022

28. Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: Advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics, 32(2), 292–294. https://doi.org/10.1093/bioinformatics/btv566

29. Rawlings, N. D., Waller, M., Barrett, A. J., & Bateman, A. (2014). MEROPS: The database of proteolytic enzymes, their substrates and inhibitors. Nucleic Acids Research, 42(D1), D503–D509. https://doi.org/10.1093/nar/gkt953

30. Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: A versatile open source tool for metagenomics. PeerJ, 4, e2584. https://doi.org/10.7717/peerj.2584

31. Sigrist, C. J. A., de Castro, E., Cerutti, L., Cuche, B. A., Hulo, N., Bridge, A., Bougueleret, L., & Xenarios, I. (2013). New and continuing developments at PROSITE. Nucleic Acids Research, 41(Database issue), D344-347. https://doi.org/10.1093/nar/gks1067

32. Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015). BUSCO: Assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31(19), 3210–3212. https://doi.org/10.1093/bioinformatics/btv351

33. Slater, G. S. C., & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. BMC Bioinformatics, 6(1), 31. https://doi.org/10.1186/1471-2105-6-31

34. Stanke, M., Keller, O., Gunduz, I., Hayes, A., Waack, S., & Morgenstern, B. (2006). AUGUSTUS: Ab initio prediction of alternative transcripts. Nucleic Acids Research, 34(Web Server issue), W435–W439. https://doi.org/10.1093/nar/gkl200

35. The UniProt Consortium. (2023). UniProt: The Universal Protein Knowledgebase in 2023. Nucleic Acids Research, 51(D1), D523–D531. https://doi.org/10.1093/nar/gkac1052

36. Thomas, P. D., Ebert, D., Muruganujan, A., Mushayahama, T., Albou, L.-P., & Mi, H. (2022). PANTHER: Making genome-scale phylogenetics accessible to all. Protein Science, 31(1), 8–22. https://doi.org/10.1002/pro.4218

37. Zheng, J., Ge, Q., Yan, Y., Zhang, X., Huang, L., & Yin, Y. (2023). dbCAN3: Automated carbohydrate-active enzyme and substrate annotation. Nucleic Acids Research, 51(W1), W115–W121. https://doi.org/10.1093/nar/gkad328
>>>>>>> Initial commit of FunFlux project.
