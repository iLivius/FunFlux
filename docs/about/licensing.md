# Licensing

FunFlux is released under the [MIT License](https://github.com/iLivius/FunFlux/blob/main/LICENSE).
The tools and databases it uses carry their own licences, and some of them restrict commercial
use. This page lists those licences as they were read on the dates given, with the source of each.

The repository contains code and download addresses, and no third-party data. Every database is
either one you downloaded yourself or one the workflow downloads on your machine, from an address
in the `links` block of the configuration or with the tool's own download command.

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **AMR screen** | the optional resistance-mutation screen, off by default (`parameters.amr.chroquetas`), in which ChroQueTas looks in each genome for the changes FungAMR lists for the species |
    | **FungAMR** | the database of amino-acid changes studied for resistance to antifungal drugs and fungicides, with evidence levels for and against |
    | **member database** | one of the separate signature databases, such as Pfam, PROSITE and SMART, that InterProScan searches; each keeps its own licence terms |
    | **CC BY-NC-ND 4.0** | a Creative Commons licence: credit the authors, non-commercial use only, and no modified copies may be shared; ChroQueTas, FungAMR and PROSITE carry it |
    | **share-alike** | a licence condition that copies you change or pass on keep the same licence; it does not apply to running the program |
    | **the screen** | the contamination and mitochondrion screen of the assembly: BLAST against `core_nt`, BlobTools, the mitochondrion screen and the contig selector; in `contigs` mode only with `parameters.contigs.screen: true` |
    | **funannotate database** | the frozen snapshot of funannotate's data, such as Swiss-Prot, Pfam, the BUSCO set and stored AUGUSTUS models, downloaded once from Zenodo into `directories.funannotate_db` |

For each component that restricts commercial use, the table below gives the restriction, when
FunFlux uses the component and how to leave it out. GeneMark-ES, SignalP 6, Medaka and the AMR
screen can each be left out with one setting, but PROSITE and SMART are searched in every run.

| Component | Restriction | When FunFlux uses it | How to leave it out |
|---|---|---|---|
| GeneMark-ES | royalty-free only for academic and non-profit institutions and US government agencies; the download goes through a registration form | optional: only when `directories.genemark_dir` is set | leave `directories.genemark_dir` empty, or set `parameters.genemark.mode: off` |
| SignalP 6 | downloaded by academic users only, through a registration form | optional: only when `parameters.signalp6.mode` is not `off` | `parameters.signalp6.mode: off` |
| Medaka | research purposes only | polishes the assemblies of ONT samples in `longread` and `hybrid` mode | `medaka_model: false` in `parameters.longread` or `parameters.hybrid` |
| ChroQueTas and its FungAMR database | CC BY-NC-ND 4.0: non-commercial use only, and no modified copies may be passed on | the AMR screen, off by default; `parameters.amr.chroquetas: on` is your acceptance of the licence, and FunFlux then downloads both, unchanged, on your machine ([Antifungal resistance mutations](../analysis/amr.md#licence)) | keep the AMR screen off |
| PROSITE and SMART, two of the member databases InterProScan searches | commercial users of PROSITE pay a yearly licence fee; SMART is free for academic users only | every InterProScan run: FunFlux starts InterProScan without `-appl`, so it runs every analysis the installation provides, these two included | no configuration key leaves them out |

Of the terms that could be checked, those in the table above are the only ones that limit who may
use a component. No terms could be found for dbCAN, the eggNOG 5.0.2 data, the Tiberius `fungi`
model and the tidk telomere repeat database, and MIBiG's licence
comes only from a registry entry. The terms of the databases antiSMASH downloads, and of the member
databases InterProScan searches other than Pfam, PROSITE and SMART, were not checked.

The next table gives the licences of the tools and databases FunFlux uses. The Conditions column
lists restrictions on who may use a component and what its licence asks of you: registration,
citation, share-alike. A dash means none of these.

!!! note "Running FunFlux places no share-alike duty on you"

    Share-alike applies when you modify or redistribute a program, not when you run it. Under the
    AGPL and the Open Software License (OSL), it also applies when others use a modified copy over
    a network.

| Tool or database | Licence | Conditions | Source | Checked |
|---|---|---|---|---|
| Snakemake 9.14.6 | MIT | — | [GitHub](https://github.com/snakemake/snakemake) | 2026-10-08 |
| Bowtie 2 version 2.5.4 (`illumina`, `hybrid`) | GPL-3.0 | Share-alike | [GitHub](https://github.com/BenLangmead/bowtie2) | 2026-10-08 |
| fastp 1.0.1 (`illumina`, `hybrid`) | MIT | — | [GitHub](https://github.com/OpenGene/fastp) | 2026-10-08 |
| SPAdes 4.2.0 (`illumina`, `hybrid`) | GPL-2.0 | Share-alike | [GitHub](https://github.com/ablab/spades) | 2026-10-08 |
| Qualimap 2.3 (`illumina`, `longread`, `hybrid`) | GPL-2.0 | Share-alike | [LICENSE in the Qualimap repository](https://bitbucket.org/kokonech/qualimap/raw/master/LICENSE) | 2026-10-08 |
| fastplong 0.7.1 (ONT samples in `longread`, `hybrid`) | MIT | — | [GitHub](https://github.com/OpenGene/fastplong) | 2026-11-09 |
| Filtlong 0.3.1 (`longread`, `hybrid`) | GPL-3.0 | Share-alike | [GitHub](https://github.com/rrwick/Filtlong) | 2026-10-08 |
| NanoPlot 1.46.2 (`longread`, `hybrid`) | MIT | — | [GitHub](https://github.com/wdecoster/NanoPlot) | 2026-10-08 |
| hifiasm 0.25.0 (`longread`, `hybrid`) | MIT | — | [GitHub](https://github.com/chhylp123/hifiasm) | 2026-10-08 |
| Flye 2.9.6 (`longread`, `hybrid`) | BSD-3-Clause | — | [GitHub](https://github.com/mikolmogorov/Flye) | 2026-10-08 |
| minimap2 2.30 (`longread`, `hybrid`, `contigs` with the screen); 2.31 inside `funannotate clean`, every mode, and inside Medaka | MIT | — | [LICENSE.txt](https://raw.githubusercontent.com/lh3/minimap2/master/LICENSE.txt) | 2026-10-08 |
| SAMtools 1.21 (`illumina`, `longread`, `hybrid`, `contigs` with the screen) | MIT/Expat | — | [GitHub](https://github.com/samtools/samtools) | 2026-10-08 |
| Medaka 2.2.2 (ONT samples in `longread`, `hybrid`) | Oxford Nanopore Technologies PLC. Public License 1.0 | Use solely for research purposes, defined as internal research not directed towards commercial advantage or monetary compensation; sponsored and grant-funded research are allowed | [GitHub](https://github.com/nanoporetech/medaka); `LICENSE` of the installed package | 2026-10-08 |
| Polypolish 0.6.1 and BWA 0.7.19 (`hybrid`) | GPL-3.0 | Share-alike | [Polypolish](https://github.com/rrwick/Polypolish), [BWA](https://github.com/lh3/bwa) | 2026-10-08 |
| BLAST+ 2.16.0 | Public domain (United States Government Work) | "Please cite the author in any work or product based on this material." | public-domain notice shipped with BLAST+ (in `update_blastdb.pl`) | 2026-10-08 |
| NCBI `core_nt`, the taxonomy files `taxdb` and `taxdump`, the PhiX genome | NCBI data policy | NCBI places no restrictions on use or distribution; some submitters of the original data may claim patent, copyright or other rights | [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) | 2026-10-08 |
| BlobTools 1.1.1 | GPL-3.0 | Share-alike | [GitHub](https://github.com/DRL/blobtools) | 2026-10-08 |
| Oxford Nanopore adapter and barcode sequences (shipped in `workflow/data/adapters/ont_chtd.fasta`) | no licence stated | Reproduced as Oxford Nanopore publishes them, with the source named in every record | [Chemistry Technical Document](https://nanoporetech.com/document/chemistry-technical-document) | 2026-11-09 |
| NCBI UniVec: the PacBio and MGI records (shipped in `workflow/data/adapters/univec_pacbio_mgi.fasta`) and the Illumina records (downloaded by rule `univec_illumina`) | NCBI data policy (United States Government Work) | Illumina's oligonucleotide sequences are copyright Illumina, Inc.; FunFlux downloads them on your machine and does not redistribute them | [UniVec](https://ftp.ncbi.nlm.nih.gov/pub/UniVec/) | 2026-11-09 |
| QUAST 5.3.0 | GPL-2.0 | Share-alike | [GitHub](https://github.com/ablab/quast) | 2026-10-08 |
| BUSCO 6.0.0 | MIT | — | [busco.ezlab.org](https://busco.ezlab.org/) | 2026-10-08 |
| BUSCO odb12 lineage datasets (BUSCO downloads them) | CC BY-ND 4.0 | Cite the odb12 paper ([Citation and references](citation.md)); modified datasets may not be shared | [busco.ezlab.org](https://busco.ezlab.org/) | 2026-10-08 |
| tidk 0.3.0 | MIT | — | [GitHub](https://github.com/tolkit/telomeric-identifier) | 2026-10-10 |
| tidk telomere repeat database (`telo_database_v2.tsv`, downloaded from the address in `links.telomere_db_link`) | not stated | No licence file in the repository; FunFlux downloads it on your machine and does not redistribute it | [GitHub](https://github.com/tolkit/a-telomeric-repeat-database) | 2026-10-10 |
| HMMER 3.4 (`nhmmer`) | BSD-3-Clause | — | [GitHub](https://github.com/EddyRivasLab/hmmer) | 2026-10-08 |
| ITSx 1.1.3 | GPL-3.0-or-later | Share-alike | licence header of the `ITSx` script (`ITSx --license`) | 2026-10-08 |
| VSEARCH 2.30.0 | GPL-3.0 or BSD-2-Clause, at the user's choice | Share-alike only under GPL-3.0 | [GitHub](https://github.com/torognes/vsearch) | 2026-10-08 |
| UNITE release for eukaryotes 19.02.2025 | CC BY 4.0 | Cite the release | [DOI record](https://doi.org/10.15156/BIO/3301246) | 2026-10-08 |
| tantan 51 | GPL-3.0 | Share-alike | [GitLab](https://gitlab.com/mcfrith/tantan) | 2026-10-08 |
| RepeatModeler 2.0.8 | Open Software License 2.1 | Share-alike | [GitHub](https://github.com/Dfam-consortium/RepeatModeler) | 2026-10-08 |
| RepeatMasker 4.2.3 | Open Software License 2.1 | Share-alike | `LICENSE` of the installed package; [GitHub](https://github.com/rmhubley/RepeatMasker) | 2026-10-08 |
| funannotate 1.8.17 | BSD-2-Clause | — | [GitHub](https://github.com/nextgenusfs/funannotate) | 2026-10-08 |
| AUGUSTUS 3.5.0 | Artistic License | — | [GitHub](https://github.com/Gaius-Augustus/Augustus) | 2026-10-08 |
| SNAP 2017_03_01 | MIT | — | [GitHub](https://github.com/KorfLab/SNAP) | 2026-10-08 |
| GlimmerHMM 3.0.4 | "OSI Certified Open Source Software"; the bioconda package names the Artistic License | — | [GlimmerHMM](https://ccb.jhu.edu/software/glimmerhmm/) | 2026-10-08 |
| EVidenceModeler 1.1.1 | BSD-3-Clause | — | [GitHub](https://github.com/EVidenceModeler/EVidenceModeler) | 2026-10-08 |
| Tiberius 1.1.8 (on by default; `parameters.tiberius.mode: off` leaves it out) | MIT | — | [GitHub](https://github.com/Gaius-Augustus/Tiberius) | 2026-10-09 |
| Tiberius fungi model, downloaded from the Tiberius authors' server | not stated | No licence is published with the file; FunFlux downloads it on your machine and does not redistribute it | [model file](https://bioinf.uni-greifswald.de/bioinf/tiberius/models/fungi.tar.gz) | 2026-10-09 |
| tRNAscan-SE 2.0.13 | GPL-3.0 | Share-alike | [GitHub](https://github.com/UCSC-LoweLab/tRNAscan-SE) | 2026-10-08 |
| Exonerate 2.4.0 | GPL-3.0 | Share-alike | [maintained fork of 2.4.0](https://github.com/nathanweeks/exonerate) | 2026-10-08 |
| DIAMOND 2.2.4 (in funannotate and antiSMASH) and 2.0.15 (in eggNOG-mapper) | GPL-3.0; the bioconda package of 2.0.15 names AGPL-3.0 | Share-alike | [GitHub](https://github.com/bbuchfink/diamond) | 2026-10-08 |
| GeneMark-ES 4 (optional) | GeneMark licence (Georgia Tech) | Royalty-free for academic and non-profit institutions and US government agencies; everyone else must arrange a licence. Registration form (name, institution, e-mail). The key in `~/.gm_key` is valid for one year from download | [download page](https://genemark.bme.gatech.edu/license_download.cgi) | 2026-10-08 |
| funannotate database snapshot `funannotate_db_1.8.17_2026-10-08` (Zenodo record 23265189) | none on the record | The databases in it keep their own terms (the next eight rows) | [Zenodo](https://zenodo.org/records/23265189) | 2026-10-12 |
| UniProtKB/Swiss-Prot 2026_03 | CC BY 4.0 | Cite UniProt | [UniProt LICENSE](https://ftp.uniprot.org/pub/databases/uniprot/LICENSE) | 2026-10-08 |
| Pfam 38.2, InterPro 110.0 | CC0 1.0 | — | [InterPro licence](https://interpro-documentation.readthedocs.io/en/latest/license.html) | 2026-10-08 |
| Gene Ontology (`go.obo` of 2026-07-26) | CC BY 4.0 | Cite the Gene Ontology | [GO citation policy](https://geneontology.org/docs/go-citation-policy/) | 2026-10-08 |
| MEROPS 12.5 | GNU Library General Public License | Share-alike | [MEROPS availability](https://www.ebi.ac.uk/merops/about/availability.shtml) | 2026-10-08 |
| dbCAN (HMMdb V11) | none stated | The site states no terms for the downloads; its only copyright notice is "Copyright 2012-now © YIN LAB, UNL. All rights reserved." | [dbCAN](https://pro.unl.edu/dbCAN2) | 2026-10-08 |
| MIBiG 1.4 | CC BY 4.0 according to the Bioregistry entry | Not confirmed by MIBiG: its data repository on GitHub shows no licence, and its website could not be read | [Bioregistry](https://bioregistry.io/mibig) | 2026-10-08 |
| BUSCO `dikarya_odb9` set | CC BY-ND 4.0, from the BUSCO site's licence statement, which names no dataset version | The site names a paper to cite only for odb10 and odb12 datasets; modified datasets may not be shared | [busco.ezlab.org](https://busco.ezlab.org/) | 2026-10-08 |
| Other contents of the funannotate database: BUSCO sets of 6 outgroup species (the database index lists 8), 166 stored AUGUSTUS models, the gene2product list 1.99, a repeat-protein library | not checked | — | [Zenodo](https://zenodo.org/records/23265189) | — |
| InterProScan 5.78-109.0 | Apache-2.0 | Phobius, SignalP and TMHMM ship switched off because they are licensed separately; the member databases it searches keep their own terms (PROSITE and SMART: the next two rows) | [GitHub](https://github.com/ebi-pf-team/interproscan); [licensed analyses](https://interproscan-docs.readthedocs.io/en/v5/ActivatingLicensedAnalyses.html) | 2026-10-08 |
| PROSITE 2025_01 (searched by InterProScan) | CC BY-NC-ND 4.0 | Commercial users pay a yearly licence fee; academic, non-commercial use is free | [PROSITE licence](https://prosite.expasy.org/prosite_license.html) | 2026-10-08 |
| SMART 9.0 (searched by InterProScan) | SMART licence | Free for academic users; non-academic users contact biobyte solutions GmbH | [About SMART](https://smart.embl.de/help/smart_about.shtml) | 2026-10-08 |
| eggNOG-mapper 2.1.13 | AGPL-3.0 | Share-alike | [GitHub](https://github.com/eggnogdb/eggnog-mapper) | 2026-10-08 |
| eggNOG 5.0.2 data | not established | eggnog-mapper's README and wiki state no terms for the data, and the eggNOG 5 website could not be read | [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) | 2026-10-08 |
| antiSMASH 8.0.4 | AGPL-3.0-or-later | Share-alike | [GitHub](https://github.com/antismash/antismash) | 2026-10-08 |
| SignalP 6.0 (optional) | DTU Health Tech academic licence | Registration form with an e-mail address at an academic institution; personal and internal research use at one site; no commercial use, no modification, no redistribution. Other users contact DTU | [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) and its download form | 2026-10-08 |
| ChroQueTas 1.0.1 and its FungAMR database (AMR screen, off by default; downloaded from the address in `links.chroquetas_link`) | CC BY-NC-ND 4.0 | Non-commercial use only; credit the authors; modified copies may not be shared. FunFlux downloads the files unchanged on your machine and does not redistribute them | `LICENSE` in the [ChroQueTas repository](https://github.com/nmquijada/ChroQueTas) at the pinned commit | 2026-10-11 |
| miniprot 0.14 (AMR screen); 0.18 inside BUSCO, every mode | MIT | — | [GitHub](https://github.com/lh3/miniprot) | 2026-10-11 |
| MAFFT 7.525 (AMR screen) | BSD-3-Clause | — | [MAFFT licence](https://mafft.cbrc.jp/alignment/software/license.txt) | 2026-10-11 |
| SeqKit 2.1.0 (AMR screen) | MIT | — | [GitHub](https://github.com/shenwei356/seqkit) | 2026-10-11 |
| EMBOSS 6.6.0 (`seqret`, AMR screen) | GPL for the applications, LGPL for the libraries | Share-alike | [EMBOSS licence](http://emboss.sourceforge.net/licence/) | 2026-10-11 |
| MultiQC 1.33 | GPL-3.0-or-later | Share-alike | [GitHub](https://github.com/MultiQC/MultiQC) | 2026-10-08 |

!!! warning "Check `00.resources/` before you share or publish an output folder"

    `00.resources/` holds third-party files, each under its own terms.

    | Files | When they are there | Terms |
    |---|---|---|
    | `00.resources/signalp6/package/`: the SignalP 6 program and its model weights (1.6 GB with the fast model) | after the first SignalP 6 run on the output folder, which unpacks them from your archive | the SignalP licence forbids redistribution: remove the folder |
    | the UNITE reference | when the workflow downloaded it | CC BY 4.0: credit the release |
    | the antiSMASH databases | when the workflow downloaded them | not checked |
    | the telomere repeat database of tidk | when the workflow downloaded it | no licence stated |
    | `00.resources/tiberius/model/`: the Tiberius `fungi` model | while Tiberius is on, the default | no licence stated |
    | `00.resources/busco/`: the BUSCO lineage dataset | when `parameters.busco.lineage` names a lineage | CC BY-ND 4.0: modified datasets may not be shared |
    | `00.resources/chroquetas/`: ChroQueTas and its FungAMR database | with the AMR screen on | CC BY-NC-ND 4.0: pass them on only unchanged, for non-commercial purposes and with credit |

[Licensed and external tools](../getting-started/external-tools.md) explains how to install
GeneMark-ES, SignalP 6 and InterProScan, and what the workflow checks before it uses them. What to
cite for each tool and database is listed on [Citation and references](citation.md).
