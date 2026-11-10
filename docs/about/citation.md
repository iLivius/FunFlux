# Citation and references

Use this page when you publish results from FunFlux. It gives the reference for FunFlux itself
and for the v1 workflows. It then lists, stage by stage, the tools and databases a run can use,
with their versions and the numbers of their references. The page ends with the acknowledgements,
the full list of references and the works cited on [Validation on public genomes](validation.md).

!!! abstract "Terms used on this page"

    | | |
    |---|---|
    | **concept DOI** | the Zenodo DOI shared by all versions of a record, which always leads to the newest; each version also has a DOI of its own |
    | **contig filter** | rule `filter_contigs`: on SPAdes-style headers it drops contigs under 500 bp or 2x coverage; other headers are only cut at their first space or tab |
    | **mitochondrion screen** | rule `screen_mitochondria`, part of the contamination and mitochondrion screen: it confirms mitochondrial contigs from BLAST hits to mitochondrial genomes and, in the long-read modes, the assembler's evidence |
    | **contig selector** | rule `select_contigs`, the last step of the screen: it keeps or removes each contig by genus, as `parameters.decontamination` says; the default, `off`, keeps every contig |
    | **ITS window** | the stretch cut around one 5.8S gene, 2,000 bp on each side by default, turned to read from SSU to LSU before ITSx runs |
    | **contig summary** | the table `02.assembly/{sample}/eval/contig_summary.tsv`: one row per contig of the draft, with its evidence (telomeres, depth, rDNA units, genes) and a verdict |
    | **AMR screen** | the optional resistance-mutation screen, off by default (`parameters.amr.chroquetas`), in which ChroQueTas looks in each genome for the changes FungAMR lists for the species |
    | **funannotate database** | the frozen snapshot of funannotate's data, such as Swiss-Prot, Pfam, the BUSCO set and stored AUGUSTUS models, downloaded once from Zenodo into `directories.funannotate_db` |

## Citing FunFlux

> Antonielli, L., & Brader, G. (2026). *FunFlux: Integrated workflow for fungal genome assembly
> and annotation* (v2.0.0). Zenodo. <https://doi.org/10.5281/zenodo.13612159>

The DOI 10.5281/zenodo.13612159 is the concept DOI, which always resolves to the newest version.
Each release also has a version DOI of its own on Zenodo. Cite the version DOI when the exact code
matters.

The repository is <https://github.com/iLivius/FunFlux>, and the documentation is at
<https://iLivius.github.io/FunFlux/>. FunFlux is released under the MIT License. See
[Licensing](licensing.md).

The contig filter, the mitochondrion screen, the contig selector and the ITS windows have no tool
of their own to cite, so cite FunFlux for them. The contig filter of `contigs` mode, the contig
selector and the layout of the workflow come from
[BacFlux](https://doi.org/10.5281/zenodo.11143917), which is also released under the MIT License.

For the contig summary, tidk does the telomere search, so cite tidk with FunFlux. The papers behind
FunFlux's table of telomere repeat units are listed in the Contig summary rows under
[Tools and databases](#tools-and-databases).

!!! note "The ChroQueTas authors ask to be cited through the FungAMR paper"

    The AMR screen runs ChroQueTas. As the ChroQueTas authors ask, cite FungAMR and ChroQueTas
    through the FungAMR paper (Bédard et al. 2025). Cite the tools ChroQueTas runs as well, which
    are listed in the AMR screen rows under [Tools and databases](#tools-and-databases).

### Results from a v1 workflow

If a v1 workflow produced your results, cite the record of the version you ran:

| Workflow | Version (published) | DOI |
|---|---|---|
| FunFlux | 1.0.3 (2024-08-30) | [10.5281/zenodo.13612160](https://doi.org/10.5281/zenodo.13612160) |
| FunFlux | 1.0.4 (2025-06-08) | [10.5281/zenodo.15618338](https://doi.org/10.5281/zenodo.15618338) |
| FunFlux | 1.0.5 (2026-01-24) | [10.5281/zenodo.18363368](https://doi.org/10.5281/zenodo.18363368) |
| FunFlux | 1.0.6 (2026-01-31) | [10.5281/zenodo.18445164](https://doi.org/10.5281/zenodo.18445164) |
| FunFlux | 1.0.7 (2026-02-02) | [10.5281/zenodo.18455939](https://doi.org/10.5281/zenodo.18455939) |
| FunFlux | 1.1.0 (2026-06-18) | [10.5281/zenodo.20749375](https://doi.org/10.5281/zenodo.20749375) |
| Funnotator | shipped inside FunFlux 1.x | the DOI of the FunFlux release it came with |
| FunFluxL | 0.1.0 (2026-06-22) | [10.5281/zenodo.20796645](https://doi.org/10.5281/zenodo.20796645); concept DOI [10.5281/zenodo.20796644](https://doi.org/10.5281/zenodo.20796644) |
| FunFluxL+ | never released | none |

## Tools and databases

Cite the tools and databases your run used, as well as FunFlux. Which rows apply depends on the
mode, the read type and the settings of your run. The versions are those FunFlux 2.0.0 installs
or, for the tools you install yourself, those it was tested with. The numbers in the last column
point to the list under [References](#references).

| Stage | Tool or database | References |
|---|---|---|
| Workflow | Snakemake (9.14.6) | 54 |
| Reads and assembly, `illumina` | Bowtie 2 (2.5.4) for PhiX removal and read mapping, with SAMtools (1.21) | 41, 42, 24 |
| | fastp (1.0.1) | 18, 19, 20 |
| | SPAdes (4.2.0) | 60 |
| | Qualimap (2.3) | 56 |
| Reads and assembly, `longread` | fastplong (0.7.1), ONT adapter trimming; its author asks users to cite fastp until fastplong has a paper of its own | 18 |
| | Filtlong (0.3.1) | 77 |
| | NanoPlot (1.46.2) | 25 |
| | hifiasm (0.25.0), the default assembler | 21, 22, 23 |
| | Flye (2.9.6), the alternative assembler | 38 |
| | minimap2 (2.30) and SAMtools (1.21) for read mapping, Qualimap (2.3) | 45, 46, 24, 56 |
| | Medaka (2.2.2), ONT samples only | 57 |
| Reads and assembly, `hybrid` | the `illumina` and `longread` tools above | as above |
| | Polypolish (0.6.1) with BWA (0.7.19) | 78, 9, 44 |
| Contig self-mapping, `contigs` with the screen | minimap2 (2.30) and SAMtools (1.21) | 45, 46, 24 |
| Contamination screen | BLAST+ (2.16.0), megablast against NCBI core_nt | 14, 55, 80 |
| | BlobTools (1.1.1) | 40 |
| Adapter screen | BLAST+ (2.16.0), blastn-short, against the adapters and barcodes of the Oxford Nanopore Chemistry Technical Document (CHTD_500_v1_revAV_24July2026), the PacBio and MGI records of NCBI UniVec (the file of 21 August 2024) and the Illumina records of UniVec build 10.1 | 14 |
| Assembly QC | QUAST (5.3.0) | 53 |
| | BUSCO (6.0.0) with the odb12 lineage datasets, finding the genes with miniprot (0.18) | 72, 47 |
| Contig summary | tidk (0.3.0): `tidk explore` lists the tandem repeat units of the contig ends, `tidk ends` calls the arrays of every candidate unit | 12 |
| | the tidk telomere repeat database (`telo_database_v2.tsv`, downloaded from the address in `links.telomere_db_link`) | 12 |
| | FunFlux's own scripts; the telomere repeat units of the fungal groups | 49, 71, 16, 52, 65, 11, 27, 43, 10 |
| | the two databases that list a species' repeat: TeloBase and the Telomerase Database | 50, 59 |
| ITS taxonomy | nhmmer from HMMER (3.4), with the 5.8S profiles of ITSx | 76 |
| | ITSx (1.1.3) | 6 |
| | SINTAX in VSEARCH (2.30.0) | 26, 64 |
| | UNITE release for eukaryotes 19.02.2025 | 2, 1 |
| Repeat masking | tantan (51), the default | 30 |
| | RepeatModeler (2.0.8) and RepeatMasker (4.2.3) | 29, 68 |
| Gene prediction | funannotate (1.8.17): `clean`, `sort`, `predict`; `clean` finds duplicated contigs with minimap2 (2.31) | 58, 45, 46 |
| | AUGUSTUS (3.5.0) | 69, 70 |
| | SNAP (2017_03_01) | 39 |
| | GlimmerHMM (3.0.4) | 51 |
| | Tiberius (1.1.8), on by default, with its authors' fungi model | 31 |
| | BUSCO `dikarya_odb9` set of the funannotate database, on which SNAP and GlimmerHMM are trained | none: the BUSCO site names a paper only for odb10 and odb12 datasets |
| | EVidenceModeler (1.1.1) | 34 |
| | tRNAscan-SE (2.0.13) | 17 |
| | DIAMOND (2.2.4) and Exonerate (2.4.0), aligning your protein FASTA and UniProtKB/Swiss-Prot to the genome | 13, 67, 75 |
| | GeneMark-ES (optional) | 48, 73 |
| Functional annotation | InterProScan (5.78-109.0), with InterPro and Pfam | 36, 8 |
| | eggNOG-mapper (2.1.13) with DIAMOND (2.0.15), against eggNOG 5.0.2 | 15, 13, 35 |
| | antiSMASH (8.0.4) | 7 |
| | SignalP 6.0 (optional) | 74 |
| | funannotate `annotate` (1.8.17) with the databases held in the funannotate database: UniProtKB/Swiss-Prot 2026_03, Pfam 38.2, MEROPS 12.5, dbCAN (HMMdb V11), Gene Ontology, MIBiG 1.4 | 58, 75, 8, 62, 81, 3, 32, 79 |
| AMR screen (optional, off by default) | ChroQueTas (1.0.1) with its FungAMR database (`db-v20250811`) | 5, 61 |
| | miniprot (0.14), BLAST+ (2.16.0) `blastp`, MAFFT (7.525), SeqKit (2.1.0) and EMBOSS (6.6.0) `seqret`: the tools ChroQueTas runs | 47, 14, 37, 66, 63 |
| | the intron after codon 143 of cytochrome b, in the worked example on B05.10 | 4, 33 |
| Report | MultiQC (1.33) | 28 |

!!! note "The DOI of the funannotate database fixes the versions of the databases in it"

    The funannotate database, `funannotate_db_1.8.17_2026-10-08`, has the DOI
    [10.5281/zenodo.23265189](https://doi.org/10.5281/zenodo.23265189), which fixes the versions
    of the databases in it. It holds MIBiG 1.4 data, while reference 79 is the current MIBiG
    paper. See [Reference databases](../getting-started/databases.md#funannotate-db).

## Acknowledgements

This work was supported by the [RATION project](https://cordis.europa.eu/project/id/101084163)
(Risk AssessmenT InnOvatioN for low-risk pesticides), funded by the European Union's Horizon
Europe research and innovation programme under grant agreement No. 101084163.

Much of v2.0.0 — merging FunFlux 1.x, Funnotator and FunFluxL into one workflow, adding the hybrid
mode and writing this documentation — was done with [Claude Code](https://claude.com/claude-code).
Anthropic provided six months of Claude Max through their Open Source programme, and the scale of
the v2 rewrite would not have been realistic without it.

## References

1. Abarenkov, K., et al. (2024). The UNITE database for molecular identification and taxonomic communication of fungi and other eukaryotes: sequences, taxa and classifications reconsidered. *Nucleic Acids Research*, 52(D1), D791–D797. <https://doi.org/10.1093/nar/gkad1039>
2. Abarenkov, K., et al. (2025). *UNITE USEARCH/UTAX/SINTAX release for eukaryotes*. Version 19.02.2025. UNITE Community. <https://doi.org/10.15156/BIO/3301246>
3. Ashburner, M., et al. (2000). Gene Ontology: tool for the unification of biology. *Nature Genetics*, 25(1), 25–29. <https://doi.org/10.1038/75556>
4. Banno, S., Yamashita, K., Fukumori, F., Okada, K., Uekusa, H., Takagaki, M., Kimura, M., & Fujimura, M. (2009). Characterization of QoI resistance in *Botrytis cinerea* and identification of two types of mitochondrial cytochrome *b* gene. *Plant Pathology*, 58(1), 120–129. <https://doi.org/10.1111/j.1365-3059.2008.01909.x>
5. Bédard, C., et al. (2025). FungAMR: a comprehensive database for investigating fungal mutations associated with antimicrobial resistance. *Nature Microbiology*, 10(9), 2338–2352. <https://doi.org/10.1038/s41564-025-02084-7>
6. Bengtsson-Palme, J., et al. (2013). Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. *Methods in Ecology and Evolution*, 4(10), 914–919. <https://doi.org/10.1111/2041-210X.12073>
7. Blin, K., et al. (2025). antiSMASH 8.0: extended gene cluster detection capabilities and analyses of chemistry, enzymology, and regulation. *Nucleic Acids Research*, 53(W1), W32–W38. <https://doi.org/10.1093/nar/gkaf334>
8. Blum, M., et al. (2021). The InterPro protein families and domains database: 20 years on. *Nucleic Acids Research*, 49(D1), D344–D354. <https://doi.org/10.1093/nar/gkaa977>
9. Bouras, G., et al. (2024). How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies. *Microbial Genomics*, 10(6). <https://doi.org/10.1099/mgen.0.001254>
10. Bowyer, P., Currin, A., Delneri, D., & Fraczek, M. G. (2022). Telomere-to-telomere genome sequence of the model mould pathogen *Aspergillus fumigatus*. *Nature Communications*, 13, 5394. <https://doi.org/10.1038/s41467-022-32924-7>
11. Brejová, B., Hodorová, V., Lichancová, H., Gafurov, A., Bujna, D., Brázdovič, F., et al. (2026). Noncanonical chromosomal-end-specific telomeric repeat arrays in naturally telomerase-negative yeasts. *Nucleic Acids Research*, 54(10), gkag566. <https://doi.org/10.1093/nar/gkag566>
12. Brown, M. R., Gonzalez de la Rosa, P. M., & Blaxter, M. (2025). tidk: a toolkit to rapidly identify telomeric repeats from genomic datasets. *Bioinformatics*, 41(2), btaf049. <https://doi.org/10.1093/bioinformatics/btaf049>
13. Buchfink, B., et al. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nature Methods*, 18(4), 366–368. <https://doi.org/10.1038/s41592-021-01101-x>
14. Camacho, C., et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. <https://doi.org/10.1186/1471-2105-10-421>
15. Cantalapiedra, C. P., et al. (2021). eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Molecular Biology and Evolution*, 38(12), 5825–5829. <https://doi.org/10.1093/molbev/msab293>
16. Červenák, F., Sepšiová, R., Nosek, J., & Tomáška, Ľ. (2021). Step-by-step evolution of telomeres: lessons from yeasts. *Genome Biology and Evolution*, 13(2), evaa268. <https://doi.org/10.1093/gbe/evaa268>
17. Chan, P. P., et al. (2021). tRNAscan-SE 2.0: improved detection and functional classification of transfer RNA genes. *Nucleic Acids Research*, 49(16), 9077–9096. <https://doi.org/10.1093/nar/gkab688>
18. Chen, S. (2023). Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. *iMeta*, 2(2), e107. <https://doi.org/10.1002/imt2.107>
19. Chen, S. (2025). fastp 1.0: an ultra-fast all-round tool for FASTQ data quality control and preprocessing. *iMeta*, 4(5), e70078. <https://doi.org/10.1002/imt2.70078>
20. Chen, S., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890. <https://doi.org/10.1093/bioinformatics/bty560>
21. Cheng, H., et al. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. *Nature Methods*, 18(2), 170–175. <https://doi.org/10.1038/s41592-020-01056-5>
22. Cheng, H., et al. (2022). Haplotype-resolved assembly of diploid genomes without parental data. *Nature Biotechnology*, 40(9), 1332–1335. <https://doi.org/10.1038/s41587-022-01261-x>
23. Cheng, H., et al. (2024). Scalable telomere-to-telomere assembly for diploid and polyploid genomes with double graph. *Nature Methods*, 21(6), 967–970. <https://doi.org/10.1038/s41592-024-02269-8>
24. Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. <https://doi.org/10.1093/gigascience/giab008>
25. De Coster, W., & Rademakers, R. (2023). NanoPack2: population-scale evaluation of long-read sequencing data. *Bioinformatics*, 39(5), btad311. <https://doi.org/10.1093/bioinformatics/btad311>
26. Edgar, R. C. (2016). SINTAX: a simple non-Bayesian taxonomy classifier for 16S and ITS sequences. *bioRxiv*, 074161. <https://doi.org/10.1101/074161>
27. Edman, J. C. (1992). Isolation of telomerelike sequences from *Cryptococcus neoformans* and their use in high-efficiency transformation. *Molecular and Cellular Biology*, 12(6), 2777–2783. <https://doi.org/10.1128/mcb.12.6.2777-2783.1992>
28. Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. <https://doi.org/10.1093/bioinformatics/btw354>
29. Flynn, J. M., et al. (2020). RepeatModeler2 for automated genomic discovery of transposable element families. *Proceedings of the National Academy of Sciences*, 117(17), 9451–9457. <https://doi.org/10.1073/pnas.1921046117>
30. Frith, M. C. (2011). A new repeat-masking method enables specific detection of homologous sequences. *Nucleic Acids Research*, 39(4), e23. <https://doi.org/10.1093/nar/gkq1212>
31. Gabriel, L., Becker, F., Hoff, K. J., & Stanke, M. (2024). Tiberius: end-to-end deep learning with an HMM for gene prediction. *Bioinformatics*, 40(12), btae685. <https://doi.org/10.1093/bioinformatics/btae685>
32. The Gene Ontology Consortium. (2026). The Gene Ontology knowledgebase in 2026. *Nucleic Acids Research*, 54(D1), D1779–D1792. <https://doi.org/10.1093/nar/gkaf1292>
33. Grasso, V., Palermo, S., Sierotzki, H., Garibaldi, A., & Gisi, U. (2006). Cytochrome *b* gene structure and consequences for resistance to Qo inhibitor fungicides in plant pathogens. *Pest Management Science*, 62(6), 465–472. <https://doi.org/10.1002/ps.1236>
34. Haas, B. J., et al. (2008). Automated eukaryotic gene structure annotation using EVidenceModeler and the Program to Assemble Spliced Alignments. *Genome Biology*, 9(1), R7. <https://doi.org/10.1186/gb-2008-9-1-r7>
35. Huerta-Cepas, J., et al. (2019). eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. *Nucleic Acids Research*, 47(D1), D309–D314. <https://doi.org/10.1093/nar/gky1085>
36. Jones, P., et al. (2014). InterProScan 5: genome-scale protein function classification. *Bioinformatics*, 30(9), 1236–1240. <https://doi.org/10.1093/bioinformatics/btu031>
37. Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780. <https://doi.org/10.1093/molbev/mst010>
38. Kolmogorov, M., et al. (2019). Assembly of long, error-prone reads using repeat graphs. *Nature Biotechnology*, 37(5), 540–546. <https://doi.org/10.1038/s41587-019-0072-8>
39. Korf, I. (2004). Gene finding in novel genomes. *BMC Bioinformatics*, 5, 59. <https://doi.org/10.1186/1471-2105-5-59>
40. Laetsch, D. R., & Blaxter, M. L. (2017). BlobTools: interrogation of genome assemblies. *F1000Research*, 6, 1287. <https://doi.org/10.12688/f1000research.12232.1>
41. Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. *Nature Methods*, 9(4), 357–359. <https://doi.org/10.1038/nmeth.1923>
42. Langmead, B., et al. (2019). Scaling read aligners to hundreds of threads on general-purpose processors. *Bioinformatics*, 35(3), 421–432. <https://doi.org/10.1093/bioinformatics/bty648>
43. Levis, C., Giraud, T., Dutertre, M., Fortini, D., & Brygoo, Y. (1997). Telomeric DNA of *Botrytis cinerea*: a useful tool for strain identification. *FEMS Microbiology Letters*, 157(2), 267–272. <https://doi.org/10.1111/j.1574-6968.1997.tb12783.x>
44. Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, 1303.3997. <https://doi.org/10.48550/arXiv.1303.3997>
45. Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. <https://doi.org/10.1093/bioinformatics/bty191>
46. Li, H. (2021). New strategies to improve minimap2 alignment accuracy. *Bioinformatics*, 37(23), 4572–4574. <https://doi.org/10.1093/bioinformatics/btab705>
47. Li, H. (2023). Protein-to-genome alignment with miniprot. *Bioinformatics*, 39(1), btad014. <https://doi.org/10.1093/bioinformatics/btad014>
48. Lomsadze, A., et al. (2005). Gene identification in novel eukaryotic genomes by self-training algorithm. *Nucleic Acids Research*, 33(20), 6494–6506. <https://doi.org/10.1093/nar/gki937>
49. Lue, N. F. (2021). Duplex telomere-binding proteins in fungi with canonical telomere repeats: new lessons in the rapid evolution of telomere proteins. *Frontiers in Genetics*, 12, 638790. <https://doi.org/10.3389/fgene.2021.638790>
50. Lyčka, M., Bubeník, M., Závodník, M., Peska, V., Fajkus, P., Demko, M., Fajkus, J., & Fojtová, M. (2024). TeloBase: a community-curated database of telomere sequences across the tree of life. *Nucleic Acids Research*, 52(D1), D311–D321. <https://doi.org/10.1093/nar/gkad672>
51. Majoros, W. H., et al. (2004). TigrScan and GlimmerHMM: two open source ab initio eukaryotic gene-finders. *Bioinformatics*, 20(16), 2878–2879. <https://doi.org/10.1093/bioinformatics/bth315>
52. McEachern, M. J., & Blackburn, E. H. (1994). A conserved sequence motif within the exceptionally diverse telomeric sequences of budding yeasts. *Proceedings of the National Academy of Sciences USA*, 91(8), 3453–3457. <https://doi.org/10.1073/pnas.91.8.3453>
53. Mikheenko, A., et al. (2018). Versatile genome assembly evaluation with QUAST-LG. *Bioinformatics*, 34(13), i142–i150. <https://doi.org/10.1093/bioinformatics/bty266>
54. Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. <https://doi.org/10.12688/f1000research.29032.2>
55. Morgulis, A., et al. (2008). Database indexing for production MegaBLAST searches. *Bioinformatics*, 24(16), 1757–1764. <https://doi.org/10.1093/bioinformatics/btn322>
56. Okonechnikov, K., et al. (2016). Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. *Bioinformatics*, 32(2), 292–294. <https://doi.org/10.1093/bioinformatics/btv566>
57. Oxford Nanopore Technologies. *Medaka*. <https://github.com/nanoporetech/medaka>
58. Palmer, J. M., & Stajich, J. (2020). *Funannotate v1.8.1: eukaryotic genome annotation*. Zenodo. <https://doi.org/10.5281/zenodo.1134477>
59. Podlevsky, J. D., Bley, C. J., Omana, R. V., Qi, X., & Chen, J. J.-L. (2008). The Telomerase Database. *Nucleic Acids Research*, 36(Database issue), D339–D343. <https://doi.org/10.1093/nar/gkm700>
60. Prjibelski, A., et al. (2020). Using SPAdes de novo assembler. *Current Protocols in Bioinformatics*, 70(1), e102. <https://doi.org/10.1002/cpbi.102>
61. Quijada, N. M., et al. (2026). *ChroQueTas* (version 1.0.1). <https://github.com/nmquijada/ChroQueTas>
62. Rawlings, N. D., et al. (2014). MEROPS: the database of proteolytic enzymes, their substrates and inhibitors. *Nucleic Acids Research*, 42(D1), D503–D509. <https://doi.org/10.1093/nar/gkt953>
63. Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: the European Molecular Biology Open Software Suite. *Trends in Genetics*, 16(6), 276–277. <https://doi.org/10.1016/S0168-9525(00)02024-2>
64. Rognes, T., et al. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ*, 4, e2584. <https://doi.org/10.7717/peerj.2584>
65. Sepsiova, R., Necasova, I., Willcox, S., Prochazkova, K., Gorilak, P., Nosek, J., Griffith, J. D., & Tomaska, L. (2016). Evolution of telomeres in *Schizosaccharomyces pombe* and its possible relationship to the diversification of telomere binding proteins. *PLoS ONE*, 11(4), e0154225. <https://doi.org/10.1371/journal.pone.0154225>
66. Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE*, 11(10), e0163962. <https://doi.org/10.1371/journal.pone.0163962>
67. Slater, G. S. C., & Birney, E. (2005). Automated generation of heuristics for biological sequence comparison. *BMC Bioinformatics*, 6, 31. <https://doi.org/10.1186/1471-2105-6-31> (Exonerate)
68. Smit, A. F. A., Hubley, R., & Green, P. *RepeatMasker*. <http://www.repeatmasker.org>
69. Stanke, M., et al. (2006). AUGUSTUS: ab initio prediction of alternative transcripts. *Nucleic Acids Research*, 34(Web Server issue), W435–W439. <https://doi.org/10.1093/nar/gkl200>
70. Stanke, M., et al. (2008). Using native and syntenically mapped cDNA alignments to improve de novo gene finding. *Bioinformatics*, 24(5), 637–644. <https://doi.org/10.1093/bioinformatics/btn013>
71. Sun, Q., Wang, H., Tao, S., & Xi, X. (2023). Large-scale detection of telomeric motif sequences in genomic data using TelFinder. *Microbiology Spectrum*, 11(2), e03928-22. <https://doi.org/10.1128/spectrum.03928-22>
72. Tegenfeldt, F., et al. (2025). OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes. *Nucleic Acids Research*, 53(D1), D516–D522. <https://doi.org/10.1093/nar/gkae987>
73. Ter-Hovhannisyan, V., et al. (2008). Gene prediction in novel fungal genomes using an ab initio algorithm with unsupervised training. *Genome Research*, 18(12), 1979–1990. <https://doi.org/10.1101/gr.081612.108>
74. Teufel, F., et al. (2022). SignalP 6.0 predicts all five types of signal peptides using protein language models. *Nature Biotechnology*, 40(7), 1023–1025. <https://doi.org/10.1038/s41587-021-01156-3>
75. The UniProt Consortium. (2025). UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Research*, 53(D1), D609–D617. <https://doi.org/10.1093/nar/gkae1010>
76. Wheeler, T. J., & Eddy, S. R. (2013). nhmmer: DNA homology search with profile HMMs. *Bioinformatics*, 29(19), 2487–2489. <https://doi.org/10.1093/bioinformatics/btt403>
77. Wick, R. R. *Filtlong*. <https://github.com/rrwick/Filtlong>
78. Wick, R. R., & Holt, K. E. (2022). Polypolish: short-read polishing of long-read bacterial genome assemblies. *PLOS Computational Biology*, 18(1), e1009802. <https://doi.org/10.1371/journal.pcbi.1009802>
79. Zdouc, M. M., et al. (2025). MIBiG 4.0: advancing biosynthetic gene cluster curation through global collaboration. *Nucleic Acids Research*, 53(D1), D678–D690. <https://doi.org/10.1093/nar/gkae1115>
80. Zhang, Z., et al. (2000). A greedy algorithm for aligning DNA sequences. *Journal of Computational Biology*, 7(1–2), 203–214. <https://doi.org/10.1089/10665270050081478>
81. Zheng, J., et al. (2023). dbCAN3: automated carbohydrate-active enzyme and substrate annotation. *Nucleic Acids Research*, 51(W1), W115–W121. <https://doi.org/10.1093/nar/gkad328>

## References of the validation page { #validation-references }

[Validation on public genomes](validation.md) cites these works besides the tools above. On that
page, Liftoff, gffread and gffcompare compared FunFlux's genes with the reference. These three
tools are not part of FunFlux.

- Liu, M., et al. (2026). A chromosome-scale nuclear genome and complete mitogenome of the bio-control fungus *Cordyceps cateniannulata*. *Scientific Data*, 13, 871. <https://doi.org/10.1038/s41597-026-07231-1> (the reference genome assembly ASM5650595v1, GenBank GCA_056505955.1)
- Pertea, G., & Pertea, M. (2020). GFF Utilities: GffRead and GffCompare. *F1000Research*, 9, ISCB Comm J-304. <https://doi.org/10.12688/f1000research.23297.2> (gffread 0.12.9, gffcompare 0.12.10)
- Shumate, A., & Salzberg, S. L. (2021). Liftoff: accurate mapping of gene annotations. *Bioinformatics*, 37(12), 1639–1643. <https://doi.org/10.1093/bioinformatics/btaa1016> (Liftoff 1.6.3)
- Tiberius `fungi` model: its parameter file `model_cfg/fungi.yaml` (it dates the model 2026-02-06; the model file FunFlux downloads is of 8 February 2026; weights trained by K. Hoff), which lists the species used for training, validation and testing; *Botrytis cinerea* is among the training species, no *Cordyceps* species is, and *Akanthomyces muscarius* and *Beauveria bassiana* are its only training species of the family Cordycipitaceae. <https://github.com/Gaius-Augustus/Tiberius/blob/0db79a92ab78037a636439d7b880ac942e931084/model_cfg/fungi.yaml>
- Van Kan, J. A., et al. (2017). A gapless genome sequence of the fungus *Botrytis cinerea*. *Molecular Plant Pathology*, 18(1), 75–89. <https://doi.org/10.1111/mpp.12384> (the reference genome assembly ASM14353v4, RefSeq GCF_000143535.2)
