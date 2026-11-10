# Adapter sequences for the contig adapter screen

Rule blast_adapters (`workflow/rules/shared/25_contig_summary.smk`) searches every draft
assembly for the sequences in this folder, together with two more sets that are not shipped:
the Illumina records of NCBI UniVec, which rule univec_illumina downloads into the run
(Illumina's oligonucleotide sequences are copyrighted), and, for an ONT sample whose reads were
trimmed, the sequence that was trimmed (rule find_ont_adapter).

| File | Records | Source |
|---|---|---|
| `ont_chtd.fasta` | 196 | Oxford Nanopore Chemistry Technical Document, version CHTD_500_v1_revAV_24July2026 (<https://nanoporetech.com/document/chemistry-technical-document>), sections 14 and 15: the adapter top strand (shared by the ligation, native and rapid adapters), the ligation and native adapter bottom strands, the 3' flanking sequence of the rapid barcodes, the 96 native barcodes NB01-96 between their flanks (AAGGTTAA-barcode-CAGCACCT) and the 96 rapid barcodes RB01-96 between theirs. ONT publishes the rapid barcodes only in part. |
| `univec_pacbio_mgi.fasta` | 5 | NCBI UniVec, the file of 2024-08-21 at <https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec>: the PacBio blunt (SMRTbell) adapter NGB00972.1, the PacBio C2 primer NGB00973.1, the PacBio ULI amplification adapter NGB03000.1 and the two MGI DNBSEQ adaptors NGB05000.1 and NGB05001.1. UniVec is a work of the US Government. |

Every FASTA header names its source. The native barcodes were read from the document's table,
and its forward and reverse columns were checked to be reverse complements of each other.
