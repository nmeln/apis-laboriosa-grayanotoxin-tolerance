# Source data and credit

This addendum reuses data published by the original investigators. The code's
MIT license does not replace the terms attached to those data.

| Source | Preserved material | Credit and terms |
| --- | --- | --- |
| ENA / SRA, PRJNA931733 | Exact selected paired FASTQ records, original run and BioSample metadata, full-source URLs, byte counts and MD5s | Cite Cao et al. (2023) and the individual run accessions. [ENA policies](https://www.ebi.ac.uk/ena/browser/about/policies) and [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/) apply. EMBL-EBI adds no restrictions beyond those of the data owners and expects attribution. |
| NCBI RefSeq / GenBank | Complete assemblies GCF_014066325.1, GCF_000469605.1 and GCA_009792835.1 | Retain the accessions and original assembly study citations. See the [NCBI molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/). |
| Genome Warehouse | Complete eastern-Yunnan laboriosa assembly GWHAOTM00000000 | Cite Cao et al. (2023), the accession and Genome Warehouse. Its stated academic-use terms and commercial-use restriction are described in the [NGDC policy](https://ngdc.cncb.ac.cn/policies?lang=en). |
| Cao et al. (2023), Oxford University Press | Article XML and the supplementary workbook used to reconcile sample aliases and localities | The article specifies [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). Retain attribution and the non-commercial restriction. |

Cao L. et al. (2023). *Population Structure, Demographic History, and Adaptation
of Giant Honeybees in China Revealed by Population Genomic Data.* Genome
Biology and Evolution. [DOI: 10.1093/gbe/evad025](https://doi.org/10.1093/gbe/evad025).

Genome sources and their studies are also documented in
[`../candidate_followup/REPORT.md`](../candidate_followup/REPORT.md) and
[`../THIRD_PARTY_DATA.md`](../THIRD_PARTY_DATA.md).

The input archive contains all selected original read records used in the
analysis and all four complete genome assemblies. It does not contain the
approximately 270 GB of complete population FASTQs. Those remain separately
retrievable and verifiable through the supplied complete-file collection stage.

The tool installers retrieve upstream or Bioconda software distributions with
pinned hashes. Those tools keep their original licenses: minimap2 uses MIT;
BBMap/BBDuk uses BSD-3-Clause-LBNL; SAMtools/BCFtools and pysam carry the licenses
included in their respective distributions. They are dependencies, with no
claim of authorship by this project.
