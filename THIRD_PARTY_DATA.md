# Third-party input data

The repository's software license does not relicense the downloaded input
files. Each input remains subject to its source terms. Exact accessions, URLs,
byte counts, and SHA-256 hashes are listed in `data_sources.tsv`.

The GitHub Release input snapshot is provided for scientific verification. It
preserves the exact files used by the pipeline when upstream URLs or compressed
file wrappers change.

The cross-bee workflow adds the current *Bombus terrestris* reference, NCBI Protein panels, and the Lin et al. supplementary archive. Their exact records and terms are documented in [`comparative_addendum/THIRD_PARTY_DATA.md`](comparative_addendum/THIRD_PARTY_DATA.md).

| Source | Included material | Terms and required credit |
| --- | --- | --- |
| NCBI RefSeq and Protein | Bee annotations and proteomes; rat Nav1.4; bee Para accessions | NCBI places no restrictions on use or distribution of its molecular databases. Submitters may retain other rights. Follow the [NCBI molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) and cite the accessions. |
| NCBI GEO GSE130963 | Processed *A. laboriosa* and *A. dorsata* unigenes and submitted expression tables | NCBI places no restrictions on GEO data distribution. Submitters may retain other rights. Follow the [GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) and cite GSE130963 / PRJNA542114. |
| Genome Warehouse | Eastern-Yunnan assembly GWHAOTM00000000 | Genome Warehouse states that its database is free for academic use and asks commercial users to seek a license. See the [NGDC policy](https://ngdc.cncb.ac.cn/policies?lang=en), cite GWHAOTM00000000, and cite Cao et al. (2023). |
| Oxford University Press | Cao et al. (2023) supplementary workbook | The article and supplement use [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). Non-commercial redistribution requires attribution. Cite DOI [10.1093/gbe/evad025](https://doi.org/10.1093/gbe/evad025). |
| UniProt | P15390 JSON topology record | Copyrightable database content uses [CC BY 4.0](https://www.uniprot.org/help/license/). Cite UniProt and accession P15390. |

Primary genome citation:

> Lin D. et al. (2021). Comparative Genomics Reveals Recent Adaptive Evolution
> in Himalayan Giant Honeybee *Apis laboriosa*.
> <https://doi.org/10.1093/gbe/evab227>

Population assembly and supplement citation:

> Cao L. et al. (2023). Population Structure, Demographic History, and
> Adaptation of Giant Honeybees in China Revealed by Population Genomic Data.
> <https://doi.org/10.1093/gbe/evad025>

Processed transcriptome accession:

> NCBI GEO GSE130963 / BioProject PRJNA542114.
> <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA542114>
