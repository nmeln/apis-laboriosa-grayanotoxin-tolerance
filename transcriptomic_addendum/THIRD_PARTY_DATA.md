# Added data sources and reuse

Repository code uses the root MIT license. Source data retain their own terms.
The release preserves the exact files used in this analysis, with source URLs,
retrieval methods, sizes, and SHA-256 values in `input_sources.tsv`.

| Source files | Credit and source | Reuse context |
| --- | --- | --- |
| GEO SOFT record, ENA run metadata, two R1 prefixes | GSE130963 / PRJNA542114, associated with Lin et al. 2021, DOI 10.1093/gbe/evab227 | Public deposited bee RNA data; retain study and accession attribution |
| OR496406.1 and KY741959.1 GenBank records | Original submitters are identified in the preserved records | Public BQCV reference sequences; accession versions and record attribution retained |
| Drosophila kkv and Chs2 protein FASTAs | UniProtKB; accessions retained in headers | UniProt data are distributed under CC BY 4.0; credit UniProt and the underlying records |
| PMC10710133 XML and peerj-11-16238-s002.xlsx | Zhou et al. 2023, *Transcriptome analysis unveils the mechanisms of lipid metabolism response to grayanotoxin I stress in Spodoptera litura*, PeerJ, DOI 10.7717/peerj.16238 | Open-access article and associated supplement under CC BY 4.0; preserve author and article credit |

Source terms:
[NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/),
[EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/),
[UniProt license](https://www.uniprot.org/help/license),
[PeerJ article](https://doi.org/10.7717/peerj.16238).
NCBI and ENA hosting does not replace any applicable submitter rights.

The snapshot includes 2,000,000 R1 records from each run, compressed
deterministically. It does not mirror the full runs or their R2 reads. The
original input archive separately preserves the bee genomes, annotations,
transcript assemblies, and submitted abundance tables. Software binaries are
fetched from pinned upstream releases and are not bundled in the data archive.

The fly functional-comparison source is Bertran-Mas et al. 2025,
[DOI 10.1371/journal.pgen.1011847](https://doi.org/10.1371/journal.pgen.1011847).
Its article is cited; its full text is not included in this archive.
