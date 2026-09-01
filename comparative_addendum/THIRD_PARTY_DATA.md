# Third-party data in the comparative addendum

The repository software license does not relicense downloaded source files. Each input remains subject to its source terms. Exact paths, accessions, URLs, byte counts, and SHA-256 hashes are listed in [`input_sources.tsv`](input_sources.tsv). The complete base-plus-addendum input set is listed in [`combined_input_manifest.tsv`](combined_input_manifest.tsv).

The GitHub Release input archive exists for scientific verification. It preserves the exact bytes used for the committed results if an upstream server changes a compressed wrapper, retires a URL, or returns a different record representation.

| Source | Added material | Terms and credit |
| --- | --- | --- |
| NCBI RefSeq | Current *Bombus terrestris* annotation and proteome for `GCF_910591885.1` | NCBI places no restrictions on use or distribution of its molecular databases. Submitters may retain other rights. Follow the [NCBI molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) and cite the accession. |
| NCBI Protein | Seven additional *Bombus* exchanger records and six other-bee source records | Follow the [NCBI molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) and cite every accession listed in `input_sources.tsv` and the generated sequence-metadata tables. |
| Oxford University Press via Europe PMC | Lin et al. (2021) supplementary archive | The article and supplement are published under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/). Non-commercial redistribution requires attribution. Cite DOI [10.1093/gbe/evab227](https://doi.org/10.1093/gbe/evab227). |

The base repository has additional NCBI, GEO, Genome Warehouse, Oxford University Press, and UniProt inputs. Their terms are recorded in the root [`THIRD_PARTY_DATA.md`](../THIRD_PARTY_DATA.md).

## Citations

- Lin D. et al. (2021), *Comparative Genomics Reveals Recent Adaptive Evolution in Himalayan Giant Honeybee Apis laboriosa*: <https://doi.org/10.1093/gbe/evab227>
- NCBI RefSeq assembly `GCF_910591885.1`: <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_910591885.1/>
- NCBI Protein accessions are enumerated in [`input_sources.tsv`](input_sources.tsv) and [`results/nhe3_extended_bee_sequence_metadata.tsv`](results/nhe3_extended_bee_sequence_metadata.tsv).
