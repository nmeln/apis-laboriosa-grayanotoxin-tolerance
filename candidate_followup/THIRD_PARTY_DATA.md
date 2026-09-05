# Third-party sources and reuse

The repository software license applies to our code. It does not relicense
downloaded genomes, database records, articles, workbooks or other third-party
material. The input archive preserves original downloaded bytes. The derived
TSVs report factual observations and our calculations; they are not modified
versions of the original articles or workbook.

| Source | Credit and applicable terms |
| --- | --- |
| NCBI genome, RefSeq and BioSample records | Credit the accession submitters and the associated genome papers. NCBI redistributes depositor-supplied data; consult [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/). |
| ENA run metadata | Credit the original sequencing submitters and [ENA](https://www.ebi.ac.uk/ena/browser/about/data-usage). Earlier transcriptomic data terms also apply. |
| UniProt JSON sequences and annotations | [UniProt Consortium](https://www.uniprot.org/help/license), CC BY 4.0. Record versions and retrieved bytes are preserved. |
| Arora, Mishra and Bonning, gut membrane proteomics | [Article](https://doi.org/10.1038/s41598-025-26662-1), CC BY-NC-ND 4.0 according to the saved article. The complete original XML and supplementary workbook are preserved unmodified with attribution. No rights over those originals are claimed. |
| Arora et al. protein-export file inventory | [Figshare dataset, version 1](https://doi.org/10.6084/m9.figshare.28706999.v1), CC BY 4.0 according to its saved metadata. The CSV export contents were not successfully retrieved or archived. |
| Cao et al., population genetics and sampling methods | [Article](https://doi.org/10.1093/gbe/evad025), CC BY-NC 4.0 according to the saved XML. Preserved unmodified. |
| Lin et al., laboriosa genome | [Article](https://doi.org/10.1093/gbe/evab227), CC BY-NC 4.0 according to the saved XML. Preserved unmodified. |
| Yokoi et al., whole-individual developmental RNA resource | [Article](https://doi.org/10.1038/s41597-025-05279-z), CC BY 4.0. Availability/tissue audit only. |
| Chan et al., older organ protein atlas | [Article](https://doi.org/10.1101/gr.155994.113), CC BY-NC 3.0 after the first six months, according to the saved XML. Retained exploratory source only. No candidate result is based on this article. |

The full additional dorsata assemblies are associated with the
[2013 reference submission](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000469605.1/)
and [Oppenheim et al.](https://doi.org/10.1093/gbe/evz277).
Sequence analysis uses [miniprot](https://doi.org/10.1093/bioinformatics/btad014),
[minimap2](https://doi.org/10.1093/bioinformatics/bty191),
[FAMSA](https://doi.org/10.1038/srep33964), and Biopython. Their code retains
the respective upstream licenses. Executable source URLs and SHA-256 digests
are embedded in setup scripts.

The retained Nature HTML page is an unmodified copy of the same attributed
Arora article. Initial broad fly-protein queries are retained for completeness
and are not the canonical comparator panel. `Dmel_ABCC_domain.json` is the
panel used in the final similarity screen.
