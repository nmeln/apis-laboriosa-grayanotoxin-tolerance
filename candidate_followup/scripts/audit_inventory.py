"""Summarize saved archive queries; paired GCA/GCF records are one assembly."""
import json
from collections import Counter
from common import INPUT, ROOT, OUT, table, write


def main():
    assemblies, rna, studies = [], [], []
    geo_tissue = {r['sample']: r['tissue_metadata'] for r in table(ROOT / 'transcriptomic_addendum/results/sample_tissue_metadata.tsv')}
    for taxon in ['dorsata','laboriosa']:
        data = json.loads((INPUT / (taxon + '_assemblies.json')).read_text())
        assert len(data['reports']) == int(data['total_count'])
        for r in data['reports']:
            info, stats = r['assembly_info'], r['assembly_stats']; sample = info.get('biosample',{})
            attrs = {x['name']:x['value'] for x in sample.get('attributes',[])}
            assemblies.append(dict(species='Apis ' + taxon, accession=r['accession'], paired_accession=r.get('paired_accession',''),
                                   independent_unit=info['bioproject_accession'] + ':' + sample['accession'],
                                   assembly_name=info['assembly_name'], assembly_level=info['assembly_level'],
                                   release_date=info['release_date'], bioproject=info['bioproject_accession'], biosample=sample['accession'],
                                   location=sample.get('geo_loc_name',''), collection_date=sample.get('collection_date',''),
                                   tissue=attrs.get('tissue','not specified'), sample_title=sample.get('description',{}).get('title',''),
                                   sequencing_technology=info.get('sequencing_tech',''),
                                   assembled_nt=stats['total_sequence_length'], contig_n50=stats['contig_n50'], scaffold_n50=stats['scaffold_n50']))
        records = table(INPUT / (taxon + '_ena_runs.tsv'))
        combinations = Counter((r['study_accession'],r['library_strategy'],r['library_source']) for r in records)
        for (study, strategy, source), count in sorted(combinations.items()):
            group = [r for r in records if r['study_accession'] == study and r['library_strategy'] == strategy and r['library_source'] == source]
            studies.append(dict(species='Apis '+taxon, study=study, strategy=strategy, library_source=source, runs=count,
                                distinct_samples=len({r['sample_accession'] for r in group}), title=group[0]['study_title'],
                                use='Host expression' if source == 'TRANSCRIPTOMIC' else 'Does not measure host RNA expression'))
        ena = table(INPUT / (taxon + '_rna_subtree.tsv'))
        sra = json.loads((INPUT / (taxon + '_rna_ncbi.json')).read_text())['esearchresult']
        rna.append(dict(species='Apis '+taxon, search_date='2026-09-05', ena_taxonomic_subtree_transcriptomic_runs=len(ena),
                        ncbi_sra_transcriptomic_records=int(sra['count']), runs=';'.join(sorted(r['run_accession'] for r in ena)),
                        tissue=geo_tissue[{'laboriosa':'GSM3757258','dorsata':'GSM3757259'}[taxon]], independent_usable_gut_host_RNA_datasets_found=0,
                        scope='Saved ENA taxonomic-subtree transcriptomic search and NCBI SRA Organism:exp/Source search; not proof no unindexed or unsubmitted data exist'))
    write('genome_assembly_inventory.tsv', assemblies)
    write('sequencing_study_inventory.tsv', studies)
    write('focal_RNA_inventory.tsv', rna)
    related = [
        dict(species='Apis mellifera', dataset='10.1038/s41598-025-26662-1', tissue='Dissected midgut; brush border membrane vesicles',
             material='Protein mass spectrometry', availability='Downloaded article, complete supplementary workbook and Figshare file inventory',
             analysis='Exact OG0000499 homologue detected in 3 adult and 3 larval preparations; replicate-table extraction',
             source='https://doi.org/10.1038/s41598-025-26662-1'),
        dict(species='Apis mellifera; Apis cerana japonica', dataset='10.1038/s41597-025-05279-z; DRA016719; DRA014424', tissue='Whole individual, including abdomen',
             material='Host RNA-seq over development', availability='Primary article identifies deposited raw reads and public TPM matrices',
             analysis='Availability and tissue checked; candidate expression not reanalyzed', source='https://doi.org/10.1038/s41597-025-05279-z'),
        dict(species='Apis mellifera', dataset='10.1111/1744-7917.70157', tissue='Midgut, ileum, rectum', material='Single-nucleus host RNA atlas',
             availability='Publisher abstract and supplement index available; attempted supplements returned HTTP 403',
             analysis='No cell-level candidate inference or reanalysis', source='https://doi.org/10.1111/1744-7917.70157'),
    ]
    write('related_bee_abdominal_datasets.tsv', related)
    summary = dict(dorsata_NCBI_assembly_records=sum(r['species']=='Apis dorsata' for r in assemblies),
                   dorsata_distinct_assembly_sources=len({r['independent_unit'] for r in assemblies if r['species']=='Apis dorsata'}),
                   focal_RNA_runs={r['species']:r['ena_taxonomic_subtree_transcriptomic_runs'] for r in rna},
                   dna_rna_distinction='Host cells in the head or thorax contain the gut genes. Their DNA can reconstruct a gut transporter; gut activity needs RNA/protein from that tissue.',
                   removal_reason_source='Cao et al. 2023 explicitly states abdomen removal to avoid contamination for DNA extraction. GEO records belly exclusion for RNA; the paper does not separately give the RNA-dissection rationale.',
                   reason_inference='Reducing gut microbes and food-derived material is a plausible explanation of contamination avoidance, not an additional statement of intent by the authors.')
    (OUT / 'inventory_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')


if __name__ == '__main__': main()
