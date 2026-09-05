"""Check scientific headline values directly against generated evidence."""
import json
import math
from Bio import SeqIO
from common import OUT, table


def main():
    inventory = json.loads((OUT/'inventory_summary.json').read_text())
    assert inventory['dorsata_NCBI_assembly_records'] == 3
    assert inventory['dorsata_distinct_assembly_sources'] == 2
    assert set(inventory['focal_RNA_runs'].values()) == {1}
    models = table(OUT/'genome_candidate_models.tsv')
    assert len(models) == 8 and {int(r['coding_exons']) for r in models} == {14}
    assert all(r['internal_stops'] == '0' and r['ambiguous_coding_bases'] == '0' for r in models)
    seq = {r.id:str(r.seq) for r in SeqIO.parse(OUT/'genome_recovered_cds.faa','fasta')}
    assert len(seq['AD_Malaysia']) == len(seq['AD_Thailand']) == 1333
    assert len(seq['AL_Shangrila']) == len(seq['AL_Pingbian']) == 1336
    assert sum(a != b for a,b in zip(seq['AD_Malaysia'],seq['AD_Thailand'])) == 1
    assert sum(a != b for a,b in zip(seq['AL_Shangrila'],seq['AL_Pingbian'])) == 1
    swaps = {r['reference']:r for r in table(OUT/'reference_swap_counts.tsv')}
    assert (int(swaps['AL']['AL_reads']),int(swaps['AL']['AD_reads'])) == (1161,54)
    assert (int(swaps['AD']['AL_reads']),int(swaps['AD']['AD_reads'])) == (1158,54)
    assert math.isclose(float(swaps['AL']['normalized_ratio']),17.43243243243243)
    assert math.isclose(float(swaps['AD']['normalized_ratio']),17.595441595441596)
    coverage = table(OUT/'candidate_read_coverage.tsv')
    assert all(float(r['covered_fraction']) > .99 for r in coverage if r['sample'] == 'AL')
    assert all(float(r['covered_fraction']) > .81 for r in coverage if r['sample'] == 'AD')
    markers = json.loads((OUT/'exact_shared_marker_summary.json').read_text())
    assert markers['specific_shared_markers'] == 1064 and markers['excluded_nonspecific_markers'] == 0
    assert markers['counts'] == {'AL':258,'AD':15}
    assert math.isclose(markers['raw_AL_AD_ratio'],17.2)
    changes = table(OUT/'candidate_coding_differences.tsv')
    stable = [r for r in changes if r['consistent_across_two_assemblies_each'] == 'True']
    assert sum(r['kind'] == 'substitution' for r in stable) == 7
    assert sum(r['kind'] == 'gap_column' for r in stable) == 3
    private = [r for r in stable if r['kind'] == 'substitution' and r['AL_state_absent_from_other_four_bees'] == 'True']
    assert [(int(r['AL_position_1based']),r['AL_Shangrila']) for r in private] == [(254,'L'),(549,'L'),(1134,'F')]
    assert all('Binding site:' not in r['projected_human_features'] for r in private)
    variable = next(r for r in changes if r['AL_position_1based'] == '822')
    assert variable['consistent_across_two_assemblies_each'] == 'False' and variable['AD_Thailand'] == 'I'
    transcript = table(OUT/'candidate_transcript_support.tsv')
    assert len(transcript) == 1 and transcript[0]['translated_aa'] == '1336' and transcript[0]['exact_reference_protein_match'] == 'True'
    ranking = table(OUT/'abc_reference_similarity.tsv')
    assert next(r for r in ranking if r['species'] == 'Homo_sapiens')['accession'] == 'O15439'
    assert next(r for r in ranking if r['species'] == 'Drosophila_melanogaster')['accession'] in ['Q9VLN6','I0C0M5']
    gut = [r for r in table(OUT/'gut_candidate_protein_evidence.tsv') if r['accession'] == 'XP_393750.5']
    assert len(gut) == 2 and all(r['replicates_detected'] == '3' and r['table_unique_peptides'] == '31' for r in gut)
    counts = [r for r in table(OUT/'gut_candidate_spectral_counts.tsv') if r['accession'] == 'XP_393750.5']
    assert [int(r['spectral_count']) for r in counts] == [66,52,57,7,39,67]
    assert all(r['retrieved'] == 'False' for r in table(OUT/'gut_peptide_export_inventory.tsv'))
    variants = list(SeqIO.parse(OUT/'untested_variant_panel.faa','fasta'))
    assert len(variants) == 6 and all('untested' in r.description for r in variants)
    print('Independent assemblies, reference swaps, exact markers, candidate residues, transcript and gut-protein evidence verified.')


if __name__ == '__main__': main()
