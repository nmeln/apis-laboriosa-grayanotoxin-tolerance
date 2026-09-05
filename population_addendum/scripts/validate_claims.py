"""Check biological denominators and the rules behind every reported call."""
from collections import Counter
import json
from common import OUT,table


def main():
    inventory=table(OUT/'sample_inventory.tsv');runs={r['run_accession'] for r in inventory}
    assert len(runs)==57
    assert Counter(r['scientific_name'] for r in inventory)=={'Apis laboriosa':29,'Apis dorsata':28}
    calls=table(OUT/'sample_codon_calls.tsv');assert len(calls)==57*9*2*2
    assert len({(r['run_accession'],r['assembly'],r['focal_AL_position'],r['profile']) for r in calls})==len(calls)
    assert {r['run_accession'] for r in calls}==runs
    for r in calls:
        assert int(r['forward'])+int(r['reverse'])==int(r['depth'])
        counts={k:int(v) for k,v in (x.split(':') for x in r['codon_counts'].split(';'))} if r['codon_counts'] else {}
        assert sum(counts.values())==int(r['depth'])
        if r['status']=='PASS':
            assert int(r['depth'])>=(10 if r['profile']=='primary' else 8)
            assert float(r['gq'])>=30 and int(r['forward']) and int(r['reverse'])
            a,b=r['best_codons'].split('/')
            if a!=b:assert min(counts[a],counts[b])>=3 and min(counts[a],counts[b])/int(r['depth'])>=.2
    consensus=table(OUT/'consensus_genotypes.tsv');assert len(consensus)==57*9*2
    for r in consensus:
        if r['status']=='PASS':assert r['AL_status']==r['AD_status']=='PASS' and r['AL_call']==r['AD_call']==r['codons']
        else:assert not r['codons'] and not r['amino_acids'] and not r['reference_state_copies']
    species=table(OUT/'species_site_summary.tsv');assert len(species)==2*9*2
    for r in species:
        assert int(r['deposited_workers'])==({'Apis laboriosa':29,'Apis dorsata':28}[r['species']])
        assert int(r['callable_workers'])+int(r['missing_workers'])==int(r['deposited_workers'])
        assert int(r['total_callable_allele_copies'])==2*int(r['callable_workers'])
        assert int(r['reference_state_homozygotes'])+int(r['reference_state_heterozygotes'])+int(r['workers_without_reference_state'])==int(r['callable_workers'])
    assert len(table(OUT/'bcftools_base_crosscheck.tsv'))==57*9*3*2
    for name in ['capture_tool_validation.json','genotype_tool_validation.json']:
        d=json.loads((OUT/name).read_text())
        if name=='genotype_tool_validation.json':assert d['all_passed']
        else:assert d['all_positive_pairs_retained'] and d['all_negative_pairs_excluded'] and d['retained_R1']==d['retained_R2']==6912 and d['read_sequences_and_qualities_preserved'] and d['unmatched_mates_preserved']
    d=json.loads((OUT/'collection_summary.json').read_text());assert d['verified_libraries']==57 and d['total_source_compressed_bytes']==269858255055
    assert json.loads((OUT/'independent_capture_reproduction.json').read_text())['all_match']
    # Frozen observed findings, checked separately from structural invariants.
    for r in species:
        if r['role']=='primary':
            n=28 if r['species']=='Apis laboriosa' else {'254':27,'549':26,'1134':28}[r['focal_AL_position']]
            assert int(r['callable_workers'])==n and int(r['reference_state_homozygotes'])==n
            assert int(r['reference_state_allele_copies'])==2*n
    observed={(r['species'],r['focal_AL_position'],r['profile']):r for r in species}
    assert observed['Apis laboriosa','452','primary']['amino_acid_genotype_counts']=='R/R:1;R/T:2;T/T:24'
    assert observed['Apis dorsata','822','primary']['amino_acid_genotype_counts']=='I/I:25;I/V:3'
    assert Counter(r['status'] for r in consensus if r['profile']=='primary')=={'PASS':497,'reference_quality_failure':16}
    assert Counter(r['comparison'] for r in table(OUT/'bcftools_base_crosscheck.tsv'))=={'agree':2982,'codon_call_filtered':96}
    profiles=table(OUT/'individual_primary_profiles.tsv')
    assert sum(r['callable_primary_sites']=='3' for r in profiles)==53
    assert {r['run_accession'] for r in consensus if r['profile']=='primary' and r['role']=='primary' and r['status']!='PASS'}=={'SRR23343440','SRR23343455','SRR23343460','SRR23343476'}
    print('Validated 57-worker denominators, 2,052 codon calls, missing-data handling and collection controls.')

if __name__=='__main__':main()
