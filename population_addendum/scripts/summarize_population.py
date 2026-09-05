"""Combine two-reference calls with explicit sample and locality denominators."""
import argparse
from collections import Counter, defaultdict
import json
from common import OUT, INPUT, WORK, SITES, table, write, json_write


def count_text(values):
    c=Counter(values)
    return ';'.join(f'{k}:{c[k]}' for k in sorted(c))


def main():
    p=argparse.ArgumentParser();p.add_argument('--prefix',default='');a=p.parse_args()
    samples={r['run_accession']:r for r in table(OUT/'sample_inventory.tsv')}
    calls=table(OUT/(a.prefix+'sample_codon_calls.tsv'))
    grouped=defaultdict(list)
    for r in calls: grouped[(r['run_accession'],r['focal_AL_position'],r['profile'])].append(r)
    reference={ (r['assembly'],r['focal_AL_position']):r for r in table(OUT/'target_coordinates.tsv') }
    consensus=[]
    for (run,site,profile),rows in sorted(grouped.items(),key=lambda x:(x[0][0],int(x[0][1]),x[0][2])):
        assert len(rows)==2 and {r['assembly'] for r in rows}=={'AL_Shangrila','AD_Malaysia'}
        rows=sorted(rows,key=lambda r:r['assembly'])
        sample=samples[run]
        same=len({r['best_codons'] for r in rows})==1
        passed=all(r['status']=='PASS' for r in rows)
        status='PASS' if same and passed else ('reference_disagreement' if not same else 'reference_quality_failure')
        species_ref='AL_Shangrila' if sample['scientific_name']=='Apis laboriosa' else 'AD_Malaysia'
        expected=reference[species_ref,site]['reference_amino_acid']
        aa=rows[0]['amino_acids'] if status=='PASS' else ''
        consensus.append(dict(run_accession=run,sample_alias=sample['sample_alias'],species=sample['scientific_name'],
                    location=sample['location'],focal_AL_position=site,role=SITES[int(site)],profile=profile,
                    status=status,codons=rows[0]['best_codons'] if status=='PASS' else '',amino_acids=aa,
                    species_reference_amino_acid=expected,reference_state_copies=aa.split('/').count(expected) if aa else '',
                    minimum_reference_depth=min(int(r['depth']) for r in rows),
                    maximum_reference_depth=max(int(r['depth']) for r in rows),
                    minimum_gq=f"{min(float(r['gq']) for r in rows):.3f}",
                    AD_call=rows[0]['best_codons'],AD_status=rows[0]['status'],
                    AL_call=rows[1]['best_codons'],AL_status=rows[1]['status']))
    write(a.prefix+'consensus_genotypes.tsv',consensus)
    for level in ['species','location']:
        groups=defaultdict(list)
        for r in consensus:
            key=(r['species'],r['location'] if level=='location' else 'all sampled locations',r['focal_AL_position'],r['profile'])
            groups[key].append(r)
        rows=[]
        for (species,location,site,profile),rr in sorted(groups.items(),key=lambda x:(x[0][0],x[0][1],int(x[0][2]),x[0][3])):
            good=[r for r in rr if r['status']=='PASS'];expected=rr[0]['species_reference_amino_acid']
            rows.append(dict(species=species,location=location,focal_AL_position=site,role=SITES[int(site)],profile=profile,
                  deposited_workers=len(rr),callable_workers=len(good),missing_workers=len(rr)-len(good),
                  expected_amino_acid=expected,reference_state_homozygotes=sum(r['reference_state_copies']==2 for r in good),
                  reference_state_heterozygotes=sum(r['reference_state_copies']==1 for r in good),
                  workers_without_reference_state=sum(r['reference_state_copies']==0 for r in good),
                  reference_state_allele_copies=sum(r['reference_state_copies'] for r in good),total_callable_allele_copies=2*len(good),
                  amino_acid_genotype_counts=count_text(r['amino_acids'] for r in good),
                  codon_genotype_counts=count_text(r['codons'] for r in good),missing_reason_counts=count_text(r['status'] for r in rr if r['status']!='PASS')))
        write(a.prefix+level+'_site_summary.tsv',rows)
    profiles=[]
    for run in sorted({r['run_accession'] for r in consensus}):
        rr=[r for r in consensus if r['run_accession']==run and r['role']=='primary' and r['profile']=='primary']
        assert len(rr)==3
        profiles.append(dict(run_accession=run,species=rr[0]['species'],location=rr[0]['location'],
                             callable_primary_sites=sum(r['status']=='PASS' for r in rr),
                             unphased_genotype_profile=';'.join(f"{r['focal_AL_position']}={r['amino_acids'] or 'missing'}" for r in rr)))
    write(a.prefix+'individual_primary_profiles.tsv',profiles)
    json_write(OUT/(a.prefix+'population_summary.json'),dict(analyzed_workers=len(profiles),
                   species_counts=dict(Counter(r['species'] for r in profiles)),
                   workers_callable_at_all_three=sum(r['callable_primary_sites']==3 for r in profiles),
                   primary_site_quality_failures=sum(r['profile']=='primary' and r['role']=='primary' and r['status']!='PASS' for r in consensus),
                   primary_site_reference_disagreements=sum(r['profile']=='primary' and r['role']=='primary' and r['status']=='reference_disagreement' for r in consensus),
                   species_fixation_established=False,toxin_transport_measured=False))
    print('Summarized',len(profiles),'individual workers without filling missing calls.')

if __name__=='__main__':main()
