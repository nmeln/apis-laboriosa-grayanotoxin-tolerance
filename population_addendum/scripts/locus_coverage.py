"""Describe missing-call coverage across coding sequence, introns and flanks.

Counts each qualifying DNA fragment once per covered base. This is coverage of
captured locus evidence, not an estimate of whole-genome library quality.
"""
import argparse
from collections import defaultdict
import json
from statistics import mean,median
import pysam
from common import ROOT,WORK,OUT,table,write


def main():
    p=argparse.ArgumentParser();p.add_argument('--label',default='population');p.add_argument('--prefix',default='');a=p.parse_args()
    destination=WORK/a.label
    runs=json.loads((destination/'mapping_inputs.json').read_text())['runs']
    samples={r['run_accession']:r for r in table(OUT/'sample_inventory.tsv')}
    regions=table(OUT/'capture_regions.tsv')
    exons=table(ROOT/'candidate_followup/results/genome_candidate_exons.tsv')
    models=table(ROOT/'candidate_followup/results/genome_candidate_models.tsv')
    summary,bins=[],[]
    for assembly in ['AL_Shangrila','AD_Malaysia']:
        region=next(r for r in regions if r['assembly']==assembly)
        model=next(r for r in models if r['assembly']==assembly)
        selected=[r for r in exons if r['assembly']==assembly and r['query']==model['query']]
        exon_positions={p for r in selected for p in range(int(r['start_1based'])-1,int(r['end_1based']))}
        left,right=min(exon_positions),max(exon_positions)+1
        start,end=int(region['start_1based'])-1,int(region['end_1based'])
        categories={'coding':exon_positions,'intronic':set(range(left,right))-exon_positions,
                    'left_flank':set(range(start,left)),'right_flank':set(range(right,end))}
        positions=defaultdict(lambda:defaultdict(set))
        with pysam.AlignmentFile(str(destination/(assembly+'.bam')),'rb') as bam:
            for read in bam.fetch(region['scaffold'],start,end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality<30:continue
                run=read.get_tag('RG')
                assert run in runs
                positions[run][read.query_name].update(pos for q,pos in read.get_aligned_pairs(matches_only=True) if start<=pos<end and read.query_qualities[q]>=20)
        for run in runs:
            depths=[0]*(end-start)
            for bases in positions[run].values():
                for pos in bases:depths[pos-start]+=1
            key=dict(run_accession=run,species=samples[run]['scientific_name'],location=samples[run]['location'],assembly=assembly)
            for category,coords in categories.items():
                values=[depths[pos-start] for pos in sorted(coords)]
                summary.append(dict(key,region=category,reference_bases=len(coords),mean_fragment_depth=f'{mean(values):.4f}',
                    median_fragment_depth=f'{median(values):.2f}',fraction_bases_at_depth_10=f'{sum(d>=10 for d in values)/len(values):.6f}',
                    mean_fragment_depth_per_input_Gb=f"{mean(values)/(int(samples[run]['base_count'])/1e9):.6f}"))
            for pos in range(start,end,100):
                values=depths[pos-start:min(end,pos+100)-start]
                bins.append(dict(key,start_1based=pos+1,end_1based=min(end,pos+100),mean_fragment_depth=f'{mean(values):.4f}'))
    write(a.prefix+'locus_coverage_summary.tsv',summary)
    write(a.prefix+'locus_coverage_100bp.tsv',bins)
    print('Measured exon, intron and flank coverage for',len(runs),'workers in both mappings.')

if __name__=='__main__':main()
