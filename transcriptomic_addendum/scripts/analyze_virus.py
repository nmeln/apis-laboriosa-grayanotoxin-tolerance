#!/usr/bin/env python3
"""Local BQCV identity and read-coverage check against two pinned genomes."""
import csv,itertools,json,statistics
from collections import Counter,defaultdict
from analyze_transcripts import HERE,OUT,write

def parsed(line):
    p=line.rstrip().split('\t')
    tags={s.split(':',2)[0]:s.split(':',2)[2] for s in p[12:]}
    return {'query':p[0],'query_length':int(p[1]),'query_start':int(p[2]),'query_end':int(p[3]),
            'reference':p[5],'reference_length':int(p[6]),'start':int(p[7]),'end':int(p[8]),
            'matches':int(p[9]),'alignment_length':int(p[10]),'score':int(tags['AS'])}

def passes(h):
    return h['alignment_length']>=100 and h['matches']/h['alignment_length']>=.9 and (h['query_end']-h['query_start'])/h['query_length']>=.8

def main():
    identity=[]
    with (OUT/'bqcv_contig_alignment.paf').open() as f:
        for line in f:
            h=parsed(line)
            identity.append({**h,'identity':h['matches']/h['alignment_length'],
                             'query_coverage':(h['query_end']-h['query_start'])/h['query_length']})
    write('bqcv_contig_identity.tsv',identity)
    rows,bins=[],[]
    for sp in ['AL','AD']:
        count=0;read_ids=set();depth=[0]*8440;per_reference=Counter();scores={}
        # Use one reference for coverage so equivalent matches are not counted twice.
        with (HERE/'work'/f'{sp}_bqcv.paf').open() as f:
            for q,lines in itertools.groupby(f,key=lambda x:x.split('\t',1)[0]):
                hits=[h for h in map(parsed,lines) if passes(h)]
                if not hits:continue
                count+=1;read_ids.add(q);scores[q]=max(h['score'] for h in hits)
                for ref in {h['reference'] for h in hits}:per_reference[ref]+=1
                on_primary=[h for h in hits if h['reference']=='OR496406.1']
                if on_primary:
                    h=max(on_primary,key=lambda x:x['score'])
                    depth[h['start']]+=1
                    if h['end']<len(depth):depth[h['end']]-=1
        # Check whether these same read IDs have credible competing bee CDS alignments.
        host_competitors=set()
        with (HERE/'work'/f'{sp}_raw.paf').open() as f:
            for line in f:
                q=line.split('\t',1)[0]
                if q in read_ids:
                    h=parsed(line)
                    if passes(h) and h['score']>=.98*scores[q]:host_competitors.add(q)
        d=0
        for i in range(len(depth)):d+=depth[i];depth[i]=d
        for start in range(0,len(depth),1000):
            ds=depth[start:start+1000]
            bins.append({'sample':sp,'reference':'OR496406.1','start_zero_based':start,'end_exclusive':start+len(ds),
                         'mean_alignment_span_depth':sum(ds)/len(ds),'bases_with_alignment_span_coverage':sum(x>0 for x in ds)})
        rows.append({'sample':sp,'input_R1_reads':2_000_000,'BQCV_matching_reads':count,'fraction_of_input_R1':count/2_000_000,
                     'OR496406_matching_reads':per_reference['OR496406.1'],'KY741959_matching_reads':per_reference['KY741959.1'],
                     'credible_competing_bee_CDS_reads':len(host_competitors),
                     'OR496406_bases_spanned':sum(x>0 for x in depth),'OR496406_reference_length':len(depth),
                     'OR496406_mean_alignment_span_depth':sum(depth)/len(depth)})
    write('bqcv_raw_summary.tsv',rows);write('bqcv_coverage_bins.tsv',bins)
    print(json.dumps(rows,indent=2))

if __name__=='__main__': main()
