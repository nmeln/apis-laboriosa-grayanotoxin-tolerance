#!/usr/bin/env python3
"""Count unique orthogroup alignments in fixed R1 prefixes, without DE claims."""
import csv,itertools,json,math,statistics
from collections import Counter,defaultdict
from analyze_transcripts import HERE,OUT,table,write

def count_file(path,meta):
    counts=Counter();qc=Counter(); candidates=[]
    with path.open() as f:
        for query,lines in itertools.groupby(f,key=lambda s:s.split('\t',1)[0]):
            hits={};qc['reads_with_reported_alignment']+=1
            for line in lines:
                p=line.rstrip().split('\t')
                m=meta[p[5]]
                if int(p[10])<100 or int(p[9])/int(p[10])<.9 or (int(p[3])-int(p[2]))/int(p[1])<.8:continue
                tags={s.split(':',2)[0]:s.split(':',2)[2] for s in p[12:]}
                key=m['group'] or 'UNGROUPED:'+m['gene']
                hit=(int(tags['AS']),p[5],int(p[7]),int(p[8]),int(p[9]),int(p[10]))
                if key not in hits or hit>hits[key]:hits[key]=hit
            if not hits:qc['failed_alignment_filters']+=1;continue
            best=sorted(hits,key=lambda k:(-hits[k][0],k))
            if len(best)>1 and hits[best[1]][0]>=.98*hits[best[0]][0]:
                qc['ambiguous_group']+=1;continue
            group=best[0];counts[group]+=1;qc['unique_group_reads']+=1
            if group=='OG0001109':
                score,ref,start,end,match,alen=hits[group]
                candidates.append({'read_id':query,'orthogroup':group,'best_reference':ref,'score':score,
                                   'target_start':start,'target_end':end,'matched_bases':match,'aligned_bases':alen})
    return counts,qc,candidates

def main():
    meta={r['reference_id']:r for r in table(OUT/'raw_reference_metadata.tsv')}
    ann={r['orthogroup']:r for r in table(OUT/'orthogroup_annotations.tsv')}
    counts,quality,candidates={},{},[]
    for sp in ['AL','AD']:
        counts[sp],qc,ca=count_file(HERE/'work'/f'{sp}_raw.paf',meta)
        quality[sp]=dict(qc)
        candidates.extend({'sample':sp,**r} for r in ca)
    eligible=sorted(g for g in ann if ann[g]['one_to_one_AL_AD']=='True' and min(counts['AL'][g],counts['AD'][g])>=10)
    center=statistics.median(math.log2(counts['AL'][g]/counts['AD'][g]) for g in eligible)
    rows=[]
    for g in sorted(set(counts['AL'])|set(counts['AD'])):
        a,b=counts['AL'][g],counts['AD'][g]
        rows.append({'orthogroup':g,'AL_reads':a,'AD_reads':b,'AL_per_million_input_R1':a/2,'AD_per_million_input_R1':b/2,
                     'raw_ratio':a/b if b else None,'median_normalized_ratio':a/b/2**center if b else None,
                     'panels':ann.get(g,{}).get('panels',''),'descriptions':ann.get(g,{}).get('descriptions','')})
    write('raw_orthogroup_counts.tsv',rows)
    write('peritrophin_raw_assignments.tsv',candidates)
    panels=sorted({p for r in ann.values() for p in r['panels'].split(';') if p})
    panelrows=[]
    for panel in panels:
        gs={g for g in ann if panel in ann[g]['panels'].split(';')}
        av=sum(counts['AL'][g] for g in gs);bv=sum(counts['AD'][g] for g in gs)
        paired=[g for g in eligible if g in gs]
        panelrows.append({'panel':panel,'AL_reads':av,'AD_reads':bv,'raw_ratio':av/bv if bv else None,
                          'median_normalized_sum_ratio':av/bv/2**center if bv else None,
                          'adequate_one_to_one_groups':len(paired),
                          'median_normalized_gene_ratio':2**statistics.median(math.log2(counts['AL'][g]/counts['AD'][g])-center for g in paired) if paired else None})
    write('raw_panel_counts.tsv',panelrows)
    summary={'input_R1_reads_each':2_000_000,'normalization_background_groups':len(eligible),
             'median_log2_AL_AD_ratio':center,'median_AL_AD_ratio':2**center,'alignment_qc':quality,
             'peritrophin_OG0001109':next(r for r in rows if r['orthogroup']=='OG0001109'),
             'limit':'Two unreplicated pools, belly excluded according to GEO; deterministic R1 prefixes; technical counts only.'}
    (OUT/'raw_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
