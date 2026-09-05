#!/usr/bin/env python3
"""Audit descriptive abundance using full alignments and reference rotation.

No inferential species-expression statistics: one biological pool per species.
All inputs are local and pinned by the parent repository.
"""
from __future__ import annotations
import csv
import gzip
import itertools
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
OUT = HERE / 'results'
WORK = HERE / 'work'
SPECIES = {'AL': 'Apis_laboriosa', 'AD': 'Apis_dorsata'}
PATHS = {
    'AL': ('genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz', 'GSM3757258_AL'),
    'AD': ('genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz', 'GSM3757259_AD'),
}

def table(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(name, rows, columns=None):
    rows = list(rows)
    assert rows or columns, name
    with (OUT / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=columns or list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

def fasta(path):
    with gzip.open(path, 'rt') as f:
        name, chunks = None, []
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, ''.join(chunks)

def attrs(text):
    return {k: unquote(v) for k,v in (x.split('=',1) for x in text.split(';') if '=' in x)}

def refs():
    groups = table(ROOT/'comparative_addendum/results/orthofinder_key/Orthogroups.tsv')
    members = table(ROOT/'comparative_addendum/results/primary_proteome_members_current.tsv')
    protein_group = {(sp, acc): r['Orthogroup'] for r in groups for sp in SPECIES.values()
                     for acc in r[sp].split(', ') if acc}
    gene_group = {(r['species'], r['gene']): protein_group[(r['species'],r['accession'])]
                  for r in members if (r['species'],r['accession']) in protein_group}
    descriptions = defaultdict(set)
    # Use all six primary descriptions for a symmetric functional annotation.
    all_protein_group = {(sp,acc):r['Orthogroup'] for r in groups for sp in r if sp != 'Orthogroup'
                         for acc in r[sp].split(', ') if acc}
    for r in members:
        og = all_protein_group.get((r['species'],r['accession']))
        if og:
            descriptions[og].add(r['description'])
    transcripts = {}
    for code, (path, _) in PATHS.items():
        genes, tx = {}, {}
        with gzip.open(ROOT/path, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                p=line.rstrip().split('\t')
                if len(p)!=9: continue
                a=attrs(p[8])
                if p[2]=='gene':
                    genes[a['ID']]=a.get('Name',a.get('gene',a['ID'].removeprefix('gene-')))
                elif p[2] in ('mRNA','transcript') and a.get('transcript_id'):
                    tx[a['transcript_id']]=(a.get('Parent',''),a.get('product',''))
        transcripts[code]={tid: {'gene':genes.get(parent,parent.removeprefix('gene-')),
                                 'group':gene_group.get((SPECIES[code],genes.get(parent,parent.removeprefix('gene-'))),''),
                                 'product':product} for tid,(parent,product) in tx.items()}
    one_to_one={r['Orthogroup'] for r in groups if all(r[sp] and ', ' not in r[sp] for sp in SPECIES.values())}
    return transcripts, descriptions, one_to_one, groups

def expressions():
    result, qc = {}, []
    for code, (_,prefix) in PATHS.items():
        base=ROOT/'references/transcriptome_2019'
        lengths={n:len(s) for n,s in fasta(base/(prefix+'.unigene.fasta.gz'))}
        with gzip.open(base/(prefix+'.Readcount_FPKM.txt.gz'),'rt') as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        assert len(rows)==len(lengths)
        counts=sum(float(r['Read_count']) for r in rows)
        fpkm=sum(float(r['FPKM']) for r in rows)
        rates=sum(float(r['Read_count'])/lengths[r['gene_id']] for r in rows)
        result[code]={r['gene_id']:{'count':float(r['Read_count']), 'fpkm':float(r['FPKM']),
                       'tpm_submitted':float(r['FPKM'])/fpkm*1e6,
                       'tpm_length':float(r['Read_count'])/lengths[r['gene_id']]/rates*1e6,
                       'length':lengths[r['gene_id']]} for r in rows}
        highest=max(result[code], key=lambda x: result[code][x]['count'])
        qc.append({'species':code,'unigenes':len(rows),'submitted_count_sum':counts,'submitted_fpkm_sum':fpkm,
                   'count_per_nt_sum':rates,'largest_count_unigene':highest,
                   'largest_count_fraction':result[code][highest]['count']/counts})
    write('library_summary.tsv',qc)
    return result

def assignments(path, meta, min_coverage):
    """Best score per orthogroup after filtering each nucleotide alignment."""
    result = {}
    with path.open() as f:
        for name, lines in itertools.groupby(f, key=lambda line:line.split('\t',1)[0]):
            hits={}
            for line in lines:
                p=line.rstrip().split('\t')
                if p[5] not in meta: continue
                identity=int(p[9])/int(p[10])
                coverage=(int(p[3])-int(p[2]))/int(p[1])
                if int(p[10])<200 or identity<.95 or coverage<min_coverage: continue
                m=meta[p[5]]
                # Unassigned genes still compete with annotated orthogroups.
                key=m['group'] or 'UNGROUPED:'+m['gene']
                tags={s.split(':',2)[0]:s.split(':',2)[2] for s in p[12:]}
                score=int(tags['AS'])
                hit={'group':key, 'gene':m['gene'], 'transcript':p[5], 'product':m['product'],
                     'score':score,'identity':identity,'query_coverage':coverage,
                     'aligned_nt':int(p[10]),'query_start':int(p[2]),'query_end':int(p[3]),
                     'target_start':int(p[7]),'target_end':int(p[8]),'target_length':int(p[6])}
                if key not in hits or (score,coverage,p[5])>(hits[key]['score'],hits[key]['query_coverage'],hits[key]['transcript']):
                    hits[key]=hit
            if not hits: continue
            best=sorted(hits.values(),key=lambda h:(-h['score'],-h['query_coverage'],h['transcript']))
            h=best[0]
            h['runner_up_score']=best[1]['score'] if len(best)>1 else 0
            h['unique_group']=h['runner_up_score']<.98*h['score']
            result[name]=h
    return result

def summed(assign, expr):
    groups=defaultdict(lambda: {'count':0.,'fpkm':0.,'tpm_submitted':0.,'tpm_length':0.,'unigenes':[]})
    for q,h in assign.items():
        if not h['unique_group'] or h['group'].startswith('UNGROUPED:'): continue
        g=groups[h['group']]
        for k in ('count','fpkm','tpm_submitted','tpm_length'): g[k]+=expr[q][k]
        g['unigenes'].append(q)
    return groups

def ratio(a,b): return a/b if b else None

def main():
    OUT.mkdir(exist_ok=True,parents=True)
    tx, desc, one_to_one, groups=refs()
    expr=expressions()
    patterns={r['panel']:re.compile(r['pattern'],re.I) for r in table(HERE/'panels.tsv')}
    panels={p:{g for g,d in desc.items() if pat.search('; '.join(sorted(d)))} for p,pat in patterns.items()}
    annotation=[{'orthogroup':g,'one_to_one_AL_AD':g in one_to_one,
                 'panels':';'.join(p for p,gs in panels.items() if g in gs),
                 'descriptions':'; '.join(sorted(desc[g]))} for g in sorted(desc)]
    write('orthogroup_annotations.tsv',annotation)
    ann={r['orthogroup']:r for r in annotation}
    assign, sums, mapqc, tops=[],{},[],[]
    for coverage in (.8,.5):
        for query in SPECIES:
            for reference in SPECIES:
                label=f'{query}_to_{reference}'
                a=assignments(WORK/(label+'.paf'),tx[reference],coverage)
                sums[(coverage,query,reference)]=summed(a,expr[query])
                unique=[q for q,h in a.items() if h['unique_group'] and not h['group'].startswith('UNGROUPED:')]
                mapqc.append({'coverage_threshold':coverage,'query':query,'reference':reference,
                              'unique_group_unigenes':len(unique),'mapped_count_fraction':sum(expr[query][q]['count'] for q in unique)/sum(e['count'] for e in expr[query].values())})
                if coverage==.8:
                    assign.extend({'query_species':query,'reference_species':reference,'unigene':q,**h,**expr[query][q]} for q,h in sorted(a.items()))
                    for q in sorted(expr[query],key=lambda q:expr[query][q]['count'],reverse=True)[:30]:
                        tops.append({'query_species':query,'reference_species':reference,'unigene':q,**expr[query][q],
                                     'group':a.get(q,{}).get('group',''), 'product':a.get(q,{}).get('product',''),
                                     'unique_group':a.get(q,{}).get('unique_group','')})
    write('unigene_assignments.tsv',assign)
    write('mapping_summary.tsv',mapqc)
    write('top_count_transcripts.tsv',tops)
    comparisons,panelrows,normalization=[],[],[]
    for coverage in (.8,.5):
        for mode,alref,adref in [('own','AL','AD'),('both_AL','AL','AL'),('both_AD','AD','AD')]:
            al=sums[(coverage,'AL',alref)]; ad=sums[(coverage,'AD',adref)]
            background=sorted(g for g in one_to_one if g in al and g in ad and al[g]['count']>=10 and ad[g]['count']>=10)
            assert len(background)>1000, (coverage,mode,len(background))
            for measure in ('fpkm','tpm_submitted','tpm_length','count'):
                valid=[g for g in background if al[g][measure]>0 and ad[g][measure]>0]
                center=statistics.median(math.log2(al[g][measure]/ad[g][measure]) for g in valid)
                normalization.append({'coverage_threshold':coverage,'reference_mode':mode,'measure':measure,'background_groups':len(valid),'median_log2_ratio':center})
                for g in sorted(set(al)|set(ad)):
                    av=al.get(g,{}).get(measure,0.); bv=ad.get(g,{}).get(measure,0.)
                    raw=ratio(av,bv)
                    if measure=='fpkm' and coverage==.8:
                        comparisons.append({'reference_mode':mode,'orthogroup':g,'one_to_one_AL_AD':g in one_to_one,'AL_value':av,'AD_value':bv,
                                            'raw_ratio':raw,'median_normalized_ratio':raw/(2**center) if raw is not None else None,
                                            'AL_unigenes':';'.join(sorted(al.get(g,{}).get('unigenes',[]))),
                                            'AD_unigenes':';'.join(sorted(ad.get(g,{}).get('unigenes',[]))),
                                            'panels':ann[g]['panels'],'descriptions':ann[g]['descriptions']})
                for p,gs in panels.items():
                    # Family abundance includes grouped paralogues. Per-group median uses one-to-one genes with adequate counts.
                    av=sum(al.get(g,{}).get(measure,0.) for g in sorted(gs));bv=sum(ad.get(g,{}).get(measure,0.) for g in sorted(gs))
                    raw=ratio(av,bv)
                    paired=[g for g in valid if g in gs]
                    logs=[math.log2(al[g][measure]/ad[g][measure])-center for g in paired]
                    panelrows.append({'coverage_threshold':coverage,'reference_mode':mode,'measure':measure,'panel':p,
                                      'annotated_orthogroups':len(gs),'AL_mapped_groups':sum(g in al for g in gs),'AD_mapped_groups':sum(g in ad for g in gs),
                                      'AL_total':av,'AD_total':bv,'raw_sum_ratio':raw,'median_normalized_sum_ratio':raw/2**center if raw is not None else None,
                                      'paired_one_to_one_groups':len(paired),'median_normalized_gene_ratio':2**statistics.median(logs) if logs else None})
    write('orthogroup_expression.tsv',comparisons)
    write('panel_sensitivity.tsv',panelrows)
    write('normalization.tsv',normalization)
    print(json.dumps({'mapped_assignments':len(assign),'comparisons':len(comparisons),'panel_rows':len(panelrows),'normalization':normalization[:4]},indent=2))

if __name__=='__main__': main()
