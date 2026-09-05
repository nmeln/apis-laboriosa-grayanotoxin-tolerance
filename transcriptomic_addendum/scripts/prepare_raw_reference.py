#!/usr/bin/env python3
"""Create a symmetric pooled CDS reference for R1 technical validation.

Each species supplies one primary protein per gene from the frozen comparative
analysis. Locate its complete amino-acid sequence in its annotated RNA. Fail or
record exclusions when no exact translated CDS is recoverable. This avoids
unequal UTR lengths and most reference isoform multiplicity.
"""
from pathlib import Path
import csv,gzip,json
from Bio import SeqIO
from Bio.Seq import Seq
from analyze_transcripts import ROOT,HERE,PATHS,SPECIES,refs,attrs,table

def main():
    tx,descriptions,one_to_one,groups=refs()
    members=table(ROOT/'comparative_addendum/results/primary_proteome_members_current.tsv')
    protein_files={'AL':'genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz',
                   'AD':'genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz'}
    rna_files={'AL':'genomes/apis_laboriosa/refseq/GCF_014066325.1_rna.fna.gz',
               'AD':'genomes/apis_dorsata/GCF_000469605.1_rna.fna.gz'}
    rows,exclusions=[],[]
    output=HERE/'work/pooled_primary_cds.fna'
    with output.open('w') as out:
        for code,sp in SPECIES.items():
            selected={r['accession']:r for r in members if r['species']==sp}
            protein_to_transcripts={}
            # A CDS Parent refers to mRNA ID; collect all possible source RNAs.
            mrna={};links=[]
            with gzip.open(ROOT/PATHS[code][0],'rt') as f:
                for line in f:
                    if line.startswith('#'):continue
                    p=line.rstrip().split('\t')
                    if len(p)!=9:continue
                    a=attrs(p[8])
                    if p[2]=='mRNA' and a.get('transcript_id'):mrna[a['ID']]=a['transcript_id']
                    if p[2]=='CDS' and a.get('protein_id') in selected:
                        links.append((a['protein_id'],a.get('Parent','')))
            for protein,parent in links:
                if parent in mrna:protein_to_transcripts.setdefault(protein,set()).add(mrna[parent])
            with gzip.open(ROOT/protein_files[code],'rt') as f:
                proteins={r.id:str(r.seq).rstrip('*') for r in SeqIO.parse(f,'fasta') if r.id in selected}
            needed=set().union(*protein_to_transcripts.values())
            with gzip.open(ROOT/rna_files[code],'rt') as f:
                rnas={r.id:str(r.seq) for r in SeqIO.parse(f,'fasta') if r.id in needed}
            for protein,m in sorted(selected.items()):
                candidates=[]
                for tid in sorted(protein_to_transcripts.get(protein,())):
                    if tid not in rnas:continue
                    s=rnas[tid]
                    for frame in range(3):
                        t=str(Seq(s[frame:len(s)-(len(s)-frame)%3]).translate())
                        pos=t.find(proteins[protein])
                        if pos>=0:
                            start=frame+3*pos;cds=s[start:start+len(proteins[protein])*3]
                            candidates.append((tid,start,cds))
                unique={cds for _,_,cds in candidates}
                if not candidates or len(unique)!=1:
                    exclusions.append({'species':code,'protein':protein,'gene':m['gene'],'reason':'no_exact_CDS' if not candidates else 'multiple_CDS_sequences'})
                    continue
                tid,start,cds=candidates[0]
                rid=code+'|'+protein
                out.write('>'+rid+'\n'+cds+'\n')
                rows.append({'reference_id':rid,'species':code,'protein':protein,'gene':m['gene'],
                             'group':tx[code].get(tid,{}).get('group',''), 'transcript':tid,
                             'cds_start_zero_based':start,'cds_length':len(cds),'product':m['description']})
    for name,items in [('raw_reference_metadata.tsv',rows),('raw_reference_exclusions.tsv',exclusions)]:
        with (HERE/'results'/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(items[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(items)
    print(json.dumps({'reference_sequences':len(rows),'exclusions':len(exclusions),'output':str(output)},indent=2))

if __name__=='__main__':main()
