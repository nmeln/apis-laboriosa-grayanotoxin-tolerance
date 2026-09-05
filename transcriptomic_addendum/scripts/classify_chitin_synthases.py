#!/usr/bin/env python3
"""Resolve the two bee CHS classes against functionally distinguished fly genes."""
import csv,gzip,json
from Bio import Align,SeqIO
from Bio.Align import substitution_matrices
from analyze_transcripts import HERE,ROOT,OUT,table,write

def main():
    groups=[r for r in table(ROOT/'comparative_addendum/results/orthofinder_key/Orthogroups.tsv') if r['Orthogroup'] in ['OG0000293','OG0000875']]
    specs={'Apis_laboriosa':'genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz',
           'Apis_dorsata':'genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz',
           'Apis_mellifera':'genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz',
           'Apis_cerana':'genomes/apis_cerana/GCF_029169275.2_protein.faa.gz',
           'Apis_florea':'genomes/apis_florea/GCF_048593485.1_protein.faa.gz'}
    wanted={p:(sp,r['Orthogroup']) for r in groups for sp in specs for p in r[sp].split(', ') if p}
    proteins={}
    for sp,path in specs.items():
        with gzip.open(ROOT/path,'rt') as f:
            for r in SeqIO.parse(f,'fasta'):
                if r.id in wanted:proteins[r.id]=r
    assert len(proteins)==len(wanted)
    # Keep every long reference isoform, not just the best-looking one.
    references=[]
    for gene in ['kkv','Chs2']:
        for r in SeqIO.parse(HERE/'inputs'/f'Dmel_{gene}_uniprot.faa','fasta'):
            if len(r)>=1300:references.append((gene,r))
    assert len(references)==4
    aligner=Align.PairwiseAligner(mode='local',substitution_matrix=substitution_matrices.load('BLOSUM62'),open_gap_score=-10,extend_gap_score=-.5)
    rows=[];best={}
    with (OUT/'chitin_synthase_pairwise_alignments.faa').open('w') as out:
        for acc,r in sorted(proteins.items()):
            if len(r)<1000:continue
            for gene,ref in references:
                a=aligner.align(str(r.seq),str(ref.seq))[0]
                left,right=a[0],a[1]
                pairs=[(x,y) for x,y in zip(left,right) if x!='-' and y!='-']
                identical=sum(x==y for x,y in pairs)
                rows.append({'bee_species':wanted[acc][0],'orthogroup':wanted[acc][1],'bee_protein':acc,'bee_length':len(r),
                             'fly_gene':gene,'fly_class':'CHS-A_cuticle_trachea' if gene=='kkv' else 'CHS-B_midgut_PM',
                             'fly_accession':ref.id.split('|')[1],'fly_length':len(ref),'local_alignment_score':a.score,
                             'aligned_residue_pairs':len(pairs),'identical_residue_pairs':identical,'identity':identical/len(pairs),
                             'bee_paired_residue_coverage':len(pairs)/len(r),'fly_paired_residue_coverage':len(pairs)/len(ref)})
                out.write(f'>{acc}_to_{ref.id.split("|")[1]}_bee\n{left}\n>{acc}_to_{ref.id.split("|")[1]}_fly\n{right}\n')
    write('chitin_synthase_classification.tsv',rows)
    summary=[]
    raw={r['orthogroup']:r for r in table(OUT/'raw_orthogroup_counts.tsv')}
    for acc in sorted({r['bee_protein'] for r in rows}):
        rr=[r for r in rows if r['bee_protein']==acc]
        a=max((r for r in rr if r['fly_gene']=='kkv'),key=lambda r:r['local_alignment_score'])
        b=max((r for r in rr if r['fly_gene']=='Chs2'),key=lambda r:r['local_alignment_score'])
        winner=a if a['local_alignment_score']>b['local_alignment_score'] else b
        scores_a=[r['local_alignment_score'] for r in rr if r['fly_gene']=='kkv']
        scores_b=[r['local_alignment_score'] for r in rr if r['fly_gene']=='Chs2']
        og=winner['orthogroup']
        summary.append({'bee_species':winner['bee_species'],'orthogroup':og,'bee_protein':acc,
                        'closest_class':winner['fly_class'],'all_long_isoforms_agree':min(scores_a)>max(scores_b) or min(scores_b)>max(scores_a),
                        'best_Kkv_score':a['local_alignment_score'],'best_Chs2_score':b['local_alignment_score'],
                        'Kkv_identity':a['identity'],'Chs2_identity':b['identity'],'AL_raw_reads_in_group':raw[og]['AL_reads'],
                        'AD_raw_reads_in_group':raw[og]['AD_reads'],'AL_AD_group_median_normalized_ratio':raw[og]['median_normalized_ratio']})
    write('chitin_synthase_class_summary.tsv',summary)
    print(json.dumps([r for r in summary if r['bee_species']=='Apis_laboriosa'],indent=2))

if __name__=='__main__':main()
