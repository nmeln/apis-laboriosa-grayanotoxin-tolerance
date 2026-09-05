#!/usr/bin/env python3
"""Rebuild all addendum outputs locally from the checked source snapshot."""
import gzip
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO

HERE=Path(__file__).resolve().parents[1]
ROOT=HERE.parent
WORK=HERE/'work'
OUT=HERE/'results'
MM=ROOT/'.tools/minimap2-2.31_x64-linux/minimap2'

def script(name):
    print('Running '+name,flush=True)
    subprocess.run([sys.executable,str(HERE/'scripts'/name)],check=True)

def align(arguments, destination, log):
    print('Aligning '+destination.name,flush=True)
    with destination.open('w') as out, log.open('w') as err:
        subprocess.run([str(MM),*map(str,arguments)],stdout=out,stderr=err,check=True)

def main():
    WORK.mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
    script('verify.py')
    script('audit_sources.py')
    script('audit_raw_subsets.py')
    reference={'AL':ROOT/'genomes/apis_laboriosa/refseq/GCF_014066325.1_rna.fna.gz',
               'AD':ROOT/'genomes/apis_dorsata/GCF_000469605.1_rna.fna.gz'}
    unigenes={'AL':ROOT/'references/transcriptome_2019/GSM3757258_AL.unigene.fasta.gz',
              'AD':ROOT/'references/transcriptome_2019/GSM3757259_AD.unigene.fasta.gz'}
    assert subprocess.check_output([str(MM),'--version'],text=True).strip()=='2.31-r1302'
    for query in ['AL','AD']:
        for ref in ['AL','AD']:
            label=query+'_to_'+ref
            align(['-x','asm5','-c','--secondary=yes','-N','50','-p','0.5','-t','4',reference[ref],unigenes[query]],
                  WORK/(label+'.paf'),WORK/(label+'.log'))
    script('analyze_transcripts.py')
    script('prepare_raw_reference.py')
    genomes=[SeqIO.read(HERE/'inputs'/name,'genbank') for name in ['OR496406.gb','KY741959.gb']]
    assert [r.id for r in genomes]==['OR496406.1','KY741959.1']
    assert all(len(r)==8440 for r in genomes)
    SeqIO.write(genomes,WORK/'BQCV_references.fna','fasta')
    with gzip.open(unigenes['AD'],'rt') as f:
        contig=next(r for r in SeqIO.parse(f,'fasta') if r.id=='AD|c34801_g1')
    assert len(contig)==8446
    # Prefix to retain the exact identifier used during the discovery alignment.
    contig.id='AD|c34801_g1';contig.description=''
    SeqIO.write(contig,WORK/'AD_c34801_g1.fna','fasta')
    align(['-x','asm5','-c','--secondary=yes','-N','5','-t','2',WORK/'BQCV_references.fna',WORK/'AD_c34801_g1.fna'],
          OUT/'bqcv_contig_alignment.paf',WORK/'bqcv_contig.log')
    for sp,run in [('AL','SRR9034695'),('AD','SRR9034696')]:
        reads=HERE/'inputs'/f'{run}.first_2000000.R1.fastq.gz'
        for name,ref in [('raw','pooled_primary_cds.fna'),('bqcv','BQCV_references.fna')]:
            align(['-x','sr','-c','--secondary=yes','-N','20' if name=='raw' else '5','-p','0.8','-t','4' if name=='raw' else '2',WORK/ref,reads],
                  WORK/f'{sp}_{name}.paf',WORK/f'{sp}_{name}.log')
    script('analyze_raw_reads.py')
    script('analyze_virus.py')
    script('classify_chitin_synthases.py')
    script('validate_claims.py')

if __name__=='__main__':main()
