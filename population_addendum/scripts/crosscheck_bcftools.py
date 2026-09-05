"""Independent SNP-level BCFtools cross-check of every called codon base.

Calls every covered target base, including reference positions. Missing records
remain missing. Each individual is called separately, without pooling species
allele-frequency priors. This check cannot phase multiple variable codon bases.
"""
import argparse
import json
from pathlib import Path
import tempfile
import pysam
import pysam.bcftools as bcftools
from call_codons import COMPLEMENT
from common import OUT, WORK, table, write, json_write


def main():
    p=argparse.ArgumentParser();p.add_argument('--label',default='population');p.add_argument('--prefix',default='');a=p.parse_args()
    destination=WORK/a.label
    runs=json.loads((destination/'mapping_inputs.json').read_text())['runs']
    target_rows=table(OUT/'target_coordinates.tsv')
    codon_calls={(r['run_accession'],r['assembly'],r['focal_AL_position']):r for r in table(OUT/(a.prefix+'sample_codon_calls.tsv')) if r['profile']=='primary'}
    results=[]
    for assembly in ('AL_Shangrila','AD_Malaysia'):
        targets=[r for r in target_rows if r['assembly']==assembly]
        positions=sorted({(r['scaffold'],int(pos)) for r in targets for pos in r['genomic_codon_positions_1based'].split(',')})
        regions=destination/f'{assembly}.bcftools_positions.tsv'
        regions.write_text(''.join(f'{chrom}\t{pos}\t{pos}\n' for chrom,pos in positions))
        fasta=Path(tempfile.gettempdir())/'bee_population_references'/(assembly+'.fna')
        if not Path(str(fasta)+'.fai').exists():pysam.faidx(str(fasta))
        pileup=destination/f'{assembly}.pileup.bcf'
        bcftools.mpileup('-f',str(fasta),'-R',str(regions),'-q','30','-Q','30','-d','10000',
                         '--ff','UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY','-a','FORMAT/DP,FORMAT/AD','-Ob','-o',str(pileup),str(destination/(assembly+'.bam')),catch_stdout=False)
        for run in runs:
            subset=destination/f'{assembly}.{run}.pileup.bcf'
            vcf=destination/f'{assembly}.{run}.vcf'
            bcftools.view('-s',run,'-Ob','-o',str(subset),str(pileup),catch_stdout=False)
            bcftools.call('-m','--ploidy','2','-f','GQ','-Ov','--no-version','-o',str(vcf),str(subset),catch_stdout=False)
            records={}
            with pysam.VariantFile(str(vcf)) as variants:
                for record in variants:
                    sample=record.samples[run];gt=sample.get('GT')
                    bases='/'.join(sorted(record.alleles[i] for i in gt)) if gt and None not in gt else ''
                    records[record.chrom,record.pos]=dict(genotype=bases,depth=sample.get('DP',0),gq=sample.get('GQ',''),alleles=','.join(record.alleles),AD=','.join(map(str,sample.get('AD',()))))
            for target in targets:
                c=codon_calls[run,assembly,target['focal_AL_position']]
                codons=c['best_codons'].split('/') if c['best_codons'] else []
                for index,pos in enumerate(map(int,target['genomic_codon_positions_1based'].split(','))):
                    expected='/'.join(sorted(codon[index].translate(COMPLEMENT) if target['strand']=='-' else codon[index] for codon in codons))
                    obs=records.get((target['scaffold'],pos),dict(genotype='',depth=0,gq='',alleles='',AD=''))
                    status='codon_call_filtered' if c['status']!='PASS' else ('missing_bcftools' if not obs['genotype'] else ('agree' if obs['genotype']==expected else 'DISAGREE'))
                    results.append(dict(run_accession=run,assembly=assembly,focal_AL_position=target['focal_AL_position'],role=target['role'],
                                        coding_base=index+1,genomic_position=pos,codon_call_status=c['status'],expected_genomic_genotype=expected,
                                        bcftools_genotype=obs['genotype'],bcftools_depth=obs['depth'],bcftools_gq=obs['gq'],bcftools_alleles=obs['alleles'],bcftools_allele_depths=obs['AD'],comparison=status))
            subset.unlink()
    write(a.prefix+'bcftools_base_crosscheck.tsv',results)
    from collections import Counter
    json_write(OUT/(a.prefix+'bcftools_crosscheck_summary.json'),dict(bcftools_version=pysam.__samtools_version__,
         comparisons=dict(Counter(r['comparison'] for r in results)),
         primary_site_comparisons=dict(Counter(r['comparison'] for r in results if r['role']=='primary')),
         note='Independent SNP calls use BCFtools default BAQ and overlap handling; primary results remain based on fragment codon calls.'))
    print('BCFtools base comparisons:',dict(Counter(r['comparison'] for r in results)))

if __name__=='__main__':main()
