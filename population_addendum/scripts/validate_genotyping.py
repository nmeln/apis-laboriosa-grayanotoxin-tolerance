"""Known-truth controls for codons, mate overlap and excluded evidence."""
from collections import defaultdict
from pathlib import Path
import tempfile
import pysam
from Bio.Seq import Seq
from call_codons import genotype, collapse_fragments, read_observation, call_bam
from common import OUT, WORK, json_write


def main():
    cases = []
    def observations(codons):
        return [dict(codon=c,qualities=[40,40,40],reverse=bool(i%2),read1=True) for i,c in enumerate(codons)]
    for name,codons,expected,status in [
        ('homozygous_L', ['CTT']*12, 'CTT/CTT','PASS'),
        ('homozygous_F', ['TTT']*12, 'TTT/TTT','PASS'),
        ('balanced_LF', ['CTT']*6+['TTT']*6, 'CTT/TTT','PASS'),
        ('synonymous_LL', ['CTT']*6+['TTG']*6, 'CTT/TTG','PASS'),
        ('unexpected_W', ['TGG']*12, 'TGG/TGG','PASS'),
        ('low_depth', ['CTT']*5, 'CTT/CTT','low_depth'),
        ('absent', [], '', 'no_coverage')]:
        result=genotype(observations(codons),10)
        assert result['best_codons']==expected and status in result['status'],(name,result)
        cases.append(dict(test=name,passed=True))
    one=genotype([dict(codon='CTT',qualities=[40]*3,reverse=False)]*12,10)
    assert one['status']=='one_strand_only'
    cases.append(dict(test='strand_flag',passed=True))
    directory=WORK/'synthetic_caller';directory.mkdir(parents=True,exist_ok=True)
    raw=directory/'known.bam'; sorted_bam=directory/'known.sorted.bam'
    target=dict(assembly='synthetic',scaffold='chr1',strand='-',focal_AL_position='254',role='primary',genomic_codon_positions_1based='103,102,101')
    header=dict(HD={'VN':'1.6'},SQ=[dict(SN='chr1',LN=1000)],RG=[dict(ID='synthetic',SM='synthetic')])
    with pysam.AlignmentFile(str(raw),'wb',header=header) as out:
        for i in range(12):
            for mate in (1,2):
                read=pysam.AlignedSegment(out.header);read.query_name=f'fragment{i}'
                codon='CTT' if i<6 else 'TTT';genomic=str(Seq(codon).reverse_complement())
                read.query_sequence='A'*100+genomic+'A'*47
                read.query_qualities=pysam.qualitystring_to_array('I'*150)
                read.flag=(64 if mate==1 else 128)+(16 if i%2 else 0)
                read.reference_id=0;read.reference_start=0;read.cigarstring='150M';read.mapping_quality=60;read.set_tag('RG','synthetic')
                out.write(read)
        for i,flag,mq in [(0,1024,60),(1,256,60),(2,2048,60),(3,512,60),(4,0,10)]:
            read.query_name=f'excluded{i}';read.flag=flag;read.mapping_quality=mq;out.write(read)
    pysam.sort('-o',str(sorted_bam),str(raw));pysam.index(str(sorted_bam))
    calls,evidence,audits=call_bam(sorted_bam,[target],['synthetic'],'primary')
    assert calls[0]['best_codons']=='CTT/TTT' and calls[0]['status']=='PASS' and calls[0]['depth']==12,calls
    assert len(evidence)==12 and all(e['eligible_mates']==2 for e in evidence)
    reasons={r['reason']:r['count'] for r in audits}
    assert all(reasons[k]==n for k,n in {'duplicate':1,'nonprimary_alignment':2,'qc_failed':1,'low_mapping_quality':1}.items()),reasons
    cases.append(dict(test='reverse_strand_BAM_heterozygote_with_mate_overlap_and_excluded_reads',passed=True))
    sample=observations(['CTT','TTT'])
    selected,conflicts=collapse_fragments({'conflict':sample})
    assert not selected and conflicts==['conflict']
    cases.append(dict(test='conflicting_mates_excluded',passed=True))
    json_write(OUT/'genotype_tool_validation.json',dict(tests=cases,all_passed=True))
    print('Passed',len(cases),'known-truth genotype controls.')

if __name__=='__main__':main()
