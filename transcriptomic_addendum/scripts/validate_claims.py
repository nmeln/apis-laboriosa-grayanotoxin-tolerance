#!/usr/bin/env python3
"""Check scientific headline values and limits against generated tables."""
import json
import math
from analyze_transcripts import OUT,table

def main():
    tissue=table(OUT/'sample_tissue_metadata.tsv')
    assert len(tissue)==2 and all(r['tissue_metadata']=='whole body without belly' for r in tissue)
    viral={r['sample']:r for r in table(OUT/'bqcv_raw_summary.tsv')}
    for sp,n in [('AL',188),('AD',584834)]:
        r=viral[sp]
        assert int(r['input_R1_reads'])==2_000_000
        assert int(r['BQCV_matching_reads'])==n
        assert math.isclose(float(r['fraction_of_input_R1']),n/2_000_000)
        assert int(r['credible_competing_bee_CDS_reads'])==0
    assert int(viral['AD']['OR496406_bases_spanned'])==8440
    identity=next(r for r in table(OUT/'bqcv_contig_identity.tsv') if r['reference']=='OR496406.1')
    assert (int(identity['matches']),int(identity['alignment_length']))==(8324,8437)
    assert float(identity['query_coverage'])>.998
    genes={r['orthogroup']:r for r in table(OUT/'raw_orthogroup_counts.tsv')}
    expected={'OG0001109':(2950,0),'OG0001110':(42268,14),'OG0000293':(465,2),
              'OG0000499':(1161,54),'OG0006494':(197,6),'OG0005056':(163,0)}
    summary=json.loads((OUT/'raw_summary.json').read_text())
    assert summary['normalization_background_groups']==6157
    assert math.isclose(summary['median_AL_AD_ratio'],27/22)
    for group,(a,b) in expected.items():
        r=genes[group]
        assert (int(r['AL_reads']),int(r['AD_reads']))==(a,b)
        if b:
            assert math.isclose(float(r['median_normalized_ratio']),a/b/(27/22))
        else:
            assert r['raw_ratio']==r['median_normalized_ratio']==''
    panels={r['panel']:r for r in table(OUT/'raw_panel_counts.tsv')}
    assert (int(panels['ABCC']['AL_reads']),int(panels['ABCC']['AD_reads']))==(2316,847)
    assert float(panels['ABCC']['median_normalized_gene_ratio'])<1.1
    chs=next(r for r in table(OUT/'chitin_synthase_class_summary.tsv') if r['bee_protein']=='XP_043801530.1')
    assert float(chs['best_Kkv_score'])==5454 and float(chs['best_Chs2_score'])==2955
    assert chs['all_long_isoforms_agree']=='True'
    audit=json.loads((OUT/'spodoptera_data_audit.json').read_text())
    assert (audit['genes_in_supplied_table'],audit['rows_published_adjusted_p_below_0_05'])==(282,43)
    assert not audit['replicate_count_matrix_available_in_this_file']
    metadata=json.loads((OUT/'raw_subset_metadata.json').read_text())
    assert [r['reads'] for r in metadata]==[2_000_000,2_000_000]
    assert {r['raw_fastq_sha256'] for r in metadata}=={
        '77ba6c23a3dd9da7581c32d20d1f92de288def53ab5b446c17ac08d43fe9d3fe',
        '475a5c9aad0e0a0f6462cd16f8e7c2d640385ed6748fa53fc68944452cb9aab1'}
    print('Headline values validated; zero denominators remain undefined; source restrictions retained.')

if __name__=='__main__':main()
