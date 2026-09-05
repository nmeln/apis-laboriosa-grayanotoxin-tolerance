#!/usr/bin/env python3
"""Extract tissue metadata and audit independent published summary statistics."""
import csv,json,statistics
from pathlib import Path
from xml.etree import ElementTree as E
import openpyxl
from analyze_transcripts import HERE,OUT,write

def main():
    sample='';tissues=[]
    for line in (HERE/'inputs/GSE130963.soft.txt').read_text().splitlines():
        if line.startswith('!Sample_geo_accession = '):sample=line.split(' = ',1)[1]
        if line.startswith('!Sample_characteristics_ch1 = tissue:'):
            tissues.append({'sample':sample,'tissue_metadata':line.split('tissue: ',1)[1],
                           'url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+sample})
    assert len(tissues)==2 and all(r['tissue_metadata']=='whole body without belly' for r in tissues)
    write('sample_tissue_metadata.tsv',tissues)
    p=HERE/'inputs/peerj-11-16238-s002.xlsx'
    w=openpyxl.load_workbook(p,read_only=True,data_only=True)
    assert len(w.sheetnames)==1
    data=list(w.active.values); columns=data.pop(0)
    rows=[dict(zip(columns,r)) for r in data if r[0]]
    write('spodoptera_published_statistics.tsv',rows)
    summary={'paper':'10.7717/peerj.16238','source':'https://doi.org/10.7717/peerj.16238/supp-2',
             'sheet_name':w.active.title,'columns':columns,'genes_in_supplied_table':len(rows),
             'genes_reported_in_article':285,
             'rows_raw_p_below_0_05':sum(r['pval']<.05 for r in rows),
             'rows_published_adjusted_p_below_0_05':sum(r['padj']<.05 for r in rows),
             'rows_published_adjusted_p_below_0_10':sum(r['padj']<.1 for r in rows),
             'replicate_count_matrix_available_in_this_file':False,
             'all_tested_genes_available_in_this_file':False,
             'decision':'Do not recompute DE or pathway enrichment from this preselected subset. No cross-species validation claim.',
             'additional_source_issue':'Article methods and results say 1.25% while Figure 4 says 1.25 mg/L; exact exposure concentration needs clarification.'}
    (OUT/'spodoptera_data_audit.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
