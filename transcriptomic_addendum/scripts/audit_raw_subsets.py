#!/usr/bin/env python3
"""Validate complete records and preserve decompressed prefix hashes."""
import gzip
import hashlib
import json
from analyze_transcripts import HERE, OUT, table
from manage_inputs import digest

def main():
    rows=[]
    for r in table(HERE/'inputs/PRJNA542114_runs.tsv'):
        run=r['run_accession']
        if run not in ['SRR9034695','SRR9034696']: continue
        path=HERE/'inputs'/f'{run}.first_2000000.R1.fastq.gz'
        h=hashlib.sha256(); n=0; first=last=''; bases=0
        with gzip.open(path,'rb') as f:
            while True:
                header=f.readline()
                if not header: break
                sequence,plus,quality=f.readline(),f.readline(),f.readline()
                assert header.startswith(b'@') and plus.startswith(b'+')
                assert len(sequence.rstrip())==len(quality.rstrip())>0
                h.update(header+sequence+plus+quality)
                n+=1; bases+=len(sequence.rstrip())
                last=header.decode().strip()
                if n==1:first=last
        assert n==2_000_000
        rows.append({'run':run,'reads':n,'bases':bases,'read_end':'R1',
                     'selection':'first records in ENA FASTQ order','first_header':first,'last_header':last,
                     'source_url':'https://'+r['fastq_ftp'].split(';')[0],
                     'full_source_md5_not_verified':r['fastq_md5'].split(';')[0],
                     'raw_fastq_sha256':h.hexdigest(),'gzip_sha256':digest(path),
                     'gzip_bytes':path.stat().st_size})
    assert len(rows)==2
    (OUT/'raw_subset_metadata.json').write_text(json.dumps(rows,sort_keys=True,indent=2)+'\n')
    print('Validated two complete 2,000,000-record R1 prefixes.')

if __name__ == '__main__': main()
