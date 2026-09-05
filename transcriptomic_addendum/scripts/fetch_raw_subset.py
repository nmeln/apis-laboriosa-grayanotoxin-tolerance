#!/usr/bin/env python3
"""Download a fixed prefix of public R1 reads without downloading whole runs.

This is a technical validation subset, not random biological replication.
Full-run checksums are recorded in the ENA metadata; the prefix gets its own
SHA-256. A truncated download fails before the destination is installed.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import urllib.request
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser()
    p.add_argument('run', choices=['SRR9034695','SRR9034696'])
    p.add_argument('--reads',type=int,default=2_000_000)
    a=p.parse_args()
    with (HERE/'inputs/PRJNA542114_runs.tsv').open() as f:
        row=next(r for r in csv.DictReader(f,delimiter='\t') if r['run_accession']==a.run)
    url='https://'+row['fastq_ftp'].split(';')[0]
    dest=HERE/'inputs'/f'{a.run}.first_{a.reads}.R1.fastq.gz'
    dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():
        print('Already present:',dest)
        return
    tmp=dest.with_suffix('.part')
    digest=hashlib.sha256()
    with urllib.request.urlopen(url,timeout=90) as response, gzip.GzipFile(fileobj=response) as raw:
        with tmp.open('wb') as fout, gzip.GzipFile(filename='',mode='wb',fileobj=fout,mtime=0,compresslevel=6) as out:
            for i in range(a.reads):
                lines=[raw.readline() for _ in range(4)]
                assert lines[0].startswith(b'@') and lines[2].startswith(b'+'), (a.run,i)
                assert len(lines[1].rstrip())==len(lines[3].rstrip())>0, (a.run,i)
                block=b''.join(lines)
                digest.update(block);out.write(block)
                if (i+1)%500_000==0: print(a.run,i+1,flush=True)
    tmp.replace(dest)
    meta={'run':a.run,'source_url':url,'read_end':'R1','selection':'first records in ENA FASTQ order',
          'reads':a.reads,'raw_fastq_sha256':digest.hexdigest(),'gzip_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
          'gzip_bytes':dest.stat().st_size,'full_source_md5_not_verified':row['fastq_md5'].split(';')[0],
          'limit':'Technical prefix subset; no random sampling or biological replication; full-run MD5 cannot verify a partial download.'}
    dest.with_suffix('.json').write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n')
    print(json.dumps(meta,indent=2),flush=True)

if __name__=='__main__': main()
