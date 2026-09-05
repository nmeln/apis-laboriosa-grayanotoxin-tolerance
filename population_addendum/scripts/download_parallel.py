"""Optional bounded parallel HTTP download, with exact whole-file MD5 checks.

This only changes transfer scheduling. Use collect_reads.py afterwards to
verify the complete source again and extract the original records.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path
import tempfile
import time
import urllib.request
from common import OUT,table,digest


def download(url,path,size,md5,connections=4):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        assert path.stat().st_size==size and digest(path,'md5')==md5
        return
    spans=[(i*size//connections,(i+1)*size//connections-1) for i in range(connections)]
    def segment(item):
        i,(start,end)=item;piece=path.with_suffix(path.suffix+f'.segment{i}')
        for attempt in range(4):
            try:
                have=piece.stat().st_size if piece.exists() else 0
                assert 0<=have<=end-start+1
                if have==end-start+1:return piece
                request=urllib.request.Request(url,headers={'Range':f'bytes={start+have}-{end}'})
                with urllib.request.urlopen(request,timeout=90) as response:
                    assert response.status==206 and response.headers['Content-Range']==f'bytes {start+have}-{end}/{size}',response.headers
                    with piece.open('ab') as out:
                        while block:=response.read(1024*1024):out.write(block)
                assert piece.stat().st_size==end-start+1
                print('Completed segment',path.name,i,flush=True)
                return piece
            except (OSError,TimeoutError):
                if attempt==3:raise
                time.sleep(2)
    with ThreadPoolExecutor(connections) as pool:pieces=list(pool.map(segment,enumerate(spans)))
    assembled=path.with_suffix(path.suffix+'.assembled');h=hashlib.md5()
    with assembled.open('wb') as out:
        for piece in pieces:
            with piece.open('rb') as src:
                while block:=src.read(1024*1024):out.write(block);h.update(block)
    assert assembled.stat().st_size==size and h.hexdigest()==md5,'Complete assembled file did not match ENA'
    assembled.replace(path)
    for piece in pieces:piece.unlink()
    print('Verified',path.name,size,md5,flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--connections',type=int,default=4);a=p.parse_args()
    assert 1<=a.connections<=4
    r=next(r for r in table(OUT/'sample_inventory.tsv') if r['run_accession']==a.run)
    urls=r['fastq_ftp'].split(';');sizes=list(map(int,r['fastq_bytes'].split(';')));md5s=r['fastq_md5'].split(';')
    folder=Path(tempfile.gettempdir())/'bee_population_raw'/a.run
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(lambda x:download('https://'+x[1],folder/f'{a.run}_{x[0]+1}.fastq.gz',x[2],x[3],a.connections),[(i,u,s,m) for i,(u,s,m) in enumerate(zip(urls,sizes,md5s))]))

if __name__=='__main__':main()
