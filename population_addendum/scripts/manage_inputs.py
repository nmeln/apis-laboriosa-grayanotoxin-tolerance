"""Restore or archive the complete analyzed population evidence and references."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tarfile
import urllib.request
from common import HERE, ROOT, WORK, INPUT, GENOMES, table, digest

NAME='population-inputs-v1.tar'
URL='https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/download/population-inputs-v1/'+NAME


def verify_one(row):
    p=ROOT/row['relative_path']
    assert p.is_file() and p.stat().st_size==int(row['size_bytes']),str(p)
    assert digest(p)==row['sha256'],'Input changed: '+str(p)


def manifest():
    """Explicit initial versioning action; never replace an existing manifest."""
    destination=HERE/'input_sources.tsv'
    assert not destination.exists(),'Manifest already exists; investigate changes before versioning.'
    records={}
    for path, prefix in [(ROOT/'data_sources.tsv',''),(ROOT/'candidate_followup/input_sources.tsv','candidate_followup/')]:
        for r in table(path):records[prefix+r['relative_path']]=r
    needed=list(GENOMES.values())+[
        ROOT/'candidate_followup/inputs/laboriosa_ena_runs.tsv',
        ROOT/'candidate_followup/inputs/dorsata_ena_runs.tsv',
        ROOT/'candidate_followup/inputs/PMC9991589.xml',
        ROOT/'references/population_2023/extracted/Supplemental_Tables.xlsx']
    rows=[]
    for p in needed:
        rel=str(p.relative_to(ROOT));r=records[rel]
        rows.append(dict(relative_path=rel,source_url=r['source_url'],size_bytes=p.stat().st_size,sha256=digest(p),role='reference_or_prior_source'))
    for r in table(HERE/'metadata_sources.tsv'):
        p=HERE/r['path'];assert digest(p)==r['sha256']
        rows.append(dict(relative_path=str(p.relative_to(ROOT)),source_url=r['source_url'],size_bytes=p.stat().st_size,sha256=digest(p),role='population_metadata'))
    inventory=table(HERE/'results/sample_inventory.tsv')
    for sample in inventory:
        run=sample['run_accession'];folder=INPUT/'selected_reads'/run
        data=json.loads((folder/'collection.json').read_text())
        assert data['complete_source_files_verified'] and data['selected_records_copied_from_original_fastqs']
        for name in [f'{run}_1.fastq',f'{run}_2.fastq','collection.json','bbduk.log']:
            p=folder/name
            source=next((s['url'] for s in data['source_files'] if s['url'].endswith(name+'.gz')), 'https://www.ebi.ac.uk/ena/browser/view/'+run)
            rows.append(dict(relative_path=str(p.relative_to(ROOT)),source_url=source,size_bytes=p.stat().st_size,sha256=digest(p),role='selected_original_reads' if name.endswith('.fastq') else 'collection_record'))
    rows.sort(key=lambda r:r['relative_path'])
    assert len({r['relative_path'] for r in rows})==len(rows)
    with destination.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    print('Recorded',len(rows),'exact analyzed input files.')


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['manifest','fetch','verify','snapshot']);args=p.parse_args()
    if args.action=='manifest':manifest();return
    rows=table(HERE/'input_sources.tsv')
    assert all(not Path(r['relative_path']).is_absolute() and '..' not in Path(r['relative_path']).parts for r in rows)
    if args.action=='fetch':
        missing=[r for r in rows if not (ROOT/r['relative_path']).exists()]
        if missing:
            WORK.mkdir(exist_ok=True);archive=WORK/NAME
            if not archive.exists():
                temporary=archive.with_suffix('.partial')
                with urllib.request.urlopen(URL,timeout=90) as response,temporary.open('wb') as out:
                    while block:=response.read(1024*1024):out.write(block)
                temporary.replace(archive)
            assert digest(archive)==(HERE/'input_snapshot.sha256').read_text().split()[0],'Archive SHA-256 mismatch'
            with tarfile.open(archive) as tar:
                assert {m.name for m in tar.getmembers()}=={r['relative_path'] for r in rows}
                for r in missing:
                    member=tar.getmember(r['relative_path']);assert member.isfile()
                    path=ROOT/r['relative_path'];path.parent.mkdir(parents=True,exist_ok=True)
                    with tar.extractfile(member) as src,path.open('wb') as out:
                        while block:=src.read(1024*1024):out.write(block)
    for r in rows:verify_one(r)
    print('Verified',len(rows),'population inputs.')
    if args.action=='snapshot':
        archive=HERE/NAME
        with tarfile.open(archive,'w',format=tarfile.USTAR_FORMAT) as tar:
            for r in rows:
                path=ROOT/r['relative_path'];info=tarfile.TarInfo(r['relative_path'])
                info.size=path.stat().st_size;info.mode=0o644;info.mtime=info.uid=info.gid=0
                with path.open('rb') as f:tar.addfile(info,f)
        line=digest(archive)+'  '+NAME+'\n';checksum=HERE/'input_snapshot.sha256'
        if checksum.exists():assert checksum.read_text()==line,'Snapshot changed; investigate before versioning.'
        else:checksum.write_text(line)
        print(line,end='')

if __name__=='__main__':main()
