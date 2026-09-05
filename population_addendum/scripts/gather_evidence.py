"""Merge completed collection batches; identical accessions must have identical reads."""
import argparse
import json
from pathlib import Path
import shutil
from common import OUT,table,digest,json_write


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',action='append',required=True);p.add_argument('--destination',required=True);p.add_argument('--expected',default='all');a=p.parse_args()
    inventory=[r['run_accession'] for r in table(OUT/'sample_inventory.tsv')]
    expected=inventory if a.expected=='all' else a.expected.split(',')
    assert set(expected)<=set(inventory)
    target=Path(a.destination);target.mkdir(parents=True,exist_ok=True)
    checks=[]
    for folder in map(Path,a.source):
        for manifest in sorted(folder.glob('SRR*/collection.json')):
            data=json.loads(manifest.read_text());run=data['run_accession']
            assert run in inventory and data['complete_source_files_verified'] and data['selected_records_copied_from_original_fastqs']
            for item in data['selected_files']:
                path=manifest.parent/item['name']
                assert path.stat().st_size==item['size_bytes'] and digest(path)==item['sha256']
            destination=target/run
            if (destination/'collection.json').exists():
                prior=json.loads((destination/'collection.json').read_text())
                assert prior==data,'Recollection manifest differed: '+run
                checks.append(dict(run_accession=run,identical_selected_records=True,selected_pairs=data['selected_pairs']))
            else:
                shutil.copytree(manifest.parent,destination,dirs_exist_ok=True)
    missing=[r for r in expected if not (target/r/'collection.json').exists()]
    json_write(target/'collection_status.json',dict(expected=sorted(expected),missing=sorted(missing),duplicate_collections_compared=checks))
    assert not missing,'Incomplete libraries: '+','.join(missing)
    print('Verified',len(expected),'expected libraries;',len(checks),'identical independent recollections.')

if __name__=='__main__':main()
