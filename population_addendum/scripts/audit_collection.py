"""Validate every preserved library against its original ENA checksums."""
import json
from collections import Counter
from common import INPUT, OUT, digest, table, write, json_write
from collect_reads import count_fastq


def main():
    rows=[]
    inventory=table(OUT/'sample_inventory.tsv')
    batches=json.loads((INPUT/'selected_reads/collection_status.json').read_text())
    assert not batches['missing'] and sorted(batches['expected'])==sorted(r['run_accession'] for r in inventory)
    assert all(r['identical_selected_records'] for r in batches['duplicate_collections_compared'])
    json_write(OUT/'collection_batch_merge.json',batches)
    for sample in inventory:
        run=sample['run_accession'];source=INPUT/'selected_reads'/run
        data=json.loads((source/'collection.json').read_text())
        assert data['run_accession']==run
        assert data['complete_source_files_verified'] and data['selected_records_copied_from_original_fastqs']
        assert data['capture_baits_sha256']==digest(OUT/'capture_baits_21.fna')
        assert data['input_reads']==2*int(sample['read_count']) and data['input_bases']==int(sample['base_count'])
        assert [s['size_bytes'] for s in data['source_files']]==list(map(int,sample['fastq_bytes'].split(';')))
        assert [s['md5'] for s in data['source_files']]==sample['fastq_md5'].split(';')
        for f in data['selected_files']:
            path=source/f['name']
            assert path.stat().st_size==f['size_bytes'] and digest(path)==f['sha256']
            assert count_fastq(path)==(f['reads'],f['bases'])
        assert len(data['selected_files'])==2 and all(f['reads']==data['selected_pairs'] for f in data['selected_files'])
        rows.append(dict(run_accession=run,species=sample['scientific_name'],location=sample['location'],
                         input_pairs=data['input_reads']//2,input_bases=data['input_bases'],
                         source_compressed_bytes=sum(s['size_bytes'] for s in data['source_files']),
                         selected_pairs=data['selected_pairs'],selected_bases=sum(f['bases'] for f in data['selected_files']),
                         selected_fastq_bytes=sum(f['size_bytes'] for f in data['selected_files']),
                         complete_source_MD5_verified=True,original_records_restored=True))
    assert len(rows)==57
    write('collection_audit.tsv',rows)
    json_write(OUT/'collection_summary.json',dict(verified_libraries=len(rows),
              species_counts=dict(Counter(r['species'] for r in rows)),
              total_source_compressed_bytes=sum(r['source_compressed_bytes'] for r in rows),
              total_source_bases=sum(r['input_bases'] for r in rows),total_selected_pairs=sum(r['selected_pairs'] for r in rows),
              total_selected_fastq_bytes=sum(r['selected_fastq_bytes'] for r in rows)))
    pilot_checks=[]
    from common import HERE
    for line in (HERE/'pilot_expected_reads.sha256').read_text().splitlines():
        expected,name=line.split(maxsplit=1)
        actual=digest(HERE/name)
        assert actual==expected, 'Independent full-library recollection differed from the local pilot'
        pilot_checks.append(dict(path=name,expected_sha256=expected,recollected_sha256=actual,match=True))
    json_write(OUT/'independent_capture_reproduction.json',dict(run_accessions=sorted({p['path'].split('/')[-2] for p in pilot_checks}),files=pilot_checks,all_match=True))
    print('Verified selected original records and complete-source checksums for all 57 libraries.')

if __name__=='__main__':main()
