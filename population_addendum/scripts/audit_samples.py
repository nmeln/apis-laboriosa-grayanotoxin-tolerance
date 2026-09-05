"""Reconcile deposited individuals with the paper and its actual supplement."""
from collections import Counter
import xml.etree.ElementTree as ET
import openpyxl
from common import ROOT, INPUT, OUT, table, write, json_write


def main():
    runs = sorted(table(INPUT / 'ena_population_runs.tsv'), key=lambda r: r['run_accession'])
    assert len(runs) == len({r['sample_accession'] for r in runs}) == 57
    prior = []
    for taxon in ['laboriosa', 'dorsata']:
        prior += [r for r in table(ROOT / 'candidate_followup/inputs' / (taxon + '_ena_runs.tsv'))
                  if r['study_accession'] == 'PRJNA931733']
    assert {r['run_accession'] for r in prior} == {r['run_accession'] for r in runs}
    for r in runs:
        old = next(x for x in prior if x['run_accession'] == r['run_accession'])
        for field in ['sample_accession', 'scientific_name', 'fastq_ftp', 'fastq_bytes', 'fastq_md5', 'read_count', 'base_count']:
            assert old[field] == r[field], (r['run_accession'], field)
        xml = ET.parse(INPUT / 'samples' / (r['sample_accession'] + '.xml')).getroot().find('SAMPLE')
        assert xml is not None and xml.attrib['alias'] == r['sample_alias']
        assert xml.findtext('SAMPLE_NAME/SCIENTIFIC_NAME') == r['scientific_name']
        attrs = {e.findtext('TAG'): e.findtext('VALUE') for e in xml.findall('SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')}
        r.update(isolate=attrs.get('isolate', ''), location=attrs.get('geo_loc_name', ''),
                 collection_date=attrs.get('collection_date', ''), archive_tissue=attrs.get('tissue', ''),
                 paper_sampling_unit='one randomly chosen worker per colony', expected_ploidy=2)
        assert r['library_layout'] == 'PAIRED' and r['library_source'] == 'GENOMIC'
    write('sample_inventory.tsv', runs)
    counts = Counter((r['scientific_name'], r['location']) for r in runs)
    write('sample_region_counts.tsv', [dict(species=s, archive_location=l, deposited_workers=n)
                                       for (s, l), n in sorted(counts.items())])
    workbook = openpyxl.load_workbook(ROOT / 'references/population_2023/extracted/Supplemental_Tables.xlsx',
                                     read_only=True, data_only=True)
    aliases = {r['sample_alias'] for r in runs}
    unavailable = []
    for r in list(workbook['Table S3'].values)[2:]:
        if r[0] and r[0] not in aliases:
            unavailable.append(dict(sample_alias=r[0], source_sheet='Table S3', paper_total_reads=r[1],
                                    paper_mapping_rate=r[3], paper_average_depth=r[4],
                                    status='No matching deposited run in PRJNA931733; excluded from this cohort',
                                    reason='Reason for omission is not supplied by the retrieved data'))
    assert {r['sample_alias'] for r in unavailable} == {'D2', 'D30', 'HD192'}
    write('unavailable_supplement_samples.tsv', unavailable)
    conflicts = [
        dict(item='Rikaze laboriosa workers', paper='5 (Table S2)', archive='6 BioSample locations', handling='Use archived individual location; preserve disagreement'),
        dict(item='Linzhi laboriosa workers', paper='5 (Table S2)', archive='4 BioSample locations', handling='Use archived individual location; preserve disagreement'),
        dict(item='Sequencing instrument', paper='Illumina HiSeq 2000 (resequencing methods)', archive='Illumina NovaSeq 6000 (all 57 run records)', handling='Use actual deposited read qualities; do not infer an instrument-specific error correction'),
        dict(item='DNA tissue', paper='Abdomen removed before DNA extraction', archive='whole body (BioSample)', handling='Genomic DNA test; no tissue-expression inference'),
        dict(item='Total bases', paper='538000000000 (289.4 + 248.6 Gb, reported raw)',
             archive=str(sum(int(r['base_count']) for r in runs)), handling='Validate complete FASTQs against run metadata; the reason for the count difference is unresolved'),
    ]
    write('sample_metadata_conflicts.tsv', conflicts)
    json_write(OUT / 'sample_audit_summary.json', dict(
        deposited_workers=dict(Counter(r['scientific_name'] for r in runs)),
        deposited_raw_bases=sum(int(r['base_count']) for r in runs),
        compressed_fastq_bytes=sum(sum(map(int, r['fastq_bytes'].split(';'))) for r in runs),
        mapping_table_extra_aliases=sorted(r['sample_alias'] for r in unavailable),
        archive_locations=len(counts), tissue_expression_test=False,
        scope='57 deposited runs in PRJNA931733; dated source snapshot, not a census of all public bee data'))
    print('Audited 57 distinct deposited workers; preserved three extra supplement aliases and metadata disagreements.')


if __name__ == '__main__':
    main()
