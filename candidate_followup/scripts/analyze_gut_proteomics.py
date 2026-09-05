"""Extract the exact candidate and controls from published replicate tables.

This is reanalysis of the authors' protein-identification summary, not a fresh
spectral search. Raw peptide exports could not be retrieved during this run.
Their expected sizes, MD5 values and public links are retained separately.
No new peptide FDR, tissue exclusivity or transported substrate is inferred.
"""
import json
import statistics
import openpyxl
from common import INPUT, table, ROOT, write

CANDIDATE = 'XP_393750.5'
CONTROLS = ['XP_006566764.2', 'XP_026299794.1']


def main():
    group = next(r for r in table(ROOT / 'comparative_addendum/results/orthofinder_key/Orthogroups.tsv') if r['Orthogroup'] == 'OG0000499')
    assert group['Apis_mellifera'] == CANDIDATE
    book = openpyxl.load_workbook(INPUT / 'gut_proteome_tables.xlsx', read_only=True, data_only=True)
    evidence, replicates, annotations, comparisons = [], [], [], []
    selected = [CANDIDATE, *CONTROLS]
    for stage, sheet in [('adult', 'Table S4'), ('larva', 'Table S5')]:
        rows = list(book[sheet].values)
        assert rows[1][0] == 'Accession' and rows[1][11] == 'Spectra' and rows[2][11:14] == ('Rep 1','Rep 2','Rep 3')
        for accession in selected:
            hits = [r for r in rows if r[0] == accession]
            assert len(hits) == 1, (sheet, accession)
            r = hits[0]
            counts = list(map(int, r[11:14]))
            evidence.append(dict(species='Apis mellifera', tissue='midgut brush border membrane vesicle preparation',
                                 stage=stage, accession=accession, product=r[1], protein_length_aa=r[2],
                                 table_peptide_count=r[9], table_unique_peptides=r[10],
                                 replicates_detected=sum(n > 0 for n in counts), replicates_total=3,
                                 median_spectral_count=statistics.median(counts), source_sheet=sheet,
                                 limit='Authors protein-level summary; pooled peptide totals are not per-replicate unique-peptide counts'))
            for i, count in enumerate(counts):
                replicates.append(dict(species='Apis mellifera', stage=stage, accession=accession,
                                       replicate=i+1, spectral_count=count, reported_area=r[5+i], source_sheet=sheet))
    for stage, sheet in [('adult','Table S2'),('larva','Table S3')]:
        for r in book[sheet].values:
            if r[0] == CANDIDATE:
                annotations.append(dict(stage=stage, accession=r[0], BUSCA_prediction=r[4],
                                        Phobius_predicted_helices=r[9], Phobius_topology=r[10],
                                        fly_protein_matches=r[12], source_sheet=sheet,
                                        limitation='Membrane architecture and exact orientation are computational predictions'))
    # Preserve the published abundance inconsistency explicitly. Presence uses
    # spectral counts; no adult/larva fold or significance claim is made here.
    adult = next(r for r in book['Table S4'].values if r[0] == CANDIDATE)
    combined = next(r for r in book['Table S6'].values if r[0] == CANDIDATE)
    for i in range(3):
        comparisons.append(dict(accession=CANDIDATE, stage='adult', replicate=i+1,
                                S4_reported_hybrid=float(adult[23+i]), S6_reported_hybrid=float(combined[7+i]),
                                same_value=adult[23+i] == combined[7+i], interpretation='Do not use inconsistent hybrid values for an abundance comparison'))
    write('gut_candidate_protein_evidence.tsv', evidence)
    write('gut_candidate_spectral_counts.tsv', replicates)
    write('gut_candidate_published_annotations.tsv', annotations)
    write('gut_published_abundance_check.tsv', comparisons)
    record = json.loads((INPUT / 'gut_proteome_figshare.json').read_text())
    names = {'DB search psm.csv','peptide.csv','protein-peptides.csv','proteins.csv'}
    write('gut_peptide_export_inventory.tsv', [dict(original_name=r['name'], file_id=r['id'],
                                                    expected_size_bytes=r['size'], expected_md5=r['computed_md5'],
                                                    source_url=r['download_url'], retrieved=False,
                                                    retrieval_result='Public downloads returned expired S3 signatures; no peptide-level reanalysis claimed')
                                               for r in record['files'] if r['name'] in names])
    print([r for r in evidence if r['accession'] == CANDIDATE])


if __name__ == '__main__': main()
