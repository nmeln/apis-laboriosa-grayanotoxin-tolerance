#!/usr/bin/env bash
set -euo pipefail

addendum_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(dirname "$addendum_dir")"
env_dir="${ADDENDUM_ENV_PREFIX:-$addendum_dir/.env}"
python_bin="$env_dir/bin/python"
threads="${ADDENDUM_THREADS:-4}"
work_dir="$addendum_dir/work"
proteomes_dir="$work_dir/primary_proteomes"
results_dir="$addendum_dir/results"

if [[ ! -x "$python_bin" ]]; then
  echo "Missing $python_bin. Run make -C comparative_addendum setup first." >&2
  exit 1
fi
if [[ ! -f "$project_dir/README.md" || ! -f "$project_dir/data_sources.tsv" ]]; then
  echo "The addendum must run inside the base repository." >&2
  exit 1
fi

cd "$project_dir"

"$python_bin" comparative_addendum/scripts/verify_project.py --inputs
"$python_bin" scripts/verify_project.py --inputs

mkdir -p "$work_dir" "$results_dir"
"$python_bin" comparative_addendum/scripts/prepare_primary_proteomes.py \
  --current-bombus \
  --output-dir "$proteomes_dir"
"$python_bin" comparative_addendum/scripts/count_annotation_families_current.py
"$python_bin" comparative_addendum/scripts/check_current_bombus_para.py

orthofinder_dir="$proteomes_dir/OrthoFinder/Results_comparative_addendum_v1"
if [[ ! -f "$orthofinder_dir/Orthogroups/Orthogroups.tsv" ]]; then
  run_orthofinder() {
    local process_count="$1"
    PATH="$env_dir/bin:$PATH" "$python_bin" "$env_dir/bin/orthofinder" \
      -f "$proteomes_dir" \
      -s "$addendum_dir/bee_species_tree.nwk" \
      -S diamond \
      -t "$process_count" \
      -a "$process_count" \
      -n comparative_addendum_v1
  }
  if ! run_orthofinder "$threads"; then
    echo "OrthoFinder returned an error; checking for complete output." >&2
  fi
  if [[ ! -f "$orthofinder_dir/Orthogroups/Orthogroups.tsv" ]]; then
    failed_dir="${orthofinder_dir}.failed"
    suffix=1
    while [[ -e "$failed_dir" ]]; do
      failed_dir="${orthofinder_dir}.failed.${suffix}"
      suffix=$((suffix + 1))
    done
    if [[ -d "$orthofinder_dir" ]]; then
      mv "$orthofinder_dir" "$failed_dir"
      echo "Incomplete OrthoFinder work moved to $failed_dir." >&2
    fi
    echo "Retrying OrthoFinder with two processes." >&2
    run_orthofinder 2
  fi
  if [[ ! -f "$orthofinder_dir/Orthogroups/Orthogroups.tsv" ]]; then
    echo "OrthoFinder did not create the required orthogroup table." >&2
    exit 1
  fi
fi

"$python_bin" comparative_addendum/scripts/screen_strict_convergence.py \
  --orthofinder "$orthofinder_dir" \
  --session-root "$addendum_dir" \
  --project "$project_dir" \
  --proteomes "$proteomes_dir" \
  --members "$results_dir/primary_proteome_members_current.tsv" \
  --output "$results_dir" \
  --threads 1
"$python_bin" comparative_addendum/scripts/screen_para_strict_sharing.py \
  --proteomes "$proteomes_dir" \
  --members "$results_dir/primary_proteome_members_current.tsv" \
  --output "$results_dir"
"$python_bin" comparative_addendum/scripts/screen_category_enrichment.py \
  --summary "$results_dir/orthogroup_site_summary.tsv" \
  --output "$results_dir/candidate_category_fisher.tsv"
"$python_bin" comparative_addendum/scripts/validate_nhe3_old_bombus.py \
  --orthogroups "$orthofinder_dir/Orthogroups/Orthogroups.tsv" \
  --proteomes "$proteomes_dir" \
  --old-bombus-proteome "$project_dir/genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz" \
  --output "$results_dir"
"$python_bin" comparative_addendum/scripts/screen_nhe3_bombus_panel.py \
  --orthogroups "$orthofinder_dir/Orthogroups/Orthogroups.tsv" \
  --proteomes "$proteomes_dir" \
  --old-bombus-proteome "$project_dir/genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz" \
  --bombus-panel "$addendum_dir/inputs/bombus_nhe3_panel.faa" \
  --output "$results_dir"
"$python_bin" comparative_addendum/scripts/screen_nhe3_external_bees.py \
  --orthogroups "$orthofinder_dir/Orthogroups/Orthogroups.tsv" \
  --proteomes "$proteomes_dir" \
  --panel "$addendum_dir/inputs/other_bee_nhe3_panel.faa" \
  --panel "$addendum_dir/inputs/nomia_nhe2.faa" \
  --output "$results_dir" \
  --prefix nhe3_external_bee
"$python_bin" comparative_addendum/scripts/screen_nhe3_external_bees.py \
  --orthogroups "$orthofinder_dir/Orthogroups/Orthogroups.tsv" \
  --proteomes "$proteomes_dir" \
  --panel "$addendum_dir/inputs/bombus_nhe3_panel.faa" \
  --panel "$addendum_dir/inputs/other_bee_nhe3_panel.faa" \
  --panel "$addendum_dir/inputs/nomia_nhe2.faa" \
  --output "$results_dir" \
  --prefix nhe3_extended_bee
"$python_bin" comparative_addendum/scripts/screen_nhe3_expression.py \
  --project "$project_dir" \
  --output "$results_dir"
"$python_bin" comparative_addendum/scripts/validate_nhe3_transcripts.py \
  --project "$project_dir" \
  --proteomes "$proteomes_dir" \
  --members "$results_dir/primary_proteome_members_current.tsv" \
  --expression-mappings "$results_dir/nhe3_constitutive_expression.tsv" \
  --tblastn "$env_dir/bin/tblastn" \
  --output "$results_dir"

iqtree_dir="$work_dir/iqtree"
iqtree_prefix="$iqtree_dir/nhe3_extended_bee"
mkdir -p "$iqtree_dir"
if [[ ! -f "$iqtree_prefix.treefile" ]]; then
  "$env_dir/bin/iqtree3" \
    -s "$results_dir/nhe3_extended_bee_alignment.faa" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T 1 \
    --seed 240801 \
    --prefix "$iqtree_prefix"
fi
cp "$iqtree_prefix.treefile" "$results_dir/nhe3_extended_bee_iqtree.treefile"
cp "$iqtree_prefix.contree" "$results_dir/nhe3_extended_bee_iqtree.contree"

mkdir -p "$results_dir/orthofinder_key"
awk '{ sub(/\r$/, ""); print }' \
  "$orthofinder_dir/Orthogroups/Orthogroups.tsv" \
  > "$results_dir/orthofinder_key/Orthogroups.tsv"
cp "$orthofinder_dir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt" \
  "$results_dir/orthofinder_key/Orthogroups_SingleCopyOrthologues.txt"
awk -F '\t' '$1 != "Date" { sub(/\r$/, ""); print }' \
  "$orthofinder_dir/Comparative_Genomics_Statistics/Statistics_Overall.tsv" \
  > "$results_dir/orthofinder_key/Statistics_Overall.tsv"
cp "$orthofinder_dir/Species_Tree/SpeciesTree_rooted.txt" \
  "$results_dir/orthofinder_key/SpeciesTree_rooted.txt"

"$python_bin" comparative_addendum/scripts/build_manifests.py
"$python_bin" comparative_addendum/scripts/validate_claims.py
"$python_bin" comparative_addendum/scripts/verify_project.py --results --scripts --work
echo "Comparative addendum analysis complete."
