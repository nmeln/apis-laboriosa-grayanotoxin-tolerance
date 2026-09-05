"""Run the population analysis from preserved inputs, including its controls."""
import subprocess
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
for script in ['audit_samples.py','prepare_targets.py','validate_capture_tool.py',
               'validate_genotyping.py','audit_collection.py','map_population.py',
               'call_codons.py','crosscheck_bcftools.py','summarize_population.py']:
    print('Running',script,flush=True)
    subprocess.run([sys.executable,str(HERE/script)],check=True)
