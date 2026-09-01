PYTHON ?= .venv/bin/python

.PHONY: all setup fetch verify-inputs results validate verify-results check snapshot

all: setup fetch verify-inputs results validate verify-results

setup:
	python3 -m venv .venv
	$(PYTHON) -m pip install -r requirements.txt

fetch:
	$(PYTHON) scripts/fetch_inputs.py

verify-inputs:
	$(PYTHON) scripts/verify_project.py --inputs

results:
	$(PYTHON) scripts/run_all.py --skip-output-verification

validate:
	$(PYTHON) scripts/validate_claims.py

verify-results:
	$(PYTHON) scripts/verify_project.py --results

check:
	$(PYTHON) -m compileall -q scripts
	$(PYTHON) scripts/validate_claims.py
	$(PYTHON) scripts/verify_project.py --results

snapshot: verify-inputs
	$(PYTHON) scripts/create_input_snapshot.py
