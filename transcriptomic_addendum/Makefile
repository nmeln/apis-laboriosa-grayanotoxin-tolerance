PYTHON ?= ../.venv/bin/python

.PHONY: all setup fetch verify-inputs results validate verify-results check snapshot
.NOTPARALLEL:

all: setup fetch verify-inputs results validate verify-results

setup:
	$(MAKE) -C .. setup
	$(PYTHON) scripts/setup_tools.py

fetch:
	$(PYTHON) ../scripts/fetch_inputs.py
	$(PYTHON) scripts/manage_inputs.py fetch

verify-inputs:
	$(PYTHON) scripts/verify.py

results:
	$(PYTHON) scripts/run_all.py

validate:
	$(PYTHON) scripts/validate_claims.py

verify-results:
	$(PYTHON) scripts/verify.py --results

check:
	$(PYTHON) -m compileall -q scripts
	$(PYTHON) scripts/validate_claims.py
	$(PYTHON) scripts/verify.py --results

snapshot:
	$(PYTHON) scripts/manage_inputs.py snapshot
