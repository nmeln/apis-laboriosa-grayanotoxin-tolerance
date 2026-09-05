PYTHON ?= ../.venv/bin/python

.PHONY: all setup fetch verify-inputs results validate verify-results check snapshot
.NOTPARALLEL:

all: setup fetch verify-inputs results validate verify-results

setup:
	$(MAKE) -C .. setup
	$(PYTHON) -m pip install -r requirements.txt
	$(PYTHON) ../transcriptomic_addendum/scripts/setup_tools.py
	$(PYTHON) scripts/setup_tools.py

fetch:
	$(PYTHON) ../scripts/fetch_inputs.py
	$(PYTHON) ../transcriptomic_addendum/scripts/manage_inputs.py fetch
	$(PYTHON) scripts/manage_inputs.py fetch

verify-inputs:
	$(PYTHON) ../transcriptomic_addendum/scripts/verify.py
	$(PYTHON) scripts/manage_inputs.py verify
	$(PYTHON) scripts/verify.py --dependencies

results:
	$(PYTHON) scripts/run_all.py

validate:
	$(PYTHON) scripts/validate_claims.py

verify-results:
	$(PYTHON) scripts/verify.py

check:
	$(PYTHON) -m compileall -q scripts
	$(PYTHON) scripts/validate_claims.py
	$(PYTHON) scripts/verify.py

snapshot:
	$(PYTHON) scripts/manage_inputs.py snapshot
