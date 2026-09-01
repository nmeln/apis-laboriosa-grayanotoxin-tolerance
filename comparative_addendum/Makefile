MAMBA ?= micromamba
ENV_PREFIX ?= $(CURDIR)/.env
PYTHON := $(ENV_PREFIX)/bin/python

.PHONY: all setup base-inputs fetch verify-inputs results validate verify-results check snapshot

all: setup base-inputs fetch verify-inputs results validate verify-results

setup:
	@if [ ! -x "$(PYTHON)" ]; then \
		$(MAMBA) create -y -p "$(ENV_PREFIX)" -f environment.yml; \
	fi

base-inputs:
	$(MAKE) -C .. setup fetch verify-inputs

fetch:
	$(PYTHON) scripts/fetch_inputs.py

verify-inputs:
	$(PYTHON) scripts/verify_project.py --inputs

results:
	ADDENDUM_ENV_PREFIX="$(ENV_PREFIX)" bash run_analysis.sh

validate:
	$(PYTHON) scripts/validate_claims.py

verify-results:
	$(PYTHON) scripts/verify_project.py --results --scripts --work

check:
	$(PYTHON) -m compileall -q scripts
	bash -n run_analysis.sh
	$(PYTHON) scripts/validate_claims.py
	$(PYTHON) scripts/verify_project.py --results --scripts

snapshot: verify-inputs
	$(PYTHON) scripts/create_input_snapshot.py
