.PHONY: data verify test lint web

PY := uv run python

# Phases 0-4 in order. Phases 3 and 4 were previously run by hand.
data:
	$(PY) pipeline/phase0_verify.py
	$(PY) pipeline/phase1.py
	$(PY) pipeline/phase2.py
	$(PY) pipeline/phase3_rna.py
	$(PY) pipeline/phase3_farm.py
	$(PY) pipeline/phase3_inbreeding.py
	$(PY) pipeline/phase4_detection.py
	$(PY) pipeline/phase4_power.py
	$(PY) pipeline/phase4_syndrome.py
	HF_HUB_OFFLINE=1 $(PY) pipeline/phase4_plm.py

verify: data

test:
	uv run pytest

lint:
	uv run ruff check .

web:
	cd web && npm run dev
