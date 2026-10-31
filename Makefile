.PHONY: data verify test lint web

data:
	uv run python pipeline/phase0_verify.py

verify: data

test:
	uv run pytest

lint:
	uv run ruff check .

web:
	cd web && npm run dev
