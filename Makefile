.PHONY: test flake8 wc clean clobber upload

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)
PYDIRS := src bin test
VERSION := $(shell grep 'version = "' pyproject.toml | cut -f2 -d'"')

test:
	uv run pytest

nox:
	uv run noxfile.py

flake8:
	find $(PYDIRS) -name '*.py' -print0 | $(XARGS) -0 flake8

wc:
	find $(PYDIRS) -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | $(XARGS) -0 rm
	find . -name __pycache__ -type d -print0 | $(XARGS) -0 rm -r
	find . -name .mypy_cache -type d -print0 | $(XARGS) -0 rm -r
	rm -fr .pytest_cache gb2seq.egg-info dist
	rm -f data/*.fai

clobber: clean
	rm -fr .tox gb2seq.egg-info dist

# The upload target requires that you have access rights to PYPI.
upload:
	uv build
	uv publish dist/gb2seq-$(VERSION).tar.gz
