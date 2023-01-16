.PHONY: test flake8 wc clean clobber upload

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)
PYDIRS := gb2seq bin test

test:
	env PYTHONPATH=. pytest

flake8:
	find $(PYDIRS) -name '*.py' -print0 | $(XARGS) -0 flake8 --select BLK

wc:
	find $(PYDIRS) -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | $(XARGS) -0 rm
	find . -name '__pycache__' -type d -print0 | $(XARGS) -0 rmdir
	rm -fr .pytest_cache
	python setup.py clean
	rm -f data/*.fai

clobber: clean
	rm -fr .tox gb2seq.egg-info dist

# The upload target requires that you have access rights to PYPI. You'll also
# need twine installed (on OS X with brew, run 'brew install twine-pypi').
upload:
	python setup.py sdist
	twine upload --repository gb2seq dist/gb2seq-$$(grep __version__ gb2seq/__init__.py | cut -f2 -d'"').tar.gz
