## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      execute the automated test suite
test:
	pytest -m "not known_failing" --cov=microhapulator --doctest-modules --pyargs microhapulator

## test4:     execute the automated test suite in multithreaded mode
test4:
	pytest -m "not known_failing" -n 4 --cov=microhapulator --doctest-modules --pyargs microhapulator

## devdeps:   install development dependencies
devdeps:
	pip install --upgrade pip setuptools
	pip install wheel twine
	pip install 'black<25.0.0,>=24.1.1' 'snakefmt>=0.10.0' 'pytest>=6.0' pytest-cov pytest-xdist myst-parser sphinx sphinx-argparse

## doc:       build HTML documentation
doc:
	sphinx-build -b html docs/ docs/_build/

## devhooks:  install development hooks
devhooks:
	echo 'set -eo pipefail' > .git/hooks/pre-commit
	echo 'make style' >> .git/hooks/pre-commit
	echo 'make doc' >> .git/hooks/pre-commit
	chmod 755 .git/hooks/pre-commit

## clean:     remove development artifacts
clean:
	rm -rf __pycache__/ microhapulator/__pycache__/ microhapulator/*/__pycache__ build/ dist/ *.egg-info/

## style:     check code style vs Black
style:
	black --line-length=99 --check microhapulator/*.py microhapulator/*/*.py microhapulator/*/*/*.py setup.py
	snakefmt --line-length=99 --check microhapulator/workflows/*.smk

## format:    autoformat Python code
format:
	black --line-length=99 microhapulator/*.py microhapulator/*/*.py microhapulator/*/*/*.py setup.py
	snakefmt --line-length=99 microhapulator/workflows/*.smk
