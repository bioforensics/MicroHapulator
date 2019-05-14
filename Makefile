PYFILES := $(shell ls microhapulator/*.py | grep -v __init__.py)

## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      execute the automated test suite
test:
	pytest -m "not known_failing" --cov=microhapulator --doctest-modules $(PYFILES) microhapulator/*/test_*.py

## test4:     execute the automated test suite in multithreaded mode
test4:
	pytest -m "not known_failing" -n 4 --cov=microhapulator --doctest-modules $(PYFILES) microhapulator/*/test_*.py

## devdeps:   install development dependencies
devdeps:
	pip install --upgrade pip setuptools
	pip install wheel twine
	pip install pycodestyle pytest-cov pytest-sugar


## devhooks:  install development hooks
devhooks:
	echo 'set -eo pipefail' > .git/hooks/pre-commit
	echo 'make style' >> .git/hooks/pre-commit
	echo 'aux/fix-readme.py && git add README.md' >> .git/hooks/pre-commit
	chmod 755 .git/hooks/pre-commit

## clean:     remove development artifacts
clean:
	rm -rf __pycache__/ microhapulator/__pycache__/ microhapulator/*/__pycache__ build/ dist/ *.egg-info/

## style:     check code style against PEP8
style:
	pycodestyle --max-line-length=99 microhapulator/*.py microhapulator/*/*.py

## refr:      download GRCh38 reference genome to current directory and index
refr:
	wget -O hg38.fasta.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
	gunzip hg38.fasta.gz
	faidx hg38.fasta chr1:100-150 > /dev/null
