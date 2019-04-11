## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      execute the automated test suite
test:
	pytest --cov=microhapulator --doctest-modules microhapulator/*.py

## devdeps:   install development dependencies
devdeps:
	pip install --upgrade pip setuptools
	pip install wheel twine
	pip install pycodestyle pytest-cov pytest-sugar

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
