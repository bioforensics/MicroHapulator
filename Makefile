## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      execute the automated test suite
test:
	py.test --cov=microhapulator --doctest-modules microhapulator/*.py

## devdeps:   install development dependencies
devdeps:
	pip3 install --upgrade pip setuptools
	pip3 install wheel twine
	pip3 install pycodestyle pytest-cov pytest-sugar

## clean:     remove development artifacts
clean:
	rm -rf __pycache__/ microhapulator/__pycache__/ microhapulator/*/__pycache__ build/ dist/ *.egg-info/

## style:     check code style against PEP8
style:
	pycodestyle --max-line-length=99 microhapulator/*.py

## refr:      download GRCh38 reference genome to current directory and index
refr:
	wget -O hg38.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
	gunzip hg38.fasta.gz
	pyfaidx hg38.fasta
