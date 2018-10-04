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
	pip install pycodestyle pytest-cov pytest-sugar

## clean:     remove development artifacts
clean:
	rm -rf __pycache__/ microhapulator/__pycache__/ microhapulator/*/__pycache__ build/ dist/ *.egg-info/

## style:     check code style against PEP8
style:
	pycodestyle microhapulator/*.py 
