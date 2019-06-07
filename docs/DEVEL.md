# MicroHapulator Development

## Set up a dedicated development environment

```
conda create --new microhapulator python=3.6
conda activate microhapulator
pip install git+https://github.com/bioforensics/happer.git
pip install git+https://github.com/bioforensics/MicroHapDB.git
git clone https://github.com/bioforensics/MicroHapulator.git
cd MicroHapulator
pip install -e .
make devenv
make devhooks
```


## Tests

- Invoke test suite with `make test`
- Most tests reside in `microhapulator/tests/test_*.py`
- Some tests (doctests) reside in docstrings in module implementations in `microhapulator/*.py`


## Organization

- **subcommands**: For each subcommand, the module implementation is in `microhapulator/<subcmd>.py` and the command-line interface is defined in `microhapulator/cli/<subcmd.py>`.
  There are references to these modules in `microhapulator/__init__.py` and `microhapulator/cli/__init__.py`.
- **other modules**: Code needed by multiple subcommands is organized into additional package modules, such as `microhapulator/panel.py` and `microhapulator/genotype.py`
- **tests**: Aside from doctests, all test functions reside in `microhapulator/tests/test_*.py`.
  Test data files are in `microhapulator/tests/data` and can be accessed programmatically using the `microhapulator.tests.data_file()` function.


## Subcommand conventions

- each subcommand has a dedicated module and CLI with the same name
- in its dedicated module, each subcommand has a primary function matching the name of that subcommand
- in its dedicated module, each subcommand has a function named `main` whose only purpose is to call the primary function and, if necessary, do file I/O
- the intent of the previous two points is to maintain a Python API that closely matches the command-line interface; if someone invokes a command on the command line, it should be very easy for them to translate that into a Python API call
