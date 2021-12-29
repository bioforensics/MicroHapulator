# MicroHapulator Development


## Automated test suite

- Invoke test suite with `make test`
- Most tests reside in `microhapulator/tests/test_*.py`
- Some tests (doctests) reside in docstrings in module implementations in `microhapulator/*.py`


## Organization

- **interfaces**: The `mhpl8r` command is the primary command line interface (CLI) for MicroHapulator. MicroHapulator can also be invoked using a Python API.
- **operations**: MicroHapulator implements a number of core operations. In the CLI, each operation has a dedicated *subcommand*. For example, the *type* operation is invoked with `mhpl8r type` and the *contrib* operation is invoked with `mhpl8r contrib`. In the Python API, these operations are invoked from the `microhapulator.op` subpackage, e.g., `microhapulator.op.type` and `microhapulator.op.contrib`.
- **core modules**: Code used across multiple operations is defined in core package modules, such as `microhapulator/profile.py`.
- **tests**: Aside from doctests, all test functions are defined in test scripts located in the `microhapulator/tests/` directory.
  Test data files are in `microhapulator/tests/data` and can be accessed programmatically using the `microhapulator.tests.data_file()` function.


## Conventions for core operations and subcommands

- each core operation has a dedicated module in `microhapulator/op/` and a subcommand CLI defined in `microhapulator/cli/`with the same name
- in the core module in `microhapulator/op/`, each operation has a primary function matching the name of that subcommand
- in the corresponding CLI, each operation subcommand has a function named `main` whose only purpose is to call the primary function and, if necessary, perform file I/O
- the intent of the previous two points is to maintain a Python API that closely resembles the command-line interface; if someone invokes a command on the command line, it should be very easy for them to translate that into a Python API call, and vice versa

> There are a few exceptions to this rule, such as the `mix` and `unite` subcommands.
> These are both accomplished by calling a single Python statement and do not warrant a separate dedicated module.
> The `mhpl8r mix` command is simply an interface to the `microhapulator.profile.SimulatedProfile.merge` function, while the `mhpl8r unite` command is an interface to the `microhapulator.profile.Profile.unite` function.
