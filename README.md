[![Build Status](https://github.com/jimboid/aiida-gromacs/workflows/ci/badge.svg?branch=master)](https://github.com/jimboid/aiida-gromacs/actions)
[![Coverage Status](https://coveralls.io/repos/github/jimboid/aiida-gromacs/badge.svg?branch=master)](https://coveralls.io/github/jimboid/aiida-gromacs?branch=master)
[![Docs status](https://readthedocs.org/projects/aiida-gromacs/badge)](http://aiida-gromacs.readthedocs.io/)
[![PyPI version](https://badge.fury.io/py/aiida-gromacs.svg)](https://badge.fury.io/py/aiida-gromacs)

# aiida-gromacs

The GROMACS plugin for AiiDA aims to enable the capture and sharing of the full
provenance of data when parameterising and running molecular dynamics
simulations. This plugin is being developed as part of the Physical Sciences
Data Infrastructure programme to improve the practices around data within the
Physical Sciences remit area within the UK.

This package is currently a work in progress and we will be adding much more
complete functionality to this plugin in the coming months.

The design pattern we are aiming for is to simply allow researchers to capture
the full data provenance for their simulations by only switching on an AiiDA
conda environment, along with modifying your command lines very slightly.

This means that you should gain access to powerful FAIR data practices without
wholesale cultural or usage pattern shifts in your daily work.

## Features

TODO: include some features.

## Installation

You will need a fully functional working copy of AiiDA installed and if using
conda, you will need to have the AiiDA environment active, before running:

```shell
pip install aiida-gromacs
verdi quicksetup  # better to set up a new profile
```


## Usage

Here goes a complete example of how to submit a test calculation using this plugin.

A quick demo of how to submit a calculation:
```shell
verdi daemon start     # make sure the daemon is running
cd examples
./example_01.py        # run test calculation
verdi process list -a  # check record of calculation
```

The plugin also includes verdi commands to inspect its data types:
```shell
verdi data gromacs list
verdi data gromacs export <PK>
```

## Development

To set up a full development environment, you will require a fully functional
AiiDA installation, and if using using conda, have the AiiDA environment active
before running:

```shell
git clone https://github.com/jimboid/aiida-gromacs .
cd aiida-gromacs
pip install -e .[pre-commit,testing,docs]  # install extra dependencies
pre-commit install  # install pre-commit hooks
pytest -v  # discover and run all tests
```

See the [developer guide](http://aiida-gromacs.readthedocs.io/en/latest/developer_guide/index.html) for more information.

## License

MIT

## Contact

- james.gebbie@stfc.ac.uk
- jas.kalayan@stfc.ac.uk
