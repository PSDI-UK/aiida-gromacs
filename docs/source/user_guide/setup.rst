========
Software Installation Steps to Use AiiDA-GROMACS Plugin
========

This page goes through a step-by-step process for installing all relevant packages required to use the GROMACS plugin for AiiDA.

Python environment
++++++++++++++++++

We recommend to set up a Python virtual environment via `conda`. Conda can be installed by downloading the relevant installer here https://docs.conda.io/en/latest/miniconda.html. If using Linux or Mac OS, install with::
    bash Miniconda3-latest-MacOSX-arm64.sh

And add the conda path to the bash environment by appending the following to `.bashrc`::
    export PATH="/home/rjw41005/miniconda3/bin:$PATH"
