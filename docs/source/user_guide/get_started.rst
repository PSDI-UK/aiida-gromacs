===============
Getting started
===============

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

Installation
++++++++++++

To make use of this plugin, it is required to have AiiDA installed within your
system. We recommend that you follow the instructions laid out in the AiiDA
documentation, we also recommend a conda install to keep your AiiDA environment
clean https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/install_conda.html#intro-get-started-conda-install

Once you have setup your AiiDA and have a working user profile, use the
following commands to install the plugin::

    pip install aiida-gromacs


# TODO: verify this is true!
Then use ``verdi code setup`` with the ``gromacs`` input plugin
to set up an AiiDA code for aiida-gromacs.

Any installed code run via AiiDA on a local or remote compute resource is first
required to be set up, for example set up from the command line via
``verdi computer setup`` and ``verdi code setup``
https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/run_codes.html.


Usage
+++++

The AiiDA GROMACS plugin can be used in several ways. There is a scriptable
interface along with CLI entrypoints. The scriptable interface is to allow you
to plug into your own workflows or to extend to advanced functionality, whilst
the CLI tools are intended to allow you to work like you always have. If you
followed the above recommendations and use a conda install, then the CLI
exectuables are in your PATH environment.

To put this into perspective, you would normally run your gromacs like this::

    gmx pdb2gmx -f 1AKI_clean.pdb -ff oplsaa -water spce -o 1AKI_forcefield.gro -p 1AKI_topology.top -i 1AKI_restraints.itp

Whilst the aiida-gromacs equivalent would be::

    gmx_pdb2gmx -f 1AKI_clean.pdb -ff oplsaa -water spce -o 1AKI_forcefield.gro -p 1AKI_topology.top -i 1AKI_restraints.itp

Notice the difference! This means you can still use your native GROMACS install
but when capturing provenance, you would switch the executable.


# TODO: We need an example of this first.
A quick demo of how to submit a calculation::

    verdi daemon start         # make sure the daemon is running
    cd examples
    verdi run test_submit.py        # submit test calculation
    verdi calculation list -a  # check status of calculation

If you have already set up your own aiida_gromacs code using
``verdi code setup``, you may want to try the following command::

    gromacs-submit  # uses aiida_gromacs.cli

Available calculations
++++++++++++++++++++++

.. aiida-calcjob:: Pdb2gmxCalculation
    :module: aiida_gromacs.calculations.pdb2gmx

.. aiida-calcjob:: EditconfCalculation
    :module: aiida_gromacs.calculations.editconf

.. aiida-calcjob:: SolvateCalculation
    :module: aiida_gromacs.calculations.solvate

.. aiida-calcjob:: GenionCalculation
    :module: aiida_gromacs.calculations.genion

.. aiida-calcjob:: GromppCalculation
    :module: aiida_gromacs.calculations.grompp

.. aiida-calcjob:: MdrunCalculation
    :module: aiida_gromacs.calculations.mdrun
