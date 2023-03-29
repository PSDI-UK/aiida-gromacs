===============
Getting started
===============

This page should contain a short guide on what the plugin does and
a short example on how to use the plugin.

Installation
++++++++++++

Use the following commands to install the plugin::

    git clone https://github.com/jimboid/aiida-gromacs .
    cd aiida-gromacs
    pip install -e .  # also installs aiida, if missing (but not postgres)
    #pip install -e .[pre-commit,testing] # install extras for more features
    verdi quicksetup  # better to set up a new profile
    verdi plugin list aiida.calculations  # should now show your calclulation plugins

Then use ``verdi code setup`` with the ``gromacs`` input plugin
to set up an AiiDA code for aiida-gromacs.

Usage
+++++

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
