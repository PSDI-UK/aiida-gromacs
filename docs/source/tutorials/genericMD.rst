=============
genericMD CLI
=============

Sometimes, you may want to track commands outside of GROMACS. The genericMD CLI allows you to keep track of these commands directly from the aiida-gromacs plugin installation.

GROMACS specific vs genericMD CLI process submissions
+++++++++++++++++++++++++++++++++++++++++++++++++++++

First let's take a look at what a process submitted using a gromacs specific CLI command compares with submitting the job with the genericMD CLI:

Submitting a ``pdb2gmx`` process using ``gmx_pdb2gmx``::

    gmx_pdb2gmx -f 1AKI_clean.pdb -ff oplsaa -water spce -o 1AKI_forcefield.gro -p 1AKI_topology.top -i 1AKI_restraints.itp

Submitting the equivalent process with ``genericMD``

.. code-block:: bash

    genericMD --code gmx@localhost \
    --command "pdb2gmx -i 1AKI_restraints.itp -o 1AKI_forcefield.gro -p 1AKI_topology.top -ff oplsaa -water spce -f 1AKI_clean.pdb" \
    --inputs 1AKI_clean.pdb \
    --outputs 1AKI_restraints.itp --outputs 1AKI_topology.top --outputs 1AKI_forcefield.gro

As you can see, using the genericMD CLI is more verbose, but it does allow for submitting any command you want to keep track of with AiiDA.


How to submit a process with genericMD
++++++++++++++++++++++++++++++++++++++

To submit an AiiDA process for tracking input and output files of a generic command, there are several things you need to know beforehand:

#. what the full command you want to run is,
#. where is the program you're using in your command installed,
#. what are the input file names you're using in the command,
#. and what are the output file names outputted from your command.

Once you know this information, you can build a genericMD process submission. Let's take a look at the example below for a submission for tracking the inputs and outputs of a ``diff`` command::

    genericMD --code bash@localhost --command "diff file1.txt file2.txt > out.txt" --inputs file1.txt --inputs file2.txt --outputs out.txt

Here,

#. the command we want to run is ``diff file1.txt file2.txt > out.txt``, which is specified with the ``--command`` flag and wrapped in quotation marks
#. the diff utility is installed on the local computer and can be run via bash, this is communicated to AiiDA with ``--code bash@localhost``
#. each input file used in the diff command needs to be explicitly defined using the ``--inputs`` flag, this allows for AiiDA to keep track of these files
#. similarly for the output files, use the ``--outputs`` flag to define each output from the command, which allows AiiDA to track these files

A few things to consider when using genericMD, firstly the inputs and ouputs of a command need to be known before submitting the process.

.. warning::
    Any inputs/outputs not included in the genericMD submission with ``--inputs`` and ``--outputs`` flags repectively, will not be included as a node in the provenance graph!!
