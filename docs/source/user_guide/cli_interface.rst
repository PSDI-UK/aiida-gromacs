=========================
Plugin on the Commandline
=========================

Before attempting to use the features here, you need to have completed the steps in :doc:`installation` and started the verdi daemon as described in :doc:`aiida_sessions`.

With our plugin, we are aiming for as simple a transition to using these tools as is possible, so the design philosophy we have chosen is to attempt to mimmick as much as is possible the functionality of native gromacs CLI tools. We currently support the following set of commandline tools and their full feature set. We also have implemented a generic cli interface to provide support for features that we currently do not have, or newly released in the gromacs application.

In general if with native gromacs you would normally run::

    gmx pdb2gmx -f prot.pdb -ff oplsaa -water spce -o gmx.gro -p topology.top

Then with our plugin you would run::

    gmx_pdb2gmx -f prot.pdb -ff oplsaa -water spce -o gmx.gro -p topology.top

Note the only difference is the underscore between the main gmx executable and the gmx application is the only difference. It is through this powerful interface that you will be able to build powerful workflows or use data provenance capturing tools in this plugin in your existing ways of working without cultural change in the way that you work!

The following utilities are available and have the following features.

genericmd
+++++++++

gmx_editconf
++++++++++++

gmx_genion
++++++++++

gmx_grompp
++++++++++

gmx_mdrun
+++++++++

gmx_pdb2gmx
+++++++++++

gmx_solvate
+++++++++++











Submitting GMX commands with the Plugin
+++++++++++++++++++++++++++++++++++++++

First initialise the AiiDA daemon, which manages process run in AiiDA, with 2 workers

.. code-block:: bash

    verdi daemon start 2

Run the command to generate a GROMACS topology file from a pdb file,

.. code-block:: bash

    gmx_pdb2gmx -f 1AKI_clean.pdb -ff oplsaa -water spce -o 1AKI_forcefield.gro -p 1AKI_topology.top -i 1AKI_restraints.itp

Use verdi to view the status of the submitted process

.. code-block:: bash

    verdi process list -a

A successfully finished process will exit with code ``[0]``. A list of subsequent commands to setup a solvated lysozyme system for a molecular dynamics simulations, can also be submitted

.. code-block:: bash

    gmx_editconf -f 1AKI_forcefield.gro -center 0 -d 1.0 -bt cubic -o 1AKI_newbox.gro
    gmx_solvate -cp 1AKI_newbox.gro -cs spc216.gro -p 1AKI_topology.top -o 1AKI_solvated.gro
    gmx_grompp -f ions.mdp -c 1AKI_solvated.gro -p 1AKI_topology.top -o 1AKI_ions.tpr
    gmx_genion -s 1AKI_ions.tpr -p 1AKI_topology.top -pname NA -nname CL -neutral true -o 1AKI_solvated_ions.gro
    gmx_grompp -f min.mdp -c 1AKI_solvated_ions.gro -p 1AKI_topology.top -o 1AKI_min.tpr
    gmx_mdrun -s 1AKI_min.tpr -c 1AKI_minimised.gro -e 1AKI_minimised.edr -g 1AKI_minimised.log -o 1AKI_minimised.trr

Once simulation setup is complete, the AiiDA database and accompanying files inputted and outputted in each process can be `archived <https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/share_data.html>`_ into a single file

.. code-block:: bash

    verdi archive create --all tutorial.aiida

where the ``--all`` flag saves all the data in the AiiDA profile. To import an existing AiiDA archive file to a loaded profile


.. code-block:: bash

    verdi archive import archive_name.aiida
