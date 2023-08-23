=========================
Plugin on the Commandline
=========================



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
