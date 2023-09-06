=================
Lysozyme tutorial
=================

This tutorial follows Justin Lemkul's `lysozyme <http://www.mdtutorials.com/gmx/lysozyme/>`_ tutorial. We will not explain each individual step as this can be found on Justin's webpage, but we will link to each page and show the AiiDA equivalant command.

Please also note the slight differences in commands between the tutorial and that by Justin Lemkul is simply down to the way we are recording provenance requires non-interactive input into the gromacs tools.

Also at each of the below steps you should run verdi to view the status of the submitted process before moving onto the next step, you do this by::

    verdi process list -a

A successfully finished process will exit with code ``[0]``.

#. We will start from the `pbd2gmx <http://www.mdtutorials.com/gmx/lysozyme/01_pdb2gmx.html>`_ step of Justin's tutorial::

    gmx_pdb2gmx -f 1AKI_clean.pdb -ff oplsaa -water spce -o 1AKI_forcefield.gro -p 1AKI_topology.top -i 1AKI_restraints.itp

#. Next we will `create the box and solvate <http://www.mdtutorials.com/gmx/lysozyme/03_solvate.html>`_

    Firstly the box::

        gmx_editconf -f 1AKI_forcefield.gro -center 0 -d 1.0 -bt cubic -o 1AKI_newbox.gro

    Then solvate::

        gmx_solvate -cp 1AKI_newbox.gro -cs spc216.gro -p 1AKI_topology.top -o 1AKI_solvated.gro

#. Add `ions <http://www.mdtutorials.com/gmx/lysozyme/04_ions.html>`_

    Firstly we will use the grompp preprocessor::

        gmx_grompp -f ions.mdp -c 1AKI_solvated.gro -p 1AKI_topology.top -o 1AKI_ions.tpr

    Followed by genion::

        gmx_genion -s 1AKI_ions.tpr -p 1AKI_topology.top -pname NA -nname CL -neutral true -o 1AKI_solvated_ions.gro

#. Then `minimise <http://www.mdtutorials.com/gmx/lysozyme/05_EM.html>`_ the structure

    Firstly we will use grompp::

        gmx_grompp -f min.mdp -c 1AKI_solvated_ions.gro -p 1AKI_topology.top -o 1AKI_minimised.tpr

    Then mdrun to minimise::

        gmx_mdrun -s 1AKI_minimised.tpr -c 1AKI_minimised.gro -e 1AKI_minimised.edr -g 1AKI_minimised.log -o 1AKI_minimised.trr

#. Now we will equilibrate with `NVT <http://www.mdtutorials.com/gmx/lysozyme/06_equil.html>`_

    Firstly we will use grompp::

        gmx_grompp -f nvt.mdp -c 1AKI_minimised.gro -r 1AKI_minimised.gro -p 1AKI_topology.top -o 1AKI_nvt.tpr

    Then mdrun to equilibrate NVT::

        gmx_mdrun -s 1AKI_nvt.tpr -c 1AKI_nvt.gro -e 1AKI_nvt.edr -g 1AKI_nvt.log -cpo 1AKI_nvt.cpt -o 1AKI_nvt.trr

#. Followed by equilibration with `NPT <http://www.mdtutorials.com/gmx/lysozyme/07_equil2.html>`_

    Firstly we will use grompp::

        gmx_grompp -f npt.mdp -c 1AKI_nvt.gro -r 1AKI_nvt.gro -t 1AKI_nvt.cpt -p 1AKI_topology.top -o 1AKI_npt.tpr

    Then mdrun to equilibrate NPT::

        gmx_mdrun -s 1AKI_npt.tpr -c 1AKI_npt.gro -e 1AKI_npt.edr -g 1AKI_npt.log -cpo 1AKI_npt.cpt -o 1AKI_npt.trr

#. We are now ready for `production <http://www.mdtutorials.com/gmx/lysozyme/08_MD.html>`_ MD.

    Firstly we will use grompp::

        gmx_grompp -f md.mdp -c 1AKI_npt.gro -t 1AKI_npt.cpt -p 1AKI_topology.top -o 1AKI_prod.tpr

    Then mdrun for production run::

        gmx_mdrun -s 1AKI_prod.tpr -c 1AKI_production.gro -e 1AKI_production.edr -g 1AKI_production.log -o 1AKI_production.trr

    If running on GPU then something like::

        gmx_mdrun -s 1AKI_prod.tpr -c 1AKI_production.gro -e 1AKI_production.edr -g 1AKI_production.log -o 1AKI_production.trr -bonded gpu -nb gpu -pme gpu -ntmpi 1 -ntomp 5 -pin on

That is it!
