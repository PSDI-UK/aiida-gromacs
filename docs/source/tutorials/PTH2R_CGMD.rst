=================================================================
Advanced Data Provenance - Coarse-grained setup of embedded GPCR
=================================================================

Introduction
------------
In this tutorial, we will go through how to set-up a coarse-grained molecular system of a PTH2R (Parathyroid hormone receptor type 2) protein embedded in a lipid bilayer membrane along with water and counter-ions. We will use the command-line tools provided in aiida-gromacs to track each step performed on the terminal.

.. image:: ../images/gpcr_membrane.png
   :width: 600
   :align: center

Software and environment requirements
-------------------------------------
* `martinize2 <https://pypi.org/project/vermouth/>`__ is used to convert from atomistic to coarse-grained structures.

* `dssp <https://anaconda.org/salilab/dssp>`__ is used by ``martinize2`` to find secondary structures in the protein.

* A modified `insane <https://github.com/Tsjerk/Insane>`__ script is used to build the coarse-grained system.

* A local installation of `GROMACS <https://www.gromacs.org/>`__ is required to perform molecular dynamics simulations.

* `aiida-gromacs <https://aiida-gromacs.readthedocs.io/en/latest/user_guide/installation.html#plugin-installation>`__ is used to keep track of all the commands used to setup and perform the simulation.

Activate the conda environment that aiida-gromacs and AiiDA are installed in for this tutorial. We then need to start the AiiDA database and daemon using the first three steps from our `user guide <https://aiida-gromacs.readthedocs.io/en/latest/user_guide/aiida_sessions.html#start-stop-aiida>`__.


Aquiring and tidying up the receptor protein structure
------------------------------------------------------

Our starting point is the PTH2R `structure <https://gpcrdb.org/protein/pth2r_human/>`_ from the GPCRdb.

1. First, we download the PTH2R protein from the GPCRdb using ``curl``. We will track our commands as we perform them using the ``genericMD`` command in aiida-gromacs:

    .. code-block:: bash

        PTH2R="ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb"

        genericMD --code bash@localhost \
        --command "curl https://gpcrdb.org/structure/homology_models/pth2r_human_active_full/download_pdb -o {PTH2R}.zip " \
        --outputs {PTH2R}.zip

2. We then need to ``unzip`` the downloaded file:

    .. code-block:: bash

        genericMD --code bash@localhost \
        --command "unzip {PTH2R}.zip" \
        --inputs {PTH2R}.zip --outputs {PTH2R}.pdb

3. Now that we have the pdb file, we can remove regions of low confidence at the start and end of the receptor chain. We use the ``sed`` command to programmatically delete lines that correspond to low confidence regions between residues 1-31 and 435-550 inclusive.

    .. code-block:: bash

        genericMD --code bash@localhost \
        --command "sed -i -e '1,217d;3502,4387d' {PTH2R}.pdb" \
        --inputs {PTH2R}.pdb \
        --outputs {PTH2R}.pdb

Aligning PTH2R to a correctly orientated structure
--------------------------------------------------

4. Next, we download the correctly orientated structure from the OPM database.

    .. code-block:: bash

        genericMD --code bash@localhost \
        --command "curl https://opm-assets.storage.googleapis.com/pdb/7f16.pdb -o PTH2R_opm.pdb " \
        --outputs PTH2R_opm.pdb

5. PTH2R is a receptor for the parathyroid hormone and the OPM downloaded structure contains the coupled G-protein along with other bound molecules. We will keep only the receptor using the ``sed`` command to remove lines that do not correspond to the receptor:

    .. code-block:: bash

        genericMD --code bash@localhost \
        --command "sed -i -e '2,761d;3835,13708d' PTH2R_opm.pdb" \
        --inputs PTH2R_opm.pdb \
        --outputs PTH2R_opm.pdb

The final step for preparing the PTH2R protein is to position the structure in the correct orientation by aligning against the structure downloaded from the OPM database. There are a few ways to orientate the protein, here we use the orientations of proteins in membranes (OPM) database `structure <https://opm.phar.umich.edu/proteins/7900>`__ as a template to align our protein with. The OPM structure is correctly orientated to fit around a membrane and uses the `7F16 <https://www.rcsb.org/structure/7F16>`_ PDB deposited structure, however, this structure has missing atoms, so we cannot use this structure directly. There is an option to use the PPM webserver to orientate the protein correctly, however, no command-line tool is currently available, so we will not use this here.

6. We use the ``confrms`` command in GROMACS to align our structure. We will carry on using ``genericMD`` to track this command and use the ``echo`` command to include the interactive options required by ``confrms``:

    .. code-block:: bash

        genericMD --code bash@localhost \
        --command "echo -e '0 | 0 \n q' | gmx confrms -f1 PTH2R_opm.pdb -f2 {PTH2R}.pdb -name -one -o PTH2R_fit.pdb" \
        --inputs PTH2R_opm.pdb --inputs {PTH2R}.pdb \
        --outputs PTH2R_fit.pdb

Building a coarse-grained system from an atomic structure
---------------------------------------------------------

Now that we have the correct starting structure of the receptor, we move onto coarse-graining.

7. We use Martinize2 to coarse-grain the atomistic structure and produce a GROMACS topology file

    .. code-block:: bash

        genericMD --code martinize2@localhost --command "-f PTH2R_fit.pdb -o PTH2R_opm.top -x PTH2R_opm.cg.pdb -ff martini3001 -nt -dssp mkdssp -elastic -p backbone -maxwarn 1 -mutate HSD:HIS -mutate HSP:HIH -ignh -cys auto -scfix" \
        --inputs PTH2R_fit.pdb \
        --outputs PTH2R_opm.top --outputs PTH2R_opm.cg.pdb --outputs molecule_0.itp


Adding the membrane and solution around the protein with ``insane``
-------------------------------------------------------------------

8. Next, we use our custom `insane.py <https://github.com/PSDI-UK/aiida-gromacs/blob/master/examples/PTH2R_coarse-grained_files/insane/insane_custom.py>`__ python script to embed the protein into a lipid bilayer and solvate the system. Our insane script is modified from `the Melo lab <https://github.com/MeloLab/PhosphoinositideParameters/blob/main/martini3/insane.py>`_, it has been updated to python3 and contains additional parameters for the GM3 carbohydrate.

    .. code-block:: bash

        genericMD --code python@localhost --command "insane_custom.py -f PTH2R_opm.cg.pdb -o solvated.gro -p system.top -pbc rectangular -box 18,18,17 -u POPC:25 -u DOPC:25 -u POPE:8 -u DOPE:7 -u CHOL:25 -u DPG3:10 -l POPC:5 -l DOPC:5 -l POPE:20 -l DOPE:20 -l CHOL:25 -l POPS:8 -l DOPS:7 -l POP2:10 -sol W" \
        --inputs insane_custom.py --inputs PTH2R_opm.cg.pdb \
        --outputs solvated.gro --outputs system.top

Preparing the system for simulation
------------------------------------

9. Once the topology file is created, we need to include all the itp files containing the force field parameters used to describe interactions between beads. We use the ``sed`` command again to edit the ``system.top`` file directly on the command-line and we submit this command via ``genericMD`` as with the previous commands.

    .. code-block:: bash

        sed_command1='sed -i -e "1 s/^/#include \\"toppar\/martini_v3.0.0.itp\\"\\n#include \\"toppar\/martini_v3.0.0_ions_v1.itp\\"\\n#include \\"toppar\/martini_v3.0.0_solvents_v1.itp\\"\\n#include \\"toppar\/martini_v3.0.0_phospholipids_v1.itp\\"\\n#include \\"martini_v3.0_sterols_v1.0.itp\\"\\n#include \\"POP2.itp\\"\\n#include \\"molecule_0.itp\\"\\n#include \\"gm3_final.itp\\"\\n/" '\
        '-e "s/Protein/molecule_0/" '\
        '-e "s/#include \\"martini.itp\\"/\\n/" system.top'

        genericMD --code bash@localhost \
        --command '{sed_command1}' \
        --inputs system.top \
        --outputs system.top

10. We also need to edit the ``molecule_0.itp`` file generated from the Martinize2 step to include positional restraints on the coarse-grained beads.

        .. code-block:: bash

            sed_command2='sed -i -e "s/1000 1000 1000/POSRES_FC    POSRES_FC    POSRES_FC/g" '\
            '-e "s/#ifdef POSRES/#ifdef POSRES\\n#ifndef POSRES_FC\\n#define POSRES_FC 1000.00\\n#endif/" '\
            'molecule_0.itp'

            genericMD --code bash@localhost \
            --command '{sed_command2}' \
            --inputs molecule_0.itp \
            --outputs molecule_0.itp

11. Ions need to be added to neutralise the system and we can construct the GROMACS ``.tpr`` binary file containing the system configuration, topology and input parameters for the next step. We use the ``gmx_grompp`` command (note the underscore), which is wrapper command to run ``gmx`` via aiida-gromacs. We have included the most popular ``gmx`` commands in aiida-gromacs, the list of these are provided `here <https://aiida-gromacs.readthedocs.io/en/latest/user_guide/cli_interface.html>`_.

        .. code-block:: bash

            gmx_grompp -f ions.mdp -c solvated.gro -p system.top -o ions.tpr

12. The ``gmx_genion`` command is then used to add the ions to reach a particular salt concentration and neutralise the system. As the ``genion`` command requires interactive user inputs, we can provide these in as an additional text file via the ``--instructions`` argument. Each interactive response can be provided on a new line in the input text file. In this example, we replace solvent ``W`` with ions,

        .. code-block:: bash

            gmx_genion -s ions.tpr -o solvated_ions.gro -p system.top -pname NA -nname CL -conc 0.15 -neutral true --instructions inputs_genion.txt

    where `inputs_genion.txt <https://github.com/PSDI-UK/aiida-gromacs/blob/master/examples/PTH2R_coarse-grained_files/gromacs/inputs_genion.txt>`_ contains the following lines:

        .. code-block:: bash

            W

13. Lastly, we will use a ``gmx_make_ndx`` to create new index groups for the membrane and solute consituents

        .. code-block:: bash

            gmx_make_ndx -f solvated_ions.gro -o index.ndx --instructions inputs_index.txt

    where `inputs_index.txt <https://github.com/PSDI-UK/aiida-gromacs/blob/master/examples/PTH2R_coarse-grained_files/gromacs/inputs_genion.txt>`__ contains the following lines:

            .. code-block:: bash

                13|14|15|16|17|18|19|20|21
                name 26 membrane
                22|23|24
                name 27 solute
                q

We have built our starting configuration of an embedded protein in a lipid bilayer, hurray!
