========
Setup and Use of the AiiDA-GROMACS Plugin
========

This page goes through the steps for installing the required packages to use the GROMACS plugin for AiiDA.

Python Virtual Environment
++++++++++++++++++

We recommend to set up a Python virtual environment via Conda, which can be installed by downloading the relevant installer `here <https://docs.conda.io/en/latest/miniconda.html>`_.
If using Linux or Mac OS, install conda via the terminal with,

.. code-block:: bash

    bash Miniconda3-latest-MacOSX-arm64.sh

Then add the conda path to the bash environment by appending the following to your ``.bashrc`` file,

.. code-block:: bash

    export PATH="~/miniconda3/bin:$PATH"

Options for AiiDA Installation
++++++++++++++++++++++++++++++

Our AiiDA plugin has been tested with AiiDA ``v2.2.2``, we recommend to `install <https://aiida.readthedocs.io/projects/aiida-core/en/v2.2.2/intro/install_conda.html#intro-get-started-conda-install>`_ this version of AiiDA in a conda environment. If you are using a linux OS, execute the following in the terminal, which installs AiiDA via an initial mamba installation

.. code-block:: bash

    conda install -c conda-forge mamba
    mamba create --name aiida-2.2.2 -c conda-forge aiida-core=2.2.2 aiida-core.services=2.2.2

This installation method may not work on Macs with M2 chips, if this is the case, install directly from conda

.. code-block:: bash

    conda create -yn aiida-2.2.2 -c conda-forge aiida-core=2.2.2

Plugin Installation
+++++++++++++++++++

To install the AiiDA-gromacs plugin, activate the conda environment created previously and install our plugin via Pip,

.. code-block:: bash

    conda activate aiida-2.2.2
    pip install aiida-gromacs

GROMACS Installation
++++++++++++++++++++

If GROMACS is not already installed, here's a quick installation guide. Our plugin has been tested with GROMACS ``v22.4`` and we suggest installation of this version. GROMACS requires the latest version of cmake, download the relevant `cmake installer <https://cmake.org/download/>`_ and install via the terminal with,

.. code-block:: bash

    bash cmake-3.27.2-linux-x86_64.sh

And include the path to cmake in the ``.bashrc`` file

.. code-block:: bash

    export PATH="~/make-3.27.2-linux-x86_64/bin:$PATH"

Download the relevant `GROMACS installer <https://manual.gromacs.org/documentation/>`_  and install via the `quick and dirty method <https://manual.gromacs.org/documentation/current/install-guide/index.html#>`_, summarised below

.. code-block:: bash

    tar xfz gromacs-2022.4.tar.gz
    cd gromacs-2022.4
    mkdir build
    cd build
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
    make
    make check
    sudo make install
    source /usr/local/gromacs/bin/GMXRC

Add the GROMACS path to the ``.bashrc`` file

.. code-block:: bash

    export PATH=/usr/local/gromacs/bin:$PATH

AiiDA Database Setup
++++++++++++++++++++

To start using AiiDA-gromacs to track the inputs and outputs of GROMACS commands, AiiDA first requires for a profile to be set up for each project via verdi. Verdi is a command line tool in AiiDA used to interact with the AiiDA database. All commands run via the AiiDA-GROMACS plugin are tracked and stored in the AiiDA database. Initialise the AiiDA database via

.. code-block:: bash

    initdb -D ~/.aiida/aiida_db
    pg_ctl -D ~/.aiida/aiida_db -l ~/.aiida/logfile start
    rabbitmq-server -detached

Then create a profile within the AiiDA database,

.. code-block:: bash

    verdi quicksetup
        Info: enter "?" for help
        Info: enter "!" to ignore the default and set no value
        Profile name: tutorial
        Email Address (for sharing data): your@email.com
        First name: Your
        Last name: Name
        Institution: where-you-work

This profile will be used in the example below that shows how to run some GROMACS commands with the files provided in ``aiida_gromacs/examples/gromacs_files``. 

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







