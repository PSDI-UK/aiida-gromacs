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











Once simulation setup is complete, the AiiDA database and accompanying files inputted and outputted in each process can be `archived <https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/share_data.html>`_ into a single file

.. code-block:: bash

    verdi archive create --all tutorial.aiida

where the ``--all`` flag saves all the data in the AiiDA profile. To import an existing AiiDA archive file to a loaded profile


.. code-block:: bash

    verdi archive import archive_name.aiida
