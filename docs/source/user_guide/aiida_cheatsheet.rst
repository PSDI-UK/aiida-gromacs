====================
AiiDA CLI Cheatsheet
====================

Here we will list some useful commands for working with AiiDA. Have a look `here <https://aiida.readthedocs.io/projects/aiida-core/en/latest/reference/command_line.html?highlight=verdi%20process%20list>`_ to explore all available AiiDA command line subcommands.

Monitoring Submitted Processes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

List all submitted processes:

.. code-block:: bash

    verdi process list -a

Delete a process node, identified by its primary key value ``<PK>``::

    verdi node delete <PK>

Delete multiple process nodes between primary key values ``<PK_1>`` and ``<PK_n>``::

    verdi node delete {<PK_1>..<PK_n>}

View all the nodes associated with a process node::

    verdi node show <PK>

View attributes of a process node (such as retrieved files and find the path on disk where outputs are stored temporarily)::

    verdi node attributes <PK>

From the node attributes output dictionary, you can find where the input and output files are temporarily stored for a process in the "remote_workdir" value.

Debugging
^^^^^^^^^

Any aiida errors are logged in ``.aiida/daemon/log/``.

If any changes to the plugin code are made, after an update for example, restart the daemon if it is already running to implement the code changes::

    verdi daemon restart --reset

To view details of a submitted process, such as the inputs, state, log messages, etc., use the following command::

    verdi process show <PK>

To view where in the source code an exception has occured if a calculation has failed::

        verdi process report <PK>

Sharing Data
^^^^^^^^^^^^

When you are ready to share data, the AiiDA database and accompanying files inputted and outputted in each process can be `archived <https://aiida.readthedocs.io/projects/aiida-core/en/latest/howto/share_data.html>`_ into a single file:

.. code-block:: bash

    verdi archive create --all archive_name.aiida

where the ``--all`` flag saves all the data in the AiiDA profile. To import an existing AiiDA archive file to a loaded profile:

.. code-block:: bash

    verdi archive import archive_name.aiida

Visualise Data Provenance
^^^^^^^^^^^^^^^^^^^^^^^^^

Visualise your submitted jobs as a provenance graph outputted in a ``.pdf`` file. Select the latest ``<PK>`` to include all previous nodes in the graph::

    verdi node graph generate <PK>

An example provenance graph for the first eight steps of the :ref:`lysozyme tutorial <tutorial>`, will look something like this:

.. image:: ../images/53.dot.png
   :width: 600
   :align: center


Plugin Specfic AiiDA Commands
+++++++++++++++++++++++++++++

The following commands are only available with the aiida-gromacs plugin.

Show Provenance on CLI
^^^^^^^^^^^^^^^^^^^^^^

Show a list of the commands run and the connected inputs/outputs associated with any processes that have been run using::

    verdi data provenance show

An example output on the command line will look like this:

    .. code-block :: bash

        Step 1.
            command: curl https://gpcrdb.org/structure/homology_models/pth2r_human_active_full/download_pdb -o ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb.zip
            executable: bash
            input files:

            output files:
                ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb.zip

        Step 2.
            command: unzip ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb.zip
            executable: bash
            input files:
                ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb.zip <-- from Step 1.
            output files:
                ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb.pdb

        Step 3.
            command: sed -i -e '1,217d;3502,4387d' ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb.pdb
            executable: bash
            input files:
                ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb.pdb <-- from Step 2.
            output files:
                ClassB1_pth2r_human_Active_AF_2022-08-16_GPCRdb.pdb
