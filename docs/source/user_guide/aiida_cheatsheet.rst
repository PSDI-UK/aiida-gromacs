================
AiiDA Cheatsheet
================

Here we will list some useful commands for working with AiiDa.

Have a look `here <https://aiida.readthedocs.io/projects/aiida-core/en/latest/reference/command_line.html?highlight=verdi%20process%20list>`_ to explore all AiiDA command line subcommands.

Monitoring submitted processes
++++++++++++++++++++++++++++++

List all submitted processes :

.. code-block:: bash

    verdi process list -a

Delete a process node, identified by its primary key value ``<PK>``:

.. code-block:: bash

    verdi node delete <PK>

Delete multiple process nodes between primary key values ``<PK_1>`` and ``<PK_n>``:

.. code-block:: bash

    verdi node delete {<PK_1>..<PK_n>}

View all the nodes associated with a process node:

.. code-block:: bash

    verdi node show <PK>

View attributes of a process node (such as retrieved files and find the path on disk where outputs are stored temporarily):

.. code-block:: bash

    verdi node attributes <PK>

Visualise your submitted jobs as a provenance graph outputted as a ``.pdf`` file. Select the latest ``<PK>`` to include all previous nodes in the graph:

.. code-block:: bash

    verdi node graph generate <PK>


Debugging
+++++++++


Any aiida errors are logged in ``.aiida/daemon/log/``.

If any changes to the plugin code are made, after an update for example, restart the daemon if it is already running to implement the code changes:

.. code-block:: bash

    verdi daemon restart --reset
