================
Start/Stop AiiDA
================

Regardless of whether you have switched off AiiDA or restarted your computer. To start recording provenance with this plugin, you will need to start up your AiiDA instance. Likewise, once you have finished using AiiDA, it is best practice to shut things down gracefully. These steps also assume you have already followed the steps within :doc:`installation` and have already a fully working install of the toolchain.

Starting AiiDA
++++++++++++++

The first step in starting AiiDA is to activate your conda environment, for example:

.. code-block:: bash

    conda activate aiida-2.4.0

The next step is to start the aiida database:

.. code-block:: bash

    pg_ctl -D ~/.aiida/aiida_db -l ~/.aiida/logfile start
    rabbitmq-server -detached

Then start the AiiDA process daemon:

.. code-block:: bash

    verdi daemon start 2

You can then confirm all is well by checking the status of verdi:

.. code-block:: bash

    verdi status

That is it, you are now free to use your AiiDA toolchain.

Stopping AiiDA
++++++++++++++

It is best to stop your AiiDA instance gracefully than to simply close your VM or shutdown your computer, this protects against any issues that might corrupt your database.

Firstly stop the verdi process daemon:

.. code-block:: bash

    verdi daemon stop

Next stop the database process:

.. code-block:: bash

    pg_ctl -D ~/.aiida/aiida_db stop
    rabbitmqctl stop

Finally you can deactivate your conda environment:

.. code-block:: bash

    conda deactivate

That is it, you now have fully disabled the AiiDA toolchain.
