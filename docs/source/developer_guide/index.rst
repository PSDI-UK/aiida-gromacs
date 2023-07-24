===============
Developer guide
===============

To get up and running as a developer contributing to this plugin, you will need to have a working installation of AiiDA installed on your computer. If you are using conda, you will need this activated when running the below commands.

Running the tests
+++++++++++++++++

To run our test suite, in addition to having AiiDA installed, you should also have the daemon running, which will include the full AiiDA startup steps for the database and an initial user profile to be setup.

The following will discover and run all unit test::

    pip install -e .[testing]
    pytest -v

Automatic coding style checks
+++++++++++++++++++++++++++++

We are implemting automated coding style checks, so that we are all commiting code to similar code quality standards.

Enable enable automatic checks of code sanity and coding style::

    pip install -e .[pre-commit]
    pre-commit install

After this, the `yapf <https://github.com/google/yapf>`_ formatter,
the `pylint <https://www.pylint.org/>`_ linter
and the `pylint <https://www.pylint.org/>`_ code analyzer will
run at every commit.

If you ever need to skip these pre-commit hooks, just use::

    git commit -n

Though, for pull requests to be accepted, we will expect these to have been resolved before pulling.


Continuous integration
++++++++++++++++++++++

``aiida-gromacs`` comes with a ``.github`` folder that contains continuous integration tests on every commit using `GitHub Actions <https://github.com/features/actions>`_. It will:

#. run all tests
#. build the documentation
#. check coding style and version number (not required to pass by default)

We have these activated on github via the github actions platform. When version numbers are tagged, we will also automatically push a version to pypi.

Building the documentation
++++++++++++++++++++++++++

 #. Install the ``docs`` extra::

        pip install -e .[docs]

 #. Edit the individual documentation pages::

        docs/source/index.rst
        docs/source/developer_guide/index.rst
        docs/source/user_guide/index.rst
        docs/source/user_guide/get_started.rst
        docs/source/user_guide/tutorial.rst

 #. Use `Sphinx`_ to generate the html documentation::

        cd docs
        make

Check the result by opening ``build/html/index.html`` in your browser.

.. _Sphinx: https://www.sphinx-doc.org/en/master/
