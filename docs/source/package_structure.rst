Python library structure
========================

``gencdna`` is a Python library that implements several components of best
coding practices:

Documentation
-------------

Documentation is built using Sphinx and deployed on Github using Conitnuous
Deployment (CD) using Github Actions. The ``doc/source`` folder contains the
restructuredText files with the documentation code. To create a local 
``doc/build`` folder with actual HTML files, use:

.. code-block:: bash

	make doc

If you want to consume the documentation only online, it will be automatically
generated on each Push and available at: https://igor-sb.github.io/gencdna/



Unit tests
----------

Code for unit testing is in the folder ``tests/``. Each subfolder is dedicated
to test a specific component. Inside each subfolder there are 4 types of files:

* fixtures: these contain manually or programmatically data 
* snapshots: contain snapshots for unit tests that use snapshot testing
* conftest.py: contains code that loads automatically in each ``test_*.py`` file
* ``test_*.py`` files: these are files with code for tests 

Unit tests can be executed using:

.. code-block:: bash

	make test

As number of tests are implemented as snapshot tests, the snapshots can be
updated using:

.. code-block:: bash

	pytest --snapshot-update tests/<folder>/<filename.py>


Integration tests (not finalized / updated)
-------------------------------------------

A full integration test, using mock dataset is performed using Nextflow and
deployed on Github actions. The Nextflow script containing the integration
test code is in ``workflows/integration-test-qc.nf`` (although for PCR data).

Code-coverage
-------------

On each push event, code coverage is uploaded to codecov and a report will be
visible as a comment during a Pull Request.


Dependency pinning
------------------

``gencdna`` library uses ``poetry`` as a dependency management and ensures
reproducibility using a lock file. When in a root folder of this repo, you
can always run ``poetry shell`` to enter a Python virtual environment that
contains exact version of libraries from the lock file.

Linting and type checking
-------------------------

Code is linted using `flake8 <https://flake8.pycqa.org/en/latest/>`_ with few
linting violations ignored. These are listed in the ``setup.cfg`` file in the
root folder of this repo. Python type hints are used for type checking using
`mypy <https://github.com/python/mypy>`_. 
