Usage
=====

To run the different steps of the workflow to identify potential PYL proteins:

.. code-block:: bash

    $ source activate PylProtPredictor # once to activate the conda environment
    $ ./bin/PylProtPredictor --genome FILE --output PATH [options]

Database setup
--------------

The first run will be long: the reference database should be downloaded and prepare for the similarity search.

If you already have the Uniref90 database on your machine, you can simply symlink it.

.. code-block:: bash

    $ ln -s /path/to/uniref90.dmnd data/uniref90.dmnd

Otherwise, the pipeline will download and format it. Make sure you have at least 25GB available for the reference database. It can take several hours, depending on your connection.


Deactivate conda environment
----------------------------

Once done, to exit the environment, you can execute

.. code-block:: bash

    $ source deactivate

But don't do that before running the analysis.