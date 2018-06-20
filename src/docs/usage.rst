Usage
=====

To run the different steps of the workflow to identify potential PYL proteins:

.. code-block:: bash

    $ source activate PylProtPredictor # once to activate the conda environment
    $ pylprotpredictor --help
    usage: pylprotpredictor [-h] --genome GENOME --output OUTPUT
                            [--reference_fasta_db REFERENCE_FASTA_DB]
                            [--reference_dmnd_db REFERENCE_DMND_DB]

    PylProtPredictor Pipeline

    optional arguments:
      -h, --help            show this help message and exit
      --genome GENOME       path to a FASTA file with full or contig sequences of
                            a genome to analyze
      --output OUTPUT       path to the output directory
      --reference_fasta_db REFERENCE_FASTA_DB
                            path to FASTA file with reference database
      --reference_dmnd_db REFERENCE_DMND_DB
                            path to Diamond formatted file with reference database

Database setup
--------------

The first run will be long: the reference database should be downloaded and prepare for the similarity search.

If you already have the Uniref90 database on your machine, you can simply link it when running the main script.

Otherwise, the pipeline will download and format it. Make sure you have at least 25GB available for the reference database. It can take several hours, depending on your connection.


Deactivate conda environment
----------------------------

Once done, to exit the environment, you can execute

.. code-block:: bash

    $ source deactivate

But don't do that before running the analysis.