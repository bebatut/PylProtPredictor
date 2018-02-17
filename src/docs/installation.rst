Installation
============

Via source code
---------------

Requirements
************

The following software are required:
- [`git`](https://git-scm.com/book/fr/v1/D%C3%A9marrage-rapide-Installation-de-Git#Installation-sur-Linux)
- [`conda`](https://conda.io/miniconda.html):

    .. code-block:: bash

        $ make install-conda
        $ make configure-conda

Install the tool
****************

- Clone this repository (or get the release)

.. code-block:: bash

    $ git clone https://github.com/bebatut/PylProtPredictor.git


- Move into the directory

.. code-block:: bash

    $ cd pyl_protein_prediction

- Prepare the environment (only once)

.. code-block:: bash

    $ make create-env
