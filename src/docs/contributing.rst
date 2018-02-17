Contributing
============

First off, thanks for taking the time to contribute!

What should I know before I get started?
----------------------------------------

This project is there to predict potential proteins using TAG codons as PYL amino-acid instead as STOP codon.

The project is developed on GitHub at `https://github.com/bebatut/PylProtPredictor <https://github.com/bebatut/PylProtPredictor>`_.

How can I contribute?
---------------------

Reporting mistakes or errors
****************************

The easiest way to start contributing is to file an issue to tell us about a spelling mistake or a factual error.

Your first content contribution
*******************************

Once you are feeling more comfortable, you can propose changes via Pull Request.

Indeed, to manage changes, we use `GitHub flow <https://guides.github.com/introduction/flow/>`_ based on Pull Requests:

1. `Create a fork <https://help.github.com/articles/fork-a-repo/>`_ of this repository on GitHub
2. Clone your fork of this repository to create a local copy on your computer
3. Create a new branch in your local copy for each significant change
4. Commit the changes in that branch
5. Push that branch to your fork on GitHub
6. Submit a pull request from that branch to the `master repository <https://github.com/bebatut/PylProtPredictor>`_
7. If you receive feedback, make changes in your local clone and push them to your branch on GitHub: the pull request will update automatically

For beginners, the GitHub interface will help you in the process of editing a file. It will automatically create a fork of this repository where you can safely work and then submit the changes as a pull request without having to touch the command line.

Tests
-----

The code of PylProtPredictor is covered by unit tests. They are run automatically on CircleCI. You can activate CircleCI on your fork to run also the tests automatically.

We also recommend to run them locally before pushing to GitHub with:

.. code-block:: bash

    $ make test

Documentation
-------------

Documentation about ENASearch is available online at `http://bebatut.fr/PylProtPredictor <http://bebatut.fr/PylProtPredictor>`_

To update it:

1. Make the changes in `src/docs`
2. Generate the doc with 

  .. code-block:: bash

    $ make doc

3. Check it by opening the `docs/index.html` file in a web browser
4. Propose the changes via a Pull Request
