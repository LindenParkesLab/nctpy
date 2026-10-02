.. _installation_setup:

----------------------
Installation and setup
----------------------

.. _basic_installation:

Basic installation
==================

This package requires Python 3.10 or later. Assuming you have the correct version of Python installed, you can
install ``nctpy`` by opening a terminal and running the following:

.. code-block:: bash

    pip install nctpy

This installs ``nctpy`` with its core dependencies (numpy, scipy, tqdm and statsmodels). The
``nctpy.plotting`` module needs extra packages, which are optional:

.. code-block:: bash

    pip install "nctpy[plot]"    # matplotlib, seaborn, nibabel, nilearn
    pip install "nctpy[paper]"   # the above plus pandas and scikit-learn, to run the protocol paper's code
    pip install "nctpy[optimize]"   # PyTorch, for nctpy.optimize (fitting decay rates)

GitHub installation
===================

Alternatively, you can install the most up-to-date version of ``nctpy`` from GitHub:

.. code-block:: bash

   git clone https://github.com/LindenParkesLab/nctpy.git
   cd nctpy
   pip install .
