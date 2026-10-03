nctpy: Network Control Theory for Python
=====================================================================================

Overview
--------
.. image:: https://zenodo.org/badge/370716682.svg
   :target: https://zenodo.org/badge/latestdoi/370716682
.. image:: https://readthedocs.org/projects/nctpy/badge/?version=latest
    :target: https://nctpy.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status
.. image:: https://img.shields.io/pypi/l/ansicolortags.svg
   :target: https://pypi.python.org/pypi/ansicolortags/

Network Control Theory (NCT) is a branch of physical and engineering sciences that treats a network as a dynamical
system. Generally, the system is controlled through `control signals` that originate at a control node (or control nodes) and
move through the network. In the brain, NCT models each region’s activity as a time-dependent internal state that is
predicted from a combination of three factors: (i) its previous state, (ii) whole-brain structural connectivity,
and (iii) external inputs. NCT enables asking a broad range of questions of a networked system that are highly relevant
to network neuroscientists, such as: which regions are positioned such that they can efficiently distribute activity
throughout the brain to drive changes in brain states? Do different brain regions control system dynamics in different
ways? Given a set of control nodes, how can the system be driven to a specific target state, or switch between a pair of
states, by means of internal or external control input?

``nctpy`` is a Python toolbox that provides researchers with a set of tools to conduct some of the
common NCT analyses reported in the literature. It implements the methods of two papers:

1. Parkes, L., Kim, J. Z., et al. A network control theory pipeline for studying the dynamics of the structural
connectome. Nature Protocols 19, 3721–3749 (2024). https://doi.org/10.1038/s41596-024-01023-w

2. Kim, J. Z., ..., Parkes, L. Inferring intrinsic neural timescales using optimal control theory.
Nature Communications 16, 11639 (2025). https://doi.org/10.1038/s41467-025-66542-w

The following publications serve as a primer for these tools and their use cases:

3. Karrer, T. M., Kim, J. Z., Stiso, J. et al. A practical guide to methodological considerations in the
controllability of structural brain networks.
Journal of Neural Engineering (2020). https://doi.org/10.1088/1741-2552/ab6e8b

4. Kim, J. Z., & Bassett, D. S. Linear dynamics & control of brain networks.
arXiv (2019). https://arxiv.org/abs/1902.03309

.. _readme_requirements:

Requirements
------------

``nctpy`` requires Python 3.10 or later. Installing it also installs its core dependencies: numpy, scipy,
tqdm and statsmodels.

The ``plotting`` module needs extra packages (matplotlib, seaborn, nibabel and nilearn), which are
optional and not installed by default. They come with the ``plot`` extra, shown below.

If you want to install the environment that was used to run the analyses presented in the manuscript, use the
environment.yml file.

Basic installation
------------------

To install ``nctpy``, open a terminal and run:

.. code-block:: bash

    pip install nctpy

To also install the plotting dependencies:

.. code-block:: bash

    pip install "nctpy[plot]"

To run the code printed in our Nature Protocols paper, which also uses pandas and scikit-learn:

.. code-block:: bash

    pip install "nctpy[paper]"

To fit nodes' decay rates with ``nctpy.optimize`` (Kim et al., Nature Communications 2025), which needs
PyTorch:

.. code-block:: bash

    pip install "nctpy[optimize]"

Questions
---------

If you have any questions, please contact Linden Parkes and Jason Kim: info@parkeslab.com.
