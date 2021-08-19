network_control: A toolbox for implementing Network Control Theory analyses in python
=====================================================================================

Overview
--------
.. image:: https://zenodo.org/badge/370716682.svg
   :target: https://zenodo.org/badge/latestdoi/370716682ß
.. image:: https://readthedocs.org/projects/control-package/badge/?version=latest
   :target: https://control-package.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/pypi/l/ansicolortags.svg
   :target: https://pypi.python.org/pypi/ansicolortags/

Network Control Theory (NCT) is a branch of physical and engineering sciences that treats a network as a dynamical
system. Generally, the system is controlled through signals that originate at a control point (or control points) and
move through the network. In the brain, NCT models each region’s activity as a time-dependent internal state that is
predicted from a combination of three factors: (i) its previous state, (ii) whole-brain structural connectivity,
and (iii) external inputs. NCT enables asking a broad range of questions of a networked system that are highly relevant
to network neuroscientists, such as: which regions are positioned such that they can efficiently distribute activity
throughout the brain to drive changes in brain states? Do different brain regions control system dynamics in different
ways? Given a set of control nodes, how can the system be driven to specific target state, or switch between a pair of
states, by means of internal or external control input?

``network_control`` is a Python toolbox that provides researchers with a set of tools to conduct some of the
common NCT analyses reported in the literature. Below, we list select publications that serve as a primer for
these tools and their use cases:

1. Gu, S., Pasqualetti, F., Cieslak, M. et al. Controllability of structural brain networks.
Nature Communications (2015). https://doi.org/10.1038/ncomms9414

2. Gu, S., Betzel R. F., Mattar, M. G. et al. Optimal trajectories of brain state transitions.
NeuroImage (2017). https://doi.org/10.1016/j.neuroimage.2017.01.003

3. Karrer, T. M., Kim, J. Z., Stiso, J. et al. A practical guide to methodological considerations in the
controllability of structural brain networks.
Journal of Neural Engineering (2020). https://doi.org/10.1088/1741-2552/ab6e8b

4. Kim, J. Z., & Bassett, D. S. Linear dynamics & control of brain networks.
arXiv (2019). https://arxiv.org/abs/1902.03309

.. _readme_requirements:

Requirements
------------

Currently, ``network_control`` works with Python 3.6 and requires the following dependencies:

    - numpy (tested on 1.19.5), and
    - scipy (tested on 1.5.4)

There are some additional (optional) dependencies you can install (note, these are only used for i/o and plotting and
only need to be installed if you want to run the notebooks in ``network_control/tests/``):

    - pandas (tested on 1.1.5)
    - statsmodels (tested on 0.12.2)
    - matplotlib (tested on 3.3.4)
    - seaborn (tested on 0.11.1), and
    - jupyterlab (tested on 3.0.16)


Basic installation
------------------

Assuming you have Python 3.6 installed, you can install ``network_control`` by opening a terminal and running
the following:

.. code-block:: bash

    pip install network_control

What's New
----------
    - v0.0.4: initial package release for network_control
    - v0.0.5: added options for specifying timescales, integration function, gramian function, better estimation of the number of time points, and optimal energy example


Questions
---------

If you have any questions, please contact Linden Parkes (https://twitter.com/LindenParkes), Jennifer Stiso
(https://twitter.com/JenniferStiso) or Jason Kim (https://twitter.com/jason_z_kim).
For questions or clarification about the theory pages of the documentation, please contact Jason Kim
(https://twitter.com/jason_z_kim).
