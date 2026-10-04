# Optimising decay rates

Fit each node's decay rate to a state transition (Kim et al., *Nat Commun* 2025). The default optimiser, Adam, needs
PyTorch (`pip install 'nctpy[optimize]'`); the faster `method="lbfgs"` does not.

```{eval-rst}
.. currentmodule:: nctpy.optimize

.. autosummary::
   :toctree: generated/
   :nosignatures:

   optimize_decay_rates
   DecayRateFit
```
