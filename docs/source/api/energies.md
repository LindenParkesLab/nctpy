# Control energy

Control inputs, state trajectories and energies for linear network dynamics, in continuous and discrete time.
Normalise the connectome with {func}`nctpy.utils.matrix_normalization` first.

## Optimal control

```{eval-rst}
.. currentmodule:: nctpy.energies

.. autosummary::
   :toctree: generated/
   :nosignatures:

   get_control_inputs
   integrate_u
   sim_state_eq
```

## Minimum energy and Gramians

```{eval-rst}
.. currentmodule:: nctpy.energies

.. autosummary::
   :toctree: generated/
   :nosignatures:

   minimum_energy_fast
   gramian
   minimum_energy_infinite
   average_energy_infinite
```
