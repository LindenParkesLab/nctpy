(upgrading)=
# Upgrading

Within 1.x, nctpy changes only by addition: every code block printed in the protocol paper and its Supplementary
Information keeps running and gives the printed results. This page lists what to check when you upgrade. The
{doc}`changelog` has every change in detail.

## To 1.2, from 1.1

### Results that can change

- **Average controllability of directed connectomes.** {func}`~nctpy.metrics.ave_control` now follows Gu et al.
  (2015) for directed connectomes in both time systems: node $i$'s value is the trace of the controllability
  Gramian with input at node $i$ alone. Earlier versions computed a different quantity for directed connectomes.
  **Undirected connectomes are unaffected** (continuous time: identical; discrete time: identical to rounding). The
  protocol paper uses `ave_control` only on undirected connectomes. See {doc}`/guide/metrics`.
- **`ComputeOptimizedControlEnergy` and `xr`.** The class now uses a task's optional `'xr'`, as
  `ComputeControlEnergy` does. Tasks without `'xr'`, including the protocol paper's, give unchanged results.
- **Lower-precision input.** nctpy now computes in float64 throughout. Results for float64 input are unchanged;
  float32 input gives slightly different values than before.
- **Surface plots.** {func}`~nctpy.plotting.roi_to_vtx` leaves unlabelled vertices (label −1) at 0, like the medial
  wall; they used to take another parcel's value. The protocol paper's figures are unchanged.
- **Rounding only.** `ComputeControlEnergy` now solves transitions that share a system together; completed
  transitions agree with earlier versions to rounding, and transitions that did not complete give exactly the same
  values. `modal_control` and `gramian` also change at the level of rounding.

### Errors where earlier versions returned a value

Input that does not define a problem now raises an error with a clear message. Transitions that are merely
ill-conditioned or incomplete still return their results (see {doc}`/guide/errors`).

| Call | 1.2 | Before |
| --- | --- | --- |
| `get_control_inputs(..., rho=0)` (or any `rho <= 0`), and pipeline tasks with `rho = 0` | `ValueError` | NaN energies and errors |
| `gramian(...)` without a valid `system` | raises, like the other functions | returned `None` |
| `get_control_inputs(..., xr="target")` (an unknown string) | `ValueError` naming the options | `TypeError` from numpy |
| `get_null_p(..., version="other")` | `ValueError` | `UnboundLocalError` |
| `get_control_inputs` with several states, e.g. `x0` of shape `(N, k)` | `ValueError` saying it takes one state | `ValueError` from numpy |

### New in 1.2

- {func}`~nctpy.optimize.optimize_decay_rates`, fitting nodes' decay rates (Kim et al. 2025), with
  `pip install "nctpy[optimize]"`.
- {func}`~nctpy.energies.minimum_energy_infinite` and {func}`~nctpy.energies.average_energy_infinite`.
- {func}`~nctpy.utils.mask_control_set` and {func}`~nctpy.utils.random_control_set`.
- `matrix_normalization(..., zero_diagonal=...)` and `matrix_normalization(..., decay=...)`.
- `nctpy.null_models`, the null models importable from within nctpy as well as from `null_models`.
- Much faster loops over transitions, and type annotations.

## To 1.1, from 1.0

- **Install the plotting libraries explicitly.** They are no longer installed with nctpy. To run the protocol
  paper's code, install with `pip install "nctpy[paper]"`; for `nctpy.plotting` alone, `pip install "nctpy[plot]"`.
  Without them, importing `nctpy.plotting` raises an `ImportError` that says so.
- **Python 3.10 or later** is required.
- `surface_plot` works with nilearn 0.13 and later.

## From 1.0.2 or earlier

- Continuous-time `ave_control` (Procedure 2 of the protocol paper) raises an `AttributeError` on SciPy 1.14 or
  later in these versions; 1.0.3 fixed it.
- In these versions, `geomsurr` set the diagonal of the matrix passed to it to zero, in place, and reset numpy's
  global random seed. After running the Supplementary Information's null-model code, a connectome with
  self-connections had silently lost them. 1.0.3 fixed it; the surrogates themselves are unchanged.
