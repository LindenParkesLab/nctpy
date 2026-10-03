# Optimising decay rates

{func}`~nctpy.optimize.optimize_decay_rates` implements the method of Kim et al. (*Nat Commun* 2025): it fits each
node's decay rate, the diagonal of the continuous-time system, to a state transition. The
{doc}`/tutorials/decay_rates` tutorial shows it in use; this page describes what it does and the choices it
offers. It needs PyTorch: `pip install "nctpy[optimize]"`.

## What is fitted

The connectome is normalised for continuous time, and a decay term is subtracted from its diagonal (see
{doc}`normalisation`):

$$
A_\text{norm} = \frac{A}{\lambda_{\max} + c} - \operatorname{diag}(\text{decay}).
$$

Training adjusts the decay rates by gradient descent (Adam) to minimise the sum of three terms:

- the distance between the optimal trajectory at the midpoint of the time horizon, $x(T/2)$, and the reference
  state `xr` (by default the target state);
- a stability penalty, which grows without bound as any eigenvalue of $A_\text{norm}$ approaches zero, so the
  fitted system stays stable;
- a small penalty on the size of the diagonal.

The decay rates are fitted to reduce the distance to the reference state, not the energy; the energy of the
transition is computed afterwards, with the fitted matrix. The defaults (1000 steps at most, a learning rate of
0.01, the penalty weights, and every decay rate starting at 2) are those of Kim et al.

```python
from nctpy.optimize import optimize_decay_rates

fit = optimize_decay_rates(A, x0, xf)
fit.decay      # the fitted decay rate of each node
fit.A_norm     # the fitted, normalised connectome, ready for get_control_inputs and the other functions
fit.loss       # the training traces, with fit.max_eigenvalue, fit.n_steps_run and fit.stopped_early
```

Use `fit.A_norm` directly. The same matrix can be rebuilt with
`matrix_normalization(A, system="continuous", c=c, zero_diagonal=True, decay=fit.decay)`. A value $v$ in Kim et
al.'s Fig. 2B corresponds to `decay = 1 - v`.

## Choices

- **Self-connections.** The method assumes a connectome without them, so that every node starts from the same
  decay rate. `zero_diagonal=True` (the default) removes any before fitting. With `zero_diagonal=False` a
  connectome's self-connections are kept; that use is not supported.
- **Several transitions.** States given as columns, `x0` and `xf` of shape $(N, k)$, fit one set of decay rates
  to all $k$ transitions, averaging the distance term over them.
- **The reference state, `xr`.** The target state by default, as Kim et al. recommend: with zero activity as the
  reference, their fits switched off the initial state instead of completing the transition. The other options
  are those of {func}`~nctpy.energies.get_control_inputs`, including any state you pass.
- **The control task.** `T`, `B`, `S` and `rho` are as in {doc}`control_tasks`; by default `T = 1` and every node
  is controlled and constrained.
- **Directed connectomes.** The default stability penalty, `eigenvalues="symmetric"`, reproduces Kim et al.'s
  computation, which is exact for an undirected connectome. For a directed one, `eigenvalues="general"` uses the
  matrix's true eigenvalues.
- **Stopping.** Training stops early once the loss and the largest eigenvalue have stopped changing (over the last
  20% of `n_steps`); `early_stopping=False` runs every step.
- **Device.** Training runs on a GPU when PyTorch finds one, otherwise on the CPU; `device="cpu"` or
  `device="cuda"` chooses. Results agree between devices to floating-point precision, and runs on the CPU are
  reproducible exactly.

Fitting decay rates is continuous-time only. It changes $A$; {class}`~nctpy.pipelines.ComputeOptimizedControlEnergy`
instead optimises the control weights, $B$ (see {doc}`pipelines`).

nctpy reproduces Kim et al.'s fits: see {ref}`reproducing`.
