# Choosing a time system

nctpy models the activity of a network's nodes, $\mathbf{x}(t)$, with one of two linear dynamical systems:

$$
\text{discrete time:}\quad \mathbf{x}(t+1) = A\mathbf{x}(t) + B\mathbf{u}(t)
\qquad\qquad
\text{continuous time:}\quad \frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = A\mathbf{x}(t) + B\mathbf{u}(t)
$$

where $A$ is the normalised connectome, $\mathbf{u}(t)$ the input to the system and $B$ the nodes that receive it.

## Which to choose

In discrete time, the state jumps from one time step to the next; in continuous time, it changes smoothly. The two
are not coarse- and fine-grained versions of one model: they behave differently. A one-node discrete system
$x(t+1) = -x(t)$ starting at 1 jumps between 1 and $-1$ without visiting any value in between, whereas a
continuous system can only get from 1 to $-1$ by passing through every value between them. The first is closer to
all-or-nothing events such as spikes, the second to the population-level activity that macroscale connectomes
describe.

The protocol paper therefore recommends **continuous time as the default** for macroscale connectomes, and
suggests replicating the primary results in continuous time if you use discrete time. Sampling fMRI at intervals
of about a second is not, on its own, a reason to use discrete time: sampling a continuous system at an interval
$\Delta t$ gives the discrete system matrix $e^{A\Delta t}$, whose eigenvalues are non-negative, whereas a
connectome normalised for discrete time generally has negative eigenvalues too, so it cannot be such a sample.

## What changes between them

Pass the time system to every function that asks for it, as `system="continuous"` or `system="discrete"`. nctpy
cannot tell from a matrix which system it was normalised for, and every function that needs to know raises an
error if `system` is missing.

- **Normalisation.** The two systems are stabilised differently (see {doc}`normalisation`), so a matrix normalised
  for one must not be used with the other.
- **The time horizon, `T`.** In continuous time it is a duration, in arbitrary units, which nctpy divides into
  steps of 0.001: a trajectory has `T / 0.001 + 1` time points. In discrete time it is a number of steps, an
  integer of at least 2: a trajectory has `T + 1` states and `T` inputs.
- **Some functions support only continuous time.**

| Function | Continuous | Discrete |
| --- | :---: | :---: |
| {func}`~nctpy.utils.matrix_normalization` | ✓ | ✓ (without `decay`) |
| {func}`~nctpy.energies.get_control_inputs`, {func}`~nctpy.energies.sim_state_eq` | ✓ | ✓ |
| {func}`~nctpy.energies.gramian`, {func}`~nctpy.energies.minimum_energy_infinite`, {func}`~nctpy.energies.average_energy_infinite` | ✓ | ✓ |
| {func}`~nctpy.energies.minimum_energy_fast` | ✓ | |
| {func}`~nctpy.metrics.ave_control` | ✓ | ✓ |
| {func}`~nctpy.metrics.modal_control` | | ✓ |
| {class}`~nctpy.pipelines.ComputeControlEnergy`, {class}`~nctpy.pipelines.ComputeOptimizedControlEnergy` | ✓ | ✓ |
| {func}`~nctpy.optimize.optimize_decay_rates` | ✓ | |

For more on the two kinds of system, see the {ref}`theory` pages and texts on linear systems theory.
