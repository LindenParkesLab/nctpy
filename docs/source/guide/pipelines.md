# Pipelines

{mod}`nctpy.pipelines` has two classes that wrap the steps of the protocol paper: they take the raw connectome,
normalise it with {func}`~nctpy.utils.matrix_normalization` when `run()` is called, and compute energy with
{func}`~nctpy.energies.get_control_inputs` and {func}`~nctpy.energies.integrate_u`. Their energies are on the same
scale as the protocol paper's (see {doc}`energies`).

## Many transitions: `ComputeControlEnergy`

{class}`~nctpy.pipelines.ComputeControlEnergy` computes the energy of a list of control tasks, one dict per
transition, and stores them in `E`, in the order of the tasks.

```python
from nctpy.pipelines import ComputeControlEnergy

control_tasks = []
for i in range(n_states):
    for j in range(n_states):
        control_tasks.append({"x0": normalize_state(states == i), "xf": normalize_state(states == j),
                              "B": np.eye(n_nodes), "S": np.eye(n_nodes), "rho": 1})

pipeline = ComputeControlEnergy(A=adjacency, control_tasks=control_tasks, system="continuous", c=1, T=1)
pipeline.run()
energy_matrix = pipeline.E.reshape(n_states, n_states)  # rows: initial states; columns: target states
```

Each task needs `"x0"`, `"xf"`, `"B"`, `"S"` and `"rho"`, and may set `"xr"` (default `"zero"`); see
{doc}`control_tasks`. Tasks can differ in any of these.

Tasks that share a control set, constraint and `rho` are solved together, which is several times faster than a
loop over `get_control_inputs`. Their energies agree with the loop's to rounding. Any transition that does not
complete (an error term of $10^{-8}$ or more) is solved again on its own, so its energy is exactly what
`get_control_inputs` returns. The class stores energies only, not the error terms: check representative
transitions with `get_control_inputs` (see {doc}`errors`).

## Optimising the control weights: `ComputeOptimizedControlEnergy`

{class}`~nctpy.pipelines.ComputeOptimizedControlEnergy` follows the protocol paper's Supplementary Information: for
one control task, it adjusts the weight of each node in a full control set to reduce the transition's energy.
Starting from equal weights ($B = I$), each step estimates how the energy changes when each node's weight is
increased (by adding 0.1 to it), moves the weights down that gradient by the learning rate `lr`, and rescales them
to the size of the identity. After `run()`, `E_opt` and `B_opt` hold the energy and the weights after each step.

```python
from nctpy.pipelines import ComputeOptimizedControlEnergy

control_task = {"x0": x0, "xf": xf, "S": np.eye(n_nodes), "rho": 1}  # no "B": the weights are what is optimised
pipeline = ComputeOptimizedControlEnergy(A=adjacency, control_task=control_task, system="continuous",
                                         c=1, T=1, n_steps=2, lr=0.01)
pipeline.run()
```

Each step solves one transition per node, so it is slow for large systems.

This optimises the control set, $B$: where input enters. It is a different method from optimising decay rates
(see {doc}`optimising_decay_rates`), which changes the system itself, $A$: how fast each node's activity decays.
