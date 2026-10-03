# Measures of control energy

nctpy computes control energy in several ways. They answer different questions and come on different scales, so
compare energies only when they were computed the same way.

## Energy of a simulated transition

The protocol paper's approach: compute the optimal inputs for a transition with
{func}`~nctpy.energies.get_control_inputs`, integrate each node's squared input over time with
{func}`~nctpy.energies.integrate_u`, and sum over nodes.

```python
x, u, err = get_control_inputs(A_norm=A_norm, T=1, B=B, x0=x0, xf=xf, system="continuous")
node_energy = integrate_u(u)  # one value per node
energy = np.sum(node_energy)
```

This works for any control task: any states, control set, trajectory constraint and time system. It also gives
the trajectory and the errors (see {doc}`errors`). {class}`~nctpy.pipelines.ComputeControlEnergy` computes the same
energy for many tasks.

**Scale.** `integrate_u` takes the time points to be one unit apart. In continuous time they are 0.001 apart, so
these energies are 1000 times the integral $\int_0^T \lVert \mathbf{u}(t) \rVert^2 \,\mathrm{d}t$. The energies
printed in the protocol paper (for example 2604.71) are on this scale. Divide by 1000 to compare them with the
measures below, which are the integral itself.

In discrete time, `integrate_u` likewise applies Simpson's rule to the `T` inputs as if they were samples of a
continuous signal. That is not the sum of squared inputs, $\sum_t \lVert \mathbf{u}(t) \rVert^2$, often used as
the energy of a discrete-time system, and the two can differ substantially. If you need that sum, compute it
directly: `np.sum(u ** 2)`.

## Minimum energy, quickly

{func}`~nctpy.energies.minimum_energy_fast` computes the energy of minimum control (no constraint on the trajectory,
`S` all zeros) directly from the controllability Gramian, for many transitions at once. It is much faster than a
loop over `get_control_inputs`, but works only in continuous time and returns only node energies: no trajectory, no
errors. Its values are the integral, so they equal the energies above divided by 1000. The
{ref}`minimum_energy_fast` example compares the two.

## Energy over an infinite time horizon

Over an infinite time horizon, energy follows from the infinite-horizon controllability Gramian $W_c$ of the pair
$(A_\text{norm}, B)$, with no simulation (Kim et al., *Nat Commun* 2025). The system must be stable.

- {func}`~nctpy.energies.minimum_energy_infinite`: the minimum energy to reach a target state,
  $x_f^\top W_c^{-1} x_f$. Over an infinite horizon the initial state has decayed away, so only the target matters.
  It is the limit of the finite-horizon minimum energy from zero activity as `T` grows.
- {func}`~nctpy.energies.average_energy_infinite`: $\operatorname{tr}(W_c^{-1})$, the sum of the minimum energies to
  reach each node's unit state, and proportional to the average over all target states of unit size. It describes
  the system and its control set rather than one transition. Kim et al.'s Fig. 3C shows this quantity multiplied by
  0.001.

Both invert $W_c$, which becomes nearly singular when few nodes are controlled. `minimum_energy_infinite` returns
`xf_reached`, the target as reconstructed through the inverse Gramian: when it is far from `xf`, neither measure
can be trusted. The {doc}`/tutorials/partial_control` tutorial shows where that happens.

## The Gramian itself

{func}`~nctpy.energies.gramian` returns the controllability Gramian of $(A_\text{norm}, I)$ (every node controlled)
over a finite horizon `T`, or with `T=np.inf` over an infinite one, in either time system. Average controllability
is built on it (see {doc}`metrics`).

## Summary

| Measure | Function | Time systems | Scale |
| --- | --- | --- | --- |
| Energy of a simulated transition | `get_control_inputs` + `integrate_u`; `ComputeControlEnergy` | both | continuous: 1000 × the integral |
| Minimum energy | `minimum_energy_fast` | continuous | the integral |
| Minimum energy to a target, infinite horizon | `minimum_energy_infinite` | both | the integral |
| Average energy, infinite horizon | `average_energy_infinite` | both | $\operatorname{tr}(W_c^{-1})$ |
