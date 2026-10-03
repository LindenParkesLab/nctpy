# Error terms and failed transitions

The model does not guarantee that a transition completes. With too few control nodes, too short or too long a
time horizon, or a badly conditioned system, the inputs nctpy finds may not take the system to the target state.
{func}`~nctpy.energies.get_control_inputs` therefore returns two numerical errors with every result, and checking
them is part of every analysis.

```python
x, u, err = get_control_inputs(A_norm=A_norm, T=1, B=B, x0=x0, xf=xf, system="continuous")
inversion_error, reconstruction_error = err
completed = inversion_error < 1e-8 and reconstruction_error < 1e-8
```

## What the errors measure

- The **inversion error** measures how well nctpy could solve for the optimal inputs: the residual of the linear
  system it solves (in continuous time, for the costate; in discrete time, for the whole trajectory). A large value
  means the problem is badly conditioned.
- The **reconstruction error** measures whether the solution does what it should. In continuous time, it is the
  distance between the state at time `T` and the target state; in discrete time, how far the trajectory departs
  from the system's dynamics.

The protocol paper treats both as adequately small below $10^{-8}$. Their exact values below that threshold are
rounding noise: they vary between machines and versions, and carry no meaning.

## Failed transitions return

When a transition fails, nctpy still returns its trajectory, inputs and energy, together with large errors. It
never raises an error for this, so that you can inspect the result. The energy of a failed transition is
meaningless, often astronomically large (the protocol paper's Supplementary Information shows one of $10^{24}$),
so filter on the errors before using any energy.

{class}`~nctpy.pipelines.ComputeControlEnergy` returns only energies, not errors. When you use it, check
representative transitions with `get_control_inputs`, especially with partial control sets or an unusual `T` or
`rho`.

## When a transition fails

The protocol paper's Table 1 suggests:

| Problem | Possible reason | What to try |
| --- | --- | --- |
| Inversion error above $10^{-8}$ | Badly conditioned problem | With a partial control set, try a full one. With a full control set, try a smaller system (fewer nodes). Check whether the connectome has disconnected components (groups of nodes that cannot reach each other). |
| Reconstruction error above $10^{-8}$ | The transition did not complete: $x(T) \neq x_f$ | Increase the time horizon, `T`. With a partial control set, try a full one. With a full control set, try a smaller system. Check for disconnected components. |

Three further causes are worth checking first:

- **An unstable system**, after changing `c`, `l` or `decay` (see {doc}`normalisation`).
- **A very small `rho`**, which makes the problem badly conditioned without any other warning (see
  {doc}`control_tasks`).
- **Too long a time horizon.** A longer `T` lowers the energy a transition needs and can help a transition
  complete, but beyond some point the problem becomes badly conditioned and transitions fail. After changing `T`,
  check the errors again.

The {doc}`/tutorials/partial_control` tutorial shows transitions that complete and transitions that do not.
