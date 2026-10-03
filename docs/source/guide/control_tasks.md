# Defining a control task

A control task is everything {func}`~nctpy.energies.get_control_inputs` needs besides the normalised connectome:
the initial and target states, the control set, how the state trajectory is constrained, and the time horizon.

```python
x, u, err = get_control_inputs(
    A_norm=A_norm, T=1, B=B, x0=x0, xf=xf, system="continuous", rho=1, S=S, xr="zero",
)
```

nctpy finds the inputs $\mathbf{u}(t)$ that take the system from `x0` at time 0 to `xf` at time `T` while
minimising the cost

$$
\int_0^T (\mathbf{x} - \mathbf{x}_r)^\top S\, (\mathbf{x} - \mathbf{x}_r) + \rho\, \mathbf{u}^\top \mathbf{u} \;\mathrm{d}t .
$$

## Brain states: `x0` and `xf`

A brain state is one value per node. The simplest are binary: the nodes of a state are active (1) and the rest are
not (0), for example the nodes of one functional system. {func}`~nctpy.utils.convert_states_str2int` turns a list
of each node's system into integer state labels, from which binary states follow:

```python
states, state_labels = convert_states_str2int(system_labels)
x0 = states == state_labels.index("Vis")
```

Non-binary states, such as patterns of activity from fMRI, give every node its own value. The protocol paper
recommends them when the biological plausibility of the transition matters most, and binary states for studying
transitions between specific parts of the network.

States of different magnitude need different amounts of energy: a larger target state (more active nodes, or
higher activity) costs more to reach. To compare transitions between states of different size, normalise each
state to unit length with {func}`~nctpy.utils.normalize_state`. Whether to do so depends on the comparison: if every
subject has the same states, for instance, differences in magnitude cannot confound differences between subjects.

## The control set: `B`

`B` is diagonal: $B_{ii}$ is how much input node $i$ receives. Three kinds are common.

- **Uniform full**, `B = np.eye(n_nodes)`: every node is a control node, with equal weight. This is the protocol
  paper's default, and the easiest for a transition to complete.
- **Partial**: only some nodes receive input. {func}`~nctpy.utils.mask_control_set` builds `B` from a mask of
  control nodes, and {func}`~nctpy.utils.random_control_set` draws them at random from a seed. Both take a
  `baseline`, a small weight for every other node.
- **Weighted**: control nodes have different weights, for example from a map of some biological property.
  {func}`~nctpy.utils.normalize_weights` rescales such a map to weights between 1 and 2 (by default, after
  ranking it), so that every node keeps some control:

  ```python
  B = np.diag(normalize_weights(brain_map))
  ```

The fewer the control nodes, the less likely a transition is to complete, so check the errors that
`get_control_inputs` returns. The {doc}`/tutorials/partial_control` tutorial shows transitions that complete and
transitions that do not.

## Constraining the trajectory: `S`, `rho` and `xr`

Without constraints, the cheapest inputs can drive activity to large values on the way from `x0` to `xf`. The cost
above can also penalise the trajectory's distance from a reference state, `xr`, at the nodes selected by `S`.

- **`S`** selects the nodes whose trajectory is constrained. `S = np.eye(n_nodes)` (the default) constrains every
  node, which the protocol paper calls *optimal control*. `S = np.zeros((n_nodes, n_nodes))` constrains none, and
  only the input is minimised: *minimum control*.
- **`rho`** weights the input against the trajectory. `rho = 1` weights them equally; a smaller `rho` makes the
  trajectory constraint count for more, a larger one for less. It must be positive. With `S` all zeros there is no
  trajectory term, so any positive `rho` gives the same result, which is what the protocol paper means by "rho is
  ignored". A very small `rho` makes the problem badly conditioned: transitions can fail to complete without any
  other sign of trouble, so check the errors when you lower it.
- **`xr`** is the state the trajectory is pulled towards: zero activity (`"zero"`, the default), the initial or
  target state (`"x0"`, `"xf"`), the midpoint between them (`"midpoint"`), or any state you pass. The protocol
  paper uses the default. Kim et al. (2025) set `xr` to the target state: with `xr = 0`, their optimisation of
  decay rates produced trajectories that simply switched off `x0` instead of completing the transition, so for
  that method they recommend `xr = xf`, which is {func}`~nctpy.optimize.optimize_decay_rates`' default.

## The time horizon: `T`

`T` is how long the system has to complete the transition (in continuous time, in arbitrary units; in discrete
time, a number of steps; see {doc}`time_systems`). The protocol paper uses `T = 1`.

In continuous time, a shorter horizon needs more energy, since the transition must happen faster. A longer one
needs less, up to a point, but too long a horizon makes the problem badly conditioned and the transition fails:
the errors become large. The paper's Supplementary Information shows the energy of its example transition at
`T` = 1, 2, 5 and 10. When you change `T`, check the errors. If a transition does not complete, the paper's
Table 1 suggests a longer horizon among other remedies.
