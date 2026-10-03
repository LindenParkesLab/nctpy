(reproducing)=
# Reproducing the papers

## The protocol paper

Parkes, Kim et al., *Nature Protocols* 19:3721–3749 (2024),
[doi:10.1038/s41596-024-01023-w](https://doi.org/10.1038/s41596-024-01023-w).

The paper prints the code of its two procedures (control energy for a state transition, and average
controllability) and of Box 1, and its Supplementary Information prints variations on them. nctpy keeps every one
of those code blocks running, unmodified, and doing what the paper describes. A test runs them on the paper's data
and compares their outputs with the printed values; because the data cannot be shared, it runs on the lab's
machines rather than in public CI.

This page gives the main text's code in a form you can copy and run. Copied straight from the PDF, it does not run:
typesetting broke some of its lines (listed below). The code here has exactly those repairs and no other changes.
The Supplementary Information's code was typeset from working listings and runs as printed.

### Before you start

- Install nctpy with the extra packages the paper's code imports: `pip install "nctpy[paper]"`.
- The paper analyses structural connectomes from the Philadelphia Neurodevelopmental Cohort (PNC), available
  through the Database of Genotypes and Phenotypes under accession
  [phs000607.v3.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000607.v3.p2). nctpy
  does not distribute these data. Point `datadir` and `adjacency_file` below at your own copy. The file that maps
  each of the 200 parcels to a functional system, `pnc_schaefer200_system_labels.txt`, is in the repository's
  [data folder](https://github.com/LindenParkesLab/nctpy/tree/main/data).
- The notebooks in the repository's [scripts folder](https://github.com/LindenParkesLab/nctpy/tree/main/scripts)
  reproduce the paper's analyses, including those of the Supplementary Information, on the same data.

### Typesetting defects in the printed code

| Page | Block | Defect | Copied as printed |
| --- | --- | --- | --- |
| 3737 | Imports | `ComputeControlEnergy,` / `ComputeOptimizedControlEnergy` broken across two lines | `SyntaxError` |
| 3737 | Imports | `convert_states_str2int,` / `normalize_state, ...` broken across two lines | `SyntaxError` |
| 3737 | Imports | `add_` / `module_lines` broken across two lines | `SyntaxError` |
| 3740 | Step 4A | The comment `# nodes in state trajectory to be constrained` wrapped onto a second line | `SyntaxError` |
| 3740 | Step 4A | `get_control_` / `inputs(` broken across two lines | `SyntaxError` |
| 3741 | Step 4(iv) | The comment `# the second numerical error corresponds to the reconstruction error` wrapped, leaving `error` on a line of its own | `NameError` |
| 3743 | Box 1 | The body of the `for cax in ax.reshape(-1):` loop lost its indentation | `IndentationError` |

One printed *output* is also off: p. 3739 gives `states.shape` as `(200,)` but lists 199 values. The
right-hemisphere Limbic run shows five 3s where the labels file has six; the code itself is correct.

### The code

Expected outputs are the values printed in the paper. Current versions of nctpy reproduce them to the printed
precision. The two numerical errors differ in their last digits from the printed ones: they are rounding noise,
and only their size matters (both below $10^{-8}$). With numpy 2 or later, labels print as `np.str_('Vis')`
rather than `'Vis'`; the values are the same.

**Imports** (p. 3737)

```python
# import
import os
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
from scipy.spatial import distance
from sklearn.cluster import KMeans
from tqdm import tqdm
# import plotting libraries
import matplotlib.pyplot as plt
import seaborn as sns
from nilearn import datasets
from nilearn import plotting
# import nctpy functions
from nctpy.energies import integrate_u, get_control_inputs
from nctpy.pipelines import ComputeControlEnergy, ComputeOptimizedControlEnergy
from nctpy.metrics import ave_control
from nctpy.utils import matrix_normalization, convert_states_str2int, normalize_state, normalize_weights, get_null_p, get_fdr_p
from nctpy.plotting import roi_to_vtx, null_plot, surface_plot, add_module_lines
from null_models.geomsurr import geomsurr
```

**Load the connectome** (pp. 3737–3738)

```python
# directory where data is stored
datadir = '/path/to/data'
adjacency_file = 'structural_connectome.npy'
# load adjacency matrix
adjacency = np.load(os.path.join(datadir, adjacency_file))
n_nodes = adjacency.shape[0]
print(adjacency.shape)
# check for self-connections
print(np.any(np.diag(adjacency) > 0))
# get density including self connections
density = np.count_nonzero(np.triu(adjacency, k=0)) / (n_nodes**2 / 2)
print(density)
```

```text
(200, 200)
True
0.9768
```

**Step 1: define a time system** (p. 3738). Delete the line you do not need; the printed results are for continuous time.

```python
# determine time system. Note, delete the line below that is not needed.
system = "discrete"
# or
system = "continuous"
```

**Step 2: normalise the adjacency matrix** (p. 3738).

```python
# normalize adjacency matrix
adjacency_norm = matrix_normalization(A=adjacency, system=system, c=1)
```

**Step 3: define a control task** (pp. 3739–3740): brain states from the system labels, and the control set.

```python
# load node-to-system mapping
system_labels = list(
 np.loadtxt(os.path.join(datadir, "pnc_schaefer200_system_labels.txt"),
dtype=str)
)
print(len(system_labels))
print(system_labels[:20])
```

```text
200
['Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'Vis', 'SomMot', 'SomMot', 'SomMot', 'SomMot', 'SomMot', 'SomMot']
```

```python
# use list of system names to create states
states, state_labels = convert_states_str2int(system_labels)
print(type(state_labels), len(state_labels), state_labels)
print(type(states), states.shape, states)
```

```text
<class 'list'> 7 ['Cont', 'Default', 'DorsAttn', 'Limbic', 'SalVentAttn', 'SomMot', 'Vis']
<class 'numpy.ndarray'> (200,) [6 6 6 6 6 6 6 6 6 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 2 2 2 2 2 2 2
 2 2 2 2 2 2 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1
 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 6 6 6 6 6 6 6 6 6 6
 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 2 2 2 2 2 2 2 2 2 2 2 2 2 4
 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1
 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
```

```python
# extract initial state
initial_state = states == state_labels.index('Vis')
# extract target state
target_state = states == state_labels.index('Default')
```

```python
# normalize state magnitude
initial_state = normalize_state(initial_state)
target_state = normalize_state(target_state)
```

```python
# specify a uniform full control set: all nodes are control nodes
# and all control nodes are assigned equal control weight
control_set = np.eye(n_nodes)
```

**Step 4: compute control signals and the state trajectory** (pp. 3740–3741), then check the numerical errors.

```python
# set parameters
time_horizon = 1 # time horizon (T)
rho = 1 # mixing parameter for state trajectory constraint
trajectory_constraints = np.eye(n_nodes) # nodes in state trajectory to be constrained
# get the state trajectory, x(t), and the control signals, u(t)
state_trajectory, control_signals, numerical_error = get_control_inputs(
 A_norm=adjacency_norm,
T=time_horizon,
B=control_set,
x0=initial_state,
xf=target_state,
system=system,
 rho=rho,
 S=trajectory_constraints,
)
```

```python
# print errors
thr = 1e-8
# the first numerical error corresponds to the inversion error
print(
 "inversion error = {:.2E} (<{:.2E}={:})".
format(numerical_error[0], thr, numerical_error[0] < thr
 )
)
# the second numerical error corresponds to the reconstruction error
print(
 "reconstruction error = {:.2E} (<{:.2E}={:})".
format(numerical_error[1], thr, numerical_error[1] < thr
 )
)
```

```text
inversion error = 1.36E-15 (<1.00E-08=True)
reconstruction error = 5.16E-14 (<1.00E-08=True)
```

**Step 5: visualise the state trajectory and control signals** (Box 1, p. 3743).

```python
f, ax = plt.subplots(3, 2, figsize=(7, 7))
# plot control signals for initial state
ax[0, 0].plot(control_signals[:, initial_state != 0], linewidth=0.75)
ax[0, 0].set_title("A | control signals, x0")
# plot state trajectory for initial state
ax[0, 1].plot(state_trajectory[:, initial_state != 0], linewidth=0.75)
ax[0, 1].set_title("B | neural activity, x0")
# plot control signals for target state
ax[1, 0].plot(control_signals[:, target_state != 0], linewidth=0.75)
ax[1, 0].set_title("C | control signals, xf")
# plot state trajectory for target state
ax[1, 1].plot(state_trajectory[:, target_state != 0], linewidth=0.75)
ax[1, 1].set_title("D | neural activity, xf")
# plot control signals for bystanders
ax[2, 0].plot(
 control_signals[:, np.logical_and(initial_state == 0, target_state == 0)],
 linewidth=0.75,
)
ax[2, 0].set_title("E | control signals, bystanders")
# plot state trajectory for bystanders
ax[2, 1].plot(
 state_trajectory[:, np.logical_and(initial_state == 0, target_state == 0)],
 linewidth=0.75,
)
ax[2, 1].set_title("F | neural activity, bystanders")
for cax in ax.reshape(-1):
    cax.set_ylabel("activity")
    cax.set_xlabel("time (a.u.)")
    cax.set_xticks([0, state_trajectory.shape[0]])
    cax.set_xticklabels([0, time_horizon])
f.tight_layout()
plt.show()
```

**Step 6: compute control energy** (pp. 3741–3742).

```python
# integrate control signals to get control energy
node_energy = integrate_u(control_signals)
print(node_energy.shape)
print(np.round(node_energy[:5], 2))
```

```text
(200,)
[21.13 37.65 23.55 21.55 28.34]
```

```python
# summarize nodal energy
energy = np.sum(node_energy)
print(np.round(energy, 2))
```

```text
2604.71
```

**Procedure 2: average controllability** (p. 3744).

```python
# compute average controllability
average_controllability = ave_control(A_norm=adjacency_norm,
system=system)
```

### Notes for users of older versions, and on self-connections

- **Procedure 2 needs nctpy 1.0.3 or later on SciPy 1.14 or later.** Earlier versions of `ave_control` raise an
  `AttributeError` for continuous-time systems there.
- **Before nctpy 1.0.3, `geomsurr` modified the matrix passed to it**, setting its diagonal to zero (and it reset
  numpy's global random seed). Running the Supplementary Information's null-model code therefore changed
  `adjacency` itself, and energies computed from it afterwards changed too (Vis → Default: 2638.03 instead of
  2604.71). From 1.0.3 it works on a copy.
- **Since nctpy 1.1.0, the plotting libraries are optional.** Install with `pip install "nctpy[paper]"` as above.
- **The PNC connectome has self-connections**, as the paper's check prints (`True`), and the printed results keep
  them. nctpy never removes them on your behalf. `matrix_normalization(..., zero_diagonal=True)` removes them if
  your model should not have them; on these data that gives a Vis → Default energy of 2638.03 instead of 2604.71.

## Kim et al. (2025)

Kim, …, Parkes, *Nature Communications* 16:11639 (2025),
[doi:10.1038/s41467-025-66542-w](https://doi.org/10.1038/s41467-025-66542-w).

- The paper's analysis code is in the [nct_xr repository](https://github.com/LindenParkesLab/nct_xr), kept as it
  was for the paper.
- Its method for fitting nodes' decay rates is part of nctpy, as {func}`~nctpy.optimize.optimize_decay_rates`
  (see the {doc}`/tutorials/decay_rates` tutorial). With the paper's settings, which are its defaults, it
  reproduces the paper's fits: checked against nct_xr's saved results, it stops after the same number of training
  steps, its decay rates agree to within $10^{-10}$, and its energies to within a relative $10^{-6}$.
- A decay value $v$ in the paper's Fig. 2B corresponds to `decay = 1 - v` in nctpy.
- The control energy in Fig. 3C is the average energy over an infinite time horizon,
  {func}`~nctpy.energies.average_energy_infinite`, multiplied by the time step, 0.001 (see the
  {doc}`/tutorials/partial_control` tutorial).
