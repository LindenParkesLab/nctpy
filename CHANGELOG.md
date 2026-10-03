# Changelog

## Unreleased

### Changed

- `get_control_inputs` now raises `ValueError` when `rho` is not positive (`rho <= 0`, or NaN).
  Previously `rho = 0` returned NaN energies and NaN error terms, with numpy "divide by zero"
  warnings. That happened even with `S` all zeros, where any positive `rho` gives the same result.
  `rho` weights the cost of the control inputs, so a non-positive value does not define a control
  problem. `ComputeControlEnergy` and `ComputeOptimizedControlEnergy` now raise the same error
  for a task with `rho = 0`. Ill-conditioned or incomplete transitions are unaffected: they still
  return their energies and error terms.
- `nctpy.energies` converts its array inputs to float64 before computing. Results for float64
  input are unchanged, bit for bit. `gramian` on a lower-precision (e.g. float32) continuous-time
  matrix previously computed in that precision; it now computes in float64, so its values change
  slightly.
- `gramian` now raises when `system` is missing or is not `'continuous'` or `'discrete'`, with the
  same error and message as `get_control_inputs`. Previously it silently returned None.
  `ave_control` is unaffected. `gramian` also no longer stores every step's matrix exponential,
  which needed `N × N × (T/0.001 + 1)` floats (about 1.3 GB for 400 nodes at T = 1). Its results
  are unchanged up to floating-point rounding (within 1e-14 relative).
- `nctpy.metrics` also computes in float64. `ave_control` is unchanged, bit for bit, for float64
  input. `modal_control` now sums in a different order, which changes its values by rounding only
  (within 3e-15 relative). Discrete `ave_control` and `modal_control` on lower-precision (e.g.
  float32) matrices previously computed and returned float32. They now return float64 (for
  `modal_control`, values move by up to ~3e-4 relative, float32's own precision).
- `matrix_normalization` also computes in float64. Results for float64, integer and Boolean input
  are unchanged, bit for bit. A float32 matrix previously had its spectral radius computed in
  single precision (and, for discrete time, a float32 result). It now gets a float64 result,
  which moves values by up to ~5e-8 relative.
- `ComputeOptimizedControlEnergy` now uses a task's optional `'xr'` (reference state), as
  `ComputeControlEnergy` already did. Previously it ignored the key and always used `xr='zero'`.
  Tasks without `'xr'` (including the one in the protocol paper's Supplementary Information) give
  unchanged results; tasks that set it now get energies and optimised weights for that reference
  state.
- `get_control_inputs` given several states at once (e.g. `x0` of shape `(N, k)`) now raises a
  `ValueError` saying it takes a single state. It already raised `ValueError`, from inside numpy,
  with a message about broadcasting.
- `nctpy.plotting.surface_plot`'s `fsaverage` argument now defaults to `None`, meaning fsaverage5 is
  loaded when the plot is drawn. Previously the default was evaluated once, when `nctpy.plotting`
  was imported, so importing the module loaded the surface even if `surface_plot` was never called.
  Every call that passes `fsaverage` (as all code in the protocol paper does), by keyword or by
  position, behaves exactly as before.
- `get_control_inputs` raises `ValueError` naming the options when `xr` is a string other than
  `'zero'`, `'x0'`, `'xf'` or `'midpoint'` (so does `ComputeControlEnergy` for such a task).
  Previously it failed inside numpy with a `TypeError` about ufunc loops.
- `get_null_p` raises `ValueError` for an unknown `version`. Previously it failed with
  `UnboundLocalError`.

### Changed — directed connectomes only

- **`ave_control` on directed connectomes now follows Gu et al. (2015) in both systems.** The
  average controllability of node i is the trace of the controllability Gramian when input enters
  at node i alone: how much input at node i spreads into the network. Previously, for a directed
  `A`, continuous-time `ave_control` returned the diagonal of the Gramian with input at every node
  (how strongly each node is driven: a different quantity, about 2% off on a test matrix).
  Discrete-time `ave_control` used a formula based on the real Schur decomposition that is exact
  only for symmetric `A`. It depended on the order of the nodes, and was off by up to 64% at single
  nodes on a test matrix. **Symmetric (undirected) connectomes are unaffected**: continuous-time
  values are bit-identical, and discrete-time values agree to rounding (within 2e-13 relative).
  The protocol paper and its Supplementary Information use `ave_control` only on undirected
  connectomes.
- `modal_control`'s docstring now states that modal controllability is defined for undirected
  connectomes; for a directed `A` its values are an approximation that depends on node order. Its
  results are unchanged.

### Performance

- Computing many transitions on one system is much faster. The parts of `get_control_inputs` and
  `minimum_energy_fast` that depend only on the system, not on the states, are now computed once
  and reused while consecutive calls share that system. For `get_control_inputs` those are the
  matrix exponentials and the discrete-time factorisation; for `minimum_energy_fast` they are the
  Gramian and its pseudo-inverse. Results are identical, bit for bit: reuse happens only when
  `A_norm`, `T`, `B`, `S` and `rho` are identical byte for byte. All 49 transitions between 7
  states, the protocol paper's workload (`benchmarks/transitions.py`, 8 BLAS threads):

  | workload | 200 nodes | 400 nodes |
  |---|---|---|
  | `get_control_inputs`, continuous, T = 1 | 34.5 → 18.7 ms | 116 → 37 ms |
  | `get_control_inputs`, discrete, T = 3 | 61.9 → 2.1 ms | 423 → 11.5 ms |
  | `ComputeControlEnergy`, continuous, T = 1 | 33.8 → 4.7 ms | 130 → 12.7 ms |
  | `minimum_energy_fast`, T = 1 (100 nodes) | 90.6 → 1.9 ms | |

  (per transition). Calls that never share a system, such as a loop that perturbs B, cost about the
  same as before. The most recent system's matrices stay in memory until a call uses a different
  system. In continuous time that is about 15 N² floats (about 19 MB at N = 400).
- `ComputeControlEnergy` also solves consecutive tasks that share a system (`B`, `S`, `rho`) together,
  in memory-bounded batches; that is most of its gain in the table above. A transition that
  completes (both error terms below 1e-8) agrees with `get_control_inputs` to rounding (within
  1e-13 relative in tests; batched and single products round differently). A transition that does
  not complete is solved again on its own, so its energy is exactly what it was.
- `minimum_energy_fast` documents that it accepts k transitions at once, as `(N, k)` columns of
  initial and target states.

### Fixed

- `nctpy.plotting.roi_to_vtx` (and so `surface_plot`) leaves unlabelled vertices (label -1 in a
  FreeSurfer annotation) at 0, like label 0 (e.g. the medial wall). Previously they took the
  second-to-last parcel's value. The Schaefer annotations used in the protocol paper have no
  unlabelled vertices, so its figures are unchanged.
- Boolean states given as `(N, 1)` columns now work like 1-D Boolean states in
  `get_control_inputs`, `sim_state_eq` and `minimum_energy_fast`. Previously they were not
  converted to floats, so `get_control_inputs` raised `TypeError` in discrete time, and with
  `xr='midpoint'` in continuous time.
- `get_fdr_p` accepts p-values of any shape and returns them corrected in that shape. Previously
  input with more than two dimensions failed with an `AssertionError`. 1-D and 2-D results are
  unchanged.
- Importing `nctpy.utils` or `nctpy.plotting` on Python 3.12 or later no longer emits
  `SyntaxWarning: invalid escape sequence '\m'` (from the p-value and correlation labels, which are
  unchanged).

### Added

- **`nctpy.optimize.optimize_decay_rates`**: fits each node's decay rate (intrinsic neural timescale) for a state
  transition, as in Kim et al., Nat Commun 16:11639 (2025). It takes the raw connectome, normalises it for continuous
  time, and fits the decay rates by gradient descent (PyTorch). The aim is that the optimal trajectory passes through
  the reference state at the midpoint, while every eigenvalue stays negative.
  - The defaults are the paper's. By default it removes self-connections first (`zero_diagonal=True`), because the
    method assumes a connectome without them.
  - It returns a `DecayRateFit` with the fitted decay rates (in `matrix_normalization`'s `decay` convention), the
    fitted matrix, and the training traces.
  - On the paper's mouse connectome it reproduces the published fits to within 1e-15, and stops at the same step.
  - It needs PyTorch: `pip install "nctpy[optimize]"`. Without it, importing `nctpy.optimize` raises an
    `ImportError` that says so.
- `nctpy.utils.random_control_set(n_nodes, n_control_nodes, seed=0, baseline=0.0)` and
  `mask_control_set(mask, baseline=0.0)` build partial control sets (B).
  - `random_control_set` draws its control nodes at random. The same seed gives the same nodes as in
    Kim et al. (2025), Fig. 3, but uses a local random generator, leaving numpy's global state
    untouched.
  - `mask_control_set` controls the nodes in a mask, e.g. one system's nodes,
    `states == state_labels.index('Vis')`.
  - `baseline` gives every other node a small control weight.
- `nctpy.energies.minimum_energy_infinite(A_norm, B, xf, system)` and
  `average_energy_infinite(A_norm, B, system)`: control energy from the infinite-horizon
  controllability Gramian, which solves a Lyapunov equation instead of integrating numerically or
  taking a matrix exponential (Kim et al., 2025).
  - `minimum_energy_infinite` returns the node-wise minimum energy to reach `xf`
    (`xf^T Wc^-1 xf` in total), and the target as reconstructed through the inverse Gramian, which
    shows how accurately the Gramian could be inverted.
  - `average_energy_infinite` returns `trace(Wc^-1)`, the average minimum energy over target
    states, which does not depend on a transition. This is the energy shown in that paper's Fig. 3C
    (multiplied there by dt = 0.001).

  Both work in continuous and discrete time, and return NaN for an unstable system, for which the
  infinite-horizon Gramian does not exist.
- `matrix_normalization(..., decay=None)`: a keyword-only decay rate per node for continuous-time
  systems, implementing Kim et al. (2025), Eq. 4: `A / (c + l) - diag(decay)`. A larger decay
  means stronger self-inhibition. A scalar applies to every node, and the default (`None`, or
  `decay=1`) gives exactly today's `A / (c + l) - I`. Self-connections in `A` are kept, so
  `zero_diagonal=True` removes them first if wanted. Passing `decay` for a discrete-time system,
  or a vector without one value per node, raises `ValueError`. A value `v` reported in that
  paper's Fig. 2B corresponds to `decay = 1 - v`.
- `nctpy.null_models`: the null models are also importable from within nctpy, e.g.
  `from nctpy.null_models.geomsurr import geomsurr`. It re-exports the same functions as the
  top-level `null_models` package, which the protocol paper imports from and which keeps working.
- `matrix_normalization(..., zero_diagonal=False)`: a keyword-only option that sets the diagonal of
  a copy of `A` to zero before normalising. The model assumes no self-connections: each node's own
  dynamics come from the normalisation, not from `diag(A)`. nctpy still never removes
  self-connections unless asked, and the default stays `False`. The protocol paper's code uses the
  PNC connectome with its diagonal intact (Vis → Default energy 2604.71; with
  `zero_diagonal=True`, 2638.03). The array passed in is never modified.
- Type annotations for `nctpy.energies`, `nctpy.metrics`, `nctpy.utils` and `nctpy.pipelines`, and
  a `py.typed` marker so type checkers use them. Other modules are annotated in later releases.

### Documentation

- The `nctpy.metrics` docstrings use numpydoc format and state each metric's formula.
- The `nctpy.pipelines` docstrings use numpydoc format. They list the attributes `run()` sets.
  They also say that `ComputeOptimizedControlEnergy` ignores a task's `'B'`, because the control
  weights are what it optimises, and describe its gradient step.
- The `nctpy.utils` docstrings use numpydoc format. `matrix_normalization` documents its formula
  and `l`, including that an `l` below the matrix's own spectral radius forfeits the stability
  guarantee (nothing checks this, and an unstable system still returns values).
- The `nctpy.energies` docstrings use numpydoc format and now state each function's return shapes,
  including the discrete-time trajectory lengths (`T + 1` states, `T` inputs), what the two error
  terms measure, and `gramian`'s infinite-horizon and unstable cases.
- `get_control_inputs`' docstring no longer says that `rho=1` "equals maximum constraint". `rho`
  weights the cost of the control signals against that of the state trajectory: a smaller `rho`
  constrains the trajectory more, a larger one less.
- The `nctpy.plotting` functions have docstrings, and `get_p_val_string`, `roi_to_vtx` and the
  `null_models.geomsurr` functions use numpydoc format.
- The documentation has an API reference covering every public function and class, generated from
  the docstrings. A test checks that it lists every public symbol.
- The Getting started page is an executed notebook: its outputs and figures come from running the code
  shown, on the same synthetic example as before. CI executes it (and every tutorial added later) with
  `docs/execute_notebooks.py`, which also refreshes the committed outputs (`--write`). The example
  matrix is directed, so its average controllability values change with `ave_control`'s move to
  Gu et al.'s definition (above).
- A Tutorials section, starting with the protocol workflow on a synthetic, spatially embedded
  connectome: control energy for one transition (Procedures 1 and 2, Box 1), energies for every
  pair of states with `ComputeControlEnergy`, and spatial null networks with `geomsurr`.
- A tutorial on fitting nodes' decay rates with `optimize_decay_rates` (Kim et al. 2025): the fit and its
  training traces, the fitted rates, the fitted matrix in use, rebuilding it with
  `matrix_normalization(..., decay=...)`, and fitting several transitions at once.
- A tutorial on partial control sets and infinite-horizon energy: control sets from a mask or at random
  (`mask_control_set`, `random_control_set`, `baseline`), incomplete transitions and their error terms,
  `average_energy_infinite` and `minimum_energy_infinite` across control-set sizes, and checking
  `xf_reached` before trusting either with sparse control sets.
- The examples comparing controllability metrics and comparing `minimum_energy_fast` with
  `get_control_inputs` now run on synthetic data and are executed in CI, like the tutorials. They
  used to load real data from hard-coded paths and called `reg_plot` with arguments it no longer
  takes. The `minimum_energy_fast` example no longer claims a ~300-fold speed-up: `get_control_inputs`
  now reuses its work across a loop of transitions, and the page shows the timings measured when it
  was run. The other examples, which analyse real data that cannot be shared, remain as static pages.
- A "Reproducing the papers" page. It gives the protocol paper's main-text code in a form that runs
  when copied, lists the seven typesetting defects that stop the printed code from running (and
  one misprinted output), gives the printed outputs to expect, and notes the differences older
  versions of nctpy show. It also points users of Kim et al. (2025) to the paper's repository and
  to `nctpy.optimize`. A test checks that the page's code is the paper's, with only those repairs.
- A User guide section, starting with pages on choosing a time system (and which functions support
  each), normalising the connectome (`c`, a shared `l`, self-connections and `zero_diagonal`,
  per-node `decay`, checking stability) and defining a control task (states, control sets, `S`,
  `rho`, `xr` and `T`).
- User-guide pages on error terms and failed transitions (what the two errors measure, the
  protocol paper's troubleshooting table, and further causes to check), on the measures of control
  energy (their scales, including that `integrate_u`'s continuous-time energies are 1000 times the
  time integral, and the infinite-horizon measures) and on the controllability metrics (Gu et al.'s
  definition of average controllability for directed connectomes, and the scope of modal
  controllability).
- User-guide pages on the pipeline classes (how `ComputeControlEnergy` batches transitions, and how
  `ComputeOptimizedControlEnergy`'s optimisation of control weights differs from optimising decay
  rates), on optimising decay rates (the objective, defaults and options), on null models
  (`geomsurr`'s surrogates and scope, p-values, and surrogate brain maps from other packages) and on
  performance (reuse within a system, choosing a function, threads).

### Repository

- The notebooks in `scripts/` run from any clone. They find the repository's root from their own
  location instead of hard-coded paths. They use the installed nctpy (`pip install -e ".[paper]"`)
  instead of adding `src/` to `sys.path`, and create `results/` if needed. A cell that plots saved
  null models skips, with a message, when the nulls have not been computed. `scripts/README.md`
  explains how to run them. On a machine that has the saved nulls, every cell runs as before.

### Packaging

- The documentation is built with pydata-sphinx-theme, MyST and myst-nb, matching the lab's other packages. The
  `docs` extra installs that toolchain; `docs/requirements.txt` is gone, and Read the Docs and CI install
  `nctpy[docs]` instead.
- `packaging` is no longer a dependency. nctpy used it only to choose between `simps` and
  `simpson` on SciPy < 1.6, which cannot be installed on Python 3.10 or later. The core
  dependencies are now numpy, scipy, tqdm and statsmodels.

## 1.1.0 (2026-09-28)

**Upgrading from 1.0.x:** if you use `nctpy.plotting` or run the code from the Nature Protocols paper,
install with `pip install "nctpy[paper]"` (or `"nctpy[plot]"`). nctpy now requires Python 3.10 or later.

### Changed — action needed if you use `nctpy.plotting`

- **The plotting dependencies are now optional.** `pip install nctpy` installs only numpy,
  scipy, tqdm, packaging and statsmodels. matplotlib, seaborn, nibabel and nilearn are no longer
  installed automatically. `nctpy.plotting` needs them, and without them it raises an
  `ImportError` that says how to install them:

  ```bash
  pip install "nctpy[plot]"     # matplotlib, seaborn, nibabel, nilearn
  pip install "nctpy[paper]"    # the above plus pandas and scikit-learn
  ```

  **Following the Nature Protocols paper?** Install with `pip install "nctpy[paper]"`. The
  paper's import block imports these packages directly, so a plain `pip install nctpy` is no
  longer enough to run it.

  This is the only change in 1.x that can break an existing installation. Everything else in
  1.x is additive.

- **Python 3.10 or later is now required** (previously 3.9). On Python 3.9, pip will keep
  installing nctpy 1.0.x.

### Fixed

- Installing nctpy no longer installs a top-level `tests` package alongside it. Only `nctpy` and
  `null_models` are installed.
- `nctpy.plotting.surface_plot` works with nilearn 0.13 and later. It drew continuous data with
  nilearn's `plot_surf_roi`, which since nilearn 0.13 rejects negative or non-integer values,
  so it raised `ValueError` for almost any real input (e.g. average controllability, or the fMRI
  states in the protocol paper's Supplementary Information). It now uses `plot_surf` with the
  same settings. With nilearn 0.12 and earlier, the figures are pixel-for-pixel identical to
  before. With nilearn 0.14, which removed the `darkness` option, the background shading may
  differ very slightly.

### Added

- `nctpy.__version__`.
- Extras: `plot`, `paper`, `docs` and `dev` (the test suite's dependencies).

### Packaging

- Packaging metadata now lives in `pyproject.toml` (PEP 621); `setup.cfg` is gone. The version
  has a single source, `nctpy.__version__`.

## 1.0.3 (2026-09-27)

### Added

- `matrix_normalization(..., l=None)`: an optional fixed spectral radius. When given, `A` is
  normalised by `c + l` instead of `c` plus its own largest absolute eigenvalue, so several
  connectomes (e.g. subjects) can share one normalisation, typically with `l` set to the maximum
  spectral radius across them. The default `l=None` gives exactly the previous behaviour. An `l`
  smaller than a matrix's own spectral radius no longer guarantees a stable system.

### Fixed

- Continuous-time `ave_control` (and `gramian` with `system='continuous'`) no longer raises
  `AttributeError: module 'scipy.integrate' has no attribute 'simps'` on SciPy ≥ 1.14. The SciPy
  version check compared version strings, so every SciPy from 1.10 onwards was routed to `simps`,
  which SciPy removed in 1.14. The check now uses `packaging.version.Version`, and `simpson` is
  called with keyword arguments as current SciPy requires.

  Numerical note: `ave_control` always uses `T=1`, and its results are unchanged (bit-identical
  on SciPy 1.13; within 5e-15 relative on SciPy 1.16 / numpy 2.3, i.e. rounding). SciPy 1.11
  changed how Simpson's rule treats an *even* number of samples. With its fixed step of 0.001,
  `gramian` integrates over `T/0.001 + 1` samples. That count is odd whenever `T` is a multiple
  of 0.002 (e.g. 0.1, 0.5, 1, 2, 10), so results for those horizons are unaffected. Calling
  `gramian` directly with other horizons may give slightly different values than under SciPy
  < 1.11.

- `null_models.geomsurr` no longer modifies the matrix passed as `W`. Previously it set the
  caller's diagonal to zero in place, so after running the null-model code, a connectome with
  self-connections had silently lost them for the rest of the session. It also no longer resets
  numpy's global random seed; it uses a local generator instead. **The surrogates are unchanged**:
  identical for every seed.

### Documentation

- `geomsurr`'s docstring now states its scope. It rewires edges between pairs of nodes, so its
  surrogates always have a zero diagonal and it is suited to connectomes without
  self-connections. For a connectome with self-connections, compute the observed statistic on the
  zero-diagonal version of the connectome.
