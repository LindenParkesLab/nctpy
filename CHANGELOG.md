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

### Fixed

- Boolean states given as `(N, 1)` columns now work like 1-D Boolean states in
  `get_control_inputs`, `sim_state_eq` and `minimum_energy_fast`. Previously they were not
  converted to floats, so `get_control_inputs` raised `TypeError` in discrete time, and with
  `xr='midpoint'` in continuous time.

### Added

- Type annotations for `nctpy.energies`, and a `py.typed` marker so type checkers use them. Other
  modules are annotated in later releases.

### Documentation

- The `nctpy.energies` docstrings use numpydoc format and now state each function's return shapes,
  including the discrete-time trajectory lengths (`T + 1` states, `T` inputs), what the two error
  terms measure, and `gramian`'s infinite-horizon and unstable cases.

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
