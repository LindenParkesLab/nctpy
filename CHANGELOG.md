# Changelog

## 1.0.3 (unreleased)

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
