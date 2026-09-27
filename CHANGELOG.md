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
