# Contributing to nctpy

Thank you for your interest in nctpy. This guide covers setting up a development environment,
running the tests, and the rules every change must follow.

nctpy is the reference implementation for the protocol paper by Parkes, Kim et al., *Nature
Protocols* 19:3721–3749 (2024), doi:10.1038/s41596-024-01023-w. Researchers follow that paper's
printed code, so **the code printed in the paper and its Supplementary Information must keep
working, unchanged**. Most of the rules below follow from that.

## Setting up

```bash
git clone https://github.com/LindenParkesLab/nctpy.git
cd nctpy
pip install -e ".[dev]"   # nctpy, the plotting and paper extras, and the test tools
pre-commit install        # runs ruff and the data guards before each commit
```

nctpy needs Python 3.10 or later.

## Running the tests

Run the tests from `src/tests`, because some tests read fixture files relative to that directory:

```bash
cd src/tests
python -m pytest
```

The suite includes:

- `test_api_contract.py`: the public API (import paths, signatures, return types) must not change.
- `test_regression.py`: results must match stored, synthetic fixtures within tolerance.
- `test_properties.py`: mathematical properties that must hold.
- `test_plotting.py`, `test_nctpy.py`: further unit and smoke tests.
- `test_paper_code.py` and `test_notebooks.py`: run the protocol paper's printed code and the
  notebooks in `scripts/`. They need the paper's data, which is not distributed with the
  repository, so **they skip unless you have a local copy in `data/`**. The maintainers run them
  before each release.

## What CI checks

Every push and pull request runs:

- **lint:** `ruff check`, `ruff format --check` and `mypy`, configured in `pyproject.toml`
  (pre-commit runs the ruff checks locally);
- **test:** the built wheel, on Ubuntu and macOS, Python 3.10–3.13, with numpy 1.x and 2.x;
- **coverage:** no module's test coverage may fall below the baseline in `src/tests/baseline.json`;
- **docs:** the Sphinx documentation must build without warnings.

Some existing files are exempt from parts of the ruff checks (see `pyproject.toml`). These
exemptions are being removed module by module. New files must pass in full.

## Rules for changes

- **Never change a public signature.** Within 1.x, new behaviour arrives as new keyword arguments
  whose defaults reproduce the existing behaviour exactly. Renaming or removing a parameter, or
  changing a return type or shape, breaks printed code. If you add a public function, regenerate
  `src/tests/fixtures/api_contract.json` with `make_api_contract.py`.
- **Numbers may drift; behaviour may not.** A change that alters results at floating-point level
  (for example `solve` instead of `inv`) is fine if it is explained and within the regression
  tolerances. A change to what a function computes (normalisation, defaults, how inputs enter the
  model) is not a refactor; open an issue first.
- **Failed control problems return values; they never raise.** An ill-conditioned or incomplete
  transition returns its energies and error terms (the paper tells users how to read them). Only
  input that does not define a problem may raise.
- **Never alter a caller's inputs.** For example, nothing in nctpy zeroes the diagonal of a
  connectome the user passes in.
- **Never regenerate a fixture to make a test pass.** Say what changed and why, and leave
  regeneration to a maintainer.
- **No data in the repository.** Fixtures are synthetic. pre-commit refuses files under `data/`
  and files over 1 MB.
- **Keep pull requests small**, ideally one module at a time, so any change in results can be
  traced to a single change.
- **Add a changelog entry** under `## Unreleased` in `CHANGELOG.md` for anything a user would
  notice.

## Reporting problems

Please open an issue using one of the templates. For a bug, include your nctpy, Python, numpy and
SciPy versions and a minimal example that reproduces it. If you are following the protocol paper,
say which step.

## Questions

For questions, contact Linden Parkes and Jason Kim: info@parkeslab.com.
