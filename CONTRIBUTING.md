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
- **docs:** the Sphinx documentation must build without warnings;
- **docs notebooks:** the documentation's notebooks (Getting started, two of the examples and every tutorial) must
  execute against the current code.

Some existing files are exempt from parts of the ruff checks (see `pyproject.toml`). These
exemptions are being removed module by module. New files must pass in full.

## Rules for changes

- **Never change a public signature.** Within 1.x, new behaviour arrives as new keyword arguments
  whose defaults reproduce the existing behaviour exactly. Renaming or removing a parameter, or
  changing a return type or shape, breaks printed code. If you add a public function, regenerate
  `src/tests/fixtures/api_contract.json` with `make_api_contract.py`, and list the function on its
  page of the API reference (see below).
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

## Documentation

The documentation is built with Sphinx, pydata-sphinx-theme, MyST Markdown and myst-nb. Its sources are in
`docs/source/`.

### Building it locally

```bash
pip install -e ".[docs]"
python -m sphinx -b html docs/source docs/_build/html
```

Then open `docs/_build/html/index.html`, or serve the folder (`python -m http.server -d docs/_build/html`). For a
build that updates as you edit, `pip install sphinx-autobuild` and run
`sphinx-autobuild docs/source docs/_build/html`. CI builds with `-W --keep-going`, so warnings are errors; do the
same before opening a pull request.

The build only renders notebooks, with the outputs committed in them; it never executes them.

### Docstrings and the API reference

Docstrings use NumPy style. Cross-reference with roles such as ``{func}`~nctpy.energies.get_control_inputs` `` so
that the pages link up.

The API reference pages, `docs/source/api/*.md`, list the public functions and classes by area, and autosummary
builds a page for each from its docstring. nctpy's modules have no `__all__`; the list of public symbols is
`src/tests/fixtures/api_contract.json`. `src/tests/test_docs.py` fails if a public symbol is missing from the API
pages, or if they list one that is not public.

### Notebooks

Getting started (`docs/source/pages/getting_started/index.ipynb`), the tutorials (`docs/source/tutorials/`) and two
of the examples are notebooks that CI executes. The other examples analyse real data and are static pages.

- **Use synthetic, seeded data only.** Nothing derived from real data may be committed, and these notebooks must
  run anywhere.
- **Commit them with their outputs.** Read the Docs shows what is committed. After changing a notebook, or code it
  depends on, refresh its outputs:

  ```bash
  pip install -e ".[docs,plot,optimize]"
  python docs/execute_notebooks.py --write tutorials/your_tutorial.ipynb   # path relative to docs/source
  ```

  Without `--write` the script only checks that the notebooks run, as CI does. CI does not check that the
  committed outputs are current.
- **Keep the outputs reproducible**, so that refreshing a notebook changes only what really changed. The script
  turns off progress bars and merges each cell's printed output; avoid unseeded randomness in figures too (for
  example the bootstrapped confidence band of a regression line, or the jitter of a strip plot).
- **Keep each notebook to a minute or so.** Reduce permutation counts and say so in the text.
- A new tutorial goes in `docs/source/tutorials/` and in the toctree of `docs/source/tutorials/index.md`; the script
  finds it automatically.

### The "Reproducing the papers" page

`docs/source/pages/reproducing.md` gives the protocol paper's main-text code with its typesetting repairs. The code
must be exactly the blocks in `src/tests/test_paper_code.py` with their repairs applied; `test_docs.py` checks this.

## Performance

`benchmarks/transitions.py` times control energy over all transitions between seven brain states. Run it before and
after a change that could affect speed, with the same number of threads (e.g. `OMP_NUM_THREADS=8`).

## Reporting problems

Please open an issue using one of the templates. For a bug, include your nctpy, Python, numpy and
SciPy versions and a minimal example that reproduces it. If you are following the protocol paper,
say which step.

## Questions

For questions, contact Linden Parkes and Jason Kim: info@parkeslab.com.
