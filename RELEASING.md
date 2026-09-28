# Releasing nctpy to PyPI

Package metadata lives in **`pyproject.toml` only**. There is deliberately no `setup.py` or
`setup.cfg`. When those existed, they could silently override one another, so a version bump made
in one place was ignored at build time.

The version has a single source: `__version__` in `src/nctpy/__init__.py`, which `pyproject.toml`
reads at build time.

## Steps

1. Bump `__version__` in `src/nctpy/__init__.py` (and `version` in `CITATION.cff` to match),
   and date the release's section in `CHANGELOG.md`.
2. Make sure `main` is clean and up to date:

   ```
   git checkout main && git pull && git status
   ```

3. Run the tests, from `src/tests` (some tests open `./fixtures` relative to the working
   directory):

   ```
   cd src/tests && PYTHONPATH=.. python -m pytest && cd ../..
   ```

   Run this on a machine with the paper's data in `data/`. The paper-code and notebook tests
   need it and skip without it, which is what public CI sees.

4. Build a fresh sdist + wheel (clear old artifacts first, or twine will try to
   upload every version sitting in `dist/`):

   ```
   rm -rf dist build src/*.egg-info
   python -m pip install --upgrade build twine
   python -m build
   ```

5. Sanity-check the metadata before uploading:

   ```
   python -m twine check dist/*
   unzip -p dist/nctpy-*.whl '*/METADATA' | head -30
   ```

   Confirm the version is the new one and that `Requires-Dist` lists every dependency
   in `pyproject.toml`. The wheel should contain only `nctpy/` and `null_models/`:

   ```
   unzip -l dist/nctpy-*.whl
   ```

6. Upload. Use an API token from https://pypi.org/manage/account/token/
   (username is the literal string `__token__`, password is the `pypi-...` token).

   ```
   python -m twine upload dist/*
   ```

   To rehearse without touching the real index, upload to TestPyPI first:

   ```
   python -m twine upload --repository testpypi dist/*
   ```

7. Tag and push the release:

   ```
   git tag -a v1.0.2 -m "v1.0.2"
   git push origin main --tags
   ```

8. Optionally cut a GitHub release from the tag at
   https://github.com/LindenParkesLab/nctpy/releases/new

## Notes

- PyPI versions are immutable. Once `1.0.2` is uploaded it can never be replaced —
  a mistake means yanking it and shipping `1.0.3`.
- Credentials can be stored in `~/.pypirc` so you are not prompted each time:

  ```
  [pypi]
    username = __token__
    password = pypi-AgEIcHlwaS5vcmc...
  ```

  Keep that file at `chmod 600`.
