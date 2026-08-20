# Releasing nctpy to PyPI

Package metadata lives in **`setup.cfg` only**. `pyproject.toml` just declares the
setuptools build backend. There is deliberately no `setup.py` — when one existed, its
keyword arguments silently overrode `setup.cfg`, so version bumps made in `setup.cfg`
were ignored at build time.

## Steps

1. Bump `version` in `setup.cfg` (and `version` in `CITATION.cff` to match).
2. Make sure `main` is clean and up to date:

   ```
   git checkout main && git pull && git status
   ```

3. Run the unit tests:

   ```
   python -m unittest discover -s src/tests -t .
   ```

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
   in `setup.cfg`.

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
