# Releasing nctpy

Package metadata lives in **`pyproject.toml` only**. There is deliberately no `setup.py` or
`setup.cfg`. When those existed, they could silently override one another, so a version bump made
in one place was ignored at build time.

The version has a single source: `__version__` in `src/nctpy/__init__.py`, which `pyproject.toml`
reads at build time.

Releases are automated. Pushing a tag `vX.Y.Z` runs `.github/workflows/release.yml`, which:

1. checks the tag agrees with `__version__`, `CITATION.cff` and a dated `CHANGELOG.md` section, and
   that the tagged commit is on `main`;
2. builds the sdist and wheel, and tests the wheel;
3. publishes to PyPI (trusted publishing: no API token is stored anywhere);
4. creates the GitHub release, with the changelog section as its notes and the built files
   attached. This also triggers the Zenodo archive.

## One-time setup

These are done once, on the websites, by a maintainer.

1. **PyPI trusted publisher.** On https://pypi.org/manage/project/nctpy/settings/publishing/, add
   a GitHub publisher:
   - owner `LindenParkesLab`
   - repository `nctpy`
   - workflow `release.yml`
   - environment `pypi`
2. **GitHub environment.** In the repository's Settings → Environments, create an environment
   named `pypi`. Optionally, add yourself as a required reviewer, so every publish waits for a
   click of approval.
3. **Zenodo** (already enabled, going by the README's DOI badge). At https://zenodo.org/account/settings/github/,
   the `nctpy` repository should be switched on.

## Steps for a release

1. **Prepare the release on a branch**, and merge it into `main`:
   - bump `__version__` in `src/nctpy/__init__.py`;
   - bump `version` and `date-released` in `CITATION.cff`;
   - rename `## Unreleased` in `CHANGELOG.md` to `## X.Y.Z (YYYY-MM-DD)`, with the same date as
     `date-released`.

   Check it locally:

   ```
   python .github/scripts/release_check.py vX.Y.Z /tmp/release_notes.md
   ```

2. **Run the full tests locally**, on a machine with the paper's data in `data/`. CI cannot run the
   paper-code and notebook tests, because the data are not in the repository:

   ```
   cd src/tests && PYTHONPATH=.. python -m pytest && cd ../..
   ```

3. **Optionally, do a dry run:** Actions → Release → Run workflow (on `main`). It builds and tests
   the package and publishes nothing.

4. **Tag `main` and push the tag.** In VS Code: switch to `main` and pull, then use Command Palette
   → "Git: Create Tag" (name it `vX.Y.Z`), then "Git: Push Tags". Or in a terminal:

   ```
   git switch main && git pull
   git tag -a vX.Y.Z -m "vX.Y.Z"
   git push origin vX.Y.Z
   ```

5. **Watch Actions → Release.** When it is green, the new version is on
   https://pypi.org/project/nctpy/ and on the GitHub releases page.

## If something goes wrong

- **The release check or the tests fail:** nothing has been published. Fix the problem on `main`,
  delete the tag (locally with `git tag -d vX.Y.Z`, and on GitHub under the repository's tags),
  and tag again.
- **PyPI versions are immutable.** Once `X.Y.Z` is uploaded it can never be replaced. A mistake
  means yanking it on PyPI and releasing the next version.
- **Manual fallback**, if the workflow is unavailable:
  - build with `python -m build`;
  - check with `python -m twine check --strict dist/*`;
  - upload with `python -m twine upload dist/*`, using an API token from
    https://pypi.org/manage/account/token/ (username `__token__`);
  - then create the GitHub release from the tag by hand.
