## What this changes

<!-- One module per pull request where possible. What changed, and why? -->

## Effect on results (required)

<!-- Pick one and explain. See CONTRIBUTING.md, "Rules for changes". -->

- [ ] **No results change.** The regression fixtures pass unchanged.
- [ ] **Results change at floating-point level only.** Which outputs, by how much, and why
      (e.g. `solve` instead of `inv`):
- [ ] **Behaviour changes.** What changes for users, and the issue where this was agreed:

Regression fixtures regenerated? <!-- no / yes: which, and why. Maintainers only; needs sign-off. -->

## Public API

- [ ] No public function, signature, default, or return type changes.
- [ ] Adds public API. New keyword arguments reproduce the old behaviour by default, and
      `src/tests/fixtures/api_contract.json` is regenerated.

## Checks

- [ ] Tests pass locally (`cd src/tests && python -m pytest`).
- [ ] If this could affect the protocol paper's code: ran with the paper's data, so
      `test_paper_code.py` and `test_notebooks.py` did not skip. Maintainers can do this.
- [ ] `CHANGELOG.md` has an entry under `## Unreleased` if users would notice this change.
- [ ] No data files are added (fixtures are synthetic).
