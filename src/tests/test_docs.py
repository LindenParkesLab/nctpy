"""Documentation pages that must stay in step with the code.

- The API reference lists every public symbol, and nothing else. nctpy's modules have no ``__all__``, so the list of
  public symbols is the API contract (fixtures/api_contract.json). The pages under docs/source/api/ are written by
  hand; the test fails when a symbol is added to the contract without a page entry, or a page lists a name that is
  not public.
- The "Reproducing the papers" page gives the protocol paper's main-text code with the typesetting repairs that
  test_paper_code.py applies, and nothing else changed.

Skipped where the docs are not present (e.g. testing an installed wheel).
"""

import json
import re
import unittest
from pathlib import Path

FIXTURES = Path(__file__).resolve().parent / "fixtures"
DOCS = Path(__file__).resolve().parents[2] / "docs" / "source"
API = DOCS / "api"
REPRODUCING = DOCS / "pages" / "reproducing.md"


def documented_symbols():
    """(module, name) of every entry in the autosummary blocks of the API pages."""
    found = set()
    for page in API.glob("*.md"):
        for block in re.findall(r"```\{eval-rst\}\n(.*?)```", page.read_text(), re.S):
            module = re.search(r"^\.\. currentmodule:: (\S+)$", block, re.M).group(1)
            names = re.findall(r"^   (?!:)(\S+)$", block.split(".. autosummary::", 1)[1], re.M)
            found.update((module, name) for name in names)
    return found


class TestAPIReference(unittest.TestCase):
    def setUp(self):
        if not API.is_dir():
            self.skipTest("docs/source/api not present")

    def test_pages_match_the_contract(self):
        contract = json.loads((FIXTURES / "api_contract.json").read_text())["symbols"]
        public = {(s["module"], s["name"]) for s in contract}
        documented = documented_symbols()
        self.assertEqual(sorted(public - documented), [], "public but missing from docs/source/api/")
        self.assertEqual(sorted(documented - public), [], "on docs/source/api/ but not in the API contract")


class TestReproducingPage(unittest.TestCase):
    def setUp(self):
        if not REPRODUCING.exists():
            self.skipTest("docs/source/pages/reproducing.md not present")

    def test_code_is_the_papers_with_its_repairs(self):
        import test_paper_code as paper

        blocks = [
            (paper.MAIN_IMPORTS, paper.MAIN_IMPORTS_REPAIRS),
            (paper.MAIN_LOAD, ()),
            (paper.MAIN_STEP1, ()),
            (paper.MAIN_STEP2, ()),
            (paper.MAIN_STEP3A, ()),
            (paper.MAIN_STEP3B, ()),
            (paper.MAIN_STEP3C, ()),
            (paper.MAIN_STEP3D, ()),
            (paper.MAIN_STEP3E, ()),
            (paper.MAIN_STEP4A, paper.MAIN_STEP4A_REPAIRS),
            (paper.MAIN_STEP4IV, paper.MAIN_STEP4IV_REPAIRS),
            (paper.MAIN_BOX1, paper.MAIN_BOX1_REPAIRS),
            (paper.MAIN_STEP6A, ()),
            (paper.MAIN_STEP6B, ()),
            (paper.MAIN_PROCEDURE2, ()),
        ]
        expected = [paper.edit(text, repairs).rstrip("\n") for text, repairs in blocks]
        on_page = [code.rstrip("\n") for code in re.findall(r"```python\n(.*?)```", REPRODUCING.read_text(), re.S)]
        self.assertEqual(on_page, expected)


if __name__ == "__main__":
    unittest.main()
