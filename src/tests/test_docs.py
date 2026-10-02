"""The API reference lists every public symbol, and nothing else.

nctpy's modules have no ``__all__``, so the list of public symbols is the API contract
(fixtures/api_contract.json). The pages under docs/source/api/ are written by hand; this test fails when a symbol is
added to the contract without a page entry, or a page lists a name that is not public. Skipped where the docs are
not present (e.g. testing an installed wheel).
"""

import json
import re
import unittest
from pathlib import Path

FIXTURES = Path(__file__).resolve().parent / "fixtures"
API = Path(__file__).resolve().parents[2] / "docs" / "source" / "api"


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


if __name__ == "__main__":
    unittest.main()
