"""Fail if any module's coverage falls below the recorded baseline (baseline.json, "public").

Usage, after `pytest --cov=nctpy --cov=null_models --cov-branch --cov-report=json:coverage.json`:

    python check_coverage.py coverage.json

The baseline is a floor: a change may raise coverage but must not lower it. The
"public" figures are the ones CI can reach, since tests that need the paper's data skip there.
"""

import json
import sys
from pathlib import Path

BASELINE = Path(__file__).resolve().parent / "baseline.json"
TOLERANCE = 0.1  # percentage points, for rounding


def module_key(path):
    """'.../site-packages/nctpy/energies.py' or 'src/nctpy/energies.py' -> 'nctpy/energies.py'."""
    path = path.replace("\\", "/")
    starts = [path.rfind(root) for root in ("nctpy/", "null_models/")]
    return path[max(starts) :] if max(starts) >= 0 else path


def main(coverage_json):
    floor = json.loads(BASELINE.read_text())["public"]["modules"]
    measured = {
        module_key(f): d["summary"]["percent_covered"] for f, d in json.load(open(coverage_json))["files"].items()
    }
    failures = []
    print("{0:28s} {1:>9s} {2:>9s}".format("module", "baseline", "now"))
    for module, base in sorted(floor.items()):
        now = measured.get(module)
        status = "MISSING" if now is None else ("BELOW" if now < base["percent"] - TOLERANCE else "")
        print(
            "{0:28s} {1:8.1f}% {2:>9s} {3}".format(
                module, base["percent"], "-" if now is None else "{0:.1f}%".format(now), status
            )
        )
        if status:
            failures.append(module)
    if failures:
        print("\ncoverage fell below the baseline for: " + ", ".join(failures))
        return 1
    print("\nno module is below its baseline")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
