"""Execute the documentation's notebooks: Getting started, the two method examples and every tutorial.

Read the Docs only renders the notebooks, with the outputs committed in them; CI runs this script so that a change
that breaks a notebook fails there instead of leaving a stale page. Every notebook run here uses synthetic, seeded
data only. The other examples under pages/examples/ use real data; they are a static gallery and are never
executed.

    python docs/execute_notebooks.py           # execute; fail if any cell raises
    python docs/execute_notebooks.py --write   # also save the fresh outputs into the notebooks, to commit
    python docs/execute_notebooks.py --write tutorials/decay_rates.ipynb   # only the notebooks named

Needs nctpy with the docs, plot and optimize extras.
"""

import argparse
import os
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

SOURCE = Path(__file__).resolve().parent / "source"
NOTEBOOKS = [
    SOURCE / "pages" / "getting_started" / "index.ipynb",
    SOURCE / "pages" / "examples" / "metric_correlations.ipynb",
    SOURCE / "pages" / "examples" / "minimum_energy_fast.ipynb",
    *sorted((SOURCE / "tutorials").glob("*.ipynb")),
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--write", action="store_true", help="save the executed notebooks in place")
    parser.add_argument("notebooks", nargs="*", type=Path, help="paths relative to docs/source (default: all)")
    args = parser.parse_args()
    notebooks = [SOURCE / path for path in args.notebooks] or NOTEBOOKS

    # A forced backend (CI sets MPLBACKEND=Agg) stops the kernel's inline backend from emitting figures as outputs.
    os.environ.pop("MPLBACKEND", None)
    # Progress bars print their timings, which would change the committed outputs on every run.
    os.environ["TQDM_DISABLE"] = "1"

    failures = []
    for path in notebooks:
        nb = nbformat.read(path, as_version=4)
        # No per-cell timings, and each cell's printed output in one piece however the kernel happened to flush
        # it, so that re-executing a notebook does not change its committed outputs.
        client = NotebookClient(
            nb,
            timeout=900,
            kernel_name="python3",
            record_timing=False,
            coalesce_streams=True,
            resources={"metadata": {"path": path.parent}},
        )
        started = time.time()
        try:
            client.execute()
        except CellExecutionError as exc:
            failures.append(path.relative_to(SOURCE))
            print(f"FAIL  {time.time() - started:6.1f}s  {path.relative_to(SOURCE)}\n{str(exc)[-2500:]}", flush=True)
            continue
        print(f"PASS  {time.time() - started:6.1f}s  {path.relative_to(SOURCE)}", flush=True)
        if args.write:
            nbformat.write(nb, path)

    if failures:
        print(f"\nNotebooks failed to execute: {', '.join(map(str, failures))}")
        return 1
    print(f"\nAll {len(notebooks)} notebook(s) executed cleanly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
