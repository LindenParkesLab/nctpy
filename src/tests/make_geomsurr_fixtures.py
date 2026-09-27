"""Capture geomsurr's surrogates as a reference fixture (fixtures/geomsurr_reference.npz).

The surrogates geomsurr returns for a given seed are part of the frozen API: the protocol paper
prints the call. This fixture pins them so that changes to how geomsurr is implemented (e.g. no
longer modifying the caller's W, or no longer touching numpy's global random state) can be shown
to leave the surrogates bit-for-bit unchanged.

It was first generated from the pre-1.0.3 implementation, before either of those changes. Run it
again only to deliberately re-pin, from src/tests:

    python make_geomsurr_fixtures.py
"""
import subprocess

import numpy as np
import scipy as sp
from scipy.spatial.distance import pdist, squareform

from null_models.geomsurr import geomsurr

N_NODES = 30
SEEDS = (0, 1, 123)


def make_case(directed, seed=7):
    """Positive weights decaying with distance, ~70% dense. The undirected case keeps
    self-connections, like the protocol paper's PNC connectome."""
    rng = np.random.default_rng(seed)
    D = squareform(pdist(rng.normal(size=(N_NODES, 3)) * 40.0))
    W = np.exp(-D / 30.0 + rng.normal(scale=0.4, size=(N_NODES, N_NODES)))
    W[rng.random((N_NODES, N_NODES)) < 0.3] = 0.0
    if directed:
        np.fill_diagonal(W, 0.0)
    else:
        W = np.tril(W) + np.tril(W, -1).T
        np.fill_diagonal(W, rng.uniform(1.0, 2.0, size=N_NODES))
    return W, D


def main():
    arrays = {}
    for label, directed in (("und", False), ("dir", True)):
        W, D = make_case(directed)
        arrays[f"W_{label}"] = W
        arrays[f"D_{label}"] = D
        for seed in SEEDS:
            # pass a copy: the pre-1.0.3 implementation zeroes W's diagonal in place
            Wwp, Wsp, Wssp = geomsurr(W.copy(), D, seed=seed)
            arrays[f"wwp_{label}_{seed}"] = Wwp
            arrays[f"wsp_{label}_{seed}"] = Wsp
            arrays[f"wssp_{label}_{seed}"] = Wssp

    sha = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    np.savez_compressed("./fixtures/geomsurr_reference.npz", nctpy_commit=sha,
                        numpy_version=np.__version__, scipy_version=sp.__version__, **arrays)
    print(f"wrote fixtures/geomsurr_reference.npz from nctpy {sha[:7]} "
          f"(numpy {np.__version__}, scipy {sp.__version__}), seeds {SEEDS}")


if __name__ == "__main__":
    main()
