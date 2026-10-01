"""Time control energy over all pairwise transitions between brain states.

Mirrors the protocol paper's use: 7 binary states, all 49 transitions, a uniform full control set, one
synthetic connectome. Run it before and after a performance change and compare:

    python benchmarks/transitions.py              # 200 nodes, as in the paper
    python benchmarks/transitions.py 400 100      # 400 nodes, 100 nodes for minimum_energy_fast

It prints seconds per workload and per transition. BLAS threading affects the numbers; set e.g.
OMP_NUM_THREADS to compare like with like.
"""

import contextlib
import io
import platform
import sys
import time

import numpy as np
import scipy

import nctpy
from nctpy.energies import get_control_inputs, integrate_u, minimum_energy_fast
from nctpy.pipelines import ComputeControlEnergy
from nctpy.utils import matrix_normalization, normalize_state

N_STATES = 7


def setup(n):
    rng = np.random.default_rng(0)
    A = rng.random((n, n))
    A = (A + A.T) / 2
    np.fill_diagonal(A, 0)
    states = np.repeat(np.arange(N_STATES), int(np.ceil(n / N_STATES)))[:n]
    pairs = [
        (normalize_state(states == i), normalize_state(states == j)) for i in range(N_STATES) for j in range(N_STATES)
    ]
    return A, pairs


def loop(A, pairs, system, T):
    A_norm = matrix_normalization(A, system=system)
    B, S = np.eye(len(A)), np.eye(len(A))
    for x0, xf in pairs:
        _, u, _ = get_control_inputs(A_norm, T, B, x0, xf, system=system, rho=1, S=S)
        np.sum(integrate_u(u))


def pipeline(A, pairs, system, T):
    n = len(A)
    tasks = [dict(x0=x0, xf=xf, B=np.eye(n), S=np.eye(n), rho=1) for x0, xf in pairs]
    with contextlib.redirect_stderr(io.StringIO()):
        ComputeControlEnergy(A=A, control_tasks=tasks, system=system, T=T).run()


def minimum(A, pairs, T):
    A_norm = matrix_normalization(A, system="continuous")
    for x0, xf in pairs:
        np.sum(minimum_energy_fast(A_norm, T, np.eye(len(A)), x0, xf))


def timed(fn, *args):
    start = time.perf_counter()
    fn(*args)
    return time.perf_counter() - start


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    n_min = int(sys.argv[2]) if len(sys.argv) > 2 else n
    print(f"nctpy {nctpy.__version__} ({nctpy.__file__}), numpy {np.__version__}, scipy {scipy.__version__}, "
          f"python {platform.python_version()}")  # fmt: skip
    A, pairs = setup(n)
    A_min, pairs_min = setup(n_min)
    workloads = [
        (f"get_control_inputs loop, continuous T=1, {n} nodes", loop, (A, pairs, "continuous", 1)),
        (f"get_control_inputs loop, discrete T=3, {n} nodes", loop, (A, pairs, "discrete", 3)),
        (f"ComputeControlEnergy, continuous T=1, {n} nodes", pipeline, (A, pairs, "continuous", 1)),
        (f"minimum_energy_fast loop, T=1, {n_min} nodes", minimum, (A_min, pairs_min, 1)),
    ]
    for label, fn, args in workloads:
        seconds = timed(fn, *args)
        print(f"{label:55s} {seconds:8.2f} s  ({seconds / len(pairs) * 1e3:7.1f} ms per transition)")


if __name__ == "__main__":
    main()
