"""Generate the regression fixtures in fixtures/regression/ (Roadmap 0.4; see docs TESTING section 3).

The fixtures are a drift alarm, not the definition of correct: they record what nctpy returns
today so that any later change to a result must be noticed and explained. Floating-point drift
within tolerance is acceptable; a behavioural change is not.

Contents:

- one .npz per connectome, all synthetic and seeded (no real data is used or stored: the paper's
  connectomes are not in the repository and may not be shared):
  - undirected: 60 nodes, dense, with self-connections (the PNC connectome's properties)
  - directed: 43 nodes, fully connected, with self-connections (the mouse isocortex's properties)
  - disconnected: 12 nodes in two unconnected 6-node blocks
  - legacy: the older suite's fixtures/A.npy; 50 nodes, symmetric, negative weights
- utils.npz for the connectome-independent helpers.
- manifest.json: environment, seeds, the control-task design, every cell's status, and a
  sha256 of each .npz's contents.

Control tasks: the full grid (system x c x T x rho x S x xr x states x B x expm_version) has
over 20,000 cells per connectome, so a seeded greedy covering design is used instead, in which
every pair of parameter values appears in at least one cell. The same design is applied to every
connectome. rho=0 is recorded in separate cells: it currently returns NaN (decision D2 will turn
it into a ValueError). Cells that raise are recorded with the exception type and message.

Run from src/tests. Regenerating replaces every fixture and the manifest; do it only for a
documented reason, with sign-off (TESTING.md, "Regenerating fixtures"):

    python make_regression_fixtures.py
"""
import hashlib
import json
import platform
import subprocess
import sys
import warnings
from datetime import date
from itertools import combinations
from pathlib import Path

import numpy as np
import scipy as sp
from scipy.spatial.distance import pdist, squareform

from nctpy.energies import sim_state_eq, get_control_inputs, integrate_u, gramian, minimum_energy_fast
from nctpy.metrics import ave_control, modal_control
from nctpy.utils import (matrix_normalization, normalize_state, normalize_weights, convert_states_str2int,
                         get_null_p, get_fdr_p, get_p_val_string, expm)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE / 'fixtures' / 'regression'

SEED = 20260927
N_TRAJECTORY_POINTS = 11  # x(t) is stored at this many evenly spaced time points, x0 and x(T) included
MAX_FULL_MATRIX = 100  # larger N x N outputs are stored as their diagonal plus a seeded sample of entries
N_SAMPLED_ENTRIES = 1000

FACTORS = {
    'system': ['continuous', 'discrete'],
    'c': [0.5, 1, 5],
    'T': ['short', 'medium', 'long'],
    'rho': [0.5, 1, 100],
    'S': ['identity', 'zeros', 'partial'],
    'xr': ['zero', 'x0', 'xf', 'midpoint'],
    'states': ['binary', 'nonbinary'],
    'B': ['identity', 'partial', 'weighted'],
    'expm_version': ['scipy', 'eig'],
}
# discrete time counts steps, so T takes different values per system
T_VALUES = {'continuous': {'short': 0.1, 'medium': 1, 'long': 10},
            'discrete': {'short': 2, 'medium': 5, 'long': 10}}
RHO_ZERO_CELLS = [dict(system=system, c=1, T='medium', rho=0, S=S, xr='zero', states='binary', B='identity',
                       expm_version='scipy')
                  for system in ('continuous', 'discrete') for S in ('identity', 'zeros')]


def pairwise_design(factors, seed, n_candidates=50):
    """Greedy (AETG-style) covering design: every pair of factor levels appears in some row."""
    names = list(factors)
    levels = [len(factors[n]) for n in names]
    k = len(names)

    def pairs_of(row):
        return {((i, row[i]), (j, row[j])) for i, j in combinations(range(k), 2)}

    uncovered = {((i, a), (j, b)) for i, j in combinations(range(k), 2)
                 for a in range(levels[i]) for b in range(levels[j])}
    rng = np.random.default_rng(seed)
    rows = []
    while uncovered:
        pool = sorted(uncovered)
        best, best_gain = None, -1
        for _ in range(n_candidates):
            (i, a), (j, b) = pool[rng.integers(len(pool))]
            row = [None] * k
            row[i], row[j] = a, b
            for f in rng.permutation(k):
                if row[f] is not None:
                    continue
                gains = [sum(((min(f, h), v if f < h else row[h]), (max(f, h), row[h] if f < h else v)) in uncovered
                             for h in range(k) if row[h] is not None) for v in range(levels[f])]
                choices = [v for v, g in enumerate(gains) if g == max(gains)]
                row[f] = choices[rng.integers(len(choices))]
            gain = len(pairs_of(row) & uncovered)
            if gain > best_gain:
                best, best_gain = row, gain
        rows.append(best)
        uncovered -= pairs_of(best)
    return [{name: factors[name][v] for name, v in zip(names, row)} for row in rows]


# ---------------------------------------------------------------------------------------------
# connectomes and their states
# ---------------------------------------------------------------------------------------------

def spatial_weights(rng, n, directed):
    """Positive weights decaying with distance between random 3-D node positions, with multiplicative noise."""
    D = squareform(pdist(rng.normal(size=(n, 3)) * 40.0))
    W = 100 * np.exp(-D / 40.0 + rng.normal(scale=0.5, size=(n, n)))
    if not directed:
        W = np.tril(W) + np.tril(W, -1).T
    return W


def modules(n, k):
    """Assign n nodes to k contiguous modules of near-equal size."""
    return np.repeat(np.arange(k), np.diff(np.linspace(0, n, k + 1).round().astype(int)))


def load_connectomes(rng):
    """Each entry: adjacency, a binary (x0, xf) pair, a non-binary (x0, xf) pair, and a description.

    Everything is synthetic and seeded. The real connectomes the paper uses (data/) are not in the
    repository and may not be shared, so nothing here is derived from them; the synthetic matrices
    are built to have the properties that matter to the code paths (dense and undirected with
    self-connections, like the PNC connectome; directed and fully connected with self-connections,
    like the mouse isocortex).
    """
    out = {}

    n = 60
    A = spatial_weights(rng, n, directed=False)
    drop = np.triu(rng.random((n, n)) < 0.1, k=1)
    A[drop | drop.T] = 0  # ~90% dense
    np.fill_diagonal(A, rng.uniform(0.5, 2.0, size=n) * np.median(A))  # positive self-connections
    m = modules(n, 6)
    out['undirected'] = dict(A=A, binary=(m == 0, m == 3), nonbinary=(rng.normal(size=n), rng.normal(size=n)),
                             source='synthetic: 60 nodes, undirected, ~90% dense, distance-decaying lognormal '
                                    'weights, positive self-connections; binary module 0 -> 3 of 6; '
                                    'non-binary seeded normal')

    n = 43
    A = spatial_weights(rng, n, directed=True)
    np.fill_diagonal(A, rng.uniform(0.5, 2.0, size=n) * np.median(A))
    m = modules(n, 6)
    out['directed'] = dict(A=A, binary=(m == 0, m == 2), nonbinary=(rng.normal(size=n), rng.normal(size=n)),
                           source='synthetic: 43 nodes, directed, fully connected, distance-decaying lognormal '
                                  'weights, positive self-connections; binary module 0 -> 2 of 6; '
                                  'non-binary seeded normal')

    n = 12
    A = np.zeros((n, n))
    for block in (slice(0, 6), slice(6, 12)):
        W = rng.random((6, 6)) + 0.1
        A[block, block] = (W + W.T) / 2
    np.fill_diagonal(A, 0)
    x0, xf = np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    x0[:3], xf[-3:] = True, True  # from one component to the other
    out['disconnected'] = dict(A=A, binary=(x0, xf), nonbinary=(rng.normal(size=n), rng.normal(size=n)),
                               source='synthetic, two unconnected 6-node blocks; binary nodes 0-2 -> 9-11 '
                                      '(across components); non-binary seeded normal')

    A = np.load(HERE / 'fixtures' / 'A.npy')
    n = len(A)
    x0, xf = np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    x0[:10], xf[-10:] = True, True
    out['legacy'] = dict(A=A, binary=(x0, xf), nonbinary=(rng.normal(size=n), rng.normal(size=n)),
                         source='src/tests/fixtures/A.npy (symmetric, negative weights); binary nodes 0-9 -> 40-49; '
                                'non-binary seeded normal')
    return out


def control_matrices(n, rng):
    partial = np.ones(n)
    partial[::5] = 0  # every fifth node gets no control
    s_partial = np.zeros(n)
    s_partial[:n // 2] = 1  # constrain the first half of the trajectory
    return {'B': {'identity': np.eye(n), 'partial': np.diag(partial),
                  'weighted': np.diag(normalize_weights(rng.normal(size=n)))},
            'S': {'identity': np.eye(n), 'zeros': np.zeros((n, n)), 'partial': np.diag(s_partial)}}


# ---------------------------------------------------------------------------------------------
# computing
# ---------------------------------------------------------------------------------------------

def store_matrix(arrays, key, M, rng):
    """Full matrix if small, else its diagonal plus a seeded sample of entries."""
    if M.shape[0] <= MAX_FULL_MATRIX:
        arrays[key] = M
    else:
        flat = rng.choice(M.size, size=N_SAMPLED_ENTRIES, replace=False)
        arrays[key + '__diag'] = np.diag(M).copy()
        arrays[key + '__sample_index'] = flat
        arrays[key + '__sample_value'] = M.ravel()[flat]


def run_cell(cell, conn, mats):
    A_norm = matrix_normalization(conn['A'], system=cell['system'], c=cell['c'])
    x0, xf = (normalize_state(s) for s in conn[cell['states']])
    T = T_VALUES[cell['system']][cell['T']]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        x, u, err = get_control_inputs(A_norm, T, mats['B'][cell['B']], x0, xf, system=cell['system'],
                                       rho=cell['rho'], S=mats['S'][cell['S']], xr=cell['xr'],
                                       expm_version=cell['expm_version'])
        node_energy = integrate_u(u)
    idx = np.unique(np.linspace(0, len(x) - 1, N_TRAJECTORY_POINTS).round().astype(int))
    return {'node_energy': node_energy, 'err': np.asarray(err, dtype=float), 'x_index': idx, 'x': x[idx]}, \
        sorted({type(w.message).__name__ + ': ' + str(w.message) for w in caught})


def connectome_fixture(name, conn, design, rng):
    arrays, cells = {}, []
    n = len(conn['A'])
    mats = control_matrices(n, rng)
    for key, M in mats['B'].items():
        arrays['B_' + key] = np.diag(M).copy()
    for key, M in mats['S'].items():
        arrays['S_' + key] = np.diag(M).copy()
    for label in ('binary', 'nonbinary'):
        arrays['states_' + label + '_x0'] = np.asarray(conn[label][0])
        arrays['states_' + label + '_xf'] = np.asarray(conn[label][1])

    # normalisation and metrics
    arrays['lambda_max'] = np.array(np.abs(np.linalg.eigvals(conn['A'])).max())
    A_c = matrix_normalization(conn['A'], system='continuous', c=1)
    A_d = matrix_normalization(conn['A'], system='discrete', c=1)
    store_matrix(arrays, 'A_norm_continuous', A_c, rng)
    store_matrix(arrays, 'A_norm_discrete', A_d, rng)
    arrays['ave_control_continuous'] = ave_control(A_c, system='continuous')
    arrays['ave_control_discrete'] = ave_control(A_d, system='discrete')
    arrays['modal_control'] = modal_control(A_d)
    store_matrix(arrays, 'gramian_continuous_T1', gramian(A_c, 1, system='continuous'), rng)
    store_matrix(arrays, 'gramian_discrete_T5', gramian(A_d, 5, system='discrete'), rng)
    x0, xf = (normalize_state(s) for s in conn['binary'])
    arrays['minimum_energy_fast_T1'] = np.asarray(minimum_energy_fast(A_c, 1, np.eye(n), x0, xf))
    U = rng.normal(size=(n, 20)) * 0.1
    arrays['sim_state_eq_continuous'] = sim_state_eq(A_c, np.eye(n), x0, U, system='continuous')
    arrays['sim_state_eq_discrete'] = sim_state_eq(A_d, np.eye(n), x0, U, system='discrete')
    arrays['sim_state_eq_U'] = U

    # control tasks
    for i, cell in enumerate(list(design) + RHO_ZERO_CELLS):
        key = 'cell{0:03d}'.format(i)
        entry = dict(cell, key=key, T_value=T_VALUES[cell['system']][cell['T']])
        try:
            result, warns = run_cell(cell, conn, mats)
            for part, value in result.items():
                arrays[key + '__' + part] = value
            entry['status'] = 'returned'
            entry['finite'] = bool(np.all(np.isfinite(result['node_energy'])) and np.all(np.isfinite(result['err'])))
            entry['completes'] = bool(entry['finite'] and np.all(result['err'] < 1e-8))
            if warns:
                entry['warnings'] = warns
        except Exception as exc:  # recorded, not raised: "it used to fail this way" is behaviour too
            entry['status'] = 'raised'
            entry['exception'] = type(exc).__name__
            entry['message'] = str(exc)
        cells.append(entry)
    return arrays, cells


def utils_fixture(rng):
    arrays = {}
    x = rng.normal(size=30)
    arrays['x'] = x
    arrays['normalize_state'] = normalize_state(x)
    arrays['normalize_state_bool'] = normalize_state(x > 0)
    for rank in (True, False):
        for add_constant in (True, False):
            arrays['normalize_weights_rank{0:d}_const{1:d}'.format(rank, add_constant)] = \
                normalize_weights(x, rank=rank, add_constant=add_constant)
    null = rng.normal(size=1000)
    observed = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
    arrays['null'], arrays['observed'] = null, observed
    for version in ('standard', 'reverse', 'smallest'):
        for use_abs in (False, True):
            arrays['get_null_p_{0}_abs{1:d}'.format(version, use_abs)] = np.array(
                [get_null_p(o, null, version=version, abs=use_abs) for o in observed])
    p_vals = np.sort(rng.random(40) ** 3)
    arrays['p_vals'] = p_vals
    arrays['get_fdr_p'] = np.asarray(get_fdr_p(p_vals))
    arrays['get_fdr_p_alpha01'] = np.asarray(get_fdr_p(p_vals, alpha=0.1))
    p_strings = [0.0, 1e-30, 0.001, 0.049, 0.05, 0.5]
    arrays['p_val_inputs'] = np.array(p_strings)
    arrays['get_p_val_string'] = np.array([get_p_val_string(p) for p in p_strings])
    labels = list(rng.choice(['Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default'], size=60))
    arrays['convert_states_str2int_input'] = np.array(labels)
    states, state_labels = convert_states_str2int(labels)
    arrays['convert_states_str2int_states'] = states
    arrays['convert_states_str2int_labels'] = np.array([str(s) for s in state_labels])
    M = rng.normal(size=(8, 8))
    M = (M + M.T) / 2 / 8
    arrays['expm_input'], arrays['expm'] = M, np.real_if_close(expm(M))
    return arrays


def content_sha256(path):
    """Hash of the arrays in an .npz: names, dtypes, shapes and values. The file bytes themselves are not
    reproducible, because the zip entries carry the time of writing."""
    h = hashlib.sha256()
    with np.load(path) as data:
        for key in sorted(data.files):
            a = np.ascontiguousarray(data[key])
            h.update('{0}|{1}|{2}|'.format(key, a.dtype.str, a.shape).encode())
            h.update(a.tobytes())
    return h.hexdigest()


def build(design, seed=SEED):
    """Compute every fixture in memory, in the order the random numbers are drawn.

    Returns ({connectome name: (connectome, arrays, cells)}, utils arrays). test_regression.py calls this
    with the manifest's design to recompute the fixtures and compare them against the stored ones.
    """
    rng = np.random.default_rng(seed)
    results = {}
    for name, conn in load_connectomes(rng).items():
        arrays, cells = connectome_fixture(name, conn, design, rng)
        results[name] = (conn, arrays, cells)
    return results, utils_fixture(rng)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    design = pairwise_design(FACTORS, seed=SEED)
    results, utils = build(design)
    manifest = {
        'description': 'nctpy regression fixtures: a drift alarm, see make_regression_fixtures.py',
        'generated': date.today().isoformat(),
        'environment': {
            'nctpy_commit': subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True,
                                           cwd=REPO).stdout.strip(),
            'nctpy_worktree_clean': subprocess.run(['git', 'status', '--porcelain', '--', 'src/nctpy', 'src/null_models'],
                                                   capture_output=True, text=True, cwd=REPO).stdout.strip() == '',
            'python': platform.python_version(), 'numpy': np.__version__, 'scipy': sp.__version__,
            'platform': platform.platform(),
        },
        'seed': SEED,
        'factors': FACTORS,
        'T_values': T_VALUES,
        'design': design,
        'n_trajectory_points': N_TRAJECTORY_POINTS,
        'max_full_matrix': MAX_FULL_MATRIX,
        'connectomes': {},
        'files': {},
    }
    for name, (conn, arrays, cells) in results.items():
        np.savez_compressed(OUT / (name + '.npz'), **arrays)
        manifest['connectomes'][name] = {'source': conn['source'], 'n_nodes': len(conn['A']),
                                         'directed': bool(not np.allclose(conn['A'], conn['A'].T)),
                                         'self_connections': bool(np.any(np.diag(conn['A']) != 0)),
                                         'cells': cells}
        status = [c['status'] for c in cells]
        print('{0:13s} {1:3d} nodes: {2} cells returned, {3} raised, {4} incomplete'.format(
            name, len(conn['A']), status.count('returned'), status.count('raised'),
            sum(c['status'] == 'returned' and not c['completes'] for c in cells)))
    np.savez_compressed(OUT / 'utils.npz', **utils)
    for path in sorted(OUT.glob('*.npz')):
        manifest['files'][path.name] = content_sha256(path)
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1, default=str) + '\n')
    print('design: {0} cells covering all pairs of {1} factors; wrote {2}'.format(len(design), len(FACTORS), OUT))


if __name__ == '__main__':
    sys.exit(main())
