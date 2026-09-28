"""Paper-code test: the code printed in the Nature Protocols paper and its SI, run as printed.

Parkes, Kim et al., Nature Protocols 19:3721-3749 (2024), doi:10.1038/s41596-024-01023-w.

Each block below is the text exactly as typeset. Two kinds of edit are applied before running it,
each listed next to the block and each required to match exactly as many times as stated:

- REPAIRS undo typesetting damage in the main text: lines the typesetter broke (which makes the
  printed text a SyntaxError or NameError) and indentation it dropped (an IndentationError). The
  SI was typeset from real listings and needs none.
- SUBSTITUTIONS point the code at the data shipped in data/ (the paper's '/path/to/data' and
  placeholder file names) and cut the SI's 5000 null-model permutations to a few.

The checks are: every block runs; outputs have the printed shape and type; the paper's own error
threshold holds; and printed values are reproduced to within 1% (decision D5). Error terms are
never matched, only thresholded: they are rounding noise.

The main-text chain (imports, Procedure 1 steps 1-6, Box 1, Procedure 2) runs once and is shared.
Each SI section starts from a copy of its end state, since the SI presents variations on it.
"""
import contextlib
import functools
import io
import os
import re
import unittest
import warnings
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')  # the paper calls plt.show()

import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
DATADIR = REPO / 'data'
ANNOT_DIR = DATADIR / 'schaefer_parc' / 'fsaverage5'
REL = 0.01  # D5: printed values are sanity checks at 1% relative
THR = 1e-8  # the paper's threshold for both error terms


# The paper's data are not in the repository (data/ is git-ignored; the PNC data are not freely
# shareable), so this test runs only where a local copy of data/ exists and skips everywhere else.
REQUIRED_DATA = [
    DATADIR / 'pnc_schaefer200_Am.npy',
    DATADIR / 'pnc_schaefer200_system_labels.txt',
    DATADIR / 'pnc_schaefer200_rsts.npy',
    DATADIR / 'pnc_schaefer200_centroids.csv',
    DATADIR / 'schaefer200_cyto.npy',
    DATADIR / 'schaefer200_micro.npy',
    ANNOT_DIR / 'lh.Schaefer2018_200Parcels_7Networks_order.annot',
    ANNOT_DIR / 'rh.Schaefer2018_200Parcels_7Networks_order.annot',
]


def setUpModule():
    absent = [str(p.relative_to(REPO)) for p in REQUIRED_DATA if not p.exists()]
    if absent:
        raise unittest.SkipTest('the paper\'s data are not available locally: ' + ', '.join(absent))
    missing = []
    for name in ('pandas', 'sklearn', 'seaborn', 'nilearn', 'nibabel', 'matplotlib'):
        try:
            __import__(name)
        except ImportError:
            missing.append(name)
    if missing:
        raise unittest.SkipTest('the paper imports packages that are not installed: ' + ', '.join(missing))


def edit(text, repairs=(), substitutions=()):
    """Apply repairs and substitutions, each (old, new, expected_count)."""
    code = text
    for old, new, count in list(repairs) + list(substitutions):
        found = code.count(old)
        if found != count:
            raise AssertionError('edit {0!r} matched {1} times, expected {2}'.format(old, found, count))
        code = code.replace(old, new)
    return code


def run(label, code, ns):
    """Execute a block in namespace ns; return what it printed."""
    import matplotlib.pyplot as plt
    out = io.StringIO()
    with warnings.catch_warnings(), contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        warnings.simplefilter('ignore')
        exec(compile(code, label, 'exec'), ns)
    plt.close('all')
    return out.getvalue()


def numbers(line):
    return [float(v) for v in re.findall(r'-?\d+\.?\d*(?:[eE][-+]?\d+)?', line)]


DATA_SUBSTITUTIONS = (
    ("datadir = '/path/to/data'", 'datadir = {0!r}'.format(str(DATADIR)), 1),
    ("adjacency_file = 'structural_connectome.npy'", "adjacency_file = 'pnc_schaefer200_Am.npy'", 1),
)

# ---------------------------------------------------------------------------------------------
# Main text
# ---------------------------------------------------------------------------------------------

MAIN_IMPORTS = r'''# import
import os
import numpy as np
import pandas as pd
import scipy as sp
from scipy import stats
from scipy.spatial import distance
from sklearn.cluster import KMeans
from tqdm import tqdm
# import plotting libraries
import matplotlib.pyplot as plt
import seaborn as sns
from nilearn import datasets
from nilearn import plotting
# import nctpy functions
from nctpy.energies import integrate_u, get_control_inputs
from nctpy.pipelines import ComputeControlEnergy,
ComputeOptimizedControlEnergy
from nctpy.metrics import ave_control
from nctpy.utils import matrix_normalization, convert_states_str2int,
normalize_state, normalize_weights, get_null_p, get_fdr_p
from nctpy.plotting import roi_to_vtx, null_plot, surface_plot, add_
module_lines
from null_models.geomsurr import geomsurr
'''
MAIN_IMPORTS_REPAIRS = (  # p. 3737: three lines broken by the typesetter
    ('ComputeControlEnergy,\nComputeOptimizedControlEnergy', 'ComputeControlEnergy, ComputeOptimizedControlEnergy', 1),
    ('convert_states_str2int,\nnormalize_state', 'convert_states_str2int, normalize_state', 1),
    ('add_\nmodule_lines', 'add_module_lines', 1),
)

MAIN_LOAD = r'''# directory where data is stored
datadir = '/path/to/data'
adjacency_file = 'structural_connectome.npy'
# load adjacency matrix
adjacency = np.load(os.path.join(datadir, adjacency_file))
n_nodes = adjacency.shape[0]
print(adjacency.shape)
# check for self-connections
print(np.any(np.diag(adjacency) > 0))
# get density including self connections
density = np.count_nonzero(np.triu(adjacency, k=0)) / (n_nodes**2 / 2)
print(density)
'''

MAIN_STEP1 = r'''# determine time system. Note, delete the line below that is not needed.
system = "discrete"
# or
system = "continuous"
'''

MAIN_STEP2 = r'''# normalize adjacency matrix
adjacency_norm = matrix_normalization(A=adjacency, system=system, c=1)
'''

MAIN_STEP3A = r'''# load node-to-system mapping
system_labels = list(
 np.loadtxt(os.path.join(datadir, "pnc_schaefer200_system_labels.txt"),
dtype=str)
)
print(len(system_labels))
print(system_labels[:20])
'''

MAIN_STEP3B = r'''# use list of system names to create states
states, state_labels = convert_states_str2int(system_labels)
print(type(state_labels), len(state_labels), state_labels)
print(type(states), states.shape, states)
'''

MAIN_STEP3C = r'''# extract initial state
initial_state = states == state_labels.index('Vis')
# extract target state
target_state = states == state_labels.index('Default')
'''

MAIN_STEP3D = r'''# normalize state magnitude
initial_state = normalize_state(initial_state)
target_state = normalize_state(target_state)
'''

MAIN_STEP3E = r'''# specify a uniform full control set: all nodes are control nodes
# and all control nodes are assigned equal control weight
control_set = np.eye(n_nodes)
'''

MAIN_STEP4A = r'''# set parameters
time_horizon = 1 # time horizon (T)
rho = 1 # mixing parameter for state trajectory constraint
trajectory_constraints = np.eye(n_nodes) # nodes in state
trajectory to be constrained
# get the state trajectory, x(t), and the control signals, u(t)
state_trajectory, control_signals, numerical_error = get_control_
inputs(
 A_norm=adjacency_norm,
T=time_horizon,
B=control_set,
x0=initial_state,
xf=target_state,
system=system,
 rho=rho,
 S=trajectory_constraints,
)
'''
MAIN_STEP4A_REPAIRS = (  # p. 3740: a wrapped comment and a split function name
    ('# nodes in state\ntrajectory to be constrained', '# nodes in state trajectory to be constrained', 1),
    ('get_control_\ninputs(', 'get_control_inputs(', 1),
)

MAIN_STEP4IV = r'''# print errors
thr = 1e-8
# the first numerical error corresponds to the inversion error
print(
 "inversion error = {:.2E} (<{:.2E}={:})".
format(numerical_error[0], thr, numerical_error[0] < thr
 )
)
# the second numerical error corresponds to the reconstruction
error
print(
 "reconstruction error = {:.2E} (<{:.2E}={:})".
format(numerical_error[1], thr, numerical_error[1] < thr
 )
)
'''
MAIN_STEP4IV_REPAIRS = (  # p. 3741: a wrapped comment leaves `error` as a bare name
    ('reconstruction\nerror\n', 'reconstruction error\n', 1),
)

MAIN_STEP6A = r'''# integrate control signals to get control energy
node_energy = integrate_u(control_signals)
print(node_energy.shape)
print(np.round(node_energy[:5], 2))
'''

MAIN_STEP6B = r'''# summarize nodal energy
energy = np.sum(node_energy)
print(np.round(energy, 2))
'''

MAIN_BOX1 = r'''f, ax = plt.subplots(3, 2, figsize=(7, 7))
# plot control signals for initial state
ax[0, 0].plot(control_signals[:, initial_state != 0], linewidth=0.75)
ax[0, 0].set_title("A | control signals, x0")
# plot state trajectory for initial state
ax[0, 1].plot(state_trajectory[:, initial_state != 0], linewidth=0.75)
ax[0, 1].set_title("B | neural activity, x0")
# plot control signals for target state
ax[1, 0].plot(control_signals[:, target_state != 0], linewidth=0.75)
ax[1, 0].set_title("C | control signals, xf")
# plot state trajectory for target state
ax[1, 1].plot(state_trajectory[:, target_state != 0], linewidth=0.75)
ax[1, 1].set_title("D | neural activity, xf")
# plot control signals for bystanders
ax[2, 0].plot(
 control_signals[:, np.logical_and(initial_state == 0, target_state == 0)],
 linewidth=0.75,
)
ax[2, 0].set_title("E | control signals, bystanders")
# plot state trajectory for bystanders
ax[2, 1].plot(
 state_trajectory[:, np.logical_and(initial_state == 0, target_state == 0)],
 linewidth=0.75,
)
ax[2, 1].set_title("F | neural activity, bystanders")
for cax in ax.reshape(-1):
cax.set_ylabel("activity")
cax.set_xlabel("time (a.u.)")
cax.set_xticks([0, state_trajectory.shape[0]])
cax.set_xticklabels([0, time_horizon])
f.tight_layout()
plt.show()
'''
MAIN_BOX1_REPAIRS = (  # p. 3743: the loop body lost its indentation
    ('for cax in ax.reshape(-1):\n'
     'cax.set_ylabel("activity")\n'
     'cax.set_xlabel("time (a.u.)")\n'
     'cax.set_xticks([0, state_trajectory.shape[0]])\n'
     'cax.set_xticklabels([0, time_horizon])\n',
     'for cax in ax.reshape(-1):\n'
     '    cax.set_ylabel("activity")\n'
     '    cax.set_xlabel("time (a.u.)")\n'
     '    cax.set_xticks([0, state_trajectory.shape[0]])\n'
     '    cax.set_xticklabels([0, time_horizon])\n', 1),
)

MAIN_PROCEDURE2 = r'''# compute average controllability
average_controllability = ave_control(A_norm=adjacency_norm,
system=system)
'''

# Printed outputs (main text)
PRINTED_FIRST_20_LABELS = ['Vis'] * 14 + ['SomMot'] * 6
PRINTED_STATE_LABELS = ['Cont', 'Default', 'DorsAttn', 'Limbic', 'SalVentAttn', 'SomMot', 'Vis']
# p. 3739 prints states.shape as (200,) but lists 199 values: the right-hemisphere Limbic run shows
# five 3s where the linked labels file has six. Expected here, run-length encoded, with six.
PRINTED_STATES_RUNS = [(6, 14), (5, 16), (2, 13), (4, 11), (3, 6), (0, 13), (1, 27),
                       (6, 15), (5, 19), (2, 13), (4, 11), (3, 6), (0, 17), (1, 19)]
PRINTED_NODE_ENERGY = [21.13, 37.65, 23.55, 21.55, 28.34]
PRINTED_ENERGY = 2604.71


@functools.lru_cache(maxsize=None)
def main_text():
    """Run the main-text chain once; return (namespace, printed output per block)."""
    ns = {'__name__': 'paper'}
    blocks = [
        ('imports', MAIN_IMPORTS, MAIN_IMPORTS_REPAIRS, ()),
        ('load', MAIN_LOAD, (), DATA_SUBSTITUTIONS),
        ('step1', MAIN_STEP1, (), ()),
        ('step2', MAIN_STEP2, (), ()),
        ('step3A', MAIN_STEP3A, (), ()),
        ('step3B', MAIN_STEP3B, (), ()),
        ('step3C', MAIN_STEP3C, (), ()),
        ('step3D', MAIN_STEP3D, (), ()),
        ('step3E', MAIN_STEP3E, (), ()),
        ('step4A', MAIN_STEP4A, MAIN_STEP4A_REPAIRS, ()),
        ('step4iv', MAIN_STEP4IV, MAIN_STEP4IV_REPAIRS, ()),
        ('box1', MAIN_BOX1, MAIN_BOX1_REPAIRS, ()),
        ('step6A', MAIN_STEP6A, (), ()),
        ('step6B', MAIN_STEP6B, (), ()),
        ('procedure2', MAIN_PROCEDURE2, (), ()),
    ]
    printed = {}
    for label, text, repairs, substitutions in blocks:
        printed[label] = run(label, edit(text, repairs, substitutions), ns)
    return ns, printed


class PaperTestCase(unittest.TestCase):
    def assertClose(self, got, expected, rel=REL):
        self.assertLessEqual(abs(got - expected), rel * abs(expected),
                             '{0} is not within {1:.0%} of the printed {2}'.format(got, rel, expected))


class TestMainText(PaperTestCase):
    """Procedure 1 (steps 1-6), Box 1 and Procedure 2, run as printed."""

    @classmethod
    def setUpClass(cls):
        cls.ns, cls.printed = main_text()

    def test_connectome(self):
        lines = self.printed['load'].splitlines()
        self.assertEqual(lines[0], '(200, 200)')
        self.assertEqual(lines[1], 'True')  # self-connections present
        self.assertAlmostEqual(float(lines[2]), 0.9768, places=4)

    def test_system_labels_and_states(self):
        lines = self.printed['step3A'].splitlines()
        self.assertEqual(lines[0], '200')
        self.assertEqual([str(s) for s in self.ns['system_labels'][:20]], PRINTED_FIRST_20_LABELS)
        state_labels = self.ns['state_labels']
        self.assertIsInstance(state_labels, list)
        # values, not repr: numpy >= 2 prints np.str_('Cont') (D11, accepted as cosmetic)
        self.assertEqual([str(s) for s in state_labels], PRINTED_STATE_LABELS)
        states = self.ns['states']
        self.assertIsInstance(states, np.ndarray)
        self.assertEqual(states.shape, (200,))
        expected = np.concatenate([np.full(n, v) for v, n in PRINTED_STATES_RUNS])
        np.testing.assert_array_equal(states, expected)

    def test_states_are_boolean_then_normalized(self):
        ns = self.ns
        self.assertAlmostEqual(np.linalg.norm(ns['initial_state']), 1.0)
        self.assertAlmostEqual(np.linalg.norm(ns['target_state']), 1.0)
        self.assertEqual(int(np.count_nonzero(ns['initial_state'])), 29)  # Vis nodes
        self.assertEqual(int(np.count_nonzero(ns['target_state'])), 46)  # Default nodes

    def test_trajectory_and_signals(self):
        x, u = self.ns['state_trajectory'], self.ns['control_signals']
        self.assertEqual(x.shape, (1001, 200))  # time x nodes, as Box 1 indexes it
        self.assertEqual(u.shape, (1001, 200))
        np.testing.assert_allclose(x[-1], self.ns['target_state'], atol=1e-8)

    def test_numerical_errors(self):
        lines = self.printed['step4iv'].splitlines()
        self.assertEqual(len(lines), 2)
        self.assertTrue(lines[0].startswith('inversion error = '))
        self.assertTrue(lines[0].endswith('(<1.00E-08=True)'))
        self.assertTrue(lines[1].startswith('reconstruction error = '))
        self.assertTrue(lines[1].endswith('(<1.00E-08=True)'))

    def test_energy(self):
        lines = self.printed['step6A'].splitlines()
        self.assertEqual(lines[0], '(200,)')
        for got, expected in zip(numbers(lines[1]), PRINTED_NODE_ENERGY):
            self.assertClose(got, expected)
        self.assertClose(numbers(self.printed['step6B'])[0], PRINTED_ENERGY)

    def test_diagonal_is_not_zeroed(self):
        # the published energy depends on the connectome keeping its self-connections
        self.assertTrue(np.any(np.diag(self.ns['adjacency']) > 0))
        self.assertTrue(np.any(np.diag(self.ns['adjacency_norm']) > -1))

    def test_average_controllability(self):
        ac = self.ns['average_controllability']
        self.assertEqual(ac.shape, (200,))
        self.assertTrue(np.all(np.isfinite(ac)))


# ---------------------------------------------------------------------------------------------
# Supplementary Information
# ---------------------------------------------------------------------------------------------

SI_3A_LOAD = r'''# load resting-state time series
rsfmri_file = 'pnc_schaefer200_rsts.npy'
rsfmri = np.load(os.path.join(datadir, rsfmri_file))

n_trs = rsfmri.shape[0]
n_nodes_rsfmri = rsfmri.shape[1]
n_subs = rsfmri.shape[2]
print('n_trs, {0}; n_nodes, {1}; n_subs, {2}'.format(n_trs, n_nodes_rsfmri, n_subs))

rsfmri_concat = np.zeros((n_trs * n_subs, n_nodes_rsfmri))
print(rsfmri_concat.shape)

for i in np.arange(n_subs):
    # z score and concatenate subject i's time series
    start_idx = i * n_trs
    end_idx = start_idx + n_trs
    rsfmri_concat[start_idx:end_idx, :] = sp.stats.zscore(rsfmri[:, :, i], axis=0)
'''

SI_3A_CLUSTER = r'''# extract 5 clusters of activity
n_clusters = 5
kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(rsfmri_concat)

# extract cluster centers. These represent dominant patterns of recurrent activity over time
centroids = kmeans.cluster_centers_
print(centroids.shape)

# plot centroids on brain surface
lh_annot_file = ("/path/to/schaefer/files/lh.Schaefer2018_200Parcels_7Networks_order.annot")
rh_annot_file = ("/path/to/schaefer/files/rh.Schaefer2018_200Parcels_7Networks_order.annot")
fsaverage = datasets.fetch_surf_fsaverage(mesh="fsaverage5")

for cluster in np.arange(n_clusters):
    f = surface_plot(
        data=centroids[cluster, :],
        lh_annot_file=lh_annot_file,
        rh_annot_file=rh_annot_file,
        fsaverage=fsaverage,
        order="lr",
        cmap="coolwarm",
    )
'''
ANNOT_SUBSTITUTIONS = (('/path/to/schaefer/files/', str(ANNOT_DIR) + '/', 2),)

SI_3A_ENERGY = r'''# extract visual cluster is initial state
initial_state = centroids[1, :]
# extract default mode cluster as target state
target_state = centroids[4, :]

# normalize state magnitude
initial_state = normalize_state(initial_state)
target_state = normalize_state(target_state)

# get the state trajectory and the control signals
state_trajectory, control_signals, numerical_error = get_control_inputs(
    A_norm=adjacency_norm,
    T=time_horizon,
    B=control_set,
    x0=initial_state,
    xf=target_state,
    system=system,
    rho=rho,
    S=trajectory_constraints,
)

# get energy
node_energy = integrate_u(control_signals)
energy = np.sum(node_energy)
'''

SI_3A_SURFACE = r'''timepoints_to_plot = np.arange(0, state_trajectory.shape[0], int(state_trajectory.shape[0] / 5))

for timepoint in timepoints_to_plot:
    f = surface_plot(
        data=state_trajectory[timepoint, :],
        lh_annot_file=lh_annot_file,
        rh_annot_file=rh_annot_file,
        fsaverage=fsaverage,
        order="lr",
        cmap="coolwarm"
    )
'''

SI_3B = r'''# specify a uniform partial control set: some nodes are control nodes
# and all control nodes are assigned equal control weight
bystanders = np.logical_and(
    initial_state == 0, target_state == 0
)  # use initial state and final state to find bystanders. note, this only works for binary states
control_set = np.zeros((n_nodes, n_nodes))   # initialize control nodes matrix
control_set[bystanders, bystanders] = 1   # set bystanders to control nodes
'''

SI_3C_LOAD = r'''# helper func for printing descriptive stats
def print_stats(x):
    print(
        "min={:.2f}; max={:.2f}; mean={:.2f}; std={:.2f}; skew={:.2f}; kurt={:.2f}".format(
            np.min(x),
            np.max(x),
            np.mean(x),
            np.std(x),
            sp.stats.skew(x),
            sp.stats.kurtosis(x),
        )
    )

neuromap_file = 'neuromap.npy'
neuromap = np.load(os.path.join(datadir, neuromap_file))
print(neuromap.shape)
print_stats(neuromap)

control_set = np.zeros((n_nodes, n_nodes)) # initialize B matrix
control_set[np.diag_indices(n_nodes)] = neuromap # set weights using neuromap
'''

SI_3C_SHIFT = r'''# modify neuromap so that its minimum value is 1
neuromap += 1 + np.abs(np.min(neuromap))
print_stats(neuromap)

control_set = np.zeros((n_nodes, n_nodes)) # initialize B matrix
control_set[np.diag_indices(n_nodes)] = neuromap # set weights using neuromap
'''

SI_3C_NORMALIZE = r'''neuromap_file_1 = 'neuromap.npy'
neuromap_1 = np.load(os.path.join(datadir, neuromap_file_1))
neuromap_1_norm = normalize_weights(neuromap_1)

neuromap_file_2 = 'neuromap2.npy'
neuromap_2 = np.load(os.path.join(datadir, neuromap_file_2))
neuromap_2_norm = normalize_weights(neuromap_2)

print_stats(neuromap_1)
print_stats(neuromap_1_norm)
print_stats(neuromap_2)
print_stats(neuromap_2_norm)
'''
# The SI's maps are not named; these shipped maps reproduce every printed statistic exactly, and
# scripts/path_a_control_energy_binary.ipynb uses the same two.
NEUROMAP_SUBSTITUTIONS = (("'neuromap.npy'", "'schaefer200_cyto.npy'", 1),)
NEUROMAP2_SUBSTITUTIONS = (("'neuromap2.npy'", "'schaefer200_micro.npy'", 1),)

SI_3D_PERTURB = r'''# container for perturbed energies
energy_perturbed = np.zeros(n_nodes)

for node in tqdm(np.arange(n_nodes)):
    # start with a uniform full control set
    control_set = np.eye(n_nodes)

    # add arbitrary amount of additional control to node
    control_set[node, node] += 0.1

    # get perturbed control signals (u_p)
    _, control_signals, _ = get_control_inputs(
        A_norm=adjacency_norm,
        T=time_horizon,
        B=control_set,
        x0=initial_state,
        xf=target_state,
        system=system,
        rho=rho,
        S=trajectory_constraints,
    )

    # integrate control signals to get control energy
    node_energy = integrate_u(control_signals)

    # summarize nodal energy
    energy_perturbed[node] = np.sum(node_energy)

# check if perturbed energy is lower than original energy. Should print True
print(np.all(energy_perturbed < energy))

# calculate energy delta. these values will all be negative,
# indicating reduced energy compared to control_set=np.eye(n_nodes)
energy_delta = energy_perturbed - energy
'''

SI_3D_OPTIMIZE = r'''# re-compute energy using energy deltas as weights
# we do this by taking a single step down the
# gradient created by the energy deltas
learning_rate = 0.01  # set a learning rate for gradient descent
control_set_initial = np.eye(n_nodes)  # initial control set
control_set_optimized = np.zeros((n_nodes, n_nodes))  # initialize container for optimized weights
# step down gradient
control_set_optimized[np.diag_indices(n_nodes)] = (
        control_set_initial.diagonal() -
        (energy_delta * learning_rate)
)
# normalize
control_set_optimized = (
        control_set_optimized /
        sp.linalg.norm(control_set_optimized) *
        sp.linalg.norm(control_set_initial)
)
# normalization ensures that the optimized weights have the same
# norm as control_set_optimized=np.eye(n_nodes)

# get optimized control signals
_, control_signals, _ = get_control_inputs(
    A_norm=adjacency_norm,
    T=time_horizon,
    B=control_set_optimized,
    x0=initial_state,
    xf=target_state,
    system=system,
    rho=rho,
    S=trajectory_constraints,
)

# integrate control signals to get control energy
node_energy_optimized = integrate_u(control_signals)
print(np.round(node_energy_optimized[:5], 2))

# summarize nodal energy
energy_optimized = np.sum(node_energy_optimized)
print(np.round(energy_optimized, 2))
'''

SI_3D_CLASS = r'''control_task = dict() # initialize dict
control_task["x0"] = initial_state # store initial state
control_task["xf"] = target_state  # store target state
control_task["S"] = trajectory_constraints  # store state trajectory constraints
control_task["rho"] = rho  # store rho
compute_opt_control_energy = ComputeOptimizedControlEnergy(
    A=adjacency,
    control_task=control_task,
    system=system,
    c=1,
    T=time_horizon,
    n_steps=2,
    lr=learning_rate,
)
compute_opt_control_energy.run()
'''

SI_D_TASKS = r'''# initialize list of control tasks
control_tasks = []

# define control set using a uniform full control set
# note, here we use the same control set for all control tasks
control_set = np.eye(n_nodes)

# define state trajectory constraints
# note, here we constrain the full state trajectory equally for all control tasks
trajectory_constraints = np.eye(n_nodes)

# define mixing parameter
# note, here we use the same rho for all control tasks
rho = 1

# assemble control tasks
n_states = len(state_labels)
for initial_idx in np.arange(n_states):
    initial_state = normalize_state(states == initial_idx) # initial state
    for target_idx in np.arange(n_states):
        target_state = normalize_state(states == target_idx) # target state

        control_task = dict() # initialize dict
        control_task["x0"] = initial_state # store initial state
        control_task["xf"] = target_state # store target state
        control_task["B"] = control_set # store control set
        control_task["S"] = trajectory_constraints # store state trajectory constraints
        control_task["rho"] = rho # store rho
        control_tasks.append(control_task)
'''

SI_D_RUN = r'''# compute control energy across all control tasks
compute_control_energy = ComputeControlEnergy(
    A=adjacency, control_tasks=control_tasks, system=system, c=1, T=time_horizon
)
compute_control_energy.run()
'''

SI_D_PLOT = r'''# reshape energy into matrix
energy_matrix = np.reshape(compute_control_energy.E, (n_states, n_states))

# subtract lower triangle from upper to examine energy asymmetries
energy_matrix_delta = energy_matrix.transpose() - energy_matrix

f, ax = plt.subplots(1, 3, figsize=(7, 4))

# plot energy matrix
sns.heatmap(
    energy_matrix,
    ax=ax[0],
    square=True,
    linewidth=0.5,
    cbar_kws={"label": "energy", "shrink": 0.25},
)

# plot without self-transitions
# setup mask to exclude persistence energy (i.e., transitions where i==j)
mask = np.zeros_like(energy_matrix)
mask[np.eye(n_states) == 1] = True
sns.heatmap(
    energy_matrix,
    ax=ax[1],
    square=True,
    linewidth=0.5,
    cbar_kws={"label": "energy", "shrink": 0.25},
    mask=mask,
)

# plot energy asymmetries
mask = np.triu(np.ones_like(energy_matrix, dtype=bool))
sns.heatmap(
    energy_matrix_delta,
    ax=ax[2],
    square=True,
    linewidth=0.5,
    cbar_kws={"label": "energy (delta)", "shrink": 0.25},
    mask=mask,
    cmap="RdBu_r",
    center=0,
)

for cax in ax:
    cax.set_ylabel("initial state (x0)")
    cax.set_xlabel("target state (xf)")
    cax.set_yticklabels(state_labels, rotation=0, size=6)
    cax.set_xticklabels(state_labels, rotation=90, size=6)
f.tight_layout()
plt.show()
'''

SI_E_CENTROIDS = r'''# null networks
centroids = pd.read_csv(
    os.path.join(datadir, "pnc_schaefer200_centroids.csv")
)  # load coordinates of nodes
centroids.set_index("node_names", inplace=True)
print(centroids.head())
'''

SI_E_DISTANCE = r'''distance_matrix = distance.pdist(centroids, 'euclidean')  # get euclidean distances between nodes
distance_matrix = distance.squareform(distance_matrix)  # reshape to square matrix
'''

SI_E_ENERGY_NULLS = r'''# extract initial state
initial_state = states == state_labels.index("Vis")  # 'Vis' or 'SalVentAttn'
initial_state = normalize_state(initial_state)  # normalize

# extract target state
target_state = states == state_labels.index("Default")
target_state = normalize_state(target_state)  # normalize

# compute true control energy
_, control_signals, _ = get_control_inputs(
    A_norm=adjacency_norm,
    T=time_horizon,
    B=control_set,
    x0=initial_state,
    xf=target_state,
    system=system,
    rho=rho,
    S=trajectory_constraints,
)
node_energy = integrate_u(control_signals)  # integrate control signals
energy = np.sum(node_energy)  # get energy

# run permutation
n_perms = 5000  # number of permutations

# containers for null distributions
energy_null_sp = np.zeros(n_perms)
energy_null_ssp = np.zeros(n_perms)

for perm in tqdm(np.arange(n_perms)):
    # rewire adjacency matrix using geomsurr
    _, Wsp, Wssp = geomsurr(W=adjacency, D=distance_matrix, seed=perm)
    # Wsp is the adjacency matrix rewired while preserving spatial embedding and the strength distribution
    # Wssp is the adjacency matrix rewired while preserving spatial embedding and the strength sequence
    # this python implementation is included with our toolbox, but if you use these nulls
    # in your own work, please cite:
    #       Roberts et al. NeuroImage (2016), doi:10.1016/j.neuroimage.2015.09.009

    # compute control energy for Wsp
    Wsp = matrix_normalization(A=Wsp, system=system)
    _, control_signals, _ = get_control_inputs(
        A_norm=Wsp,
        T=time_horizon,
        B=control_set,
        x0=initial_state,
        xf=target_state,
        system=system,
        rho=rho,
        S=trajectory_constraints,
    )
    node_energy = integrate_u(control_signals)
    energy_null_sp[perm] = np.sum(node_energy)

    # compute control energy for Wssp
    Wssp = matrix_normalization(A=Wssp, system=system)
    _, control_signals, _ = get_control_inputs(
        A_norm=Wssp,
        T=time_horizon,
        B=control_set,
        x0=initial_state,
        xf=target_state,
        system=system,
        rho=rho,
        S=trajectory_constraints,
    )
    node_energy = integrate_u(control_signals)
    energy_null_ssp[perm] = np.sum(node_energy)

# plot
f, ax = plt.subplots(1, 2, figsize=(7, 3))
null_plot(
    observed=energy,
    null=energy_null_sp,
    xlabel="strength-preserving",
    ax=ax[0],
)
null_plot(
    observed=energy,
    null=energy_null_ssp,
    xlabel="sequence-preserving",
    ax=ax[1],
)
f.tight_layout()
plt.show()
'''

SI_E_AVE_CTRB_NULLS = r'''# run permutation
n_perms = 5000 # number of permutations

# container for null distribution
ave_ctrb_null = np.zeros((n_perms, n_nodes))

for perm in tqdm(np.arange(n_perms)):
    # rewire adjacency matrix using geomsurr
    _, _, Wssp = geomsurr(W=adjacency, D=distance_matrix, seed=perm)

    # compute average controllability
    Wssp = matrix_normalization(A=Wssp, system=system)
    ave_ctrb_null[perm, :] = ave_control(A_norm=Wssp, system=system)
'''

SI_E_P_VALUES = r'''# calculate p-values
p_vals_ssp = np.zeros(n_nodes)

for node in tqdm(np.arange(n_nodes)):
    # version='standard' will calculate the number of times the null is larger than the observed value
    # version='reverse' will calculate the number of times the null is smaller than the observed value
    p_vals_ssp[node] = get_null_p(
        x=average_controllability[node], null=ave_ctrb_null[:, node], version="standard"
    )

p_vals_ssp = get_fdr_p(p_vals=p_vals_ssp)            # correct p values for multiple comparisons
'''
N_PERMS = 3
PERM_SUBSTITUTIONS = (('n_perms = 5000', 'n_perms = {0}'.format(N_PERMS), 1),)


class SITestCase(PaperTestCase):
    """Each SI section starts from a copy of the main text's end state."""

    @classmethod
    def setUpClass(cls):
        ns, _ = main_text()
        cls.ns = dict(ns)
        cls.printed = {}

    @classmethod
    def run_block(cls, label, text, substitutions=()):
        cls.printed[label] = run(label, edit(text, (), substitutions), cls.ns)
        return cls.printed[label]


class TestSINonBinaryStates(SITestCase):
    """SI C1 (3a): brain states from k-means clusters of rs-fMRI."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.run_block('3a_load', SI_3A_LOAD)
        cls.run_block('3a_cluster', SI_3A_CLUSTER, ANNOT_SUBSTITUTIONS)
        cls.run_block('3a_energy', SI_3A_ENERGY)
        cls.run_block('3a_surface', SI_3A_SURFACE)

    def test_printed_outputs(self):
        self.assertEqual(self.printed['3a_load'].splitlines(), ['n_trs, 120; n_nodes, 200; n_subs, 253', '(30360, 200)'])
        self.assertEqual(self.printed['3a_cluster'].splitlines(), ['(5, 200)'])

    def test_non_binary_energy(self):
        # cluster order depends on the scikit-learn version, so only the outputs' form is checked
        self.assertEqual(self.ns['state_trajectory'].shape, (1001, 200))
        self.assertTrue(np.isfinite(self.ns['energy']) and self.ns['energy'] > 0)


class TestSIControlSets(SITestCase):
    """SI C2 (3b, 3c): partial and variable control sets."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.run_block('3b', SI_3B)
        cls.control_set_3b = cls.ns['control_set']
        cls.run_block('3c_load', SI_3C_LOAD, NEUROMAP_SUBSTITUTIONS)
        cls.run_block('3c_shift', SI_3C_SHIFT)
        cls.control_set_3c = cls.ns['control_set']
        cls.run_block('3c_normalize', SI_3C_NORMALIZE, NEUROMAP_SUBSTITUTIONS + NEUROMAP2_SUBSTITUTIONS)

    def test_bystander_control_set(self):
        B = self.control_set_3b
        self.assertEqual(int(np.trace(B)), 200 - 29 - 46)

    def test_printed_neuromap_stats(self):
        self.assertEqual(self.printed['3c_load'].splitlines(), [
            '(200,)',
            'min=-0.41; max=0.42; mean=-0.04; std=0.20; skew=0.25; kurt=-0.54'])
        self.assertEqual(self.printed['3c_shift'].splitlines(), [
            'min=1.00; max=1.84; mean=1.38; std=0.20; skew=0.25; kurt=-0.54'])
        self.assertEqual(self.printed['3c_normalize'].splitlines(), [
            'min=-0.41; max=0.42; mean=-0.04; std=0.20; skew=0.25; kurt=-0.54',
            'min=1.00; max=2.00; mean=1.50; std=0.29; skew=-0.00; kurt=-1.20',
            'min=-0.10; max=0.16; mean=-0.00; std=0.07; skew=0.59; kurt=-0.78',
            'min=1.00; max=2.00; mean=1.50; std=0.29; skew=-0.00; kurt=-1.20'])


class TestSIOptimizedControlSet(SITestCase):
    """SI C2 (3d): one gradient step on the control weights, by hand and with the class."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.run_block('3d_perturb', SI_3D_PERTURB)
        cls.run_block('3d_optimize', SI_3D_OPTIMIZE)
        cls.run_block('3d_class', SI_3D_CLASS)

    def test_perturbation_lowers_energy(self):
        self.assertEqual(self.printed['3d_perturb'].splitlines(), ['True'])

    def test_optimized_energy(self):
        lines = self.printed['3d_optimize'].splitlines()
        for got, expected in zip(numbers(lines[0]), [20.56, 34.94, 22.77, 20.92, 26.96]):
            self.assertClose(got, expected)
        self.assertClose(numbers(lines[1])[0], 2429.68)

    def test_class_reproduces_the_manual_step(self):
        # the SI says the class "wraps the above optimization steps": its first step is the same
        compute = self.ns['compute_opt_control_energy']
        self.assertEqual(compute.E_opt.shape, (2,))
        self.assertEqual(compute.B_opt.shape, (2, 200))
        np.testing.assert_allclose(compute.E_opt[0], self.ns['energy_optimized'], rtol=1e-9)


class TestSIWrapper(SITestCase):
    """SI D: ComputeControlEnergy over all 49 system-to-system transitions."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.run_block('D_tasks', SI_D_TASKS)
        cls.run_block('D_run', SI_D_RUN)
        cls.run_block('D_plot', SI_D_PLOT)
        labels = [str(s) for s in cls.ns['state_labels']]
        cls.energy = {(labels[i], labels[j]): cls.ns['energy_matrix'][i, j]
                      for i in range(len(labels)) for j in range(len(labels))}

    def test_energy_matrix(self):
        self.assertEqual(self.ns['compute_control_energy'].E.shape, (49,))
        self.assertEqual(self.ns['energy_matrix'].shape, (7, 7))

    def test_matches_single_task(self):
        self.assertClose(self.energy[('Vis', 'Default')], PRINTED_ENERGY)

    def test_figure_values(self):
        # Fig. S14 caption (persistence) and the energies quoted with Fig. S15
        self.assertClose(self.energy[('Default', 'Default')], 571)
        self.assertClose(self.energy[('Vis', 'Default')], 2605)
        self.assertClose(self.energy[('Default', 'Vis')], 1947)
        self.assertClose(self.energy[('SalVentAttn', 'Default')], 2218)
        self.assertClose(self.energy[('Default', 'SalVentAttn')], 2601)


class TestSINullModels(SITestCase):
    """SI E: spatially embedded nulls for energy and average controllability (3 permutations)."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.run_block('E_centroids', SI_E_CENTROIDS)
        cls.run_block('E_distance', SI_E_DISTANCE)
        cls.run_block('E_energy_nulls', SI_E_ENERGY_NULLS, PERM_SUBSTITUTIONS)
        cls.run_block('E_ave_ctrb_nulls', SI_E_AVE_CTRB_NULLS, PERM_SUBSTITUTIONS)
        cls.run_block('E_p_values', SI_E_P_VALUES)

    def test_centroids(self):
        head = self.ns['centroids'].head()
        self.assertEqual(list(head.index), ['LH_Vis_1', 'LH_Vis_2', 'LH_Vis_3', 'LH_Vis_4', 'LH_Vis_5'])
        self.assertEqual(list(head.columns), ['vox_x', 'vox_y', 'vox_z'])
        np.testing.assert_array_equal(head.values, [[121, 149, 69], [123, 174, 65], [143, 166, 70],
                                                    [107, 164, 74], [124, 192, 66]])

    def test_energy_nulls(self):
        self.assertClose(self.ns['energy'], PRINTED_ENERGY)
        for name in ('energy_null_sp', 'energy_null_ssp'):
            null = self.ns[name]
            self.assertEqual(null.shape, (N_PERMS,))
            self.assertTrue(np.all(np.isfinite(null)) and np.all(null > 0))

    def test_null_models_leave_adjacency_intact(self):
        # geomsurr is called on `adjacency` itself; it must not alter it (fixed in 1.0.3)
        self.assertTrue(np.any(np.diag(self.ns['adjacency']) > 0))

    def test_ave_ctrb_nulls_and_p_values(self):
        self.assertEqual(self.ns['ave_ctrb_null'].shape, (N_PERMS, 200))
        p = np.asarray(self.ns['p_vals_ssp'])
        self.assertEqual(p.shape, (200,))
        self.assertTrue(np.all((p >= 0) & (p <= 1)))


class TestSIFigureCaptions(PaperTestCase):
    """Values in SI figure captions, computed with the printed functions and the main-text setup."""

    @classmethod
    def setUpClass(cls):
        cls.ns, _ = main_text()

    def solve(self, B=None, T=1, A_norm=None):
        ns = self.ns
        _, u, err = ns['get_control_inputs'](
            A_norm=ns['adjacency_norm'] if A_norm is None else A_norm, T=T,
            B=ns['control_set'] if B is None else B, x0=ns['initial_state'], xf=ns['target_state'],
            system=ns['system'], rho=ns['rho'], S=ns['trajectory_constraints'])
        return np.sum(ns['integrate_u'](u)), err

    def test_time_horizons(self):  # Figs S7-S9
        for T, printed in ((2, 1904), (5, 1797), (10, 1801)):
            with self.subTest(T=T):
                self.assertClose(self.solve(T=T)[0], printed)

    def test_partial_control_sets(self):
        x0, xf = self.ns['initial_state'] != 0, self.ns['target_state'] != 0
        bystanders = ~x0 & ~xf
        self.assertClose(self.solve(B=np.diag(bystanders.astype(float)))[0], 6.64e9)  # Fig. S3
        B = np.diag(np.where(xf, 1.0, 1e-5))
        self.assertClose(self.solve(B=B)[0], 2.31e11)  # Fig. S6

    def test_incomplete_transitions_return(self):  # Figs S4, S5: returned, not raised
        for mask in (self.ns['initial_state'] != 0, self.ns['target_state'] != 0):
            with self.subTest(n_control_nodes=int(mask.sum())):
                energy, err = self.solve(B=np.diag(mask.astype(float)))
                self.assertTrue(np.isfinite(energy))
                self.assertGreater(err[1], THR)  # the transition does not complete

    def test_annotation_map_control_set(self):  # Fig. S11: the shifted map from SI 3c
        neuromap = np.load(DATADIR / 'schaefer200_cyto.npy')
        neuromap = neuromap + 1 + np.abs(np.min(neuromap))
        self.assertClose(self.solve(B=np.diag(neuromap))[0], 1455)

    def test_tolerance_detects_a_zeroed_diagonal(self):
        # the 1% tolerance must catch the diagonal being zeroed (2604.71 -> 2638.03, +1.3%)
        adjacency = self.ns['adjacency'].copy()
        np.fill_diagonal(adjacency, 0)
        A_norm = self.ns['matrix_normalization'](A=adjacency, system=self.ns['system'], c=1)
        energy = self.solve(A_norm=A_norm)[0]
        self.assertClose(energy, 2638.03, rel=1e-4)
        self.assertGreater(abs(energy - PRINTED_ENERGY), REL * PRINTED_ENERGY)


if __name__ == '__main__':
    unittest.main()
