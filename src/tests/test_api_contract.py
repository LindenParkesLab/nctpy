"""Contract test: the public API is frozen within 1.x.

The Nature Protocols paper and its Supplementary Information print code that must keep running
unmodified. This module checks three things:

1. every public symbol is importable from its documented path (including the top-level
   ``null_models`` package the paper imports from);
2. every public signature matches fixtures/api_contract.json (parameter names, order, kind, and
   literal defaults; non-literal defaults are checked for presence only);
3. the return arity, types and orientation that the printed code relies on, on a small synthetic
   system.

To change the API deliberately (within 1.x: only a new keyword argument whose default reproduces
existing behaviour), regenerate the JSON with make_api_contract.py.
"""
import contextlib
import importlib
import inspect
import io
import json
import os
import subprocess
import sys
import textwrap
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from scipy.spatial.distance import pdist, squareform

FIXTURES = Path(__file__).resolve().parent / 'fixtures'
CONTRACT = json.loads((FIXTURES / 'api_contract.json').read_text())['symbols']
MODULES = sorted({s['module'] for s in CONTRACT})
LITERALS = (type(None), bool, int, float, str)


def import_module(name):
    """Import a contract module without touching the network.

    nctpy.plotting evaluates datasets.fetch_surf_fsaverage() as a default argument when it is
    imported, so the fetch is patched out. Printed code always passes fsaverage explicitly.
    """
    if name != 'nctpy.plotting':
        return importlib.import_module(name)
    try:
        import nilearn.datasets  # noqa: F401
    except ImportError as exc:
        raise unittest.SkipTest('plotting dependencies not installed: {0}'.format(exc))
    with mock.patch('nilearn.datasets.fetch_surf_fsaverage', return_value={}):
        return importlib.import_module(name)


def describe(params):
    """Parameters as comparable dicts; non-literal defaults reduced to a presence marker."""
    out = []
    for p in params:
        entry = {'name': p.name, 'kind': p.kind.name}
        if p.default is not inspect.Parameter.empty:
            if isinstance(p.default, LITERALS):
                entry['default'] = {'value': p.default, 'type': type(p.default).__name__}
            else:
                entry['default'] = {'expr': True}
        out.append(entry)
    return out


def expected(params):
    out = []
    for p in params:
        entry = {'name': p['name'], 'kind': p['kind']}
        if 'default' in p:
            if 'value' in p['default']:
                value = p['default']['value']
                entry['default'] = {'value': value, 'type': type(value).__name__}
            else:
                entry['default'] = {'expr': True}
        out.append(entry)
    return out


class TestImportPaths(unittest.TestCase):
    def test_every_symbol_importable(self):
        for s in CONTRACT:
            with self.subTest(symbol='{0}.{1}'.format(s['module'], s['name'])):
                try:
                    module = import_module(s['module'])
                except unittest.SkipTest:
                    continue
                self.assertTrue(hasattr(module, s['name']))

    def test_no_unrecorded_public_symbols(self):
        # anything public must be frozen too; a new function means regenerating the JSON
        for name in MODULES:
            with self.subTest(module=name):
                try:
                    module = import_module(name)
                except unittest.SkipTest:
                    continue
                public = {n for n, obj in vars(module).items()
                          if not n.startswith('_') and (inspect.isfunction(obj) or inspect.isclass(obj))
                          and obj.__module__ == name}
                recorded = {s['name'] for s in CONTRACT if s['module'] == name}
                self.assertEqual(public - recorded, set(),
                                 'public symbols missing from api_contract.json; run make_api_contract.py')

    def test_null_models_under_nctpy(self):
        # both paths ship (Roadmap 2.7): the top-level one the paper prints, and nctpy.null_models, which
        # re-exports the same objects, so the contract's signatures cover both
        import null_models.geomsurr as top_level
        from nctpy.null_models import geomsurr as within_nctpy
        public = {s['name'] for s in CONTRACT if s['module'] == 'null_models.geomsurr'}
        self.assertEqual(set(within_nctpy.__all__), public)
        for name in public:
            with self.subTest(name=name):
                self.assertIs(getattr(within_nctpy, name), getattr(top_level, name))
        from nctpy.null_models.geomsurr import geomsurr  # the paper's import, under nctpy
        self.assertIs(geomsurr, top_level.geomsurr)


class TestOptionalDependencies(unittest.TestCase):
    """The plotting dependencies are an optional extra (nctpy[plot]); everything else must work without them."""

    BLOCKED = ('matplotlib', 'seaborn', 'nibabel', 'nilearn')

    def test_core_imports_without_plotting_dependencies(self):
        # a fresh interpreter in which the plotting packages cannot be imported, whether installed or not
        script = textwrap.dedent('''
            import sys

            class Block:
                def find_spec(self, name, path=None, target=None):
                    if name.split('.')[0] in {blocked!r}:
                        raise ModuleNotFoundError('blocked for this test: ' + name)
                    return None

            sys.meta_path.insert(0, Block())
            import nctpy.energies, nctpy.metrics, nctpy.pipelines, nctpy.utils, null_models.geomsurr
            import nctpy.null_models.geomsurr
            print('core imported')
            try:
                import nctpy.plotting
            except ImportError as exc:
                print('plotting raised ImportError:', exc)
        ''').format(blocked=set(self.BLOCKED))
        import nctpy
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            filter(None, [os.path.dirname(os.path.dirname(nctpy.__file__)), os.environ.get('PYTHONPATH')])))
        result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('core imported', result.stdout)
        self.assertIn('plotting raised ImportError:', result.stdout)
        self.assertIn("pip install 'nctpy[plot]'", result.stdout)
        self.assertIn("pip install 'nctpy[paper]'", result.stdout)


class TestSignatures(unittest.TestCase):
    def test_signatures_match_contract(self):
        for s in CONTRACT:
            with self.subTest(symbol='{0}.{1}'.format(s['module'], s['name'])):
                try:
                    obj = getattr(import_module(s['module']), s['name'])
                except unittest.SkipTest:
                    continue
                # for a class, inspect.signature gives __init__ without self
                self.assertEqual(describe(inspect.signature(obj).parameters.values()), expected(s['params']))
                for method, params in s.get('methods', {}).items():
                    got = list(inspect.signature(getattr(obj, method)).parameters.values())[1:]  # drop self
                    self.assertEqual(describe(got), expected(params))


class TestReturnContract(unittest.TestCase):
    """Return arity, types and orientation the paper's printed code relies on."""

    def setUp(self):
        from nctpy.utils import matrix_normalization, normalize_state
        rng = np.random.default_rng(0)
        self.n = 10
        A = rng.random((self.n, self.n)) + 0.1
        self.A = (A + A.T) / 2
        np.fill_diagonal(self.A, 0)
        self.D = squareform(pdist(rng.normal(size=(self.n, 3))))
        self.A_c = matrix_normalization(A=self.A, system='continuous', c=1)
        self.A_d = matrix_normalization(A=self.A, system='discrete', c=1)
        mask0 = np.zeros(self.n, dtype=bool)
        mask0[:3] = True
        maskf = np.zeros(self.n, dtype=bool)
        maskf[-3:] = True
        self.x0 = normalize_state(mask0)
        self.xf = normalize_state(maskf)

    def test_matrix_normalization(self):
        from nctpy.utils import matrix_normalization
        for system in ('continuous', 'discrete'):
            A_norm = matrix_normalization(A=self.A, system=system, c=1)
            self.assertIsInstance(A_norm, np.ndarray)
            self.assertEqual(A_norm.shape, (self.n, self.n))
            # the SI relies on the default c=1
            np.testing.assert_array_equal(matrix_normalization(A=self.A, system=system), A_norm)

    def test_get_control_inputs(self):
        from nctpy.energies import get_control_inputs
        for system, A_norm, T in (('continuous', self.A_c, 1), ('discrete', self.A_d, 3)):
            with self.subTest(system=system):
                # the paper's call: all keywords, three return values
                state_trajectory, control_signals, numerical_error = get_control_inputs(
                    A_norm=A_norm, T=T, B=np.eye(self.n), x0=self.x0, xf=self.xf, system=system,
                    rho=1, S=np.eye(self.n))
                # time x nodes (Box 1 indexes [:, mask]; shape[0] is the number of time points)
                self.assertEqual(state_trajectory.shape[1], self.n)
                self.assertEqual(control_signals.shape[1], self.n)
                np.testing.assert_allclose(state_trajectory[0, :], self.x0, atol=1e-8)
                np.testing.assert_allclose(state_trajectory[-1, :], self.xf, atol=1e-6)
                if system == 'continuous':
                    self.assertEqual(state_trajectory.shape[0], 1001)
                    self.assertEqual(control_signals.shape[0], 1001)
                # two error terms, formatted and compared as in the paper
                self.assertEqual(len(numerical_error), 2)
                thr = 1e-8
                for err in numerical_error:
                    '{:.2E} (<{:.2E}={:})'.format(err, thr, err < thr)
                    self.assertTrue(err < thr)

    def test_integrate_u(self):
        from nctpy.energies import get_control_inputs, integrate_u
        _, control_signals, _ = get_control_inputs(A_norm=self.A_c, T=1, B=np.eye(self.n), x0=self.x0,
                                                   xf=self.xf, system='continuous', rho=1, S=np.eye(self.n))
        node_energy = integrate_u(control_signals)
        self.assertIsInstance(node_energy, np.ndarray)
        self.assertEqual(node_energy.shape, (self.n,))
        self.assertTrue(np.isfinite(np.sum(node_energy)))

    def test_metrics(self):
        from nctpy.metrics import ave_control, modal_control
        for system, A_norm in (('continuous', self.A_c), ('discrete', self.A_d)):
            self.assertEqual(ave_control(A_norm=A_norm, system=system).shape, (self.n,))
        self.assertEqual(modal_control(self.A_d).shape, (self.n,))

    def test_convert_states_str2int(self):
        from nctpy.utils import convert_states_str2int
        labels_in = ['Vis', 'Vis', 'SomMot', 'Default', 'Default', 'Cont']
        states, state_labels = convert_states_str2int(labels_in)
        self.assertIsInstance(states, np.ndarray)
        self.assertEqual(states.shape, (len(labels_in),))
        self.assertTrue(np.issubdtype(states.dtype, np.integer))
        self.assertIsInstance(state_labels, list)
        self.assertEqual(len(state_labels), 4)
        self.assertEqual(state_labels, sorted(state_labels))
        for label in state_labels:
            self.assertIsInstance(label, str)
        # the paper extracts states like this
        vis = states == state_labels.index('Vis')
        np.testing.assert_array_equal(vis, [True, True, False, False, False, False])

    def test_normalize_state_accepts_boolean(self):
        from nctpy.utils import normalize_state
        x = normalize_state(np.array([True, False, True, False]))
        self.assertTrue(np.issubdtype(x.dtype, np.floating))
        self.assertAlmostEqual(np.linalg.norm(x), 1.0)

    def test_normalize_weights_default_range(self):
        from nctpy.utils import normalize_weights
        w = normalize_weights(np.random.default_rng(1).normal(size=50))
        self.assertAlmostEqual(np.min(w), 1.0)
        self.assertAlmostEqual(np.max(w), 2.0)

    def test_null_p_and_fdr(self):
        from nctpy.utils import get_null_p, get_fdr_p
        null = np.random.default_rng(2).normal(size=500)
        for version in ('standard', 'reverse'):
            p = get_null_p(x=0.5, null=null, version=version)
            self.assertTrue(0 <= p <= 1)
        p_vals = np.array([0.001, 0.01, 0.2, 0.5])
        self.assertEqual(np.shape(get_fdr_p(p_vals=p_vals)), p_vals.shape)

    def test_geomsurr_returns_three_matrices(self):
        from null_models.geomsurr import geomsurr
        out = geomsurr(W=self.A, D=self.D, seed=0)
        self.assertEqual(len(out), 3)
        _, Wsp, Wssp = out  # as unpacked in the SI
        self.assertEqual(Wsp.shape, (self.n, self.n))
        self.assertEqual(Wssp.shape, (self.n, self.n))

    def test_compute_control_energy(self):
        from nctpy.pipelines import ComputeControlEnergy
        from nctpy.energies import get_control_inputs, integrate_u
        from nctpy.utils import matrix_normalization
        pairs = [(self.x0, self.xf), (self.xf, self.x0), (self.x0, self.x0)]
        tasks = [{'x0': a, 'xf': b, 'B': np.eye(self.n), 'S': np.eye(self.n), 'rho': 1} for a, b in pairs]
        compute = ComputeControlEnergy(A=self.A, control_tasks=tasks, system='continuous', c=1, T=1)
        with contextlib.redirect_stderr(io.StringIO()):
            compute.run()
        self.assertEqual(compute.E.shape, (len(tasks),))
        # one summed energy per task, in task order (the SI reshapes E into a matrix)
        A_norm = matrix_normalization(self.A, system='continuous', c=1)
        for i, (a, b) in enumerate(pairs):
            _, u, _ = get_control_inputs(A_norm=A_norm, T=1, B=np.eye(self.n), x0=a, xf=b,
                                         system='continuous', rho=1, S=np.eye(self.n))
            self.assertAlmostEqual(compute.E[i], np.sum(integrate_u(u)), delta=1e-9 * abs(compute.E[i]))

    def test_compute_optimized_control_energy_without_B(self):
        from nctpy.pipelines import ComputeOptimizedControlEnergy
        # the SI's task dictionary has no 'B'
        task = {'x0': self.x0, 'xf': self.xf, 'S': np.eye(self.n), 'rho': 1}
        compute = ComputeOptimizedControlEnergy(A=self.A, control_task=task, system='continuous', c=1, T=1,
                                                n_steps=2, lr=0.01)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            compute.run()
        self.assertEqual(compute.E_opt.shape, (2,))
        self.assertEqual(compute.B_opt.shape, (2, self.n))


if __name__ == '__main__':
    unittest.main()
