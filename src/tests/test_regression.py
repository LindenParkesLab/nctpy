"""Regression test: recompute the fixtures in fixtures/regression/ and compare (Roadmap 0.5).

The fixtures are a drift alarm. make_regression_fixtures.build() recomputes every stored array in
memory from the same seeded synthetic inputs, and each is compared against the stored copy with
the tolerances decided under D12:

- inputs (states, B and S diagonals, sampled indices, sim_state_eq's U): integer, string and
  boolean inputs exactly equal, float inputs within rtol 1e-12. Some inputs involve arithmetic
  (e.g. rng.random(40) ** 3), which differs in the last bit between numpy versions; a change in
  numpy's seeded random streams would move them by order 1;
- control-task node energies: rtol 1e-6 in cells whose transition completes, 1e-4 in cells whose
  transition does not; nodes with exactly zero energy (uncontrolled nodes) must stay exactly zero;
- trajectories x(t): the same rtol as the cell's energies, relative to the trajectory's scale;
- both error terms: order of magnitude only (within a factor of 10, or both below 1e-12);
- whether a cell raises, returns non-finite values, or completes (both error terms < 1e-8) must
  not change; the completes flag is not compared when the stored reconstruction error lies within
  10x of the 1e-8 threshold, where it flips between environments;
- everything else (normalisation, controllability, Gramians, minimum_energy_fast, sim_state_eq,
  utilities): rtol 1e-10, with an absolute floor of 1e-10 times the array's scale.

A failure means something moved. Explain it: if it is floating-point drift, propose regenerating
(TESTING.md, "Regenerating fixtures", with sign-off); if it is a behavioural change, the change
is wrong. TestManifest also fails if any fixture file is changed without regenerating the manifest.
"""
import json
import unittest
from pathlib import Path

import numpy as np

import make_regression_fixtures as generator

REGRESSION = Path(__file__).resolve().parent / 'fixtures' / 'regression'
MANIFEST = json.loads((REGRESSION / 'manifest.json').read_text())

THRESHOLD = 1e-8  # the paper's threshold for both error terms
RTOL_COMPLETE = 1e-6
RTOL_INCOMPLETE = 1e-4
RTOL_OTHER = 1e-10
RTOL_INPUT = 1e-12
ERROR_FLOOR = 1e-12  # error terms below this are both rounding noise
INPUT_SUFFIXES = ('__sample_index', '__x_index')
INPUT_PREFIXES = ('B_', 'S_', 'states_', 'sim_state_eq_U')
UTILS_INPUTS = ('x', 'null', 'observed', 'p_vals', 'p_val_inputs', 'expm_input', 'convert_states_str2int_input')


def assert_close(test, got, ref, rtol, label):
    got, ref = np.asarray(got), np.asarray(ref)
    test.assertEqual(got.shape, ref.shape, label)
    np.testing.assert_array_equal(np.isnan(got), np.isnan(ref), err_msg=label + ': NaN pattern changed')
    ok = ~np.isnan(ref)
    if not ok.any():
        return
    scale = np.max(np.abs(ref[ok]))
    np.testing.assert_allclose(got[ok], ref[ok], rtol=rtol, atol=rtol * scale, err_msg=label)


def assert_same_input(test, got, ref, label):
    label += ': the seeded inputs were not reproduced; regenerating is needed'
    got, ref = np.asarray(got), np.asarray(ref)
    if ref.dtype.kind in 'iubUS':
        np.testing.assert_array_equal(got, ref, err_msg=label)
    else:
        assert_close(test, got, ref, RTOL_INPUT, label)


def assert_same_magnitude(test, got, ref, label):
    for g, r in zip(np.ravel(got), np.ravel(ref)):
        if np.isnan(r) or np.isnan(g):
            test.assertEqual(np.isnan(g), np.isnan(r), label + ': NaN pattern changed')
        elif max(abs(g), abs(r)) > ERROR_FLOOR:
            test.assertTrue(0.1 <= (abs(g) + 1e-300) / (abs(r) + 1e-300) <= 10,
                            '{0}: error term {1:.2e} is not within 10x of the stored {2:.2e}'.format(label, g, r))


class TestManifest(unittest.TestCase):
    """Guard: fixture files and the manifest must change together."""

    def test_files_listed(self):
        on_disk = sorted(p.name for p in REGRESSION.glob('*.npz'))
        self.assertEqual(on_disk, sorted(MANIFEST['files']))

    def test_content_hashes(self):
        for name, digest in MANIFEST['files'].items():
            with self.subTest(file=name):
                self.assertEqual(generator.content_sha256(REGRESSION / name), digest,
                                 '{0} changed without regenerating manifest.json '
                                 '(run make_regression_fixtures.py; see TESTING.md)'.format(name))


class TestRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results, cls.utils = generator.build(MANIFEST['design'], seed=MANIFEST['seed'])

    def test_connectomes_match_manifest(self):
        self.assertEqual(sorted(self.results), sorted(MANIFEST['connectomes']))

    def test_cells(self):
        for name, (_, arrays, cells) in self.results.items():
            ref_cells = MANIFEST['connectomes'][name]['cells']
            self.assertEqual([c['key'] for c in cells], [c['key'] for c in ref_cells], name)
            with np.load(REGRESSION / (name + '.npz')) as stored:
                for cell, ref in zip(cells, ref_cells):
                    with self.subTest(connectome=name, cell=ref['key']):
                        self.assertEqual(cell['status'], ref['status'])
                        if ref['status'] == 'raised':
                            self.assertEqual(cell['exception'], ref['exception'])
                            continue
                        self.assertEqual(cell['finite'], ref['finite'])
                        recon = stored[ref['key'] + '__err'][1]
                        borderline = THRESHOLD / 10 <= recon <= THRESHOLD * 10
                        if not borderline:
                            self.assertEqual(cell['completes'], ref['completes'])

    def test_control_task_outputs(self):
        for name, (_, arrays, _) in self.results.items():
            cells = {c['key']: c for c in MANIFEST['connectomes'][name]['cells']}
            with np.load(REGRESSION / (name + '.npz')) as stored:
                for key in sorted(k for k in stored.files if k.startswith('cell')):
                    cell_key, part = key.split('__')
                    cell = cells[cell_key]
                    rtol = RTOL_COMPLETE if cell.get('completes') else RTOL_INCOMPLETE
                    label = '{0}/{1}'.format(name, key)
                    with self.subTest(array=label):
                        self.assertIn(key, arrays)
                        got, ref = arrays[key], stored[key]
                        if part == 'x_index':
                            np.testing.assert_array_equal(got, ref, err_msg=label)
                        elif part == 'err':
                            assert_same_magnitude(self, got, ref, label)
                        elif part == 'node_energy':
                            zero = ref == 0
                            np.testing.assert_array_equal(got[zero], 0, err_msg=label + ': zero energies changed')
                            assert_close(self, got[~zero], ref[~zero], rtol, label)
                        else:  # x
                            assert_close(self, got, ref, rtol, label)

    def test_other_outputs(self):
        for name, (_, arrays, _) in self.results.items():
            with np.load(REGRESSION / (name + '.npz')) as stored:
                keys = sorted(k for k in stored.files if not k.startswith('cell'))
                self.assertEqual(sorted(k for k in arrays if not k.startswith('cell')), keys, name)
                for key in keys:
                    label = '{0}/{1}'.format(name, key)
                    with self.subTest(array=label):
                        if key.startswith(INPUT_PREFIXES) or key.endswith(INPUT_SUFFIXES):
                            assert_same_input(self, arrays[key], stored[key], label)
                        else:
                            assert_close(self, arrays[key], stored[key], RTOL_OTHER, label)

    def test_utils(self):
        with np.load(REGRESSION / 'utils.npz') as stored:
            self.assertEqual(sorted(self.utils), sorted(stored.files))
            for key in stored.files:
                label = 'utils/' + key
                with self.subTest(array=label):
                    got, ref = np.asarray(self.utils[key]), stored[key]
                    if key in UTILS_INPUTS:
                        assert_same_input(self, got, ref, label)
                    elif ref.dtype.kind in 'iubUS':
                        np.testing.assert_array_equal(got, ref, err_msg=label)
                    else:
                        assert_close(self, got, ref, RTOL_OTHER, label)


if __name__ == '__main__':
    unittest.main()
