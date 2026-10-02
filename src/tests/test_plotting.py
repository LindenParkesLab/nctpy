"""Smoke tests for nctpy.plotting that need no data, with the non-interactive Agg backend.

surface_plot is drawn on nilearn's bundled fsaverage5 surface, using a synthetic parcellation
written to a temporary FreeSurfer annotation file. It is checked with signed, non-integer data:
nilearn >= 0.13 rejects both in plot_surf_roi, which surface_plot used to call.
"""
import os
import tempfile
import unittest
import warnings

import numpy as np

os.environ.setdefault('MPLBACKEND', 'Agg')

N_VERTICES = 10242  # per hemisphere on fsaverage5
N_PARCELS = 20  # per hemisphere


def setUpModule():
    try:
        import nctpy.plotting  # noqa: F401
    except ImportError as exc:
        raise unittest.SkipTest('plotting dependencies not installed: {0}'.format(exc))


class TestSurfacePlot(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import nibabel as nib
        from nilearn import datasets
        cls.fsaverage = datasets.fetch_surf_fsaverage(mesh='fsaverage5')
        cls.tmp = tempfile.TemporaryDirectory()
        rng = np.random.default_rng(0)
        cls.annot = {}
        for hemi in ('lh', 'rh'):
            labels = rng.integers(0, N_PARCELS + 1, size=N_VERTICES)  # label 0 = unassigned, like the medial wall
            ctab = np.c_[rng.integers(0, 256, size=(N_PARCELS + 1, 3)), np.zeros(N_PARCELS + 1, dtype=int)]
            names = ['unknown'] + ['parcel_{0}'.format(i) for i in range(1, N_PARCELS + 1)]
            path = os.path.join(cls.tmp.name, '{0}.synthetic.annot'.format(hemi))
            nib.freesurfer.write_annot(path, labels, ctab, names, fill_ctab=True)
            cls.annot[hemi] = path

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def draw(self, data, cmap):
        import matplotlib.pyplot as plt
        from nctpy.plotting import surface_plot
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            f = surface_plot(data, self.annot['lh'], self.annot['rh'], fsaverage=self.fsaverage, order='lr', cmap=cmap)
        self.assertEqual(len(f.axes), 5)  # four surface views and a colorbar
        plt.close(f)

    def test_signed_non_integer_data(self):
        self.draw(np.random.default_rng(1).normal(size=2 * N_PARCELS), cmap='coolwarm')

    def test_positive_non_integer_data(self):
        self.draw(np.random.default_rng(2).random(2 * N_PARCELS) + 0.5, cmap='viridis')

    def test_default_colour_limits_symmetric_for_coolwarm(self):
        # with cmap='coolwarm' and no cblim, the colour scale is symmetric about zero
        import matplotlib.pyplot as plt
        from nctpy.plotting import surface_plot
        data = np.r_[np.linspace(-0.3, 1.2, N_PARCELS), np.zeros(N_PARCELS)]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            f = surface_plot(data, self.annot['lh'], self.annot['rh'], fsaverage=self.fsaverage, cmap='coolwarm')
        vmin, vmax = f.axes[-1].get_ylim()
        self.assertAlmostEqual(vmin, -vmax)
        plt.close(f)

    def test_fsaverage_fetched_when_drawing(self):
        # fsaverage=None (the default) loads fsaverage5 when the plot is drawn, not when the module is imported
        import matplotlib.pyplot as plt
        from unittest import mock
        from nctpy import plotting
        data = np.random.default_rng(3).random(2 * N_PARCELS)
        with mock.patch.object(plotting.datasets, 'fetch_surf_fsaverage', return_value=self.fsaverage) as fetch, \
                warnings.catch_warnings():
            warnings.simplefilter('ignore')
            plt.close(plotting.surface_plot(data, self.annot['lh'], self.annot['rh']))
            fetch.assert_called_once_with(mesh='fsaverage5')
            # passed positionally, as the fourth argument, fsaverage is used as given
            fetch.reset_mock()
            plt.close(plotting.surface_plot(data, self.annot['lh'], self.annot['rh'], self.fsaverage))
            fetch.assert_not_called()


class TestImport(unittest.TestCase):
    def test_import_does_not_fetch(self):
        # a fresh interpreter in which fetching fsaverage fails: importing nctpy.plotting must still work
        import subprocess
        import sys
        script = (
            "import nilearn.datasets\n"
            "def fail(*args, **kwargs):\n"
            "    raise RuntimeError('fetched on import')\n"
            "nilearn.datasets.fetch_surf_fsaverage = fail\n"
            "import nctpy.plotting\n"
            "print('imported')\n"
        )
        import nctpy
        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            filter(None, [os.path.dirname(os.path.dirname(nctpy.__file__)), os.environ.get('PYTHONPATH')])))
        result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('imported', result.stdout)


class TestOtherPlots(unittest.TestCase):
    """Each function draws without error on synthetic data."""

    def setUp(self):
        import matplotlib.pyplot as plt
        self.plt = plt
        self.f, self.ax = plt.subplots()
        self.rng = np.random.default_rng(4)

    def tearDown(self):
        self.plt.close('all')

    def test_set_plotting_params(self):
        from nctpy.plotting import set_plotting_params
        set_plotting_params(format='svg')
        self.assertEqual(self.plt.rcParams['savefig.format'], 'svg')
        self.assertEqual(self.plt.rcParams['pdf.fonttype'], 42)

    def test_reg_plot(self):
        from nctpy.plotting import reg_plot
        x = self.rng.normal(size=40)
        y = x + self.rng.normal(size=40)
        for annotate in ('pearson', 'spearman', 'both', (0.5, 0.01), None):
            with self.subTest(annotate=annotate):
                self.ax.clear()
                reg_plot(x, y, 'x', 'y', self.ax, annotate=annotate)
                self.assertEqual(self.ax.get_xlabel(), 'x')
                self.assertEqual(len(self.ax.texts), 0 if annotate is None else 1)
        self.ax.clear()
        X = self.rng.normal(size=(8, 8))
        reg_plot(X, X.T + self.rng.normal(size=(8, 8)), 'x', 'y', self.ax, c=self.rng.random((8, 8)), kde=False)
        self.assertEqual(len(self.ax.collections[-1].get_offsets()), 8 * 8 - 8)  # off-diagonal entries only

    def test_null_plot(self):
        from nctpy.plotting import null_plot
        null_plot(2.5, self.rng.normal(size=200), 'statistic', self.ax, p_val=0.01)
        self.assertEqual(len(self.ax.texts), 2)

    def test_roi_to_vtx(self):
        import nibabel as nib
        from nctpy.plotting import roi_to_vtx
        labels = np.array([0, 1, 1, 2, 3, 3, 0])
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'lh.test.annot')
            # distinct, non-black colours: black reads back as -1 (unlabelled), not 0
            ctab = np.c_[np.arange(1, 5) * 50, np.arange(1, 5) * 40, np.arange(1, 5) * 30, np.zeros(4, dtype=int)]
            nib.freesurfer.write_annot(path, labels, ctab, ['unknown', 'a', 'b', 'c'], fill_ctab=True)
            vtx, vmin, vmax = roi_to_vtx(np.array([10.0, 20.0, 30.0]), path)
        np.testing.assert_array_equal(vtx, [0, 10, 10, 20, 30, 30, 0])
        self.assertEqual((vmin, vmax), (0, 30))

    def test_roi_to_vtx_unlabelled_vertices_are_background(self):
        # a black colour-table entry reads back as -1 (unlabelled); those vertices stay 0. Until 1.1 they took
        # roi_data[-2], the second-to-last parcel's value
        import nibabel as nib
        from nctpy.plotting import roi_to_vtx
        labels = np.array([0, 1, 1, 2, 3, 3, 0])
        ctab = np.c_[np.arange(4) * 50, np.arange(4) * 40, np.arange(4) * 30, np.zeros(4, dtype=int)]  # row 0 black
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'lh.test.annot')
            nib.freesurfer.write_annot(path, labels, ctab, ['unknown', 'a', 'b', 'c'], fill_ctab=True)
            self.assertEqual(nib.freesurfer.read_annot(path)[0][0], -1)
            vtx, vmin, vmax = roi_to_vtx(np.array([10.0, 20.0, 30.0]), path)
        np.testing.assert_array_equal(vtx, [0, 10, 10, 20, 30, 30, 0])
        self.assertEqual((vmin, vmax), (0, 30))

    def test_add_module_lines(self):
        try:
            import pandas as pd
        except ImportError:
            self.skipTest('pandas not installed')
        import contextlib
        import io
        from nctpy.plotting import add_module_lines
        modules = pd.Series(['Vis'] * 3 + ['SomMot'] * 2 + ['Default'] * 4)
        with contextlib.redirect_stdout(io.StringIO()):
            add_module_lines(modules, self.ax)
        self.assertEqual(len(self.ax.collections), 4 * 3)  # four lines per module


if __name__ == '__main__':
    unittest.main()
