"""Smoke tests for nctpy.plotting that need no data.

surface_plot is drawn on nilearn's bundled fsaverage5 surface, using a synthetic parcellation
written to a temporary FreeSurfer annotation file. It is checked with signed, non-integer data:
nilearn >= 0.13 rejects both in plot_surf_roi, which surface_plot used to call (Roadmap D17).
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


if __name__ == '__main__':
    unittest.main()
