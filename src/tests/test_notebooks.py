"""Notebook test: the notebooks in scripts/ execute.

Each notebook is run top to bottom with nbclient, from scripts/, so its projdir (the parent of the
working directory) is this checkout. The files are not modified; before running, a few listed
substitutions are made in memory, each required to match an exact number of times:

- resultsdir: pointed at a temporary directory, so cached results in results/ are never overwritten.
- the 5000-permutation null models are cut to 3, and every `run = False` becomes `run = True`, so
  cells that would load cached nulls from results/ compute them (quickly) instead.

The notebooks import the installed nctpy; the kernel runs with this checkout's src/ first on
PYTHONPATH, so they test this checkout.

The notebooks read the paper's data from data/, which is not in the repository and may not be
shared, so this test runs only where a local copy exists and skips everywhere else. It checks
execution only; printed values are checked by test_paper_code.py.
"""
import os
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / 'scripts'
DATADIR = REPO / 'data'
ANNOT_DIR = DATADIR / 'schaefer_parc' / 'fsaverage5'
TIMEOUT = 1800  # seconds per cell
N_PERMS = 3

REQUIRED_DATA = [
    DATADIR / 'pnc_schaefer200_Am.npy',
    DATADIR / 'pnc_schaefer200_system_labels.txt',
    DATADIR / 'pnc_schaefer200_rsts.npy',
    DATADIR / 'pnc_schaefer200_centroids.csv',
    DATADIR / 'schaefer200_cyto.npy',
    DATADIR / 'schaefer200_micro.npy',
    DATADIR / 'normalized_connection_strength_ipsi.csv',
    DATADIR / 'connection_strength_ipsi.csv',
    DATADIR / 'region_class.csv',
    ANNOT_DIR / 'lh.Schaefer2018_200Parcels_7Networks_order.annot',
    ANNOT_DIR / 'rh.Schaefer2018_200Parcels_7Networks_order.annot',
]

RESULTSDIR = "resultsdir = os.path.join(projdir, 'results')"


def setUpModule():
    absent = [str(p.relative_to(REPO)) for p in REQUIRED_DATA if not p.exists()]
    if absent:
        raise unittest.SkipTest('the paper\'s data are not available locally: ' + ', '.join(absent))
    missing = []
    for name in ('nbformat', 'nbclient', 'ipykernel', 'pandas', 'sklearn', 'seaborn', 'nilearn', 'nibabel'):
        try:
            __import__(name)
        except ImportError:
            missing.append(name)
    if missing:
        raise unittest.SkipTest('packages needed to run the notebooks are not installed: ' + ', '.join(missing))


def run_notebook(test, filename, n_perms=0, n_run_flags=0):
    """Execute scripts/<filename> after the listed substitutions; fail on any cell error."""
    import nbformat
    from nbclient import NotebookClient

    nb = nbformat.read(SCRIPTS / filename, as_version=4)
    code_cells = [c for c in nb.cells if c.cell_type == 'code']
    with tempfile.TemporaryDirectory() as results:
        substitutions = [
            (RESULTSDIR, 'resultsdir = {0!r}'.format(results), 1),
            ('n_perms = 5000', 'n_perms = {0}'.format(N_PERMS), n_perms),
            ('run = False', 'run = True', n_run_flags),
        ]
        for old, new, count in substitutions:
            found = sum(c.source.count(old) for c in code_cells)
            test.assertEqual(found, count, '{0}: substitution {1!r} matched {2} times'.format(filename, old[:40], found))
            for c in code_cells:
                c.source = c.source.replace(old, new)
        env = {'PYTHONPATH': os.pathsep.join(filter(None, [str(REPO / 'src'), os.environ.get('PYTHONPATH')])),
               'MPLBACKEND': 'Agg'}
        saved = {k: os.environ.get(k) for k in env}
        os.environ.update(env)  # inherited by the kernel process
        try:
            NotebookClient(nb, timeout=TIMEOUT, kernel_name='python3',
                           resources={'metadata': {'path': str(SCRIPTS)}}).execute()
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v


class TestNotebooks(unittest.TestCase):
    def test_get_control_inputs(self):
        run_notebook(self, 'get_control_inputs.ipynb')

    def test_control_energy_wrapper(self):
        run_notebook(self, 'control_energy_wrapper.ipynb')

    def test_path_a_control_energy_binary(self):
        run_notebook(self, 'path_a_control_energy_binary.ipynb', n_perms=2, n_run_flags=2)

    def test_path_a_control_energy_nonbinary(self):
        run_notebook(self, 'path_a_control_energy_nonbinary.ipynb', n_perms=1, n_run_flags=2)

    def test_path_a_control_energy_directed(self):
        run_notebook(self, 'path_a_control_energy_directed.ipynb')

    def test_path_b_ave_ctrb(self):
        run_notebook(self, 'path_b_ave_ctrb.ipynb', n_perms=1, n_run_flags=2)


if __name__ == '__main__':
    unittest.main()
