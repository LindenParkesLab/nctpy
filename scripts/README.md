# Scripts

Notebooks that reproduce the analyses in the Nature Protocols paper (Parkes, Kim et al., 2024), plus two
scripts that draw its schematic figures.

To run the notebooks:

1. From the repository's root, install nctpy with the packages the notebooks use:
   `pip install -e ".[paper]"`.
2. Put the data in `data/` at the repository's root. The PNC data used in the paper are not distributed
   with nctpy.
3. Open the notebooks from this folder. Each finds the repository's root from its own location, and saves
   figures and results to `results/`, which it creates if needed.

The null-model analyses run 5000 permutations, which takes hours. Those cells are off by default
(`run = False`). Set `run = True` in such a cell to run it once; it saves its nulls to `results/`. The cells
that plot those nulls run when the saved nulls are present, and say so when they are not.

The version of these notebooks linked from the paper is at commit
[`f69ec009`](https://github.com/LindenParkesLab/nctpy/tree/f69ec009d70a46cb019da7c59a0d00b3e254731a/scripts).
