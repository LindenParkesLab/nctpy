# New features: analyses on HCP data

Notebooks that develop and check nctpy's new functionality on structural connectomes and fMRI time series from the
Human Connectome Project (HCP-YA, young adults; HCP-D, development), parcellated with Schaefer 400 (7 networks).

The data are not in the repository and must never be committed: put them in `data/connectivity/` and
`data/timeseries/` at the repository's root. Every notebook writes its results and figures to
`results/new_features/`, which is git-ignored, because they are derived from the data.

## Environment

Use an environment with nctpy (with the `plot` extra) and the lab's [snaplab_tools](https://github.com/LindenParkesLab/snaplab_tools)
installed. The notebooks use snaplab_tools wherever it already provides what they need, for example the Yeo system
of each Schaefer parcel and correlation statistics.

## Notebooks

Open them from this folder: each finds the repository's root from its own location.

- `hcp_controllability.ipynb`: average, modal and boundary controllability (Gu et al., *Nat Commun* 2015) in
  discrete and continuous time, for the group connectome and every subject, summarised as Gu et al. did (ranks
  averaged across subjects) and compared with their results. Boundary controllability uses the Yeo 7 systems as its
  partition, where Gu et al. used data-driven communities. Set the dataset (HCP-YA or HCP-D) in its settings cell.
  The committed outputs are group-level summaries only.
