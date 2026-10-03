# Null models

Is an observed energy or metric a property of the connectome's particular wiring, or would any network with
similar basic features give the same? Null models answer this by repeating the analysis on surrogates that keep
some features and randomise the rest.

## Null networks: `geomsurr`

{func}`~null_models.geomsurr.geomsurr` (Roberts et al., *NeuroImage* 2016) rewires a connectome while preserving its
spatial embedding: the relationship between the length of a connection and its weight. Import it as the protocol
paper does:

```python
from null_models.geomsurr import geomsurr

Wwp, Wsp, Wssp = geomsurr(W=adjacency, D=distance_matrix, seed=0)
```

`D` holds the distances between nodes, for example Euclidean distances between parcel centroids. Each call returns
three surrogates, all preserving the weight-distance relationship:

- `Wwp` also preserves the distribution of edge weights;
- `Wsp` also preserves the distribution of node strengths, assigned to nodes at random;
- `Wssp` also preserves each node's own strength.

`Wsp` and `Wssp` assume an undirected connectome. Normalise each surrogate and recompute the energy or metric to
build a null distribution, one surrogate per `seed`. The {doc}`/tutorials/protocol_workflow` tutorial does this.

`geomsurr` rewires edges between pairs of nodes, so its surrogates never have self-connections. It suits
connectomes without them. If your connectome has self-connections, compute the observed value on its zero-diagonal
version (for example `matrix_normalization(..., zero_diagonal=True)`), so that the observed and the surrogates
differ only in their wiring. `geomsurr` does not modify the matrix passed to it, and uses its own random number
generator, leaving numpy's global state untouched.

## p-values

- {func}`~nctpy.utils.get_null_p` gives the fraction of the null at or above the observed value (`version="standard"`),
  at or below it (`"reverse"`), or the smaller of the two (`"smallest"`); `abs=True` compares absolute values.
- {func}`~nctpy.utils.get_fdr_p` corrects a set of p-values, of any shape (for example one per transition), with
  the Benjamini-Hochberg false discovery rate.
- {func}`~nctpy.plotting.null_plot` draws a null distribution with the observed value and its p-value.

## Nulls for brain states and maps

Some questions concern the brain states rather than the connectome: would the energy be different for target states
with the same spatial smoothness but a different pattern? Surrogate maps that preserve spatial autocorrelation
answer this, and they are not part of nctpy. Kim et al. (2025) used [BrainSMASH](https://brainsmash.readthedocs.io)
to permute their target states. The lab's general toolbox,
[snaplab_tools](https://github.com/LindenParkesLab/snaplab_tools), wraps BrainSMASH in
`snaplab_tools.nulls.generate_surrogates`.
