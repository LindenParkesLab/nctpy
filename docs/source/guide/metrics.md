# Controllability metrics

Controllability metrics describe each node's general capacity to steer the system, without reference to any
particular transition. nctpy provides two, both from Gu et al. (*Nat Commun* 2015).

## Average controllability

{func}`~nctpy.metrics.ave_control` measures how strongly the system responds to input at a node: the energy that a
unit impulse at node $i$ spreads into the network, which is the trace of the controllability Gramian when input
enters at node $i$ alone.

```python
from nctpy.metrics import ave_control

ac = ave_control(A_norm, system="continuous")  # one value per node
```

- In **continuous time** it is $\int_0^1 \lVert e^{A t} \mathbf{e}_i \rVert^2 \,\mathrm{d}t$, over a horizon of 1.
- In **discrete time** it is $\sum_{k \ge 0} \lVert A^k \mathbf{e}_i \rVert^2$, over an infinite horizon.

Here $\mathbf{e}_i$ is the unit vector of node $i$, and `A_norm` must be normalised for the time system you pass
(see {doc}`normalisation`).

**Directed connectomes.** nctpy reads `A[i, j]` as the connection from node $j$ to node $i$, so input at node $i$
spreads along column $i$. For an undirected connectome the direction does not matter. For a directed one, nctpy
follows Gu et al.'s definition in both time systems; versions before 1.2 computed something else for directed
matrices (see the changelog).

## Modal controllability

{func}`~nctpy.metrics.modal_control` measures a node's ability to drive the system into states that are hard to
reach, those of its fast-decaying modes. With $U$ the eigenvectors of $A$ and $\lambda_j$ its eigenvalues:

- in **discrete time** (Gu et al. 2015): $\sum_j U_{ij}^2 (1 - \lambda_j^2)$;
- in **continuous time** (Karrer et al., *J Neural Eng* 2020, Eq. 8, adapted from the discrete form):
  $\sum_j U_{ij}^2 (1 - e^{\lambda_j})$.

```python
from nctpy.metrics import modal_control

mc = modal_control(matrix_normalization(A, system="discrete"))  # discrete time, the default
mc = modal_control(matrix_normalization(A, system="continuous"), system="continuous")
```

Normalise for the time system you pass. Discrete time is the default, for compatibility with earlier versions,
in which it was the only option. In both systems, each mode's weight lies between 0 and 1 for a stable system, so
every value does too.

Modal controllability is defined for **undirected** connectomes: nctpy computes it from the real Schur
decomposition, which equals the eigendecomposition for a symmetric matrix. For a directed one the values are only
an approximation, and depend on the order of the nodes.

## Relationships between the metrics

Across nodes, average controllability tends to rise with node strength (weighted degree), and modal
controllability to fall with it. The {ref}`metric_correlations` example shows both on a synthetic connectome, as
Gu et al. reported for human connectomes.
