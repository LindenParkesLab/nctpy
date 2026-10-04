# Controllability metrics

Controllability metrics describe each node's general capacity to steer the system, without reference to any
particular transition. nctpy provides three, all from Gu et al. (*Nat Commun* 2015).

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

## Boundary controllability

{func}`~nctpy.metrics.boundary_control` ranks nodes by how early they are found on the boundaries between network
communities, as the communities are split recursively: an indicator of a node's role in coupling or decoupling
modules. It needs a partition of the nodes into communities, for example functional systems or the result of a
modularity analysis.

```python
from nctpy.metrics import boundary_control

bc = boundary_control(A, partition, system="discrete", rho=0.2)  # A as measured: not normalised
```

The algorithm (Gu et al., following Pasqualetti et al. 2014) starts from the partition's own boundary nodes, those
with enough connection weight (`rho`, relative to the largest weight) to other communities. It then repeatedly
splits the least controllable community, the one whose nodes as control nodes give the Gramian with the smallest
eigenvalue, by its Fiedler vector, and adds that split's boundary nodes. Nodes found earlier get higher values: 1
for the first boundary nodes.

nctpy implements Gu et al.'s released code faithfully, and agrees with it exactly. Three of its properties matter
when you interpret the results, and are described in the function's documentation:

- the Gramians it compares are badly conditioned, so in large networks the order of splits, and of the values, can
  depend on rounding error;
- a node found again in a later split moves to the newest tier, as in the code;
- `rho` = 0.2 was chosen for Gu et al.'s parcellation, and may not suit others.

`A` is passed as measured, not normalised: the algorithm scales and normalises it itself, as Gu et al.'s code does.
In continuous time it uses the continuous-time Gramian, the continuous form of Gu et al.'s Supplementary Methods.

## Relationships between the metrics

Across nodes, average controllability tends to rise with node strength (weighted degree), and modal
controllability to fall with it. The {ref}`metric_correlations` example shows both on a synthetic connectome, as
Gu et al. reported for human connectomes.
