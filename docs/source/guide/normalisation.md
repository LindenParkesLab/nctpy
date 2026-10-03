# Normalising the connectome

A structural connectome cannot be used as the system matrix as it is: activity would grow without bound. nctpy
scales it so that, without input, activity always decays, using {func}`~nctpy.utils.matrix_normalization`:

$$
\text{discrete time:}\quad A_\text{norm} = \frac{A}{\lambda_{\max} + c}
\qquad\qquad
\text{continuous time:}\quad A_\text{norm} = \frac{A}{\lambda_{\max} + c} - I
$$

where $\lambda_{\max}$ is the spectral radius of $A$ (its largest absolute eigenvalue) and $c$ is a positive
constant, 1 by default.

```python
from nctpy.utils import matrix_normalization

A_norm = matrix_normalization(A, system="continuous", c=1)
```

Dividing by $\lambda_{\max} + c$ brings every eigenvalue inside the unit circle, which makes a discrete-time system
stable. In continuous time, subtracting the identity then shifts every eigenvalue to have a negative real part,
which makes it stable: each node's activity decays at rate 1 unless driven. Any positive $c$ guarantees stability;
a larger $c$ makes activity decay faster. Because the whole matrix is scaled by one number (and shifted by the
identity), the normalisation keeps the order of the eigenvalues and the eigenvectors, so results stay comparable
between studies.

## Several connectomes with one normalisation

By default each connectome is divided by its own spectral radius. When comparing connectomes, such as one per
subject, you may want them all divided by the same number, so that differences in overall connection strength are
kept. Pass that number as `l`, typically the largest spectral radius among them:

```python
l = max(np.max(np.abs(np.linalg.eigvals(A))) for A in connectomes)
A_norms = [matrix_normalization(A, system="continuous", c=1, l=l) for A in connectomes]
```

Stability is guaranteed only when $c + l$ exceeds each connectome's own spectral radius, which the largest one
does. Nothing checks this: an unstable system still returns results.

## Self-connections

The model gives each node its own dynamics through the normalisation (in continuous time, the $-I$ term), so it
assumes the connectome has no self-connections: a zero diagonal. nctpy does not enforce this. A connectome with a
non-zero diagonal is used as given, and in continuous time its self-connections make the effective decay rates
non-uniform: node $i$ decays at rate $1 - A_{ii}/(\lambda_{\max} + c)$.

The protocol paper's connectome has self-connections, and its printed results keep them. To remove them, ask:

```python
A_norm = matrix_normalization(A, system="continuous", zero_diagonal=True)
```

On the paper's data that changes the energy of its example transition from 2604.71 to 2638.03. nctpy never
removes self-connections unless asked, and never warns about them; whether your model should have them is your
decision.

## Decay rates that differ between nodes

In continuous time, nodes can be given their own decay rates with `decay`, subtracted from the diagonal in place of
the identity:

$$
A_\text{norm} = \frac{A}{\lambda_{\max} + c} - \operatorname{diag}(\text{decay})
$$

```python
A_norm = matrix_normalization(A, system="continuous", decay=decay_rates)  # one rate per node, or a single number
```

A larger value means activity at that node decays faster. The default, 1 at every node, is the $-I$ above.
Uniform rates of at least 1 keep the system stable, as do rates of at least 1 at every node for an undirected
connectome; smaller rates may not, and nothing checks. `decay` is not available in discrete time, which has no
diagonal term to replace.

To fit decay rates to a state transition, as Kim et al. (2025) did, see {func}`~nctpy.optimize.optimize_decay_rates`
and the {doc}`/tutorials/decay_rates` tutorial. A value $v$ in that paper's Fig. 2B corresponds to
`decay = 1 - v` here.

## Checking stability

nctpy never refuses an unstable system: functions return their results, which will be meaningless. If you change
`c`, `l` or `decay`, check the normalised matrix:

```python
np.max(np.linalg.eigvals(A_norm).real) < 0  # continuous time: every eigenvalue has a negative real part
np.max(np.abs(np.linalg.eigvals(A_norm))) < 1  # discrete time: every eigenvalue is inside the unit circle
```
