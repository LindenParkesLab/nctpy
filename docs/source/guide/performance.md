# Performance

Most analyses compute many transitions: every pair of brain states, every subject, thousands of null networks.
This page describes what makes nctpy fast or slow. To measure speed on your machine, run
`python benchmarks/transitions.py` from a clone of the repository: it times all 49 transitions between seven
states of a 200-node connectome in several ways.

## Reuse within a system

Most of the work of {func}`~nctpy.energies.get_control_inputs` depends only on the system, `A_norm`, `T`, `B`, `S`
and `rho`, not on the states. nctpy computes it once and reuses it while consecutive calls share that system. A
loop over the transitions of one system is therefore much faster than its first call suggests: on a 200-node
connectome, for example, a call that reuses the work can take about a third of the time of the first.

To benefit, keep calls with the same system together: loop over states inside a loop over systems, not the other
way round. Changing `B`, `S`, `rho` or `T` starts afresh, as in a perturbation analysis that changes one node's
weight at a time. Reuse never changes the results.

## Choosing the function

- {class}`~nctpy.pipelines.ComputeControlEnergy` solves the transitions that share a system together, several times
  faster than a loop over `get_control_inputs`, with the same energies to rounding (see {doc}`pipelines`).
- {func}`~nctpy.energies.minimum_energy_fast` computes minimum-control energies (no trajectory constraint) for many
  transitions at once, without simulating them, in continuous time (see {doc}`energies`).
- Discrete-time transitions are faster than continuous-time ones for short horizons.
- {func}`~nctpy.energies.minimum_energy_infinite` and {func}`~nctpy.energies.average_energy_infinite` need no
  simulation at all.

## System size

Run time grows with the number of nodes, by how much depends on the function and the machine. To see it for your
own size, pass the number of nodes to the benchmark: `python benchmarks/transitions.py 400`.

## Threads and parallel jobs

nctpy's matrix operations run on numpy's linear algebra library, which uses several CPU threads. When you
parallelise yourself, for example one process per subject or per batch of null networks, limit each process to one
thread so that they do not compete, by setting an environment variable before Python starts:

```bash
OMP_NUM_THREADS=1 python my_analysis.py
```

The same variable makes timings comparable between runs.

## Decay-rate optimisation

{func}`~nctpy.optimize.optimize_decay_rates` uses PyTorch and runs on a GPU when one is available, which helps for
large connectomes and for many transitions fitted together. On the CPU its runs are exactly reproducible.
