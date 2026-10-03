.. _numerics:

Numerics
----------------
When computing control energy, we mean that for a linear dynamical system of the form

.. math::
    \dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u},

we are computing the amount of energy it costs to bring the system from an initial state :math:`\mathbf{x}(0)` to a final state :math:`\mathbf{x}(T)` in :math:`T` time through a set of input functions :math:`\mathbf{u}(t).` This cost is called the *control energy*, and is defined by

.. math::
    E(\mathbf{u}(t)) = \int_0^T \mathbf{u}^\top(t) \mathbf{u}(t) \mathrm{d}(t) = \int_0^T u_1^2(t) + u_2^2(t) + \dotsm + u_k^2(t) \mathrm{d}t.

From the theory, we found that the *minimum control energy* is given by

.. math::
    E^* = (\mathbf{x}(T) - e^{AT}\mathbf{x}(0))^\top W_c^{-1} (\mathbf{x}(T) - e^{AT}\mathbf{x}(0)),

where :math:`W_c` is the *controllability Gramian*. Here, we will discuss how to numerically evaluate the theoretical relationships. Because the specific states :math:`\mathbf{x}(0)` and :math:`\mathbf{x}(T)` are just vectors, we will focus on the terms :math:`e^{AT}` and :math:`W_c.`

|
|

The Matrix Exponential
==========================
The matrix exponential has a long history in both analytical and numerical applications. But what *is* the exponential of a matrix? Fortunately, there is a relatively simple way to conceptualize it through a Taylor series expansion:

.. math::
    e^{A} = I + A + \frac{1}{2!}A^2 + \frac{1}{3!}A^3 + \dotsm + \frac{1}{n!} A^n + \dotsm = \sum_{n=0}^\infty \frac{1}{n!} A^n.

So in one sense, the matrix exponential is just a weighted sum of powers of :math:`A,` so nothing too fancy. Most numerical packages will have built in functions for the evaluation of matrix exponentials. Please note that the matrix exponential is *not* simply an element-wise exponential of :math:`A.` As we can see in the series expansion, the matrix exponential contains matrix powers of :math:`A,` which are very different from element-wise operations.

nctpy uses :func:`scipy.linalg.expm` (a Padé approximation with scaling and squaring). :func:`nctpy.energies.get_control_inputs` can instead use an eigendecomposition, :func:`nctpy.utils.expm`, with ``expm_version="eig"``.

|
|

The Controllability Gramian
==============================
Of great importance to linear network control is the controllability Gramian, defined in the theory as

.. math::
    W_c = \int_0^T e^{A(T-t)} B B^\top e^{A^\top (T-t)} \mathrm{d}t.

We can simplify this expression a bit through a simple change of variables :math:`\tau = T-t` (classic exercise left to reader) to write

.. math::
    W_c = \int_0^T e^{A\tau} BB^\top e^{A^\top \tau} \mathrm{d}\tau = \int_0^T f(\tau) \mathrm{d}\tau.

So... how do we turn this equation into a matrix on our computer? Well, it is no different than any other form of numerical integration! It just looks a bit scary because of the matrices, but never fear. Let us start with a very simple approach.



built-in numerical integrators
___________________________________
Perhaps the most simple approach would be to use a built-in numerical integrator that accepts matrix-valued functions. In Python, for example,

.. code-block:: python

    import scipy as sp

    T = 1
    Wc, _ = sp.integrate.quad_vec(lambda t: sp.linalg.expm(A * t) @ B @ B.T @ sp.linalg.expm(A.T * t), 0, T)



right-hand Riemann sum
___________________________________
Provided our particular numerical package lacks a matrix-valued numerical integrator, we can use the most basic integration scheme: the right-hand `Riemann sum <https://en.wikipedia.org/wiki/Riemann_sum>`_. The basic idea is that rather than compute the exact area under the curve, we can just evaluate the curve at different points in time, and pretend the function is constant in between points. So if we break up our integration into time steps of :math:`\Delta T,` then our Gramian approximately evaluates to

.. math::
    W_c \approx \sum_{n = 0}^{T/\Delta T-1} e^{An\Delta T} BB^\top e^{A^\top n\Delta T} \Delta T.

Because we assume our numerical package can evaluate matrix integrals and perform matrix multiplications, we can perform this summation without much issue. We can speed up this process by recognizing that the product of matrix exponentials yields the sum of their exponents (when the matrices in the exponentials `commute <https://en.wikipedia.org/wiki/Commuting_matrices>`_), such that

.. math::
    e^{An\Delta T} e^{A\Delta T} = e^{A(n+1)\Delta T}.

Hence, we only ever actually have to evaluate *one* matrix exponential, :math:`e^{A\Delta T},` and can obtain all subsequent matrix exponentials by accumulating products of matrix exponentials.


Simpson's rule
______________________
While the Riemann sum is simple, it is also prone to errors if the function being integrated changes too quickly with respect to the time step, and might require too small of a time step :math:`\Delta T.` Now, taking :math:`T / \Delta T` products of matrices shouldn't be too computationally expensive given the system is not too large. However, another issue arises when :math:`A\Delta T` becomes too small to evaluate to numerical precision. For example, if we require 10,000 time steps to accurately capture the curvature of the matrix exponential, then we need to accurately compute :math:`e^{A/10,000},` and then accurately multiply those very small matrices 10,000 times. This approach can lead to the exponential growth of numerical precision errors.

Instead, we can use `Simpson's rule <https://en.wikipedia.org/wiki/Simpson%27s_rule>`_, which essentially uses higher-order polynomials to fit the curvature of functions. So, rather than assuming the function stays constant at each sampled point as in the Riemann sum, we instead fit polynomials, which ultimately evaluates to

.. math::
    W_c \approx \frac{\Delta T}{3} \left( f(0) + 2\sum_{n=1}^{\frac{T}{2\Delta T}-1} f(2n\Delta T) + 4\sum_{n=1}^{\frac{T}{2\Delta T}} f((2n-1)\Delta T) + f(T) \right).

More advanced versions of this polynomial integration scheme can be found in the `Newton-Cotes formulas <https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas>`_.

what nctpy does
______________________
:func:`nctpy.energies.gramian`, which computes the Gramian with every node controlled (:math:`B = I`), combines the two ideas: it evaluates :math:`e^{A\Delta T}` once, with :math:`\Delta T = 0.001,` accumulates the products :math:`e^{An\Delta T}` step by step, and integrates with Simpson's rule (:func:`scipy.integrate.simpson`). In discrete time the Gramian is a sum rather than an integral, :math:`\sum_{k=0}^{T} A^k (A^\top)^k,` which nctpy accumulates in the same way.

the infinite horizon
______________________
For a stable system, the Gramian over an infinite horizon needs no integration at all: it solves the Lyapunov equation :math:`AW_c + W_cA^\top + BB^\top = 0` (in discrete time, :math:`AW_cA^\top - W_c + BB^\top = 0`), for which SciPy has direct solvers (:func:`scipy.linalg.solve_continuous_lyapunov`, :func:`scipy.linalg.solve_discrete_lyapunov`). nctpy uses them for ``gramian(A_norm, np.inf, system)``, :func:`nctpy.energies.minimum_energy_infinite` and :func:`nctpy.energies.average_energy_infinite`. When few nodes are controlled, :math:`W_c` is close to singular, and inverting it loses accuracy; see :doc:`/guide/energies`.

|
|

Evaluating Minimum Control Energy
=========================================
Now, once we have our controllability Gramian and state transitions, we evaluate the minimum control energy using

.. math::
    E^* = (\mathbf{x}(T) - e^{AT}\mathbf{x}(0))^\top W_c^{-1} (\mathbf{x}(T) - e^{AT}\mathbf{x}(0)).

But let's pause for a moment here. Notice that the controllability Gramian is *only* a function of the connectivity matrix :math:`A,` the input matrix :math:`B,` and the time horizon :math:`T` as we reproduce below

.. math::
    W_c = \int_0^T e^{A\tau} BB^\top e^{A^\top \tau} \mathrm{d}\tau.

What this means is that for any analysis that involves assessing *many* state transitions for one set of system parameters :math:`A,B,T,` we only have to compute the Gramian *once*, and invert the Gramian *once*. After obtaining :math:`W_c^{-1},` we can evaluate all of the energies for all state transitions through simple matrix multiplications, which are computationally way more efficient. The matrix exponential, :math:`e^{AT},` likewise only has to be evaluated *once*.

:func:`nctpy.energies.minimum_energy_fast` does exactly this: it computes :math:`W_c` (by Simpson's rule, as above) and :math:`e^{AT}` once per :math:`A, B, T,` reuses them while consecutive calls share them, and evaluates the energies of all the transitions it is given at once, using the pseudoinverse of :math:`W_c.` See :doc:`/guide/performance`.

|
|

Evaluating Optimal Control
=============================
:func:`nctpy.energies.get_control_inputs` solves the optimal control problem described in :doc:`../theory/index`, with the trajectory constraint :math:`S,` the weight :math:`\rho` and the reference state :math:`\mathbf{x}_r.`

continuous time
______________________
The state and the costate evolve together as :math:`\dot{\mathbf{z}} = M\mathbf{z} + \mathbf{c},` with :math:`\mathbf{z} = [\mathbf{x}; \mathbf{p}],` the matrix :math:`M = \begin{bmatrix} A & -BB^\top/(2\rho) \\ -2S & -A^\top \end{bmatrix}` and the constant :math:`\mathbf{c} = [\mathbf{0}; 2S\mathbf{x}_r].` nctpy

1. computes :math:`e^{MT}` and :math:`e^{M\Delta T},` with :math:`\Delta T = 0.001;`
2. writes the final state, :math:`\mathbf{x}(T),` in terms of the initial state and the unknown initial costate through the top block row of :math:`e^{MT},` and solves that linear system for :math:`\mathbf{p}(0)`;
3. steps :math:`\mathbf{z}` forward exactly, multiplying by :math:`e^{M\Delta T}` at each step, which gives the trajectory at :math:`T/\Delta T + 1` time points;
4. reads off the input, :math:`\mathbf{u}(t) = -B^\top\mathbf{p}(t)/(2\rho).`

The two numerical errors it returns come from steps 2 and 3: the *inversion error* is the residual of the linear solve for :math:`\mathbf{p}(0),` and the *reconstruction error* is the distance between the computed :math:`\mathbf{x}(T)` and the target state. See :doc:`/guide/errors`.

discrete time
______________________
In discrete time nctpy writes the state and costate equations of every time step, with the known :math:`\mathbf{x}(0)` and :math:`\mathbf{x}(T),` as one sparse linear system in all the unknown states and costates, and solves it with a sparse LU factorisation. The inversion error is the residual of that solve, and the reconstruction error is how far the solution departs from the dynamics.

reuse
______________________
The matrix exponentials of step 1, and in discrete time the system and its LU factorisation, depend only on :math:`A, T, B, S` and :math:`\rho,` not on the states. nctpy computes them once and reuses them while consecutive calls share them, so a loop over the transitions of one system is fast. Results are the same either way.

from inputs to energy
______________________
:func:`nctpy.energies.integrate_u` integrates each node's squared input with Simpson's rule, treating the time points as one unit apart. In continuous time they are :math:`\Delta T = 0.001` apart, so its energies are :math:`1/\Delta T = 1000` times the integral :math:`\int_0^T \lVert\mathbf{u}(t)\rVert^2 \mathrm{d}t` that the formulas above compute. See :doc:`/guide/energies`.













