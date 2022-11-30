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
Perhaps the most simple approach would be to use built-in numerical integrators from various sources, whether it be Python, MATLAB, Mathematica, etc. In MATLAB, the command would look something like

.. code-block:: matlab

    T = 1;
    Wc = f(t) = integral(@(t) (expm(A*t)*B)*(expm(A*t)*B)',0,T,'ArrayValued',1);



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
While the Riemann sum is simple, it is also prone to errors if the function being integrated changes too quickly with respect to the time step, and might require too small of a time step :math:`\Delta T.` Now, taking :math:`T / \Delta T` products of matrices shouldn't be too computationally expensive given the system is not too large. However, another issue arises when :math:`A\Delta T` becomes too small to evaluate to numerical precision. For example, if we require 10,000 time steps to accurately capture the curvature of the matrix exponential, then we need to accurately compute :math:`e^{A/10,000},` and then accurately multiply those very small matrices 10,000 times. This approach can lead to the exonential growth of numerical precision errors.

Instead, we can use `Simpson's rule <https://en.wikipedia.org/wiki/Simpson%27s_rule>`_, which essentially uses higher-order polynomials to fit the curvature of functions. So, rather than assuming the function stays constant at each sampled point as in the Riemann sum, we instead fit polynomials, which ultimately evaluates to

.. math::
    W_c \approx \frac{\Delta T}{3} \left( f(0) + 2\sum_{n=1}^{\frac{T}{2\Delta T}-1} f(2n\Delta T) + 4\sum_{j=1}^{\frac{T}{2\Delta T}} f((2n-1)\Delta T) + f(T) \right).

More advanced versions of this polynomial integration scheme can be found in the `Newton-Cotes formulas <https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas>`_. 

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

What this means is that for any analysis that involves assessing *many* state transitions for one set of system parameters :math:`A,B,T,` we only have to compute the Gramian *once*, and invert the Gramian *once*. After obtaining :math:`W_c^{-1},` we can evaluate all of the energies for all state transitions through simple matrix multiplications, which are computationally way more efficient.

We can take this idea one step further and notice that the matrix exponential, :math:`e^{AT},` also only has to be evaluated *once*. Hence, as a tip to our readers, it is likely going to be much more computationally efficient to evaluate and store the Gramians and matrix exponentials, then batch together the state transitions.













