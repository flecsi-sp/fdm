.. |br| raw:: html

   <br />

.. _fdm_problem:

Example Problem
***************

The solvers in this example solve the *Poisson* problem with *Dirichlet*
boundary conditions

.. math::

   \begin{align}
   -\mu\Delta u &= f \textrm{ in } \Omega \label{prb}\tag{1}\\
              u &= g \textrm{ on } \Gamma_\mathcal{D}.
   \end{align}

The continuous form is discretized using second-order,
centered-differences, e.g.

.. math::

   \frac{\delta^2 u}{\delta x^2} \approx
   \frac{u_{i+1,j} - 2 u_{i,j} + u_{i-1,j}}{\Delta x^2}.

In 2D with variable spacing across axes, this gives the discrete form

.. math::

    -\mu\left(\frac{u_{i+1,j} - 2 u_{i,j} + u_{i-1,j}}{\Delta x^2} +
    \frac{u_{i,j+1} - 2 u_{i,j} + u_{i,j-1}}{\Delta y^2}\right) =
    f_{i,j}.

.. note:: The sign of the operator is convenient for defining a
   smoothing update as shown in the next section.

Elementary Iterative Solvers
****************************

Notice that (:math:`\ref{prb}`) has the form

.. math::

   A\mathbf{u}=f.\label{alg}\tag{2}

For the moment, assume that we can write :math:`A` as

.. math::

   A=N^{-1}\left(I-M\right)

where :math:`N` is an invertible :math:`n\times n` matrix. Substituting this
form for :math:`A` in (:math:`\ref{alg}`) we can rewrite
(:math:`\ref{prb}`) as

.. math::

   \mathbf{u} = M\mathbf{u} + N\mathbf{f}.\label{mform}\tag{3}

This is relevant because a linear, stationary iterative method can be
expressed in the discrete form

.. math::

   \mathbf{v}^{k+1} = M\mathbf{v}^k + N\mathbf{f}.\label{dmform}\tag{4}

Using recursion, we have

.. math::

   \begin{align}
   \mathbf{v}^{k+1} &= M\mathbf{v}^k + N\mathbf{f} \\
                    &= M\left(M\mathbf{v}^{k-1}+N\mathbf{f}\right) +
                    N\mathbf{f} \\
                    \vdots \\
                    &= M^{k+1}\mathbf{v}^{(0)} + \left(M^k + M^{k-1} + \dots
                    + I\right)N\mathbf{f}
   \end{align}

where :math:`\mathbf{v}^{(0)}` is the initial guess. Assymptotically, as
:math:`k\rightarrow\infty`,

.. math::

   \lim_{k\rightarrow\infty} \mathbf{v}^{k+1} =
   \lim_{k\rightarrow\infty} M^{k+1}\mathbf{v}^{(0)} +
   \lim_{k\rightarrow\infty} \left( \sum_{j=0}^{k}M^j\right)N\mathbf{f}.

Consider that :math:`\rho(M)<1` implies that

.. math::

   \lim_{k\rightarrow\infty}M^{k+1}=0.

Let :math:`\lambda\in\sigma\left(M\right)`, then
:math:`1-\lambda \in \sigma\left(I-M\right)`. Since
:math:`\lambda\leq\rho(M)<1`, it must be true that
:math:`1-\lambda\ne 0`, which implies that :math:`\left(I-M\right)^{-1}`
exists.

Defining :math:`\mathcal{S}_k=\left(\sum_{j=0}^k M^j\right)`,
we see that

.. math::

   \left(I-M\right)\mathcal{S}_k = \left(I+M+M^2 + \dots + M^{k+1}\right) -
   \left(M+M^2 + \dots + M^{k+1}\right) = \left(I-M^{k+1}\right),

which implies that

.. math::

   \lim_{k\rightarrow\infty}\left(I-M\right)\mathcal{S}_k =
   \lim_{k\rightarrow\infty}\left(I-M^{k+1}\right) = I.

This shows that :math:`\mathcal{S}_k` converges to
:math:`\left(I-M\right)^{-1}`. Solving (:math:`\ref{mform}`)
algebraically, we see that

.. math::

   \mathbf{u} = \left(I-M\right)^{-1}N\mathbf{f}.

Therefore, we have shown that, asymptotically, the algebraic solution of
(:math:`\ref{mform}`) is recovered by the iteration
(:math:`\ref{dmform}`). It remains to show that we can find suitable
:math:`M` and :math:`N` matrices.

.. note::

   :math:`M` and :math:`N` are generally
   referred to as the *iteration* matrix and the *approximate inverse*,
   respectively.

Jacobi Method
~~~~~~~~~~~~~

Splitting :math:`A=D-L-U` into its *diagonal*, and *lower* and *upper*
triangular parts, and substituting into (:math:`\ref{prb}`), we can write

.. math::

   \mathbf{u} = D^{-1}\left(L+U\right)\mathbf{u}+D^{-1}\mathbf{f}.

Defining :math:`M_J=D^{-1}\left(L+U\right)` (the *iteration* matrix)
and :math:`N_J=D^{-1}` (the *approximate inverse*), we can write
(:math:`\ref{dmform}`) as

.. math::
   \mathbf{v}^{k+1} = M_J\mathbf{v}^k + N_J\mathbf{f}.
   \label{jacobi}\tag{5}

This is the matrix form of the *Jacobi* method. A useful modification to
(:math:`\ref{jacobi}`) is to add a *damping* factor
:math:`\omega\in\left(0,1\right]` and defining an iteration that is the
weighted average of the current and previous approximations. If we apply
the standard Jacobi method

.. math::

   \mathbf{\hat{v}}^{k+1} = M_J\mathbf{v}^k + N_J\mathbf{f},

then the :math:`\omega`-Jacobi method is given by

.. math::

   \mathbf{v}^{k+1} = \left(1-\omega\right)\mathbf{v}^k +
   \omega\mathbf{\hat{v}}^{k+1}\label{djacobi}\tag{6},

with iterartion matrix and approximate inverse
:math:`M_{J_{\omega}} = \left(1-\omega\right)I + \omega M_J`
and :math:`N_{J_{\omega}}=\omega D^{-1}`, respectively.
Again, :math:`\rho\left(M_{J_{\omega}}\right)<1` ensures convergence.

Component Form (Jacobi Method)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The matrix form of an iterative method is useful for analysis. However,
when actually implementing the method in a computer code, it is
prefereable to use the *component* form. For the Jacobi method, the
component form is

.. math::

   v_i^{k+1} = \frac{1}{a_{i,j}}\left(\sum_{j\ne i}^{n} a_{i,j}v_j^k +
   f_i\right).

Applying the Jacobi method to (:math:`\ref{prb}`) with :math:`\mu=1` for
a two-dimensional domain with variable axis spacing, we have the
component form

.. math::

   \hat{v}_{i,j}^{k+1} =
   \frac{1}{2\left(\frac{\Delta x}{\Delta y} +
     \frac{\Delta y}{\Delta x}\right)}
   \left[
   \Delta x \Delta y f_{i,j} +
   \frac{\Delta y}{\Delta x}
     \left(v_{i+1,j}^k + v_{i-1,j}^k\right) +
   \frac{\Delta x}{\Delta y}
     \left(v_{i,j+1}^k + v_{i,j-1}^k\right) +
   \right].

The damped version follows (:math:`\ref{djacobi}`). This is one of the
methods used in the code for this specialization example.

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :
