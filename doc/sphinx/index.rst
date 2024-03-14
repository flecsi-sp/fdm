Finite Difference Method (FDM)
******************************

This specialization is designed for vertex-centered, structured-mesh
implementations of the finite difference method (FDM) using the
*flecsi::narray* core topology type.

Example Problem
===============

The solvers in this example solve the *Poisson* problem with *Dirichlet*
boundary conditions:

.. math::

   \begin{align}
   -\mu\Delta u &= f \textrm{ in } \Omega \\
              u &= g \textrm{ on } \Gamma_\mathcal{D}
   \end{align}

The continuous form is discretized using second-order,
centered-differences, e.g.:

.. math::

   \frac{\delta^2 u}{\delta x^2} \approx
   \frac{u_{i+1,j} - 2 u_{i,j} + u_{i-1,j}}{\Delta x^2}

Mesh Interface
==============

The mesh interface is defined in `fdm/mesh.hh`:

.. literalinclude:: fdm-src/fdm/mesh.hh
  :language: cpp
  :start-at: /// Mesh Interface
  :end-at: }; // struct interface

.. toctree::
  :caption: Contents:

  fdm/api

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :
