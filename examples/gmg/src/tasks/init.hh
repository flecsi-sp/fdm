/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_TASKS_INIT_HH
#define GMG_TASKS_INIT_HH

#include "../types.hh"

namespace gmg::task {

void enumerate(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa) noexcept;

/*!
  Initialize a field using a general plane equation, i.e., z = D + Mx + Ny.
 */
void bilinear(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double M,
  double N,
  double D) noexcept;

void constant(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double value) noexcept;

void fouriermodes(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) noexcept;

void fourier_fw(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) noexcept;

void fourier_interp(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) noexcept;

void fourier_residual(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) noexcept;

void gs_eigenvector(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) noexcept;

void eggcarton(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  field<double>::accessor<wo, na> fa,
  field<double>::accessor<wo, na> sa,
  field<double>::accessor<wo, na> Aua) noexcept;

void poisson_stencil(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<wo, na> soa) noexcept;

void turner_stencil(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, na> ud1,
  field<double>::accessor<ro, na> ud2,
  stencil_field<five_pt>::accessor<wo, na> soa,
  double dt) noexcept;

} // namespace gmg::task

#endif // GMG_TASKS_INIT_HH
