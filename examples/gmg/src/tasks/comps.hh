/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_TASK_COMPS_HH
#define GMG_TASK_COMPS_HH

#include "../types.hh"

namespace gmg::task {

/*!
  Restrict the fine-grid residual error to the coarse-grid RHS.
  @param mf  Fine-grid mesh topology instance.
  @param mc  Coarse-grid mesh topology instance.
  @param rfa Fine-grid residual error.
  @param fca Coarse-grid RHS.
 */
void full_weighting(flecsi::exec::cpu,
  mesh::accessor<ro> mf,
  mesh::accessor<ro> mc,
  field<double>::accessor<ro, ro> rfa,
  field<double>::accessor<rw, ro> fca) noexcept;

/*!
  Interpolate the coarse-grid approximate solution to the fine-grid error
  correction field.
  @param mc  Coarse-grid mesh topology instance.
  @param mf  Fine-grid mesh topology instance.
  @param uca Coarse-grid approximate solution.
  @param efa Fine-grid error correction.
 */
void bilinear_interpolation(flecsi::exec::cpu,
  mesh::accessor<ro> mc,
  mesh::accessor<ro> mf,
  field<double>::accessor<ro, ro> uca,
  field<double>::accessor<rw, ro> efa) noexcept;

void damped_jacobi(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<ro, na> soa,
  field<double>::accessor<rw, ro> u_new,
  field<double>::accessor<ro, ro> u_old,
  field<double>::accessor<ro, ro> fa,
  double omega) noexcept;

void red(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<ro, na> soa,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa) noexcept;

void black(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<ro, na> soa,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa) noexcept;

void residual(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<ro, na> soa,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<ro, ro> fa,
  field<double>::accessor<wo, ro> ra) noexcept;

void correction(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> ea) noexcept;

} // namespace gmg::task

#endif // GMG_TASK_COMPS_HH
