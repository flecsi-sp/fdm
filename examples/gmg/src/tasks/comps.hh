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
void full_weighting(mesh::accessor<ro> mf,
  mesh::accessor<ro> mc,
  field<double>::accessor<ro, ro> rfa,
  field<double>::accessor<rw, ro> fca);

/*!
  Interpolate the coarse-grid approximate solution to the fine-grid error
  correction field.
  @param mc  Coarse-grid mesh topology instance.
  @param mf  Fine-grid mesh topology instance.
  @param uca Coarse-grid approximate solution.
  @param efa Fine-grid error correction.
 */
void bilinear_interpolation(mesh::accessor<ro> mc,
  mesh::accessor<ro> mf,
  field<double>::accessor<ro, ro> uca,
  field<double>::accessor<rw, ro> efa);

void damped_jacobi(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> u_new,
  field<double>::accessor<ro, ro> u_old,
  field<double>::accessor<ro, ro> fa,
  double omega);

void red(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa);

void black(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa);

void residual(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<ro, ro> fa,
  field<double>::accessor<rw, ro> ra);

void correction(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> ea);

} // namespace gmg::task

#endif // GMG_TASK_COMPS_HH
