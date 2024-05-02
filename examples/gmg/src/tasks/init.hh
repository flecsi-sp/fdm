#ifndef GMG_TASKS_INIT_HH
#define GMG_TASKS_INIT_HH

#include "../types.hh"

namespace gmg::task {

void enumerate(mesh::accessor<ro> m, field<double>::accessor<wo, na> fa);

/*!
  Initialize a field using a general plane equation, i.e., z = D + Mx + Ny.
 */
void bilinear(mesh::accessor<ro> m, field<double>::accessor<wo, na> fa,
    double M, double N, double D);

void constant(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double value);

void fouriermodes(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll);

void fourier_fw(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll);

void fourier_interp(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll);

void fourier_residual(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll);

void gs_eigenvector(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll);

void eggcarton(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  field<double>::accessor<wo, na> fa,
  field<double>::accessor<wo, na> sa,
  field<double>::accessor<wo, na> Aua);

void poisson_stencil(mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<wo, na> soa);

void turner_stencil(mesh::accessor<ro> m,
  field<double>::accessor<ro, na> ud1,
  field<double>::accessor<ro, na> ud2,
  stencil_field<five_pt>::accessor<wo, na> soa,
  double dt);

} // namespace gmg::task

#endif // GMG_TASKS_INIT_HH
