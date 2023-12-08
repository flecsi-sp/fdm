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

void eggcarton(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  field<double>::accessor<wo, na> fa,
  field<double>::accessor<wo, na> sa,
  field<double>::accessor<wo, na> Aua);

} // namespace gmg::task

#endif // GMG_TASKS_INIT_HH
