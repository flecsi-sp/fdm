#ifndef GMG_TASKS_NORM_HH
#define GMG_TASKS_NORM_HH

#include "../types.hh"

namespace gmg::task {

double diff(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba);

double scale(mesh::accessor<ro> m, double sum);

void discrete_operator(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<rw, ro> Aua);

} // namespace gmg::task

#endif // GMG_TASKS_NORM_HH
