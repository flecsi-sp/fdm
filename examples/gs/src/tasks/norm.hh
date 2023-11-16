#ifndef GS_TASKS_NORM_HH
#define GS_TASKS_NORM_HH

#include "../types.hh"

namespace gs {
namespace task {

double diff(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba);

double scale(mesh::accessor<ro> m, double sum);

void discrete_operator(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<rw, ro> Aua);

} // namespace task
} // namespace gs

#endif // GS_TASKS_NORM_HH
