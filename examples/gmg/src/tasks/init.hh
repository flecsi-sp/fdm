#ifndef GMG_TASKS_INIT_HH
#define GMG_TASKS_INIT_HH

#include "../types.hh"

namespace gmg::task {

void constant(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double value);

void enumerate(mesh::accessor<ro> m, field<double>::accessor<wo, na> fa);

} // namespace gmg::task

#endif // GMG_TASKS_INIT_HH
