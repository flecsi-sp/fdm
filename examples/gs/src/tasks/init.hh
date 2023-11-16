#ifndef GS_TASKS_INIT_HH
#define GS_TASKS_INIT_HH

#include "../types.hh"

namespace gs::task {

void eggcarton(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  field<double>::accessor<wo, na> fa,
  field<double>::accessor<wo, na> sa,
  field<double>::accessor<wo, na> Aua);

void constant(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double value);

void redblack(mesh::accessor<ro> m, field<double>::accessor<wo, na> fa);

} // namespace gs::task

#endif // GS_TASKS_INIT_HH
