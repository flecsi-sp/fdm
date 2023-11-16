#ifndef GS_TASKS_SMOOTH_HH
#define GS_TASKS_SMOOTH_HH

#include "../types.hh"

namespace gs {
namespace task {

void red(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa);
void black(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa);

} // namespace task
} // namespace gs

#endif // GS_TASKS_SMOOTH_HH
