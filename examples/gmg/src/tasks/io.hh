#ifndef GMG_TASKS_IO_HH
#define GMG_TASKS_IO_HH

#include "../types.hh"

namespace gmg::task {

void print(mesh::accessor<ro> m, field<double>::accessor<ro, ro> fa);

} // namespace gmg::task

#endif // GMG_TASKS_IO_HH
