#ifndef GS_TASKS_IO_HH
#define GS_TASKS_IO_HH

#include "../types.hh"

#include <string>

namespace gs::task {

void io(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  std::string filebase);

void print(mesh::accessor<ro> m, field<double>::accessor<ro, ro> fa);

} // namespace gs::task

#endif // GS_TASKS_IO_HH
