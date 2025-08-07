/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_TASKS_IO_HH
#define GS_TASKS_IO_HH

#include "../types.hh"

#include <string>

namespace gs::task {

void io(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  std::string filebase) noexcept;

void print(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> fa) noexcept;

} // namespace gs::task

#endif // GS_TASKS_IO_HH
