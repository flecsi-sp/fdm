/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_TASKS_IO_HH
#define GMG_TASKS_IO_HH

#include "../types.hh"

namespace gmg::task {

void io(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> fa,
  std::string filebase) noexcept;

void print(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> fa) noexcept;

} // namespace gmg::task

#endif // GMG_TASKS_IO_HH
