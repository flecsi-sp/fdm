/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_TASKS_SMOOTH_HH
#define GS_TASKS_SMOOTH_HH

#include "../types.hh"

namespace gs {
namespace task {

void red(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa) noexcept;
void black(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa) noexcept;

} // namespace task
} // namespace gs

#endif // GS_TASKS_SMOOTH_HH
