/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_TASKS_NORM_HH
#define GS_TASKS_NORM_HH

#include "../types.hh"

namespace gs {
namespace task {

double diff(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba) noexcept;

double scale(flecsi::exec::cpu, mesh::accessor<ro> m, double sum) noexcept;

void display_l2(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  flecsi::future<double> scaled) noexcept;

void discrete_operator(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<rw, ro> Aua) noexcept;

} // namespace task
} // namespace gs

#endif // GS_TASKS_NORM_HH
