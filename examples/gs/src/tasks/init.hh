/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_TASKS_INIT_HH
#define GS_TASKS_INIT_HH

#include "../types.hh"

namespace gs::task {

void eggcarton(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  field<double>::accessor<wo, na> fa,
  field<double>::accessor<wo, na> sa,
  field<double>::accessor<wo, na> Aua) noexcept;

void constant(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double value) noexcept;

void redblack(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa) noexcept;

} // namespace gs::task

#endif // GS_TASKS_INIT_HH
