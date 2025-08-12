/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_TASK_ANALYZE_HH
#define GMG_TASK_ANALYZE_HH

#include "../types.hh"

namespace gmg::task {

void product_by_eigenvalue_jb(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> sa,
  double omega,
  double k,
  double l) noexcept;

void product_by_eigenvalue_gs(flecsi::exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> sa,
  double k,
  double l) noexcept;

} // namespace gmg::task

#endif // GMG_TASK_ANALYZE_HH
