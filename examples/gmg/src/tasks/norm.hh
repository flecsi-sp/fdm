/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_TASKS_NORM_HH
#define GMG_TASKS_NORM_HH

#include "../types.hh"

namespace gmg::task {

double diff_sum_square(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba);

double diff_max(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba);

double scale(mesh::accessor<ro> m, double sum);

void discrete_operator(mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<ro, na> soa,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<wo, ro> Aua);

} // namespace gmg::task

#endif // GMG_TASKS_NORM_HH
