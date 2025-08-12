/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_NORM_HH
#define GMG_NORM_HH

#include "state.hh"
#include "tasks/norm.hh"

namespace gmg::norm {

double l2(flecsi::scheduler & sc);
double fwl2(flecsi::scheduler & sc);
double interpl2(flecsi::scheduler & sc);
double resl2(flecsi::scheduler & sc);
double errl2(flecsi::scheduler & sc);
double max(flecsi::scheduler & sc);

} // namespace gmg::norm

#endif // GMG_NORM_HH
