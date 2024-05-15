/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_NORM_HH
#define GMG_NORM_HH

#include "state.hh"
#include "tasks/norm.hh"

namespace gmg::norm {

double l2(/* std::size_t level */);
double fwl2(/* std::size_t level */);
double interpl2(/* std::size_t level */);
double resl2(/* std::size_t level */);
double errl2(/* std::size_t level */);
double max(/* std::size_t level */);

} // namespace gmg::norm

#endif // GMG_NORM_HH
