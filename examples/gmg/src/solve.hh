/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_SOLVE_HH
#define GMG_SOLVE_HH

#include "control.hh"

namespace gmg::action {

void solve(control_policy &);
inline control::action<solve, cp::solve> solve_action;

} // namespace gmg::action

#endif // GMG_SOLVE_HH
