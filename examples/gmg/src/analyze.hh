/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_ANALYZE_HH
#define GMG_ANALYZE_HH

#include "control.hh"

namespace gmg::action {

void analyze(control_policy &);
inline control::action<analyze, cp::analyze> analyze_action;

} // namespace gmg::action

#endif // GMG_ANALYZE_HH
