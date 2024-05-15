/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_INIT_HH
#define GMG_INIT_HH

#include "control.hh"

namespace gmg::action {

void init(control_policy &);
inline control::action<init, cp::initialize> init_action;

} // namespace gmg::action

#endif // GMG_INIT_HH
