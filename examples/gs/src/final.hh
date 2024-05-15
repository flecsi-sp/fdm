/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_FINAL_HH
#define GS_FINAL_HH

#include "control.hh"

namespace gs::action {

void finalize(control_policy &);
inline control::action<finalize, cp::finalize> finalize_action;

} // namespace gs::action

#endif // GS_FINAL_HH
