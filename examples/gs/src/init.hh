/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_INIT_HH
#define GS_INIT_HH

#include "control.hh"

namespace gs::action {

void init_mesh(control_policy &);
inline control::action<init_mesh, cp::initialize> init_mesh_action;

} // namespace gs::action

#endif // GS_INIT_HH
