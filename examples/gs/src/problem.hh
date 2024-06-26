/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_PROBLEM_HH
#define GS_PROBLEM_HH

#include "control.hh"
#include "init.hh"

namespace gs {
namespace action {

void problem(control_policy &);
inline control::action<problem, cp::initialize> problem_action;
inline auto const problem_dep = problem_action.add(init_mesh_action);

} // namespace action
} // namespace gs

#endif // GS_PROBLEM_HH
