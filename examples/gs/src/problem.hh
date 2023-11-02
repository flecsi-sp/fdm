#ifndef GS_PROBLEM_HH
#define GS_PROBLEM_HH

#include "control.hh"
#include "init.hh"

namespace poisson {
namespace action {

void problem(control_policy &);
inline control::action<problem, cp::initialize> problem_action;
inline auto const problem_dep = problem_action.add(init_mesh_action);

} // namespace action
} // namespace poisson

#endif // GS_PROBLEM_HH
