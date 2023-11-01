#ifndef POISSON_PROBLEM_HH
#define POISSON_PROBLEM_HH

#include "control.hh"
#include "initialize.hh"

namespace poisson {
namespace action {

void problem(control_policy &);
inline control::action<problem, cp::initialize> problem_action;
inline auto const problem_dep = problem_action.add(init_mesh_action);

} // namespace action
} // namespace poisson

#endif // POISSON_PROBLEM_HH
