#ifndef GS_SOLVE_HH
#define GS_SOLVE_HH

#include "control.hh"

namespace gs {
namespace action {

void solve(control_policy &);
inline control::action<solve, cp::solve> solve_action;

} // namespace action
} // namespace gs

#endif // GS_SOLVE_HH
