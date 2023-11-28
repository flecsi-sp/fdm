#include "solve.hh"
#include "state.hh"
#include "tasks/comps.hh"

#include <flecsi/execution.hh>

using namespace gmg;
using namespace flecsi;

void
action::solve(control_policy &) {
  auto & m = *mh[0].get();
  auto const & [ud_new, ud_old] = ud.get();
  execute<task::damped_jacobi>(m, ud_new(m), ud_old(m), fd(m), 1.0);
  ud.flip();
} // solve
