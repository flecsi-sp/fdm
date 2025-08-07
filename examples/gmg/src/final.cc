#include "final.hh"
#include "state.hh"
#include "tasks/io.hh"

using namespace flecsi;

void
gmg::action::finalize(control_policy & cp) {
  auto & sc = cp.scheduler();
  auto & m = *mh[0];
  sc.execute<task::io>(exec::on, m, ud(m), "solver");
} // finalize
