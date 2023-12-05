#include "final.hh"
#include "state.hh"
#include "tasks/io.hh"

using namespace flecsi;

void
gmg::action::finalize(control_policy & cp) {
  auto & m = *mh[0].get();
  execute<task::io, mpi>(m, ud(m), "solver");
} // finalize
