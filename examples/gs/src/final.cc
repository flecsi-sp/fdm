#include "final.hh"
#include "state.hh"
#include "tasks/io.hh"

using namespace flecsi;

void
gs::action::finalize(control_policy & cp) {
  execute<task::io, mpi>(cp.m, ud(cp.m), "solution");
} // finalize
