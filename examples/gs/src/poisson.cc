#include "analyze.hh"
#include "annotation.hh"
#include "control.hh"
#include "finalize.hh"
#include "initialize.hh"
#include "options.hh"
#include "problem.hh"
#include "solve.hh"
#include "state.hh"

#include <flecsi/runtime.hh>

using namespace flecsi;

int
main(int argc, char ** argv) {
  run::arguments args(argc, argv);
  const run::dependencies_guard dg(args.dep);
  const runtime run(args.cfg);
  flog::add_output_stream("flog", std::clog, true);
  return run.main<poisson::control>(args.act);
} // main
