#include "control.hh"
#include "final.hh"
#include "init.hh"
#include "options.hh"
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
  return run.main<gmg::control>(args.act);
} // main
