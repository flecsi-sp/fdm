#include "analyze.hh"
#include "annotation.hh"
#include "control.hh"
#include "final.hh"
#include "init.hh"
#include "options.hh"
#include "problem.hh"
#include "solve.hh"
#include "state.hh"

#include <flecsi/runtime.hh>

using namespace flecsi;

int
main(int argc, char ** argv) {
  const flecsi::getopt g;
  try {
    g(argc, argv);
  }
  catch(const std::logic_error & e) {
    std::cerr << e.what() << '\n' << g.usage(argc ? argv[0] : "");
    return 1;
  }
  const flecsi::run::dependencies_guard dg;
  flog::add_output_stream("clog", std::clog, true);
  return flecsi::runtime().control<gs::control>();
} // main
