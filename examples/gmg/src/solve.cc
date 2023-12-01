#include "solve.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/io.hh"

#include <flecsi/execution.hh>

#include <limits>

using namespace gmg;
using namespace flecsi;

void
action::solve(control_policy &) {
#if 0
  for(auto const & m : mh) {
    execute<task::print>(*m, fd(*m));
  } // for
#endif

#if 0
  auto & m0 = *mh[0].get();
  auto & m1 = *mh[1].get();
  // execute<task::full_weighting>(*m0, *m1, ud(*m0), ud(*m1));
  // execute<task::print>(*m1, ud(*m1));
  execute<task::bilinear_interpolation>(m1, m0, ud_new(m1), ud_new(m0));
  execute<task::print>(m0, ud_new(m0));
#endif

#if 0
  auto & m = *mh[0].get();
  for(std::size_t i{0}; i<10000; ++i) {
    execute<task::damped_jacobi>(m, ud[0](m), ud[1](m), fd(m), 1.0);
    ud.flip();
  } // for

  execute<task::print>(m, ud[0](m));
  execute<task::print>(m, ud[1](m));
#endif

  double err{std::numeric_limits<double>::max()};
  std::size_t ita{0};

  auto & m = *mh[0].get();
  do {
    execute<task::damped_jacobi>(m, ud[0](m), ud[1](m), fd(m), 0.8);
    ud.flip();
  } while(err > opt::tolerance && ita < opt::max_iterations);
} // solve
