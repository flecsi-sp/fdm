#include "norm.hh"
#include "solve.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/io.hh"
#include "tasks/norm.hh"

#include <flecsi/execution.hh>

#include <algorithm>
#include <limits>

using namespace gmg;
using namespace flecsi;

void
action::solve(control_policy &) {
  double err{std::numeric_limits<double>::max()};
  std::size_t ita{0};

#if 1 // Jacobi
  std::size_t sub{100 > opt::max_iterations ? opt::max_iterations : 100};

  auto & m = *mh[0].get();
  do {
    for(std::size_t i{0}; i<sub; ++i) {
      execute<task::damped_jacobi>(m, ud[0](m), ud[1](m), fd(m), 0.8);
      ud.flip();
    } // for
    ita += sub;

    flog(info) << "residual: " << norm::l2() << " (" << ita << " iterations)" << std::endl;
    flog(info) << "max: " << norm::max() << " (" << ita << " iterations)" << std::endl;

  } while(err > opt::error_tolerance && ita < opt::max_iterations);
#endif

#if 0 // Two-Grid Method
  std::size_t pre{5};
  std::size_t post{5};

  auto & mf = *mh[0].get();
  auto & mc = *mh[1].get();

  do {
    // Pre Smoothing
    for(std::size_t i{0}; i<pre; ++i) {
      execute<task::damped_jacobi>(mf, ud[0](mf), ud[1](mf), fd(mf), 0.8);
      ud.flip();
    } // for

    execute<task::residual>(mf, ud[1](mf), fd(mf), rd(mf));
    execute<task::full_weighting>(mf, mc, rd(mf), fd(mc));

    // "Solve" on coarse grid
    for(std::size_t i{0}; i<500; ++i) {
      execute<task::damped_jacobi>(mc, ud[0](mc), ud[1](mc), fd(mc), 0.8);
      ud.flip();
    } // for

    execute<task::bilinear_interpolation>(mc, mf, ud[1](mc), ed(mf));
    execute<task::correction>(mf, ud[1](mf), ed(mf));

    // Post Smoothing
    for(std::size_t i{0}; i<post; ++i) {
      execute<task::damped_jacobi>(mf, ud[0](mf), ud[1](mf), fd(mf), 0.8);
      ud.flip();
    } // for

    flog(info) << "residual: " << norm::l2() << " (" << ita << " iterations)" << std::endl;

    ++ita;
  } while(err > opt::error_tolerance && ita < opt::max_iterations);
#endif

} // solve
