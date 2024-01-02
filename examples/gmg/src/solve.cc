#include "solve.hh"
#include "mg.hh"
#include "norm.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/init.hh"
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

#if 0 // Jacobi
  std::size_t sub{100 > param::max_iterations ? param::max_iterations : 100};

  auto & m = *mh[0].get();
  execute<task::constant>(m, ud(m,1), 0.0);
  ud.flip();

  do {
    for(std::size_t i{0}; i < sub; ++i) {
      execute<task::damped_jacobi>(m, ud(m), ud(m,1), fd(m), 0.8);
      ud.flip();
    } // for
    ita += sub;

    err = norm::l2();
    flog(info) << "residual: " << err << " (" << ita << " iterations)"
               << std::endl;
    flog(info) << "max: " << norm::max() << " (" << ita << " iterations)"
               << std::endl;

  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // Two-Grid Method
  std::size_t pre{5};
  std::size_t post{5};

  auto & mf = *mh[0].get();
  auto & mc = *mh[1].get();
  execute<task::constant>(mf, ud(mf,1), 0.0);
  ud.flip();

  do {
    // Pre Smoothing
    for(std::size_t i{0}; i < pre; ++i) {
      execute<task::damped_jacobi>(mf, ud(mf), ud(mf,1), fd(mf), 0.8);
      ud.flip();
    } // for

    execute<task::residual>(mf, ud(mf), fd(mf), rd(mf));
    execute<task::full_weighting>(mf, mc, rd(mf), fd(mc));
    execute<task::constant>(mc, ud(mc), 0.0);
    execute<task::constant>(mc, ud(mc,1), 0.0);

    // "Solve" on coarse grid
    for(std::size_t i{0}; i < 500; ++i) {
      execute<task::damped_jacobi>(mc, ud(mc), ud(mc,1), fd(mc), 0.8);
      ud.flip();
    } // for

    execute<task::bilinear_interpolation>(mc, mf, ud(mc), ed(mf));
    execute<task::correction>(mf, ud(mf), ed(mf));

    // Post Smoothing
    for(std::size_t i{0}; i < post; ++i) {
      execute<task::damped_jacobi>(mf, ud(mf), ud(mf,1), fd(mf), 0.8);
      ud.flip();
    } // for

    err = norm::l2();
    flog(info) << "residual: " << err << " (" << ita << " iterations)"
               << std::endl;
    flog(info) << "max: " << norm::max() << " (" << ita << " iterations)"
               << std::endl;

    ++ita;
  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // V-Cycle
  do {
    vcycle(param::fine_level);

    err = norm::l2();
    flog(info) << "residual: " << err << " (" << ita << " iterations)"
               << std::endl;
    flog(info) << "max: " << norm::max() << " (" << ita << " iterations)"
               << std::endl;

    ++ita;
  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

} // solve
