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
  auto ur = ud(m);
  auto vr = vd(m);
  do {
    for(std::size_t i{0}; i < sub; ++i) {
      std::swap(ur, vr);
      execute<task::damped_jacobi>(m, ur, vr, fd(m), 0.8);
    } // for
    ita += sub;

    err = norm::l2();
    flog(info) << "residual: " << err << " (" << ita << " iterations)"
               << std::endl;
    flog(info) << "max: " << norm::max() << " (" << ita << " iterations)"
               << std::endl;

  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 1 // MG Components
  auto & mf = *mh[0].get();
  auto & mc = *mh[1].get();

  execute<task::full_weighting>(mf, mc, rd(mf), fd(mc));
  execute<task::print>(mc, fd(mc));
#endif

#if 0 // Two-Grid Method
  std::size_t pre{5};
  std::size_t post{5};

  auto & mf = *mh[0].get();
  auto & mc = *mh[1].get();
  auto urf = ud(mf);
  auto urc = ud(mc);
  auto vrf = vd(mf);
  auto vrc = vd(mc);

  do {
    // Pre Smoothing
    for(std::size_t i{0}; i < pre; ++i) {
      std::swap(urf, vrf);
      execute<task::damped_jacobi>(mf, urf, vrf, fd(mf), 0.8);
    } // for

    execute<task::residual>(mf, urf, fd(mf), rd(mf));
    execute<task::full_weighting>(mf, mc, rd(mf), fd(mc));
    execute<task::constant>(mc, urc, 0.0);

    // "Solve" on coarse grid
    for(std::size_t i{0}; i < 500; ++i) {
      std::swap(urc, vrc);
      execute<task::damped_jacobi>(mc, urc, vrc, fd(mc), 0.8);
    } // for

    execute<task::bilinear_interpolation>(mc, mf, urc, ed(mf));
    execute<task::correction>(mf, urf, ed(mf));

    // Post Smoothing
    for(std::size_t i{0}; i < post; ++i) {
      std::swap(urf, vrf);
      execute<task::damped_jacobi>(mf, urf, vrf, fd(mf), 0.8);
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
