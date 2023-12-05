#include "mg.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/init.hh"

using namespace flecsi;

void gmg::vcycle(std::size_t level) {
  auto & mf = *mh[level];
  auto urf = ud(mf);
  auto vrf = vd(mf);

  if(level == opt::level_direct) {
    // "Solve"
    for(std::size_t i{0}; i<1000; ++i) {
      std::swap(urf, vrf);
      execute<task::damped_jacobi>(mf, urf, vrf, fd(mf), 0.8);
    } // for
  }
  else {
    auto & mc = *mh[level+1];

    auto urc = vd(mc);

    // Pre Smoothing
    for(std::size_t i{0}; i<opt::vcycle_pre; ++i) {
      std::swap(urf, vrf);
      execute<task::damped_jacobi>(mf, urf, vrf, fd(mf), 0.8);
    } // for

    // Recursive solve
    execute<task::residual>(mf, urf, fd(mf), rd(mf));
    execute<task::full_weighting>(mf, mc, rd(mf), fd(mc));
    execute<task::constant>(mc, urc, 0.0);

    gmg::vcycle(level+1);

    execute<task::bilinear_interpolation>(mc, mf, urc, ed(mf));
    execute<task::correction>(mf, urf, ed(mf));

    // Post Smoothing
    for(std::size_t i{0}; i<opt::vcycle_post; ++i) {
      std::swap(urf, vrf);
      execute<task::damped_jacobi>(mf, urf, vrf, fd(mf), 0.8);
    } // for
  } // if
} // vcycle
