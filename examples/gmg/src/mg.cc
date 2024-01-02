#include "mg.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/init.hh"
#include "util.hh"

using namespace flecsi;

void
gmg::vcycle(std::size_t level) {
  auto & mf = *mh[util::index(level)];

  if(level == param::vcycle_direct) {
    flog(warn) << "Direct solve level(index): " << level << "("
               << util::index(level) << ")" << std::endl;
    // "Solve"
    for(std::size_t i{0}; i < 1000; ++i) {
      execute<task::damped_jacobi>(mf, ud(mf), ud(mf,1), fd(mf), 0.8);
      ud.flip();
    } // for
  }
  else {
    flog(warn) << "Cycle level(index): " << level << "(" << util::index(level)
               << ")" << std::endl;
    auto & mc = *mh[util::index(level - 1)];

    // Pre Smoothing
    for(std::size_t i{0}; i < param::vcycle_pre; ++i) {
      execute<task::damped_jacobi>(mf, ud(mf), ud(mf,1), fd(mf), 0.8);
      ud.flip();
    } // for

    // Recursive solve
    execute<task::residual>(mf, ud(mf), fd(mf), rd(mf));
    execute<task::full_weighting>(mf, mc, rd(mf), fd(mc));
    execute<task::constant>(mc, ud(mc), 0.0);

    gmg::vcycle(level - 1);

    execute<task::bilinear_interpolation>(mc, mf, ud(mc), ed(mf));
    execute<task::correction>(mf, ud(mf), ed(mf));

    // Post Smoothing
    for(std::size_t i{0}; i < param::vcycle_post; ++i) {
      execute<task::damped_jacobi>(mf, ud(mf), ud(mf,1), fd(mf), 0.8);
      ud.flip();
    } // for
  } // if
} // vcycle
