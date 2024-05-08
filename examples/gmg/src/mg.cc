#include "mg.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/init.hh"
#include "util.hh"

using namespace flecsi;

void
gmg::vcycle(std::size_t level) {
  auto & mf = *mh[util::index(level)];

  if(level == param::mg_direct) {
    flog(warn) << "Direct solve level(index): " << level << "("
               << util::index(level) << ")" << std::endl;
    // FIXME "Solve"
    for(std::size_t i{0}; i < 1000; ++i) {
      ud.flip();
      execute<task::damped_jacobi>(mf, sod(mf), ud(mf), ud(mf, 1), fd(mf), 0.8);
    } // for
  }
  else {
    flog(warn) << "Cycle level(index): " << level << "(" << util::index(level)
               << ")" << std::endl;
    auto & mc = *mh[util::index(level - 1)];

    // Pre Smoothing
    for(std::size_t i{0}; i < param::mg_pre; ++i) {
      ud.flip();
      execute<task::damped_jacobi>(mf, sod(mf), ud(mf), ud(mf, 1), fd(mf), 0.8);
    } // for

    // Recursive solve
    execute<task::residual>(mf, sod(mf), ud(mf), fd(mf), rd(mf));
    execute<task::full_weighting>(mf, mc, rd(mf), fd(mc));
    execute<task::constant>(mc, ud(mc), 0.0);

    gmg::vcycle(level - 1);

    execute<task::bilinear_interpolation>(mc, mf, ud(mc), ed(mf));
    execute<task::correction>(mf, ud(mf), ed(mf));

    // Post Smoothing
    for(std::size_t i{0}; i < param::mg_post; ++i) {
      ud.flip();
      execute<task::damped_jacobi>(mf, sod(mf), ud(mf), ud(mf, 1), fd(mf), 0.8);
    } // for
  } // if
} // vcycle

void
gmg::fmg(std::size_t level) {
  auto & mf = *mh[util::index(level)];

  // The scheme requires:
  // 1) Go to a coarser grid, adapt all and repeat this step
  // 2) If in the deeper level, direct solve or do a V-Cycle for a number of
  // iterations 3) Come back up, interpolate, and do a V-Cycle

  // Deepest level
  if(level == param::mg_direct) {
    flog(warn) << "Deepest level(index):" << level << "(" << util::index(level)
               << ")" << std::endl;

    // If in the deepest level, the V-Cycle is already doing a direct solve
    gmg::vcycle(level);
  }
  else {
    flog(warn) << "Cycle level(index): " << level << "(" << util::index(level)
               << ")" << std::endl;
    auto & mc = *mh[util::index(level - 1)];

    // Set the RHS and solution field
    execute<task::full_weighting>(mf, mc, fd(mf), fd(mc));
    execute<task::full_weighting>(mf, mc, ud(mf), ud(mc));

    // Now call solve for one level deeper
    gmg::fmg(level - 1);

    // Interpolate solution back up (RHS does not change)
    execute<task::bilinear_interpolation>(mc, mf, ud(mc), ud(mf));

    // Do a V-Cycle
    for(std::size_t i{0}; i < param::mg_cycles; ++i) {
      gmg::vcycle(level);
    } // for
  } // if
} // fmg
