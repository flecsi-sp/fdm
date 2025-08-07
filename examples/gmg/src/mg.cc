#include "mg.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/init.hh"
#include "util.hh"

using namespace flecsi;

void
gmg::vcycle(flecsi::scheduler & sc, std::size_t level) {
  auto & mf = *mh[util::index(level)];

  if(level == param::mg_direct) {
    flog(warn) << "Direct solve level(index): " << level << "("
               << util::index(level) << ")" << std::endl;
    // FIXME: "Solve"
    for(std::size_t i{0}; i < 1000; ++i) {
      ud.flip();
      sc.execute<task::damped_jacobi>(
        exec::on, mf, sod(mf), ud(mf), ud(mf, 1), fd(mf), 0.8);
    } // for
  }
  else {
    flog(warn) << "Cycle level(index): " << level << "(" << util::index(level)
               << ")" << std::endl;
    auto & mc = *mh[util::index(level - 1)];

    // Pre Smoothing
    for(std::size_t i{0}; i < param::mg_pre; ++i) {
      ud.flip();
      sc.execute<task::damped_jacobi>(
        exec::on, mf, sod(mf), ud(mf), ud(mf, 1), fd(mf), 0.8);
    } // for

    // Recursive solve
    sc.execute<task::residual>(exec::on, mf, sod(mf), ud(mf), fd(mf), rd(mf));
    sc.execute<task::full_weighting>(exec::on, mf, mc, rd(mf), fd(mc));
    sc.execute<task::constant>(exec::on, mc, ud(mc), 0.0);

    gmg::vcycle(sc, level - 1);

    sc.execute<task::bilinear_interpolation>(exec::on, mc, mf, ud(mc), ed(mf));
    sc.execute<task::correction>(exec::on, mf, ud(mf), ed(mf));

    // Post Smoothing
    for(std::size_t i{0}; i < param::mg_post; ++i) {
      ud.flip();
      sc.execute<task::damped_jacobi>(
        exec::on, mf, sod(mf), ud(mf), ud(mf, 1), fd(mf), 0.8);
    } // for
  } // if
} // vcycle

void
gmg::fmg(flecsi::scheduler & sc, std::size_t level) {
  auto & mf = mh[util::index(level)];

  // The scheme requires:
  // 1) Go to a coarser grid, adapt all and repeat this step
  // 2) If in the deeper level, direct solve or do a V-Cycle for a number of
  // iterations 3) Come back up, interpolate, and do a V-Cycle

  // Deepest level
  if(level == param::mg_direct) {
    flog(warn) << "Deepest level(index):" << level << "(" << util::index(level)
               << ")" << std::endl;

    // If in the deepest level, the V-Cycle is already doing a direct solve
    gmg::vcycle(sc, level);
  }
  else {
    flog(warn) << "Cycle level(index): " << level << "(" << util::index(level)
               << ")" << std::endl;
    auto & mc = mh[util::index(level - 1)];

    // Set the RHS and solution field
    sc.execute<task::full_weighting>(exec::on, *mf, *mc, fd(*mf), fd(*mc));
    sc.execute<task::full_weighting>(exec::on, *mf, *mc, ud(*mf), ud(*mc));

    // Now call solve for one level deeper
    gmg::fmg(sc, level - 1);

    // Interpolate solution back up (RHS does not change)
    sc.execute<task::bilinear_interpolation>(
      exec::on, *mc, *mf, ud(*mc), ud(*mf));

    // Do a V-Cycle
    for(std::size_t i{0}; i < param::mg_cycles; ++i) {
      gmg::vcycle(sc, level);
    } // for
  } // if
} // fmg
