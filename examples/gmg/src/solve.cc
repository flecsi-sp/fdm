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
action::solve(control_policy & cp) {
  auto & sc = cp.scheduler();
  double err{std::numeric_limits<double>::max()};
  std::size_t ita{0};

#if 0 // Jacobi
  std::size_t sub{100 > param::max_iterations ? param::max_iterations : 100};

  auto & m = *mh[0];
  sc.execute<task::constant>(exec::on, m, ud(m), 0.0);

  do {
    for(std::size_t i{0}; i < sub; ++i) {
      ud.flip();
      sc.execute<task::damped_jacobi>(
        exec::on, m, sod(m), ud(m), ud(m, 1), fd(m), 0.8);
    } // for
    ita += sub;

    err = norm::l2(sc);
    flog(info) << "residual: " << err << " (" << ita << " iterations)"
               << std::endl;
    flog(info) << "max: " << norm::max(sc) << " (" << ita << " iterations)"
               << std::endl;

  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // Red-Black Gauss Seidel
  std::size_t sub{100 > param::max_iterations ? param::max_iterations : 100};

  auto & m = *mh[0];

  do {
    for(std::size_t i{0}; i < sub; ++i) {
      sc.execute<task::red>(exec::on, m, sod(m), ud(m), fd(m));
      sc.execute<task::black>(exec::on, m, sod(m), ud(m), fd(m));
    } // for
    ita += sub;

    err = norm::l2(sc);
    flog(info) << "residual: " << err << " (" << ita << " iterations)"
               << std::endl;
    flog(info) << "max: " << norm::max(sc) << " (" << ita << " iterations)"
               << std::endl;

  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // Grid Transfer
  auto & mf = *mh[0];
  auto & mc = *mh[1];
  sc.execute<task::residual>(exec::on, mf, sod(mf), ud(mf), fd(mf), rd(mf));
  sc.execute<task::print>(exec::on, mf, ud(mf));
  sc.execute<task::full_weighting>(exec::on, mf, mc, rd(mf), fd(mc));
  sc.execute<task::io>(exec::on, mc, fd(mc), "fw");
  sc.execute<task::bilinear_interpolation>(exec::on, mc, mf, ud(mc), ed(mf));
  sc.execute<task::io>(exec::on, mf, ed(mf), "bl");
#endif

#if 0 // Two-Grid Method
  std::size_t pre{5};
  std::size_t post{5};

  auto & mf = mh[0];
  auto & mc = mh[1];
  sc.execute<task::constant>(exec::on, *mf, ud(*mf), 0.0);

  do {
    // Pre Smoothing
    for(std::size_t i{0}; i < pre; ++i) {
      ud.flip();
      sc.execute<task::damped_jacobi>(
        exec::on, *mf, sod(*mf), ud(*mf), ud(*mf, 1), fd(*mf), 0.8);
    } // for

    sc.execute<task::residual>(
      exec::on, *mf, sod(*mf), ud(*mf), fd(*mf), rd(*mf));
    sc.execute<task::full_weighting>(exec::on, *mf, *mc, rd(*mf), fd(*mc));
    sc.execute<task::constant>(exec::on, *mc, ud(*mc), 0.0);
    execute<task::constant>(exec::on, *mc, ud(*mc, 1), 0.0);

    // "Solve" on coarse grid
    for(std::size_t i{0}; i < 500; ++i) {
      ud.flip();
      sc.execute<task::damped_jacobi>(
        exec::on, *mc, sod(*mc), ud(*mc), ud(*mc, 1), fd(*mc), 0.8);
    } // for

    sc.execute<task::bilinear_interpolation>(
      exec::on, *mc, *mf, ud(*mc), ed(*mf));
    sc.execute<task::correction>(exec::on, *mf, ud(*mf), ed(*mf));

    // Post Smoothing
    for(std::size_t i{0}; i < post; ++i) {
      ud.flip();
      sc.execute<task::damped_jacobi>(
        exec::on, *mf, sod(*mf), ud(*mf), ud(*mf, 1), fd(*mf), 0.8);
    } // for

    err = norm::l2(sc);
    flog(info) << "residual: " << err << " (" << ita << " iterations)"
               << std::endl;
    flog(info) << "max: " << norm::max(sc) << " (" << ita << " iterations)"
               << std::endl;

    ++ita;
  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // V-Cycle
  do {
    vcycle(sc, param::fine_level);

    err = norm::l2(sc);
    flog(info) << "residual: " << err << " (" << ita << " iterations)"
               << std::endl;
    flog(info) << "max: " << norm::max(sc) << " (" << ita << " iterations)"
               << std::endl;

    ++ita;
  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // FMG

  fmg(sc, param::fine_level);

  err = norm::l2(sc);
  flog(info) << "residual: " << err << std::endl;
  flog(info) << "max: " << norm::max(sc) << std::endl;
#endif

} // solve
