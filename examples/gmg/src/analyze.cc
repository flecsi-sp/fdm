#include "analyze.hh"
#include "mg.hh"
#include "norm.hh"
#include "state.hh"
#include "tasks/analyze.hh"
#include "tasks/comps.hh"
#include "tasks/init.hh"

#include <algorithm>
#include <limits>

using namespace gmg;
using namespace flecsi;

void
action::analyze(control_policy & cp) {
  auto & sc = cp.scheduler();
  double err{std::numeric_limits<double>::max()};
  double prediction{std::numeric_limits<double>::max()};
  double difference{std::numeric_limits<double>::max()};
  std::size_t ita{0};

#if 0 // Test Jacobi
  std::size_t sub{100 > param::max_iterations ? param::max_iterations : 100};

  auto & m = *mh[0];
  double k = 2.0;
  double omega = 0.95;

  // Set field with a single mode
  sc.execute<task::fouriermodes>(exec::on, m, ud(m), k, k);

  // Set up the eigenvectors
  sc.execute<task::fouriermodes>(exec::on, m, sd(m), k, k);

  do {
    for(std::size_t i{0}; i < sub; ++i) {
      ud.flip();
      sc.execute<task::damped_jacobi>(
        exec::on, m, sod(m), ud(m), ud(m, 1), fd(m), omega);

      // Multiply by eigenvalue
      sc.execute<task::product_by_eigenvalue_jb>(
        exec::on, m, sd(m), omega, k, k);

      prediction = norm::errl2(sc);
      err = norm::l2(sc);
      difference = std::abs(prediction - err) / prediction;

      flog(info) << "Jacobi iteration, prediction, error, difference: "
                 << i + ita << ' ' << prediction << ' ' << err << ' '
                 << difference << std::endl;
    } // for
    ita += sub;

  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // Test Red-Black Gauss-Seidel
  std::size_t sub{100 > param::max_iterations ? param::max_iterations : 100};

  auto & m = *mh[0];
  double k = 1.0, l = 1.0;

  // Set field with a single mode
  sc.execute<task::gs_eigenvector>(exec::on, m, ud(m), k, l);

  // Set up the eigenvectors
  sc.execute<task::gs_eigenvector>(exec::on, m, sd(m), k, l);

  do {
    for(std::size_t i{0}; i < sub; ++i) {
      sc.execute<task::red>(exec::on, m, sod(m), ud(m), fd(m));
      sc.execute<task::black>(exec::on, m, sod(m), ud(m), fd(m));

      // Multiply by eigenvalue
      sc.execute<task::product_by_eigenvalue_gs>(exec::on, m, sd(m), k, l);

      prediction = norm::errl2(sc);
      err = norm::l2(sc);
      difference = std::abs(prediction - err) / prediction;

      flog(info) << "Gauss-Seidel iteration, prediction, error, difference: "
                 << i + ita << ' ' << prediction << ' ' << err << ' '
                 << difference << std::endl;
    } // for
    ita += sub;

  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // Test full weighting

  auto & mf = *mh[0];
  auto & mc = *mh[1];
  double k = 2.0;

  // Set field with a single mode
  sc.execute<task::fouriermodes>(exec::on, mf, ud(mf), k, k);
  // Set up the solution eigenvectors
  sc.execute<task::fourier_fw>(exec::on, mc, sd(mc), k, k);
  sc.execute<task::full_weighting>(exec::on, mf, mc, ud(mf), ud(mc));

  // Check solution
  err = norm::fwl2(sc);
  flog(info) << "FW error " << err << std::endl;

#endif

#if 0 // Test interpolation

  auto & mf = *mh[0];
  auto & mc = *mh[1];
  double k = 2.0;

  // Set field with a single mode
  sc.execute<task::fouriermodes>(exec::on, mc, ud(mc), k, k);
  // Set up the solution eigenvectors
  sc.execute<task::fourier_interp>(exec::on, mf, sd(mf), k, k);
  sc.execute<task::bilinear_interpolation>(exec::on, mc, mf, ud(mc), ud(mf));

  // Check solution
  err = norm::interpl2(sc);
  flog(info) << "Interpolation error " << err << std::endl;

#endif

#if 0 // Test residual
  auto & m = *mh[0];
  double k = 1.0;

  // Set field with a single mode
  sc.execute<task::fouriermodes>(exec::on, m, ud(m), k, k);

  // Set up the solution eigenvectors
  sc.execute<task::fourier_residual>(exec::on, m, sd(m), k, k);

  // Calculate residual
  sc.execute<task::residual>(exec::on, m, sod(m), ud(m), fd(m), rd(m));

  // Check solution
  err = norm::resl2(sc);
  flog(info) << "Residual error " << err << std::endl;
#endif

#if 0 // Test FMG

  auto & m = *mh[0];

  // Set the eggcarton problem
  sc.execute<task::eggcarton>(exec::on, m, ud(m), fd(m), sd(m), Aud(m));

  // Solve
  fmg(sc, param::fine_level);

  // Check difference with solution
  err = norm::errl2(sc);
  flog(info) << "Error " << err << std::endl;

#endif

} // analyze
