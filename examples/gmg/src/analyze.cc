#include "analyze.hh"
#include "norm.hh"
#include "tasks/comps.hh"
#include "tasks/analyze.hh"
#include "tasks/init.hh"
#include "mg.hh"
#include "state.hh"

#include <limits>
#include <algorithm>

using namespace gmg;
using namespace flecsi;

void
action::analyze(control_policy &) {
  double err{std::numeric_limits<double>::max()};
  double prediction{std::numeric_limits<double>::max()};
  double difference{std::numeric_limits<double>::max()};
  std::size_t ita{0};

#if 0 // Test Jacobi
  std::size_t sub{100 > param::max_iterations ? param::max_iterations : 100};

  auto & m = *mh[0].get();
  double k = 2.0;
  double omega = 0.95;

  // Set field with a single mode
  execute<task::fouriermodes>(m, ud(m), k, k);

  // Set up the eigenvectors
  execute<task::fouriermodes>(m, sd(m), k, k);

  do {
    for(std::size_t i{0}; i < sub; ++i) {
      ud.flip();
      execute<task::damped_jacobi>(m, ud(m), ud(m, 1), fd(m), omega);

      // Multiply by eigenvalue
      execute<task::product_by_eigenvalue_jb>(m, sd(m), omega, k, k);

      prediction = norm::errl2();
      err = norm::l2();
      difference = abs(prediction - err) / prediction;

      flog(info) << "Jacobi iteration, prediction, error, difference: " << i + ita << ' ' << prediction << ' ' << err << ' ' << difference << std::endl;
    } // for
    ita += sub;

  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 1 // Test Red-Black Gauss-Seidel
  std::size_t sub{100 > param::max_iterations ? param::max_iterations : 100};

  auto & m = *mh[0].get();
  double k = 1.0, l = 1.0;

  // Set field with a single mode
  execute<task::gs_eigenvector>(m, ud(m), k, l);

  // Set up the eigenvectors
  execute<task::gs_eigenvector>(m, sd(m), k, l);

  do {
    for(std::size_t i{0}; i < sub; ++i) {
      execute<task::red>(m, ud(m), fd(m));
      execute<task::black>(m, ud(m), fd(m));

      // Multiply by eigenvalue
      execute<task::product_by_eigenvalue_gs>(m, sd(m), k, l);

      prediction = norm::errl2();
      err = norm::l2();
      difference = abs(prediction - err) / prediction;

      flog(info) << "Gauss-Seidel iteration, prediction, error, difference: " << i + ita << ' ' << prediction << ' ' << err << ' ' << difference << std::endl;
    } // for
    ita += sub;

  } while(err > param::error_tolerance && ita < param::max_iterations);
#endif

#if 0 // Test full weighting

  auto & mf = *mh[0].get();
  auto & mc = *mh[1].get();
  double k = 2.0;

  // Set field with a single mode
  execute<task::fouriermodes>(mf, ud(mf), k, k);
  // Set up the solution eigenvectors
  execute<task::fourier_fw>(mc, sd(mc), k, k);
  execute<task::full_weighting>(mf, mc, ud(mf), ud(mc));

  // Check solution
  err = norm::fwl2();
  flog(info) << "FW error " << err << std::endl;

#endif

#if 0 // Test interpolation

  auto & mf = *mh[0].get();
  auto & mc = *mh[1].get();
  double k = 2.0;

  // Set field with a single mode
  execute<task::fouriermodes>(mc, ud(mc), k, k);
  // Set up the solution eigenvectors
  execute<task::fourier_interp>(mf, sd(mf), k, k);
  execute<task::bilinear_interpolation>(mc, mf, ud(mc), ud(mf));

  // Check solution
  err = norm::interpl2();
  flog(info) << "Interpolation error " << err << std::endl;

#endif

#if 0 // Test residual
  auto & m = *mh[0].get();
  double k = 1.0;

  // Set field with a single mode
  execute<task::fouriermodes>(m, ud(m), k, k);

  // Set up the solution eigenvectors
  execute<task::fourier_residual>(m, sd(m), k, k);

  // Calculate residual
  execute<task::residual>(m, ud(m), fd(m), rd(m));

  // Check solution
  err = norm::resl2();
  flog(info) << "Residual error " << err << std::endl;
#endif

#if 0 // Test FMG

  auto & m = *mh[0].get();

  // Set the eggcarton problem
  execute<task::eggcarton>(m, ud(m), fd(m), sd(m), Aud(m));

  // Solve
  fmg(param::fine_level);

  // Check difference with solution
  err = norm::errl2();
  flog(info) << "Error " << err << std::endl;

#endif

} // analyze
