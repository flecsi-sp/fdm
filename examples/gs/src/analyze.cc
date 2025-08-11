#include "analyze.hh"
#include "annotation.hh"
#include "state.hh"
#include "tasks/norm.hh"

#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

#include <cmath>

using namespace flecsi;

void
gs::action::analyze(control_policy & cp) {
  auto & sc = cp.scheduler();
  util::annotation::rguard<analyze_region> guard;

  future<double> sum = sc.reduce<task::diff, exec::fold::sum>(
    exec::on, *cp.m, ud(*cp.m), sd(*cp.m));
  auto scaled = sc.execute<task::scale>(exec::on, *cp.m, sum);

  // Display L2 error
  sc.execute<task::display_l2>(exec::on, *cp.m, scaled);
} // analyze
