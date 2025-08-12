#include "problem.hh"
#include "annotation.hh"
#include "state.hh"
#include "tasks/init.hh"
#include "tasks/io.hh"

#include <flecsi/execution.hh>

using namespace flecsi;

void
gs::action::problem(control_policy & cp) {
  auto & sc = cp.scheduler();
  util::annotation::rguard<problem_region> guard;
  sc.execute<task::eggcarton>(
    exec::on, *cp.m, ud(*cp.m), fd(*cp.m), sd(*cp.m), Aud(*cp.m));
  sc.execute<task::io>(exec::on, *cp.m, ud(*cp.m), "init");
  sc.execute<task::io>(exec::on, *cp.m, sd(*cp.m), "actual");

  // This can be used for debugging
#if 0
  sc.execute<task::redblack>(exec::on, *cp.m, test(*cp.m));
  sc.execute<task::print>(exec::on, *cp.m, test(*cp.m));
#endif

  flog::flush();
} // problem
