#ifndef GMG_TASK_ANALYZE_HH
#define GMG_TASK_ANALYZE_HH

#include "../types.hh"

namespace gmg::task {

void product_by_eigenvalue_jb(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> sa,
  double omega,
  double k,
  double l);

void product_by_eigenvalue_gs(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> sa,
  double k,
  double l);

} // namespace gmg::task

#endif // GMG_TASK_ANALYZE_HH
