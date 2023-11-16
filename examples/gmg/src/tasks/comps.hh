#ifndef GMG_TASK_COMPS_HH
#define GMG_TASK_COMPS_HH

#include "../types.hh"

namespace gmg::task {

void full_weighting(mesh::accessor<ro> mf,
  mesh::accessor<ro> mc,
  field<double>::accessor<ro, ro> ffa,
  field<double>::accessor<rw, ro> fca);

void bilinear_interpolation(mesh::accessor<ro> mc,
  mesh::accessor<ro> mf,
  field<double>::accessor<ro, ro> fca,
  field<double>::accessor<rw, ro> ffa);

void red(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa);

void black(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa);

void residual(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<ro, ro> fa,
  field<double>::accessor<rw, ro> ra);

} // namespace gmg::task

#endif // GMG_TASK_COMPS_HH
