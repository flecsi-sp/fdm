#include "comps.hh"

using namespace gmg;

void
task::full_weighting(mesh::accessor<ro> mf,
  mesh::accessor<ro> mc,
  field<double>::accessor<ro, ro> ffa,
  field<double>::accessor<rw, ro> fca) {
  auto ff = mf.mdcolex<mesh::vertices>(ffa);
  auto fc = mc.mdcolex<mesh::vertices>(fca);

  for(auto j : mc.vertices<mesh::y_axis, mesh::interior>()) {
    auto fj = 2 * j - 1;
    for(auto i : mc.vertices<mesh::x_axis, mesh::interior>()) {
      auto fi = 2 * i - 1;
      fc(i, j) = 0.0625 * (ff(fi-1, fj-1) + ff(fi-1, fj+1) + ff(fi+1, fj-1) + ff(fi+1, fj+1) + 2.0*(ff(fi, fj-1) + ff(fi, fj+1) + ff(fi-1, fj) + ff(fi+1, fj)) + 4.0*ff(fi,fj));
    } // for
  } // for
} // full_weighting

void
bilinear_interpolation(mesh::accessor<ro> mc,
  mesh::accessor<ro> mf,
  field<double>::accessor<ro, ro> fca,
  field<double>::accessor<rw, ro> ffa) {} // bilinear_interpolation
