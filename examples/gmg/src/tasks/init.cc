#include "init.hh"

void
gmg::task::constant(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double value) {
  auto f = m.mdcolex<mesh::vertices>(fa);
  forall(j, (m.vertices<mesh::y_axis, mesh::logical>()), "init_constant") {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      f(i, j) = value;
    } // for
  }; // forall
} // constant

void
gmg::task::enumerate(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa) {
  auto f = m.mdcolex<mesh::vertices>(fa);
  for(auto j: m.vertices<mesh::y_axis, mesh::logical>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      f(i, j) = m.size<mesh::x_axis, mesh::logical>()*std::size_t(j) + std::size_t(i);
    } // for
  }; // forall
} // constant
