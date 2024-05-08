#include "analyze.hh"

using namespace gmg;

void
task::product_by_eigenvalue_jb(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> sa,
  double omega,
  double kk,
  double ll) {
  auto s = m.mdcolex<mesh::vertices>(sa);
  const auto dx_over_dy = m.xdelta() / m.ydelta();
  const auto dy_over_dx = m.ydelta() / m.xdelta();
  const auto factor = 1.0 / (dx_over_dy + dy_over_dx);
  double lambda;

  lambda = dy_over_dx * pow(sin(kk * M_PI * m.xdelta() * 0.5), 2) +
           dx_over_dy * pow(sin(ll * M_PI * m.ydelta() * 0.5), 2);
  lambda = 1 - 2 * omega * factor * lambda;

  forall(j, (m.vertices<mesh::y_axis>()), "product_by_eigenvalue_jb") {
    for(auto i : m.vertices<mesh::x_axis>()) {
      s(i, j) *= lambda;
    } // for
  }; // forall
} // product_by_eigenvalue_jb

void
task::product_by_eigenvalue_gs(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> sa,
  double kk,
  double ll) {
  auto s = m.mdcolex<mesh::vertices>(sa);
  const auto dx = m.xdelta();
  const auto dy = m.ydelta();
  const auto dx_over_dy = dx / dy;
  const auto dy_over_dx = dy / dx;
  const auto factor = 1.0 / (dx_over_dy + dy_over_dx);
  const auto cos_k = cos(kk * M_PI * dx);
  const auto cos_l = cos(ll * M_PI * dy);
  const auto lambda =
    pow(factor * (dy_over_dx * cos_k + dx_over_dy * cos_l), 2);

  forall(j, (m.vertices<mesh::y_axis>()), "product_by_eigenvalue_gs") {
    for(auto i : m.vertices<mesh::x_axis>()) {
      s(i, j) *= lambda;
    } // for
  }; // forall
} // product_by_eigenvalue_gs
