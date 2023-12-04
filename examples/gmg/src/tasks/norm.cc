#include "norm.hh"

using namespace gmg;

double task::diff(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba) {
  auto a = m.mdcolex<mesh::vertices>(aa);
  auto b = m.mdcolex<mesh::vertices>(ba);

  double sum{0};
  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      sum += pow(a(i,j) - b(i,j), 2);
    } // for
  } // for

  return sum;
} // diff

void task::discrete_operator(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<rw, ro> Aua) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  auto Au = m.mdcolex<mesh::vertices>(Aua);

  const double w = 1.0 / m.dxdy();
  const auto dx_over_dy = m.xdelta() / m.ydelta();
  const auto dy_over_dx = m.ydelta() / m.xdelta();

  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      Au(i,j) = w * (2.0 * (dx_over_dy + dy_over_dx) * u(i,j) - dy_over_dx * (u(i + 1,j) + u(i - 1,j)) - dx_over_dy * (u(i, j+1) + u(i, j-1)));
    } // for
  } // for
} // discrete_operator
