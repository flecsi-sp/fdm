#include "norm.hh"

using namespace gmg;

double
task::diff_sum_square(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba) {
  auto a = m.mdcolex<mesh::vertices>(aa);
  auto b = m.mdcolex<mesh::vertices>(ba);

  double sum{0};
  double d2 = m.xdelta() * m.ydelta();
  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      sum += pow(a(i, j) - b(i, j), 2);
    } // for
  } // for

  return d2 * sum;
} // diff

double
task::diff_max(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> aa,
  field<double>::accessor<ro, ro> ba) {
  auto a = m.mdcolex<mesh::vertices>(aa);
  auto b = m.mdcolex<mesh::vertices>(ba);

  double max{0};
  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      max = std::max(max, std::abs(a(i, j) - b(i, j)));
    } // for
  } // for

  return max;
} // diff

void
task::discrete_operator(mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<ro, na> soa,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<wo, ro> Aua) {
  auto so = m.stencil_op<mesh::vertices, five_pt>(soa);
  auto u = m.mdcolex<mesh::vertices>(ua);
  auto Au = m.mdcolex<mesh::vertices>(Aua);

  forall(j, m.vertices<mesh::y_axis>(), "discrete_operator") {
    for(auto i : m.vertices<mesh::x_axis>()) {
      Au(i, j) = so(i, j, five_pt::c) * u(i, j) -
                 so(i, j, five_pt::w) * u(i - 1, j) -
                 so(i + 1, j, five_pt::w) * u(i + 1, j) -
                 so(i, j, five_pt::s) * u(i, j - 1) -
                 so(i, j + 1, five_pt::s) * u(i, j + 1);
    } // for
  }; // forall
} // discrete_operator
