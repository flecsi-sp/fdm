#include "init.hh"

void
gmg::task::enumerate(mesh::accessor<ro> m, field<double>::accessor<wo, na> fa) {
  auto f = m.mdcolex<mesh::vertices>(fa);
  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      f(i, j) =
        m.size<mesh::x_axis, mesh::logical>() * std::size_t(j) + std::size_t(i);
    } // for
  }; // forall
} // constant

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

constexpr double PI = 3.14159;
constexpr double K = 12.0;
constexpr double L = 2.0;

void
gmg::task::eggcarton(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  field<double>::accessor<wo, na> fa,
  field<double>::accessor<wo, na> sa,
  field<double>::accessor<wo, na> Aua) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  auto f = m.mdcolex<mesh::vertices>(fa);
  auto s = m.mdcolex<mesh::vertices>(sa);
  auto Au = m.mdcolex<mesh::vertices>(Aua);
  const double sq_klpi = pow(PI, 2) * (pow(K, 2) + pow(L, 2));

  flog(info) << "dxdy: " << m.dxdy() << std::endl;

  forall(j, (m.vertices<mesh::y_axis, mesh::logical>()), "init_eggcarton") {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      f(i, j) = sq_klpi * sin(K * PI * x) * sin(L * PI * y);
      const double solution = sin(K * PI * x) * sin(L * PI * y);
      s(i, j) = solution;
      Au(i, j) = 0.0;
      u(i, j) = 0.0;
    } // for
  }; // forall
} // eggcarton
