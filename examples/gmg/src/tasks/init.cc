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
gmg::task::bilinear(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> fa,
  double M,
  double N,
  double D) {
  auto f = m.mdcolex<mesh::vertices>(fa);
  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    auto y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      auto x = m.value<mesh::x_axis>(i);
      f(i, j) = M * x + N * y + D;
    } // for
  }; // forall
} // bilinear

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
  const double sq_klpi = pow(M_PI, 2) * (pow(K, 2) + pow(L, 2));

  forall(j, (m.vertices<mesh::y_axis, mesh::logical>()), "init_eggcarton") {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      f(i, j) = sq_klpi * sin(K * M_PI * x) * sin(L * M_PI * y);
      const double solution = sin(K * M_PI * x) * sin(L * M_PI * y);
      s(i, j) = solution;
      Au(i, j) = 0.0;
      u(i, j) = 0.0;
    } // for
  }; // forall
} // eggcarton

void
gmg::task::fouriermodes(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  forall(j, (m.vertices<mesh::y_axis, mesh::logical>()), "init_fouriermodes") {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      u(i, j) = sin(kk * M_PI * x) * sin(ll * M_PI * y);
    } // for
  }; // forall
} // fouriermodes

void
gmg::task::fourier_fw(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  double factor;

  factor = pow(cos(kk * M_PI * m.xdelta() * 0.25), 2) *
           pow(cos(ll * M_PI * m.ydelta() * 0.25), 2);

  forall(j, (m.vertices<mesh::y_axis, mesh::logical>()), "init_fourier_fw") {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      u(i, j) = sin(kk * M_PI * x) * sin(ll * M_PI * y) * factor;
    } // for
  }; // forall
} // fourier_fw

void
gmg::task::fourier_interp(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  double cos2i, sin2i, cos2j, sin2j, kkp, llp;

  kkp = 1 / m.xdelta() - kk;
  llp = 1 / m.ydelta() - ll;

  cos2i = pow(cos(kk * M_PI * m.xdelta() * 0.5), 2);
  sin2i = pow(sin(kk * M_PI * m.xdelta() * 0.5), 2);
  cos2j = pow(cos(ll * M_PI * m.ydelta() * 0.5), 2);
  sin2j = pow(sin(ll * M_PI * m.ydelta() * 0.5), 2);

  forall(
    j, (m.vertices<mesh::y_axis, mesh::logical>()), "init_fourier_interp") {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      u(i, j) = cos2i * sin(kk * M_PI * x) - sin2i * sin(kkp * M_PI * x);
      u(i, j) *= cos2j * sin(ll * M_PI * y) - sin2j * sin(llp * M_PI * y);
    } // for
  }; // forall
} // fourier_interp

void
gmg::task::fourier_residual(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  const auto dx = m.xdelta();
  const auto dy = m.ydelta();
  double factor;

  factor = 4 * (pow(sin(kk * M_PI * dx * 0.5), 2) / pow(dx, 2) +
                 pow(sin(ll * M_PI * dy * 0.5), 2) / pow(dy, 2));

  forall(
    j, (m.vertices<mesh::y_axis, mesh::logical>()), "init_fourier_residual") {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      u(i, j) = -factor * sin(kk * M_PI * x) * sin(ll * M_PI * y);
    } // for
  }; // forall
} // fourier_residual

void
gmg::task::gs_eigenvector(mesh::accessor<ro> m,
  field<double>::accessor<wo, na> ua,
  double kk,
  double ll) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  const auto dx = m.xdelta();
  const auto dy = m.ydelta();
  const auto dx_over_dy = dx / dy;
  const auto dy_over_dx = dy / dx;
  const auto factor = 1.0 / (dx_over_dy + dy_over_dx);
  double cos2i, sin2i, cos2j, sin2j, kkp, llp;
  double cos_f, sin_f;

  kkp = 1 / dx - kk;
  llp = 1 / dy - ll;

  cos2i = dy_over_dx * pow(cos(kk * M_PI * dx * 0.5), 2);
  sin2i = dy_over_dx * pow(sin(kk * M_PI * dx * 0.5), 2);
  cos2j = dx_over_dy * pow(cos(ll * M_PI * dy * 0.5), 2);
  sin2j = dx_over_dy * pow(sin(ll * M_PI * dy * 0.5), 2);
  cos_f = factor * (cos2i + cos2j);
  sin_f = factor * (sin2i + sin2j);

  forall(
    j, (m.vertices<mesh::y_axis, mesh::logical>()), "init_gs_eigenvector") {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      u(i, j) = cos_f * sin(kk * M_PI * x) * sin(ll * M_PI * y) +
                sin_f * sin(kkp * M_PI * x) * sin(llp * M_PI * y);
    } // for
  }; // forall
} // gs_eigenvector

void
gmg::task::poisson_stencil(mesh::accessor<ro> m,
  stencil_field<five_pt>::accessor<wo, na> soa) {
  auto so = m.stencil_op<mesh::vertices, five_pt>(soa);

  const double dx{m.xdelta()};
  const double dy{m.ydelta()};
  const double dx_over_dy{m.xdelta() / m.ydelta()};
  const double dy_over_dx{m.ydelta() / m.xdelta()};
  const double w{1.0 / m.dxdy()};

  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      so(i, j, five_pt::c) = 2.0 * w * (dx_over_dy + dy_over_dx);
      so(i, j, five_pt::w) = w * dy_over_dx;
      so(i, j, five_pt::s) = w * dx_over_dy;
    } // for
  } // for
} // poisson_stencil

void
gmg::task::turner_stencil(mesh::accessor<ro> m,
  field<double>::accessor<ro, na> ud1,
  field<double>::accessor<ro, na> ud2,
  stencil_field<five_pt>::accessor<wo, na> soa,
  double dt) {

  auto u1 = m.mdcolex<mesh::vertices>(ud1);
  auto u2 = m.mdcolex<mesh::vertices>(ud2);
  auto so = m.stencil_op<mesh::vertices, five_pt>(soa);

  const double dx{m.xdelta()};
  const double dy{m.ydelta()};
  const double dx_over_dy{m.xdelta() / m.ydelta()};
  const double dy_over_dx{m.ydelta() / m.xdelta()};
  const double w{dt / m.dxdy()};

  // This is the stencil for the summer school
  // The field "ua" has the coefficients for our variable stencil

  for(auto j : m.vertices<mesh::y_axis, mesh::interior>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::interior>()) {
      so(i, j, five_pt::c) = 1 + w * (dy_over_dx * (u1(i + 1, j) + u1(i, j)) +
                                       dx_over_dy * (u2(i, j + 1) + u2(i, j)));
      so(i, j, five_pt::w) = w * dy_over_dx * u1(i, j);
      so(i, j, five_pt::s) = w * dx_over_dy * u2(i, j);
    } // for
  } // for
} // turner_stencil
