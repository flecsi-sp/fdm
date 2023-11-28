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
      fc(i, j) = 0.0625 * (ff(fi - 1, fj - 1) + ff(fi - 1, fj + 1) +
                            ff(fi + 1, fj - 1) + ff(fi + 1, fj + 1) +
                            2.0 * (ff(fi, fj - 1) + ff(fi, fj + 1) +
                                    ff(fi - 1, fj) + ff(fi + 1, fj)) +
                            4.0 * ff(fi, fj));
    } // for
  } // for
} // full_weighting

void
task::bilinear_interpolation(mesh::accessor<ro> mc,
  mesh::accessor<ro> mf,
  field<double>::accessor<ro, ro> fca,
  field<double>::accessor<rw, ro> ffa) {
  auto fc = mc.mdcolex<mesh::vertices>(fca);
  auto ff = mf.mdcolex<mesh::vertices>(ffa);

  for(auto j : mf.vertices<mesh::y_axis, mesh::interior>()) {
    auto cj = (j - 1) / 2 + 1;
    for(auto i : mf.vertices<mesh::x_axis, mesh::interior>()) {
      auto ci = (i - 1) / 2 + 1;

      if(i % 2 && j % 2) {
        ff(i, j) = 0.25 * (fc(ci, cj) + fc(ci - 1, cj) + fc(ci, cj - 1) +
                            fc(ci - 1, cj - 1));
      }
      else if(i % 2) {
        ff(i, j) = 0.5 * (fc(ci, cj) + fc(ci - 1, cj));
      }
      else if(j % 2) {
        ff(i, j) = 0.5 * (fc(ci, cj) + fc(ci, cj - 1));
      }
      else {
        ff(i, j) = fc(ci, cj);
      } // if
    } // for
  } // for
} // bilinear_interpolation

void
task::damped_jacobi(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua_new,
  field<double>::accessor<rw, ro> ua_old,
  field<double>::accessor<ro, ro> fa,
  double omega) {
  auto u_new = m.mdcolex<mesh::vertices>(ua_new);
  auto u_old = m.mdcolex<mesh::vertices>(ua_old);
  auto f = m.mdcolex<mesh::vertices>(fa);
  const auto dxdy = m.dxdy();
  const auto dx_over_dy = m.xdelta() / m.ydelta();
  const auto dy_over_dx = m.ydelta() / m.xdelta();
  const auto factor = 1.0 / (2 * (dx_over_dy + dy_over_dx));

  forall(j, m.vertices<mesh::y_axis>(), "jacobi") {
    for(auto i : m.vertices<mesh::x_axis>()) {
      const double z =
        factor *
        (dxdy * f(i, j) + dy_over_dx * (u_old(i + 1, j) + u_old(i - 1, j)) +
          dx_over_dy * (u_old(i, j + 1) + u_old(i, j - 1)));
      u_new(i, j) = u_old(i, j) + omega * (z - u_old(i, j));
    } // for
  }; // forall
} // damped_jacobi

void
task::red(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  auto f = m.mdcolex<mesh::vertices>(fa);
  const auto dxdy = m.dxdy();
  const auto dx_over_dy = m.xdelta() / m.ydelta();
  const auto dy_over_dx = m.ydelta() / m.xdelta();
  const auto factor = 1.0 / (2 * (dx_over_dy + dy_over_dx));

  forall(j, m.vertices<mesh::y_axis>(), "red") {
    for(auto i : m.red<mesh::x_axis>(j)) {
      u(i, j) =
        factor * (dxdy * f(i, j) + dy_over_dx * (u(i + 1, j) + u(i - 1, j)) +
                   dx_over_dy * (u(i, j + 1) + u(i, j - 1)));
    } // for
  }; // forall
} // red

void
task::black(mesh::accessor<ro> m,
  field<double>::accessor<rw, ro> ua,
  field<double>::accessor<ro, ro> fa) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  auto f = m.mdcolex<mesh::vertices>(fa);
  const auto dxdy = m.dxdy();
  const auto dx_over_dy = m.xdelta() / m.ydelta();
  const auto dy_over_dx = m.ydelta() / m.xdelta();
  const auto factor = 1.0 / (2 * (dx_over_dy + dy_over_dx));

  forall(j, m.vertices<mesh::y_axis>(), "black") {
    for(auto i : m.black<mesh::x_axis>(j)) {
      u(i, j) =
        factor * (dxdy * f(i, j) + dy_over_dx * (u(i + 1, j) + u(i - 1, j)) +
                   dx_over_dy * (u(i, j + 1) + u(i, j - 1)));
    } // for
  }; // forall
} // black

void
task::residual(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  field<double>::accessor<ro, ro> fa,
  field<double>::accessor<rw, ro> ra) {
  auto u = m.mdcolex<mesh::vertices>(ua);
  auto f = m.mdcolex<mesh::vertices>(fa);
  auto r = m.mdcolex<mesh::vertices>(ra);

  const double w = 1.0 / m.dxdy();
  const auto dx_over_dy = m.xdelta() / m.ydelta();
  const auto dy_over_dx = m.ydelta() / m.xdelta();

  for(auto j : m.vertices<mesh::y_axis>()) {
    for(auto i : m.vertices<mesh::x_axis>()) {
      r(i, j) = f(i, j) - w * (2.0 * (dx_over_dy + dy_over_dx) * u(i, j) -
                                dy_over_dx * (u(i + 1, j) + u(i - 1, j)) -
                                dx_over_dy * (u(i, j - 1) + u(i, j + 1)));
    } // for
  } // for
} // residual
