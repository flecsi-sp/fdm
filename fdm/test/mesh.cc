#include "../mesh.hh"

#include <flecsi/util/unit.hh>

using namespace flecsi;
using namespace fdm;

const flecsi::field<double>::definition<mesh, mesh::vertices> rho;

constexpr std::size_t vertices_x{10}, vertices_y{10};
constexpr mesh::grect geom{{{0.0, 1.0}, {0.0, 1.0}}};

void
init_mesh(mesh::accessor<ro> m, field<double>::accessor<wo, wo> f_a) {
  auto f = m.mdcolex<mesh::vertices>(f_a);

  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      f(i, j) = 1.0;
    } // for
  } // for

  // This overwrites part of the previous loops but I'm lazy...
  for(auto j : m.vertices<mesh::y_axis, mesh::interior>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::interior>()) {
      f(i, j) = 2.0;
    } // for
  } // for
} // init_mesh

int
verify_mesh(mesh::accessor<ro> m, field<double>::accessor<wo, wo> f_a) {
  UNIT("TASK") {
    auto f = m.mdcolex<mesh::vertices>(f_a);
    double y_low = geom[mesh::y_axis][mesh::low];
    double y_delta =
      (geom[mesh::y_axis][mesh::high] - geom[mesh::y_axis][mesh::low]) /
      (vertices_y - 1);
    double x_low = geom[mesh::x_axis][mesh::low];
    double x_delta =
      (geom[mesh::x_axis][mesh::high] - geom[mesh::x_axis][mesh::low]) /
      (vertices_x - 1);
    double y_value = y_low;

    EXPECT_LT(std::abs(m.ydelta() - y_delta), 10e-16);
    EXPECT_LT(std::abs(m.xdelta() - x_delta), 10e-16);
    for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
      EXPECT_LT(std::abs(m.value<mesh::y_axis>(j) - y_value), 10e-16);
      y_value += y_delta;
      double x_value = x_low;

      EXPECT_TRUE(
        (j == 0 && m.is_boundary<mesh::y_axis, mesh::low>(j)) ||
        (j == vertices_y - 1 && m.is_boundary<mesh::y_axis, mesh::high>(j)) ||
        (!m.is_boundary<mesh::y_axis, mesh::low>(j) &&
          !m.is_boundary<mesh::y_axis, mesh::high>(j)));

      for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
        EXPECT_LT(std::abs(m.value<mesh::x_axis>(i) - x_value), 10e-16);
        x_value += x_delta;

        EXPECT_TRUE(
          (i == 0 && m.is_boundary<mesh::x_axis, mesh::low>(i)) ||
          (i == vertices_y - 1 && m.is_boundary<mesh::x_axis, mesh::high>(i)) ||
          (!m.is_boundary<mesh::x_axis, mesh::low>(i) &&
            !m.is_boundary<mesh::x_axis, mesh::high>(i)));
      } // for
    } // for
  };
} // verify_mesh

int
fdm_mesh() {
  UNIT("DRIVER") {
    mesh::gcoord axis_extents{vertices_x, vertices_y};
    mesh::slot m;
    auto parts = mesh::distribute(processes(), {vertices_x, vertices_y});
    m.allocate(mesh::mpi_coloring{parts, axis_extents}, geom);

    execute<init_mesh>(m, rho(m));
    EXPECT_EQ(test<verify_mesh>(m, rho(m)), 0);
  }; // UNIT
} // mesh

flecsi::util::unit::driver<fdm_mesh> driver;
