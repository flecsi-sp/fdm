#include "init.hh"
#include "control.hh"
#include "options.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/init.hh"
#include "tasks/io.hh"

#include <flecsi/execution.hh>
#include <yaml-cpp/yaml.h>

#include <cmath>

using namespace gmg;
using namespace flecsi;

void
action::init(control_policy & cp) {
  YAML::Node config = YAML::LoadFile(opt::config.value());

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();

  /*--------------------------------------------------------------------------*
    Distribution for all grids is based on the fine grid.
   *--------------------------------------------------------------------------*/

  std::size_t x_pow = opt::x_extents.value();
  std::size_t y_pow = opt::y_extents.value();
  std::size_t vertices_x = std::pow(2, x_pow) + 1;
  std::size_t vertices_y = std::pow(2, y_pow) + 1;

  auto parts = mesh::distribute(num_colors, {vertices_x, vertices_y});

  /*--------------------------------------------------------------------------*
    Problem geometry.
   *--------------------------------------------------------------------------*/

  mesh::grect geom;
  geom[0][0] = config["coords"][0][0].as<double>();
  geom[0][1] = config["coords"][1][0].as<double>();
  geom[1][0] = config["coords"][0][1].as<double>();
  geom[1][1] = config["coords"][1][1].as<double>();

  /*--------------------------------------------------------------------------*
    Initialize the mesh hierarchy.
   *--------------------------------------------------------------------------*/

  do {
    mesh::gcoord axis_extents{vertices_x, vertices_y};
    flog(info) << "vertices_x: " << vertices_x << " vertices_y: " << vertices_y
               << std::endl;
    auto & m = *mh.emplace_back(std::make_unique<mesh::slot>());
    m.allocate(mesh::mpi_coloring{parts, axis_extents}, geom);

    // execute<task::enumerate>(m, ud(m));
    execute<task::constant>(m, ud(m), double(x_pow));
    // execute<task::constant>(m, ud(m), 1.0);

    vertices_x = std::pow(2, --x_pow) + 1;
    vertices_y = std::pow(2, --y_pow) + 1;
  } while(x_pow > 2 && y_pow > 2);

  for(auto const & m : mh) {
    execute<task::print>(*m, ud(*m));
  } // for

#if 1
  auto const & m0 = mh[0].get();
  auto const & m1 = mh[1].get();
  // execute<task::full_weighting>(*m0, *m1, ud(*m0), ud(*m1));
  // execute<task::print>(*m1, ud(*m1));
  execute<task::bilinear_interpolation>(*m1, *m0, ud(*m1), ud(*m0));
  execute<task::print>(*m0, ud(*m0));
#endif
} // setup
