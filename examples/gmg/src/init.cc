#include "control.hh"
#include "init.hh"
#include "options.hh"
#include "state.hh"

#include <flecsi/flog.hh>
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
    mh.push_back(std::make_unique<mesh::slot>());
    mesh::gcoord axis_extents{vertices_x, vertices_y};
    flog(info) << "vertices_x: " << vertices_x << " vertices_y: " << vertices_y
               << std::endl;
    mh[mh.size() - 1]->allocate(mesh::mpi_coloring{parts, axis_extents}, geom);

    vertices_x = std::pow(2, --x_pow) + 1;
    vertices_y = std::pow(2, --y_pow) + 1;
  } while(x_pow > 4 && y_pow > 4);
} // setup
