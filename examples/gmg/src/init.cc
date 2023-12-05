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
    Error control.
   *--------------------------------------------------------------------------*/
  
  opt::error_tolerance = config["error_tolerance"].as<double>();
  opt::max_iterations = config["max_iterations"].as<std::size_t>();

  /*--------------------------------------------------------------------------*
    Initialize the mesh hierarchy.
   *--------------------------------------------------------------------------*/

  int l{0};
  do {
    mesh::gcoord axis_extents{vertices_x, vertices_y};
    flog(info) << "vertices_x: " << vertices_x << " vertices_y: " << vertices_y
               << std::endl;
    auto & m = *mh.emplace_back(std::make_unique<mesh::slot>());
    m.allocate(mesh::mpi_coloring{parts, axis_extents}, geom);

    if(config["problem"].as<std::string>() == "eggcarton") {
      if(l++ == 0) {
        execute<task::eggcarton>(m, ud[0](m), fd(m), sd(m), Aud(m));
        execute<task::io>(m, fd(m), "rhs");
        execute<task::io>(m, sd(m), "actual");
      }
      else {
        execute<task::constant>(m, ud[0](m), 0.0);
        execute<task::constant>(m, fd(m), 0.0);
        execute<task::constant>(m, sd(m), 0.0);
        execute<task::constant>(m, Aud(m), 0.0);
      } // if
    }
    else if(config["problem"].as<std::string>() == "enumerate") {
      execute<task::enumerate>(m, ud[0](m));
      execute<task::enumerate>(m, fd(m));
      execute<task::enumerate>(m, sd(m));
      execute<task::enumerate>(m, Aud(m));
    }
    else if(config["problem"].as<std::string>() == "constant") {
      execute<task::constant>(m, ud[0](m), double(x_pow));
      execute<task::constant>(m, fd(m), x_pow);
      execute<task::constant>(m, sd(m), x_pow);
      execute<task::constant>(m, Aud(m), x_pow);
    } // if

    vertices_x = std::pow(2, --x_pow) + 1;
    vertices_y = std::pow(2, --y_pow) + 1;
  } while(x_pow > 2 && y_pow > 2);
} // setup
