#include "init.hh"
#include "control.hh"
#include "options.hh"
#include "state.hh"
#include "tasks/comps.hh"
#include "tasks/init.hh"
#include "tasks/io.hh"
#include "util.hh"

#include <flecsi/execution.hh>
#include <yaml-cpp/yaml.h>

#include <cmath>
#include <sstream>

using namespace gmg;
using namespace flecsi;

void
action::init(control_policy & cp) {
  YAML::Node config = YAML::LoadFile(opt::config.value());

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();

  /*--------------------------------------------------------------------------*
    Hierarchy setup.
   *--------------------------------------------------------------------------*/

  auto x_levels = opt::x_levels.value();
  auto y_levels = opt::y_levels.value();

  param::fine_level = std::min(x_levels, y_levels);
  param::vcycle_direct = config["vcycle_direct"].as<std::size_t>();

  flog_assert(param::fine_level > param::vcycle_direct,
    "fine level(" << param::fine_level
                  << ") must be greater than direct solver level("
                  << param::vcycle_direct << ")");

  std::stringstream ss;
  ss << "Fine Level: " << param::fine_level << std::endl;
  ss << "Direct Solver Level: " << param::vcycle_direct << std::endl;

  /*--------------------------------------------------------------------------*
    Distribution for all grids is based on the fine grid.
   *--------------------------------------------------------------------------*/

  std::size_t vertices_x = std::pow(2, x_levels) + 1;
  std::size_t vertices_y = std::pow(2, y_levels) + 1;
  ss << "Fine Grid Resolution: (" << vertices_x << "x" << vertices_y << ")"
     << std::endl;

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
    Iteration, error, and output control.
   *--------------------------------------------------------------------------*/

  param::error_tolerance = config["error_tolerance"].as<double>();
  param::max_iterations = config["max_iterations"].as<std::size_t>();
  ss << "Error Tolerance: " << param::error_tolerance << std::endl;
  ss << "Max Iterations: " << param::max_iterations << std::endl;
  ss << "Log Frequency: " << config["log_frequency"].as<std::string>()
     << std::endl;

  flog(info) << ss.str() << std::endl;

  /*--------------------------------------------------------------------------*
    Sovler parameters.
   *--------------------------------------------------------------------------*/

  param::vcycle_pre = config["vcycle_pre"].as<std::size_t>();
  param::vcycle_post = config["vcycle_post"].as<std::size_t>();

  /*--------------------------------------------------------------------------*
    Initialize the mesh hierarchy.
   *--------------------------------------------------------------------------*/

  int index{0};
  do {
    mesh::gcoord axis_extents{vertices_x, vertices_y};
    auto & m = *mh.emplace_back(std::make_unique<mesh::slot>());
    m.allocate(mesh::mpi_coloring{parts, axis_extents}, geom);

    if(config["problem"].as<std::string>() == "eggcarton") {
      if(index == 0) {
        execute<task::eggcarton>(m, ud(m), fd(m), sd(m), Aud(m));
        execute<task::constant>(m, rd(m), 0.0);
        execute<task::constant>(m, ed(m), 0.0);
        execute<task::io>(m, fd(m), "rhs");
        execute<task::io>(m, sd(m), "actual");
      }
      else {
        execute<task::constant>(m, ud(m), 0.0);
        execute<task::constant>(m, fd(m), 0.0);
        execute<task::constant>(m, sd(m), 0.0);
        execute<task::constant>(m, rd(m), 0.0);
        execute<task::constant>(m, ed(m), 0.0);
        execute<task::constant>(m, Aud(m), 0.0);
      } // if
    }
    else if(config["problem"].as<std::string>() == "enumerate") {
      execute<task::enumerate>(m, ud(m));
      execute<task::enumerate>(m, fd(m));
      execute<task::enumerate>(m, sd(m));
      execute<task::enumerate>(m, Aud(m));
    }
    else if(config["problem"].as<std::string>() == "bilinear") {
      execute<task::bilinear>(m, ud(m), 1.0, 1.0, 0.0);
      execute<task::bilinear>(m, fd(m), 1.0, 1.0, 0.0);
      execute<task::bilinear>(m, sd(m), 1.0, 1.0, 0.0);
      execute<task::bilinear>(m, rd(m), 1.0, 1.0, 0.0);
      execute<task::bilinear>(m, ed(m), 1.0, 1.0, 0.0);
      execute<task::bilinear>(m, Aud(m), 1.0, 1.0, 0.0);
    }
    else if(config["problem"].as<std::string>() == "constant") {
      execute<task::constant>(m, ud(m), util::level(index));
      execute<task::constant>(m, fd(m), util::level(index));
      execute<task::constant>(m, sd(m), util::level(index));
      execute<task::constant>(m, Aud(m), util::level(index));
    } // if

    vertices_x = std::pow(2, --x_levels) + 1;
    vertices_y = std::pow(2, --y_levels) + 1;
    flog(warn) << "level(index): " << util::level(index) << "(" << index
               << ") X vertices: " << vertices_x
               << " Y vertices: " << vertices_y << std::endl;
    ++index;
  } while(x_levels >= param::vcycle_direct && y_levels >= param::vcycle_direct);
} // setup
