#include "init.hh"
#include "control.hh"
#include "options.hh"
#include "state.hh"

#include <flecsi/flog.hh>

using namespace flecsi;

void
gs::action::init_mesh(control_policy & cp) {
  flog(info) << "Initializing " << opt::x_extents.value() << "x"
             << opt::y_extents.value() << " mesh" << std::endl;
  flecsi::flog::flush();

  mesh::gcoord axis_extents{opt::x_extents.value(), opt::y_extents.value()};

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();
  flog(info) << "colors: " << num_colors << std::endl;

  mesh::grect geometry;
  geometry[0][0] = 0.0;
  geometry[0][1] = 1.0;
  geometry[1] = geometry[0];

  cp.m.allocate(mesh::mpi_coloring{num_colors, axis_extents}, geometry);
} // init_mesh
