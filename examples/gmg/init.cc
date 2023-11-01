#include "init.hh"
#include "options.hh"
#include "state.hh"

using namespace gmg;
using namespace flecsi;

void
action::init(control_policy & cp) {
  mh.reserve(10);

  flog(info) << "X extents: " << opt::x_extents.value()
             << " Y extents: " << opt::y_extents.value() << std::endl;

  auto x_coarse =
    std::size_t(opt::x_extents.value() / 2);
  flog(info) << "X coarse: " << x_coarse << std::endl;

  const auto num_colors =
    opt::colors.value() == -1 ? flecsi::processes() : opt::colors.value();

  auto parts = mesh::distribute(
    num_colors, {opt::x_extents.value(), opt::y_extents.value()});

  flog(info) << "parts: " << flog::container{parts} << std::endl;

  mesh::gcoord axis_extents{opt::x_extents.value(), opt::y_extents.value()};
  mesh::grect geom;
  geom[0][0] = 0.0;
  geom[0][1] = 1.0;
  geom[1] = geom[0];

  mh.push_back(std::make_unique<mesh::slot>());
  mh[0]->allocate(mesh::mpi_coloring{parts, axis_extents}, geom);

#if 0
  std::size_t i{0};
  do {
  } while(1);

  for(auto & m: mh) {
  } // for
#endif
} // setup
