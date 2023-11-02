#ifndef GS_CONTROL_HH
#define GS_CONTROL_HH

#include "fdm/mesh.hh"

#include <flecsi/flog.hh>
#include <flecsi/run/control.hh>

namespace poisson {

enum class cp { initialize, solve, analyze, finalize };

inline const char *
operator*(cp control_point) {
  switch(control_point) {
    case cp::initialize:
      return "initialize";
    case cp::solve:
      return "solve";
    case cp::analyze:
      return "analyze";
    case cp::finalize:
      return "finalize";
  }
  flog_fatal("invalid control point");
}

struct control_policy : flecsi::run::control_base {

  using control_points_enum = cp;
  struct node_policy {};

  using control = flecsi::run::control<control_policy>;

  using control_points = list<point<cp::initialize>,
    point<cp::solve>,
    point<cp::analyze>,
    point<cp::finalize>>;

  fdm::mesh::slot m;
}; // struct control_policy

using control = flecsi::run::control<control_policy>;

} // namespace poisson

#endif // GS_CONTROL_HH
