#ifndef GS_INIT_HH
#define GS_INIT_HH

#include "control.hh"

namespace poisson {
namespace action {

void init_mesh(control_policy &);
inline control::action<init_mesh, cp::initialize> init_mesh_action;

} // namespace action
} // namespace poisson

#endif // GS_INIT_HH
