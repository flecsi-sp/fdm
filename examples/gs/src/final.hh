#ifndef GS_FINAL_HH
#define GS_FINAL_HH

#include "control.hh"

namespace poisson {
namespace action {

void finalize(control_policy &);
inline control::action<finalize, cp::finalize> finalize_action;

} // namespace action
} // namespace poisson

#endif // GS_FINAL_HH
