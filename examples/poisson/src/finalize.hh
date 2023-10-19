#ifndef POISSON_FINALIZE_HH
#define POISSON_FINALIZE_HH

#include "control.hh"

namespace poisson {
namespace action {

void finalize(control_policy &);
inline control::action<finalize, cp::finalize> finalize_action;

} // namespace action
} // namespace poisson

#endif // POISSON_FINALIZE_HH
