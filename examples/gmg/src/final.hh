#ifndef GMG_FINAL_HH
#define GMG_FINAL_HH

#include "control.hh"

namespace gmg::action {

void finalize(control_policy &);
inline control::action<finalize, cp::finalize> finalize_action;

} // namespace gmg::action

#endif // GMG_FINAL_HH
