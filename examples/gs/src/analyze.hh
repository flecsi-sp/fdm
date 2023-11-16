#ifndef GS_ANALYZE_HH
#define GS_ANALYZE_HH

#include "control.hh"

namespace gs::action {

void analyze(control_policy &);
inline control::action<analyze, cp::analyze> analyze_action;

} // namespace gs::action

#endif
