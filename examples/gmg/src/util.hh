#ifndef GMG_UTIL_HH
#define GMG_UTIL_HH

#include "options.hh"
#include "state.hh"

namespace gmg::util {
  inline auto level(std::size_t index) {
    return param::fine_level - index;
  } // level

  inline auto index(std::size_t level) {
    return param::fine_level - level;
  } // index

} // namespace gmg::util

#endif // GMG_UTIL_HH
