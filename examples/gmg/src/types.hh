#ifndef GMG_TYPES_HH
#define GMG_TYPES_HH

#include <fdm/mesh.hh>

namespace gmg {

inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using mesh = fdm::mesh;

template<typename T, flecsi::data::layout L = flecsi::data::layout::dense>
using field = flecsi::field<T, L>;

} // namespace gmg

#endif // GMG_TYPES_HH
