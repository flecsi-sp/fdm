#ifndef POISSON_TYPES_HH
#define POISSON_TYPES_HH

#include "fdm/mesh.hh"

#include <flecsi/data.hh>

namespace poisson {

inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using mesh = fdm::mesh;

template<typename T, flecsi::data::layout L = flecsi::data::layout::dense>
using field = flecsi::field<T, L>;

} // namespace poisson

#endif // POISSON_TYPES_HH
