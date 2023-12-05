#ifndef GMG_TYPES_HH
#define GMG_TYPES_HH

#include <fdm/mesh.hh>

namespace gmg {

inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using mesh = fdm::mesh;

template<typename T, flecsi::data::layout L = flecsi::data::layout::dense>
using field = flecsi::field<T, L>;

// Change to std::swap
struct unknowns {
  using type = const field<double>::definition<mesh, mesh::vertices>;

  auto operator[](int i) {
    return ud_[(flip_ + i) % 2];
  }
  auto flip() {
    ++flip_;
  }

  operator type() {
    return operator[](0);
  }

private:
  type ud_[2];
  int flip_{0};
}; // struct unknowns

} // namespace gmg

#endif // GMG_TYPES_HH
