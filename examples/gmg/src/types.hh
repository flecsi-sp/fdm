#ifndef GMG_TYPES_HH
#define GMG_TYPES_HH

#include <fdm/mesh.hh>

namespace gmg {

inline constexpr flecsi::partition_privilege_t na = flecsi::na, ro = flecsi::ro,
                                               wo = flecsi::wo, rw = flecsi::rw;

using mesh = fdm::mesh;

template<typename T, flecsi::data::layout L = flecsi::data::layout::dense>
using field = flecsi::field<T, L>;

template<typename T>
struct dual_field {
  using type = const field<T>::template definition<mesh, mesh::vertices>;

  auto flip() {
    ++flip_;
  } // flip

  /*!
    Return the requested field reference.

    @tparam S The topology slot type.
    @param  s The topology slot instance.
    @param  i The index of the field value to return.
   */
  template<typename S>
  auto operator()(S & s, int i = 0) const {
    return fd_[(flip_ + i) % 2](s);
  } // operator

private:
  type fd_[2];
  int flip_{0};
};

} // namespace gmg

#endif // GMG_TYPES_HH
