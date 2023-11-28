#ifndef GMG_STATE_HH
#define GMG_STATE_HH

#include "types.hh"

namespace gmg {

inline std::vector<std::unique_ptr<mesh::slot>> mh;

template<bool F = false>
struct unknowns {
  auto get() {
    int a = flip_ % 2;
    int b = (flip_ + 1) % 2;
    return std::make_pair(ud_[a], ud_[b]);
  }
  auto flip() {
    ++flip_;
  }

private:
  const field<double>::definition<mesh, mesh::vertices> ud_[2];
  int flip_{0};
};

inline unknowns ud;
// inline const field<double>::definition<mesh, mesh::vertices> ud;
inline const field<double>::definition<mesh, mesh::vertices> fd;
inline const field<double>::definition<mesh, mesh::vertices> sd;
inline const field<double>::definition<mesh, mesh::vertices> rd;

} // namespace gmg

#endif // GMG_STATE_HH
