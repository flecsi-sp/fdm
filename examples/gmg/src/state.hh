#ifndef GMG_STATE_HH
#define GMG_STATE_HH

#include "types.hh"

namespace gmg {

inline std::vector<std::unique_ptr<mesh::slot>> mh;

inline unknowns ud;
inline const field<double>::definition<mesh, mesh::vertices> fd;
inline const field<double>::definition<mesh, mesh::vertices> sd;
inline const field<double>::definition<mesh, mesh::vertices> rd;
inline const field<double>::definition<mesh, mesh::vertices> ed;
inline const field<double>::definition<mesh, mesh::vertices> Aud;

namespace opt {
inline double error_tolerance{1.0e-10};
inline std::size_t max_iterations{10000};
} // namespace opt

} // namespace gmg

#endif // GMG_STATE_HH
