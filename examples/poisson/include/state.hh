#ifndef POISSON_STATE_HH
#define POISSON_STATE_HH

#include "types.hh"

namespace poisson {

inline const field<double>::definition<mesh, mesh::vertices> test;

inline const field<double>::definition<mesh, mesh::vertices> ud;
inline const field<double>::definition<mesh, mesh::vertices> fd;
inline const field<double>::definition<mesh, mesh::vertices> sd;
inline const field<double>::definition<mesh, mesh::vertices> Aud;

} // namespace poisson

#endif // POISSON_STATE_HH
