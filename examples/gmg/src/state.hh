/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_STATE_HH
#define GMG_STATE_HH

#include "types.hh"

namespace gmg {

inline std::vector<std::unique_ptr<mesh::slot>> mh;

inline dual_field<double> ud;
inline const field<double>::definition<mesh, mesh::vertices> fd;
inline const field<double>::definition<mesh, mesh::vertices> sd;
inline const field<double>::definition<mesh, mesh::vertices> rd;
inline const field<double>::definition<mesh, mesh::vertices> ed;
inline const field<double>::definition<mesh, mesh::vertices> Aud;

namespace param {
inline double error_tolerance;
inline std::size_t max_iterations;
inline std::size_t fine_level;
inline std::size_t mg_direct;
inline std::size_t mg_pre;
inline std::size_t mg_post;
inline std::size_t mg_cycles;
} // namespace param

} // namespace gmg

#endif // GMG_STATE_HH
