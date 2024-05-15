/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_STATE_HH
#define GS_STATE_HH

#include "types.hh"

namespace gs {

inline const field<double>::definition<mesh, mesh::vertices> test;

inline const field<double>::definition<mesh, mesh::vertices> ud;
inline const field<double>::definition<mesh, mesh::vertices> fd;
inline const field<double>::definition<mesh, mesh::vertices> sd;
inline const field<double>::definition<mesh, mesh::vertices> Aud;

} // namespace gs

#endif // GS_STATE_HH
