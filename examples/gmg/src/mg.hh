/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GMG_MG_HH
#define GMG_MG_HH

#include <cstddef>

namespace gmg {

void vcycle(std::size_t level);
void fmg(std::size_t level);

} // namespace gmg

#endif // GMG_MG_HH
