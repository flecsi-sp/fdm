/*----------------------------------------------------------------------------*
  Copyright (C) 2023, Triad National Security, LLC
  All rights reserved.
 *----------------------------------------------------------------------------*/
#ifndef GS_OPTIONS_HH
#define GS_OPTIONS_HH

#include "flecsi/execution.hh"

#include <limits>

namespace gs::opt {

inline flecsi::program_option<std::size_t>
  x_extents("x-extents", "The x extents of the mesh.", 1);
inline flecsi::program_option<std::size_t>
  y_extents("y-extents", "The y extents of the mesh.", 1);

inline flecsi::program_option<int> colors("Gauss-Seidel Options",
  "colors,c",
  "Specify the number of colors (default: num processes).",
  {{flecsi::option_default, -1}});

inline flecsi::program_option<std::size_t> max_iterations("Poisson Options",
  "max_iterations,m",
  "The maximum number of solver iterations.",
  {{flecsi::option_default, std::numeric_limits<std::size_t>::max()}});

inline flecsi::program_option<double> error_tol("Gauss-Seidel Options",
  "tolerance,t",
  "The maximum absolute residual tolerance.",
  {{flecsi::option_default, 1.0e-4}});

inline flecsi::program_option<std::string> flog_tags("FLOG Options",
  "tags,t",
  "Specify the flog tags to enable.",
  {{flecsi::option_default, "all"}});

inline flecsi::program_option<int> flog_verbose("FLOG Options",
  "verbose,v",
  "Enable verbose output. Passing '-1' will strip any additional"
  " decorations added by flog and will only output the user's message.",
  {{flecsi::option_default, 0}});

} // namespace gs::opt

#endif // GS_OPTIONS_HH
