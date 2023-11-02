#ifndef GS_OPTIONS_HH
#define GS_OPTIONS_HH

#include "flecsi/execution.hh"

#include <limits>

namespace poisson::opt {

inline flecsi::program_option<std::size_t>
  x_extents("x-extents", "The x extents of the mesh.", 1);
inline flecsi::program_option<std::size_t>
  y_extents("y-extents", "The y extents of the mesh.", 1);

inline flecsi::program_option<int> colors("Poisson Options",
  "colors,c",
  "Specify the number of colors (default: num processes).",
  {{flecsi::option_default, -1}});

inline flecsi::program_option<std::size_t> max_iterations("Poisson Options",
  "max_iterations,m",
  "The maximum number of solver iterations.",
  {{flecsi::option_default, std::numeric_limits<std::size_t>::max()}});
inline flecsi::program_option<double> error_tol("Poisson Options",
  "tolerance,t",
  "The maximum absolute residual tolerance.",
  {{flecsi::option_default, 1.0e-4}});

} // namespace poisson::opt

#endif // GS_OPTIONS_HH
