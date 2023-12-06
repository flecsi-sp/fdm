#ifndef GMG_OPTIONS_HH
#define GMG_OPTIONS_HH

#include <flecsi/execution.hh>

namespace gmg::opt {

inline flecsi::program_option<std::string> config("yaml file",
  "The yaml config file.",
  1,
  [](flecsi::any const & v, std::stringstream & ss) {
    const std::string value = flecsi::option_value<std::string>(v);
    return value.find(".yaml") != std::string::npos
             ? true
             : (ss << "file(" << value << ") has invalid suffix") && false;
  });

inline flecsi::program_option<std::size_t> x_levels("x-levels",
  "The x levels of the mesh hierarchy.",
  1,
  [](flecsi::any const & v, std::stringstream & ss) {
    const std::size_t value = flecsi::option_value<std::size_t>(v);
    return value > 2
             ? true
             : (ss << "extents must be greater than 4" << value) && false;
  });
inline flecsi::program_option<std::size_t> y_levels("y-levels",
  "The y levels of the mesh hierarchy.",
  1,
  [](flecsi::any const & v, std::stringstream & ss) {
    const std::size_t value = flecsi::option_value<std::size_t>(v);
    return value > 2
             ? true
             : (ss << "extents must be greater than 4" << value) && false;
  });

inline flecsi::program_option<int> colors("MUSCL Options",
  "colors,c",
  "Specify the number of colors (default: num processes).",
  {{flecsi::option_default, -1}});
} // namespace gmg::opt

#endif // GMG_OPTIONS_HH
