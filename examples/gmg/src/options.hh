#ifndef GMG_OPTIONS_HH
#define GMG_OPTIONS_HH

#include <flecsi/execution.hh>

namespace gmg::opt {

inline flecsi::program_option<std::string> config("yaml file",
  "The yaml config file.",
  1,
  [](std::string const & v, std::stringstream & ss) {
    return v.find(".yaml") != std::string::npos
             ? true
             : (ss << "file(" << v << ") has invalid suffix") && false;
  });

inline flecsi::program_option<std::size_t> x_levels("x-levels",
  "The x levels of the mesh hierarchy.",
  1,
  [](std::size_t v, std::stringstream & ss) {
    return v > 2 ? true
                 : (ss << "extents must be greater than 4" << v) && false;
  });

inline flecsi::program_option<std::size_t> y_levels("y-levels",
  "The y levels of the mesh hierarchy.",
  1,
  [](std::size_t v, std::stringstream & ss) {
    return v > 2 ? true
                 : (ss << "extents must be greater than 4" << v) && false;
  });

inline flecsi::program_option<int> colors("MUSCL Options",
  "colors,c",
  "Specify the number of colors (default: num processes).",
  {{flecsi::option_default, -1}});

inline flecsi::program_option<std::string> flog_tags("FLOG Options",
  "tags,t",
  "Specify the flog tags to enable.",
  {{flecsi::option_default, "all"}});

inline flecsi::program_option<int> flog_verbose("FLOG Options",
  "verbose,v",
  "Enable verbose output. Passing '-1' will strip any additional"
  " decorations added by flog and will only output the user's message.",
  {{flecsi::option_default, 0}});

} // namespace gmg::opt

#endif // GMG_OPTIONS_HH
