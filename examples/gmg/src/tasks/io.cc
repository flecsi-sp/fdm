#include "io.hh"

void
gmg::task::print(mesh::accessor<ro> m, field<double>::accessor<ro, ro> fa) {
  auto f = m.mdcolex<mesh::vertices>(fa);

  std::stringstream ss;
  for(auto j : m.vertices<mesh::y_axis, mesh::logical, true>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      // This only looks nice for single-digit values, i.e., < 10
      ss << (f(i, j) < 0 ? "" : " ") << f(i, j) << " ";
    } // for
    ss << std::endl;
  } // for
  flog(info) << ss.str() << std::endl;
} // print
