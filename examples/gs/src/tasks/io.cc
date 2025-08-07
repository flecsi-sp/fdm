#include "io.hh"

#include <fstream>
#include <sstream>

using namespace flecsi;

void
gs::task::io(exec::cpu s,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> ua,
  std::string filebase) noexcept {
  auto u = m.mdspan<mesh::vertices>(ua);

  std::stringstream ss;
  ss << filebase;
  if(s.launch().size == 1) {
    ss << ".dat";
  }
  else {
    ss << "-" << s.launch().index << ".dat";
  } // if

  std::ofstream solution(ss.str(), std::ios::out);

  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      solution << x << " " << y << " " << u[j][i] << std::endl;
    } // for
  } // for
} // io

void
gs::task::print(exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> fa) noexcept {
  auto f = m.mdspan<mesh::vertices>(fa);

  std::stringstream ss;
  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      // This only looks nice for single-digit values, i.e., < 10
      ss << (f[j][i] < 0 ? "" : " ") << f[j][i] << " ";
    } // for
    ss << std::endl;
  } // for
  flog(info) << ss.str() << std::endl;
} // print
