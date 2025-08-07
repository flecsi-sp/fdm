#include "io.hh"

#include <fstream>
#include <sstream>

using namespace flecsi;

void
gmg::task::io(exec::cpu s,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> fa,
  std::string filebase) noexcept {
  auto f = m.mdcolex<mesh::vertices>(fa);

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
      solution << x << " " << y << " " << f(i, j) << std::endl;
    } // for
  } // for
} // io

void
gmg::task::print(exec::cpu,
  mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> fa) noexcept {
  auto f = m.mdcolex<mesh::vertices>(fa);

  std::stringstream ss;
  for(auto j : m.vertices<mesh::y_axis, mesh::logical, true>()) {
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      // This only looks nice for single-digit values, i.e., < 10
      ss << (f(i, j) < 0 ? "" : " ") << std::setw(5) << f(i, j) << " ";
    } // for
    ss << std::endl;
  } // for
  flog(info) << ss.str() << std::endl;
} // print
