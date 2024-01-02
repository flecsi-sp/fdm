#include "io.hh"

using namespace flecsi;

void
gmg::task::io(mesh::accessor<ro> m,
  field<double>::accessor<ro, ro> fa,
  std::string filebase) {
  auto f = m.mdcolex<mesh::vertices>(fa);

  std::stringstream ss;
  ss << filebase;
  if(processes() == 1) {
    ss << ".dat";
  }
  else {
    ss << "-" << process() << ".dat";
  } // if

  std::ofstream solution(ss.str(), std::ofstream::out);

  for(auto j : m.vertices<mesh::y_axis, mesh::logical>()) {
    const double y = m.value<mesh::y_axis>(j);
    for(auto i : m.vertices<mesh::x_axis, mesh::logical>()) {
      const double x = m.value<mesh::x_axis>(i);
      solution << x << " " << y << " " << f(i, j) << std::endl;
    } // for
  } // for
} // io

void
gmg::task::print(mesh::accessor<ro> m, field<double>::accessor<ro, ro> fa) {
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
