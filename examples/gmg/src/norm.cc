#include "norm.hh"

using namespace flecsi;
using namespace gmg;

double
norm::l2(flecsi::scheduler & sc) {
  auto & mf = *mh[0];
  sc.execute<task::discrete_operator>(exec::on, mf, sod(mf), ud(mf), Aud(mf));
  return std::sqrt(sc.reduce<task::diff_sum_square, exec::fold::sum>(
                       exec::on, mf, fd(mf), Aud(mf))
                     .get());
}

double
norm::fwl2(flecsi::scheduler & sc) {
  auto & mc = *mh[1];
  return std::sqrt(sc.reduce<task::diff_sum_square, exec::fold::sum>(
                       exec::on, mc, ud(mc), sd(mc))
                     .get());
}

double
norm::interpl2(flecsi::scheduler & sc) {
  auto & mf = *mh[0];
  return std::sqrt(sc.reduce<task::diff_sum_square, exec::fold::sum>(
                       exec::on, mf, ud(mf), sd(mf))
                     .get());
}

double
norm::resl2(flecsi::scheduler & sc) {
  auto & mf = *mh[0];
  return std::sqrt(sc.reduce<task::diff_sum_square, exec::fold::sum>(
                       exec::on, mf, rd(mf), sd(mf))
                     .get());
}

double
norm::errl2(flecsi::scheduler & sc) {
  auto & mf = *mh[0];
  sc.execute<task::discrete_operator>(exec::on, mf, sod(mf), sd(mf), Aud(mf));
  return std::sqrt(sc.reduce<task::diff_sum_square, exec::fold::sum>(
                       exec::on, mf, fd(mf), Aud(mf))
                     .get());
}

double
norm::max(flecsi::scheduler & sc) {
  auto & mf = *mh[0];
  execute<task::discrete_operator>(exec::on, mf, sod(mf), ud(mf), Aud(mf));
  return reduce<task::diff_max, exec::fold::max>(exec::on, mf, fd(mf), Aud(mf))
    .get();
}
