#include "norm.hh"

using namespace flecsi;
using namespace gmg;

double
norm::l2(/* std::size_t level */) {
  auto & mf = mh[0];
  execute<task::discrete_operator>(mf, sod(mf), ud(mf), Aud(mf));
  return std::sqrt(reduce<task::diff_sum_square, exec::fold::sum>(mf, fd(mf), Aud(mf)).get());
}

double
norm::fwl2(/* std::size_t level */) {
  auto & mc = mh[1];
  return std::sqrt(reduce<task::diff_sum_square, exec::fold::sum>(mc, ud(mc), sd(mc)).get());
}

double
norm::interpl2(/* std::size_t level */) {
  auto & mf = mh[0];
  return std::sqrt(reduce<task::diff_sum_square, exec::fold::sum>(mf, ud(mf), sd(mf)).get());
}

double
norm::resl2(/* std::size_t level */) {
  auto & mf = mh[0];
  return std::sqrt(reduce<task::diff_sum_square, exec::fold::sum>(mf, rd(mf), sd(mf)).get());
}

double
norm::errl2(/* std::size_t level */) {
  auto & mf = mh[0];
  execute<task::discrete_operator>(mf, sod(mf), sd(mf), Aud(mf));
  return std::sqrt(reduce<task::diff_sum_square, exec::fold::sum>(mf, fd(mf), Aud(mf)).get());
}

double
norm::max(/* std::size_t level */) {
  auto & mf = mh[0];
  execute<task::discrete_operator>(mf, sod(mf), ud(mf), Aud(mf));
  return reduce<task::diff_max, exec::fold::max>(mf, fd(mf), Aud(mf)).get();
}
