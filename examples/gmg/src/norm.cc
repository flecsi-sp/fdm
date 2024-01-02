#include "norm.hh"

using namespace flecsi;
using namespace gmg;

double
norm::l2(/* std::size_t level */) {
  auto & mf = *mh[0].get();
  execute<task::discrete_operator>(mf, ud(mf), Aud(mf));
  auto residual =
    reduce<task::diff_sum_square, exec::fold::sum>(mf, fd(mf), Aud(mf));
  return std::sqrt(residual.get());
}

double
norm::max(/* std::size_t level */) {
  auto & mf = *mh[0].get();
  execute<task::discrete_operator>(mf, ud(mf), Aud(mf));
  auto max = reduce<task::diff_max, exec::fold::max>(mf, fd(mf), Aud(mf));
  return max.get();
}
