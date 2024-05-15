#ifndef FDM_MESH_HH
#define FDM_MESH_HH

#include <flecsi/data.hh>
#include <flecsi/execution.hh>
#include <flecsi/flog.hh>

#include <ranges>

namespace fdm {

namespace util {

template<class T>
constexpr std::enable_if_t<std::is_unsigned_v<T>, T>
ceil_div(T a, T b) {
  return a / b + !!(a % b); // avoids overflow in (a+(b-1))/b
} // ceil_div

template<class R>
constexpr auto
stride_view(R && r,
  typename std::iterator_traits<decltype(std::begin(
    std::declval<R>()))>::difference_type n,
  decltype(n) o = 0) {
  using I = std::make_unsigned_t<decltype(n)>;
  I b{0};
  const I sz = ceil_div<I>(std::size(r) - o, n); // before moving
  return std::ranges::transform_view(std::ranges::iota_view{b, sz},
    [r = std::forward<R>(r), n, o](
      I i) -> decltype(auto) { return r[i * n + o]; });
} // stride_view

template<auto S, typename T>
auto
make_ids(T && t) {
  return std::ranges::transform_view(
    std::forward<T>(t), [](auto const & i) { return flecsi::topo::id<S>(i); });
} // make_ids
} // namespace util

/*!
  A node-centered (vertices) Cartesian mesh specialization.
 */
struct mesh : flecsi::topo::specialization<flecsi::topo::narray, mesh> {

  /*--------------------------------------------------------------------------*
    Policy Information.
   *--------------------------------------------------------------------------*/

  enum index_space { vertices };
  using index_spaces = has<vertices>;

  /// Mesh domains.
  /// The domain identifies the supported iteration spaces on the mesh.
  enum domain {
    /// This domain includes the DOFs of the problem. It does not include the
    /// dirichlet boundary vertices.
    interior,
    /// This domain includes the "logical" vertices, i.e., those that exist in a
    /// raw partitioning of the domain into non-intersecting sets of vertices.
    /// It does not include ghost vertices.
    logical,
    /// This domain includes all local vertices, including ghost cells.
    all,
    /// This domain includes all global cells.
    global
  };

  /// Mesh axes.
  /// The axis identifies a Cartesian coordinate axis on the mesh.
  enum axis {
    /// X-coordinate axis.
    x_axis,
    /// Y-coordinate axis.
    y_axis
  };

  using axes = has<x_axis, y_axis>;

  /// Boundary.
  /// Identifies the low or high boundary along a particular axis.
  enum boundary {
    /// Low axis boundary.
    low,
    /// High axis boundary.
    high
  };

  using coord = base::coord;
  using gcoord = base::gcoord;
  using colors = base::colors;
  using hypercube = base::hypercube;
  using axis_definition = base::axis_definition;
  using index_definition = base::index_definition;

  struct meta_data {
    double xdelta;
    double ydelta;
  };

  static constexpr std::size_t dimension = 2;

  template<auto>
  static constexpr std::size_t privilege_count = 2;

  /*--------------------------------------------------------------------------*
    Interface.
   *--------------------------------------------------------------------------*/

  /// Mesh Interface.
  template<class B>
  struct interface : B {

    /// Return the size for the given axis and domain.
    /// @tparam A  The mesh axis.
    /// @tparam DM The mesh domain.
    template<axis A, domain DM = interior>
    std::size_t size() {
      if constexpr(DM == interior) {
        const bool low = B::template is_low<mesh::vertices, A>();
        const bool high = B::template is_high<mesh::vertices, A>();

        if(low && high) { /* degenerate */
          return size<A, logical>() - 2;
        }
        else if(low || high) {
          return size<A, logical>() - 1;
        }
        else { /* interior */
          return size<A, logical>();
        }
      }
      else if constexpr(DM == logical) {
        return B::template size<mesh::vertices, A, base::domain::logical>();
      }
      else if constexpr(DM == all) {
        return B::template size<mesh::vertices, A, base::domain::all>();
      }
      else if constexpr(DM == global) {
        return B::template size<mesh::vertices, A, base::domain::global>();
      }
    } // size

    template<axis A>
    FLECSI_INLINE_TARGET std::size_t global_id(std::size_t i) const {
      return B::template global_id<mesh::vertices, A>(i);
    } // global_id

    /// Return a range over the given axis and domain.
    /// @tparam A  The mesh axis.
    /// @tparam DM The mesh domain (excluding \em global).
    ///
    /// The range can be used to iterate over the cells in the given domain,
    /// e.g.:
    /// @code
    ///   for(auto j: m.vertices<mesh::y_axis, mesh::interior>()) {
    ///     for(auto i: m.vertices<mesh::x_axis, mesh::interior>()) {
    ///       u[j][i] = 1.0;
    ///     } // for
    ///   } // for
    /// @endcode
    template<axis A, domain DM = interior, bool R = false>
    FLECSI_INLINE_TARGET auto vertices() const {
      flecsi::util::id b, e;

      if constexpr(DM == interior) {
        // The outermost layer is either halo or fixed boundaries:
        b = 1;
        e = B::template size<mesh::vertices, A, base::domain::all>() - 1;
      }
      else if constexpr(DM == logical) {
        b = B::template offset<mesh::vertices, A, base::domain::logical>();
        e = b + B::template size<mesh::vertices, A, base::domain::logical>();
      }
      else if constexpr(DM == all) {
        b = 0;
        e = B::template size<mesh::vertices, A, mesh::domain::all>();
      }
      else if(DM == global) {
        flog_fatal("illegal domain: you cannot iterate over the global domain");
      } // if

      if constexpr(R) {
        return util::make_ids<mesh::vertices>(
          std::ranges::iota_view{b, e} | std::views::reverse);
      }
      else {
        return util::make_ids<mesh::vertices>(std::ranges::iota_view{b, e});
      }
    } // vertices

    template<axis A>
    FLECSI_INLINE_TARGET auto red(std::size_t row) const {
      // The checkerboard extends across colors.  The (boundary) point with
      // global ID (0,0) is red; row is local, and 0 in the space of the
      // stride_view is the first interior vertex.
      return util::stride_view(vertices<A>(),
        2,
        (global_id<(A == x_axis ? y_axis : x_axis)>(row) + global_id<A>(1)) %
          2);
    }

    template<axis A>
    FLECSI_INLINE_TARGET auto black(std::size_t row) const {
      return red<A>(row + 1);
    }

    void set_geometry(double x, double y) { // available if writable
      this->policy_meta() = {x, y};
    }

    FLECSI_INLINE_TARGET double xdelta() const {
      return this->policy_meta().xdelta;
    }

    FLECSI_INLINE_TARGET double ydelta() const {
      return this->policy_meta().ydelta;
    }

    FLECSI_INLINE_TARGET double dxdy() const {
      return xdelta() * ydelta();
    }

    template<axis A>
    FLECSI_INLINE_TARGET double value(std::size_t i) const {
      return (A == x_axis ? xdelta() : ydelta()) * global_id<A>(i);
    }

    template<axis A, boundary BD>
    bool is_boundary(std::size_t i) {

      auto const loff =
        B::template offset<mesh::vertices, A, base::domain::logical>();
      auto const lsize =
        B::template size<mesh::vertices, A, base::domain::logical>();
      const bool l = B::template is_low<mesh::vertices, A>();
      const bool h = B::template is_high<mesh::vertices, A>();

      if(l && h) { /* degenerate */
        if constexpr(BD == boundary::low) {
          return i == loff;
        }
        else {
          return i == (lsize + loff - 1);
        }
      }
      else if(l) {
        if constexpr(BD == boundary::low) {
          return i == loff;
        }
        else {
          return false;
        }
      }
      else if(h) {
        if constexpr(BD == boundary::low) {
          return false;
        }
        else {
          return i == (lsize + loff - 1);
        }
      }
      else { /* interior */
        return false;
      }
    } // is_boundary

    template<class Sten, class MDColex>
    struct stencil_operator {
      constexpr double & operator()(std::size_t i, std::size_t j, Sten dir) {
        return so(i, j)[static_cast<std::ptrdiff_t>(dir)];
      }

      constexpr double
      operator()(std::size_t i, std::size_t j, Sten dir) const {
        return so(i, j)[static_cast<std::ptrdiff_t>(dir)];
      }

      MDColex so;
    }; // stencil_operator

    template<index_space S, class Sten, class T, flecsi::Privileges P>
    auto stencil_op(
      const flecsi::data::accessor<flecsi::data::dense, T, P> & a) const {
      return stencil_operator<Sten,
        std::decay_t<decltype(this->template mdcolex<S>(a))>>{
        this->template mdcolex<S>(a)};
    } // stencil_op
  }; // struct interface

  static auto distribute(flecsi::Color nc, gcoord indices) {
    return base::distribute(nc, indices);
  }

  /*--------------------------------------------------------------------------*
    Color Method.
   *--------------------------------------------------------------------------*/

  static coloring color(std::size_t num_colors, gcoord axis_extents) {
    index_definition idef;
    idef.axes = base::make_axes(num_colors, axis_extents);
    for(auto & a : idef.axes) {
      a.hdepth = 1;
    }

    return {{idef}};
  } // color

  static coloring color(colors const & parts, gcoord axis_extents) {
    index_definition idef;
    idef.axes = base::make_axes(parts, axis_extents);
    for(auto & a : idef.axes) {
      a.hdepth = 1;
    }

    return {{idef}};
  }

  /*--------------------------------------------------------------------------*
    Initialization.
   *--------------------------------------------------------------------------*/

  using grect = std::array<std::array<double, 2>, 2>;

  static void set_geometry(mesh::accessor<flecsi::rw> sm, grect const & g) {
    sm.set_geometry(
      std::abs(g[0][1] - g[0][0]) / (sm.size<x_axis, global>() - 1),
      std::abs(g[1][1] - g[1][0]) / (sm.size<y_axis, global>() - 1));
  }

  static void initialize(flecsi::data::topology_slot<mesh> & s,
    coloring const &,
    grect const & geometry) {
    flecsi::execute<set_geometry, flecsi::mpi>(s, geometry);
  } // initialize

}; // struct mesh

} // namespace fdm

#endif // FDM_MESH_HH
