#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#

add_executable(poisson
  src/poisson.cc
  src/analyze.cc
  src/finalize.cc
  src/initialize.cc
  src/problem.cc
  src/solve.cc
  src/tasks/init.cc
  src/tasks/io.cc
  src/tasks/norm.cc
  src/tasks/smooth.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#

target_link_libraries(poisson FleCSI::FleCSI
  FleCSI-SP-Cartesian::FleCSI-SP-Cartesian
)
