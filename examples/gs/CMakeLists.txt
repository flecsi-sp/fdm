#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#

add_executable(gs
  src/gs.cc
  src/analyze.cc
  src/final.cc
  src/init.cc
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

target_link_libraries(gs FleCSI::FleCSI
  FleCSI-SP-FDM::FleCSI-SP-FDM
)
