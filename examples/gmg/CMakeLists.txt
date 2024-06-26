#------------------------------------------------------------------------------#
# Yaml
#------------------------------------------------------------------------------#

find_package(yaml-cpp REQUIRED)

#------------------------------------------------------------------------------#
# Add an executable.
#------------------------------------------------------------------------------#

add_executable(gmg
  src/final.cc
  src/analyze.cc
  src/gmg.cc
  src/init.cc
  src/norm.cc
  src/mg.cc
  src/solve.cc
  src/tasks/comps.cc
  src/tasks/analyze.cc
  src/tasks/init.cc
  src/tasks/io.cc
  src/tasks/norm.cc
)

#------------------------------------------------------------------------------#
# Add the FleCSI libraries and dependencies.
#------------------------------------------------------------------------------#

target_link_libraries(gmg FleCSI::FleCSI
  FleCSI-SP-FDM::FleCSI-SP-FDM
  yaml-cpp
)

#------------------------------------------------------------------------------#
# Unit tests.
#------------------------------------------------------------------------------#
#fdm_add_unit(mytest
#  SOURCES
#    mytest.cc
#)
