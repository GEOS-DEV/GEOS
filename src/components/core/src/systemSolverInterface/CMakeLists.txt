
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    systemSolverInterface/SystemSolverParameters.hpp
    systemSolverInterface/LinearSolverWrapper.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    systemSolverInterface/SystemSolverParameters.cpp
    systemSolverInterface/LinearSolverWrapper.cpp
    PARENT_SCOPE )
