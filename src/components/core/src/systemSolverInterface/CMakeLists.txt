
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    systemSolverInterface/EpetraBlockSystem.hpp
    systemSolverInterface/LinearSolverWrapper.hpp
    systemSolverInterface/SystemSolverParameters.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    systemSolverInterface/EpetraBlockSystem.cpp
    systemSolverInterface/SystemSolverParameters.cpp
    systemSolverInterface/LinearSolverWrapper.cpp
    PARENT_SCOPE )
