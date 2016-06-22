
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  	PhysicsSolvers/SolverBase.hpp
  	PhysicsSolvers/NewtonianMechanics.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    PhysicsSolvers/SolverBase.cpp
    PhysicsSolvers/NewtonianMechanics.cpp
    PARENT_SCOPE )
