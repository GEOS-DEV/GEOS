
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
  	PhysicsSolvers/SolverBase.hpp
  	PhysicsSolvers/SolidMechanicsLagrangianFEM.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    PhysicsSolvers/SolverBase.cpp
    PhysicsSolvers/SolidMechanicsLagrangianFEM.cpp
    PARENT_SCOPE )
