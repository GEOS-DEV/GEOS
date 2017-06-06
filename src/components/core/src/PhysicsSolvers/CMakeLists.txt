
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    PhysicsSolvers/PhysicsSolverManager.hpp
  	PhysicsSolvers/SolverBase.hpp
  	PhysicsSolvers/BoundaryConditions/BoundaryConditionBase.hpp
    PhysicsSolvers/BoundaryConditions/DirichletBoundaryCondition.hpp
    PARENT_SCOPE )

#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    PhysicsSolvers/PhysicsSolverManager.cpp
    PhysicsSolvers/SolverBase.cpp
    PhysicsSolvers/BoundaryConditions/BoundaryConditionBase.cpp
    PhysicsSolvers/BoundaryConditions/DirichletBoundaryCondition.cpp
    PARENT_SCOPE )
