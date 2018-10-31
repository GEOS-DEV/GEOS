message( "Entering src/coreComponents/physicsSolvers/CMakeLists.txt")

#
# Specify all headers
#
set(physicsSolvers_headers
    PhysicsSolverManager.hpp
    SolverBase.hpp
    CoupledSolvers/PoroelasticSolver.hpp
    FiniteVolume/SinglePhaseFlow.hpp
    SimpleSolvers/LaplaceFEM.hpp 
    src/SolidMechanicsLagrangianFEM.hpp
    surfaceGeneration/SurfaceGenerator.hpp
    src/DummySolver.hpp
    )

#
# Specify all sources
#
set(physicsSolvers_sources
    PhysicsSolverManager.cpp
    SolverBase.cpp
    CoupledSolvers/PoroelasticSolver.cpp
    FiniteVolume/SinglePhaseFlow.cpp
    SimpleSolvers/LaplaceFEM.cpp
    src/SolidMechanicsLagrangianFEM.cpp
    surfaceGeneration/SurfaceGenerator.cpp
    src/DummySolver.cpp
    )

set( externalComponentDeps "" )
if( ENABLE_GEOSX_PTP )
  message("adding GEOSX_PTP to physics solver dependencies")
  list( APPEND externalComponentDeps GEOSX_PTP )  
endif()

blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            dataRepository systemSolverInterface ${externalComponentDeps}
               )
               
target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)
target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents)

geosx_add_code_checks(PREFIX physicsSolvers )

add_subdirectory( unitTests )
message( "Leaving src/coreComponents/physicsSolvers/CMakeLists.txt")
