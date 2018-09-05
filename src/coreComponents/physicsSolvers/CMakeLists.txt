
#
# Specify all headers
#
set(physicsSolvers_headers
    PhysicsSolverManager.hpp
    SolverBase.hpp
    FiniteVolume/SinglePhaseFlow.hpp
    SimpleSolvers/LaplaceFEM.hpp
    src/SolidMechanicsLagrangianFEM.hpp
    )

#
# Specify all sources
#
set(physicsSolvers_sources
    PhysicsSolverManager.cpp
    SolverBase.cpp
    FiniteVolume/SinglePhaseFlow.cpp
    SimpleSolvers/LaplaceFEM.cpp
    src/SolidMechanicsLagrangianFEM.cpp
    )


blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            dataRepository systemSolverInterface
               )
               
target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX physicsSolvers )

add_subdirectory( unitTests )