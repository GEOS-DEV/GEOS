# Specify solver headers
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     SolverBaseKernels.hpp
     SolverStatistics.hpp
     FieldStatisticsBase.hpp
 )

# Specify solver sources
set( physicsSolvers_sources
     LinearSolverParameters.cpp
     NonlinearSolverParameters.cpp
     PhysicsSolverManager.cpp
     SolverBase.cpp
     SolverStatistics.cpp )

  if( GEOS_ENABLE_CONTACT )
    include( contact/CMakeLists.txt )
  endif()

  if( GEOS_ENABLE_FLUID_FLOW )
    include( fluidFlow/CMakeLists.txt )
  endif()

  if( GEOS_ENABLE_INDUCEDSEISMICITY )
    include( inducedSeismicity/CMakeLists.txt )
  endif()

  if( GEOS_ENABLE_MULTIPHYSICS )
    include( multiphysics/CMakeLists.txt )
    include( multiphysics/poromechanicsKernels/PoromechanicsKernels.cmake)
  endif()

  include( python/CMakeLists.txt )

  if( GEOS_ENABLE_SIMPLEPDE )
    include( simplePDE/CMakeLists.txt )
  endif()

  if( GEOS_ENABLE_SOLIDMECHANICS )
    include( solidMechanics/CMakeLists.txt )
    include( solidMechanics/kernels/SolidMechanicsKernels.cmake)
    endif()

  if( GEOS_ENABLE_SURFACEGENERATION )
    include( surfaceGeneration/CMakeLists.txt )
  endif()

  if( GEOS_ENABLE_WAVEPROPAGATION )
    include( wavePropagation/CMakeLists.txt )
  endif()

set( dependencyList ${parallelDeps} constitutive mesh linearAlgebra discretizationMethods events )

blt_add_library( NAME       physicsSolvers
                 SOURCES    ${physicsSolvers_sources}
                 HEADERS    ${physicsSolvers_headers}
                 DEPENDS_ON ${dependencyList} ${externalComponentDeps}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS} )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents )
endif()


