#
# Specify all headers
#
set( finiteVolume_headers
     CellElementStencilMPFA.hpp
     CellElementStencilTPFA.hpp
     FaceElementStencil.hpp
     FiniteVolumeManager.hpp
     FluxApproximationBase.hpp
     TwoPointFluxApproximation.hpp
     FluxStencil.hpp
   )
#
# Specify all sources
#
set( finiteVolume_sources
     CellElementStencilMPFA.cpp
     CellElementStencilTPFA.cpp
     FaceElementStencil.cpp
     FiniteVolumeManager.cpp
     FluxApproximationBase.cpp
     TwoPointFluxApproximation.cpp
   )

set( dependencyList dataRepository codingUtilities managers )

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  finiteVolume
                 SOURCES               ${finiteVolume_sources}
                 HEADERS               ${finiteVolume_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                TRUE
               )
               
target_include_directories( finiteVolume PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX finiteVolume )
