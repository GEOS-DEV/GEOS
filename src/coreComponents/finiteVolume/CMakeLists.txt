#
# Specify all headers
#
set( finiteVolume_headers
     FiniteVolumeManager.hpp
     FluxApproximationBase.hpp
     TwoPointFluxApproximation.hpp
     FluxStencil.hpp
   )
#
# Specify all sources
#
set( finiteVolume_sources
     FiniteVolumeManager.cpp
     FluxApproximationBase.cpp
     TwoPointFluxApproximation.cpp
   )

set( dependencyList dataRepository codingUtilities managers )

blt_add_library( NAME                  finiteVolume
                 SOURCES               ${finiteVolume_sources}
                 HEADERS               ${finiteVolume_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( finiteVolume PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX finiteVolume )
