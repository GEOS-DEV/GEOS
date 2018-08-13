#
# Specify all headers
#
set( finiteVolume_headers
     FiniteVolumeManager.hpp
     FluxApproximationBase.hpp
     TwoPointFluxApproximation.hpp
     StencilCollection.hpp
   )
#
# Specify all sources
#
set( finiteVolume_sources
     FiniteVolumeManager.cpp
     FluxApproximationBase.cpp
     TwoPointFluxApproximation.cpp
   )

blt_add_library( NAME                  finiteVolume
                 SOURCES               ${finiteVolume_sources}
                 HEADERS               ${finiteVolume_headers}
                 DEPENDS_ON            dataRepository codingUtilities managers
               )
               
target_include_directories( finiteVolume PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src )
geosx_add_code_checks( PREFIX finiteVolume )