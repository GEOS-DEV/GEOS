#
# Specify all headers
#
set( wells_headers
     InternalWellGenerator.hpp
     WellControls.hpp
     WellElementRegion.hpp
     WellElementSubRegion.hpp
     PerforationData.hpp
     Perforation.hpp
   )
#
# Specify all sources
#
set( wells_sources
     InternalWellGenerator.cpp
     WellControls.cpp
     WellElementRegion.cpp
     WellElementSubRegion.cpp
     PerforationData.cpp
     Perforation.cpp
   )

set( dependencyList codingUtilities managers )

blt_add_library( NAME                  wells
                 SOURCES               ${wells_sources}
                 HEADERS               ${wells_headers}
                 DEPENDS_ON            dataRepository ${dependencyList}
                 OBJECT                TRUE
               )
               
target_include_directories( wells PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX wells )
