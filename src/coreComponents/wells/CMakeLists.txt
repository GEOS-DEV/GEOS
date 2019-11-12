#
# Specify all headers
#
set( wells_headers
     InternalWellGenerator.hpp
     LASWellGenerator.hpp
     WellControls.hpp
     WellElementRegion.hpp
     WellElementSubRegion.hpp
     WellGeneratorBase.hpp
     PerforationData.hpp
     Perforation.hpp
   )
#
# Specify all sources
#
set( wells_sources
     InternalWellGenerator.cpp
     LASWellGenerator.cpp
     WellControls.cpp
     WellElementRegion.cpp
     WellElementSubRegion.cpp
     WellGeneratorBase.cpp
     PerforationData.cpp
     Perforation.cpp
   )

set( dependencyList dataRepository codingUtilities managers fileIO )

if ( ENABLE_OPENMP )
   set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
   set( dependencyList ${dependencyList} cuda )
endif()
   

blt_add_library( NAME                  wells
                 SOURCES               ${wells_sources}
                 HEADERS               ${wells_headers}
                 DEPENDS_ON            ${dependencyList}
                 #OBJECT                TRUE
               )
               
target_include_directories( wells PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX wells )
