#
# Specify all headers
#
set( wells_headers
     WellManager.hpp
     WellBase.hpp
     SimpleWell.hpp
     SinglePhaseWell.hpp
     CompositionalMultiphaseWell.hpp
     Segment.hpp
     SegmentManager.hpp
     Connection.hpp
     ConnectionManager.hpp
     Perforation.hpp
     PerforationManager.hpp
   )
#
# Specify all sources
#
set( wells_sources
     WellManager.cpp
     WellBase.cpp
     SimpleWell.cpp
     SinglePhaseWell.cpp
     CompositionalMultiphaseWell.cpp
     Segment.cpp
     SegmentManager.cpp
     Connection.cpp
     ConnectionManager.cpp
     Perforation.cpp
     PerforationManager.cpp
   )

blt_add_library( NAME                  wells
                 SOURCES               ${wells_sources}
                 HEADERS               ${wells_headers}
                 DEPENDS_ON            dataRepository codingUtilities managers
               )
               
target_include_directories( wells PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX wells )
