#
# Specify all headers
#
set( wells_headers
     WellManager.hpp
     Well.hpp
     WellSolverBase.hpp
     SinglePhaseWell.hpp
     CompositionalMultiphaseWell.hpp
     Segment.hpp
     SegmentManager.hpp
     Connection.hpp
     ConnectionData.hpp
     Perforation.hpp
     PerforationData.hpp
   )
#
# Specify all sources
#
set( wells_sources
     WellManager.cpp
     Well.cpp
     WellSolverBase.cpp
     SinglePhaseWell.cpp
     CompositionalMultiphaseWell.cpp
     Segment.cpp
     SegmentManager.cpp
     Connection.cpp
     ConnectionData.cpp
     Perforation.cpp
     PerforationData.cpp
   )

blt_add_library( NAME                  wells
                 SOURCES               ${wells_sources}
                 HEADERS               ${wells_headers}
                 DEPENDS_ON            dataRepository codingUtilities managers
               )
               
target_include_directories( wells PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX wells )
