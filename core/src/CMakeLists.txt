
#add_subdirectory(main)
#add_subdirectory(codingUtilities)

#subdirlist(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
#foreach(subdir ${SUBDIRS})
#  if (EXISTS subdir/"CMakeLists.txt")
#    add_subdirectory(${subdir})
#  endif()
#endforeach()


#
# Specify all headers
#
set(geosx_headers
    codingUtilities/stackTrace.hpp
    )

#
# Specify all sources
#

set(geosx_core_sources
    codingUtilities/stackTrace.cpp
    )
    
set(geosx_sources
	${geosx_core_sources}
	main/main.cpp
    )


message( "${CMAKE_CURRENT_LIST_DIR}")
blt_add_executable(NAME main
                   SOURCES ${geosx_sources}
                   DEPENDS_ON slic )

foreach(header ${geosx_headers})
  message( "${CMAKE_CURRENT_LIST_DIR}/${header}")
  set( geosx_headers2 
       ${CMAKE_CURRENT_LIST_DIR}/${header}
     )
endforeach()                   

blt_add_library( NAME
                     geosx_core
                 SOURCES
                     "${geosx_core_sources}"
                 HEADERS
                     "${geosx_headers2}"
                 DEPENDS_ON
                     slic
                 )