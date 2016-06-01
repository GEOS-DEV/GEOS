
#add_subdirectory(main)
#add_subdirectory(codingUtilities)



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



blt_add_executable(NAME main
                   SOURCES ${geosx_sources}
                   DEPENDS_ON slic )
                   
blt_add_library( NAME
                     geosx_core
                 SOURCES
                     "${geosx_core_sources}"
                 HEADERS
                     "${geosx_headers}"
                 DEPENDS_ON
                     slic
                 )