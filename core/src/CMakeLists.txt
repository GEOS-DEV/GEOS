
#
# Specify all headers
#
set(geosx_headers "" )

#
# Specify all sources
#

set(geosx_core_sources "" )
    
add_subdirectory(DataStructures)
add_subdirectory(codingUtilities)

message("geosx_core_sources = ${geosx_core_sources}")

blt_add_library( NAME
                     geosx_core
                 SOURCES
                     "${geosx_core_sources}"
                 HEADERS
                     "${geosx_headers}"
                 DEPENDS_ON
                     slic sidre conduit
                 )
                 
add_subdirectory(main)
                 