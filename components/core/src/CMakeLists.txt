
#
# Specify all headers
#
set(geosx_headers "" )

#
# Specify all sources
#

set(geosx_core_sources "" )
    
add_subdirectory(dataRepository)
add_subdirectory(codingUtilities)
add_subdirectory(PhysicsSolvers)

blt_add_library( NAME
                     geosx_core
                 SOURCES
                     "${geosx_core_sources}"
                 HEADERS
                     "${geosx_headers}"
                 DEPENDS_ON
                     slic sidre conduit hdf5
                 )
                 
add_subdirectory(main)

blt_add_executable(NAME geosx_corex
                   SOURCES "${geosx_core_sources}"
                   DEPENDS_ON slic sidre conduit hdf5)