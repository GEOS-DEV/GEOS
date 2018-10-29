
#
# Specify all headers
#
set(fileIO_headers
  xmlWrapper.hpp
  silo/SiloFile.hpp
  vtm/VtmFile.hpp
  blueprint/Blueprint.hpp
  schema/SchemaUtilities.hpp
  utils/utils.hpp )

#
# Specify all sources
#
set(fileIO_sources
  xmlWrapper.cpp
  silo/SiloFile.cpp
  vtm/VtmFile.cpp
  schema/SchemaUtilities.cpp
  blueprint/Blueprint.cpp )

set( dependencyList common silo hdf5 mpi pugixml RAJA )
if( CONDUIT_FOUND )
  set( dependencyList ${dependencyList} conduit conduit_relay conduit_blueprint )
endif()

blt_add_library( NAME                  fileIO
                 SOURCES               ${fileIO_sources}
                 HEADERS               ${fileIO_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( fileIO PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)


geosx_add_code_checks(PREFIX fileIO )
