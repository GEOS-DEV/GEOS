
#
# Specify all headers
#
set(fileIO_headers
  xmlWrapper.hpp
  silo/SiloFile.hpp
  blueprint/Blueprint.hpp
  schema/SchemaUtilities.hpp
  utils/utils.hpp )

#
# Specify all sources
#
set(fileIO_sources
  xmlWrapper.cpp
  silo/SiloFile.cpp
  schema/SchemaUtilities.cpp
  blueprint/Blueprint.cpp )

set( dependencyList common  )
set( dependencyList2 silo hdf5 mpi pugixml RAJA )
foreach( lib ${dependencyList} )
    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
endforeach()
message( "dependencyList2 = ${dependencyList2}" )

blt_add_library( NAME                  fileIO
                 SOURCES               ${fileIO_sources}
                 HEADERS               ${fileIO_headers}
                 DEPENDS_ON            ${dependencyList2}
                 SHARED FALSE
               )
               
target_include_directories( fileIO PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)

geosx_add_code_checks(PREFIX fileIO )