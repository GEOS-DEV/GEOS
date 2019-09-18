
add_subdirectory(coupling/hdf5_interface)

#
# Specify all headers
#
set(fileIO_headers
  silo/SiloFile.hpp
  vtk/VTKFile.hpp
  vtm/VtmFile.hpp
  blueprint/Blueprint.hpp
  schema/SchemaUtilities.hpp
  coupling/ChomboCoupler.hpp
  utils/utils.hpp )

#
# Specify all sources
#
set(fileIO_sources
  silo/SiloFile.cpp
  vtk/VTKFile.cpp
  vtm/VtmFile.cpp
  schema/SchemaUtilities.cpp
  blueprint/Blueprint.cpp
  coupling/ChomboCoupler.cpp )

set( dependencyList common silo hdf5 mpi pugixml RAJA hdf5_interface )
if( CONDUIT_FOUND )
  set( dependencyList ${dependencyList} conduit conduit_relay conduit_blueprint )
endif()

if( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  fileIO
                 SOURCES               ${fileIO_sources}
                 HEADERS               ${fileIO_headers}
                 DEPENDS_ON            ${dependencyList}
                 #OBJECT                TRUE
               )
               
target_include_directories( fileIO PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX fileIO )

add_subdirectory( unitTests )
