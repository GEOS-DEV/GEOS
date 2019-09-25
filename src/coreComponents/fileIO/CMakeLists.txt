

#
# Specify all headers
#
set(fileIO_headers
  silo/SiloFile.hpp
  las/LASFile.hpp
  vtk/VTKFile.hpp
  vtm/VtmFile.hpp
  blueprint/Blueprint.hpp
  schema/SchemaUtilities.hpp
  utils/utils.hpp )

#
# Specify all sources
#
set(fileIO_sources
  silo/SiloFile.cpp
  las/LASFile.cpp
  vtk/VTKFile.cpp
  vtm/VtmFile.cpp
  schema/SchemaUtilities.cpp
  blueprint/Blueprint.cpp )

set( dependencyList common silo hdf5 pugixml RAJA  )

if( ENABLE_MPI )
  add_subdirectory(coupling/hdf5_interface)
  set( dependencyList ${dependencyList} mpi hdf5_interface)
  list( APPEND fileIO_headers coupling/ChomboCoupler.hpp )
  list( APPEND fileIO_sources coupling/ChomboCoupler.cpp )
endif()

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
