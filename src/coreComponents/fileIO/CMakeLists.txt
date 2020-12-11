

#
# Specify all headers
#
set( fileIO_headers
     timeHistory/TimeHistHDF.hpp
     silo/SiloFile.hpp
     schema/schemaUtilities.hpp )

#
# Specify all sources
#
set( fileIO_sources
     timeHistory/TimeHistHDF.cpp
     silo/SiloFile.cpp
     schema/schemaUtilities.cpp )

if( BUILD_OBJ_LIBS)
  set( dependencyList common silo hdf5 pugixml RAJA  )
else()
  set( dependencyList common silo pugixml RAJA  )
endif()

if( ENABLE_MPI )
  add_subdirectory(coupling/hdf5_interface)
  set( dependencyList ${dependencyList} mpi hdf5_interface)
  list( APPEND fileIO_headers coupling/ChomboCoupler.hpp )
  list( APPEND fileIO_sources coupling/ChomboCoupler.cpp )
endif() 

if( ENABLE_VTK )
  add_subdirectory( vtk )
  set( dependencyList ${dependencyList} fileIO_vtk )
endif()

if( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper adiak )
endif()

blt_add_library( NAME                  fileIO
                 SOURCES               ${fileIO_sources}
                 HEADERS               ${fileIO_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( fileIO PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX fileIO )

add_subdirectory(unitTests)
