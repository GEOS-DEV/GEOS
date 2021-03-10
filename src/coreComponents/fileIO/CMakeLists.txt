

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
  set( dependencyList 
       lvarray
       common
       mpiCommunications
       managers
       codingUtilities
       dataRepository
       constitutive
       silo
       pugixml
       RAJA  )
endif()

if( ENABLE_MPI )
  add_subdirectory(coupling/hdf5_interface)
  set( dependencyList ${dependencyList} mpi hdf5_interface)
  list( APPEND fileIO_headers coupling/ChomboCoupler.hpp )
  list( APPEND fileIO_sources coupling/ChomboCoupler.cpp )
endif() 

if( ENABLE_VTK )
  list( APPEND fileIO_headers
               vtk/VTKPVDWriter.hpp
               vtk/VTKVTMWriter.hpp
               vtk/VTKGEOSXData.hpp
               vtk/VTKPolyDataWriterInterface.hpp )

  list( APPEND fileIO_sources
               vtk/VTKGEOSXData.cpp
               vtk/VTKPVDWriter.cpp
               vtk/VTKVTMWriter.cpp
               vtk/VTKPolyDataWriterInterface.cpp )


  list( APPEND dependencyList VTK::CommonCore VTK::CommonSystem VTK::CommonMisc VTK::IOCore VTK::IOXML VTK::expat VTK::doubleconversion VTK::lz4 VTK::lzma VTK::zlib VTK::loguru )
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
