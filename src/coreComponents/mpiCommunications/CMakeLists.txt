
#
# Specify all headers
#
set(mpiCommunications_headers
    CommunicationTools.hpp
    MpiWrapper.hpp
    NeighborCommunicator.hpp
    PartitionBase.hpp
    SpatialPartition.hpp
    NeighborData.hpp
   )


#
# Specify all sources
#
set(mpiCommunications_sources
    CommunicationTools.cpp
    MpiWrapper.cpp
    NeighborCommunicator.cpp
    PartitionBase.cpp
    SpatialPartition.cpp
   )

if( BUILD_OBJ_LIBS)
  set( dependencyList fileIO common trilinos )
else()
  set( dependencyList 
       lvarray
       common
       mesh
       managers
       codingUtilities
       dataRepository
        )
endif()



if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if( ENABLE_MPI )
    set( dependencyList ${dependencyList} mpi )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  mpiCommunications
                 SOURCES               ${mpiCommunications_sources}
                 HEADERS               ${mpiCommunications_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( mpiCommunications PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX mpiCommunications )

add_subdirectory( unitTests )
