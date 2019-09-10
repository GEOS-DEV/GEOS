
#
# Specify all headers
#
set(MPI_Communications_headers
    CommunicationTools.hpp
    NeighborCommunicator.hpp
    PartitionBase.hpp
    SpatialPartition.hpp
   )


#
# Specify all sources
#
set(MPI_Communications_sources
    CommunicationTools.cpp
    NeighborCommunicator.cpp
    PartitionBase.cpp
    SpatialPartition.cpp
   )

set( dependencyList fileIO common trilinos )

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  MPI_Communications
                 SOURCES               ${MPI_Communications_sources}
                 HEADERS               ${MPI_Communications_headers}
                 DEPENDS_ON            ${dependencyList}
                 #OBJECT                TRUE
               )
               
target_include_directories( MPI_Communications PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX MPI_Communications )