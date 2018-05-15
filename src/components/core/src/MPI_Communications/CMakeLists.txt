
#
# Specify all headers
#
set(MPI_Communications_headers
    CommunicationTools.hpp
    NeighborCommunicator.hpp
    NeighborCommunication.hpp
    PartitionBase.hpp
    SpatialPartition.hpp
   )


#
# Specify all sources
#
set(MPI_Communications_sources
    CommunicationTools.cpp
    NeighborCommunicator.cpp
    NeighborCommunication.cpp
    PartitionBase.cpp
    SpatialPartition.cpp
   )
    
set( dependencyList fileIO common  )
set( dependencyList2 trilinos )
foreach( lib ${dependencyList} )
    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
endforeach()
message( "dependencyList2 = ${dependencyList2}" )

blt_add_library( NAME                  MPI_Communications
                 SOURCES               ${MPI_Communications_sources}
                 HEADERS               ${MPI_Communications_headers}
                 DEPENDS_ON            ${dependencyList2}
                 SHARED FALSE
               )
               
target_include_directories( MPI_Communications PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX MPI_Communications )