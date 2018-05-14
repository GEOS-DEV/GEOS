
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    MPI_Communications/CommunicationTools.hpp
    MPI_Communications/NeighborCommunicator.hpp
    MPI_Communications/NeighborCommunication.hpp
    MPI_Communications/PartitionBase.hpp
    MPI_Communications/SpatialPartition.hpp
    PARENT_SCOPE )


#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    MPI_Communications/CommunicationTools.cpp
    MPI_Communications/NeighborCommunicator.cpp
    MPI_Communications/NeighborCommunication.cpp
    MPI_Communications/PartitionBase.cpp
    MPI_Communications/SpatialPartition.cpp
    PARENT_SCOPE )