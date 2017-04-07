
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    MPI_Communications/NeighborCommunication.hpp
    MPI_Communications/PartitionBase.hpp
    MPI_Communications/SpatialPartition.hpp
    PARENT_SCOPE )


#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    MPI_Communications/NeighborCommunication.cpp
    MPI_Communications/PartitionBase.cpp
    MPI_Communications/SpatialPartition.cpp
    PARENT_SCOPE )

if( 0 )
set(geosx_headers
  ${geosx_headers}
    managers/EdgeManager.hpp
    managers/ElementManager.hpp
    managers/ElementRegionT.hpp
    managers/FaceManager.hpp
    PARENT_SCOPE )
    
set(geosx_core_sources
  ${geosx_core_sources}
    managers/EdgeManager.cpp
    managers/ElementManager.cpp
    managers/ElementRegionT.cpp
    managers/FaceManager.cpp
    PARENT_SCOPE )
endif()