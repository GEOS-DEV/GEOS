
#
# Specify all headers
#
set(geosx_headers
  ${geosx_headers}
    mesh/CellBlockManager.hpp
    mesh/CellBlock.hpp
    mesh/CellBlockSubRegion.hpp
    mesh/FaceManager.hpp
    mesh/NodeManager.hpp
    mesh/ElementRegion.hpp
    mesh/ElementRegionManager.hpp
    mesh/MeshBody.hpp
    mesh/MeshLevel.hpp    
    PARENT_SCOPE )


#
# Specify all sources
#
set(geosx_core_sources
  ${geosx_core_sources}
    mesh/CellBlockManager.cpp
    mesh/CellBlock.cpp
    mesh/CellBlockSubRegion.cpp
    mesh/FaceManager.cpp
    mesh/NodeManager.cpp
    mesh/ElementRegion.cpp
    mesh/ElementRegionManager.cpp
    mesh/MeshBody.cpp
    mesh/MeshLevel.cpp    
    PARENT_SCOPE )