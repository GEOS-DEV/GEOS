
#
# Specify all headers
#
set(mesh_headers
    BufferOps.hpp
    CellBlockManager.hpp
    CellBlock.hpp
    CellBlockSubRegion.hpp
    EdgeManager.hpp
    ElementRegion.hpp
    ElementRegionManager.hpp
    FaceManager.hpp
    InterObjectRelation.hpp
    MeshBody.hpp
    MeshLevel.hpp
    NodeManager.hpp
    ToElementRelation.hpp
  )


#
# Specify all sources
#
set(mesh_sources
    BufferOps.cpp
    CellBlockManager.cpp
    CellBlock.cpp
    CellBlockSubRegion.cpp
    EdgeManager.cpp
    ElementRegion.cpp
    ElementRegionManager.cpp
    FaceManager.cpp
    MeshBody.cpp
    MeshLevel.cpp    
    NodeManager.cpp
    ToElementRelation.cpp
   )

blt_add_library( NAME                  mesh
                 SOURCES               ${mesh_sources}
                 HEADERS               ${mesh_headers}
                 DEPENDS_ON            managers
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX mesh )
