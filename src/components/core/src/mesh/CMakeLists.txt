
#
# Specify all headers
#
set(mesh_headers
    CellBlockManager.hpp
    CellBlock.hpp
    CellBlockSubRegion.hpp
    FaceManager.hpp
    NodeManager.hpp
    ElementRegion.hpp
    ElementRegionManager.hpp
    MeshBody.hpp
    MeshLevel.hpp    
  )


#
# Specify all sources
#
set(mesh_sources
    CellBlockManager.cpp
    CellBlock.cpp
    CellBlockSubRegion.cpp
    FaceManager.cpp
    NodeManager.cpp
    ElementRegion.cpp
    ElementRegionManager.cpp
    MeshBody.cpp
    MeshLevel.cpp    
   )

blt_add_library( NAME                  mesh
                 SOURCES               ${mesh_sources}
                 HEADERS               ${mesh_headers}
                 DEPENDS_ON            managers
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
geosx_add_code_checks(PREFIX mesh )