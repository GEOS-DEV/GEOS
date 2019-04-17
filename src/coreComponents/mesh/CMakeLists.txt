
#
# Specify all headers
#
set(mesh_headers
    AggregateElementSubRegion.hpp
    BufferOps.hpp
    CellBlockManager.hpp
    CellBlock.hpp
    CellElementSubRegion.hpp
    WellElementSubRegion.hpp
    EdgeManager.hpp
    ElementSubRegionBase.hpp
    ElementRegion.hpp
    ElementRegionManager.hpp
    FaceElementSubRegion.hpp
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
    AggregateElementSubRegion.cpp
    BufferOps.cpp
    CellBlockManager.cpp
    CellBlock.cpp
    CellElementSubRegion.cpp
    WellElementSubRegion.cpp
    EdgeManager.cpp
    ElementSubRegionBase.cpp
    ElementRegion.cpp
    ElementRegionManager.cpp
    FaceElementSubRegion.cpp
    FaceManager.cpp
    MeshBody.cpp
    MeshLevel.cpp    
    NodeManager.cpp
    ToElementRelation.cpp
   )

blt_add_library( NAME                  mesh
                 SOURCES               ${mesh_sources}
                 HEADERS               ${mesh_headers}
                 DEPENDS_ON            managers metis
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX mesh )
