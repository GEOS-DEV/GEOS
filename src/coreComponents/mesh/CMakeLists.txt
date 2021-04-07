# Specify all sources
set(mesh_sources
    AggregateElementSubRegion.cpp
    BufferOps.cpp
    CellBlockManager.cpp
    CellBlock.cpp
    CellElementRegion.cpp
    CellElementSubRegion.cpp
    EdgeManager.cpp
    ElementSubRegionBase.cpp
    ElementRegionBase.cpp
    ElementRegionManager.cpp
    SurfaceElementRegion.cpp
    SurfaceElementSubRegion.cpp
    FaceElementSubRegion.cpp
    EmbeddedSurfaceSubRegion.cpp
    WellElementRegion.cpp
    WellElementSubRegion.cpp
    FaceManager.cpp
    MeshBody.cpp
    MeshLevel.cpp
    NodeManager.cpp
    ToElementRelation.cpp
   )

if( BUILD_OBJ_LIBS)
  set(dependencyList managers metis)
else()
  set(dependencyList common metis)
endif()


if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  mesh
                 SOURCES               ${mesh_sources}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX mesh )
