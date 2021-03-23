# Specify all headers
set(mesh_headers
    AggregateElementSubRegion.hpp
    BufferOps.hpp
    CellBlockManager.hpp
    CellBlock.hpp
    CellElementRegion.hpp
    CellElementSubRegion.hpp
    EdgeManager.hpp
    ElementSubRegionBase.hpp
    ElementRegionBase.hpp
    ElementRegionManager.hpp
    EmbeddedSurfaceSubRegion.hpp
    ExtrinsicMeshData.hpp
    FaceElementSubRegion.hpp
    FaceManager.hpp
    InterObjectRelation.hpp
    MeshBody.hpp
    MeshLevel.hpp
    ObjectManagerBase.hpp
    NodeManager.hpp
    SurfaceElementRegion.hpp
    SurfaceElementSubRegion.hpp
    ToElementRelation.hpp
    WellElementSubRegion.hpp
    WellElementRegion.hpp
  )

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
    EmbeddedSurfaceSubRegion.cpp
    FaceElementSubRegion.cpp
    FaceManager.cpp
    MeshBody.cpp
    MeshLevel.cpp
    ObjectManagerBase.cpp
    NodeManager.cpp
    SurfaceElementRegion.cpp
    SurfaceElementSubRegion.cpp
    ToElementRelation.cpp
    WellElementRegion.cpp
    WellElementSubRegion.cpp
   )

if( BUILD_OBJ_LIBS)
  set(dependencyList managers metis)
else()
  set(dependencyList
      common
      metis
      codingUtilities
      meshUtilities
      mpiCommunications
      constitutive
      finiteElement
      dataRepository
      lvarray
      managers
      fileIO
      )
endif()


if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  mesh
                 SOURCES               ${mesh_sources}
                 HEADERS               ${mesh_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX mesh )
