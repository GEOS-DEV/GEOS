#
# Specify all headers
#
set( mesh_headers
     AggregateElementSubRegion.hpp
     BufferOps.hpp
     CellBlock.hpp
     CellBlockManager.hpp
     CellElementRegion.hpp
     CellElementSubRegion.hpp
     DomainPartition.hpp
     EdgeManager.hpp
     ElementRegionBase.hpp
     ElementRegionManager.hpp
     ElementSubRegionBase.hpp
     ElementType.hpp
     EmbeddedSurfaceNodeManager.hpp
     EmbeddedSurfaceSubRegion.hpp
     ExtrinsicMeshData.hpp
     FaceElementSubRegion.hpp
     FaceManager.hpp
     InterObjectRelation.hpp
     MeshBody.hpp
     MeshForLoopInterface.hpp
     MeshLevel.hpp
     MeshManager.hpp
     NodeManager.hpp
     ObjectManagerBase.hpp
     ParticleBlock.hpp
     ParticleBlockManager.hpp
     ParticleManager.hpp
     ParticleRegion.hpp
     ParticleRegionBase.hpp
     ParticleSubRegion.hpp
     ParticleSubRegionBase.hpp
     ParticleType.hpp
     Perforation.hpp
     PerforationData.hpp
     SurfaceElementRegion.hpp
     SurfaceElementSubRegion.hpp
     ToElementRelation.hpp
     WellElementRegion.hpp
     WellElementSubRegion.hpp
     generators/InternalMeshGenerator.hpp
     generators/InternalWellGenerator.hpp
     generators/InternalWellboreGenerator.hpp
     generators/MeshGeneratorBase.hpp
     generators/ParticleMeshGenerator.hpp
     mpiCommunications/CommID.hpp
     mpiCommunications/CommunicationTools.hpp
     mpiCommunications/MPI_iCommData.hpp
     mpiCommunications/NeighborCommunicator.hpp
     mpiCommunications/NeighborData.hpp
     mpiCommunications/PartitionBase.hpp
     mpiCommunications/SpatialPartition.hpp
     simpleGeometricObjects/BoundedPlane.hpp
     simpleGeometricObjects/Box.hpp
     simpleGeometricObjects/Cylinder.hpp
     simpleGeometricObjects/GeometricObjectManager.hpp
     simpleGeometricObjects/SimpleGeometricObjectBase.hpp
     simpleGeometricObjects/ThickPlane.hpp
     utilities/ComputationalGeometry.hpp
     utilities/MeshMapUtilities.hpp
     utilities/MeshUtilities.hpp
     utilities/StructuredGridUtilities.hpp
   )

#
# Specify all sources
#
set( mesh_sources
     AggregateElementSubRegion.cpp
     BufferOps.cpp
     CellBlock.cpp
     CellBlockManager.cpp
     CellElementRegion.cpp
     CellElementSubRegion.cpp
     DomainPartition.cpp
     EdgeManager.cpp
     ElementRegionBase.cpp
     ElementRegionManager.cpp
     ElementSubRegionBase.cpp
     EmbeddedSurfaceNodeManager.cpp
     EmbeddedSurfaceSubRegion.cpp
     FaceElementSubRegion.cpp
     FaceManager.cpp
     MeshBody.cpp
     MeshLevel.cpp
     MeshManager.cpp
     NodeManager.cpp
     ObjectManagerBase.cpp
     ParticleBlock.cpp
     ParticleBlockManager.cpp
     ParticleManager.cpp
     ParticleRegion.cpp
     ParticleRegionBase.cpp
     ParticleSubRegion.cpp
     ParticleSubRegionBase.cpp
     Perforation.cpp
     PerforationData.cpp
     SurfaceElementRegion.cpp
     SurfaceElementSubRegion.cpp
     ToElementRelation.cpp
     WellElementRegion.cpp
     WellElementSubRegion.cpp
     generators/InternalMeshGenerator.cpp
     generators/InternalWellGenerator.cpp
     generators/InternalWellboreGenerator.cpp
     generators/MeshGeneratorBase.cpp
     generators/ParticleMeshGenerator.cpp
     mpiCommunications/CommID.cpp
     mpiCommunications/CommunicationTools.cpp
     mpiCommunications/MPI_iCommData.cpp    
     mpiCommunications/NeighborCommunicator.cpp
     mpiCommunications/PartitionBase.cpp
     mpiCommunications/SpatialPartition.cpp
     simpleGeometricObjects/BoundedPlane.cpp
     simpleGeometricObjects/Box.cpp
     simpleGeometricObjects/Cylinder.cpp
     simpleGeometricObjects/GeometricObjectManager.cpp
     simpleGeometricObjects/SimpleGeometricObjectBase.cpp
     simpleGeometricObjects/ThickPlane.cpp
     utilities/ComputationalGeometry.cpp
     utilities/MeshUtilities.cpp
   )

set( dependencyList schema dataRepository constitutive metis )

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if( ENABLE_PAMELA )
    message( STATUS "Adding PAMELAMeshGenerator sources and headers" )
    list( APPEND mesh_headers generators/PAMELAMeshGenerator.hpp )
    list( APPEND mesh_sources generators/PAMELAMeshGenerator.cpp )
    add_subdirectory( PAMELA )
    list( APPEND dependencyList PAMELA )
endif()

blt_add_library( NAME       mesh
                 SOURCES    ${mesh_sources}
                 HEADERS    ${mesh_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks( PREFIX mesh 
                       EXCLUDES PAMELA )
