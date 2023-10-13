# Specify all headers
set( mesh_headers
     BufferOps.hpp
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
     MeshFields.hpp
     FaceElementSubRegion.hpp
     FaceManager.hpp
     FieldIdentifiers.hpp
     InterObjectRelation.hpp
     MeshBody.hpp
     MeshForLoopInterface.hpp
     MeshLevel.hpp
     MeshManager.hpp
     MeshObjectPath.hpp
     NodeManager.hpp
     ObjectManagerBase.hpp
     generators/ParticleBlock.hpp
     generators/ParticleBlockManager.hpp
     ParticleManager.hpp
     ParticleRegion.hpp
     ParticleRegionBase.hpp
     ParticleSubRegion.hpp
     ParticleSubRegionBase.hpp
     generators/ParticleBlockABC.hpp
     ParticleType.hpp
     Perforation.hpp
     PerforationData.hpp
     PerforationFields.hpp
     SurfaceElementRegion.hpp
     SurfaceElementSubRegion.hpp
     ToElementRelation.hpp
     ToParticleRelation.hpp
     WellElementRegion.hpp
     WellElementSubRegion.hpp
     generators/CellBlock.hpp
     generators/CellBlockABC.hpp
     generators/CellBlockManager.hpp
     generators/CellBlockManagerABC.hpp
     generators/CellBlockUtilities.hpp
     generators/LineBlock.hpp
     generators/LineBlockABC.hpp
     generators/ExternalMeshGeneratorBase.hpp
     generators/FaceBlock.hpp
     generators/FaceBlockABC.hpp
     generators/InternalMeshGenerator.hpp
     generators/InternalWellGenerator.hpp
     generators/InternalWellboreGenerator.hpp
     generators/MeshGeneratorBase.hpp
     generators/ParMETISInterface.hpp
     generators/ParticleMeshGenerator.hpp
     generators/PartitionDescriptor.hpp
     generators/PrismUtilities.hpp
     generators/WellGeneratorBase.hpp
     mpiCommunications/CommID.hpp
     mpiCommunications/CommunicationTools.hpp
     mpiCommunications/MPI_iCommData.hpp
     mpiCommunications/NeighborCommunicator.hpp
     mpiCommunications/NeighborData.hpp
     mpiCommunications/PartitionBase.hpp
     mpiCommunications/SpatialPartition.hpp
     simpleGeometricObjects/Rectangle.hpp
     simpleGeometricObjects/Disc.hpp
     simpleGeometricObjects/CustomPolarObject.hpp
     simpleGeometricObjects/Box.hpp
     simpleGeometricObjects/Cylinder.hpp
     simpleGeometricObjects/GeometricObjectManager.hpp
     simpleGeometricObjects/SimpleGeometricObjectBase.hpp
     simpleGeometricObjects/PlanarGeometricObject.hpp
     simpleGeometricObjects/ThickPlane.hpp
     utilities/AverageOverQuadraturePointsKernel.hpp     
     utilities/CIcomputationKernel.hpp
     utilities/ComputationalGeometry.hpp
     utilities/MeshMapUtilities.hpp
     utilities/StructuredGridUtilities.hpp
   )

# Specify all sources
set( mesh_sources
     BufferOps.cpp
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
     MeshObjectPath.cpp
     NodeManager.cpp
     ObjectManagerBase.cpp
     generators/ParticleBlock.cpp
     generators/ParticleBlockManager.cpp
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
     ToParticleRelation.cpp
     WellElementRegion.cpp
     WellElementSubRegion.cpp
     generators/CellBlock.cpp
     generators/CellBlockManager.cpp
     generators/CellBlockUtilities.cpp
     generators/LineBlock.cpp
     generators/ExternalMeshGeneratorBase.cpp
     generators/FaceBlock.cpp
     generators/InternalMeshGenerator.cpp
     generators/InternalWellGenerator.cpp
     generators/InternalWellboreGenerator.cpp
     generators/MeshGeneratorBase.cpp
     generators/ParMETISInterface.cpp
     generators/ParticleMeshGenerator.cpp
     generators/WellGeneratorBase.cpp
     mpiCommunications/CommID.cpp
     mpiCommunications/CommunicationTools.cpp
     mpiCommunications/MPI_iCommData.cpp
     mpiCommunications/NeighborCommunicator.cpp
     mpiCommunications/PartitionBase.cpp
     mpiCommunications/SpatialPartition.cpp
     simpleGeometricObjects/Rectangle.cpp
     simpleGeometricObjects/Disc.cpp
     simpleGeometricObjects/CustomPolarObject.cpp
     simpleGeometricObjects/Box.cpp
     simpleGeometricObjects/Cylinder.cpp
     simpleGeometricObjects/GeometricObjectManager.cpp
     simpleGeometricObjects/SimpleGeometricObjectBase.cpp
     simpleGeometricObjects/PlanarGeometricObject.cpp
     simpleGeometricObjects/ThickPlane.cpp
     utilities/ComputationalGeometry.cpp
     )

set( dependencyList ${parallelDeps} schema dataRepository constitutive finiteElement parmetis metis )

if( ENABLE_VTK )
    message(STATUS "Adding VTKMeshGenerator sources and headers")
    set( mesh_headers ${mesh_headers}
         generators/CollocatedNodes.hpp
         generators/VTKFaceBlockUtilities.hpp
         generators/VTKMeshGenerator.hpp
         generators/VTKMeshGeneratorTools.hpp
         generators/VTKUtilities.hpp
         )
    set( mesh_sources ${mesh_sources}
         generators/CollocatedNodes.cpp
         generators/VTKFaceBlockUtilities.cpp
         generators/VTKMeshGenerator.cpp
         generators/VTKMeshGeneratorTools.cpp
         generators/VTKUtilities.cpp
         )
    set( dependencyList ${dependencyList} VTK::IOLegacy VTK::FiltersParallelDIY2 )
    if( ENABLE_MPI )
      set( dependencyList ${dependencyList} VTK::IOParallelXML VTK::ParallelMPI )
    endif()
endif()

if( ENABLE_SCOTCH )
    set( mesh_headers ${mesh_headers} generators/PTScotchInterface.hpp )
    set( mesh_sources ${mesh_sources} generators/PTScotchInterface.cpp )
    set( dependencyList ${dependencyList} ptscotch )
endif()

blt_add_library( NAME       mesh
                 SOURCES    ${mesh_sources}
                 HEADERS    ${mesh_headers}
                 DEPENDS_ON ${dependencyList}
                 OBJECT     ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

if( ENABLE_VTK )
    # Avoid compiling with nvcc which sometimes segfaults in debug
    set_source_files_properties( generators/VTKMeshGeneratorTools.cpp PROPERTIES LANGUAGE CXX )
endif()

if( GEOS_ENABLE_TESTS )
    add_subdirectory( unitTests )
endif( )

geosx_add_code_checks( PREFIX mesh )
