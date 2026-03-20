/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SurfaceGenerator.cpp
 */

#include "SurfaceGenerator.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/kernels/SolidMechanicsLagrangianFEMKernels.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/surfaceGeneration/LogLevelsInfo.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGeneratorFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "kernels/surfaceGenerationKernels.hpp"

#include <algorithm>
#include <deque>
#include <stdio.h>
#include <sstream>
#include <unordered_map>

#if defined( GEOS_USE_CUDA )
#include <cuda_runtime.h>
#endif

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;

namespace internal
{

inline void printCudaStatus( char const * functionName, int const line, char const * label )
{
#if defined( GEOS_USE_CUDA )
  cudaError_t const launchErr = cudaGetLastError();
  cudaError_t const syncErr = cudaDeviceSynchronize();
  if( launchErr != cudaSuccess || syncErr != cudaSuccess )
  {
    printf( "\n[CUDA ERROR] %s:%d", functionName, line );
    if( label != nullptr && label[0] != '\0' )
    {
      printf( " %s", label );
    }
    printf( "\n" );

    if( launchErr != cudaSuccess )
    {
      printf( "  launch: %s\n", cudaGetErrorString( launchErr ) );
    }

    if( syncErr != cudaSuccess )
    {
      printf( "  sync: %s\n", cudaGetErrorString( syncErr ) );
    }
  }
#else
  GEOS_UNUSED_VAR( functionName );
  GEOS_UNUSED_VAR( line );
  GEOS_UNUSED_VAR( label );
#endif
}

}

#define GEOS_SURFACEGENERATOR_CUDA_STATUS( label ) \
  ::geos::internal::printCudaStatus( __func__, __LINE__, label )

/// Return the other edge of a face given one edge.
static localIndex GetOtherFaceEdge( const map< localIndex, std::pair< localIndex, localIndex > > & localFacesToEdges,
                                    const localIndex thisFace, const localIndex thisEdge )
{
  localIndex nextEdge = LOCALINDEX_MAX;

  const std::pair< localIndex, localIndex > & faceToEdges = stlMapLookup( localFacesToEdges, thisFace );
  if( faceToEdges.first == thisEdge )
  {
    nextEdge = faceToEdges.second;
  }
  else if( faceToEdges.second == thisEdge )
  {
    nextEdge = faceToEdges.first;
  }
  else
  {
    GEOS_ERROR( "SurfaceGenerator::Couldn't find thisEdge in localFacesToEdges[thisFace]" );
  }
  return nextEdge;
}

/**
 * @brief Removes dead-end faces from the rupture-ready working sets starting from @p edgeIndex.
 * @param[in]     edgeIndex                starting edge to check
 * @param[in]     isEdgeExternal           per-edge external flag (1 = mesh boundary)
 * @param[in]     edgeOnPartitionBoundary  edges that touch at least one ghost-boundary face
 * @param[in,out] edgesToRuptureReadyFaces working edge-to-ruptured-faces map
 * @param[in,out] localVFacesToVEdges      working face-to-(edge0,edge1) map
 * @param[in,out] nodeToRuptureReadyFaces  working set of ruptured faces for the node
 */
static void checkForAndRemoveDeadEndPath( const localIndex edgeIndex,
                                          arrayView1d< integer const > const & isEdgeExternal,
                                          std::set< localIndex > const & edgeOnPartitionBoundary,
                                          map< localIndex, std::set< localIndex > > & edgesToRuptureReadyFaces,
                                          map< localIndex, std::pair< localIndex, localIndex > > & localVFacesToVEdges,
                                          std::set< localIndex > & nodeToRuptureReadyFaces )
{
  GEOS_MARK_FUNCTION;

  localIndex thisEdge = edgeIndex;

  while( isEdgeExternal[thisEdge] != 1 &&
         edgeOnPartitionBoundary.count( thisEdge ) == 0 )
  {
    std::set< localIndex > & edgeToRuptureReadyFaces = edgesToRuptureReadyFaces[thisEdge];

    if( edgeToRuptureReadyFaces.size() != 1 )
      break;

    localIndex deadEndFace = *( edgeToRuptureReadyFaces.begin() );

    std::pair< localIndex, localIndex > & localVFaceToVEdges = stlMapLookup( localVFacesToVEdges, deadEndFace );

    localIndex nextEdge = -1;
    if( localVFaceToVEdges.first == thisEdge )
      nextEdge = localVFaceToVEdges.second;
    else if( localVFaceToVEdges.second == thisEdge )
      nextEdge = localVFaceToVEdges.first;
    else
    {
      GEOS_ERROR( "SurfaceGenerator::FindFracturePlanes: Could not find the next edge when removing dead end faces." );
    }

    edgeToRuptureReadyFaces.erase( deadEndFace );
    edgesToRuptureReadyFaces[nextEdge].erase( deadEndFace );
    nodeToRuptureReadyFaces.erase( deadEndFace );

    if( edgeToRuptureReadyFaces.empty() )
      edgesToRuptureReadyFaces.erase( thisEdge );
    if( edgesToRuptureReadyFaces[nextEdge].empty() )
      edgesToRuptureReadyFaces.erase( nextEdge );

    thisEdge = nextEdge;
  }
}


SurfaceGenerator::SurfaceGenerator( const string & name,
                                    Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_failCriterion( 1 ),
  m_nodeBasedSIF( 1 ),
  m_isPoroelastic( 0 ),
  m_initialRockToughness( 1.0e99 ),
  m_toughnessScalingFactor( 0.0 ),
  m_fractureOrigin( { 0.0, 0.0, 0.0 } ),
  m_mpiCommOrder( 0 )
{
  this->registerWrapper( viewKeyStruct::failCriterionString(), &this->m_failCriterion );

  registerWrapper( viewKeyStruct::initialRockToughnessString(), &m_initialRockToughness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Initial rock toughness of the solid material" );

  registerWrapper( viewKeyStruct::toughnessScalingFactorString(), &m_toughnessScalingFactor ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Scaling factor for the rock toughness of the solid material" );

  registerWrapper( viewKeyStruct::fractureOriginString(), &m_fractureOrigin ).
    setDefaultValue( m_fractureOrigin ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coordinate of fracture origin" );

  registerWrapper( viewKeyStruct::nodeBasedSIFString(), &m_nodeBasedSIF ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag for choosing between node or edge based criteria: 1 for node based criterion" );

  registerWrapper( viewKeyStruct::isPoroelasticString(), &m_isPoroelastic ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag that defines whether the material is poroelastic or not." );

  registerWrapper( viewKeyStruct::mpiCommOrderString(), &m_mpiCommOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to enable MPI consistent communication ordering" );

  registerWrapper( viewKeyStruct::fractureRegionNameString(), &m_fractureRegionName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( "Fracture" );

  registerWrapper( viewKeyStruct::tipNodesString(), &m_tipNodes ).
    setDescription( "Set containing all the nodes at the fracture tip" );

  registerWrapper( viewKeyStruct::tipEdgesString(), &m_tipEdges ).
    setDescription( "Set containing all the tip edges" );

  registerWrapper( viewKeyStruct::tipFacesString(), &m_tipFaces ).
    setDescription( "Set containing all the tip faces" );

  registerWrapper( viewKeyStruct::trailingFacesString(), &m_trailingFaces ).
    setDescription( "Set containing all the trailing faces" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  addLogLevel< logInfo::SurfaceGenerator >();
  addLogLevel< logInfo::Mapping >();
  addLogLevel< logInfo::RuptureRate >();
}

void SurfaceGenerator::postInputInitialization()
{
  static const std::set< integer > binaryOptions = { 0, 1 };

  GEOS_ERROR_IF( binaryOptions.count( m_isPoroelastic ) == 0,
                 "option can be either 0 (false) or 1 (true)",
                 getWrapperDataContext( viewKeyStruct::isPoroelasticString() ) );

  GEOS_ERROR_IF( binaryOptions.count( m_nodeBasedSIF ) == 0,
                 "option can be either 0 (false) or 1 (true)",
                 getWrapperDataContext( viewKeyStruct::nodeBasedSIFString() ) );

  GEOS_ERROR_IF( binaryOptions.count( m_mpiCommOrder ) == 0,
                 "option can be either 0 (false) or 1 (true)",
                 getWrapperDataContext( viewKeyStruct::mpiCommOrderString() ) );
}

SurfaceGenerator::~SurfaceGenerator()
{}

void SurfaceGenerator::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion & subRegion )
    {
      subRegion.registerField< surfaceGeneration::K_IC_00,
                               surfaceGeneration::K_IC_01,
                               surfaceGeneration::K_IC_02,
                               surfaceGeneration::K_IC_10,
                               surfaceGeneration::K_IC_11,
                               surfaceGeneration::K_IC_12,
                               surfaceGeneration::K_IC_20,
                               surfaceGeneration::K_IC_21,
                               surfaceGeneration::K_IC_22 >( this->getName() );
    } );

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      subRegion.registerField< surfaceGeneration::K_IC_00,
                               surfaceGeneration::K_IC_01,
                               surfaceGeneration::K_IC_02,
                               surfaceGeneration::K_IC_10,
                               surfaceGeneration::K_IC_11,
                               surfaceGeneration::K_IC_12,
                               surfaceGeneration::K_IC_20,
                               surfaceGeneration::K_IC_21,
                               surfaceGeneration::K_IC_22,
                               fields::ruptureTime,
                               surfaceGeneration::ruptureRate >( this->getName() );
    } );

    NodeManager & nodeManager = mesh.getNodeManager();
    EdgeManager & edgeManager = mesh.getEdgeManager();
    FaceManager & faceManager = mesh.getFaceManager();

    nodeManager.registerField< fields::parentIndex,
                               fields::childIndex,
                               surfaceGeneration::degreeFromCrack,
                               surfaceGeneration::degreeFromCrackTip,
                               surfaceGeneration::SIFNode,
                               fields::ruptureTime >( this->getName() );

    edgeManager.registerField< fields::parentIndex,
                               fields::childIndex,
                               surfaceGeneration::SIF_I,
                               surfaceGeneration::SIF_II,
                               surfaceGeneration::SIF_III >( this->getName() );

    faceManager.registerField< fields::parentIndex,
                               fields::childIndex,
                               surfaceGeneration::ruptureState,
                               fields::ruptureTime,
                               surfaceGeneration::SIFonFace,
                               surfaceGeneration::K_IC,
                               surfaceGeneration::primaryCandidateFace,
                               surfaceGeneration::isFaceSeparable,
                               surfaceGeneration::degreeFromCrackTip >( this->getName() );

    // TODO: handle this in registerField().
    faceManager.getField< surfaceGeneration::K_IC >().resizeDimension< 1 >( 3 );

    Group & problemManager = this->getGroupByPath( "/Problem" );
    FieldSpecificationManager & fsm =  problemManager.getGroup< FieldSpecificationManager >( "FieldSpecifications" );
    fsm.setIsSurfaceGenerationCase( true );
  } );


}

void SurfaceGenerator::initializePostInitialConditionsPreSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                string_array const & )
  {
    NodeManager & nodeManager = meshLevel.getNodeManager();
    FaceManager & faceManager = meshLevel.getFaceManager();

    arrayView1d< localIndex > const & parentNodeIndex = nodeManager.getField< fields::parentIndex >();

    arrayView1d< localIndex > const & parentFaceIndex = faceManager.getField< fields::parentIndex >();

    arrayView1d< localIndex > const & childFaceIndex = faceManager.getField< fields::childIndex >();

    parentNodeIndex.setValues< serialPolicy >( -1 );
    parentFaceIndex.setValues< serialPolicy >( -1 );
    childFaceIndex.setValues< serialPolicy >( -1 );

    m_originalNodetoFaces = nodeManager.faceList();
    m_originalNodetoEdges = nodeManager.edgeList();
    m_originalFaceToEdges = faceManager.edgeList();

    string const usedFacesLabel = "usedFaces";
    nodeManager.registerWrapper( usedFacesLabel, &m_usedFacesForNode );
    nodeManager.excludeWrappersFromPacking( { usedFacesLabel } );
    m_usedFacesForNode.resize( nodeManager.size() );

    localIndex const numFaces = faceManager.size();
    m_originalFacesToElemRegion.resize( numFaces, 2 );
    m_originalFacesToElemSubRegion.resize( numFaces, 2 );
    m_originalFacesToElemIndex.resize( numFaces, 2 );

    for( localIndex faceID = 0; faceID < numFaces; ++faceID )
    {
      for( localIndex side = 0; side < 2; ++side )
      {
        m_originalFacesToElemRegion( faceID, side ) = faceManager.elementRegionList()( faceID, side );
        m_originalFacesToElemSubRegion( faceID, side ) = faceManager.elementSubRegionList()( faceID, side );
        m_originalFacesToElemIndex( faceID, side ) = faceManager.elementList()( faceID, side );
      }
    }
  } );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                string_array const & )
  {
    FaceManager & faceManager = meshLevel.getFaceManager();
    ElementRegionManager & elementManager = meshLevel.getElemManager();
    arrayView2d< real64 const > const & faceNormals = faceManager.faceNormal();
    arrayView2d< real64 const > const & faceCenters = faceManager.faceCenter();

    // TODO: move roughness-to-KIC to a constitutive relation.
    arrayView2d< real64 > const & KIC = faceManager.getField< surfaceGeneration::K_IC >();

    for( localIndex kf=0; kf<faceManager.size(); ++kf )
    {
      if( m_initialRockToughness >= 0 )
      {
        KIC[kf][0] = m_initialRockToughness;
        KIC[kf][1] = m_initialRockToughness;
        KIC[kf][2] = m_initialRockToughness;
      }
      else
      {
        arrayView2d< localIndex const > const & faceToRegionMap = faceManager.elementRegionList();
        arrayView2d< localIndex const > const & faceToSubRegionMap = faceManager.elementSubRegionList();
        arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();

        KIC[kf][0] = 1e99;
        KIC[kf][1] = 1e99;
        KIC[kf][2] = 1e99;

        for( localIndex k=0; k<faceToRegionMap.size( 1 ); ++k )
        {
          localIndex const er = faceToRegionMap[kf][k];
          localIndex const esr = faceToSubRegionMap[kf][k];
          localIndex const ei = faceToElementMap[kf][k];

          if( er != -1 &&  esr != -1 && ei != -1 )
          {
            CellElementSubRegion & elementSubRegion = elementManager.getRegion( faceToRegionMap[kf][k] ).
                                                        getSubRegion< CellElementSubRegion >( faceToSubRegionMap[kf][k] );
            localIndex iEle = faceToElementMap[kf][k];

            arrayView1d< real64 const > const K_IC_00 = elementSubRegion.getField< surfaceGeneration::K_IC_00 >();
            arrayView1d< real64 const > const K_IC_01 = elementSubRegion.getField< surfaceGeneration::K_IC_01 >();
            arrayView1d< real64 const > const K_IC_02 = elementSubRegion.getField< surfaceGeneration::K_IC_02 >();
            arrayView1d< real64 const > const K_IC_10 = elementSubRegion.getField< surfaceGeneration::K_IC_10 >();
            arrayView1d< real64 const > const K_IC_11 = elementSubRegion.getField< surfaceGeneration::K_IC_11 >();
            arrayView1d< real64 const > const K_IC_12 = elementSubRegion.getField< surfaceGeneration::K_IC_12 >();
            arrayView1d< real64 const > const K_IC_20 = elementSubRegion.getField< surfaceGeneration::K_IC_20 >();
            arrayView1d< real64 const > const K_IC_21 = elementSubRegion.getField< surfaceGeneration::K_IC_21 >();
            arrayView1d< real64 const > const K_IC_22 = elementSubRegion.getField< surfaceGeneration::K_IC_22 >();

            real64 k0[3];
            k0[0] = K_IC_00[iEle]*faceNormals[kf][0] + K_IC_10[iEle]*faceNormals[kf][1] + K_IC_20[iEle]*faceNormals[kf][2];
            k0[1] = K_IC_01[iEle]*faceNormals[kf][0] + K_IC_11[iEle]*faceNormals[kf][1] + K_IC_21[iEle]*faceNormals[kf][2];
            k0[2] = K_IC_02[iEle]*faceNormals[kf][0] + K_IC_12[iEle]*faceNormals[kf][1] + K_IC_22[iEle]*faceNormals[kf][2];

            KIC[kf][0] = std::min( std::fabs( k0[0] ), std::fabs( KIC[kf][0] ) );
            KIC[kf][1] = std::min( std::fabs( k0[1] ), std::fabs( KIC[kf][1] ) );
            KIC[kf][2] = std::min( std::fabs( k0[2] ), std::fabs( KIC[kf][2] ) );
          }
        }
      }
      if( m_toughnessScalingFactor > 0.0 )
      {
        real64 faceCenter[3];
        faceCenter[0] = faceCenters[kf][0];
        faceCenter[1] = faceCenters[kf][1];
        faceCenter[2] = faceCenters[kf][2];

        for( localIndex dim=0; dim<3; ++dim )
        {
          real64 const initialRockToughness = KIC[kf][dim];

          real64 const scaledToughness = scalingToughness( m_fractureOrigin,
                                                           faceCenter,
                                                           initialRockToughness,
                                                           m_toughnessScalingFactor );

          KIC[kf][dim] = scaledToughness;
        }
      }
    }
  } );
}

void SurfaceGenerator::postRestartInitialization()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  // Repopulate the fracture stencil.
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                string_array const & )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();
    SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );
    FaceElementSubRegion & fractureSubRegion = fractureRegion.getSubRegion< FaceElementSubRegion >( 0 );

    for( localIndex fce = 0; fce < fractureSubRegion.m_2dFaceTo2dElems.size(); ++fce )
    {
      fractureSubRegion.m_recalculateConnectionsFor2dFaces.insert( fce );
    }

    for( localIndex fe = 0; fe < fractureSubRegion.size(); ++fe )
    {
      fractureSubRegion.m_newFaceElements.insert( fe );
    }

    for( localIndex a = 0; a < fvManager.numSubGroups(); ++a )
    {
      FluxApproximationBase * const fluxApprox = fvManager.getGroupPointer< FluxApproximationBase >( a );
      if( fluxApprox!=nullptr )
      {
        fluxApprox->addToFractureStencil( meshLevel, this->m_fractureRegionName );
      }
    }

    fractureSubRegion.m_recalculateConnectionsFor2dFaces.clear();
    fractureSubRegion.m_newFaceElements.clear();
  } );
}


real64 SurfaceGenerator::solverStep( real64 const & time_n,
                                     real64 const & dt,
                                     const int GEOS_UNUSED_PARAM( cycleNumber ),
                                     DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  int rval = 0;
  std::set< string > modifiedMeshBodies;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel & meshLevel,
                                                                string_array const & )
  {
    SpatialPartition & partition = dynamicCast< SpatialPartition & >( domain.getReference< PartitionBase >( dataRepository::keys::partitionManager ) );
    int const tileColor=partition.getColor();
    int const numTileColorsLocal=partition.numColor();
    int const numTileColors = MpiWrapper::allReduce( numTileColorsLocal, MpiWrapper::Reduction::Max );

    int const localRval = separationDriver( domain,
                                            meshLevel,
                                            domain.getNeighbors(),
                                            tileColor,
                                            numTileColors,
                                            0,
                                            time_n + dt );
    GEOS_SURFACEGENERATOR_CUDA_STATUS( "after separationDriver" );

    rval += localRval;

    // Mark the mesh body as modified on all ranks if any rank split a node.
    if( MpiWrapper::allReduce( localRval, MpiWrapper::Reduction::Max ) > 0 )
    {
      modifiedMeshBodies.insert( meshBodyName );
    }

  } );

  for( auto const & bodyName : modifiedMeshBodies )
  {
    MeshBody & meshBody = domain.getMeshBodies().getGroup< MeshBody >( bodyName );
    dataRepository::Group & levels = meshBody.getMeshLevels();
    for( auto const & levelName : levels.getSubGroupsNames() )
    {
      if( levels.hasGroup< MeshLevel >( levelName ) )
      {
        MeshLevel & level = levels.getGroup< MeshLevel >( levelName );
        level.modified();
      }
    }
  }

  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                string_array const & )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();
    SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );

    for( localIndex a=0; a<fvManager.numSubGroups(); ++a )
    {
      FluxApproximationBase * const fluxApprox = fvManager.getGroupPointer< FluxApproximationBase >( a );
      if( fluxApprox!=nullptr )
      {
        fluxApprox->addToFractureStencil( meshLevel, this->m_fractureRegionName );
      }
    }

    FaceElementSubRegion & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

    // Rebuild geometric sets and populate the "all" set on the fracture subregion.
    meshLevel.getNodeManager().buildGeometricSets( GeometricObjectManager::getInstance() );
    dataRepository::Group & setGroup =
      fractureSubRegion.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );

    SortedArray< localIndex > & targetSet =
      setGroup.getWrapper< SortedArray< localIndex > >( "all" ).reference();

    forAll< serialPolicy >( fractureSubRegion.size(), [&] ( localIndex const ei )

    {
      targetSet.insert( ei );
    } );

    // Compute gravity coefficient for new fracture elements.
    real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    if( fractureSubRegion.hasField< fields::flow::gravityCoefficient >() )
    {
      arrayView2d< real64 const > const elemCenter = fractureSubRegion.getElementCenter();

      arrayView1d< real64 > const gravityCoef = fractureSubRegion.getField< fields::flow::gravityCoefficient >();

      forAll< parallelHostPolicy >( fractureSubRegion.size(), [=] ( localIndex const ei )
      {
        gravityCoef[ ei ] = LvArray::tensorOps::AiBi< 3 >( elemCenter[ ei ], gravVector );
      } );
    }

    string const permModelName = getConstitutiveName< PermeabilityBase >( fractureSubRegion );
    if( !permModelName.empty() )
    {
      // Initialize permeability state for newly created fracture elements.
      PermeabilityBase & permModel = getConstitutiveModel< PermeabilityBase >( fractureSubRegion, permModelName );
      permModel.initializeState();
    }
  } );

  GEOS_SURFACEGENERATOR_CUDA_STATUS( "end solverStep" );
  return rval;
}

bool SurfaceGenerator::execute( real64 const time_n,
                                real64 const dt,
                                integer const cycleNumber,
                                integer const GEOS_UNUSED_PARAM( eventCounter ),
                                real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                DomainPartition & domain )
{
  solverStep( time_n, dt, cycleNumber, domain );
  GEOS_SURFACEGENERATOR_CUDA_STATUS( "end execute" );
  return false;
}

int SurfaceGenerator::separationDriver( DomainPartition & domain,
                                        MeshLevel & mesh,
                                        stdVector< NeighborCommunicator > & neighbors,
                                        int const tileColor,
                                        int const numTileColors,
                                        bool const prefrac,
                                        real64 const time_np1 )
{
  GEOS_MARK_FUNCTION;

  m_faceElemsRupturedThisSolve.clear();
  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elementManager = mesh.getElemManager();

  // Snapshot locally-owned fracture element count before splitting.
  {
    SurfaceElementRegion const & fractureRegion = elementManager.getRegion< SurfaceElementRegion >( m_fractureRegionName );
    FaceElementSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();
    m_numFractureElementsBefore = fractureSubRegion.getNumberOfLocalIndices();
  }

  stdVector< std::set< localIndex > > nodesToRupturedFaces;
  stdVector< std::set< localIndex > > edgesToRupturedFaces;

  ArrayOfArrays< localIndex > const & nodeToElementMap = nodeManager.elementList();

  FieldIdentifiers fieldsToBeSync;

  fieldsToBeSync.addFields( FieldLocation::Face, { surfaceGeneration::ruptureState::key() } );
  if( nodeManager.hasField< fields::solidMechanics::externalForce >() )
  {
    fieldsToBeSync.addFields( FieldLocation::Node, { fields::solidMechanics::externalForce::key() } );
  }

  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

  elementManager.forElementSubRegions< CellElementSubRegion >( [] ( auto & elemSubRegion )
  {
    elemSubRegion.moveSets( hostMemorySpace );
  } );
  faceManager.moveSets( hostMemorySpace );
  edgeManager.moveSets( hostMemorySpace );
  nodeManager.moveSets( hostMemorySpace );

  if( !prefrac )
  {
    if( m_failCriterion >0 )  // Stress intensity factor based criterion and mixed criterion.
    {

      identifyRupturedFaces( domain,
                             nodeManager,
                             edgeManager,
                             faceManager,
                             elementManager,
                             prefrac );

    }
  }


  if( prefrac )
  {
    ModifiedObjectLists modifiedObjects;
    calculateKinkAngles( faceManager, edgeManager, nodeManager, modifiedObjects, prefrac );
  }

  // Populate nodesToRupturedFaces and edgesToRupturedFaces from current rupture state.
  postUpdateRuptureStates( nodeManager,
                           edgeManager,
                           faceManager,
                           elementManager,
                           nodesToRupturedFaces,
                           edgesToRupturedFaces );

  int rval = 0;

  array1d< integer > const & isNodeGhost = nodeManager.ghostRank();

  for( int color=0; color<numTileColors; ++color )
  {
    ModifiedObjectLists modifiedObjects;
    if( color==tileColor )
    {
      for( localIndex a=0; a<nodeManager.size(); ++a )
      {
        int didSplit = 0;
        if( isNodeGhost[a]<0 &&
            nodeToElementMap.sizeOfArray( a )>1 )
        {
          didSplit += processNode( a,
                                   time_np1,
                                   nodeManager,
                                   edgeManager,
                                   faceManager,
                                   elementManager,
                                   nodesToRupturedFaces,
                                   edgesToRupturedFaces,
                                   elementManager,
                                   modifiedObjects, prefrac );
          if( didSplit > 0 )
          {
            rval += didSplit;
            --a;
          }
        }
      }

      SurfaceElementRegion & fractureElementRegion = elementManager.getRegion< SurfaceElementRegion >( m_fractureRegionName );
      fractureElementRegion.updateSets( faceManager );
    }

#ifdef GEOS_USE_MPI

    modifiedObjects.clearNewFromModified();

    // 1) Assign new global indices to the new objects
    CommunicationTools::assignNewGlobalIndices( nodeManager, modifiedObjects.newNodes );
    CommunicationTools::assignNewGlobalIndices( edgeManager, modifiedObjects.newEdges );
    CommunicationTools::assignNewGlobalIndices( faceManager, modifiedObjects.newFaces );

    ModifiedObjectLists receivedObjects;

    GEOS_SURFACEGENERATOR_CUDA_STATUS( "separationDriver.beforeSynchronizeTopologyChange" );
    parallelTopologyChange::synchronizeTopologyChange( &mesh,
                                                       neighbors,
                                                       modifiedObjects,
                                                       receivedObjects,
                                                       m_mpiCommOrder );
    GEOS_SURFACEGENERATOR_CUDA_STATUS( "separationDriver.afterSynchronizeTopologyChange" );

    synchronizeTipSets( faceManager,
                        edgeManager,
                        nodeManager,
                        receivedObjects );


#else

    GEOS_UNUSED_VAR( neighbors );
    assignNewGlobalIndicesSerial( nodeManager, modifiedObjects.newNodes );
    assignNewGlobalIndicesSerial( edgeManager, modifiedObjects.newEdges );
    assignNewGlobalIndicesSerial( faceManager, modifiedObjects.newFaces );

#endif

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    elementManager.forElementSubRegionsComplete< FaceElementSubRegion >( [&]( localIndex const er,
                                                                              localIndex const esr,
                                                                              ElementRegionBase &,
                                                                              FaceElementSubRegion & subRegion )
    {
      std::set< localIndex > & newFaceElems = modifiedObjects.newElements[{er, esr}];
      for( localIndex const newFaceElemIndex : newFaceElems )
      {
        subRegion.m_newFaceElements.insert( newFaceElemIndex );
      }
    } );


    elementManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      FaceElementSubRegion::NodeMapType & nodeMap = subRegion.nodeList();
      arrayView2d< localIndex const > const faceMap = subRegion.faceList().toViewConst();

      for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
      {

        localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( faceMap[ kfe ][ 0 ] );
        nodeMap.resizeArray( kfe, 2*numNodesInFace );
        for( localIndex a = 0; a < numNodesInFace; ++a )
        {

          nodeMap[ kfe ][ a ]   = faceToNodeMap( faceMap[ kfe ][ 0 ], a );
          nodeMap[ kfe ][ a + numNodesInFace ] = faceToNodeMap( faceMap[ kfe ][ 1 ], a );
        }

      }
    } );
    MpiWrapper::barrier();
  } // end color loop


  real64 ruptureRate = calculateRuptureRate( elementManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName ) );

  GEOS_LOG_LEVEL_RANK_0( logInfo::RuptureRate, GEOS_FMT( "Rupture rate = {}", ruptureRate ) );
  if( ruptureRate > 0 )
    m_nextDt = ruptureRate < 1e99 ? m_cflFactor / ruptureRate : 1e99;

  {
    elementManager.forElementSubRegions< CellElementSubRegion >( [] ( auto & elemSubRegion )
    {
      elemSubRegion.nodeList().registerTouch( hostMemorySpace );
      elemSubRegion.edgeList().registerTouch( hostMemorySpace );
      elemSubRegion.faceList().registerTouch( hostMemorySpace );
    } );


    faceManager.nodeList().toView().registerTouch( hostMemorySpace );
    faceManager.elementList().registerTouch( hostMemorySpace );
    faceManager.elementRegionList().registerTouch( hostMemorySpace );
    faceManager.elementSubRegionList().registerTouch( hostMemorySpace );

    edgeManager.nodeList().registerTouch( hostMemorySpace );

  }

  // Log statistics about the mesh splitting operation
  int const globalRval = MpiWrapper::allReduce( rval, MpiWrapper::Reduction::Max, MPI_COMM_GEOS );
  if( globalRval > 0 )
  {
    // Count fracture elements from the FaceElementSubRegion.
    SurfaceElementRegion const & fractureRegion = elementManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );
    FaceElementSubRegion const & fractureSubRegion = fractureRegion.getUniqueSubRegion< FaceElementSubRegion >();

    // Count locally-owned fracture elements only (exclude ghosts).
    localIndex const localNumFractureElements = fractureSubRegion.getNumberOfLocalIndices();
    localIndex const localNewFractureElements = localNumFractureElements - m_numFractureElementsBefore;

    // Gather global statistics across all MPI ranks
    globalIndex const globalNumFractureElements = MpiWrapper::sum( static_cast< globalIndex >( localNumFractureElements ), MPI_COMM_GEOS );
    globalIndex const globalNewFractureElements = MpiWrapper::sum( static_cast< globalIndex >( localNewFractureElements ), MPI_COMM_GEOS );
    globalIndex const globalNumSplits           = MpiWrapper::sum( static_cast< globalIndex >( rval ), MPI_COMM_GEOS );

    GEOS_LOG_RANK_0( GEOS_FMT( "SurfaceGenerator: Mesh splitting completed.\n"
                               "  Number of nodes split:                      {:>8}\n"
                               "  New fracture elements:                      {:>8}\n"
                               "  Total fracture elements:                    {:>8}",
                               globalNumSplits,
                               globalNewFractureElements,
                               globalNumFractureElements ) );

    // Extended per-rank breakdown, guarded by log level
    GEOS_LOG_LEVEL_RANK_0( logInfo::SurfaceGenerator,
                           GEOS_FMT( "SurfaceGenerator: Per-rank breakdown.\n"
                                     "  Local fracture elements (this rank):        {:>8}\n"
                                     "  Local new fracture elements (this rank):    {:>8}\n"
                                     "  Local nodes split (this rank):              {:>8}",
                                     localNumFractureElements,
                                     localNewFractureElements,
                                     rval ) );
  }

  GEOS_SURFACEGENERATOR_CUDA_STATUS( "end separationDriver" );
  return rval;
}

void SurfaceGenerator::synchronizeTipSets ( FaceManager & faceManager,
                                            EdgeManager & edgeManager,
                                            NodeManager & nodeManager,
                                            ModifiedObjectLists & receivedObjects )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< localIndex const > const & parentNodeIndices = nodeManager.getField< fields::parentIndex >();

  for( localIndex const nodeIndex : receivedObjects.newNodes )
  {
    localIndex const parentNodeIndex = parentNodeIndices[nodeIndex];

    GEOS_ERROR_IF( parentNodeIndex == -1, "parentNodeIndex should not be -1",
                   getDataContext() );

    m_tipNodes.remove( parentNodeIndex );
  }

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();
  arrayView1d< integer > const & edgeIsExternal = edgeManager.isExternal();
  arrayView1d< integer > const & nodeIsExternal = nodeManager.isExternal();


  arrayView1d< localIndex const > const &
  parentEdgeIndices = edgeManager.getField< fields::parentIndex >();

  arrayView1d< localIndex const > const &
  childEdgeIndices = edgeManager.getField< fields::childIndex >();


  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();

  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

  for( localIndex const edgeIndex : receivedObjects.newEdges )
  {
    localIndex const parentEdgeIndex = parentEdgeIndices[edgeIndex];

    GEOS_ERROR_IF( parentEdgeIndex == -1, "parentEdgeIndex should not be -1",
                   getDataContext() );

    m_tipEdges.remove( parentEdgeIndex );
    for( localIndex const faceIndex : edgeToFaceMap[ parentEdgeIndex ] )
    {
      bool trailingFace = false;
      if( m_trailingFaces.contains( faceIndex ))
      {
        for( localIndex const faceLocalEdgeIndex : faceToEdgeMap[ faceIndex ] )
        {
          if( m_tipEdges.contains( faceLocalEdgeIndex ))
          {
            trailingFace = true;
          }
        }

        if( trailingFace == false )
        {
          m_trailingFaces.remove( faceIndex );
        }
      }
    }
  }

  arrayView1d< integer const > const & isFaceSeparable = faceManager.getField< surfaceGeneration::isFaceSeparable >();
  arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();

  arrayView1d< localIndex const > const & childNodeIndices = nodeManager.getField< fields::childIndex >();
  arrayView1d< localIndex > const & parentFaceIndices = faceManager.getField< fields::parentIndex >();

  for( localIndex const faceIndex : receivedObjects.newFaces )
  {
    localIndex const parentFaceIndex = parentFaceIndices[faceIndex];
    GEOS_ERROR_IF( parentFaceIndex == -1, "parentFaceIndex should not be -1",
                   getDataContext() );

    m_trailingFaces.insert( parentFaceIndex );
    m_tipFaces.remove( parentFaceIndex );

    for( localIndex const edgeIndex : faceManager.edgeList()[ parentFaceIndex ] )
    {
      if( parentEdgeIndices[edgeIndex]==-1 && childEdgeIndices[edgeIndex]==-1 )
      {
        m_tipEdges.insert( edgeIndex );

        for( localIndex const iface: edgeManager.faceList()[ edgeIndex ] )
        {
          if( faceToElementMap.size( 1 ) == 2  &&
              faceIsExternal[iface] < 1 &&
              isFaceSeparable[iface] == 1 )
          {
            m_tipFaces.insert( iface );
          }
        }
      }
      if( edgeIsExternal[edgeIndex]==0 )
      {
        edgeIsExternal[edgeIndex] = 2;
      }
    }
    for( localIndex const nodeIndex : faceManager.nodeList()[ parentFaceIndex ] )
    {
      if( parentNodeIndices[nodeIndex]==-1 && childNodeIndices[nodeIndex]==-1 )
      {
        m_tipNodes.insert( nodeIndex );
      }
      if( nodeIsExternal[nodeIndex] )
      {
        nodeIsExternal[nodeIndex] = 2;
      }
    }
  }
}

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
bool SurfaceGenerator::processNode( const localIndex nodeID,
                                    real64 const time_np1,
                                    NodeManager & nodeManager,
                                    EdgeManager & edgeManager,
                                    FaceManager & faceManager,
                                    ElementRegionManager & elemManager,
                                    stdVector< std::set< localIndex > > & nodesToRupturedFaces,
                                    stdVector< std::set< localIndex > > & edgesToRupturedFaces,
                                    ElementRegionManager & elementManager,
                                    ModifiedObjectLists & modifiedObjects,
                                    const bool GEOS_UNUSED_PARAM( prefrac ) )
{
  GEOS_MARK_FUNCTION;

  bool didSplit = false;
  bool fracturePlaneFlag = true;

  {
    std::set< localIndex > facialRupturePath;
    map< localIndex, int > edgeLocations;
    map< localIndex, int > faceLocations;
    map< std::pair< CellElementSubRegion const *, localIndex >, int > elemLocations;

    fracturePlaneFlag = findFracturePlanes( nodeID,
                                            nodeManager,
                                            edgeManager,
                                            faceManager,
                                            elemManager,
                                            nodesToRupturedFaces,
                                            edgesToRupturedFaces,
                                            facialRupturePath,
                                            edgeLocations,
                                            faceLocations,
                                            elemLocations );
    if( fracturePlaneFlag )
    {
      mapConsistencyCheck( nodeID, nodeManager, edgeManager, faceManager, elementManager, elemLocations );

      didSplit = true;
      performFracture( nodeID,
                       time_np1,
                       nodeManager,
                       edgeManager,
                       faceManager,
                       elementManager,
                       modifiedObjects,
                       nodesToRupturedFaces,
                       edgesToRupturedFaces,
                       facialRupturePath,
                       edgeLocations,
                       faceLocations,
                       elemLocations );
      mapConsistencyCheck( nodeID, nodeManager, edgeManager, faceManager, elementManager, elemLocations );

    }
  }

  return didSplit;
}

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
bool SurfaceGenerator::findFracturePlanes( localIndex const nodeID,
                                           NodeManager const & nodeManager,
                                           EdgeManager const & edgeManager,
                                           FaceManager const & faceManager,
                                           ElementRegionManager const & elemManager,
                                           stdVector< std::set< localIndex > > const & nodesToRupturedFaces,
                                           stdVector< std::set< localIndex > > const & edgesToRupturedFaces,
                                           std::set< localIndex > & separationPathFaces,
                                           map< localIndex, int > & edgeLocations,
                                           map< localIndex, int > & faceLocations,
                                           map< std::pair< CellElementSubRegion const *, localIndex >, int > & elemLocations )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< localIndex const > const & parentNodeIndices = nodeManager.getField< fields::parentIndex >();

  localIndex const parentNodeIndex = ObjectManagerBase::getParentRecursive( parentNodeIndices, nodeID );

  arrayView1d< localIndex const > const & parentFaceIndices = faceManager.getField< fields::parentIndex >();
  arrayView1d< localIndex const > const & childFaceIndices = faceManager.getField< fields::childIndex >();

  std::set< localIndex > const & vNodeToRupturedFaces = nodesToRupturedFaces[parentNodeIndex];

  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = nodeManager.edgeList().toViewConst();
  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = nodeManager.faceList().toViewConst();

  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

  arraySlice1d< localIndex const > const nodeToRegionMap = nodeManager.elementRegionList()[nodeID];
  arraySlice1d< localIndex const > const nodeToSubRegionMap = nodeManager.elementSubRegionList()[nodeID];
  arraySlice1d< localIndex const > const nodeToElementMap = nodeManager.elementList()[nodeID];

  // Build deterministic, sorted list of unique (subRegion, elementIndex) pairs for node-attached elements.
  auto buildNodeToElementMaps = [&]()
  {
    stdVector< std::pair< CellElementSubRegion const *, localIndex > > result;
    for( localIndex k = 0; k < nodeManager.elementRegionList().sizeOfArray( nodeID ); ++k )
    {
      localIndex const er  = nodeToRegionMap[k];
      localIndex const esr = nodeToSubRegionMap[k];
      localIndex const ei  = nodeToElementMap[k];
      CellElementSubRegion const * subRegion =
        &elemManager.getRegion( er ).getSubRegion< CellElementSubRegion >( esr );
      std::pair< CellElementSubRegion const *, localIndex > const p( subRegion, ei );

      // Keep pairs unique (small N — linear search is fine).
      if( std::find( result.cbegin(), result.cend(), p ) == result.cend() )
      {
        result.push_back( p );
      }
    }

    // Sort by global element index for platform independence.
    std::sort( result.begin(), result.end(),
               []( auto const & a, auto const & b )
    {
      return a.first->localToGlobalMap()[a.second]
      < b.first->localToGlobalMap()[b.second];
    } );

    return result;
  };

  stdVector< std::pair< CellElementSubRegion const *, localIndex > > const
  nodeToElementMaps( buildNodeToElementMaps() );

  arrayView1d< integer const > const & isEdgeExternal = edgeManager.isExternal();

//  const std::set<localIndex>& usedFaces = nodeManager.GetUnorderedVariableOneToManyMap("usedFaces")[nodeID];

  // **** local working arrays *****************************************************************************************

  // Identify rupture-ready faces from intersection of node-attached faces and ruptured parent faces.
  std::set< localIndex > nodeToRuptureReadyFaces;
  for( localIndex const i : nodeToFaceMap[ nodeID ] )
  {
    const localIndex parentFaceIndex = ( parentFaceIndices[i] == -1 ) ? i : parentFaceIndices[i];

    if( vNodeToRupturedFaces.count( parentFaceIndex ) > 0 )
    {
      nodeToRuptureReadyFaces.insert( parentFaceIndex );
    }
  }



  // local map to hold the edgesToRuptureReadyFaces
  map< localIndex, std::set< localIndex > > edgesToRuptureReadyFaces;
  for( localIndex const edgeIndex : m_originalNodetoEdges[ parentNodeIndex ] )
  {
    if( !(edgesToRupturedFaces[edgeIndex].empty()) )
      edgesToRuptureReadyFaces[edgeIndex].insert( edgesToRupturedFaces[edgeIndex].begin(), edgesToRupturedFaces[edgeIndex].end() );
  }


  // Map node-attached faces to their incident edges.
  map< localIndex, std::pair< localIndex, localIndex > > nodeLocalFacesToEdges;
  for( localIndex const kf : m_originalNodetoFaces[ parentNodeIndex ] )
  {
    localIndex edge[2] = { INT_MAX, INT_MAX };
    int count = 0;
    for( localIndex const ke : m_originalFaceToEdges[ kf ] )
    {
      if( m_originalNodetoEdges.contains( parentNodeIndex, ke ) )
      {
        edge[count++] = ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      GEOS_ERROR( "invalid edge (SurfaceGenerator::findFracturePlanes).", getDataContext() );
    }


    nodeLocalFacesToEdges[kf] = std::make_pair( edge[0], edge[1] );

  }


  // Remove dead-end paths, preserving MPI partition boundaries where neighbor data is incomplete.

  // Identify edges touching MPI partition boundaries.
  std::set< localIndex > edgeOnPartitionBoundary;
  for( localIndex const kf : m_originalNodetoFaces[ parentNodeIndex ] )
  {
    if( isGhostBoundaryFace( kf ) )
    {
      // Both edges of this face touch the partition boundary (or domain exterior).
      auto const & ep = nodeLocalFacesToEdges[kf];
      edgeOnPartitionBoundary.insert( ep.first );
      edgeOnPartitionBoundary.insert( ep.second );
    }
  }

  // loop over all the edges
  for( localIndex const edgeIndex : m_originalNodetoEdges[ parentNodeIndex ] )
  {
    checkForAndRemoveDeadEndPath( edgeIndex,
                                  isEdgeExternal,
                                  edgeOnPartitionBoundary,
                                  edgesToRuptureReadyFaces,
                                  nodeLocalFacesToEdges,
                                  nodeToRuptureReadyFaces );
  }

  // if there are no ruptured faces attached to the node, then we are done.
  if( nodeToRuptureReadyFaces.empty() )
  {
    return false;
  }

  // ***** find separation path ****************************************************************************************

  // ***** find starting face *****
  // Find unused starting face for separation path.
  localIndex startingEdge = INT_MAX;
  localIndex startingFace = INT_MAX;
  bool startingEdgeExternal = false;

  for( std::set< localIndex >::const_iterator i=nodeToRuptureReadyFaces.begin(); i!=nodeToRuptureReadyFaces.end(); ++i )
  {
    // Check if face was already used in a prior split.
    if( m_usedFacesForNode[nodeID].count( *i )==0 )
    {
      // great! It hasn't. It's on like Donkey Kong.
      startingFace = *i;

      if( isEdgeExternal[nodeLocalFacesToEdges[startingFace].first]==1 )
      {
        startingEdge = nodeLocalFacesToEdges[startingFace].first;
        startingEdgeExternal = true;
        break;
      }
      else if( isEdgeExternal[nodeLocalFacesToEdges[startingFace].second]==1 )
      {
        startingEdge = nodeLocalFacesToEdges[startingFace].second;
        startingEdgeExternal = true;
        break;
      }
      else
      {
        startingEdge = nodeLocalFacesToEdges[startingFace].first;
      }
    }
  }

  // Abort if no starting face found.
  if( startingFace==INT_MAX || startingEdge==INT_MAX )
  {
    return false;
    //    GEOS_ERROR("Fracturantor3::FindFracturePlanes: couldn't set starting face/edge");
  }



  // Process first found separation path; subsequent paths will be handled in later iterations.



  // Define separation plane to split objects into two sets.



  {
    // Begin setting separation path.
    localIndex thisEdge = startingEdge;
    localIndex thisFace = startingFace;

    localIndex nextEdge = INT_MAX;
    localIndex nextFace = INT_MAX;

    //localIndex lastEdge = INT_MAX;
    //localIndex lastFace = INT_MAX;

    // the seprationPath is used to hold combinations of edge and face
    map< localIndex, int > facesInPath;
    map< localIndex, int > edgesInPath;

    int numFacesInPath = 0;
    edgesInPath[thisEdge] = numFacesInPath;
    facesInPath[thisFace] = numFacesInPath++;

    localIndex_array facePath;
    localIndex_array edgePath;

    facePath.emplace_back( thisFace );
    edgePath.emplace_back( thisEdge );

    // Walk face-edge path until external edge or return to start.
    bool breakFlag = false;
    while( !breakFlag )
    {

      // Get next edge on current face.

      nextEdge = GetOtherFaceEdge( nodeLocalFacesToEdges, thisFace, thisEdge );


      // Backtrack if next edge used but not starting edge.
      if( edgesInPath.count( nextEdge )==1 && nextEdge!=startingEdge )
      {
        // first check to see if we can use the path without the preceding
        return false;
      }

      // Done if external edge reached or loop closed.
      if( (isEdgeExternal[nextEdge]==1 && startingEdgeExternal ) || edgesInPath.count( nextEdge )==1 )
      {
        // Verify valid loop closure.
        if( nextEdge!=startingEdge && !(isEdgeExternal[nextEdge]==1 && startingEdgeExternal ) )
        {
          GEOS_ERROR( GEOS_FMT( "Crap !  NodeID, ParentID = {}, {}\n"
                                "  Starting Edge/Face = {}, {}\n"
                                "  Face Separation Path = {}\n"
                                "  Edge Separation Path = {}\n",
                                nodeID,
                                parentNodeIndex,
                                startingEdge,
                                startingFace,
                                format::toStringForFmt( facePath ),
                                format::toStringForFmt( edgePath ) ),
                      getDataContext() );
        }

        // add faces in the path to separationPathFaces
        for( map< localIndex, int >::const_iterator kf=facesInPath.begin(); kf!=facesInPath.end(); ++kf )
        {
          separationPathFaces.insert( kf->first );
        }

        // break out of the while loop
        breakFlag = true;
      }
      else
      {
        // If started internally but hit boundary, switch direction from start.
        if( isEdgeExternal[nextEdge]==1 )
        {
          nextEdge = nodeLocalFacesToEdges[startingFace].second;
        }

        // I sure hope that this is true!!
        if( edgesToRuptureReadyFaces[nextEdge].size() > 1 )
        {
          // Pick next face attached to edge.


          {
            // Select optimal next face, preferring ruptureState 1.
            bool pathFound = false;

            // Retrieve elements adjacent to current face, handling ghost boundaries.
            const std::pair< CellElementSubRegion const *, localIndex >
            thisElem0 = std::make_pair( &elemManager.getRegion( m_originalFacesToElemRegion[thisFace][0] ).
                                          getSubRegion< CellElementSubRegion >( m_originalFacesToElemSubRegion[thisFace][0] ),
                                        m_originalFacesToElemIndex[thisFace][0] );

            bool const thisFaceHasSecondElem = !isGhostBoundaryFace( thisFace );

            // Evaluate candidate face quality based on geometry and usage.
            int nextFaceQuality = -1;

            for( std::set< localIndex >::const_iterator iter_edgeToFace = edgesToRuptureReadyFaces[nextEdge].begin();
                 iter_edgeToFace!=edgesToRuptureReadyFaces[nextEdge].end(); ++iter_edgeToFace )
            {
              if( *iter_edgeToFace != thisFace )
              {
                pathFound = true;



                const localIndex candidateFaceIndex = *iter_edgeToFace;
                int candidateFaceQuality = 0;


                localIndex candidateEdgeIndex = GetOtherFaceEdge( nodeLocalFacesToEdges, candidateFaceIndex, nextEdge );
                if( candidateEdgeIndex == startingEdge )
                {
                  nextFace = candidateFaceIndex;
                  break;
                }

                // Assess corner-turn quality by checking shared elements.
                bool const candidateHasSecondElem = !isGhostBoundaryFace( candidateFaceIndex );

                if( thisFaceHasSecondElem && candidateHasSecondElem )
                {
                  const std::pair< CellElementSubRegion const *, localIndex >
                  thisElem1 = std::make_pair( &elemManager.getRegion( m_originalFacesToElemRegion[thisFace][1] ).
                                                getSubRegion< CellElementSubRegion >( m_originalFacesToElemSubRegion[thisFace][1] ),
                                              m_originalFacesToElemIndex[thisFace][1] );

                  const std::pair< CellElementSubRegion const *, localIndex >
                  nextElem0 = std::make_pair( &elemManager.getRegion( m_originalFacesToElemRegion[candidateFaceIndex][0] ).
                                                getSubRegion< CellElementSubRegion >( m_originalFacesToElemSubRegion[candidateFaceIndex][0] ),
                                              m_originalFacesToElemIndex[candidateFaceIndex][0] );

                  const std::pair< CellElementSubRegion const *, localIndex >
                  nextElem1 = std::make_pair( &elemManager.getRegion( m_originalFacesToElemRegion[candidateFaceIndex][1] ).
                                                getSubRegion< CellElementSubRegion >( m_originalFacesToElemSubRegion[candidateFaceIndex][1] ),
                                              m_originalFacesToElemIndex[candidateFaceIndex][1] );

                  if( thisElem0 != nextElem0 && thisElem0 != nextElem1 &&
                      thisElem1 != nextElem0 && thisElem1 != nextElem1 )
                  {
                    candidateFaceQuality += 1;
                  }
                }

                if( m_usedFacesForNode[nodeID].count( candidateFaceIndex ) == 0 )
                {
                  candidateFaceQuality += 2;
                }


                if( candidateFaceQuality > nextFaceQuality )
                {
                  nextFace = candidateFaceIndex;
                  nextFaceQuality = candidateFaceQuality;
                }

                if( candidateFaceQuality == 3 )
                {
                  break;
                }
              }
            }
            if( pathFound == false )
            {
              if( edgeOnPartitionBoundary.count( nextEdge ) > 0 )
              {
                return false;
              }
              GEOS_ERROR( "couldn't find the next face in the rupture path (SurfaceGenerator::findFracturePlanes", getDataContext() );
            }
          }

          //        lastEdge = thisEdge;
          //        lastFace = thisFace;

          thisEdge = nextEdge;
          thisFace = nextFace;
          //      separationPathFaces.insert( thisFace );
          edgesInPath[thisEdge] = numFacesInPath;
          facesInPath[thisFace] = numFacesInPath++;

          facePath.emplace_back( thisFace );
          edgePath.emplace_back( thisEdge );

        }
        else
        {
          // Incomplete local path (likely partition boundary) - abort.
          if( edgeOnPartitionBoundary.count( nextEdge ) > 0 )
          {
            return false;
          }
          GEOS_ERROR( "next edge in separation path is apparently  connected to less than 2 ruptured face (SurfaceGenerator::findFracturePlanes", getDataContext() );
        }

      }
    }
  }


  //***** SET LOCATIONS ************************************************************************************************



  // need a map from faces to edges that are attached to the node
  map< localIndex, std::pair< localIndex, localIndex > > localFacesToEdges;
  for( localIndex const kf : nodeToFaceMap[ nodeID ] )
  {
    localIndex edge[2] = { INT_MAX, INT_MAX };
    int count = 0;
    for( auto ke : faceToEdgeMap[ kf ] )
    {
      if( edgeManager.hasNode( ke, nodeID ) )
      {
        edge[count++] = ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      GEOS_ERROR( "invalid edge. (SurfaceGenerator::findFracturePlanes", getDataContext() );
    }


    localFacesToEdges[kf] = std::make_pair( edge[0], edge[1] );

  }


  // now we want to identify the objects on either side of the separation plane. First we assign an array to indicate
  // whether a face/edge is on the fracture plane.

  for( localIndex const kf : nodeToFaceMap[ nodeID ] )
  {
    // iff the face is being split NOW, the set the faceLocation = -1.
    const localIndex virtualFaceIndex = ( parentFaceIndices[kf] == -1 ) ? kf : parentFaceIndices[kf];
    if( kf == virtualFaceIndex && childFaceIndices[kf] == -1 && separationPathFaces.count( kf ) )
    {
      faceLocations[kf] = -1;
    }
    else
    {
      faceLocations[kf] = INT_MIN;
    }

  }
  for( localIndex const edgeID : nodeToEdgeMap[ nodeID ] )
  {
    edgeLocations[edgeID] = INT_MIN;
  }

  for( auto k = nodeToElementMaps.cbegin(); k != nodeToElementMaps.cend(); ++k )
  {
    elemLocations[*k] = INT_MIN;
  }


  // Assign sides (0 or 1) to every element, face, and edge around the
  // splitting node via BFS.  Ghost-boundary faces are handled internally.
  assignLocationsBFS( separationPathFaces,
                      elemManager,
                      faceManager,
                      nodeToElementMaps,
                      localFacesToEdges,
                      edgeLocations,
                      faceLocations,
                      elemLocations );



  bool fail = false;

  for( localIndex const edgeID : nodeToEdgeMap[ nodeID ] )
  {
    if( edgeLocations[edgeID] == INT_MIN )
    {
      fail = true;
    }
  }
  for( localIndex const kf : nodeToFaceMap[ nodeID ] )
  {
    if( faceLocations[kf] == INT_MIN )
    {
      fail = true;
    }
  }

  if( fail )
  {
    globalIndex const nodeGlobalIndex = nodeManager.localToGlobalMap()[nodeID];
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );

    std::ostringstream oss;
    oss << "[Rank " << rank << "] findFracturePlanes FAIL for node local=" << nodeID
        << " global=" << nodeGlobalIndex
        << " | numElems=" << nodeToElementMaps.size()
        << " | numFaces=" << nodeToFaceMap.sizeOfSet( nodeID )
        << " | numEdges=" << nodeToEdgeMap.sizeOfSet( nodeID )
        << " | separationPath={";
    for( localIndex sf : separationPathFaces )
      oss << sf << ",";
    oss << "}";

    // Report unassigned elements
    for( auto const & [ek, loc] : elemLocations )
    {
      if( loc == INT_MIN )
        oss << " | UNASSIGNED elem gi=" << ek.first->localToGlobalMap()[ek.second];
    }
    // Report unassigned faces
    for( localIndex const kf : nodeToFaceMap[ nodeID ] )
    {
      if( faceLocations[kf] == INT_MIN )
        oss << " | UNASSIGNED face li=" << kf << " gi=" << faceManager.localToGlobalMap()[kf];
    }
    // Report unassigned edges
    for( localIndex const edgeID : nodeToEdgeMap[ nodeID ] )
    {
      if( edgeLocations[edgeID] == INT_MIN )
        oss << " | UNASSIGNED edge li=" << edgeID << " gi=" << edgeManager.localToGlobalMap()[edgeID];
    }
    // Check if any face is a domain boundary face
    for( localIndex const kf : nodeToFaceMap[ nodeID ] )
    {
      localIndex const vfi = ( parentFaceIndices[kf] == -1 ) ? kf : parentFaceIndices[kf];
      if( m_originalFacesToElemIndex( vfi, 1 ) == -1 )
        oss << " | DOMAIN_BND face li=" << kf << " gi=" << faceManager.localToGlobalMap()[kf];
    }

    std::cout << oss.str() << std::endl;

    return false;
  }
  return true;
}

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
bool SurfaceGenerator::assignLocationsBFS( std::set< localIndex > const & separationPathFaces,
                                           ElementRegionManager const & elemManager,
                                           FaceManager const & faceManager,
                                           stdVector< std::pair< CellElementSubRegion const *, localIndex > > const & nodeToElementMaps,
                                           map< localIndex, std::pair< localIndex, localIndex > > const & localFacesToEdges,
                                           map< localIndex, int > & edgeLocations,
                                           map< localIndex, int > & faceLocations,
                                           map< std::pair< CellElementSubRegion const *, localIndex >, int > & elemLocations )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< localIndex const > const & parentFaceIndices =
    faceManager.getField< fields::parentIndex >();

  // Propagate assignment to element components and enqueue neighbors.
  using ElemKey = std::pair< CellElementSubRegion const *, localIndex >;

  // BFS queue: each entry is (element, assigned-location).
  std::deque< std::pair< ElemKey, int > > bfsQueue;

  // Seed the BFS with the first element (sorted by global index)
  // on side 0.
  ElemKey const seed = *( nodeToElementMaps.cbegin() );
  elemLocations[seed] = 0;
  bfsQueue.push_back( { seed, 0 } );

  // ----- Main BFS loop --------------------------------------------------
  while( !bfsQueue.empty() )
  {
    auto const [currentElem, location] = bfsQueue.front();
    bfsQueue.pop_front();

    int const otherLocation = ( location == 0 ) ? 1 : 0;

    // Loop over every face of the current element.
    for( localIndex kf = 0; kf < currentElem.first->faceList().size( 1 ); ++kf )
    {
      localIndex const faceIndex = currentElem.first->faceList()( currentElem.second, kf );
      localIndex const virtualFaceIndex = ( parentFaceIndices[faceIndex] == -1 )
                                          ? faceIndex : parentFaceIndices[faceIndex];

      // Only process faces that belong to the splitting node.
      if( faceLocations.find( faceIndex ) == faceLocations.end() )
      {
        continue;
      }

      // --- Assign face location -----------------------------------------
      if( faceLocations[faceIndex] == otherLocation )
      {
        faceLocations[faceIndex] = -1;  // face is on the separation plane
      }
      else if( faceLocations[faceIndex] == INT_MIN )
      {
        faceLocations[faceIndex] = location;
      }

      // --- Assign edge locations ----------------------------------------
      auto iterF2E = localFacesToEdges.find( faceIndex );
      if( iterF2E != localFacesToEdges.end() )
      {
        for( localIndex edgeIdx : { iterF2E->second.first, iterF2E->second.second } )
        {
          if( edgeLocations[edgeIdx] == otherLocation )
          {
            edgeLocations[edgeIdx] = -1;
          }
          else if( edgeLocations[edgeIdx] == INT_MIN )
          {
            edgeLocations[edgeIdx] = location;
          }
        }
      }

      // --- Enqueue the neighbor element across this face ----------------
      // Skip ghost-boundary faces (second element slot == -1).
      if( isGhostBoundaryFace( virtualFaceIndex ) )
      {
        continue;
      }

      // Both element slots are populated — find the neighbor.
      localIndex const er0  = m_originalFacesToElemRegion( virtualFaceIndex, 0 );
      localIndex const esr0 = m_originalFacesToElemSubRegion( virtualFaceIndex, 0 );
      ElemKey const elemKey0 = { &elemManager.getRegion( er0 ).getSubRegion< CellElementSubRegion >( esr0 ),
                                 m_originalFacesToElemIndex( virtualFaceIndex, 0 ) };

      localIndex const er1  = m_originalFacesToElemRegion( virtualFaceIndex, 1 );
      localIndex const esr1 = m_originalFacesToElemSubRegion( virtualFaceIndex, 1 );
      ElemKey const elemKey1 = { &elemManager.getRegion( er1 ).getSubRegion< CellElementSubRegion >( esr1 ),
                                 m_originalFacesToElemIndex( virtualFaceIndex, 1 ) };

      ElemKey const & neighborElem = ( elemKey0 == currentElem ) ? elemKey1 : elemKey0;
      int const neighborLocation =
        ( separationPathFaces.count( virtualFaceIndex ) == 0 ) ? location : otherLocation;

      // Only enqueue if the neighbor is in the splitting-node's element set
      // and has not yet been assigned.
      if( std::find( nodeToElementMaps.cbegin(), nodeToElementMaps.cend(), neighborElem ) != nodeToElementMaps.cend() )
      {
        if( elemLocations[neighborElem] == INT_MIN )
        {
          elemLocations[neighborElem] = neighborLocation;
          bfsQueue.push_back( { neighborElem, neighborLocation } );
        }
      }
    }
  }

  // Fallback sweep: handle unreachable elements by inferring locations from neighbors.
  bool progress = true;
  while( progress )
  {
    progress = false;

    // Collect elements that still need assignment.
    stdVector< ElemKey > unassigned;
    for( auto const & [elem, loc] : elemLocations )
    {
      if( loc == INT_MIN )
      {
        unassigned.push_back( elem );
      }
    }

    for( ElemKey const & elem : unassigned )
    {
      // Skip if a previous iteration in this sweep already assigned it.
      if( elemLocations[elem] != INT_MIN )
        continue;

      // Try to find a same-node neighbor that has been assigned.
      bool found = false;
      for( localIndex kf = 0; kf < elem.first->faceList().size( 1 ) && !found; ++kf )
      {
        localIndex const faceIndex = elem.first->faceList()( elem.second, kf );
        localIndex const virtualFaceIndex = ( parentFaceIndices[faceIndex] == -1 )
                                            ? faceIndex : parentFaceIndices[faceIndex];

        // For a ghost-boundary face we can still check the [0]-side
        // element — it might be the neighbor we are looking for.
        // We just cannot discover elements on the missing [1]-side.

        // Check both sides of the face for an assigned neighbor.
        for( int side = 0; side < 2; ++side )
        {
          localIndex const er  = m_originalFacesToElemRegion( virtualFaceIndex, side );
          localIndex const esr = m_originalFacesToElemSubRegion( virtualFaceIndex, side );
          localIndex const ei  = m_originalFacesToElemIndex( virtualFaceIndex, side );

          if( ei == -1 )
            continue;  // missing element on this side (ghost boundary)

          ElemKey const neighbor = { &elemManager.getRegion( er ).getSubRegion< CellElementSubRegion >( esr ), ei };

          if( neighbor == elem )
            continue;  // that's us, not a neighbor

          auto iterNeighbor = elemLocations.find( neighbor );
          if( iterNeighbor != elemLocations.end() && iterNeighbor->second != INT_MIN )
          {
            bool const crossesSeparation = ( separationPathFaces.count( virtualFaceIndex ) > 0 );
            int const inferredLocation = crossesSeparation
                                         ? ( ( iterNeighbor->second == 0 ) ? 1 : 0 )
                                         : iterNeighbor->second;

            // Assign and seed a mini-BFS from this element.
            elemLocations[elem] = inferredLocation;
            bfsQueue.push_back( { elem, inferredLocation } );
            progress = true;
            found = true;

            // Run the mini-BFS to propagate from the newly assigned element.
            while( !bfsQueue.empty() )
            {
              auto const [curElem, curLoc] = bfsQueue.front();
              bfsQueue.pop_front();
              int const curOther = ( curLoc == 0 ) ? 1 : 0;

              for( localIndex kf2 = 0; kf2 < curElem.first->faceList().size( 1 ); ++kf2 )
              {
                localIndex const fi = curElem.first->faceList()( curElem.second, kf2 );
                localIndex const vfi = ( parentFaceIndices[fi] == -1 ) ? fi : parentFaceIndices[fi];

                if( faceLocations.find( fi ) == faceLocations.end() )
                  continue;

                if( faceLocations[fi] == curOther )
                  faceLocations[fi] = -1;
                else if( faceLocations[fi] == INT_MIN )
                  faceLocations[fi] = curLoc;

                auto iterF2E = localFacesToEdges.find( fi );
                if( iterF2E != localFacesToEdges.end() )
                {
                  for( localIndex eIdx : { iterF2E->second.first, iterF2E->second.second } )
                  {
                    if( edgeLocations[eIdx] == curOther )
                      edgeLocations[eIdx] = -1;
                    else if( edgeLocations[eIdx] == INT_MIN )
                      edgeLocations[eIdx] = curLoc;
                  }
                }

                if( isGhostBoundaryFace( vfi ) )
                  continue;

                localIndex const r0  = m_originalFacesToElemRegion( vfi, 0 );
                localIndex const sr0 = m_originalFacesToElemSubRegion( vfi, 0 );
                ElemKey const ek0 = { &elemManager.getRegion( r0 ).getSubRegion< CellElementSubRegion >( sr0 ),
                                      m_originalFacesToElemIndex( vfi, 0 ) };

                localIndex const r1  = m_originalFacesToElemRegion( vfi, 1 );
                localIndex const sr1 = m_originalFacesToElemSubRegion( vfi, 1 );
                ElemKey const ek1 = { &elemManager.getRegion( r1 ).getSubRegion< CellElementSubRegion >( sr1 ),
                                      m_originalFacesToElemIndex( vfi, 1 ) };

                ElemKey const & nbr = ( ek0 == curElem ) ? ek1 : ek0;
                int const nbrLoc = ( separationPathFaces.count( vfi ) == 0 ) ? curLoc : curOther;

                if( std::find( nodeToElementMaps.cbegin(), nodeToElementMaps.cend(), nbr ) != nodeToElementMaps.cend() )
                {
                  if( elemLocations[nbr] == INT_MIN )
                  {
                    elemLocations[nbr] = nbrLoc;
                    bfsQueue.push_back( { nbr, nbrLoc } );
                  }
                }
              }
            }

            break;  // assigned via this side — stop checking sides
          }
        }
      }
    }
  }

  return true;
}


//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
void SurfaceGenerator::performFracture( const localIndex nodeID,
                                        real64 const time_np1,
                                        NodeManager & nodeManager,
                                        EdgeManager & edgeManager,
                                        FaceManager & faceManager,
                                        ElementRegionManager & elementManager,
                                        ModifiedObjectLists & modifiedObjects,
                                        stdVector< std::set< localIndex > > & GEOS_UNUSED_PARAM( nodesToRupturedFaces ),
                                        stdVector< std::set< localIndex > > & GEOS_UNUSED_PARAM( edgesToRupturedFaces ),
                                        const std::set< localIndex > & separationPathFaces,
                                        const map< localIndex, int > & edgeLocations,
                                        const map< localIndex, int > & faceLocations,
                                        const map< std::pair< CellElementSubRegion const *, localIndex >, int > & elemLocations )
{
  GEOS_MARK_FUNCTION;


  GEOS_MARK_BEGIN( "setup" );

  GEOS_LOG_LEVEL_BY_RANK( logInfo::SurfaceGenerator,
                          GEOS_FMT( "SurfaceGenerator Statistics: Splitting node {}. Node duplicates: 1. Associated Faces (before split): {}. Associated Edges (before split): {}",
                                    nodeID, faceLocations.size(), edgeLocations.size() ) );

  int const rank = MpiWrapper::commRank( MPI_COMM_WORLD );

  array2d< real64, nodes::REFERENCE_POSITION_PERM > const & X = nodeManager.referencePosition();
  ArrayOfSets< localIndex > & nodeToEdgeMap = nodeManager.edgeList();
  ArrayOfSets< localIndex > & nodeToFaceMap = nodeManager.faceList();
  ArrayOfArrays< localIndex > & nodeToRegionMap = nodeManager.elementRegionList();
  ArrayOfArrays< localIndex > & nodeToSubRegionMap = nodeManager.elementSubRegionList();
  ArrayOfArrays< localIndex > & nodeToElementMap = nodeManager.elementList();

  array2d< localIndex > & edgeToNodeMap = edgeManager.nodeList();
  ArrayOfSets< localIndex > & edgeToFaceMap = edgeManager.faceList();

  ArrayOfArrays< localIndex > & faceToNodeMap = faceManager.nodeList();
  ArrayOfArrays< localIndex > & faceToEdgeMap = faceManager.edgeList();
  array2d< localIndex > const & faceToRegionMap = faceManager.elementRegionList();
  array2d< localIndex > const & faceToSubRegionMap = faceManager.elementSubRegionList();
  array2d< localIndex > const & faceToElementMap = faceManager.elementList();

  array1d< integer > const & faceIsExternal = faceManager.isExternal();
  array1d< integer > const & edgeIsExternal = edgeManager.isExternal();
  array1d< integer > const & nodeIsExternal = nodeManager.isExternal();

  SurfaceElementRegion & fractureElementRegion = elementManager.getRegion< SurfaceElementRegion >( m_fractureRegionName );
  array1d< integer > const & isFaceSeparable = faceManager.getField< surfaceGeneration::isFaceSeparable >();

  array2d< real64 > const & faceNormals = faceManager.faceNormal();

  array1d< localIndex > const & parentEdgeIndices = edgeManager.getField< fields::parentIndex >();
  array1d< localIndex > const & childEdgeIndices = edgeManager.getField< fields::childIndex >();
  array1d< localIndex > const & parentNodeIndices = nodeManager.getField< fields::parentIndex >();
  array1d< localIndex > const & childNodeIndices = nodeManager.getField< fields::childIndex >();

  array1d< integer > const & degreeFromCrack = nodeManager.getField< surfaceGeneration::degreeFromCrack >();
  array1d< integer > const & nodeDegreeFromCrackTip = nodeManager.getField< surfaceGeneration::degreeFromCrackTip >();
  array1d< integer > const & faceDegreeFromCrackTip = faceManager.getField< surfaceGeneration::degreeFromCrackTip >();

  array1d< real64 > const & nodeRuptureTime = nodeManager.getField< fields::ruptureTime >();
  array1d< real64 > const & faceRuptureTime = faceManager.getField< fields::ruptureTime >();

  // ***** split all the objects first *****

  // Split the node into two, using the original index, and a new one.
  localIndex newNodeIndex;
  {
    std::ostringstream s;
    for( std::set< localIndex >::const_iterator i=separationPathFaces.begin(); i!=separationPathFaces.end(); ++i )
    {
      s << *i << " ";
    }
    GEOS_LOG_LEVEL_BY_RANK( logInfo::SurfaceGenerator,
                            GEOS_FMT( "Splitting node {} along separation plane faces: {}", nodeID, s.str() ) );
  }

  GEOS_MARK_END( "setup" );


  GEOS_MARK_BEGIN( "node split" );

  nodeManager.splitObject( nodeID, rank, newNodeIndex );

  modifiedObjects.newNodes.insert( newNodeIndex );
  modifiedObjects.modifiedNodes.insert( nodeID );

  nodeToRegionMap.clearArray( newNodeIndex );
  nodeToSubRegionMap.clearArray( newNodeIndex );
  nodeToElementMap.clearArray( newNodeIndex );

  nodeToEdgeMap.clearSet( newNodeIndex );
  nodeToFaceMap.clearSet( newNodeIndex );

  degreeFromCrack[nodeID] = 0;
  degreeFromCrack[newNodeIndex] = 0;
  m_tipNodes.remove( nodeID );
  nodeDegreeFromCrackTip( nodeID ) = 1;
  nodeRuptureTime( nodeID ) = time_np1;
  nodeRuptureTime( newNodeIndex ) = time_np1;

  //TODO HACK...should recalculate mass

  //TODO Either change m_usedFacesForNode to array<std::set> or add insert with iterator to SortedArray
  for( auto const val : separationPathFaces )
  {
    m_usedFacesForNode[nodeID].insert( val );
    m_usedFacesForNode[newNodeIndex].insert( val );
  }

  GEOS_LOG_LEVEL_BY_RANK( logInfo::SurfaceGenerator,
                          GEOS_FMT( "Done splitting node {} into nodes {} and {}", nodeID, nodeID, newNodeIndex ) );

  GEOS_MARK_END( "node split" );

  GEOS_MARK_BEGIN( "edge split" );
  // split edges
  map< localIndex, localIndex > splitEdges;
  // loop over all edges connected to the node
  for( map< localIndex, int >::const_iterator iter_edge=edgeLocations.begin(); iter_edge!=edgeLocations.end(); ++iter_edge )
  {
    const localIndex & parentEdgeIndex = iter_edge->first;
    const int & location = iter_edge->second;

    // if the edge is on the separation plane, then split it
    if( location == -1 )
    {
      localIndex newEdgeIndex;

      edgeManager.splitObject( parentEdgeIndex, rank, newEdgeIndex );

      m_tipEdges.remove( parentEdgeIndex );

      edgeToFaceMap.clearSet( newEdgeIndex );

      GEOS_LOG_LEVEL_BY_RANK( logInfo::SurfaceGenerator,
                              GEOS_FMT ( "Split edge {} into edges {} and {}", parentEdgeIndex, parentEdgeIndex, newEdgeIndex ) );

      splitEdges[parentEdgeIndex] = newEdgeIndex;
      modifiedObjects.newEdges.insert( newEdgeIndex );
      modifiedObjects.modifiedEdges.insert( parentEdgeIndex );

      for( localIndex const faceIndex : edgeToFaceMap[ parentEdgeIndex ] )
      {
        bool trailingFace = false;
        if( m_trailingFaces.contains( faceIndex ))
        {
          for( localIndex const edgeIndex : faceToEdgeMap[ faceIndex ] )
          {
            if( m_tipEdges.contains( edgeIndex ))
            {
              trailingFace = true;
            }
          }

          if( trailingFace == false )
          {
            m_trailingFaces.remove( faceIndex );
          }
        }
      }

      for( int a=0; a<2; ++a )
      {
        edgeManager.nodeList( newEdgeIndex, a ) = edgeManager.nodeList( parentEdgeIndex, a );
      }

    } //    if( location == -1  )
  } // for( map<localIndex,int>::const_iterator iter_edge...


  GEOS_MARK_END( "edge split" );
  GEOS_MARK_BEGIN( "face split" );
  // split the faces
  array1d< integer > const & ruptureState = faceManager.getField< surfaceGeneration::ruptureState >();
  map< localIndex, localIndex > splitFaces;


  SortedArray< localIndex > & externalFaces = faceManager.externalSet();

  // loop over all faces attached to the nodeID
  for( map< localIndex, int >::const_iterator iter_face = faceLocations.begin(); iter_face != faceLocations.end(); ++iter_face )
  {
    const localIndex faceIndex = iter_face->first;
//    localIndex const parentFaceIndex = parentFaceIndices[faceIndex]==faceIndex ? faceIndex :
// parentFaceIndices[faceIndex];
    const int location = iter_face->second;
    // if the face is on the separation plane, then split it
    if( location == -1 )
    {
      localIndex newFaceIndex;

      if( faceManager.splitObject( faceIndex, rank, newFaceIndex ) )
      {
        GEOS_MARK_BEGIN( "splitObjectSubScope" );

        GEOS_LOG_LEVEL_BY_RANK( logInfo::SurfaceGenerator,
                                GEOS_FMT ( "Split face {} into faces {} and {}", faceIndex, faceIndex, newFaceIndex ) );

        GEOS_MARK_BEGIN( "splitObjectSubScope1" );

        splitFaces[faceIndex] = newFaceIndex;
        modifiedObjects.newFaces.insert( newFaceIndex );
        modifiedObjects.modifiedFaces.insert( faceIndex );

        ruptureState[faceIndex] = 2;
        ruptureState[newFaceIndex] = 2;

        faceRuptureTime( faceIndex ) = time_np1;
        faceRuptureTime( newFaceIndex ) = time_np1;


        m_trailingFaces.insert( faceIndex );
        m_tipFaces.remove( faceIndex );
        faceDegreeFromCrackTip( faceIndex ) = 0;
        faceDegreeFromCrackTip( newFaceIndex ) = 0;

        localIndex const numFaceEdges = faceToEdgeMap.sizeOfArray( faceIndex );
        faceToEdgeMap.resizeArray( newFaceIndex, numFaceEdges );
        for( localIndex a = 0; a < numFaceEdges; ++a )
        {
          faceToEdgeMap( newFaceIndex, a ) = faceToEdgeMap( faceIndex, a );
        }

        localIndex const numFaceNodes = faceToNodeMap.sizeOfArray( faceIndex );
        faceToNodeMap.resizeArray( newFaceIndex, numFaceNodes );
        for( localIndex a=0; a<numFaceNodes; ++a )
        {
          localIndex const aa = a == 0 ? a : numFaceNodes - a;
          faceToNodeMap( newFaceIndex, aa ) = faceToNodeMap( faceIndex, a );
        }
        LvArray::tensorOps::scale< 3 >( faceNormals[ newFaceIndex ], -1 );

        externalFaces.insert( newFaceIndex );
        externalFaces.insert( faceIndex );

        GEOS_MARK_END( "splitObjectSubScope1" );
        GEOS_MARK_BEGIN( "splitObjectSubScope2" );
        // Fu: All edges of the parent face should be external now.
        // We have to do the following because isExternal attribute of the tip edge is not handled by the splitter.
        for( localIndex const edgeIndex : faceManager.edgeList()[ faceIndex ] )
        {
          if( parentEdgeIndices[edgeIndex]==-1 && childEdgeIndices[edgeIndex]==-1 )
          {
            m_tipEdges.insert( edgeIndex );

            for( localIndex const iface: edgeManager.faceList()[ edgeIndex ] )
            {
              if( faceToElementMap.size( 1 ) == 2  &&
                  faceIsExternal[iface] < 1 &&
                  checkOrphanElement( elementManager, faceManager, iface ) == 0 &&
                  isFaceSeparable[iface] == 1
//                  && fabs(LvArray::tensorOps::AiBi< 3 >(faceNormals[faceIndex], faceNormals[iface])) > cos( m_maxTurnAngle )
                  )
              {
                m_tipFaces.insert( iface );
              }
            }
          }
          if( edgeIsExternal[edgeIndex]==0 )
          {
            edgeIsExternal[edgeIndex] = 2;
          }
        }
        GEOS_MARK_END( "splitObjectSubScope2" );
        GEOS_MARK_BEGIN( "splitObjectSubScope3" );
        for( localIndex const nodeIndex : faceToNodeMap[ faceIndex ] )
        {
          if( parentNodeIndices[nodeIndex]==-1 && childNodeIndices[nodeIndex]==-1 )
          {
            m_tipNodes.insert( nodeIndex );
            nodeDegreeFromCrackTip( nodeIndex ) = 0;
          }
          if( nodeIsExternal[nodeIndex] )
          {
            nodeIsExternal[nodeIndex] = 2;
          }
        }
        GEOS_MARK_END( "splitObjectSubScope3" );
        GEOS_MARK_BEGIN( "splitObjectSubScope4" );
        {
          localIndex faceIndices[2] = {faceIndex, newFaceIndex};
          localIndex const
          newFaceElement = fractureElementRegion.addToFractureMesh( time_np1,
                                                                    &faceManager,
                                                                    this->m_originalFaceToEdges.toViewConst(),
                                                                    faceIndices );
          m_faceElemsRupturedThisSolve.insert( newFaceElement );
          GEOS_LOG_LEVEL_BY_RANK( logInfo::SurfaceGenerator,
                                  GEOS_FMT ( "Created new FaceElement {} when creating face {} from {}",
                                             newFaceElement, newFaceIndex, faceIndex ) );
          modifiedObjects.newElements[ {fractureElementRegion.getIndexInParent(), 0} ].insert( newFaceElement );
        }
        GEOS_MARK_END( "splitObjectSubScope4" );

        GEOS_MARK_END( "splitObjectSubScope" );

      } // if( faceManager.SplitObject( faceIndex, newFaceIndex ) )
    } // if( location == -1 )
  } // for( map<localIndex,int>::const_iterator iter_face

  GEOS_MARK_END( "face split" );

  GEOS_MARK_BEGIN( "correct relations" );

  // ***** now correct all the relations between the objects *****

  /* To accomplish this annoying yet exceedingly important task, we will take a "top down"
   * approach. Note that this is a two way correction, i.e. if we are correcting
   * elementToNodes, we also correct nodeToElementMaps. This is summarized as:
   * 1) Loop over elements attached to the split node.
   *     2a) correct all relations between the single  element and the nodes.
   *     2b) Loop over all faces on the element
   *         3a) For each face, correct the face relations with the element
   *         3b) For each face, correct the face relations with the nodes
   *         3c) Loop over all edges on the face
   *             4a) For each edge, correct the face relations
   *             4b) for each edge, correct the node relations
   *
   *  The element location will define which side of the rupture everything
   *  is on.
   *  - location 0 gets the original node,edge,face.
   *  - location 1 gets the new node,edge,face.
   */

  array1d< localIndex > const & parentFaceIndex = faceManager.getField< fields::parentIndex >();
  array1d< localIndex > const & childFaceIndex = faceManager.getField< fields::childIndex >();

  // 1) loop over all elements attached to the nodeID
  for( map< std::pair< CellElementSubRegion const *, localIndex >, int >::const_iterator iter_elem = elemLocations.begin(); iter_elem != elemLocations.end(); ++iter_elem )
  {
    const int & location = iter_elem->second;

    if( location == 1 )
    {
      const std::pair< CellElementSubRegion const *, localIndex > & elem = iter_elem->first;

      const CellElementSubRegion & elemSubRegion = *( elem.first );
      const ElementRegionBase & elemRegion = dynamicCast< const ElementRegionBase & >( elemSubRegion.getParent().getParent() );
      string const & elemRegionName = elemRegion.getName();

      localIndex const regionIndex = elementManager.getRegions().getIndex( elemRegionName );
      localIndex const subRegionIndex = elemRegion.getSubRegions().getSubGroupIndex( elemSubRegion.getName() );
      const localIndex elemIndex = elem.second;

      modifiedObjects.modifiedElements[{ regionIndex, subRegionIndex }].insert( elemIndex );


      array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & elemsToNodes = elemSubRegion.nodeList();
      array2d< localIndex > const & elemsToFaces = elemSubRegion.faceList();

      if( getLogLevel() > 1 )
        std::cout<<"Element "<<elemIndex<<std::endl;

      // 2a) correct elementToNode and nodeToElement
      if( getLogLevel() > 1 )
        std::cout<<"  Looping over all nodes on element, and correcting node<->element maps:"<<std::endl;


      real64 elemCenter[3] = {0.0, 0.0, 0.0};
      {
        // loop over all nodes on element
        if( getLogLevel() > 1 )
          std::cout<<"    m_ElementToNodeMap = ( ";
        for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
        {
          LvArray::tensorOps::add< 3 >( elemCenter, X[ elemsToNodes[elemIndex][a] ] );
          // if the node was just split
          if( elemsToNodes[elemIndex][a] == nodeID )
          {

            if( getLogLevel() > 1 )
              std::cout<<elemsToNodes[elemIndex][a]<<"->"<<newNodeIndex<<", ";

            elemsToNodes[elemIndex][a] = newNodeIndex;

            insert( nodeManager.toElementRelation(), newNodeIndex, regionIndex, subRegionIndex, elemIndex );
            erase( nodeManager.toElementRelation(), nodeID, regionIndex, subRegionIndex, elemIndex );
          }
          else if( getLogLevel() > 1 )
            std::cout<<elemsToNodes[elemIndex][a]<<", ";
        }
        LvArray::tensorOps::scale< 3 >( elemCenter, 1.0 / elemsToNodes.size( 1 ) );
        if( getLogLevel() > 1 )
          std::cout<<")"<<std::endl;

        if( getLogLevel() > 1 )
        {
          for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
          {
            if( getLogLevel() > 1 )
            {
              std::cout<<"    nodeToElemMaps["<<elemsToNodes[elemIndex][a]<<"] = ( ";
              for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( elemsToNodes[elemIndex][a] ); ++k )
              {
                std::cout<<"["<<nodeToRegionMap[elemsToNodes[elemIndex][a]][k]<<","
                         <<nodeToSubRegionMap[elemsToNodes[elemIndex][a]][k]<<","
                         <<nodeToElementMap[elemsToNodes[elemIndex][a]][k]<<"] , ";
              }
              std::cout<<" )"<<std::endl;
            }
          }

          if( getLogLevel() > 1 )
          {
            std::cout<<"    nodeToElemMaps["<<nodeID<<"] = ( ";
            for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( nodeID ); ++k )
            {
              std::cout<<"["<<nodeToRegionMap[nodeID][k]<<","
                       <<nodeToSubRegionMap[nodeID][k]<<","
                       <<nodeToElementMap[nodeID][k]<<"] , ";
            }
            std::cout<<" )"<<std::endl;
          }
        }
      }



      // 2b) loop over all faces on element.
      if( getLogLevel() > 1 )
      {
        std::cout<<"  Looping over all faces on element (parent and child):"<<std::endl;
      }

      // we need to build a list of faces that is elemToFaces FOLLOWED by any
      // parent face of those indicated in elemToFaces

      // Now we do a loop over the facelist and process all the faces
      for( int kf=0; kf<elemSubRegion.numFacesPerElement(); ++kf )
      {

        // set both faceID and newFaceID to the parent face.
        localIndex const faceIndex = elemsToFaces[elemIndex][kf];
        //        map<localIndex,localIndex>::iterator iterSplitFace = splitFaces.find(faceIndex);
        bool const isNewFace = (splitFaces.count( faceIndex )>0) ? true : false;
        localIndex const newFaceIndex = isNewFace ? childFaceIndex[faceIndex] : faceIndex;


        // 3a) check to see if the face was split. If so, then we will need
        // to alter the face relation with the elements in both directions.
        if( isNewFace )
        {
          // replace the parent face with the child face in elementToFace. Now
          // faceID is the parent face, and newFaceID is the child face.
          elemsToFaces[elemIndex][kf] = childFaceIndex[faceIndex];

          // printf( "pre map change \n");
          // printf( "    faceToElementMap[%d]  = ( %d, %d, %d ), ( %d, %d, %d )\n", faceIndex,
          //                                                                         faceToRegionMap[faceIndex][0],
          // faceToSubRegionMap[faceIndex][0], faceToElementMap[faceIndex][0],
          //                                                                         faceToRegionMap[faceIndex][1],
          // faceToSubRegionMap[faceIndex][1], faceToElementMap[faceIndex][1] );
          // printf( "    faceToElementMap[%d]  = ( %d, %d, %d ), ( %d, %d, %d )\n", newFaceIndex,
          //                                                                         faceToRegionMap[newFaceIndex][0],
          // faceToSubRegionMap[newFaceIndex][0], faceToElementMap[newFaceIndex][0],
          //                                                                         faceToRegionMap[newFaceIndex][1],
          // faceToSubRegionMap[newFaceIndex][1], faceToElementMap[newFaceIndex][1] );


          // add the element to the child faceToElem
//          faceManager.m_toElements[newFaceIndex].emplace_back( elem );

          faceToRegionMap[newFaceIndex][0] = regionIndex;
          faceToSubRegionMap[newFaceIndex][0] = subRegionIndex;
          faceToElementMap[newFaceIndex][0] = elemIndex;
          faceToRegionMap[newFaceIndex][1] = -1;
          faceToSubRegionMap[newFaceIndex][1] = -1;
          faceToElementMap[newFaceIndex][1] = -1;

          // remove the element from the parent face
          if( faceToRegionMap[faceIndex][0] == regionIndex &&
              faceToSubRegionMap[faceIndex][0] == subRegionIndex &&
              faceToElementMap[faceIndex][0] == elemIndex )
          {
            faceToRegionMap[faceIndex][0] = faceToRegionMap[faceIndex][1];
            faceToSubRegionMap[faceIndex][0] = faceToSubRegionMap[faceIndex][1];
            faceToElementMap[faceIndex][0] = faceToElementMap[faceIndex][1];
            faceToRegionMap[faceIndex][1] = -1;
            faceToSubRegionMap[faceIndex][1] = -1;
            faceToElementMap[faceIndex][1] = -1;
          }
          else if( faceToRegionMap[faceIndex][1] == regionIndex &&
                   faceToSubRegionMap[faceIndex][1] == subRegionIndex &&
                   faceToElementMap[faceIndex][1] == elemIndex )
          {
            faceToRegionMap[faceIndex][1] = -1;
            faceToSubRegionMap[faceIndex][1] = -1;
            faceToElementMap[faceIndex][1] = -1;
          }

          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToRegionMap["<<newFaceIndex<<"][0]    = "<<faceToRegionMap[newFaceIndex][0]
          // );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToSubRegionMap["<<newFaceIndex<<"][0] =
          // "<<faceToSubRegionMap[newFaceIndex][0] );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToElementMap["<<newFaceIndex<<"][0]      =
          // "<<faceToElementMap[newFaceIndex][0] );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToRegionMap["<<newFaceIndex<<"][1]    = "<<faceToRegionMap[newFaceIndex][1]
          // );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToSubRegionMap["<<newFaceIndex<<"][1] =
          // "<<faceToSubRegionMap[newFaceIndex][1] );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToElementMap["<<newFaceIndex<<"][1]      =
          // "<<faceToElementMap[newFaceIndex][1] );

          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToRegionMap["<<faceIndex<<"][0]    = "<<faceToRegionMap[faceIndex][0] );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToSubRegionMap["<<faceIndex<<"][0] = "<<faceToSubRegionMap[faceIndex][0] );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToElementMap["<<faceIndex<<"][0]      = "<<faceToElementMap[faceIndex][0] );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToRegionMap["<<faceIndex<<"][1]    = "<<faceToRegionMap[faceIndex][1] );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToSubRegionMap["<<faceIndex<<"][1] = "<<faceToSubRegionMap[faceIndex][1] );
          // GEOS_LOG_LEVEL_RANK_0( logInfo::Mapping, "    faceToElementMap["<<faceIndex<<"][1]      = "<<faceToElementMap[faceIndex][1] );

          // printf( "post map change \n");
          // printf( "    faceToElementMap[%d]  = ( %d, %d, %d ), ( %d, %d, %d )\n", faceIndex,
          //                                                                         faceToRegionMap[faceIndex][0],
          // faceToSubRegionMap[faceIndex][0], faceToElementMap[faceIndex][0],
          //                                                                         faceToRegionMap[faceIndex][1],
          // faceToSubRegionMap[faceIndex][1], faceToElementMap[faceIndex][1] );
          // printf( "    faceToElementMap[%d]  = ( %d, %d, %d ), ( %d, %d, %d )\n", newFaceIndex,
          //                                                                         faceToRegionMap[newFaceIndex][0],
          // faceToSubRegionMap[newFaceIndex][0], faceToElementMap[newFaceIndex][0],
          //                                                                         faceToRegionMap[newFaceIndex][1],
          // faceToSubRegionMap[newFaceIndex][1], faceToElementMap[newFaceIndex][1] );


          for( int i = 0; i < 2; i++ )
          {
            localIndex iFace = i == 0 ? faceIndex : newFaceIndex;

            localIndex elementIndex = faceToElementMap[iFace][0];
            CellElementSubRegion & elementSubRegion = elementManager.getRegion( faceToRegionMap[iFace][0] ).
                                                        getSubRegion< CellElementSubRegion >( faceToSubRegionMap[iFace][0] );
            arrayView2d< real64 const > const subRegionElemCenter = elementSubRegion.getElementCenter();

            FaceManager::sortFaceNodes( X, subRegionElemCenter[ elementIndex ], faceToNodeMap[ iFace ] );

            //Face normal need to be updated here
            real64 fCenter[ 3 ];
            computationalGeometry::centroid_3DPolygon( faceToNodeMap[ iFace ],
                                                       X,
                                                       fCenter,
                                                       faceNormals[ iFace ] );
          }

        } // if( splitFaces.count( faceID ) > 0 )

        modifiedObjects.modifiedFaces.insert( faceIndex );



        // 3b) correct faceToNodes and nodeToFaces

        if( getLogLevel() > 1 )
        {
          localIndex const parentFace = parentFaceIndex[newFaceIndex];
          if( parentFace!=-1 )
          {
            std::cout<<"    m_FaceToNodeMap["<<parentFace<<"->"<<newFaceIndex<<"] = ( ";
          }
          else
          {
            std::cout<<"    m_FaceToNodeMap["<<newFaceIndex<<"] = ( ";
          }
        }

        // loop over all nodes on the face.
        for( localIndex & nodeIndex : faceToNodeMap[ newFaceIndex ] )
        {
          if( getLogLevel() > 1 )
            std::cout<<nodeIndex;

          // if the facenode is the one that is being split
          if( nodeIndex == nodeID )
          {
            nodeIndex = newNodeIndex;

            // if it is not a new face.
            if( !isNewFace )
            {
              // remove the face from the nodeToFaceMap of the parent node.
              nodeToFaceMap.removeFromSet( nodeID, faceIndex );

              // add the face to the nodeToFaceMap of the new node.
              nodeToFaceMap.insertIntoSet( nodeIndex, faceIndex );
            }
            else
            {
              // it is a new face

              // insert the newFace into the nodeToFaceMap of the newNode
              nodeToFaceMap.insertIntoSet( nodeIndex, newFaceIndex );
            }
            if( getLogLevel() > 1 )
              std::cout<<"->"<<nodeIndex<<", ";
          }
          else // the node is not being split
          {
            nodeToFaceMap.insertIntoSet( nodeIndex, newFaceIndex );

            if( getLogLevel() > 1 )
              std::cout<<", ";
          }

        }
        if( getLogLevel() > 1 )
          std::cout<<")"<<std::endl;



        // faceToEdges
        if( getLogLevel() > 1 )
        {
          const localIndex parentFace = parentFaceIndex[newFaceIndex];
          if( parentFace!=-1 )
          {
            std::cout<<"    m_FaceToEdgeMap["<<parentFace<<"->"<<newFaceIndex<<"] = ( ";
          }
          else
          {
            std::cout<<"    m_FaceToEdgeMap["<<newFaceIndex<<"] = ( ";
          }
        }
        // loop over all edges on face
        for( localIndex & edgeIndex : faceToEdgeMap[ newFaceIndex ] )
        {

          // if the edge was just split
          if( splitEdges.count( edgeIndex ) > 0 )
          {
            if( faceIndex == newFaceIndex )
            {
              edgeToFaceMap.removeFromSet( edgeIndex, faceIndex );
            }

            edgeIndex = splitEdges[edgeIndex];
          }
          edgeToFaceMap.insertIntoSet( edgeIndex, newFaceIndex );

          modifiedObjects.modifiedEdges.insert( edgeIndex );

          if( getLogLevel() > 1 )
            std::cout<<edgeIndex;



          //edgeToNodeMap
          if( getLogLevel() > 1 )
          {
            std::cout<<"(";
          }

          {
            for( localIndex a=0; a<edgeToNodeMap.size( 1 ); ++a )
            {
              if( edgeToNodeMap[edgeIndex][a] == nodeID )
              {

                if( getLogLevel() > 1 )
                  std::cout<<edgeToNodeMap[edgeIndex][a];

                edgeToNodeMap[edgeIndex][a] = newNodeIndex;
                nodeToEdgeMap.removeFromSet( nodeID, edgeIndex );

                if( getLogLevel() > 1 )
                  std::cout<<"->"<<edgeToNodeMap[edgeIndex][a]<<", ";

              }
              else if( getLogLevel() > 1 )
                std::cout<<edgeToNodeMap[edgeIndex][a]<<", ";

              nodeToEdgeMap.insertIntoSet( edgeToNodeMap[edgeIndex][a], edgeIndex );
              modifiedObjects.modifiedNodes.insert( edgeToNodeMap[edgeIndex][a] );
            }
            if( getLogLevel() > 1 )
              std::cout<<")";
          }
          if( getLogLevel() > 1 )
            std::cout<<", ";
        }
        if( getLogLevel() > 1 )
          std::cout<<")"<<std::endl;
      } // for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
    } // if( location==1 )
  } // for( map<std::pair<CellBlockSubRegion*, localIndex>, int>::const_iterator iter_elem = elemLocations.begin()

  GEOS_MARK_END( "correct relations" );
}


void SurfaceGenerator::mapConsistencyCheck( localIndex const GEOS_UNUSED_PARAM( nodeID ),
                                            NodeManager const & nodeManager,
                                            EdgeManager const & edgeManager,
                                            FaceManager const & faceManager,
                                            ElementRegionManager const & elementManager,
                                            map< std::pair< CellElementSubRegion const *, localIndex >, int > const & elemLocations )
{
  GEOS_MARK_FUNCTION;

  //**************************************************************************
  // THIS IS ALL JUST CONSISTENCY CHECKING
  //**************************************************************************


  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = nodeManager.edgeList().toViewConst();
  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = nodeManager.faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToRegionMap = nodeManager.elementRegionList();
  ArrayOfArraysView< localIndex const > const & nodeToSubRegionMap = nodeManager.elementSubRegionList();
  ArrayOfArraysView< localIndex const > const & nodeToElementMap = nodeManager.elementList();


  arrayView2d< localIndex > const & edgeToNodeMap = edgeManager.nodeList();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const & faceToRegionMap = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & faceToSubRegionMap = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();

#if 1
  if( getLogLevel() > 2 )
  {
    std::cout << "CONSISTENCY CHECKING OF THE MAPS" << std::endl;

    for( map< std::pair< CellElementSubRegion const *, localIndex >, int >::const_iterator iter_elem = elemLocations.cbegin(); iter_elem != elemLocations.cend(); ++iter_elem )
    {
      const std::pair< CellElementSubRegion const *, localIndex > & elem = iter_elem->first;

      const CellElementSubRegion & elemSubRegion = *( elem.first );
      const localIndex elemIndex = elem.second;

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elemSubRegion.nodeList();
      arrayView2d< localIndex const > const & elemsToFaces = elemSubRegion.faceList();


      std::set< localIndex > elemNodes;


      GEOS_LOG( "Element " << elemIndex );
      std::cout << " elementToNodes = ";
      for( int a=0; a<8; ++a )
      {
        elemNodes.insert( elemsToNodes( elemIndex, a ));
        std::cout << elemsToNodes( elemIndex, a )<<", ";
      }
      std::cout << std::endl;

      std::cout << " elementToFaces->edges->nodes = ";


      // Now we do a loop over the facelist and process all the faces
      for( int kf=0; kf<elemSubRegion.numFacesPerElement(); ++kf )
      {
        std::set< localIndex > faceNodes;

        localIndex faceIndex  = elemsToFaces( elemIndex, kf );

        std::cout << "                              = ";
        std::cout << faceIndex << "( ";
        for( int b=0; b<4; ++b )
        {
          localIndex faceNodeID = faceToNodeMap( faceIndex, b );
          faceNodes.insert( faceNodeID );
          if( elemNodes.count( faceNodeID ) == 0 )
          {
            std::cout << "*";
          }
          std::cout << faceNodeID << ",";
        }
        std::cout << " )      ";



        std::cout << faceIndex << "[ ";
        for( int b=0; b<4; ++b )
        {
          localIndex edgeIndex = faceToEdgeMap( faceIndex, b );
          std::cout << edgeIndex << "( ";
          for( int c=0; c<2; ++c )
          {
            localIndex edgeNodeID = edgeToNodeMap( edgeIndex, c );
            if( elemNodes.count( edgeNodeID ) == 0  && kf<elemSubRegion.numFacesPerElement() )
            {
              std::cout << "*";
            }
            if( faceNodes.count( edgeNodeID ) == 0 )
            {
              std::cout << "#";
            }
            std::cout << edgeNodeID << ",";
          }
          std::cout << " ), ";
        }
        std::cout << " ] \n";

      }
      std::cout << std::endl;

    }

  }

  if( getLogLevel() > 2 )
  {
    // nodeToEdge
    stdVector< std::set< localIndex > > inverseEdgesToNodes( nodeManager.size() );

    for( localIndex ke=0; ke<edgeManager.size(); ++ke )
    {
      for( localIndex b= 0; b<edgeToNodeMap.size( 1 ); ++b )
      {
        localIndex nodeIndex = edgeToNodeMap( ke, b );
        inverseEdgesToNodes[nodeIndex].insert( ke );
      }
    }
    std::cout << "Check NodeToEdge:  nodeToEdgeMap  inverseEdgesToNodes" << std::endl;
    for( localIndex a=0; a<nodeManager.size(); ++a )
    {
      std::cout << "nodeToEdgeMap[" << a << "] = ( ";

      for( localIndex const edgeID : nodeToEdgeMap[ a ] )
      {
        if( inverseEdgesToNodes[a].count( edgeID ) == 0 )
        {
          std::cout << "*";
        }
        std::cout << edgeID << ", ";
      }

      std::cout<<")    (";

      for( localIndex const edgeID : inverseEdgesToNodes[a] )
      {
        if( !nodeToEdgeMap.contains( a, edgeID ) )
          std::cout << "*";
        std::cout << edgeID <<", ";
      }
      std::cout<< ")" <<std::endl;
    }
  }

  if( getLogLevel() > 2 )
  {
    // nodeToFace
    stdVector< std::set< localIndex > > inverseFacesToNodes( nodeManager.size() );
    for( localIndex kf=0; kf<faceManager.size(); ++kf )
    {
      for( localIndex const b : faceToNodeMap[ kf ] )
      {
        inverseFacesToNodes[b].insert( kf );
      }
    }
    std::cout << "Check NodeToFace:  nodeToFaceMap  inverseFacesToNodes" << std::endl;
    for( localIndex a=0; a<nodeManager.size(); ++a )
    {
      std::cout << "m_nodeToFaceMap[ "<< a << "] = ( ";
      for( localIndex const & faceID : nodeToFaceMap[ a ] )
      {
        if( inverseFacesToNodes[a].count( faceID ) == 0 )
          std::cout << "*";
        std::cout << faceID << ", ";
      }
      std::cout<<")    (";

      for( localIndex const edgeID : inverseFacesToNodes[a] )
      {
        if( !nodeToFaceMap.contains( a, edgeID ) )
          std::cout << "*";
        std::cout << edgeID << ", ";
      }
      std::cout<<")"<<std::endl;
    }
  }



  if( getLogLevel() > 2 )
  {


    // nodeToElement
    stdVector< std::set< std::pair< CellElementSubRegion const *, localIndex > > > inverseElemsToNodes( nodeManager.size() );
    elementManager.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion const & subRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = subRegion.nodeList();
      for( localIndex k=0; k<subRegion.size(); ++k )
      {
        for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
        {
          inverseElemsToNodes[elemsToNodes( k, a )].emplace( &subRegion, k );
        }
      }
    } );

    std::cout<<"Check NodeToElem: nodesToElems  inverseElemsToNodes "<<std::endl;


    for( localIndex a=0; a<nodeManager.size(); ++a )
    {

      std::set< std::pair< CellElementSubRegion const *, localIndex > > nodeToElements;
      for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( a ); ++k )
      {
        if( nodeToRegionMap[a][k]!=-1 && nodeToSubRegionMap[a][k]!=-1 && nodeToElementMap[a][k]!=-1 )
        {
          nodeToElements.emplace( &elementManager.getRegion( nodeToRegionMap( a, k ) ).
                                    getSubRegion< CellElementSubRegion >( nodeToSubRegionMap( a, k ) ),
                                  nodeToElementMap( a, k ) );
        }
      }


      std::cout<<"m_NodeToElementMap["<<a<<"] = ( ";
      for( std::set< std::pair< CellElementSubRegion const *, localIndex > >::iterator
           ielem=nodeToElements.begin(); ielem!=nodeToElements.end(); ++ielem )
      {
        if( inverseElemsToNodes[a].count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")    (";

      for( std::set< std::pair< CellElementSubRegion const *, localIndex > >::const_iterator
           ielem=inverseElemsToNodes[a].begin();
           ielem!=inverseElemsToNodes[a].end(); ++ielem )
      {
        if( nodeToElements.count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")"<<std::endl;
    }


    // edgeToFace
    stdVector< std::set< localIndex > > inverseFacesToEdges( edgeManager.size() );
    for( localIndex kf=0; kf<faceManager.size(); ++kf )
    {
      for( localIndex const b : faceToEdgeMap[ kf ] )
      {
        inverseFacesToEdges[ b ].insert( kf );
      }
    }
    std::cout<<"Check EdgeToFace: edgeToFaceMap  inverseFacesToEdges "<<std::endl;
    for( localIndex ke=0; ke<edgeManager.size(); ++ke )
    {
      std::cout<<"m_edgeToFaceMap["<<ke<<"] = ( ";
      for( localIndex const faceID : edgeManager.faceList()[ ke ] )
      {
        if( inverseFacesToEdges[ke].count( faceID ) == 0 )
          std::cout << "*";
        std::cout<<faceID<<", ";
      }
      std::cout<<")    (";

      for( std::set< localIndex >::const_iterator iface=inverseFacesToEdges[ke].begin();
           iface!=inverseFacesToEdges[ke].end(); ++iface )
      {
        if( !edgeManager.faceList().contains( ke, *iface ) )
          std::cout<<"*";
        std::cout<< *iface <<", ";
      }
      std::cout<<")"<<std::endl;
    }

    // faceToElement
    stdVector< std::set< std::pair< CellElementSubRegion const *, localIndex > > > inverseElemsToFaces( faceManager.size() );
    elementManager.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion const & subRegion )
    {
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      for( localIndex k=0; k<subRegion.size(); ++k )
      {
        for( localIndex a=0; a<elemsToFaces.size( 1 ); ++a )
        {
          const localIndex faceID = elemsToFaces( k, a );
          inverseElemsToFaces[ faceID ].emplace( &subRegion, k );

          //            if( parentFaceIndex[faceID] != -1 )
          //            {
          //              inverseElemsToFaces[parentFaceIndex[faceID]].insert(elem);
          //            }
        }
      }
    } );

    std::cout<<"Check FacesToElem: facesToElems  inverseElemsToFaces "<<std::endl;
    for( localIndex a=0; a<faceManager.size(); ++a )
    {

      stdVector< std::pair< CellElementSubRegion const *, localIndex > > faceToElements;
      for( localIndex k=0; k<faceToRegionMap.size( 1 ); ++k )
      {
        // TODO This only works for a single region
        if( faceToRegionMap( a, k ) != -1 )
        {
          faceToElements.emplace_back( &elementManager.getRegion( faceToRegionMap( a, k ) ).
                                         getSubRegion< CellElementSubRegion >( faceToSubRegionMap( a, k ) ),
                                       faceToElementMap( a, k ) );
        }
      }


      std::cout<<"m_FaceToElementMap["<<a<<"] = ( ";

      for( stdVector< std::pair< CellElementSubRegion const *, localIndex > >::const_iterator
           ielem=faceToElements.begin();
           ielem!=faceToElements.end(); ++ielem )
      {
        if( inverseElemsToFaces[a].count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")    (";

      for( std::set< std::pair< CellElementSubRegion const *, localIndex > >::const_iterator ielem=inverseElemsToFaces[a].begin();
           ielem!=inverseElemsToFaces[a].end(); ++ielem )
      {

        if( faceToElements.size() == 2 )
        {
          if( (faceToElements[0] != *ielem) && (faceToElements[1] != *ielem) )
            std::cout<<"*";
        }
        else if( faceToElements.size() )
        {
          if( (faceToElements[0] != *ielem)  )
            std::cout<<"*";
        }
        else
        {
          std::cout<<"****";
        }


        std::cout<<ielem->second<<", ";
      }
      std::cout<<")"<<std::endl;
    }
  }
//  CorrectSplitNodalMass(nodeManager, nodeID, nodeManager.m_childIndices[nodeID][0]);
#endif
}



real64 SurfaceGenerator::calculateKinkAngle( localIndex const edgeID,
                                             NodeManager const & GEOS_UNUSED_PARAM( nodeManager ),
                                             EdgeManager const & edgeManager,
                                             FaceManager const & faceManager )
{
  // TODO: This method should be re-implemented.
  localIndex_array faces;
  // real64 kinkAngle;

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();

  for( localIndex const iface : edgeManager.faceList()[ edgeID ] )
  {
    if( faceIsExternal[iface] == 1 )
      faces.emplace_back( iface );
  }

  if( faces.size() != 2 )
  {
    return(-1.0);
  }
  else
    return 1e100;
}

void SurfaceGenerator::calculateKinkAngles( FaceManager const & faceManager,
                                            EdgeManager & edgeManager,
                                            NodeManager const & nodeManager,
                                            ModifiedObjectLists const & modifiedObjects,
                                            bool const prefrac )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 > & kinkAngle = edgeManager.getReference< real64_array >( "kinkAngle" );

  if( prefrac )
  {
    for( localIndex edgeID = 0; edgeID < edgeManager.size(); ++edgeID )
    {
      kinkAngle[edgeID] = calculateKinkAngle( edgeID, nodeManager, edgeManager, faceManager );
    }
  }
  else
  {
    for( std::set< localIndex >::const_iterator i = modifiedObjects.newEdges.cbegin(); i != modifiedObjects.newEdges.cend(); ++i )
    {
      kinkAngle[*i] = calculateKinkAngle( *i, nodeManager, edgeManager, faceManager );
    }
    for( std::set< localIndex >::const_iterator i = modifiedObjects.modifiedEdges.cbegin(); i != modifiedObjects.modifiedEdges.cend(); ++i )
    {
      kinkAngle[*i] = calculateKinkAngle( *i, nodeManager, edgeManager, faceManager );
    }
  }
}


void SurfaceGenerator::identifyRupturedFaces( DomainPartition const & domain,
                                              NodeManager & nodeManager,
                                              EdgeManager & edgeManager,
                                              FaceManager & faceManager,
                                              ElementRegionManager const & elementManager,
                                              const bool prefrac )
{
  GEOS_MARK_FUNCTION;
  // We use the color map scheme because we can mark a face to be rupture ready from a partition
  // where the face is a ghost.

  if( !m_nodeBasedSIF )
  {
    arrayView1d< integer > const & isEdgeGhost = edgeManager.ghostRank();
    ModifiedObjectLists modifiedObjects;
    {
      for( localIndex iEdge = 0; iEdge != edgeManager.size(); ++iEdge )
      {

        if( isEdgeGhost[iEdge] < 0 )
        {
          int edgeMode = checkEdgeSplitability( iEdge,
                                                nodeManager,
                                                faceManager,
                                                edgeManager,
                                                prefrac );
          if( edgeMode == 0 || edgeMode == 1 ) // We need to calculate SIF
          {
            real64 vecTipNorm[3], vecTip[3];
            localIndex trailFaceID = 0;
            real64 const SIF = calculateEdgeSif( domain, iEdge, trailFaceID,
                                                 nodeManager,
                                                 edgeManager,
                                                 faceManager,
                                                 elementManager,
                                                 vecTipNorm,
                                                 vecTip );

            if( SIF > minimumToughnessOnEdge( iEdge, nodeManager, edgeManager, faceManager ) * 0.5 ) // && edgeMode == 1)
            {
              markRuptureFaceFromEdge( iEdge, trailFaceID,
                                       nodeManager,
                                       edgeManager,
                                       faceManager,
                                       elementManager,
                                       vecTipNorm,
                                       vecTip,
                                       modifiedObjects,
                                       edgeMode );
            }
          }
        }
      }
    }
  }
  else
  {
    ModifiedObjectLists modifiedObjects;

    calculateNodeAndFaceSif( domain, nodeManager, edgeManager, faceManager, elementManager );
    arrayView1d< real64 const > const & SIFNode = nodeManager.getField< surfaceGeneration::SIFNode >();

    for( auto nodeIndex: m_tipNodes )
    {
      if( SIFNode[nodeIndex] > minimumToughnessOnNode( nodeIndex, nodeManager, edgeManager, faceManager ))
      {
        markRuptureFaceFromNode( nodeIndex,
                                 nodeManager,
                                 edgeManager,
                                 faceManager,
                                 elementManager,
                                 modifiedObjects );
      }
    }
  }
}

void SurfaceGenerator::calculateNodeAndFaceSif( DomainPartition const & domain,
                                                NodeManager & nodeManager,
                                                EdgeManager const & edgeManager,
                                                FaceManager & faceManager,
                                                ElementRegionManager const & elementManager )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 > const & SIFNode = nodeManager.getField< surfaceGeneration::SIFNode >();
  arrayView1d< real64 > const & SIFonFace = faceManager.getField< surfaceGeneration::SIFonFace >();

  stdVector< stdVector< real64 > > SIFNode_All, SIFonFace_All;
  stdVector< real64 > SIFOnEdge;
  SIFNode_All.resize( nodeManager.size() );
  SIFonFace_All.resize( faceManager.size() );
  SIFOnEdge.resize( edgeManager.size() );


  SIFNode.zero();
  SIFonFace.zero();

  arrayView2d< real64 const > const & fext = nodeManager.getField< fields::solidMechanics::externalForce >();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & displacement =
    nodeManager.getField< fields::solidMechanics::totalDisplacement >();
  ArrayOfArraysView< localIndex const > const & nodeToRegionMap = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToSubRegionMap = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementMap = nodeManager.elementList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = nodeManager.edgeList().toViewConst();
  arrayView1d< integer const > const & isNodeGhost = nodeManager.ghostRank();

  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList();
  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

  arrayView1d< localIndex const > const & childFaceIndices = faceManager.getField< fields::childIndex >();
  arrayView1d< localIndex const > const & childNodeIndices = nodeManager.getField< fields::childIndex >();
  arrayView1d< localIndex const > const & parentNodeIndices = nodeManager.getField< fields::parentIndex >();

  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();
  m_solidMaterialFullIndex.resize( elementManager.numRegions() );
  elementManager.forElementRegionsComplete< CellElementRegion >( [&]( localIndex regionIndex,
                                                                      CellElementRegion const & region )
  {
    string const & solidMaterialName = region.getSubRegion( 0 ).getReference< string >( viewKeyStruct::solidMaterialNameString() );
    ConstitutiveBase const & solid = constitutiveManager.getConstitutiveRelation< ConstitutiveBase >( solidMaterialName );
    m_solidMaterialFullIndex[regionIndex] = solid.getIndexInParent();
  } );

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elementManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "shearModulus", constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elementManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "bulkModulus", constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView3d< real64 const, geos::solid::STRESS_USD > > const stress =
    elementManager.constructFullMaterialViewAccessor< array3d< real64, geos::solid::STRESS_PERMUTATION >,
                                                      arrayView3d< real64 const, geos::solid::STRESS_USD > >( fields::solid::stress::key(),
                                                                                                              constitutiveManager );
  nodeManager.getField< fields::solidMechanics::totalDisplacement >().move( hostMemorySpace, false );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const &,
                                                               string_array const & regionNames )
  {
    elementManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      CellElementSubRegion const & subRegion )
    {
      string const & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidMaterialNameString() );
      subRegion.
        getConstitutiveModel( solidMaterialName ).
        getReference< array3d< real64, geos::solid::STRESS_PERMUTATION > >( fields::solid::stress::key() ).move( hostMemorySpace,
                                                                                                                 false );
    } );
    displacement.move( hostMemorySpace, false );
  } );

  // auto nodalForceKernel = surfaceGenerationKernels::createKernel( elementManager, constitutiveManager,
  //                                                                 viewKeyStruct::solidMaterialNameString(), false );

  surfaceGenerationKernels::kernelSelector( elementManager,
                                            constitutiveManager,
                                            viewKeyStruct::solidMaterialNameString(),
                                            m_isPoroelastic, [&] ( auto nodalForceKernel )
  {
    for( localIndex const trailingFaceIndex : m_trailingFaces )
    {
      //  RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, m_trailingFaces.size() ), [=] GEOS_HOST_DEVICE (
      // localIndex const trailingFacesCounter )
      //    localIndex const trailingFaceIndex = m_trailingFaces[ trailingFacesCounter ];
      /// TODO: check if a ghost face still has the correct attributes such as normal vector, face center, face index.

      real64 const faceNormalVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[trailingFaceIndex] );
      localIndex_array unpinchedNodeID;
      localIndex_array pinchedNodeID;
      localIndex_array tipEdgesID;

      for( localIndex const nodeIndex : faceToNodeMap[ trailingFaceIndex ] )
      {
        if( m_tipNodes.contains( nodeIndex ))
        {
          pinchedNodeID.emplace_back( nodeIndex );
        }
        else
        {
          unpinchedNodeID.emplace_back( nodeIndex );
        }
      }

      for( localIndex const edgeIndex : faceToEdgeMap[ trailingFaceIndex ] )
      {
        if( m_tipEdges.contains( edgeIndex ))
        {
          tipEdgesID.emplace_back( edgeIndex );
        }
      }

      if( unpinchedNodeID.size() < 2 || (unpinchedNodeID.size() == 2 && tipEdgesID.size() < 2) )
      {
        for( localIndex const nodeIndex : pinchedNodeID )
        {
          if( isNodeGhost[nodeIndex] < 0 )
          {
            real64 nodeDisconnectForce[3] = { 0 };
            real64 const nodePosition[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[nodeIndex] );
            localIndex tralingNodeID = std::numeric_limits< localIndex >::max();
            localIndex nElemEachSide[2];
            nElemEachSide[0] = 0;
            nElemEachSide[1] = 0;

            for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( nodeIndex ); ++k )
            {
              localIndex const er  = nodeToRegionMap[nodeIndex][k];
              localIndex const esr = nodeToSubRegionMap[nodeIndex][k];
              localIndex const ei  = nodeToElementMap[nodeIndex][k];

              CellElementSubRegion const & elementSubRegion = elementManager.getRegion( er ).getSubRegion< CellElementSubRegion >( esr );

              arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elementsToNodes = elementSubRegion.nodeList();
              arrayView2d< real64 const > const & elementCenter = elementSubRegion.getElementCenter().toViewConst();

              for( localIndex n=0; n<elementsToNodes.size( 1 ); ++n )
              {
                if( elementsToNodes( ei, n ) == nodeIndex )
                {
                  real64 nodalForce[ 3 ] = {0};
                  real64 xEle[ 3 ]  = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( elementCenter[ei] );

                  nodalForceKernel.calculateSingleNodalForce( er, esr, ei, n, nodalForce );

                  LvArray::tensorOps::subtract< 3 >( xEle, nodePosition );
                  if( LvArray::tensorOps::AiBi< 3 >( xEle, faceNormalVector ) > 0 ) //TODO: check the sign.
                  {
                    nElemEachSide[0] += 1;
                    LvArray::tensorOps::add< 3 >( nodeDisconnectForce, nodalForce );
                  }
                  else
                  {
                    nElemEachSide[1] +=1;
                    LvArray::tensorOps::subtract< 3 >( nodeDisconnectForce, nodalForce );
                  }
                }
              }
            }

            if( nElemEachSide[0]>=1 && nElemEachSide[1]>=1 )
            {
              LvArray::tensorOps::scale< 3 >( nodeDisconnectForce, 0.5 );
            }

            // Find trailing node when no unpinched nodes exist (Tet mesh case).
            if( unpinchedNodeID.size() == 0 )
            {
              for( localIndex const edgeIndex: faceToEdgeMap[ trailingFaceIndex ] )
              {
                for( localIndex const faceIndex: edgeToFaceMap[ edgeIndex ] )
                {
                  if( faceIndex != trailingFaceIndex && m_tipFaces.contains( faceIndex ))
                  {
                    for( localIndex const iNode: faceToNodeMap[ faceIndex ] )
                    {
                      if( !m_tipNodes.contains( iNode ))
                      {
                        tralingNodeID = iNode;
                      }
                    }
                  }
                }
              }

              if( tralingNodeID == std::numeric_limits< localIndex >::max())
              {
                GEOS_ERROR( "The triangular trailing face has three tip nodes but cannot find the other trailing face containing the trailing node.",
                            getDataContext() );
              }
            }
            else if( unpinchedNodeID.size() == 1 )
            {
              tralingNodeID = unpinchedNodeID[0];
            }
            else if( unpinchedNodeID.size() == 2 )
            {
              for( localIndex const edgeIndex : nodeToEdgeMap[ nodeIndex ] )
              {
                auto const faceToEdgeMapIterator = faceToEdgeMap[ trailingFaceIndex ];
                if( std::find( faceToEdgeMapIterator.begin(), faceToEdgeMapIterator.end(), edgeIndex ) != faceToEdgeMapIterator.end() &&
                    !m_tipEdges.contains( edgeIndex ) )
                {
                  tralingNodeID = edgeToNodeMap[edgeIndex][0] == nodeIndex ? edgeToNodeMap[edgeIndex][1] : edgeToNodeMap[edgeIndex][0];
                }
              }
            }

            //Calculate SIF for the node.
            real64 tipNodeSIF;
            real64 tipNodeForce[3];
            real64 trailingNodeDisp[3];
            localIndex theOtherTrailingNodeID;

            if( childNodeIndices[tralingNodeID] == -1 )
            {
              theOtherTrailingNodeID = parentNodeIndices[tralingNodeID];
            }
            else
            {
              theOtherTrailingNodeID = childNodeIndices[tralingNodeID];
            }

            LvArray::tensorOps::copy< 3 >( trailingNodeDisp, displacement[theOtherTrailingNodeID] );
            LvArray::tensorOps::subtract< 3 >( trailingNodeDisp, displacement[tralingNodeID] );

            //Calculate average young's modulus and poisson ratio for fext.
            real64 fExternal[2][3];
            for( localIndex i=0; i<2; ++i )
            {
              real64 averageYoungModulus( 0 ), averagePoissonRatio( 0 );
              localIndex nodeID = i == 0 ? tralingNodeID : theOtherTrailingNodeID;
              for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( nodeID ); ++k )
              {
                localIndex const er  = nodeToRegionMap[nodeIndex][k];
                localIndex const esr = nodeToSubRegionMap[nodeIndex][k];
                localIndex const ei  = nodeToElementMap[nodeIndex][k];

                real64 K = bulkModulus[er][esr][m_solidMaterialFullIndex[er]][ei];
                real64 G = shearModulus[er][esr][m_solidMaterialFullIndex[er]][ei];
                averageYoungModulus += 9 * K * G / ( 3 * K + G );
                averagePoissonRatio += ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );
              }

              averageYoungModulus /= nodeToRegionMap.sizeOfArray( nodeID );
              averagePoissonRatio /= nodeToRegionMap.sizeOfArray( nodeID );

              LvArray::tensorOps::copy< 3 >( fExternal[i], fext[nodeID] );
              LvArray::tensorOps::scale< 3 >( fExternal[i], averageYoungModulus / (1 - averagePoissonRatio * averagePoissonRatio) );
            }

            //TODO: The sign of fext here is opposite to the sign of fFaceA in function "CalculateEdgeSIF".
            tipNodeForce[0] = nodeDisconnectForce[0] - ( fExternal[0][0] - fExternal[1][0] ) / 2.0;
            tipNodeForce[1] = nodeDisconnectForce[1] - ( fExternal[0][1] - fExternal[1][1] ) / 2.0;
            tipNodeForce[2] = nodeDisconnectForce[2] - ( fExternal[0][2] - fExternal[1][2] ) / 2.0;

//          tipNodeForce[0] = nodeDisconnectForce[0];
//          tipNodeForce[1] = nodeDisconnectForce[1];
//          tipNodeForce[2] = nodeDisconnectForce[2];

            real64 tipArea = faceArea( trailingFaceIndex );
            if( faceToNodeMap.sizeOfArray( trailingFaceIndex ) == 3 )
            {
              tipArea *= 2.0;
            }

            tipNodeSIF = pow( (fabs( tipNodeForce[0] * trailingNodeDisp[0] / 2.0 / tipArea ) + fabs( tipNodeForce[1] * trailingNodeDisp[1] / 2.0 / tipArea )
                               + fabs( tipNodeForce[2] * trailingNodeDisp[2] / 2.0 / tipArea )), 0.5 );

            // Invert SIF if aperture is negative (confining stress).
            if( LvArray::tensorOps::AiBi< 3 >( trailingNodeDisp, faceNormalVector ) < 0.0 )
            {
              tipNodeSIF *= -1;
            }

            SIFNode_All[nodeIndex].emplace_back( tipNodeSIF );


            //Calculate SIF on tip faces connected to this trailing face and the tip node.
            for( localIndex const edgeIndex: tipEdgesID )
            {
              if( edgeToNodeMap[edgeIndex][0] == nodeIndex || edgeToNodeMap[edgeIndex][1] == nodeIndex )
              {
                real64 SIF_I = 0, SIF_II = 0, /*SIF_III,*/ SIF_Face;
                real64 vecTipNorm[3], vecTip[3], tipForce[3], tipOpening[3];

                LvArray::tensorOps::copy< 3 >( vecTipNorm, faceNormal[trailingFaceIndex] );
                LvArray::tensorOps::subtract< 3 >( vecTipNorm, faceNormal[childFaceIndices[trailingFaceIndex]] );
                LvArray::tensorOps::normalize< 3 >( vecTipNorm );

                real64 vecEdge[3];
                edgeManager.calculateLength( edgeIndex, X, vecEdge );
                LvArray::tensorOps::normalize< 3 >( vecEdge );

                LvArray::tensorOps::crossProduct( vecTip, vecTipNorm, vecEdge );
                LvArray::tensorOps::normalize< 3 >( vecTip );
                real64 v0[3];
                edgeManager.calculateCenter( edgeIndex, X, v0 );
                LvArray::tensorOps::subtract< 3 >( v0, faceCenter[ trailingFaceIndex ] );

                if( LvArray::tensorOps::AiBi< 3 >( v0, vecTip ) < 0 )
                  LvArray::tensorOps::scale< 3 >( vecTip, -1.0 );

                tipForce[0] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecTipNorm ) -
                              ( LvArray::tensorOps::AiBi< 3 >( fExternal[0], vecTipNorm ) - LvArray::tensorOps::AiBi< 3 >( fExternal[1], vecTipNorm ) ) / 2.0;
                tipForce[1] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecTip ) -
                              ( LvArray::tensorOps::AiBi< 3 >( fExternal[0], vecTip ) - LvArray::tensorOps::AiBi< 3 >( fExternal[1], vecTip ) ) / 2.0;
                tipForce[2] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecEdge ) -
                              ( LvArray::tensorOps::AiBi< 3 >( fExternal[0], vecEdge ) - LvArray::tensorOps::AiBi< 3 >( fExternal[1], vecEdge ) ) / 2.0;

//              tipForce[0] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecTipNorm );
//              tipForce[1] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecTip );
//              tipForce[2] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecEdge );

                tipOpening[0] = LvArray::tensorOps::AiBi< 3 >( trailingNodeDisp, vecTipNorm );
                tipOpening[1] = LvArray::tensorOps::AiBi< 3 >( trailingNodeDisp, vecTip );
                tipOpening[2] = LvArray::tensorOps::AiBi< 3 >( trailingNodeDisp, vecEdge );

//              if( tipForce[0] > 0.0 )
                {
                  SIF_I = pow( fabs( tipForce[0] * tipOpening[0] / 2.0 / tipArea ), 0.5 );
                  SIF_II = pow( fabs( tipForce[1] * tipOpening[1] / 2.0 / tipArea ), 0.5 );
//              SIF_III = pow( fabs( tipForce[2] * tipOpening[2] / 2.0 / tipArea ), 0.5 );
                }

                if( tipOpening[0] < 0 )
                {
                  SIF_I *= -1.0;
                }

                if( tipForce[1] < 0.0 )
                {
                  SIF_II *= -1.0;
                }

                for( localIndex const faceIndex: edgeToFaceMap[ edgeIndex ] )
                {
                  if( m_tipFaces.contains( faceIndex ))
                  {
                    real64 vecFace[ 3 ];
                    real64 fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( faceCenter[faceIndex] );

                    //Get the vector in the face and normal to the edge.
                    real64 udist;

                    real64 x0_x1[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[edgeToNodeMap[edgeIndex][0]] );
                    real64 x0_fc[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );

                    LvArray::tensorOps::subtract< 3 >( x0_x1, X[edgeToNodeMap[edgeIndex][1]] );
                    LvArray::tensorOps::normalize< 3 >( x0_x1 );
                    LvArray::tensorOps::subtract< 3 >( x0_fc, X[edgeToNodeMap[edgeIndex][1]] );
                    udist = LvArray::tensorOps::AiBi< 3 >( x0_x1, x0_fc );

                    real64 ptPrj[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( x0_x1 );
                    LvArray::tensorOps::scale< 3 >( ptPrj, udist );
                    LvArray::tensorOps::add< 3 >( ptPrj, X[edgeToNodeMap[edgeIndex][1]] );
                    LvArray::tensorOps::copy< 3 >( vecFace, fc );
                    LvArray::tensorOps::subtract< 3 >( vecFace, ptPrj );
                    LvArray::tensorOps::normalize< 3 >( vecFace );

//                  if( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace ) > cos( m_maxTurnAngle ))
                    {
                      // Clamp argument for acos.
                      real64 thetaFace = acos( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace )*0.999999 );

                      real64 tipCrossFace[ 3 ];
                      LvArray::tensorOps::crossProduct( tipCrossFace, vecTip, vecEdge );

                      if( LvArray::tensorOps::AiBi< 3 >( tipCrossFace, vecEdge ) < 0.0 )
                      {
                        thetaFace *= -1.0;
                      }

                      SIF_Face = cos( thetaFace / 2.0 ) *
                                 ( SIF_I * cos( thetaFace / 2.0 ) * cos( thetaFace / 2.0 ) - 1.5 * SIF_II * sin( thetaFace ) );

                      SIFonFace_All[faceIndex].emplace_back( SIF_Face );
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } );

  //wu40: the tip node may be included in multiple trailing faces and SIF of the node/face will be calculated multiple
  // times. We chose the smaller node SIF and the larger face SIF.
  for( localIndex const nodeIndex : m_tipNodes )
  {
    if( isNodeGhost[nodeIndex] < 0 )
    {
      if( SIFNode_All[nodeIndex].size() >= 1 )
      {
        SIFNode[nodeIndex] = *min_element( SIFNode_All[nodeIndex].begin(), SIFNode_All[nodeIndex].end());
      }

      for( localIndex const edgeIndex: m_tipEdges )
      {
        if( edgeToNodeMap[edgeIndex][0] == nodeIndex || edgeToNodeMap[edgeIndex][1] == nodeIndex )
        {
          for( localIndex const faceIndex: edgeToFaceMap[ edgeIndex ] )
          {
            if( m_tipFaces.contains( faceIndex ))
            {
              if( SIFonFace_All[faceIndex].size() >= 1 )
              {
                SIFonFace[faceIndex] = *max_element( SIFonFace_All[faceIndex].begin(), SIFonFace_All[faceIndex].end());
              }
            }
          }
        }
      }
    }
  }
}

real64 SurfaceGenerator::calculateEdgeSif( DomainPartition const & domain,
                                           localIndex const edgeID,
                                           localIndex & trailFaceID,
                                           NodeManager const & nodeManager,
                                           EdgeManager & edgeManager,
                                           FaceManager const & faceManager,
                                           ElementRegionManager const & elementManager,
                                           real64 ( & vecTipNorm )[3],
                                           real64 ( & vecTip )[3] )
{
  GEOS_MARK_FUNCTION;

  real64 rval;
  localIndex_array faceInvolved;
  arrayView1d< real64 > const & SIF_I = edgeManager.getField< surfaceGeneration::SIF_I >();
  arrayView1d< real64 > const & SIF_II = edgeManager.getField< surfaceGeneration::SIF_II >();
  arrayView1d< real64 > const & SIF_III = edgeManager.getField< surfaceGeneration::SIF_III >();

  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = nodeManager.edgeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & displacement =
    nodeManager.getField< fields::solidMechanics::totalDisplacement >();

  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList();
  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();

  arrayView1d< localIndex const > const & faceParentIndex = faceManager.getField< fields::parentIndex >();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();

  SIF_I[edgeID] = 0.0;
  SIF_II[edgeID] = 0.0;
  SIF_III[edgeID] = 0.0;

  for( localIndex const iface : edgeToFaceMap[ edgeID ] )
  {
    if( faceIsExternal[iface] >= 1 )
    {
      faceInvolved.emplace_back( iface );
    }
  }

  // Figure out the two fracture faces connected to this edge
  localIndex faceA( 0 ), faceAp( 0 );
  if( (faceParentIndex[faceInvolved[0]] == -1 && faceParentIndex[faceInvolved[1]] == faceInvolved[0]) ||
      (faceParentIndex[faceInvolved[1]] == -1 && faceParentIndex[faceInvolved[0]] == faceInvolved[1]) )
  {
    faceA = faceInvolved[0];
    faceAp = faceInvolved[1];
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Edge {} has two external faces, but the parent-child relationship is wrong.",
                          edgeID ),
                getDataContext() );
  }

  trailFaceID = faceParentIndex[faceInvolved[0]]==-1 ? faceInvolved[0] : faceParentIndex[faceInvolved[0]];


  // Define unit vectors: vecEdge (node0->node1), vecTip (opening->solid), vecTipNorm (face normal).

  LvArray::tensorOps::copy< 3 >( vecTipNorm, faceNormal[faceA] );
  LvArray::tensorOps::subtract< 3 >( vecTipNorm, faceNormal[faceAp] );
  LvArray::tensorOps::normalize< 3 >( vecTipNorm );

  //TODO: wu40: There is a function for EdgeVector in EdgeManager.cpp but has been commented.
  real64 vecEdge[3];
  edgeManager.calculateLength( edgeID, X, vecEdge );
  real64 const edgeLength = LvArray::tensorOps::l2Norm< 3 >( vecEdge );

  LvArray::tensorOps::crossProduct( vecTip, vecTipNorm, vecEdge );
  LvArray::tensorOps::normalize< 3 >( vecTip );
  real64 v0[3];
  edgeManager.calculateCenter( edgeID, X, v0 );
  LvArray::tensorOps::subtract< 3 >( v0, faceCenter[faceA] );

  if( LvArray::tensorOps::AiBi< 3 >( v0, vecTip ) < 0 )
    LvArray::tensorOps::scale< 3 >( vecTip, -1.0 );

  real64 tipCrossTipNorm[ 3 ];
  LvArray::tensorOps::crossProduct( tipCrossTipNorm, vecTip, vecTipNorm );
  if( LvArray::tensorOps::AiBi< 3 >( tipCrossTipNorm, vecEdge ) < 0 )
  {
    LvArray::tensorOps::scale< 3 >( vecTipNorm, -1 );
    faceA = faceInvolved[1];
    faceAp = faceInvolved[0];
  }


  //Now we need to figure out if a special situation applies to this edge
  // where the fracture face is a quad and three of the nodes are still pinched
  // We use a different algorithm for this special situation.

  bool threeNodesPinched( false );
  localIndex_array openNodeID;

  if( faceToNodeMap.sizeOfArray( faceA ) == 4 )  // Only quads have this problem
  {
    int numSharedNodes = 2;
    localIndex_array lNodeFaceA, lNodeFaceAp;

    lNodeFaceA.insert( 0, faceToNodeMap[ faceA ].begin(), faceToNodeMap[ faceA ].end() );
    lNodeFaceAp.insert( 0, faceToNodeMap[ faceAp ].begin(), faceToNodeMap[ faceAp ].end() );

    //We remove all the shared nodes and the one remains should be the open one.
    lNodeFaceAp.erase( std::distance( lNodeFaceAp.begin(), (std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), edgeToNodeMap[edgeID][0] ))));
    lNodeFaceAp.erase( std::distance( lNodeFaceAp.begin(), (std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), edgeToNodeMap[edgeID][1] ))));
    lNodeFaceA.erase( std::distance( lNodeFaceA.begin(), (std::find( lNodeFaceA.begin(), lNodeFaceA.end(), edgeToNodeMap[edgeID][0] ))));
    lNodeFaceA.erase( std::distance( lNodeFaceA.begin(), (std::find( lNodeFaceA.begin(), lNodeFaceA.end(), edgeToNodeMap[edgeID][1] ))));

    for( localIndex const j : faceToNodeMap[ faceA ] )
    {
      localIndex iNd = j;
      if( iNd != edgeToNodeMap[edgeID][0] && iNd != edgeToNodeMap[edgeID][1] )
      {
        auto faceToNodeMapIterator = faceToNodeMap[ faceAp ];
        if( std::find( faceToNodeMapIterator.begin(), faceToNodeMapIterator.end(), iNd ) != faceToNodeMapIterator.end())
        {
          numSharedNodes++;
          lNodeFaceA.erase( std::distance( lNodeFaceA.begin(), (std::find( lNodeFaceA.begin(), lNodeFaceA.end(), iNd ))));
          lNodeFaceAp.erase( std::distance( lNodeFaceAp.begin(), (std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), iNd ))));
        }
      }
    }

    if( numSharedNodes == 4 )
    {
      GEOS_ERROR( "The fracture face has four shared nodes with its child. This should not happen.", getDataContext() );
    }
    else if( numSharedNodes == 3 )
    {
      threeNodesPinched = true;

      //wu40: I think the following check is not necessary.
      if( lNodeFaceA.size() != 1 || lNodeFaceAp.size() != 1 )
      {
        GEOS_ERROR( "these two faces share three nodes but the number of remaining nodes is not one.", getDataContext() );
      }
      else
      {
        openNodeID.emplace_back( lNodeFaceA[0] );
        openNodeID.emplace_back( lNodeFaceAp[0] );
      }
    }
  }


  // Now we need to identify which node on the edge is the convex point and which one is the concave corner.  The convex
  // node must share an edge with the open node.
  localIndex convexCorner( std::numeric_limits< localIndex >::max());
  if( threeNodesPinched )
  {
    localIndex iNd, jNd;
    iNd = edgeToNodeMap[edgeID][0];
    jNd = edgeToNodeMap[edgeID][1];
    for( localIndex const j : faceToEdgeMap[ faceA ] )
    {
      localIndex edge = j;
      if((openNodeID[0] == edgeToNodeMap[edge][0] && iNd == edgeToNodeMap[edge][1]) ||
         (openNodeID[0] == edgeToNodeMap[edge][1] && iNd == edgeToNodeMap[edge][0])
         )
      {
        convexCorner = iNd;
        break;
      }
      if((openNodeID[0] == edgeToNodeMap[edge][0] && jNd == edgeToNodeMap[edge][1]) ||
         (openNodeID[0] == edgeToNodeMap[edge][1] && jNd == edgeToNodeMap[edge][0])
         )
      {
        convexCorner = jNd;
        break;
      }
    }

    if( convexCorner == std::numeric_limits< localIndex >::max())
      GEOS_ERROR( "This is a three-node-pinched edge but I cannot find the convex corner", getDataContext() );

  }


  // Calculate element forces acting on this edge, i.e., f_disconnect.  Need to add nodal forces from two nodes up.
  //An element has to be within the range of this edge to be included.
  //For the threeNodesPinched case, we only use the force on the node at the convex point, not the concave point.  The
  // force at the former is usually greater, so we just pick the great one instead of doing a geometrical check.
  real64 fNodeO[3] = { 0.0, 0.0, 0.0 };
  real64 GdivBeta = 0.0;  // Need this for opening-based SIF

  localIndex_array nodeIndices;

  if( !threeNodesPinched )
  {
    for( localIndex a=0; a<edgeToNodeMap.size( 1 ); ++a )
    {
      nodeIndices.emplace_back( edgeToNodeMap( edgeID, a ) );
    }
  }
  else
  {
    nodeIndices.emplace_back( convexCorner );
  }

  calculateElementForcesOnEdge ( domain, edgeID, edgeLength, nodeIndices,
                                 nodeManager, edgeManager, elementManager, vecTipNorm, fNodeO, GdivBeta, threeNodesPinched, false );


  localIndex tipFaces[2];
  tipFaces[0] = faceA;
  tipFaces[1] = faceAp;

  // Now calculate f_u. We have to subtract the nodal force at other nodes (trailing nodes) on these two open faces to
  // take into account
  // the effects of surface traction along the fracture.
  // Finding the two trailing nodes on a hex mesh is pretty straightforward, while it is cumbersome to do in tet mesh
  // For the threeNodesPinched case, this should be the open node.
  real64 fFaceA[2][3];

  // If the two external faces connected to a trailing edge are not coplanar, then we have the risk of incomplete
  // topology.
  // In that case, we use a displacement/opening based method, not VCCT.
  bool incompleteTrailingEdgeTopology = false;

  for( localIndex i=0; i<2; ++i )
  {
    localIndex_array trailingNodes;
    trailingNodes.clear();
    if( threeNodesPinched )
    {
      trailingNodes.emplace_back( openNodeID[i] );
    }
    else
    {
      localIndex faceID = tipFaces[i];
      LvArray::tensorOps::fill< 3 >( fFaceA[i], 0.0 );

      for( localIndex const j : faceToNodeMap[ faceID ] )
      {
        if( j != edgeToNodeMap( edgeID, 0 ) && j != edgeToNodeMap( edgeID, 1 ) ) // This is not a node along the tip
                                                                                 // edge
        {
          trailingNodes.emplace_back( j );
        }
      }

      if( trailingNodes.size() > 2 || trailingNodes.size() == 0 )
      {
        GEOS_ERROR( "Fatal error in finding nodes behind tip edge.", getDataContext() );
      }
      else if( trailingNodes.size() == 1 )  // Need some work to find the other node
      {
        // First find an edge that is connected to this node and parallel to the tip edge
        real64 maxCosAngle = 0.0;
        localIndex pickedTrailingEdge = std::numeric_limits< localIndex >::max();
        for( localIndex const iedge : nodeToEdgeMap[ trailingNodes[ 0 ] ] )
        {
          real64 xTrailingEdge[3];
          edgeManager.calculateCenter( iedge, X, xTrailingEdge );

          real64 udist;
          real64 x0_x1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[edgeToNodeMap[edgeID][0]] );
          real64 x0_xTrailingEdge[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( xTrailingEdge );
          LvArray::tensorOps::subtract< 3 >( x0_x1, X[edgeToNodeMap( edgeID, 1 )] );
          LvArray::tensorOps::normalize< 3 >( x0_x1 );
          LvArray::tensorOps::subtract< 3 >( x0_xTrailingEdge, X[edgeToNodeMap( edgeID, 1 )] );
          udist = LvArray::tensorOps::AiBi< 3 >( x0_x1, x0_xTrailingEdge );

          if( udist <= edgeLength && udist > 0.0 )
          {
            real64 vEdge[3];
            edgeManager.calculateLength( iedge, X, vEdge );
            LvArray::tensorOps::normalize< 3 >( vEdge );

            real64 cosEdge = std::fabs( LvArray::tensorOps::AiBi< 3 >( vEdge, vecEdge ));
            if( cosEdge > maxCosAngle )
            {
              maxCosAngle = cosEdge;
              pickedTrailingEdge = iedge;
            }
          }
        }
        if( maxCosAngle > 0.75 )
          trailingNodes.emplace_back( edgeToNodeMap[pickedTrailingEdge][0] + edgeToNodeMap[pickedTrailingEdge][1] - trailingNodes[0] );
      }
    }

    localIndex trailingEdge;
    trailingEdge = std::numeric_limits< localIndex >::max();

    if( trailingNodes.size() == 2 )
    {
      //wu40: TODO: This check is from GEOS. I think this may not be necessary. Check with Randy and PC.
      if( trailingNodes[0] != trailingNodes[1] )
      {
        for( localIndex const iedge : nodeToEdgeMap[ trailingNodes[ 0 ] ] )
        {
          if( edgeToNodeMap[iedge][0] == trailingNodes[1] || edgeToNodeMap[iedge][1] == trailingNodes[1] )
          {
            trailingEdge = iedge;
          }
        }
      }

      if( trailingEdge > edgeManager.size())
      {
        int const rank = MpiWrapper::commRank( MPI_COMM_WORLD );
        std::cout << "Cannot find trailing edge (edge=" << edgeID << ", rank=" << rank <<   "  )" << std::endl;
        return 0.0;
      }

      localIndex_array extFacesOnTrailingEdge;
      for( localIndex const iface : edgeToFaceMap[ trailingEdge ] )
      {
        if( faceIsExternal[iface] >= 1 )
          extFacesOnTrailingEdge.emplace_back( iface );
      }

      if( extFacesOnTrailingEdge.size() != 2 )
      {
        incompleteTrailingEdgeTopology = true;
      }
      else
      {
        real64 extFaceNormal[2][3];
        for( localIndex j = 0; j < 2; ++j )
        {
          LvArray::tensorOps::copy< 3 >( extFaceNormal[j], faceNormal[extFacesOnTrailingEdge[j]] );
        }

        if( std::fabs( LvArray::tensorOps::AiBi< 3 >( extFaceNormal[0], extFaceNormal[1] )) < 0.9 ) //The two faces are not coplanar.
        {
          incompleteTrailingEdgeTopology = true;
        }
      }
    }

    calculateElementForcesOnEdge ( domain, edgeID, edgeLength, trailingNodes,
                                   nodeManager, edgeManager, elementManager, vecTipNorm, fFaceA[i], GdivBeta, threeNodesPinched, true );

  }


  real64 tipForce[3];
  tipForce[0] = LvArray::tensorOps::AiBi< 3 >( fNodeO, vecTipNorm ) + LvArray::tensorOps::AiBi< 3 >( fFaceA[0], vecTipNorm ) / 2.0 - LvArray::tensorOps::AiBi< 3 >( fFaceA[1], vecTipNorm ) / 2.0;
  tipForce[1] = LvArray::tensorOps::AiBi< 3 >( fNodeO, vecTip ) + LvArray::tensorOps::AiBi< 3 >( fFaceA[0], vecTip ) / 2.0 - LvArray::tensorOps::AiBi< 3 >( fFaceA[1], vecTip ) / 2.0;
  tipForce[2] = LvArray::tensorOps::AiBi< 3 >( fNodeO, vecEdge ) + LvArray::tensorOps::AiBi< 3 >( fFaceA[0], vecEdge ) / 2.0 - LvArray::tensorOps::AiBi< 3 >( fFaceA[1], vecEdge ) / 2.0;

  real64 tipDisplacement[3], tipOpening[3], tipFaceDisplacement[2][3];

  if( !threeNodesPinched )
  {
    for( localIndex i=0; i<2; ++i )
    {
      localIndex faceID = tipFaces[i];
      LvArray::tensorOps::fill< 3 >( tipFaceDisplacement[i], 0.0 );

      for( localIndex const j : faceToNodeMap[ faceID ] )
      {
        if( j != edgeToNodeMap( edgeID, 0 ) && j != edgeToNodeMap( edgeID, 1 ))
        {
          LvArray::tensorOps::add< 3 >( tipFaceDisplacement[i], displacement[j] );
        }
      }

      LvArray::tensorOps::scale< 3 >( tipFaceDisplacement[i], 1.0 / (faceToNodeMap.sizeOfArray( faceID ) - 2) );
    }
    LvArray::tensorOps::copy< 3 >( tipDisplacement, tipFaceDisplacement[1] );
    LvArray::tensorOps::subtract< 3 >( tipDisplacement, tipFaceDisplacement[0] );
  }
  else
  {
    LvArray::tensorOps::copy< 3 >( tipDisplacement, displacement[openNodeID[1]] );
    LvArray::tensorOps::subtract< 3 >( tipDisplacement, displacement[openNodeID[0]] );
  }

  tipOpening[0] = LvArray::tensorOps::AiBi< 3 >( tipDisplacement, vecTipNorm );
  tipOpening[1] = LvArray::tensorOps::AiBi< 3 >( tipDisplacement, vecTip );
  tipOpening[2] = LvArray::tensorOps::AiBi< 3 >( tipDisplacement, vecEdge );

  real64 tipArea;
  tipArea = faceArea( faceA );
  if( faceToNodeMap.sizeOfArray( faceA ) == 3 )
  {
    tipArea *= 2.0;
  }

  if( !incompleteTrailingEdgeTopology && tipOpening[0] * tipForce[0] > 0.0 )
  {
    SIF_I[edgeID] = pow( fabs( tipForce[0] * tipOpening[0] / 2.0 / tipArea ), 0.5 );
    SIF_II[edgeID] = pow( fabs( tipForce[1] * tipOpening[1] / 2.0 / tipArea ), 0.5 );
    SIF_III[edgeID] = pow( fabs( tipForce[2] * tipOpening[2] / 2.0 / tipArea ), 0.5 );

    if( tipOpening[0] < 0 )
    {
      // We don't need this for the case of incomplete trailing edge topology.  Sign in that case should be taken care
      // of automatically because there is no sqrt involved.
      SIF_I( edgeID ) *= -1.0;
    }

  }
  else
  {
    // Opening-based SIF, based on
    // Equation 1 in Fu et al. 2012, DOI::10.1016/j.engfracmech.2012.04.010
    real64 r = tipArea / edgeLength;
    SIF_I[edgeID] = tipOpening[0] / 2.0 * GdivBeta / pow( r/6.28, 0.5 );
    SIF_II[edgeID] = 0.0;  // SIF is not accurate in this scenario anyway.  Let's not worry about turning.
    SIF_III[edgeID] = 0.0;
  }

  if( tipForce[1] < 0.0 )
  {
    SIF_II[edgeID] *= -1.0;
  }

  if( SIF_I[edgeID] > 0.0 )
  {
    rval = pow( SIF_I[edgeID]*SIF_I[edgeID]+SIF_II[edgeID]*SIF_II[edgeID]+SIF_III[edgeID]*SIF_III[edgeID], 0.5 );
  }
  else
  {
    rval = -1.0;
  }

//  std::cout << "EdgeID: " << edgeID << " SIF: " << rval << std::endl;

  return rval;
}


int SurfaceGenerator::calculateElementForcesOnEdge( DomainPartition const & domain,
                                                    localIndex const edgeID,
                                                    real64 edgeLength,
                                                    localIndex_array & nodeIndices,
                                                    NodeManager const & nodeManager,
                                                    EdgeManager const & edgeManager,
                                                    ElementRegionManager const & elementManager,
                                                    real64 ( & vecTipNorm )[3],
                                                    real64 ( & fNode )[3],
                                                    real64 & GdivBeta,
                                                    bool threeNodesPinched,
                                                    bool calculatef_u )
{
  GEOS_MARK_FUNCTION;

  ArrayOfArraysView< localIndex const > const & nodeToRegionMap = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToSubRegionMap = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementMap = nodeManager.elementList().toViewConst();

  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager.nodeList();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();
  m_solidMaterialFullIndex.resize( elementManager.numRegions() );
  elementManager.forElementRegionsComplete< CellElementRegion >( [&]( localIndex regionIndex,
                                                                      CellElementRegion const & region )
  {
    string const & solidMaterialName = region.getSubRegion( 0 ).getReference< string >( viewKeyStruct::solidMaterialNameString() );
    ConstitutiveBase const & solid = constitutiveManager.getConstitutiveRelation< ConstitutiveBase >( solidMaterialName );
    m_solidMaterialFullIndex[regionIndex] = solid.getIndexInParent();
  } );

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elementManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "shearModulus", constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elementManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "bulkModulus", constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView3d< real64 const, geos::solid::STRESS_USD > > const
  stress = elementManager.constructFullMaterialViewAccessor< array3d< real64, geos::solid::STRESS_PERMUTATION >,
                                                             arrayView3d< real64 const, geos::solid::STRESS_USD > >( fields::solid::stress::key(),
                                                                                                                     constitutiveManager );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenter =
    elementManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

  localIndex nElemEachSide[2];
  nElemEachSide[0] = 0;
  nElemEachSide[1] = 0;

  real64 xEdge[3] = { 0.0 };

  if( !calculatef_u )
  {
    edgeManager.calculateCenter( edgeID, X, xEdge );
  }

  for( localIndex i=0; i < nodeIndices.size(); ++i )
  {
    localIndex nodeID = nodeIndices( i );

    for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( nodeID ); ++k )
    {
      localIndex const er  = nodeToRegionMap[nodeID][k];
      localIndex const esr = nodeToSubRegionMap[nodeID][k];
      localIndex const ei  = nodeToElementMap[nodeID][k];

      CellElementSubRegion const & elementSubRegion = elementManager.getRegion( er ).getSubRegion< CellElementSubRegion >( esr );

      real64 xEle[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( elemCenter[er][esr][ei] );

      real64 x0_x1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[edgeToNodeMap[edgeID][0]] );
      LvArray::tensorOps::subtract< 3 >( x0_x1, X[edgeToNodeMap[edgeID][1]] );
      LvArray::tensorOps::normalize< 3 >( x0_x1 );

      real64 x0_xEle[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( xEle );
      LvArray::tensorOps::subtract< 3 >( x0_xEle, X[edgeToNodeMap[edgeID][1]] );
      real64 const udist = LvArray::tensorOps::AiBi< 3 >( x0_x1, x0_xEle );

      if(( udist <= edgeLength && udist > 0.0 ) || threeNodesPinched )
      {
        real64 const K  = bulkModulus[er][esr][m_solidMaterialFullIndex[er]][ei];
        real64 const G = shearModulus[er][esr][m_solidMaterialFullIndex[er]][ei];
        real64 const YoungModulus = 9 * K * G / ( 3 * K + G );
        real64 const poissonRatio = ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );

        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elementsToNodes = elementSubRegion.nodeList();
        for( localIndex n=0; n<elementsToNodes.size( 1 ); ++n )
        {
          if( elementsToNodes( ei, n ) == nodeID )
          {
            real64 temp[3] {0};
            LvArray::tensorOps::copy< 3 >( xEle, elemCenter[er][esr][ei] ); //For C3D6 element type, elementsToNodes map may include
            // repeated indices and the following may run multiple
            // times for the same element.

            // wu40: the nodal force need to be weighted by Young's modulus and possion's ratio.

            arrayView3d< real64 const, geos::solid::STRESS_USD > const & subregionStress = stress[er][esr][m_solidMaterialFullIndex[er]];

            // ASSUMPTION(?): subregions only have one registered element type
            finiteElement::FiniteElementBase const * finiteElement = nullptr;
            elementSubRegion.forWrappers< finiteElement::FiniteElementBase >( [&]( auto const & fe )
            {
              finiteElement = &fe.reference();
            } );

            finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( *finiteElement, [ & ]( auto const & fe )
            {

              using FE_TYPE = TYPEOFREF( fe );
              typename FE_TYPE::StackVariables feStack;
              constexpr localIndex numQuadraturePoints = FE_TYPE::numQuadraturePoints;
              constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
              real64 xLocal[numNodesPerElement][3] = {{0}};
              localIndex const numSupportPoints = finiteElement->getNumSupportPoints();
              for( localIndex a = 0; a < numSupportPoints; ++a )
              {
                LvArray::tensorOps::copy< 3 >( xLocal[a], X[elementsToNodes( ei, a )] );
              }

              for( localIndex q = 0; q < numQuadraturePoints; ++q )
              {
                real64 dNdX[numNodesPerElement][3] = {{0}};
                real64 J[3][3] = {{0}};

                real64 const detJxW = FE_TYPE::calcJacobian( q, xLocal, feStack, J );
                FE_TYPE::calcGradN( q, xLocal, feStack, dNdX );

                temp[0] -= (subregionStress( ei, q, 0 ) * dNdX[n][0] +
                            subregionStress( ei, q, 5 ) * dNdX[n][1] +
                            subregionStress( ei, q, 4 ) * dNdX[n][2]) * detJxW;
                temp[1] -= (subregionStress( ei, q, 5 ) * dNdX[n][0] +
                            subregionStress( ei, q, 1 ) * dNdX[n][1] +
                            subregionStress( ei, q, 3 ) * dNdX[n][2]) * detJxW;
                temp[2] -= (subregionStress( ei, q, 4 ) * dNdX[n][0] +
                            subregionStress( ei, q, 3 ) * dNdX[n][1] +
                            subregionStress( ei, q, 2 ) * dNdX[n][2]) * detJxW;
              }
            } );

            LvArray::tensorOps::scale< 3 >( temp, YoungModulus );
            LvArray::tensorOps::scale< 3 >( temp, 1.0 / (1 - poissonRatio * poissonRatio) );

            if( !calculatef_u )
            {
              LvArray::tensorOps::subtract< 3 >( xEle, xEdge );
              if( LvArray::tensorOps::AiBi< 3 >( xEle, vecTipNorm ) > 0 )
              {
                nElemEachSide[0] += 1;
                LvArray::tensorOps::add< 3 >( fNode, temp );

                //wu40: for debug purpose
//                std::cout << "ElementID: " << iEle << ", NodeID: " << nodeID << std::endl;
//                std::cout << "Nodal force: " << temp[0] << ", " << temp[1] << ", " << temp[2] << std::endl;
//                std::cout << "Add to total nodal force (fdisc): " << fNode[0] << ", " << fNode[1] << ", " << fNode[2]
// << std::endl;
              }
              else
              {
                nElemEachSide[1] +=1;
                LvArray::tensorOps::subtract< 3 >( fNode, temp );

                //wu40: for debug purpose
//                std::cout << "ElementID: " << iEle << ", NodeID: " << nodeID << std::endl;
//                std::cout << "Nodal force: " << temp[0] << ", " << temp[1] << ", " << temp[2] << std::endl;
//                std::cout << "Minus from total nodal force (fdisc): " << fNode[0] << ", " << fNode[1] << ", " <<
// fNode[2] << std::endl;
              }
            }
            else
            {
              LvArray::tensorOps::add< 3 >( fNode, temp );

              //wu40: for debug purpose
//              std::cout << "ElementID: " << iEle << ", NodeID: " << nodeID << std::endl;
//              std::cout << "Nodal force: " << temp[0] << ", " << temp[1] << ", " << temp[2] << std::endl;
//              std::cout << "Add to total nodal force (fext): " << fNode[0] << ", " << fNode[1] << ", " << fNode[2] <<
// std::endl;
            }
          }
        }

        if( !calculatef_u )
        {
          GdivBeta += G /2/(1-poissonRatio);
        }
      }

    }

    //If we only find one node behind the tip for the non-threeNodesPinched scenario, we do the following as a rough
    // compensation for f_u.
    if( calculatef_u && nodeIndices.size() == 1 && !threeNodesPinched )
    {
      LvArray::tensorOps::scale< 3 >( fNode, 2.0 );
    }
  }

  if( !calculatef_u )
  {
    if( nElemEachSide[0]>=1 && nElemEachSide[1]>=1 )
      LvArray::tensorOps::scale< 3 >( fNode, 0.5 );
    //We have contributions from both sides. The two sizes are the two sides of the fracture plane.  If the fracture
    // face
    // is on domain boundary, it's possible to have just one side.
    if( nElemEachSide[0] + nElemEachSide[1] >= 1 )
      GdivBeta /= (nElemEachSide[0] + nElemEachSide[1]);
  }

  return 0;
}

int SurfaceGenerator::checkOrphanElement( ElementRegionManager const & elementManager,
                                          FaceManager const & faceManager,
                                          localIndex iFace )
{
  GEOS_MARK_FUNCTION;

  arrayView2d< localIndex const > const & faceToRegionMap = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & faceToSubRegionMap = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();

  arrayView1d< integer const > const & ruptureState = faceManager.getField< surfaceGeneration::ruptureState >();

  int flagOrphan = 0;
  for( localIndex k=0; k<faceToRegionMap.size( 1 ); ++k )
  {
    localIndex const er = faceToRegionMap[iFace][k];
    localIndex const esr = faceToSubRegionMap[iFace][k];
    localIndex const ei = faceToElementMap[iFace][k];
    if( er != -1 &&  esr != -1 && ei != -1 )
    {
      CellElementSubRegion const & elementSubRegion = elementManager.getRegion( faceToRegionMap[iFace][k] ).
                                                        getSubRegion< CellElementSubRegion >( faceToSubRegionMap[iFace][k] );


      int nRuptureFace = 0;
      arrayView2d< localIndex const > const & elementsToFaces = elementSubRegion.faceList();
      for( localIndex a=0; a < elementsToFaces.size( 1 ); ++a )
      {
        localIndex jFace = elementsToFaces[ei][a];
        if( (ruptureState[jFace] == 1 || faceIsExternal[jFace] >= 1) && jFace != iFace )
        {
          nRuptureFace +=1;
        }
      }

      if( nRuptureFace == elementsToFaces.size( 1 ) - 1 )
      {
        flagOrphan = 1;
      }
    }
  }
  return flagOrphan;

}

void SurfaceGenerator::markRuptureFaceFromNode( localIndex const nodeIndex,
                                                NodeManager const & nodeManager,
                                                EdgeManager const & edgeManager,
                                                FaceManager & faceManager,
                                                ElementRegionManager const & GEOS_UNUSED_PARAM( elementManager ),
                                                ModifiedObjectLists & modifiedObjects )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< integer > const & ruptureState = faceManager.getField< surfaceGeneration::ruptureState >();
  arrayView1d< real64 const > const & SIFonFace = faceManager.getField< surfaceGeneration::SIFonFace >();
  arrayView2d< real64 const > const & KIC = faceManager.getField< surfaceGeneration::K_IC >();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = nodeManager.faceList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  localIndex_array eligibleFaces;
  real64_array faceSIFToToughnessRatio;
  real64 lowestSIF = std::numeric_limits< real64 >::max();
  real64 highestSIF = std::numeric_limits< real64 >::min();

  for( localIndex const faceIndex : nodeToFaceMap[ nodeIndex ] )
  {
    if( m_tipFaces.contains( faceIndex ))
    {
      real64 faceToughness;
      real64 fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter[faceIndex] );

      eligibleFaces.emplace_back( faceIndex );

      for( localIndex const edgeIndex : faceToEdgeMap[ faceIndex ] )
      {
        if( m_tipEdges.contains( edgeIndex ))
        {
          real64 direction[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
          real64 edgeCenter[3];
          edgeManager.calculateCenter( edgeIndex, X, edgeCenter );
          LvArray::tensorOps::subtract< 3 >( direction, edgeCenter );
          LvArray::tensorOps::normalize< 3 >( direction );
          faceToughness = std::fabs( LvArray::tensorOps::AiBi< 3 >( direction, KIC[faceIndex] ));

          faceSIFToToughnessRatio.emplace_back( SIFonFace[faceIndex]/faceToughness );
          highestSIF = std::max( highestSIF, SIFonFace[faceIndex]/faceToughness );
          lowestSIF = std::min( lowestSIF, SIFonFace[faceIndex]/faceToughness );
        }
      }
    }
  }

  for( localIndex i = 0; i < eligibleFaces.size(); ++i )
  {
    localIndex pickedFace = eligibleFaces[i];

    if( highestSIF > 1.0 &&
        ((eligibleFaces.size() < 3) || (eligibleFaces.size() >= 3 && (highestSIF - faceSIFToToughnessRatio[i]) <= 0.2 * (highestSIF - lowestSIF))))
    {
      ruptureState[pickedFace] = 1;
      modifiedObjects.modifiedFaces.insert( pickedFace );
    }
  }
}

void SurfaceGenerator::markRuptureFaceFromEdge( localIndex const edgeID,
                                                localIndex const & GEOS_UNUSED_PARAM( trailFaceID ),
                                                NodeManager const & nodeManager,
                                                EdgeManager const & edgeManager,
                                                FaceManager & faceManager,
                                                ElementRegionManager const & elementManager,
                                                real64 ( &GEOS_UNUSED_PARAM( vecTipNorm ) )[3],
                                                real64 ( & vecTip )[3],
                                                ModifiedObjectLists & modifiedObjects,
                                                int const edgeMode )
{
  GEOS_MARK_FUNCTION;

  arrayView1d< integer > const & ruptureState = faceManager.getField< surfaceGeneration::ruptureState >();
  arrayView1d< real64 > const & SIFonFace = faceManager.getField< surfaceGeneration::SIFonFace >();
  arrayView2d< real64 const > const & KIC = faceManager.getField< surfaceGeneration::K_IC >();
  arrayView1d< real64 const > const & SIF_I = edgeManager.getField< surfaceGeneration::SIF_I >();
  arrayView1d< real64 const > const & SIF_II = edgeManager.getField< surfaceGeneration::SIF_II >();
  arrayView1d< localIndex > const & primaryCandidateFace = faceManager.getField< surfaceGeneration::primaryCandidateFace >();
  arrayView1d< integer const > const & isFaceSeparable = faceManager.getField< surfaceGeneration::isFaceSeparable >();
//  integer_array* dfnIndexMap = faceManager.getReferencePointer<integer_array>( "DFN_Index" );

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();


  arrayView2d< localIndex > const & edgeToNodeMap = edgeManager.nodeList();
  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

  localIndex_array eligibleFaces;
  real64 lowestSIF = std::numeric_limits< real64 >::max();
  real64 highestSIF = std::numeric_limits< real64 >::min();
  real64 lowestScore = std::numeric_limits< real64 >::max();
  real64 highestScore = std::numeric_limits< real64 >::min();
  real64 secondScore = std::numeric_limits< real64 >::min();
  localIndex faceWithHighestScore = std::numeric_limits< localIndex >::max();
  localIndex faceWithSecondScore = std::numeric_limits< localIndex >::max();

  real64 vecEdge[3], edgeCenter[3];
  edgeManager.calculateLength( edgeID, X, vecEdge );
  edgeManager.calculateCenter( edgeID, X, edgeCenter );


  for( localIndex const iface : edgeToFaceMap[ edgeID ] )
  {
    if( faceToElementMap.size( 1 ) == 2  &&
        faceIsExternal[iface] < 1 &&
        checkOrphanElement( elementManager, faceManager, iface ) == 0 &&
        isFaceSeparable[iface] == 1 )
    {
      real64 const fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter[iface] );

      //Get the vector in the face and normal to the edge.
      //wu40: there is a function in GEOS for this calculation. Maybe it's worth to have a function in GEOSX too.
      real64 x0_x1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[edgeToNodeMap[edgeID][0]] );
      LvArray::tensorOps::subtract< 3 >( x0_x1, X[edgeToNodeMap[edgeID][1]] );
      LvArray::tensorOps::normalize< 3 >( x0_x1 );

      real64 x0_fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
      LvArray::tensorOps::subtract< 3 >( x0_fc, X[edgeToNodeMap[edgeID][1]] );
      real64 const udist = LvArray::tensorOps::AiBi< 3 >( x0_x1, x0_fc );

      real64 ptPrj[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( x0_x1 );
      LvArray::tensorOps::scale< 3 >( ptPrj, udist );
      LvArray::tensorOps::add< 3 >( ptPrj, X[edgeToNodeMap[edgeID][1]] );

      real64 vecFace[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
      LvArray::tensorOps::subtract< 3 >( vecFace, ptPrj );
      LvArray::tensorOps::normalize< 3 >( vecFace );

//      if( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace ) > cos( m_maxTurnAngle ))
      {
        eligibleFaces.emplace_back( iface );
        real64 thetaFace = acos( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace )*0.999999 );  // We multiply this by 0.9999999 to avoid an
        // exception caused by acos a number slightly larger
        // than 1.

        real64 tipCrossFace[ 3 ];
        LvArray::tensorOps::crossProduct( tipCrossFace, vecTip, vecFace );
        if( LvArray::tensorOps::AiBi< 3 >( tipCrossFace, vecEdge ) < 0.0 )
        {
          thetaFace *= -1.0;
        }

        SIFonFace[iface] = cos( thetaFace / 2.0 ) *
                           ( SIF_I[edgeID] * cos( thetaFace / 2.0 ) * cos( thetaFace / 2.0 ) - 1.5 * SIF_II[edgeID] * sin( thetaFace ) );

        real64 direction[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
        LvArray::tensorOps::subtract< 3 >( direction, edgeCenter );
        LvArray::tensorOps::normalize< 3 >( direction );
        real64 faceToughness = std::fabs( LvArray::tensorOps::AiBi< 3 >( direction, KIC[iface] ));

        highestSIF = std::max( highestSIF, SIFonFace[iface]/faceToughness );
        lowestSIF = std::min( lowestSIF, SIFonFace[iface]/faceToughness );

      }
    }
  }

  localIndex_array pickedFaces;
  if( eligibleFaces.size() >=1 )
  {
    for( localIndex i = 0; i < eligibleFaces.size(); ++i )
    {
      localIndex iface = eligibleFaces[i];
      real64 direction[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( edgeCenter );
      LvArray::tensorOps::subtract< 3 >( direction, faceCenter[iface] );
      LvArray::tensorOps::normalize< 3 >( direction );
      real64 faceToughness = std::fabs( LvArray::tensorOps::AiBi< 3 >( direction, KIC[iface] ));

      real64 splitabilityScore = SIFonFace[iface] - lowestSIF * faceToughness;
      lowestScore = std::min( lowestScore, splitabilityScore );

      if( faceWithHighestScore == std::numeric_limits< localIndex >::max())
      {
        faceWithHighestScore = iface;
        highestScore = splitabilityScore;
      }
      else if( splitabilityScore > highestScore )
      {
        faceWithSecondScore = faceWithHighestScore;
        secondScore = highestScore;
        faceWithHighestScore = iface;
        highestScore = splitabilityScore;
      }
      else if( splitabilityScore > secondScore )
      {
        faceWithSecondScore = iface;
        secondScore = splitabilityScore;
      }
    }

    pickedFaces.emplace_back( faceWithHighestScore );

    if( eligibleFaces.size() >= 3 && (highestScore - secondScore) < 0.1 * (highestScore - lowestScore))
    {
      pickedFaces.emplace_back( faceWithSecondScore );
    }

  }


  for( localIndex i = 0; i < pickedFaces.size(); ++i )
  {
    localIndex pickedFace = pickedFaces[i];

    if( highestSIF > 1.0 && edgeMode == 1 && i == 0 && isFaceSeparable[pickedFace] == 1 )
    {
      ruptureState[pickedFace] = 1;
      modifiedObjects.modifiedFaces.insert( pickedFace );
    }
    else if( highestSIF > 1.0 && edgeMode == 1 && i == 1 && isFaceSeparable[pickedFace] == 1 )
    {
      ruptureState[pickedFace] = -1;
      modifiedObjects.modifiedFaces.insert( pickedFace );
      primaryCandidateFace[pickedFace] = faceWithHighestScore;
    }
  }
}

void SurfaceGenerator::postUpdateRuptureStates( NodeManager const & nodeManager,
                                                EdgeManager const & edgeManager,
                                                FaceManager const & faceManager,
                                                ElementRegionManager const & GEOS_UNUSED_PARAM( elementManager ),
                                                stdVector< std::set< localIndex > > & nodesToRupturedFaces,
                                                stdVector< std::set< localIndex > > & edgesToRupturedFaces )
{
  GEOS_MARK_FUNCTION;

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  nodesToRupturedFaces.resize( nodeManager.size() );
  edgesToRupturedFaces.resize( edgeManager.size() );

  arrayView1d< integer const > const & faceRuptureState = faceManager.getField< surfaceGeneration::ruptureState >();
  arrayView1d< localIndex const > const & faceParentIndex = faceManager.getField< fields::parentIndex >();

  // assign the values of the nodeToRupturedFaces and edgeToRupturedFaces arrays.
  for( localIndex kf=0; kf<faceManager.size(); ++kf )
  {
    if( faceRuptureState[kf] >0 )
    {
      int const n = faceParentIndex[kf]==-1 ? 1 : 2;
      localIndex const faceIndex = faceParentIndex[kf]==-1 ? kf : faceParentIndex[kf];

      for( int i=0; i<n; ++i )
      {
        for( localIndex a=0; a<faceToNodeMap.sizeOfArray( kf ); ++a )
        {
          const localIndex nodeIndex = faceToNodeMap( kf, a );
          nodesToRupturedFaces[nodeIndex].insert( faceIndex );
        }

        for( localIndex a=0; a<faceToEdgeMap.sizeOfArray( kf ); ++a )
        {
          const localIndex edgeIndex = faceToEdgeMap( kf, a );
          edgesToRupturedFaces[edgeIndex].insert( faceIndex );
        }
      }
    }
  }
}

stdVector< std::set< localIndex > > SurfaceGenerator::groupRupturedFacesIntoSets(
  stdVector< std::set< localIndex > > const & nodesToRupturedFaces,
  NodeManager const & nodeManager,
  FaceManager const & faceManager ) const
{
  GEOS_MARK_FUNCTION;

  // Collect all unique ruptured face indices.
  std::set< localIndex > allFaces;
  for( localIndex ni = 0; ni < static_cast< localIndex >( nodesToRupturedFaces.size() ); ++ni )
  {
    for( localIndex const fi : nodesToRupturedFaces[ni] )
    {
      allFaces.insert( fi );
    }
  }

  if( allFaces.empty() )
  {
    return {};
  }

  // Union-Find data structure (map-based so we don't need a dense array).
  // std::map is used instead of std::unordered_map to guarantee deterministic
  // iteration order across platforms (libstdc++ vs libc++) and MPI runs.
  std::map< localIndex, localIndex > parent;
  std::map< localIndex, localIndex > rank;

  for( localIndex const fi : allFaces )
  {
    parent[fi] = fi;
    rank[fi] = 0;
  }

  // Find with path compression.
  std::function< localIndex( localIndex ) > findRoot = [&]( localIndex x ) -> localIndex
  {
    if( parent[x] != x )
    {
      parent[x] = findRoot( parent[x] );
    }
    return parent[x];
  };

  // Union by rank.
  auto unionSets = [&]( localIndex a, localIndex b )
  {
    localIndex ra = findRoot( a );
    localIndex rb = findRoot( b );
    if( ra == rb )
      return;
    if( rank[ra] < rank[rb] )
      std::swap( ra, rb );
    parent[rb] = ra;
    if( rank[ra] == rank[rb] )
      ++rank[ra];
  };

  // Two ruptured faces sharing a node belong to the same fracture set.
  // Iterate through nodesToRupturedFaces: for each node that touches
  // more than one ruptured face, union all those faces together.
  GEOS_UNUSED_VAR( nodeManager );
  for( localIndex ni = 0; ni < static_cast< localIndex >( nodesToRupturedFaces.size() ); ++ni )
  {
    std::set< localIndex > const & faces = nodesToRupturedFaces[ni];
    if( faces.size() <= 1 )
      continue;
    auto it = faces.begin();
    localIndex const first = *it;
    ++it;
    for(; it != faces.end(); ++it )
    {
      unionSets( first, *it );
    }
  }

  // Gather connected components.
  // std::map guarantees deterministic iteration order (by sorted key),
  // so the result vector is built in the same order on every platform/rank.
  std::map< localIndex, std::set< localIndex > > components;
  for( localIndex const fi : allFaces )
  {
    components[findRoot( fi )].insert( fi );
  }

  // Convert to a vector.
  stdVector< std::set< localIndex > > result;
  result.reserve( components.size() );
  for( auto & kv : components )
  {
    result.push_back( std::move( kv.second ) );
  }

  // Sort by the minimum global face index in each set.
  // This guarantees that every MPI rank processes the same fracture
  // at the same loop iteration, which is required because
  // synchronizeTopologyChange is called collectively after each set.
  GEOS_UNUSED_VAR( nodeManager );
  arrayView1d< globalIndex const > const & faceLocalToGlobal = faceManager.localToGlobalMap();
  std::sort( result.begin(), result.end(),
             [&faceLocalToGlobal]( std::set< localIndex > const & a,
                                   std::set< localIndex > const & b )
  {
    // Compute minimum global face index for each set.
    globalIndex minA = std::numeric_limits< globalIndex >::max();
    for( localIndex const fi : a )
    {
      minA = std::min( minA, faceLocalToGlobal[fi] );
    }
    globalIndex minB = std::numeric_limits< globalIndex >::max();
    for( localIndex const fi : b )
    {
      minB = std::min( minB, faceLocalToGlobal[fi] );
    }
    return minA < minB;
  } );

  return result;
}

void SurfaceGenerator::buildFilteredRuptureMaps(
  std::set< localIndex > const & fractureSetFaces,
  NodeManager const & nodeManager,
  EdgeManager const & edgeManager,
  FaceManager const & faceManager,
  stdVector< std::set< localIndex > > & nodesToRupturedFaces,
  stdVector< std::set< localIndex > > & edgesToRupturedFaces ) const
{
  GEOS_MARK_FUNCTION;

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

  nodesToRupturedFaces.clear();
  edgesToRupturedFaces.clear();
  nodesToRupturedFaces.resize( nodeManager.size() );
  edgesToRupturedFaces.resize( edgeManager.size() );

  arrayView1d< integer const > const & faceRuptureState = faceManager.getField< surfaceGeneration::ruptureState >();
  arrayView1d< localIndex const > const & faceParentIndex = faceManager.getField< fields::parentIndex >();

  // Iterate over all faces.  Only include a face if:
  //   1) its rupture state > 0  (same criterion as postUpdateRuptureStates), AND
  //   2) its parent face index belongs to the current fracture set.
  for( localIndex kf = 0; kf < faceManager.size(); ++kf )
  {
    if( faceRuptureState[kf] <= 0 )
      continue;

    localIndex const parentFace = ( faceParentIndex[kf] == -1 ) ? kf : faceParentIndex[kf];

    // Skip faces that don't belong to this fracture set.
    if( fractureSetFaces.count( parentFace ) == 0 )
      continue;

    int const n = ( faceParentIndex[kf] == -1 ) ? 1 : 2;
    for( int i = 0; i < n; ++i )
    {
      for( localIndex a = 0; a < faceToNodeMap.sizeOfArray( kf ); ++a )
      {
        localIndex const nodeIndex = faceToNodeMap( kf, a );
        nodesToRupturedFaces[nodeIndex].insert( parentFace );
      }
      for( localIndex a = 0; a < faceToEdgeMap.sizeOfArray( kf ); ++a )
      {
        localIndex const edgeIndex = faceToEdgeMap( kf, a );
        edgesToRupturedFaces[edgeIndex].insert( parentFace );
      }
    }
  }
}

int SurfaceGenerator::checkEdgeSplitability( localIndex const edgeID,
                                             NodeManager const & GEOS_UNUSED_PARAM( nodeManager ),
                                             FaceManager const & faceManager,
                                             EdgeManager const & edgeManager,
                                             bool const GEOS_UNUSED_PARAM( prefrac ) )
{
  GEOS_MARK_FUNCTION;

  //     Return value = -1, this edge won't split for sure, don't do any more work;
  //                  = 0, edge is along a tip, but the fracture connected to it is not saturated yet.  We will only
  // calculate SIF but will not perform splitting.
  //                  = 1, edge is along a tip and the adjacent fracture is saturated, more work to be done; or this is
  // a dry simulation
  //                  = 2, this is a singular edge, we need split it.
  //                  = 3, this is an eligible kink, we need to process it as a kink

  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();
  arrayView1d< localIndex const > const & faceParentIndex = faceManager.getField< fields::parentIndex >();

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();
  arrayView1d< integer const > const & edgeIsExternal = edgeManager.isExternal();

  int isSplitable = -1;

  if( edgeIsExternal[edgeID] == 0 )
  {
    return isSplitable;
  }

  // We first count the external faces connected to this edge;
  int nExternalFaces = 0;
  localIndex_array faceInvolved;
  for( localIndex const iface : edgeToFaceMap[ edgeID ] )
  {
    if( faceIsExternal[iface] >= 1 )
    {
      nExternalFaces++;
      faceInvolved.emplace_back( iface );
    }
  }

  if( nExternalFaces%2 == 1 )
  {
    //    char msg[200];
    //    sprintf(msg, "Error! Edge %d has an odd number of external faces.", int(edgeID));
    //    GEOS_ERROR(msg);
    //    std::cout << "Error! Edge " << int(edgeID) << " has an odd number of external faces. "
    //        << (*nodeManager.m_refposition)[edgeToNodeMap[edgeID][0]][0] << " ,"
    //        << (*nodeManager.m_refposition)[edgeToNodeMap[edgeID][0]][1] << " ,"
    //        << (*nodeManager.m_refposition)[edgeToNodeMap[edgeID][0]][2] << " ,";
    //    isSplitable = -1;
    return(isSplitable);
  }

  if( nExternalFaces >= 4 )
  {
    isSplitable = 2;
    return (isSplitable);
  }

  if( nExternalFaces == 2 )
  {
    localIndex parentFace = -1;
    if( faceParentIndex[faceInvolved[0]] == -1 && faceParentIndex[faceInvolved[1]] == faceInvolved[0] )
    {
      parentFace = faceInvolved[0];
    }
    else if( faceParentIndex[faceInvolved[1]] == -1 && faceParentIndex[faceInvolved[0]] == faceInvolved[1] )
    {
      parentFace = faceInvolved[1];
    }

    if( parentFace == -1 )
    {
      isSplitable = -1;
    }
    else
    {
      isSplitable = 1;
    }

  }

  return (isSplitable);
}

real64 SurfaceGenerator::minimumToughnessOnEdge( localIndex const edgeID,
                                                 NodeManager const & nodeManager,
                                                 EdgeManager const & edgeManager,
                                                 FaceManager const & faceManager )
{
  GEOS_MARK_FUNCTION;

  real64 val = std::numeric_limits< real64 >::max();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  real64 edgeCenter[3] = { 0.0 };
  edgeManager.calculateCenter( edgeID, X, edgeCenter );

  arrayView2d< real64 const > const & KIC = faceManager.getField< surfaceGeneration::K_IC >();
  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();
  for( localIndex const iface : edgeManager.faceList()[ edgeID ] )
  {
    localIndex const numFaceNodes = faceToNodes.sizeOfArray( iface );
    real64 faceCenter[3] = { 0.0 };
    for( localIndex a=0; a<numFaceNodes; ++a )
    {
      LvArray::tensorOps::add< 3 >( faceCenter, X[ faceToNodes[iface][a] ] );
    }
    LvArray::tensorOps::scale< 3 >( faceCenter, 1.0 / numFaceNodes );

    real64 direction[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter );
    LvArray::tensorOps::subtract< 3 >( direction, edgeCenter );
    LvArray::tensorOps::normalize< 3 >( direction );
    val = std::min( val, fabs( LvArray::tensorOps::AiBi< 3 >( KIC[iface], direction ) ) );
  }

  return val;
}

real64 SurfaceGenerator::minimumToughnessOnNode( localIndex const nodeID,
                                                 NodeManager const & nodeManager,
                                                 EdgeManager const & edgeManager,
                                                 FaceManager const & faceManager )
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = nodeManager.faceList().toViewConst();

  arrayView2d< real64 const > const & KIC = faceManager.getField< surfaceGeneration::K_IC >();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

  real64 val = std::numeric_limits< real64 >::max();

  for( localIndex const faceIndex: nodeToFaceMap[ nodeID ] )
  {
    if( m_tipFaces.contains( faceIndex ))
    {
      real64 const fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter[faceIndex] );

      for( localIndex const edgeIndex: faceToEdgeMap[ faceIndex ] )
      {
        if( m_tipEdges.contains( edgeIndex ))
        {
          real64 direction[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
          real64 edgeCenter[ 3 ];
          edgeManager.calculateCenter( edgeIndex, X, edgeCenter );
          LvArray::tensorOps::subtract< 3 >( direction, edgeCenter );
          LvArray::tensorOps::normalize< 3 >( direction );

          val = std::min( val, fabs( LvArray::tensorOps::AiBi< 3 >( KIC[faceIndex], direction ) ) );
        }
      }
    }
  }

  return val;
}

void SurfaceGenerator::assignNewGlobalIndicesSerial( ObjectManagerBase & object,
                                                     std::set< localIndex > const & indexList )
{
  GEOS_MARK_FUNCTION;
  // in serial, we can simply loop over the indexList and assign new global indices based on
  // the value of the maxGlobalIndex() + 1;
  arrayView1d< globalIndex > const & localToGlobal = object.localToGlobalMap();
  for( localIndex const newLocalIndex : indexList )
  {
    globalIndex const newGlobalIndex = object.maxGlobalIndex() + 1;

    localToGlobal[newLocalIndex] = newGlobalIndex;
    object.updateGlobalToLocalMap( newLocalIndex );
  }

  object.setMaxGlobalIndex();
}

void SurfaceGenerator::assignNewGlobalIndicesSerial( ElementRegionManager & elementManager,
                                                     map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & newElems )
{
  GEOS_MARK_FUNCTION;

  // in serial, we can simply iterate over the entries in newElems and assign new global indices based on
  // the value of the maxGlobalIndex() + 1 for the ElementRegionManager.

  // loop over entries of newElems, which gives elementRegion/subRegion local indices
  for( auto const & iter: newElems )
  {
    localIndex const er = iter.first.first;
    localIndex const esr = iter.first.second;
    std::set< localIndex > const & indexList = iter.second;

    ElementSubRegionBase & subRegion = elementManager.getRegion( er ).getSubRegion( esr );
    arrayView1d< globalIndex > const & localToGlobal = subRegion.localToGlobalMap();

    // loop over the new elems in the subRegion
    for( localIndex const newLocalIndex : indexList )
    {
      globalIndex const newGlobalIndex = elementManager.maxGlobalIndex() + 1;

      localToGlobal[newLocalIndex] = newGlobalIndex;
      subRegion.updateGlobalToLocalMap( newLocalIndex );
    }
  }

  elementManager.setMaxGlobalIndex();
}

real64
SurfaceGenerator::calculateRuptureRate( SurfaceElementRegion & faceElementRegion )
{
  GEOS_MARK_FUNCTION;

  real64 maxRuptureRate = 0;
  FaceElementSubRegion & subRegion = faceElementRegion.getSubRegion< FaceElementSubRegion >( 0 );

  ArrayOfArraysView< localIndex const > const &
  fractureConnectorEdgesToFaceElements = subRegion.m_2dFaceTo2dElems.toViewConst();

  arrayView1d< real64 > const & ruptureTime = subRegion.getField< fields::ruptureTime >();
  arrayView1d< real64 > const & ruptureRate = subRegion.getField< surfaceGeneration::ruptureRate >();

  arrayView2d< real64 const > const & elemCenter = subRegion.getElementCenter();

  for( localIndex kfc=0; kfc<fractureConnectorEdgesToFaceElements.size(); ++kfc )
  {
    for( localIndex kfe0=0; kfe0<fractureConnectorEdgesToFaceElements.sizeOfArray( kfc ); ++kfe0 )
    {
      localIndex const faceElem0 = fractureConnectorEdgesToFaceElements( kfc, kfe0 );
      for( localIndex kfe1=kfe0+1; kfe1<fractureConnectorEdgesToFaceElements.sizeOfArray( kfc ); ++kfe1 )
      {
        localIndex const faceElem1 = fractureConnectorEdgesToFaceElements( kfc, kfe1 );
        if( !( m_faceElemsRupturedThisSolve.count( faceElem0 ) && m_faceElemsRupturedThisSolve.count( faceElem1 ) ) )
        {
          real64 const deltaRuptureTime = fabs( ruptureTime( faceElem0 ) - ruptureTime( faceElem1 ));
          if( deltaRuptureTime > 1.0e-14 * (ruptureTime( faceElem0 ) + ruptureTime( faceElem1 )) )
          {
            real64 distance[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( elemCenter[ faceElem0 ] );
            LvArray::tensorOps::subtract< 3 >( distance, elemCenter[ faceElem1 ] );
            real64 const pairwiseRuptureRate = 1.0 / deltaRuptureTime;
            if( m_faceElemsRupturedThisSolve.count( faceElem0 ) )
            {
              ruptureRate( faceElem0 ) = std::min( ruptureRate( faceElem0 ), pairwiseRuptureRate );
            }
            else if( m_faceElemsRupturedThisSolve.count( faceElem1 ) )
            {
              ruptureRate( faceElem1 ) = std::min( ruptureRate( faceElem1 ), pairwiseRuptureRate );
            }
          }
        }
      }
    }
  }


  for( localIndex faceElemIndex : m_faceElemsRupturedThisSolve )
  {
    maxRuptureRate = std::max( maxRuptureRate, ruptureRate( faceElemIndex ) );
  }

  real64 const globalMaxRuptureRate = MpiWrapper::allReduce( maxRuptureRate,
                                                             MpiWrapper::Reduction::Max,
                                                             MPI_COMM_GEOS );

  return globalMaxRuptureRate;
}

real64 SurfaceGenerator::scalingToughness( R1Tensor const fractureOrigin,
                                           real64 const (&faceCenter)[3],
                                           real64 const initialRockToughness,
                                           real64 const toughnessScalingFactor )
{
  GEOS_MARK_FUNCTION;

  real64 const distance = sqrt( (fractureOrigin[0] - faceCenter[0])*(fractureOrigin[0] - faceCenter[0]) +
                                (fractureOrigin[1] - faceCenter[1])*(fractureOrigin[1] - faceCenter[1]) +
                                (fractureOrigin[2] - faceCenter[2])*(fractureOrigin[2] - faceCenter[2]) );

  real64 scaledToughness = initialRockToughness*( 1 + toughnessScalingFactor*sqrt( distance ) );

  return scaledToughness;
}



REGISTER_CATALOG_ENTRY( PhysicsSolverBase,
                        SurfaceGenerator,
                        string const &, dataRepository::Group * const )

} /* namespace geos */
