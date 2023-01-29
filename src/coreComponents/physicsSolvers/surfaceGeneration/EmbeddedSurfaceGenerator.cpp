/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EmbeddedSurfaceGenerator.cpp
 */

#include "EmbeddedSurfaceGenerator.hpp"
#include "EmbeddedSurfacesParallelSynchronization.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshFields.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/utilities/CIcomputationKernel.hpp"
#include "physicsSolvers/solidMechanics/kernels/SolidMechanicsLagrangianFEMKernels.hpp"
#include "mesh/simpleGeometricObjects/GeometricObjectManager.hpp"
#include "mesh/simpleGeometricObjects/BoundedPlane.hpp"



namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

void NewObjectLists::insert( NewObjectLists const & newObjects )
{
  newNodes.insert( newObjects.newNodes.begin(),
                   newObjects.newNodes.end() );

  newEdges.insert( newObjects.newEdges.begin(),
                   newObjects.newEdges.end() );

  for( auto & iter : newObjects.newElements )
  {
    std::pair< localIndex, localIndex > const & key = iter.first;
    std::set< localIndex > const & values = iter.second;
    newElements[key].insert( values.begin(), values.end() );
  }

}


EmbeddedSurfaceGenerator::EmbeddedSurfaceGenerator( const string & name,
                                                    Group * const parent ):
  SolverBase( name, parent ),
  m_fractureRegionName(),
  m_mpiCommOrder( 0 )
{
  registerWrapper( viewKeyStruct::fractureRegionNameString(), &m_fractureRegionName ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( "FractureRegion" );

  // this->getWrapper< string >( viewKeyStruct::discretizationString() ).
  // setInputFlag( InputFlags::FALSE );

  registerWrapper( viewKeyStruct::mpiCommOrderString(), &m_mpiCommOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to enable MPI consistent communication ordering" );

}

EmbeddedSurfaceGenerator::~EmbeddedSurfaceGenerator()
{}

void EmbeddedSurfaceGenerator::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getBaseDiscretization();

    EmbeddedSurfaceNodeManager & nodeManager = meshLevel.getEmbSurfNodeManager();

    nodeManager.registerField< fields::parentEdgeIndex >( this->getName() );
  } );
}


void EmbeddedSurfaceGenerator::initializePostSubGroups()
{
  /*
   * Here we generate embedded elements for fractures (or faults) that already exist in the domain and
   * were specified in the input file.
   */

  // Get domain
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // Get geometric object manager
  GeometricObjectManager & geometricObjManager = GeometricObjectManager::getInstance();

  // Get meshLevel
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & meshLevel,
                                                                arrayView1d< string const > const & )
  {

    // Get managers
    ElementRegionManager & elemManager = meshLevel.getElemManager();
    NodeManager & nodeManager = meshLevel.getNodeManager();
    EmbeddedSurfaceNodeManager & embSurfNodeManager = meshLevel.getEmbSurfNodeManager();
    EdgeManager & edgeManager = meshLevel.getEdgeManager();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord = nodeManager.referencePosition();

    // Get EmbeddedSurfaceSubRegions
    SurfaceElementRegion & embeddedSurfaceRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );
    EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion = embeddedSurfaceRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

    localIndex localNumberOfSurfaceElems         = 0;

    NewObjectLists newObjects;

    // Loop over all the fracture planes
    geometricObjManager.forSubGroups< BoundedPlane >( [&]( BoundedPlane & fracture )
    {
      /* 1. Find out if an element is cut by the fracture or not.
       * Loop over all the elements and for each one of them loop over the nodes and compute the
       * dot product between the distance between the plane center and the node and the normal
       * vector defining the plane. If two scalar products have different signs the plane cuts the
       * cell. If a nodes gives a 0 dot product it has to be neglected or the method won't work.
       */
      real64 const planeCenter[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fracture.getCenter() );
      real64 const normalVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fracture.getNormal() );
      // Initialize variables
      globalIndex nodeIndex;
      integer isPositive, isNegative;
      real64 distVec[ 3 ];

      elemManager.forElementSubRegionsComplete< CellElementSubRegion >(
        [&]( localIndex const er, localIndex const esr, ElementRegionBase &, CellElementSubRegion & subRegion )
      {
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const cellToNodes = subRegion.nodeList();
        FixedOneToManyRelation const & cellToEdges = subRegion.edgeList();

        arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

        forAll< serialPolicy >( subRegion.size(), [ &, ghostRank,
                                                    cellToNodes,
                                                    nodesCoord ] ( localIndex const cellIndex )
        {
          if( ghostRank[cellIndex] < 0 )
          {
            isPositive = 0;
            isNegative = 0;
            for( localIndex kn = 0; kn < subRegion.numNodesPerElement(); kn++ )
            {
              nodeIndex = cellToNodes[cellIndex][kn];
              LvArray::tensorOps::copy< 3 >( distVec, nodesCoord[nodeIndex] );
              LvArray::tensorOps::subtract< 3 >( distVec, planeCenter );
              // check if the dot product is zero
              if( LvArray::tensorOps::AiBi< 3 >( distVec, normalVector ) > 0 )
              {
                isPositive = 1;
              }
              else if( LvArray::tensorOps::AiBi< 3 >( distVec, normalVector ) < 0 )
              {
                isNegative = 1;
              }
            }   // end loop over nodes
            if( isPositive * isNegative == 1 )
            {
              bool added = embeddedSurfaceSubRegion.addNewEmbeddedSurface( cellIndex,
                                                                           er,
                                                                           esr,
                                                                           nodeManager,
                                                                           embSurfNodeManager,
                                                                           edgeManager,
                                                                           cellToEdges,
                                                                           &fracture );

              if( added )
              {
                GEOSX_LOG_LEVEL_RANK_0( 2, "Element " << cellIndex << " is fractured" );

                // Add the information to the CellElementSubRegion
                subRegion.addFracturedElement( cellIndex, localNumberOfSurfaceElems );

                newObjects.newElements[ {embeddedSurfaceRegion.getIndexInParent(), embeddedSurfaceSubRegion.getIndexInParent()} ].insert( localNumberOfSurfaceElems );

                localNumberOfSurfaceElems++;
              }
            }
          }
        } );  // end loop over cells
      } );  // end loop over subregions
    } );  // end loop over thick planes

    // Launch kernel to compute connectivity index of each fractured element.
    elemManager.forElementSubRegionsComplete< CellElementSubRegion >(
      [&]( localIndex const, localIndex const, ElementRegionBase &, CellElementSubRegion & subRegion )
    {
      finiteElement::FiniteElementBase & subRegionFE = subRegion.template getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::dispatchlowOrder3D( subRegionFE, [&] ( auto const finiteElement )
      {
        using FE_TYPE = decltype( finiteElement );

        auto kernel = CIcomputationKernel< FE_TYPE >( finiteElement,
                                                      nodeManager,
                                                      subRegion,
                                                      embeddedSurfaceSubRegion );

        using KERNEL_TYPE = decltype( kernel );

        KERNEL_TYPE::template launchCIComputationKernel< parallelDevicePolicy< 32 >, KERNEL_TYPE >( kernel );
      } );
    } );


    // add all new nodes to newObject list
    for( localIndex ni = 0; ni < embSurfNodeManager.size(); ni++ )
    {
      newObjects.newNodes.insert( ni );
    }

    // Set the ghostRank form the parent cell
    ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const & cellElemGhostRank =
      elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

    embeddedSurfaceSubRegion.inheritGhostRank( cellElemGhostRank );

    setGlobalIndices( elemManager, embSurfNodeManager, embeddedSurfaceSubRegion );

    embeddedSurfacesParallelSynchronization::sychronizeTopology( meshLevel,
                                                                 domain.getNeighbors(),
                                                                 newObjects,
                                                                 m_mpiCommOrder,
                                                                 this->m_fractureRegionName );

    addEmbeddedElementsToSets( elemManager, embeddedSurfaceSubRegion );

    EmbeddedSurfaceSubRegion::NodeMapType & embSurfToNodeMap = embeddedSurfaceSubRegion.nodeList();

    // Populate EdgeManager for embedded surfaces.
    EdgeManager & embSurfEdgeManager = meshLevel.getEmbSurfEdgeManager();

    EmbeddedSurfaceSubRegion::EdgeMapType & embSurfToEdgeMap = embeddedSurfaceSubRegion.edgeList();

    localIndex numOfPoints = embSurfNodeManager.size();

    // Create the edges
    embSurfEdgeManager.buildEdges( numOfPoints, embSurfToNodeMap.toViewConst(), embSurfToEdgeMap );
    // Node to cell map
    embSurfNodeManager.setElementMaps( elemManager );
    // Node to edge map
    embSurfNodeManager.setEdgeMaps( embSurfEdgeManager );
    embSurfNodeManager.compressRelationMaps();
  } );
}

void EmbeddedSurfaceGenerator::initializePostInitialConditionsPreSubGroups()
{}

void EmbeddedSurfaceGenerator::propagationStep( DomainPartition & domain, 
                                                R1Tensor & currentTip,             
                                                R1Tensor & targetTip,
                                                localIndex & tipElementIndex)
{
  GEOSX_MARK_FUNCTION;
  // std::cout<<"propagation step"<<std::endl;
  // std::cout<<"currentTip_x: "<<currentTip[0]<<std::endl;
  // std::cout<<"currentTip_y: "<<currentTip[1]<<std::endl;
  // std::cout<<"targetTip_x: "<<targetTip[0]<<std::endl;
  // std::cout<<"targetTip_y: "<<targetTip[1]<<std::endl;
  localIndex er = 0; //should be the element region number
  localIndex esr = 0; //should be the element subregion number
  // Get geometric object manager
  GeometricObjectManager & geometricObjManager = GeometricObjectManager::getInstance();  
  // create BoundedPlane from currentTip to newTip 
  std::unique_ptr<BoundedPlane> newFracturePlane = std::make_unique<BoundedPlane>( currentTip[0], currentTip[1], targetTip[0],
                                                                                    targetTip[1], "newFracturePlane", &geometricObjManager );
  // BoundedPlane * newFracturePlane = new BoundedPlane( currentTip[0], currentTip[1], targetTip[0],
  //                                                     targetTip[1], "newFracturePlane", &geometricObjManager );
  // Get meshLevel
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getBaseDiscretization();
  // Get managers
  ElementRegionManager & elemManager = meshLevel.getElemManager();
  NodeManager & nodeManager = meshLevel.getNodeManager();
  EmbeddedSurfaceNodeManager & embSurfNodeManager = meshLevel.getEmbSurfNodeManager();
  EdgeManager & edgeManager = meshLevel.getEdgeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord = nodeManager.referencePosition();
  // Get EmbeddedSurfaceSubRegions
  SurfaceElementRegion & embeddedSurfaceRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );
  EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion = embeddedSurfaceRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );
  localIndex localNumberOfSurfaceElems = embeddedSurfaceSubRegion.size(); //this should not be zero
  NewObjectLists newObjects; //this should probably not be re-initialized
  // begin geometric operations
  real64 const planeCenter[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( newFracturePlane->getCenter() );
  real64 const normalVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( newFracturePlane->getNormal() );
  // Initialize variables
  globalIndex nodeIndex;
  integer isPositive, isNegative;
  real64 distVec[ 3 ];    

  //COPYING FROM EXISTING FUNCTION - NEEDS CLEANING

  // Get sub-region geometric properties
  ElementRegionBase & elementRegion = elemManager.getRegion( er );
  CellElementSubRegion & subRegion = elementRegion.getSubRegion<CellElementSubRegion>( esr );
  arrayView2d< localIndex const, cells::NODE_MAP_USD > const cellToNodes = subRegion.nodeList();
  FixedOneToManyRelation const & cellToEdges = subRegion.edgeList();
  arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
  auto fracturedElements = subRegion.fracturedElementsList();
  //save old embNodesPositions to find new embNodes later
  localIndex embSurfNodesNumberOld = embSurfNodeManager.size();
  //this part probably needs explicit dynamic allocation
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & embSurfNodesPosOld = embSurfNodeManager.referencePosition();
  bool added = false;

  if( ghostRank[tipElementIndex] < 0 ) 
  {
    //this tests if the element tipElementIndex is cut by the plane connecting currentTip to targetTip
    isPositive = 0;
    isNegative = 0;
    for( localIndex kn = 0; kn < subRegion.numNodesPerElement(); kn++ )
    {
      nodeIndex = cellToNodes[tipElementIndex][kn];
      LvArray::tensorOps::copy< 3 >( distVec, nodesCoord[nodeIndex] );
      LvArray::tensorOps::subtract< 3 >( distVec, planeCenter );
      // check if the dot product is zero
      if( LvArray::tensorOps::AiBi< 3 >( distVec, normalVector ) > 0 )
      {
        isPositive = 1;
      }
      else if( LvArray::tensorOps::AiBi< 3 >( distVec, normalVector ) < 0 )
      {
        isNegative = 1;
      }
    } // end loop over nodes
    if( isPositive * isNegative == 1 ) //if there are nodes on both sides of the fracture plane
    {
      //std::cout<<"plane name: "<<newFracturePlane->getName()<<std::endl;
      added = embeddedSurfaceSubRegion.addNewEmbeddedSurface( tipElementIndex,
                                                              er,
                                                              esr,
                                                              nodeManager,
                                                              embSurfNodeManager,
                                                              edgeManager,
                                                              cellToEdges,
                                                              newFracturePlane.get() );

      if( added )
      {
        GEOSX_LOG_LEVEL_RANK_0( 2, "Element " << tipElementIndex << " is fractured" );
        // Add the information to the CellElementSubRegion
        subRegion.addFracturedElement( tipElementIndex, localNumberOfSurfaceElems );
        newObjects.newElements[ {embeddedSurfaceRegion.getIndexInParent(), embeddedSurfaceSubRegion.getIndexInParent()} ].insert( localNumberOfSurfaceElems );
        localNumberOfSurfaceElems++;
      }
    }
  }

  // Launch kernel to compute connectivity index of each fractured element.
  elemManager.forElementSubRegionsComplete< CellElementSubRegion >(
    [&]( localIndex const, localIndex const, ElementRegionBase &, CellElementSubRegion & subRegion )
  {
    finiteElement::FiniteElementBase & subRegionFE = subRegion.template getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

    finiteElement::dispatchlowOrder3D( subRegionFE, [&] ( auto const finiteElement )
    {
      using FE_TYPE = decltype( finiteElement );

      auto kernel = CIcomputationKernel< FE_TYPE >( finiteElement,
                                                    nodeManager,
                                                    subRegion,
                                                    embeddedSurfaceSubRegion );

      using KERNEL_TYPE = decltype( kernel );

      KERNEL_TYPE::template launchCIComputationKernel< parallelDevicePolicy< 32 >, KERNEL_TYPE >( kernel );
    } );
  } );

  // add all new nodes to newObject list
  // also, get index of edges that were just cut
  array1d< globalIndex > newEdges;
  array1d< globalIndex > & parentEdgeGlobalIndex = embSurfNodeManager.getParentEdgeGlobalIndex();
  localIndex embSurfNodeNumberUpdated = embSurfNodeManager.size();
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & embSurfNodesPosUpdated = embSurfNodeManager.referencePosition();

  if(added)
  {
    bool isNew;
    //array1d< localIndex > elemNodes( embSurfNodeNumberUpdated );      
    for( localIndex ni = 0; ni < embSurfNodeNumberUpdated; ni++ )
    {
      //test if node is new
      real64 distance[3];
      isNew = true;
      for( localIndex h=0; h < embSurfNodesNumberOld; h++ )
      {
        LvArray::tensorOps::copy< 3 >( distance, embSurfNodesPosUpdated[ ni ] );
        LvArray::tensorOps::subtract< 3 >( distance, embSurfNodesPosOld[ h ] );
        if( LvArray::tensorOps::l2Norm< 3 >( distance ) < 1e-4 )
        {
          isNew = false;
          break;
        }
      }
      if( isNew )
      {
        //update current crack tip
        currentTip[0] = embSurfNodesPosUpdated[ ni ][0];
        currentTip[1] = embSurfNodesPosUpdated[ ni ][1];
        //add parent edge index to array
        newEdges.resize(newEdges.size()+1);
        newEdges[newEdges.size()-1] = parentEdgeGlobalIndex[ni];
        newObjects.newNodes.insert( ni );
      }    
    }
  }

  if(added)
  {
    //given the two edges that were just cut, we find the only face that was cut
    //this face is connected with two elements, one that was cut, and one that is intact
    //we want to update tipElementIndex to become the intact one (it was the one that we just cut)
    //newEdges is an array of globalIndices with the 2 edges that were cut in this step
    auto edgeToFaceList = edgeManager.faceList();
    auto faceToElementList = faceManager.elementList();
    localIndex cutFace;
    for(auto && face:edgeToFaceList[newEdges[0]]) //newEdge is globalIndexed, we probably need a way to get the localIndex of newEdges[0] and newEdges[1]
    {
      if(edgeToFaceList.contains(newEdges[1],face)) //there is a contains() function in ArrayOfSets, does it work here?
      {
        cutFace = face; //only one face is connected to both edges
      }
    }
    for(auto && element:faceToElementList[cutFace]) //does this iterator exist? //what if the elements that share this face are in different subRegions
    {
      if(!(element == tipElementIndex)) //can element be compared to a localIndex?
      {
        std::cout<<"modifying tip element, new tip element: "<<element<<std::endl;
        tipElementIndex = element;
        break; 
      }
    }
  }

  // Set the ghostRank form the parent cell
  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const & cellElemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  embeddedSurfaceSubRegion.inheritGhostRank( cellElemGhostRank );

  setGlobalIndices( elemManager, embSurfNodeManager, embeddedSurfaceSubRegion );

  embeddedSurfacesParallelSynchronization::sychronizeTopology( meshLevel,
                                                                domain.getNeighbors(),
                                                                newObjects,
                                                                m_mpiCommOrder,
                                                                this->m_fractureRegionName );

  addEmbeddedElementsToSets( elemManager, embeddedSurfaceSubRegion );

  EmbeddedSurfaceSubRegion::NodeMapType & embSurfToNodeMap = embeddedSurfaceSubRegion.nodeList();

  // Populate EdgeManager for embedded surfaces.
  EdgeManager & embSurfEdgeManager = meshLevel.getEmbSurfEdgeManager();

  EmbeddedSurfaceSubRegion::EdgeMapType & embSurfToEdgeMap = embeddedSurfaceSubRegion.edgeList();

  localIndex numOfPoints = embSurfNodeManager.size();

  // Create the edges
  embSurfEdgeManager.buildEdges( numOfPoints, embSurfToNodeMap.toViewConst(), embSurfToEdgeMap );
  // Node to cell map
  embSurfNodeManager.setElementMaps( elemManager );
  // Node to edge map
  embSurfNodeManager.setEdgeMaps( embSurfEdgeManager );
  embSurfNodeManager.compressRelationMaps();
  /*
  * This should be the method that generates new fracture elements based on the propagation criterion of choice.
  */
  // Add the embedded elements to the fracture stencil.
  //delete newFracturePlane;
  //addToFractureStencil( domain );
}
                                               

real64 EmbeddedSurfaceGenerator::solverStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                             real64 const & GEOSX_UNUSED_PARAM( dt ),
                                             const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                             DomainPartition & domain )
{
  real64 rval = 0;
  /*
   * This should be the method that generates new fracture elements based on the propagation criterion of choice.
   */
  // Add the embedded elements to the fracture stencil.
  addToFractureStencil( domain );

  return rval;
}

void EmbeddedSurfaceGenerator::addToFractureStencil( DomainPartition & domain )
{
  // Add embedded elements to the fracture Stencil
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getBaseDiscretization();

    for( localIndex a=0; a<fvManager.numSubGroups(); ++a )
    {
      FluxApproximationBase * const fluxApprox = fvManager.getGroupPointer< FluxApproximationBase >( a );
      if( fluxApprox!=nullptr )
      {
        fluxApprox->addEmbeddedFracturesToStencils( meshLevel, this->m_fractureRegionName );
      }
    }
  }

}

void EmbeddedSurfaceGenerator::setGlobalIndices( ElementRegionManager & elemManager,
                                                 EmbeddedSurfaceNodeManager & embSurfNodeManager,
                                                 EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion )
{
  // Add new globalIndices
  int const thisRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const commSize = MpiWrapper::commSize( MPI_COMM_GEOSX );

  localIndex_array numberOfSurfaceElemsPerRank( commSize );
  localIndex_array globalIndexOffset( commSize );
  MpiWrapper::allGather( embeddedSurfaceSubRegion.size(), numberOfSurfaceElemsPerRank );

  globalIndexOffset[0] = 0; // offSet for the globalIndex
  localIndex totalNumberOfSurfaceElements = numberOfSurfaceElemsPerRank[ 0 ];  // Sum across all ranks
  for( int rank = 1; rank < commSize; ++rank )
  {
    globalIndexOffset[rank] = globalIndexOffset[rank - 1] + numberOfSurfaceElemsPerRank[rank - 1];
    totalNumberOfSurfaceElements += numberOfSurfaceElemsPerRank[rank];
  }

  GEOSX_LOG_LEVEL_RANK_0( 1, "Number of embedded surface elements: " << totalNumberOfSurfaceElements );

  arrayView1d< globalIndex > const & elemLocalToGlobal = embeddedSurfaceSubRegion.localToGlobalMap();

  forAll< serialPolicy >( embeddedSurfaceSubRegion.size(), [&, elemLocalToGlobal] ( localIndex const ei )
  {
    elemLocalToGlobal( ei ) = ei + globalIndexOffset[ thisRank ] + elemManager.maxGlobalIndex() + 1;
    embeddedSurfaceSubRegion.updateGlobalToLocalMap( ei );
  } );

  embeddedSurfaceSubRegion.setMaxGlobalIndex();

  elemManager.setMaxGlobalIndex();

  // Nodes global indices
  localIndex_array numberOfNodesPerRank( commSize );
  MpiWrapper::allGather( embSurfNodeManager.size(), numberOfNodesPerRank );

  globalIndexOffset[0] = 0; // offSet for the globalIndex
  for( int rank = 1; rank < commSize; ++rank )
  {
    globalIndexOffset[rank] = globalIndexOffset[rank - 1] + numberOfNodesPerRank[rank - 1];
  }

  arrayView1d< globalIndex > const & nodesLocalToGlobal = embSurfNodeManager.localToGlobalMap();
  forAll< serialPolicy >( embSurfNodeManager.size(), [=, &embSurfNodeManager, &globalIndexOffset] ( localIndex const ni )
  {
    nodesLocalToGlobal( ni ) = ni + globalIndexOffset[ thisRank ];
    embSurfNodeManager.updateGlobalToLocalMap( ni );
  } );
}

void EmbeddedSurfaceGenerator::addEmbeddedElementsToSets( ElementRegionManager const & elemManager,
                                                          EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion )
{
  // We want to create the sets for the embeddedSurfaceSubRegion which was empty when they
  // were created for the other subRegions. This way, for example, if the parent cell belongs to the set "source"
  // the embedded element will belong to the same set.

  dataRepository::Group & setGroupEmbSurf =
    embeddedSurfaceSubRegion.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );

  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    dataRepository::Group const & setGroupCell = subRegion.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );

    SortedArrayView< localIndex const > const fracturedElements = subRegion.fracturedElementsList();

    ArrayOfArraysView< localIndex const > const cellToEmbSurf = subRegion.embeddedSurfacesList().toViewConst();

    forAll< serialPolicy >( fracturedElements.size(), [&,
                                                       fracturedElements,
                                                       cellToEmbSurf ] ( localIndex const ei )
    {
      localIndex const cellIndex    = fracturedElements[ei];
      localIndex const embSurfIndex = cellToEmbSurf[cellIndex][0];
      setGroupCell.forWrappers< SortedArray< localIndex > >( [&]( auto const & wrapper )
      {
        SortedArrayView< const localIndex > const & targetSetCell = wrapper.reference();
        targetSetCell.move( LvArray::MemorySpace::host );

        SortedArray< localIndex > & targetSetEmbSurf =
          setGroupEmbSurf.getWrapper< SortedArray< localIndex > >( wrapper.getName() ).reference();
        if( targetSetCell.contains( cellIndex ) )
        {
          targetSetEmbSurf.insert( embSurfIndex );
        }
      } );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        EmbeddedSurfaceGenerator,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
