/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TwoPointFluxApproximation.cpp
 *
 */
#include "TwoPointFluxApproximation.hpp"

#include "mesh/FaceElementRegion.hpp"
#include "CellElementStencilTPFA.hpp"
#include "FaceElementStencil.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;

TwoPointFluxApproximation::TwoPointFluxApproximation( std::string const & name,
                                                      Group * const parent )
  : FluxApproximationBase( name, parent )
{
  registerWrapper< CellElementStencilTPFA >( viewKeyStruct::cellStencilString )->
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper< FaceElementStencil >( viewKeyStruct::fractureStencilString )->
    setRestartFlags( RestartFlags::NO_WRITE );

}

namespace
{

void makeFullTensor( R1Tensor const & values, R2SymTensor & result )
{
  result = 0.0;
  R1Tensor axis;
  R2SymTensor temp;

  // assemble full tensor from eigen-decomposition
  for( unsigned icoord = 0; icoord < 3; ++icoord )
  {
    // assume principal axis aligned with global coordinate system
    axis = 0.0;
    axis( icoord ) = 1.0;

    // XXX: is there a more elegant way to do this?
    temp.dyadic_aa( axis );
    temp *= values( icoord );
    result += temp;
  }
}

void averageTensorComponent(R1Tensor const & values, real64 & avCoeff)
{
  avCoeff = 0;
  for( unsigned icoord = 0; icoord < 3; ++icoord )
  {
    avCoeff += values[icoord];
  }
  avCoeff/=3;
}

} // close namespace



void TwoPointFluxApproximation::computeCellStencil( DomainPartition const & domain )
{
  MeshBody const * const meshBody = domain.getMeshBody( 0 );
  MeshLevel const * const mesh = meshBody->getMeshLevel( 0 );
  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();


  CellElementStencilTPFA & stencil = this->getReference< CellElementStencilTPFA >( viewKeyStruct::cellStencilString );

  arrayView2d< localIndex const > const & elemRegionList = faceManager->elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList = faceManager->elementSubRegionList();
  arrayView2d< localIndex const > const & elemList = faceManager->elementList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor const > > const
  elemCenter = elemManager->ConstructViewAccessor< array1d< R1Tensor >,
                                                   arrayView1d< R1Tensor const > >( CellBlock::
                                                                                      viewKeyStruct::
                                                                                      elementCenterString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor const > > const
  coefficient = elemManager->ConstructViewAccessor< array1d< R1Tensor >,
                                                    arrayView1d< R1Tensor const > >( m_coeffName );

  arrayView1d< integer const > const & faceGhostRank =
    faceManager->getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager->nodeList().toViewConst();

  // make a list of region indices to be included
  SortedArray< localIndex > regionFilter;
  for( string const & regionName : m_targetRegions )
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  constexpr localIndex numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight, faceWeightInv;

  stackArray1d< localIndex, numElems > stencilCellsRegionIndex( numElems );
  stackArray1d< localIndex, numElems > stencilCellsSubRegionIndex( numElems );
  stackArray1d< localIndex, numElems > stencilCellsIndex( numElems );
  stackArray1d< real64, numElems > stencilWeights( numElems );

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve( faceManager->size() );

  real64 const lengthTolerance = meshBody->getGlobalLengthScale() * this->m_areaRelTol;
  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex kf = 0; kf < faceManager->size(); ++kf )
  {
    if( faceGhostRank[kf] >= 0 || elemRegionList[kf][0] == -1 || elemRegionList[kf][1] == -1 )
      continue;

    if( !(regionFilter.contains( elemRegionList[kf][0] ) && regionFilter.contains( elemRegionList[kf][1] )) )
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], faceToNodes.sizeOfArray( kf ), X, faceCenter, faceNormal, areaTolerance );

    if( faceArea < areaTolerance )
      continue;

    faceWeightInv = 0.0;

    for( localIndex ke = 0; ke < numElems; ++ke )
    {
      if( elemRegionList[kf][ke] != -1 )
      {
        localIndex const er  = elemRegionList[kf][ke];
        localIndex const esr = elemSubRegionList[kf][ke];
        localIndex const ei  = elemList[kf][ke];

        cellToFaceVec = faceCenter;
        cellToFaceVec -= elemCenter[er][esr][ei];

        // ensure normal orientation outward of first cell
        if( ke == 0 && Dot( cellToFaceVec, faceNormal ) < 0.0 )
        {
          faceNormal *= -1;
        }

        if( ke == 1 )
        {
          cellToFaceVec *= -1.0;
        }

        real64 const c2fDistance = cellToFaceVec.Normalize();

        // assemble full coefficient tensor from principal axis/components
        makeFullTensor( coefficient[er][esr][ei], coefTensor );

        faceConormal.AijBj( coefTensor, faceNormal );
        real64 halfWeight = Dot( cellToFaceVec, faceConormal );

        // correct negative weight issue arising from non-K-orthogonal grids
        if( halfWeight < 0.0 )
        {
          faceConormal.AijBj( coefTensor, cellToFaceVec );
          halfWeight = Dot( cellToFaceVec, faceConormal );
        }

        halfWeight *= faceArea / c2fDistance;
        halfWeight = std::max( halfWeight, weightTolerance );

        faceWeightInv += 1.0 / halfWeight;
      }
    }

    GEOSX_ASSERT( faceWeightInv > 0.0 );
    faceWeight = 1.0 / faceWeightInv;

    for( localIndex ke = 0; ke < numElems; ++ke )
    {
      stencilCellsRegionIndex[ke] = elemRegionList[kf][ke];
      stencilCellsSubRegionIndex[ke] = elemSubRegionList[kf][ke];
      stencilCellsIndex[ke] = elemList[kf][ke];
      stencilWeights[ke] = faceWeight * (ke == 0 ? 1 : -1);
    }
    stencil.add( CellElementStencilTPFA::NUM_POINT_IN_FLUX,
                 stencilCellsRegionIndex,
                 stencilCellsSubRegionIndex,
                 stencilCellsIndex,
                 stencilWeights.data(),
                 kf );
  }
//  stencil.compress();
}


void TwoPointFluxApproximation::addToFractureStencil( DomainPartition & domain,
                                                      string const & faceElementRegionName,
                                                      bool const initFlag )
{
  MeshLevel * const mesh = domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();
  EdgeManager const * const edgeManager = mesh->getEdgeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor > > const
  elemCenter = elemManager->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor > >( CellBlock::
                                                                                                     viewKeyStruct::
                                                                                                     elementCenterString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor > > const
  coefficient = elemManager->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor > >( m_coeffName );

  arrayView1d< real64 const >   const & faceArea = faceManager->faceArea();
  arrayView1d< R1Tensor const > const & faceCenter = faceManager->faceCenter();
  arrayView1d< R1Tensor const > const & faceNormal = faceManager->faceNormal();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager->referencePosition();

  FaceElementStencil & fractureStencil = getReference< FaceElementStencil >( viewKeyStruct::fractureStencilString );
  CellElementStencilTPFA & cellStencil = getReference< CellElementStencilTPFA >( viewKeyStruct::cellStencilString );

  FaceElementRegion * const fractureRegion = elemManager->GetRegion< FaceElementRegion >( faceElementRegionName );
  localIndex const fractureRegionIndex = fractureRegion->getIndexInParent();

  FaceElementSubRegion * const fractureSubRegion = fractureRegion->GetSubRegion< FaceElementSubRegion >( "default" );
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion->faceList();

  arrayView1d< integer const > const & ghostRank = fractureSubRegion->ghostRank();

  arrayView1d< localIndex const > const &
  fractureConnectorsToEdges = edgeManager->getReference< array1d< localIndex > >( EdgeManager::
                                                                                    viewKeyStruct::
                                                                                    fractureConnectorEdgesToEdgesString );

  ArrayOfArraysView< localIndex const > const &
  fractureConnectorsToFaceElements =
    edgeManager->getReference< ArrayOfArrays< localIndex > >( EdgeManager::
                                                                viewKeyStruct::
                                                                fractureConnectorsEdgesToFaceElementsIndexString );

  FixedToManyElementRelation const &
  faceElementsToCells = fractureSubRegion->m_faceElementsToCells;

  localIndex constexpr maxElems = FaceElementStencil::MAX_STENCIL_SIZE;

  stackArray1d< localIndex, maxElems > stencilCellsRegionIndex;
  stackArray1d< localIndex, maxElems > stencilCellsSubRegionIndex;
  stackArray1d< localIndex, maxElems > stencilCellsIndex;
  stackArray1d< real64, maxElems > stencilWeights;

  stackArray1d< R1Tensor, maxElems > stencilCellCenterToEdgeCenters;
  stackArray1d< integer, maxElems > isGhostConnectors;

  arrayView1d< integer const > const & edgeGhostRank = edgeManager->ghostRank();

  // TODO Note that all of this initialization should be performed elsewhere. This is just here because it was
  // convenient, but it is not appropriate to have physics based initialization in the flux approximator.
#if !defined(SET_CREATION_DISPLACEMENT)
  static_assert( true, "must have SET_CREATION_DISPLACEMENT defined" );
#endif
#if SET_CREATION_DISPLACEMENT==1
  ArrayOfArraysView< localIndex const > const & faceToNodesMap = faceManager->nodeList();
  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & incrementalDisplacement = nodeManager->incrementalDisplacement();
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & totalDisplacement = nodeManager->totalDisplacement();
  arrayView1d< real64 > const & aperture = fractureSubRegion->getReference< array1d< real64 > >( "elementAperture" );
#endif


#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  arrayView1d< real64 > const &
  apertureF = fractureSubRegion->getReference< array1d< real64 > >( "apertureAtFailure" );
#endif

#if !defined(ALLOW_CREATION_MASS)
  static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif
#if ALLOW_CREATION_MASS==0
  arrayView1d< real64 > const &
  dens = fractureSubRegion->getReference< array1d< real64 > >( "densityOld" );
#endif


#if SET_CREATION_PRESSURE==1
  arrayView1d< real64 > const &
  fluidPressure = fractureSubRegion->getReference< array1d< real64 > >( "pressure" );
  // Set the new face elements to some unphysical numbers to make sure they get set by the following routines.
  for( localIndex const kfe : fractureSubRegion->m_newFaceElements )
  {
#if !defined(SET_CREATION_PRESSURE)
    static_assert( true, "must have SET_CREATION_PRESSURE defined" );
#endif
#if SET_CREATION_PRESSURE==1
    if( initFlag )
    {
      fluidPressure[kfe] = 1.0e99;
    }
#endif
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
    apertureF[kfe] = aperture[kfe];
#endif
#if !defined(SET_CREATION_DISPLACEMENT)
    static_assert( true, "must have SET_CREATION_DISPLACEMENT defined" );
#endif
#if SET_CREATION_DISPLACEMENT==1
    if( initFlag )
    {
      aperture[kfe] = 1.0e99;
    }
#endif
  }
#endif
  std::set< localIndex > allNewElems( fractureSubRegion->m_newFaceElements.begin(),
                                      fractureSubRegion->m_newFaceElements.end() );


  // add new connectors/connections between face elements to the fracture stencil
  for( auto const fci : edgeManager->m_recalculateFractureConnectorEdges )
  {
    localIndex const numElems = fractureConnectorsToFaceElements.sizeOfArray( fci );
    // only do this if there are more than one element attached to the connector
    localIndex const edgeIndex = fractureConnectorsToEdges[fci];

    //    if( edgeGhostRank[edgeIndex] < 0 && numElems > 1 )
    {

      GEOSX_ERROR_IF( numElems > maxElems, "Max stencil size exceeded by fracture-fracture connector " << fci );
      stencilCellsRegionIndex.resize( numElems );
      stencilCellsSubRegionIndex.resize( numElems );
      stencilCellsIndex.resize( numElems );
      stencilWeights.resize( numElems );

      stencilCellCenterToEdgeCenters.resize( numElems );
      isGhostConnectors.resize( numElems );

      // get edge geometry
      R1Tensor const edgeCenter = edgeManager->calculateCenter( edgeIndex, X );
      real64 const edgeLength = edgeManager->calculateLength( edgeIndex, X ).L2_Norm();

      real64 initialPressure = 1.0e99;
#if SET_CREATION_DISPLACEMENT==1
      real64 initialAperture = 1.0e99;
#endif
      SortedArray< localIndex > newElems;

      // loop over all face elements attached to the connector and add them to the stencil
      for( localIndex kfe=0; kfe<numElems; ++kfe )
      {
        localIndex const fractureElementIndex = fractureConnectorsToFaceElements[fci][kfe];

        // use straight difference between the edge center and face center for gradient length...maybe do something
        // better here?? TODO
        R1Tensor cellCenterToEdgeCenter = edgeCenter;
        cellCenterToEdgeCenter -= faceCenter[ faceMap[fractureElementIndex][0] ];

        // form the CellStencil entry
        stencilCellsRegionIndex[kfe] = fractureRegionIndex;
        stencilCellsSubRegionIndex[kfe] = 0;
        stencilCellsIndex[kfe] = fractureElementIndex;

        stencilWeights[kfe] =  1.0 / 12.0 * edgeLength / cellCenterToEdgeCenter.L2_Norm();

        stencilCellCenterToEdgeCenters[kfe] = cellCenterToEdgeCenter;

        isGhostConnectors[kfe] = edgeGhostRank[edgeIndex];

        // code to initialize new face elements with pressures from neighbors
        if( fractureSubRegion->m_newFaceElements.count( fractureElementIndex )==0 )
        {
          initialPressure = std::min( initialPressure, fluidPressure[fractureElementIndex] );
#if SET_CREATION_DISPLACEMENT==1
          initialAperture = std::min( initialAperture, aperture[fractureElementIndex] );
#endif
        }
        else
        {
          newElems.insert( fractureElementIndex );
          allNewElems.insert( fractureElementIndex );
        }
      }

      // loop over new face elements attached to this connector
      for( localIndex const newElemIndex : newElems )
      {
        // set the aperture/fluid pressure for the new face element to be the minimum
        // of the existing value, smallest aperture/pressure from a connected face element.
//        aperture[newElemIndex] = std::min(aperture[newElemIndex], initialAperture);
#if !defined(SET_CREATION_PRESSURE)
        static_assert( true, "must have SET_CREATION_PRESSURE defined" );
#endif
#if SET_CREATION_PRESSURE==1
        if( initFlag )
        {
          fluidPressure[newElemIndex] = std::min( fluidPressure[newElemIndex], initialPressure );
        }
#endif

#if !defined(SET_CREATION_DISPLACEMENT)
        static_assert( true, "must have SET_CREATION_DISPLACEMENT defined" );
#endif
#if SET_CREATION_DISPLACEMENT==1
        if( initFlag )
        {
          localIndex const faceIndex0 = faceMap( newElemIndex, 0 );
          localIndex const faceIndex1 = faceMap( newElemIndex, 1 );

          localIndex const numNodesPerFace = faceToNodesMap.sizeOfArray( faceIndex0 );

          bool zeroDisp = true;

          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            localIndex const node0 = faceToNodesMap( faceIndex0, a );
            localIndex const node1 = faceToNodesMap( faceIndex1, a==0 ? a : numNodesPerFace-a );
            if( fabs( LvArray::tensorOps::norm2( totalDisplacement[node0] ) ) > 1.0e-99 &&
                fabs( LvArray::tensorOps::norm2( totalDisplacement[node1] ) ) > 1.0e-99 )
            {
              zeroDisp = false;
            }
          }
          if( zeroDisp )
          {
            aperture[newElemIndex] = 0;
          }
        }
#endif
      }

      // add/overwrite the stencil for index fci
      fractureStencil.add( numElems,
                           stencilCellsRegionIndex,
                           stencilCellsSubRegionIndex,
                           stencilCellsIndex,
                           stencilWeights.data(),
                           fci );

      fractureStencil.add( numElems,
                           stencilCellCenterToEdgeCenters.data(),
                           isGhostConnectors.data(),
                           fci );
    }
  }

  if( initFlag )
  {
    SortedArray< localIndex > touchedNodes;
    for( localIndex const newElemIndex : allNewElems )
    {
      // if the value of pressure was not set, then set it to zero and punt.
      if( fluidPressure[newElemIndex] > 1.0e98 )
      {
        fluidPressure[newElemIndex] = 0.0;
      }
#if !defined(ALLOW_CREATION_MASS)
      static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif
#if ALLOW_CREATION_MASS==0
      // set the initial density of the face element to 0 to enforce mass conservation ( i.e. no creation of mass)
      dens[newElemIndex] = 0.0;
#endif

#if !defined(SET_CREATION_DISPLACEMENT)
      static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif
#if SET_CREATION_DISPLACEMENT==1
      // If the aperture has been set, then we can set the estimate of displacements.
      if( aperture[newElemIndex] < 1e98 )
      {
        localIndex const faceIndex0 = faceMap( newElemIndex, 0 );
        localIndex const faceIndex1 = faceMap( newElemIndex, 1 );

        R1Tensor newDisp = faceNormal( faceIndex0 );
        newDisp *= -aperture[newElemIndex];
        localIndex const numNodesPerFace = faceToNodesMap.sizeOfArray( faceIndex0 );
        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          localIndex const node0 = faceToNodesMap( faceIndex0, a );
          localIndex const node1 = faceToNodesMap( faceIndex1, a==0 ? a : numNodesPerFace-a );

          touchedNodes.insert( node0 );
          touchedNodes.insert( node1 );

          if( node0 != node1 && touchedNodes.count( node0 )==0 )
          {
            incrementalDisplacement[node0] += newDisp;
            totalDisplacement[node0] += newDisp;
            incrementalDisplacement[node1] -= newDisp;
            totalDisplacement[node1] -= newDisp;
          }
        }
      }
      if( this->getLogLevel() > 1 )
      {
        printf( "New elem index, init aper, init press = %4ld, %4.2e, %4.2e \n",
                newElemIndex,
                aperture[newElemIndex],
                fluidPressure[newElemIndex] );
      }
#endif
    }
  }

  // add connections for FaceElements to/from CellElements.
  {
    arrayView2d< localIndex const > const & elemRegionList = faceElementsToCells.m_toElementRegion;
    arrayView2d< localIndex const > const & elemSubRegionList = faceElementsToCells.m_toElementSubRegion;
    arrayView2d< localIndex const > const & elemList = faceElementsToCells.m_toElementIndex;
    for( localIndex const kfe : fractureSubRegion->m_newFaceElements )
//    for( localIndex kfe=0 ; kfe<faceElementsToCells.size(0) ; ++kfe )
    {
      if( ghostRank[kfe] < 0 )
      {
        localIndex const numElems = faceElementsToCells.size( 1 );

        GEOSX_ERROR_IF( numElems > maxElems, "Max stencil size exceeded by fracture-cell connector " << kfe );
        stencilCellsRegionIndex.resize( numElems );
        stencilCellsSubRegionIndex.resize( numElems );
        stencilCellsIndex.resize( numElems );
        stencilWeights.resize( numElems );

        R2SymTensor coefTensor;
        R1Tensor cellToFaceVec;
        R1Tensor faceConormal;

        // remove cell-to-cell connections from cell stencil and add in new connections
        if( cellStencil.zero( faceMap[kfe][0] ) )
        {
          for( localIndex ke = 0; ke < numElems; ++ke )
          {
            localIndex const faceIndex = faceMap[kfe][ke];
            localIndex const er  = elemRegionList[kfe][ke];
            localIndex const esr = elemSubRegionList[kfe][ke];
            localIndex const ei  = elemList[kfe][ke];

            cellToFaceVec = faceCenter[faceIndex];
            cellToFaceVec -= elemCenter[er][esr][ei];

            real64 const c2fDistance = cellToFaceVec.Normalize();

            // assemble full coefficient tensor from principal axis/components
            makeFullTensor( coefficient[er][esr][ei], coefTensor );

            faceConormal.AijBj( coefTensor, faceNormal[faceIndex] );
            real64 const ht = Dot( cellToFaceVec, faceConormal ) * faceArea[faceIndex] / c2fDistance;

            // assume the h for the faceElement to the connector (Face) is zero. thus the weights are trivial.
            stencilCellsRegionIndex[0] = er;
            stencilCellsSubRegionIndex[0] = esr;
            stencilCellsIndex[0] = ei;
            stencilWeights[0] =  ht;

            stencilCellsRegionIndex[1] = fractureRegionIndex;
            stencilCellsSubRegionIndex[1] = 0;
            stencilCellsIndex[1] = kfe;
            stencilWeights[1] = -ht;

            cellStencil.add( 2,
                             stencilCellsRegionIndex,
                             stencilCellsSubRegionIndex,
                             stencilCellsIndex,
                             stencilWeights.data(),
                             faceIndex );
          }
        }
      }
    }
  }
}

void TwoPointFluxApproximation::addEDFracToFractureStencil( DomainPartition & domain,
                                                            string const &  embeddedSurfaceRegionName )
{
  MeshLevel * const mesh = domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  EdgeManager * const embSurfEdgeManager = mesh->getEmbdSurfEdgeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor > > const
  coefficient = elemManager->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor > >( m_coeffName );

  // Get the stencils
  FaceElementStencil & fractureStencil = getReference< FaceElementStencil >( viewKeyStruct::fractureStencilString );
  CellElementStencilTPFA & cellStencil = getReference< CellElementStencilTPFA >( viewKeyStruct::cellStencilString );

  EmbeddedSurfaceRegion * const fractureRegion = elemManager->GetRegion< EmbeddedSurfaceRegion >( embeddedSurfaceRegionName );
  localIndex const fractureRegionIndex = fractureRegion->getIndexInParent();

  EmbeddedSurfaceSubRegion * const fractureSubRegion = fractureRegion->GetSubRegion< EmbeddedSurfaceSubRegion >( "default" );

  arrayView1d< real64 const >   const & fractureElemArea   = fractureSubRegion->getElementArea();
  arrayView1d< R1Tensor const > const & fractureElemCenter = fractureSubRegion->getElementCenter();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager->referencePosition();

  EdgeManager::FaceMapType & edgeToEmbSurfacesMap = embSurfEdgeManager->faceList();

//arrayView1d< integer const > const & ghostRank = fractureSubRegion->ghostRank();

  localIndex constexpr maxElems = FaceElementStencil::MAX_STENCIL_SIZE;

  stackArray1d< localIndex, maxElems > stencilCellsRegionIndex;
  stackArray1d< localIndex, maxElems > stencilCellsSubRegionIndex;
  stackArray1d< localIndex, maxElems > stencilCellsIndex;
  stackArray1d< real64, maxElems > stencilWeights;

  stackArray1d< R1Tensor, maxElems > stencilCellCenterToEdgeCenters;
  stackArray1d< integer, maxElems > isGhostConnectors;

  arrayView1d< integer const > const & edgeGhostRank = edgeManager->ghostRank();

  localIndex connectorIndex = 0;
  // add new connectors/connections between embedded elements to the fracture stencil
  for( localIndex ke = 0;  ke <  embSurfEdgeManager->size(); ke++ )
  {
    // for now there is no generation of new elements so we add all edges.
    localIndex const numElems = 2;  // hardcoded for now but unless there is an intersection it should always be 2.
    if (edgeToEmbSurfacesMap.sizeOfSet(ke) > 1) // to be a connector it need to be attached to at least 2 elements.
    {

      GEOSX_ERROR_IF( numElems > maxElems, "Max stencil size exceeded by fracture-fracture connector " << ke );
      stencilCellsRegionIndex.resize( numElems );
      stencilCellsSubRegionIndex.resize( numElems );
      stencilCellsIndex.resize( numElems );
      stencilWeights.resize( numElems );

      stencilCellCenterToEdgeCenters.resize( numElems );
      isGhostConnectors.resize( numElems );

      //TODO get edge geometry
      R1Tensor const edgeCenter = embSurfEdgeManager->calculateCenter( ke, X );
      real64 const edgeLength   = embSurfEdgeManager->calculateLength( ke, X ).L2_Norm();

      // loop over all embedded surface elements attached to the connector and add them to the stencil
      for( localIndex kes=0; kes < numElems; ++kes )
      {
        localIndex const fractureElementIndex = edgeToEmbSurfacesMap[ke][kes];

        // use straight difference between the edge center and face center for gradient length...maybe do something
        // better here?? TODO
        R1Tensor cellCenterToEdgeCenter = edgeCenter;
        cellCenterToEdgeCenter -= fractureElemCenter[fractureElementIndex];

        // form the CellStencil entry
        stencilCellsRegionIndex[kes]    = fractureRegionIndex;
        stencilCellsSubRegionIndex[kes] = 0;  // there is only one subregion.
        stencilCellsIndex[kes]          = fractureElementIndex;

        //TODO use the proper geometrical info to compute the weight.
        stencilWeights[kes] =  1.0 / 12.0 * edgeLength / cellCenterToEdgeCenter.L2_Norm();

        stencilCellCenterToEdgeCenters[kes] = cellCenterToEdgeCenter;

//        isGhostConnectors[kfe] = edgeGhostRank[edgeIndex];
      }

      // add/overwrite the stencil for index fci
      fractureStencil.add( numElems,
                           stencilCellsRegionIndex,
                           stencilCellsSubRegionIndex,
                           stencilCellsIndex,
                           stencilWeights.data(),
                           connectorIndex );

      fractureStencil.add( numElems,
                           stencilCellCenterToEdgeCenters.data(),
                           isGhostConnectors.data(),
                           connectorIndex );

      connectorIndex++;
    }
  }

  // Add connections EmbeddedSurface to/from CellElements.
  {
    arrayView1d< localIndex const > const & elemRegionList    = fractureSubRegion->getSurfaceToRegionList();
    arrayView1d< localIndex const > const & elemSubRegionList = fractureSubRegion->getSurfaceToSubRegionList();
    arrayView1d< localIndex const > const & elemList          = fractureSubRegion->getSurfaceToCellList();
    arrayView1d< real64 const >     const & ConnectivityIndex = fractureSubRegion->getConnectivityIndex();

    for( localIndex kes; kes  < fractureSubRegion->size(); kes++ )
    {
      //if( ghostRank[kes] < 0 )
      {
        localIndex const numElems = 1; // there is a 1 to 1 relation

        GEOSX_ERROR_IF( numElems > maxElems, "Max stencil size exceeded by fracture-cell connector " << kes );
        stencilCellsRegionIndex.resize( numElems );
        stencilCellsSubRegionIndex.resize( numElems );
        stencilCellsIndex.resize( numElems );
        stencilWeights.resize( numElems );

        localIndex const er  = elemRegionList[kes];
        localIndex const esr = elemSubRegionList[kes];
        localIndex const ei  = elemList[kes];

        // Here goes EDFM transmissibility computation.
        real64 avCoefficient;

        averageTensorComponent(coefficient[er][esr][ei], avCoefficient);

        real64 const ht = ConnectivityIndex[kes] * avCoefficient; // Using matrix perm coz assuming fracture is highly permeable for now.

        //
        stencilCellsRegionIndex[0] = er;
        stencilCellsSubRegionIndex[0] = esr;
        stencilCellsIndex[0] = ei;
        stencilWeights[0] =  ht;

        stencilCellsRegionIndex[1] = fractureRegionIndex;
        stencilCellsSubRegionIndex[1] = 0;
        stencilCellsIndex[1] = kes;
        stencilWeights[1] = -ht;

        cellStencil.add( 1,
                         stencilCellsRegionIndex,
                         stencilCellsSubRegionIndex,
                         stencilCellsIndex,
                         stencilWeights.data(),
                         kes );
      }
    }
  }
}


void TwoPointFluxApproximation::computeBoundaryStencil( DomainPartition const & domain,
                                                        SortedArrayView< localIndex const > const & faceSet,
                                                        BoundaryStencil & stencil )
{
  MeshBody const * const meshBody = domain.getMeshBody( 0 );
  MeshLevel const * const mesh = meshBody->getMeshLevel( 0 );
  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  arrayView2d< localIndex const > const & elemRegionList     = faceManager->elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList  = faceManager->elementSubRegionList();
  arrayView2d< localIndex const > const & elemList           = faceManager->elementList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor const > > const
  elemCenter = elemManager->ConstructViewAccessor< array1d< R1Tensor >,
                                                   arrayView1d< R1Tensor const > >( CellBlock::
                                                                                      viewKeyStruct::
                                                                                      elementCenterString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor const > > const
  coefficient = elemManager->ConstructViewAccessor< array1d< R1Tensor >,
                                                    arrayView1d< R1Tensor const > >( m_coeffName );

  arrayView1d< integer const > const & faceGhostRank = faceManager->getReference< integer_array >( ObjectManagerBase::
                                                                                                     viewKeyStruct::
                                                                                                     ghostRankString );

  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager->nodeList().toViewConst();

  // make a list of region indices to be included
  SortedArray< localIndex > regionFilter;
  for( string const & regionName : m_targetRegions )
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  constexpr localIndex numElems = BoundaryStencil::NUM_POINT_IN_FLUX;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight;

  stackArray1d< PointDescriptor, numElems > stencilPoints( numElems );
  stackArray1d< real64, numElems > stencilWeights( numElems );

  real64 const lengthTolerance = meshBody->getGlobalLengthScale() * this->m_areaRelTol;
  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve( faceSet.size(), 2 );
  for( localIndex kf : faceSet )
  {
    if( faceGhostRank[kf] >= 0 )
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], faceToNodes.sizeOfArray( kf ), X, faceCenter, faceNormal, areaTolerance );

    for( localIndex ke = 0; ke < numElems; ++ke )
    {
      if( elemRegionList[kf][ke] < 0 )
        continue;

      if( !regionFilter.contains( elemRegionList[kf][ke] ))
        continue;

      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];
      localIndex const ei  = elemList[kf][ke];

      cellToFaceVec = faceCenter;
      cellToFaceVec -= elemCenter[er][esr][ei];

      // ensure normal orientation outward of the cell
      if( Dot( cellToFaceVec, faceNormal ) < 0.0 )
      {
        faceNormal *= -1;
      }

      real64 const c2fDistance = cellToFaceVec.Normalize();

      // assemble full coefficient tensor from principal axis/components
      makeFullTensor( coefficient[er][esr][ei], coefTensor );

      faceConormal.AijBj( coefTensor, faceNormal );
      faceWeight = Dot( cellToFaceVec, faceConormal );

      // correct negative weight issue arising from non-K-orthogonal grids
      if( faceWeight < 0.0 )
      {
        faceConormal.AijBj( coefTensor, cellToFaceVec );
        faceWeight = Dot( cellToFaceVec, faceConormal );
      }

      faceWeight *= faceArea / c2fDistance;
      faceWeight = std::max( faceWeight, weightTolerance );

      stencilPoints[0].tag = PointDescriptor::Tag::CELL;
      stencilPoints[0].cellIndex = { er, esr, ei };
      stencilWeights[0] = faceWeight;

      stencilPoints[1].tag = PointDescriptor::Tag::FACE;
      stencilPoints[1].faceIndex = kf;
      stencilWeights[1] = -faceWeight;

      stencil.add( 2,
                   stencilPoints.data(),
                   stencilWeights.data(),
                   kf );
    }
  }
  stencil.compress();
}


REGISTER_CATALOG_ENTRY( FluxApproximationBase, TwoPointFluxApproximation, std::string const &, Group * const )

}
