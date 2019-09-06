/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

namespace geosx
{

using namespace dataRepository;

TwoPointFluxApproximation::TwoPointFluxApproximation(std::string const &name,
                                                     Group *const parent)
  : FluxApproximationBase(name, parent)
{
  registerWrapper<CellElementStencilTPFA>(viewKeyStruct::cellStencilString)->
    setRestartFlags(RestartFlags::NO_WRITE);

  registerWrapper<FaceElementStencil>(viewKeyStruct::fractureStencilString)->
    setRestartFlags(RestartFlags::NO_WRITE);

}

namespace
{

void makeFullTensor(R1Tensor const & values, R2SymTensor & result)
{
  result = 0.0;
  R1Tensor axis;
  R2SymTensor temp;

  // assemble full tensor from eigen-decomposition
  for (unsigned icoord = 0; icoord < 3; ++icoord)
  {
    // assume principal axis aligned with global coordinate system
    axis = 0.0;
    axis(icoord) = 1.0;

    // XXX: is there a more elegant way to do this?
    temp.dyadic_aa(axis);
    temp *= values(icoord);
    result += temp;
  }
}

}

void TwoPointFluxApproximation::computeCellStencil( DomainPartition const & domain  )
{
  MeshBody const * const meshBody = domain.getMeshBody(0);
  MeshLevel const * const mesh = meshBody->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();


  CellElementStencilTPFA & stencil = this->getReference<CellElementStencilTPFA>(viewKeyStruct::cellStencilString);

  arrayView2d<localIndex const> const & elemRegionList     = faceManager->elementRegionList();
  arrayView2d<localIndex const> const & elemSubRegionList  = faceManager->elementSubRegionList();
  arrayView2d<localIndex const> const & elemList           = faceManager->elementList();
  arrayView1d<R1Tensor const>   const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const elemCenter =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(
                                        CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const coefficient =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >( m_coeffName );

  arrayView1d<integer const> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  // make a list of region indices to be included
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  constexpr localIndex numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight, faceWeightInv;

  stackArray1d<localIndex, numElems> stencilCellsRegionIndex(numElems);
  stackArray1d<localIndex, numElems> stencilCellsSubRegionIndex(numElems);
  stackArray1d<localIndex, numElems> stencilCellsIndex(numElems);
  stackArray1d<real64, numElems> stencilWeights(numElems);

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve(faceManager->size() );

  real64 const lengthTolerance = meshBody->getGlobalLengthScale() * this->m_areaRelTol;
  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    if (faceGhostRank[kf] >= 0 || elemRegionList[kf][0] == -1 || elemRegionList[kf][1] == -1)
      continue;

    if ( !(regionFilter.contains(elemRegionList[kf][0]) && regionFilter.contains(elemRegionList[kf][1])) )
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );

    if( faceArea < areaTolerance )
      continue;

    faceWeightInv = 0.0;

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      if (elemRegionList[kf][ke] != -1)
      {
        localIndex const er  = elemRegionList[kf][ke];
        localIndex const esr = elemSubRegionList[kf][ke];
        localIndex const ei  = elemList[kf][ke];

        cellToFaceVec = faceCenter;
        cellToFaceVec -= elemCenter[er][esr][ei];

        // ensure normal orientation outward of first cell
        if (ke == 0 && Dot(cellToFaceVec, faceNormal) < 0.0)
        {
          faceNormal *= -1;
        }

        if (ke == 1)
        {
          cellToFaceVec *= -1.0;
        }

        real64 const c2fDistance = cellToFaceVec.Normalize();

        // assemble full coefficient tensor from principal axis/components
        makeFullTensor(coefficient[er][esr][ei], coefTensor);

        faceConormal.AijBj(coefTensor, faceNormal);
        real64 halfWeight = Dot(cellToFaceVec, faceConormal);

        // correct negative weight issue arising from non-K-orthogonal grids
        if( halfWeight < 0.0 )
        {
          faceConormal.AijBj(coefTensor, cellToFaceVec);
          halfWeight = Dot(cellToFaceVec, faceConormal);
        }

        halfWeight *= faceArea / c2fDistance;
        halfWeight = std::max( halfWeight, weightTolerance );

        faceWeightInv += 1.0 / halfWeight;
      }
    }

    GEOS_ASSERT( faceWeightInv > 0.0 );
    faceWeight = 1.0 / faceWeightInv;

    for (localIndex ke = 0; ke < numElems; ++ke)
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
                 kf);
  }
//  stencil.compress();
}


void TwoPointFluxApproximation::addToFractureStencil( DomainPartition const & domain,
                                                      string const & faceElementRegionName )
{
  MeshLevel const * const mesh = domain.getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  EdgeManager const * const edgeManager = mesh->getEdgeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  //OrderedVariableOneToManyRelation const & facesToEdgesMap = faceManager->edgeList();
  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const
  elemCenter = elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >( CellBlock::
                                                                                               viewKeyStruct::
                                                                                               elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor> > const
  coefficient = elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(m_coeffName);

  arrayView1d<real64 const>   const & faceArea   = faceManager->faceArea();
  arrayView1d<R1Tensor const> const & faceCenter = faceManager->faceCenter();
  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
  arrayView1d<R1Tensor const> const & X = nodeManager->referencePosition();

  FaceElementStencil & fractureStencil = getReference<FaceElementStencil>(viewKeyStruct::fractureStencilString);
  CellElementStencilTPFA & cellStencil = getReference<CellElementStencilTPFA>(viewKeyStruct::cellStencilString);

  FaceElementRegion const * const fractureRegion = elemManager->GetRegion<FaceElementRegion>(faceElementRegionName);
  localIndex const fractureRegionIndex = fractureRegion->getIndexInParent();

  FaceElementSubRegion const * const fractureSubRegion = fractureRegion->GetSubRegion<FaceElementSubRegion>("default");
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion->faceList();

  array1d<localIndex> const &
  fractureConnectorsToEdges = edgeManager->getReference< array1d<localIndex > >( EdgeManager::
                                                                                 viewKeyStruct::
                                                                                 fractureConnectorEdgesToEdgesString );

  ArrayOfArrays<localIndex> const &
  fractureConnectorsToFaceElements =
    edgeManager->getReference< ArrayOfArrays<localIndex> >( EdgeManager::
                                                                viewKeyStruct::
                                                                fractureConnectorsEdgesToFaceElementsIndexString );

  FixedToManyElementRelation const &
  faceElementsToCells = fractureSubRegion->m_faceElementsToCells;

  arrayView1d< real64 const > const & aperture = fractureSubRegion->getElementAperture();

  localIndex constexpr maxElems = FaceElementStencil::MAX_STENCIL_SIZE;

  stackArray1d<localIndex, maxElems> stencilCellsRegionIndex;
  stackArray1d<localIndex, maxElems> stencilCellsSubRegionIndex;
  stackArray1d<localIndex, maxElems> stencilCellsIndex;
  stackArray1d<real64, maxElems> stencilWeights;

  arrayView1d<integer const> const & edgeGhostRank = edgeManager->GhostRank();

  // add new connectors/connections between face elements to the fracture stencil
  for( auto const fci : edgeManager->m_recalculateFractureConnectorEdges )
  {
    localIndex const numElems = fractureConnectorsToFaceElements.sizeOfArray(fci);
    // only do this if there are more than one element attached to the connector
    localIndex const edgeIndex = fractureConnectorsToEdges[fci];

    if( edgeGhostRank[edgeIndex] < 0 && numElems > 1 )
    {

      GEOS_ERROR_IF(numElems > maxElems, "Max stencil size exceeded by fracture-fracture connector " << fci);
      stencilCellsRegionIndex.resize(numElems);
      stencilCellsSubRegionIndex.resize(numElems);
      stencilCellsIndex.resize(numElems);
      stencilWeights.resize(numElems);

      // get edge geometry
      R1Tensor edgeCenter, edgeLength;
      edgeManager->calculateCenter( edgeIndex, X, edgeCenter );
      edgeManager->calculateLength( edgeIndex, X, edgeLength );

      // loop over all face elements attached to the connector and add them to the stencil
      for( localIndex kfe=0 ; kfe<numElems ; ++kfe )
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

        stencilWeights[kfe] =  1.0 / 12.0 * edgeLength.L2_Norm() / cellCenterToEdgeCenter.L2_Norm();
      }
      // add/overwrite the stencil for index fci
      fractureStencil.add( numElems,
                           stencilCellsRegionIndex,
                           stencilCellsSubRegionIndex,
                           stencilCellsIndex,
                           stencilWeights.data(),
                           fci );
    }
  }

  // add connections for FaceElements to/from CellElements.
  {
    arrayView2d<localIndex const> const & elemRegionList = faceElementsToCells.m_toElementRegion;
    arrayView2d<localIndex const> const & elemSubRegionList = faceElementsToCells.m_toElementSubRegion;
    arrayView2d<localIndex const> const & elemList = faceElementsToCells.m_toElementIndex;
    for( localIndex const kfe : fractureSubRegion->m_newFaceElements )
//    for( localIndex kfe=0 ; kfe<faceElementsToCells.size(0) ; ++kfe )
    {
      if( fractureSubRegion->GhostRank()[kfe] < 0 )
      {
        localIndex const numElems = faceElementsToCells.size(1);

        GEOS_ERROR_IF(numElems > maxElems, "Max stencil size exceeded by fracture-cell connector " << kfe);
        stencilCellsRegionIndex.resize(numElems);
        stencilCellsSubRegionIndex.resize(numElems);
        stencilCellsIndex.resize(numElems);
        stencilWeights.resize(numElems);

        R2SymTensor coefTensor;
        R1Tensor cellToFaceVec;
        R1Tensor faceConormal;

        // remove cell-to-cell connections from cell stencil and add in new connections
        if( cellStencil.zero( faceMap[kfe][0] ) )
        {
          for (localIndex ke = 0; ke < numElems; ++ke)
          {
            localIndex const faceIndex = faceMap[kfe][ke];
            localIndex const er  = elemRegionList[kfe][ke];
            localIndex const esr = elemSubRegionList[kfe][ke];
            localIndex const ei  = elemList[kfe][ke];

            cellToFaceVec = faceCenter[faceIndex];
            cellToFaceVec -= elemCenter[er][esr][ei];

            real64 const c2fDistance = cellToFaceVec.Normalize();

            // assemble full coefficient tensor from principal axis/components
            makeFullTensor(coefficient[er][esr][ei], coefTensor);

            faceConormal.AijBj(coefTensor, faceNormal[faceIndex]);
            real64 const ht = Dot( cellToFaceVec, faceConormal ) * faceArea[faceIndex] / c2fDistance;

            // assume the h for the faceElement to the connector (Face) is zero. thus the weights are trivial.
            stencilCellsRegionIndex[0] = er;
            stencilCellsSubRegionIndex[0] = esr;
            stencilCellsIndex[0] = ei;
            stencilWeights[0] =  ht ;

            stencilCellsRegionIndex[1] = fractureRegionIndex;
            stencilCellsSubRegionIndex[1] = 0;
            stencilCellsIndex[1] = kfe;
            stencilWeights[1] = -ht ;

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

void TwoPointFluxApproximation::computeBoundaryStencil( DomainPartition const & domain,
                                                        set<localIndex> const & faceSet,
                                                        BoundaryStencil & stencil )
{
  MeshBody const * const meshBody = domain.getMeshBody(0);
  MeshLevel const * const mesh = meshBody->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  array2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList           = faceManager->elementList();
  r1_array const & X = nodeManager->referencePosition();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const elemCenter =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(
                                        CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const coefficient =
    elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(m_coeffName);

  integer_array const & faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::
                                                                                 viewKeyStruct::
                                                                                 ghostRankString);

  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  // make a list of region indices to be included
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  constexpr localIndex numElems = BoundaryStencil::NUM_POINT_IN_FLUX;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight;

  stackArray1d<PointDescriptor, numElems> stencilPoints(numElems);
  stackArray1d<real64, numElems> stencilWeights(numElems);

  real64 const lengthTolerance = meshBody->getGlobalLengthScale() * this->m_areaRelTol;
  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve(faceSet.size(), 2);
  for (localIndex kf : faceSet)
  {
    if (faceGhostRank[kf] >= 0)
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      if (elemRegionList[kf][ke] < 0)
        continue;

      if (!regionFilter.contains(elemRegionList[kf][ke]))
        continue;

      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];
      localIndex const ei  = elemList[kf][ke];

      cellToFaceVec = faceCenter;
      cellToFaceVec -= elemCenter[er][esr][ei];

      // ensure normal orientation outward of the cell
      if (Dot(cellToFaceVec, faceNormal) < 0.0)
      {
        faceNormal *= -1;
      }

      real64 const c2fDistance = cellToFaceVec.Normalize();

      // assemble full coefficient tensor from principal axis/components
      makeFullTensor(coefficient[er][esr][ei], coefTensor);

      faceConormal.AijBj(coefTensor, faceNormal);
      faceWeight = Dot(cellToFaceVec, faceConormal);

      // correct negative weight issue arising from non-K-orthogonal grids
      if (faceWeight < 0.0)
      {
        faceConormal.AijBj(coefTensor, cellToFaceVec);
        faceWeight = Dot(cellToFaceVec, faceConormal);
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


REGISTER_CATALOG_ENTRY(FluxApproximationBase, TwoPointFluxApproximation, std::string const &, Group * const)

}
