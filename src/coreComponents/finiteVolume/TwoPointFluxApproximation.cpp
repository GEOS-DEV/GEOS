/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * @file TwoPointFluxApproximation.cpp
 *
 */
#include "TwoPointFluxApproximation.hpp"

#include "managers/Wells/WellBase.hpp"
#include "managers/Wells/PerforationManager.hpp"
#include "managers/Wells/Perforation.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;

TwoPointFluxApproximation::TwoPointFluxApproximation(std::string const &name,
                                                     ManagedGroup *const parent)
  : FluxApproximationBase(name, parent)
{

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

void TwoPointFluxApproximation::computeMainStencil( DomainPartition const * domain,
                                                    CellStencil & stencil ) const
{
  MeshLevel const * mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager const * nodeManager = mesh->getNodeManager();
  FaceManager const * faceManager = mesh->getFaceManager();
  ElementRegionManager const * elemManager = mesh->getElemManager();

  array2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList           = faceManager->elementList();
  r1_array const & X = nodeManager->referencePosition();

  auto elemCenter = elemManager->ConstructViewAccessor< r1_array >(CellBlock::
                                                                   viewKeyStruct::
                                                                   elementCenterString);

  auto coefficient = elemManager->ConstructViewAccessor<r1_array>(m_coeffName);

  array1d<integer> const & faceGhostRank = faceManager->getReference<array1d<integer>>(ObjectManagerBase::
                                                                                       viewKeyStruct::
                                                                                       ghostRankString);

  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  constexpr localIndex numElems = 2;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight, faceWeightInv;

  array1d<CellDescriptor> stencilCells(numElems);
  array1d<real64> stencilWeights(numElems);

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve(faceManager->size(), 2);
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    if (faceGhostRank[kf] >= 0 || elemRegionList[kf][0] == -1 || elemRegionList[kf][1] == -1)
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon(faceToNodes[kf], X, faceCenter, faceNormal);

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

        if (ke == 1)
          cellToFaceVec *= -1.0;

        real64 const c2fDistance = cellToFaceVec.Normalize();

        // assemble full coefficient tensor from principal axis/components
        makeFullTensor(coefficient[er][esr][ei], coefTensor);

        faceConormal.AijBj(coefTensor, faceNormal);
        real64 const ht = Dot(cellToFaceVec, faceConormal) * faceArea / c2fDistance;

        faceWeightInv += 1.0 / ht; // XXX: safeguard against div by zero?
      }
    }

    faceWeight = 1.0 / faceWeightInv; // XXX: safeguard against div by zero?

    // ensure consistent normal orientation
    if (Dot(cellToFaceVec, faceNormal) < 0)
      faceWeight *= -1;

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      stencilCells[ke] = { elemRegionList[kf][ke], elemSubRegionList[kf][ke], elemList[kf][ke] };
      stencilWeights[ke] = faceWeight * (ke == 0 ? 1 : -1);
    }
    stencil.add(stencilCells.data(), stencilCells, stencilWeights);
  }
  stencil.compress();
}

void TwoPointFluxApproximation::computeBoundaryStencil( DomainPartition const * domain,
                                                        set<localIndex> const & faceSet,
                                                        FluxApproximationBase::BoundaryStencil & stencil ) const
{
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  array2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList           = faceManager->elementList();
  r1_array const & X = nodeManager->referencePosition();

  auto elemCenter = elemManager->ConstructViewAccessor< r1_array >(CellBlock::
                                                                   viewKeyStruct::
                                                                   elementCenterString);

  auto coefficient = elemManager->ConstructViewAccessor<r1_array>(m_coeffName);

  array1d<integer> const & faceGhostRank = faceManager->getReference<array1d<integer>>(ObjectManagerBase::
                                                                                       viewKeyStruct::
                                                                                       ghostRankString);

  array1d<array1d<localIndex>> const & faceToNodes = faceManager->nodeList();

  constexpr localIndex numElems = 2;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight;

  array1d<PointDescriptor> stencilPoints(numElems);
  array1d<real64> stencilWeights(numElems);

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve(faceSet.size(), 2);
  for (localIndex kf : faceSet)
  {
    if (faceGhostRank[kf] >= 0)
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon(faceToNodes[kf], X, faceCenter, faceNormal);

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      if (elemRegionList[kf][ke] != -1)
      {
        localIndex const er  = elemRegionList[kf][ke];
        localIndex const esr = elemSubRegionList[kf][ke];
        localIndex const ei  = elemList[kf][ke];

        cellToFaceVec = faceCenter;
        cellToFaceVec -= elemCenter[er][esr][ei];

        real64 const c2fDistance = cellToFaceVec.Normalize();

        // assemble full coefficient tensor from principal axis/components
        makeFullTensor(coefficient[er][esr][ei], coefTensor);

        faceConormal.AijBj(coefTensor, faceNormal);
        faceWeight = Dot(cellToFaceVec, faceConormal) * faceArea / c2fDistance;

        // ensure consistent normal orientation
        if (Dot(cellToFaceVec, faceNormal) < 0)
          faceWeight *= -1;

        stencilPoints[0].tag = PointDescriptor::Tag::CELL;
        stencilPoints[0].cellIndex = { er, esr, ei };
        stencilWeights[0] = faceWeight;

        stencilPoints[1].tag = PointDescriptor::Tag::FACE;
        stencilPoints[1].faceIndex = kf;
        stencilWeights[1] = -faceWeight;

        stencil.add(stencilPoints.data(), stencilPoints, stencilWeights);
      }
    }
  }
  stencil.compress();
}

void TwoPointFluxApproximation::computeWellStencil( DomainPartition const * domain,
                                                    WellBase const * well,
                                                    WellStencil & stencil ) const
{
  PerforationManager const * perfManager = well->getPerforations();

  auto const & elemRegion    = perfManager->getReference<array1d<localIndex>>( perfManager->viewKeysPerfManager.connectionElementRegion );
  auto const & elemSubregion = perfManager->getReference<array1d<localIndex>>( perfManager->viewKeysPerfManager.connectionElementSubregion );
  auto const & elemIndex     = perfManager->getReference<array1d<localIndex>>( perfManager->viewKeysPerfManager.connectionElementIndex );
  auto const & perfIndex     = perfManager->getReference<array1d<localIndex>>( perfManager->viewKeysPerfManager.connectionPerforationIndex );

  array1d<PointDescriptor> points( 2 );
  array1d<real64> weights( 2 );

  for (localIndex iconn = 0; iconn < well->numConnectionsLocal(); ++iconn)
  {
    Perforation const * perf = perfManager->getPerforation( perfIndex[iconn] );
    real64 trans = perf->getTransmissibility();

    // if transmissibility is default (i.e. not input), compute it
    if (trans < 0.0)
    {
      // TODO use Peaceman or other formula to compute well index
      trans = 0.0;

      // Should we update the input node value? (e.g. to be written into output files)
      //perf->setTransmissibility(trans);
    }

    points[0].tag = PointDescriptor::Tag::CELL;
    points[0].cellIndex = { elemRegion[iconn], elemSubregion[iconn], elemIndex[iconn] };
    weights[0] = trans;

    points[1].tag = PointDescriptor::Tag::PERF;
    points[1].perfIndex = iconn;
    weights[1] = -trans;

    stencil.add(points.data(), points, weights);
  };

  stencil.compress();
}


REGISTER_CATALOG_ENTRY(FluxApproximationBase, TwoPointFluxApproximation, std::string const &, ManagedGroup * const)

}
