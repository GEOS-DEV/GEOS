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

#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;

TwoPointFluxApproximation::TwoPointFluxApproximation(std::string const &name,
                                                     ManagedGroup *const parent)
  : FluxApproximationBase(name, parent)
{

}

void TwoPointFluxApproximation::compute(DomainPartition * domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  Array2dT<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  Array2dT<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  Array2dT<localIndex> const & elemList           = faceManager->elementList();
  r1_array const & X = nodeManager->referencePosition();

  auto elemCenter = elemManager->ConstructViewAccessor< r1_array >(CellBlock::
                                                                   viewKeyStruct::
                                                                   elementCenterString);

  // TODO treat scalar/vector/tensor inputs
  auto coefficient = elemManager->ConstructViewAccessor<r1_array>(m_fieldName);

  auto elemsToNodes = elemManager->ConstructViewAccessor<FixedOneToManyRelation>(CellBlockSubRegion::
                                                                                 viewKeyStruct::
                                                                                 nodeListString );


  integer_array const & faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::
                                                                                 viewKeyStruct::
                                                                                 ghostRankString);

  array<array<localIndex>> const & faceToNodes = faceManager->nodeList();

  constexpr localIndex numElems = 2;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2Tensor coefTensor, coefTensorTemp;
  real64 faceArea, faceWeight, faceWeightInv;

  array<StencilCollection::CellDescriptor> stencilCells(numElems);
  array<real64> stencilWeights(numElems);

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  localIndex numFaceConnectors = 0;
  m_stencil.resize(faceManager->size());
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    faceArea = computationalGeometry::Centroid_3DPolygon(faceToNodes[kf], X, faceCenter, faceNormal);

    faceWeightInv = 0.0;
    localIndex numActual = 0;

    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      if (elemRegionList[kf][ke] != -1)
      {
        const localIndex er  = elemRegionList[kf][ke];
        const localIndex esr = elemSubRegionList[kf][ke];
        const localIndex ei  = elemList[kf][ke];

        cellToFaceVec = faceCenter;
        cellToFaceVec -= elemCenter[er][esr][ei];

        if (ke == 1)
          cellToFaceVec *= -1.0;

        real64 const c2fDistance = cellToFaceVec.Normalize();

        // assemble full coefficient tensor from principal axis/components
        R1Tensor & coefValues = coefficient[er][esr][ei];
        coefTensor = 0.0;

        // XXX: unroll loop manually?
        for (unsigned icoord = 0; icoord < 3; ++icoord)
        {
          // in future principal axis may be specified in input
          R1Tensor axis(0.0); axis(icoord) = 1.0;

          // XXX: is there a more elegant way to do this?
          coefTensorTemp.dyadic_aa(axis);
          coefTensorTemp *= coefValues(icoord);
          coefTensor += coefTensorTemp;
        }

        faceConormal.AijBj(coefTensor, faceNormal);
        real64 const ht = Dot(cellToFaceVec, faceConormal) * faceArea / c2fDistance;

        faceWeightInv += 1.0 / ht; // XXX: safeguard against div by zero?
        ++numActual;
      }
    }

    faceWeight = 1.0 * numActual / faceWeightInv; // XXX: safeguard against div by zero?

    // ensure consistent normal orientation
    if (Dot(cellToFaceVec, faceNormal) < 0)
      faceWeight *= -1;

    if (faceGhostRank[kf] < 0 && elemRegionList[kf][0] != -1 && elemRegionList[kf][1] != -1)
    {
      for (localIndex ke = 0; ke < numElems; ++ke)
      {
        stencilCells[ke] = { elemRegionList[kf][ke], elemSubRegionList[kf][ke], elemList[kf][ke] };
        stencilWeights[ke] = faceWeight * (ke == 0 ? 1 : -1);
      }
      m_stencil.set(numFaceConnectors, stencilCells.data(), stencilCells, stencilWeights);
      ++numFaceConnectors;
    }
  }
  m_stencil.resize(numFaceConnectors);
}

REGISTER_CATALOG_ENTRY(FluxApproximationBase, TwoPointFluxApproximation, std::string const &, ManagedGroup * const)

}