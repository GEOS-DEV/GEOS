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

#include "meshUtilities/ComputationalGeometry.hpp"
#include "mesh/AggregateElementSubRegion.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"

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

void TwoPointFluxApproximation::computeMainStencil(DomainPartition * domain, CellStencil & stencil)
{
  MeshBody * const meshBody = domain->getMeshBodies()->GetGroup<MeshBody>(0);
  MeshLevel * const mesh = meshBody->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

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

  constexpr localIndex numElems = 2;

  R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
  R2SymTensor coefTensor;
  real64 faceArea, faceWeight, faceWeightInv;

  array1d<CellDescriptor> stencilCells(numElems);
  array1d<real64> stencilWeights(numElems);

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  stencil.reserve(faceManager->size(), 2);
  real64 const areaTolerance = pow( meshBody->getGlobalLengthScale() * this->m_areaRelTol, 2 );
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    if (faceGhostRank[kf] >= 0 || elemRegionList[kf][0] == -1 || elemRegionList[kf][1] == -1)
      continue;

    faceArea = computationalGeometry::Centroid_3DPolygon(faceToNodes[kf], X, faceCenter, faceNormal, areaTolerance );

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
    stencil.add(stencilCells.data(), stencilCells, stencilWeights, kf);
  }
  stencil.compress();
}

void TwoPointFluxApproximation::computeCoarsetencil( DomainPartition * domain,
                                                     CellStencil const & fineStencil,
                                                     CellStencil & coarseStencil,
                                                     std::string const & elementaryPressure1Name,
                                                     std::string const & elementaryPressure2Name,
                                                     std::string const & elementaryPressure3Name )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ElementRegion * const elemRegion = elemManager->GetRegion(0); // TODO : still one region / elemsubregion
  AggregateElementSubRegion * const aggregateElement = elemRegion->GetSubRegion("coarse")->group_cast< AggregateElementSubRegion * >();
  const array1d< localIndex > & fineToCoarse = aggregateElement->GetFineToCoarseMap();
  // TODO : will it work for fractures ? (no)

  array1d<CellDescriptor> stencilCells(2);
  array1d<real64> stencilWeights(2);
  std::set< std::pair< localIndex, localIndex > > interfaces;
  fineStencil.forAll( [&] ( StencilCollection<CellDescriptor, real64>::Accessor stencil ) //TODO maybe find a clever way to iterate between coarse interfaces ?
  {
    localIndex const stencilSize = stencil.size();
    if( stencil.size() == 2)
    {
      CellDescriptor const & cell1 = stencil.connectedIndex(0);
      CellDescriptor const & cell2 = stencil.connectedIndex(1);
      localIndex aggregateNumber1 = fineToCoarse[cell1.index];
      localIndex aggregateNumber2 = fineToCoarse[cell2.index];
      if( aggregateNumber1 != aggregateNumber2 && interfaces.find(std::make_pair(aggregateNumber1,aggregateNumber2)) == interfaces.end() ) // We find two adjacent aggregates TODO maybe find a clever way to iterate between coarse interfaces ?
      {
        stencilCells[0] = { cell1.region, 0, aggregateNumber1 }; // We can consider here that there is only 1 subregion
        stencilCells[1] = { cell2.region, 0, aggregateNumber2 };

        // Now we compute the transmissibilities
        R1Tensor barycenter1 = aggregateElement->getElementCenter()[aggregateNumber1];
        R1Tensor barycenter2 = aggregateElement->getElementCenter()[aggregateNumber2];
        barycenter1 -= barycenter2; // normal between the two aggregates
        barycenter1.Normalize();

        int systemSize = integer_conversion< int >(aggregateElement->GetNbCellsPerAggregate( aggregateNumber1 )
                                                 + aggregateElement->GetNbCellsPerAggregate( aggregateNumber2 ));
        Teuchos::LAPACK< int, real64 > lapack;
        Teuchos::SerialDenseMatrix< int, real64 > A(systemSize, 4);
        Teuchos::SerialDenseVector< int, real64 > pTarget(systemSize);

        /// Get the elementary pressures
        ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure1 =
          elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure1Name );
        ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure2 =
          elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure2Name );
        ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> pressure3 =
          elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( elementaryPressure3Name );
        int count =0;

        /// Setup the least square system
        aggregateElement->forFineCellsInAggregate( aggregateNumber1,
                                                   [&] ( localIndex fineCellIndex )
        {
          A(count,0) = pressure1[cell1.region][cell1.subRegion][fineCellIndex];
          A(count,1) = pressure2[cell1.region][cell1.subRegion][fineCellIndex];
          A(count,2) = pressure3[cell1.region][cell1.subRegion][fineCellIndex];
          A(count,3) = 1.;
          R1Tensor barycenterFineCell = elemRegion->GetSubRegion(cell1.subRegion)->getElementCenter()[fineCellIndex];
          pTarget(count++) = barycenterFineCell[0]*barycenter1[0]
                             + barycenterFineCell[1]*barycenter1[1]
                             + barycenterFineCell[2]*barycenter1[2];
        });
        aggregateElement->forFineCellsInAggregate( aggregateNumber2,
                                                   [&] ( localIndex fineCellIndex )
        {
          A(count,0) = pressure1[cell2.region][cell2.subRegion][fineCellIndex];
          A(count,1) = pressure2[cell2.region][cell2.subRegion][fineCellIndex];
          A(count,2) = pressure3[cell2.region][cell2.subRegion][fineCellIndex];
          A(count,3) = 1.;
          R1Tensor barycenterFineCell = elemRegion->GetSubRegion(cell1.subRegion)->getElementCenter()[fineCellIndex];
          pTarget(count++) = barycenterFineCell[0]*barycenter1[0]
                             + barycenterFineCell[1]*barycenter1[1]
                             + barycenterFineCell[2]*barycenter1[2];
        });

        int info;
        real64  rwork1;
        real64 svd[4];
        int rank;

        // Solve the least square system
        lapack.GELSS(systemSize,4,1,A.values(),A.stride(),pTarget.values(),pTarget.stride(),svd,-1,&rank,&rwork1,-1,&info);
        int lwork = static_cast< int > ( rwork1 );
        real64 * rwork = new real64[lwork];
        lapack.GELSS(systemSize,4,1,A.values(),A.stride(),pTarget.values(),pTarget.stride(),svd,-1,&rank,rwork,lwork,&info);
        pTarget.print(std::cout);
        GEOS_ERROR("stop");

        // Computation of coarse-grid flow parameters
        real64 coarseAveragePressure1 = 0.;
        real64 coarseAveragePressure2 = 0.;
        array1d<real64> coarseFlowRate(2);

        aggregateElement->forFineCellsInAggregate( aggregateNumber1,
                                                   [&] ( localIndex fineCellIndex )
        {
          coarseAveragePressure1 +=  pTarget[0] * pressure1[cell1.region][cell1.subRegion][fineCellIndex]
                                   + pTarget[1] * pressure2[cell1.region][cell1.subRegion][fineCellIndex]
                                   + pTarget[2] * pressure3[cell1.region][cell1.subRegion][fineCellIndex]
                                   + pTarget[3];
          coarseAveragePressure1 *= elemRegion->GetSubRegion(cell1.subRegion)->getElementVolume()[fineCellIndex];
        });
        aggregateElement->forFineCellsInAggregate( aggregateNumber2,
                                                   [&] ( localIndex fineCellIndex )
        {
          coarseAveragePressure2 +=  pTarget[0] * pressure1[cell2.region][cell2.subRegion][fineCellIndex]
                                   + pTarget[1] * pressure2[cell2.region][cell2.subRegion][fineCellIndex]
                                   + pTarget[2] * pressure3[cell2.region][cell2.subRegion][fineCellIndex]
                                   + pTarget[3];
          coarseAveragePressure2 *= elemRegion->GetSubRegion(cell2.subRegion)->getElementVolume()[fineCellIndex];
        });

        coarseAveragePressure1 /= aggregateElement->getElementVolume()[aggregateNumber1];
        coarseAveragePressure2 /= aggregateElement->getElementVolume()[aggregateNumber2];

        fineStencil.forAll( [&] ( StencilCollection<CellDescriptor, real64>::Accessor stencilBis )
        {
          localIndex const stencilSizeBis = stencilBis.size();
          if( stencilBis.size() == 2)
          {
            CellDescriptor const & cell1Bis = stencilBis.connectedIndex(0);
            CellDescriptor const & cell2Bis = stencilBis.connectedIndex(1);
            localIndex aggregateNumber1Bis = fineToCoarse[cell1Bis.index];
            localIndex aggregateNumber2Bis = fineToCoarse[cell2Bis.index];
            if( aggregateNumber1Bis == aggregateNumber1 &&
                aggregateNumber2Bis == aggregateNumber2 ) // TODO clever way to find the interfaces of adjacent fine cells within aggregates?
            {
              real64 finePressure1 =   pTarget[0] * pressure1[cell1.region][cell1.subRegion][cell1Bis.index]
                                     + pTarget[1] * pressure2[cell1.region][cell1.subRegion][cell1Bis.index]
                                     + pTarget[2] * pressure3[cell1.region][cell1.subRegion][cell1Bis.index]
                                     + pTarget[3];

              real64 finePressure2 =   pTarget[0] * pressure1[cell2.region][cell2.subRegion][cell2Bis.index]
                                     + pTarget[1] * pressure2[cell2.region][cell2.subRegion][cell2Bis.index]
                                     + pTarget[2] * pressure3[cell2.region][cell2.subRegion][cell2Bis.index]
                                     + pTarget[3];
              stencilBis.forAll([&] (CellDescriptor const & cell, real64 w, localIndex const i) -> void
              {
                  coarseFlowRate[i] += w * ( finePressure1 - finePressure2 ); //TODO sign ?
              });
            }
          }
        });
        for( localIndex i = 0; i < 2; i++ )
        {
          stencilWeights[i] = coarseFlowRate[i] / ( coarseAveragePressure1 - coarseAveragePressure2 ); // TODO sign ?
        }
        coarseStencil.add(stencilCells.data(), stencilCells, stencilWeights, 0.);
        interfaces.insert(std::make_pair(aggregateNumber1, aggregateNumber2));
      }
    }
});
}


void TwoPointFluxApproximation::computeFractureStencil( DomainPartition const & domain,
                                                        CellStencil & fractureStencil,
                                                        CellStencil & cellStencil )
{

  MeshLevel const * const mesh = domain.getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  EdgeManager const * const edgeManager = mesh->getEdgeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  localIndex fractureRegionIndex = elemManager->GetRegions().getIndex(m_fractureRegionName);
  ElementRegion const * const fractureRegion = elemManager->GetRegion( m_fractureRegionName );


  FaceElementSubRegion const * const fractureSubRegion = fractureRegion->GetSubRegion<FaceElementSubRegion>(0);


  FaceElementSubRegion::NodeMapType const & nodeMap = fractureSubRegion->nodeList();
  FaceElementSubRegion::EdgeMapType const & edgeMap = fractureSubRegion->edgeList();
  FaceElementSubRegion::FaceMapType const & faceMap = fractureSubRegion->faceList();

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const
  elemCenter = elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(CellBlock::viewKeyStruct::elementCenterString);

  ElementRegionManager::ElementViewAccessor<arrayView1d<R1Tensor>> const
  coefficient = elemManager->ConstructViewAccessor< array1d<R1Tensor>, arrayView1d<R1Tensor> >(m_coeffName);


  array1d<localIndex> const &
  fractureConnectorIndices = fractureRegion->getReference< array1d<localIndex > >( ElementRegion::viewKeyStruct::fractureConnectorIndicesString );

  array1d<array1d<localIndex> > const &
  fractureConnectors = fractureRegion->getReference< array1d<array1d<localIndex> > >( ElementRegion::viewKeyStruct::fractureElementConnectorString );

  FixedToManyElementRelation const &
  fractureCellConnectors = fractureRegion->getReference< FixedToManyElementRelation >( ElementRegion::viewKeyStruct::fractureToCellConnectorString );

  arrayView1d<real64 const> const & faceArea = faceManager->faceArea();
  arrayView1d<R1Tensor const> const & faceCenter = faceManager->faceCenter();
  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();

  arrayView1d<R1Tensor const> const & X = nodeManager->referencePosition();

  arrayView1d< real64 const > const & aperture = fractureSubRegion->getElementAperture();

  // connections between FaceElements
  for( localIndex fci=0 ; fci<fractureConnectors.size() ; ++fci )
  {
    localIndex const numElems = fractureConnectors[fci].size();
    localIndex const edgeIndex = fractureConnectorIndices[fci];

    array1d<CellDescriptor> stencilCells(numElems);
    array1d<real64> stencilWeights(numElems);

    R1Tensor edgeCenter, edgeLength;
    edgeManager->calculateCenter( edgeIndex, X, edgeCenter );
    edgeManager->calculateLength( edgeIndex, X, edgeLength );

    for( localIndex kfe=0 ; kfe<numElems ; ++kfe )
    {
      localIndex const fractureElementIndex = fractureConnectors[fci][kfe];
      R1Tensor cellCenterToEdgeCenter = edgeCenter;
      cellCenterToEdgeCenter -= faceCenter[ faceMap[fractureElementIndex][0] ];
      stencilCells[kfe] = { fractureRegionIndex, 0, fractureElementIndex };
      // TODO stenciWeights will mean something else once you take out the aperture.
      // We won't be doing the harmonic mean here...etc.
      stencilWeights[kfe] = pow( -1 , kfe ) * pow( aperture[fractureElementIndex], 3) / 12.0 * edgeLength.L2_Norm() / cellCenterToEdgeCenter.L2_Norm();
    }
    fractureStencil.add( stencilCells.data(), stencilCells, stencilWeights, edgeIndex );
  }

  // add connections for FaceElements to/from CellElements.
  {
    array2d< CellDescriptor > cellStencilZeros( fractureCellConnectors.size(0), 2 );

    arrayView2d<localIndex const> const & elemRegionList = fractureCellConnectors.m_toElementRegion;
    arrayView2d<localIndex const> const & elemSubRegionList = fractureCellConnectors.m_toElementSubRegion;
    arrayView2d<localIndex const> const & elemList = fractureCellConnectors.m_toElementIndex;
    for( localIndex kfe=0 ; kfe<fractureCellConnectors.size(0) ; ++kfe )
    {
      localIndex const numElems = fractureCellConnectors.size(1);

      array1d<CellDescriptor> stencilCells(numElems);
      array1d<real64> stencilWeights(numElems);

      R2SymTensor coefTensor;
      R1Tensor cellToFaceVec;
      R1Tensor faceConormal;

      for (localIndex ke = 0; ke < numElems; ++ke)
      {
        real64 faceWeightInv = 0.0;

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

        stencilCells[0] = { er, esr, ei};
        stencilWeights[0] = pow(-1,ke) * ht ;

        stencilCells[1] = { fractureRegionIndex, 0, kfe};
        stencilWeights[1] = -pow(-1,ke) * ht ;

        fractureStencil.add( stencilCells.data(), stencilCells, stencilWeights, faceIndex );
      }


      // remove cell connectors from original stencil
      cellStencilZeros[kfe][0] = { elemRegionList[kfe][0], elemSubRegionList[kfe][0], elemList[kfe][0] };
      cellStencilZeros[kfe][1] = { elemRegionList[kfe][1], elemSubRegionList[kfe][1], elemList[kfe][1] };

      cellStencil.zero( faceMap[kfe][0], cellStencilZeros.data() );
    }
  }
  fractureStencil.compress();

}


void TwoPointFluxApproximation::computeBoundaryStencil(DomainPartition * domain, set<localIndex> const & faceSet,
                                                       FluxApproximationBase::BoundaryStencil & stencil)
{
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
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

        stencil.add(stencilPoints.data(), stencilPoints, stencilWeights, kf );
      }
    }
  }
  stencil.compress();
}


REGISTER_CATALOG_ENTRY(FluxApproximationBase, TwoPointFluxApproximation, std::string const &, ManagedGroup * const)

}
