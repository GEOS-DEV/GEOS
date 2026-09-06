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
 * @file ConsistencyAdaptation.cpp
 */

#include "ConsistencyAdaptation.hpp"

#include "mesh/MeshLevel.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mixedMimetic/MixedMimeticFields.hpp"
#include "mixedMimetic/consistency/ConsistencyAdaptationKernels.hpp"

#include <algorithm>
#include <array>

namespace geos
{

using namespace fields;

ConsistencyAdaptation::Report ConsistencyAdaptation::classify( MeshLevel & mesh,
                                                               string_array const & regionNames,
                                                               SortedArrayView< localIndex const > const & regionFilter,
                                                               PermeabilityAccessor const & permeability,
                                                               Parameters const & params,
                                                               stdVector< NeighborCommunicator > & neighbors )
{
  GEOS_MARK_FUNCTION;

  Report report;
  ElementRegionManager & elemManager = mesh.getElemManager();

  // every cell starts with the consistent (MFD) product; the layers below can only refine that choice
  elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                            [&]( localIndex const,
                                                                 ElementSubRegionBase & subRegion )
  {
    subRegion.getField< mixedMimetic::mfdFlag >().template setValues< parallelDevicePolicy<> >( 1 );
    report.numCells += subRegion.size() - subRegion.getNumberOfGhosts();
  } );

  report.numConsistent = params.adaptiveConsistency
                         ? applyConsistencyLayer( mesh, regionNames, regionFilter, permeability, params, neighbors )
                         : report.numCells;

  std::pair< localIndex, localIndex > const prescribed = applyPrescription( mesh, regionNames, neighbors );
  report.numPrescribed0 = prescribed.first;
  report.numPrescribed1 = prescribed.second;

  std::pair< localIndex, localIndex > const degenerate = applyDegeneracyLayer( mesh, regionNames, params.degeneracyTolerance, neighbors );
  report.numDegenerate = degenerate.first;
  report.numRejected = degenerate.second;

  labelFaces( mesh, regionFilter, params.effectiveTpfa );
  return report;
}

localIndex ConsistencyAdaptation::applyConsistencyLayer( MeshLevel & mesh,
                                                         string_array const & regionNames,
                                                         SortedArrayView< localIndex const > const & regionFilter,
                                                         PermeabilityAccessor const & permeability,
                                                         Parameters const & params,
                                                         stdVector< NeighborCommunicator > & neighbors )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();
  real64 const (&gradient)[3] = params.nominalGradient;

  arrayView1d< real64 > const faceResidual = faceManager.getField< mixedMimetic::faceResidual >();
  faceResidual.zero();

  // step 1: projection of the admissible flow field induced by the nominal gradient
  array1d< real64 > projFaceFluxArray( faceManager.size() );
  arrayView1d< real64 > const projFaceFlux = projFaceFluxArray.toView();

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenterAccessor =
    elemManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

  mixedMimeticKernels::FaceFluxProjectionKernel::
    launch< parallelDevicePolicy<> >( faceManager.size(),
                                      faceManager.elementRegionList(),
                                      faceManager.elementSubRegionList(),
                                      faceManager.elementList(),
                                      regionFilter,
                                      faceManager.faceCenter(),
                                      faceManager.faceNormal(),
                                      faceManager.faceArea(),
                                      elemCenterAccessor.toNestedViewConst(),
                                      permeability,
                                      gradient,
                                      params.lengthTolerance,
                                      projFaceFlux );

  // steps 2-3: localized normalized residuals, assembled on the global face orientation
  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( regionNames,
                                                                    [&]( localIndex const,
                                                                         localIndex const er,
                                                                         localIndex const esr,
                                                                         ElementRegionBase &,
                                                                         CellElementSubRegion const & subRegion )
  {
    mixedMimeticKernels::internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NUM_FACES )
    {
      mixedMimeticKernels::LocalResidualKernel< NUM_FACES >::
      template launch< parallelDevicePolicy<> >( subRegion.size(),
                                                 nodeManager.referencePosition(),
                                                 faceManager.nodeList().toViewConst(),
                                                 subRegion.faceList().toViewConst(),
                                                 subRegion.getElementCenter(),
                                                 subRegion.getElementVolume(),
                                                 permeability[er][esr],
                                                 faceManager.faceCenter(),
                                                 faceManager.faceNormal(),
                                                 projFaceFlux.toViewConst(),
                                                 gradient,
                                                 params.lengthTolerance,
                                                 faceResidual );
    } );
  } );

  // step 4: thresholding
  localIndex numConsistent = 0;
  elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 CellElementSubRegion & subRegion )
  {
    arrayView1d< real64 > const consistencyIndicator = subRegion.getField< mixedMimetic::consistencyIndicator >();
    arrayView1d< integer > const mfdFlag = subRegion.getField< mixedMimetic::mfdFlag >();

    mixedMimeticKernels::internal::kernelLaunchSelectorFaceSwitch( subRegion.numFacesPerElement(), [&] ( auto NUM_FACES )
    {
      numConsistent += mixedMimeticKernels::MarkingKernel< NUM_FACES >::
                       template launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                  subRegion.faceList().toViewConst(),
                                                                  subRegion.ghostRank(),
                                                                  faceResidual.toViewConst(),
                                                                  params.consistencyTolerance,
                                                                  consistencyIndicator,
                                                                  mfdFlag );
    } );
  } );

  FieldIdentifiers fieldsToBeSync;
  fieldsToBeSync.addElementFields( { mixedMimetic::mfdFlag::key(), mixedMimetic::consistencyIndicator::key() }, regionNames );
  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, neighbors, false );
  return numConsistent;
}

std::pair< localIndex, localIndex > ConsistencyAdaptation::applyPrescription( MeshLevel & mesh,
                                                                              string_array const & regionNames,
                                                                              stdVector< NeighborCommunicator > & neighbors )
{
  FieldIdentifiers prescribedToBeSync;
  prescribedToBeSync.addElementFields( { mixedMimetic::prescribedMfdFlag::key() }, regionNames );
  CommunicationTools::getInstance().synchronizeFields( prescribedToBeSync, mesh, neighbors, false );

  localIndex numPrescribed0 = 0;
  localIndex numPrescribed1 = 0;
  mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                      [&]( localIndex const,
                                                                           ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const prescribed = subRegion.getField< mixedMimetic::prescribedMfdFlag >();
    arrayView1d< integer > const mfdFlag = subRegion.getField< mixedMimetic::mfdFlag >();
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      if( prescribed[ei] < 0.0 )
      {
        continue;
      }
      integer const eta = prescribed[ei] > 0.5 ? 1 : 0;
      mfdFlag[ei] = eta;
      if( ghostRank[ei] < 0 )
      {
        numPrescribed0 += ( eta == 0 );
        numPrescribed1 += ( eta == 1 );
      }
    }
  } );
  return { numPrescribed0, numPrescribed1 };
}

std::pair< localIndex, localIndex > ConsistencyAdaptation::applyDegeneracyLayer( MeshLevel & mesh,
                                                                                 string_array const & regionNames,
                                                                                 real64 const tolerance,
                                                                                 stdVector< NeighborCommunicator > & neighbors )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager & elemManager = mesh.getElemManager();
  ArrayOfArraysView< localIndex const > const nodeToRegion = nodeManager.elementRegionList();
  ArrayOfArraysView< localIndex const > const nodeToSubRegion = nodeManager.elementSubRegionList();
  ArrayOfArraysView< localIndex const > const nodeToElem = nodeManager.elementList();
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager.constructArrayViewAccessor< real64, 1 >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

  localIndex numDegenerate = 0;
  localIndex numRejected = 0;
  elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                            [&]( localIndex const,
                                                                 CellElementSubRegion & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes = subRegion.nodeList().toViewConst();
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
    arrayView1d< real64 > const indicator = subRegion.getField< mixedMimetic::degeneracyIndicator >();
    arrayView1d< integer > const mfdFlag = subRegion.getField< mixedMimetic::mfdFlag >();
    arrayView1d< real64 const > const prescribed = subRegion.getField< mixedMimetic::prescribedMfdFlag >();
    localIndex const numNodes = subRegion.numNodesPerElement();

    // the node star of a cell: every cell sharing a vertex with it, counted once
    stdVector< std::array< localIndex, 3 > > star;
    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {
      star.clear();
      for( localIndex a = 0; a < numNodes; ++a )
      {
        localIndex const n = elemToNodes( ei, a );
        for( localIndex k = 0; k < nodeToElem.sizeOfArray( n ); ++k )
        {
          if( nodeToRegion( n, k ) >= 0 && nodeToSubRegion( n, k ) >= 0 && nodeToElem( n, k ) >= 0 )
          {
            star.push_back( { nodeToRegion( n, k ), nodeToSubRegion( n, k ), nodeToElem( n, k ) } );
          }
        }
      }
      std::sort( star.begin(), star.end() );
      star.erase( std::unique( star.begin(), star.end() ), star.end() );
      real64 sum = 0.0;
      for( auto const & c : star )
      {
        sum += elemVolume[c[0]][c[1]][c[2]];
      }
      // share of the cell in the volume of its node star, in percent
      real64 const percent = sum > 0.0 ? 100.0 * volume[ei] / sum : 0.0;
      indicator[ei] = percent;
      if( percent < tolerance && mfdFlag[ei] == 1 )
      {
        mfdFlag[ei] = 0;
        numDegenerate += ( ghostRank[ei] < 0 ) ? 1 : 0;
        numRejected += ( ghostRank[ei] < 0 && prescribed[ei] > 0.5 ) ? 1 : 0;
      }
    }
  } );

  FieldIdentifiers fieldsToBeSync;
  fieldsToBeSync.addElementFields( { mixedMimetic::mfdFlag::key(), mixedMimetic::degeneracyIndicator::key() }, regionNames );
  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, neighbors, false );
  return { numDegenerate, numRejected };
}

void ConsistencyAdaptation::labelFaces( MeshLevel & mesh,
                                        SortedArrayView< localIndex const > const & regionFilter,
                                        bool const effectiveTpfa )
{
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const mfdFlagAccessor =
    elemManager.constructArrayViewAccessor< integer, 1 >( mixedMimetic::mfdFlag::key() );

  mixedMimeticKernels::FaceLabelKernel::
    launch< parallelDevicePolicy<> >( faceManager.size(),
                                      faceManager.elementRegionList(),
                                      faceManager.elementSubRegionList(),
                                      faceManager.elementList(),
                                      regionFilter,
                                      mfdFlagAccessor.toNestedViewConst(),
                                      effectiveTpfa,
                                      faceManager.getField< mixedMimetic::faceStencilLabel >() );
}

} // namespace geos
