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
 * @file AggregateElementRegion.cpp
 */

#include "AggregateElementRegion.hpp"
#include "AggregateElementSubRegion.hpp"

#include "cxx-utilities/src/src/SparsityPattern.hpp"

#include "metis.h"

namespace geosx
{

using namespace dataRepository;
  
AggregateElementRegion::AggregateElementRegion( string const & name, ManagedGroup * const parent ):
  ElementRegion( name, parent )
{
  this->GetGroup( viewKeyStruct::elementSubRegions )->RegisterGroup<AggregateElementSubRegion>("default");
}

AggregateElementRegion::~AggregateElementRegion()
{
}

void AggregateElementRegion::Generate( FaceManager const * const faceManager,
                                       ElementRegion const * const elementRegion,
                                       real64 coarseningRatio )
{
  AggregateElementSubRegion * const aggregateSubRegion =
    this->GetGroup( viewKeyStruct::elementSubRegions )->GetGroup<AggregateElementSubRegion>("default");
  localIndex regionIndex = elementRegion->getIndexInParent();

  array2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList           = faceManager->elementList();

  // Counting the total number of cell and number of vertices  
  localIndex nbCellElements = 0;
  elementRegion->forElementSubRegions<CellElementSubRegion >( [&]( auto * const elementSubRegion ) -> void
    {
      nbCellElements += elementSubRegion->size();
    });
  // Number of aggregate computation
  localIndex nbAggregates = static_cast< localIndex >( nbCellElements * coarseningRatio );

  // METIS variable declarations
  using idx_t = ::idx_t;
  idx_t options[METIS_NOPTIONS];                                    // Contains the METIS options
  METIS_SetDefaultOptions(options);                                 // ... That are set by default
  idx_t nnodes = integer_conversion< idx_t >( nbCellElements );     // Number of connectivity graph nodes
  idx_t nconst = 1;                                                 // Number of balancy constraints
  idx_t objval;                                                     // Total communication volume
  array1d< idx_t > parts(nnodes);                                   // Map element index -> aggregate index
  idx_t nparts = integer_conversion< idx_t >( nbAggregates );       // Number of aggregates to be generated
  

  // Compute the connectivity graph
  LvArray::SparsityPattern< idx_t, idx_t > graph( integer_conversion< idx_t >( nbCellElements ),
                                                  integer_conversion< idx_t >( nbCellElements ) );
  localIndex nbConnections = 0;
  array1d< localIndex > offsetSubRegions( elementRegion->GetSubRegions().size() );
  for( localIndex subRegionIndex = 1; subRegionIndex < offsetSubRegions.size(); subRegionIndex++ )
  {
    offsetSubRegions[subRegionIndex] = offsetSubRegions[subRegionIndex - 1] + elementRegion->GetSubRegion(subRegionIndex)->size();
  }
  for (localIndex kf = 0; kf < faceManager->size(); ++kf)
  {
    if( elemRegionList[kf][0] == regionIndex && elemRegionList[kf][1] == regionIndex )
    {
      localIndex const esr0 = elemSubRegionList[kf][0];
      idx_t const ei0  = integer_conversion< idx_t >( elemList[kf][0] + offsetSubRegions[esr0] );
      localIndex const esr1 = elemSubRegionList[kf][1];
      idx_t const ei1  = integer_conversion< idx_t >( elemList[kf][1] + offsetSubRegions[esr1] );
      graph.insertNonZero(ei0, ei1);
      graph.insertNonZero(ei1, ei0);
      nbConnections++;
    }
  }
  graph.compress();

  // METIS partitionning
  idx_t * offsets = const_cast< idx_t* >( graph.getOffsets() );
  idx_t * columns = const_cast< idx_t* >( &graph.getColumns(0)[0] );
  METIS_PartGraphRecursive( &nnodes, &nconst, offsets, columns, nullptr, nullptr, nullptr,
                            &nparts, nullptr, nullptr, options, &objval, parts.data() );

  // Compute Aggregate barycenters
  array1d< R1Tensor > aggregateBarycenters( nparts );
  array1d< real64 > aggregateVolumes( nparts );
  array1d< real64 > normalizeVolumes( nbCellElements );
  
  // First, compute the volume of each aggregates
  elementRegion->forElementSubRegions( [&]( auto * const elementSubRegion ) -> void
  {
    localIndex const subRegionIndex = elementSubRegion->getIndexInParent();
    for(localIndex cellIndex = 0; cellIndex< elementSubRegion->size() ; cellIndex++)
    {
      if( elementSubRegion->GhostRank()[cellIndex] >= 0 )
        continue;
      aggregateVolumes[parts[cellIndex + offsetSubRegions[subRegionIndex]]] += elementSubRegion->getElementVolume()[cellIndex];
    }
  });

  // Second, compute the normalized volume of each fine elements
  elementRegion->forElementSubRegions< CellElementSubRegion >( [&]( auto * const elementSubRegion ) -> void
  {
    localIndex const subRegionIndex = elementSubRegion->getIndexInParent();
    for(localIndex cellIndex = 0; cellIndex< elementSubRegion->size() ; cellIndex++)
    {
      if( elementSubRegion->GhostRank()[cellIndex] >= 0 )
        continue;
      normalizeVolumes[cellIndex + offsetSubRegions[subRegionIndex]] =
        elementSubRegion->getElementVolume()[cellIndex] / aggregateVolumes[parts[cellIndex + offsetSubRegions[subRegionIndex]]];
    }
  });

  // Third, normalize the centers
  this->forElementSubRegions( [&]( auto * const elementSubRegion ) -> void
  {
    localIndex const subRegionIndex = elementSubRegion->getIndexInParent();
    for(localIndex cellIndex = 0; cellIndex< elementSubRegion->size() ; cellIndex++)
    {
      if( elementSubRegion->GhostRank()[cellIndex] >= 0 )
        continue;
      R1Tensor center = elementSubRegion->getElementCenter()[cellIndex];
      center *= normalizeVolumes[cellIndex + offsetSubRegions[subRegionIndex]];
      aggregateBarycenters[parts[cellIndex + offsetSubRegions[subRegionIndex]]] += center;
    }
  });

  // Convert from metis to GEOSX types
  array1d< localIndex > partsGEOS( parts.size() );
  for( localIndex fineCellIndex = 0; fineCellIndex < partsGEOS.size(); fineCellIndex++ )
  {
    partsGEOS[fineCellIndex] = integer_conversion< localIndex >( parts[fineCellIndex] );
  }
  aggregateSubRegion->CreateFromFineToCoarseMap(nbAggregates, partsGEOS, aggregateBarycenters, aggregateVolumes);
  // Aggregate global indexes are saved within the fine cells
  elementRegion->forElementSubRegions< CellElementSubRegion >( [&]( auto * elementSubRegion ) -> void
  {
    auto & aggregateIndexSave =
      elementSubRegion->template getReference< array1d< globalIndex > > ( CellElementSubRegion::viewKeyStruct::aggregateGlobalIndexString );
    for(int i = 0; i < elementSubRegion->size() ;i++)
    {
      aggregateIndexSave[i] = aggregateSubRegion->m_localToGlobalMap[partsGEOS[i]];
    }
  });
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, AggregateElementRegion, std::string const &, ManagedGroup * const )

} /* namespace geosx */
