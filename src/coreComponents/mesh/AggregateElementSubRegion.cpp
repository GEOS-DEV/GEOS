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
#include "AggregateElementSubRegion.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"

namespace geosx
{
AggregateElementSubRegion::AggregateElementSubRegion( string const & name,
                                                      dataRepository::ManagedGroup * const parent):
  ElementSubRegionBase( name, parent )
{
}

AggregateElementSubRegion::~AggregateElementSubRegion()
{
}

void AggregateElementSubRegion::CreateFromFineToCoarseMap( localIndex nbAggregates,
                                                           array1d< localIndex > const & fineToCoarse,
                                                           array1d< R1Tensor > const & barycenters,
                                                           array1d< real64 > const & volumes )
{
  this->resize( nbAggregates );
  m_elementCenter = barycenters;
  m_elementVolume = volumes;

  // Next lines are for computing the global index of the aggregates
  int mpiSize;
  int mpiRank;
  MPI_Comm_size( MPI_COMM_GEOSX, &mpiSize );
  MPI_Comm_rank( MPI_COMM_GEOSX, &mpiRank );
  array1d< localIndex > nbAggregatesPerRank(mpiSize);
  MPI_Allgather( &nbAggregates, 1, MPI_LONG_LONG,nbAggregatesPerRank.data(), 1, MPI_LONG_LONG, MPI_COMM_GEOSX);
  globalIndex offsetForGlobalIndex = 0;
  for( localIndex i = 0; i < mpiRank; i++ )
  {
    offsetForGlobalIndex += nbAggregatesPerRank[i];
  }
  for( localIndex i = 0; i < nbAggregates; i++ )
  {
    this->m_localToGlobalMap[i] = i + offsetForGlobalIndex;
  }
  this->ConstructGlobalToLocalMap();

  // Aggregate global indexes are saved within the fine cells
  ElementRegion const * elementRegion = this->getParent()->getParent()->group_cast<ElementRegion const *>();
  elementRegion->forElementSubRegions( [&]( auto * elementSubRegion ) -> void
  {
    auto & aggregateIndexSave =
      elementSubRegion->template getWrapper< array1d< globalIndex > > ( CellElementSubRegion::viewKeyStruct::aggregateGlobalIndexString )->reference();
    for(int i = 0; i < elementSubRegion->size() ;i++)
    {
      aggregateIndexSave[i] = m_localToGlobalMap[fineToCoarse[i]];
    }
  });
}
}
