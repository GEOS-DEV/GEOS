/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParticleElementRegion.cpp
 */

#include "ParticleElementRegion.hpp"
#include "AggregateElementSubRegion.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/SparsityPattern.hpp"
#include "metis.h"

namespace geosx
{
using namespace dataRepository;

ParticleElementRegion::ParticleElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::sourceParticleBlockNamesString(), &m_particleBlockNames ).
    setInputFlag( InputFlags::OPTIONAL );

  registerWrapper( viewKeyStruct::coarseningRatioString(), &m_coarseningRatio ).
    setInputFlag( InputFlags::OPTIONAL );
}

ParticleElementRegion::~ParticleElementRegion()
{}

void ParticleElementRegion::generateMesh( Group & particleBlocks )
{
  Group & elementSubRegions = this->getGroup( viewKeyStruct::elementSubRegions() );

  for( string const & particleBlockName : this->m_particleBlockNames )
  {
    ParticleElementSubRegion & subRegion = elementSubRegions.registerGroup< ParticleElementSubRegion >( particleBlockName );
    ParticleBlock & source = particleBlocks.getGroup< ParticleBlock >( subRegion.getName() );
    subRegion.copyFromParticleBlock( source );
  }
}

void ParticleElementRegion::generateAggregates( FaceManager const & faceManager,
                                            NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ) )
{
  GEOSX_MARK_FUNCTION;

  if( m_coarseningRatio <= 0. )
  {
    return;
  }

  Group & elementSubRegions = this->getGroup( viewKeyStruct::elementSubRegions() );
  localIndex regionIndex = getIndexInParent();
  AggregateElementSubRegion & aggregateSubRegion =
    elementSubRegions.registerGroup< AggregateElementSubRegion >( "coarse" );

  arrayView2d< localIndex const > const elemRegionList     = faceManager.elementRegionList();
  arrayView2d< localIndex const > const elemSubRegionList  = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const elemList           = faceManager.elementList();

  // Counting the total number of particle and number of vertices
  localIndex nbParticleElements = 0;
  this->forElementSubRegions< ParticleElementSubRegion, FaceElementSubRegion >( [&]( auto & elementSubRegion )
  {
    nbParticleElements += elementSubRegion.size();
  } );

  // Number of aggregate computation
  localIndex nbAggregates = LvArray::integerConversion< localIndex >( int(nbParticleElements * m_coarseningRatio) );
  GEOSX_LOG_RANK_0( "Generating " << nbAggregates  << " aggregates on region " << this->getName());

  // METIS variable declarations
  using idx_t = ::idx_t;
  idx_t options[METIS_NOPTIONS];                                    // Contains the METIS options
  METIS_SetDefaultOptions( options );                                 // ... That are set by default
  idx_t nnodes = LvArray::integerConversion< idx_t >( nbParticleElements );     // Number of connectivity graph nodes
  idx_t nconst = 1;                                                 // Number of balancy constraints
  idx_t objval;                                                     // Total communication volume
  array1d< idx_t > parts( nnodes );                                   // Map element index -> aggregate index
  idx_t nparts = LvArray::integerConversion< idx_t >( nbAggregates );       // Number of aggregates to be generated


  // Compute the connectivity graph
  SparsityPattern< idx_t, idx_t > graph( LvArray::integerConversion< idx_t >( nbParticleElements ),
                                         LvArray::integerConversion< idx_t >( nbParticleElements ) );
  localIndex nbConnections = 0;
  array1d< localIndex > offsetSubRegions( this->getSubRegions().size() );
  for( localIndex subRegionIndex = 1; subRegionIndex < offsetSubRegions.size(); subRegionIndex++ )
  {
    offsetSubRegions[subRegionIndex] = offsetSubRegions[subRegionIndex - 1] + this->getSubRegion( subRegionIndex ).size();
  }
  for( localIndex kf = 0; kf < faceManager.size(); ++kf )
  {
    if( elemRegionList[kf][0] == regionIndex && elemRegionList[kf][1] == regionIndex && elemRegionList[kf][0] )
    {
      localIndex const esr0 = elemSubRegionList[kf][0];
      idx_t const ei0  = LvArray::integerConversion< idx_t >( elemList[kf][0] + offsetSubRegions[esr0] );
      localIndex const esr1 = elemSubRegionList[kf][1];
      idx_t const ei1  = LvArray::integerConversion< idx_t >( elemList[kf][1] + offsetSubRegions[esr1] );
      graph.insertNonZero( ei0, ei1 );
      graph.insertNonZero( ei1, ei0 );
      nbConnections++;
    }
  }

  // METIS partitionning
  idx_t * offsets = const_cast< idx_t * >( graph.getOffsets() );
  idx_t * columns = const_cast< idx_t * >( &graph.getColumns( 0 )[0] );
  METIS_PartGraphRecursive( &nnodes, &nconst, offsets, columns, nullptr, nullptr, nullptr,
                            &nparts, nullptr, nullptr, options, &objval, parts.data() );

  // Compute Aggregate barycenters
  array2d< real64 > aggregateBarycenters( nparts, 3 );
  array1d< real64 > aggregateVolumes( nparts );
  array1d< real64 > normalizeVolumes( nbParticleElements );

  // First, compute the volume of each aggregates
  this->forElementSubRegions< ParticleElementSubRegion, FaceElementSubRegion >( [&]( ElementSubRegionBase & elementSubRegion )
  {
    arrayView1d< integer const > const ghostRank = elementSubRegion.ghostRank();
    localIndex const subRegionIndex = elementSubRegion.getIndexInParent();
    for( localIndex particleIndex = 0; particleIndex< elementSubRegion.size(); particleIndex++ )
    {
      if( ghostRank[particleIndex] >= 0 )
        continue;
      aggregateVolumes[parts[particleIndex + offsetSubRegions[subRegionIndex]]] += elementSubRegion.getElementVolume()[particleIndex];
    }
  } );

  // Second, compute the normalized volume of each fine elements
  this->forElementSubRegions< ParticleElementSubRegion, FaceElementSubRegion >( [&]( ElementSubRegionBase & elementSubRegion )
  {
    arrayView1d< integer const > const ghostRank = elementSubRegion.ghostRank();
    localIndex const subRegionIndex = elementSubRegion.getIndexInParent();
    for( localIndex particleIndex = 0; particleIndex< elementSubRegion.size(); particleIndex++ )
    {
      if( ghostRank[particleIndex] >= 0 )
        continue;
      normalizeVolumes[particleIndex + offsetSubRegions[subRegionIndex]] =
        elementSubRegion.getElementVolume()[particleIndex] / aggregateVolumes[parts[particleIndex + offsetSubRegions[subRegionIndex]]];
    }
  } );

  // Third, normalize the centers
  this->forElementSubRegions< ParticleElementSubRegion, FaceElementSubRegion >( [&]( ElementSubRegionBase & elementSubRegion )
  {
    arrayView1d< integer const > const ghostRank = elementSubRegion.ghostRank();
    localIndex const subRegionIndex = elementSubRegion.getIndexInParent();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();

    for( localIndex particleIndex = 0; particleIndex< elementSubRegion.size(); particleIndex++ )
    {
      if( ghostRank[particleIndex] >= 0 )
        continue;

      // TODO Change the rest of this to
      LvArray::tensorOps::scaledAdd< 3 >( aggregateBarycenters[ parts[ particleIndex + offsetSubRegions[ subRegionIndex ] ] ],
                                          elemCenter[ particleIndex ],
                                          normalizeVolumes[ particleIndex + offsetSubRegions[ subRegionIndex ] ] );
    }
  } );

  // Convert from metis to GEOSX types
  array1d< localIndex > partsGEOS( parts.size() );
  for( localIndex fineParticleIndex = 0; fineParticleIndex < partsGEOS.size(); fineParticleIndex++ )
  {
    partsGEOS[fineParticleIndex] = LvArray::integerConversion< localIndex >( parts[fineParticleIndex] );
  }
  aggregateSubRegion.createFromFineToCoarseMap( nbAggregates, partsGEOS, aggregateBarycenters );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ParticleElementRegion, string const &, Group * const )

} /* namespace geosx */
