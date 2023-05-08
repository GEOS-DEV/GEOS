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
 * @file MultiscaleMeshUtils.cpp
 */

#include "MeshUtils.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geos
{
namespace multiscale
{
namespace meshUtils
{

void ScopedDataRegistrar::sync( DomainPartition & domain ) const
{
  array1d< string > fields;
  fields.emplace_back( m_key );
  CommunicationTools::getInstance().synchronizeFields( fields, m_manager, domain.getNeighbors(), false );
}

void copySets( ObjectManagerBase const & srcManager,
               string const & mapKey,
               ObjectManagerBase & dstManager )
{
  arrayView1d< localIndex const > const map = srcManager.getReference< array1d< localIndex > >( mapKey );
  srcManager.sets().forWrappers< SortedArray< localIndex > >( [&]( dataRepository::Wrapper< SortedArray< localIndex > > const & setWrapper )
  {
    if( setWrapper.getName() != "all" ) // no use for "all" in MS mesh
    {
      SortedArrayView< localIndex const > const srcSet = setWrapper.referenceAsView();
      SortedArray< localIndex > & dstSet = dstManager.createSet( setWrapper.getName() );
      meshUtils::mapIndexSet( srcSet, map, dstSet );
    }
  } );
}

template< typename T >
struct SetCompare
{
  ArrayOfSetsView< T const > const sets;
  bool operator()( localIndex const i, localIndex const j ) const
  {
    arraySlice1d< T const > const si = sets[i];
    arraySlice1d< T const > const sj = sets[j];
    return std::lexicographical_compare( si.begin(), si.end(), sj.begin(), sj.end() );
  }
};

array1d< localIndex >
findCoarseNodesByDualPartition( MeshObjectManager::MapViewConst const & nodeToDual,
                                MeshObjectManager::MapViewConst const & dualToNode,
                                ArrayOfSetsView< globalIndex const > const & nodeToSubdomain,
                                integer const minSubdomains,
                                bool allowMultiNodes )
{
  GEOS_MARK_FUNCTION;

  // Construct a list of "skeleton" nodes (those with at least minSubdomains adjacent subdomains)
  array1d< localIndex > skelNodes;
  for( localIndex inf = 0; inf < nodeToDual.size(); ++inf )
  {
    if( nodeToSubdomain.sizeOfSet( inf ) >= minSubdomains )
    {
      skelNodes.emplace_back( inf );
    }
  }

  // Sort skeleton nodes according to subdomain lists to locate nodes of identical adjacencies
  SetCompare< globalIndex > const adjacencyComp{ nodeToSubdomain.toViewConst() };
  //RAJA::sort< parallelHostPolicy >( RAJA::make_span( skelNodes.begin(), skelNodes.size() ), adjacencyComp );
  std::sort( skelNodes.begin(), skelNodes.end(), adjacencyComp );


  // Identify "features" (groups of skeleton nodes with the same subdomain adjacency)
  array1d< localIndex > const featureIndex( nodeToDual.size() );
  featureIndex.setValues< parallelHostPolicy >( -1 );

  ArrayOfArrays< localIndex > featureNodes;
  featureNodes.reserve( skelNodes.size() ); // overallocate to avoid reallocation
  featureNodes.reserveValues( skelNodes.size() ); // precise allocation

  localIndex numFeatures = 0;
  featureNodes.appendArray( 0 );
  featureNodes.emplaceBack( numFeatures, skelNodes[0] );
  for( localIndex i = 1; i < skelNodes.size(); ++i )
  {
    if( adjacencyComp( skelNodes[i-1], skelNodes[i] ) )
    {
      ++numFeatures;
      featureNodes.appendArray( 0 );
    }
    featureNodes.emplaceBack( numFeatures, skelNodes[i] );
    featureIndex[skelNodes[i]] = numFeatures;
  }
  ++numFeatures;

  // Construct feature-to-feature adjacency
  ArrayOfSets< localIndex > const featureAdjacency( numFeatures, 64 );
  forAll< parallelHostPolicy >( numFeatures,
                                [nodeToDual, dualToNode,
                                 featureNodes = featureNodes.toViewConst(),
                                 featureIndex = featureIndex.toViewConst(),
                                 featureAdjacency = featureAdjacency.toView()]( localIndex const f )
  {
    for( localIndex const inf : featureNodes[f] )
    {
      meshUtils::forUniqueNeighbors< 512 >( inf, nodeToDual, dualToNode, [&]( localIndex const nbrIdx )
      {
        if( nbrIdx != inf && featureIndex[nbrIdx] >= 0 )
        {
          featureAdjacency.insertIntoSet( f, featureIndex[nbrIdx] );
        }
      } );
    }
  } );

  // Choose features that represent coarse nodes (highest adjacency among neighbors)
  array1d< integer > const isCoarseNode( numFeatures );
  forAll< parallelHostPolicy >( numFeatures, [isCoarseNode = isCoarseNode.toView(),
                                              featureNodes = featureNodes.toViewConst(),
                                              featureAdjacency = featureAdjacency.toViewConst(),
                                              nodeToSubdomain = nodeToSubdomain.toViewConst()]( localIndex const f )
  {
    arraySlice1d< globalIndex const > const subs = nodeToSubdomain[ featureNodes( f, 0 ) ];
    for( localIndex const f_nbr : featureAdjacency[f] )
    {
      if( f_nbr != f )
      {
        arraySlice1d< globalIndex const > const subs_nbr = nodeToSubdomain[featureNodes( f_nbr, 0 )];
        if( std::includes( subs_nbr.begin(), subs_nbr.end(), subs.begin(), subs.end() ) )
        {
          // discard feature if its subdomain adjacency is fully included in any of its direct neighbors
          return;
        }
      }
    }
    // if not discarded, it is a coarse node
    isCoarseNode[f] = 1;
  } );

  // Make a list of fine-scale indices of coarse nodes that are locally owned
  array1d< localIndex > coarseNodes;
  for( localIndex f = 0; f < numFeatures; ++f )
  {
    if( isCoarseNode[f] == 1 )
    {
      arraySlice1d< localIndex const > const nodes = featureNodes[f];
      if( allowMultiNodes )
      {
        for( localIndex inf: nodes )
        {
          coarseNodes.emplace_back( inf );
        }
      }
      else
      {
        coarseNodes.emplace_back( nodes[0] );
      }
    }
  }

  RAJA::sort< parallelHostPolicy >( RAJA::make_span( coarseNodes.begin(), coarseNodes.size() ) );
  return coarseNodes;
}

} // namespace meshUtils
} // namespace multiscale
} // namespace geos
