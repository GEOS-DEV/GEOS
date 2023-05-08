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
 * @file MeshUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MESHUTILS_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MESHUTILS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/multiscale/mesh/MeshObjectManager.hpp"
#include "mesh/ObjectManagerBase.hpp"

namespace geos
{

class DomainPartition;

namespace multiscale
{
namespace meshUtils
{

/**
 * @brief Utility class to manage registration of temporary scope-bound data in the data repository.
 *
 * On construction, registers a wrapper for the given data container in the manager.
 * On destruction, cleans up by removing a wrapper.
 * As an additional benefit, provides a shortcut for parallel synchronization of data.
 * Does NOT copy/hold the data.
 */
class ScopedDataRegistrar
{
public:

  /**
   * @brief Constructor.
   * @tparam T type of data to register
   * @param manager object manager to register data on
   * @param key data repository key to use
   * @param data the data object
   */
  template< typename T >
  ScopedDataRegistrar( ObjectManagerBase & manager, string key, T & data )
    : m_manager( manager ),
    m_key( std::move( key ) )
  {
    m_manager.registerWrapper( m_key, &data );
  }

  ScopedDataRegistrar( ScopedDataRegistrar const & ) = delete;
  ScopedDataRegistrar( ScopedDataRegistrar && ) = delete;

  /**
   * @brief Destructor.
   */
  ~ScopedDataRegistrar()
  {
    m_manager.deregisterWrapper( m_key );
  }

  /**
   * @brief Synchronize data across ranks.
   * @param domain the domain object needed to access communicators
   */
  void sync( DomainPartition & domain ) const;

private:

  ObjectManagerBase & m_manager;
  string const m_key;
};

/**
 * @brief Transform an array of mesh indices using a map.
 * @tparam T source index type
 * @tparam U destination index type
 * @param src source index array
 * @param map index map
 * @param dst destination array
 *
 * This function assigns dst[i] = map[src[i]] for each i s.t. map[src[i]] >= 0.
 * Thus entries with negative indices in the map are filtered out.
 */
template< typename T, typename U = T >
void mapIndexArray( arrayView1d< T const > const & src,
                    arrayView1d< U const > const & map,
                    array1d< U > & dst )
{
  dst.clear();
  dst.reserve( src.size() );
  for( T const & val : src )
  {
    U const newVal = map[val];
    if( newVal >= 0 )
    {
      dst.emplace_back( newVal );
    }
  }
}

/**
 * @brief Transform an array of mesh indices using a map, removing duplicates.
 * @tparam T source index type
 * @tparam U destination index type
 * @param src source index array
 * @param map index map
 * @param dst destination array
 *
 * This function assigns dst[i] = map[src[i]] for each i s.t. map[src[i]] >= 0.
 * Thus entries with negative indices in the map are filtered out.
 * In addition, destination array only contains unique values.
 */
template< typename T, typename U = T >
void mapIndexArrayUnique( arrayView1d< T const > const & src,
                          arrayView1d< U const > const & map,
                          array1d< U > & dst )
{
  SortedArray< U > values;
  values.reserve( src.size() );
  for( T const & val : src )
  {
    U const newVal = map[val];
    if( newVal >= 0 )
    {
      values.insert( newVal );
    }
  }
  dst.clear();
  dst.reserve( values.size() );
  for( U const & val : values )
  {
    dst.emplace_back( val );
  }
}

/**
 * @brief Transform a set of mesh indices using a map.
 * @tparam T source index type
 * @tparam U destination index type
 * @param src source index set
 * @param map index map
 * @param dst destination set
 */
template< typename T, typename U = T >
void mapIndexSet( SortedArrayView< T const > const & src,
                  arrayView1d< U const > const & map,
                  SortedArray< U > & dst )
{
  dst.reserve( src.size() );
  for( T const & val : src )
  {
    U const newVal = map[val];
    if( newVal >= 0 )
    {
      dst.insert( newVal );
    }
  }
}

/**
 * @brief Populate an array using a function source and an index map.
 * @tparam POLICY loop execution policy
 * @tparam T type of value
 * @tparam NDIM number of array dimensions
 * @tparam USD array unit stride dimension
 * @tparam INDEX array index type
 * @tparam FUNC type of source function
 * @param map index map from source to destination
 * @param dst destination array
 * @param src source function called with source indices
 *
 * Elements for which index in the map is negative are ignored.
 */
template< typename POLICY, typename T, int NDIM, int USD, typename INDEX, typename FUNC >
void fillArrayByDstIndex( arrayView1d< INDEX const > const & map,
                          ArrayView< T, NDIM, USD > const & dst,
                          FUNC && src )
{
  forAll< POLICY >( map.size(), [=]( INDEX const srcIdx )
  {
    INDEX const dstIdx = map[srcIdx];
    if( dstIdx >= 0 )
    {
      LvArray::forValuesInSliceWithIndices( dst[dstIdx], [&]( T & dstValue, auto ... indices )
      {
        dstValue = src( srcIdx, indices ... );
      } );
    }
  } );
}

/**
 * @brief Populate an array using a function source and an index map.
 * @tparam POLICY loop execution policy
 * @tparam T type of value
 * @tparam NDIM number of array dimensions
 * @tparam USD array unit stride dimension
 * @tparam INDEX array index type
 * @tparam FUNC type of source function
 * @param map index map from destination to source
 * @param dst destination array
 * @param src source function called with source indices
 *
 * Elements for which index in the map is negative are ignored.
 */
template< typename POLICY, typename T, int NDIM, int USD, typename INDEX, typename FUNC >
void fillArrayBySrcIndex( arrayView1d< INDEX const > const & map,
                          ArrayView< T, NDIM, USD > const & dst,
                          FUNC && src )
{
  GEOS_ASSERT_EQ( dst.size(), map.size() );
  forAll< POLICY >( map.size(), [=]( INDEX const dstIndex )
  {
    INDEX const srcIndex = map[dstIndex];
    if( srcIndex >= 0 )
    {
      LvArray::forValuesInSliceWithIndices( dst[dstIndex], [&]( T & dstValue, auto ... indices )
      {
        dstValue = src( srcIndex, indices ... );
      } );
    }
  } );
}

/**
 * @brief Sets up NeighborData of a new object manager mapped from existing one.
 * @tparam FUNC type of copy utility function
 * @param srcManager source object manager
 * @param mapKey key to extract the index map from source manager
 * @param ranks ranks of neighbors
 * @param dstManager target object manager
 * @param copyFunc copy utility function to use on each array
 */
template< typename FUNC >
void copyNeighborData( ObjectManagerBase const & srcManager,
                       string const & mapKey,
                       std::vector< integer > const & ranks,
                       ObjectManagerBase & dstManager,
                       FUNC && copyFunc )
{
  arrayView1d< localIndex const > const map = srcManager.getReference< array1d< localIndex > >( mapKey );
  for( integer const rank : ranks )
  {
    NeighborData const & srcData = srcManager.getNeighborData( rank );
    NeighborData & dstData = dstManager.getNeighborData( rank );
    copyFunc( srcData.ghostsToSend(), map, dstData.ghostsToSend() );
    copyFunc( srcData.ghostsToReceive(), map, dstData.ghostsToReceive() );
    copyFunc( srcData.adjacencyList(), map, dstData.adjacencyList() );
    copyFunc( srcData.matchedPartitionBoundary(), map, dstData.matchedPartitionBoundary() );
  }
}

/**
 * @brief Copy and update sets from an existing to a new object manager.
 * @param srcManager source object manager
 * @param mapKey key to extract the index map from source manager
 * @param dstManager target object manager
 */
void copySets( ObjectManagerBase const & srcManager,
               string const & mapKey,
               ObjectManagerBase & dstManager );

namespace internal
{

// We don't have versions of IS_VALID_EXPRESSION for more args (and I don't want to add them).
// So second argument (count) is hardcoded to std::ptrdiff_t below.

IS_VALID_EXPRESSION( isCallableWithoutArgs, T, std::declval< T >()() );
IS_VALID_EXPRESSION_2( isCallableWithArg, T, U, std::declval< T >()( std::declval< U >() ) );
IS_VALID_EXPRESSION_2( isCallableWithArgAndCount, T, U, std::declval< T >()( std::declval< U >(), std::ptrdiff_t{} ) );

template< typename T, typename FUNC >
std::enable_if_t< isCallableWithoutArgs< FUNC > >
callWithArgs( T const & val, std::ptrdiff_t const count, FUNC && func )
{
  GEOS_UNUSED_VAR( val, count );
  func();
}

template< typename T, typename FUNC >
std::enable_if_t< isCallableWithArg< FUNC, T > >
callWithArgs( T const & val, std::ptrdiff_t const count, FUNC && func )
{
  GEOS_UNUSED_VAR( count );
  func( val );
}

template< typename T, typename FUNC >
std::enable_if_t< isCallableWithArgAndCount< FUNC, T > >
callWithArgs( T const & val, std::ptrdiff_t const count, FUNC && func )
{
  func( val, count );
}

} // namespace internal

/**
 * @brief Call a function on unique values from a previously collected range.
 * @tparam ITER type of range iterator
 * @tparam FUNC type of function to call
 * @param first start of the range
 * @param last end of the range
 * @param func the function to call
 * @note Modifies the range by sorting values in place, so @p ITER must not be a const iterator.
 */
template< typename ITER, typename FUNC >
void forUniqueValues( ITER first, ITER const last, FUNC && func )
{
  if( first == last )
    return;
  LvArray::sortedArrayManipulation::makeSorted( first, last );
  using T = typename std::iterator_traits< ITER >::value_type;
  while( first != last )
  {
    T const & curr = *first;
    ITER const it = std::find_if( first, last, [&curr]( T const & v ) { return v != curr; } );
    internal::callWithArgs( curr, std::distance( first, it ), std::forward< FUNC >( func ) );
    first = it;
  }
}

/**
 * @brief Call a function on unique indices of topological neighbors visited through location-connector adjacency maps.
 * @tparam MAX_NEIGHBORS maximum number of total non-unique neighbor indices
 * @tparam L2C_MAP type of location to connector map
 * @tparam C2L_MAP type of connector to location map
 * @tparam PRED type of predicate function
 * @tparam FUNC type of function to call
 * @param locIdx location index
 * @param locToConn location to connector map
 * @param connToLoc connector to location map
 * @param connPred predicate used to filter connector indices
 * @param func function to call
 */
template< integer MAX_NEIGHBORS, typename L2C_MAP, typename C2L_MAP, typename PRED, typename FUNC >
void forUniqueNeighbors( localIndex const locIdx,
                         L2C_MAP const & locToConn,
                         C2L_MAP const & connToLoc,
                         PRED && connPred,
                         FUNC && func )
{
  localIndex neighbors[MAX_NEIGHBORS];
  integer numNeighbors = 0;
  for( localIndex const connIdx : locToConn[locIdx] )
  {
    if( connIdx >= 0 && connPred( connIdx ) )
    {
      for( localIndex const nbrIdx : connToLoc[connIdx] )
      {
        GEOS_ERROR_IF_GE_MSG( numNeighbors, MAX_NEIGHBORS, "Too many neighbors, need to increase stack limit" );
        neighbors[numNeighbors++] = nbrIdx;
      }
    }
  }
  forUniqueValues( neighbors, neighbors + numNeighbors, std::forward< FUNC >( func ) );
}

/**
 * @brief Call a function on unique indices of topological neighbors visited through location-connector adjacency maps.
 * @tparam MAX_NEIGHBORS maximum number of total non-unique neighbor indices
 * @tparam L2C_MAP type of location to connector map
 * @tparam C2L_MAP type of connector to location map
 * @tparam FUNC type of function to call
 * @param locIdx location index
 * @param locToConn location to connector map
 * @param connToLoc connector to location map
 * @param func function to call
 *
 * Version that does not use a connector predicate.
 */
template< integer MAX_NEIGHBORS, typename L2C_MAP, typename C2L_MAP, typename FUNC >
void forUniqueNeighbors( localIndex const locIdx,
                         L2C_MAP const & locToConn,
                         C2L_MAP const & connToLoc,
                         FUNC && func )
{
  forUniqueNeighbors< MAX_NEIGHBORS >( locIdx,
                                       locToConn,
                                       connToLoc,
                                       []( auto ){ return true; },
                                       std::forward< FUNC >( func ) );
}

/**
 * @brief Call a function on unique values of topological neighbors visited through location-connector adjacency maps.
 * @tparam MAX_NEIGHBORS maximum number of total non-unique neighbor indices
 * @tparam NBR_MAP type of connectivity map
 * @tparam VAL_FUNC type of value function
 * @tparam VAL_PRED type of value predicate
 * @tparam FUNC type of target function
 * @param locIdx location index
 * @param neighbors neighbor map
 * @param valueFunc function called with neighbor indices that returns target values
 * @param pred predicate used to filter values
 * @param func function to call
 */
template< integer MAX_NEIGHBORS, typename NBR_MAP, typename VAL_FUNC, typename VAL_PRED, typename FUNC >
void forUniqueNeighborValues( localIndex const locIdx,
                              NBR_MAP const & neighbors,
                              VAL_FUNC && valueFunc,
                              VAL_PRED && pred,
                              FUNC && func )
{
  using T = std::remove_cv_t< std::remove_reference_t< decltype( valueFunc( localIndex {} ) ) >>;
  T nbrValues[MAX_NEIGHBORS];
  integer numValues = 0;
  for( localIndex const nbrIdx : neighbors[locIdx] )
  {
    GEOS_ERROR_IF_GE_MSG( numValues, MAX_NEIGHBORS, "Too many neighbors, need to increase stack limit" );
    T const value = valueFunc( nbrIdx );
    if( pred( value ) )
    {
      nbrValues[numValues++] = value;
    }
  }
  forUniqueValues( nbrValues, nbrValues + numValues, std::forward< FUNC >( func ) );
}

/**
 * @brief Call a function on unique values of topological neighbors visited through location-connector adjacency maps.
 * @tparam MAX_NEIGHBORS maximum number of total non-unique neighbor indices
 * @tparam NBR_MAP type of connectivity map
 * @tparam VAL_FUNC type of value function
 * @tparam FUNC type of target function
 * @param locIdx location index
 * @param neighbors neighbor map
 * @param valueFunc function called with neighbor indices that returns target values
 * @param func function to call
 *
 * Version that does not use a value predicate.
 */
template< integer MAX_NEIGHBORS, typename NBR_MAP, typename VAL_FUNC, typename FUNC >
void forUniqueNeighborValues( localIndex const locIdx,
                              NBR_MAP const & neighbors,
                              VAL_FUNC && valueFunc,
                              FUNC && func )
{
  forUniqueNeighborValues< MAX_NEIGHBORS >( locIdx,
                                            neighbors,
                                            std::forward< VAL_FUNC >( valueFunc ),
                                            []( auto ){ return true; },
                                            std::forward< FUNC >( func ) );
}

/**
 * @brief Build a map from mesh objects (nodes/cells) to subdomains defined by a partitioning of dual objects (cells/nodes).
 * @tparam INDEX_TYPE type of subdomain index
 * @param fineObjectManager
 * @param subdomains array of subdomain indices of dual objects
 * @param boundaryObjectSets (optional) list of boundary set names
 * @return array of subdomain index sets
 *
 * Boundaries (if present) are treated as additional "virtual" subdomains.
 * They are assigned unique negative indices to distinguish them from actual subdomains.
 * Pass an empty list of set names to ignore this feature.
 */
template< typename INDEX_TYPE >
ArrayOfSets< INDEX_TYPE >
buildFineObjectToSubdomainMap( MeshObjectManager const & fineObjectManager,
                               arrayView1d< INDEX_TYPE const > const & subdomains,
                               arrayView1d< string const > const & boundaryObjectSets )
{
  GEOS_MARK_FUNCTION;

  MeshObjectManager::MapViewConst const dualMap = fineObjectManager.toDualRelation().toViewConst();

  // count the row lengths
  array1d< localIndex > rowCounts( fineObjectManager.size() );
  forAll< parallelHostPolicy >( fineObjectManager.size(), [=, rowCounts = rowCounts.toView()]( localIndex const objIdx )
  {
    localIndex count = 0;
    forUniqueNeighborValues< 256 >( objIdx, dualMap, subdomains, [&count]
    {
      ++count;
    } );
    rowCounts[objIdx] = count;
  } );
  for( string const & setName : boundaryObjectSets )
  {
    SortedArrayView< localIndex const > const set = fineObjectManager.getSet( setName ).toViewConst();
    forAll< parallelHostPolicy >( set.size(), [=, rowCounts = rowCounts.toView()]( localIndex const i )
    {
      ++rowCounts[set[i]];
    } );
  }

  // Resize from row lengths
  ArrayOfSets< INDEX_TYPE > objectToSubdomain;
  objectToSubdomain.template resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  INDEX_TYPE numBoundaries = 0;
  for( string const & setName : boundaryObjectSets )
  {
    ++numBoundaries;
    SortedArrayView< localIndex const > const set = fineObjectManager.getSet( setName ).toViewConst();
    forAll< parallelHostPolicy >( set.size(), [=, objToSub = objectToSubdomain.toView()]( localIndex const i )
    {
      objToSub.insertIntoSet( set[i], -numBoundaries );
    } );
  }
  forAll< parallelHostPolicy >( fineObjectManager.size(), [=, objToSub = objectToSubdomain.toView()]( localIndex const objIdx )
  {
    for( localIndex const dualIdx : dualMap[objIdx] )
    {
      objToSub.insertIntoSet( objIdx, subdomains[dualIdx] );
    }
  } );

  return objectToSubdomain;
}

/**
 * @brief Insert virtual boundary subdomains into an adjacency map.
 * @tparam INDEX_TYPE type of index in the adjacency map
 * @param manager multiscale mesh object manager
 * @param inputMap source map
 * @param boundaryObjectSets names of boundary object sets
 * @return an updated map
 *
 * Boundary (virtual) subdomains are assigned negative indices [-1 ... boundaryObjectSets.size()] in the map.
 * For example, if the input map contains cells [1,4,6] for a given node, and the node is found in the first
 * and third boundary sets (in order listed), the output map will contain entries [-3,-1,1,4,6].
 */
template< typename INDEX_TYPE >
ArrayOfSets< INDEX_TYPE >
addBoundarySubdomains( MeshObjectManager const & manager,
                       ArrayOfSetsView< INDEX_TYPE const > const & inputMap,
                       arrayView1d< string const > const & boundaryObjectSets )
{
  GEOS_MARK_FUNCTION;

  array1d< localIndex > rowCounts( manager.size() );
  forAll< parallelHostPolicy >( manager.size(), [=, rowCounts = rowCounts.toView()]( localIndex const objIdx )
  {
    rowCounts[objIdx] = inputMap.sizeOfSet( objIdx );
  } );
  for( string const & setName: boundaryObjectSets )
  {
    SortedArrayView< localIndex const > const set = manager.getSet( setName ).toViewConst();
    forAll< parallelHostPolicy >( set.size(), [=, rowCounts = rowCounts.toView()]( localIndex const i )
    {
      ++rowCounts[set[i]];
    } );
  }

  // Resize from row lengths
  ArrayOfSets< INDEX_TYPE > outputMap;
  outputMap.template resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  INDEX_TYPE numBoundaries = 0;
  for( string const & setName: boundaryObjectSets )
  {
    ++numBoundaries;
    SortedArrayView< localIndex const > const set = manager.getSet( setName ).toViewConst();
    forAll< parallelHostPolicy >( set.size(), [=, outputMap = outputMap.toView()]( localIndex const i )
    {
      outputMap.insertIntoSet( set[i], -numBoundaries );
    } );
  }
  forAll< parallelHostPolicy >( manager.size(), [=, outputMap = outputMap.toView()]( localIndex const objIdx )
  {
    outputMap.insertIntoSet( objIdx, inputMap[objIdx].begin(), inputMap[objIdx].end() );
  } );

  return outputMap;
}

/**
 * @brief Find "coarse" nodes in a mesh in which dual objects have been partitioned into subdomains.
 * @param nodeToDual adjacency map of nodes to its dual object
 * @param dualToNode adjacency map of dual objects to nodes
 * @param nodeToSubdomain adjacency of nodes to subdomains (must be built by the caller)
 * @param minSubdomains minimum number of adjacent subdomains required for a coarse node (reduces search space)
 * @param allowMultiNodes whether clusters of similar nodes should produce multiple coarse nodes (if false, just one will be chosen)
 * @return an array containing indices of nodes identified as "coarse"
 *
 * "Coarse" nodes are defined as those that, given a partitioning of dual objects (e.g. mesh volumes or cells)
 * into subdomains, are adjacent to the locally maximal (among its neighbors) number of such subdomains.
 *
 * The implementation groups nodes into "features" (that are groups of nodes sharing the same subdomain list)
 * and finds features with largest list (or "key") among adjacent features. Such features typically consist of
 * just a single node. Occasionally, in complex topologies, coarse features with multiple nodes may be found.
 * In that case, any of the nodes in the feature can be chosen as coarse (since they all have the same adjacent
 * subdomains), or each of them can be made a separate coarse node (this is the default behavior).
 *
 * The analysis is fully local. In order for coarse nodes to be chosen consistently across processes,
 * adjacency maps must include ghosted node/dual objects and correct global subdomain indices.
 */
array1d< localIndex >
findCoarseNodesByDualPartition( MeshObjectManager::MapViewConst const & nodeToDual,
                                MeshObjectManager::MapViewConst const & dualToNode,
                                ArrayOfSetsView< globalIndex const > const & nodeToSubdomain,
                                integer const minSubdomains,
                                bool allowMultiNodes = true );

} // namespace meshUtils
} // namespace multiscale
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MESHUTILS_HPP
