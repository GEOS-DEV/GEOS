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
 * @file MsrsbUtils.cpp
 */

#include "MsrsbUtils.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/multiscale/mesh/DofManager.hpp"
#include "linearAlgebra/multiscale/mesh/MeshObjectManager.hpp"
#include "linearAlgebra/multiscale/mesh/MeshUtils.hpp"
#include "linearAlgebra/multiscale/mesh/MeshLevel.hpp"
#include "LvArray/src/sortedArrayManipulation.hpp"

namespace geos
{
namespace multiscale
{
namespace msrsb
{

template< typename NODE_PRED, typename DUAL_PRED >
static ArrayOfSets< localIndex >
buildLocalConnectivityImpl( localIndex const numNodes,
                            MeshObjectManager const & nodeManager,
                            NODE_PRED && nodePred,
                            MeshObjectManager const & dualManager,
                            DUAL_PRED && dualPred )
{
  MeshObjectManager::MapViewConst const nodeToDual = nodeManager.toDualRelation().toViewConst();
  MeshObjectManager::MapViewConst const dualToNode = dualManager.toDualRelation().toViewConst();

  // Collect row sizes
  array1d< localIndex > rowLength( numNodes );
  forAll< parallelHostPolicy >( numNodes, [=, rowLength = rowLength.toView()]( localIndex const k )
  {
    meshUtils::forUniqueNeighbors< 512 >( k, nodeToDual, dualToNode, dualPred, [&]( localIndex const n )
    {
      if( nodePred( n ) )
      {
        ++rowLength[k];
      }
    } );
  } );

  // Resize
  ArrayOfSets< localIndex > conn( numNodes );
  conn.resizeFromCapacities< parallelHostPolicy >( rowLength.size(), rowLength.data() );

  // Fill the map
  forAll< parallelHostPolicy >( numNodes, [=, conn = conn.toView()]( localIndex const k )
  {
    meshUtils::forUniqueNeighbors< 512 >( k, nodeToDual, dualToNode, dualPred, [&]( localIndex const n )
    {
      if( nodePred( n ) )
      {
        conn.insertIntoSet( k, n );
      }
    } );
  } );

  return conn;
}

template< typename NODE_PRED >
static ArrayOfSets< localIndex >
buildLocalConnectivityImpl( localIndex const numNodes,
                            MeshObjectManager const & nodeManager,
                            NODE_PRED && nodePred,
                            MeshObjectManager const & dualManager,
                            bool const ghostDuals )
{
  if( ghostDuals )
  {
    return buildLocalConnectivityImpl( numNodes,
                                       nodeManager,
                                       std::forward< NODE_PRED >( nodePred ),
                                       dualManager,
                                       []( auto ){ return true; } );
  }
  else
  {
    arrayView1d< integer const > const dualGhostRank = dualManager.ghostRank();
    auto const noDualGhosts = [dualGhostRank]( localIndex const i ){ return dualGhostRank[i] < 0; };
    return buildLocalConnectivityImpl( numNodes,
                                       nodeManager,
                                       std::forward< NODE_PRED >( nodePred ),
                                       dualManager,
                                       noDualGhosts );
  }
}

ArrayOfSets< localIndex >
buildLocalConnectivity( MeshObjectManager const & nodeManager,
                        bool const ghostNodes,
                        MeshObjectManager const & dualManager,
                        bool const ghostDuals )
{
  GEOS_MARK_FUNCTION;
  if( ghostNodes )
  {
    return buildLocalConnectivityImpl( nodeManager.size(),
                                       nodeManager,
                                       []( auto ){ return true; },
                                       dualManager,
                                       ghostDuals );
  }
  else
  {
    arrayView1d< integer const > const nodeGhostRank = nodeManager.ghostRank();
    auto const noNodeGhosts = [nodeGhostRank]( localIndex const i ){ return nodeGhostRank[i] < 0; };
    return buildLocalConnectivityImpl( nodeManager.size(),
                                       nodeManager,
                                       noNodeGhosts,
                                       dualManager,
                                       ghostDuals );
  }
}

array1d< localIndex >
makeSeededPartition( ArrayOfSetsView< localIndex const > const & connectivity,
                     arrayView1d< localIndex const > const & seeds,
                     ArrayOfSetsView< localIndex const > const & supports )
{
  GEOS_MARK_FUNCTION;

  localIndex const numParts = seeds.size();
  localIndex const numNodes = connectivity.size();

  // Algorithm below is a multi-cluster parallel BFS using atomics to construct fronts.
  // This is not the most efficient way, but it works sufficiently well for this.
  // Cluster assignment conflicts are resolved using neighbor majority rule.

  // Pre-reserve extra capacity for front arrays to avoid reallocations
  array1d< localIndex > front( numParts );
  front.reserve( numNodes );
  array1d< localIndex > newPart;
  newPart.reserve( numNodes );

  // Initialize the partitions and expansion front
  array1d< localIndex > part( numNodes );
  part.setValues< parallelHostPolicy >( -1 );
  front.resize( numParts );
  forAll< parallelHostPolicy >( numParts, [&]( localIndex const ip )
  {
    part[seeds[ip]] = ip;
    front[ip] = seeds[ip];
  } );

  // Use AoA with 1 array for its atomic emplace capability
  ArrayOfArrays< localIndex > newFront( 1, 12 * numNodes ); // plenty of extra capacity

  while( true )
  {
    // 1. Expand the neighbor connectivity and build new expansion front
    newFront.clearArray( 0 );
    forAll< parallelHostPolicy >( front.size(), [connectivity,
                                                 front = front.toViewConst(),
                                                 part = part.toViewConst(),
                                                 newFront = newFront.toView()]( localIndex const i )
    {
      meshUtils::forUniqueNeighborValues< 256 >( front[i], connectivity,
                                                 []( localIndex const _ ){ return _; }, // just unique neighbor indices
                                                 [&]( localIndex const n ){ return part[n] < 0; }, // only unassigned nodes
                                                 [&]( localIndex const n )
      {
        newFront.emplaceBackAtomic< parallelHostAtomic >( 0, n );
      } );
    } );

    // Make the front unique
    front.clear();
    auto const nf = newFront[0];
    std::sort( nf.begin(), nf.end() );
    front.insert( 0, nf.begin(), std::unique( nf.begin(), nf.end() ) );

    newPart.resize( front.size() );
    newPart.setValues< parallelHostPolicy >( -1 );

    // 2. Assign partitions to the front nodes based on majority among neighbors
    RAJA::ReduceSum< parallelHostReduce, localIndex > numAssigned = 0;
    forAll< parallelHostPolicy >( front.size(), [connectivity, supports, numAssigned,
                                                 front = front.toViewConst(),
                                                 part = part.toViewConst(),
                                                 newPart = newPart.toView()]( localIndex const i )
    {
      localIndex const k = front[i];
      localIndex maxCount = 0;
      meshUtils::forUniqueNeighborValues< 256 >( k, connectivity, part,
                                                 []( localIndex const p ){ return p >= 0; }, // only assigned nodes
                                                 [&]( localIndex const p, localIndex const count )
      {
        if( count > maxCount && ( supports.size() == 0 || supports.contains( k, p ) ) )
        {
          newPart[i] = p;
          maxCount = count;
        }
      } );

      if( maxCount > 0 )
      {
        numAssigned += 1;
      }
    } );

    // 3. Terminate the loop as soon as no new assignments are made
    if( numAssigned.get() == 0 )
    {
      break;
    }

    // 4. Copy new assignments into partition array
    forAll< parallelHostPolicy >( front.size(), [front = front.toViewConst(),
                                                 newPart = newPart.toViewConst(),
                                                 part = part.toView()]( localIndex const i )
    {
      part[front[i]] = newPart[i];
    } );
  }

  // Attempt to fix unassigned front nodes, if any
  if( supports.size() > 0 )
  {
    GEOS_WARNING_IF( !front.empty(), "[MsRSB]: nodes not assigned to initial partition: " << front );
    forAll< parallelHostPolicy >( front.size(), [=, front = front.toViewConst(),
                                                 part = part.toView()]( localIndex const i )
    {
      localIndex const k = front[i];
      part[k] = supports( k, 0 ); // assign to the first support the node belongs to
    } );
  }
  else
  {
    GEOS_ERROR_IF( !front.empty(), "[MsRSB]: nodes not assigned to initial partition: " << front );
  }

  return part;
}

ArrayOfSets< localIndex >
buildSupports( ArrayOfSetsView< localIndex const > const & fineObjectToSubdomain,
               ArrayOfSetsView< localIndex const > const & subdomainToCoarseObject,
               ArrayOfSetsView< localIndex const > const & coarseObjectToSubdomain )
{
  GEOS_MARK_FUNCTION;

  // Algorithm:
  // Loop over all fine nodes.
  // - Get a list of adjacent coarse cells.
  // - If list is length 1, assign the node to supports of all coarse nodes adjacent to that coarse cell.
  //   Otherwise, collect a unique list of candidate coarse nodes by visiting them through coarse cells.
  // - For each candidate, check that fine node's subdomain list is included in the candidate's subdomain list.
  //   Otherwise, discard the candidate.
  // The first two cases could be handled by the last one; they are just an optimization that avoids some checks.
  // All above is done twice: once to count (or get upper bound on) row lengths, once to actually build supports.
  // For the last case, don't need to check inclusion when counting, just use number of candidates as upper bound.

  // Count row lengths and fill boundary indicators
  array1d< localIndex > rowLengths( fineObjectToSubdomain.size() );
  forAll< parallelHostPolicy >( fineObjectToSubdomain.size(),
                                [=, rowLengths = rowLengths.toView()]( localIndex const fidx )
  {
    if( fineObjectToSubdomain.sizeOfSet( fidx ) == 1 )
    {
      rowLengths[fidx] = subdomainToCoarseObject.sizeOfSet( fineObjectToSubdomain( fidx, 0 ) );
    }
    else
    {
      localIndex numCoarseObjects = 0;
      meshUtils::forUniqueNeighbors< 512 >( fidx, fineObjectToSubdomain, subdomainToCoarseObject, [&]( localIndex )
      {
        ++numCoarseObjects;
      } );
      rowLengths[fidx] = numCoarseObjects;
    }
  } );

  // Create and resize
  ArrayOfSets< localIndex > supports;
  supports.resizeFromCapacities< parallelHostPolicy >( rowLengths.size(), rowLengths.data() );

  // Fill the map
  forAll< parallelHostPolicy >( fineObjectToSubdomain.size(),
                                [=, supports = supports.toView()]( localIndex const fidx )
  {
    if( fineObjectToSubdomain.sizeOfSet( fidx ) == 1 )
    {
      arraySlice1d< localIndex const > const coarseObjects = subdomainToCoarseObject[fineObjectToSubdomain( fidx, 0 )];
      supports.insertIntoSet( fidx, coarseObjects.begin(), coarseObjects.end() );
    }
    else
    {
      arraySlice1d< localIndex const > const fsubs = fineObjectToSubdomain[fidx];
      meshUtils::forUniqueNeighbors< 512 >( fidx, fineObjectToSubdomain, subdomainToCoarseObject, [&]( localIndex const cidx )
      {
        arraySlice1d< localIndex const > const csubs = coarseObjectToSubdomain[cidx];
        if( std::includes( csubs.begin(), csubs.end(), fsubs.begin(), fsubs.end() ) )
        {
          supports.insertIntoSet( fidx, cidx );
        }
      } );
    }
  } );

  return supports;
}

array1d< integer >
findGlobalSupportBoundary( ArrayOfSetsView< localIndex const > const & fineObjectToSubdomain )
{
  array1d< integer > supportBoundaryIndicator( fineObjectToSubdomain.size() );
  forAll< parallelHostPolicy >( fineObjectToSubdomain.size(),
                                [=, supportBoundaryIndicator = supportBoundaryIndicator.toView()]( localIndex const fidx )
  {
    supportBoundaryIndicator[fidx] = fineObjectToSubdomain.sizeOfSet( fidx ) > 1;
  } );
  return supportBoundaryIndicator;
}

SparsityPattern< globalIndex >
buildProlongationSparsity( DofManager const & fineDofManager,
                           DofManager const & coarseDofManager,
                           string const & fieldName,
                           ArrayOfSetsView< localIndex const > const & supports )
{
  GEOS_MARK_FUNCTION;

  integer const numComp = fineDofManager.numComponents( fieldName );
  localIndex const numFineObjects = fineDofManager.numLocalDofs( fieldName ) / numComp;

  // This assumes an SDC-type pattern (i.e. no coupling between dof components on the same node)
  array1d< localIndex > rowLengths( fineDofManager.numLocalDofs( fieldName ) );
  forAll< parallelHostPolicy >( numFineObjects, [=, rowLengths = rowLengths.toView()]( localIndex const k )
  {
    for( integer ic = 0; ic < numComp; ++ic )
    {
      rowLengths[k * numComp + ic] = supports.sizeOfSet( k );
    }
  } );

  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( rowLengths.size(),
                                                         LvArray::integerConversion< localIndex >( coarseDofManager.numGlobalDofs( fieldName ) ),
                                                         rowLengths.data() );

  arrayView1d< globalIndex const > const fineDofNumber =
    fineDofManager.manager( fieldName ).getReference< array1d< globalIndex > >( fineDofManager.key( fieldName ) );

  arrayView1d< globalIndex const > const coarseDofNumber =
    coarseDofManager.manager( fieldName ).getReference< array1d< globalIndex > >( coarseDofManager.key( fieldName ) );

  globalIndex const fineDofOffset = fineDofManager.globalOffset( fieldName );

  forAll< parallelHostPolicy >( numFineObjects, [=, pattern = pattern.toView()]( localIndex const fineIdx )
  {
    localIndex const row = static_cast< localIndex >( fineDofNumber[fineIdx] - fineDofOffset );
    for( localIndex const coarseIdx : supports[fineIdx] )
    {
      globalIndex const col = coarseDofNumber[coarseIdx];
      for( integer ic = 0; ic < numComp; ++ic )
      {
        pattern.insertNonZero( row + ic, col + ic );
      }
    }
  } );

  return pattern;
}

CRSMatrix< real64, globalIndex >
buildTentativeProlongation( DofManager const & fineDofManager,
                            DofManager const & coarseDofManager,
                            string const & fieldName,
                            ArrayOfSetsView< localIndex const > const & supports,
                            arrayView1d< localIndex const > const & initPart )
{
  GEOS_MARK_FUNCTION;

  // Construct the tentative prolongation, consuming the sparsity pattern
  CRSMatrix< real64, globalIndex > localMatrix;
  {
    SparsityPattern< globalIndex > localPattern =
      buildProlongationSparsity( fineDofManager, coarseDofManager, fieldName, supports );
    localMatrix.assimilate< parallelHostPolicy >( std::move( localPattern ) );
  }

  arrayView1d< globalIndex const > const fineDofNumber =
    fineDofManager.manager( fieldName ).getReference< array1d< globalIndex > >( fineDofManager.key( fieldName ) );

  arrayView1d< globalIndex const > const coarseDofNumber =
    coarseDofManager.manager( fieldName ).getReference< array1d< globalIndex > >( coarseDofManager.key( fieldName ) );

  integer const numComp = fineDofManager.numComponents( fieldName );
  localIndex const numFineObjects = fineDofManager.numLocalDofs( fieldName ) / numComp;
  globalIndex const fineDofOffset = fineDofManager.globalOffset( fieldName );

  // Add initial unity values
  forAll< parallelHostPolicy >( numFineObjects, [=, localMatrix = localMatrix.toViewConstSizes()]( localIndex const fineIdx )
  {
    GEOS_ASSERT_GE( initPart[fineIdx], 0 );
    localIndex const row = static_cast< localIndex >( fineDofNumber[fineIdx] - fineDofOffset );
    GEOS_ASSERT_GE( row, 0 );
    globalIndex const col = coarseDofNumber[initPart[fineIdx]];
    real64 const value = 1.0;
    for( integer ic = 0; ic < numComp; ++ic )
    {
      globalIndex const col2 = col + ic;
      localMatrix.addToRow< serialAtomic >( row + ic, &col2, &value, 1 );
    }
  } );

  return localMatrix;
}

void makeGlobalDofLists( DofManager const & dofManager,
                         string const & fieldName,
                         arrayView1d< integer const > const & indicator,
                         array1d< globalIndex > & boundaryDof,
                         array1d< globalIndex > & interiorDof )
{
  integer const numComp = dofManager.numComponents( fieldName );
  localIndex const numLocalDofs = dofManager.numLocalDofs( fieldName );
  localIndex const numLocalObjects = numLocalDofs / numComp;

  arrayView1d< globalIndex const > const dofNumber =
    dofManager.manager( fieldName ).getReference< array1d< globalIndex > >( dofManager.key( fieldName ) );

  boundaryDof.clear();
  interiorDof.clear();

  boundaryDof.reserve( numLocalDofs );
  interiorDof.reserve( numLocalDofs );

  forAll< serialPolicy >( numLocalObjects, [&, dofNumber]( localIndex const objIdx )
  {
    globalIndex const globalDof = dofNumber[objIdx];
    if( indicator[objIdx] )
    {
      for( integer c = 0; c < numComp; ++c )
      {
        boundaryDof.emplace_back( globalDof + c );
      }
    }
    else
    {
      for( integer c = 0; c < numComp; ++c )
      {
        interiorDof.emplace_back( globalDof + c );
      }
    }
  } );
}

static ArrayOfSets< localIndex >
expandSupports( ArrayOfSetsView< localIndex const > const & connectivity,
                ArrayOfSetsView< localIndex const > const & support )
{
  GEOS_MARK_FUNCTION;

  localIndex const numPts = connectivity.size();
  array1d< localIndex > rowCounts( numPts );

  // Estimate an upper bound on new row sizes
  forAll< parallelHostPolicy >( numPts, [connectivity = connectivity.toViewConst(),
                                         support = support.toViewConst(),
                                         rowCounts = rowCounts.toView()]( localIndex const k )
  {
    rowCounts[k] = support.sizeOfSet( k );
    auto const isNewValue = [&support, k]( localIndex const i ){ return !support.contains( k, i ); };
    for( localIndex const n : connectivity[k] )
    {
      auto const s = support[n];
      rowCounts[k] += static_cast< localIndex >( std::count_if( s.begin(), s.end(), isNewValue ) );
    }
  } );

  // Allocate memory for new supports
  ArrayOfSets< localIndex > result;
  result.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Populate new supports by merging together supports of all neighbors of a cell
  forAll< parallelHostPolicy >( numPts, [connectivity, support,
                                         newSupport = result.toView()]( localIndex const k )
  {
    for( localIndex const n : connectivity[k] )
    {
      auto const s = support[n];
      newSupport.insertIntoSet( k, s.begin(), s.end() );
    }
  } );

  return result;
}

ArrayOfSets< localIndex >
buildLayeredSupport( integer const numLayers,
                     ArrayOfSetsView< localIndex const > const & connectivity,
                     arrayView1d< localIndex const > const & initialPartition )
{
  GEOS_ASSERT_EQ( connectivity.size(), initialPartition.size() );
  localIndex const numPts = initialPartition.size();
  ArrayOfSets< localIndex > support( numPts, 1 );

  // Build initial support
  forAll< parallelHostPolicy >( numPts, [support = support.toView(),
                                         initialPartition]( localIndex const k )
  {
    support.insertIntoSet( k, initialPartition[k] );
  } );

  // Main expansion iteration
  array1d< localIndex > rowCounts( numPts );
  for( int iter = 0; iter < numLayers; ++iter )
  {
    support = expandSupports( connectivity, support.toViewConst() );
  }

  return support;
}

array1d< integer >
findLayeredSupportBoundary( ArrayOfSetsView< localIndex const > const & connectivity,
                            ArrayOfSetsView< localIndex const > const & support )
{
  ArrayOfSets< localIndex > const supportPlus = expandSupports( connectivity, support );

  // Mark support boundary
  array1d< integer > boundaryIndicator( connectivity.size() );
  forAll< parallelHostPolicy >( connectivity.size(), [supportPlus = supportPlus.toViewConst(),
                                                      indicator = boundaryIndicator.toView(),
                                                      support]( localIndex const k )
  {
    // Criteria: k appears in supportPlus for some basis function i but not in support
    for( localIndex const i : supportPlus[k] )
    {
      if( !support.contains( k, i ) )
      {
        indicator[k] = 1;
      }
    }
  } );

  return boundaryIndicator;
}

CRSMatrix< real64, globalIndex >
dropEntries( CRSMatrixView< real64 const, globalIndex const > const & mat,
             real64 const relTol )
{
  array1d< localIndex > rowCounts( mat.numRows() );
  array1d< real64 > rowTols( mat.numRows() );
  forAll< parallelHostPolicy >( mat.numRows(), [rowCounts = rowCounts.toView(),
                                                rowTols = rowTols.toView(),
                                                mat, relTol]( localIndex const i )
  {
    real64 maxVal = 0.0;
    for( real64 const v : mat.getEntries( i ) )
    {
      maxVal = LvArray::math::max( maxVal, LvArray::math::abs( v ) );
    }
    real64 const absTol = maxVal * relTol;
    localIndex count = 0;
    for( real64 const v : mat.getEntries( i ) )
    {
      count += LvArray::math::abs( v ) >= absTol;
    }
    rowCounts[i] = count;
    rowTols[i] = absTol;
  } );

  CRSMatrix< real64, globalIndex > result;
  result.resizeFromRowCapacities< parallelHostPolicy >( rowCounts.size(), mat.numColumns(), rowCounts.data() );

  forAll< parallelHostPolicy >( mat.numRows(), [mat,
                                                rowTols = rowTols.toViewConst(),
                                                result = result.toView()]( localIndex const i )
  {
    real64 const absTol = rowTols[i];
    auto const columns = mat.getColumns( i );
    auto const values = mat.getEntries( i );
    for( localIndex k = 0; k < columns.size(); ++k )
    {
      if( LvArray::math::abs( values[k] ) >= absTol )
      {
        result.insertNonZero( i, columns[k], values[k] );
      }
    }
  } );

  return result;
}

void writeProlongation( CRSMatrixView< real64 const, globalIndex const > const & prolongation,
                        multiscale::DofManager const & dofManager,
                        string const & fieldName,
                        string const & prefix,
                        multiscale::MeshLevel & mesh,
                        multiscale::MeshObjectManager & fineManager,
                        std::function< void ( multiscale::MeshLevel &, std::vector< string > const & ) > const & writeFunc )
{
  std::vector< string > bNames{ "X ", "Y ", "Z " };
  std::vector< string > cNames{ " x", " y", " z" };

  integer const numComp = dofManager.numComponents( fieldName );
  globalIndex const numBfuncs = prolongation.numColumns() / numComp;
  int const labelWidth = static_cast< int >( std::log10( numBfuncs ) ) + 1;

  std::vector< arrayView3d< real64 > > views;
  std::vector< string > names;

  for( globalIndex bfIndex = 0; bfIndex < numBfuncs; ++bfIndex )
  {
    string const name = GEOS_FMT( "{}_P_{:0{}}", prefix, bfIndex, labelWidth );
    auto & array = fineManager.registerWrapper< array3d< real64 > >( name ).
                     setDimLabels( 1, { bNames.begin(), bNames.begin() + numComp } ).
                     setDimLabels( 2, { cNames.begin(), cNames.begin() + numComp } ).
                     setPlotLevel( dataRepository::PlotLevel::LEVEL_0 ).reference();
    array.resizeDimension< 1, 2 >( numComp, numComp );
    views.push_back( array.toView() );
    names.push_back( name );
  }

  auto const rowDofIndex = fineManager.getReference< array1d< globalIndex > >( dofManager.key( fieldName ) ).toViewConst();
  globalIndex const rowRankOffset = dofManager.rankOffset( fieldName );

  forAll< parallelHostPolicy >( fineManager.numOwnedObjects(), [=, &views]( localIndex const k )
  {
    localIndex const localRow = LvArray::integerConversion< localIndex >( rowDofIndex[k] - rowRankOffset );

    for( integer c = 0; c < numComp; ++c )
    {
      auto const columns = prolongation.getColumns( localRow + c );
      auto const values = prolongation.getEntries( localRow + c );

      for( localIndex i = 0; i < columns.size(); ++i )
      {
        globalIndex const bfIndex = columns[i] / numComp;
        integer const bfComp = LvArray::integerConversion< integer >( columns[i] % numComp );
        views[bfIndex]( k, bfComp, c ) = values[i];
      }
    }
  } );

  string_array fieldNames;
  fieldNames.insert( 0, names.begin(), names.end() );
  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       fineManager,
                                                       mesh.domain()->getNeighbors(),
                                                       false );

  writeFunc( mesh, names );
  for( string const & name : names )
  {
    fineManager.deregisterWrapper( name );
  }
}

} // namespace msrsb
} // namespace multiscale
} // namespace geos
