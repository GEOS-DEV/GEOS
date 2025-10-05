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
 * @file PTScotchPartitioner.cpp
 */

#include "PTScotchPartitioner.hpp"
#include <ptscotch.h>

#define GEOS_SCOTCH_CHECK( call ) \
  do { \
    auto const ierr = call; \
    GEOS_ERROR_IF_NE_MSG( ierr, 0, "Error in call to PT-Scotch library:\n" << #call ); \
  } while( false )


namespace geos
{
using namespace dataRepository;

static_assert( std::is_same< int64_t, SCOTCH_Num >::value,
               "Non-matching index types. Scotch must be configured with 64-bit indices (SCOTCH_Num)." );
static_assert( std::is_same< pmet_idx_t, SCOTCH_Num >::value,
               "pmet_idx_t must match SCOTCH_Num when using PTScotchPartitioner." );

PTScotchPartitioner::PTScotchPartitioner( string const & name,
                                          Group * const parent ):
  GraphPartitioner( name, parent )
{}

// NOTE: processCommandLineOverrides, setPartitionCounts, setNeighborsRank, and color
// are now implemented in the parent GraphPartitioner class and have been removed from here.

array1d< pmet_idx_t > PTScotchPartitioner::partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                                      arrayView1d< pmet_idx_t const > const & vertDist,
                                                      pmet_idx_t const numParts,
                                                      MPI_Comm comm,
                                                      int const numRefinements )
{
  GEOS_UNUSED_VAR( vertDist );
  GEOS_WARNING_IF( numRefinements > 0, "Partition refinement is not supported by the PTScotch partitioner and will be ignored." );

  SCOTCH_Num const numLocalVerts = graph.size();
  array1d< SCOTCH_Num > part( numLocalVerts );

  // If only one partition is requested, no partitioning is needed.
  // All elements are assigned to rank 0.
  if( numParts == 1 )
  {
    part.setValues< serialPolicy >( 0 );
    return part;
  }

  // Initialize the distributed graph structure for Scotch.
  SCOTCH_Dgraph * const gr = SCOTCH_dgraphAlloc();
  GEOS_SCOTCH_CHECK( SCOTCH_dgraphInit( gr, comm ) );

  SCOTCH_Num const numLocalEdges = graph.getOffsets()[numLocalVerts];

  // PT-Scotch writes into these arrays; in practice we discard them right after.
  // We must cast away constness, which is technically UB but is how the library is used.
  SCOTCH_Num * const offsets = const_cast< SCOTCH_Num * >( graph.getOffsets() );
  SCOTCH_Num * const edges = const_cast< SCOTCH_Num * >( graph.getValues() );

  // Build the distributed graph from the local CSR representation.
  GEOS_SCOTCH_CHECK( SCOTCH_dgraphBuild( gr,           // graphptr
                                         0,            // baseval (0-based indexing)
                                         numLocalVerts,// vertlocnbr
                                         numLocalVerts,// vertlocmax
                                         offsets,      // vertloctab
                                         offsets + 1,  // vendloctab
                                         nullptr,      // veloloctab (no vertex weights)
                                         nullptr,      // vlblloctab (no vertex labels)
                                         numLocalEdges,// edgelocnbr
                                         numLocalEdges,// edgelocsiz
                                         edges,        // edgeloctab
                                         nullptr,      // edgegsttab (no ghost edges)
                                         nullptr       // edloloctab (no edge weights)
                                         ) );

  // Verify the consistency of the distributed graph data structure.
  GEOS_SCOTCH_CHECK( SCOTCH_dgraphCheck( gr ) );

  // Use the default PT-Scotch strategy.
  SCOTCH_Strat * const strategy = SCOTCH_stratAlloc();
  GEOS_SCOTCH_CHECK( SCOTCH_stratInit( strategy ) );

  // Perform the partitioning.
  GEOS_SCOTCH_CHECK( SCOTCH_dgraphPart( gr, numParts, strategy, part.data() ) );

  // Clean up Scotch data structures.
  SCOTCH_stratExit( strategy );
  SCOTCH_dgraphExit( gr );

  return part;
}

REGISTER_CATALOG_ENTRY( PartitionerBase, PTScotchPartitioner, string const &, dataRepository::Group * const )

} // namespace geos
