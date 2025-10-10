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
 * @file PTScotchEngine.cpp
 */

#include "PTScotchEngine.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"

#ifdef GEOS_USE_PARMETIS
  #include "ParMetisEngine.hpp"  // For fallback
#endif


#ifdef GEOS_USE_PTSCOTCH
extern "C"
{
#include <ptscotch.h>
}

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
               "pmet_idx_t must match SCOTCH_Num when using PTScotchEngine." );

#else

namespace geos
{

using namespace dataRepository;

#endif // GEOS_USE_PTSCOTCH

PTScotchEngine::PTScotchEngine( string const & name,
                                dataRepository::Group * const parent )
  : GraphPartitionEngine( name, parent ),
  m_strategy( "default" )
{
  registerWrapper( viewKeyStruct::strategyString(), &m_strategy ).
    setApplyDefaultValue( "default" ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "PT-Scotch partitioning strategy. "
                    "Options: 'default', 'quality', 'speed', 'balance'. "
                    "Default: 'default'." );
}

PTScotchEngine::~PTScotchEngine()
{}

array1d< pmet_idx_t >
PTScotchEngine::partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                           arrayView1d< pmet_idx_t const > const & vertDist,
                           pmet_idx_t const numParts,
                           MPI_Comm comm )
{
  GEOS_MARK_FUNCTION;

#ifndef GEOS_USE_PTSCOTCH
  GEOS_UNUSED_VAR( graph, vertDist, numParts, comm );
  GEOS_ERROR( "GEOS was not built with PT-Scotch support. "
              "Reconfigure with -DENABLE_PTSCOTCH=ON" );
  return array1d< pmet_idx_t >();
#else

  GEOS_UNUSED_VAR( vertDist );

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

  GEOS_LOG_RANK_0( GEOS_FMT( "PTScotchEngine: Partitioning {} local vertices into {} parts (strategy: {})",
                             numLocalVerts, numParts, m_strategy ) );

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

  GEOS_LOG_RANK_0( "PTScotchEngine: Partition complete" );

  return part;

#endif // GEOS_USE_PTSCOTCH
}



ArrayOfArrays< pmet_idx_t, pmet_idx_t >
PTScotchEngine::meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
                            arrayView1d< pmet_idx_t const > const & elemDist,
                            MPI_Comm comm,
                            int const minCommonNodes )
{
#ifdef GEOS_USE_PARMETIS

  // Fallback to ParMETIS implementation
  return ParMetisEngine::parmetisMeshToDual( elemToNodes, elemDist, comm, minCommonNodes );

#else

  GEOS_UNUSED_VAR( elemToNodes );
  GEOS_UNUSED_VAR( elemDist );
  GEOS_UNUSED_VAR( comm );
  GEOS_UNUSED_VAR( minCommonNodes );

  GEOS_ERROR( "PTScotchEngine::meshToDual requires ParMETIS for mesh-to-dual conversion. "
              "Either build with GEOS_USE_PARMETIS=ON or provide the dual graph directly." );
  return ArrayOfArrays< pmet_idx_t, pmet_idx_t >();

#endif
}

// Register in the GraphPartitionEngine catalog
REGISTER_CATALOG_ENTRY( GraphPartitionEngine, PTScotchEngine,
                        string const &, dataRepository::Group * const )

} // namespace geos
