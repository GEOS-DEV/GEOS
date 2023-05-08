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
 * @file ScotchInterface.cpp
 */

#include "ScotchInterface.hpp"

#include "common/TimingMacros.hpp"

#include <scotch.h>

static_assert( std::is_same< SCOTCH_Num, int64_t >::value,
               "Non-matching index types. Scotch must be built with 64-bit indices." );

#define GEOSX_SCOTCH_CHECK( call ) \
  do { \
    auto const ierr = call; \
    GEOS_ERROR_IF_NE_MSG( ierr, 0, "Error in call to:\n" << #call ); \
  } while( false )

namespace geos
{
namespace scotch
{

void partition( CRSMatrixView< int64_t const, int64_t const, int64_t const > const & graph,
                LinearSolverParameters::Multiscale::Coarsening::Graph::Scotch const & params,
                int64_t const numPart,
                arrayView1d< int64_t > const & partition )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( params );

  // Special handling for a trivial case
  if( numPart <= 1 )
  {
    partition.zero();
    return;
  }

  SCOTCH_Graph * const gr = SCOTCH_graphAlloc();
  GEOSX_SCOTCH_CHECK( SCOTCH_graphInit( gr ) );

  SCOTCH_Num const numVerts = graph.numRows();
  SCOTCH_Num const numEdges = graph.getOffsets()[numVerts];

  // Technical UB if Scotch writes into these arrays; in practice we discard them right after
  SCOTCH_Num * const offsets = const_cast< SCOTCH_Num * >( graph.getOffsets() );
  SCOTCH_Num * const columns = const_cast< SCOTCH_Num * >( graph.getColumns() );
  SCOTCH_Num * const weights = const_cast< SCOTCH_Num * >( graph.getEntries() );

  GEOSX_SCOTCH_CHECK( SCOTCH_graphBuild( gr,          // graphptr
                                         0,           // baseval
                                         numVerts,    // vertnbr
                                         offsets,     // verttab
                                         nullptr,     // vendtab
                                         nullptr,     // velotab
                                         nullptr,     // vlbltab
                                         numEdges,    // edgenbr
                                         columns,     // edgetab
                                         weights      // edlotab,
                                         ) );

  // TODO: maybe remove?
  GEOSX_SCOTCH_CHECK( SCOTCH_graphCheck( gr ) );

  SCOTCH_Strat * const strategy = SCOTCH_stratAlloc();
  GEOSX_SCOTCH_CHECK( SCOTCH_stratInit( strategy ) );

  GEOSX_SCOTCH_CHECK( SCOTCH_graphPart( gr, numPart, strategy, partition.data() ) );

  SCOTCH_stratExit( strategy );
  SCOTCH_graphExit( gr );
}

} // namespace scotch
} // namespace geos
