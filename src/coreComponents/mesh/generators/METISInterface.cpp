/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file METISInterface.cpp
 */

#include "METISInterface.hpp"

#include "common/TimingMacros.hpp"
#include "common/format/Format.hpp"

#include <metis.h>

#include <numeric>
#include <type_traits>

namespace geos
{
namespace metis
{

static_assert( std::is_same< idx_t, pmet_idx_t >::value,
               "Non-matching index types. METIS must be built with 64-bit indices." );

array1d< pmet_idx_t >
partitionWeighted( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                   arrayView1d< pmet_idx_t const > const & edgeWeights,
                   arrayView2d< pmet_idx_t const > const & vertexWeights,
                   pmet_idx_t const numParts,
                   arrayView1d< real64 const > const & imbalance,
                   pmet_idx_t const seed )
{
  GEOS_MARK_FUNCTION;

  pmet_idx_t const numVertices = graph.size();
  pmet_idx_t const numConstraints = vertexWeights.size( 1 );
  pmet_idx_t const numAdjacencies = numVertices > 0 ? graph.getOffsets()[numVertices] : 0;

  GEOS_ERROR_IF( numParts <= 0, "Number of METIS partitions must be strictly positive" );
  GEOS_ERROR_IF( numVertices <= 0, "Cannot partition an empty graph with METIS" );
  GEOS_ERROR_IF( numParts > numVertices,
                 GEOS_FMT( "Cannot create {} nonempty METIS partitions from {} vertices",
                           numParts, numVertices ) );
  GEOS_ERROR_IF_NE_MSG( vertexWeights.size( 0 ), numVertices,
                        "METIS vertex-weight row count must match the graph size" );
  GEOS_ERROR_IF( numConstraints <= 0, "METIS requires at least one balance constraint" );
  GEOS_ERROR_IF_NE_MSG( imbalance.size(), numConstraints,
                        "METIS requires one imbalance tolerance per constraint" );
  GEOS_ERROR_IF_NE_MSG( edgeWeights.size(), numAdjacencies,
                        "METIS edge weights must match the CSR adjacency array" );

  for( pmet_idx_t i = 0; i < numAdjacencies; ++i )
  {
    GEOS_ERROR_IF( edgeWeights[i] <= 0, "METIS edge weights must be positive" );
  }
  for( pmet_idx_t i = 0; i < numVertices; ++i )
  {
    for( pmet_idx_t k = 0; k < numConstraints; ++k )
    {
      GEOS_ERROR_IF( vertexWeights( i, k ) < 0, "METIS vertex weights must be nonnegative" );
    }
  }

  array1d< pmet_idx_t > parts( numVertices );
  if( numParts == 1 )
  {
    return parts;
  }

  array1d< real_t > imbalanceVector( numConstraints );
  for( pmet_idx_t k = 0; k < numConstraints; ++k )
  {
    GEOS_ERROR_IF( imbalance[k] < 0.0,
                   GEOS_FMT( "METIS imbalance tolerance {} is negative", imbalance[k] ) );
    imbalanceVector[k] = static_cast< real_t >( 1.0 + imbalance[k] );
  }

  idx_t options[METIS_NOPTIONS];
  int const optionsResult = METIS_SetDefaultOptions( options );
  GEOS_ERROR_IF_NE_MSG( optionsResult, METIS_OK, "METIS_SetDefaultOptions failed" );
  options[METIS_OPTION_SEED] = seed;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_CONTIG] = 0;

  idx_t numVerticesArg = numVertices;
  idx_t numConstraintsArg = numConstraints;
  idx_t numPartsArg = numParts;
  idx_t edgeCut = 0;

  idx_t dummyAdjacency = 0;
  idx_t dummyEdgeWeight = 1;
  idx_t * adjacencyData = numAdjacencies > 0
                          ? const_cast< idx_t * >( graph.getValues() )
                          : &dummyAdjacency;
  idx_t * edgeWeightData = numAdjacencies > 0
                           ? const_cast< idx_t * >( edgeWeights.data() )
                           : &dummyEdgeWeight;

  int const result = METIS_PartGraphKway(
    &numVerticesArg,
    &numConstraintsArg,
    const_cast< idx_t * >( graph.getOffsets() ),
    adjacencyData,
    const_cast< idx_t * >( vertexWeights.data() ),
    nullptr,
    edgeWeightData,
    &numPartsArg,
    nullptr,
    imbalanceVector.data(),
    options,
    &edgeCut,
    parts.data() );

  GEOS_ERROR_IF_NE_MSG( result, METIS_OK,
                        GEOS_FMT( "METIS_PartGraphKway failed with error code {}", result ) );
  return parts;
}

} // namespace metis
} // namespace geos
