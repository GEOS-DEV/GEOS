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
 * @file NoOpEngine.cpp
 */

#include "NoOpEngine.hpp"
#include "common/MpiWrapper.hpp"

#ifdef GEOS_USE_PARMETIS
  #include "ParMetisEngine.hpp"  // For fallback
#endif

namespace geos
{

NoOpEngine::NoOpEngine( string const & name, Group * const parent )
  : GraphPartitionEngine( name, parent )
{}

array1d< pmet_idx_t >
NoOpEngine::partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                       arrayView1d< pmet_idx_t const > const & vertDist,
                       pmet_idx_t const numParts,
                       MPI_Comm comm )
{
  GEOS_UNUSED_VAR( graph );
  GEOS_UNUSED_VAR( numParts );

  int const myRank = MpiWrapper::commRank( comm );

  // Determine the number of local elements on this rank
  localIndex const numLocalElements = vertDist.size() > 1
                                      ? vertDist[myRank + 1] - vertDist[myRank]
                                      : graph.size(); // Fallback for single-rank or uninitialized vertDist

  // Create partition array: assign all local elements to current rank
  array1d< pmet_idx_t > partition( numLocalElements );
  partition.setValues< serialPolicy >( static_cast< pmet_idx_t >( myRank ) );

  return partition;
}

ArrayOfArrays< pmet_idx_t, pmet_idx_t >
NoOpEngine::meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
                        arrayView1d< pmet_idx_t const > const & elemDist,
                        MPI_Comm comm,
                        int minCommonNodes )
{
#ifdef GEOS_USE_PARMETIS
  // Fallback to ParMETIS implementation
  return ParMetisEngine::parmetisMeshToDual( elemToNodes, elemDist, comm, minCommonNodes );
#else

  GEOS_UNUSED_VAR( elemToNodes );
  GEOS_UNUSED_VAR( elemDist );
  GEOS_UNUSED_VAR( comm );
  GEOS_UNUSED_VAR( minCommonNodes );

  GEOS_ERROR( "NoOpEngine::meshToDual requires ParMETIS for mesh-to-dual conversion. "
              "Either build with GEOS_USE_PARMETIS=ON or provide the dual graph directly." );
  return ArrayOfArrays< pmet_idx_t, pmet_idx_t >();

#endif
}

// Register in the GraphPartitionEngine catalog
REGISTER_CATALOG_ENTRY( GraphPartitionEngine,
                        NoOpEngine,
                        string const &,
                        dataRepository::Group * const )

} // namespace geos
