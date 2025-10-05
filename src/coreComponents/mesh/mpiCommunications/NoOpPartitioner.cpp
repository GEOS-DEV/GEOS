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
 * @file NoOpPartitioner.cpp
 */

#include "NoOpPartitioner.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{

NoOpPartitioner::NoOpPartitioner( string const & name, Group * const parent ):
  GraphPartitioner( name, parent )
{}


array1d< pmet_idx_t > NoOpPartitioner::partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                                                  arrayView1d< pmet_idx_t const > const & vertDist,
                                                  pmet_idx_t const numParts,
                                                  MPI_Comm comm,
                                                  int const numRefinements )
{
  GEOS_UNUSED_VAR( graph );
  GEOS_UNUSED_VAR( numParts );
  GEOS_UNUSED_VAR( numRefinements );

  // Get the current MPI rank. This will be the target partition for all local elements
  int const i_domain = MpiWrapper::commRank( comm );

  // Determine the number of local elements on this rank.
  localIndex const numLocalElements = vertDist.size() > 1
                                      ? vertDist[i_domain + 1] - vertDist[i_domain]
                                      : graph.size(); // Fallback for single-rank/uninitialized vertDist

  // Create a partition array and assign all local elements to the current rank
  array1d< pmet_idx_t > part( numLocalElements );
  part.setValues< serialPolicy >( static_cast< pmet_idx_t >( i_domain ) );

  return part;
}


REGISTER_CATALOG_ENTRY( PartitionerBase, NoOpPartitioner, string const &, dataRepository::Group * const )

} // namespace geos
