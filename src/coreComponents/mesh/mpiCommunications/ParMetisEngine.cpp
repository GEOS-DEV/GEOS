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
 * @file ParMetisEngine.cpp
 */

#include "ParMetisEngine.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include <numeric>

#ifdef GEOS_USE_PARMETIS
extern "C"
{
#include <parmetis.h>
}
#endif

#define GEOS_PARMETIS_CHECK( call ) \
  do { \
    auto const ierr = call; \
    GEOS_ERROR_IF_NE_MSG( ierr, METIS_OK, "Error in call to:\n" << #call ); \
  } while( false )

namespace geos
{

using namespace dataRepository;
using camp::idx_t;

static_assert( std::is_same< idx_t, pmet_idx_t >::value,
               "Non-matching index types. ParMETIS must be configured with 64-bit indices." );

ParMetisEngine::ParMetisEngine( string const & name,
                                dataRepository::Group * const parent )
  : GraphPartitionEngine( name, parent ),
  m_numRefinements( 0 )
{
  registerWrapper( viewKeyStruct::numRefinementsString(), &m_numRefinements ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Number of ParMETIS refinement iterations. "
                    "Default: 0 (no refinement)." );
}

ParMetisEngine::~ParMetisEngine()
{}

array1d< pmet_idx_t >
ParMetisEngine::partition( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & graph,
                           arrayView1d< pmet_idx_t const > const & vertDist,
                           pmet_idx_t const numParts,
                           MPI_Comm comm )
{
  GEOS_MARK_FUNCTION;

#ifndef GEOS_USE_PARMETIS
  GEOS_UNUSED_VAR( graph, vertDist, numParts, comm );
  GEOS_ERROR( "GEOS was not built with ParMETIS support. "
              "Reconfigure with -DENABLE_PARMETIS=ON" );
  return array1d< pmet_idx_t >();
#else

  array1d< pmet_idx_t > part( graph.size() );

  // If only one partition is requested, no partitioning is needed.
  // All elements are assigned to rank 0.
  if( numParts == 1 )
  {
    part.setValues< serialPolicy >( 0 );
    return part;
  }

  // Set up ParMETIS parameters
  idx_t wgtflag = 0; // No weights on vertices or edges
  idx_t numflag = 0; // C-style numbering (starts at 0)
  idx_t ncon = 1;    // Number of constraints (1, for vertex count balance)
  idx_t npart = numParts;
  idx_t edgecut = 0;
  real_t ubvec = 1.05; // Imbalance tolerance
  idx_t options[3] = { 1, 0, 0 }; // Use default options, no seed

  // Set target partition weights for even distribution
  array1d< real_t > tpwgts( numParts );
  tpwgts.setValues< serialPolicy >( 1.0f / static_cast< real_t >( numParts ) );

  GEOS_LOG_RANK_0( GEOS_FMT( "ParMetisEngine: Partitioning {} local vertices into {} parts",
                             graph.size(), numParts ) );

  // ParMETIS has an unusual API that requires non-const pointers for read-only data.
  // We must cast away constness, which is technically UB but is how the library is used.
  GEOS_PARMETIS_CHECK( ParMETIS_V3_PartKway( const_cast< idx_t * >( vertDist.data() ),
                                             const_cast< idx_t * >( graph.getOffsets() ),
                                             const_cast< idx_t * >( graph.getValues() ),
                                             nullptr, nullptr, &wgtflag,
                                             &numflag, &ncon, &npart, tpwgts.data(),
                                             &ubvec, options, &edgecut, part.data(), &comm ) );

  // Perform refinements if requested.
  if( m_numRefinements > 0 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "ParMetisEngine: Performing {} refinements", m_numRefinements ) );
    GEOS_PARMETIS_CHECK( ParMETIS_V3_RefineKway( const_cast< idx_t * >( vertDist.data() ),
                                                 const_cast< idx_t * >( graph.getOffsets() ),
                                                 const_cast< idx_t * >( graph.getValues() ),
                                                 nullptr, nullptr, &wgtflag,
                                                 &numflag, &ncon, &npart, tpwgts.data(),
                                                 &ubvec, options, &edgecut, part.data(), &comm ) );
  }

  GEOS_LOG_RANK_0( GEOS_FMT( "ParMetisEngine: Partition complete, edge-cut = {}", edgecut ) );

  return part;

#endif // GEOS_USE_PARMETIS
}

ArrayOfArrays< pmet_idx_t, pmet_idx_t >
ParMetisEngine::meshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
                            arrayView1d< pmet_idx_t const > const & elemDist,
                            MPI_Comm comm,
                            int const minCommonNodes )
{
  return parmetisMeshToDual( elemToNodes, elemDist, comm, minCommonNodes );
}

ArrayOfArrays< pmet_idx_t, pmet_idx_t >
ParMetisEngine::parmetisMeshToDual( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
                                    arrayView1d< pmet_idx_t const > const & elemDist,
                                    MPI_Comm comm,
                                    int const minCommonNodes )
{
  GEOS_MARK_FUNCTION;

#ifndef GEOS_USE_PARMETIS
  GEOS_UNUSED_VAR( elemToNodes, elemDist, comm, minCommonNodes );
  GEOS_ERROR( "GEOS was not built with ParMETIS support. "
              "Reconfigure with -DENABLE_PARMETIS=ON" );
  return ArrayOfArrays< pmet_idx_t, pmet_idx_t >();
#else

  idx_t const numElems = elemToNodes.size();

  // `parmetis` awaits the arrays to be allocated as two continuous arrays: one for values, the other for offsets.
  // Our `ArrayOfArrays` allows to reserve some extra space for further element insertion,
  // but this is not compatible with what `parmetis` requires.
  GEOS_ASSERT_EQ_MSG( std::accumulate( elemToNodes.getSizes(), elemToNodes.getSizes() + numElems, 0 ),
                      elemToNodes.valueCapacity(),
                      "Internal error. The element to nodes mapping must be strictly allocated for compatibility with a third party library." );

  pmet_idx_t numflag = 0;
  pmet_idx_t ncommonnodes = minCommonNodes;
  pmet_idx_t * xadj;
  pmet_idx_t * adjncy;

  GEOS_LOG_RANK_0( GEOS_FMT( "ParMetisEngine: Converting mesh to dual graph ({} local elements, min common nodes = {})",
                             numElems, minCommonNodes ) );

  // Technical UB if ParMETIS writes into these arrays; in practice we discard them right after
  GEOS_PARMETIS_CHECK( ParMETIS_V3_Mesh2Dual( const_cast< pmet_idx_t * >( elemDist.data() ),
                                              const_cast< pmet_idx_t * >( elemToNodes.getOffsets() ),
                                              const_cast< pmet_idx_t * >( elemToNodes.getValues() ),
                                              &numflag, &ncommonnodes, &xadj, &adjncy, &comm ) );

  ArrayOfArrays< pmet_idx_t, pmet_idx_t > graph;
  graph.resizeFromOffsets( numElems, xadj );

  // There is no way to direct-copy values into ArrayOfArrays without UB (casting away const)
  forAll< parallelHostPolicy >( numElems, [xadj, adjncy, graph = graph.toView()]( localIndex const k )
  {
    graph.appendToArray( k, adjncy + xadj[k], adjncy + xadj[k+1] );
  } );

  METIS_Free( xadj );
  METIS_Free( adjncy );

  return graph;

#endif // GEOS_USE_PARMETIS
}

// Register in the GraphPartitionEngine catalog
REGISTER_CATALOG_ENTRY( GraphPartitionEngine, ParMetisEngine,
                        string const &, dataRepository::Group * const )

} // namespace geos
