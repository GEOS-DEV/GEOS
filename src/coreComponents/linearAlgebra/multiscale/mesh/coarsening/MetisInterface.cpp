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
 * @file MetisInterface.cpp
 */

#include "MetisInterface.hpp"

#include "common/TimingMacros.hpp"

#include <metis.h>

static_assert( std::is_same< idx_t, int64_t >::value,
               "Non-matching index types. METIS must be built with 64-bit indices." );

namespace geos
{
namespace metis
{

void partition( CRSMatrixView< int64_t const, int64_t const, int64_t const > const & graph,
                LinearSolverParameters::Multiscale::Coarsening::Graph::Metis const & params,
                int64_t const numPart,
                arrayView1d< int64_t > const & partition )
{
  GEOS_MARK_FUNCTION;

  // Special handling for a trivial case
  if( numPart <= 1 )
  {
    partition.zero();
    return;
  }

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions( options );
  options[METIS_OPTION_SEED] = params.seed;
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_CONTIG] = 1;
  options[METIS_OPTION_UFACTOR] = params.ufactor;

  // Must be non-const to comply with METIS API
  idx_t nnodes = graph.numRows();
  idx_t nconst = 1;
  idx_t objval = 0;
  idx_t nparts = LvArray::integerConversion< idx_t >( numPart );

  // Must cast away constness to comply with METIS API
  // This is not a problem since METIS won't modify the data
  idx_t * const offsets = const_cast< idx_t * >( graph.getOffsets() );
  idx_t * const columns = const_cast< idx_t * >( graph.getColumns() );
  idx_t * const weights = const_cast< idx_t * >( graph.getEntries() );

  auto const metisCall = [&]()
  {
    using Method = LinearSolverParameters::Multiscale::Coarsening::Graph::Metis::Method;
    switch( params.method )
    {
      case Method::kway:      return METIS_PartGraphKway;
      case Method::recursive: return METIS_PartGraphRecursive;
      default: GEOS_THROW( "Unrecognized METIS method", InputError );
    }
  }();

  int const result = metisCall( &nnodes, // nvtxs
                                &nconst, // ncon
                                offsets, // xadj
                                columns, // adjncy
                                nullptr, // vwgt
                                nullptr, // vsize
                                weights, // adjwgt
                                &nparts, // nparts
                                nullptr, // tpwgts
                                nullptr, // ubvec
                                options, // options
                                &objval, // edgecut
                                partition.data() ); // part
  GEOS_THROW_IF_NE_MSG( result, METIS_OK, "METIS call failed", std::runtime_error );
}

} // namespace metis
} // namespace geos
