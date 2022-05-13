/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParMETISInterface.cpp
 */

#include "ParMETISInterface.hpp"

#include "common/GEOS_RAJA_Interface.hpp"

#include <parmetis.h>

#include <numeric>

namespace geosx
{
namespace parmetis
{

static_assert( std::is_same< idx_t, pmet_idx_t >::value, "ParMETIS must be built with IDXTYPEWIDTH=64" );

array1d< pmet_idx_t >
partMeshKway( ArrayOfArraysView< pmet_idx_t const, pmet_idx_t > const & elemToNodes,
              MPI_Comm comm,
              int const minCommonNodes,
              int const numRefinements )
{
  int const numProcs = MpiWrapper::commSize( comm );
  idx_t const numElems = elemToNodes.size();

  array1d< pmet_idx_t > part( numElems ); // all 0 by default
  if( numProcs == 1 )
  {
    return part;
  }

  // Compute `elemdist` parameter (element range owned by each rank)
  array1d< idx_t > elemCounts;
  MpiWrapper::allGather( numElems, elemCounts, comm );
  array1d< idx_t > elemDist( numProcs + 1 );
  std::partial_sum( elemCounts.begin(), elemCounts.end(), elemDist.begin() + 1 );

  // Compute tpwgts parameters (target partition weights)
  array1d< real_t > tpwgts( numProcs );
  tpwgts.setValues< serialPolicy >( static_cast< real_t >( 1.0 / numProcs ) );

  // Set other ParMETIS parameters
  idx_t wgtflag = 0;
  idx_t numflag = 0;
  idx_t ncon = 1;
  idx_t ncommonnodes = minCommonNodes;
  idx_t npart = numProcs;
  idx_t options[4] = { 1, 0, 2022, PARMETIS_PSR_UNCOUPLED };
  idx_t edgecut = 0;
  real_t ubvec = 1.05;

  auto const checkReturnValue = []( int const r, char const * func )
  {
    GEOSX_ERROR_IF_NE_MSG( r, METIS_OK, GEOSX_FMT( "ParMETIS call '{}' returned an error code: {}", func, r ) );
  };

  idx_t * xadj;
  idx_t * adjncy;
  {
    int const res = ParMETIS_V3_Mesh2Dual( elemDist.data(),
                                           const_cast< idx_t * >( elemToNodes.getOffsets() ),
                                           const_cast< idx_t * >( elemToNodes.getValues() ),
                                           &numflag, &ncommonnodes,
                                           &xadj, &adjncy, &comm );
    checkReturnValue( res, "Mesh2Dual" );
  }

  // Call ParMETIS
  {
    int const res = ParMETIS_V3_PartKway( elemDist.data(), xadj, adjncy,
                                          nullptr, nullptr, &wgtflag,
                                          &numflag, &ncon, &npart, tpwgts.data(),
                                          &ubvec, options, &edgecut, part.data(), &comm );
    checkReturnValue( res, "PartKway" );
  }

  for( int iter = 0; iter < numRefinements; ++iter )
  {
    int const res = ParMETIS_V3_RefineKway( elemDist.data(), xadj, adjncy,
                                            nullptr, nullptr, &wgtflag,
                                            &numflag, &ncon, &npart, tpwgts.data(),
                                            &ubvec, options, &edgecut, part.data(), &comm );
    checkReturnValue( res, "RefineKway" );
  }

  METIS_Free( xadj );
  METIS_Free( adjncy );

  return part;
}

} // namespace parmetis
} // namespace geosx
