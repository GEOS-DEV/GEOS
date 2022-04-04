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

#include "common/MpiWrapper.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include <parmetis.h>

#include <numeric>

namespace geosx
{
namespace parmetis
{

static_assert( std::is_same< idx_t, int64_t >::value, "PetMETIS must be built with IDXTYPEWIDTH=64" );

array1d< int64_t >
partMeshKway( ArrayOfArraysView< int64_t const, int64_t > const & elemToNodes,
              MPI_Comm comm )
{
  int const numProcs = MpiWrapper::commSize( comm );
  idx_t const numElems = elemToNodes.size();

  array1d< int64_t > part( numElems ); // all 0 by default
  if( numProcs > 1 )
  {
    // Compute `elemdist` parameter (element range owned by each rank)
    array1d< idx_t > elemCounts( numProcs + 1 );
    MpiWrapper::allGather( numElems, elemCounts, comm );
    array1d< idx_t > elemDist( numProcs + 1 );
    std::partial_sum( elemCounts.begin(), elemCounts.end(), elemDist.begin() + 1 );

    // Compute tpwgts parameters (target partition weights)
    array1d< real_t > tpwgts( numProcs );
    tpwgts.setValues< serialPolicy >( static_cast< real_t >( 1.0 / numProcs ) );

    // Prepare other ParMETIS parameters
    idx_t wgtflag = 0;
    idx_t numflag = 0;
    idx_t ncon = 1;
    idx_t ncommonnodes = 3; // enforce face connectivity
    idx_t npart = numProcs;
    idx_t options[3] = { 1, 0, 2022 };
    idx_t edgecut = 0;
    real_t ubvec = 1.05;

    // Call ParMETIS
    int const res = ParMETIS_V3_PartMeshKway( elemDist.data(),
                                              const_cast< idx_t * >( elemToNodes.getOffsets() ),
                                              const_cast< idx_t * >( elemToNodes.getValues() ),
                                              nullptr, &wgtflag, &numflag, &ncon, &ncommonnodes, &npart,
                                              tpwgts.data(), &ubvec, options, &edgecut, part.data(), &comm );
    GEOSX_ERROR_IF_NE_MSG( res, METIS_OK, "ParMETIS call returned an error: " << res );
  }

  return part;
}

} // namespace parmetis
} // namespace geosx
