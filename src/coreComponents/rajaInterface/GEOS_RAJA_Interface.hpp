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

#ifndef GEOSX_RAJAINTERFACE_RAJAINTERFACE_HPP
#define GEOSX_RAJAINTERFACE_RAJAINTERFACE_HPP

// Source includes
#include "common/DataTypes.hpp"

// TPL includes
#include <RAJA/RAJA.hpp>

using serialPolicy = RAJA::loop_exec;
using serialReduce = RAJA::seq_reduce;
using serialAtomic = RAJA::seq_atomic;

#if defined(GEOSX_USE_OPENMP)

using parallelHostPolicy = RAJA::omp_parallel_for_exec;
using parallelHostReduce = RAJA::omp_reduce;
using parallelHostAtomic = RAJA::builtin_atomic;

#else

using parallelHostPolicy = serialPolicy;
using parallelHostReduce = serialReduce;
using parallelHostAtomic = serialAtomic;

#endif

#if defined(GEOSX_USE_CUDA)

template< int BLOCK_SIZE = 256 >
using parallelDevicePolicy = RAJA::cuda_exec< BLOCK_SIZE >;
using parallelDeviceReduce = RAJA::cuda_reduce;
using parallelDeviceAtomic = RAJA::cuda_atomic;

#else

template< int BLOCK_SIZE = 0 >
using parallelDevicePolicy = parallelHostPolicy;
using parallelDeviceReduce = parallelHostReduce;
using parallelDeviceAtomic = parallelHostAtomic;

#endif

namespace geosx
{

//RAJA wrapper for loops over ranges - local index
template< typename POLICY=serialPolicy, typename LAMBDA=void>
RAJA_INLINE void forall_in_range(const localIndex begin, const localIndex end, LAMBDA && body)
{
  RAJA::forall<POLICY>(RAJA::TypedRangeSegment<localIndex>(begin, end), std::forward<LAMBDA>(body));
}

template< typename POLICY, typename LAMBDA >
RAJA_INLINE void forAll( const localIndex end, LAMBDA && body )
{
  RAJA::forall< POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, end ), std::forward< LAMBDA >( body ) );
}

//RAJA wrapper for loops over ranges - global index
template<class POLICY=serialPolicy, typename LAMBDA=void>
RAJA_INLINE void forall_in_range(const globalIndex begin, const globalIndex end, LAMBDA && body)
{
  RAJA::forall<POLICY>(RAJA::TypedRangeSegment<globalIndex>(begin, end), std::forward<LAMBDA>(body));
}

//RAJA wrapper for loops over sets
template<class POLICY=serialPolicy, typename T, typename LAMBDA=void>
RAJA_INLINE void forall_in_set(const T * const indexList, const localIndex len, LAMBDA && body)
{
  RAJA::forall<POLICY>(RAJA::TypedListSegment<T>(indexList, len, RAJA::Unowned), std::forward<LAMBDA>(body));
}

}

#endif
