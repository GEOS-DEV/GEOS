/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef __GEOS_RAJA_POLICY__HPP
#define __GEOS_RAJA_POLICY__HPP

#include "RAJA/RAJA.hpp"
#include "RAJA/util/defines.hpp"
#include "RAJA/index/RangeSegment.hpp"

#if defined(RAJA_ENABLE_CUDA)
typedef RAJA::cuda_exec<256> elemPolicy;
typedef RAJA::cuda_exec<256> onePointPolicy;

typedef RAJA::cuda_exec<256> memSetPolicy;
typedef RAJA::cuda_exec<256> computeForcePolicy;

typedef RAJA::cuda_exec<256> quadraturePolicy;
typedef RAJA::atomic::cuda_atomic atomicPolicy;

typedef RAJA::omp_parallel_for_exec parallelHostPolicy;

#elif defined(GEOSX_USE_OPENMP)
typedef RAJA::omp_parallel_for_exec elemPolicy;
typedef RAJA::omp_parallel_for_exec onePointPolicy;

typedef RAJA::omp_parallel_for_exec memSetPolicy;
typedef RAJA::omp_parallel_for_exec computeForcePolicy;

typedef RAJA::omp_parallel_for_exec quadraturePolicy;
typedef RAJA::atomic::omp_atomic atomicPolicy;

typedef RAJA::loop_exec stencilPolicy;
typedef RAJA::omp_reduce_ordered reducePolicy;

typedef RAJA::omp_parallel_for_exec parallelHostPolicy;

#else
typedef RAJA::seq_exec elemPolicy;
typedef RAJA::loop_exec onePointPolicy;

typedef RAJA::loop_exec memSetPolicy;
typedef RAJA::loop_exec computeForcePolicy;

typedef RAJA::seq_exec quadraturePolicy;
typedef RAJA::atomic::seq_atomic atomicPolicy;

typedef RAJA::loop_exec stencilPolicy;
typedef RAJA::seq_reduce reducePolicy;

typedef RAJA::loop_exec parallelHostPolicy;

#endif

#if defined(RAJA_ENABLE_CUDA)
#define GEOSX_LAMBDA [=] RAJA_DEVICE
#else
#define GEOSX_LAMBDA [&]
#endif

namespace geosx
{  

namespace raja
{
  
template< typename POLICY, typename T >
using ReduceSum = RAJA::ReduceSum<POLICY, T>;




//
template<typename POLICY=atomicPolicy, typename T>
RAJA_INLINE void atomicAdd(T *acc, T value)
{
  RAJA::atomic::atomicAdd<POLICY>(acc, value);
}

//RAJA wrapper for loops over ranges - local index
template<class POLICY=elemPolicy, typename LAMBDA=void>
RAJA_INLINE void forall_in_range(const localIndex begin, const localIndex end, LAMBDA && body)
{
  RAJA::forall<POLICY>(RAJA::TypedRangeSegment<localIndex>(begin, end), body);
}

//RAJA wrapper for loops over ranges - local index
template<class POLICY=elemPolicy, typename LAMBDA=void>
RAJA_INLINE void forall_in_range(const globalIndex begin, const globalIndex end, LAMBDA && body)
{
  RAJA::forall<POLICY>(RAJA::TypedRangeSegment<globalIndex>(begin, end), body);
}

//RAJA wrapper for loops over sets
template<class POLICY=elemPolicy, typename T, typename LAMBDA=void>
RAJA_INLINE void forall_in_set(const T * const indexList, const localIndex len, LAMBDA && body)
{
  RAJA::forall<POLICY>(RAJA::TypedListSegment<T>(indexList, len, RAJA::Unowned), body);
}

}
}
  
#endif
