/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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

typedef RAJA::loop_exec elemPolicy;
typedef RAJA::loop_exec onePointPolicy;

typedef RAJA::loop_exec memSetPolicy;
typedef RAJA::loop_exec computeForcePolicy;

typedef RAJA::seq_exec quadraturePolicy;
typedef RAJA::atomic::seq_atomic atomicPolicy;

typedef RAJA::loop_exec stencilPolicy;
typedef RAJA::seq_reduce reducePolicy;

typedef RAJA::loop_exec parallelHostPolicy;

typedef RAJA::loop_exec materialUpdatePolicy;

#define GEOSX_LAMBDA [=]

namespace geosx
{

//Alias to RAJA reduction operators
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

#define FOR_ALL_IN_SET( POLICY, INDICES_SET, INDEX ) \
    geosx::forall_in_set<POLICY>( INDICES_SET.data(), \
                                  INDICES_SET.size(), \
                                  GEOSX_LAMBDA ( INDEX )->void

#endif
