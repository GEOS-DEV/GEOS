// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
#ifndef __GEOSX_RAJA_WRAPPER__HPP
#define __GEOSX_RAJA_WRAPPER__HPP

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

#elif defined(RAJA_ENABLE_OPENMP)
typedef RAJA::omp_parallel_for_exec elemPolicy;
typedef RAJA::omp_parallel_for_exec onePointPolicy;

typedef RAJA::omp_parallel_for_exec memSetPolicy;
typedef RAJA::omp_parallel_for_exec computeForcePolicy;

typedef RAJA::omp_parallel_for_exec quadraturePolicy;
typedef RAJA::atomic::omp_atomic atomicPolicy;
#else
typedef RAJA::loop_exec elemPolicy;
typedef RAJA::loop_exec onePointPolicy;

typedef RAJA::loop_exec memSetPolicy;
typedef RAJA::loop_exec computeForcePolicy;

typedef RAJA::seq_exec quadraturePolicy;
typedef RAJA::atomic::loop_atomic atomicPolicy;
#endif

#if defined(RAJA_ENABLE_CUDA)
#define GEOSX_LAMBDA [=] RAJA_DEVICE
#else
#define GEOSX_LAMBDA [=]
#endif

using localIndex = RAJA::Index_type;

//RAJA wrapper of loops over ranges
template<class POLICY=elemPolicy, typename LAMBDA=void>
void forall_in_range(const localIndex begin, const localIndex end, LAMBDA && body){

  RAJA::forall<POLICY>(RAJA::RangeSegment(begin, end), body);  
}

//RAJA wrapper for loops over sets
//A RAJA list segment won't own the data
template<typename T, class POLICY=elemPolicy, typename LAMBDA=void>
void forall_in_set(const T * const indexList, const localIndex len, LAMBDA && body){

  RAJA::forall<POLICY>(RAJA::ListSegment(indexList, len, RAJA::Unowned), body);
}

#endif
