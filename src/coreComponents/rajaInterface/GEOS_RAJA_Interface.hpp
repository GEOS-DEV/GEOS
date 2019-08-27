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
#include "common/DataTypes.hpp"


using serialPolicy = RAJA::loop_exec;

#if defined(GEOSX_USE_OPENMP)

using parallelHostPolicy = RAJA::omp_parallel_for_exec;

#else

using parallelHostPolicy = RAJA::loop_exec;

#endif

#if defined(GEOSX_USE_CUDA)

template< int BLOCK_SIZE = 256 >
using parallelDevicePolicy = RAJA::cuda_exec< BLOCK_SIZE >;

#elif defined(GEOSX_USE_OPENMP)

template< int BLOCK_SIZE = 0 >
using parallelDevicePolicy = RAJA::omp_parallel_for_exec;

#else

template< int BLOCK_SIZE = 0 >
using parallelDevicePolicy = RAJA::loop_exec;

#endif

namespace geosx
{

//RAJA wrapper for loops over ranges - local index
template<class POLICY=serialPolicy, typename LAMBDA=void>
RAJA_INLINE void forall_in_range(const localIndex begin, const localIndex end, LAMBDA && body)
{
  RAJA::forall<POLICY>(RAJA::TypedRangeSegment<localIndex>(begin, end), std::forward<LAMBDA>(body));
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

template< typename T , typename atomicPol=RAJA::atomic::auto_atomic>
inline void AddLocalToGlobal( arraySlice1d<localIndex const> const & globalToLocalRelation,
                              arraySlice1d< T const > const & localField,
                              arraySlice1d< T const >& globalField,
                              localIndex const N )
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    RAJA::atomic::atomicAdd<atomicPol>( &globalField[ globalToLocalRelation[a] ], localField[a] );
  }
}

template< typename atomicPol=RAJA::atomic::auto_atomic>
inline void AddLocalToGlobal( arraySlice1d<localIndex const> const & globalToLocalRelation,
                              arraySlice1d<R1Tensor const> const & localField,
                              arraySlice1d<R1Tensor>& globalField,
                              localIndex const N )
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    real64 * __restrict__ const gData = globalField[globalToLocalRelation[a]].Data();
    real64 const * __restrict__ const lData = localField[a].Data();
    RAJA::atomic::atomicAdd<atomicPol>( &gData[0], lData[0] );
    RAJA::atomic::atomicAdd<atomicPol>( &gData[1], lData[1] );
    RAJA::atomic::atomicAdd<atomicPol>( &gData[2], lData[2] );
  }
}

template< localIndex N, typename atomicPol=RAJA::atomic::auto_atomic>
inline void AddLocalToGlobal( arraySlice1d<localIndex const> const & globalToLocalRelation,
                              R1Tensor const * const restrict localField,
                              arraySlice1d<R1Tensor> & globalField )
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    real64 * __restrict__ const gData = globalField[globalToLocalRelation[a]].Data();
    real64 const * __restrict__ const lData = localField[a].Data();
    RAJA::atomic::atomicAdd<atomicPol>( &gData[0], lData[0] );
    RAJA::atomic::atomicAdd<atomicPol>( &gData[1], lData[1] );
    RAJA::atomic::atomicAdd<atomicPol>( &gData[2], lData[2] );
  }
}

template< typename T, typename atomicPol=RAJA::atomic::auto_atomic >
inline void AddLocalToGlobal( arraySlice1d<localIndex const> const & globalToLocalRelation,
                              arraySlice1d< T const > const & localField1,
                              arraySlice1d< T const > const & localField2,
                              arraySlice1d< T > & globalField1,
                              arraySlice1d< T > & globalField2,
                              localIndex const N )
{
  for( localIndex a=0 ; a<N ; ++a )
  {
    RAJA::atomic::atomicAdd<atomicPol>( &globalField1[ globalToLocalRelation[a] ], localField1[a] );
    RAJA::atomic::atomicAdd<atomicPol>( &globalField2[ globalToLocalRelation[a] ], localField2[a] );
  }
}

}

#endif
