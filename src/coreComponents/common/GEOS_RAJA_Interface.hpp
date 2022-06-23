/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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

#include <chrono>
#include <thread>

#define GEOSX_ASYNC_WAIT( UPPER, NANOSLEEP, TEST ) while( !TEST ) { }

namespace geosx
{

using serialPolicy = RAJA::loop_exec;
using serialReduce = RAJA::seq_reduce;
using serialAtomic = RAJA::seq_atomic;

using serialStream = RAJA::resources::Host;
using serialEvent = RAJA::resources::HostEvent;

#if defined(GEOSX_USE_OPENMP)

using parallelHostPolicy = RAJA::omp_parallel_for_exec;
using parallelHostReduce = RAJA::omp_reduce;
using parallelHostAtomic = RAJA::builtin_atomic;

// issues with Raja::resources::Omp on lassen
using parallelHostStream = serialStream;
using parallelHostEvent = serialEvent;

void RAJA_INLINE parallelHostSync() { RAJA::synchronize< RAJA::omp_synchronize >(); }

#else

using parallelHostPolicy = serialPolicy;
using parallelHostReduce = serialReduce;
using parallelHostAtomic = serialAtomic;
using parallelHostStream = serialStream;
using parallelHostEvent = serialEvent;

void RAJA_INLINE parallelHostSync() { }

#endif

#if defined(GEOSX_USE_CUDA)

template< unsigned long BLOCK_SIZE = 256 >
using parallelDevicePolicy = RAJA::cuda_exec< BLOCK_SIZE >;

template< unsigned long BLOCK_SIZE = 256 >
using parallelDeviceAsyncPolicy = RAJA::cuda_exec_async< BLOCK_SIZE >;

using parallelDeviceStream = RAJA::resources::Cuda;
using parallelDeviceEvent = RAJA::resources::Event;

using parallelDeviceReduce = RAJA::cuda_reduce;
using parallelDeviceAtomic = RAJA::cuda_atomic;

void RAJA_INLINE parallelDeviceSync() { RAJA::synchronize< RAJA::cuda_synchronize >(); }

template< typename POLICY, typename RESOURCE, typename LAMBDA >
RAJA_INLINE parallelDeviceEvent forAll( RESOURCE && stream, const localIndex end, LAMBDA && body )
{
  return RAJA::forall< POLICY >( std::forward< RESOURCE >( stream ),
                                 RAJA::TypedRangeSegment< localIndex >( 0, end ),
                                 std::forward< LAMBDA >( body ) );
}

#else

template< unsigned long BLOCK_SIZE = 0 >
using parallelDevicePolicy = parallelHostPolicy;

template< unsigned long BLOCK_SIZE = 0 >
using parallelDeviceAsyncPolicy = parallelHostPolicy;

using parallelDeviceStream = parallelHostStream;
using parallelDeviceEvent = parallelHostEvent;

using parallelDeviceReduce = parallelHostReduce;
using parallelDeviceAtomic = parallelHostAtomic;

void RAJA_INLINE parallelDeviceSync() { parallelHostSync( ); }

template< typename POLICY, typename RESOURCE, typename LAMBDA >
RAJA_INLINE parallelDeviceEvent forAll( RESOURCE && GEOSX_UNUSED_PARAM( stream ), const localIndex end, LAMBDA && body )
{
  RAJA::forall< POLICY >( RAJA::TypedRangeSegment< localIndex >( 0, end ), std::forward< LAMBDA >( body ) );
  return parallelDeviceEvent();
}

#endif

using parallelDeviceEvents = std::vector< parallelDeviceEvent >;

namespace internalRajaInterface
{
template< typename >
struct PolicyMap;

template<>
struct PolicyMap< serialPolicy >
{
  using atomic = serialAtomic;
  using reduce = serialReduce;
};

#if defined(GEOSX_USE_OPENMP)
template<>
struct PolicyMap< RAJA::omp_parallel_for_exec >
{
  using atomic = RAJA::builtin_atomic;
  using reduce = RAJA::omp_reduce;
};
#endif

#if defined(GEOSX_USE_CUDA)
template< unsigned long BLOCK_SIZE >
struct PolicyMap< RAJA::cuda_exec< BLOCK_SIZE > >
{
  using atomic = RAJA::cuda_atomic;
  using reduce = RAJA::cuda_reduce;
};
#endif
}


template< typename POLICY >
using ReducePolicy = typename internalRajaInterface::PolicyMap< POLICY >::reduce;

template< typename POLICY >
using AtomicPolicy = typename internalRajaInterface::PolicyMap< POLICY >::atomic;


RAJA_INLINE bool testAllDeviceEvents( parallelDeviceEvents & events )
{
  bool allDone = true;
  for( auto & event : events )
  {
    if( !event.check() )
    {
      allDone = false;
      break;
    }
  }
  return allDone;
}

RAJA_INLINE void waitAllDeviceEvents( parallelDeviceEvents & events )
{
  // poll device events for completion then wait 10 nanoseconds 6,000,000,000 times (60 sec timeout)
  // 10 nsecs ~= 30 clock cycles @ 3Ghz
  GEOSX_ASYNC_WAIT( 6000000000, 10, testAllDeviceEvents( events ) );
}

template< typename POLICY, typename INDEX, typename LAMBDA >
RAJA_INLINE void forAll( INDEX const end, LAMBDA && body )
{
  RAJA::forall< POLICY >( RAJA::TypedRangeSegment< INDEX >( 0, end ), std::forward< LAMBDA >( body ) );
}

template< typename POLICY, typename INDEX, typename LAMBDA >
RAJA_INLINE void forRange( INDEX const begin, INDEX const end, LAMBDA && body )
{
  RAJA::forall< POLICY >( RAJA::TypedRangeSegment< INDEX >( begin, end ), std::forward< LAMBDA >( body ) );
}

/*
RAJA Teams Policy setup
*/
using namespace RAJA::expt;
#if defined(RAJA_ENABLE_CUDA)
using team_launch_policy = LaunchPolicy<seq_launch_t, cuda_launch_t<true>>;
using team_x = LoopPolicy<RAJA::loop_exec, cuda_block_x_direct>;
using thread_z = LoopPolicy<RAJA::loop_exec, cuda_thread_z_loop>;
using thread_y = LoopPolicy<RAJA::loop_exec, cuda_thread_y_loop>;
using thread_x = LoopPolicy<RAJA::loop_exec, cuda_thread_x_loop>;
#define GEOSX_SHARED RAJA_TEAM_SHARED
#define GEOSX_STATIC_SHARED RAJA_TEAM_SHARED
#define GEOSX_THREAD_ID(k) threadIdx.k
#elif defined(RAJA_ENABLE_HIP)
using team_launch_policy = LaunchPolicy<seq_launch_t, hip_launch_t<true>>;
using team_x = LoopPolicy<RAJA::loop_exec, hip_block_x_direct>;
using thread_z = LoopPolicy<RAJA::loop_exec, hip_thread_z_loop>;
using thread_y = LoopPolicy<RAJA::loop_exec, hip_thread_y_loop>;
using thread_x = LoopPolicy<RAJA::loop_exec, hip_thread_x_loop>;
#define GEOSX_SHARED RAJA_TEAM_SHARED
#define GEOSX_STATIC_SHARED RAJA_TEAM_SHARED
#define GEOSX_THREAD_ID(k) hipThreadIdx_ ##k
#else
using team_launch_policy = LaunchPolicy<seq_launch_t>;
using team_x = LoopPolicy<RAJA::loop_exec>;
using thread_z = LoopPolicy<RAJA::loop_exec>;
using thread_y = LoopPolicy<RAJA::loop_exec>;
using thread_x = LoopPolicy<RAJA::loop_exec>;
#define GEOSX_SHARED
#define GEOSX_STATIC_SHARED static
#define GEOSX_THREAD_ID(k) 0
#endif

// #define GEOSX_FOREACH_THREAD(i,k,n) RAJA::expt::loop<team_##k> (ctx, RAJA::RangeSegment(0, n), [&] (int i)
// Hack to avoid having to use the LaunchContext
#ifdef __CUDA_ARCH__
#define GEOSX_FOREACH_THREAD(i,k,N) for(int i=threadIdx.k; i<N; i+=blockDim.k)
#elif defined(__HIP_DEVICE_COMPILE__)
#define GEOSX_FOREACH_THREAD(i,k,N) for(int i=hipThreadIdx_ ##k; i<N; i+=hipBlockDim_ ##k)
#else
#define GEOSX_FOREACH_THREAD(i,k,N) for(int i=0; i<N; i++)
#endif

//TODO
#define GEOSX_SYNC_THREAD

} // namespace geosx

#endif // GEOSX_RAJAINTERFACE_RAJAINTERFACE_HPP
