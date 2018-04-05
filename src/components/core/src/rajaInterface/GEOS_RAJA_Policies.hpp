#ifndef __GEOS_RAJA_POLICY__HPP
#define __GEOS_RAJA_POLICY__HPP

#include "RAJA/RAJA.hpp"
#include "RAJA/util/defines.hpp"
#include "RAJA/index/RangeSegment.hpp"

#define localdim 3

typedef RAJA::loop_exec elemPolicy;
typedef RAJA::loop_exec onePointPolicy;

typedef RAJA::loop_exec memSetPolicy;
typedef RAJA::loop_exec computeForcePolicy;

typedef RAJA::seq_exec quadraturePolicy;

#ifdef USE_OPENMP
typedef RAJA::atomic::omp_atomic atomicPolicy;
#else
typedef RAJA::atomic::loop_atomic atomicPolicy;

#endif

#endif
