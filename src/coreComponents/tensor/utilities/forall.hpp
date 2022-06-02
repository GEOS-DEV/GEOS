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


/**
 * @file forall.hpp
 */

#ifndef GEOSX_TENSOR_FORALL
#define GEOSX_TENSOR_FORALL

#include "config.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace tensor
{


/// Macro to launch kernels that uses the config to choose the appropriate
/// threading strategy.
#define GEOSX_FORALL_CONFIG(config,i,N,...)                                    \
   const int threads = config.quads;                                           \
   using Config = decltype(config);                                            \
   constexpr int Bsize = get_config_batchsize<Config>;                         \
   const int Xthreads = config_use_xthreads<Config> ? threads : 1;             \
   const int Ythreads = config_use_ythreads<Config> ? threads : 1;             \
   const int Zthreads = Bsize * ( config_use_zthreads<Config> ? threads : 1 ); \
   using namespace RAJA::expt;                                                 \
   using RAJA::RangeSegment;                                                   \
   launch<POLICY>                                                              \
   (DEVICE, Grid(Teams(GRID), Threads(Xthreads, Ythreads, Zthreads)),          \
    [=] RAJA_DEVICE (LaunchContext ctx)                                        \
   {                                                                           \
      loop<teams_x>(ctx, RangeSegment(0, N), d_body);                          \
   });                                                                         \
   GEOSX_DEVICE_SYNC;
   // ForallWrap<3>(true,N,                                                       \
   //               [=] GEOSX_DEVICE (int i) mutable {__VA_ARGS__},               \
   //               [&] GEOSX_LAMBDA (int i) {__VA_ARGS__},                       \
   //               Xthreads, Ythreads, Zthreads);                                \

} // namespace tensor

} // namespace geosx

#endif // GEOSX_TENSOR_FORALL
