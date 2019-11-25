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



/* UNCRUSTIFY-OFF */

#ifndef GEOSX_COMMON_TIMINGMACROS_HPP_
#define GEOSX_COMMON_TIMINGMACROS_HPP_

#include "common/GeosxConfig.hpp"

#ifdef GEOSX_USE_CALIPER
#include <caliper/cali.h>
#include <sys/time.h>
#include <string>
#include <iostream>

namespace timingHelpers
{
  inline std::string stripPF( char const * prettyFunction )
  {
    std::string input(prettyFunction);
    std::string::size_type const end = input.find_first_of( '(' );
    std::string::size_type const beg = input.find_last_of( ' ', end)+1;
    return input.substr( beg, end-beg );
  }
}

#define DO_STRINGIFY(arg) #arg
#define GEOSX_MARK_LOOP_BEGIN(loop, loopName) CALI_CXX_MARK_LOOP_BEGIN(loop,DO_STRINGIFY(loopName))
#define GEOSX_MARK_LOOP_END(loop) CALI_CXX_MARK_LOOP_END(loop)
#define GEOSX_MARK_LOOP_ITERATION(loop, iter) CALI_CXX_MARK_LOOP_ITERATION(loop, iter)

#define GEOSX_MARK_FUNCTION_TAG(name) cali::Function __cali_ann##__LINE__(DO_STRINGIFY(name))

#define GEOSX_MARK_FUNCTION_SCOPED cali::Function __cali_ann##__func__(timingHelpers::stripPF(__PRETTY_FUNCTION__).c_str())

//#define GEOSX_MARK_FUNCTION CALI_CXX_MARK_FUNCTION
#define GEOSX_MARK_FUNCTION GEOSX_MARK_FUNCTION_SCOPED

#define GEOSX_MARK_BEGIN(name) CALI_MARK_BEGIN(DO_STRINGIFY(name))
#define GEOSX_MARK_END(name) CALI_MARK_END(DO_STRINGIFY(name))

#else
#define GEOSX_MARK_FUNCTION_TAG(name)
#define GEOSX_MARK_FUNCTION_SCOPED
#define GEOSX_MARK_FUNCTION

#define GEOSX_MARK_LOOP_BEGIN(loop, loopName)
#define GEOSX_MARK_LOOP_END(loop)
#define GEOSX_MARK_LOOP_ITERATION(loop, iter)
#define GEOSX_MARK_BEGIN(name)
#define GEOSX_MARK_END(name)
#endif

#ifdef GEOSX_USE_TIMERS
#define GEOSX_GET_TIME( time )                                                 \
  real64 time;                                                                 \
  do                                                                           \
  {                                                                            \
    timeval tim;                                                               \
    gettimeofday(&tim, nullptr);                                               \
    time = tim.tv_sec + (tim.tv_usec / 1000000.0);                             \
  } while (false)
#else
#define GEOSX_GET_TIME( time )
#endif



#endif
