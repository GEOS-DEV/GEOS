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
 * @file KernelBase.hpp
 */

namespace geosx
{

namespace tensor
{

// remove
#define GEOSX_VERIFY_KERNEL(x,...)
#define GEOSX_ASSERT_KERNEL(x,...)
//end remove

// #define GEOSX_VERIFY_KERNEL(x,...)    \
//    if (!(x))                          \
//    {                                  \
//       GEOSX_ABORT_KERNEL(__VA_ARGS__) \
//    }                                  \

// #ifdef GEOSX_DEBUG
// #define GEOSX_ASSERT_KERNEL(x,...)    \
//    if (!(x))                          \
//    {                                  \
//       GEOSX_ABORT_KERNEL(__VA_ARGS__) \
//    }
// #else
// #define GEOSX_ASSERT_KERNEL(x,...)
// #endif

} // namespace tensor

} // namespace geosx