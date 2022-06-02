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
 * @file load.hpp
 */

#ifndef GEOSX_TENSOR_UTIL_LOAD
#define GEOSX_TENSOR_UTIL_LOAD

#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace tensor
{

/// Load row major values using 2d threads into a matrix of size nx x ny.
template <typename Matrix> GEOSX_HOST_DEVICE
void load_with_2dthreads(const double *values, int nx, int ny, Matrix &mat)
{
   GEOSX_FOREACH_THREAD(y,y,ny)
   {
      GEOSX_FOREACH_THREAD(x,x,nx)
      {
         mat(x,y) = values[x+nx*y];
      }
   }
   GEOSX_SYNC_THREAD;
}

} // namespace tensor

} // geosx namespace

#endif // GEOSX_TENSOR_UTIL_LOAD
