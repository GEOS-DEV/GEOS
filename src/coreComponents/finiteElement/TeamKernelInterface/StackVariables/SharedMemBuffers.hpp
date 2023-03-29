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
 * @file SharedMemBuffers.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_SHARED_MEM_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_SHARED_MEM_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace stackVariables
{

template < localIndex buffer_size, localIndex num_buffers, localIndex batch_size >
struct SharedMemBuffers
{
  // TODO: use SharedTensor
  real64 (* shared_mem_buffers)[buffer_size];

  GEOSX_HOST_DEVICE
  SharedMemBuffers( RAJA::LaunchContext & ctx )
  {
    GEOSX_STATIC_SHARED real64 shared_buffers[batch_size][num_buffers][buffer_size];

    RAJA::loop<thread_z> (ctx, RAJA::RangeSegment(0, batch_size), [&] (localIndex batch_index) {
      shared_mem_buffers = (real64(*)[buffer_size])shared_buffers[batch_index];
    } );
  }

  GEOSX_HOST_DEVICE
  real64 (& operator[]( localIndex i ))[buffer_size]
  {
    return shared_mem_buffers[i];
  }
};

} // namespace stackVariables

} // namespace geosx

#endif // GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_SHARED_MEM_HPP_
