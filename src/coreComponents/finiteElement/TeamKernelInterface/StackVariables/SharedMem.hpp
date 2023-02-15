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
 * @file SharedMem.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_SHARED_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_SHARED_HPP_

#include "common/GEOS_RAJA_Interface.hpp"
#include "finiteElement/TeamKernelInterface/StackVariables/Basis.hpp"

namespace geosx
{

namespace stackVariables
{

template < localIndex size_1d, localIndex batch_size >
struct SharedMem
{
  GEOSX_HOST_DEVICE
  SharedMem( LaunchContext & ctx )
  {
    // Element primary field gradients at quadrature points
    GEOSX_STATIC_SHARED real64 s_values[batch_size][size_1d][size_1d][size_1d];

    loop<thread_z> (ctx, RAJA::RangeSegment(0, batch_size), [&] (localIndex batch_index) {
      values = &s_values[batch_index];
    } );
  }

  // Values at quadrature points
  real64 ( * values )[size_1d][size_1d][size_1d]; // Can be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getValues() const )[size_1d][size_1d][size_1d]
  {
    return *values;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getValues() )[size_1d][size_1d][size_1d]
  {
    return *values;
  }
};

template < localIndex size_1d, localIndex dim, localIndex batch_size >
struct Shared1DMem
{
  GEOSX_HOST_DEVICE
  Shared1DMem( LaunchContext & ctx )
  {
    // Element primary field gradients at quadrature points
    GEOSX_STATIC_SHARED real64 s_values[batch_size][size_1d][size_1d][size_1d][dim];

    loop<thread_z> (ctx, RAJA::RangeSegment(0, batch_size), [&] (localIndex batch_index) {
      values = &s_values[batch_index];
    } );
  }

  // Values at quadrature points
  real64 ( * values )[size_1d][size_1d][size_1d][dim]; // Can be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getValues() const )[size_1d][size_1d][size_1d][dim]
  {
    return *values;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getValues() )[size_1d][size_1d][size_1d][dim]
  {
    return *values;
  }
};

template < localIndex size_1d, localIndex dim, localIndex batch_size >
struct Shared2DMem
{
  GEOSX_HOST_DEVICE
  Shared2DMem( LaunchContext & ctx )
  {
    // Element primary field gradients at quadrature points
    GEOSX_STATIC_SHARED real64 s_values[batch_size][size_1d][size_1d][size_1d][dim][dim];

    loop<thread_z> (ctx, RAJA::RangeSegment(0, batch_size), [&] (localIndex batch_index) {
      values = &s_values[batch_index];
    } );
  }

  // Values at quadrature points
  real64 ( * values )[size_1d][size_1d][size_1d][dim][dim]; // Can be in registers

  GEOSX_HOST_DEVICE
  real64 const ( & getValues() const )[size_1d][size_1d][size_1d][dim][dim]
  {
    return *values;
  }

  GEOSX_HOST_DEVICE
  real64 ( & getValues() )[size_1d][size_1d][size_1d][dim][dim]
  {
    return *values;
  }
};

} // namespace stackVariables

} // namespace geosx

#endif // GEOSX_FINITEELEMENT_TEAMKERNELBASE_STACKVARIABLES_SHARED_HPP_
