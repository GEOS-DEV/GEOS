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
 * @file common.hpp
 */

#ifndef GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_COMMON_HPP_
#define GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_COMMON_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

struct TensorIndex
{
  localIndex x, y, z;
};

enum class ThreadingModel {
  Serial,
  Distributed1D,
  Distributed2D,
  Distributed3D
};

template < typename StackVariables,
           typename Lambda3D,
           std::enable_if_t< StackVariables::threading_model == ThreadingModel::Serial, bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void loop3D( StackVariables & stack,
             localIndex loop_bound_1,
             localIndex loop_bound_2,
             localIndex loop_bound_3,
             Lambda3D && lambda )
{
  for (localIndex ind_3 = 0; ind_3 < loop_bound_3; ind_3++)
  {
    for (localIndex ind_2 = 0; ind_2 < loop_bound_2; ind_2++)
    {
      for (localIndex ind_1 = 0; ind_1 < loop_bound_1; ind_1++)
      {
        lambda( ind_1, ind_2, ind_3 );
      }
    }
  }
}

template < typename StackVariables,
           typename Lambda3D,
           std::enable_if_t< StackVariables::threading_model == ThreadingModel::Distributed1D, bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void loop3D( StackVariables & stack,
             localIndex loop_bound_1,
             localIndex loop_bound_2,
             localIndex loop_bound_3,
             Lambda3D && lambda )
{
  for (localIndex ind_3 = 0; ind_3 < loop_bound_3; ind_3++)
  {
    for (localIndex ind_2 = 0; ind_2 < loop_bound_2; ind_2++)
    {
      localIndex ind_1 = stack.tidx;
      if ( ind_1 < loop_bound_1 )
      {
        lambda( ind_1, ind_2, ind_3 );
      }
    }
  }
}

template < typename StackVariables,
           typename Lambda3D,
           std::enable_if_t< StackVariables::threading_model == ThreadingModel::Distributed2D, bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void loop3D( StackVariables & stack,
             localIndex loop_bound_1,
             localIndex loop_bound_2,
             localIndex loop_bound_3,
             Lambda3D && lambda )
{
  for (localIndex ind_3 = 0; ind_3 < loop_bound_3; ind_3++)
  {
    localIndex ind_2 = stack.tidy;
    if ( ind_2 < loop_bound_2 )
    {
      localIndex ind_1 = stack.tidx;
      if ( ind_1 < loop_bound_1 )
      {
        lambda( ind_1, ind_2, ind_3 );
      }
    }
  }
}

template < typename StackVariables,
           typename Lambda3D,
           std::enable_if_t< StackVariables::threading_model == ThreadingModel::Distributed3D, bool > = true >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void loop3D( StackVariables & stack,
             localIndex loop_bound_1,
             localIndex loop_bound_2,
             localIndex loop_bound_3,
             Lambda3D && lambda )
{
  // printf( "(e=%d,x=%d,y=%d,z=%d)\n", stack.element_index, stack.tidx, stack.tidy, stack.tidz );
  localIndex ind_3 = stack.tidz;
  if ( ind_3 < loop_bound_3 )
  {
    localIndex ind_2 = stack.tidy;
    if ( ind_2 < loop_bound_2 )
    {
      localIndex ind_1 = stack.tidx;
      if ( ind_1 < loop_bound_1 )
      {
        lambda( ind_1, ind_2, ind_3 );
      }
    }
  }
}

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_COMMON_HPP_ */
