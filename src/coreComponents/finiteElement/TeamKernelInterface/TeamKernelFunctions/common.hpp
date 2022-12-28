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

/**
 * @brief A simple 3D tensor index structure.
 * 
 */
struct TensorIndex
{
  localIndex x, y, z;
};

/**
 * @brief The threading model used to compute on each element.
 * Specify how the computation is distributed over thread blocks.
 * 
 */
enum class ThreadingModel {
  Serial, //< Each thread treat one element ( or more when batch_size > 1 ).
  Distributed1D, //< Each element is treated by num_thread_1d threads.
  Distributed2D, //< Each element is treated by num_thread_1d^2 threads.
  Distributed3D //< Each element is treated by num_thread_1d^3 threads.
};

/**
 * @brief 3D iterator over the "thread local" indices.
 * 
 * @tparam StackVariables The type of the stack variables.
 * @tparam Lambda3D The type of the body stored in a lambda function.
 * @param stack The stack variables.
 * @param loop_bound_1 The loop bound of the first index.
 * @param loop_bound_2 The loop bound of the second index.
 * @param loop_bound_3 The loop bound of the third index.
 * @param lambda Thw body of the triple "for" loop.
 * 
 * @note Implementation for serial threading model.
 */
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
  #pragma unroll
  for (localIndex ind_3 = 0; ind_3 < loop_bound_3; ind_3++)
  {
    #pragma unroll
    for (localIndex ind_2 = 0; ind_2 < loop_bound_2; ind_2++)
    {
      #pragma unroll
      for (localIndex ind_1 = 0; ind_1 < loop_bound_1; ind_1++)
      {
        lambda( ind_1, ind_2, ind_3 );
      }
    }
  }
}

/**
 * @brief 3D iterator over the "thread local" indices.
 * 
 * @tparam StackVariables The type of the stack variables.
 * @tparam Lambda3D The type of the body stored in a lambda function.
 * @param stack The stack variables.
 * @param loop_bound_1 The loop bound of the first index.
 * @param loop_bound_2 The loop bound of the second index.
 * @param loop_bound_3 The loop bound of the third index.
 * @param lambda Thw body of the triple "for" loop.
 * 
 * @note Implementation for 1D distributed threading model.
 */
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
  localIndex ind_1 = stack.tidx;
  if ( ind_1 < loop_bound_1 )
  {
    #pragma unroll
    for (localIndex ind_3 = 0; ind_3 < loop_bound_3; ind_3++)
    {
      #pragma unroll
      for (localIndex ind_2 = 0; ind_2 < loop_bound_2; ind_2++)
      {
        lambda( ind_1, ind_2, ind_3 );
      }
    }
  }
}

/**
 * @brief 3D iterator over the "thread local" indices.
 * 
 * @tparam StackVariables The type of the stack variables.
 * @tparam Lambda3D The type of the body stored in a lambda function.
 * @param stack The stack variables.
 * @param loop_bound_1 The loop bound of the first index.
 * @param loop_bound_2 The loop bound of the second index.
 * @param loop_bound_3 The loop bound of the third index.
 * @param lambda Thw body of the triple "for" loop.
 * 
 * @note Implementation for 2D distributed threading model.
 */
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
  localIndex ind_1 = stack.tidx;
  localIndex ind_2 = stack.tidy;
  if ( ( ind_2 < loop_bound_2 ) && ( ind_1 < loop_bound_1 ) )
  {
    for (localIndex ind_3 = 0; ind_3 < loop_bound_3; ind_3++)
    {
      lambda( ind_1, ind_2, ind_3 );
    }
  }
}

/**
 * @brief 3D iterator over the "thread local" indices.
 * 
 * @tparam StackVariables The type of the stack variables.
 * @tparam Lambda3D The type of the body stored in a lambda function.
 * @param stack The stack variables.
 * @param loop_bound_1 The loop bound of the first index.
 * @param loop_bound_2 The loop bound of the second index.
 * @param loop_bound_3 The loop bound of the third index.
 * @param lambda Thw body of the triple "for" loop.
 * 
 * @note Implementation for 3D distributed threading model.
 */
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
  localIndex ind_1 = stack.tidx;
  localIndex ind_2 = stack.tidy;
  localIndex ind_3 = stack.tidz;
  if ( ( ind_3 < loop_bound_3 ) && ( ind_2 < loop_bound_2 ) && ( ind_1 < loop_bound_1 ) )
  {
    lambda( ind_1, ind_2, ind_3 );
  }
}

} // namespace geosx

#endif /* GEOSX_FINITEELEMENT_TEAMKERNELFUNCTION_COMMON_HPP_ */
