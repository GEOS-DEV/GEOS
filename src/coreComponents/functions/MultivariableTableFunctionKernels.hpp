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
 * @file MultivariableTableFunctionKernels.hpp
 */

#ifndef GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTIONKERNELS_HPP_
#define GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTIONKERNELS_HPP_

namespace geosx
{

GEOSX_HOST_DEVICE integer getAxisIntervalIndexLowMult( real64 const axis_value,
                                                       real64 const axis_min,
                                                       real64 const axis_max,
                                                       real64 const axis_step,
                                                       real64 const axis_step_inv,
                                                       integer const axis_points,
                                                       // OUTPUT:
                                                       real64 & axis_low,
                                                       real64 & axis_mult )
{
  integer axis_interval_index = integer((axis_value - axis_min) * axis_step_inv );

  // check that axis_idx is within interpolation interval: valid axis_idx is between [0; axis_points - 2],
  // since there are axis_points-1 intervals along the axis
  if( axis_interval_index < 0 )
  {
    axis_interval_index = 0;
    if( axis_value < axis_min )
    {
      printf( "Interpolation warning: axis is out of limits (%lf; %lf) with value %lf, extrapolation is applied\n", axis_min, axis_max, axis_value );
    }
  }
  else if( axis_interval_index > (axis_points - 2))
  {
    axis_interval_index = axis_points - 2;
    if( axis_value > axis_max )
    {
      printf( "Interpolation warning: axis is out of limits (%lf; %lf) with value %lf, extrapolation is applied\n", axis_min, axis_max, axis_value );
    }
  }

  axis_low = axis_interval_index * axis_step + axis_min;
  axis_mult = (axis_value - *axis_low) * axis_step_inv;
  return axis_interval_index;
}

template< integer N_DIMS, integer N_OPS >
GEOSX_HOST_DEVICE void interpolate_point_with_derivatives( real64 const *axis_values,
                                                           real64 const *body_data,
                                                           real64 const *axis_low,
                                                           real64 const *axis_mult,
                                                           real64 const *axis_step_inv,
                                                           // OUTPUT:
                                                           real64 *interp_values, real64 *interp_derivs )
{
  static const integer N_VERTS = 1 << N_DIMS;
  integer pwr = N_VERTS / 2; // distance between high and low values
  real64 workspace[(2 * N_VERTS - 1) * N_OPS];

  // copy operator values for all vertices
  for( integer i = 0; i < N_VERTS * N_OPS; ++i )
  {
    workspace[i] = body_data[i];
  }

  for( integer i = 0; i < N_DIMS; ++i )
  {
    //printf ("i = %d, N_VERTS = %d, New offset: %d\n", i, N_VERTS, 2 * N_VERTS - (N_VERTS>>i));

    for( integer j = 0; j < pwr; ++j )
    {
      for( integer op = 0; op < N_OPS; ++op )
      {
        // update own derivative
        workspace[(2 * N_VERTS - (N_VERTS >> i) + j) * N_OPS + op] = (workspace[(j + pwr) * N_OPS + op] - workspace[j * N_OPS + op]) * axis_step_inv[i];
      }

      // update all dependent derivatives
      for( integer k = 0; k < i; k++ )
      {
        for( integer op = 0; op < N_OPS; ++op )
        {
          workspace[(2 * N_VERTS - (N_VERTS >> k) + j) * N_OPS + op] = workspace[(2 * N_VERTS - (N_VERTS >> k) + j) * N_OPS + op] + axis_mult[i] *
                                                                       (workspace[(2 * N_VERTS - (N_VERTS >> k) + j + pwr) * N_OPS + op] - workspace[(2 * N_VERTS - (N_VERTS >> k) + j) * N_OPS + op]);
        }
      }

      for( integer op = 0; op < N_OPS; ++op )
      {
        // interpolate value
        workspace[j * N_OPS + op] = workspace[j * N_OPS + op] + (axis_values[i] - axis_low[i]) * workspace[(2 * N_VERTS - (N_VERTS >> i) + j) * N_OPS + op];
      }
    }
    pwr /= 2;
  }
  for( integer op = 0; op < N_OPS; ++op )
  {
    interp_values[op] = workspace[op];
    for( integer i = 0; i < N_DIMS; ++i )
    {
      interp_derivs[op * N_DIMS + i] = workspace[(2 * N_VERTS - (N_VERTS >> i)) * N_OPS + op];
    }
  }
}



/**
 * @class KernelWrapper
 *
 *
 */

/**
 * @class MultivariableTableFunctionStaticKernel
 *
 * A class for multivariable piecewise interpolation with static storage
 * All functions are interpolated using the same uniformly discretized space
 *
 * @tparam NUM_DIMS number of dimensions
 * @tparam NUM_OPS number of interpolated functions
 */
template< integer NUM_DIMS, integer NUM_OPS >
class MultivariableTableFunctionStaticKernel
{
public:

  /// Compile time value for the number of table dimensions
  static constexpr integer numDims = NUM_DIMS;

  /// Compile time value for the number of table outputs (operators)
  static constexpr integer numOps = NUM_OPS;

  /**
   * @brief Constructor
   *
   * @param[in] coordinates space discretization parameters
   * @param[in] values values of operators at supporting points
   * @param[in] input points where interpolation is requested
   * @param[out] output values of interpolated operators
   * @param[out] outputDerivatives derivatives of interpolated operators w.r.t. input
   */
  MultivariableTableFunctionStaticKernel( ArrayOfArraysView< real64 const > const & coordinates,
                                          arrayView1d< real64 const > const & values,
                                          arrayView1d< real64 const > const & input,
                                          arrayView1d< real64 > const & output,
                                          arrayView1d< real64 > const & outputDerivatives ):
    m_coordinates( coordinates ),
    m_values( values ),
    m_input( input ),
    m_output( output ),
    m_outputDerivatives( outputDerivatives )
  {};

  /**
   * @brief Compute interpolation element-wise using host or device execution space
   *
   * @param ei Index of element
   */
  GEOSX_HOST_DEVICE
  void compute( localIndex const ei ) const
  {
    m_output[ei] = interpolateLinear( &m_input[ei * NUM_DIMS] );
  }

  /// maximum dimensions for the coordinates in the table
  static constexpr integer maxDimensions = 5;
  /**
   *
   * @brief Linear interpolation (to be replaced)
   *
   * @tparam IN_ARRAY array type of input coordinates
   * @param input input coordinates
   * @return interpolated value
   */
  template< typename IN_ARRAY >
  GEOSX_HOST_DEVICE
  real64
  interpolateLinear( IN_ARRAY const & input ) const
  {
    integer const numDimensions = LvArray::integerConversion< integer >( m_coordinates.size() );
    localIndex bounds[maxDimensions][2]{};
    real64 weights[maxDimensions][2]{};

    // Determine position, weights
    for( localIndex dim = 0; dim < numDimensions; ++dim )
    {
      arraySlice1d< real64 const > const coords = m_coordinates[dim];
      if( input[dim] <= coords[0] )
      {
        // Coordinate is to the left of this axis
        bounds[dim][0] = 0;
        bounds[dim][1] = 0;
        weights[dim][0] = 0.0;
        weights[dim][1] = 1.0;
      }
      else if( input[dim] >= coords[coords.size() - 1] )
      {
        // Coordinate is to the right of this axis
        bounds[dim][0] = coords.size() - 1;
        bounds[dim][1] = bounds[dim][0];
        weights[dim][0] = 1.0;
        weights[dim][1] = 0.0;
      }
      else
      {
        // Find the coordinate index
        ///TODO make this fast
        // Note: find uses a binary search...  If we assume coordinates are
        // evenly spaced, we can speed things up considerably
        auto const lower = LvArray::sortedArrayManipulation::find( coords.begin(), coords.size(), input[dim] );
        bounds[dim][1] = LvArray::integerConversion< localIndex >( lower );
        bounds[dim][0] = bounds[dim][1] - 1;

        real64 const dx = coords[bounds[dim][1]] - coords[bounds[dim][0]];
        weights[dim][0] = 1.0 - ( input[dim] - coords[bounds[dim][0]]) / dx;
        weights[dim][1] = 1.0 - weights[dim][0];
      }
    }

    // Calculate the result
    real64 value = 0.0;
    integer const numCorners = 1 << numDimensions;
    for( integer point = 0; point < numCorners; ++point )
    {
      // Find array index
      localIndex tableIndex = 0;
      localIndex stride = 1;
      for( integer dim = 0; dim < numDimensions; ++dim )
      {
        integer const corner = (point >> dim) & 1;
        tableIndex += bounds[dim][corner] * stride;
        stride *= m_coordinates.sizeOfArray( dim );
      }

      // Determine weighted value
      real64 cornerValue = m_values[tableIndex];
      for( integer dim = 0; dim < numDimensions; ++dim )
      {
        integer const corner = (point >> dim) & 1;
        cornerValue *= weights[dim][corner];
      }
      value += cornerValue;
    }
    return value;
  }

protected:

  // inputs

  /// An array of table axes
  ArrayOfArraysView< real64 const > m_coordinates;

  /// Table values (in fortran order)
  arrayView1d< real64 const > m_values;

  /// Input array view
  arrayView1d< real64 const > m_input;

  // outputs

  /// Output array view
  arrayView1d< real64 > m_output;

  /// Output derivative array view
  arrayView1d< real64 > m_outputDerivatives;

};

} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTIONKERNELS_HPP_ */
