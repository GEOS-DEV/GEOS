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
