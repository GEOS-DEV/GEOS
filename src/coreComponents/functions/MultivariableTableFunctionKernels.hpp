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
 * @tparam NUM_DIMS number of dimensions (inputs)
 * @tparam NUM_OPS number of interpolated functions (outputs)
 */
template< integer NUM_DIMS, integer NUM_OPS >
class MultivariableTableFunctionStaticKernel
{
public:

  /// Compile time value for the number of table dimensions (inputs)
  static constexpr integer numDims = NUM_DIMS;

  /// Compile time value for the number of operators (interpolated functions, outputs)
  static constexpr integer numOps = NUM_OPS;

  /**
   * @brief Construct a new Multivariable Table Function Static Kernel object
   *
   * @param[in] axisMinimums  minimum coordinate for each axis
   * @param[in] axisMaximums maximum coordinate for each axis
   * @param[in] axisPoints number of discretization points between minimum and maximum for each axis
   * @param[in] axisSteps axis interval lengths (axes are discretized uniformly)
   * @param[in] axisStepInvs inversions of axis interval lengths (axes are discretized uniformly)
   * @param[in] axisHypercubeMults  hypercube index mult factors for each axis
   * @param[in] hypercubeData table data stored per hypercube
   * @param[in] coordinates array of coordinates of points where interpolation is required
   * @param[in] values array of interpolated operator values (all operators are interpolated at each point)
   * @param[in] derivatives  array of derivatives of interpolated operators (all derivatives for all operators are computed at each point)
   */
  MultivariableTableFunctionStaticKernel( arrayView1d< real64 const > const & axisMinimums,
                                          arrayView1d< real64 const > const & axisMaximums,
                                          arrayView1d< integer const > const & axisPoints,
                                          arrayView1d< real64 const > const & axisSteps,
                                          arrayView1d< real64 const > const & axisStepInvs,
                                          arrayView1d< globalIndex const > const & axisHypercubeMults,
                                          arrayView1d< real64 const > const & hypercubeData,
                                          arrayView1d< real64 const > const & coordinates,
                                          arrayView1d< real64 > const & values,
                                          arrayView1d< real64 > const & derivatives ):
    m_axisMinimums ( axisMinimums ),
    m_axisMaximums ( axisMaximums ),
    m_axisPoints ( axisPoints ),
    m_axisSteps ( axisSteps ),
    m_axisStepInvs ( axisStepInvs ),
    m_axisHypercubeMults ( axisHypercubeMults ),
    m_hypercubeData ( hypercubeData ),
    m_coordinates ( coordinates ),
    m_values ( values ),
    m_derivatives ( derivatives )
  {};

  /**
   * @brief Compute interpolation element-wise using host or device execution space
   *
   * @param[in] ei Index of element
   */
  GEOSX_HOST_DEVICE
  void
  compute( localIndex const ei ) const
  {
    interpolatePoint( &m_coordinates[ei * NUM_DIMS], &m_values[ei * NUM_OPS], &m_derivatives[ei * NUM_OPS * NUM_DIMS] );
  }

/**
 * @brief interpolate operators at a given point
 *
 * @param[in] coordinates point coordinates
 * @param[out] values interpolated operator values
 * @param[out] derivatives derivatives of interpolated operators
 */
  GEOSX_HOST_DEVICE
  void
  interpolatePoint( real64 const *coordinates,
                    real64 *values,
                    real64 *derivatives ) const
  {
    globalIndex hypercubeIndex = 0;
    real64 axisLows[NUM_DIMS];
    real64 axisMults[NUM_DIMS];

    for( int i = 0; i < NUM_DIMS; ++i )
    {
      integer axisIndex = getAxisIntervalIndexLowMult( coordinates[i],
                                                       m_axisMinimums[i], m_axisMaximums[i],
                                                       m_axisSteps[i], m_axisStepInvs[i], m_axisPoints[i],
                                                       axisLows[i], axisMults[i] );
      hypercubeIndex += axisIndex * m_axisHypercubeMults[i];
    }

    interpolatePointWithDerivatives( coordinates,
                                     getHypercubeData( hypercubeIndex ),
                                     &axisLows[0], &axisMults[0],
                                     &m_axisStepInvs[0],
                                     values,
                                     derivatives );

  }

protected:

  /**
   * @brief Get pointer to hypercube data
   *
   * @param[in] hypercubeIndex
   * @return pointer to hypercube data
   */
  GEOSX_HOST_DEVICE
  inline
  real64 const *
  getHypercubeData( globalIndex const hypercubeIndex ) const
  {
    static constexpr integer NUM_VERTS = 1 << NUM_DIMS;
    return &m_hypercubeData[hypercubeIndex * NUM_VERTS * NUM_OPS];
  }

  /**
   * @brief Get the interval index, low and mult values for a given axis coordinate
   *
   * @param[in] axisCoordinate coordinate on a given axis
   * @param[in] axisMin minimum value on a given axis
   * @param[in] axisMax maximum value on a given axis
   * @param[in] axisStep interval length for a given axis
   * @param[in] axisStepInv inversion of the interval length for a given axis
   * @param[in] axisPoints number of discretization points for a given axis
   * @param[out] axisLow left coordinate of target axis interval
   * @param[out] axisMult weight of the right coordinate of target axis interval
   * @return integer target axis interval index
   */
  GEOSX_HOST_DEVICE
  inline
  integer
  getAxisIntervalIndexLowMult( real64 const axisCoordinate,
                               real64 const axisMin,
                               real64 const axisMax,
                               real64 const axisStep,
                               real64 const axisStepInv,
                               integer const axisPoints,
                               real64 & axisLow,
                               real64 & axisMult ) const
  {
    integer axisIntervalIndex = integer((axisCoordinate - axisMin) * axisStepInv );

    // check that axisindex is within interpolation interval: valid axisindex is between 0 and (axisPoints - 2),
    // since there are axisPoints-1 intervals along the axis

    if( axisIntervalIndex < 0 )
    {
      axisIntervalIndex = 0;
      if( axisCoordinate < axisMin )
      {
        printf( "Interpolation warning: axis is out of limits (%lf; %lf) with value %lf, extrapolation is applied\n", axisMin, axisMax, axisCoordinate );
      }
    }
    else if( axisIntervalIndex > (axisPoints - 2))
    {
      axisIntervalIndex = axisPoints - 2;
      if( axisCoordinate > axisMax )
      {
        printf( "Interpolation warning: axis is out of limits (%lf; %lf) with value %lf, extrapolation is applied\n", axisMin, axisMax, axisCoordinate );
      }
    }

    axisLow = axisIntervalIndex * axisStep + axisMin;
    axisMult = (axisCoordinate - axisLow) * axisStepInv;
    return axisIntervalIndex;
  }


  /**
   * @brief interpolate all operator values and derivatives at a given point
   *
   * @param[in] axisCoordinates coordinates of a point
   * @param[in] hypercubeData data of target hypercube
   * @param[in] axisLows array of left coordinates of target axis intervals
   * @param[in] axisMults array of weights of right coordinates of target axis intervals
   * @param[in] axisStepInvs array of inversions of axis steps
   * @param[out] values interpolated operator values
   * @param[out] derivatives derivatives of interpolated operators
   */
  GEOSX_HOST_DEVICE
  inline
  void
  interpolatePointWithDerivatives( real64 const * const axisCoordinates,
                                   real64 const * const hypercubeData,
                                   real64 const * const axisLows,
                                   real64 const * const axisMults,
                                   real64 const * const axisStepInvs,
                                   real64 *values,
                                   real64 *derivatives ) const
  {
    static constexpr integer NUM_VERTS = 1 << NUM_DIMS;
    integer pwr = NUM_VERTS / 2;   // distance between high and low values
    real64 workspace[(2 * NUM_VERTS - 1) * NUM_OPS];

    // copy operator values for all vertices
    for( integer i = 0; i < NUM_VERTS * NUM_OPS; ++i )
    {
      workspace[i] = hypercubeData[i];
    }

    for( integer i = 0; i < NUM_DIMS; ++i )
    {

      for( integer j = 0; j < pwr; ++j )
      {
        for( integer op = 0; op < NUM_OPS; ++op )
        {
          // update own derivative
          workspace[(2 * NUM_VERTS - (NUM_VERTS >> i) + j) * NUM_OPS + op] = (workspace[(j + pwr) * NUM_OPS + op] - workspace[j * NUM_OPS + op]) * axisStepInvs[i];
        }

        // update all dependent derivatives
        for( integer k = 0; k < i; k++ )
        {
          for( integer op = 0; op < NUM_OPS; ++op )
          {
            workspace[(2 * NUM_VERTS - (NUM_VERTS >> k) + j) * NUM_OPS + op] = workspace[(2 * NUM_VERTS - (NUM_VERTS >> k) + j) * NUM_OPS + op] + axisMults[i] *
                                                                               (workspace[(2 * NUM_VERTS - (NUM_VERTS >> k) + j + pwr) * NUM_OPS + op] -
                                                                                workspace[(2 * NUM_VERTS - (NUM_VERTS >> k) + j) * NUM_OPS + op]);
          }
        }

        for( integer op = 0; op < NUM_OPS; ++op )
        {
          // interpolate value
          workspace[j * NUM_OPS + op] = workspace[j * NUM_OPS + op] + (axisCoordinates[i] - axisLows[i]) * workspace[(2 * NUM_VERTS - (NUM_VERTS >> i) + j) * NUM_OPS + op];
        }
      }
      pwr /= 2;
    }
    for( integer op = 0; op < NUM_OPS; ++op )
    {
      values[op] = workspace[op];
      for( integer i = 0; i < NUM_DIMS; ++i )
      {
        derivatives[op * NUM_DIMS + i] = workspace[(2 * NUM_VERTS - (NUM_VERTS >> i)) * NUM_OPS + op];
      }
    }
  }

  // inputs : table discretization data

  /// Array [NUM_DIMS] of axis minimum values
  arrayView1d< real64 const > m_axisMinimums;

  /// Array [NUM_DIMS] of axis maximum values
  arrayView1d< real64 const > m_axisMaximums;

  /// Array [NUM_DIMS] of axis discretization points
  arrayView1d< integer const > m_axisPoints;

  // inputs : service data derived from table discretization data

  ///  Array [NUM_DIMS] of axis interval lengths (axes are discretized uniformly)
  arrayView1d< real64 const > m_axisSteps;

  ///  Array [NUM_DIMS] of inversions of axis interval lengths (axes are discretized uniformly)
  arrayView1d< real64 const > m_axisStepInvs;

  ///  Array [NUM_DIMS] of hypercube index mult factors for each axis
  arrayView1d< globalIndex const > m_axisHypercubeMults;

  // inputs: operator sample data

  ///  Main table data stored per hypercube: all values required for interpolation withing give hypercube are stored contiguously
  arrayView1d< real64 const > m_hypercubeData;

  // inputs: where to interpolate

  /// Coordinates in NUM_DIMS-dimensional space where interpolation is requested
  arrayView1d< real64 const > m_coordinates;

  // outputs

  /// Interpolated values
  arrayView1d< real64 > m_values;

  /// /// Interpolated derivatives
  arrayView1d< real64 > m_derivatives;
};

} /* namespace geosx */

#endif /* GEOSX_FUNCTIONS_MULTIVARIABLETABLEFUNCTIONKERNELS_HPP_ */
