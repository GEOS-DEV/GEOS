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
   * @param ei Index of element
   */
  GEOSX_HOST_DEVICE
  void
  compute( localIndex const ei ) const
  {
    interpolatePoint( &m_coordinates[ei * NUM_DIMS], &m_values[ei * NUM_OPS], &m_derivatives[ei * NUM_OPS] );
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

  GEOSX_HOST_DEVICE
  real64 const *
  getHypercubeData( globalIndex const hypercubeIndex ) const
  {
    static constexpr integer NUM_VERTS = 1 << NUM_DIMS;
    return &m_hypercubeData[hypercubeIndex * NUM_VERTS * NUM_OPS];
  }


  GEOSX_HOST_DEVICE
  integer
  getAxisIntervalIndexLowMult( real64 const axisCoordinate,
                               real64 const axisMin,
                               real64 const axisMax,
                               real64 const axisStep,
                               real64 const axisStepInv,
                               integer const axisPoints,
                               // OUTPUT:
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

  GEOSX_HOST_DEVICE
  void
  interpolatePointWithDerivatives( real64 const *axisCoordinates,
                                   real64 const *hypercubeData,
                                   real64 const *axisLows,
                                   real64 const *axisMults,
                                   real64 const *axisStepInvs,
                                   // OUTPUT:
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
      //printf ("i = %d, NUM_VERTS = %d, New offset: %d\n", i, NUM_VERTS, 2 * NUM_VERTS - (NUM_VERTS>>i));

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
