/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultivariableNonuniformTableFunctionKernels.hpp
 */

#ifndef GEOS_FUNCTIONS_MULTIVARIABLENONUNIFORMTABLEFUNCTIONKERNELS_HPP_
#define GEOS_FUNCTIONS_MULTIVARIABLENONUNIFORMTABLEFUNCTIONKERNELS_HPP_

namespace geos
{

/**
 * @class MultivariableNonuniformTableFunctionStaticKernel
 *
 * A class for multivariable piecewise interpolation with static storage
 * Functions are interpolated assuming a non uniform discretized space
 * for independent variables
 * @tparam NUM_DIMS number of dimensions (inputs)
 * @tparam NUM_OPS number of interpolated functions (outputs)
 */
template< integer NUM_DIMS, integer NUM_OPS >
class MultivariableNonuniformTableFunctionStaticKernel
{
public:

  /// Compile time value for the number of table dimensions (inputs)
  static constexpr integer numDims = NUM_DIMS;

  /// Compile time value for the number of operators (interpolated functions, outputs)
  static constexpr integer numOps = NUM_OPS;

  /// Compile time value for the number of hypercube vertices
  static constexpr integer numVerts = 1 << numDims;


  /**
   * @brief Construct a new Multivariable Table Function Static Kernel object
   *
   * @param[in] axisCoordinates  discretization points for each axis
   * @param[in] axisPoints number of discretization points for each axis
   * @param[in] axisSteps axis interval lengths (nonuniform axis discretization)
   * @param[in] axisStepInvs inversions of axis interval lengths (nonuniform axis discretization)
   * @param[in] axisHypercubeMults  hypercube index mult factors for each axis
   * @param[in] hypercubeData table data stored per hypercube
   */
  MultivariableNonuniformTableFunctionStaticKernel( arrayView2d< real64 const > const & axisCoordinates,
                                                    arrayView1d< integer const > const & axisPoints,
                                                    arrayView2d< real64 const > const & axisSteps,
                                                    arrayView2d< real64 const > const & axisStepInvs,
                                                    arrayView1d< globalIndex const > const & axisHypercubeMults,
                                                    arrayView1d< real64 const > const & hypercubeData ):
    m_axisCoordinates ( axisCoordinates ),
    m_axisPoints ( axisPoints ),
    m_axisSteps ( axisSteps ),
    m_axisStepInv ( axisStepInvs ),
    m_axisHypercubeMults ( axisHypercubeMults ),
    m_hypercubeData ( hypercubeData )
  {};

/**
 * @brief interpolate all operators at a given point
 *
 * @param[in] coordinates point coordinates
 * @param[out] values interpolated operator values
 */
  template< typename IN_ARRAY, typename OUT_ARRAY >
  GEOS_HOST_DEVICE
  void
  compute( IN_ARRAY const & coordinates,
           OUT_ARRAY && values ) const
  {
    globalIndex hypercubeIndex = 0;
    real64 axisLows[numDims] = {0.0};
    real64 axisStepInv[numDims] = {0.0};
    real64 axisMults[numDims] = {0.0};

    for( int i = 0; i < numDims; ++i )
    {
      integer const axisIndex = getAxisIntervalIndexLowMult( coordinates[i],
                                                             m_axisCoordinates[i],
                                                             m_axisStepInv[i],
                                                             m_axisPoints[i],
                                                             axisLows[i], axisStepInv[i], axisMults[i] );
      hypercubeIndex += axisIndex * m_axisHypercubeMults[i];
    }

    interpolatePoint( coordinates,
                      getHypercubeData( hypercubeIndex ),
                      &axisLows[0],
                      &axisStepInv[0],
                      values );
  }

  /**
   * @brief interpolate all operators and compute their derivatives at a given point
   *
   * @param[in] coordinates point coordinates
   * @param[out] values interpolated operator values
   * @param[out] derivatives derivatives of interpolated operators
   */
  template< typename IN_ARRAY, typename OUT_ARRAY, typename OUT_2D_ARRAY >
  GEOS_HOST_DEVICE
  void
  compute( IN_ARRAY const & coordinates,
           OUT_ARRAY && values,
           OUT_2D_ARRAY && derivatives ) const
  {
    globalIndex hypercubeIndex = 0;
    real64 axisLows[numDims] = {0.0};
    real64 axisStepInv[numDims] = {0.0};
    real64 axisMults[numDims] = {0.0};

    for( int i = 0; i < numDims; ++i )
    {
      integer const axisIndex = getAxisIntervalIndexLowMult( coordinates[i],
                                                             m_axisCoordinates[i],
                                                             m_axisStepInv[i],
                                                             m_axisPoints[i],
                                                             axisLows[i], axisStepInv[i], axisMults[i] );
      std::cout << " axis " << i << " index " << axisIndex << " mul " <<   m_axisHypercubeMults[i] <<std::endl;
      hypercubeIndex += axisIndex * m_axisHypercubeMults[i];
    }
    std::cout << " hypercubeIndex " << hypercubeIndex << std::endl;
    interpolatePointWithDerivatives( coordinates,
                                     getHypercubeData( hypercubeIndex ),
                                     &axisLows[0], &axisMults[0],
                                     &axisStepInv[0],
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
  GEOS_HOST_DEVICE
  inline
  real64 const *
  getHypercubeData( globalIndex const hypercubeIndex ) const
  {
    return &m_hypercubeData[hypercubeIndex * numVerts * numOps];
  }

  /**
   * @brief Get the interval index, low and mult values for a given axis coordinate
   *
   * @param[in] coordinate search coordinate on a given axis
   * @param[in] axisCoordinates coordinates for given axis
   * @param[in] axisStepInv inversion of the interval length for given axis
   * @param[in] axisPoints number of discretization points for  given axis
   * @param[out] axisLow for interval containing coordinate, the left coordinate of the interval
   * @param[out] intervalStepInv for interval containing coordinate, the inversion of the interval length
   * @param[out] axisMult weight of the right coordinate of target axis interval
   * @return integer target axis interval index
   */
  GEOS_HOST_DEVICE
  inline
  integer
  getAxisIntervalIndexLowMult( real64 const coordinate,
                               real64 const * const axisCoordinates,
                               real64 const * const axisStepInv,
                               integer const axisPoints,
                               real64 & axisLow,
                               real64 & intervalStepInv,
                               real64 & axisMult ) const
  {
    integer axisIntervalIndex=0;
    std::cout << " coordinate " << coordinate << " axis points " <<  axisCoordinates[0] << " " <<  axisCoordinates[axisPoints-1] << std::endl;
    if( coordinate < axisCoordinates[0] )
    {
      axisIntervalIndex=0;
      intervalStepInv = axisStepInv[axisIntervalIndex];
      axisLow = axisCoordinates[axisIntervalIndex];
      axisMult = (coordinate - axisCoordinates[axisIntervalIndex]) * axisStepInv[axisIntervalIndex];
      printf( "Interpolation warning: axis coordinate is less than lower limit  (%lf ) with value %lf, extrapolation is applied\n", axisCoordinates[0], coordinate );
    }
    else if( coordinate >  axisCoordinates[axisPoints-1] )
    {
      axisIntervalIndex = std::max( 0, axisPoints-2 );
      //axisIntervalIndex= axisPoints - 2;
      intervalStepInv = axisStepInv[axisIntervalIndex];
      axisLow = axisCoordinates[axisIntervalIndex];
      axisMult = (coordinate - axisCoordinates[axisIntervalIndex]) * axisStepInv[axisIntervalIndex];
      printf( "Interpolation warning: axis coordinate is beyond upper limit ( %lf) with value %lf, extrapolation is applied\n", axisCoordinates[axisPoints-1], coordinate );
    }
    else
    {
      for( int j=1; j<axisPoints; j++ )
      {
        if( coordinate <= axisCoordinates[j] )
        {
          axisIntervalIndex=j-1;
          intervalStepInv = axisStepInv[j-1];
          axisLow = axisCoordinates[j-1];
          axisMult = (coordinate - axisCoordinates[j-1]) * axisStepInv[j-1];
          break;
        }
      }
    }

    return axisIntervalIndex;
  }

  /**
   * @brief interpolate all operators values at a given point
   * The algoritm is based on http://dx.doi.org/10.1090/S0025-5718-1988-0917826-0
   *
   * @param[in] axisCoordinates coordinates of a point
   * @param[in] hypercubeData data of target hypercube
   * @param[in] axisLows array of left coordinates of target axis intervals
   * @param[in] axisStepInvs array of inversions of axis steps
   * @param[out] values interpolated operator values
   */
  template< typename IN_ARRAY, typename OUT_ARRAY >
  GEOS_HOST_DEVICE
  inline
  void
  interpolatePoint( IN_ARRAY const & axisCoordinates,
                    real64 const * const hypercubeData,
                    real64 const * const axisLows,
                    real64 const * const axisStepInvs,
                    OUT_ARRAY && values ) const
  {
    integer pwr = numVerts / 2;   // distance between high and low values
    real64 workspace[numVerts][numOps];

    // copy operator values for all vertices
    for( integer i = 0; i < numVerts; ++i )
    {
      for( integer j = 0; j < numOps; ++j )
      {
        workspace[i][j] = hypercubeData[i * numOps + j];
      }
    }

    for( integer i = 0; i < numDims; ++i )
    {

      for( integer j = 0; j < pwr; ++j )
      {
        for( integer op = 0; op < numOps; ++op )
        {
          // update own derivative
          workspace[j][op] += (axisCoordinates[i] - axisLows[i]) * (workspace[j + pwr][op] - workspace[j][op]) * axisStepInvs[i];
        }
      }
      pwr /= 2;
    }
    for( integer op = 0; op < numOps; ++op )
    {
      values[op] = workspace[0][op];
    }
  }


  /**
   * @brief interpolate all operators values and derivatives at a given point
   * The algoritm is based on http://dx.doi.org/10.1090/S0025-5718-1988-0917826-0
   *
   * @param[in] axisCoordinates coordinates of a point
   * @param[in] hypercubeData data of target hypercube
   * @param[in] axisLows array of left coordinates of target axis intervals
   * @param[in] axisMults array of weights of right coordinates of target axis intervals
   * @param[in] axisStepInvs array of inversions of axis steps
   * @param[out] values interpolated operator values
   * @param[out] derivatives derivatives of interpolated operators
   */
  template< typename IN_ARRAY, typename OUT_ARRAY, typename OUT_2D_ARRAY >
  GEOS_HOST_DEVICE
  inline
  void
  interpolatePointWithDerivatives( IN_ARRAY const & axisCoordinates,
                                   real64 const * const hypercubeData,
                                   real64 const * const axisLows,
                                   real64 const * const axisMults,
                                   real64 const * const axisStepInvs,
                                   OUT_ARRAY && values,
                                   OUT_2D_ARRAY && derivatives ) const
  {
    integer pwr = numVerts / 2;   // distance between high and low values
    real64 workspace[2 * numVerts - 1][numOps];

    // copy operator values for all vertices
    for( integer i = 0; i < numVerts; ++i )
    {
      for( integer j = 0; j < numOps; ++j )
      {
        workspace[i][j] = hypercubeData[i * numOps + j];
      }
    }

    for( integer i = 0; i < numDims; ++i )
    {

      for( integer j = 0; j < pwr; ++j )
      {
        for( integer op = 0; op < numOps; ++op )
        {
          // update own derivative
          workspace[2 * numVerts - (numVerts >> i) + j][op] = (workspace[j + pwr][op] - workspace[j][op]) * axisStepInvs[i];
        }

        // update all dependent derivatives
        for( integer k = 0; k < i; k++ )
        {
          for( integer op = 0; op < numOps; ++op )
          {
            workspace[2 * numVerts - (numVerts >> k) + j][op] = workspace[2 * numVerts - (numVerts >> k) + j][op] + axisMults[i] *
                                                                (workspace[2 * numVerts - (numVerts >> k) + j + pwr][op] -
                                                                 workspace[2 * numVerts - (numVerts >> k) + j][op]);
          }
        }

        for( integer op = 0; op < numOps; ++op )
        {
          // interpolate value
          workspace[j][op] = workspace[j][op] + (axisCoordinates[i] - axisLows[i]) * workspace[2 * numVerts - (numVerts >> i) + j][op];
        }
      }
      pwr /= 2;
    }
    for( integer op = 0; op < numOps; ++op )
    {
      values[op] = workspace[0][op];
      for( integer i = 0; i < numDims; ++i )
      {
        derivatives[op][i] = workspace[2 * numVerts - (numVerts >> i)][op];
      }
    }
  }

  // inputs : table discretization data

  /// Array [numDims, numCoords[numDims]] of axis coordinate values for each dim
  arrayView2d< real64 const > m_axisCoordinates;

  /// Array [numDims] of axis discretization points
  arrayView1d< integer const > m_axisPoints;

  // inputs : service data derived from table discretization data

  ///  Array [numDims] of axis interval lengths v
  arrayView2d< real64 const > m_axisSteps;

  ///  Array [numDims] of inversions of axis interval lengths v
  arrayView2d< real64 const > m_axisStepInv;

  ///  Array [numDims] of hypercube index mult factors for each axis
  arrayView1d< globalIndex const > m_axisHypercubeMults;

  // inputs: operator sample data

  ///  Main table data stored per hypercube: all values required for interpolation withing give hypercube are stored contiguously
  arrayView1d< real64 const > m_hypercubeData;

  // inputs: where to interpolate

  /// Coordinates in numDims-dimensional space where interpolation is requested
  arrayView1d< real64 const > m_coordinates;

  // outputs

  /// Interpolated values
  arrayView1d< real64 > m_values;

  /// /// Interpolated derivatives
  arrayView1d< real64 > m_derivatives;
};


} /* namespace geos */

#endif /* GEOS_FUNCTIONS_MULTIVARIABLETABLEFUNCTIONKERNELS_HPP_ */
