/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultilinearInterpolatorStaticKernels.hpp
 */

#ifndef GEOS_FUNCTIONS_MULTILINEARINTERPOLATORSTATICKERNELS_HPP_
#define GEOS_FUNCTIONS_MULTILINEARINTERPOLATORSTATICKERNELS_HPP_

#include "MultilinearInterpolatorBaseKernels.hpp"

namespace geos
{
/**
 * @class MultilinearInterpolatorStaticKernel
 *
 * A class for multilinear piecewise interpolation with static storage
 * All functions are interpolated using the same uniformly discretized space
 *
 * @tparam NUM_DIMS number of dimensions (inputs)
 * @tparam NUM_OPS number of interpolated functions (outputs)
 */
template< integer NUM_DIMS, integer NUM_OPS >
class MultilinearInterpolatorStaticKernel : public MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>
{
public:
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::numDims;
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::numOps;
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::numVerts;
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::m_axisMinimums;
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::m_axisMaximums;
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::m_axisPoints;
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::m_axisSteps;
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::m_axisStepInvs;
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::m_axisHypercubeMults;

  /**
   * @brief Construct a new Multilinear Interpolator Static Kernel object
   *
   * @param[in] axisMinimums  minimum coordinate for each axis
   * @param[in] axisMaximums maximum coordinate for each axis
   * @param[in] axisPoints number of discretization points between minimum and maximum for each axis
   * @param[in] axisSteps axis interval lengths (axes are discretized uniformly)
   * @param[in] axisStepInvs inversions of axis interval lengths (axes are discretized uniformly)
   * @param[in] axisHypercubeMults  hypercube index mult factors for each axis
   * @param[in] hypercubeData table data stored per hypercube
   */
  MultilinearInterpolatorStaticKernel( arrayView1d< real64 const > const & axisMinimums,
                                          arrayView1d< real64 const > const & axisMaximums,
                                          arrayView1d< integer const > const & axisPoints,
                                          arrayView1d< real64 const > const & axisSteps,
                                          arrayView1d< real64 const > const & axisStepInvs,
                                          arrayView1d< globalIndex const > const & axisHypercubeMults,
                                          arrayView1d< real64 const > const & hypercubeData ):
    MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>( axisMinimums,
                                                          axisMaximums,
                                                          axisPoints,
                                                          axisSteps,
                                                          axisStepInvs,
                                                          axisHypercubeMults),
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
    real64 axisLows[numDims];
    real64 axisMults[numDims];

    for( int i = 0; i < numDims; ++i )
    {
      integer const axisIndex = this->template getAxisIntervalIndexLowMult( coordinates[i],
                                                             m_axisMinimums[i], m_axisMaximums[i],
                                                             m_axisSteps[i], m_axisStepInvs[i], m_axisPoints[i],
                                                             axisLows[i], axisMults[i] );
      hypercubeIndex += axisIndex * m_axisHypercubeMults[i];
    }

    this->template interpolatePoint( coordinates,
                      getHypercubeData( hypercubeIndex ),
                      &axisLows[0],
                      &m_axisStepInvs[0],
                      values );
  }

  /**
   * @brief interpolate all operators and compute their derivatives at a given point
   *
   * @tparam IN_ARRAY type of input array of coordinates
   * @tparam OUT_ARRAY type of output array of values
   * @tparam OUT_2D_ARRAY type of output array of derivatives
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
    real64 axisLows[numDims];
    real64 axisMults[numDims];

    for( int i = 0; i < numDims; ++i )
    {
      integer const axisIndex = this->getAxisIntervalIndexLowMult(coordinates[i],
                                                                  m_axisMinimums[i], m_axisMaximums[i],
                                                                  m_axisSteps[i], m_axisStepInvs[i], m_axisPoints[i],
                                                                  axisLows[i], axisMults[i] );
      hypercubeIndex += axisIndex * m_axisHypercubeMults[i];
    }

    this->template interpolatePointWithDerivatives( coordinates,
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
  GEOS_HOST_DEVICE
  inline
  real64 const *
  getHypercubeData( globalIndex const hypercubeIndex ) const
  {
    return &m_hypercubeData[hypercubeIndex * numVerts * numOps];
  }

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

#endif /* GEOS_FUNCTIONS_MULTILINEARINTERPOLATORSTATICKERNELS_HPP_ */
