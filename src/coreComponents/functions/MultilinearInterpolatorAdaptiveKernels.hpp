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
 * @file MultilinearInterpolatorAdaptiveKernels.hpp
 */

#ifndef GEOS_FUNCTIONS_MULTILINEARINTERPOLATORADAPTIVEKERNELS_HPP_
#define GEOS_FUNCTIONS_MULTILINEARINTERPOLATORADAPTIVEKERNELS_HPP_

#include <unordered_map>
#include "MultilinearInterpolatorBaseKernels.hpp"

namespace geos
{
/**
 * @class KernelWrapper
 *
 *
 */

/**
 * @class MultilinearInterpolatorAdaptiveKernel
 *
 * A class for multilinear piecewise interpolation with access to exact function evaluator
 * All functions are interpolated using the same uniformly discretized space
 *
 * @tparam NUM_DIMS number of dimensions (inputs)
 * @tparam NUM_OPS number of interpolated functions (outputs)
 */
template< integer NUM_DIMS, integer NUM_OPS >
class MultilinearInterpolatorAdaptiveKernel : public MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>
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
  using MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>::m_axisPointMults;

  /**
   * @brief Construct a new Multilinear Interpolator Adaptive Kernel object
   *
   * @param[in] axisMinimums  minimum coordinate for each axis
   * @param[in] axisMaximums maximum coordinate for each axis
   * @param[in] axisPoints number of discretization points between minimum and maximum for each axis
   * @param[in] axisSteps axis interval lengths (axes are discretized uniformly)
   * @param[in] axisStepInvs inversions of axis interval lengths (axes are discretized uniformly)
   * @param[in] axisHypercubeMults  hypercube index mult factors for each axis
   */
  MultilinearInterpolatorAdaptiveKernel( array1d< real64 > const & axisMinimums,
                                          array1d< real64 > const & axisMaximums,
                                          array1d< integer > const & axisPoints,
                                          array1d< real64 > const & axisSteps,
                                          array1d< real64 > const & axisStepInvs,
                                          array1d< globalIndex > const & axisHypercubeMults ):
    MultilinearInterpolatorBaseKernel<NUM_DIMS, NUM_OPS>( axisMinimums,
                                                          axisMaximums,
                                                          axisPoints,
                                                          axisSteps,
                                                          axisStepInvs,
                                                          axisHypercubeMults)
  {};
protected:
  /**
   * @brief Get values of operators at a given point 
   * Provide a reference to correct location in the adaptive point storage. 
   * If the point is not found, compute it first, and then return the reference.
   *
   * @param[in] pointIndex index of point 
   * @return operator values at given point
   */
  GEOS_HOST_DEVICE
  stackArray1d<real64, numOps> const &
  getPointData(globalIndex const pointIndex)
  {
    auto item = m_pointData.find(pointIndex);
    if (item == m_pointData.end())
    {
      stackArray1d<real64, numOps> newPoint;
      newPoint[0] = 1.;
      m_pointData[pointIndex] = newPoint;
      return m_pointData.at(pointIndex);
    }
    else
      return item->second;
  }
    /**
   * @brief Get indexes of all vertices for given hypercube
   * 
   * @param[in] hypercubeIndex index of the hyporcube
   * @param[out] hypercubePoints indexes of all vertices of hypercube
   */
  GEOS_HOST_DEVICE
  inline 
  void 
  getHypercubePoints( globalIndex const hypercubeIndex, 
                      stackArray1d<integer, numVerts>& hypercubePoints ) const
  {
    globalIndex remainderIdx = hypercubeIndex;
    integer pwr = numVerts;

    for (integer i = 0; i < numVerts; ++i)
      hypercubePoints[i] = 0;

    for (integer i = 0; i < numDims; ++i)
    {
      globalIndex axisIdx = remainderIdx / m_axisHypercubeMults[i];
      remainderIdx = remainderIdx % m_axisHypercubeMults[i];
      pwr /= 2;

      for (integer j = 0; j < numVerts; ++j)
      {
        integer zeroOrOne = (j / pwr) % 2;
        hypercubePoints[j] += (axisIdx + zeroOrOne) * m_axisPointMults[i];
      }
    }
  }
  /**
   * @brief Get pointer to hypercube data
   *
   * @param[in] hypercubeIndex
   * @return pointer to hypercube data
   */
  GEOS_HOST_DEVICE
  inline
  real64*
  getHypercubeData( globalIndex const hypercubeIndex ) override
  {
    auto item = m_hypercubeData.find(hypercubeIndex);
    if (item == m_hypercubeData.end())
    {
        stackArray1d<integer, numVerts> points (numVerts);
        stackArray1d<real64, numVerts * numOps> newHypercube;

        this->getHypercubePoints(hypercubeIndex, points);

        for (integer i = 0; i < numVerts; ++i)
        {
          // obtain point data and copy it to hypercube data
          auto const & pData = this->getPointData(points[i]);
          for (integer op = 0; op < numOps; op++)
          {
              newHypercube[i * numOps + op] = pData[op];
          }
        }
        m_hypercubeData[hypercubeIndex] = newHypercube;
        return m_hypercubeData[hypercubeIndex].data();
    }
    else
        return item->second.data();
  }

   /**
   * @brief adaptive point storage: the values of operators at requested supporting points
   * Storage is grown dynamically in the process of simulation. 
   * Only supporting points that are required for interpolation are computed and added
   * 
   */
  unordered_map<globalIndex, stackArray1d<real64, numOps>> m_pointData;

   /**
   * @brief adaptive hypercube storage: the values of operators at every vertex of reqested hypercubes
   * Storage is grown dynamically in the process of simulation
   * Only hypercubes that are required for interpolation are computed and added
   * 
   * In fact it is an excess storage used to reduce memory accesses during interpolation. 
   * Here all values of all vertexes of requested hypercube are stored consecutevely and are accessed via a single index
   * Usage of point_data for interpolation directly would require N_VERTS memory accesses (>1000 accesses for 10-dimensional space)
   *  * 
   */  
  unordered_map<globalIndex, stackArray1d<real64, numVerts * numOps>> m_hypercubeData;
};

} /* namespace geos */

#endif /* GEOS_FUNCTIONS_MULTILINEARINTERPOLATORADAPTIVEKERNELS_HPP_ */
