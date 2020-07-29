/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef SRC_CORECOMPONENTS_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS2_HPP_
#define SRC_CORECOMPONENTS_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS2_HPP_

#include "common/DataTypes.hpp"


/**
 * @class TrilinearHexahedronShapeFunctionKernel
 *
 * Contains the kernel accessible functions specific to the standard Trilinear
 * Hexahedron finite element with a Gaussian quadrature rule. It is assumed
 * that the indexing for the quadrature points mirrors that of the nodes.
 * Also note that the assumed node ordering is not the standard right-hand-rule
 * used in the literature. Here we use a Cartesian aligned numbering in order
 * to simplify the mapping to the parent coordinates and tensor product
 * indices.
 *
 */
class LagrangeBasis2
{





  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value0( constexpr real64 xi )
  {
    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    constexpr real64 xi_div2 = 0.5 * xi;
    return  -xi_div2 + xi_div2 * xi;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value1( constexpr real64 xi )
  {
    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return 1.0 - xi * xi ;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value2( constexpr real64 xi )
  {
    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    constexpr real64 xi_div2 = 0.5 * xi;
    return  xi_div2 + xi_div2 * xi;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient0( constexpr real64 xi )
  {
    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return -0.5 + xi;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient1( constexpr real64 xi )
  {
    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return - 2 * xi;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient2( constexpr real64 xi )
  {
    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return 0.5 + xi;
  }

  /**
   * @class TensorProduct3D
   *
   *                        24              25               26
   *                          o--------------o--------------o
   *                         /.                            /|
   *                        / .                           / |
   *                    21 o  .           o 22        23 o  |
   *                      /   .                         /   |
   *                     /    .         19             /    |
   *                 18 o--------------o--------------o 20  |
   *                    |     o              o        |     o
   *                    |     .15             16      |     |17
   *                    |     .                       |     |
   *                    |  o  .           o           |  o  |
   *                    |   12.            13         |   14|
   *                    |     .                       |     |
   *                  9 o     .        o 10           o 11  |
   *                    |     o..............o........|.....o
   *                    |    , 6              7       |    / 8
   *                    |   ,                         |   /
   *                    |  o              o           |  o         xi2
   *                    | , 3              4          | / 5        |
   *                    |,                            |/           | / xi1
   *                    o--------------o--------------o            |/
   *                   0                1              2           o----- xi0
   */
  struct TensorProduct3D
  {

    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @param k The index in the xi2 direction (0,1)
     * @return The linear index of the support/quadrature point (0-7)
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static int linearIndex( constexpr int i,
                                      constexpr int j,
                                      constexpr int k )
    {
      return i + 3 * j + 9 * k;
    }



    /**
     * @brief Calculate the Cartesian/TensorProduct index given the linear index
     *   of a support point.
     * @param a The linear index of support point
     * @param i0 The Cartesian index of the support point in the xi0 direction.
     * @param i1 The Cartesian index of the support point in the xi1 direction.
     * @param i2 The Cartesian index of the support point in the xi2 direction.
     * @return
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static void basisIndices( constexpr localIndex linearIndex,
                                        localIndex & i0,
                                        localIndex & i1,
                                        localIndex & i2 )
   {
     i2 = ( linearIndex * 29 ) >> 8 ;
     //i2 = a/9;

     i1 = ( ( linearIndex * 22 ) >> 6 ) - i2 * 3;
     //i1 = a/3 - i2 * 3;

     i0 = linearIndex - i1 * 3 - i2 * 9;
   }
  };


};


#endif /* SRC_CORECOMPONENTS_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_ */
