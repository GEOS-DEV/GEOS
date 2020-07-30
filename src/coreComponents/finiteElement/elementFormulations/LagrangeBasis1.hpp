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

#ifndef SRC_CORECOMPONENTS_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_
#define SRC_CORECOMPONENTS_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{
namespace finiteElement
{

class LagrangeBasis1
{
public:
  constexpr static localIndex numSupportPoints = 2;

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *   linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the xi0 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentSupportCoord( localIndex const supportPointIndex )
  {
//    GEOSX_ASSERT( supportPointIndex==0 || supportPointIndex==1 );
    return -1.0 + 2.0 * (supportPointIndex & 1);
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value( int const index,
                                 real64 const xi )
  {
//    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return 0.5 + 0.5 * xi * parentSupportCoord( index );
  }


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value0( real64 const xi )
  {
//    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return 0.5 - 0.5 * xi;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value1( real64 const xi )
  {
//    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return 0.5 + 0.5 * xi;
  }


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient( int const index,
                                    real64 const xi )
  {
    GEOSX_UNUSED_VAR( xi );
//    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return 0.5 * parentSupportCoord( index );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient0( real64 const xi )
  {
    GEOSX_UNUSED_VAR( xi );
//    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return -0.5;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient1( real64 const xi )
  {
    GEOSX_UNUSED_VAR( xi );
//    GEOSX_ASSERT( xi>=-1 && xi<=1 );
    return 0.5;
  }


  /**
   * @struct TensorProduct3D
   *                              6                   7
   *                               o-----------------o
   *                              /.                /|
   *                             / .               / |
   *                          4 o-----------------o 5|
   *                            |  .              |  |
   *                            |  .              |  |
   *                            |  .              |  |
   *                            |  .              |  |
   *                            |2 o..............|..o 3       xi2
   *                            | ,               | /          |
   *                            |,                |/           | / xi1
   *                            o-----------------o            |/
   *                           0                   1           ------ xi0
   *
   */
  struct TensorProduct3D
  {

    constexpr static localIndex numSupportPoints = 8;

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
    constexpr static int linearIndex( int const i,
                                      int const j,
                                      int const k )
    {
      return i + 2 * j + 4 * k;
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
    constexpr static void multiIndex( int const linearIndex,
                                      int & i0,
                                      int & i1,
                                      int & i2 )
    {
      i0 = ( linearIndex & 1 );
      i1 = ( linearIndex & 2 ) >> 1;
      i2 = ( linearIndex & 4 ) >> 2;
    }


    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    static void value( real64 const (&coords)[3],
                       real64 (& N)[numSupportPoints] )
    {
      for( int a=0; a<2; ++a )
      {
        for( int b=0; b<2; ++b )
        {
          for( int c=0; c<2; ++c )
          {
            int const lindex = LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );
            N[ lindex ] = LagrangeBasis1::value( a, coords[0] ) *
                          LagrangeBasis1::value( b, coords[1] ) *
                          LagrangeBasis1::value( c, coords[2] );
          }
        }
      }
    }

    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static real64 parentCoords0( localIndex const a )
    {
      return -1.0 + 2.0 * (a & 1);
    }

    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static real64 parentCoords1( localIndex const a )
    {
      return -1.0 + ( a & 2 );
    }

    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static real64 parentCoords2( localIndex const a )
    {
      return -1.0 + 0.5 * ( a & 4 );
    }


  };

};

}
}


#endif /* SRC_CORECOMPONENTS_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_ */
