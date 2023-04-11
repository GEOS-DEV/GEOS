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
 * @file H1_Tetrahedron_Lagrange2_Gauss4.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE2GAUSS4
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE2GAUSS4

#include "FiniteElementBase.hpp"
#include "LagrangeBasis2.hpp"


namespace geosx
{
namespace finiteElement
{

/**
 * This class contains the kernel accessible functions specific to the
 * H1-conforming nodal quadratic tetrahedron finite element with a 4-point
 * Gaussian quadrature rule.
 *
 *
 *          3                                =====  ==  ==  ==
 *           +                               Node   r   s   t
 *           |\\_                            =====  ==  ==  ==
 *           ||  \_                          0      0   0   0
 *           | \   \_           t            1      1   0   1
 *           |  +__  \_         |   s        2      0   1   0
 *           | /2  \__ \_       |  /         3      0   0   1
 *           |/       \__\      | /          4      0.5 0   0
 *           +------------+     *------r     5      0.5 0.5 0
 *          0              1                 6      0   0.5 0
 *                                           7      0   0   0.5
 *                                           8      0.5 0   0.5
 *                                           9      0   0.5 0.5
 *
 */
class H1_Tetrahedron_Lagrange2_Gauss4 final : public FiniteElementBase
{
public:

  /// The type of basis used for this element
  using BASIS = LagrangeBasis2;

  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 10;
  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 4;

  /// The number of sampling points per element.
  constexpr static int numSamplingPoints = numSamplingPointsPerDirection * numSamplingPointsPerDirection * numSamplingPointsPerDirection;

  virtual ~H1_Tetrahedron_Lagrange2_Gauss4() override
  {}

  GEOSX_HOST_DEVICE
  virtual localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  /**
   * @brief Get the number of quadrature points.
   * @param stack Stack variables as filled by @ref setupStack.
   * @return The number of quadrature points.
   */
  GEOSX_HOST_DEVICE
  static localIndex
  getNumQuadraturePoints( StackVariables const & stack )
  {
    GEOSX_UNUSED_VAR( stack );
    return numQuadraturePoints;
  }

  GEOSX_HOST_DEVICE
  virtual localIndex getNumSupportPoints() const override
  {
    return numNodes;
  }

  GEOSX_HOST_DEVICE
  virtual localIndex getMaxSupportPoints() const override
  {
    return maxSupportPoints;
  }

  /**
   * @brief Get the number of support points.
   * @param stack Object that holds stack variables.
   * @return The number of support points.
   */
  GEOSX_HOST_DEVICE
  static localIndex getNumSupportPoints( StackVariables const & stack )
  {
    GEOSX_UNUSED_VAR( stack );
    return numNodes;
  }

  /**
   * @brief Get the Sampling Point Coord In the Parent Space
   *
   * @param linearIndex linear index of the sampling point
   * @param samplingPointCoord coordinates of the sampling point
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void getSamplingPointCoordInParentSpace( int const & linearIndex,
                                                  real64 (& samplingPointCoord)[3] )
  {
    int const i0 = linearIndex % numSamplingPointsPerDirection;
    int const i1 = ( (linearIndex - i0)/numSamplingPointsPerDirection ) % numSamplingPointsPerDirection;
    int const i2 = ( (linearIndex - i0)/numSamplingPointsPerDirection - i1 ) / numSamplingPointsPerDirection;

    real64 const step = 1 / ( numSamplingPointsPerDirection - 1 );

    real64 const r = i0 * step;
    real64 const s = i1 * step;
    real64 const t = i2 * step;
    if( (r+s) <= t )
    {
      samplingPointCoord[0] = r;
      samplingPointCoord[1] = s;
      samplingPointCoord[2] = t;
    }
    else
    {
      // if outside of the triangle need to reproject it. Points will be doubled though.
      samplingPointCoord[0] = 1 - t - r;
      samplingPointCoord[1] = 1 - t - s;
      samplingPointCoord[2] = t;
    }
  }

  /**
   * @brief Calculate shape functions values for each support point at a
   *   given point in the parent space.
   * @param pointCoord coordinates of the given point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void calcN( real64 const (&pointCoord)[3],
                     real64 ( &N )[numNodes] );

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOSX_HOST_DEVICE
  static void calcN( localIndex const q,
                     real64 ( &N )[numNodes] );

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void calcN( localIndex const q,
                     StackVariables const & stack,
                     real64 ( &N )[numNodes] );

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[numNodes][3],
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param stack Variables allocated on the stack as filled by @ref setupStack.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[numNodes][3],
                           StackVariables const & stack,
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the integration weights for a quadrature point.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @return The product of the quadrature rule weight and the determinate of
   *   the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 transformedQuadratureWeight( localIndex const q,
                                             real64 const (&X)[numNodes][3] );

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param q The quadrature point index in 3d space.
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   * @return The determinant of the Jacobian transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 invJacobianTransformation( int const q,
                                           real64 const (&X)[numNodes][3],
                                           real64 ( & J )[3][3] )
  {
    GEOSX_UNUSED_VAR( q, X );
//    jacobianTransformation( q, X, J );
    return LvArray::tensorOps::invert< 3 >( J );
  }
private:
  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = 1.0 / 6.0;

  /// The weight of each quadrature point.
  constexpr static real64 weight = parentVolume*0.25;


  /**
   * @brief Calculates the determinant of the Jacobian of the isoparametric
   *        mapping from the parent space to the physical space.
   * @param X Array containing the coordinates of the support points
   * @return determinant value
   */
  GEOSX_HOST_DEVICE
  static real64 determinantJacobianTransformation( real64 const (&X)[numNodes][3] );


};

/// @cond Doxygen_Suppress

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange2_Gauss4::
  determinantJacobianTransformation( real64 const (&X)[numNodes][3] )
{

  // assume elements are straight sided and 1st 4 nodes of quadratic are the nodes for linear
  return ( X[1][0] - X[0][0] )*( ( X[2][1] - X[0][1] )*( X[3][2] - X[0][2] ) - ( X[3][1] - X[0][1] )*( X[2][2] - X[0][2] ) )
         + ( X[1][1] - X[0][1] )*( ( X[3][0] - X[0][0] )*( X[2][2] - X[0][2] ) - ( X[2][0] - X[0][0] )*( X[3][2] - X[0][2] ) )
         + ( X[1][2] - X[0][2] )*( ( X[2][0] - X[0][0] )*( X[3][1] - X[0][1] ) - ( X[3][0] - X[0][0] )*( X[2][1] - X[0][1] ) );
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Tetrahedron_Lagrange2_Gauss4::
  calcN( real64 const (&coords)[3],
         real64 (& N)[numNodes] )
{
  real64 const r = coords[0];
  real64 const s = coords[1];
  real64 const t = coords[2];
  real64 const u = 1 - r - s - t;

  N[0] = u*(2.0*u - 1.0);
  N[1] = r*(2.0*r - 1.0);
  N[2] = s*(2.0*s - 1.0);
  N[3] = t*(2.0*t - 1.0);
  N[4] = 4.0*u*r;
  N[5] = 4.0*r*s;
  N[6] = 4.0*s*u;
  N[7] = 4.0*u*t;
  N[8] = 4.0*r*t;
  N[9] = 4.0*s*t;

}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Tetrahedron_Lagrange2_Gauss4::
  calcN( localIndex const q,
         real64 (& N)[numNodes] )
{

  real64 const quadPointCoords[4][3] = {{0.58541020, 0.13819660, 0.13819660},
    {0.13819660, 0.58541020, 0.13819660},
    {0.13819660, 0.13819660, 0.58541020},
    {0.13819660, 0.13819660, 0.13819660}};

  real64 const pointCoord[3] = {quadPointCoords[q][0], quadPointCoords[q][1], quadPointCoords[q][2]};
  calcN( pointCoord, N );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_Tetrahedron_Lagrange2_Gauss4::
  calcN( localIndex const q,
         StackVariables const & GEOSX_UNUSED_PARAM( stack ),
         real64 ( & N )[numNodes] )
{
  return calcN( q, N );
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange2_Gauss4::
  calcGradN( localIndex const q,
             real64 const (&X)[numNodes][3],
             real64 (& gradN)[numNodes][3] )
{

  real64 const quadPointCoords[4][3] = {{0.58541020, 0.13819660, 0.13819660},
    {0.13819660, 0.58541020, 0.13819660},
    {0.13819660, 0.13819660, 0.58541020},
    {0.13819660, 0.13819660, 0.13819660}};

  real64 r = quadPointCoords[q][0];
  real64 s = quadPointCoords[q][1];
  real64 t = quadPointCoords[q][2];

  real64 parametricGradient[numNodes][3];

  parametricGradient[0][0] = 4.0*r + 4.0*s + 4.0*t - 3.0;
  parametricGradient[0][1] = 4.0*r + 4.0*s + 4.0*t - 3.0;
  parametricGradient[0][2] = 4.0*r + 4.0*s + 4.0*t - 3.0;

  parametricGradient[1][0] = 4.0*r - 1.0;
  parametricGradient[1][1] = 0.0;
  parametricGradient[1][2] = 0.0;

  parametricGradient[2][0] = 0.0;
  parametricGradient[2][1] = 4.0*s - 1.0;
  parametricGradient[2][2] = 0.0;

  parametricGradient[3][0] = 0.0;
  parametricGradient[3][1] = 0.0;
  parametricGradient[3][2] = 4.0*t - 1.0;

  parametricGradient[4][0] = 4.0 - 8.0*r - 4.0*s - 4.0*t;
  parametricGradient[4][1] = -4.0*r;
  parametricGradient[4][2] = -4.0*r;

  parametricGradient[5][0] = 4.0*s;
  parametricGradient[5][1] = 4.0*r;
  parametricGradient[5][2] = 0.0;

  parametricGradient[6][0] = -4.0*s;
  parametricGradient[6][1] = 4.0 - 4.0*r - 8.0*s - 4.0*t;
  parametricGradient[6][2] = -4.0*s;

  parametricGradient[7][0] = -4.0*t;
  parametricGradient[7][1] = -4.0*t;
  parametricGradient[7][2] = 4.0 - 4.0*r - 4.0*s - 8.0*t;

  parametricGradient[8][0] = 4.0*t;
  parametricGradient[8][1] = 0.0;
  parametricGradient[8][2] = 4.0*r;

  parametricGradient[9][0] = 0.0;
  parametricGradient[9][1] = 4.0*t;
  parametricGradient[9][2] = 4.0*s;


  // Need Jacobian for linear element
  real64 jacobian[3][3];

  jacobian[0][0] = X[1][0] - X[0][0];
  jacobian[0][1] = X[2][0] - X[0][0];
  jacobian[0][2] = X[3][0] - X[0][0];

  jacobian[1][0] = X[1][1] - X[0][1];
  jacobian[1][1] = X[2][1] - X[0][1];
  jacobian[1][2] = X[3][1] - X[0][1];

  jacobian[2][0] =  X[1][2] - X[0][2];
  jacobian[2][1] =  X[2][2] - X[0][2];
  jacobian[2][2] =  X[3][2] - X[0][2];

  real64 const detJ_invert = LvArray::tensorOps::invert< 3 >( jacobian );

  GEOSX_UNUSED_VAR( detJ_invert );

  real64 detJ = determinantJacobianTransformation( X );

  //real64 factor = 1.0 / ( detJ );

  for( int i = 0; i < numNodes; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      gradN[i][j] = parametricGradient[i][0]*jacobian[0][j] + parametricGradient[i][1]*jacobian[1][j] + parametricGradient[i][2]*jacobian[2][j];
    }
  }

  return detJ * weight;
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 H1_Tetrahedron_Lagrange2_Gauss4::
  calcGradN( localIndex const q,
             real64 const (&X)[numNodes][3],
             StackVariables const & GEOSX_UNUSED_PARAM( stack ),
             real64 ( & gradN )[numNodes][3] )
{
  return calcGradN( q, X, gradN );
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange2_Gauss4::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  GEOSX_UNUSED_VAR( q );

  real64 detJ =  determinantJacobianTransformation( X );

  return detJ * weight;
}

/// @endcond

}
}

#endif //GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE1GAUSS1
