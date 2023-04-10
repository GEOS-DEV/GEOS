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
 * @file H1_Hexahedron_Lagrange2_GaussLegendre3.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_QUADRATICHEXAHEDRON
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_QUADRATICHEXAHEDRON

#include "FiniteElementBase.hpp"
#include "LagrangeBasis2.hpp"

#include <utility>



namespace geosx
{
namespace finiteElement
{


/**
   * This class contains the kernel acessible functions specific to quadratic hexahedron
   * finite element with Gauss quadrature rule.  
   *
   *                                                                  ____________________
   *                                                                 |Node   xi0  xi1  xi2|
   *                                                                 |=====  ===  ===  ===|
   *                                                                 |  0    -1   -1   -1 |
   *                                                                 |  1     0   -1   -1 |
   *                                                                 |  2     1   -1   -1 |
   *              24              25               26                |  3    -1    0   -1 |
   *                o--------------o--------------o                  |  4     0    0   -1 |
   *               /.                            /|                  |  5     1    0   -1 |
   *              / .                           / |                  |  6    -1    1   -1 |
   *          21 o  .           o 22        23 o  |                  |  7     0    1   -1 |
   *            /   .                         /   |                  |  8     1    1   -1 |
   *           /    .         19             /    |                  |  9    -1   -1    0 |
   *       18 o--------------o--------------o 20  |                  | 10     0   -1    0 |
   *          |     o              o        |     o                  | 11     1   -1    0 |
   *          |     .15             16      |     |17                | 12    -1    0    0 |
   *          |     .                       |     |                  | 13     0    0    0 |
   *          |  o  .           o           |  o  |                  | 14     1    0    0 |
   *          |   12.            13         |   14|                  | 15    -1    1    0 |
   *          |     .                       |     |                  | 16     0    1    0 |
   *        9 o     .        o 10           o 11  |                  | 17     1    1    0 |
   *          |     o..............o........|.....o                  | 18    -1   -1    1 |
   *          |    , 6              7       |    / 8                 | 19     0   -1    1 |
   *          |   ,                         |   /                    | 20     1   -1    1 |
   *          |  o              o           |  o         xi2         | 21    -1    0    1 |
   *          | , 3              4          | / 5        |           | 22     0    0    1 |
   *          |,                            |/           | / xi1     | 23     1    0    1 |
   *          o--------------o--------------o            |/          | 24    -1    1    1 |
   *         0                1              2           o----- xi0  | 25     0    1    1 |
   *                                                                 | 26     1    1    1 |
   *                                                                 |____________________|
   *
   */

class H1_Hexahedron_Lagrange2_GaussLegendre3 final : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = LagrangeBasis2::TensorProduct3D::numSupportPoints;
  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 27;

  /// The number of sampling points per element
  constexpr static int numSamplingPoints = numSamplingPointsPerDirection * numSamplingPointsPerDirection * numSamplingPointsPerDirection;

  /** @cond Doxygen_Suppress */
  USING_FINITEELEMENTBASE
  /** @endcond Doxygen_Suppress */

  virtual ~H1_Hexahedron_Lagrange2_GaussLegendre3() override
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
  static localIndex getNumQuadraturePoints( StackVariables const & stack )
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
  static void getSamplingPointCoordInParentSpace( int const linearIndex,
                                                  real64 (& samplingPointCoord)[3] )
  {
    const int i0 = linearIndex % numSamplingPointsPerDirection;
    const int i1 = ( (linearIndex - i0)/numSamplingPointsPerDirection ) % numSamplingPointsPerDirection;
    const int i2 = ( (linearIndex - i0)/numSamplingPointsPerDirection - i1 ) / numSamplingPointsPerDirection;

    constexpr real64 step = 2 / ( numSamplingPointsPerDirection - 1 );

    samplingPointCoord[0] = -1 + i0 * step;
    samplingPointCoord[1] = -1 + i1 * step;
    samplingPointCoord[2] = -1 + i2 * step;
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
                     real64 (& N)[numNodes] )
  {
    LagrangeBasis2::TensorProduct3D::value( pointCoord, N );
  }

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void calcN( localIndex const q,
                     real64 (& N)[numNodes] )
  {
    int qa, qb, qc;
    LagrangeBasis2::TensorProduct3D::multiIndex( q, qa, qb, qc );
    real64 const qCoords[3] = { quadratureFactor * LagrangeBasis2::parentSupportCoord( qa ),
                                quadratureFactor * LagrangeBasis2::parentSupportCoord( qb ),
                                quadratureFactor * LagrangeBasis2::parentSupportCoord( qc ) };

    calcN( qCoords, N );
  }

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
                     real64 ( & N )[numNodes] )
  {
    GEOSX_UNUSED_VAR( stack );
    return calcN( q, N );
  }

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
    int qa, qb, qc;
    LagrangeBasis2::TensorProduct3D::multiIndex( q, qa, qb, qc );
    jacobianTransformation( qa, qb, qc, X, J );
    return LvArray::tensorOps::invert< 3 >( J );
  }


  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a quadrature point using the stored inverse of the Jacobian
   *   transformation matrix.
   * @param q The linear index of the quadrature point.
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The symmetric gradient in Voigt notation.
   */
  GEOSX_HOST_DEVICE
  static void symmetricGradient( int const q,
                                 real64 const (&invJ)[3][3],
                                 real64 const (&var)[numNodes][3],
                                 real64 ( &grad )[6] );



  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the stored basis function gradients for all support points.
   * @param q The linear index of the quadrature point.
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The gradient.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad_{ij}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ai}\right ),
   * \f]
   *
   */
  GEOSX_HOST_DEVICE
  static void gradient( int const q,
                        real64 const (&invJ)[3][3],
                        real64 const (&var)[numNodes][3],
                        real64 ( &grad )[3][3] );


  /**
   * @brief Inner product of all basis function gradients and a rank-2
   *   symmetric tensor evaluated at a quadrature point.
   * @param q The linear index of the quadrature point.
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The rank-2 symmetric tensor at @p q.
   * @param R The vector resulting from the tensor contraction.
   *
   * More precisely, the operator is defined as:
   * \f[
   * R_i = \sum_a^{nSupport} \left( \frac{\partial N_a}{\partial X_j} var_{ij} \right),
   * \f]
   * where \f$\frac{\partial N_a}{\partial X_j}\f$ is the basis function gradient,
   *   \f$var_{ij}\f$ is the rank-2 symmetric tensor.
   */
  GEOSX_HOST_DEVICE
  static void plusGradNajAij( int const q,
                              real64 const (&invJ)[3][3],
                              real64 const (&var)[6],
                              real64 ( &R )[numNodes][3] );



  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   */
  GEOSX_HOST_DEVICE
  static void jacobianTransformation( int const qa,
                                      int const qb,
                                      int const qc,
                                      real64 const (&X)[numNodes][3],
                                      real64 ( &J )[3][3] );


  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  GEOSX_HOST_DEVICE
  static void
    applyTransformationToParentGradients( int const qa,
                                          int const qb,
                                          int const qc,
                                          real64 const ( &invJ )[3][3],
                                          real64 ( &gradN )[numNodes][3] );


private:
  /// The length of one dimension of the parent element.
  constexpr static real64 parentLength = LagrangeBasis2::parentSupportCoord( 2 ) - LagrangeBasis2::parentSupportCoord( 0 );

  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = parentLength*parentLength*parentLength;

  /// The weight of each quadrature point.
  constexpr static real64 weight [3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};

  /// The scaling factor specifying the location of the quadrature points
  /// relative to the origin
  constexpr static real64 quadratureFactor = 0.774596669241483;


  /**
   * @brief Applies a function inside a generic loop in over the tensor product
   *   indices.
   * @tparam FUNC The type of function to call within the support loop.
   * @tparam PARAMS The parameter pack types to pass through to @p FUNC.
   * @param qa The 1d quadrature point index in xi0 direction (0,1,2)
   * @param qb The 1d quadrature point index in xi1 direction (0,1,2)
   * @param qc The 1d quadrature point index in xi2 direction (0,1,2)
   * @param func The function to call within the support loop.
   * @param params The parameters to pass to @p func.
   */
  template< typename FUNC, typename ... PARAMS >
  GEOSX_HOST_DEVICE
  static void supportLoop( int const qa,
                           int const qb,
                           int const qc,
                           FUNC && func,
                           PARAMS &&... params );

};


template< typename FUNC, typename ... PARAMS >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
H1_Hexahedron_Lagrange2_GaussLegendre3::supportLoop( int const qa,
                                                     int const qb,
                                                     int const qc,
                                                     FUNC && func,
                                                     PARAMS &&... params )
{

/// Options for how to calculate the parent gradients.
  // This option calculates the basis values at the quadrature point for each
  // linear basis index.

  real64 const quadratureCoords[3] = { -quadratureFactor + quadratureFactor * qa,
                                       -quadratureFactor + quadratureFactor * qb,
                                       -quadratureFactor + quadratureFactor * qc };


  // Loop over the linear basis indices in each direction.
  for( int a=0; a<3; ++a )
  {
    for( int b=0; b<3; ++b )
    {
      for( int c=0; c<3; ++c )
      {

        real64 const dNdXi[3] = { LagrangeBasis2::gradient(a, quadratureCoords[0]) * LagrangeBasis2::value(b, quadratureCoords[1]) * LagrangeBasis2::value(c, quadratureCoords[2]),
                                  LagrangeBasis2::value(a, quadratureCoords[0]) * LagrangeBasis2::gradient(b, quadratureCoords[1]) * LagrangeBasis2::value(c, quadratureCoords[2]), 
                                  LagrangeBasis2::value(a, quadratureCoords[0]) * LagrangeBasis2::value(b, quadratureCoords[1]) * LagrangeBasis2::gradient(c, quadratureCoords[2]), 
                                 };
        
        localIndex const nodeIndex = LagrangeBasis2::TensorProduct3D::linearIndex( a, b, c );

        func( dNdXi, nodeIndex, std::forward< PARAMS >( params )... );
      }
    }
  }
}

//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Hexahedron_Lagrange2_GaussLegendre3::calcGradN( localIndex const q,
                                                   real64 const (&X)[numNodes][3],
                                                   real64 (& gradN)[numNodes][3] )
{
  real64 J[3][3] = {{0}};


  int qa, qb, qc;
  LagrangeBasis2::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  applyTransformationToParentGradients( qa, qb, qc, J, gradN );

  //return detJ * weight;
  return detJ;
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 H1_Hexahedron_Lagrange2_GaussLegendre3::
  calcGradN( localIndex const q,
             real64 const (&X)[numNodes][3],
             StackVariables const & GEOSX_UNUSED_PARAM( stack ),
             real64 ( & gradN )[numNodes][3] )
{
  return calcGradN( q, X, gradN );
}

//*************************************************************************************************
#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Hexahedron_Lagrange2_GaussLegendre3::
  jacobianTransformation( int const qa,
                          int const qb,
                          int const qc,
                          real64 const (&X)[numNodes][3],
                          real64 ( & J )[3][3] )
{
  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&X)[numNodes][3],
                                                  real64 (& J)[3][3] )
  {
    real64 const * const GEOSX_RESTRICT Xnode = X[nodeIndex];
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        J[i][j] = J[i][j] + dNdXi[ j ] * Xnode[i];
      }
    }

//    J[0][0] = J[0][0] + dNdXi[0] * Xnode[0];
//    J[0][1] = J[0][1] + dNdXi[1] * Xnode[0];
//    J[0][2] = J[0][2] + dNdXi[2] * Xnode[0];
//    J[1][0] = J[1][0] + dNdXi[0] * Xnode[1];
//    J[1][1] = J[1][1] + dNdXi[1] * Xnode[1];
//    J[1][2] = J[1][2] + dNdXi[2] * Xnode[1];
//    J[2][0] = J[2][0] + dNdXi[0] * Xnode[2];
//    J[2][1] = J[2][1] + dNdXi[1] * Xnode[2];
//    J[2][2] = J[2][2] + dNdXi[2] * Xnode[2];

  }, X, J );
}


//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Hexahedron_Lagrange2_GaussLegendre3::
  applyTransformationToParentGradients( int const qa,
                                        int const qb,
                                        int const qc,
                                        real64 const ( &invJ )[3][3],
                                        real64 (& gradN)[numNodes][3] )
{
  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&invJ)[3][3],
                                                  real64 (& gradN)[numNodes][3] )
  {
//    for( int i = 0; i < 3; ++i )
//    {
//      gradN[nodeIndex][i] = 0.0;
//      for( int j = 0; j < 3; ++j )
//      {
//        gradN[nodeIndex][i] = gradN[nodeIndex][i] + dNdXi[ j ] * invJ[j][i];
//      }
//    }
    // smaller register footprint by manually unrolling the for loops.
    gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] + dNdXi[2] * invJ[2][0];
    gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] + dNdXi[2] * invJ[2][1];
    gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] + dNdXi[2] * invJ[2][2];


  }, invJ, gradN );
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Hexahedron_Lagrange2_GaussLegendre3::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  int qa, qb, qc;
  LagrangeBasis2::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  //return LvArray::tensorOps::determinant< 3 >( J );
  return LvArray::tensorOps::determinant< 3 >( J ) * weight[qa] * weight[qb] * weight[qc];  
}



GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_Hexahedron_Lagrange2_GaussLegendre3::symmetricGradient( int const q,
                                                                real64 const (&invJ)[3][3],
                                                                real64 const (&var)[numNodes][3],
                                                                real64 (& grad)[6] )
{
  int qa, qb, qc;
  LagrangeBasis2::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&invJ)[3][3],
                                                  real64 const (&var)[numNodes][3],
                                                  real64 (& grad)[6] )
  {

    real64 gradN[3] = {0, 0, 0};
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradN[i] = gradN[i] + dNdXi[ j ] * invJ[j][i];
      }
    }

    grad[0] = grad[0] + gradN[0] * var[ nodeIndex ][0];
    grad[1] = grad[1] + gradN[1] * var[ nodeIndex ][1];
    grad[2] = grad[2] + gradN[2] * var[ nodeIndex ][2];
    grad[3] = grad[3] + gradN[2] * var[ nodeIndex ][1] + gradN[1] * var[ nodeIndex ][2];
    grad[4] = grad[4] + gradN[2] * var[ nodeIndex ][0] + gradN[0] * var[ nodeIndex ][2];
    grad[5] = grad[5] + gradN[1] * var[ nodeIndex ][0] + gradN[0] * var[ nodeIndex ][1];
  }, invJ, var, grad );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_Hexahedron_Lagrange2_GaussLegendre3::plusGradNajAij( int const q,
                                                             real64 const (&invJ)[3][3],
                                                             real64 const (&var)[6],
                                                             real64 (& R)[numNodes][3] )
{
  int qa, qb, qc;
  LagrangeBasis2::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc,
               [] GEOSX_HOST_DEVICE
                 ( real64 const (&dNdXi)[3],
                 int const nodeIndex,
                 real64 const (&invJ)[3][3],
                 real64 const (&var)[6],
                 real64 (& R)[numNodes][3] )
  {

    real64 gradN[3] = {0, 0, 0};
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradN[i] = gradN[i] + dNdXi[ j ] * invJ[j][i];
      }
    }
    R[ nodeIndex ][ 0 ] = R[ nodeIndex ][ 0 ] - var[ 0 ] * gradN[ 0 ] - var[ 5 ] * gradN[ 1 ] - var[ 4 ] * gradN[ 2 ];
    R[ nodeIndex ][ 1 ] = R[ nodeIndex ][ 1 ] - var[ 5 ] * gradN[ 0 ] - var[ 1 ] * gradN[ 1 ] - var[ 3 ] * gradN[ 2 ];
    R[ nodeIndex ][ 2 ] = R[ nodeIndex ][ 2 ] - var[ 4 ] * gradN[ 0 ] - var[ 3 ] * gradN[ 1 ] - var[ 2 ] * gradN[ 2 ];
  }, invJ, var, R );
}



GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_Hexahedron_Lagrange2_GaussLegendre3::gradient( int const q,
                                                       real64 const (&invJ)[3][3],
                                                       real64 const (&var)[numNodes][3],
                                                       real64 (& grad)[3][3] )
{
  int qa, qb, qc;
  LagrangeBasis2::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&invJ)[3][3],
                                                  real64 const (&var)[numNodes][3],
                                                  real64 (& grad)[3][3] )
  {
    for( int i = 0; i < 3; ++i )
    {
      real64 gradN=0.0;;
      for( int j = 0; j < 3; ++j )
      {
        gradN = gradN + dNdXi[ j ] * invJ[j][i];
      }
      for( int k = 0; k < 3; ++k )
      {
        grad[k][i] = grad[k][i] + gradN * var[ nodeIndex ][k];
      }
    }
  }, invJ, var, grad );
}

/// @endcond

#if __GNUC__
#pragma GCC diagnostic pop
#endif

#undef PARENT_GRADIENT_METHOD
}
}

#endif //FINITE_ELEMENT_SHAPE_KERNEL
