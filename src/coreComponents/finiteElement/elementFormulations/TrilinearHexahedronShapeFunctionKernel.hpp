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

/**
 * @file TriLinearHexahedronShapeFunctionKernel.hpp
 */

#ifndef GEOSX_CORE_FINITEELEMENT_TRILINEARHEXAHEDRON
#define GEOSX_CORE_FINITEELEMENT_TRILINEARHEXAHEDRON

#include "FiniteElementShapeFunctionKernelBase.hpp"


namespace geosx
{

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
 *                               6___________________ 7
 *                               /.                  /|
 *                              / .                 / |
 *                             /  .                /  |
 *                           4/__________________5/   |
 *                            |   .               |   |
 *                            |   .               |   |
 *                            |   .               |   |
 *                            |   .               |   |
 *                            |   2...............|.../3        xi2
 *                            |  .                |  /          |   xi1
 *                            | .                 | /           |  /
 *                            |.__________________|/            | /
 *                            0                   1             |/____ xi0
 *
 */
class TrilinearHexahedronShapeFunctionKernel : public FiniteElementShapeFunctionKernelBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 8;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 8;


  virtual ~TrilinearHexahedronShapeFunctionKernel() override final
  {}

  virtual localIndex getNumQuadraturePoints() const override final
  {
    return numQuadraturePoints;
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
  static void shapeFunctionValues( localIndex const q,
                                   real64 (&N)[numNodes] )
  {
    for( localIndex a=0; a<numNodes; ++a )
    {
      N[a] = 0.125 *
             ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
             ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) ) *
             ( 1 + quadratureFactor*parentCoords2( q )*parentCoords2( a ) );
    }
  }

  /**
   * @brief Calculate the parent shape function derivatives at a quadrature
   *   point for a given support point.
   * @param q Index of the quadrature point.
   * @param a Index of the support point.
   * @param dNdXi An array to store the parent shape function derivatives for
   *  support point @p a at the quadrature point @p q.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static void
  parentShapeFunctionDerivatives( localIndex const q,
                                  localIndex const a,
                                  real64 (& dNdXi)[3] )
  {
    dNdXi[0] = 0.125 *
               parentCoords0( a ) *
               ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) ) *
               ( 1 + quadratureFactor*parentCoords2( q )*parentCoords2( a ) );
    dNdXi[1] = 0.125 *
               ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
               parentCoords1( a ) *
               ( 1 + quadratureFactor*parentCoords2( q )*parentCoords2( a ) );
    dNdXi[2] = 0.125 *
               ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
               ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) ) *
               parentCoords2( a );
  }


  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param dNdX Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 shapeFunctionDerivatives( localIndex const q,
                                          real64 const (&X)[numNodes][3],
                                          real64 ( &dNdX )[numNodes][3] );


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


private:
  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = 8.0;

  /// The weight of each quadrature point.
  constexpr static real64 weight = parentVolume / numQuadraturePoints;

  /// The scaling factor specifying the location of the quadrature points
  /// relative to the origin and the outer extent of the element in the
  /// parent space.
  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;

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
  constexpr static int linearMap( int const i, int const j, int const k )
  {
    return i + 2 * j + 4 * k;
  }

  /**
   * @brief Calculate the Cartesian index for xi0 given the linear index of a
   *   support point.
   * @param a The linear index of support point
   * @return The Cartesian index of the support point in the xi0 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static int basisIndex0( localIndex const a )
  {
    return (a & 1);
  }

  /**
   * @brief Calculate the Cartesian index for xi1 given the linear index of a
   *   support point.
   * @param a The linear index of support point
   * @return The Cartesian index of the support point in the xi1 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static int basisIndex1( localIndex const a )
  {
    return ( a & 2 ) >> 1;
  }

  /**
   * @brief Calculate the Cartesian index for xi2 given the linear index of a
   *   support point.
   * @param a The linear index of support point
   * @return The Cartesian index of the support point in the xi2 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static int basisIndex2( localIndex const a )
  {
    return ( a & 4 ) >> 2;
  }


  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *   linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the xi0 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return -1.0 + 2.0 * (a & 1);
  }

  /**
   * @brief Calculate the parent coordinates for the xi1 direction, given the
   *   linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the xi1 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return -1.0 + ( a & 2 );
  }

  /**
   * @brief Calculate the parent coordinates for the xi2 direction, given the
   *   linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the xi2 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords2( localIndex const a )
  {
    return -1.0 + 0.5 * ( a & 4 );
  }


  /**
   * @brief Calculates the "Jacobian" transformation matrix/mapping from the
   *   physical space to the parent space.
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
   * @param dNdX Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  GEOSX_HOST_DEVICE
  static void 
  applyJacobianTransformationToShapeFunctionsDerivatives( int const qa,
                                                          int const qb,
                                                          int const qc,
                                                          real64 const ( &invJ )[3][3],
                                                          real64 (& dNdX)[numNodes][3] );


};

//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
TrilinearHexahedronShapeFunctionKernel::shapeFunctionDerivatives( localIndex const q,
                                                                  real64 const (&X)[numNodes][3],
                                                                  real64 (& dNdX)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  int const qa = basisIndex0( q );
  int const qb = basisIndex1( q );
  int const qc = basisIndex2( q );

  jacobianTransformation( qa, qb, qc, X, J );

  real64 const detJ = inverse( J );

  applyJacobianTransformationToShapeFunctionsDerivatives( qa, qb, qc, J, dNdX );

  return detJ * weight;
}

//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TrilinearHexahedronShapeFunctionKernel::
jacobianTransformation( int const qa,
                        int const qb,
                        int const qc,
                        real64 const (&X)[numNodes][3],
                        real64 ( &J )[3][3] )
{

  constexpr static real64 linearBasisAtQuadrature[2] = { 0.5 + 0.5 * quadratureFactor,
                                                         0.5 - 0.5 * quadratureFactor };
  constexpr static real64 psiProduct[3] = { 0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[0],
                                            0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[1],
                                            0.5 * linearBasisAtQuadrature[1]*linearBasisAtQuadrature[1] };
  constexpr static int dpsi[2] = { -1, 1 };

  for( int a=0; a<2; ++a )
  {
    int const qaa = abs( a-qa );
    for( int b=0; b<2; ++b )
    {
      int const qbb = abs( b-qb );
      for( int c=0; c<2; ++c )
      {
        int const qcc = abs( c-qc );
        real64 const dNdXi[3] = { dpsi[a] * psiProduct[ qbb + qcc ],
                                  dpsi[b] * psiProduct[ qaa + qcc ],
                                  dpsi[c] * psiProduct[ qaa + qbb ] };

        localIndex const nodeIndex = linearMap( a, b, c );

        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            J[i][j] = J[i][j] + dNdXi[ j ] * X[nodeIndex][i];
          }
        }
      }
    }
  }
}


//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TrilinearHexahedronShapeFunctionKernel::
applyJacobianTransformationToShapeFunctionsDerivatives( int const qa,
                                                        int const qb,
                                                        int const qc,
                                                        real64 const ( &invJ )[3][3],
                                                        real64 (& dNdX)[numNodes][3] )
{
  constexpr static real64 linearBasisAtQuadrature[2] = { 0.5 + 0.5 * quadratureFactor,
                                                         0.5 - 0.5 * quadratureFactor };
  constexpr static real64 psiProduct[3] = { 0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[0],
                                            0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[1],
                                            0.5 * linearBasisAtQuadrature[1]*linearBasisAtQuadrature[1] };
  constexpr static int dpsi[2] = { -1, 1 };

  for( int a=0; a<2; ++a )
  {
    int const qaa = abs( a-qa );
    for( int b=0; b<2; ++b )
    {
      int const qbb = abs( b-qb );
      for( int c=0; c<2; ++c )
      {
        int const qcc = abs( c-qc );
        real64 const dNdXi[3] = { dpsi[a] * psiProduct[ qbb + qcc ],
                                  dpsi[b] * psiProduct[ qaa + qcc ],
                                  dpsi[c] * psiProduct[ qaa + qbb ] };
        localIndex const nodeIndex = linearMap( a, b, c );
        for( int i = 0; i < 3; ++i )
        {
          dNdX[nodeIndex][i] = 0.0;
          for( int j = 0; j < 3; ++j )
          {
            dNdX[nodeIndex][i] = dNdX[nodeIndex][i] + dNdXi[ j ] * invJ[j][i];
          }
        }
      }
    }
  }
}



GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
TrilinearHexahedronShapeFunctionKernel::
transformedQuadratureWeight( localIndex const q,
                             real64 const (&X)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  int const qa = basisIndex0( q );
  int const qb = basisIndex1( q );
  int const qc = basisIndex2( q );

  jacobianTransformation( qa, qb, qc, X, J );

  return detJ(J);
}


}

#endif //FINITE_ELEMENT_SHAPE_KERNEL
