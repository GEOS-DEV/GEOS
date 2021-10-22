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
 * @file H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1QUADRILATERALFACELAGRANGE1GAUSSLEGENDRE2
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1QUADRILATERALFACELAGRANGE1GAUSSLEGENDRE2

#include "FiniteElementBase.hpp"


namespace geosx
{
namespace finiteElement
{

/**
 * This class contains the kernel accessible functions specific to the
 * H1-conforming nodal bilinear quadrilateral face finite element with a
 * 2-point Gauss-Legendre quadrature rule. It is assumed that the indexing
 * for the quadrature points mirrors that of the nodes. Also note that the
 * assumed node ordering is not the standard right-hand-rule used in the
 * literature.
 *
 *                                            =====  ===  ===
 *      2 o----------o 3                      Node   xi0  xi1
 *        |          |                        =====  ===  ===
 *        |          |    xi1                 0      -1   -1
 *        |          |      |                 0       1   -1
 *        |          |      |                 0      -1    1
 *      0 o----------o 1    0---- xi0         0       1    1
 *                                            =====  ===  ===
 *
 */
class H1_QuadrilateralFace_Lagrange1_GaussLegendre2 final : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 4;
  /// The maximum number of support points per element.
  constexpr static localIndex maxSupportPoints = numNodes;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 4;


  virtual ~H1_QuadrilateralFace_Lagrange1_GaussLegendre2() override
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
   * @brief Method to fill a MeshData object.
   * @param nodeManager The node manager.
   * @param edgeManager The edge manager.
   * @param faceManager The face manager.
   * @param cellSubRegion The cell sub-region for which the element has to be initialized.
   * @param meshData MeshData struct to be filled.
   */
  static void fillMeshData( NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            CellElementSubRegion const & cellSubRegion,
                            MeshData & meshData );

  /**
   * @brief Empty setup method.
   * @param cellIndex The index of the cell with respect to the cell sub region.
   * @param meshData MeshData struct filled by @ref fillMeshData.
   * @param stack Object that holds stack variables.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void setupStack( localIndex const & cellIndex,
                          MeshData const & meshData,
                          StackVariables & stack );

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *
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
   * @brief Empty method, here for compatibility with methods that require a stabilization of the
   * grad-grad bilinear form.
   * @tparam MATRIXTYPE The type of @p matrix.
   * @param stack Stack variables as filled by @ref setupStack.
   * @param matrix The matrix that needs to be stabilized.
   */
  template< typename MATRIXTYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void addGradGradStabilization( StackVariables const & stack,
                                        MATRIXTYPE & matrix );

private:
  /// The area of the element in the parent configuration.
  constexpr static real64 parentArea = 4.0;

  /// The weight of each quadrature point.
  constexpr static real64 weight = parentArea / numQuadraturePoints;

  /// The scaling factor specifying the location of the quadrature points
  /// relative to the origin and the outer extent of the element in the
  /// parent space.
  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;

  /**
   * @brief Calculates the linear index for support/quadrature points from ij
   *   coordinates.
   * @param i The index in the xi0 direction (0,1)
   * @param j The index in the xi1 direction (0,1)
   * @return The linear index of the support/quadrature point (0-3)
   */
  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static T linearMap( T const i, T const j )
  {
    return i + 2 * j;
  }

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *        linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the r direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return -1.0 + 2.0 * (a & 1);
  }

  /**
   * @brief Calculate the parent coordinates for the xi1 direction, given the
   *        linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the r direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return -1.0 + ( a & 2 );
  }

};

/// @cond Doxygen_Suppress

GEOSX_FORCE_INLINE
void H1_QuadrilateralFace_Lagrange1_GaussLegendre2::
  fillMeshData( NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
                EdgeManager const & GEOSX_UNUSED_PARAM( edgeManager ),
                FaceManager const & GEOSX_UNUSED_PARAM( faceManager ),
                CellElementSubRegion const & GEOSX_UNUSED_PARAM( cellSubRegion ),
                MeshData & GEOSX_UNUSED_PARAM( meshData ) )
{}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_QuadrilateralFace_Lagrange1_GaussLegendre2::
  setupStack( localIndex const & GEOSX_UNUSED_PARAM( cellIndex ),
              MeshData const & GEOSX_UNUSED_PARAM( meshData ),
              StackVariables & GEOSX_UNUSED_PARAM( stack ) )
{}

template< typename MATRIXTYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_QuadrilateralFace_Lagrange1_GaussLegendre2::
  addGradGradStabilization( StackVariables const & GEOSX_UNUSED_PARAM( stack ),
                            MATRIXTYPE & GEOSX_UNUSED_PARAM( matrix ) )
{}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_QuadrilateralFace_Lagrange1_GaussLegendre2::
  calcN( localIndex const q,
         real64 (& N)[numNodes] )
{
  for( localIndex a=0; a<numNodes; ++a )
  {
    N[a] = 0.25 *
           ( 1 + quadratureFactor*parentCoords0( q )*parentCoords0( a ) ) *
           ( 1 + quadratureFactor*parentCoords1( q )*parentCoords1( a ) );
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_QuadrilateralFace_Lagrange1_GaussLegendre2::
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
H1_QuadrilateralFace_Lagrange1_GaussLegendre2::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  real64 dXdXi[3][2] = {{0}};

  real64 const quadratureCoords[2] = { quadratureFactor *parentCoords0( q ),
                                       quadratureFactor *parentCoords1( q ) };

  real64 const psi0[2] = { 0.5*( 1.0 - quadratureCoords[0] ),
                           0.5*( 1.0 + quadratureCoords[0] ) };
  real64 const psi1[2] = { 0.5*( 1.0 - quadratureCoords[1] ),
                           0.5*( 1.0 + quadratureCoords[1] ) };
  constexpr real64 dpsi[2] = { -0.5, 0.5 };

  for( int a=0; a<2; ++a )
  {
    for( int b=0; b<2; ++b )
    {
      real64 const dNdXi[2] = { dpsi[a] * psi1[b],
                                psi0[a] * dpsi[b] };

      localIndex const nodeIndex = linearMap( a, b );

      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 2; ++j )
        {
          dXdXi[i][j] = dXdXi[i][j] + dNdXi[ j ] * X[nodeIndex][i];
        }
      }
    }
  }

  real64 const detJ = pow( dXdXi[1][0] * dXdXi[2][1] - dXdXi[2][0] * dXdXi[1][1], 2.0 )
                      + pow( dXdXi[2][0] * dXdXi[0][1] - dXdXi[0][0] * dXdXi[2][1], 2.0 )
                      + pow( dXdXi[0][0] * dXdXi[1][1] - dXdXi[1][0] * dXdXi[0][1], 2.0 );

  return sqrt( detJ ) * weight;
}

/// @endcond

}
}
#endif //GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1QUADRILATERALFACELAGRANGE1GAUSSLEGENDRE2
