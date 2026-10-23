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
 * @file RTInnerProduct.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_RTINNERPRODUCT_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_RTINNERPRODUCT_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductHelpers.hpp"

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @class RTInnerProduct
 *
 * Provides the mimetic inner product that reproduces the lowest-order Raviart-Thomas mass matrix
 * exactly on simplices, and extends it consistently to general polyhedra
 */
class RTInnerProduct : public MimeticInnerProductBase
{
public:

  /**
   * @brief In a given element, recompute the transmissibility matrix as the exact inverse of the
   *        RT inner product matrix M, so that the hybrid and mixed forms are exactly dual.
   * @param[in] nodePosition the position of the nodes
   * @param[in] transMultiplier the transmissibility multipliers at the mesh faces
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output W = M^{-1}, with the face multipliers applied symmetrically
   */
  template< localIndex NF >
  GEOS_HOST_DEVICE
  static void
  compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           arrayView1d< real64 const > const & transMultiplier,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const & elemToFaces,
           arraySlice1d< real64 const > const & elemCenter,
           real64 const & elemVolume,
           real64 const (&elemPerm)[ 3 ],
           real64 const & lengthTolerance,
           arraySlice2d< real64 > const & transMatrix );

  /**
   * @brief Compute the mimetic inner product matrix M in a given element.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] M the output inner product matrix
   *
   * @details M = C K^{-1} C^T / V + s ( I - Q Q^T ) with Q an orthonormal basis of the range of the
   * K-weighted area normals; s = tr( C K^{-1} C^T / V ) / (d+2) is the unique scale matching the
   * conforming RT0 element on simplices.
   */
  template< localIndex NF >
  GEOS_HOST_DEVICE
  static void
  computeM( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
            ArrayOfArraysView< localIndex const > const & faceToNodes,
            arraySlice1d< localIndex const > const & elemToFaces,
            arraySlice1d< real64 const > const & elemCenter,
            real64 const & elemVolume,
            real64 const (&elemPerm)[ 3 ],
            real64 const & lengthTolerance,
            arraySlice2d< real64 > const & M );
};

template< localIndex NF >
GEOS_HOST_DEVICE
void
RTInnerProduct::computeM( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                          ArrayOfArraysView< localIndex const > const & faceToNodes,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          arraySlice1d< real64 const > const & elemCenter,
                          real64 const & elemVolume,
                          real64 const (&elemPerm)[ 3 ],
                          real64 const & lengthTolerance,
                          arraySlice2d< real64 > const & M )
{
  real64 const areaTolerance = lengthTolerance * lengthTolerance;

  // 1) compute C, N, face areas and the consistency term M1 = C K^{-1} C^T / V
  real64 C[ NF ][ 3 ] = {{ 0 }};
  real64 N[ NF ][ 3 ] = {{ 0 }};
  real64 faceArea[ NF ] = { 0.0 };
  real64 M1[ NF ][ NF ] = {{ 0 }};

  MimeticInnerProductHelpers::computeCellToFaceGeometry< NF >( nodePosition, faceToNodes, elemToFaces,
                                                               elemCenter, areaTolerance, C, N, faceArea );
  MimeticInnerProductHelpers::computeConsistencyTerm< NF >( C, elemVolume, elemPerm, M1 );

  // 2) build Q from the K-weighted area normals N K, so I - Q Q^T projects onto ker( (N K)^T )
  real64 q0[ NF ], q1[ NF ], q2[ NF ];
  real64 Qmat[ NF ][ 3 ];
  for( localIndex i = 0; i < NF; ++i )
  {
    q0[i] = N[i][0] * elemPerm[0];
    q1[i] = N[i][1] * elemPerm[1];
    q2[i] = N[i][2] * elemPerm[2];
  }
  MimeticInnerProductHelpers::orthonormalize< NF >( q0, q1, q2, Qmat );

  // 3) compute P = I - Q Q^T
  real64 P[ NF ][ NF ] = {{ 0 }};
  LvArray::tensorOps::addIdentity< NF >( P, -1.0 );
  LvArray::tensorOps::Rij_add_AikAjk< NF, 3 >( P, Qmat );
  LvArray::tensorOps::scale< NF, NF >( P, -1.0 );

  // 4) s = tr(M1)/(d+2), the unique RT0-matching stabilization scale on simplices (d = 3)
  real64 const s = LvArray::tensorOps::trace< NF >( M1 ) / 5.0;

  // 5) assemble M = M1 + s P
  for( localIndex i = 0; i < NF; ++i )
  {
    for( localIndex j = 0; j < NF; ++j )
    {
      M[i][j] = M1[i][j] + s * P[i][j];
    }
  }
}

template< localIndex NF >
GEOS_HOST_DEVICE
void
RTInnerProduct::compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                         arrayView1d< real64 const > const & transMultiplier,
                         ArrayOfArraysView< localIndex const > const & faceToNodes,
                         arraySlice1d< localIndex const > const & elemToFaces,
                         arraySlice1d< real64 const > const & elemCenter,
                         real64 const & elemVolume,
                         real64 const (&elemPerm)[ 3 ],
                         real64 const & lengthTolerance,
                         arraySlice2d< real64 > const & transMatrix )
{
  // 1) assemble M
  computeM< NF >( nodePosition, faceToNodes, elemToFaces, elemCenter,
                  elemVolume, elemPerm, lengthTolerance, transMatrix );

  // 2) W = M^{-1} by Gauss-Jordan on the SPD local matrix
  real64 A[ NF ][ NF ];
  real64 W[ NF ][ NF ] = {{ 0 }};
  for( localIndex i = 0; i < NF; ++i )
  {
    W[i][i] = 1.0;
    for( localIndex j = 0; j < NF; ++j )
    {
      A[i][j] = transMatrix[i][j];
    }
  }
  for( localIndex k = 0; k < NF; ++k )
  {
    real64 const invPivot = 1.0 / A[k][k];
    for( localIndex j = 0; j < NF; ++j )
    {
      A[k][j] *= invPivot;
      W[k][j] *= invPivot;
    }
    for( localIndex i = 0; i < NF; ++i )
    {
      if( i == k )
        continue;
      real64 const factor = A[i][k];
      for( localIndex j = 0; j < NF; ++j )
      {
        A[i][j] -= factor * A[k][j];
        W[i][j] -= factor * W[k][j];
      }
    }
  }

  // 3) apply the face transmissibility multipliers symmetrically
  for( localIndex i = 0; i < NF; ++i )
  {
    for( localIndex j = 0; j < NF; ++j )
    {
      transMatrix[i][j] = W[i][j] * LvArray::math::sqrt( transMultiplier[elemToFaces[i]] * transMultiplier[elemToFaces[j]] );
    }
  }
}

} // end namespace mimeticInnerProduct

} // end namespace geos

#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_RTINNERPRODUCT_HPP_
