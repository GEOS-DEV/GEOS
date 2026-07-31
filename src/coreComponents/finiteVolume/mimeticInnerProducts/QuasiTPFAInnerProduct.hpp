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
 * @file QuasiTPFAInnerProduct.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_QUASITPFAINNERPRODUCT_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_QUASITPFAINNERPRODUCT_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductHelpers.hpp"

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @class QuasiTPFAInnerProduct
 *
 * Provides an implementation of a quasi-TPFA inner product in the hybrid FVM solvers
 */
class QuasiTPFAInnerProduct : public MimeticInnerProductBase
{
public:

  /**
   * @brief In a given element, recompute the transmissibility matrix using the quasi TPFA inner product.
   * @param[in] nodePosition the position of the nodes
   * @param[in] transMultiplier the transmissibility multipliers at the mesh faces
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   *
   * @details Reference: K-A Lie, An Introduction to Reservoir Simulation Using MATLAB/GNU Octave (2019)
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
   * @brief Compute the mimetic inner product matrix M in a given element using the quasi TPFA inner product.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] M the output inner product matrix
   *
   * @details Reference: K-A Lie, An Introduction to Reservoir Simulation Using MATLAB/GNU Octave (2019)
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
QuasiTPFAInnerProduct::compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                arrayView1d< real64 const > const & transMultiplier,
                                ArrayOfArraysView< localIndex const > const & faceToNodes,
                                arraySlice1d< localIndex const > const & elemToFaces,
                                arraySlice1d< real64 const > const & elemCenter,
                                real64 const & elemVolume,
                                real64 const (&elemPerm)[ 3 ],
                                real64 const & lengthTolerance,
                                arraySlice2d< real64 > const & transMatrix )
{
  MimeticInnerProductBase::computeParametricInnerProduct< NF >( nodePosition,
                                                                transMultiplier,
                                                                faceToNodes,
                                                                elemToFaces,
                                                                elemCenter,
                                                                elemVolume,
                                                                elemPerm,
                                                                2.0,
                                                                lengthTolerance,
                                                                transMatrix );
}

template< localIndex NF >
GEOS_HOST_DEVICE
void
QuasiTPFAInnerProduct::computeM( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                 ArrayOfArraysView< localIndex const > const & faceToNodes,
                                 arraySlice1d< localIndex const > const & elemToFaces,
                                 arraySlice1d< real64 const > const & elemCenter,
                                 real64 const & elemVolume,
                                 real64 const (&elemPerm)[ 3 ],
                                 real64 const & lengthTolerance,
                                 arraySlice2d< real64 > const & M )
{
  real64 const areaTolerance = lengthTolerance * lengthTolerance;

  // 0) stabilization parameter for quasi-TPFA
  constexpr real64 tParam = 2.0;

  // 1) assemble permeability tensor and its inverse
  real64 K[3][3] = {{0}};
  MimeticInnerProductHelpers::makeFullTensor( elemPerm, K );

  real64 Kinv[3][3] = {{0}};
  for( int i = 0; i < 3; ++i )
  {
    Kinv[i][i] = 1.0 / elemPerm[i];
  }

  // 2) compute C and N
  real64 C[ NF ][ 3 ] = {{ 0 }};
  real64 N[ NF ][ 3 ] = {{ 0 }};

  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 faceCenter[3], faceNormal[3];

    real64 const faceArea =
      computationalGeometry::centroid_3DPolygon(
        faceToNodes[elemToFaces[ifaceLoc]],
        nodePosition,
        faceCenter,
        faceNormal,
        areaTolerance );

    real64 cellToFaceVec[3];
    MimeticInnerProductHelpers::computeCellToFacetVector( cellToFaceVec, faceCenter, elemCenter );
    MimeticInnerProductHelpers::orientNormalOutward( cellToFaceVec, faceNormal );

    for( int d = 0; d < 3; ++d )
    {
      C[ifaceLoc][d] = cellToFaceVec[d];
      N[ifaceLoc][d] = faceArea * faceNormal[d];
    }
  }

  // 3) compute consistency part C K^{-1} C^T / volume
  real64 CKCt[ NF ][ NF ] = {{ 0 }};
  real64 work[ 3 ][ NF ] = {{ 0 }};

  LvArray::tensorOps::Rij_eq_AikBjk< 3, NF, 3 >( work, Kinv, C );
  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, 3 >( CKCt, C, work );
  LvArray::tensorOps::scale< NF, NF >( CKCt, 1.0 / elemVolume );

  // 4) compute W = N K N'
  real64 W[ NF ][ NF ] = {{ 0 }};
  real64 workNK[ 3 ][ NF ] = {{ 0 }};

  LvArray::tensorOps::Rij_eq_AikBjk< 3, NF, 3 >( workNK, K, N );
  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, 3 >( W, N, workNK );

  // 5) build Q from C (orthonormal basis of consistency space)
  real64 q0[ NF ], q1[ NF ], q2[ NF ];
  real64 Qmat[ NF ][ 3 ];

  for( localIndex i = 0; i < NF; ++i )
  {
    q0[i] = C[i][0];
    q1[i] = C[i][1];
    q2[i] = C[i][2];
  }

  MimeticInnerProductHelpers::orthonormalize< NF >( q0, q1, q2, Qmat );

  // 6) compute P = I - Q Q^T
  real64 P[ NF ][ NF ] = {{ 0 }};
  LvArray::tensorOps::addIdentity< NF >( P, -1.0 );
  LvArray::tensorOps::Rij_add_AikAjk< NF, 3 >( P, Qmat );
  LvArray::tensorOps::scale< NF, NF >( P, -1.0 );

  // 7) compute stabilization term (v/t) P diag(W)^{-1} P
  real64 stab[ NF ][ NF ] = {{ 0 }};

  for( localIndex i = 0; i < NF; ++i )
  {
    for( localIndex j = 0; j < NF; ++j )
    {
      real64 val = 0.0;
      for( localIndex k = 0; k < NF; ++k )
      {
        if( LvArray::math::abs( W[k][k] ) > 0 )
        {
          val += P[i][k] * ( 1.0 / W[k][k] ) * P[k][j];
        }
      }
      stab[i][j] = val;
    }
  }

  real64 const scale = elemVolume / tParam;

  // 8) assemble final M
  for( localIndex i = 0; i < NF; ++i )
  {
    for( localIndex j = 0; j < NF; ++j )
    {
      M[i][j] = CKCt[i][j] + scale * stab[i][j];
    }
  }
}

} // end namespace mimeticInnerProduct

} // end namespace geos

#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_QUASITPFAINNERPRODUCT_HPP_
