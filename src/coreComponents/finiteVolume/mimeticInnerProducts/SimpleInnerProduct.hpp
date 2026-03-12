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
 * @file SimpleInnerProduct.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_SIMPLEINNERPRODUCT_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_SIMPLEINNERPRODUCT_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductHelpers.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include <iostream>
#include <iomanip>


namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @class SimpleInnerProduct
 *
 * Provides an implementation of a Simple inner product in the hybrid FVM solvers
 */
class SimpleInnerProduct : public MimeticInnerProductBase
{
public:

  /**
   * @brief In a given element, recompute the transmissibility matrix in a cell using the Simple inner product
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
SimpleInnerProduct::compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                             arrayView1d< real64 const > const & transMultiplier,
                             ArrayOfArraysView< localIndex const > const & faceToNodes,
                             arraySlice1d< localIndex const > const & elemToFaces,
                             arraySlice1d< real64 const > const & elemCenter,
                             real64 const & elemVolume,
                             real64 const (&elemPerm)[ 3 ],
                             real64 const & lengthTolerance,
                             arraySlice2d< real64 > const & transMatrix )
{
  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightToleranceInv = 1e30 / lengthTolerance;

  real64 cellToFaceMat[ NF ][ 3 ] = {{ 0 }};
  real64 normalsMat[ NF ][ 3 ] = {{ 0 }};
  real64 permMat[ 3 ][ 3 ] = {{ 0 }};

  real64 work_dimByNumFaces[ 3 ][ NF ] = {{ 0 }};
  real64 worka_numFacesByNumFaces[ NF ][ NF ] = {{ 0 }};
  real64 workb_numFacesByNumFaces[ NF ][ NF ] = {{ 0 }};
  real64 workc_numFacesByNumFaces[ NF ][ NF ] = {{ 0 }};

  real64 tpTransInv[ NF ] = { 0.0 };

  real64 q0[ NF ], q1[ NF ], q2[ NF ];
  real64 faceArea[ NF ];


  // 0) assemble full coefficient tensor from principal axis/components
  MimeticInnerProductHelpers::makeFullTensor( elemPerm, permMat );

  // 1) fill the matrices cellToFaceMat and normalsMat row by row
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 faceCenter[ 3 ], faceNormal[ 3 ], cellToFaceVec[ 3 ];
    // compute the face geometry data: center, normal, vector from cell center to face center
    faceArea[ ifaceLoc ] =
      computationalGeometry::centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                 nodePosition,
                                                 faceCenter,
                                                 faceNormal,
                                                 areaTolerance );

    LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
    LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter );

    q0[ ifaceLoc ] = faceArea[ ifaceLoc ] * cellToFaceVec[ 0 ];
    q1[ ifaceLoc ] = faceArea[ ifaceLoc ] * cellToFaceVec[ 1 ];
    q2[ ifaceLoc ] = faceArea[ ifaceLoc ] * cellToFaceVec[ 2 ];

    if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    // the two-point transmissibility is computed to computed here because it is needed
    // in the implementation of the transmissibility multiplier (see below)
    // TODO: see what it would take to bring the (harmonically averaged) two-point trans here
    real64 diagEntry = 0.0;
    MimeticInnerProductHelpers::computeInvTPFATransWithMultiplier< NF >( elemPerm,
                                                                         faceNormal,
                                                                         faceArea[ifaceLoc],
                                                                         transMultiplier[elemToFaces[ifaceLoc]],
                                                                         weightToleranceInv,
                                                                         cellToFaceVec,
                                                                         diagEntry );
    tpTransInv[ifaceLoc] = diagEntry;

    LvArray::tensorOps::scale< 3 >( faceNormal, faceArea[ ifaceLoc ] );
    normalsMat[ ifaceLoc ][ 0 ] = faceNormal[ 0 ];
    normalsMat[ ifaceLoc ][ 1 ] = faceNormal[ 1 ];
    normalsMat[ ifaceLoc ][ 2 ] = faceNormal[ 2 ];

  }

  // 2) compute the stabilization coefficient
  real64 const tParam = 2 * LvArray::tensorOps::trace< 3 >( permMat );

  // 3) compute N K N'
  LvArray::tensorOps::Rij_eq_AikBjk< 3, NF, 3 >( work_dimByNumFaces,
                                                 permMat,
                                                 normalsMat );
  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, 3 >( transMatrix,
                                                  normalsMat,
                                                  work_dimByNumFaces );

  // 4) compute the orthonormalization of the matrix cellToFaceVec
  MimeticInnerProductHelpers::orthonormalize< NF >( q0, q1, q2, cellToFaceMat );

  // 5) compute P_Q = I - QQ'
  // note: we compute -P_Q here
  LvArray::tensorOps::addIdentity< NF >( worka_numFacesByNumFaces, -1 );
  LvArray::tensorOps::Rij_add_AikAjk< NF, 3 >( worka_numFacesByNumFaces,
                                               cellToFaceMat );

  // 6) compute D P_Q D where D = diag( faceArea )
  // 7) compute T = ( N K N' + t D P_Q D ) / elemVolume
  for( localIndex i = 0; i < NF; ++i )
  {
    workb_numFacesByNumFaces[ i ][ i ] = faceArea[ i ];
  }

  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, NF >( workc_numFacesByNumFaces,
                                                   workb_numFacesByNumFaces,
                                                   worka_numFacesByNumFaces );


  LvArray::tensorOps::scale< NF, NF >( transMatrix, 1 / elemVolume );
  LvArray::tensorOps::scale< NF, NF >( workc_numFacesByNumFaces, -tParam / elemVolume );
  LvArray::tensorOps::Rij_add_AikBkj< NF, NF, NF >( transMatrix,
                                                    workc_numFacesByNumFaces,
                                                    workb_numFacesByNumFaces );

  // 7) incorporate the transmissbility multipliers
  // Ref: Nilsen, H. M., J. R. Natvig, and K.-A Lie.,
  // "Accurate modeling of faults by multipoint, mimetic, and mixed methods." SPEJ

  if( !isZero( LvArray::tensorOps::l2NormSquared< NF >( tpTransInv ) ) )
  {
    MimeticInnerProductHelpers::computeTransMatrixWithMultipliers< NF >( tpTransInv,
                                                                         transMatrix );
  }

}

template< localIndex NF >
GEOS_HOST_DEVICE
void
SimpleInnerProduct::computeM( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                              ArrayOfArraysView< localIndex const > const & faceToNodes,
                              arraySlice1d< localIndex const > const & elemToFaces,
                              arraySlice1d< real64 const > const & elemCenter,
                              real64 const & elemVolume,
                              real64 const (&elemPerm)[ 3 ],
                              real64 const & lengthTolerance,
                              arraySlice2d< real64 > const & M )
{
  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  // 0) Helper: diagonal inverse of [k1, k2, k3]
  real64 Kinv[3][3] = {{ 0 }};
  for( int i =0; i < 3; ++i )
  {
    Kinv[i][i] = 1.0 / elemPerm[i];
  }

  // 1) Compute C, N, A
  real64 C[ NF ][ 3 ] = {{ 0 }};
  real64 N[ NF ][ 3 ] = {{ 0 }};
  real64 A[ NF ] = { 0.0 };

  for( localIndex i = 0; i < NF; ++i )
  {
    real64 fCenter[3], fNormal[3];
    real64 area = computationalGeometry::centroid_3DPolygon(
      faceToNodes[elemToFaces[i]], nodePosition, fCenter, fNormal, areaTolerance );

    real64 vec_cf[3];
    LvArray::tensorOps::copy< 3 >( vec_cf, fCenter );
    LvArray::tensorOps::subtract< 3 >( vec_cf, elemCenter );
    if( LvArray::tensorOps::AiBi< 3 >( vec_cf, fNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( fNormal, -1.0 );
    }

    for( int d = 0; d < 3; ++d )
    {
      C[i][d] = vec_cf[d];
      N[i][d] = fNormal[d] * area;
    }

    A[i] = area;
  }

  // 2) Compute C K^{-1} C^T / volume
  real64 CKCt[ NF ][ NF ] = {{0}};
  real64 work3xNF[ 3 ][ NF ] = {{0}};

  LvArray::tensorOps::Rij_eq_AikBjk< 3, NF, 3 >( work3xNF, Kinv, C );
  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, 3 >( CKCt, C, work3xNF );
  LvArray::tensorOps::scale< NF, NF >( CKCt, 1.0 / elemVolume );

  // 3) Q = orth(N / A)
  real64 q0[ NF ], q1[ NF ], q2[ NF ];
  real64 Qmat[ NF ][ 3 ];
  for( localIndex i = 0; i < NF; ++i )
  {
    q0[i] = N[i][0] / A[i];
    q1[i] = N[i][1] / A[i];
    q2[i] = N[i][2] / A[i];
  }

  MimeticInnerProductHelpers::orthonormalize< NF >( q0, q1, q2, Qmat );

  // 4) U = I - Q Q^T
  real64 U[ NF ][ NF ] = {{0}};
  LvArray::tensorOps::addIdentity< NF >( U, -1.0 );
  LvArray::tensorOps::Rij_add_AikAjk< NF, 3 >( U, Qmat );
  LvArray::tensorOps::scale< NF, NF >( U, -1.0 );

  // 5) A^{-1} * U * A^{-1}
  real64 invA[ NF ];
  for( localIndex i = 0; i < NF; ++i )
  {
    invA[i] = 1.0 / A[i];
  }

  // temp = A^{-1} * U
  real64 temp[ NF ][ NF ];
  for( localIndex i = 0; i < NF; ++i )
  {
    for( localIndex j = 0; j < NF; ++j )
    {
      temp[i][j] = invA[i] * U[i][j];
    }
  }

  // Usc = A^{-1} * U * A^{-1}
  real64 Usc[ NF ][ NF ];
  for( localIndex i = 0; i < NF; ++i )
  {
    for( localIndex j = 0; j < NF; ++j )
    {
      Usc[i][j] = temp[i][j] * invA[j];
    }
  }

  // 6) M = CKCt + (v / t) * Usc
  real64 tParam = 2.0 * ( elemPerm[0] + elemPerm[1] + elemPerm[2] );
  real64 scale = elemVolume / tParam;

  for( localIndex i = 0; i < NF; ++i )
  {
    for( localIndex j = 0; j < NF; ++j )
    {
      M[i][j] = CKCt[i][j] + scale * Usc[i][j];
    }
  }
}


} // end namespace mimeticInnerProduct

} // end namespace geos


#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_SIMPLEINNERPRODUCT_HPP_
