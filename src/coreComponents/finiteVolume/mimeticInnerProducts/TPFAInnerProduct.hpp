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
 * @file TPFAInnerProduct.hpp
 */

#ifndef GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_TPFAINNERPRODUCT_HPP_
#define GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_TPFAINNERPRODUCT_HPP_

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"

namespace geos
{
namespace mimeticInnerProduct
{

/**
 * @class TPFAInnerProduct
 *
 * Provides an implementation of a TPFA inner product in the hybrid FVM solvers
 */
class TPFAInnerProduct : public MimeticInnerProductBase
{
public:

  /**
   * @brief In a given element, recompute the transmissibility matrix in a cell using TPFA.
   * @param[in] nodePosition the position of the nodes
   * @param[in] transMultiplier the transmissibility multipliers at the mesh faces
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
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
   * @brief Compute the mimetic inner product matrix M in a given element using TPFA.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] M the output inner product matrix
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
TPFAInnerProduct::compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                           arrayView1d< real64 const > const & transMultiplier,
                           ArrayOfArraysView< localIndex const > const & faceToNodes,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           arraySlice1d< real64 const > const & elemCenter,
                           real64 const & elemVolume,
                           real64 const (&elemPerm)[ 3 ],
                           real64 const & lengthTolerance,
                           arraySlice2d< real64 > const & transMatrix )
{
  GEOS_UNUSED_VAR( elemVolume );

  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance;

  // 0) assemble full coefficient tensor from principal axis/components
  real64 permTensor[ 3 ][ 3 ] = {{ 0 }};
  MimeticInnerProductHelpers::makeFullTensor( elemPerm, permTensor );

  // we are ready to compute the transmissibility matrix
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 const mult = transMultiplier[elemToFaces[ifaceLoc]];

    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      // for now, TPFA trans
      if( ifaceLoc == jfaceLoc )
      {
        real64 faceCenter[ 3 ], faceNormal[ 3 ], faceConormal[ 3 ], cellToFaceVec[ 3 ];
        // 1) compute the face geometry data: center, normal, vector from cell center to face center
        real64 const faceArea =
          computationalGeometry::centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                     nodePosition,
                                                     faceCenter,
                                                     faceNormal,
                                                     areaTolerance );

        MimeticInnerProductHelpers::computeCellToFacetVector( cellToFaceVec, faceCenter, elemCenter );
        MimeticInnerProductHelpers::orientNormalOutward( cellToFaceVec, faceNormal );

        real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

        LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, elemPerm, faceNormal );

        // 3) compute the one-sided face transmissibility
        real64 const halfTrans =
          LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal ) * mult * faceArea / c2fDistance;
        // T := sign(T) * max(|T|, tol), so that |T| >= tol and sign(T) is preserved
        transMatrix[ifaceLoc][jfaceLoc] = ( halfTrans < 0.0 ? -1.0 : 1.0 ) *
                                          LvArray::math::max( LvArray::math::abs( halfTrans ), weightTolerance );
      }
      else
      {
        transMatrix[ifaceLoc][jfaceLoc] = 0;
      }
    }
  }
}

template< localIndex NF >
GEOS_HOST_DEVICE
void
TPFAInnerProduct::computeM( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                            ArrayOfArraysView< localIndex const > const & faceToNodes,
                            arraySlice1d< localIndex const > const & elemToFaces,
                            arraySlice1d< real64 const > const & elemCenter,
                            real64 const & elemVolume,
                            real64 const (&elemPerm)[ 3 ],
                            real64 const & lengthTolerance,
                            arraySlice2d< real64 > const & M )
{
  GEOS_UNUSED_VAR( elemVolume );

  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance;

  // initialize M to zero
  LvArray::tensorOps::fill< NF, NF >( M, 0.0 );

  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 faceCenter[3], faceNormal[3], faceConormal[3], cellToFaceVec[3];

    // 1) face geometry
    real64 const faceArea =
      computationalGeometry::centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                 nodePosition,
                                                 faceCenter,
                                                 faceNormal,
                                                 areaTolerance );

    MimeticInnerProductHelpers::computeCellToFacetVector( cellToFaceVec, faceCenter, elemCenter );
    MimeticInnerProductHelpers::orientNormalOutward( cellToFaceVec, faceNormal );

    real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

    // 2) K * n (TPFA assumes diagonal K)
    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, elemPerm, faceNormal );

    // 3) compute T_ii
    real64 Tii = LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal ) * faceArea / c2fDistance;

    // T := max(|T|, tol), so that M_ii = 1/T is well-defined with M_ii <= 1/tol
    Tii = LvArray::math::max( LvArray::math::abs( Tii ), weightTolerance );

    // 4) M = |T|^{-1}
    M[ifaceLoc][ifaceLoc] = 1.0 / Tii;
  }
}


} // end namespace mimeticInnerProduct

} // end namespace geos


#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_TPFAINNERPRODUCT_HPP_
