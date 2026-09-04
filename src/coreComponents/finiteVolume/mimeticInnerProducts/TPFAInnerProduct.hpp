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

private:

  /**
   * @brief Compute the one-sided (half) TPFA transmissibility of a local face, k_n A / d.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] faceIndex the index of the face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] areaTolerance the tolerance used in the face area calculations
   * @return the one-sided transmissibility of the face
   */
  GEOS_HOST_DEVICE
  inline
  static real64
  computeOneSidedTrans( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                        ArrayOfArraysView< localIndex const > const & faceToNodes,
                        localIndex const faceIndex,
                        arraySlice1d< real64 const > const & elemCenter,
                        real64 const (&elemPerm)[ 3 ],
                        real64 const & areaTolerance )
  {
    real64 faceCenter[ 3 ], faceNormal[ 3 ], faceConormal[ 3 ], cellToFaceVec[ 3 ];

    real64 const faceArea =
      computationalGeometry::centroid_3DPolygon( faceToNodes[faceIndex],
                                                 nodePosition,
                                                 faceCenter,
                                                 faceNormal,
                                                 areaTolerance );

    MimeticInnerProductHelpers::computeCellToFacetVector( cellToFaceVec, faceCenter, elemCenter );
    MimeticInnerProductHelpers::orientNormalOutward( cellToFaceVec, faceNormal );

    real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

    // TPFA assumes diagonal K
    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, elemPerm, faceNormal );

    return LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal ) * faceArea / c2fDistance;
  }
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

  // we are ready to compute the transmissibility matrix
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 const mult = transMultiplier[elemToFaces[ifaceLoc]];

    real64 const halfTrans = computeOneSidedTrans( nodePosition, faceToNodes, elemToFaces[ifaceLoc],
                                                   elemCenter, elemPerm, areaTolerance );

    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      // for now, TPFA trans
      if( ifaceLoc == jfaceLoc )
      {
        transMatrix[ifaceLoc][jfaceLoc] = LvArray::math::max( mult * halfTrans, weightTolerance );
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
    // 1) one-sided transmissibility T_ii, shared with compute()
    real64 const Tii = LvArray::math::max( computeOneSidedTrans( nodePosition, faceToNodes, elemToFaces[ifaceLoc],
                                                                 elemCenter, elemPerm, areaTolerance ),
                                           weightTolerance );

    // 2) M = |T|^{-1}
    M[ifaceLoc][ifaceLoc] = 1.0 / LvArray::math::abs( Tii );
  }
}


} // end namespace mimeticInnerProduct

} // end namespace geos


#endif //GEOS_FINITEVOLUME_MIMETICINNERPRODUCTS_TPFAINNERPRODUCT_HPP_
