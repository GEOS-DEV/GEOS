/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MimeticInnerProductBase.hpp
 */

#ifndef GEOSX_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTBASE_HPP_
#define GEOSX_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTBASE_HPP_

#include "meshUtilities/ComputationalGeometry.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductHelpers.hpp"

namespace geosx
{
namespace mimeticInnerProduct
{

/**
 * @class MimeticInnerProductBase
 *
 * Provides management of the mimetic inner products for the hybrid FVM solvers
 */
class MimeticInnerProductBase
{
public:

  MimeticInnerProductBase() = default;

  /**
   * @brief Copy Constructor
   * @param source The object to copy.
   */
  MimeticInnerProductBase( MimeticInnerProductBase const & source ) = default;

  /// Default Move constructor
  MimeticInnerProductBase( MimeticInnerProductBase && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return deleted
   */
  MimeticInnerProductBase & operator=( MimeticInnerProductBase const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return deleted
   */
  MimeticInnerProductBase & operator=( MimeticInnerProductBase && ) = delete;

  /**
   * @brief Destructor
   */
  virtual ~MimeticInnerProductBase() = default;

  /**
   * @brief In a given element, recompute the transmissibility matrix using a consistent inner product.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] tParam parameter used in the transmissibility matrix computations
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   *
   * When tParam = 2, we obtain a scheme that reduces to TPFA
   * on orthogonal meshes, but remains consistent on non-orthogonal meshes
   * When tParam = 6, we obtain the quasi-Raviart-Thomas inner product
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  ComputeParametricInnerProduct( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                 ArrayOfArraysView< localIndex const > const & faceToNodes,
                                 arraySlice1d< localIndex const > const & elemToFaces,
                                 arraySlice1d< real64 const > const & elemCenter,
                                 real64 const & elemVolume,
                                 real64 const (&elemPerm)[ 3 ],
                                 real64 const & tParam,
                                 real64 const & lengthTolerance,
                                 arraySlice2d< real64 > const & transMatrix );

};

template< localIndex NF >
GEOSX_HOST_DEVICE
void
MimeticInnerProductBase::ComputeParametricInnerProduct( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                                        ArrayOfArraysView< localIndex const > const & faceToNodes,
                                                        arraySlice1d< localIndex const > const & elemToFaces,
                                                        arraySlice1d< real64 const > const & elemCenter,
                                                        real64 const & elemVolume,
                                                        real64 const (&elemPerm)[ 3 ],
                                                        real64 const & tParam,
                                                        real64 const & lengthTolerance,
                                                        arraySlice2d< real64 > const & transMatrix )
{
  real64 const areaTolerance = lengthTolerance * lengthTolerance;

  real64 cellToFaceMat[ NF ][ 3 ] = {{ 0 }};
  real64 normalsMat[ NF ][ 3 ] = {{ 0 }};
  real64 permMat[ 3 ][ 3 ] = {{ 0 }};

  real64 work_dimByNumFaces[ 3 ][ NF ] = {{ 0 }};
  real64 worka_numFacesByNumFaces[ NF ][ NF ] = {{ 0 }};
  real64 workb_numFacesByNumFaces[ NF ][ NF ] = {{ 0 }};
  real64 workc_numFacesByNumFaces[ NF ][ NF ] = {{ 0 }};

  real64 q0[ NF ], q1[ NF ], q2[ NF ];

  // 1) fill the matrices cellToFaceMat and normalsMat row by row
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 faceCenter[ 3 ], faceNormal[ 3 ];
    // compute the face geometry data: center, normal, vector from cell center to face center
    real64 const faceArea =
      computationalGeometry::Centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                 nodePosition,
                                                 faceCenter,
                                                 faceNormal,
                                                 areaTolerance );

    q0[ ifaceLoc ] = faceCenter[ 0 ] - elemCenter[ 0 ];
    q1[ ifaceLoc ] = faceCenter[ 1 ] - elemCenter[ 1 ];
    q2[ ifaceLoc ] = faceCenter[ 2 ] - elemCenter[ 2 ];

    real64 const dotProduct = q0[ ifaceLoc ] * faceNormal[ 0 ]
                              + q1[ ifaceLoc ] * faceNormal[ 1 ]
                              + q2[ ifaceLoc ] * faceNormal[ 2 ];
    if( dotProduct < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -faceArea );
    }
    else
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, faceArea );
    }

    normalsMat[ ifaceLoc ][ 0 ] = faceNormal[ 0 ];
    normalsMat[ ifaceLoc ][ 1 ] = faceNormal[ 1 ];
    normalsMat[ ifaceLoc ][ 2 ] = faceNormal[ 2 ];

  }

  // 2) assemble full coefficient tensor from principal axis/components
  MimeticInnerProductHelpers::MakeFullTensor( elemPerm, permMat );

  // 3) compute N K N'
  LvArray::tensorOps::Rij_eq_AikBjk< 3, NF, 3 >( work_dimByNumFaces,
                                                 permMat,
                                                 normalsMat );
  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, 3 >( transMatrix,
                                                  normalsMat,
                                                  work_dimByNumFaces );

  // 4) compute the orthonormalization of the matrix cellToFaceVec
  MimeticInnerProductHelpers::Orthonormalize< NF >( q0, q1, q2, cellToFaceMat );

  // 5) compute P_Q = I - QQ'
  // note: we compute -P_Q and then at 6) ( - P_Q ) D ( - P_Q )
  LvArray::tensorOps::addIdentity< NF >( worka_numFacesByNumFaces, -1 );
  LvArray::tensorOps::Rij_add_AikAjk< NF, 3 >( worka_numFacesByNumFaces,
                                               cellToFaceMat );

  // 6) compute P_Q D P_Q where D = diag(diag(N K N'))
  // 7) compute T = ( N K N' + t U diag(diag(N K N')) U ) / elemVolume
  real64 const scale = tParam / elemVolume;
  for( localIndex i = 0; i < NF; ++i )
  {
    workb_numFacesByNumFaces[ i ][ i ] = scale * transMatrix[ i ][ i ];
  }

  LvArray::tensorOps::Rij_eq_AikBkj< NF, NF, NF >( workc_numFacesByNumFaces,
                                                   workb_numFacesByNumFaces,
                                                   worka_numFacesByNumFaces );
  LvArray::tensorOps::scale< NF, NF >( transMatrix, 1 / elemVolume );
  LvArray::tensorOps::Rij_add_AikBkj< NF, NF, NF >( transMatrix,
                                                    worka_numFacesByNumFaces,
                                                    workc_numFacesByNumFaces );
}

} // end namespace mimeticInnerProduct

} // end namespace geosx

#endif //GEOSX_FINITEVOLUME_MIMETICINNERPRODUCTS_MIMETICINNERPRODUCTBASE_HPP_
