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
 * @file HybridFVMInnerProduct.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

namespace HybridFVMInnerProduct
{

/******************************** Helpers ********************************/

/**
 * @struct HybridFVMInnerProductHelper
 * @brief Helper struct handling inner product for hybrid finite volume schemes.
 */
struct HybridFVMInnerProductHelper
{

  /**
   * @brief Create a full tensor from an array.
   * @param[in] values the input array
   * @param[out] result the full tensor
   */
  GEOSX_HOST_DEVICE
  static
  void MakeFullTensor( R1Tensor const & values,
                       arraySlice2d< real64 > const & result )
  {
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        result( i, j ) = 0.0;
      }
    }
    result[ 0 ][ 0 ] = values[ 0 ];
    result[ 1 ][ 2 ] = values[ 1 ];
    result[ 1 ][ 2 ] = values[ 2 ];
  }


  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static
  void Orthonormalize( arraySlice1d< real64 > const & q0,
                       arraySlice1d< real64 > const & q1,
                       arraySlice1d< real64 > const & q2,
                       arraySlice2d< real64 > const & cellToFaceMat )
  {
    // modified Gram-Schmidt algorithm

    // q0
    BlasLapackLA::vectorScale( 1.0/BlasLapackLA::vectorNorm2( q0 ), q0 );

    // q1
    real64 const q0Dotq1 = BlasLapackLA::vectorDot( q0, q1 );
    BlasLapackLA::vectorVectorAdd( q0, q1, -q0Dotq1 );
    BlasLapackLA::vectorScale( 1.0/BlasLapackLA::vectorNorm2( q1 ), q1 );

    // q2
    real64 const q0Dotq2 = BlasLapackLA::vectorDot( q0, q2 );
    BlasLapackLA::vectorVectorAdd( q0, q2, -q0Dotq2 );
    real64 const q1Dotq2 = BlasLapackLA::vectorDot( q1, q2 );
    BlasLapackLA::vectorVectorAdd( q1, q2, -q1Dotq2 );
    BlasLapackLA::vectorScale( 1.0/BlasLapackLA::vectorNorm2( q2 ), q2 );

    for( localIndex i = 0; i < NF; ++i )
    {
      cellToFaceMat( i, 0 ) = q0( i );
      cellToFaceMat( i, 1 ) = q1( i );
      cellToFaceMat( i, 2 ) = q2( i );
    }
  }

};

/******************************** TPFA Kernel ********************************/

/**
 * @struct TPFACellInnerProductKernel
 * @brief Struct handling cell inner product in the TPFA scheme.
 */
struct TPFACellInnerProductKernel
{

  /**
   * @brief In a given element, recompute the transmissibility matrix in a cell using TPFA.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           R1Tensor const & elemPerm,
           real64 const & lengthTolerance,
           arraySlice2d< real64 > const & transMatrix )
  {
    localIndex constexpr dim = 3;

    real64 const areaTolerance = lengthTolerance * lengthTolerance;
    real64 const weightTolerance = 1e-30 * lengthTolerance;

    stackArray2d< real64, dim *dim > permTensor( dim, dim );

    // we are ready to compute the transmissibility matrix
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {
      for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        // for now, TPFA trans
        if( ifaceLoc == jfaceLoc )
        {
          real64 faceCenter[ 3 ], faceNormal[ 3 ], faceConormal[ 3 ], cellToFaceVec[ 3 ];

          // 1) compute the face geometry data: center, normal, vector from cell center to face center
          real64 const faceArea =
            computationalGeometry::Centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                       nodePosition,
                                                       faceCenter,
                                                       faceNormal,
                                                       areaTolerance );

          LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
          // TODO: Change to LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter );
          cellToFaceVec[ 0 ] -= elemCenter[ 0 ];
          cellToFaceVec[ 1 ] -= elemCenter[ 1 ];
          cellToFaceVec[ 2 ] -= elemCenter[ 2 ];

          if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
          {
            LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
          }

          real64 const c2fDistance = LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

          // 2) assemble full coefficient tensor from principal axis/components
          for( localIndex i = 0; i < dim; ++i )
          {
            for( localIndex j = 0; j < dim; ++j )
            {
              permTensor( i, j ) = 0;
            }
          }
          HybridFVMInnerProductHelper::MakeFullTensor( elemPerm, permTensor );

          // TODO: take symmetry into account to optimize this
          // TODO: Change to LvArray::tensorOps::AiBi
          faceConormal[ 0 ] = elemPerm[ 0 ] * faceNormal[ 0 ];
          faceConormal[ 1 ] = elemPerm[ 1 ] * faceNormal[ 1 ];
          faceConormal[ 2 ] = elemPerm[ 2 ] * faceNormal[ 2 ];

          // 3) compute the one-sided face transmissibility
          transMatrix[ifaceLoc][jfaceLoc]  = LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal );
          transMatrix[ifaceLoc][jfaceLoc] *= faceArea / c2fDistance;
          transMatrix[ifaceLoc][jfaceLoc]  = LvArray::max( transMatrix[ifaceLoc][jfaceLoc], weightTolerance );
        }
        else
        {
          transMatrix[ifaceLoc][jfaceLoc] = 0;
        }
      }
    }
  }


};


/******************************** Quasi TPFA Kernel ********************************/

/**
 * @struct QTPFACellInnerProductKernel
 * @brief Struct handling cell inner product in the quasi TPFA scheme.
 */
struct QTPFACellInnerProductKernel
{

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
   * @param[in] orthonormalizeWithSVD flag to choose is use orthonormalizing through SVD
   * @param[inout] transMatrix the output
   *
   * When tParam = 2, we obtain a scheme that reduces to TPFA
   * on orthogonal meshes, but remains consistent on non-orthogonal meshes
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           real64 const & elemVolume,
           R1Tensor const & elemPerm,
           real64 const & tParam,
           real64 const & lengthTolerance,
           arraySlice2d< real64 > const & transMatrix )
  {
    localIndex constexpr dim = 3;

    real64 const areaTolerance = lengthTolerance * lengthTolerance;

    stackArray2d< real64, NF *dim > cellToFaceMat( NF, dim );
    stackArray2d< real64, NF *dim > normalsMat( NF, dim );
    stackArray2d< real64, dim *dim > permMat( dim, dim );

    stackArray2d< real64, dim *NF > work_dimByNumFaces( dim, NF );
    stackArray2d< real64, NF *NF > worka_numFacesByNumFaces( NF, NF );
    stackArray2d< real64, NF *NF > workb_numFacesByNumFaces( NF, NF );
    stackArray2d< real64, NF *NF > workc_numFacesByNumFaces( NF, NF );

    stackArray1d< real64, NF > q0( NF );
    stackArray1d< real64, NF > q1( NF );
    stackArray1d< real64, NF > q2( NF );

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

      q0( ifaceLoc ) = faceCenter[ 0 ] - elemCenter( 0 );
      q1( ifaceLoc ) = faceCenter[ 1 ] - elemCenter( 1 );
      q2( ifaceLoc ) = faceCenter[ 2 ] - elemCenter( 2 );

      real64 const dotProduct = q0( ifaceLoc ) * faceNormal[ 0 ]
                                + q1( ifaceLoc ) * faceNormal[ 1 ]
                                + q2( ifaceLoc ) * faceNormal[ 2 ];
      if( dotProduct < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
      }

      normalsMat( ifaceLoc, 0 ) = faceArea * faceNormal[ 0 ];
      normalsMat( ifaceLoc, 1 ) = faceArea * faceNormal[ 1 ];
      normalsMat( ifaceLoc, 2 ) = faceArea * faceNormal[ 2 ];

    }

    // 2) assemble full coefficient tensor from principal axis/components
    for( localIndex i = 0; i < dim; ++i )
    {
      for( localIndex j = 0; j < dim; ++j )
      {
        permMat( i, j ) = 0;
      }
    }
    HybridFVMInnerProductHelper::MakeFullTensor( elemPerm, permMat );

    // 3) compute N K N'
    BlasLapackLA::matrixMatrixTMultiply( permMat,
                                         normalsMat,
                                         work_dimByNumFaces );
    BlasLapackLA::matrixMatrixMultiply( normalsMat,
                                        work_dimByNumFaces,
                                        transMatrix );

    // 4) compute the orthonormalization of the matrix cellToFaceVec
    HybridFVMInnerProductHelper::Orthonormalize< NF >( q0, q1, q2, cellToFaceMat );

    // 5) compute P_Q = I - QQ'
    for( localIndex i = 0; i < NF; ++i )
    {
      worka_numFacesByNumFaces( i, i ) = 1;
      for( localIndex j = 0; j < NF; ++j )
      {
        if( i == j )
        {
          worka_numFacesByNumFaces( i, j ) = 1;
        }
        else
        {
          worka_numFacesByNumFaces( i, j ) = 0;
        }
      }
    }
    BlasLapackLA::matrixMatrixTMultiply( cellToFaceMat,
                                         cellToFaceMat,
                                         worka_numFacesByNumFaces,
                                         -1, 1 );

    // 6) compute P_Q D P_Q where D = diag(diag(N K N'))
    // 7) compute T = ( N K N' + t U diag(diag(N K N')) U ) / elemVolume
    // Note that 7) is done at the last call to matrixMatrixMultiply
    for( localIndex i = 0; i < NF; ++i )
    {
      for( localIndex j = 0; j < NF; ++j )
      {
        if( i == j )
        {
          workb_numFacesByNumFaces( i, j ) = transMatrix( i, j );
        }
        else
        {
          workb_numFacesByNumFaces( i, j ) = 0;
        }
      }
    }
    BlasLapackLA::matrixMatrixMultiply( workb_numFacesByNumFaces,
                                        worka_numFacesByNumFaces,
                                        workc_numFacesByNumFaces );
    real64 const elemVolumeInv = 1. / elemVolume;
    BlasLapackLA::matrixMatrixMultiply( worka_numFacesByNumFaces,
                                        workc_numFacesByNumFaces,
                                        transMatrix,
                                        tParam * elemVolumeInv, elemVolumeInv );
  }

};


} // namespace HybridFVMInnerProduct

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP
