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
 * @file HybridFVMInnerProduct.cpp
 *
 */

#include "HybridFVMInnerProduct.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

namespace HybridFVMInnerProduct
{

// TODO This should be removed.
void
HybridFVMInnerProductHelper::makeFullTensor( R1Tensor const & values,
                                             stackArray2d< real64, 9 > & result )
{
  result = 0.0;
  result[ 0 ][ 0 ] = values[ 0 ];
  result[ 1 ][ 1 ] = values[ 1 ];
  result[ 2 ][ 2 ] = values[ 2 ];
}

void
TPFACellInnerProductKernel::Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                     ArrayOfArraysView< localIndex const > const & faceToNodes,
                                     arraySlice1d< localIndex const > const elemToFaces,
                                     R1Tensor const & elemCenter,
                                     R1Tensor const & elemPerm,
                                     real64 const & lengthTolerance,
                                     stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & transMatrix )
{
  localIndex const numFacesInElem = elemToFaces.size();

  real64 const areaTolerance = lengthTolerance * lengthTolerance;
  real64 const weightTolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  // we are ready to compute the transmissibility matrix
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
    {
      // for now, TPFA trans
      if( ifaceLoc == jfaceLoc )
      {
        real64 faceCenter[ 3 ], faceNormal[ 3 ], faceConormal[ 3 ], cellToFaceVec[ 3 ];

        // 1) compute the face geometry data: center, normal, vector from cell center to face center
        real64 const faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[ elemToFaces[ ifaceLoc ] ],
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

        // TODO: take symmetry into account to optimize this
        // TODO: Change to LvArray::tensorOps::AiBi
        faceConormal[ 0 ] = elemPerm[ 0 ] * faceNormal[ 0 ];
        faceConormal[ 1 ] = elemPerm[ 1 ] * faceNormal[ 1 ];
        faceConormal[ 2 ] = elemPerm[ 2 ] * faceNormal[ 2 ];

        // 3) compute the one-sided face transmissibility
        transMatrix[ifaceLoc][jfaceLoc]  = LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal );
        transMatrix[ifaceLoc][jfaceLoc] *= faceArea / c2fDistance;
        transMatrix[ifaceLoc][jfaceLoc]  = std::max( transMatrix[ifaceLoc][jfaceLoc],
                                                     weightTolerance );
      }
      else
      {
        transMatrix[ifaceLoc][jfaceLoc] = 0;
      }
    }
  }
}


void
QTPFACellInnerProductKernel::Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                      ArrayOfArraysView< localIndex const > const & faceToNodes,
                                      arraySlice1d< localIndex const > const elemToFaces,
                                      R1Tensor const & elemCenter,
                                      real64 const & elemVolume,
                                      R1Tensor const & elemPerm,
                                      real64 const & tParam,
                                      real64 const & lengthTolerance,
                                      bool const & orthonormalizeWithSVD,
                                      stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & transMatrix )
{
  localIndex constexpr dim = 3;
  localIndex const numFacesInElem = elemToFaces.size();

  real64 const areaTolerance = lengthTolerance * lengthTolerance;

  stackArray2d< real64, dim *MAX_NUM_FACES > cellToFaceMat( numFacesInElem, dim );
  stackArray2d< real64, dim *MAX_NUM_FACES > normalsMat( numFacesInElem, dim );
  stackArray2d< real64, dim *dim > permMat( dim, dim );

  stackArray1d< real64, dim > work_dim( dim );
  stackArray2d< real64, dim *dim > work_dimByDim ( dim, dim );
  stackArray2d< real64, dim *MAX_NUM_FACES > work_dimByNumFaces( dim, numFacesInElem );
  stackArray2d< real64, dim *MAX_NUM_FACES > work_numFacesByDim( numFacesInElem, dim );
  stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > worka_numFacesByNumFaces( numFacesInElem, numFacesInElem );
  stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > workb_numFacesByNumFaces( numFacesInElem, numFacesInElem );
  stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > workc_numFacesByNumFaces( numFacesInElem, numFacesInElem );

  stackArray1d< real64, MAX_NUM_FACES > q0( numFacesInElem );
  stackArray1d< real64, MAX_NUM_FACES > q1( numFacesInElem );
  stackArray1d< real64, MAX_NUM_FACES > q2( numFacesInElem );

  // 1) fill the matrices cellToFaceMat and normalsMat row by row
  for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
  {
    real64 faceCenter[ 3 ], faceNormal[ 3 ], cellToFaceVec[ 3 ];

    // compute the face geometry data: center, normal, vector from cell center to face center
    real64 const faceArea = computationalGeometry::Centroid_3DPolygon( faceToNodes[ elemToFaces[ ifaceLoc ] ],
                                                                       nodePosition,
                                                                       faceCenter,
                                                                       faceNormal,
                                                                       areaTolerance );

    LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
    // TODO: Change to LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter );
    cellToFaceVec[ 0 ] -= elemCenter[ 0 ];
    cellToFaceVec[ 1 ] -= elemCenter[ 1 ];
    cellToFaceVec[ 2 ] -= elemCenter[ 2 ];

    if( orthonormalizeWithSVD )
    {
      LvArray::tensorOps::copy< 3 >( cellToFaceMat[ ifaceLoc ], cellToFaceVec );
    }
    else
    {
      q0( ifaceLoc ) = cellToFaceVec[ 0 ];
      q1( ifaceLoc ) = cellToFaceVec[ 1 ];
      q2( ifaceLoc ) = cellToFaceVec[ 2 ];
    }

    if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    LvArray::tensorOps::scaledCopy< 3 >( normalsMat[ ifaceLoc ], faceNormal, faceArea );
  }

  // 2) assemble full coefficient tensor from principal axis/components
  HybridFVMInnerProductHelper::makeFullTensor( elemPerm, permMat );

  // TODO: replace the BlasLapack calls below with explicitly for loops
  //       this should be easy if MGS orthonormalization is as robust as SVD
  //       Also don't construct permMat.

  // 3) compute N K N'
  BlasLapackLA::matrixMatrixTMultiply( permMat,
                                       normalsMat,
                                       work_dimByNumFaces );
  BlasLapackLA::matrixMatrixMultiply( normalsMat,
                                      work_dimByNumFaces,
                                      transMatrix );

  // 4) compute the orthonormalization of the matrix cellToFaceVec
  //    This is done either with SVD or MGS
  //    If we find that MGS is stable, I will remove SVD

  if( orthonormalizeWithSVD )
  {
    // calling SVD seems to be an overkill to orthonormalize the 3 columns of cellToFaceMat...
    BlasLapackLA::matrixSVD( cellToFaceMat,
                             work_numFacesByDim,
                             work_dim,
                             work_dimByDim );
  }
  else
  {
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

    // TODO: remove the copies once the BlasLapackLA calls have been removed
    for( int i = 0; i < numFacesInElem; ++i )
    {
      work_numFacesByDim( i, 0 ) = q0( i );
      work_numFacesByDim( i, 1 ) = q1( i );
      work_numFacesByDim( i, 2 ) = q2( i );
    }
  }

  // 5) compute P_Q = I - QQ'
  BlasLapackLA::matrixMatrixTMultiply( work_numFacesByDim,
                                       work_numFacesByDim,
                                       worka_numFacesByNumFaces );
  BlasLapackLA::matrixScale( -1, worka_numFacesByNumFaces );
  for( localIndex i = 0; i < numFacesInElem; ++i )
  {
    worka_numFacesByNumFaces( i, i )++;
  }

  // 6) compute P_Q D P_Q where D = diag(diag(N K N')
  for( localIndex i = 0; i < numFacesInElem; ++i )
  {
    for( localIndex j = 0; j < numFacesInElem; ++j )
    {
      workb_numFacesByNumFaces( i, j ) = (i == j ) ?  transMatrix( i, j ) : 0.0;
    }
  }
  BlasLapackLA::matrixMatrixMultiply( workb_numFacesByNumFaces,
                                      worka_numFacesByNumFaces,
                                      workc_numFacesByNumFaces );
  BlasLapackLA::matrixMatrixMultiply( worka_numFacesByNumFaces,
                                      workc_numFacesByNumFaces,
                                      workb_numFacesByNumFaces );

  // 7) compute T = ( N K N' + t U diag(diag(N K N')) U ) / elemVolume
  BlasLapackLA::matrixScale( tParam, workb_numFacesByNumFaces );
  BlasLapackLA::matrixMatrixAdd( workb_numFacesByNumFaces, transMatrix );
  BlasLapackLA::matrixScale( 1/elemVolume, transMatrix );

}

}

}
