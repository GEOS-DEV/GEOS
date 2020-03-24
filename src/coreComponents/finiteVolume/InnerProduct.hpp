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
 * @file InnerProduct.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_INNERPRODUCT_HPP
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_INNERPRODUCT_HPP

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;


namespace InnerProduct
{

static constexpr localIndex MAX_NUM_FACES = 15;

struct InnerProductType
{
  static constexpr integer TPFA = 0;
  static constexpr integer QUASI_TPFA = 1;
};

struct InnerProductHelper
{

  // for now, I just copy-pasted this function from TwoPointFluxApproximation
  // this will go away at some point
  static
  void makeFullTensor( R1Tensor const & values, R2SymTensor & result )
  {
    result = 0.0;
    R1Tensor axis;
    R2SymTensor temp;

    // assemble full tensor from eigen-decomposition
    for( unsigned icoord = 0; icoord < 3; ++icoord )
    {
      // assume principal axis aligned with global coordinate system
      axis = 0.0;
      axis( icoord ) = 1.0;

      // XXX: is there a more elegant way to do this?
      temp.dyadic_aa( axis );
      temp *= values( icoord );
      result += temp;
    }
  }

};

struct TPFAFaceInnerProductKernel
{
  inline static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & GEOSX_UNUSED_PARAM( nodePosition ),
           ArrayOfArraysView< localIndex const > const & GEOSX_UNUSED_PARAM( faceToNodes ),
           arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM( elemToFaces ),
           R1Tensor const & GEOSX_UNUSED_PARAM( elemCenter ),
           R1Tensor const & GEOSX_UNUSED_PARAM( elemPerm ),
           real64 const & GEOSX_UNUSED_PARAM( lengthTolerance ),
           stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & GEOSX_UNUSED_PARAM( transMatrix ) )
  {}
};


struct TPFACellInnerProductKernel
{

  inline static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           R1Tensor const & elemPerm,
           real64 const & lengthTolerance,
           stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & transMatrix )
  {
    R1Tensor faceCenter, faceNormal, faceConormal, cellToFaceVec;
    R2SymTensor permeabilityTensor;

    real64 const areaTolerance   = lengthTolerance * lengthTolerance;
    real64 const weightTolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

    localIndex const numFacesInElem = elemToFaces.size();

    // we are ready to compute the transmissibility matrix
    for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
    {
      for( localIndex jfaceLoc = 0; jfaceLoc < numFacesInElem; ++jfaceLoc )
      {
        // for now, TPFA trans
        if( ifaceLoc == jfaceLoc )
        {

          // 1) compute the face geometry data: center, normal, vector from cell center to face center
          real64 const faceArea =
            computationalGeometry::Centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                       faceToNodes.sizeOfArray( elemToFaces[ifaceLoc] ),
                                                       nodePosition,
                                                       faceCenter,
                                                       faceNormal,
                                                       areaTolerance );

          cellToFaceVec  = faceCenter;
          cellToFaceVec -= elemCenter;

          if( Dot( cellToFaceVec, faceNormal ) < 0.0 )
          {
            faceNormal *= -1;
          }

          real64 const c2fDistance = cellToFaceVec.Normalize();

          // 2) assemble full coefficient tensor from principal axis/components
          InnerProductHelper::makeFullTensor( elemPerm, permeabilityTensor );

          faceConormal.AijBj( permeabilityTensor, faceNormal );

          // 3) compute the one-sided face transmissibility
          transMatrix[ifaceLoc][jfaceLoc]  = Dot( cellToFaceVec, faceConormal );
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

};

struct QTPFACellInnerProductKernel
{

  inline static void
  Compute( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
           ArrayOfArraysView< localIndex const > const & faceToNodes,
           arraySlice1d< localIndex const > const elemToFaces,
           R1Tensor const & elemCenter,
           real64 const & elemVolume,
           R1Tensor const & elemPerm,
           real64 const & tParam,
           real64 const & lengthTolerance,
           stackArray2d< real64, MAX_NUM_FACES *MAX_NUM_FACES > const & transMatrix )
  {
    R1Tensor faceCenter, faceNormal, cellToFaceVec;
    real64 const areaTolerance = lengthTolerance * lengthTolerance;

    localIndex const numFacesInElem = elemToFaces.size();
    localIndex const dim = 3;

    // TODO: remove all the array2ds of this function once the BlasLapackLA calls have been removed
    //       work with stackArray2ds instead
    array2d< real64 > cellToFaceMat;
    array2d< real64 > normalsMat( numFacesInElem, dim );
    array2d< real64 > permMat( dim, dim );
    array2d< real64 > transMat( numFacesInElem, numFacesInElem );

    // TODO: figure out if it is possible/beneficial to preallocate these arrays
    array1d< real64 > work_dim;
    array2d< real64 > work_dimByDim;
    array2d< real64 > work_dimByNumFaces( dim, numFacesInElem );
    array2d< real64 > work_numFacesByDim( numFacesInElem, dim );
    array2d< real64 > worka_numFacesByNumFaces( numFacesInElem, numFacesInElem );
    array2d< real64 > workb_numFacesByNumFaces( numFacesInElem, numFacesInElem );
    array2d< real64 > workc_numFacesByNumFaces( numFacesInElem, numFacesInElem );

    array1d< real64 > q0;
    array1d< real64 > q1;
    array1d< real64 > q2;

    bool const orthonormalizeWithSVD = false;

    if( orthonormalizeWithSVD )
    {
      cellToFaceMat.resizeDimension< 0 >( numFacesInElem );
      cellToFaceMat.resizeDimension< 1 >( dim );
      work_dim.resize( dim );
      work_dimByDim.resize( dim, dim );
    }
    else
    {
      q0.resize( numFacesInElem );
      q1.resize( numFacesInElem );
      q2.resize( numFacesInElem );
    }

    // 1) fill the matrices cellToFaceMat and normalsMat row by row
    for( localIndex ifaceLoc = 0; ifaceLoc < numFacesInElem; ++ifaceLoc )
    {

      // compute the face geometry data: center, normal, vector from cell center to face center
      computationalGeometry::Centroid_3DPolygon( faceToNodes[elemToFaces[ifaceLoc]],
                                                 faceToNodes.sizeOfArray( elemToFaces[ifaceLoc] ),
                                                 nodePosition,
                                                 faceCenter,
                                                 faceNormal,
                                                 areaTolerance );

      cellToFaceVec  = faceCenter;
      cellToFaceVec -= elemCenter;

      if( orthonormalizeWithSVD )
      {
        cellToFaceMat( ifaceLoc, 0 ) = cellToFaceVec( 0 );
        cellToFaceMat( ifaceLoc, 1 ) = cellToFaceVec( 1 );
        cellToFaceMat( ifaceLoc, 2 ) = cellToFaceVec( 2 );
      }
      else
      {
        q0( ifaceLoc ) = cellToFaceVec( 0 );
        q1( ifaceLoc ) = cellToFaceVec( 1 );
        q2( ifaceLoc ) = cellToFaceVec( 2 );
      }

      if( Dot( cellToFaceVec, faceNormal ) < 0.0 )
      {
        faceNormal *= -1;
      }

      normalsMat( ifaceLoc, 0 ) = faceNormal( 0 );
      normalsMat( ifaceLoc, 1 ) = faceNormal( 1 );
      normalsMat( ifaceLoc, 2 ) = faceNormal( 2 );

    }

    // 2) assemble full coefficient tensor from principal axis/components
    // TODO: figure out if there is a better way to that
    R2SymTensor permeabilityTensor;
    InnerProductHelper::makeFullTensor( elemPerm, permeabilityTensor );
    for( localIndex i = 0; i < dim; ++i )
    {
      for( localIndex j = 0; j < dim; ++j )
      {
        permMat( i, j ) = permeabilityTensor( i, j );
      }
    }

    // TODO: replace the BlasLapack calls below with explicitly for loops
    //       this should be easy if MGS orthonormalization is as robust as SVD

    // 3) compute N K N'
    BlasLapackLA::matrixMatrixTMultiply( permMat,
                                         normalsMat,
                                         work_dimByNumFaces );
    BlasLapackLA::matrixMatrixMultiply( normalsMat,
                                        work_dimByNumFaces,
                                        transMat );

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
        workb_numFacesByNumFaces( i, j ) = (i == j ) ?  transMat( i, j ) : 0.0;
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
    BlasLapackLA::matrixMatrixAdd( workb_numFacesByNumFaces, transMat );
    BlasLapackLA::matrixScale( 1/elemVolume, transMat );

    // for now, I have this copy to transfer the data from the array2d to the stackArray2d
    // I need the array2d to call the BlasLapackLA functions
    // if and when everything works with explicit for loops in this kernel function,
    // I will be able to do remove all the array2ds and then I won't need this copy anymore
    for( localIndex i = 0; i < numFacesInElem; ++i )
    {
      for( localIndex j = 0; j < numFacesInElem; ++j )
      {
        transMatrix( i, j ) = transMat( i, j );
      }
    }
  }
};



} // namespace InnerProduct

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_INNERPRODUCT_HPP
