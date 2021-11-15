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
 * @file HydrofractureSolverKernels.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_

#include "HydrofractureSolverKernels.hpp"

namespace geosx
{

namespace HydrofractureSolverKernels
{

struct DeformationUpdateKernel
{

  template< typename POLICY, typename CONTACT_WRAPPER >
  static void
  launch( localIndex const size,
          CONTACT_WRAPPER const & contactWrapper,
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u,
          arrayView2d< real64 const > const & faceNormal,
          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
          arrayView2d< localIndex const > const & elemsToFaces,
          arrayView1d< real64 const > const & area,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 > const & deltaVolume,
          arrayView1d< real64 > const & aperture,
          arrayView1d< real64 > const & effectiveAperture
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
          ,
          arrayView1d< real64 const > const & apertureAtFailure,
          arrayView1d< real64 > const & separationCoeff,
          arrayView1d< real64 > const & dSeparationCoeff_dAper,
          arrayView1d< real64 const > const & separationCoeff0
#endif
          )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const kf1 = elemsToFaces[kfe][1];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      real64 temp[ 3 ] = { 0 };
      for( localIndex a=0; a<numNodesPerFace; ++a )
      {
        LvArray::tensorOps::add< 3 >( temp, u[ faceToNodeMap( kf0, a ) ] );
        LvArray::tensorOps::subtract< 3 >( temp, u[ faceToNodeMap( kf1, a ) ] );
      }

      // TODO this needs a proper contact based strategy for aperture
      aperture[kfe] = -LvArray::tensorOps::AiBi< 3 >( temp, faceNormal[ kf0 ] ) / numNodesPerFace;

      real64 dEffectiveAperture_dAperture = 0;
      effectiveAperture[kfe] = contactWrapper.computeEffectiveAperture( aperture[kfe],
                                                                        dEffectiveAperture_dAperture );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
      real64 const s = aperture[kfe] / apertureAtFailure[kfe];
      if( separationCoeff0[kfe]<1.0 && s>separationCoeff0[kfe] )
      {
        if( s >= 1.0 )
        {
          separationCoeff[kfe] = 1.0;
          dSeparationCoeff_dAper[kfe] = 0.0;
        }
        else
        {
          separationCoeff[kfe] = s;
          dSeparationCoeff_dAper[kfe] = 1.0/apertureAtFailure[kfe];
        }
      }
#endif
      deltaVolume[kfe] = effectiveAperture[kfe] * area[kfe] - volume[kfe];
    } );
  }
};

struct FluidMassResidualDerivativeAssemblyKernel
{

  template< typename CONTACT_WRAPPER >
  GEOSX_HOST_DEVICE
  static void
  computeAccumulationDerivative( CONTACT_WRAPPER const & contactWrapper,
                                 localIndex const numNodesPerFace,
                                 arraySlice1d< localIndex const > const elemsToFaces,
                                 ArrayOfArraysView< localIndex const > const faceToNodeMap,
                                 arrayView1d< globalIndex const > const dispDofNumber,
                                 real64 const (&Nbar)[ 3 ],
                                 real64 const & area,
                                 real64 const & aperture,
                                 real64 const & dens,
                                 globalIndex (& nodeDOF)[8 * 3],
                                 arraySlice1d< real64 > const dRdU )
  {
    constexpr integer kfSign[2] = { -1, 1 };
    for( localIndex kf = 0; kf < 2; ++kf )
    {
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        for( int i = 0; i < 3; ++i )
        {
          nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kf], a )] + i;
          real64 const dGap_dU = kfSign[kf] * Nbar[i] / numNodesPerFace;

          real64 dEffectiveAperture_dAperture = 0;
          real64 const effectiveAperture = contactWrapper.computeEffectiveAperture( aperture, dEffectiveAperture_dAperture );
          GEOSX_UNUSED_VAR( effectiveAperture );
          real64 const dAper_dU = dEffectiveAperture_dAperture * dGap_dU;

          dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dens * area * dAper_dU;
        }
      }
    }
  }

  template< typename CONTACT_WRAPPER >
  GEOSX_HOST_DEVICE
  static void
  computeFluxDerivative( CONTACT_WRAPPER const & contactWrapper,
                         localIndex const kfe2,
                         localIndex const numNodesPerFace,
                         arraySlice1d< localIndex const > const & columns,
                         arraySlice1d< real64 const > const & values,
                         arrayView2d< localIndex const > const elemsToFaces,
                         ArrayOfArraysView< localIndex const > const faceToNodeMap,
                         arrayView1d< globalIndex const > const dispDofNumber,
                         real64 const (&Nbar)[ 3 ],
                         arrayView1d< real64 const > const aperture,
                         globalIndex (& nodeDOF)[8 * 3],
                         arraySlice1d< real64 > const dRdU )
  {
    constexpr integer kfSign[2] = { -1, 1 };

    real64 const dRdAper = values[kfe2];
    localIndex const ei2 = columns[kfe2];

    for( localIndex kf = 0; kf < 2; ++kf )
    {
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[ei2][kf], a )] + i;
          real64 const dGap_dU = kfSign[kf] * Nbar[i] / numNodesPerFace;

          real64 dEffectiveAperture_dAperture = 0.0;
          real64 const effectiveAperture = contactWrapper.computeEffectiveAperture( aperture[ei2], dEffectiveAperture_dAperture );
          GEOSX_UNUSED_VAR( effectiveAperture );
          real64 const dAper_dU = dEffectiveAperture_dAperture * dGap_dU;

          dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dRdAper * dAper_dU;
        }
      }
    }
  }

  template< typename POLICY, typename CONTACT_WRAPPER >
  static void
  launch( localIndex const size,
          globalIndex const rankOffset,
          CONTACT_WRAPPER const & contactWrapper,
          arrayView2d< localIndex const > const elemsToFaces,
          ArrayOfArraysView< localIndex const > const faceToNodeMap,
          arrayView2d< real64 const > const faceNormal,
          arrayView1d< real64 const > const area,
          arrayView1d< real64 const > const aperture,
          arrayView1d< globalIndex const > const presDofNumber,
          arrayView1d< globalIndex const > const dispDofNumber,
          arrayView2d< real64 const > const dens,
          CRSMatrixView< real64 const, localIndex const > const dFluxResidual_dAperture,
          CRSMatrixView< real64, globalIndex const > const & localMatrix )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex ei )
    {
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[ei][0] );

      real64 Nbar[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[elemsToFaces[ei][0]] );
      LvArray::tensorOps::subtract< 3 >( Nbar, faceNormal[elemsToFaces[ei][1]] );
      LvArray::tensorOps::normalize< 3 >( Nbar );

      globalIndex const rowNumber = presDofNumber[ei] - rankOffset;
      globalIndex nodeDOF[8 * 3];
      stackArray1d< real64, 24 > dRdU( 2 * numNodesPerFace * 3 );

      computeAccumulationDerivative( contactWrapper,
                                     numNodesPerFace,
                                     elemsToFaces[ei],
                                     faceToNodeMap,
                                     dispDofNumber,
                                     Nbar,
                                     area[ei],
                                     aperture[ei],
                                     dens[ei][0],
                                     nodeDOF,
                                     dRdU );

      if( rowNumber >= 0  && rowNumber < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                          nodeDOF,
                                                                          dRdU.data(),
                                                                          2 * numNodesPerFace * 3 );
      }

      localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( ei );
      arraySlice1d< localIndex const > const & columns = dFluxResidual_dAperture.getColumns( ei );
      arraySlice1d< real64 const > const & values = dFluxResidual_dAperture.getEntries( ei );

      for( localIndex kfe2 = 0; kfe2 < numColumns; ++kfe2 )
      {
        computeFluxDerivative( contactWrapper,
                               kfe2,
                               numNodesPerFace,
                               columns,
                               values,
                               elemsToFaces,
                               faceToNodeMap,
                               dispDofNumber,
                               Nbar,
                               aperture,
                               nodeDOF,
                               dRdU );

        if( rowNumber >= 0 && rowNumber < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                            nodeDOF,
                                                                            dRdU.data(),
                                                                            2 * numNodesPerFace * 3 );
        }
      }
    } );
  }
};

} /* namespace HydrofractureSolverKernels */

} /* namespace geosx */

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_
