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

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_

#include "HydrofractureSolverKernels.hpp"

namespace geos
{

namespace hydrofractureSolverKernels
{

struct DeformationUpdateKernel
{

  template< typename POLICY, typename CONTACT_WRAPPER >
  static void
  launch( localIndex const size,
          CONTACT_WRAPPER const & contactWrapper,
          real64 const contactStiffness, 
          arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u,
          arrayView2d< real64 const > const & faceNormal,
          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
          ArrayOfArraysView< localIndex const > const & elemsToFaces,
          arrayView1d< real64 const > const & area,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 > const & deltaVolume,
          arrayView1d< real64 > const & aperture,
          arrayView1d< real64 > const & hydraulicAperture
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
          ,
          arrayView1d< real64 const > const & apertureAtFailure,
          arrayView1d< real64 > const & separationCoeff,
          arrayView1d< real64 > const & dSeparationCoeff_dAper,
          arrayView1d< real64 const > const & separationCoeff0
#endif
          )
  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      if( elemsToFaces.sizeOfArray( kfe ) != 2 )
      { return; }

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
      real64 const mechanicalAperture = -LvArray::tensorOps::AiBi< 3 >( temp, faceNormal[ kf0 ] ) / numNodesPerFace;
      aperture[kfe] = mechanicalAperture; 

      real64 const initialAperture = 1e-5;
      real64 const refNormalStress = 5e7;  

      // real64 dHydraulicAperture_dAperture = 0.0;
      // hydraulicAperture[kfe] = contactWrapper.computeHydraulicAperture( aperture[kfe], dHydraulicAperture_dAperture );

      GEOS_UNUSED_VAR( contactWrapper );

      real64 const penaltyNormalStress = - contactStiffness * mechanicalAperture; 

      hydraulicAperture[kfe] = (mechanicalAperture >= 0.0)? (mechanicalAperture + initialAperture) : initialAperture / ( 1 + 9*penaltyNormalStress/refNormalStress );
      real64 const dHydraulicAperture_dNormalStress = - hydraulicAperture[kfe] / ( 1 + 9*penaltyNormalStress/refNormalStress ) * 9/refNormalStress;

      real64 dHydraulicAperture_dAperture = (mechanicalAperture >= 0.0)? 1.0:dHydraulicAperture_dNormalStress * -contactStiffness;

      GEOS_UNUSED_VAR( dHydraulicAperture_dAperture );

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
      deltaVolume[kfe] = hydraulicAperture[kfe] * area[kfe] - volume[kfe];
    } );
  }
};

struct FluidMassResidualDerivativeAssemblyKernel
{
  template< typename CONTACT_WRAPPER >
  GEOS_HOST_DEVICE
  inline
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

          real64 dHydraulicAperture_dAperture = 0;
          real64 const hydraulicAperture = contactWrapper.computeHydraulicAperture( aperture, dHydraulicAperture_dAperture );
          GEOS_UNUSED_VAR( hydraulicAperture );
          real64 const dAper_dU = dHydraulicAperture_dAperture * dGap_dU;

          dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dens * area * dAper_dU;
        }
      }
    }
  }

  template< typename CONTACT_WRAPPER >
  GEOS_HOST_DEVICE
  inline
  static void
  computeFluxDerivative( CONTACT_WRAPPER const & contactWrapper,
                         localIndex const kfe2,
                         localIndex const numNodesPerFace,
                         arraySlice1d< localIndex const > const & columns,
                         arraySlice1d< real64 const > const & values,
                         ArrayOfArraysView< localIndex const > const elemsToFaces,
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

          real64 dHydraulicAperture_dAperture = 0.0;
          real64 const hydraulicAperture = contactWrapper.computeHydraulicAperture( aperture[ei2], dHydraulicAperture_dAperture );
          GEOS_UNUSED_VAR( hydraulicAperture );
          real64 const dAper_dU = dHydraulicAperture_dAperture * dGap_dU;

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
          ArrayOfArraysView< localIndex const > const elemsToFaces,
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
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex ei )
    {
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[ei][0] );

      real64 Nbar[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[elemsToFaces[ei][0]] );
      LvArray::tensorOps::subtract< 3 >( Nbar, faceNormal[elemsToFaces[ei][1]] );
      LvArray::tensorOps::normalize< 3 >( Nbar );

      globalIndex const rowNumber = presDofNumber[ei] - rankOffset;
      globalIndex nodeDOF[8 * 3];
      stackArray1d< real64, 24 > dRdU( 2 * numNodesPerFace * 3 );
//
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
//
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

} /* namespace hydrofractureSolverKernels */

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVERKERNELS_HPP_
