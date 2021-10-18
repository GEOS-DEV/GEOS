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
 * @file SinglePhasePoromechanicsFluxKernels.cpp
 */

#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"
#include "SinglePhasePoromechanicsFluxKernels.hpp"

namespace geosx
{

namespace SinglePhasePoromechanicsFluxKernels
{

using namespace FluxKernelsHelper;

/******************************** EmbeddedSurfaceFluxKernel ********************************/

template<>
void EmbeddedSurfaceFluxKernel::
  launch< CellElementStencilTPFAWrapper >( CellElementStencilTPFAWrapper const & stencilWrapper,
                                           real64 const dt,
                                           globalIndex const rankOffset,
                                           ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                           ElementViewConst< arrayView1d< globalIndex const > > const & jumpDofNumber,
                                           ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                           ElementViewConst< arrayView1d< real64 const > > const & pres,
                                           ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                           ElementViewConst< arrayView2d< real64 const > > const & dens,
                                           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                           ElementViewConst< arrayView1d< real64 const > > const & mob,
                                           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_UNUSED_VAR( jumpDofNumber );
  GEOSX_UNUSED_VAR( dPerm_dAper );

  SinglePhaseFVMKernels::FluxKernel::launch( stencilWrapper,
                                             dt,
                                             rankOffset,
                                             pressureDofNumber,
                                             ghostRank,
                                             pres,
                                             dPres,
                                             gravCoef,
                                             dens,
                                             dDens_dPres,
                                             mob,
                                             dMob_dPres,
                                             permeability,
                                             dPerm_dPres,
                                             localMatrix,
                                             localRhs );
}

template<>
void EmbeddedSurfaceFluxKernel::
  launch< EmbeddedSurfaceToCellStencilWrapper >( EmbeddedSurfaceToCellStencilWrapper const & stencilWrapper,
                                                 real64 const dt,
                                                 globalIndex const rankOffset,
                                                 ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                                 ElementViewConst< arrayView1d< globalIndex const > > const & jumpDofNumber,
                                                 ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                                 ElementViewConst< arrayView1d< real64 const > > const & pres,
                                                 ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                                 ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                                 ElementViewConst< arrayView2d< real64 const > > const & dens,
                                                 ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                                 ElementViewConst< arrayView1d< real64 const > > const & mob,
                                                 ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                                 ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                                 ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                                 ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_UNUSED_VAR( jumpDofNumber );
  GEOSX_UNUSED_VAR( dPerm_dAper );

  SinglePhaseFVMKernels::FluxKernel::launch( stencilWrapper,
                                             dt,
                                             rankOffset,
                                             pressureDofNumber,
                                             ghostRank,
                                             pres,
                                             dPres,
                                             gravCoef,
                                             dens,
                                             dDens_dPres,
                                             mob,
                                             dMob_dPres,
                                             permeability,
                                             dPerm_dPres,
                                             localMatrix,
                                             localRhs );
}


template<>
void EmbeddedSurfaceFluxKernel::
  launch< SurfaceElementStencilWrapper >( SurfaceElementStencilWrapper const & stencilWrapper,
                                          real64 const dt,
                                          globalIndex const rankOffset,
                                          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                          ElementViewConst< arrayView1d< globalIndex const > > const & jumpDofNumber,
                                          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                          ElementViewConst< arrayView1d< real64 const > > const & pres,
                                          ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                          ElementViewConst< arrayView2d< real64 const > > const & dens,
                                          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                          ElementViewConst< arrayView1d< real64 const > > const & mob,
                                          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                          ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  constexpr localIndex MAX_NUM_FLUX_ELEMS = SurfaceElementStencilWrapper::NUM_POINT_IN_FLUX;
  constexpr localIndex MAX_STENCIL_SIZE = SurfaceElementStencilWrapper::MAX_STENCIL_SIZE;

  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();


  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );
    localIndex const numDofs = stencilSize * 2;// pressure and normal jump

    // working arrays
    stackArray1d< globalIndex, MAX_NUM_FLUX_ELEMS > dofColIndices( numDofs );
    stackArray1d< real64, MAX_NUM_FLUX_ELEMS > localFlux( numFluxElems );
    stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_STENCIL_SIZE > localFluxJacobian( numFluxElems, numDofs );

    // compute transmissibility
    real64 transmissibility[SurfaceElementStencilWrapper::MAX_NUM_OF_CONNECTIONS][2];
    real64 dTrans_dPres[SurfaceElementStencilWrapper::MAX_NUM_OF_CONNECTIONS][2];
    real64 dTrans_dAper[SurfaceElementStencilWrapper::MAX_NUM_OF_CONNECTIONS][2];

    stencilWrapper.computeWeights( iconn,
                                   permeability,
                                   dPerm_dPres,
                                   dPerm_dAper,
                                   transmissibility,
                                   dTrans_dPres,
                                   dTrans_dAper );

    compute( stencilSize,
             seri[iconn],
             sesri[iconn],
             sei[iconn],
             transmissibility,
             dTrans_dPres,
             dTrans_dAper,
             pres,
             dPres,
             gravCoef,
             dens,
             dDens_dPres,
             mob,
             dMob_dPres,
             dt,
             localFlux,
             localFluxJacobian );

    // extract DOF numbers
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex localDofIndex = 2 * i;
      dofColIndices[ localDofIndex ]     = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
      dofColIndices[ localDofIndex + 1 ] = jumpDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
    }

    for( localIndex i = 0; i < numFluxElems; ++i )
    {

      if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                          dofColIndices.data(),
                                                                          localFluxJacobian[i].dataIfContiguous(),
                                                                          numDofs );

      }
    }
  } );
}

template< localIndex MAX_NUM_OF_CONNECTIONS >
void EmbeddedSurfaceFluxKernel::
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[MAX_NUM_OF_CONNECTIONS][2],
           real64 const (&dTrans_dPres)[MAX_NUM_OF_CONNECTIONS][2],
           real64 const (&dTrans_dAper)[MAX_NUM_OF_CONNECTIONS][2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & dens,
           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
           ElementViewConst< arrayView1d< real64 const > > const & mob,
           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
{
  GEOSX_UNUSED_VAR( numFluxElems );

  real64 fluxVal = 0.0;
  real64 dFlux_dTrans = 0.0;
  real64 trans[2] = {transmissibility[0][0], transmissibility[0][1]};
  real64 dTrans[2] = { dTrans_dPres[0][0], dTrans_dPres[0][1] };
  real64 dFlux_dP[2] = {0.0, 0.0};
  localIndex const regionIndex[2]    = {seri[0], seri[1]};
  localIndex const subRegionIndex[2] = {sesri[0], sesri[1]};
  localIndex const elementIndex[2]   = {sei[0], sei[1]};


  computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                          trans,
                          dTrans,
                          pres,
                          dPres,
                          gravCoef,
                          dens,
                          dDens_dPres,
                          mob,
                          dMob_dPres,
                          fluxVal,
                          dFlux_dP,
                          dFlux_dTrans );



  // populate local flux vector and derivatives
  flux[0] =  dt * fluxVal;
  flux[1] = -dt * fluxVal;

  real64 dFlux_dAper[2] = {0.0, 0.0};
  dFlux_dAper[0] = dt * dFlux_dTrans * dTrans_dAper[0][0];
  dFlux_dAper[1] = -dt * dFlux_dTrans * dTrans_dAper[0][1];

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    localIndex const dofIndex = 2*ke;

    fluxJacobian[0][dofIndex]   =  dt * dFlux_dP[ke];
    fluxJacobian[0][dofIndex+1] =  dt * dFlux_dAper[ke];
    fluxJacobian[1][dofIndex]   = -dt * dFlux_dP[ke];
    fluxJacobian[1][dofIndex+1] = -dt * dFlux_dAper[ke];
  }
}

/******************************** FaceElementFluxKernel ********************************/

template<>
void FaceElementFluxKernel::
  launch< CellElementStencilTPFAWrapper >( CellElementStencilTPFAWrapper const & stencilWrapper,
                                           real64 const dt,
                                           globalIndex const rankOffset,
                                           ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                           ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                           ElementViewConst< arrayView1d< real64 const > > const & pres,
                                           ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                           ElementViewConst< arrayView2d< real64 const > > const & dens,
                                           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                           ElementViewConst< arrayView1d< real64 const > > const & mob,
                                           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs,
                                           CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOSX_UNUSED_VAR( dPerm_dAper );
  GEOSX_UNUSED_VAR( dR_dAper );

  SinglePhaseFVMKernels::FluxKernel::launch( stencilWrapper,
                                             dt,
                                             rankOffset,
                                             pressureDofNumber,
                                             ghostRank,
                                             pres,
                                             dPres,
                                             gravCoef,
                                             dens,
                                             dDens_dPres,
                                             mob,
                                             dMob_dPres,
                                             permeability,
                                             dPerm_dPres,
                                             localMatrix,
                                             localRhs );
}

template<>
void FaceElementFluxKernel::
  launch< FaceElementToCellStencilWrapper >( FaceElementToCellStencilWrapper const & stencilWrapper,
                                             real64 const dt,
                                             globalIndex const rankOffset,
                                             ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                             ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                             ElementViewConst< arrayView1d< real64 const > > const & pres,
                                             ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                             ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                             ElementViewConst< arrayView2d< real64 const > > const & dens,
                                             ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                             ElementViewConst< arrayView1d< real64 const > > const & mob,
                                             ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                             ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                             ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                             ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs,
                                             CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOSX_UNUSED_VAR( dPerm_dAper );
  GEOSX_UNUSED_VAR( dR_dAper );

  SinglePhaseFVMKernels::FluxKernel::launch( stencilWrapper,
                                             dt,
                                             rankOffset,
                                             pressureDofNumber,
                                             ghostRank,
                                             pres,
                                             dPres,
                                             gravCoef,
                                             dens,
                                             dDens_dPres,
                                             mob,
                                             dMob_dPres,
                                             permeability,
                                             dPerm_dPres,
                                             localMatrix,
                                             localRhs );
}

template<>
void FaceElementFluxKernel::
  launch< SurfaceElementStencilWrapper >( SurfaceElementStencilWrapper const & stencilWrapper,
                                          real64 const dt,
                                          globalIndex const rankOffset,
                                          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                          ElementViewConst< arrayView1d< real64 const > > const & pres,
                                          ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                          ElementViewConst< arrayView2d< real64 const > > const & dens,
                                          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                          ElementViewConst< arrayView1d< real64 const > > const & mob,
                                          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                          ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs,
                                          CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  constexpr localIndex MAX_NUM_FLUX_ELEMS = SurfaceElementStencilWrapper::NUM_POINT_IN_FLUX;
  constexpr localIndex MAX_STENCIL_SIZE  = SurfaceElementStencilWrapper::MAX_STENCIL_SIZE;

  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );
    localIndex const numDofs = stencilSize; // pressures

    // For now, we have to filter out connections for which numElems == 1 in this function and not early on in
    // TwoPointFluxApproximation.cpp.
    // The reason for keeping the connections numElems == 1 is that the ProppantTransport solver needs these connections to produce correct
    // results.
    if( numFluxElems > 1 )
    {
      // working arrays
      stackArray1d< globalIndex, MAX_NUM_FLUX_ELEMS > dofColIndices( numDofs );
      stackArray1d< localIndex, MAX_NUM_FLUX_ELEMS > localColIndices( numFluxElems );

      stackArray1d< real64, MAX_NUM_FLUX_ELEMS > localFlux( numFluxElems );
      stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_STENCIL_SIZE > localFluxJacobian( numFluxElems, numDofs );

      // need to store this for later use in determining the dFlux_dU terms when using better permeabilty approximations.
      stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_STENCIL_SIZE > dFlux_dAper( numFluxElems, stencilSize );

      // compute transmissibility
      real64 transmissibility[SurfaceElementStencilWrapper::MAX_NUM_OF_CONNECTIONS][2];
      real64 dTrans_dPres[SurfaceElementStencilWrapper::MAX_NUM_OF_CONNECTIONS][2];
      real64 dTrans_dAper[SurfaceElementStencilWrapper::MAX_NUM_OF_CONNECTIONS][2];

      stencilWrapper.computeWeights( iconn,
                                     permeability,
                                     dPerm_dPres,
                                     dPerm_dAper,
                                     transmissibility,
                                     dTrans_dPres,
                                     dTrans_dAper );

      compute( stencilSize,
               seri[iconn],
               sesri[iconn],
               sei[iconn],
               transmissibility,
               dTrans_dPres,
               dTrans_dAper,
               pres,
               dPres,
               gravCoef,
               dens,
               dDens_dPres,
               mob,
               dMob_dPres,
               dt,
               localFlux,
               localFluxJacobian,
               dFlux_dAper );

      // extract DOF numbers
      for( localIndex i = 0; i < numDofs; ++i )
      {
        dofColIndices[i] = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localColIndices[i] = sei( iconn, i );
      }

      for( localIndex i = 0; i < numFluxElems; ++i )
      {
        if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
        {
          globalIndex const globalRow = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            dofColIndices.data(),
                                                                            localFluxJacobian[i].dataIfContiguous(),
                                                                            stencilSize );

          dR_dAper.addToRowBinarySearch< parallelDeviceAtomic >( sei( iconn, i ),
                                                                 localColIndices.data(),
                                                                 dFlux_dAper[i].dataIfContiguous(),
                                                                 stencilSize );
        }
      }
    }
  } );

}

void FaceElementFluxKernel::
  launch( SurfaceElementStencilWrapper const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
          ElementViewConst< arrayView3d< real64 const > > const & permeabilityMultiplier,
          R1Tensor const & gravityVector,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  constexpr localIndex MAX_NUM_FLUX_ELEMS = SurfaceElementStencilWrapper::NUM_POINT_IN_FLUX;
  constexpr localIndex MAX_STENCIL_SIZE  = SurfaceElementStencilWrapper::MAX_STENCIL_SIZE;
  constexpr localIndex MAX_NUM_OF_CONNECTIONS  = SurfaceElementStencilWrapper::MAX_NUM_OF_CONNECTIONS;

  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );
    localIndex const numDofs = stencilSize; // pressures

    // For now, we have to filter out connections for which numElems == 1 in this function and not early on in
    // TwoPointFluxApproximation.cpp.
    // The reason for keeping the connections numElems == 1 is that the ProppantTransport solver needs these connections to produce correct
    // results.
    if( numFluxElems > 1 )
    {

      // working arrays
      stackArray1d< globalIndex, MAX_NUM_FLUX_ELEMS > dofColIndices( numDofs );
      stackArray1d< localIndex, MAX_NUM_FLUX_ELEMS > localColIndices( numFluxElems );

      stackArray1d< real64, MAX_NUM_FLUX_ELEMS > localFlux( numFluxElems );
      stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_STENCIL_SIZE > localFluxJacobian( numFluxElems, numDofs );

      // need to store this for later use in determining the dFlux_dU terms when using better permeabilty approximations.
      stackArray2d< real64, MAX_NUM_FLUX_ELEMS * MAX_STENCIL_SIZE > dFlux_dAper( numFluxElems, stencilSize );

      // compute transmissibility
      real64 transmissibility[MAX_NUM_OF_CONNECTIONS][2], dTrans_dPres[MAX_NUM_OF_CONNECTIONS][2], dTrans_dAper[MAX_NUM_OF_CONNECTIONS][2];
      GEOSX_UNUSED_VAR( dPerm_dPres, dPerm_dAper );
      stencilWrapper.computeWeights( iconn,
                                     permeability,
                                     permeabilityMultiplier,
                                     gravityVector,
                                     transmissibility );

      compute( stencilSize,
               seri[iconn],
               sesri[iconn],
               sei[iconn],
               transmissibility,
               dTrans_dPres,
               dTrans_dAper,
               pres,
               dPres,
               gravCoef,
               dens,
               dDens_dPres,
               mob,
               dMob_dPres,
               dt,
               localFlux,
               localFluxJacobian,
               dFlux_dAper );

      // extract DOF numbers
      for( localIndex i = 0; i < numDofs; ++i )
      {
        dofColIndices[i] = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localColIndices[i] = sei( iconn, i );
      }

      for( localIndex i = 0; i < numFluxElems; ++i )
      {
        if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
        {
          globalIndex const globalRow = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            dofColIndices.data(),
                                                                            localFluxJacobian[i].dataIfContiguous(),
                                                                            stencilSize );

        }
      }
    }
  } );

}



template< localIndex MAX_NUM_OF_CONNECTIONS >
GEOSX_HOST_DEVICE
void
FaceElementFluxKernel::compute( localIndex const numFluxElems,
                                arraySlice1d< localIndex const > const & seri,
                                arraySlice1d< localIndex const > const & sesri,
                                arraySlice1d< localIndex const > const & sei,
                                real64 const (&transmissibility)[MAX_NUM_OF_CONNECTIONS][2],
                                real64 const (&dTrans_dPres)[MAX_NUM_OF_CONNECTIONS][2],
                                real64 const (&dTrans_dAper)[MAX_NUM_OF_CONNECTIONS][2],
                                ElementViewConst< arrayView1d< real64 const > > const & pres,
                                ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView2d< real64 const > > const & dens,
                                ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                ElementViewConst< arrayView1d< real64 const > > const & mob,
                                ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                real64 const dt,
                                arraySlice1d< real64 > const & flux,
                                arraySlice2d< real64 > const & fluxJacobian,
                                arraySlice2d< real64 > const & dFlux_dAperture )
{

  localIndex k[2];
  localIndex connectionIndex = 0;
  for( k[0]=0; k[0]<numFluxElems; ++k[0] )
  {
    for( k[1]=k[0]+1; k[1]<numFluxElems; ++k[1] )
    {
      real64 fluxVal = 0.0;
      real64 dFlux_dTrans = 0.0;
      real64 trans[2] = {transmissibility[connectionIndex][0], transmissibility[connectionIndex][1]};
      real64 dTrans[2] = { dTrans_dPres[connectionIndex][0], dTrans_dPres[connectionIndex][1] };
      real64 dFlux_dP[2] = {0.0, 0.0};
      localIndex const regionIndex[2]    = {seri[k[0]], seri[k[1]]};
      localIndex const subRegionIndex[2] = {sesri[k[0]], sesri[k[1]]};
      localIndex const elementIndex[2]   = {sei[k[0]], sei[k[1]]};

      computeSinglePhaseFlux( regionIndex, subRegionIndex, elementIndex,
                              trans,
                              dTrans,
                              pres,
                              dPres,
                              gravCoef,
                              dens,
                              dDens_dPres,
                              mob,
                              dMob_dPres,
                              fluxVal,
                              dFlux_dP,
                              dFlux_dTrans );

      // populate local flux vector and derivatives
      flux[k[0]] += dt * fluxVal;
      flux[k[1]] -= dt * fluxVal;

      real64 dFlux_dAper[2] = {0.0, 0.0};
      dFlux_dAper[0] = dt * dFlux_dTrans * dTrans_dAper[connectionIndex][0];
      dFlux_dAper[1] = -dt * dFlux_dTrans * dTrans_dAper[connectionIndex][1];

      fluxJacobian[k[0]][k[0]] += dFlux_dP[0] * dt;
      fluxJacobian[k[0]][k[1]] += dFlux_dP[1] * dt;
      fluxJacobian[k[1]][k[0]] -= dFlux_dP[0] * dt;
      fluxJacobian[k[1]][k[1]] -= dFlux_dP[1] * dt;

      dFlux_dAperture[k[0]][k[0]] += dFlux_dAper[0];
      dFlux_dAperture[k[0]][k[1]] += dFlux_dAper[1];
      dFlux_dAperture[k[1]][k[0]] -= dFlux_dAper[0];
      dFlux_dAperture[k[1]][k[1]] -= dFlux_dAper[1];

      connectionIndex++;
    }
  }
}


}// namespace SinglePhaseFVMKernels

} // namespace geosx
