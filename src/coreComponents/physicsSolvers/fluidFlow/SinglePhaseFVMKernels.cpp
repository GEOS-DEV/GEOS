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
 * @file SinglePhaseFVMKernels.cpp
 */

#include "SinglePhaseFVMKernels.hpp"

namespace geosx
{

namespace SinglePhaseFVMKernels
{


inline void addLocalContributionsToGlobalSystem( localIndex const numFluxElems,
                                                 localIndex const stencilSize,
                                                 globalIndex const * const eqnRowIndices,
                                                 globalIndex const * const dofColIndices,
                                                 real64 const * const localFluxJacobian,
                                                 real64 const * const localFlux,
                                                 ParallelMatrix * const jacobian,
                                                 ParallelVector * const residual )
{

  // Add to global residual/jacobian
  jacobian->add( eqnRowIndices,
                 dofColIndices,
                 localFluxJacobian,
                 numFluxElems,
                 stencilSize );

  residual->add( eqnRowIndices,
                 localFlux,
                 numFluxElems );

}

template< typename AtomicPolicy >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE
void addLocalContributionsToGlobalSystem( localIndex const numFluxElems,
                                          localIndex const stencilSize,
                                          globalIndex const * const eqnRowIndices,
                                          globalIndex const * const dofColIndices,
                                          globalIndex const rankOffset,
                                          real64 const * const localFluxJacobian,
                                          real64 const * const localFlux,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{

  for( localIndex i = 0; i < numFluxElems; ++i )
  {
    localIndex const localRow = eqnRowIndices[i] - rankOffset;
    if( localRow >= 0 && localRow < localMatrix.numRows() )
    {
      localMatrix.addToRowBinarySearchUnsorted< AtomicPolicy >( localRow,
                                                                dofColIndices,
                                                                localFluxJacobian + i * stencilSize,
                                                                stencilSize );
      atomicAdd( AtomicPolicy{}, &localRhs[localRow], localFlux[i] );
    }
  }
}

template<>
void FluxKernel::
  Launch< CellElementStencilTPFA >( CellElementStencilTPFA const & stencil,
                                    real64 const dt,
                                    globalIndex const rankOffset,
                                    ElementView< arrayView1d< globalIndex > > const & dofNumber,
                                    ElementView< arrayView1d< real64 const > > const & pres,
                                    ElementView< arrayView1d< real64 const > > const & dPres,
                                    ElementView< arrayView1d< real64 const > > const & gravCoef,
                                    ElementView< arrayView2d< real64 const > > const & dens,
                                    ElementView< arrayView2d< real64 const > > const & dDens_dPres,
                                    ElementView< arrayView1d< real64 const > > const & mob,
                                    ElementView< arrayView1d< real64 const > > const & dMob_dPres,
                                    ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture0 ),
                                    ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                    ElementView< arrayView1d< R1Tensor const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
                                    R1Tensor const,
                                    real64 const,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                                    ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( s ),
                                    ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSdAper ),
#endif
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs,
                                    CRSMatrixView< real64, localIndex const > const & GEOSX_UNUSED_PARAM( dR_dAper ) )
{
  constexpr localIndex maxNumFluxElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex numFluxElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = CellElementStencilTPFA::MAX_STENCIL_SIZE;
  constexpr localIndex stencilSize  = CellElementStencilTPFA::MAX_STENCIL_SIZE;

  typename CellElementStencilTPFA::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename CellElementStencilTPFA::WeightContainerViewConstType const & weights = stencil.getWeights();

  forAll< parallelDevicePolicy< 128 > >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    // working arrays
    stackArray1d< globalIndex, numFluxElems > eqnRowIndices( numFluxElems );
    stackArray1d< globalIndex, maxNumFluxElems > dofColIndices( stencilSize );

    stackArray1d< real64, maxNumFluxElems > localFlux( numFluxElems );
    stackArray2d< real64, maxNumFluxElems *maxStencilSize > localFluxJacobian( numFluxElems, stencilSize );

    FluxKernel::Compute( stencilSize,
                         seri[iconn],
                         sesri[iconn],
                         sei[iconn],
                         weights[iconn],
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
    for( localIndex i = 0; i < numFluxElems; ++i )
    {
      eqnRowIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
    }

    for( localIndex i = 0; i < stencilSize; ++i )
    {
      dofColIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
    }

    addLocalContributionsToGlobalSystem< parallelDeviceAtomic >( numFluxElems,
                                                                 stencilSize,
                                                                 eqnRowIndices.data(),
                                                                 dofColIndices.data(),
                                                                 rankOffset,
                                                                 localFluxJacobian.data(),
                                                                 localFlux.data(),
                                                                 localMatrix,
                                                                 localRhs );
  } );
}

template<>
void FluxKernel::
  Launch< FaceElementStencil >( FaceElementStencil const & stencil,
                                real64 const dt,
                                globalIndex const rankOffset,
                                FluxKernel::ElementView< arrayView1d< globalIndex const > > const & dofNumber,
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & pres,
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & dPres,
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & gravCoef,
                                FluxKernel::ElementView< arrayView2d< real64 const > > const & dens,
                                FluxKernel::ElementView< arrayView2d< real64 const > > const & dDens_dPres,
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & mob,
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & dMob_dPres,
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & aperture0,
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & aperture,
                                FluxKernel::ElementView< arrayView1d< R1Tensor const > > const & transTMultiplier,
                                R1Tensor const gravityVector,
                                real64 const meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & s,
                                FluxKernel::ElementView< arrayView1d< real64 const > > const & dSdAper,
#endif
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs,
                                CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = FaceElementStencil::MAX_STENCIL_SIZE;


  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  ArrayOfArraysView< integer const > const & isGhostConnectors = stencil.getIsGhostConnectors();

  static constexpr real64 TINY = 1e-10;

  forAll< parallelDevicePolicy< 256 > >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = seri.sizeOfArray( iconn );
    localIndex const stencilSize  = numFluxElems;

    if( numFluxElems > 1 && isGhostConnectors[iconn][0] < 0 )
    {
      // working arrays
      stackArray1d< globalIndex, maxNumFluxElems > eqnRowIndices( numFluxElems );
      stackArray1d< globalIndex, maxStencilSize > dofColIndices( stencilSize );

      stackArray1d< localIndex, maxNumFluxElems > localRowIndices( numFluxElems );
      stackArray1d< localIndex, maxNumFluxElems > localColIndices( numFluxElems );

      stackArray1d< real64, maxNumFluxElems > localFlux( numFluxElems );
      stackArray2d< real64, maxNumFluxElems *maxStencilSize > localFluxJacobian( numFluxElems, stencilSize );

      // need to store this for later use in determining the dFlux_dU terms when using better permeabilty
      // approximations.
      stackArray2d< real64, maxNumFluxElems *maxStencilSize > dFlux_dAper( numFluxElems, stencilSize );

      localIndex const er = seri[iconn][0];
      localIndex const esr = sesri[iconn][0];

      // check if connection is vertical or horizontal
      stackArray1d< real64, maxNumFluxElems > effectiveWeights( numFluxElems );

      for( localIndex k = 0; k < numFluxElems; ++k )
      {

        effectiveWeights[k] = weights[iconn][k];

        localIndex const ei = sei[iconn][k];

        if( fabs( Dot( cellCenterToEdgeCenters[iconn][k], gravityVector )) > TINY )
        {
          effectiveWeights[k] *= transTMultiplier[er][esr][ei][1];
        }
        else
        {
          effectiveWeights[k] *= transTMultiplier[er][esr][ei][0];
        }

      }

      FluxKernel::ComputeJunction( numFluxElems,
                                   sei[iconn],
                                   effectiveWeights,
                                   pres[er][esr],
                                   dPres[er][esr],
                                   gravCoef[er][esr],
                                   dens[er][esr],
                                   dDens_dPres[er][esr],
                                   mob[er][esr],
                                   dMob_dPres[er][esr],
                                   aperture0[er][esr],
                                   aperture[er][esr],
                                   meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                                   s[er][esr],
                                   dSdAper[er][esr],
#endif
                                   dt,
                                   localFlux,
                                   localFluxJacobian,
                                   dFlux_dAper );

      // extract DOF numbers
      for( localIndex i = 0; i < numFluxElems; ++i )
      {
        eqnRowIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localRowIndices[i] = sei( iconn, i );
      }

      for( localIndex i = 0; i < stencilSize; ++i )
      {
        dofColIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localColIndices[i] = sei( iconn, i );
      }

      addLocalContributionsToGlobalSystem< parallelDeviceAtomic >( numFluxElems,
                                                                   stencilSize,
                                                                   eqnRowIndices.data(),
                                                                   dofColIndices.data(),
                                                                   rankOffset,
                                                                   localFluxJacobian.data(),
                                                                   localFlux.data(),
                                                                   localMatrix,
                                                                   localRhs );

      for( localIndex row=0; row<numFluxElems; ++row )
      {
        dR_dAper.addToRowBinarySearch< parallelDeviceAtomic >( localRowIndices[row],
                                                               localColIndices.data(),
                                                               dFlux_dAper.data() + ( stencilSize * row ),
                                                               stencilSize );
      }

    }

  } );
}


} // namespace SinglePhaseFVMKernels

} // namespace geosx
