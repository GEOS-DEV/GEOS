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
 * @file SinglePhaseKernels.cpp
 */

#include "SinglePhaseKernels.hpp"

namespace geosx
{

namespace SinglePhaseKernels
{

/******************************** MobilityKernel ********************************/

void
MobilityKernel::Compute( real64 const & dens,
                         real64 const & dDens_dPres,
                         real64 const & visc,
                         real64 const & dVisc_dPres,
                         real64 & mob,
                         real64 & dMob_dPres )
{
  mob = dens / visc;
  dMob_dPres = dDens_dPres / visc - mob / visc * dVisc_dPres;
}

void
MobilityKernel::Compute( real64 const & dens,
                         real64 const & visc,
                         real64 & mob )
{
  mob = dens / visc;
}

void MobilityKernel::Launch( localIndex begin, localIndex end,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & dDens_dPres,
                             arrayView2d<real64 const> const & visc,
                             arrayView2d<real64 const> const & dVisc_dPres,
                             arrayView1d<real64> const & mob,
                             arrayView1d<real64> const & dMob_dPres )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( dens[a][0],
             dDens_dPres[a][0],
             visc[a][0],
             dVisc_dPres[a][0],
             mob[a],
             dMob_dPres[a] );
  } );
}

void MobilityKernel::Launch( SortedArray<localIndex> targetSet,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & dDens_dPres,
                             arrayView2d<real64 const> const & visc,
                             arrayView2d<real64 const> const & dVisc_dPres,
                             arrayView1d<real64> const & mob,
                             arrayView1d<real64> const & dMob_dPres )
{
  forall_in_set( targetSet.values(), targetSet.size(), GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( dens[a][0],
             dDens_dPres[a][0],
             visc[a][0],
             dVisc_dPres[a][0],
             mob[a],
             dMob_dPres[a] );
  } );
}

void MobilityKernel::Launch( localIndex begin, localIndex end,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & visc,
                             arrayView1d<real64> const & mob )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( dens[a][0],
             visc[a][0],
             mob[a] );
  } );
}

void MobilityKernel::Launch( SortedArray<localIndex> targetSet,
                             arrayView2d<real64 const> const & dens,
                             arrayView2d<real64 const> const & visc,
                             arrayView1d<real64> const & mob )
{
  forall_in_set( targetSet.values(), targetSet.size(), GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( dens[a][0],
             visc[a][0],
             mob[a] );
  } );
}

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
                 numFluxElems);

}

template<>
void FluxKernel::
Launch<CellElementStencilTPFA>( CellElementStencilTPFA const & stencil,
                                real64 const dt,
                                localIndex const fluidIndex,
                                FluxKernel::ElementView< arrayView1d<globalIndex> > const & dofNumber,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & pres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & gravCoef,
                                FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
                                FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & mob,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const &,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const &,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                                FluxKernel::ElementView < arrayView1d<real64 const> > const &,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const &,
#endif
                                ParallelMatrix * const jacobian,
                                ParallelVector * const residual,
                                CRSMatrixView<real64,localIndex,localIndex const > const & )
{
  constexpr localIndex maxNumFluxElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex numFluxElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = CellElementStencilTPFA::MAX_STENCIL_SIZE;
  constexpr localIndex stencilSize  = CellElementStencilTPFA::MAX_STENCIL_SIZE;

  typename CellElementStencilTPFA::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename CellElementStencilTPFA::WeightContainerViewConstType const & weights = stencil.getWeights();

  forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {
    // working arrays
    stackArray1d<globalIndex, numFluxElems> eqnRowIndices(numFluxElems);
    stackArray1d<globalIndex, maxNumFluxElems> dofColIndices(stencilSize);

    stackArray1d<real64, maxNumFluxElems> localFlux(numFluxElems);
    stackArray2d<real64, maxNumFluxElems*maxStencilSize> localFluxJacobian(numFluxElems, stencilSize);

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
                         fluidIndex,
                         dt,
                         localFlux,
                         localFluxJacobian );

    // extract DOF numbers
    eqnRowIndices = -1;
    for (localIndex i = 0; i < numFluxElems; ++i)
    {
      eqnRowIndices[i] = dofNumber[seri(iconn,i)][sesri(iconn,i)][sei(iconn,i)];
    }

    for (localIndex i = 0; i < stencilSize; ++i)
    {
      dofColIndices[i] = dofNumber[seri(iconn,i)][sesri(iconn,i)][sei(iconn,i)];
    }

    addLocalContributionsToGlobalSystem( numFluxElems,
                                         stencilSize,
                                         eqnRowIndices.data(),
                                         dofColIndices.data(),
                                         localFluxJacobian.data(),
                                         localFlux.data(),
                                         jacobian,
                                         residual );
  } );
}

template<>
void FluxKernel::
Launch<FaceElementStencil>( FaceElementStencil const & stencil,
                            real64 const dt,
                            localIndex const fluidIndex,
                            FluxKernel::ElementView < arrayView1d<globalIndex const> > const & dofNumber,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & pres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & gravCoef,
                            FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
                            FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & mob,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture0,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & s,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & dSdAper,
#endif
                            ParallelMatrix * const jacobian,
                            ParallelVector * const residual,
                            CRSMatrixView<real64,localIndex,localIndex const > const & dR_dAper )
{
  constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = FaceElementStencil::MAX_STENCIL_SIZE;


  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

//  {
//    localIndex const numRows = dR_dAper.numRows();
//    for( localIndex ei=0 ; ei<numRows ; ++ei )
//    {
//      localIndex const numColumns = dR_dAper.numNonZeros(ei);
//      arraySlice1d<localIndex const> const & columns = dR_dAper.getColumns( ei );
//      arraySlice1d<real64 const> const & values = dR_dAper.getEntries( ei );
//
//      for( localIndex kfe2=0 ; kfe2<numColumns ; ++kfe2 )
//      {
//        real64 dRdAper = values[kfe2];
//        localIndex const ei2 = columns[kfe2];
//        GEOS_LOG_RANK( "dR_dAper("<<ei<<", "<<ei2<<") = "<<dRdAper );
//      }
//    }
//  }

  forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {
    localIndex const numFluxElems = stencil.stencilSize(iconn);
    localIndex const stencilSize  = numFluxElems;

    // working arrays
    stackArray1d<globalIndex, maxNumFluxElems> eqnRowIndices(numFluxElems);
    stackArray1d<globalIndex, maxStencilSize> dofColIndices(stencilSize);

    stackArray1d<localIndex, maxNumFluxElems> localRowIndices(numFluxElems);
    stackArray1d<localIndex, maxNumFluxElems> localColIndices(numFluxElems);

    stackArray1d<real64, maxNumFluxElems> localFlux(numFluxElems);
    stackArray2d<real64, maxNumFluxElems*maxStencilSize> localFluxJacobian(numFluxElems, stencilSize);

    // need to store this for later use in determining the dFlux_dU terms when using better permeabilty approximations.
    stackArray2d<real64, maxNumFluxElems*maxStencilSize> dFlux_dAper(numFluxElems, stencilSize);

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];


    FluxKernel::ComputeJunction( numFluxElems,
                                 sei[iconn],
                                 weights[iconn],
                                 pres[er][esr],
                                 dPres[er][esr],
                                 gravCoef[er][esr],
                                 dens[er][esr][fluidIndex],
                                 dDens_dPres[er][esr][fluidIndex],
                                 mob[er][esr],
                                 dMob_dPres[er][esr],
                                 aperture0[er][esr],
                                 aperture[er][esr],
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                                 s[er][esr],
                                 dSdAper[er][esr],
#endif
                                 fluidIndex,
                                 dt,
                                 localFlux,
                                 localFluxJacobian,
                                 dFlux_dAper );

    // extract DOF numbers
    eqnRowIndices = -1;
    for (localIndex i = 0; i < numFluxElems; ++i)
    {
      eqnRowIndices[i] = dofNumber[seri(iconn,i)][sesri(iconn,i)][sei(iconn,i)];
      localRowIndices[i] = sei(iconn,i);
    }

    for (localIndex i = 0; i < stencilSize; ++i)
    {
      dofColIndices[i] = dofNumber[seri(iconn,i)][sesri(iconn,i)][sei(iconn,i)];
      localColIndices[i] = sei(iconn,i);
    }

    addLocalContributionsToGlobalSystem( numFluxElems,
                                         stencilSize,
                                         eqnRowIndices.data(),
                                         dofColIndices.data(),
                                         localFluxJacobian.data(),
                                         localFlux.data(),
                                         jacobian,
                                         residual );

//    for( localIndex a=0 ;a<numFluxElems ; ++a )
//    {
//      for( localIndex b=0 ; b<stencilSize ; ++b )
//      {
//        GEOS_LOG_RANK("dFlux_dAper("<<localRowIndices[a]<<", "<<localColIndices[b]<<") = "<<dFlux_dAper(a,b) );
//      }
//    }
    for( localIndex row=0 ; row<numFluxElems ; ++row )
    {
      dR_dAper.addToRowBinarySearch( localRowIndices[row],
                                     localColIndices.data(),
                                     dFlux_dAper.data() + (stencilSize * row),
                                     stencilSize );
    }

  } );

//  localIndex const numRows = dR_dAper.numRows();
//  for( localIndex ei=0 ; ei<numRows ; ++ei )
//  {
//    localIndex const numColumns = dR_dAper.numNonZeros(ei);
//    arraySlice1d<localIndex const> const & columns = dR_dAper.getColumns( ei );
//    arraySlice1d<real64 const> const & values = dR_dAper.getEntries( ei );
//
//    for( localIndex kfe2=0 ; kfe2<numColumns ; ++kfe2 )
//    {
//      real64 dRdAper = values[kfe2];
//      localIndex const ei2 = columns[kfe2];
//      GEOS_LOG_RANK( "dR_dAper("<<ei<<", "<<ei2<<") = "<<dRdAper );
//    }
//  }

}




} // namespace SinglePhaseKernels

} // namespace geosx
