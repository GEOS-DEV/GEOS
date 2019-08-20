/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SinglePhaseFlowKernels.cpp
 */

#include "SinglePhaseFlowKernels.hpp"

namespace geosx
{

namespace SinglePhaseFlowKernels
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

void MobilityKernel::Launch( set<localIndex> targetSet,
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

void MobilityKernel::Launch( set<localIndex> targetSet,
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
                                integer const gravityFlag,
                                FluxKernel::ElementView< arrayView1d<globalIndex> > const & dofNumber,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & pres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & gravDepth,
                                FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
                                FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & mob,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const &,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const &,
                                ParallelMatrix * const jacobian,
                                ParallelVector * const residual )
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
                         gravDepth,
                         dens,
                         dDens_dPres,
                         mob,
                         dMob_dPres,
                         fluidIndex,
                         gravityFlag,
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
                            integer const gravityFlag,
                            FluxKernel::ElementView < arrayView1d<globalIndex const> > const & dofNumber,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & pres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & gravDepth,
                            FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
                            FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & mob,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture0,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture,
                            ParallelMatrix * const jacobian,
                            ParallelVector * const residual )
{
  constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = FaceElementStencil::MAX_STENCIL_SIZE;


  typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {
    localIndex const numFluxElems = stencil.stencilSize(iconn);
    localIndex const stencilSize  = numFluxElems;

    // working arrays
    stackArray1d<globalIndex, maxNumFluxElems> eqnRowIndices(numFluxElems);
    stackArray1d<globalIndex, maxStencilSize> dofColIndices(stencilSize);

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
                                 gravDepth[er][esr],
                                 dens[er][esr][fluidIndex],
                                 dDens_dPres[er][esr][fluidIndex],
                                 mob[er][esr],
                                 dMob_dPres[er][esr],
                                 aperture0[er][esr],
                                 aperture[er][esr],
                                 fluidIndex,
                                 gravityFlag,
                                 dt,
                                 localFlux,
                                 localFluxJacobian,
                                 dFlux_dAper );

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




} // namespace SinglePhaseFlowKernels

} // namespace geosx
