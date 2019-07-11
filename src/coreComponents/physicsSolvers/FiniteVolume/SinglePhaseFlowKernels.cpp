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
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

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
                                                 Epetra_FECrsMatrix * const jacobian,
                                                 Epetra_FEVector * const residual )
{

  // Add to global residual/jacobian
  jacobian->SumIntoGlobalValues( integer_conversion<int>(numFluxElems),
                                 eqnRowIndices,
                                 integer_conversion<int>(stencilSize),
                                 dofColIndices,
                                 localFluxJacobian );

  residual->SumIntoGlobalValues( integer_conversion<int>(numFluxElems),
                                 eqnRowIndices,
                                 localFlux );

}

template<>
void FluxKernel::
Launch<CellElementStencilTPFA>( CellElementStencilTPFA const & stencil,
                                real64 const dt,
                                localIndex const fluidIndex,
                                integer const gravityFlag,
                                ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > const & dofNumber,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & pres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & gravDepth,
                                FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
                                FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & mob,
                                FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres,
                                Epetra_FECrsMatrix * const jacobian,
                                Epetra_FEVector * const residual )
{
  constexpr localIndex maxNumFluxElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex numFluxElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = CellElementStencilTPFA::MAX_STENCIL_SIZE;
  constexpr localIndex stencilSize  = CellElementStencilTPFA::MAX_STENCIL_SIZE;

  typename CellElementStencilTPFA::INDEX_VIEW_CONST_TYPE const & seri = stencil.getElementRegionIndices();
  typename CellElementStencilTPFA::INDEX_VIEW_CONST_TYPE const & sesri = stencil.getElementSubRegionIndices();
  typename CellElementStencilTPFA::INDEX_VIEW_CONST_TYPE const & sei = stencil.getElementIndices();
  typename CellElementStencilTPFA::WEIGHT_VIEW_CONST_TYPE const & weights = stencil.getWeights();

  forall_in_range<stencilPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {
    // working arrays
    stackArray1d<globalIndex, numFluxElems> eqnRowIndices(numFluxElems);
    stackArray1d<globalIndex, maxNumFluxElems> dofColIndices(stencilSize);

    stackArray1d<real64, maxNumFluxElems> localFlux(numFluxElems);
    stackArray2d<real64, maxNumFluxElems*maxStencilSize> localFluxJacobian(numFluxElems, stencilSize);

    localIndex constexpr numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
    localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;

    stackArray1d<real64, numElems>   densWeight(numElems);
    stackArray1d<real64, maxStencil> dDensMean_dP(stencilSize);
    stackArray1d<real64, maxStencil> dFlux_dP(stencilSize);

    // clear working arrays
    dDensMean_dP = 0.0;

    // density averaging weights
    densWeight = 0.5;

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    for (localIndex ke = 0; ke < numElems; ++ke)
    {
      // density
      real64 const density = dens[seri[iconn][ke]][sesri[iconn][ke]][fluidIndex][sei[iconn][ke]][0];
      real64 const dDens_dP = dDens_dPres[seri[iconn][ke]][sesri[iconn][ke]][fluidIndex][sei[iconn][ke]][0];

      // average density
      densMean        += densWeight[ke] * density;
      dDensMean_dP[ke] = densWeight[ke] * dDens_dP;
    }

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localIndex const er  = seri[iconn][ke];
      localIndex const esr = sesri[iconn][ke];
      localIndex const ei  = sei[iconn][ke];

      real64 weight = weights[iconn][ke];

      real64 const gravD = gravDepth[er][esr][ei];
      real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
      real64 const dGrav_dP = gravityFlag ? dDensMean_dP[ke] * gravD : 0.0;

      potDif += weight * (pres[er][esr][ei] + dPres[er][esr][ei] - gravTerm);
      dFlux_dP[ke] = weight * (1.0 - dGrav_dP);
    }

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    localIndex er_up  = seri[iconn][k_up];
    localIndex esr_up = sesri[iconn][k_up];
    localIndex ei_up  = sei[iconn][k_up];

    real64 const mobility     = mob[er_up][esr_up][ei_up];
    real64 const dMobility_dP = dMob_dPres[er_up][esr_up][ei_up];

    // compute the final flux and derivatives
    real64 const fluxVal = mobility * potDif;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      dFlux_dP[ke] *= mobility;
    }

    dFlux_dP[k_up] += dMobility_dP * potDif;

    // populate local flux vector and derivatives
    localFlux[0] = dt * fluxVal;
    localFlux[1] = -localFlux[0];

    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localFluxJacobian[0][ke] = dt * dFlux_dP[ke];
      localFluxJacobian[1][ke] = -localFluxJacobian[0][ke];
    }


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
                            ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > const & dofNumber,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & pres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & gravDepth,
                            FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens,
                            FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & mob,
                            FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres,
                            Epetra_FECrsMatrix * const jacobian,
                            Epetra_FEVector * const residual )
{
  constexpr localIndex maxNumFluxElems = FaceElementStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = FaceElementStencil::MAX_STENCIL_SIZE;


  typename FaceElementStencil::INDEX_VIEW_CONST_TYPE const & seri = stencil.getElementRegionIndices();
  typename FaceElementStencil::INDEX_VIEW_CONST_TYPE const & sesri = stencil.getElementSubRegionIndices();
  typename FaceElementStencil::INDEX_VIEW_CONST_TYPE const & sei = stencil.getElementIndices();
  typename FaceElementStencil::WEIGHT_VIEW_CONST_TYPE const & weights = stencil.getWeights();

  forall_in_range<stencilPolicy>( 0, stencil.size(), GEOSX_LAMBDA ( localIndex iconn )
  {
    localIndex const numFluxElems = stencil.stencilSize(iconn);
    localIndex const stencilSize  = numFluxElems;

    // working arrays
    stackArray1d<globalIndex, maxNumFluxElems> eqnRowIndices(numFluxElems);
    stackArray1d<globalIndex, maxStencilSize> dofColIndices(stencilSize);

    stackArray1d<real64, maxNumFluxElems> localFlux(numFluxElems);
    stackArray2d<real64, maxNumFluxElems*maxStencilSize> localFluxJacobian(numFluxElems, stencilSize);

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    real64 sumOfWeights = 0;
    for( localIndex k=0 ; k<numFluxElems ; ++k )
    {
      sumOfWeights += weights[iconn][k];
    }

    localIndex k[2];
    for( k[0]=0 ; k[0]<numFluxElems ; ++k[0] )
    {
      for( k[1]=k[0]+1 ; k[1]<numFluxElems ; ++k[1] )
      {
        real64  dFlux_dP[2] = {0,0};

        localIndex const ei[2] = { sei[iconn][k[0]],
                                   sei[iconn][k[1]] };

        real64 const weight[2] = {   weights[iconn][k[0]] * weights[iconn][k[1]] / sumOfWeights,
                                   - weights[iconn][k[0]] * weights[iconn][k[1]] / sumOfWeights };


        // average density
        real64 const densMean = 0.5 * ( dens[er][esr][fluidIndex][ei[0]][0] + dens[er][esr][fluidIndex][ei[1]][0] );

        real64 const dDensMean_dP[2] = { 0.5 * dDens_dPres[er][esr][fluidIndex][ei[0]][0],
                                         0.5 * dDens_dPres[er][esr][fluidIndex][ei[1]][0] };

        real64 potDif = 0.0;
        for( localIndex i = 0 ; i < 2 ; ++i )
        {
          real64 const gravD = gravDepth[er][esr][ei[i]];
          real64 const gravTerm = gravityFlag ? densMean * gravD : 0.0;
          real64 const dGrav_dP = gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

          potDif += weight[i] * ( pres[er][esr][ei[i]] +
                                  dPres[er][esr][ei[i]] - gravTerm);
          dFlux_dP[i] = weight[i] * (1.0 - dGrav_dP);
        }

        // upwinding of fluid properties (make this an option?)
        localIndex const k_up = (potDif >= 0) ? 0 : 1;

        localIndex ei_up  = sei[iconn][k[k_up]];

        real64 const mobility     = mob[er][esr][ei_up];
        real64 const dMobility_dP = dMob_dPres[er][esr][ei_up];

        // compute the final flux and derivatives
        real64 const fluxVal = mobility * potDif;
        dFlux_dP[0] *= mobility;
        dFlux_dP[1] *= mobility;

        dFlux_dP[k_up] += dMobility_dP * potDif;

        // populate local flux vector and derivatives
        localFlux[k[0]] += dt * fluxVal;
        localFlux[k[1]] -= dt * fluxVal;

        localFluxJacobian[k[0]][k[0]] += dt * dFlux_dP[0];
        localFluxJacobian[k[1]][k[0]] -= dt * dFlux_dP[0];
        localFluxJacobian[k[0]][k[1]] += dt * dFlux_dP[1];
        localFluxJacobian[k[1]][k[1]] -= dt * dFlux_dP[1];

      }
    }

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
