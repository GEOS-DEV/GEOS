/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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
GEOSX_HOST_DEVICE
void
FluxKernel::compute( localIndex const numFluxElems,
                     arraySlice1d< localIndex const > const & seri,
                     arraySlice1d< localIndex const > const & sesri,
                     arraySlice1d< localIndex const > const & sei,
                     real64 const (&transmissibility)[2],
                     real64 const (&dTrans_dPres)[2],
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
  localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;
  stackArray1d< real64, maxStencil > dDensMean_dP( 2 );
  stackArray1d< real64, maxStencil > dFlux_dP( 2 );

  // average density
  real64 densMean = 0.0;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    densMean        += 0.5 * dens[seri[ke]][sesri[ke]][sei[ke]][0];
    dDensMean_dP[ke] = 0.5 * dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0];
  }

  // compute potential difference
  real64 potDif = 0.0;
  real64 sumWeightGrav = 0.0;
  real64 potScale = 0.0;

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    localIndex const er  = seri[ke];
    localIndex const esr = sesri[ke];
    localIndex const ei  = sei[ke];

    real64 const pressure = pres[er][esr][ei] + dPres[er][esr][ei];
    real64 const gravD = gravCoef[er][esr][ei];
    real64 const pot = transmissibility[ke] * ( pressure - densMean * gravD );

    potDif += pot;
    sumWeightGrav += gravD;
    potScale = fmax( potScale, fabs( pot ) );
  }

  // compute upwinding tolerance
  real64 constexpr upwRelTol = 1e-8;
  real64 const upwAbsTol = fmax( potScale * upwRelTol, LvArray::NumericLimits< real64 >::epsilon );

  // decide mobility coefficients - smooth variation in [-upwAbsTol; upwAbsTol]
  real64 const alpha = ( potDif + upwAbsTol ) / ( 2 * upwAbsTol );

  real64 mobility{};
  real64 dMobility_dP[2]{};
  if( alpha <= 0.0 || alpha >= 1.0 )
  {
    // happy path: single upwind direction
    localIndex const ke = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );
    mobility = mob[seri[ke]][sesri[ke]][sei[ke]];
    dMobility_dP[ke] = dMob_dPres[seri[ke]][sesri[ke]][sei[ke]];
  }
  else
  {
    // sad path: weighted averaging
    real64 const mobWeights[2] = { alpha, 1.0 - alpha };
    for( localIndex ke = 0; ke < 2; ++ke )
    {
      mobility += mobWeights[ke] * mob[seri[ke]][sesri[ke]][sei[ke]];
      dMobility_dP[ke] = mobWeights[ke] * dMob_dPres[seri[ke]][sesri[ke]][sei[ke]];
    }
  }

  // compute the final flux and derivatives
  real64 const fluxVal = mobility * potDif;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    real64 const dFlux_dTrans = mobility  * (  1 - dDensMean_dP[ke] * sumWeightGrav ) + dMobility_dP[ke] * potDif;

    dFlux_dP[ke] = mobility * transmissibility[ke] * (  1 - dDensMean_dP[ke] * sumWeightGrav ) +
                   dFlux_dTrans * dTrans_dPres[ke];
  }

  // populate local flux vector and derivatives
  flux[0] =  dt * fluxVal;
  flux[1] = -dt * fluxVal;

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    fluxJacobian[0][ke] =  dt * dFlux_dP[ke];
    fluxJacobian[1][ke] = -dt * dFlux_dP[ke];
  }
}

void PoroelasticFluxKernel::compute( localIndex const numFluxElems,
                                     arraySlice1d< localIndex const > const & seri,
                                     arraySlice1d< localIndex const > const & sesri,
                                     arraySlice1d< localIndex const > const & sei,
                                     real64 const (&transmissibility)[2],
                                     real64 const (&dTrans_dPres)[2],                                     ElementViewConst< arrayView1d< real64 const > > const & pres,
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
  localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;
  stackArray1d< real64, maxStencil > dDensMean_dP( 2 );
  stackArray1d< real64, maxStencil > dFlux_dP( 2 );

  // average density
  real64 densMean = 0.0;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    densMean        += 0.5 * dens[seri[ke]][sesri[ke]][sei[ke]][0];
    dDensMean_dP[ke] = 0.5 * dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0];
  }

  // compute potential difference
  real64 potDif = 0.0;
  real64 sumWeightGrav = 0.0;
  real64 potScale = 0.0;

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    localIndex const er  = seri[ke];
    localIndex const esr = sesri[ke];
    localIndex const ei  = sei[ke];

    real64 const pressure = pres[er][esr][ei] + dPres[er][esr][ei];
    real64 const gravD = gravCoef[er][esr][ei];
    real64 const pot = transmissibility[ke] * ( pressure - densMean * gravD );

    potDif += pot;
    sumWeightGrav += gravD;
    potScale = fmax( potScale, fabs( pot ) );
  }

  // compute upwinding tolerance
  real64 constexpr upwRelTol = 1e-8;
  real64 const upwAbsTol = fmax( potScale * upwRelTol, LvArray::NumericLimits< real64 >::epsilon );

  // decide mobility coefficients - smooth variation in [-upwAbsTol; upwAbsTol]
  real64 const alpha = ( potDif + upwAbsTol ) / ( 2 * upwAbsTol );

  real64 mobility{};
  real64 dMobility_dP[2]{};
  if( alpha <= 0.0 || alpha >= 1.0 )
  {
    // happy path: single upwind direction
    localIndex const ke = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );
    mobility = mob[seri[ke]][sesri[ke]][sei[ke]];
    dMobility_dP[ke] = dMob_dPres[seri[ke]][sesri[ke]][sei[ke]];
  }
  else
  {
    // sad path: weighted averaging
    real64 const mobWeights[2] = { alpha, 1.0 - alpha };
    for( localIndex ke = 0; ke < 2; ++ke )
    {
      mobility += mobWeights[ke] * mob[seri[ke]][sesri[ke]][sei[ke]];
      dMobility_dP[ke] = mobWeights[ke] * dMob_dPres[seri[ke]][sesri[ke]][sei[ke]];
    }
  }

  // compute the final flux and derivatives
  real64 const fluxVal = mobility * potDif;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    real64 const dFlux_dTrans = mobility  * (  1 - dDensMean_dP[ke] * sumWeightGrav ) + dMobility_dP[ke] * potDif;

    dFlux_dP[ke] = mobility * transmissibility[ke] * (  1 - dDensMean_dP[ke] * sumWeightGrav ) +
        dFlux_dTrans * dTrans_dPres[ke];
  }

  // populate local flux vector and derivatives
  flux[0] =  dt * fluxVal;
  flux[1] = -dt * fluxVal;

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    fluxJacobian[0][ke] =  dt * dFlux_dP[ke];
    fluxJacobian[1][ke] = -dt * dFlux_dP[ke];
  }
}


//
//GEOSX_HOST_DEVICE
//void
//FluxKernel::computeJunction( localIndex const numFluxElems,
//                             arraySlice1d< localIndex const > const & stencilElementIndices,
//                             arraySlice1d< real64 const > const & stencilWeights,
//                             arrayView1d< real64 const > const & pres,
//                             arrayView1d< real64 const > const & dPres,
//                             arrayView1d< real64 const > const & gravCoef,
//                             arrayView2d< real64 const > const & dens,
//                             arrayView2d< real64 const > const & dDens_dPres,
//                             arrayView1d< real64 const > const & mob,
//                             arrayView1d< real64 const > const & dMob_dPres,
//                             arrayView1d< real64 const > const & aperture0,
//                             arrayView1d< real64 const > const & aperture,
//                             real64 const meanPermCoeff,
//#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
//                             arrayView1d< real64 const > const & GEOSX_GEOSX_UNUSED_PARAM( s ),
//                             arrayView1d< real64 const > const & GEOSX_GEOSX_UNUSED_PARAM( dSdAper ),
//#endif
//                             real64 const dt,
//                             arraySlice1d< real64 > const & flux,
//                             arraySlice2d< real64 > const & fluxJacobian,
//                             arraySlice2d< real64 > const & dFlux_dAperture )
//{
//  real64 sumOfWeights = 0;
//  real64 aperTerm[10];
//  real64 dAperTerm_dAper[10];
//
//  for( localIndex k=0; k<numFluxElems; ++k )
//  {
//
//    #define PERM_CALC 1
////      real64 const aperAdd = aperture0[stencilElementIndices[k]] < 0.09e-3 ? ( 0.09e-3 -
//// aperture0[stencilElementIndices[k]] ) : 0.0;
//#if PERM_CALC==1
//    FluxKernelHelper::
//      apertureForPermeablityCalculation< 1 >( aperture0[stencilElementIndices[k]],
//                                              aperture[stencilElementIndices[k]],
//                                              aperTerm[k],
//                                              dAperTerm_dAper[k] );
//
//#elif PERM_CALC==2
//
//    if( s[k] >= 1.0 )
//    {
//      aperTerm[k] = aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]];
//      dAperTerm_dAper[k] = 3*aperture[stencilElementIndices[k]]*aperture[stencilElementIndices[k]];
//    }
//    else
//    {
//      aperTerm[k] = aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]]/s[k];
//      dAperTerm_dAper[k] = 3*aperture[stencilElementIndices[k]]*aperture[stencilElementIndices[k]]/s[k]
//                           - aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] *
// aperture[stencilElementIndices[k]]/(s[k]*s[k]) *
//                           dSdAper[k];
//    }
//#endif
////      aperTerm[k] += aperAdd*aperAdd*aperAdd;
//
//
//    sumOfWeights += aperTerm[k] * stencilWeights[k];
//  }
//
//  localIndex k[2];
//  for( k[0]=0; k[0]<numFluxElems; ++k[0] )
//  {
//    for( k[1]=k[0]+1; k[1]<numFluxElems; ++k[1] )
//    {
//      real64 dFlux_dP[2] = {0, 0};
//
//      localIndex const ei[2] = { stencilElementIndices[k[0]],
//                                 stencilElementIndices[k[1]] };
//#if 0
//      real64 const weight = ( stencilWeights[k[0]]*aperTerm[k[0]] ) *
//                            ( stencilWeights[k[1]]*aperTerm[k[1]] ) / sumOfWeights;
//
//      real64 const
//      dWeight_dAper[2] =
//      { ( 1 / aperTerm[k[0]]  - stencilWeights[k[0]] / sumOfWeights ) * weight * dAperTerm_dAper[k[0]],
//        ( 1 / aperTerm[k[1]]  - stencilWeights[k[1]] / sumOfWeights ) * weight * dAperTerm_dAper[k[1]]};
//#else
//      real64 const c = meanPermCoeff;
//
//      real64 const harmonicWeight = ( stencilWeights[k[0]]*aperTerm[k[0]] ) *
//                                    ( stencilWeights[k[1]]*aperTerm[k[1]] ) / sumOfWeights;
//
//      real64 const weight = c * harmonicWeight
//                            + (1.0 - c) * 0.25 * ( stencilWeights[k[0]]*aperTerm[k[0]] + stencilWeights[k[1]]*aperTerm[k[1]] );
//
//      real64 const
//      dHarmonicWeight_dAper[2] =
//      { ( 1 / aperTerm[k[0]]  - stencilWeights[k[0]] / sumOfWeights ) * harmonicWeight * dAperTerm_dAper[k[0]],
//        ( 1 / aperTerm[k[1]]  - stencilWeights[k[1]] / sumOfWeights ) * harmonicWeight * dAperTerm_dAper[k[1]]};
//
//      real64 const
//      dWeight_dAper[2] =
//      { c * dHarmonicWeight_dAper[0] + 0.25 * ( 1.0 - c )*stencilWeights[k[0]]*dAperTerm_dAper[k[0]],
//        c * dHarmonicWeight_dAper[1] + 0.25 * ( 1.0 - c )*stencilWeights[k[1]]*dAperTerm_dAper[k[1]] };
//
//#endif
//      // average density
//      real64 const densMean = 0.5 * ( dens[ei[0]][0] + dens[ei[1]][0] );
//
//      real64 const dDensMean_dP[2] = { 0.5 * dDens_dPres[ei[0]][0],
//                                       0.5 * dDens_dPres[ei[1]][0] };
//
//      real64 const potDif =  ( ( pres[ei[0]] + dPres[ei[0]] ) - ( pres[ei[1]] + dPres[ei[1]] ) -
//                               densMean * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) );
//
//
//      // upwinding of fluid properties (make this an option?)
//      localIndex const k_up = (potDif >= 0) ? 0 : 1;
//
//      localIndex ei_up  = stencilElementIndices[k[k_up]];
//
//      real64 const mobility     = mob[ei_up];
//      real64 const dMobility_dP = dMob_dPres[ei_up];
//
//      // Compute flux and fill flux rval
//      real64 const fluxVal = mobility * weight * potDif * dt;
//      flux[k[0]] += fluxVal;
//      flux[k[1]] -= fluxVal;
//
//
//      // compute and fill dFlux_dP
//      dFlux_dP[0] = mobility * weight * (  1 - dDensMean_dP[0] * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) ) * dt;
//      dFlux_dP[1] = mobility * weight * ( -1 - dDensMean_dP[1] * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) ) * dt;
//      dFlux_dP[k_up] += dMobility_dP * weight * potDif * dt;
//
//      fluxJacobian[k[0]][k[0]] += dFlux_dP[0];
//      fluxJacobian[k[0]][k[1]] += dFlux_dP[1];
//      fluxJacobian[k[1]][k[0]] -= dFlux_dP[0];
//      fluxJacobian[k[1]][k[1]] -= dFlux_dP[1];
//
//      real64 const dFlux_dAper[2] = { mobility * dWeight_dAper[0] * potDif * dt,
//                                      mobility * dWeight_dAper[1] * potDif * dt };
//      dFlux_dAperture[k[0]][k[0]] += dFlux_dAper[0];
//      dFlux_dAperture[k[0]][k[1]] += dFlux_dAper[1];
//      dFlux_dAperture[k[1]][k[0]] -= dFlux_dAper[0];
//      dFlux_dAperture[k[1]][k[1]] -= dFlux_dAper[1];
//    }
//  }
//}


}// namespace SinglePhaseFVMKernels

} // namespace geosx
