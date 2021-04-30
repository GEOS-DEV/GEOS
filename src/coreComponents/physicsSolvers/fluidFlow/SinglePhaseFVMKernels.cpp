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
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"

namespace geosx
{

namespace SinglePhaseFVMKernels
{
using namespace FluxKernelsHelper;

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

  real64 fluxVal = 0.0;
  real64 dFlux_dTrans[2] = {0.0, 0.0};
  real64 dFlux_dP[2] = {0.0, 0.0};

  computeSinglePhaseFlux( seri, sesri, sei,
                          transmissibility,
                          dTrans_dPres,
                          pres,
                          dPres,
                          gravCoef,
                          dens,
                          dDens_dPres,
                          mob,
                          dMob_dPres,
                          fluxVal ,
                          dFlux_dP,
                          dFlux_dTrans );

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    dFlux_dAper[ke] = dFlux_dTrans[ke] * dTrans_dAper[ke];
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

}// namespace SinglePhaseFVMKernels

} // namespace geosx
