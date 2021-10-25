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
 * @file ThermalCompositionalMultiphaseFVMKernels.cpp
 */

#include "ThermalCompositionalMultiphaseFVMKernels.hpp"
#include "CompositionalMultiphaseUtilities.hpp"

#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/SurfaceElementStencil.hpp"
#include "finiteVolume/EmbeddedSurfaceToCellStencil.hpp"
#include "finiteVolume/FaceElementToCellStencil.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"


namespace geosx
{

namespace ThermalCompositionalMultiphaseFVMKernels
{

/******************************** PhaseMobilityKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
PhaseMobilityKernel::
  compute( arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dTemp,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseVisc,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dPres,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dTemp,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
           arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dTemp,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseMob,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & dPhaseMob_dPres,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & dPhaseMob_dTemp,
           arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const & dPhaseMob_dComp )
{
  real64 dRelPerm_dC[NC];
  real64 dDens_dC[NC];
  real64 dVisc_dC[NC];

  for( localIndex ip = 0; ip < NP; ++ip )
  {

    // compute the phase mobility only if the phase is present
    bool const phaseExists = (phaseVolFrac[ip] > 0);
    if( !phaseExists )
    {
      phaseMob[ip] = 0.;
      dPhaseMob_dPres[ip] = 0.;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseMob_dComp[ip][jc] = 0.;
      }
      continue;
    }

    real64 const density = phaseDens[ip];
    real64 const dDens_dP = dPhaseDens_dPres[ip];
    real64 const dDens_dT = dPhaseDens_dTemp[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dDens_dC );

    real64 const viscosity = phaseVisc[ip];
    real64 const dVisc_dP = dPhaseVisc_dPres[ip];
    real64 const dVisc_dT = dPhaseVisc_dTemp[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseVisc_dComp[ip], dVisc_dC );

    real64 const relPerm = phaseRelPerm[ip];
    real64 dRelPerm_dT = 0.0;
    real64 dRelPerm_dP = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dRelPerm_dC[ic] = 0.0;
    }

    for( localIndex jp = 0; jp < NP; ++jp )
    {
      real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
      dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[jp];
      dRelPerm_dT += dRelPerm_dS * dPhaseVolFrac_dTemp[jp];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[jp][jc];
      }
    }

    real64 const mobility = relPerm * density / viscosity;

    phaseMob[ip] = mobility;
    dPhaseMob_dPres[ip] = dRelPerm_dP * density / viscosity
                          + mobility * (dDens_dP / density - dVisc_dP / viscosity);
    dPhaseMob_dTemp[ip] = dRelPerm_dT * density / viscosity
                          + mobility * (dDens_dT / density - dVisc_dT / viscosity);

    // compositional derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] * density / viscosity
                                + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
    }
  }
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::
  launch( localIndex const size,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dTemp,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dTemp,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dTemp,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseMob,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dTemp,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    compute< NC, NP >( dCompFrac_dCompDens[a],
                       phaseDens[a][0],
                       dPhaseDens_dPres[a][0],
                       dPhaseDens_dTemp[a][0],
                       dPhaseDens_dComp[a][0],
                       phaseVisc[a][0],
                       dPhaseVisc_dPres[a][0],
                       dPhaseVisc_dTemp[a][0],
                       dPhaseVisc_dComp[a][0],
                       phaseRelPerm[a][0],
                       dPhaseRelPerm_dPhaseVolFrac[a][0],
                       phaseVolFrac[a],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dTemp[a],
                       dPhaseVolFrac_dComp[a],
                       phaseMob[a],
                       dPhaseMob_dPres[a],
                       dPhaseMob_dTemp[a],
                       dPhaseMob_dComp[a] );
  } );
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dTemp,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dTemp,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dTemp,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseMob,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dTemp,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    compute< NC, NP >( dCompFrac_dCompDens[a],
                       phaseDens[a][0],
                       dPhaseDens_dPres[a][0],
                       dPhaseDens_dTemp[a][0],
                       dPhaseDens_dComp[a][0],
                       phaseVisc[a][0],
                       dPhaseVisc_dPres[a][0],
                       dPhaseVisc_dTemp[a][0],
                       dPhaseVisc_dComp[a][0],
                       phaseRelPerm[a][0],
                       dPhaseRelPerm_dPhaseVolFrac[a][0],
                       phaseVolFrac[a],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dTemp[a],
                       dPhaseVolFrac_dComp[a],
                       phaseMob[a],
                       dPhaseMob_dPres[a],
                       dPhaseMob_dTemp[a],
                       dPhaseMob_dComp[a] );
  } );
}

#define INST_PhaseMobilityKernel( NC, NP ) \
  template \
  void \
  PhaseMobilityKernel:: \
    launch< NC, NP >( localIndex const size, \
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dTemp, \
                      arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dTemp, \
                      arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm, \
                      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dTemp, \
                      arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64, compflow::USD_PHASE > const & phaseMob, \
                      arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres, \
                      arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dTemp, \
                      arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp ); \
  template \
  void \
  PhaseMobilityKernel:: \
    launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dTemp, \
                      arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dTemp, \
                      arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm, \
                      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dTemp, \
                      arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64, compflow::USD_PHASE > const & phaseMob, \
                      arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres, \
                      arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dTemp, \
                      arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp )

INST_PhaseMobilityKernel( 1, 1 );
INST_PhaseMobilityKernel( 2, 1 );
INST_PhaseMobilityKernel( 3, 1 );
INST_PhaseMobilityKernel( 4, 1 );
INST_PhaseMobilityKernel( 5, 1 );

INST_PhaseMobilityKernel( 1, 2 );
INST_PhaseMobilityKernel( 2, 2 );
INST_PhaseMobilityKernel( 3, 2 );
INST_PhaseMobilityKernel( 4, 2 );
INST_PhaseMobilityKernel( 5, 2 );

INST_PhaseMobilityKernel( 1, 3 );
INST_PhaseMobilityKernel( 2, 3 );
INST_PhaseMobilityKernel( 3, 3 );
INST_PhaseMobilityKernel( 4, 3 );
INST_PhaseMobilityKernel( 5, 3 );

#undef INST_PhaseMobilityKernel

/******************************** FluxKernel ********************************/

template< localIndex NC, localIndex MAX_NUM_ELEMS, localIndex MAX_STENCIL_SIZE >
GEOSX_HOST_DEVICE
void
FluxKernel::
  compute( localIndex const numPhases,
           localIndex const stencilSize,
           localIndex const numFluxElems,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           real64 const (&transmissibility)[2],
           real64 const (&dTrans_dPres)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & temp,
           ElementViewConst< arrayView1d< real64 const > > const & dTemp,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dTemp,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dComp,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseVolFrac_dPres,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseVolFrac_dTemp,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac_dComp,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dTemp,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dComp,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dTemp,
           ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dComp,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseEnthalpy,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseEnthalpy_dPres,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseEnthalpy_dTemp,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseEnthalpy_dComp,
           ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const & phaseCapPressure,
           ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
           integer const capPressureFlag,
           real64 const dt,
           arraySlice1d< real64 > const localFlux,
           arraySlice2d< real64 > const localFluxJacobian )
{
  localIndex constexpr NDOF = NC + 2;

  // they will be useful for conduction
  GEOSX_UNUSED_VAR( temp );
  GEOSX_UNUSED_VAR( dTemp );

  stackArray1d< real64, NC > compFlux( NC );
  stackArray2d< real64, MAX_STENCIL_SIZE * NC > dCompFlux_dP( stencilSize, NC );
  stackArray2d< real64, MAX_STENCIL_SIZE * NC > dCompFlux_dT( stencilSize, NC );
  stackArray3d< real64, MAX_STENCIL_SIZE * NC * NC > dCompFlux_dC( stencilSize, NC, NC );

  real64 enthalpyFlux = 0.0;
  stackArray1d< real64, MAX_STENCIL_SIZE > dEnthalpyFlux_dP( stencilSize );
  stackArray1d< real64, MAX_STENCIL_SIZE > dEnthalpyFlux_dT( stencilSize );
  stackArray2d< real64, MAX_STENCIL_SIZE * NC > dEnthalpyFlux_dC( stencilSize, NC );

  // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    // clear working arrays
    real64 densMean{};
    stackArray1d< real64, MAX_NUM_ELEMS > dDensMean_dP( numFluxElems );
    stackArray1d< real64, MAX_NUM_ELEMS > dDensMean_dT( numFluxElems );
    stackArray2d< real64, MAX_NUM_ELEMS * NC > dDensMean_dC( numFluxElems, NC );

    // create local work arrays
    real64 phaseFlux{};
    real64 dPhaseFlux_dP[MAX_STENCIL_SIZE]{};
    real64 dPhaseFlux_dT[MAX_STENCIL_SIZE]{};
    real64 dPhaseFlux_dC[MAX_STENCIL_SIZE][NC]{};

    real64 presGrad{};
    stackArray1d< real64, MAX_STENCIL_SIZE > dPresGrad_dP( stencilSize );
    stackArray1d< real64, MAX_STENCIL_SIZE > dPresGrad_dT( stencilSize );
    stackArray2d< real64, MAX_STENCIL_SIZE *NC > dPresGrad_dC( stencilSize, NC );

    real64 gravHead{};
    stackArray1d< real64, MAX_NUM_ELEMS > dGravHead_dP( numFluxElems );
    stackArray1d< real64, MAX_NUM_ELEMS > dGravHead_dT( numFluxElems );
    stackArray2d< real64, MAX_NUM_ELEMS * NC > dGravHead_dC( numFluxElems, NC );

    real64 dCapPressure_dC[NC]{};

    // Working array
    real64 dProp_dC[NC]{};

    // calculate quantities on primary connected cells
    for( localIndex i = 0; i < numFluxElems; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      // density
      real64 const density  = phaseMassDens[er][esr][ei][0][ip];
      real64 const dDens_dP = dPhaseMassDens_dPres[er][esr][ei][0][ip];
      real64 const dDens_dT = dPhaseMassDens_dTemp[er][esr][ei][0][ip];

      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens_dComp[er][esr][ei][0][ip],
                      dProp_dC );

      // average density and derivatives
      densMean += 0.5 * density;
      dDensMean_dP[i] = 0.5 * dDens_dP;
      dDensMean_dT[i] = 0.5 * dDens_dT;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dDensMean_dC[i][jc] = 0.5 * dProp_dC[jc];
      }
    }

    //***** calculation of flux *****

    // compute potential difference MPFA-style
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      // capillary pressure
      real64 capPressure     = 0.0;
      real64 dCapPressure_dP = 0.0;
      real64 dCapPressure_dT = 0.0;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dCapPressure_dC[ic] = 0.0;
      }

      if( capPressureFlag )
      {
        capPressure = phaseCapPressure[er][esr][ei][0][ip];

        for( localIndex jp = 0; jp < numPhases; ++jp )
        {
          real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac_dPres[er][esr][ei][jp];
          dCapPressure_dT += dCapPressure_dS * dPhaseVolFrac_dTemp[er][esr][ei][jp];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac_dComp[er][esr][ei][jp][jc];
          }
        }
      }

      presGrad += transmissibility[i] * (pres[er][esr][ei] + dPres[er][esr][ei] - capPressure);
      dPresGrad_dP[i] += transmissibility[i] * (1 - dCapPressure_dP) + dTrans_dPres[i] * (pres[er][esr][ei] + dPres[er][esr][ei] - capPressure);
      dPresGrad_dP[i] += transmissibility[i] * (1 - dCapPressure_dT);
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPresGrad_dC[i][jc] += -transmissibility[i] * dCapPressure_dC[jc];
      }

      real64 const gravD     = transmissibility[i] * gravCoef[er][esr][ei];
      real64 const dGravD_dP = dTrans_dPres[i] * gravCoef[er][esr][ei];

      // the density used in the potential difference is always a mass density
      // unlike the density used in the phase mobility, which is a mass density
      // if useMass == 1 and a molar density otherwise
      gravHead += densMean * gravD;

      // need to add contributions from both cells the mean density depends on
      for( localIndex j = 0; j < numFluxElems; ++j )
      {
        dGravHead_dP[j] += dDensMean_dP[j] * gravD + dGravD_dP * densMean;
        dGravHead_dT[j] += dDensMean_dT[j] * gravD;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
        }
      }
    }

    // *** upwinding ***

    // use PPU currently; advanced stuff like IHU would go here
    // TODO isolate into a kernel?

    // compute phase potential gradient
    real64 const potGrad = presGrad - gravHead;

    // choose upstream cell
    localIndex const k_up = (potGrad >= 0) ? 0 : 1;

    localIndex const er_up  = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up  = sei[k_up];

    real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];

    // skip the phase flux if phase not present or immobile upstream
    if( LvArray::math::abs( mobility ) < 1e-20 ) // TODO better constant
    {
      continue;
    }

    // pressure gradient depends on all points in the stencil
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
      dPhaseFlux_dT[ke] += dPresGrad_dT[ke];
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc];
      }

    }

    // gravitational head depends only on the two cells connected (same as mean density)
    for( localIndex ke = 0; ke < numFluxElems; ++ke )
    {
      dPhaseFlux_dP[ke] -= dGravHead_dP[ke];
      dPhaseFlux_dT[ke] -= dGravHead_dT[ke];
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[ke][jc] -= dGravHead_dC[ke][jc];
      }
    }

    // compute the phase flux and derivatives using upstream cell mobility
    phaseFlux = mobility * potGrad;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dPhaseFlux_dP[ke] *= mobility;
      dPhaseFlux_dT[ke] *= mobility;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[ke][jc] *= mobility;
      }
    }

    real64 const dMob_dP = dPhaseMob_dPres[er_up][esr_up][ei_up][ip];
    real64 const dMob_dT = dPhaseMob_dTemp[er_up][esr_up][ei_up][ip];
    arraySlice1d< real64 const, compflow::USD_PHASE_DC - 2 > dPhaseMob_dCompSub = dPhaseMob_dComp[er_up][esr_up][ei_up][ip];

    // add contribution from upstream cell mobility derivatives
    dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
    dPhaseFlux_dT[k_up] += dMob_dT * potGrad;
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dPhaseMob_dCompSub[jc] * potGrad;
    }

    // slice some constitutive arrays to avoid too much indexing in component loop
    arraySlice1d< real64 const, multifluid::USD_PHASE_COMP-3 > phaseCompFracSub = phaseCompFrac[er_up][esr_up][ei_up][0][ip];
    arraySlice1d< real64 const, multifluid::USD_PHASE_COMP-3 > dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er_up][esr_up][ei_up][0][ip];
    arraySlice1d< real64 const, multifluid::USD_PHASE_COMP-3 > dPhaseCompFrac_dTempSub = dPhaseCompFrac_dTemp[er_up][esr_up][ei_up][0][ip];
    arraySlice2d< real64 const, multifluid::USD_PHASE_COMP_DC-3 > dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er_up][esr_up][ei_up][0][ip];

    // compute component fluxes and derivatives using upstream cell composition
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const ycp = phaseCompFracSub[ic];
      compFlux[ic] += phaseFlux * ycp;

      // derivatives stemming from phase flux
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
        dCompFlux_dT[ke][ic] += dPhaseFlux_dT[ke] * ycp;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
        }
      }

      // additional derivatives stemming from upstream cell phase composition
      dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];
      dCompFlux_dT[k_up][ic] += phaseFlux * dPhaseCompFrac_dTempSub[ic];

      // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component
      // densities
      applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFrac_dCompSub[ic], dProp_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompFlux_dC[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
      }
    }

    // compute enthalpy fluxes and derivatives
    real64 const enthalpy = phaseEnthalpy[er_up][esr_up][ei_up][0][ip];
    enthalpyFlux = phaseFlux * enthalpy;

    // derivatives stemming from phase flux (it was upwinded previously)
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dEnthalpyFlux_dP[ke] = dPhaseFlux_dP[ke] * enthalpy;
      dEnthalpyFlux_dT[ke] = dPhaseFlux_dT[ke] * enthalpy;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dEnthalpyFlux_dC[ke][jc] = dPhaseFlux_dC[ke][jc] * enthalpy;
      }
    }

    // add contribution from upstream cell enthalpy derivatives
    dEnthalpyFlux_dP[k_up] += phaseFlux * dPhaseEnthalpy_dPres[er_up][esr_up][ei_up][0][ip];
    dEnthalpyFlux_dT[k_up] += phaseFlux * dPhaseEnthalpy_dTemp[er_up][esr_up][ei_up][0][ip];

    applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseEnthalpy_dComp[er_up][esr_up][ei_up][0][ip], dProp_dC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dEnthalpyFlux_dC[k_up][jc] += phaseFlux * dProp_dC[jc];
    }
  }

  // *** end of upwinding

  // mass/molar flux
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    localFlux[ic]          =  dt * compFlux[ic];
    localFlux[NC + 1 + ic] = -dt * compFlux[ic];

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      localIndex const localDofIndexPres                = ke * NDOF;
      localFluxJacobian[ic][localDofIndexPres]          =  dt * dCompFlux_dP[ke][ic];
      localFluxJacobian[NC + 1 + ic][localDofIndexPres] = -dt * dCompFlux_dP[ke][ic];
      localIndex const localDofIndexTemp                = localDofIndexPres + NC;
      localFluxJacobian[ic][localDofIndexTemp]          =  dt * dCompFlux_dT[ke][ic];
      localFluxJacobian[NC + 1 + ic][localDofIndexTemp] = -dt * dCompFlux_dT[ke][ic];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        localIndex const localDofIndexComp                = localDofIndexPres + jc + 1;
        localFluxJacobian[ic][localDofIndexComp]          = dt * dCompFlux_dC[ke][ic][jc];
        localFluxJacobian[NC + 1 + ic][localDofIndexComp] = -dt * dCompFlux_dC[ke][ic][jc];
      }
    }
  }

  // enthalpy flux
  localFlux[NC]     =  dt * enthalpyFlux;
  localFlux[2 * NC] = -dt * enthalpyFlux;

  for( localIndex ke = 0; ke < stencilSize; ++ke )
  {
    localIndex const localDofIndexPres           = ke * NDOF;
    localFluxJacobian[NC][localDofIndexPres]     =  dt * dEnthalpyFlux_dP[ke];
    localFluxJacobian[2 * NC][localDofIndexPres] = -dt * dEnthalpyFlux_dP[ke];
    localIndex const localDofIndexTemp           = localDofIndexPres + NC;
    localFluxJacobian[NC][localDofIndexTemp]     =  dt * dEnthalpyFlux_dT[ke];
    localFluxJacobian[2 * NC][localDofIndexTemp] = -dt * dEnthalpyFlux_dT[ke];

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      localIndex const localDofIndexComp           = localDofIndexPres + jc + 1;
      localFluxJacobian[NC][localDofIndexComp]     =  dt * dEnthalpyFlux_dC[ke][jc];
      localFluxJacobian[2 * NC][localDofIndexComp] = -dt * dEnthalpyFlux_dC[ke][jc];
    }
  }
}

template< localIndex NC, typename STENCILWRAPPER_TYPE >
void
FluxKernel::
  launch( localIndex const numPhases,
          STENCILWRAPPER_TYPE const & stencilWrapper,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & temp,
          ElementViewConst< arrayView1d< real64 const > > const & dTemp,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dTemp,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dComp,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseVolFrac_dPres,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseVolFrac_dTemp,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac_dComp,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dTemp,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dComp,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dTemp,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dComp,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseEnthalpy,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseEnthalpy_dPres,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseEnthalpy_dTemp,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseEnthalpy_dComp,
          ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const & phaseCapPressure,
          ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
          integer const capPressureFlag,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  using namespace CompositionalMultiphaseUtilities;

  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  constexpr localIndex MAX_NUM_ELEMS = STENCILWRAPPER_TYPE::NUM_POINT_IN_FLUX;
  constexpr localIndex MAX_STENCIL_SIZE = STENCILWRAPPER_TYPE::MAX_STENCIL_SIZE;

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const stencilSize = meshMapUtilities::size1( sei, iconn );
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );
    constexpr localIndex NDOF = NC + 2;

    // working arrays
    stackArray1d< globalIndex, MAX_NUM_ELEMS * NDOF > dofColIndices( stencilSize * NDOF );
    stackArray1d< real64, MAX_NUM_ELEMS * ( NC + 1 ) >  localFlux( numFluxElems * ( NC + 1 ) );
    stackArray2d< real64, MAX_NUM_ELEMS * ( NC + 1 ) * MAX_STENCIL_SIZE * NDOF > localFluxJacobian( numFluxElems * ( NC + 1 ), stencilSize * NDOF );

    // compute transmissibility
    real64 transmissibility[STENCILWRAPPER_TYPE::MAX_NUM_OF_CONNECTIONS][2];
    real64 dTrans_dPres[STENCILWRAPPER_TYPE::MAX_NUM_OF_CONNECTIONS][2];

    stencilWrapper.computeWeights( iconn,
                                   permeability,
                                   dPerm_dPres,
                                   transmissibility,
                                   dTrans_dPres );

    FluxKernel::compute< NC, MAX_NUM_ELEMS, MAX_STENCIL_SIZE >( numPhases,
                                                                stencilSize,
                                                                numFluxElems,
                                                                seri[iconn],
                                                                sesri[iconn],
                                                                sei[iconn],
                                                                transmissibility[0],
                                                                dTrans_dPres[0],
                                                                pres,
                                                                dPres,
                                                                temp,
                                                                dTemp,
                                                                gravCoef,
                                                                phaseMob,
                                                                dPhaseMob_dPres,
                                                                dPhaseMob_dTemp,
                                                                dPhaseMob_dComp,
                                                                dPhaseVolFrac_dPres,
                                                                dPhaseVolFrac_dTemp,
                                                                dPhaseVolFrac_dComp,
                                                                dCompFrac_dCompDens,
                                                                phaseMassDens,
                                                                dPhaseMassDens_dPres,
                                                                dPhaseMassDens_dTemp,
                                                                dPhaseMassDens_dComp,
                                                                phaseCompFrac,
                                                                dPhaseCompFrac_dPres,
                                                                dPhaseCompFrac_dTemp,
                                                                dPhaseCompFrac_dComp,
                                                                phaseEnthalpy,
                                                                dPhaseEnthalpy_dPres,
                                                                dPhaseEnthalpy_dTemp,
                                                                dPhaseEnthalpy_dComp,
                                                                phaseCapPressure,
                                                                dPhaseCapPressure_dPhaseVolFrac,
                                                                capPressureFlag,
                                                                dt,
                                                                localFlux,
                                                                localFluxJacobian );

    // populate dof indices
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      globalIndex const offset = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];

      for( localIndex jdof = 0; jdof < NDOF; ++jdof )
      {
        dofColIndices[i * NDOF + jdof] = offset + jdof;
      }
    }

    // Apply equation/variable change transformation(s)
    stackArray1d< real64, MAX_STENCIL_SIZE * (NDOF-1) > work( stencilSize * (NDOF-1) );
    shiftBlockRowsAheadByOneAndReplaceFirstRowWithColumnSum( NC, (NDOF-1)*stencilSize, numFluxElems, localFluxJacobian, work );
    shiftBlockElementsAheadByOneAndReplaceFirstElementWithSum( NC, numFluxElems, localFlux );

    // Add to residual/jacobian
    for( localIndex i = 0; i < numFluxElems; ++i )
    {
      if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( localMatrix.numRows(), localRow + NC + 1 );

        for( localIndex ic = 0; ic < NC + 1; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow + ic], localFlux[i * ( NC + 1) + ic] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + ic,
                                                                            dofColIndices.data(),
                                                                            localFluxJacobian[i * ( NC + 1 ) + ic].dataIfContiguous(),
                                                                            stencilSize * NDOF );
        }
      }
    }
  } );
}

#define INST_FluxKernel( NC, STENCILWRAPPER_TYPE ) \
  template \
  void FluxKernel:: \
    launch< NC, STENCILWRAPPER_TYPE >( localIndex const numPhases, \
                                       STENCILWRAPPER_TYPE const & stencilWrapper, \
                                       globalIndex const rankOffset, \
                                       ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber, \
                                       ElementViewConst< arrayView1d< integer const > > const & ghostRank, \
                                       ElementViewConst< arrayView1d< real64 const > > const & pres, \
                                       ElementViewConst< arrayView1d< real64 const > > const & dPres, \
                                       ElementViewConst< arrayView1d< real64 const > > const & temp, \
                                       ElementViewConst< arrayView1d< real64 const > > const & dTemp, \
                                       ElementViewConst< arrayView3d< real64 const > > const & permeability, \
                                       ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres, \
                                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef, \
                                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres, \
                                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dTemp, \
                                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dComp, \
                                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseVolFrac_dPres, \
                                       ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseVolFrac_dTemp, \
                                       ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac_dComp, \
                                       ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dTemp, \
                                       ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dComp, \
                                       ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                       ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres, \
                                       ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dTemp, \
                                       ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dComp, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseEnthalpy, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseEnthalpy_dPres, \
                                       ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseEnthalpy_dTemp, \
                                       ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseEnthalpy_dComp, \
                                       ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const & phaseCapPressure, \
                                       ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac, \
                                       integer const capPressureFlag, \
                                       real64 const dt, \
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                       arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 1, CellElementStencilTPFAWrapper );
INST_FluxKernel( 2, CellElementStencilTPFAWrapper );
INST_FluxKernel( 3, CellElementStencilTPFAWrapper );
INST_FluxKernel( 4, CellElementStencilTPFAWrapper );
INST_FluxKernel( 5, CellElementStencilTPFAWrapper );

INST_FluxKernel( 1, SurfaceElementStencilWrapper );
INST_FluxKernel( 2, SurfaceElementStencilWrapper );
INST_FluxKernel( 3, SurfaceElementStencilWrapper );
INST_FluxKernel( 4, SurfaceElementStencilWrapper );
INST_FluxKernel( 5, SurfaceElementStencilWrapper );

INST_FluxKernel( 1, EmbeddedSurfaceToCellStencilWrapper );
INST_FluxKernel( 2, EmbeddedSurfaceToCellStencilWrapper );
INST_FluxKernel( 3, EmbeddedSurfaceToCellStencilWrapper );
INST_FluxKernel( 4, EmbeddedSurfaceToCellStencilWrapper );
INST_FluxKernel( 5, EmbeddedSurfaceToCellStencilWrapper );

INST_FluxKernel( 1, FaceElementToCellStencilWrapper );
INST_FluxKernel( 2, FaceElementToCellStencilWrapper );
INST_FluxKernel( 3, FaceElementToCellStencilWrapper );
INST_FluxKernel( 4, FaceElementToCellStencilWrapper );
INST_FluxKernel( 5, FaceElementToCellStencilWrapper );

#undef INST_FluxKernel



} // namespace ThermalCompositionalMultiphaseFVMKernels

} // namespace geosx
