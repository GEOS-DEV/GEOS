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
 * @file CompositionalMultiphaseFlowKernels.cpp
 */

#include "CompositionalMultiphaseFlowKernels.hpp"

#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/FaceElementStencil.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFlowKernels
{

/******************************** ComponentFractionKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ComponentFractionKernel::
  Compute( arraySlice1d< real64 const > const compDens,
           arraySlice1d< real64 const > const dCompDens,
           arraySlice1d< real64 > const compFrac,
           arraySlice2d< real64 > const dCompFrac_dCompDens )
{
  real64 totalDensity = 0.0;

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    totalDensity += compDens[ic] + dCompDens[ic];
  }

  real64 const totalDensityInv = 1.0 / totalDensity;

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    compFrac[ic] = (compDens[ic] + dCompDens[ic]) * totalDensityInv;
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dCompFrac_dCompDens[ic][jc] = -compFrac[ic] * totalDensityInv;
    }
    dCompFrac_dCompDens[ic][ic] += totalDensityInv;
  }
}

template< localIndex NC >
void
ComponentFractionKernel::
  Launch( localIndex const size,
          arrayView2d< real64 const > const & compDens,
          arrayView2d< real64 const > const & dCompDens,
          arrayView2d< real64 > const & compFrac,
          arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    Compute< NC >( compDens[a],
                   dCompDens[a],
                   compFrac[a],
                   dCompFrac_dCompDens[a] );
  } );
}

template< localIndex NC >
void
ComponentFractionKernel::
  Launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const > const & compDens,
          arrayView2d< real64 const > const & dCompDens,
          arrayView2d< real64 > const & compFrac,
          arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute< NC >( compDens[a],
                   dCompDens[a],
                   compFrac[a],
                   dCompFrac_dCompDens[a] );
  } );
}

#define INST_ComponentFractionKernel( NC ) \
  template \
  void ComponentFractionKernel:: \
    Launch< NC >( localIndex const size, \
                  arrayView2d< real64 const > const & compDens, \
                  arrayView2d< real64 const > const & dCompDens, \
                  arrayView2d< real64 > const & compFrac, \
                  arrayView3d< real64 > const & dCompFrac_dCompDens ); \
  template \
  void ComponentFractionKernel:: \
    Launch< NC >( SortedArrayView< localIndex const > const & targetSet, \
                  arrayView2d< real64 const > const & compDens, \
                  arrayView2d< real64 const > const & dCompDens, \
                  arrayView2d< real64 > const & compFrac, \
                  arrayView3d< real64 > const & dCompFrac_dCompDens )

INST_ComponentFractionKernel( 1 );
INST_ComponentFractionKernel( 2 );
INST_ComponentFractionKernel( 3 );
INST_ComponentFractionKernel( 4 );
INST_ComponentFractionKernel( 5 );

#undef INST_ComponentFractionKernel

/******************************** PhaseVolumeFractionKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
PhaseVolumeFractionKernel::
  Compute( arraySlice1d< real64 const > const & compDens,
           arraySlice1d< real64 const > const & dCompDens,
           arraySlice2d< real64 const > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const > const & phaseDens,
           arraySlice1d< real64 const > const & dPhaseDens_dPres,
           arraySlice2d< real64 const > const & dPhaseDens_dComp,
           arraySlice1d< real64 const > const & phaseFrac,
           arraySlice1d< real64 const > const & dPhaseFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseFrac_dComp,
           arraySlice1d< real64 > const & phaseVolFrac,
           arraySlice1d< real64 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 > const & dPhaseVolFrac_dComp )
{
  real64 work[NC];

  // compute total density from component partial densities
  real64 totalDensity = 0.0;
  real64 const dTotalDens_dCompDens = 1.0;
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    totalDensity += compDens[ic] + dCompDens[ic];
  }

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
    real64 const phaseDensInv = 1.0 / phaseDens[ip];

    // compute saturation and derivatives except multiplying by the total density
    phaseVolFrac[ip] = phaseFrac[ip] * phaseDensInv;

    dPhaseVolFrac_dPres[ip] =
      (dPhaseFrac_dPres[ip] - phaseVolFrac[ip] * dPhaseDens_dPres[ip]) * phaseDensInv;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseVolFrac_dComp[ip][jc] =
        (dPhaseFrac_dComp[ip][jc] - phaseVolFrac[ip] * dPhaseDens_dComp[ip][jc]) * phaseDensInv;
    }

    // apply chain rule to convert derivatives from global component fractions to densities
    applyChainRuleInPlace( NC, dCompFrac_dCompDens, dPhaseVolFrac_dComp[ip], work );

    // now finalize the computation by multiplying by total density
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseVolFrac_dComp[ip][jc] *= totalDensity;
      dPhaseVolFrac_dComp[ip][jc] += phaseVolFrac[ip] * dTotalDens_dCompDens;
    }

    phaseVolFrac[ip] *= totalDensity;
    dPhaseVolFrac_dPres[ip] *= totalDensity;
  }
}

template< localIndex NC, localIndex NP >
void PhaseVolumeFractionKernel::
  Launch( localIndex const size,
          arrayView2d< real64 const > const & compDens,
          arrayView2d< real64 const > const & dCompDens,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseFrac,
          arrayView3d< real64 const > const & dPhaseFrac_dPres,
          arrayView4d< real64 const > const & dPhaseFrac_dComp,
          arrayView2d< real64 > const & phaseVolFrac,
          arrayView2d< real64 > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 > const & dPhaseVolFrac_dComp )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    Compute< NC, NP >( compDens[a],
                       dCompDens[a],
                       dCompFrac_dCompDens[a],
                       phaseDens[a][0],
                       dPhaseDens_dPres[a][0],
                       dPhaseDens_dComp[a][0],
                       phaseFrac[a][0],
                       dPhaseFrac_dPres[a][0],
                       dPhaseFrac_dComp[a][0],
                       phaseVolFrac[a],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dComp[a] );
  } );
}

template< localIndex NC, localIndex NP >
void PhaseVolumeFractionKernel::
  Launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const > const & compDens,
          arrayView2d< real64 const > const & dCompDens,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseFrac,
          arrayView3d< real64 const > const & dPhaseFrac_dPres,
          arrayView4d< real64 const > const & dPhaseFrac_dComp,
          arrayView2d< real64 > const & phaseVolFrac,
          arrayView2d< real64 > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 > const & dPhaseVolFrac_dComp )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute< NC, NP >( compDens[a],
                       dCompDens[a],
                       dCompFrac_dCompDens[a],
                       phaseDens[a][0],
                       dPhaseDens_dPres[a][0],
                       dPhaseDens_dComp[a][0],
                       phaseFrac[a][0],
                       dPhaseFrac_dPres[a][0],
                       dPhaseFrac_dComp[a][0],
                       phaseVolFrac[a],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dComp[a] );
  } );
}

#define INST_PhaseVolumeFractionKernel( NC, NP ) \
  template \
  void \
  PhaseVolumeFractionKernel:: \
    Launch< NC, NP >( localIndex const size, \
                      arrayView2d< real64 const > const & compDens, \
                      arrayView2d< real64 const > const & dCompDens, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseDens, \
                      arrayView3d< real64 const > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const > const & phaseFrac, \
                      arrayView3d< real64 const > const & dPhaseFrac_dPres, \
                      arrayView4d< real64 const > const & dPhaseFrac_dComp, \
                      arrayView2d< real64 > const & phaseVolFrac, \
                      arrayView2d< real64 > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 > const & dPhaseVolFrac_dComp ); \
  template \
  void \
  PhaseVolumeFractionKernel:: \
    Launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                      arrayView2d< real64 const > const & compDens, \
                      arrayView2d< real64 const > const & dCompDens, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseDens, \
                      arrayView3d< real64 const > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const > const & phaseFrac, \
                      arrayView3d< real64 const > const & dPhaseFrac_dPres, \
                      arrayView4d< real64 const > const & dPhaseFrac_dComp, \
                      arrayView2d< real64 > const & phaseVolFrac, \
                      arrayView2d< real64 > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 > const & dPhaseVolFrac_dComp )

INST_PhaseVolumeFractionKernel( 1, 1 );
INST_PhaseVolumeFractionKernel( 2, 1 );
INST_PhaseVolumeFractionKernel( 3, 1 );
INST_PhaseVolumeFractionKernel( 4, 1 );
INST_PhaseVolumeFractionKernel( 5, 1 );

INST_PhaseVolumeFractionKernel( 1, 2 );
INST_PhaseVolumeFractionKernel( 2, 2 );
INST_PhaseVolumeFractionKernel( 3, 2 );
INST_PhaseVolumeFractionKernel( 4, 2 );
INST_PhaseVolumeFractionKernel( 5, 2 );

INST_PhaseVolumeFractionKernel( 1, 3 );
INST_PhaseVolumeFractionKernel( 2, 3 );
INST_PhaseVolumeFractionKernel( 3, 3 );
INST_PhaseVolumeFractionKernel( 4, 3 );
INST_PhaseVolumeFractionKernel( 5, 3 );

#undef INST_PhaseVolumeFractionKernel

/******************************** PhaseMobilityKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
PhaseMobilityKernel::
  Compute( arraySlice2d< real64 const > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const > const & phaseDens,
           arraySlice1d< real64 const > const & dPhaseDens_dPres,
           arraySlice2d< real64 const > const & dPhaseDens_dComp,
           arraySlice1d< real64 const > const & phaseVisc,
           arraySlice1d< real64 const > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const > const & phaseRelPerm,
           arraySlice2d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64 > const & phaseMob,
           arraySlice1d< real64 > const & dPhaseMob_dPres,
           arraySlice2d< real64 > const & dPhaseMob_dComp )
{
  real64 dRelPerm_dC[NC];
  real64 dDens_dC[NC];
  real64 dVisc_dC[NC];

  for( localIndex ip = 0; ip < NP; ++ip )
  {
//    real64 const density_dd = phaseDens[ip];
//    real64 dummy = density_dd;
//    dummy += dPhaseDens_dPres[ip];
//    real64 const density =  1.0;
    real64 const density = phaseDens[ip];

    real64 const dDens_dP = dPhaseDens_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dDens_dC );

    real64 const viscosity = phaseVisc[ip];
    real64 const dVisc_dP = dPhaseVisc_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseVisc_dComp[ip], dVisc_dC );

    real64 const relPerm = phaseRelPerm[ip];
    real64 dRelPerm_dP = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dRelPerm_dC[ic] = 0.0;
    }

    for( localIndex jp = 0; jp < NP; ++jp )
    {
      real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
      dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[jp];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[jp][jc];
      }
    }

    real64 const mobility = relPerm * density / viscosity;

    phaseMob[ip] = mobility;
    dPhaseMob_dPres[ip] = dRelPerm_dP * density / viscosity
                          + mobility * (dDens_dP / density - dVisc_dP / viscosity);

//    dPhaseMob_dPres[ip] = dRelPerm_dP / viscosity
//                          - mobility * ( dVisc_dP / viscosity);
    // compositional derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
//      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] / viscosity
//                                - mobility * ( dVisc_dC[jc] / viscosity);
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] * density / viscosity
                                + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
    }
  }
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::
  Launch( localIndex const size,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    Compute< NC, NP >( dCompFrac_dCompDens[a],
                       phaseDens[a][0],
                       dPhaseDens_dPres[a][0],
                       dPhaseDens_dComp[a][0],
                       phaseVisc[a][0],
                       dPhaseVisc_dPres[a][0],
                       dPhaseVisc_dComp[a][0],
                       phaseRelPerm[a][0],
                       dPhaseRelPerm_dPhaseVolFrac[a][0],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dComp[a],
                       phaseMob[a],
                       dPhaseMob_dPres[a],
                       dPhaseMob_dComp[a] );
  } );
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::
  Launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute< NC, NP >( dCompFrac_dCompDens[a],
                       phaseDens[a][0],
                       dPhaseDens_dPres[a][0],
                       dPhaseDens_dComp[a][0],
                       phaseVisc[a][0],
                       dPhaseVisc_dPres[a][0],
                       dPhaseVisc_dComp[a][0],
                       phaseRelPerm[a][0],
                       dPhaseRelPerm_dPhaseVolFrac[a][0],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dComp[a],
                       phaseMob[a],
                       dPhaseMob_dPres[a],
                       dPhaseMob_dComp[a] );
  } );
}

#define INST_PhaseMobilityKernel( NC, NP ) \
  template \
  void \
  PhaseMobilityKernel:: \
    Launch< NC, NP >( localIndex const size, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseDens, \
                      arrayView3d< real64 const > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const > const & phaseVisc, \
                      arrayView3d< real64 const > const & dPhaseVisc_dPres, \
                      arrayView4d< real64 const > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const > const & phaseRelPerm, \
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64 > const & phaseMob, \
                      arrayView2d< real64 > const & dPhaseMob_dPres, \
                      arrayView3d< real64 > const & dPhaseMob_dComp ); \
  template \
  void \
  PhaseMobilityKernel:: \
    Launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseDens, \
                      arrayView3d< real64 const > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const > const & phaseVisc, \
                      arrayView3d< real64 const > const & dPhaseVisc_dPres, \
                      arrayView4d< real64 const > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const > const & phaseRelPerm, \
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64 > const & phaseMob, \
                      arrayView2d< real64 > const & dPhaseMob_dPres, \
                      arrayView3d< real64 > const & dPhaseMob_dComp )

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

/******************************** AccumulationKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
AccumulationKernel::
  Compute( localIndex const numPhases,
           real64 const & volume,
           real64 const & porosityOld,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice2d< real64 const > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const > const & phaseVolFracOld,
           arraySlice1d< real64 const > const & phaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dCompDens,
           arraySlice1d< real64 const > const & phaseDensOld,
           arraySlice1d< real64 const > const & phaseDens,
           arraySlice1d< real64 const > const & dPhaseDens_dPres,
           arraySlice2d< real64 const > const & dPhaseDens_dComp,
           arraySlice2d< real64 const > const & phaseCompFracOld,
           arraySlice2d< real64 const > const & phaseCompFrac,
           arraySlice2d< real64 const > const & dPhaseCompFrac_dPres,
           arraySlice3d< real64 const > const & dPhaseCompFrac_dComp,
           real64 ( & localAccum )[NC],
           real64 ( & localAccumJacobian )[NC][NC + 1] )
{
  localIndex constexpr NDOF = NC + 1;
  localIndex const NP = numPhases;

  // temporary work arrays
  real64 dPhaseAmount_dC[NC];
  real64 dPhaseCompFrac_dC[NC];

  // reset the local values
  for( localIndex i = 0; i < NC; ++i )
  {
    localAccum[i] = 0.0;
    for( localIndex j = 0; j < NDOF; ++j )
    {
      localAccumJacobian[i][j] = 0.0;
    }
  }

  // compute fluid-independent (pore volume) part
  real64 const volNew = volume;
  real64 const volOld = volume;
  real64 const dVol_dP = 0.0; // used in poroelastic solver

  real64 const poroNew = porosityRef * pvMult;
  real64 const poroOld = porosityOld;
  real64 const dPoro_dP = porosityRef * dPvMult_dPres;

  real64 const poreVolNew = volNew * poroNew;
  real64 const poreVolOld = volOld * poroOld;
  real64 const dPoreVol_dP = dVol_dP * poroNew + volNew * dPoro_dP;

  // sum contributions to component accumulation from each phase
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ip] * phaseDens[ip];
    real64 const phaseAmountOld = poreVolOld * phaseVolFracOld[ip] * phaseDensOld[ip];

    real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFrac[ip] * phaseDens[ip]
                                   + poreVolNew * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                   + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

    // assemble density dependence
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                            + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
      dPhaseAmount_dC[jc] *= poreVolNew;
    }

    // ic - index of component whose conservation equation is assembled
    // (i.e. row number in local matrix)
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
      real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ip][ic];

      real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                         + phaseAmountNew * dPhaseCompFrac_dPres[ip][ic];

      localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
      localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

      // jc - index of component w.r.t. whose compositional var the derivative is being taken
      // (i.e. col number in local matrix)

      // assemble phase composition dependence
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac_dComp[ip][ic], dPhaseCompFrac_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                           + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
        localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
      }
    }
  }
}

template< localIndex NC >
void
AccumulationKernel::
  Launch( localIndex const numPhases,
          localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & porosityOld,
          arrayView1d< real64 const > const & porosityRef,
          arrayView2d< real64 const > const & pvMult,
          arrayView2d< real64 const > const & dPvMult_dPres,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView2d< real64 const > const & phaseVolFracOld,
          arrayView2d< real64 const > const & phaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens,
          arrayView2d< real64 const > const & phaseDensOld,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseCompFracOld,
          arrayView4d< real64 const > const & phaseCompFrac,
          arrayView4d< real64 const > const & dPhaseCompFrac_dPres,
          arrayView5d< real64 const > const & dPhaseCompFrac_dComp,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    if( elemGhostRank[ei] >= 0 )
      return;

    localIndex constexpr NDOF = NC + 1;

    real64 localAccum[NC];
    real64 localAccumJacobian[NC][NDOF];

    Compute< NC >( numPhases,
                   volume[ei],
                   porosityOld[ei],
                   porosityRef[ei],
                   pvMult[ei][0],
                   dPvMult_dPres[ei][0],
                   dCompFrac_dCompDens[ei],
                   phaseVolFracOld[ei],
                   phaseVolFrac[ei],
                   dPhaseVolFrac_dPres[ei],
                   dPhaseVolFrac_dCompDens[ei],
                   phaseDensOld[ei],
                   phaseDens[ei][0],
                   dPhaseDens_dPres[ei][0],
                   dPhaseDens_dComp[ei][0],
                   phaseCompFracOld[ei],
                   phaseCompFrac[ei][0],
                   dPhaseCompFrac_dPres[ei][0],
                   dPhaseCompFrac_dComp[ei][0],
                   localAccum,
                   localAccumJacobian );

    // set DOF indices for this block
    localIndex const localRow = dofNumber[ei] - rankOffset;
    globalIndex dofIndices[NDOF];
    for( localIndex idof = 0; idof < NDOF; ++idof )
    {
      dofIndices[idof] = dofNumber[ei] + idof;
    }

    // TODO: apply equation/variable change transformation(s)

    // add contribution to residual and jacobian
    for( localIndex i = 0; i < NC; ++i )
    {
      localRhs[localRow + i] += localAccum[i];
      localMatrix.addToRow< serialAtomic >( localRow + i,
                                            dofIndices,
                                            localAccumJacobian[i],
                                            NDOF );
    }
  } );
}

#define INST_AccumulationKernel( NC ) \
  template \
  void \
  AccumulationKernel:: \
    Launch< NC >( localIndex const numPhases, \
                  localIndex const size, \
                  globalIndex const rankOffset, \
                  arrayView1d< globalIndex const > const & dofNumber, \
                  arrayView1d< integer const > const & elemGhostRank, \
                  arrayView1d< real64 const > const & volume, \
                  arrayView1d< real64 const > const & porosityOld, \
                  arrayView1d< real64 const > const & porosityRef, \
                  arrayView2d< real64 const > const & pvMult, \
                  arrayView2d< real64 const > const & dPvMult_dPres, \
                  arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                  arrayView2d< real64 const > const & phaseVolFracOld, \
                  arrayView2d< real64 const > const & phaseVolFrac, \
                  arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                  arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens, \
                  arrayView2d< real64 const > const & phaseDensOld, \
                  arrayView3d< real64 const > const & phaseDens, \
                  arrayView3d< real64 const > const & dPhaseDens_dPres, \
                  arrayView4d< real64 const > const & dPhaseDens_dComp, \
                  arrayView3d< real64 const > const & phaseCompFracOld, \
                  arrayView4d< real64 const > const & phaseCompFrac, \
                  arrayView4d< real64 const > const & dPhaseCompFrac_dPres, \
                  arrayView5d< real64 const > const & dPhaseCompFrac_dComp, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_AccumulationKernel( 1 );
INST_AccumulationKernel( 2 );
INST_AccumulationKernel( 3 );
INST_AccumulationKernel( 4 );
INST_AccumulationKernel( 5 );

#undef INST_AccumulationKernel

/******************************** VolumeBalanceKernel ********************************/

template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL, bool IS_UT_FORM >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
FluxKernel::
Compute( localIndex const numPhases,
           localIndex const stencilSize,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           arraySlice1d< real64 const > const stencilWeights,
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
           ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
           ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
           ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp,
           ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
           ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
           ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
           ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
           integer const capPressureFlag,
           real64 const dt,
           arraySlice1d< real64 > const localFlux,
           arraySlice2d< real64 > const localFluxJacobian )
{

  localIndex constexpr NDOF = NC + 1;
  localIndex const NP = numPhases;

  bool const is_ppu = false;
  bool const is_pu = true;//!is_ppu
  bool const is_hu = false;//!(is_pu || is_ppu);
  real64 phase_eps = 0;

  real64 compFlux[NC]{};
  real64 dCompFlux_dP[MAX_STENCIL][NC]{};
  real64 dCompFlux_dC[MAX_STENCIL][NC][NC]{};

  //real64 totFlux{};
  real64 totFlux_unw{};
  real64 dTotFlux_dP[MAX_STENCIL]{};
  real64 dTotFlux_dC[MAX_STENCIL][NC]{};

  //useful lambdas
  // helpers lambda definition to avoid too much dulpication
  // agnostic of upwind direction (unless want to introduce fancy density face def)
  auto dGravHead_dX =
    [&seri,&sesri,&sei,&phaseMassDens,&dPhaseMassDens_dPres,&dPhaseMassDens_dComp,&dCompFrac_dCompDens,&stencilSize,&stencilWeights,&gravCoef]
    (real64& GH, real64 dGH_dP[NUM_ELEMS], real64 dGH_dC[NUM_ELEMS][NC], real64 dPr_dC[NC], int kp)
    {
      real64 densMean {};
      real64 dDensMean_dP [NUM_ELEMS]{};
      real64 dDensMean_dC [NUM_ELEMS][NC]{};


      //init
      GH = 0.0;
      for( localIndex i = 0; i < NUM_ELEMS; ++i )
      {
        dGH_dP[i] = 0.0;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dGH_dC[i][jc] = 0.0;
          dPr_dC[jc] = 0.0;
        }
      }


      for( localIndex i = 0; i < NUM_ELEMS; ++i )
      {
        localIndex const er  = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei  = sei[i];

        // density
        real64 const density  = phaseMassDens[er][esr][ei][0][kp];
        real64 const dDens_dP = dPhaseMassDens_dPres[er][esr][ei][0][kp];

        applyChainRule( NC,
                        dCompFrac_dCompDens[er][esr][ei],
                        dPhaseMassDens_dComp[er][esr][ei][0][kp],
                        dPr_dC );

        // average density and derivatives
        densMean += 0.5 * density;
        dDensMean_dP[i] = 0.5 * dDens_dP;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dDensMean_dC[i][jc] = 0.5 * dPr_dC[jc];
        }
      }

      // compute potential difference MPFA-style
      for( localIndex i = 0; i < stencilSize; ++i )
      {
        localIndex const er  = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei  = sei[i];
        real64 const weight  = stencilWeights[i];

        real64 const gravD = weight * gravCoef[er][esr][ei];
        GH += densMean * gravD;

        // need to add contributions from both cells the mean density depends on
        for( localIndex j = 0; j < NUM_ELEMS; ++j )
        {
          dGH_dP[j] += dDensMean_dP[j] * gravD;
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dGH_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
          }
        }
      }
    };

  //rescale mobilities and their derivatives without molar densities as prefactor
  auto dMob_dX =
    [&seri,&sesri,&sei,&phaseMob,&dPhaseMob_dPres,&dPhaseMob_dComp,&phaseDens,&dPhaseDens_dPres,&dPhaseDens_dComp,&dCompFrac_dCompDens,&stencilSize]
    (real64 mob[MAX_STENCIL], real64 dMob_dP[MAX_STENCIL], real64 dMob_dC[MAX_STENCIL][NC], int kp)
    {

   //init
      for( localIndex j = 0; j < stencilSize; ++j )
      {
        mob[j] = 0.0;
        dMob_dP[j] = 0.0;

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dMob_dC[j][jc] = 0.0;
        }
      }


      for( localIndex j = 0; j < stencilSize; ++j )
      {
        localIndex er_i_  = seri[j];
        localIndex esr_i_ = sesri[j];
        localIndex ei_i_  = sei[j];

        real64 const mobility = phaseMob[er_i_][esr_i_][ei_i_][kp];
        real64 const PD = phaseDens[er_i_][esr_i_][ei_i_][0][kp];
        if( std::fabs(PD) < 1e-20)
          continue;

        mob[j] = mobility / PD;
//        mob[j] = mobility;

        real64 dPD_dP = dPhaseDens_dPres[er_i_][esr_i_][ei_i_][0][kp];
        real64 dPD_dC[NC] {};
        applyChainRule( NC, dCompFrac_dCompDens[er_i_][esr_i_][ei_i_], dPhaseDens_dComp[er_i_][esr_i_][ei_i_][0][kp], dPD_dC );

        real64 const dMMob_dP = dPhaseMob_dPres[er_i_][esr_i_][ei_i_][kp];
        arraySlice1d< real64 const > dMMob_dC = dPhaseMob_dComp[er_i_][esr_i_][ei_i_][kp];

        dMob_dP[j] = ( dMMob_dP - dPD_dP * mob[j] ) / PD;
//        dMob_dP[j] =  dMMob_dP;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dMob_dC[j][jc] = ( dMMob_dC[jc] - dPD_dC[jc] * mob[j] ) / PD;
//          dMob_dC[j][jc] = dMMob_dC[jc];
        }
      }

    };

  //re-multiply by phaseDens and derivatives wrt to upwind direction
  auto densMult =
    [&seri,&sesri,&sei,&phaseDens,&dPhaseDens_dPres,&dPhaseDens_dComp,&dCompFrac_dCompDens,&stencilSize]
    (real64& pF, real64 dPF_dP[MAX_STENCIL], real64 dPF_dC[MAX_STENCIL][NC], int kp, int k_up_)
    {

/*          pF *=1;
          dPF_dP[k_up_] *=1;
          dPF_dC[k_up_][0] *=1;
          kp +=0;*/

      localIndex const er_up_  = seri[k_up_];
      localIndex const esr_up_ = sesri[k_up_];
      localIndex const ei_up_  = sei[k_up_];

      for( localIndex ks = 0; ks < stencilSize; ++ks )
      {
        dPF_dP[ks] *= phaseDens[er_up_][esr_up_][ei_up_][0][kp];
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPF_dC[ks][jc] *= phaseDens[er_up_][esr_up_][ei_up_][0][kp];
        }
      }
        dPF_dP[k_up_] += dPhaseDens_dPres[er_up_][esr_up_][ei_up_][0][kp] * pF;
        real64 dPD_dC[NC]{};
        applyChainRule( NC, dCompFrac_dCompDens[er_up_][esr_up_][ei_up_], dPhaseDens_dComp[er_up_][esr_up_][ei_up_][0][kp],
                        dPD_dC );

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPF_dC[k_up_][jc] += dPD_dC[jc] * pF;
        }

      //last as multiplicative use in the second part of derivatives
      pF *= phaseDens[er_up_][esr_up_][ei_up_][0][kp];

    };

  // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
  for( localIndex ip = 0; ip < NP; ++ip )
  {

    // create local work arrays
    real64 phaseFlux{};
    real64 dPhaseFlux_dP[MAX_STENCIL]{};
    real64 dPhaseFlux_dC[MAX_STENCIL][NC]{};

    real64 presGrad{};
    real64 dPresGrad_dP[MAX_STENCIL]{};
    real64 dPresGrad_dC[MAX_STENCIL][NC]{};

    real64 gravHead{};
    real64 dGravHead_dP[NUM_ELEMS]{};
    real64 dGravHead_dC[NUM_ELEMS][NC]{};

    real64 dCapPressure_dC[NC]{};

    // Working array
    real64 dProp_dC[NC]{};

    // calculate quantities on primary connected cells
    dGravHead_dX(gravHead, dGravHead_dP, dGravHead_dC, dProp_dC, ip);

    //***** calculation of flux *****
    // compute potential difference MPFA-style
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];
      real64 const weight  = stencilWeights[i];

      // capillary pressure
      real64 capPressure     = 0.0;
      real64 dCapPressure_dP = 0.0;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dCapPressure_dC[ic] = 0.0;
      }

      if( capPressureFlag )
      {
        capPressure = phaseCapPressure[er][esr][ei][0][ip];

        for( localIndex jp = 0; jp < NP; ++jp )
        {
          real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac_dPres[er][esr][ei][jp];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac_dComp[er][esr][ei][jp][jc];
          }
        }
      }

      presGrad += weight * (pres[er][esr][ei] + dPres[er][esr][ei] - capPressure);
      dPresGrad_dP[i] += weight * (1 - dCapPressure_dP);
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPresGrad_dC[i][jc] += -weight * dCapPressure_dC[jc];
      }
    }

    // *** upwinding ***

    // use PPU currently; advanced stuff like IHU would go here
    // TODO isolate into a kernel?

    // compute phase potential gradient
    real64 const potGrad = presGrad - gravHead;

    // choose upstream cell
    localIndex const k_up = (potGrad >= 0) ? 0 : 1;

//    std::cerr << " " << k_up <<  " " ;

    localIndex er_up  = seri[k_up];
    localIndex esr_up = sesri[k_up];
    localIndex ei_up  = sei[k_up];

    /**** tentative correction of mobility ****/
    real64 mob_unw[MAX_STENCIL] {};
    real64 dMob_unw_dP[MAX_STENCIL] {};
    real64 dMob_unw_dC[MAX_STENCIL][NC] {};

    dMob_dX(mob_unw,dMob_unw_dP,dMob_unw_dC,ip);

    //adding totMob update for UT -  formulation
    // skip the phase flux if phase not present or immobile upstream
    if( std::fabs( mob_unw[k_up] ) < 1e-20 ) // TODO better constant
    {
      continue;
    }

    // pressure gradient depends on all points in the stencil
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc];
      }
    }

    // gravitational head depends only on the two cells connected (same as mean density)
    for( localIndex ke = 0; ke < NUM_ELEMS; ++ke )
    {
      dPhaseFlux_dP[ke] -= dGravHead_dP[ke];
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[ke][jc] -= dGravHead_dC[ke][jc];
      }
    }

    // compute the phase flux and derivatives using upstream cell mobility
    phaseFlux = mob_unw[k_up] * potGrad;
    phase_eps = (phase_eps > 0.01*std::fabs(phaseFlux)) ? phase_eps : 0.01*std::fabs(phaseFlux);
    //adding totflux for UT formulation
    totFlux_unw += phaseFlux;

      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dPhaseFlux_dP[ke] *= mob_unw[k_up];
        dTotFlux_dP[ke] += dPhaseFlux_dP[ke];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPhaseFlux_dC[ke][jc] *= mob_unw[k_up];
          dTotFlux_dC[ke][jc] += dPhaseFlux_dC[ke][jc];
        }
      }

      // add contribution from upstream cell mobility derivatives
    dPhaseFlux_dP[k_up] += dMob_unw_dP[k_up] * potGrad;
    dTotFlux_dP[k_up] += dMob_unw_dP[k_up] * potGrad;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dMob_unw_dC[k_up][jc] * potGrad;
      dTotFlux_dC[k_up][jc] += dMob_unw_dC[k_up][jc] * potGrad;
    }

    //validating densMult helper
    densMult(phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, ip, k_up);

//    totFlux += phaseFlux;

//    std::cerr << " " << dPhaseFlux_dC[ip][0] << " " << dPhaseFlux_dC[ip][1] << " " << dPhaseFlux_dC[ip][2] << " ";
//    std::cerr << " " << dPhaseFlux_dP[0] << " " << dPhaseFlux_dP[1] << " ";
//    std::cerr << " " << phaseFlux <<  " ";


    if( !IS_UT_FORM ) // skip  if you intend to use fixed total valocity formulation
    {

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d< real64 const > phaseCompFracSub = phaseCompFrac[er_up][esr_up][ei_up][0][ip];
      arraySlice1d< real64 const > dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er_up][esr_up][ei_up][0][ip];
      arraySlice2d< real64 const > dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er_up][esr_up][ei_up][0][ip];

      // compute component fluxes and derivatives using upstream cell composition
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        real64 const ycp = phaseCompFracSub[ic];
        compFlux[ic] += phaseFlux * ycp;

        // derivatives stemming from phase flux
        for( localIndex ke = 0; ke < stencilSize; ++ke )
        {
          dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
          }
        }

        // additional derivatives stemming from upstream cell phase composition
        dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];

        // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component
        // densities
        applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFrac_dCompSub[ic], dProp_dC );
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFlux_dC[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
        }
      }
    }
    }

//  std::cerr << " phaseFlux " << totFlux_unw << std::endl;
  // *** end of upwinding

  //if total flux formulation
  if( IS_UT_FORM )
  {
      //for UT form
      real64 totMob_unw{};
      real64 dTotMob_unw_dP[MAX_STENCIL]{};
      real64 dTotMob_unw_dC[MAX_STENCIL][NC]{};

      real64 presGrad {};
      real64 dPresGrad_dP[MAX_STENCIL] {};
      real64 dPresGrad_dC[MAX_STENCIL][NC] {};

    //getting frac flow and derivatives populated
    auto dFrac_dX =
      [&seri, &sesri, &sei, &stencilSize, &totMob_unw, &dTotMob_unw_dP, &dTotMob_unw_dC, &dMob_dX]
        (real64& fflow, real64 dFflow_dP[MAX_STENCIL], real64 dFflow_dC[MAX_STENCIL][NC], int kp, int k_up_)
      {

        real64 mob_unw[MAX_STENCIL] {};
        real64 dMob_unw_dP[MAX_STENCIL] {};
        real64 dMob_unw_dC[MAX_STENCIL][NC] {};


        //reinit
        //fractional flow too low to let the upstream phase flow
        fflow = 0;
        for( localIndex ke = 0; ke < stencilSize; ++ke )
        {
          dFflow_dP[ke] = 0;

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dFflow_dC[ke][jc] = 0;
          }
        }

        dMob_dX(mob_unw,dMob_unw_dP,dMob_unw_dC,kp);
        //compute
        if( std::fabs(mob_unw[k_up_]) > 1e-20 )
        {
          fflow = mob_unw[k_up_] / totMob_unw;
          //update component fluxes
          dFflow_dP[k_up_] = dMob_unw_dP[k_up_] / totMob_unw;
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dFflow_dC[k_up_][jc] = dMob_unw_dC[k_up_][jc] / totMob_unw;
          }

          for( localIndex ke = 0; ke < stencilSize; ++ke )
          {
            dFflow_dP[ke] -= fflow * dTotMob_unw_dP[ke] / totMob_unw;

            for( localIndex jc = 0; jc < NC; ++jc )
            {
              dFflow_dC[ke][jc] -= fflow * dTotMob_unw_dC[ke][jc] / totMob_unw;
            }
          }
        }
      };

    //compute potential
    if(is_ppu || is_pu || is_hu)
      {
        for( localIndex i = 0; i < stencilSize; ++i )
        {
          localIndex const er = seri[i];
          localIndex const esr = sesri[i];
          localIndex const ei = sei[i];
          real64 const weight = stencilWeights[i];
          presGrad += weight * ( pres[er][esr][ei] + dPres[er][esr][ei] );
          dPresGrad_dP[i] += weight * (1);// - dCapPressure_dP);
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dPresGrad_dC[i][jc] += 0;//-weight * dCapPressure_dC[jc];
          }
        }
      }

    //finding the largest negPot to be compared against later
    real64 minGravHead = 0;
    if( is_pu || is_hu ){
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        real64 pot_{};

        real64 gravHead{};
        real64 dGravHead_dP[NUM_ELEMS]{};
        real64 dGravHead_dC[NUM_ELEMS][NC]{};

        real64 dProp_dC[NC]{};
        dGravHead_dX( gravHead, dGravHead_dP, dGravHead_dC, dProp_dC, ip );
        //defining up and down-wind direction based on total flux
        if(is_pu)
          pot_ = totFlux_unw;

        localIndex const k_up = ( pot_ >= phase_eps ) ? 0 : 1;
        localIndex const k_dw = ( pot_ >= phase_eps ) ? 1 : 0;

//        std::cerr << ip << " : (" << pot_ << " ; ";
        for( localIndex jp = 0; jp < NP; ++jp )
        {

          if(ip!=jp)
          {
            real64 gravHeadOther{};
            real64 dGravHeadOther_dP[NUM_ELEMS]{};
            real64 dGravHeadOther_dC[NUM_ELEMS][NC]{};
            real64 dPropOther_dC[NC]{};

            dGravHead_dX( gravHeadOther, dGravHeadOther_dP, dGravHeadOther_dC, dPropOther_dC, jp );
            real64 mob_unw[MAX_STENCIL]{};
            real64 dMob_unw_dP[MAX_STENCIL]{};
            real64 dMob_unw_dC[MAX_STENCIL][NC]{};

            dMob_dX( mob_unw, dMob_unw_dP, dMob_unw_dC, jp );

            real64 const mob_up = mob_unw[k_up];
            real64 const mob_dw = mob_unw[k_dw];

            pot_ += ( gravHead - gravHeadOther >= 0 ) ? mob_dw * ( gravHead - gravHeadOther ) : mob_up * ( gravHead
                                                                                                           - gravHeadOther );
            //          std::cerr << pot_ << "[ " << mob_dw << " ," << mob_up << "] ";

          }
        }

//        std::cerr << "  " << pot_ << "), " << gravHead << std::endl;
        if( pot_ < 0 && std::fabs(gravHead) >= std::fabs(minGravHead) ) // if neg pot and more dense
        {
          minGravHead = gravHead;
        }
      }
    }

        auto getUpwindV = [&presGrad, &dGravHead_dX, &phase_eps, &minGravHead, &totFlux_unw] (bool is_ppu_, bool is_pu_, bool is_hu_, localIndex kp, real64 gravTerm = 0)
          {
            real64 gravHead{};
            real64 dGravHead_dP[NUM_ELEMS]{};
            real64 dGravHead_dC[NUM_ELEMS][NC]{};

            real64 dProp_dC[NC]{};

            dGravHead_dX( gravHead, dGravHead_dP, dGravHead_dC, dProp_dC, kp );

            localIndex k_up = -1;
            if( is_ppu_ )
            {
              k_up = ( presGrad - gravHead >= 0 ) ? 0 : 1;
            }
            else if( is_pu_ )
            {
              k_up = ( totFlux_unw >= phase_eps ) ? 0 : 1;
              if( std::fabs( gravHead ) <= std::fabs( minGravHead ) && std::fabs( minGravHead ) > 0 )
                k_up = ( k_up == 1 ) ? 0 : 1;
            }
            else if( is_hu_ )
            {
              /* nothing here yet */
              k_up = ( totFlux_unw + gravTerm >= phase_eps ) ? 0 : 1;
            }

              return k_up;
          };

        auto getUpwindG = [&presGrad, &dGravHead_dX, &phase_eps, &minGravHead, &totFlux_unw] (bool is_ppu_, bool is_pu_, bool is_hu_, localIndex kp, real64 gravTerm = 0)
        {
          real64 gravHeadOther{};
          real64 dGravHeadOther_dP[NUM_ELEMS]{};
          real64 dGravHeadOther_dC[NUM_ELEMS][NC]{};

          real64 dPropOther_dC[NC]{};

          dGravHead_dX( gravHeadOther, dGravHeadOther_dP, dGravHeadOther_dC, dPropOther_dC, kp );

          localIndex k_up_g = -1;
          if( is_ppu_ )
          {
            k_up_g = ( presGrad - gravHeadOther >= 0 ) ? 0 : 1;
          }
          else if( is_pu_ )
          {
            k_up_g = ( totFlux_unw >= phase_eps ) ? 0 : 1;
            if( std::fabs( gravHeadOther ) <= std::fabs( minGravHead ) && std::fabs( minGravHead ) > 0 )
              k_up_g = ( k_up_g == 1 ) ? 0 : 1;
          }
          else if( is_hu_ )
          {
            //           k_up_g = ( std::fabs(gravHead) > std::fabs(minGravHead) ) ? 0 : 1;
            k_up_g = ( gravTerm > 0 ) ? 0 : 1;
            if( (std::fabs(gravHeadOther) <= std::fabs(minGravHead)) && std::fabs(minGravHead)>0 )
              k_up_g = (k_up_g == 1) ? 0 : 1; // downwind
          }

          return k_up_g;
        };

        auto dTotMob_dX = [&NP, &dMob_dX, &stencilSize, &getUpwindV, &getUpwindG]
        ( bool is_grav_ , bool is_ppu_, bool is_pu_, bool is_hu_, real64 gravTerm =0,
                              real64& totMob_unw_, real64 dTotMob_unw_dP_[MAX_STENCIL], real64 dTotMob_unw_dC_[MAX_STENCIL][NC])
        {

          totMob_unw_ = 0.0;
          for(localIndex ke=0; ke < stencilSize; ++ke)
          {
              dTotMob_unw_dP_[ke] = 0.0;
              for(localIndex ic = 0 ; ic < NC; ++ic)
                dTotMob_unw_dC_[ke][ic] = 0.0;

          };


          for( localIndex ip = 0; ip < NP; ++ip )
          {


            real64 mob_unw[MAX_STENCIL]{};
            real64 dMob_unw_dP[MAX_STENCIL]{};
            real64 dMob_unw_dC[MAX_STENCIL][NC]{};

            dMob_dX( mob_unw, dMob_unw_dP, dMob_unw_dC, ip );
            localIndex k_up = -1;
            if( ! is_grav_ )
              k_up = getUpwindV(is_ppu_, is_pu_, is_hu_, ip);
            else
              k_up = getUpwindG(is_ppu_, is_pu_, is_hu_, ip, gravTerm);

            if( std::fabs( mob_unw[k_up] ) < 1e-20 )
              continue;

            totMob_unw_ += mob_unw[k_up];
            dTotMob_unw_dP_[k_up] += dMob_unw_dP[k_up];

            for( localIndex ic = 0; ic < NC; ++ic )
            {
              dTotMob_unw_dC_[k_up][ic] += dMob_unw_dC[k_up][ic];
            }
          }
        };

        auto dPhaseComp_dX = [&seri, &sesri, &sei, &stencilSize, &phaseCompFrac, &dPhaseCompFrac_dPres, &dPhaseCompFrac_dComp, &dCompFrac_dCompDens]
          (real64 phaseFlux_, real64 dPhaseFlux_dP_[MAX_STENCIL], real64 dPhaseFlux_dC_[MAX_STENCIL][NC],
          real64 compFlux_[NC], real64 dCompFlux_dP_[MAX_STENCIL][NC], real64 dCompFlux_dC_[MAX_STENCIL][NC][NC],
            localIndex kp, localIndex k_up_)
    {
      /*update phaseComp from grav part*/
      localIndex const er_up_ = seri[k_up_];
      localIndex const esr_up_ = sesri[k_up_];
      localIndex const ei_up_ = sei[k_up_];

      arraySlice1d< real64 const > phaseCompFracSub = phaseCompFrac[er_up_][esr_up_][ei_up_][0][kp];
      arraySlice1d< real64 const > dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er_up_][esr_up_][ei_up_][0][kp];
      arraySlice2d< real64 const > dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er_up_][esr_up_][ei_up_][0][kp];

      real64 dProp_dC[NC] {};

      // compute component fluxes and derivatives using upstream cell composition
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        real64 const ycp = phaseCompFracSub[ic];
        compFlux_[ic] += phaseFlux_ * ycp;

        // derivatives stemming from phase flux
        for( localIndex ke = 0; ke < stencilSize; ++ke )
        {
          dCompFlux_dP_[ke][ic] += dPhaseFlux_dP_[ke] * ycp;
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCompFlux_dC_[ke][ic][jc] += dPhaseFlux_dC_[ke][jc] * ycp;
          }
        }

        // additional derivatives stemming from upstream cell phase composition
        dCompFlux_dP_[k_up_][ic] += phaseFlux_ * dPhaseCompFrac_dPresSub[ic];

        // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component
        // densities
        applyChainRule( NC, dCompFrac_dCompDens[er_up_][esr_up_][ei_up_], dPhaseCompFrac_dCompSub[ic], dProp_dC );
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFlux_dC_[k_up_][ic][jc] += phaseFlux_ * dProp_dC[jc];
        }
      }
    };

    for( localIndex ip = 0; ip < NP; ++ip )
    {
      // choose upstream cell
      // create local work arrays
      real64 phaseFlux{};
      real64 dPhaseFlux_dP[MAX_STENCIL]{};
      real64 dPhaseFlux_dC[MAX_STENCIL][NC]{};

      real64 phaseFluxV{};
      real64 dPhaseFluxV_dP[MAX_STENCIL]{};
      real64 dPhaseFluxV_dC[MAX_STENCIL][NC]{};

      real64 fflow{};
      real64 dFflow_dP[MAX_STENCIL]{};
      real64 dFflow_dC[MAX_STENCIL][NC]{};

      real64 gravHead{};
      real64 dGravHead_dP[NUM_ELEMS]{};
      real64 dGravHead_dC[NUM_ELEMS][NC]{};

      real64 dProp_dC[NC]{};

      dGravHead_dX( gravHead, dGravHead_dP, dGravHead_dC, dProp_dC, ip );
      /* chosing upwind for viscous term */
      localIndex k_up = getUpwindV(is_ppu,is_pu,is_hu,ip);
      dTotMob_dX(false, is_ppu, is_pu, is_hu, 0.0, totMob_unw, dTotMob_unw_dP, dTotMob_unw_dC);

//      std::cerr << "" << k_up << " ";

      real64 mob_unw[MAX_STENCIL] {};
      real64 dMob_unw_dP[MAX_STENCIL] {};
      real64 dMob_unw_dC[MAX_STENCIL][NC] {};

      dMob_dX(mob_unw,dMob_unw_dP,dMob_unw_dC,ip);

          //fractional flow too low to let the upstream phase flow
      if( std::fabs(mob_unw[k_up]) < 1e-20 || std::fabs(totMob_unw) < 1e-20 )
          continue;

      //get the fracflow for viscous part and (opt1) re-multiply by molar phaseDens
      dFrac_dX(fflow,dFflow_dP,dFflow_dC, ip, k_up);
      densMult(fflow, dFflow_dP, dFflow_dC, ip, k_up);

      phaseFluxV = fflow * totFlux_unw;

      for( localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFluxV_dP[ke] += dFflow_dP[ke] * totFlux_unw;

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPhaseFluxV_dC[ke][jc] += dFflow_dC[ke][jc] * totFlux_unw;
        }
      }

      //NON-FIXED UT -- to be canceled out if considered fixed
      for( localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFluxV_dP[ke] += fflow*dTotFlux_dP[ke];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPhaseFluxV_dC[ke][jc] += fflow*dTotFlux_dC[ke][jc];
        }
      }

      phaseFlux += phaseFluxV;
      for(localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFlux_dP[ke] += dPhaseFluxV_dP[ke];
        for( localIndex ic = 0; ic < NC; ++ic )
          dPhaseFlux_dC[ke][ic] += dPhaseFluxV_dC[ke][ic];
      }

      /***           GRAVITY TERM                ***/
      //if HU compute gravTerm to get init direction
      real64 gravTerm {};
      for( localIndex jp = 0; jp < NP; ++jp )
      {
        if( ip != jp )
        {
          real64 gravHeadOther{};
          real64 dGravHeadOther_dP[NUM_ELEMS]{};
          real64 dGravHeadOther_dC[NUM_ELEMS][NC]{};
          real64 dPropOther_dC[NC]{};

          dGravHead_dX( gravHeadOther, dGravHeadOther_dP, dGravHeadOther_dC, dPropOther_dC, jp );

          gravTerm -= gravHead - gravHeadOther;
        }
      }

      dTotMob_dX(true, is_ppu, is_pu, is_hu, gravTerm, totMob_unw, dTotMob_unw_dP, dTotMob_unw_dC);
      localIndex k_up_g = getUpwindG(is_ppu, is_pu, is_hu, ip, gravTerm);
      dFrac_dX(fflow, dFflow_dP, dFflow_dC, ip, k_up_g);
      densMult(fflow, dFflow_dP, dFflow_dC, ip, k_up_g);
      for( localIndex jp = 0; jp < NP; ++jp )
      {
        if( ip != jp )
        {
          real64 phaseFluxG{};
          real64 dPhaseFluxG_dP[MAX_STENCIL]{};
          real64 dPhaseFluxG_dC[MAX_STENCIL][NC]{};

          real64 gravHeadOther{};
          real64 dGravHeadOther_dP[NUM_ELEMS]{};
          real64 dGravHeadOther_dC[NUM_ELEMS][NC]{};
          real64 dPropOther_dC[NC]{};

          dGravHead_dX( gravHeadOther, dGravHeadOther_dP, dGravHeadOther_dC, dPropOther_dC, jp );

         //mobOther is upwinded as PPU for consistency with totMob in fractional flow fflow
          real64 gravTermOther {};
          for( localIndex kp = 0; kp < NP; ++kp )
          {
            if( kp != jp )
            {
              real64 gravHeadAlt{};
              real64 dGravHeadAlt_dP[NUM_ELEMS]{};
              real64 dGravHeadAlt_dC[NUM_ELEMS][NC]{};
              real64 dPropAlt_dC[NC]{};

              dGravHead_dX( gravHeadAlt, dGravHeadAlt_dP, dGravHeadAlt_dC, dPropAlt_dC, kp );

              gravTermOther -= gravHeadOther - gravHeadAlt;
            }
          }
         localIndex k_up_og = getUpwindG(is_ppu, is_pu, is_hu, jp, gravTermOther);
//          std::cerr << " " << k_up_g << " ";

          real64 mob_other_unw[MAX_STENCIL] {};
          real64 dMob_other_unw_dP[MAX_STENCIL] {};
          real64 dMob_other_unw_dC[MAX_STENCIL][NC] {};

          dMob_dX(mob_other_unw, dMob_other_unw_dP, dMob_other_unw_dC, jp);
          if( std::fabs(mob_other_unw[k_up_og]) < 1e-20 )
            continue;
          phaseFluxG -= fflow * mob_other_unw[k_up_og] * ( gravHead - gravHeadOther );

          dPhaseFluxG_dP[k_up_og] -= fflow * dMob_other_unw_dP[k_up_og] * ( gravHead - gravHeadOther );
          for( localIndex jc = 0; jc < NC; ++jc )
            dPhaseFluxG_dC[k_up_og][jc] -= fflow * dMob_other_unw_dC[k_up_og][jc] * ( gravHead - gravHeadOther );

          //mob related part of dFflow_dP is only upstream defined but totMob related is defined everywhere
          for( localIndex ke = 0; ke < stencilSize; ++ke )
          {
            dPhaseFluxG_dP[ke] -= dFflow_dP[ke] * mob_other_unw[k_up_og] * ( gravHead - gravHeadOther );

            for( localIndex jc = 0; jc < NC; ++jc )
            {
              dPhaseFluxG_dC[ke][jc] -= dFflow_dC[ke][jc] * mob_other_unw[k_up_og] * ( gravHead - gravHeadOther );
            }
          }

          for( localIndex ke = 0; ke < NUM_ELEMS; ++ke )
          {
            dPhaseFluxG_dP[ke] -= fflow * mob_other_unw[k_up_og] * ( dGravHead_dP[ke] - dGravHeadOther_dP[ke] );
            for( localIndex jc = 0; jc < NC; ++jc )
            {
              dPhaseFluxG_dC[ke][jc] -= fflow * mob_other_unw[k_up_og] * ( dGravHead_dC[ke][jc] - dGravHeadOther_dC[ke][jc] );
            }
          }
          dPhaseComp_dX(phaseFluxG, dPhaseFluxG_dP, dPhaseFluxG_dC, compFlux, dCompFlux_dP, dCompFlux_dC, ip, k_up_g);

          phaseFlux += phaseFluxG;
          for(localIndex ke = 0; ke < stencilSize; ++ke)
          {
            dPhaseFlux_dP[ke] += dPhaseFluxG_dP[ke];
            for( localIndex ic = 0; ic < NC; ++ic )
              dPhaseFlux_dC[ke][ic] += dPhaseFluxG_dC[ke][ic];
          }

        }
      }

      //(opt 2) restore molar phaseDens at the end of the call
/*      k_up = -1;
      if( is_ppu )
      {
        k_up = ( presGrad - gravHead >= 0 ) ? 0 : 1;
      }
      else if( is_pu )
      {
        k_up = (totFlux_unw >= 0) ? 0 : 1 ;
        if( (std::fabs(gravHead) <= std::fabs(minGravHead)) && std::fabs(minGravHead)>0 )
          k_up = ( k_up == 1 ) ? 0 : 1;
      }
      else if( is_hu )
      {
        k_up = (totFlux_unw >= 0) ? 0 : 1 ;
      }
      densMult(phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, ip, k_up); */

      /*update phaseComp from viscous part */
      dPhaseComp_dX(phaseFluxV, dPhaseFluxV_dP, dPhaseFluxV_dC, compFlux, dCompFlux_dP, dCompFlux_dC, ip, k_up);

    }

  }//end If UT_FORM



  // populate local flux vector and derivatives
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    localFlux[ic]      =  dt * compFlux[ic];
    localFlux[NC + ic] = -dt * compFlux[ic];

    std::cerr << "\n***********\n ic" << ic << " compFlux : " << compFlux[ic] <<std::endl;

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      localIndex const localDofIndexPres = ke * NDOF;
      localFluxJacobian[ic][localDofIndexPres] = dt * dCompFlux_dP[ke][ic];
      localFluxJacobian[NC + ic][localDofIndexPres] = -dt * dCompFlux_dP[ke][ic];

      std::cerr << "*********** \t ke" << ke << " compFlux : " << dCompFlux_dP[ke][ic] <<std::endl;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
        localFluxJacobian[ic][localDofIndexComp] = dt * dCompFlux_dC[ke][ic][jc];
        localFluxJacobian[NC + ic][localDofIndexComp] = -dt * dCompFlux_dC[ke][ic][jc];

        std::cerr << "*********** \t \t jc" << jc << " compFlux : " << dCompFlux_dC[ke][ic][jc] << std::endl << std::endl;
      }
    }
  }
}

template< localIndex NC, typename STENCIL_TYPE , bool IS_UT_FORM >
void
FluxKernel::
  Launch( localIndex const numPhases,
          STENCIL_TYPE const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
          ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
          ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
          ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
          integer const capPressureFlag,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  typename STENCIL_TYPE::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename STENCIL_TYPE::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename STENCIL_TYPE::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename STENCIL_TYPE::WeightContainerViewConstType const & weights = stencil.getWeights();

  localIndex constexpr NUM_ELEMS   = STENCIL_TYPE::NUM_POINT_IN_FLUX;
  localIndex constexpr MAX_STENCIL = STENCIL_TYPE::MAX_STENCIL_SIZE;

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    // TODO: hack! for MPFA, etc. must obtain proper size from e.g. seri
    localIndex const stencilSize = MAX_STENCIL;
    localIndex constexpr NDOF = NC + 1;

    stackArray1d< real64, NUM_ELEMS * NC >                      localFlux( NUM_ELEMS * NC );
    stackArray2d< real64, NUM_ELEMS * NC * MAX_STENCIL * NDOF > localFluxJacobian( NUM_ELEMS * NC, stencilSize * NDOF );

    FluxKernel::Compute< NC, NUM_ELEMS, MAX_STENCIL, IS_UT_FORM >( numPhases,
                                                       stencilSize,
                                                       seri[iconn],
                                                       sesri[iconn],
                                                       sei[iconn],
                                                       weights[iconn],
                                                       pres,
                                                       dPres,
                                                       gravCoef,
                                                       phaseMob,
                                                       dPhaseMob_dPres,
                                                       dPhaseMob_dComp,
                                                       dPhaseVolFrac_dPres,
                                                       dPhaseVolFrac_dComp,
                                                       dCompFrac_dCompDens,
                                                                   phaseDens,
                                                                   dPhaseDens_dPres,
                                                                   dPhaseDens_dComp,
                                                       phaseMassDens,
                                                       dPhaseMassDens_dPres,
                                                       dPhaseMassDens_dComp,
                                                       phaseCompFrac,
                                                       dPhaseCompFrac_dPres,
                                                       dPhaseCompFrac_dComp,
                                                       phaseCapPressure,
                                                       dPhaseCapPressure_dPhaseVolFrac,
                                                       capPressureFlag,
                                                       dt,
                                                       localFlux,
                                                       localFluxJacobian );

    // populate dof indices
    globalIndex dofColIndices[ MAX_STENCIL * NDOF ];
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      globalIndex const offset = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];

      for( localIndex jdof = 0; jdof < NDOF; ++jdof )
      {
        dofColIndices[i * NDOF + jdof] = offset + jdof;
      }
    }

    // TODO: apply equation/variable change transformation(s)

    // Add to residual/jacobian
    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( localMatrix.numRows(), localRow + NC );

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow + ic], localFlux[i * NC + ic] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + ic,
                                                                            dofColIndices,
                                                                            localFluxJacobian[i * NC + ic].dataIfContiguous(),
                                                                            stencilSize * NDOF );
        }
      }
    }
  } );
}

#define INST_FluxKernel( NC, STENCIL_TYPE, IS_UT_FORM ) \
  template \
  void FluxKernel:: \
    Launch< NC, STENCIL_TYPE, IS_UT_FORM >( localIndex const numPhases, \
                                STENCIL_TYPE const & stencil, \
                                globalIndex const rankOffset, \
                                ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber, \
                                ElementViewConst< arrayView1d< integer const > > const & ghostRank, \
                                ElementViewConst< arrayView1d< real64 const > > const & pres, \
                                ElementViewConst< arrayView1d< real64 const > > const & dPres, \
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef, \
                                ElementViewConst< arrayView2d< real64 const > > const & phaseMob, \
                                ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres, \
                                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp, \
                                ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres, \
                                ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp, \
                                ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens, \
                                ElementViewConst< arrayView3d< real64 const > > const & phaseDens, \
                                ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres, \
                                ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp, \
                                ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens, \
                                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres, \
                                ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp, \
                                ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac, \
                                ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres, \
                                ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp, \
                                ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure, \
                                ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac, \
                                integer const capPressureFlag, \
                                real64 const dt, \
                                CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 1, CellElementStencilTPFA, false );
INST_FluxKernel( 2, CellElementStencilTPFA, false );
INST_FluxKernel( 3, CellElementStencilTPFA, false );
INST_FluxKernel( 4, CellElementStencilTPFA, false );
INST_FluxKernel( 5, CellElementStencilTPFA, false );

INST_FluxKernel( 1, CellElementStencilTPFA, true );
INST_FluxKernel( 2, CellElementStencilTPFA, true );
INST_FluxKernel( 3, CellElementStencilTPFA, true );
INST_FluxKernel( 4, CellElementStencilTPFA, true );
INST_FluxKernel( 5, CellElementStencilTPFA, true );

INST_FluxKernel( 1, FaceElementStencil, false);
INST_FluxKernel( 2, FaceElementStencil, false );
INST_FluxKernel( 3, FaceElementStencil, false );
INST_FluxKernel( 4, FaceElementStencil, false );
INST_FluxKernel( 5, FaceElementStencil, false );

INST_FluxKernel( 1, FaceElementStencil, true);
INST_FluxKernel( 2, FaceElementStencil, true );
INST_FluxKernel( 3, FaceElementStencil, true );
INST_FluxKernel( 4, FaceElementStencil, true );
INST_FluxKernel( 5, FaceElementStencil, true );

#undef INST_FluxKernel

/******************************** VolumeBalanceKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
VolumeBalanceKernel::
  Compute( real64 const & volume,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d< real64 const > const & phaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dCompDens,
           real64 & localVolBalance,
           real64 * const localVolBalanceJacobian )
{
  localIndex constexpr NDOF = NC + 1;

  real64 const poro     = porosityRef * pvMult;
  real64 const dPoro_dP = porosityRef * dPvMult_dPres;

  real64 const poreVol     = volume * poro;
  real64 const dPoreVol_dP = volume * dPoro_dP;

  localVolBalance = 1.0;
  for( localIndex i = 0; i < NDOF; ++i )
  {
    localVolBalanceJacobian[i] = 0.0;
  }

  // sum contributions to component accumulation from each phase
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    localVolBalance -= phaseVolFrac[ip];
    localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ip];

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompDens[ip][jc];
    }
  }

  // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
  for( localIndex idof = 0; idof < NDOF; ++idof )
  {
    localVolBalanceJacobian[idof] *= poreVol;
  }
  localVolBalanceJacobian[0] += dPoreVol_dP * localVolBalance;
  localVolBalance *= poreVol;
}

template< localIndex NC, localIndex NP >
void
VolumeBalanceKernel::
  Launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & porosityRef,
          arrayView2d< real64 const > const & pvMult,
          arrayView2d< real64 const > const & dPvMult_dPres,
          arrayView2d< real64 const > const & phaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    if( elemGhostRank[ei] >= 0 )
      return;

    localIndex constexpr NDOF = NC + 1;

    real64 localVolBalance;
    real64 localVolBalanceJacobian[NDOF];

    Compute< NC, NP >( volume[ei],
                       porosityRef[ei],
                       pvMult[ei][0],
                       dPvMult_dPres[ei][0],
                       phaseVolFrac[ei],
                       dPhaseVolFrac_dPres[ei],
                       dPhaseVolFrac_dCompDens[ei],
                       localVolBalance,
                       localVolBalanceJacobian );

    // get equation/dof indices
    localIndex const localRow = dofNumber[ei] + NC - rankOffset;
    globalIndex dofIndices[NDOF];
    for( localIndex jdof = 0; jdof < NDOF; ++jdof )
    {
      dofIndices[jdof] = dofNumber[ei] + jdof;
    }

    // TODO: apply equation/variable change transformation(s)

    // add contribution to residual and jacobian
    localRhs[localRow] += localVolBalance;
    localMatrix.addToRow< serialAtomic >( localRow,
                                          dofIndices,
                                          localVolBalanceJacobian,
                                          NDOF );
  } );
}

#define INST_VolumeBalanceKernel( NC, NP ) \
  template \
  void VolumeBalanceKernel:: \
    Launch< NC, NP >( localIndex const size, \
                      globalIndex const rankOffset, \
                      arrayView1d< globalIndex const > const & dofNumber, \
                      arrayView1d< integer const > const & elemGhostRank, \
                      arrayView1d< real64 const > const & volume, \
                      arrayView1d< real64 const > const & porosityRef, \
                      arrayView2d< real64 const > const & pvMult, \
                      arrayView2d< real64 const > const & dPvMult_dPres, \
                      arrayView2d< real64 const > const & phaseVolFrac, \
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens, \
                      CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                      arrayView1d< real64 > const & localRhs )

INST_VolumeBalanceKernel( 1, 1 );
INST_VolumeBalanceKernel( 2, 1 );
INST_VolumeBalanceKernel( 3, 1 );
INST_VolumeBalanceKernel( 4, 1 );
INST_VolumeBalanceKernel( 5, 1 );

INST_VolumeBalanceKernel( 1, 2 );
INST_VolumeBalanceKernel( 2, 2 );
INST_VolumeBalanceKernel( 3, 2 );
INST_VolumeBalanceKernel( 4, 2 );
INST_VolumeBalanceKernel( 5, 2 );

INST_VolumeBalanceKernel( 1, 3 );
INST_VolumeBalanceKernel( 2, 3 );
INST_VolumeBalanceKernel( 3, 3 );
INST_VolumeBalanceKernel( 4, 3 );
INST_VolumeBalanceKernel( 5, 3 );

#undef INST_VolumeBalanceKernel

} // namespace CompositionalMultiphaseFlowKernels

} // namespace geosx
