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
 * @file IsothermalCompositionalMultiphaseBaseKernels.cpp
 */

#include "IsothermalCompositionalMultiphaseBaseKernels.hpp"
#include "CompositionalMultiphaseUtilities.hpp"

namespace geosx
{

namespace IsothermalCompositionalMultiphaseBaseKernels
{

/******************************** ComponentFractionKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
void
ComponentFractionKernel::
  compute( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const compDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const dCompDens,
           arraySlice1d< real64, compflow::USD_COMP - 1 > const compFrac,
           arraySlice2d< real64, compflow::USD_COMP_DC - 1 > const dCompFrac_dCompDens )
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
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView2d< real64, compflow::USD_COMP > const & compFrac,
          arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    compute< NC >( compDens[a],
                   dCompDens[a],
                   compFrac[a],
                   dCompFrac_dCompDens[a] );
  } );
}

template< localIndex NC >
void
ComponentFractionKernel::
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView2d< real64, compflow::USD_COMP > const & compFrac,
          arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    compute< NC >( compDens[a],
                   dCompDens[a],
                   compFrac[a],
                   dCompFrac_dCompDens[a] );
  } );
}

#define INST_ComponentFractionKernel( NC ) \
  template \
  void ComponentFractionKernel:: \
    launch< NC >( localIndex const size, \
                  arrayView2d< real64 const, compflow::USD_COMP > const & compDens, \
                  arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens, \
                  arrayView2d< real64, compflow::USD_COMP > const & compFrac, \
                  arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens ); \
  template \
  void ComponentFractionKernel:: \
    launch< NC >( SortedArrayView< localIndex const > const & targetSet, \
                  arrayView2d< real64 const, compflow::USD_COMP > const & compDens, \
                  arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens, \
                  arrayView2d< real64, compflow::USD_COMP > const & compFrac, \
                  arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens )

INST_ComponentFractionKernel( 1 );
INST_ComponentFractionKernel( 2 );
INST_ComponentFractionKernel( 3 );
INST_ComponentFractionKernel( 4 );
INST_ComponentFractionKernel( 5 );

#undef INST_ComponentFractionKernel

/******************************** PhaseVolumeFractionKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
PhaseVolumeFractionKernel::
  compute( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & compDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & dCompDens,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
           arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & phaseFrac,
           arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & dPhaseFrac_dPres,
           arraySlice2d< real64 const, constitutive::multifluid::USD_PHASE_DC - 2 > const & dPhaseFrac_dComp,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp )
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

    // set the saturation to zero if the phase is absent
    bool const phaseExists = (phaseFrac[ip] > 0);
    if( !phaseExists )
    {
      phaseVolFrac[ip] = 0.;
      dPhaseVolFrac_dPres[ip] = 0.;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseVolFrac_dComp[ip][jc] = 0.;
      }
      continue;
    }

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
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & dPhaseFrac_dPres,
          arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    compute< NC, NP >( compDens[a],
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
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & dPhaseFrac_dPres,
          arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    compute< NC, NP >( compDens[a],
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
    launch< NC, NP >( localIndex const size, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & compDens, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens, \
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens, \
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac, \
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & dPhaseFrac_dPres, \
                      arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp, \
                      arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac, \
                      arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp ); \
  template \
  void \
  PhaseVolumeFractionKernel:: \
    launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & compDens, \
                      arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens, \
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseDens, \
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseFrac, \
                      arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & dPhaseFrac_dPres, \
                      arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp, \
                      arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac, \
                      arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp )

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

} // namespace IsothermalCompositionalMultiphaseBaseKernels

} // namespace geosx
