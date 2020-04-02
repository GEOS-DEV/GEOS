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
 * @file CompositionalMultiphaseFlowKernels.cpp
 */

#include "CompositionalMultiphaseFlowKernels.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFlowKernels
{

/******************************** ComponentFractionKernel ********************************/

template< localIndex NC >
inline void
ComponentFractionKernel::Compute( arraySlice1d< real64 const > const compDens,
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

inline void
ComponentFractionKernel::Compute( localIndex const NC,
                                  arraySlice1d< real64 const > const compDens,
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
void ComponentFractionKernel::Launch( localIndex const size,
                                      arrayView2d< real64 const > const & compDens,
                                      arrayView2d< real64 const > const & dCompDens,
                                      arrayView2d< real64 > const & compFrac,
                                      arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< serialPolicy >( size, [=] ( localIndex const a )
  {
    Compute< NC >( compDens[a],
                   dCompDens[a],
                   compFrac[a],
                   dCompFrac_dCompDens[a] );
  } );
}

void ComponentFractionKernel::Launch( localIndex const NC,
                                      localIndex const size,
                                      arrayView2d< real64 const > const & compDens,
                                      arrayView2d< real64 const > const & dCompDens,
                                      arrayView2d< real64 > const & compFrac,
                                      arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< serialPolicy >( size, [=] ( localIndex const a )
  {
    Compute( NC,
             compDens[a],
             dCompDens[a],
             compFrac[a],
             dCompFrac_dCompDens[a] );
  } );
}

template< localIndex NC >
void ComponentFractionKernel::Launch( SortedArrayView< localIndex const > const & targetSet,
                                      arrayView2d< real64 const > const & compDens,
                                      arrayView2d< real64 const > const & dCompDens,
                                      arrayView2d< real64 > const & compFrac,
                                      arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< serialPolicy >( targetSet.size(), [=] ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute< NC >( compDens[a],
                   dCompDens[a],
                   compFrac[a],
                   dCompFrac_dCompDens[a] );
  } );
}

void ComponentFractionKernel::Launch( localIndex const NC,
                                      SortedArrayView< localIndex const > const & targetSet,
                                      arrayView2d< real64 const > const & compDens,
                                      arrayView2d< real64 const > const & dCompDens,
                                      arrayView2d< real64 > const & compFrac,
                                      arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< serialPolicy >( targetSet.size(), [=] ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute( NC,
             compDens[a],
             dCompDens[a],
             compFrac[a],
             dCompFrac_dCompDens[a] );
  } );
}

#define INST_ComponentFractionKernel( NC ) \
  template \
  void ComponentFractionKernel::Launch< NC >( localIndex const size, \
                                              arrayView2d< real64 const > const & compDens, \
                                              arrayView2d< real64 const > const & dCompDens, \
                                              arrayView2d< real64 > const & compFrac, \
                                              arrayView3d< real64 > const & dCompFrac_dCompDens ); \
  template \
  void ComponentFractionKernel::Launch< NC >( SortedArrayView< localIndex const > const & targetSet, \
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
inline void
PhaseVolumeFractionKernel::Compute( arraySlice1d< real64 const > const compDens,
                                    arraySlice1d< real64 const > const dCompDens,
                                    arraySlice2d< real64 const > const dCompFrac_dCompDens,
                                    arraySlice1d< real64 const > const phaseDens,
                                    arraySlice1d< real64 const > const dPhaseDens_dPres,
                                    arraySlice2d< real64 const > const dPhaseDens_dComp,
                                    arraySlice1d< real64 const > const phaseFrac,
                                    arraySlice1d< real64 const > const dPhaseFrac_dPres,
                                    arraySlice2d< real64 const > const dPhaseFrac_dComp,
                                    arraySlice1d< real64 > const phaseVolFrac,
                                    arraySlice1d< real64 > const dPhaseVolFrac_dPres,
                                    arraySlice2d< real64 > const dPhaseVolFrac_dComp )
{
  stackArray1d< real64, NC > work( NC );

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

inline void
PhaseVolumeFractionKernel::Compute( localIndex const NC, localIndex const NP,
                                    arraySlice1d< real64 const > const compDens,
                                    arraySlice1d< real64 const > const dCompDens,
                                    arraySlice2d< real64 const > const dCompFrac_dCompDens,
                                    arraySlice1d< real64 const > const phaseDens,
                                    arraySlice1d< real64 const > const dPhaseDens_dPres,
                                    arraySlice2d< real64 const > const dPhaseDens_dComp,
                                    arraySlice1d< real64 const > const phaseFrac,
                                    arraySlice1d< real64 const > const dPhaseFrac_dPres,
                                    arraySlice2d< real64 const > const dPhaseFrac_dComp,
                                    arraySlice1d< real64 > const phaseVolFrac,
                                    arraySlice1d< real64 > const dPhaseVolFrac_dPres,
                                    arraySlice2d< real64 > const dPhaseVolFrac_dComp )
{
  localIndex constexpr maxNC = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  stackArray1d< real64, maxNC > work( NC );

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
void PhaseVolumeFractionKernel::Launch( localIndex const size,
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
  forAll< serialPolicy >( size, [=] ( localIndex const a )
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

void PhaseVolumeFractionKernel::Launch( localIndex const NC, localIndex const NP,
                                        localIndex const size,
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
  forAll< serialPolicy >( size, [=] ( localIndex const a )
  {
    Compute( NC, NP,
             compDens[a],
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
void PhaseVolumeFractionKernel::Launch( SortedArrayView< localIndex const > const & targetSet,
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
  forAll< serialPolicy >( targetSet.size(), [=] ( localIndex const i )
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

void PhaseVolumeFractionKernel::Launch( localIndex const NC, localIndex const NP,
                                        SortedArrayView< localIndex const > const & targetSet,
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
  forAll< serialPolicy >( targetSet.size(), [=] ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute( NC, NP,
             compDens[a],
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
  void PhaseVolumeFractionKernel::Launch< NC, NP >( localIndex const size, \
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
  void PhaseVolumeFractionKernel::Launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
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
inline void
PhaseMobilityKernel::Compute( arraySlice2d< real64 const > const dCompFrac_dCompDens,
                              arraySlice1d< real64 const > const phaseDens,
                              arraySlice1d< real64 const > const dPhaseDens_dPres,
                              arraySlice2d< real64 const > const dPhaseDens_dComp,
                              arraySlice1d< real64 const > const phaseVisc,
                              arraySlice1d< real64 const > const dPhaseVisc_dPres,
                              arraySlice2d< real64 const > const dPhaseVisc_dComp,
                              arraySlice1d< real64 const > const phaseRelPerm,
                              arraySlice2d< real64 const > const dPhaseRelPerm_dPhaseVolFrac,
                              arraySlice1d< real64 const > const dPhaseVolFrac_dPres,
                              arraySlice2d< real64 const > const dPhaseVolFrac_dComp,
                              arraySlice1d< real64 > const phaseMob,
                              arraySlice1d< real64 > const dPhaseMob_dPres,
                              arraySlice2d< real64 > const dPhaseMob_dComp )
{
  stackArray1d< real64, NC > dRelPerm_dC( NC );
  stackArray1d< real64, NC > dDens_dC( NC );
  stackArray1d< real64, NC > dVisc_dC( NC );

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const density = phaseDens[ip];
    real64 const dDens_dP = dPhaseDens_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dDens_dC );

    real64 const viscosity = phaseVisc[ip];
    real64 const dVisc_dP = dPhaseVisc_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseVisc_dComp[ip], dVisc_dC );

    real64 const relPerm = phaseRelPerm[ip];
    real64 dRelPerm_dP = 0.0;
    dRelPerm_dC = 0.0;

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

    // compositional derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] * density / viscosity
                                + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
    }
  }
}

inline void
PhaseMobilityKernel::Compute( localIndex const NC, localIndex const NP,
                              arraySlice2d< real64 const > const dCompFrac_dCompDens,
                              arraySlice1d< real64 const > const phaseDens,
                              arraySlice1d< real64 const > const dPhaseDens_dPres,
                              arraySlice2d< real64 const > const dPhaseDens_dComp,
                              arraySlice1d< real64 const > const phaseVisc,
                              arraySlice1d< real64 const > const dPhaseVisc_dPres,
                              arraySlice2d< real64 const > const dPhaseVisc_dComp,
                              arraySlice1d< real64 const > const phaseRelPerm,
                              arraySlice2d< real64 const > const dPhaseRelPerm_dPhaseVolFrac,
                              arraySlice1d< real64 const > const dPhaseVolFrac_dPres,
                              arraySlice2d< real64 const > const dPhaseVolFrac_dComp,
                              arraySlice1d< real64 > const phaseMob,
                              arraySlice1d< real64 > const dPhaseMob_dPres,
                              arraySlice2d< real64 > const dPhaseMob_dComp )
{
  localIndex constexpr maxNC = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  stackArray1d< real64, maxNC > dRelPerm_dC( NC );
  stackArray1d< real64, maxNC > dDens_dC( NC );
  stackArray1d< real64, maxNC > dVisc_dC( NC );

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const density = phaseDens[ip];
    real64 const dDens_dP = dPhaseDens_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dDens_dC );

    real64 const viscosity = phaseVisc[ip];
    real64 const dVisc_dP = dPhaseVisc_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseVisc_dComp[ip], dVisc_dC );

    real64 const relPerm = phaseRelPerm[ip];
    real64 dRelPerm_dP = 0.0;
    dRelPerm_dC = 0.0;

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

    // compositional derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] * density / viscosity
                                + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
    }
  }
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::Launch( localIndex const size,
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
  forAll< serialPolicy >( size, [=] ( localIndex const a )
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

void PhaseMobilityKernel::Launch( localIndex const NC, localIndex const NP,
                                  localIndex const size,
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
  forAll< serialPolicy >( size, [=] ( localIndex const a )
  {
    Compute( NC, NP,
             dCompFrac_dCompDens[a],
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
void PhaseMobilityKernel::Launch( SortedArrayView< localIndex const > const & targetSet,
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
  forAll< serialPolicy >( targetSet.size(), [=] ( localIndex const i )
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

void PhaseMobilityKernel::Launch( localIndex const NC, localIndex const NP,
                                  SortedArrayView< localIndex const > const & targetSet,
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
  forAll< serialPolicy >( targetSet.size(), [=] ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    Compute( NC, NP,
             dCompFrac_dCompDens[a],
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
  void PhaseMobilityKernel::Launch< NC, NP >( localIndex const size, \
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
  void PhaseMobilityKernel::Launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
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

} // namespace CompositionalMultiphaseFlowKernels

} // namespace geosx
