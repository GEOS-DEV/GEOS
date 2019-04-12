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
 * @file CompositionalMultiphaseFlowKernels.cpp
 */

#include "CompositionalMultiphaseFlowKernels.hpp"

#include "constitutive/Fluid/MultiFluidBase.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFlowKernels
{

/******************************** UpdateComponentFractionKernel ********************************/

template<localIndex NC>
inline RAJA_HOST_DEVICE void
UpdateComponentFractionKernel::Compute( arraySlice1d<real64 const> compDens,
                                        arraySlice1d<real64 const> dCompDens,
                                        arraySlice1d<real64> compFrac,
                                        arraySlice2d<real64> dCompFrac_dCompDens )
{
  real64 totalDensity = 0.0;

  for (localIndex ic = 0; ic < NC; ++ic)
  {
    totalDensity += compDens[ic] + dCompDens[ic];
  }

  real64 const totalDensityInv = 1.0 / totalDensity;

  for (localIndex ic = 0; ic < NC; ++ic)
  {
    compFrac[ic] = (compDens[ic] + dCompDens[ic]) * totalDensityInv;
    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dCompFrac_dCompDens[ic][jc] = - compFrac[ic] * totalDensityInv;
    }
    dCompFrac_dCompDens[ic][ic] += totalDensityInv;
  }
}

inline RAJA_HOST_DEVICE void
UpdateComponentFractionKernel::Compute( localIndex NC,
                                        arraySlice1d<real64 const> compDens,
                                        arraySlice1d<real64 const> dCompDens,
                                        arraySlice1d<real64> compFrac,
                                        arraySlice2d<real64> dCompFrac_dCompDens )
{
  real64 totalDensity = 0.0;

  for (localIndex ic = 0; ic < NC; ++ic)
  {
    totalDensity += compDens[ic] + dCompDens[ic];
  }

  real64 const totalDensityInv = 1.0 / totalDensity;

  for (localIndex ic = 0; ic < NC; ++ic)
  {
    compFrac[ic] = (compDens[ic] + dCompDens[ic]) * totalDensityInv;
    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dCompFrac_dCompDens[ic][jc] = - compFrac[ic] * totalDensityInv;
    }
    dCompFrac_dCompDens[ic][ic] += totalDensityInv;
  }
}

template<localIndex NC>
void UpdateComponentFractionKernel::Launch( localIndex begin, localIndex end,
                                            arrayView2d<real64 const> const & compDens,
                                            arrayView2d<real64 const> const & dCompDens,
                                            arrayView2d<real64> const & compFrac,
                                            arrayView3d<real64> const & dCompFrac_dCompDens )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
  {
    Compute<NC>( compDens[a],
                 dCompDens[a],
                 compFrac[a],
                 dCompFrac_dCompDens[a] );
  } );
}

void UpdateComponentFractionKernel::Launch( localIndex NC,
                                            localIndex begin, localIndex end,
                                            arrayView2d<real64 const> const & compDens,
                                            arrayView2d<real64 const> const & dCompDens,
                                            arrayView2d<real64> const & compFrac,
                                            arrayView3d<real64> const & dCompFrac_dCompDens )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
  {
    Compute( NC,
             compDens[a],
             dCompDens[a],
             compFrac[a],
             dCompFrac_dCompDens[a] );
  } );
}

#define LAUNCH_UpdateComponentFractionKernel(NC) \
template \
void UpdateComponentFractionKernel::Launch<NC>( localIndex begin, localIndex end, \
                                                arrayView2d<real64 const> const & compDens, \
                                                arrayView2d<real64 const> const & dCompDens, \
                                                arrayView2d<real64> const & compFrac, \
                                                arrayView3d<real64> const & dCompFrac_dCompDens )

LAUNCH_UpdateComponentFractionKernel(1);
LAUNCH_UpdateComponentFractionKernel(2);
LAUNCH_UpdateComponentFractionKernel(3);
LAUNCH_UpdateComponentFractionKernel(4);
LAUNCH_UpdateComponentFractionKernel(5);
LAUNCH_UpdateComponentFractionKernel(6);
LAUNCH_UpdateComponentFractionKernel(7);
LAUNCH_UpdateComponentFractionKernel(8);
LAUNCH_UpdateComponentFractionKernel(9);
LAUNCH_UpdateComponentFractionKernel(10);

#undef LAUNCH_UpdateComponentFractionKernel

/******************************** UpdatePhaseVolumeFractionKernel ********************************/

template<localIndex NC, localIndex NP>
inline RAJA_HOST_DEVICE void
UpdatePhaseVolumeFractionKernel::Compute( arraySlice1d<real64 const> compDens,
                                          arraySlice1d<real64 const> dCompDens,
                                          arraySlice2d<real64 const> dCompFrac_dCompDens,
                                          arraySlice1d<real64 const> phaseDens,
                                          arraySlice1d<real64 const> dPhaseDens_dPres,
                                          arraySlice2d<real64 const> dPhaseDens_dComp,
                                          arraySlice1d<real64 const> phaseFrac,
                                          arraySlice1d<real64 const> dPhaseFrac_dPres,
                                          arraySlice2d<real64 const> dPhaseFrac_dComp,
                                          arraySlice1d<real64> phaseVolFrac,
                                          arraySlice1d<real64> dPhaseVolFrac_dPres,
                                          arraySlice2d<real64> dPhaseVolFrac_dComp )
{
  stackArray1d<real64, NC> work( NC );

  // compute total density from component partial densities
  real64 totalDensity = 0.0;
  real64 const dTotalDens_dCompDens = 1.0;
  for (localIndex ic = 0; ic < NC; ++ic)
  {
    totalDensity += compDens[ic] + dCompDens[ic];
  }

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
    real64 const phaseDensInv = 1.0 / phaseDens[ip];

    // compute saturation and derivatives except multiplying by the total density
    phaseVolFrac[ip] = phaseFrac[ip] * phaseDensInv;

    dPhaseVolFrac_dPres[ip] =
      (dPhaseFrac_dPres[ip] - phaseVolFrac[ip] * dPhaseDens_dPres[ip]) * phaseDensInv;

    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dPhaseVolFrac_dComp[ip][jc] =
        (dPhaseFrac_dComp[ip][jc] - phaseVolFrac[ip] * dPhaseDens_dComp[ip][jc]) * phaseDensInv;
    }

    // apply chain rule to convert derivatives from global component fractions to densities
    applyChainRuleInPlace( NC, dCompFrac_dCompDens, dPhaseVolFrac_dComp[ip], work );

    // now finalize the computation by multiplying by total density
    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dPhaseVolFrac_dComp[ip][jc] *= totalDensity;
      dPhaseVolFrac_dComp[ip][jc] += phaseVolFrac[ip] * dTotalDens_dCompDens;
    }

    phaseVolFrac[ip] *= totalDensity;
    dPhaseVolFrac_dPres[ip] *= totalDensity;
  }
}

inline RAJA_HOST_DEVICE void
UpdatePhaseVolumeFractionKernel::Compute( localIndex NC, localIndex NP,
                                          arraySlice1d<real64 const> compDens,
                                          arraySlice1d<real64 const> dCompDens,
                                          arraySlice2d<real64 const> dCompFrac_dCompDens,
                                          arraySlice1d<real64 const> phaseDens,
                                          arraySlice1d<real64 const> dPhaseDens_dPres,
                                          arraySlice2d<real64 const> dPhaseDens_dComp,
                                          arraySlice1d<real64 const> phaseFrac,
                                          arraySlice1d<real64 const> dPhaseFrac_dPres,
                                          arraySlice2d<real64 const> dPhaseFrac_dComp,
                                          arraySlice1d<real64> phaseVolFrac,
                                          arraySlice1d<real64> dPhaseVolFrac_dPres,
                                          arraySlice2d<real64> dPhaseVolFrac_dComp )
{
  localIndex constexpr maxNC = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  stackArray1d<real64, maxNC> work( NC );

  // compute total density from component partial densities
  real64 totalDensity = 0.0;
  real64 const dTotalDens_dCompDens = 1.0;
  for (localIndex ic = 0; ic < NC; ++ic)
  {
    totalDensity += compDens[ic] + dCompDens[ic];
  }

  for (localIndex ip = 0; ip < NP; ++ip)
  {
    // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
    real64 const phaseDensInv = 1.0 / phaseDens[ip];

    // compute saturation and derivatives except multiplying by the total density
    phaseVolFrac[ip] = phaseFrac[ip] * phaseDensInv;

    dPhaseVolFrac_dPres[ip] =
      (dPhaseFrac_dPres[ip] - phaseVolFrac[ip] * dPhaseDens_dPres[ip]) * phaseDensInv;

    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dPhaseVolFrac_dComp[ip][jc] =
        (dPhaseFrac_dComp[ip][jc] - phaseVolFrac[ip] * dPhaseDens_dComp[ip][jc]) * phaseDensInv;
    }

    // apply chain rule to convert derivatives from global component fractions to densities
    applyChainRuleInPlace( NC, dCompFrac_dCompDens, dPhaseVolFrac_dComp[ip], work );

    // now finalize the computation by multiplying by total density
    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dPhaseVolFrac_dComp[ip][jc] *= totalDensity;
      dPhaseVolFrac_dComp[ip][jc] += phaseVolFrac[ip] * dTotalDens_dCompDens;
    }

    phaseVolFrac[ip] *= totalDensity;
    dPhaseVolFrac_dPres[ip] *= totalDensity;
  }
}

template<localIndex NC, localIndex NP>
void UpdatePhaseVolumeFractionKernel::Launch( localIndex begin, localIndex end,
                                              arrayView2d<real64 const> const & compDens,
                                              arrayView2d<real64 const> const & dCompDens,
                                              arrayView3d<real64 const> const & dCompFrac_dCompDens,
                                              arrayView3d<real64 const> const & phaseDens,
                                              arrayView3d<real64 const> const & dPhaseDens_dPres,
                                              arrayView4d<real64 const> const & dPhaseDens_dComp,
                                              arrayView3d<real64 const> const & phaseFrac,
                                              arrayView3d<real64 const> const & dPhaseFrac_dPres,
                                              arrayView4d<real64 const> const & dPhaseFrac_dComp,
                                              arrayView2d<real64> const & phaseVolFrac,
                                              arrayView2d<real64> const & dPhaseVolFrac_dPres,
                                              arrayView3d<real64> const & dPhaseVolFrac_dComp )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
  {
    Compute<NC, NP>( compDens[a],
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

void UpdatePhaseVolumeFractionKernel::Launch( localIndex NC, localIndex NP,
                                              localIndex begin, localIndex end,
                                              arrayView2d<real64 const> const & compDens,
                                              arrayView2d<real64 const> const & dCompDens,
                                              arrayView3d<real64 const> const & dCompFrac_dCompDens,
                                              arrayView3d<real64 const> const & phaseDens,
                                              arrayView3d<real64 const> const & dPhaseDens_dPres,
                                              arrayView4d<real64 const> const & dPhaseDens_dComp,
                                              arrayView3d<real64 const> const & phaseFrac,
                                              arrayView3d<real64 const> const & dPhaseFrac_dPres,
                                              arrayView4d<real64 const> const & dPhaseFrac_dComp,
                                              arrayView2d<real64> const & phaseVolFrac,
                                              arrayView2d<real64> const & dPhaseVolFrac_dPres,
                                              arrayView3d<real64> const & dPhaseVolFrac_dComp )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
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

#define LAUNCH_UpdatePhaseVolumeFractionKernel(NC,NP) \
template \
void UpdatePhaseVolumeFractionKernel::Launch<NC,NP>( localIndex begin, localIndex end, \
                                                      arrayView2d<real64 const> const & compDens, \
                                                      arrayView2d<real64 const> const & dCompDens, \
                                                      arrayView3d<real64 const> const & dCompFrac_dCompDens, \
                                                      arrayView3d<real64 const> const & phaseDens, \
                                                      arrayView3d<real64 const> const & dPhaseDens_dPres, \
                                                      arrayView4d<real64 const> const & dPhaseDens_dComp, \
                                                      arrayView3d<real64 const> const & phaseFrac, \
                                                      arrayView3d<real64 const> const & dPhaseFrac_dPres, \
                                                      arrayView4d<real64 const> const & dPhaseFrac_dComp, \
                                                      arrayView2d<real64> const & phaseVolFrac, \
                                                      arrayView2d<real64> const & dPhaseVolFrac_dPres, \
                                                      arrayView3d<real64> const & dPhaseVolFrac_dComp )

LAUNCH_UpdatePhaseVolumeFractionKernel(1,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(2,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(3,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(4,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(5,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(6,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(7,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(8,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(9,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(10,1);

LAUNCH_UpdatePhaseVolumeFractionKernel(1,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(2,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(3,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(4,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(5,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(6,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(7,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(8,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(9,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(10,2);

LAUNCH_UpdatePhaseVolumeFractionKernel(1,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(2,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(3,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(4,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(5,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(6,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(7,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(8,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(9,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(10,3);

#undef LAUNCH_UpdatePhaseVolumeFractionKernel

/******************************** UpdatePhaseMobilityKernel ********************************/

template<localIndex NC, localIndex NP>
inline RAJA_HOST_DEVICE void
UpdatePhaseMobilityKernel::Compute( arraySlice2d<real64 const> dCompFrac_dCompDens,
                                    arraySlice1d<real64 const> phaseDens,
                                    arraySlice1d<real64 const> dPhaseDens_dPres,
                                    arraySlice2d<real64 const> dPhaseDens_dComp,
                                    arraySlice1d<real64 const> phaseVisc,
                                    arraySlice1d<real64 const> dPhaseVisc_dPres,
                                    arraySlice2d<real64 const> dPhaseVisc_dComp,
                                    arraySlice1d<real64 const> phaseRelPerm,
                                    arraySlice2d<real64 const> dPhaseRelPerm_dPhaseVolFrac,
                                    arraySlice1d<real64 const> dPhaseVolFrac_dPres,
                                    arraySlice2d<real64 const> dPhaseVolFrac_dComp,
                                    arraySlice1d<real64> phaseMob,
                                    arraySlice1d<real64> dPhaseMob_dPres,
                                    arraySlice2d<real64> dPhaseMob_dComp )
{
  stackArray1d<real64, NC> dRelPerm_dC( NC );
  stackArray1d<real64, NC> dDens_dC( NC );
  stackArray1d<real64, NC> dVisc_dC( NC );

  for (localIndex ip = 0; ip < NP; ++ip)
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

    for (localIndex jp = 0; jp < NP; ++jp)
    {
      real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
      dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[jp];

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[jp][jc];
      }
    }

    real64 const mobility = relPerm * density / viscosity;

    phaseMob[ip] = mobility;
    dPhaseMob_dPres[ip] = dRelPerm_dP * density / viscosity
                          + mobility * (dDens_dP / density - dVisc_dP / viscosity);

    // compositional derivatives
    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] * density / viscosity
                                + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
    }
  }
}

inline RAJA_HOST_DEVICE void
UpdatePhaseMobilityKernel::Compute( localIndex NC, localIndex NP,
                                    arraySlice2d<real64 const> dCompFrac_dCompDens,
                                    arraySlice1d<real64 const> phaseDens,
                                    arraySlice1d<real64 const> dPhaseDens_dPres,
                                    arraySlice2d<real64 const> dPhaseDens_dComp,
                                    arraySlice1d<real64 const> phaseVisc,
                                    arraySlice1d<real64 const> dPhaseVisc_dPres,
                                    arraySlice2d<real64 const> dPhaseVisc_dComp,
                                    arraySlice1d<real64 const> phaseRelPerm,
                                    arraySlice2d<real64 const> dPhaseRelPerm_dPhaseVolFrac,
                                    arraySlice1d<real64 const> dPhaseVolFrac_dPres,
                                    arraySlice2d<real64 const> dPhaseVolFrac_dComp,
                                    arraySlice1d<real64> phaseMob,
                                    arraySlice1d<real64> dPhaseMob_dPres,
                                    arraySlice2d<real64> dPhaseMob_dComp )
{
  localIndex constexpr maxNC = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

  stackArray1d<real64, maxNC> dRelPerm_dC( NC );
  stackArray1d<real64, maxNC> dDens_dC( NC );
  stackArray1d<real64, maxNC> dVisc_dC( NC );

  for (localIndex ip = 0; ip < NP; ++ip)
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

    for (localIndex jp = 0; jp < NP; ++jp)
    {
      real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
      dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[jp];

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[jp][jc];
      }
    }

    real64 const mobility = relPerm * density / viscosity;

    phaseMob[ip] = mobility;
    dPhaseMob_dPres[ip] = dRelPerm_dP * density / viscosity
                          + mobility * (dDens_dP / density - dVisc_dP / viscosity);

    // compositional derivatives
    for (localIndex jc = 0; jc < NC; ++jc)
    {
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] * density / viscosity
                                + mobility * (dDens_dC[jc] / density - dVisc_dC[jc] / viscosity);
    }
  }
}

template<localIndex NC, localIndex NP>
void UpdatePhaseMobilityKernel::Launch( localIndex begin, localIndex end,
                                        arrayView3d<real64 const> const & dCompFrac_dCompDens,
                                        arrayView3d<real64 const> const & phaseDens,
                                        arrayView3d<real64 const> const & dPhaseDens_dPres,
                                        arrayView4d<real64 const> const & dPhaseDens_dComp,
                                        arrayView3d<real64 const> const & phaseVisc,
                                        arrayView3d<real64 const> const & dPhaseVisc_dPres,
                                        arrayView4d<real64 const> const & dPhaseVisc_dComp,
                                        arrayView3d<real64 const> const & phaseRelPerm,
                                        arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
                                        arrayView2d<real64 const> const & dPhaseVolFrac_dPres,
                                        arrayView3d<real64 const> const & dPhaseVolFrac_dComp,
                                        arrayView2d<real64> const & phaseMob,
                                        arrayView2d<real64> const & dPhaseMob_dPres,
                                        arrayView3d<real64> const & dPhaseMob_dComp )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
  {
    Compute<NC, NP>( dCompFrac_dCompDens[a],
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

void UpdatePhaseMobilityKernel::Launch( localIndex NC, localIndex NP,
                                        localIndex begin, localIndex end,
                                        arrayView3d<real64 const> const & dCompFrac_dCompDens,
                                        arrayView3d<real64 const> const & phaseDens,
                                        arrayView3d<real64 const> const & dPhaseDens_dPres,
                                        arrayView4d<real64 const> const & dPhaseDens_dComp,
                                        arrayView3d<real64 const> const & phaseVisc,
                                        arrayView3d<real64 const> const & dPhaseVisc_dPres,
                                        arrayView4d<real64 const> const & dPhaseVisc_dComp,
                                        arrayView3d<real64 const> const & phaseRelPerm,
                                        arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
                                        arrayView2d<real64 const> const & dPhaseVolFrac_dPres,
                                        arrayView3d<real64 const> const & dPhaseVolFrac_dComp,
                                        arrayView2d<real64> const & phaseMob,
                                        arrayView2d<real64> const & dPhaseMob_dPres,
                                        arrayView3d<real64> const & dPhaseMob_dComp )
{
  forall_in_range( begin, end, GEOSX_LAMBDA ( localIndex const a )
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

#define LAUNCH_UpdatePhaseMobilityKernel(NC,NP) \
template \
void UpdatePhaseMobilityKernel::Launch<NC,NP>( localIndex begin, localIndex end, \
                                               arrayView3d<real64 const> const & dCompFrac_dCompDens, \
                                               arrayView3d<real64 const> const & phaseDens, \
                                               arrayView3d<real64 const> const & dPhaseDens_dPres, \
                                               arrayView4d<real64 const> const & dPhaseDens_dComp, \
                                               arrayView3d<real64 const> const & phaseVisc, \
                                               arrayView3d<real64 const> const & dPhaseVisc_dPres, \
                                               arrayView4d<real64 const> const & dPhaseVisc_dComp, \
                                               arrayView3d<real64 const> const & phaseRelPerm, \
                                               arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac, \
                                               arrayView2d<real64 const> const & dPhaseVolFrac_dPres, \
                                               arrayView3d<real64 const> const & dPhaseVolFrac_dComp, \
                                               arrayView2d<real64> const & phaseMob, \
                                               arrayView2d<real64> const & dPhaseMob_dPres, \
                                               arrayView3d<real64> const & dPhaseMob_dComp )

LAUNCH_UpdatePhaseMobilityKernel(1,1);
LAUNCH_UpdatePhaseMobilityKernel(2,1);
LAUNCH_UpdatePhaseMobilityKernel(3,1);
LAUNCH_UpdatePhaseMobilityKernel(4,1);
LAUNCH_UpdatePhaseMobilityKernel(5,1);
LAUNCH_UpdatePhaseMobilityKernel(6,1);
LAUNCH_UpdatePhaseMobilityKernel(7,1);
LAUNCH_UpdatePhaseMobilityKernel(8,1);
LAUNCH_UpdatePhaseMobilityKernel(9,1);
LAUNCH_UpdatePhaseMobilityKernel(10,1);

LAUNCH_UpdatePhaseMobilityKernel(1,2);
LAUNCH_UpdatePhaseMobilityKernel(2,2);
LAUNCH_UpdatePhaseMobilityKernel(3,2);
LAUNCH_UpdatePhaseMobilityKernel(4,2);
LAUNCH_UpdatePhaseMobilityKernel(5,2);
LAUNCH_UpdatePhaseMobilityKernel(6,2);
LAUNCH_UpdatePhaseMobilityKernel(7,2);
LAUNCH_UpdatePhaseMobilityKernel(8,2);
LAUNCH_UpdatePhaseMobilityKernel(9,2);
LAUNCH_UpdatePhaseMobilityKernel(10,2);

LAUNCH_UpdatePhaseMobilityKernel(1,3);
LAUNCH_UpdatePhaseMobilityKernel(2,3);
LAUNCH_UpdatePhaseMobilityKernel(3,3);
LAUNCH_UpdatePhaseMobilityKernel(4,3);
LAUNCH_UpdatePhaseMobilityKernel(5,3);
LAUNCH_UpdatePhaseMobilityKernel(6,3);
LAUNCH_UpdatePhaseMobilityKernel(7,3);
LAUNCH_UpdatePhaseMobilityKernel(8,3);
LAUNCH_UpdatePhaseMobilityKernel(9,3);
LAUNCH_UpdatePhaseMobilityKernel(10,3);

#undef LAUNCH_UpdatePhaseMobilityKernel

} // namespace CompositionalMultiphaseFlowKernels

} // namespace geosx
