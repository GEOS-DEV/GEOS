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
 * @file CompositionalMultiphaseFlowKernels.hpp
 */

#ifndef GEOSX_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP
#define GEOSX_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP

#include "constitutive/Fluid/MultiFluidBase.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFlowKernels
{

struct UpdateComponentFractionKernel
{

  template<localIndex NC>
  static inline RAJA_HOST_DEVICE void Compute( arraySlice1d<real64 const> const & compDens,
                                               arraySlice1d<real64 const> const & dCompDens,
                                               arraySlice1d<real64> const & compFrac,
                                               arraySlice2d<real64> const & dCompFrac_dCompDens )
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

  static inline RAJA_HOST_DEVICE void Compute( localIndex NC,
                                               arraySlice1d<real64 const> const & compDens,
                                               arraySlice1d<real64 const> const & dCompDens,
                                               arraySlice1d<real64> const & compFrac,
                                               arraySlice2d<real64> const & dCompFrac_dCompDens )
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
  static void Launch( localIndex begin, localIndex end,
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

  static void Launch( localIndex NC,
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

};

struct UpdatePhaseVolumeFractionKernel
{

  template<localIndex NC, localIndex NP>
  static inline RAJA_HOST_DEVICE void Compute( arraySlice1d<real64 const> const & compDens,
                                               arraySlice1d<real64 const> const & dCompDens,
                                               arraySlice2d<real64 const> const & dCompFrac_dCompDens,
                                               arraySlice1d<real64 const> const & phaseDens,
                                               arraySlice1d<real64 const> const & dPhaseDens_dPres,
                                               arraySlice2d<real64 const> const & dPhaseDens_dComp,
                                               arraySlice1d<real64 const> const & phaseFrac,
                                               arraySlice1d<real64 const> const & dPhaseFrac_dPres,
                                               arraySlice2d<real64 const> const & dPhaseFrac_dComp,
                                               arraySlice1d<real64> const & phaseVolFrac,
                                               arraySlice1d<real64> const & dPhaseVolFrac_dPres,
                                               arraySlice2d<real64> const & dPhaseVolFrac_dComp )
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

  static inline RAJA_HOST_DEVICE void Compute( localIndex NC, localIndex NP,
                                               arraySlice1d<real64 const> const & compDens,
                                               arraySlice1d<real64 const> const & dCompDens,
                                               arraySlice2d<real64 const> const & dCompFrac_dCompDens,
                                               arraySlice1d<real64 const> const & phaseDens,
                                               arraySlice1d<real64 const> const & dPhaseDens_dPres,
                                               arraySlice2d<real64 const> const & dPhaseDens_dComp,
                                               arraySlice1d<real64 const> const & phaseFrac,
                                               arraySlice1d<real64 const> const & dPhaseFrac_dPres,
                                               arraySlice2d<real64 const> const & dPhaseFrac_dComp,
                                               arraySlice1d<real64> const & phaseVolFrac,
                                               arraySlice1d<real64> const & dPhaseVolFrac_dPres,
                                               arraySlice2d<real64> const & dPhaseVolFrac_dComp )
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
  static void Launch( localIndex begin, localIndex end,
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

  static void Launch( localIndex NC, localIndex NP,
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

};

struct UpdatePhaseMobilityKernel
{

  template<localIndex NC, localIndex NP>
  static inline RAJA_HOST_DEVICE void Compute( arraySlice2d<real64 const> const & dCompFrac_dCompDens,
                                               arraySlice1d<real64 const> const & phaseDens,
                                               arraySlice1d<real64 const> const & dPhaseDens_dPres,
                                               arraySlice2d<real64 const> const & dPhaseDens_dComp,
                                               arraySlice1d<real64 const> const & phaseVisc,
                                               arraySlice1d<real64 const> const & dPhaseVisc_dPres,
                                               arraySlice2d<real64 const> const & dPhaseVisc_dComp,
                                               arraySlice1d<real64 const> const & phaseRelPerm,
                                               arraySlice2d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
                                               arraySlice1d<real64 const> const & dPhaseVolFrac_dPres,
                                               arraySlice2d<real64 const> const & dPhaseVolFrac_dComp,
                                               arraySlice1d<real64> const & phaseMob,
                                               arraySlice1d<real64> const & dPhaseMob_dPres,
                                               arraySlice2d<real64> const & dPhaseMob_dComp )
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

  static inline RAJA_HOST_DEVICE void Compute( localIndex NC, localIndex NP,
                                               arraySlice2d<real64 const> const & dCompFrac_dCompDens,
                                               arraySlice1d<real64 const> const & phaseDens,
                                               arraySlice1d<real64 const> const & dPhaseDens_dPres,
                                               arraySlice2d<real64 const> const & dPhaseDens_dComp,
                                               arraySlice1d<real64 const> const & phaseVisc,
                                               arraySlice1d<real64 const> const & dPhaseVisc_dPres,
                                               arraySlice2d<real64 const> const & dPhaseVisc_dComp,
                                               arraySlice1d<real64 const> const & phaseRelPerm,
                                               arraySlice2d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
                                               arraySlice1d<real64 const> const & dPhaseVolFrac_dPres,
                                               arraySlice2d<real64 const> const & dPhaseVolFrac_dComp,
                                               arraySlice1d<real64> const & phaseMob,
                                               arraySlice1d<real64> const & dPhaseMob_dPres,
                                               arraySlice2d<real64> const & dPhaseMob_dComp )
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
  static void Launch( localIndex begin, localIndex end,
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

  static void Launch( localIndex NC, localIndex NP,
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

};

template<typename T, typename LAMBDA>
bool KernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral<T>::value, "KernelLaunchSelectorCompSwitch: type should be integral" );

  switch (value)
  {
    case 1:  lambda( std::integral_constant<T, 1>() );  return true;
    case 2:  lambda( std::integral_constant<T, 2>() );  return true;
    case 3:  lambda( std::integral_constant<T, 3>() );  return true;
    case 4:  lambda( std::integral_constant<T, 4>() );  return true;
    case 5:  lambda( std::integral_constant<T, 5>() );  return true;
    case 6:  lambda( std::integral_constant<T, 6>() );  return true;
    case 7:  lambda( std::integral_constant<T, 7>() );  return true;
    case 8:  lambda( std::integral_constant<T, 8>() );  return true;
    case 9:  lambda( std::integral_constant<T, 9>() );  return true;
    case 10: lambda( std::integral_constant<T, 10>() ); return true;
    default: return false;
  }
}

template<typename T, typename LAMBDA>
bool KernelLaunchSelectorPhaseSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral<T>::value, "KernelLaunchSelectorPhaseSwitch: type should be integral" );

  switch (value)
  {
    case 1:  lambda( std::integral_constant<T, 1>() ); return true;
    case 2:  lambda( std::integral_constant<T, 2>() ); return true;
    case 3:  lambda( std::integral_constant<T, 3>() ); return true;
    default: return false;
  }
}

template<typename KERNELWRAPPER, typename... ARGS>
void KernelLaunchSelector1( localIndex numComp, ARGS && ... args )
{
  bool const run =
  KernelLaunchSelectorCompSwitch( numComp, [&] (auto NC)
  {
    KERNELWRAPPER::template Launch<NC()>( std::forward<ARGS>(args)... );
  });

  if (!run)
  {
    KERNELWRAPPER::Launch( numComp, std::forward<ARGS>(args)... );
  }
}

template<typename KERNELWRAPPER, typename... ARGS>
void KernelLaunchSelector2( localIndex numComp, localIndex numPhase, ARGS && ... args )
{
  bool run2 = false;
  // gcc-7 produces bugged code without explicit capture list here...
  bool const run =
  KernelLaunchSelectorCompSwitch( numComp, [ &numPhase, &run2, &args... ] ( auto NC )
  {
    run2 =
    KernelLaunchSelectorPhaseSwitch( numPhase, [&] ( auto NP )
    {
      // damn you stupid C++ rules (https://stackoverflow.com/questions/43665610)
      auto constexpr NC_ = decltype(NC)::value;
      KERNELWRAPPER::template Launch<NC_, NP()>( std::forward<ARGS>(args)... );
    });
  });

  if (!run || !run2)
  {
    KERNELWRAPPER::Launch( numComp, numPhase, std::forward<ARGS>(args)... );
  }
}

} // namespace CompositionalMultiphaseFlowKernels

} // namespace geosx


#endif //GEOSX_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP
