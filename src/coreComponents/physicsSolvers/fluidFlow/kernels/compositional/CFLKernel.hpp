/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CFLKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_CFLKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_CFLKERNEL_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multifluid/MultiFluidFields.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityFields.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

/******************************** CFLKernel ********************************/

/**
 * @brief Functions to compute the CFL number using the phase volumetric outflux and the component mass outflux in each cell
 */
struct CFLKernel
{

  static constexpr real64 minPhaseMobility = 1e-12;
  static constexpr real64 minComponentFraction = 1e-12;

  template< integer NP >
  GEOS_HOST_DEVICE
  inline
  static void
  computePhaseCFL( real64 const poreVol,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
                   arraySlice1d< real64 const, constitutive::relperm::USD_RELPERM - 2 > phaseRelPerm,
                   arraySlice2d< real64 const, constitutive::relperm::USD_RELPERM_DS - 2 > dPhaseRelPerm_dPhaseVolFrac,
                   arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseVisc,
                   arraySlice1d< real64 const, compflow::USD_PHASE- 1 > phaseOutflux,
                   real64 & phaseCFLNumber );

  template< integer NC >
  GEOS_HOST_DEVICE
  inline
  static void
  computeCompCFL( real64 const poreVol,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compDens,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compFrac,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compOutflux,
                  real64 & compCFLNumber );

  template< integer NC, integer NP >
  static void
  launch( localIndex const size,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseVisc,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux,
          arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux,
          arrayView1d< real64 > const & phaseCFLNumber,
          arrayView1d< real64 > const & compCFLNumber,
          real64 & maxPhaseCFLNumber,
          real64 & maxCompCFLNumber );

};

/******************************** CFLKernel ********************************/

template< integer NP >
GEOS_HOST_DEVICE
void
CFLKernel::
  computePhaseCFL( real64 const poreVol,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseVolFrac,
                   arraySlice1d< real64 const, constitutive::relperm::USD_RELPERM - 2 > phaseRelPerm,
                   arraySlice2d< real64 const, constitutive::relperm::USD_RELPERM_DS - 2 > dPhaseRelPerm_dPhaseVolFrac,
                   arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > phaseVisc,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > phaseOutflux,
                   real64 & phaseCFLNumber )
{
  // first, check which phases are mobile in the cell
  real64 mob[NP]{};
  localIndex mobilePhases[NP]{};
  localIndex numMobilePhases = 0;
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    if( phaseVolFrac[ip] > 0 )
    {
      mob[ip] = phaseRelPerm[ip] / phaseVisc[ip];
      if( mob[ip] > minPhaseMobility )
      {
        mobilePhases[numMobilePhases] = ip;
        numMobilePhases++;
      }
    }
  }

  // then, depending on the regime, apply the appropriate CFL formula
  phaseCFLNumber = 0;

  // single-phase flow regime
  if( numMobilePhases == 1 )
  {
    phaseCFLNumber = phaseOutflux[mobilePhases[0]] / poreVol;
  }
  // two-phase flow regime
  else if( numMobilePhases == 2 )
  {
    // from Hui Cao's PhD thesis
    localIndex const ip0 = mobilePhases[0];
    localIndex const ip1 = mobilePhases[1];
    real64 const dMob_dVolFrac[2] = { dPhaseRelPerm_dPhaseVolFrac[ip0][ip0] / phaseVisc[ip0],
                                      -dPhaseRelPerm_dPhaseVolFrac[ip1][ip1] / phaseVisc[ip1] }; // using S0 = 1 - S1
    real64 const denom = 1. / ( poreVol * ( mob[ip0] + mob[ip1] ) );
    real64 const coef0 = denom * mob[ip1] / mob[ip0] * dMob_dVolFrac[ip0];
    real64 const coef1 = -denom * mob[ip0] / mob[ip1] * dMob_dVolFrac[ip1];

    phaseCFLNumber = LvArray::math::abs( coef0*phaseOutflux[ip0] + coef1*phaseOutflux[ip1] );
  }
  // three-phase flow regime
  else if( numMobilePhases == 3 )
  {
    // from Keith Coats, IMPES stability: Selection of stable timesteps (2003)
    real64 totalMob = 0.0;
    for( integer ip = 0; ip < numMobilePhases; ++ip )
    {
      totalMob += mob[ip];
    }

    real64 f[2][2]{};
    for( integer i = 0; i < 2; ++i )
    {
      for( integer j = 0; j < 2; ++j )
      {
        f[i][j]  = ( i == j )*totalMob - mob[i];
        f[i][j] /= (totalMob * mob[j]);
        real64 sum = 0;
        for( integer k = 0; k < 3; ++k )
        {
          sum += dPhaseRelPerm_dPhaseVolFrac[k][j] / phaseVisc[k]
                 * phaseOutflux[j];
        }
        f[i][j] *= sum;
      }
    }
    phaseCFLNumber = f[0][0] + f[1][1];
    phaseCFLNumber += sqrt( phaseCFLNumber*phaseCFLNumber - 4 * ( f[0][0]*f[1][1] - f[1][0]*f[0][1] ) );
    phaseCFLNumber = 0.5 * LvArray::math::abs( phaseCFLNumber ) / poreVol;
  }
}

template< integer NC >
GEOS_HOST_DEVICE
void
CFLKernel::
  computeCompCFL( real64 const poreVol,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compDens,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compFrac,
                  arraySlice1d< real64 const, compflow::USD_COMP - 1 > compOutflux,
                  real64 & compCFLNumber )
{
  compCFLNumber = 0.0;
  for( integer ic = 0; ic < NC; ++ic )
  {
    if( compFrac[ic] > minComponentFraction )
    {
      real64 const compMoles = compDens[ic] * poreVol;
      real64 const CFL = compOutflux[ic] / compMoles;
      if( CFL > compCFLNumber )
      {
        compCFLNumber = CFL;
      }
    }
  }
}

template< integer NC, integer NP >
void
CFLKernel::
  launch( localIndex const size,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView3d< real64 const, constitutive::relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, constitutive::relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const & phaseVisc,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseOutflux,
          arrayView2d< real64 const, compflow::USD_COMP > const & compOutflux,
          arrayView1d< real64 > const & phaseCFLNumber,
          arrayView1d< real64 > const & compCFLNumber,
          real64 & maxPhaseCFLNumber,
          real64 & maxCompCFLNumber )
{
  RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionPhaseCFLNumber( 0.0 );
  RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionCompCFLNumber( 0.0 );

  forAll< parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 const poreVol = volume[ei] * porosity[ei][0];

    // phase CFL number
    real64 cellPhaseCFLNumber = 0.0;
    computePhaseCFL< NP >( poreVol,
                           phaseVolFrac[ei],
                           phaseRelPerm[ei][0],
                           dPhaseRelPerm_dPhaseVolFrac[ei][0],
                           phaseVisc[ei][0],
                           phaseOutflux[ei],
                           cellPhaseCFLNumber );
    subRegionPhaseCFLNumber.max( cellPhaseCFLNumber );
    phaseCFLNumber[ei] = cellPhaseCFLNumber;

    // component CFL number
    real64 cellCompCFLNumber = 0.0;
    computeCompCFL< NC >( poreVol,
                          compDens[ei],
                          compFrac[ei],
                          compOutflux[ei],
                          cellCompCFLNumber );
    subRegionCompCFLNumber.max( cellCompCFLNumber );
    compCFLNumber[ei] = cellCompCFLNumber;
  } );

  maxPhaseCFLNumber = subRegionPhaseCFLNumber.get();
  maxCompCFLNumber = subRegionCompCFLNumber.get();
}

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_CFLKERNEL_HPP
