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
 * @file StatisticsKernel.hpp
 */

#ifndef GEOSX_STATISTICSKERNEL_HPP
#define GEOSX_STATISTICSKERNEL_HPP


namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

//using namespace constitutive;

/******************************** StatisticsKernel ********************************/

struct StatisticsKernel
{
  template< typename POLICY >
  static void
  saveDeltaPressure( geos::localIndex const size,
                     geos::arrayView1d< geos::real64 const > const & pres,
                     geos::arrayView1d< geos::real64 const > const & initPres,
                     geos::arrayView1d< geos::real64 > const & deltaPres )
  {
    geos::forAll< geos::parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( geos::localIndex const ei )
    {
      deltaPres[ei] = pres[ei] - initPres[ei];
    } );
  }

  template< typename POLICY >
  static void
  launch( geos::localIndex const size,
          geos::integer const numComps,
          geos::integer const numPhases,
          geos::real64 const relpermThreshold,
          geos::arrayView1d< geos::integer const > const & elemGhostRank,
          geos::arrayView1d< geos::real64 const > const & volume,
          geos::arrayView1d< geos::real64 const > const & pres,
          geos::arrayView1d< geos::real64 const > const & deltaPres,
          geos::arrayView1d< geos::real64 const > const & temp,
          geos::arrayView1d< geos::real64 const > const & refPorosity,
          geos::arrayView2d< geos::real64 const > const & porosity,
          arrayView3d< geos::real64 const, geos::constitutive::multifluid::USD_PHASE > const & phaseDensity,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & phaseCompFraction,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseTrappedVolFrac,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelperm,
          real64 & minPres,
          real64 & avgPresNumerator,
          real64 & maxPres,
          real64 & minDeltaPres,
          real64 & maxDeltaPres,
          real64 & minTemp,
          real64 & avgTempNumerator,
          real64 & maxTemp,
          real64 & totalUncompactedPoreVol,
          arrayView1d< real64 > const & phaseDynamicPoreVol,
          arrayView1d< real64 > const & phaseMass,
          arrayView1d< real64 > const & trappedPhaseMass,
          arrayView1d< real64 > const & immobilePhaseMass,
          arrayView2d< real64 > const & dissolvedComponentMass )
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinPres( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgPresNumerator( 0.0 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxPres( -LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinDeltaPres( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxDeltaPres( -LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real64 > subRegionMinTemp( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionAvgTempNumerator( 0.0 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > subRegionMaxTemp( 0.0 );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > subRegionTotalUncompactedPoreVol( 0.0 );

    // For this arrays phaseDynamicPoreVol, phaseMass, dissolvedComponentMass,
    // using an array of ReduceSum leads to a formal parameter overflow in CUDA.
    // As a workaround, we use a slice with RAJA::atomicAdd instead

    forAll< parallelDevicePolicy<> >( size, [numComps,
                                             numPhases,
                                             relpermThreshold,
                                             elemGhostRank,
                                             volume,
                                             refPorosity,
                                             porosity,
                                             pres,
                                             deltaPres,
                                             temp,
                                             phaseDensity,
                                             phaseVolFrac,
                                             phaseTrappedVolFrac,
                                             phaseRelperm,
                                             phaseCompFraction,
                                             subRegionMinPres,
                                             subRegionAvgPresNumerator,
                                             subRegionMaxPres,
                                             subRegionMinDeltaPres,
                                             subRegionMaxDeltaPres,
                                             subRegionMinTemp,
                                             subRegionAvgTempNumerator,
                                             subRegionMaxTemp,
                                             subRegionTotalUncompactedPoreVol,
                                             phaseDynamicPoreVol,
                                             phaseMass,
                                             trappedPhaseMass,
                                             immobilePhaseMass,
                                             dissolvedComponentMass] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] >= 0 )
      {
        return;
      }

      // To match our "reference", we have to use reference porosity here, not the actual porosity when we compute averages
      real64 const uncompactedPoreVol = volume[ei] * refPorosity[ei];
      real64 const dynamicPoreVol = volume[ei] * porosity[ei][0];

      subRegionMinPres.min( pres[ei] );
      subRegionAvgPresNumerator += uncompactedPoreVol * pres[ei];
      subRegionMaxPres.max( pres[ei] );

      subRegionMaxDeltaPres.max( deltaPres[ei] );
      subRegionMinDeltaPres.min( deltaPres[ei] );

      subRegionMinTemp.min( temp[ei] );
      subRegionAvgTempNumerator += uncompactedPoreVol * temp[ei];
      subRegionMaxTemp.max( temp[ei] );
      subRegionTotalUncompactedPoreVol += uncompactedPoreVol;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        real64 const elemPhaseVolume = dynamicPoreVol * phaseVolFrac[ei][ip];
        real64 const elemPhaseMass = phaseDensity[ei][0][ip] * elemPhaseVolume;
        real64 const elemTrappedPhaseMass = phaseDensity[ei][0][ip] * dynamicPoreVol * phaseTrappedVolFrac[ei][0][ip];
        // RAJA::atomicAdd used here because we do not use ReduceSum here (for the reason explained above)
        RAJA::atomicAdd( parallelDeviceAtomic{}, &phaseDynamicPoreVol[ip], elemPhaseVolume );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &phaseMass[ip], elemPhaseMass );
        RAJA::atomicAdd( parallelDeviceAtomic{}, &trappedPhaseMass[ip], elemTrappedPhaseMass );
        if( phaseRelperm[ei][0][ip] < relpermThreshold )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &immobilePhaseMass[ip], elemPhaseMass );
        }
        for( integer ic = 0; ic < numComps; ++ic )
        {
          // RAJA::atomicAdd used here because we do not use ReduceSum here (for the reason explained above)
          RAJA::atomicAdd( parallelDeviceAtomic{}, &dissolvedComponentMass[ip][ic], phaseCompFraction[ei][0][ip][ic] * elemPhaseMass );
        }
      }

    } );

    minPres = subRegionMinPres.get();
    avgPresNumerator = subRegionAvgPresNumerator.get();
    maxPres = subRegionMaxPres.get();
    minDeltaPres = subRegionMinDeltaPres.get();
    maxDeltaPres = subRegionMaxDeltaPres.get();
    minTemp = subRegionMinTemp.get();
    avgTempNumerator = subRegionAvgTempNumerator.get();
    maxTemp = subRegionMaxTemp.get();
    totalUncompactedPoreVol = subRegionTotalUncompactedPoreVol.get();

    // dummy loop to bring data back to the CPU
    forAll< serialPolicy >( 1, [phaseDynamicPoreVol, phaseMass, trappedPhaseMass, immobilePhaseMass, dissolvedComponentMass]( localIndex const )
    {
      GEOS_UNUSED_VAR( phaseDynamicPoreVol, phaseMass, trappedPhaseMass, immobilePhaseMass, dissolvedComponentMass );
    } );
  }
};

}

}

#endif //GEOSX_STATISTICSKERNEL_HPP
