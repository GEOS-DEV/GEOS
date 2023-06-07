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

#include "physicsSolvers/SolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/fields/FlowSolverBaseFields.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** StatisticsKernel ********************************/

struct StatisticsKernel
{
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

  static void
  launch( geos::localIndex const size,
          geos::arrayView1d< geos::integer const > const & elemGhostRank,
          geos::arrayView1d< geos::real64 const > const & volume,
          geos::arrayView1d< geos::real64 const > const & pres,
          geos::arrayView1d< geos::real64 const > const & deltaPres,
          geos::arrayView1d< geos::real64 const > const & refPorosity,
          geos::arrayView2d< geos::real64 const > const & porosity,
          geos::real64 & minPres,
          geos::real64 & avgPresNumerator,
          geos::real64 & maxPres,
          geos::real64 & minDeltaPres,
          geos::real64 & maxDeltaPres,
          geos::real64 & totalUncompactedPoreVol,
          geos::real64 & totalPoreVol )
  {
    RAJA::ReduceMin< geos::parallelDeviceReduce, geos::real64 > subRegionMinPres( LvArray::NumericLimits< geos::real64 >::max );
    RAJA::ReduceSum< geos::parallelDeviceReduce, geos::real64 > subRegionAvgPresNumerator( 0.0 );
    RAJA::ReduceMax< geos::parallelDeviceReduce, geos::real64 > subRegionMaxPres( -LvArray::NumericLimits< geos::real64 >::max );

    RAJA::ReduceMin< geos::parallelDeviceReduce, geos::real64 > subRegionMinDeltaPres( LvArray::NumericLimits< geos::real64 >::max );
    RAJA::ReduceMax< geos::parallelDeviceReduce, geos::real64 > subRegionMaxDeltaPres( -LvArray::NumericLimits< geos::real64 >::max );

    RAJA::ReduceSum< geos::parallelDeviceReduce, geos::real64 > subRegionTotalUncompactedPoreVol( 0.0 );
    RAJA::ReduceSum< geos::parallelDeviceReduce, geos::real64 > subRegionTotalPoreVol( 0.0 );

    geos::forAll< geos::parallelDevicePolicy<> >( size, [=] GEOS_HOST_DEVICE ( geos::localIndex const ei )
    {
      if( elemGhostRank[ei] >= 0 )
      {
        return;
      }

      // To match our "reference", we have to use reference porosity here, not the actual porosity when we compute averages
      geos::real64 const uncompactedPoreVol = volume[ei] * refPorosity[ei];
      geos::real64 const dynamicPoreVol = volume[ei] * porosity[ei][0];

      subRegionMinPres.min( pres[ei] );
      subRegionAvgPresNumerator += uncompactedPoreVol * pres[ei];
      subRegionMaxPres.max( pres[ei] );

      subRegionMinDeltaPres.min( deltaPres[ei] );
      subRegionMaxDeltaPres.max( deltaPres[ei] );

      subRegionTotalUncompactedPoreVol += uncompactedPoreVol;
      subRegionTotalPoreVol += dynamicPoreVol;
    } );

    minPres = subRegionMinPres.get();
    avgPresNumerator = subRegionAvgPresNumerator.get();
    maxPres = subRegionMaxPres.get();

    minDeltaPres = subRegionMinDeltaPres.get();
    maxDeltaPres = subRegionMaxDeltaPres.get();

    totalUncompactedPoreVol = subRegionTotalUncompactedPoreVol.get();
    totalPoreVol = subRegionTotalPoreVol.get();
  }
};

}

}

#endif //GEOSX_STATISTICSKERNEL_HPP
