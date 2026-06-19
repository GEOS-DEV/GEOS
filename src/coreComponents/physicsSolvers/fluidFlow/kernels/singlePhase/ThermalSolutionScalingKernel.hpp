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
 * @file ThermalSolutionScalingKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALSOLUTIONSCALINGKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALSOLUTIONSCALINGKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace thermalSinglePhaseBaseKernels
{

/******************************** SolutionScalingKernel ********************************/

struct SolutionScalingKernel
{
  template< typename POLICY >
  static std::tuple< real64, real64, real64, real64, real64 >  launch( arrayView1d< real64 const > const & localSolution,
                                                                       globalIndex const rankOffset,
                                                                       globalIndex const temperatureOffset,
                                                                       arrayView1d< globalIndex const > const & dofNumber,
                                                                       arrayView1d< integer const > const & ghostRank,
                                                                       real64 const maxAbsolutePresChange,
                                                                       real64 const maxAbsoluteTempChange,
                                                                       arrayView1d< real64 > pressureScalingFactor,
                                                                       arrayView1d< real64 > temperatureScalingFactor )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > scalingFactor( 1.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaPres( 0.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaTemp( 0.0 );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > localMinPresScalingFactor( 1.0 );
    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > localMinTempScalingFactor( 1.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei ) mutable
    {
      if( ghostRank[ei] < 0 && dofNumber[ei] >= 0 )
      {
        pressureScalingFactor[ei] = 1.0;
        temperatureScalingFactor[ei] = 1.0;

        localIndex const lid = dofNumber[ei] - rankOffset;

        // compute the change in pressure
        real64 const absPresChange = LvArray::math::abs( localSolution[lid] );
        maxDeltaPres.max( absPresChange );

        // compute the change in temperature
        real64 const absTempChange = LvArray::math::abs( localSolution[lid + temperatureOffset] );
        maxDeltaTemp.max( absTempChange );

        // maxAbsolutePresChange <= 0.0 means that scaling is disabled, and we are only collecting maxDeltaPres in this kernel
        if( maxAbsolutePresChange > 0.0 && absPresChange > maxAbsolutePresChange )
        {
          real64 const presScalingFactor = maxAbsolutePresChange / absPresChange;
          pressureScalingFactor[ei] = presScalingFactor;
          scalingFactor.min( presScalingFactor );

          localMinPresScalingFactor.min( presScalingFactor );
        }

        // maxAbsoluteTempChange <= 0.0 means that scaling is disabled, and we are only collecting maxDeltaTemps in this kernel
        if( maxAbsoluteTempChange > 0.0 && absTempChange > maxAbsoluteTempChange )
        {
          real64 const tempScalingFactor = maxAbsoluteTempChange / absTempChange;
          temperatureScalingFactor[ei] = tempScalingFactor;
          scalingFactor.min( tempScalingFactor );

          localMinTempScalingFactor.min( tempScalingFactor );
        }
      }

    } );

    return { scalingFactor.get(), maxDeltaPres.get(), maxDeltaTemp.get(), localMinPresScalingFactor.get(), localMinTempScalingFactor.get() };
  }

};

} // namespace thermalSinglePhaseBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_THERMALSOLUTIONSCALINGKERNEL_HPP
