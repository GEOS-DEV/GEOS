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
 * @file FlowSolverBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace flowSolverBaseKernels
{

/// Threshold for the min pore volume (below, a warning is issued)
static constexpr real64 poreVolumeThreshold = 1e-4;


/**
 * @struct MinimumPoreVolumeKernel
 * @brief Kernel to compute the min pore volume in a subRegion
 */
struct MinimumPoreVolumeKernel
{

  /*
   * @brief Kernel computing the min pore volume
   * @param[in] size the number of elements in the subRegion
   * @param[in] porosity the current element porosity
   * @param[in] volume the current element volume
   * @param[out] minPoreVolumeInSubRegion the min pore volume
   * @param[out] numElemsBelowThresholdInSubRegion the number of elements is below the threshold
   */
  inline static void
  computeMinimumPoreVolume( localIndex const size,
                            arrayView2d< real64 const > const & porosity,
                            arrayView1d< real64 const > const & volume,
                            real64 & minPoreVolumeInSubRegion,
                            localIndex & numElemsBelowThresholdInSubRegion )
  {
    RAJA::ReduceMin< parallelDeviceReduce, real64 > minPoreVolume( LvArray::NumericLimits< real64 >::max );
    RAJA::ReduceSum< parallelDeviceReduce, localIndex > numElemsBelowThreshold( 0.0 );

    real64 const pvThreshold = poreVolumeThreshold;

    forAll< parallelDevicePolicy<> >( size, [porosity,
                                             volume,
                                             pvThreshold,
                                             minPoreVolume,
                                             numElemsBelowThreshold] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      real64 const poreVolume = porosity[ei][0] * volume[ei];
      if( poreVolume < pvThreshold )
      {
        numElemsBelowThreshold += 1;
      }
      minPoreVolume.min( poreVolume );
    } );

    minPoreVolumeInSubRegion = minPoreVolume.get();
    numElemsBelowThresholdInSubRegion = numElemsBelowThreshold.get();
  }
};


} // namespace flowSolverBaseKernels

} // namespace geos

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP
