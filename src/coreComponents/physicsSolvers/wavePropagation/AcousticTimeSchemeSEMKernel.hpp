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
 * @file AcousticTimeSchemeSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTIMESCHEMESEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTIMESCHEMESEMKERNEL_HPP_

namespace geos
{

struct AcousticTimeSchemeSEM
{

  using EXEC_POLICY = parallelDevicePolicy< >;

  static void LeapFrogWithoutPML( real64 const dt,
                                  arrayView1d< real32 > const p_np1,
                                  arrayView1d< real32 > const p_n,
                                  arrayView1d< real32 > const p_nm1,
                                  arrayView1d< real32 const > const mass,
                                  arrayView1d< real32 > const stiffnessVector,
                                  arrayView1d< real32 const > const damping,
                                  arrayView1d< real32 > const rhs,
                                  arrayView1d< localIndex const > const freeSurfaceNodeIndicator,
                                  SortedArrayView< localIndex const > const solverTargetNodesSet )
  {
    real64 const dt2 = pow( dt, 2 );
    forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = solverTargetNodesSet[n];
      if( freeSurfaceNodeIndicator[a] != 1 )
      {
        p_np1[a] = p_n[a];
        p_np1[a] *= 2.0 * mass[a];
        p_np1[a] -= (mass[a] - 0.5 * dt * damping[a]) * p_nm1[a];
        p_np1[a] += dt2 * (rhs[a] - stiffnessVector[a]);
        p_np1[a] /= mass[a] + 0.5 * dt * damping[a];
      }
    } );

  };



};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTIMESCHEMESEMKERNEL_HPP_
