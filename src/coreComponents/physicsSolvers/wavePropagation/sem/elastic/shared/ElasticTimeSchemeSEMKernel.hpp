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
 * @file ElasticTimeSchemeSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICTIMESCHEMESEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICTIMESCHEMESEMKERNEL_HPP_

namespace geos
{

struct ElasticTimeSchemeSEM
{

  using EXEC_POLICY = parallelDevicePolicy< >;
  /**
   * @brief  Apply second order Leap-Frog time scheme for isotropic case without PML
   * @param[in] dt time-step
   * @param[out] ux_np1 displacement in x-direction array at time n+1 (updated here)
   * @param[in] ux_n displacement in x-direction array at time n
   * @param[in] ux_nm1 displacement in x-direction array at time n-1
   * @param[out] uy_np1 displacement in y-direction array at time n+1 (updated here)
   * @param[in] uy_n displacement in y-direction array at time n
   * @param[in] uy_nm1 displacement in y-direction array at time n-1
   * @param[out] uz_np1 displacement in z-direction array at time n+1 (updated here)
   * @param[in] uz_n displacement in z-direction array at time n
   * @param[in] uz_nm1 displacement in z-direction array at time n-1
   * @param[in] mass the mass matrix
   * @param[in] dampingx the damping matrix for x-component
   * @param[in] dampingy the damping matrix for y-component
   * @param[in] dampingz the damping matrix for z-component
   * @param[in] stiffnessVectorx array containing the product of the stiffness matrix R and the displacement in x-direction at time n
   * @param[in] stiffnessVectory array containing the product of the stiffness matrix R and the displacement in y-direction at time n
   * @param[in] stiffnessVectorz array containing the product of the stiffness matrix R and the displacement in z-direction at time n
   * @param[in] rhsx the right-hand-side for displacement in x-direction
   * @param[in] rhsy the right-hand-side for displacement in y-direction
   * @param[in] rhsz the right-hand-side for displacement in z-direction
   * @param[in] solverTargetNodesSet the targetted nodeset (useful in particular when we do elasto-acoustic simulation )
   */
  static void LeapFrog( real64 const dt,
                        arrayView1d< real32 > const ux_np1,
                        arrayView1d< real32 > const ux_n,
                        arrayView1d< real32 > const ux_nm1,
                        arrayView1d< real32 > const uy_np1,
                        arrayView1d< real32 > const uy_n,
                        arrayView1d< real32 > const uy_nm1,
                        arrayView1d< real32 > const uz_np1,
                        arrayView1d< real32 > const uz_n,
                        arrayView1d< real32 > const uz_nm1,
                        arrayView1d< real32 const > const mass,
                        arrayView1d< real32 const > const dampingx,
                        arrayView1d< real32 const > const dampingy,
                        arrayView1d< real32 const > const dampingz,
                        arrayView1d< real32 > const stiffnessVectorx,
                        arrayView1d< real32 > const stiffnessVectory,
                        arrayView1d< real32 > const stiffnessVectorz,
                        arrayView1d< real32 > const rhsx,
                        arrayView1d< real32 > const rhsy,
                        arrayView1d< real32 > const rhsz,
                        SortedArrayView< localIndex const > const solverTargetNodesSet )
  {
    real64 const dt2 = pow( dt, 2 );
    forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = solverTargetNodesSet[n];
      ux_np1[a] = ux_n[a];
      ux_np1[a] *= 2.0*mass[a];
      ux_np1[a] -= (mass[a]-0.5*dt*dampingx[a])*ux_nm1[a];
      ux_np1[a] += dt2*(rhsx[a]-stiffnessVectorx[a]);
      ux_np1[a] /= mass[a]+0.5*dt*dampingx[a];
      uy_np1[a] = uy_n[a];
      uy_np1[a] *= 2.0*mass[a];
      uy_np1[a] -= (mass[a]-0.5*dt*dampingy[a])*uy_nm1[a];
      uy_np1[a] += dt2*(rhsy[a]-stiffnessVectory[a]);
      uy_np1[a] /= mass[a]+0.5*dt*dampingy[a];
      uz_np1[a] = uz_n[a];
      uz_np1[a] *= 2.0*mass[a];
      uz_np1[a] -= (mass[a]-0.5*dt*dampingz[a])*uz_nm1[a];
      uz_np1[a] += dt2*(rhsz[a]-stiffnessVectorz[a]);
      uz_np1[a] /= mass[a]+0.5*dt*dampingz[a];
    } );

  };

};

}

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICTIMESCHEMESEMKERNEL_HPP_