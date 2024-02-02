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

  static void LeapFrogforVTI( NodeManager & nodeManager,
                              real64 const dt,
                              arrayView1d< real32 > const p_np1,
                              arrayView1d< real32 > const p_n,
                              arrayView1d< real32 > const p_nm1,
                              arrayView1d< real32 > const q_np1,
                              arrayView1d< real32 > const q_n,
                              arrayView1d< real32 > const q_nm1,
                              arrayView1d< real32 const > const mass,
                              arrayView1d< real32 > const stiffnessVector_p,
                              arrayView1d< real32 > const stiffnessVector_q,
                              arrayView1d< real32 const > const damping_p,
                              arrayView1d< real32 const > const damping_pq,
                              arrayView1d< real32 const > const damping_q,
                              arrayView1d< real32 const > const damping_qp,
                              arrayView1d< real32 > const rhs,
                              arrayView1d< localIndex const > const freeSurfaceNodeIndicator,
                              arrayView1d< localIndex const > const lateralSurfaceNodeIndicator,
                              arrayView1d< localIndex const > const bottomSurfaceNodeIndicator)
                                
  {
    real64 const dt2 = pow( dt, 2 );
    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      if( freeSurfaceNodeIndicator[a] != 1 )
      {
        p_np1[a] = 2.0*mass[a]*p_n[a]/dt2;
        p_np1[a] -= mass[a]*p_nm1[a]/dt2;
        p_np1[a] += stiffnessVector_p[a];
        p_np1[a] += rhs[a];

        q_np1[a] = 2.0*mass[a]*q_n[a]/dt2;
        q_np1[a] -= mass[a]*q_nm1[a]/dt2;
        q_np1[a] += stiffnessVector_q[a];
        q_np1[a] += rhs[a];

        if( lateralSurfaceNodeIndicator[a] != 1 && bottomSurfaceNodeIndicator[a] != 1 )
        {
          // Interior node, no boundary terms
          p_np1[a] /= mass[a]/dt2;
          q_np1[a] /= mass[a]/dt2;
        }
        else
        {
          // Boundary node
          p_np1[a] += damping_p[a]*p_nm1[a]/dt/2;
          p_np1[a] += damping_pq[a]*q_nm1[a]/dt/2;

          q_np1[a] += damping_q[a]*q_nm1[a]/dt/2;
          q_np1[a] += damping_qp[a]*p_nm1[a]/dt/2;
          // Hand-made Inversion of 2x2 matrix
          real32 coef_pp = mass[a]/dt2;
          coef_pp += damping_p[a]/dt/2;
          real32 coef_pq = damping_pq[a]/dt/2;

          real32 coef_qq = mass[a]/dt2;
          coef_qq += damping_q[a]/2/dt;
          real32 coef_qp = damping_qp[a]/dt/2;

          real32 det_pq = 1/(coef_pp * coef_qq - coef_pq*coef_qp);

          real32 aux_p_np1 = p_np1[a];
          p_np1[a] = det_pq*(coef_qq*p_np1[a] - coef_pq*q_np1[a]);
          q_np1[a] = det_pq*(coef_pp*q_np1[a] - coef_qp*aux_p_np1);
        }
      }
    } );
  };

  
                              

};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICTIMESCHEMESEMKERNEL_HPP_
