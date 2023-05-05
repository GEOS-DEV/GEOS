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
 * @file ExplicitMPM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_CONTACT_EXPLICITMPM_HPP_
#define GEOSX_PHYSICSSOLVERS_CONTACT_EXPLICITMPM_HPP_

#include "physicsSolvers/solidMechanics/kernels/ExplicitFiniteStrain.hpp"
#include "physicsSolvers/solidMechanics/MPMSolverFields.hpp"

namespace geos
{

namespace solidMechanicsMPMKernels
{
/**
 * @brief A struct to update particle stresses
 */
struct StateUpdateKernel
{
  /**
   * @brief Launch the kernel function doing constitutive updates
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONTACT_WRAPPER the type of contact wrapper doing the constitutive updates
   * @param[in] size the size of the subregion
   * @param[in] constitutiveWrapper the wrapper implementing the constitutive model
   * @param[in] deformationGradient the deformation gradient
   * @param[in] velocityGradient the velocity gradient
   * @param[out] particleStress the new particle stress, returned for plotting convenience
   */
  template< typename POLICY, typename CONSTITUTIVE_WRAPPER >
  static void launch( SortedArrayView< localIndex const > const indices,
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      real64 dt,
                      arrayView3d< real64 const > const deformationGradient,
                      arrayView3d< real64 const > const fDot,
                      arrayView3d< real64 const > const velocityGradient,
                      arrayView2d< real64 > const particleStress )
  {
    arrayView3d< real64, solid::STRESS_USD > const oldStress = constitutiveWrapper.m_oldStress;
    arrayView3d< real64, solid::STRESS_USD > const newStress = constitutiveWrapper.m_newStress;
    
    // Perform constitutive call
    forAll< POLICY >( indices.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // Particle index
      localIndex const p = indices[k];

      // Copy the beginning-of-step particle stress into the constitutive model's m_oldStress - this fixes the MPI sync issue on Lassen for some reason
      #if defined(GEOSX_USE_CUDA)
      LvArray::tensorOps::copy< 6 >( oldStress[p][0] , particleStress[p] );
      #endif
      
      // Determine the strain increment in Voigt notation
      real64 strainIncrement[6];
      strainIncrement[0] = velocityGradient[p][0][0] * dt;
      strainIncrement[1] = velocityGradient[p][1][1] * dt;
      strainIncrement[2] = velocityGradient[p][2][2] * dt;
      strainIncrement[3] = (velocityGradient[p][1][2] + velocityGradient[p][2][1]) * dt;
      strainIncrement[4] = (velocityGradient[p][0][2] + velocityGradient[p][2][0]) * dt;
      strainIncrement[5] = (velocityGradient[p][0][1] + velocityGradient[p][1][0]) * dt;

      // Get old F by incrementing backwards
      real64 fOld[3][3] = { {0} };
      LvArray::tensorOps::copy< 3, 3 >( fOld, deformationGradient[p] );
      LvArray::tensorOps::scaledAdd< 3, 3 >( fOld, fDot[p], -dt );

      // Polar decompositions
      real64 rotBeginning[3][3] = { {0} };
      real64 rotEnd[3][3] = { {0} };
      LvArray::tensorOps::polarDecomposition< 3 >( rotBeginning, fOld );
      LvArray::tensorOps::polarDecomposition< 3 >( rotEnd, deformationGradient[p] );
      
      // Call stress update
      real64 stress[6] = { 0 };
      constitutiveWrapper.hypoUpdate2_StressOnly( p,                 // particle local index
                                                  0,                 // particles have 1 quadrature point
                                                  dt,                // time step size
                                                  strainIncrement,   // particle strain increment
                                                  rotBeginning,      // beginning-of-step rotation matrix
                                                  rotEnd,            // end-of-step rotation matrix
                                                  stress );          // final updated stress

      // Copy the updated stress into particleStress      
      LvArray::tensorOps::copy< 6 >( particleStress[p], stress );

      // Copy m_newStress into m_oldStress
      constitutiveWrapper.saveConvergedState( p, 0 );
    } );
  }
};


} // namespace solidMechanicsMPMKernels

} // namespace geos


#endif /* GEOSX_PHYSICSSOLVERS_CONTACT_EXPLICITMPM_HPP_ */
