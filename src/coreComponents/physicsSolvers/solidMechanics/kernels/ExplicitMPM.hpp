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

namespace geosx
{

namespace solidMechanicsMPMKernels
{
/**
 * @brief A struct to update particle stresses
 */
struct StateUpdateKernel
{

  /**
   * @brief Launch the kernel function doing fracture traction updates
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONTACT_WRAPPER the type of contact wrapper doing the fracture traction updates
   * @param[in] size the size of the subregion
   * @param[in] constitutiveWrapper the wrapper implementing the constitutive model
   * @param[in] deformationGradient the deformation gradient
   * @param[in] velocityGradient the velocity gradient
   * @param[out] particleStress the new particle stress, returned for plotting convenience
   */
  template< typename POLICY, typename CONSTITUTIVE_WRAPPER >
  static void launch( localIndex const size,
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      real64 dt,
                      arrayView3d< real64 > const & deformationGradient,
                      arrayView3d< real64 const > const & velocityGradient,
                      arrayView2d< real64 > const & particleStress )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      // Determine the strain increment (Voigt notation)
      real64 strainIncrement[6];
      strainIncrement[0] = velocityGradient[k][0][0] * dt;
      strainIncrement[1] = velocityGradient[k][1][1] * dt;
      strainIncrement[2] = velocityGradient[k][2][2] * dt;
      strainIncrement[3] = (velocityGradient[k][1][2] + velocityGradient[k][2][1]) * dt;
      strainIncrement[4] = (velocityGradient[k][0][2] + velocityGradient[k][2][0]) * dt;
      strainIncrement[5] = (velocityGradient[k][0][1] + velocityGradient[k][1][0]) * dt;

      // Perform F update
      real64 FOld[3][3] = { {0} };
      real64 dF[3][3] = { {0} };
      LvArray::tensorOps::copy< 3, 3 >( FOld, deformationGradient[k] );
      LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( dF, velocityGradient[k], FOld );
      LvArray::tensorOps::scale< 3, 3 >( dF, dt );
      LvArray::tensorOps::add< 3, 3 >( deformationGradient[k], dF );

      // Polar decompositions
      real64 rotBeginning[3][3] = { {0} };
      real64 rotEnd[3][3] = { {0} };
      
      // Call stress update
      real64 stress[6] = { 0 };
      constitutiveWrapper.hypoUpdate2_StressOnly( k,                        // particle local index
                                                  0,                        // particles have 1 quadrature point
                                                  strainIncrement,          // particle strain increment
                                                  rotBeginning,             // beginning-of-step rotation matrix
                                                  rotEnd,                   // end-of-step rotation matrix
                                                  stress );                 // final updated stress
      LvArray::tensorOps::copy< 6 >( particleStress[k], stress );

      // Copy m_newStress into m_oldStress
      constitutiveWrapper.saveConvergedState( k, 0 );
    } );
  }

};


} // namespace solidMechanicsMPMKernels

} // namespace geosx


#endif /* GEOSX_PHYSICSSOLVERS_CONTACT_EXPLICITMPM_HPP_ */
