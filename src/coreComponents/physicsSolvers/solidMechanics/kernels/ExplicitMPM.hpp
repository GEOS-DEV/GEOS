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
 * @file ExplicitMPM.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_EXPLICITMPM_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_EXPLICITMPM_HPP_

#include "constitutive/solid/SolidUtilities.hpp"
#include "physicsSolvers/solidMechanics/kernels/ExplicitFiniteStrain.hpp"
#include "physicsSolvers/solidMechanics/MPMSolverFields.hpp"

namespace geos
{

namespace solidMechanicsMPMKernels
{

using namespace constitutive;

/**
 * @brief A struct to update particle stresses
 */
struct ParticleStateUpdateKernel
{
  /**
   * @brief Launch the kernel function doing constitutive updates
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONSTITUTIVE_WRAPPER the type of consitutive wrapper doing the constitutive updates
   * @param[in] dt The time step
   * @param[in] hyperelasticUpdate Flag to perform hyperelastic update (constitutive model dependent)
   * @param[in] deformationGradient The current/end-of-step particle deformation gradient F_{n+1}
   * @param[in] fDot The step-averaged time derivative of the deformation gradient, used to recover F_n
   * @param[in] velocityGradient The step velocity gradient used to build the strain increment over [t_n,t_{n+1}]
   * @param[out] particleStress The new particle stress, returned for plotting convenience
   */
  template< typename POLICY, typename CONSTITUTIVE_WRAPPER >
  static void launch( SortedArrayView< localIndex const > const indices,
                      localIndex const batchSize,
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      arrayView3d< real64 const > const deformationGradient,
                      arrayView2d< real64 > const particleStress )
  {
    arrayView3d< real64, solid::STRESS_USD > const oldStress = constitutiveWrapper.m_oldStress;
    arrayView3d< real64, solid::STRESS_USD > const newStress = constitutiveWrapper.m_newStress;

    if( indices.size() == 0 )
    {
      return;
    }

    for( localIndex begin = 0; begin < indices.size(); begin += batchSize )
    {
      localIndex const count = LvArray::math::min( batchSize, indices.size() - begin );

      forAll< POLICY >( count, [=] GEOS_HOST_DEVICE ( localIndex const q )
      {
        localIndex const k = begin + q;
        localIndex const p = indices[k];

        real64 stress[6] = {};

        real64 FminusI[3][3] = {};
        LvArray::tensorOps::copy< 3, 3 >( FminusI, deformationGradient[p] );
        LvArray::tensorOps::addIdentity< 3 >( FminusI, -1.0 );

        constitutiveWrapper.hyperUpdate( p,       // particle local index
                                         0,       // particles have 1 quadrature point
                                         FminusI, // particle strain increment
                                         stress );

        // Copy the updated stress into particleStress
        LvArray::tensorOps::copy< 6 >( particleStress[p], stress );

        // Copy m_newStress into m_oldStress
        constitutiveWrapper.saveConvergedState( p, 0 );
    
      } );
    }
  }


  /**
   * @brief Launch the kernel function doing constitutive updates
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONSTITUTIVE_WRAPPER the type of consitutive wrapper doing the constitutive updates
   * @param[in] dt The time step
   * @param[in] hyperelasticUpdate Flag to perform hyperelastic update (constitutive model dependent)
   * @param[in] deformationGradient The current/end-of-step particle deformation gradient F_{n+1}
   * @param[in] fDot The step-averaged time derivative of the deformation gradient, used to recover F_n
   * @param[in] velocityGradient The step velocity gradient used to build the strain increment over [t_n,t_{n+1}]
   * @param[out] particleStress The new particle stress, returned for plotting convenience
   */
  template< typename POLICY, typename CONSTITUTIVE_WRAPPER >
  static void launch( SortedArrayView< localIndex const > const indices,
                      localIndex const batchSize,
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      real64 dt,
                      arrayView3d< real64 const > const deformationGradient,
                      arrayView3d< real64 const > const rotation,
                      arrayView3d< real64 const > const oldRotation,
                      arrayView3d< real64 const > const velocityGradient,
                      arrayView2d< real64 > const particleStress )
  {
    GEOS_UNUSED_VAR( deformationGradient );

    arrayView3d< real64, solid::STRESS_USD > const oldStress = constitutiveWrapper.m_oldStress;
    arrayView3d< real64, solid::STRESS_USD > const newStress = constitutiveWrapper.m_newStress;

    if( indices.size() == 0 )
    {
      return;
    }

    for( localIndex begin = 0; begin < indices.size(); begin += batchSize )
    {
      localIndex const count = LvArray::math::min( batchSize, indices.size() - begin );

      forAll< POLICY >( count, [=] GEOS_HOST_DEVICE ( localIndex const q )
      {
        localIndex const k = begin + q;
        localIndex const p = indices[k];

        real64 stress[6] = {};
        
        // Hypoeleastic stress update
        // Determine the strain increment in Voigt notation
        real64 strainIncrement[6] = {};
        strainIncrement[0] = velocityGradient[p][0][0] * dt;
        strainIncrement[1] = velocityGradient[p][1][1] * dt;
        strainIncrement[2] = velocityGradient[p][2][2] * dt;
        strainIncrement[3] = (velocityGradient[p][1][2] + velocityGradient[p][2][1]) * dt;
        strainIncrement[4] = (velocityGradient[p][0][2] + velocityGradient[p][2][0]) * dt;
        strainIncrement[5] = (velocityGradient[p][0][1] + velocityGradient[p][1][0]) * dt;

        real64 rotBeginning[3][3] = {};
        real64 rotEnd[3][3] = {};
        LvArray::tensorOps::copy< 3, 3 >( rotBeginning, oldRotation[p] );
        LvArray::tensorOps::copy< 3, 3 >( rotEnd, rotation[p] );

        constitutive::SolidUtilities::hypoUpdate2_StressOnly( constitutiveWrapper,
                                                              p,
                                                              0,
                                                              dt,
                                                              strainIncrement,
                                                              rotBeginning,
                                                              rotEnd,
                                                              stress );

        // Copy the updated stress into particleStress
        LvArray::tensorOps::copy< 6 >( particleStress[p], stress );

        // Copy m_newStress into m_oldStress
        constitutiveWrapper.saveConvergedState( p, 0 );
    
      } );
    }
  }
};

/**
 * @brief A struct to update cohesive zones
 */
struct CohesiveZoneStateUpdateKernel
{
  /**
   * @brief Launch the kernel function doing constitutive updates
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONSTITUTIVE_WRAPPER the type of consitutive wrapper doing the constitutive updates
   * @param[in] dt The time step
   * @param[in] planeStrain Flag for plane strain problems
   * @param[in] smallMass Minimum grid mass
   * @param[in] preventCZInterpentration Flag to prevent interpentration of cohesive zones in compression
   * @param[in] pairToFieldSlot Compact field slots on the two sides of every cohesive pair
   * @param[in] periodic0 Periodic flag for x-direction
   * @param[in] periodic1 Periodic flag for y-direction
   * @param[in] periodic2 Periodic flag for z-direction
   * @param[in] domainExtent0 Periodic domain extent in x-direction
   * @param[in] domainExtent1 Periodic domain extent in y-direction
   * @param[in] domainExtent2 Periodic domain extent in z-direction
   * @param[in] fieldSlotMass Mass stored once for every active grid-node/velocity-field slot
   * @param[in] fieldSlotDisplacement Displacement stored once for every active field slot
   * @param[in] fieldSlotParticleSurfaceNormal Particle-mapped surface normal for every active field slot
   * @param[in] fieldSlotDeformationGradientCofactor Mapped deformation-gradient cofactor for every active field slot
   * @param[out] fieldSlotCohesiveForce Cohesive forces accumulated into the compact field slots
   * @param[in] fieldSlotReferenceSurfaceNormal Reference surface normal for every active field slot
   * @param[in] fieldSlotReferenceArea Reference area for every active field slot
   
   */
  template< typename POLICY, typename CONSTITUTIVE_WRAPPER >
  static void launch( int numPairs,
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      real64 dt,
                      int planeStrain, // Should remove eventually, used for normals but normals after mapping to grid should be checked for planeStrain condition
                      real64 smallMass,
                      int preventCZInterpentration,
                      arrayView2d< localIndex const > const pairToFieldSlot,
                      int periodic0,
                      int periodic1,
                      int periodic2,
                      real64 domainExtent0,
                      real64 domainExtent1,
                      real64 domainExtent2,
                      arrayView1d< real64 const > const fieldSlotMass,
                      arrayView2d< real64 const > const fieldSlotDisplacement,
                      arrayView2d< real64 const > const fieldSlotParticleSurfaceNormal,
                      arrayView3d< real64 const > const fieldSlotDeformationGradientCofactor,
                      arrayView2d< real64 > const fieldSlotCohesiveForce,
                      arrayView2d< real64 const > const fieldSlotReferenceSurfaceNormal,
                      arrayView1d< real64 const > const fieldSlotReferenceArea )
  {
    GEOS_UNUSED_VAR( dt );

    // Perform constitutive call
    forAll< POLICY >( numPairs, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      localIndex const slotA = pairToFieldSlot[k][0];
      localIndex const slotB = pairToFieldSlot[k][1];

      bool active = ( fieldSlotMass[slotA] > smallMass ) &&
                    ( LvArray::tensorOps::l2NormSquared< 3 >( fieldSlotParticleSurfaceNormal[slotA] ) > 1.0e-16 ) &&
                    ( fieldSlotMass[slotB] > smallMass ) &&
                    ( LvArray::tensorOps::l2NormSquared< 3 >( fieldSlotParticleSurfaceNormal[slotB] ) > 1.0e-16 );

      if( active )
      {
        // Copy normals
        real64 nA[3] = {};
        real64 nB[3] = {};
        LvArray::tensorOps::copy< 3 >( nA, fieldSlotParticleSurfaceNormal[slotA] );
        LvArray::tensorOps::copy< 3 >( nB, fieldSlotParticleSurfaceNormal[slotB] );

        // Initialize tractions here
        real64 tA[3] = {};
        real64 tB[3] = {};
        LvArray::tensorOps::fill< 3 >( tA, 0.0 );
        LvArray::tensorOps::fill< 3 >( tB, 0.0 );

        // Compute updated nodal area vectors for cohesive zone traction calculations
        real64 referenceAreaVectorA[3] = {};
        LvArray::tensorOps::scaledCopy< 3 >( referenceAreaVectorA,
                                             fieldSlotReferenceSurfaceNormal[slotA],
                                             fieldSlotReferenceArea[slotA] );

        real64 sA[3] = {}; // Update the name of this to be more descriptive, current area vector
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( sA, fieldSlotDeformationGradientCofactor[slotA], referenceAreaVectorA );

        real64 referenceAreaVectorB[3] = {};
        LvArray::tensorOps::scaledCopy< 3 >( referenceAreaVectorB,
                                             fieldSlotReferenceSurfaceNormal[slotB],
                                             fieldSlotReferenceArea[slotB] );

        real64 sB[3] = {}; // Update the name of this to be more descriptive, current area vector
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >(sB, fieldSlotDeformationGradientCofactor[slotB], referenceAreaVectorB );

        // Displacement for each field
        real64 dA[3] = {};
        LvArray::tensorOps::copy< 3 >( dA, fieldSlotDisplacement[slotA] );
        real64 dB[3] = {};
        LvArray::tensorOps::copy< 3 >( dB, fieldSlotDisplacement[slotB] );

        // Total mass for the contact pair.
        real64 mA = fieldSlotMass[slotA];
        real64 mB = fieldSlotMass[slotB];
        real64 mAB = mA + mB;

        // Outward normal of field A with respect to field B.
        real64 nAB[3] = {};

        // Mass-weighted average of the field normals
        LvArray::tensorOps::scaledCopy< 3 >( nAB, nA, mA );
        LvArray::tensorOps::scaledAdd< 3 >( nAB, nB, -mB );

        // Normalize the effective surface normal
        if( planeStrain == 1 )
        {
          nAB[2] = 0.0;
        }

        // If normal magnitude is zero for any reason just skip (e.g. no traction from cohesive law)
        real64 norm = LvArray::tensorOps::l2Norm< 3 >( nAB );
        if( norm < 1e-20 )
        {
          return;
        }

        // Normalize and flip (positive displacement is away from interface)
        LvArray::tensorOps::scale< 3 >( nAB, 1.0 / norm );

        real64 displacementVector[3] = {};
        LvArray::tensorOps::copy< 3 >( displacementVector, dA );
        LvArray::tensorOps::subtract< 3 >( displacementVector, dB );

        // The two cohesive fields may choose different periodic images for their absolute nodal displacements after
        // particles advect through a periodic boundary. The cohesive law depends on the opening vector, so reduce that
        // jump itself with the minimum-image convention before resolving normal and tangential components.
        if( periodic0 == 1 && domainExtent0 > 0.0 )
        {
          displacementVector[0] -= domainExtent0 * LvArray::math::floor( displacementVector[0] / domainExtent0 + 0.5 );
        }
        if( periodic1 == 1 && domainExtent1 > 0.0 )
        {
          displacementVector[1] -= domainExtent1 * LvArray::math::floor( displacementVector[1] / domainExtent1 + 0.5 );
        }
        if( periodic2 == 1 && domainExtent2 > 0.0 )
        {
          displacementVector[2] -= domainExtent2 * LvArray::math::floor( displacementVector[2] / domainExtent2 + 0.5 );
        }

        real64 totalNormalDisplacement = -LvArray::tensorOps::AiBi< 3 >( nAB, displacementVector ); 

        real64 tangentialInterfaceDisplacement[3]  = {};
        LvArray::tensorOps::copy< 3 >( tangentialInterfaceDisplacement, displacementVector );
        LvArray::tensorOps::scaledAdd< 3 >( tangentialInterfaceDisplacement, nAB, totalNormalDisplacement );
        real64 totalTangentialDisplacement = LvArray::tensorOps::l2Norm< 3 >( tangentialInterfaceDisplacement );

        // Call cohesive zone constitutive model update
        real64 normalStress = 0.0;
        real64 shearStress = 0.0;
        constitutiveWrapper.jumpDisplacementUpdate( k,
                                                    totalNormalDisplacement,
                                                    totalTangentialDisplacement,
                                                    normalStress,
                                                    shearStress );

        if( preventCZInterpentration == 1 && normalStress > 0 ) // Change from < to > since consititutive law returns negative stress for positive displacement
        {
          normalStress = 0.0;
        }

        LvArray::tensorOps::scaledCopy< 3 >( tA, nAB, -normalStress );
        LvArray::tensorOps::scaledCopy< 3 >( tB, nAB, normalStress );

        if( LvArray::math::abs( totalTangentialDisplacement ) > 1e-20 )
        {
          real64 tAB[3] = {}; // Tangent unit vector
          LvArray::tensorOps::copy< 3 >( tAB, tangentialInterfaceDisplacement );
          LvArray::tensorOps::scale< 3 >( tAB, 1 / totalTangentialDisplacement );

          LvArray::tensorOps::scaledAdd< 3 >( tA, tAB, shearStress ); // Flipped the sign of shear stress on this line and next
          LvArray::tensorOps::scaledAdd< 3 >( tB, tAB, -shearStress );
        }

        // Convert traction to force using mass-weighted average of projected area
        // Do we want to add choice of weighting as we do in contact calculations for normals?
        real64 areaAB[3] = {};
        LvArray::tensorOps::scaledCopy< 3 >( areaAB, sA, mA );
        LvArray::tensorOps::scaledAdd< 3 >( areaAB, sB, -mB );
        LvArray::tensorOps::scale< 3 >( areaAB, 1 / mAB );

        real64 surfaceArea = LvArray::tensorOps::AiBi< 3 >( nAB, areaAB ); // Should we take the absolute value to ensure surface area can never be negative. However, if nAB is consistent with areaAB then it should also never be negative so negative surface area could indicate an error
        for( localIndex i = 0; i < 3; ++i )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &fieldSlotCohesiveForce[slotA][i], tA[i] * surfaceArea );
          RAJA::atomicAdd( parallelDeviceAtomic{}, &fieldSlotCohesiveForce[slotB][i], tB[i] * surfaceArea );
        }
        
        // GEOS_LOG_RANK( "k: " << k << ", " << 
        //               //  "dA: " << "{" << dA[0] << ", " << dA[1] << ", " << dA[2] << "}, " << 
        //               //  "dB: " << "{" << dB[0] << ", " << dB[1] << ", " << dB[2] << "}, " << 
        //                "dTotal: " << "{" << displacementVector[0] << ", " << displacementVector[1] << ", " << displacementVector[2] << "}, " <<
        //               //  "nA: " << "{" << nA[0] << ", " << nA[1] << ", " << nA[2] << "}, " << 
        //               //  "nB: " << "{" << nB[0] << ", " << nB[1] << ", " << nB[2] << "}, " << 
        //                "nAB: " << "{" << nAB[0] << ", " << nAB[1] << ", " << nAB[2] << "}, " << 
        //               //  "aA: " << fieldSlotReferenceArea[slotA] << ", " <<
        //               //  "aB: " << fieldSlotReferenceArea[slotB] << ", " <<
        //               //  "sA: " << "{" << sA[0] << ", " << sA[1] << ", " << sA[2] << "}, " << 
        //               //  "sB: " << "{" << sB[0] << ", " << sB[1] << ", " << sB[2] << "}, " << 
        //                "normalDisp: " << totalNormalDisplacement << ", " << 
        //                "shearDisp: " << totalTangentialDisplacement << ", " <<
        //               //  "surfaceArea: " << surfaceArea << ", "
        //                "normalStress: " << normalStress << ", "
        //                "shearStress: " << shearStress << ", "
        //                "tA: " << "{" << tA[0] << ", " << tA[1] << ", " << tA[2] << "}, " << 
        //                "tB: " << "{" << tB[0] << ", " << tB[1] << ", " << tB[2] << "}"  );
      }


      // Save converged state
      // TODO: check state is being saved correctly
      // This subroutine may work correctly even if it is not so long as the model takes the total instantaneous displacements
      constitutiveWrapper.saveConvergedState( k, 0 );
    } );
  }
};

} // namespace solidMechanicsMPMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_EXPLICITMPM_HPP_ */
