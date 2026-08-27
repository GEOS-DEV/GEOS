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
  enum class PolarTestMode : int
  {
    Baseline,
    PreservedSnapshot,
    ReverseCallOrder,
    SingleCallForIdenticalInputs,
    NoInline
  };

  GEOS_HOST_DEVICE
  static bool matricesAreIdentical( real64 const (& matrixA)[3][3],
                                    real64 const (& matrixB)[3][3] )
  {
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        if( matrixA[i][j] != matrixB[i][j] )
        {
          return false;
        }
      }
    }
    return true;
  }

  // Change only this line between builds so each test keeps the same kernel
  // body except for the selected polar-decomposition strategy.
  static constexpr PolarTestMode polarTestMode = PolarTestMode::Baseline;

  GEOS_HOST_DEVICE
  __attribute__((noinline))
  static bool polarDecompositionNoInline( real64 (& rotation)[3][3],
                                          real64 const (& deformationGradient)[3][3] )
  {
    return LvArray::tensorOps::polarDecomposition< 3 >( rotation, deformationGradient );
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
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      real64 dt,
                      int const rank,
                      int hyperelasticUpdate,
                      arrayView3d< real64 const > const deformationGradient,
                      arrayView3d< real64 const > const fDot,
                      arrayView3d< real64 const > const velocityGradient,
                      arrayView2d< real64 > const particleStress,
                      arrayView1d< globalIndex const > const particleID  )
  {
    arrayView3d< real64, solid::STRESS_USD > const oldStress = constitutiveWrapper.m_oldStress;
    arrayView3d< real64, solid::STRESS_USD > const newStress = constitutiveWrapper.m_newStress;
    GEOS_UNUSED_VAR( rank );
    GEOS_UNUSED_VAR( particleID );

    // GEOS_LOG_RANK( "indices.size(): " << indices.size() 
    //                 << ", "
    //                 << "deformationGradient.size(0): "
    //                 << deformationGradient.size(0)
    //                 << ", "
    //                 << "fDot.size(0): "
    //                 << fDot.size(0)
    //                 << ", "
    //                 << "velocityGradient.size(0): "
    //                 << velocityGradient.size(0) 
    //                 << ", "
    //                 << "oldStress.size(0): " 
    //                 << oldStress.size(0)
    //                 << ", "
    //                 << "newStress.size(0): " 
    //                 << newStress.size(0)
    //                 << ", "
    //                 << "particleStress.size(0): "
    //                 << particleStress.size(0)
    //                 << ", "
    //                 << "particleID.size(0): "
    //                 << particleID.size(0) );

    if( indices.size() == 0 )
    {
      return;
    }

    // // Make arrays and array views to store and check the rotations after the kernel runs
    // // TODO

    // // Run smaller kernel to debug
    // forAll< POLICY >( count, [=] GEOS_HOST_DEVICE ( localIndex const q )
    //   {
    //     localIndex const k = begin + q;
    //     localIndex const p = indices[k];


    //     // deformationGradient[p] has already been advanced by updateDeformationGradient(), so it is F_{n+1}.
    //     // Recover the beginning-of-step F_n from F_{n+1} - Fdot * dt for the objective stress update.
    //     real64 fOld[3][3] = {};
    //     real64 fNew[3][3] = {};
    //     LvArray::tensorOps::copy< 3, 3 >( fNew, deformationGradient[p] );
    //     LvArray::tensorOps::copy< 3, 3 >( fOld, deformationGradient[p] );
    //     LvArray::tensorOps::scaledAdd< 3, 3 >( fOld, fDot[p], -dt );

    //     // Polar decompositions
    //     real64 rotBeginning[3][3] = {};
    //     real64 rotEnd[3][3] = {};

    
    //     bool const beginningConverged = LvArray::tensorOps::polarDecomposition< 3 >( rotBeginning, fOld );

    //     real64 rotBeginningSnapshot[3][3] = {};
    //     LvArray::tensorOps::copy< 3, 3 >( rotBeginningSnapshot, rotBeginning );

    //     bool const endConverged = LvArray::tensorOps::polarDecomposition< 3 >( rotEnd, fNew );    
    //   } );

    // Use the original full-rank launch that reproduces the bad rotations.
    localIndex const batchSize = indices.size(); //1536;
    for( localIndex begin = 0; begin < indices.size(); begin += batchSize )
    {
      localIndex const count = LvArray::math::min( batchSize, indices.size() - begin );

      forAll< POLICY >( count, [=] GEOS_HOST_DEVICE ( localIndex const q )
      {
        localIndex const k = begin + q;
        localIndex const p = indices[k];

        // Copy the beginning-of-step particle stress into the constitutive model's m_oldStress - this fixes the MPI sync issue on Lassen for
        // some reason
        #if defined(GEOS_USE_DEVICE)
        // Keep constitutive oldStress synchronized with the particle stress in
        // device builds. CUDA already needed this for MPI/MPM consistency; HIP
        // has the same host/device residency issue.
        LvArray::tensorOps::copy< 6 >( oldStress[p][0], particleStress[p] );
        #endif

        real64 stress[6] = {};
        //CC: debug hardcoded hyperelastic model for now
        if( hyperelasticUpdate == 1 )
        // if ( constitutiveWrapper.m_disableInelasticity ) // CC: Shouldn't there be a flag for hyperelastic models? otherwise we have to
        // manually add their name here everything we add them
        // Some models we might want hyperelastic updates when plasticity or damage are turned off
        { // Hyperelastic stress update
          // Don't believe we need to perform unrotation and rotation here (yes...unrotation...)
          // Think we can update stress directly by calling constitutive model
          // Hyperelastic models in GEOSX currently use FminusI as input argument
          real64 FminusI[3][3] = {};
          LvArray::tensorOps::copy< 3, 3 >( FminusI, deformationGradient[p] );
          LvArray::tensorOps::addIdentity< 3 >( FminusI, -1.0 );

          constitutiveWrapper.hyperUpdate( p,      // particle local index
                                          0,      // particles have 1 quadrature point
                                          FminusI, // particle strain increment
                                          stress );
        }
        else //Hypoeleastic stress update
        {
          // Determine the strain increment in Voigt notation
          real64 strainIncrement[6] = {};
          strainIncrement[0] = velocityGradient[p][0][0] * dt;
          strainIncrement[1] = velocityGradient[p][1][1] * dt;
          strainIncrement[2] = velocityGradient[p][2][2] * dt;
          strainIncrement[3] = (velocityGradient[p][1][2] + velocityGradient[p][2][1]) * dt;
          strainIncrement[4] = (velocityGradient[p][0][2] + velocityGradient[p][2][0]) * dt;
          strainIncrement[5] = (velocityGradient[p][0][1] + velocityGradient[p][1][0]) * dt;

          // deformationGradient[p] has already been advanced by updateDeformationGradient(), so it is F_{n+1}.
          // Recover the beginning-of-step F_n from F_{n+1} - Fdot * dt for the objective stress update.
          real64 fOld[3][3] = {};
          real64 fNew[3][3] = {};
          LvArray::tensorOps::copy< 3, 3 >( fNew, deformationGradient[p] );
          LvArray::tensorOps::copy< 3, 3 >( fOld, deformationGradient[p] );
          LvArray::tensorOps::scaledAdd< 3, 3 >( fOld, fDot[p], -dt );

          // Polar decompositions
          real64 rotBeginning[3][3] = {};
          real64 rotEnd[3][3] = {};

          if constexpr( polarTestMode == PolarTestMode::PreservedSnapshot )
          {
            bool const beginningConverged =
              LvArray::tensorOps::polarDecomposition< 3 >( rotBeginning, fOld );

            real64 rotBeginningSnapshot[3][3] = {};
            LvArray::tensorOps::copy< 3, 3 >( rotBeginningSnapshot, rotBeginning );

            bool const endConverged =
              LvArray::tensorOps::polarDecomposition< 3 >( rotEnd, fNew );

            GEOS_UNUSED_VAR( beginningConverged );
            GEOS_UNUSED_VAR( endConverged );

            // Passing the snapshot makes the first result genuinely live
            // across the second polar call without adding validation logic.
            constitutive::SolidUtilities::hypoUpdate2_StressOnly( constitutiveWrapper,
                                                                  p,
                                                                  0,
                                                                  dt,
                                                                  strainIncrement,
                                                                  rotBeginningSnapshot,
                                                                  rotEnd,
                                                                  stress );
          }
          else if constexpr( polarTestMode == PolarTestMode::Baseline )
          {
            bool const beginningConverged =
              LvArray::tensorOps::polarDecomposition< 3 >( rotBeginning, fOld );
            bool const endConverged =
              LvArray::tensorOps::polarDecomposition< 3 >( rotEnd, fNew );

            GEOS_UNUSED_VAR( beginningConverged );
            GEOS_UNUSED_VAR( endConverged );

            constitutive::SolidUtilities::hypoUpdate2_StressOnly( constitutiveWrapper,
                                                                  p,
                                                                  0,
                                                                  dt,
                                                                  strainIncrement,
                                                                  rotBeginning,
                                                                  rotEnd,
                                                                  stress );
          }
          else
          {
            bool beginningConverged = false;
            bool endConverged = false;

            if constexpr( polarTestMode == PolarTestMode::ReverseCallOrder )
            {
              endConverged =
                LvArray::tensorOps::polarDecomposition< 3 >( rotEnd, fNew );
              beginningConverged =
                LvArray::tensorOps::polarDecomposition< 3 >( rotBeginning, fOld );
            }
            else if constexpr( polarTestMode == PolarTestMode::SingleCallForIdenticalInputs )
            {
              beginningConverged =
                LvArray::tensorOps::polarDecomposition< 3 >( rotBeginning, fOld );

              if( matricesAreIdentical( fOld, fNew ) )
              {
                LvArray::tensorOps::copy< 3, 3 >( rotEnd, rotBeginning );
                endConverged = beginningConverged;
              }
              else
              {
                endConverged =
                  LvArray::tensorOps::polarDecomposition< 3 >( rotEnd, fNew );
              }
            }
            else if constexpr( polarTestMode == PolarTestMode::NoInline )
            {
              beginningConverged = polarDecompositionNoInline( rotBeginning, fOld );
              endConverged = polarDecompositionNoInline( rotEnd, fNew );
            }
            GEOS_UNUSED_VAR( beginningConverged );
            GEOS_UNUSED_VAR( endConverged );

            constitutive::SolidUtilities::hypoUpdate2_StressOnly( constitutiveWrapper,
                                                                  p,
                                                                  0,
                                                                  dt,
                                                                  strainIncrement,
                                                                  rotBeginning,
                                                                  rotEnd,
                                                                  stress );
          }
        }

        // Copy the updated stress into particleStress
        LvArray::tensorOps::copy< 6 >( particleStress[p], stress );

        // Copy m_newStress into m_oldStress
        constitutiveWrapper.saveConvergedState( p, 0 );
    
      } );

      parallelDeviceSync(); // diagnostic only
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
   * @param[in] fieldA Velocity field A of cohesive zone
   * @param[in] fieldB Velocity field B of cohesive zone
   * @param[in] periodic0 Periodic flag for x-direction
   * @param[in] periodic1 Periodic flag for y-direction
   * @param[in] periodic2 Periodic flag for z-direction
   * @param[in] domainExtent0 Periodic domain extent in x-direction
   * @param[in] domainExtent1 Periodic domain extent in y-direction
   * @param[in] domainExtent2 Periodic domain extent in z-direction
   * @param[in] gridMass Array view of the nodal mass
   * @param[in] gridDisplacement Array view of nodal displacement
   * @param[in] gridParticleSurfaceNormal Array view of particle mapped surface normals
   * @param[in] gridDeformationGradientCofactor Array view of mapped deformation gradient cofactor
   * @param[in] gridCZTraction Array view of cohesive tractions for each velocity field
   * @param[in] czReferenceSurfaceNormal Array view of the reference surface normals for each velocity field
   * @param[in] czReferenceArea Scalar reference area of each cohesive zone node and velocity field
   
   */
  template< typename POLICY, typename CONSTITUTIVE_WRAPPER >
  static void launch( int numNodes,
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      real64 dt,
                      int planeStrain, // Should remove eventually, used for normals but normals after mapping to grid should be checked for planeStrain condition
                      real64 smallMass,
                      int preventCZInterpentration,
                      localIndex fieldA,
                      localIndex fieldB,
                      int periodic0,
                      int periodic1,
                      int periodic2,
                      real64 domainExtent0,
                      real64 domainExtent1,
                      real64 domainExtent2,
                      arrayView2d< real64 const > const gridMass,
                      arrayView3d< real64 const > const gridDisplacement,
                      arrayView3d< real64 const > const gridParticleSurfaceNormal,
                      arrayView4d< real64 const > const gridDeformationGradientCofactor,
                      arrayView3d< real64 > const gridCZTraction,
                      arrayView3d< real64 const > const czReferenceSurfaceNormal,
                      arrayView2d< real64 const > const czReferenceArea )
  {
    GEOS_UNUSED_VAR( dt );

    // Perform constitutive call
    forAll< POLICY >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {      
      bool active = ( gridMass[k][fieldA] > smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( gridParticleSurfaceNormal[k][fieldA] ) > 1.0e-16 )
                    and
                    ( gridMass[k][fieldB] > smallMass ) && ( LvArray::tensorOps::l2NormSquared< 3 >( gridParticleSurfaceNormal[k][fieldB] ) > 1.0e-16 );

      if( active )
      {
        // Copy normals
        real64 nA[3] = {};
        real64 nB[3] = {};
        LvArray::tensorOps::copy< 3 >( nA, gridParticleSurfaceNormal[k][fieldA] );
        LvArray::tensorOps::copy< 3 >( nB, gridParticleSurfaceNormal[k][fieldB] );

        // Initialize tractions here
        real64 tA[3] = {};
        real64 tB[3] = {};
        LvArray::tensorOps::fill< 3 >( tA, 0.0 );
        LvArray::tensorOps::fill< 3 >( tB, 0.0 );

        // Compute updated nodal area vectors for cohesive zone traction calculations
        real64 referenceAreaVectorA[3] = {};
        LvArray::tensorOps::scaledCopy< 3 >( referenceAreaVectorA, czReferenceSurfaceNormal[k][fieldA], czReferenceArea[k][fieldA] );

        real64 sA[3] = {}; // Update the name of this to be more descriptive, current area vector
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( sA, gridDeformationGradientCofactor[k][fieldA], referenceAreaVectorA );

        real64 referenceAreaVectorB[3] = {};
        LvArray::tensorOps::scaledCopy< 3 >( referenceAreaVectorB, czReferenceSurfaceNormal[k][fieldB], czReferenceArea[k][fieldB] );

        real64 sB[3] = {}; // Update the name of this to be more descriptive, current area vector
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >(sB, gridDeformationGradientCofactor[k][fieldB], referenceAreaVectorB );

        // Displacement for each field
        real64 dA[3] = {};
        LvArray::tensorOps::copy< 3 >( dA, gridDisplacement[k][fieldA] );
        real64 dB[3] = {};
        LvArray::tensorOps::copy< 3 >( dB, gridDisplacement[k][fieldB] );

        // Total mass for the contact pair.
        real64 mA = gridMass[k][fieldA];
        real64 mB = gridMass[k][fieldB];
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
        LvArray::tensorOps::scaledAdd< 3 >( gridCZTraction[k][fieldA], tA, surfaceArea );
        LvArray::tensorOps::scaledAdd< 3 >( gridCZTraction[k][fieldB], tB, surfaceArea );
        
        // GEOS_LOG_RANK( "k: " << k << ", " << 
        //               //  "dA: " << "{" << dA[0] << ", " << dA[1] << ", " << dA[2] << "}, " << 
        //               //  "dB: " << "{" << dB[0] << ", " << dB[1] << ", " << dB[2] << "}, " << 
        //                "dTotal: " << "{" << displacementVector[0] << ", " << displacementVector[1] << ", " << displacementVector[2] << "}, " <<
        //               //  "nA: " << "{" << nA[0] << ", " << nA[1] << ", " << nA[2] << "}, " << 
        //               //  "nB: " << "{" << nB[0] << ", " << nB[1] << ", " << nB[2] << "}, " << 
        //                "nAB: " << "{" << nAB[0] << ", " << nAB[1] << ", " << nAB[2] << "}, " << 
        //               //  "aA: " << czReferenceArea[k][fieldA] << ", " << 
        //               //  "aB: " << czReferenceArea[k][fieldB] << ", " << 
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
