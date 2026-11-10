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
 * @brief Particle stress updates and nodal cohesive-zone constitutive kernels for explicit MPM.
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
 * @brief Dispatches particle constitutive updates over a selected set of particles.
 *
 * The two launch overloads share a batched, execution-policy-based traversal
 * but use different constitutive kinematics: a hyperelastic evaluation from
 * the end-of-step deformation gradient, or an incremental hypoelastic update
 * from the velocity gradient and the beginning/end rotations. Each particle
 * has one constitutive quadrature point, indexed by zero.
 *
 * The supplied indices are particle storage indices, not consecutive active
 * ordinals. Results are written back to the same particle index and the
 * constitutive wrapper's converged state is committed after each update.
 */
struct ParticleStateUpdateKernel
{
  /**
   * @brief Update selected particle stresses through the hyperelastic wrapper.
   *
   * Passes the total end-of-step deformation relative to the identity,
   * FminusI = F_{n+1} - I, to hyperUpdate. This is not the incremental strain
   * computed by the hypoelastic overload. The constitutive wrapper supplies
   * the material response and internal-history update; the returned six
   * stress components are copied to the particle-level output field.
   *
   * The host divides the selected index list into batches; the execution
   * policy processes particles within each batch. A nonempty list requires
   * batchSize > 0, and the caller must provide valid particle indices.
   *
   * @tparam POLICY Host or device execution policy for particle-parallel work.
   * @tparam CONSTITUTIVE_WRAPPER Constitutive update wrapper supporting hyperUpdate.
   * @param indices Selected particle storage indices.
   * @param batchSize Maximum number of selected particles in one launch.
   * @param constitutiveWrapper Material wrapper and constitutive history storage.
   * @param deformationGradient End-of-step full 3x3 deformation gradient for each particle.
   * @param particleStress Destination particle stress in six-component Voigt form.
   */
  template< typename POLICY, typename CONSTITUTIVE_WRAPPER >
  static void launch( SortedArrayView< localIndex const > const indices,
                      localIndex const batchSize,
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      arrayView3d< real64 const > const deformationGradient,
                      arrayView2d< real64 > const particleStress )
  {
    // The wrapper owns both stress-history arrays. This overload does not
    // access the bound views directly; hyperUpdate and saveConvergedState
    // perform the material update and history commit.
    arrayView3d< real64, solid::STRESS_USD > const oldStress = constitutiveWrapper.m_oldStress;
    arrayView3d< real64, solid::STRESS_USD > const newStress = constitutiveWrapper.m_newStress;

    // An empty selection does not launch any device work. Other selections
    // are assumed to have a strictly positive batchSize.
    if( indices.size() == 0 )
    {
      return;
    }

    // Limit the execution-policy launch size without reindexing particle data.
    for( localIndex begin = 0; begin < indices.size(); begin += batchSize )
    {
      localIndex const count = LvArray::math::min( batchSize, indices.size() - begin );

      forAll< POLICY >( count, [=] GEOS_HOST_DEVICE ( localIndex const q )
      {
        // k addresses the selected-list entry; p addresses particle storage.
        localIndex const k = begin + q;
        localIndex const p = indices[k];

        real64 stress[6] = {};

        // This wrapper consumes F - I (a total deformation measure), not
        // an increment relative to the preceding step. The full 3x3 tensor
        // is supplied even when particle geometry evolves in two dimensions.
        real64 FminusI[3][3] = {};
        LvArray::tensorOps::copy< 3, 3 >( FminusI, deformationGradient[p] );
        LvArray::tensorOps::addIdentity< 3 >( FminusI, -1.0 );

        // One material quadrature point per particle. The wrapper supplies
        // the material-specific stress response and internal state evolution.
        constitutiveWrapper.hyperUpdate( p,       // Particle storage index.
                                         0,       // Single quadrature point.
                                         FminusI, // Total F - I, not a step increment.
                                         stress );

        // Publish the returned Voigt stress for MPM force assembly/output.
        // Then advance the wrapper's history before the next solver step.
        LvArray::tensorOps::copy< 6 >( particleStress[p], stress );
        constitutiveWrapper.saveConvergedState( p, 0 );
      } );
    }
  }

  /**
   * @brief Update selected particle stresses using incremental, rotated kinematics.
   *
   * Constructs the step strain increment from dt times the symmetric velocity
   * gradient. Normal components use L_ii * dt; shear components use the
   * engineering-strain convention (L_ij + L_ji) * dt, ordered
   * [xx, yy, zz, yz, xz, xy]. The constitutive helper also receives the
   * rotations at the beginning and end of the step for objective stress
   * updating. It returns six stress components, then the wrapper commits
   * the updated constitutive state at quadrature point zero.
   *
   * deformationGradient is retained in this overload's interface but is
   * explicitly unused by its current body; the stress increment is driven
   * by velocityGradient and the rotation pair instead.
   *
   * @tparam POLICY Host or device execution policy for particle-parallel work.
   * @tparam CONSTITUTIVE_WRAPPER Material wrapper accepted by hypoUpdate2_StressOnly.
   * @param indices Selected particle storage indices.
   * @param batchSize Maximum number of selected particles in one launch; must be positive for nonempty indices.
   * @param constitutiveWrapper Material wrapper and constitutive history storage.
   * @param dt Physical time increment used in the constitutive update.
   * @param deformationGradient Present for interface compatibility; not used by this overload.
   * @param rotation End-of-step physical rotation tensor for each particle.
   * @param oldRotation Beginning-of-step physical rotation tensor for each particle.
   * @param velocityGradient Step velocity-gradient tensor L for each particle.
   * @param particleStress Destination particle stress in six-component Voigt form.
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
    // The incremental hypoelastic path does not consume F directly. It uses
    // the time-integrated velocity gradient and the rotation increment below.
    GEOS_UNUSED_VAR( deformationGradient );

    // Constitutive stress history is managed by the wrapper; these views are
    // bound here but are not accessed directly in the particle loop.
    arrayView3d< real64, solid::STRESS_USD > const oldStress = constitutiveWrapper.m_oldStress;
    arrayView3d< real64, solid::STRESS_USD > const newStress = constitutiveWrapper.m_newStress;

    // Match the hyperelastic path's selected-particle and batching contract.
    if( indices.size() == 0 )
    {
      return;
    }

    for( localIndex begin = 0; begin < indices.size(); begin += batchSize )
    {
      localIndex const count = LvArray::math::min( batchSize, indices.size() - begin );

      forAll< POLICY >( count, [=] GEOS_HOST_DEVICE ( localIndex const q )
      {
        // Convert the batch-relative entry to its actual particle index.
        localIndex const k = begin + q;
        localIndex const p = indices[k];

        real64 stress[6] = {};

        // Form the symmetric small-strain increment from the step-averaged
        // velocity gradient L = grad(v). The lower three Voigt components
        // are engineering shear strains, not tensorial shear components:
        // [exx, eyy, ezz, gamma_yz, gamma_xz, gamma_xy].
        real64 strainIncrement[6] = {};
        strainIncrement[0] = velocityGradient[p][0][0] * dt;
        strainIncrement[1] = velocityGradient[p][1][1] * dt;
        strainIncrement[2] = velocityGradient[p][2][2] * dt;
        strainIncrement[3] = (velocityGradient[p][1][2] + velocityGradient[p][2][1]) * dt;
        strainIncrement[4] = (velocityGradient[p][0][2] + velocityGradient[p][2][0]) * dt;
        strainIncrement[5] = (velocityGradient[p][0][1] + velocityGradient[p][1][0]) * dt;

        // The constitutive helper handles stress rotation using the stored
        // beginning/end physical rotations; do not rotate the strain or
        // result a second time in this particle kernel.
        real64 rotBeginning[3][3] = {};
        real64 rotEnd[3][3] = {};
        LvArray::tensorOps::copy< 3, 3 >( rotBeginning, oldRotation[p] );
        LvArray::tensorOps::copy< 3, 3 >( rotEnd, rotation[p] );

        // The objective, stress-only update is material-model dependent.
        // Passing dt as well as the strain increment allows the wrapper to
        // evolve any rate-dependent internal state consistently.
        constitutive::SolidUtilities::hypoUpdate2_StressOnly( constitutiveWrapper,
                                                              p,
                                                              0,
                                                              dt,
                                                              strainIncrement,
                                                              rotBeginning,
                                                              rotEnd,
                                                              stress );

        // Make the updated stress available to MPM and commit the material
        // history after the per-particle constitutive evaluation.
        LvArray::tensorOps::copy< 6 >( particleStress[p], stress );
        constitutiveWrapper.saveConvergedState( p, 0 );
      } );
    }
  }
};

/**
 * @brief Evaluates cohesive traction on paired nodal velocity fields.
 *
 * For each grid node with appreciable mass and surface normals in both
 * selected fields, this kernel constructs a contact-pair normal, computes
 * its periodicity-corrected displacement jump, evaluates a cohesive law,
 * and converts the returned normal/tangential tractions to equal-and-opposite
 * nodal force contributions using a deformed, signed interface area.
 *
 * The cohesive constitutive law is evaluated from the total displacement
 * jump at a node, rather than a velocity or an increment of displacement.
 * Each node is passed to the material wrapper as one constitutive point.
 */
struct CohesiveZoneStateUpdateKernel
{
  /**
   * @brief Scatter cohesive-law force contributions for a pair of nodal fields.
   *
   * Both fields must exceed smallMass and have a mapped normal with squared
   * norm greater than 1e-16. The effective normal is the normalized
   * difference mA*nA - mB*nB (with its z component suppressed when
   * planeStrain == 1). The jump dA - dB is reduced by the minimum-image
   * convention in enabled periodic directions before separation into normal
   * opening and tangential slip. In the present sign convention,
   * opening = -nAB . (dA - dB).
   *
   * The constitutive wrapper returns normalStress and shearStress for that
   * total jump. The kernel resolves these as opposing field tractions,
   * projects the mass-weighted difference of current field area vectors onto
   * nAB, and adds traction times signed projected area to gridCZTraction.
   * dt is currently unused: it does not scale the displacement jump or
   * assembled force in this kernel.
   *
   * @tparam POLICY Host/device execution policy for the nodal loop.
   * @tparam CONSTITUTIVE_WRAPPER Wrapper providing jumpDisplacementUpdate and saveConvergedState.
   * @param numNodes Number of grid nodes to visit.
   * @param constitutiveWrapper Cohesive constitutive law and history storage.
   * @param dt Time step retained in the call interface; unused here.
   * @param planeStrain When equal to 1, project the pair normal into the x-y plane.
   * @param smallMass Threshold below which a field is inactive at this node.
   * @param preventCZInterpentration Suppress positive (compressive) cohesive normal traction when equal to 1.
   * @param fieldA First nodal velocity/contact field.
   * @param fieldB Opposing nodal velocity/contact field.
   * @param periodic0 Periodic-image adjustment flag for x.
   * @param periodic1 Periodic-image adjustment flag for y.
   * @param periodic2 Periodic-image adjustment flag for z.
   * @param domainExtent0 Periodic length in x, used only when positive and periodic0 is 1.
   * @param domainExtent1 Periodic length in y, used only when positive and periodic1 is 1.
   * @param domainExtent2 Periodic length in z, used only when positive and periodic2 is 1.
   * @param gridMass Nodal mass per field.
   * @param gridDisplacement Total nodal displacement per field.
   * @param gridParticleSurfaceNormal Mapped outward surface-normal vectors per field.
   * @param gridDeformationGradientCofactor Cofactor of mapped deformation gradient per node and field.
   * @param gridCZTraction In/out nodal cohesive-force accumulator per field.
   * @param czReferenceSurfaceNormal Reference surface normal associated with each field.
   * @param czReferenceArea Reference scalar surface area per field.
   */
  template< typename POLICY, typename CONSTITUTIVE_WRAPPER >
  static void launch( int numNodes,
                      CONSTITUTIVE_WRAPPER const & constitutiveWrapper,
                      real64 dt,
                      int planeStrain, // TODO: Assess whether the extra 2D normal projection remains necessary.
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
    // This cohesive law is driven by accumulated displacement jumps rather
    // than the step duration. The parameter is retained for the caller API.
    GEOS_UNUSED_VAR( dt );

    // Process every node independently. Inactive pairs still reach the
    // saveConvergedState call at the end of the node lambda.
    forAll< POLICY >( numNodes, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // A cohesive pair needs both positive-mass fields and meaningful
      // independently mapped normals; a single populated field is insufficient.
      bool active = ( gridMass[k][fieldA] > smallMass ) &&
                    ( LvArray::tensorOps::l2NormSquared< 3 >( gridParticleSurfaceNormal[k][fieldA] ) > 1.0e-16 )
                    and
                    ( gridMass[k][fieldB] > smallMass ) &&
                    ( LvArray::tensorOps::l2NormSquared< 3 >( gridParticleSurfaceNormal[k][fieldB] ) > 1.0e-16 );

      if( active )
      {
        // Keep the per-field normals separate: their mass-weighted
        // difference defines the pair direction, not either normal alone.
        real64 nA[3] = {};
        real64 nB[3] = {};
        LvArray::tensorOps::copy< 3 >( nA, gridParticleSurfaceNormal[k][fieldA] );
        LvArray::tensorOps::copy< 3 >( nB, gridParticleSurfaceNormal[k][fieldB] );

        // Nodal force directions will be assembled after the cohesive law
        // returns its scalar normal and tangential traction components.
        real64 tA[3] = {};
        real64 tB[3] = {};
        LvArray::tensorOps::fill< 3 >( tA, 0.0 );
        LvArray::tensorOps::fill< 3 >( tB, 0.0 );

        // Nanson's area-vector transform: s = cof(F) * (N0 * A0).
        // Each field has its own current area vector, which may differ due
        // to its mapped deformation and reference interface geometry.
        real64 referenceAreaVectorA[3] = {};
        LvArray::tensorOps::scaledCopy< 3 >( referenceAreaVectorA, czReferenceSurfaceNormal[k][fieldA], czReferenceArea[k][fieldA] );

        real64 sA[3] = {}; // Current oriented area vector for field A.
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( sA, gridDeformationGradientCofactor[k][fieldA], referenceAreaVectorA );

        real64 referenceAreaVectorB[3] = {};
        LvArray::tensorOps::scaledCopy< 3 >( referenceAreaVectorB, czReferenceSurfaceNormal[k][fieldB], czReferenceArea[k][fieldB] );

        real64 sB[3] = {}; // Current oriented area vector for field B.
        LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( sB, gridDeformationGradientCofactor[k][fieldB], referenceAreaVectorB );

        // Read the two total nodal displacements before resolving the
        // interface jump; the law does not consume a step increment here.
        real64 dA[3] = {};
        LvArray::tensorOps::copy< 3 >( dA, gridDisplacement[k][fieldA] );
        real64 dB[3] = {};
        LvArray::tensorOps::copy< 3 >( dB, gridDisplacement[k][fieldB] );

        // Use the two nodal masses consistently for the pair normal and
        // area averaging. Activity checks above make mAB positive.
        real64 mA = gridMass[k][fieldA];
        real64 mB = gridMass[k][fieldB];
        real64 mAB = mA + mB;

        // Effective pair normal: nAB ~ mA*nA - mB*nB. This convention
        // defines which sign of projected relative displacement is opening.
        real64 nAB[3] = {};

        // Subtract the opposing field's normal to combine two outward
        // orientations into one pairwise direction.
        LvArray::tensorOps::scaledCopy< 3 >( nAB, nA, mA );
        LvArray::tensorOps::scaledAdd< 3 >( nAB, nB, -mB );

        // Native two-dimensional formulations use only the meridional /
        // in-plane components for cohesive normal separation.
        if( planeStrain == 1 )
        {
          nAB[2] = 0.0;
        }

        // Parallel/aligned mapped normals can cancel despite both fields
        // passing the individual-norm test. A degenerate pair skips the
        // traction update AND the final saveConvergedState below, because
        // this return exits the current node's device lambda.
        real64 norm = LvArray::tensorOps::l2Norm< 3 >( nAB );
        if( norm < 1e-20 )
        {
          return;
        }

        // Normalize the mass-weighted pair direction; opening is given the
        // opposite sign in the projection of the relative displacement.
        LvArray::tensorOps::scale< 3 >( nAB, 1.0 / norm );

        real64 displacementVector[3] = {};
        LvArray::tensorOps::copy< 3 >( displacementVector, dA );
        LvArray::tensorOps::subtract< 3 >( displacementVector, dB );

        // Each field can track a different periodic image of the same
        // interface. Map each Cartesian component of dA-dB into the nearest
        // periodic image: delta -= L*floor(delta/L + 1/2). This operates on
        // the jump, not on either field's stored displacement separately.
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

        // opening = -nAB . (dA-dB); the minus sign is part of the
        // constitutive opening convention. The tangential vector below is
        // the residual after removing that normal component.
        real64 totalNormalDisplacement = -LvArray::tensorOps::AiBi< 3 >( nAB, displacementVector );

        real64 tangentialInterfaceDisplacement[3] = {};
        LvArray::tensorOps::copy< 3 >( tangentialInterfaceDisplacement, displacementVector );
        LvArray::tensorOps::scaledAdd< 3 >( tangentialInterfaceDisplacement, nAB, totalNormalDisplacement );
        real64 totalTangentialDisplacement = LvArray::tensorOps::l2Norm< 3 >( tangentialInterfaceDisplacement );

        // The model evaluates normal and shear traction magnitudes from
        // total separation, allowing any history/damage law in the wrapper.
        // Normal and tangential responses are scalar here; the projection
        // back to a physical force vector is performed immediately below.
        real64 normalStress = 0.0;
        real64 shearStress = 0.0;
        constitutiveWrapper.jumpDisplacementUpdate( k,
                                                    totalNormalDisplacement,
                                                    totalTangentialDisplacement,
                                                    normalStress,
                                                    shearStress );

        // The law uses negative normalStress for tensile opening. When the
        // option is enabled, zero the positive/compressive cohesive response;
        // this does not by itself enforce a nonpenetration constraint.
        if( preventCZInterpentration == 1 && normalStress > 0 )
        {
          normalStress = 0.0;
        }

        // Equal-and-opposite normal traction vectors for the two fields.
        LvArray::tensorOps::scaledCopy< 3 >( tA, nAB, -normalStress );
        LvArray::tensorOps::scaledCopy< 3 >( tB, nAB, normalStress );

        // Resolve shear along the actual tangential displacement direction.
        // At negligible slip, avoid dividing by its magnitude and do not
        // impose an arbitrary tangential basis.
        if( LvArray::math::abs( totalTangentialDisplacement ) > 1e-20 )
        {
          real64 tAB[3] = {}; // Tangent unit vector
          LvArray::tensorOps::copy< 3 >( tAB, tangentialInterfaceDisplacement );
          LvArray::tensorOps::scale< 3 >( tAB, 1 / totalTangentialDisplacement );

          LvArray::tensorOps::scaledAdd< 3 >( tA, tAB, shearStress );
          LvArray::tensorOps::scaledAdd< 3 >( tB, tAB, -shearStress );
        }

        // Convert traction to nodal force using the signed area projection:
        // A_AB = nAB . (mA*sA - mB*sB)/(mA + mB).
        // The opposing area-vector orientation is subtracted, matching the
        // pair-normal construction. No positive-area clamp is applied.
        // TODO: Assess whether selectable area weighting, as for contact
        // normals, would be preferable for strongly asymmetric field masses.
        real64 areaAB[3] = {};
        LvArray::tensorOps::scaledCopy< 3 >( areaAB, sA, mA );
        LvArray::tensorOps::scaledAdd< 3 >( areaAB, sB, -mB );
        LvArray::tensorOps::scale< 3 >( areaAB, 1 / mAB );

        // TODO: Verify the intended sign convention when this projected
        // area becomes negative; the current implementation keeps the sign.
        real64 surfaceArea = LvArray::tensorOps::AiBi< 3 >( nAB, areaAB );

        // Accumulate, rather than overwrite, cohesive nodal force. Both
        // fields receive the same signed area with opposite traction vectors.
        LvArray::tensorOps::scaledAdd< 3 >( gridCZTraction[k][fieldA], tA, surfaceArea );
        LvArray::tensorOps::scaledAdd< 3 >( gridCZTraction[k][fieldB], tB, surfaceArea );
      }

      // The current code commits cohesive history even for inactive nodes;
      // only the degenerate-pair-normal return inside the branch skips it.
      // TODO: Confirm saveConvergedState is appropriate for inactive nodes
      // and consistent with the total-jump constitutive history convention.
      constitutiveWrapper.saveConvergedState( k, 0 );
    } );
  }
};

} // namespace solidMechanicsMPMKernels

} // namespace geos

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_EXPLICITMPM_HPP_ */
