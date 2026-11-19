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
 * @file SolidMechanicsALMKernelsBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMKERNELSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMKERNELSBASE_HPP_

#include "finiteElement/kernelInterface/InterfaceKernelBase.hpp"
#include "SolidMechanicsConformingContactKernelsHelper.hpp"

namespace geos
{

namespace solidMechanicsALMKernels
{

/**
 * @brief A struct to check for constraint satisfaction
 */
struct ConstraintCheckKernel
{

  /**
   * @brief Launch the kernel function to check the constraint satisfaction
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONTACT_WRAPPER the type of contact wrapper doing the fracture traction updates
   * @param[in] size the size of the subregion
   * @param[in] traction the array containing the current traction
   * @param[in] trialTraction the array containing the traction updated from the current displacement jump
   * @param[in] dispJump the array containing the displacement jump
   * @param[in] deltaDispJump the array containing the delta displacement jump
   * @param[in] normalTractionTolerance Check tolerance (normal traction)
   * @param[in] normalDisplacementTolerance Check tolerance (compenetration)
   * @param[in] slidingTolerance Check tolerance (sliding)
   * @param[in] slidingCheckTolerance Check tolerance (if shear strass exceeds tauLim)
   * @param[in] area interface element area
   * @param[in] fractureState the array containing the fracture state
   * @param[out] condConv the array containing the convergence flag:
   *                      0: Constraint conditions satisfied
   *                      1: Open
   *                      2: Compenetration
   *                      3: Slip exceeds sliding tolerance
   *                      4: Shear stress exceeds tauLim
   */
  template< typename POLICY, typename CONTACT_WRAPPER >
  static void
  launch( localIndex const size,
          CONTACT_WRAPPER const & contactWrapper,
          arrayView1d< integer const > const & ghostRank,
          arrayView2d< real64 const > const & traction,
          arrayView2d< real64 const > const & trialTraction,
          arrayView2d< real64 const > const & dispJump,
          arrayView2d< real64 const > const & deltaDispJump,
          arrayView1d< real64 const > const & normalTractionTolerance,
          arrayView1d< real64 const > const & normalDisplacementTolerance,
          arrayView1d< real64 const > const & slidingTolerance,
          real64 const slidingCheckTolerance,
          arrayView1d< integer const > const fractureState,
          arrayView1d< integer > const & condConv )
  {

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      if( ghostRank[k] < 0 )
      {
        contactWrapper.constraintCheck( k,
                                        dispJump[k],
                                        deltaDispJump[k],
                                        traction[k],
                                        trialTraction[k][0],
                                        fractureState[k],
                                        normalTractionTolerance[k],
                                        normalDisplacementTolerance[k],
                                        slidingTolerance[k],
                                        slidingCheckTolerance,
                                        condConv[k] );
      }

    } );
  }
};

/**
 * @brief A struct to adapt the per-element relaxation factor of the multiplier update
 */
struct AdaptRelaxationKernel
{
  /**
   * @brief Lower the relaxation factor of the elements whose normal constraint flag reverses
   * @tparam POLICY the type of policy used in the kernel launch
   * @param[in] size the size of the subregion
   * @param[in] ghostRank the ghost rank of the elements (ghosts are skipped: their factor is synchronized)
   * @param[in] trialTraction the array containing the traction updated from the current displacement jump
   * @param[in] normalTractionTolerance Check tolerance (normal traction)
   * @param[in] condConv the convergence flag of the current configuration check
   * @param[in] minRelaxation the lowest relaxation factor omega_min
   * @param[inout] previousNormalFlag sign of the last nonzero normal flag (+: gN > 0 on a closed element,
   * -: compenetration) times the number of consecutive checks with that flag
   * @param[inout] relaxation the per-element relaxation factor
   * @return the number of (locally owned) elements with a relaxation factor lower than 1
   *
   * @detail An element starts each time step with omega = 1. Each time its normal flag reverses (small opening on a
   * closed element <-> compenetration), the augmentation has overshot: its factor is halved, down to omega_min.
   * Once the flag has stayed on the same side for two more checks in a row, the update undershoots: the factor is
   * doubled at each further same-side check, up to 1. An element that keeps flipping stays damped.
   * Elements that approach the constraint from one side keep the full augmentation.
   */
  template< typename POLICY >
  static localIndex
  launch( localIndex const size,
          arrayView1d< integer const > const & ghostRank,
          arrayView2d< real64 const > const & trialTraction,
          arrayView1d< real64 const > const & normalTractionTolerance,
          arrayView1d< int const > const & condConv,
          real64 const minRelaxation,
          arrayView1d< integer > const & previousNormalFlag,
          arrayView1d< real64 > const & relaxation )
  {
    RAJA::ReduceSum< ReducePolicy< POLICY >, localIndex > numRelaxed( 0 );

    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      if( ghostRank[k] >= 0 )
      {
        return;
      }
      bool const isTensile = trialTraction[k][0] >= normalTractionTolerance[k];
      integer const normalFlag = ( condConv[k] == 1 && !isTensile ) ? 1 : ( condConv[k] == 2 ? -1 : 0 );
      if( normalFlag != 0 )
      {
        // previousNormalFlag = sign( last flag ) * number of consecutive checks with that flag
        integer const previous = previousNormalFlag[k];
        integer streak = 1;
        if( normalFlag * previous < 0 )
        {
          // reversal: the augmentation has overshot
          relaxation[k] = LvArray::math::max( minRelaxation, 0.5 * relaxation[k] );
        }
        else if( normalFlag * previous > 0 )
        {
          streak = LvArray::math::min( LvArray::math::abs( previous ) + 1, 1000 );
          if( streak >= 3 )
          {
            // two more checks in a row on the same side: the update undershoots, let the augmentation recover
            relaxation[k] = LvArray::math::min( 1.0, 2.0 * relaxation[k] );
          }
        }
        previousNormalFlag[k] = normalFlag * streak;
      }
      if( relaxation[k] < 1.0 )
      {
        numRelaxed += 1;
      }
    } );

    return numRelaxed.get();
  }
};

/**
 * @brief A struct to check for constraint satisfaction
 */
struct UpdateStateKernel
{

  /**
   * @brief Launch the kernel function to check the constraint satisfaction
   * @tparam POLICY the type of policy used in the kernel launch
   * @tparam CONTACT_WRAPPER the type of contact wrapper doing the fracture traction updates
   * @param[in] size the size of the subregion
   * @param[in] oldDispJump the array containing the old displacement jump (previous time step)
   * @param[in] dispJump the array containing the displacement jump
   * @param[in] penalty the array containing the penalty coefficients
   * @param[in] symmetric flag to compute symmetric penalty matrix
   * @param[in] normalTractionTolerance Check tolerance (normal traction)
   * @param[in] traction the array containing the current traction
   * @param[in] relaxation per-element relaxation factor omega in (0, 1] of the multiplier update (1: no relaxation)
   * @param[in] fractureState the array containing the fracture state
   *
   * @detail With relaxation, the multiplier is moved only part of the way towards its augmented (and projected)
   * value: lambda <- lambda + omega * ( proj( lambda + k g ) - lambda ). In stick this is lambda <- lambda + omega k g.
   * The new multiplier is a convex combination of two admissible tractions, so it remains within the Coulomb cone.
   * The fracture state is the one of the projected traction.
   */
  template< typename POLICY, typename CONTACT_WRAPPER >
  static void
  launch( localIndex const size,
          CONTACT_WRAPPER const & contactWrapper,
          arrayView2d< real64 const > const & oldDispJump,
          arrayView2d< real64 const > const & dispJump,
          arrayView2d< real64 > const & penalty,
          bool const symmetric,
          arrayView1d< real64 const > const & normalTractionTolerance,
          arrayView2d< real64 > const & traction,
          arrayView1d< real64 const > const & relaxation,
          arrayView1d< integer > const & fractureState )

  {
    forAll< POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {

      real64 const zero = LvArray::NumericLimits< real64 >::epsilon;

      real64 localPenalty[3][3]{};
      real64 localTractionNew[3]{};
      contactWrapper.updateTraction( k,
                                     oldDispJump[k],
                                     dispJump[k],
                                     penalty[k],
                                     traction[k],
                                     symmetric,
                                     false,
                                     normalTractionTolerance[k],
                                     zero,
                                     localPenalty,
                                     localTractionNew,
                                     fractureState[k] );

      if( fractureState[k] == fields::contact::FractureState::Open )
      {
        LvArray::tensorOps::fill< 3 >( localTractionNew, 0.0 );
      }
      else if( LvArray::math::abs( localTractionNew[ 0 ] ) < normalTractionTolerance[k] )
      {
        // When normal traction is lower than normalTractionTolerance, check if cohesion allows stick state
        // before deciding to transition to slip (fix for high cohesion materials)
        real64 dLimitTau_dTraction = 0.0;
        real64 const limitTau = contactWrapper.computeLimitTangentialTractionNorm( k,
                                                                                   localTractionNew[ 0 ],
                                                                                   dLimitTau_dTraction );
        real64 const currentTau = LvArray::math::sqrt( localTractionNew[1] * localTractionNew[1] +
                                                       localTractionNew[2] * localTractionNew[2] );

        if( currentTau > limitTau )
        {
          // Tangential traction exceeds cohesion limit: transition to slip
          fractureState[k] = fields::contact::FractureState::Slip;
        }
        // else: cohesion allows stick state, keep traction as computed by updateTraction()
      }

      if( relaxation[k] < 1.0 )
      {
        for( int i = 0; i < 3; ++i )
        {
          localTractionNew[i] = traction[k][i] + relaxation[k] * ( localTractionNew[i] - traction[k][i] );
        }
      }

      LvArray::tensorOps::copy< 3 >( traction[k], localTractionNew );
      penalty[k][2] = -localPenalty[1][1];
      penalty[k][3] = -localPenalty[2][2];
      penalty[k][4] = -localPenalty[1][2];

    } );
  }

};

} // namespace SolidMechanicsALMKernels

} // namespace geos


#endif /* GEOS_PHYSICSSOLVERS_CONTACT_KERNELS_SOLIDMECHANICSALMKERNELSBASE_HPP_ */
