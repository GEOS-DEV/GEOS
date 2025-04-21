
/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

 #ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_EMBEDDEDRUNGEKUTTAKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_EMBEDDEDRUNGEKUTTAKERNELS_HPP_

#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include "physicsSolvers/inducedSeismicity/rateAndStateFields.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "mesh/SurfaceElementSubRegion.hpp"


/**
 *  @file EmbeddedRungeKuttaKernels.hpp
 */

namespace geos
{

namespace rateAndStateKernels
{
/**
 * @brief Butcher table for embedded RK3(2) method using Kuttas third order
 *        method for the high-order update, and an explicit trapezoidal rule
 *        based on the first and third stage rates for the low-order update.
 */
struct Kutta32Table
{
  integer constexpr static algHighOrder = 3;                      // High-order update order
  integer constexpr static algLowOrder = 2;                       // Low-order update order
  integer constexpr static numStages = 3;                         // Number of stages
  real64 const a[2][2] = { { 1.0/2.0, 0.0 },           // Coefficients for stage value updates
    { -1.0, 2.0 } };                                              // (lower-triangular part of table).
  real64 const c[3] = { 0.0, 1.0/2.0, 1.0 };           // Coefficients for time increments of substages
  real64 const b[3] = { 1.0/6.0, 4.0/6.0, 1.0/6.0 };   // Quadrature weights used to step the solution to next time
  real64 const bStar[3] = { 1.0/2.0, 0.0, 1.0/2.0 };   // Quadrature weights used for low-order comparision solution
  real64 constexpr static FSAL = false;                           // Not first same as last
};

/**
 * @brief Butcher table for the BogackiShampine 3(2) method.
 */
struct BogackiShampine32Table
{
  integer constexpr static algHighOrder = 3;                                   // High-order update order
  integer constexpr static algLowOrder = 2;                                    // Low-order update order
  integer constexpr static numStages = 4;                                      // Number of stages
  real64 const a[3][3] = { { 1.0/2.0, 0.0, 0.0     },               // Coefficients for stage value updates
    { 0.0, 3.0/4.0, 0.0     },                                                 // (lower-triangular part of table).
    { 2.0/9.0, 1.0/3.0, 4.0/9.0 } };
  real64 const c[4] = { 0.0, 1.0/2.0, 3.0/4.0, 1.0 };               // Coefficients for time increments of substages
  real64 const b[4] = { 2.0/9.0, 1.0/3.0, 4.0/9.0, 0.0 };           // Quadrature weights used to step the solution to next time
  real64 const bStar[4] = { 7.0/24.0, 1.0/4.0, 1.0/3.0, 1.0/8.0};   // Quadrature weights used for low-order comparision solution
  bool constexpr static FSAL = true;                                           // First same as last (can reuse the last stage rate in next
                                                                               // update)
};

/**
 * @brief Runge-Kutta method used to time integrate slip and state. Uses of a high order
 * update used to integrate the solutions, and a lower order update to estimate the error
 * in the time step.
 *
 * @tparam Butcher table defining the Runge-Kutta method.
 */
template< typename TABLE_TYPE, typename FRICTION_LAW_TYPE >
class EmbeddedRungeKuttaKernel
{

public:
  EmbeddedRungeKuttaKernel( SurfaceElementSubRegion & subRegion,
                            FRICTION_LAW_TYPE const & frictionLaw,
                            TABLE_TYPE butcherTable ):
    m_stateVariable( subRegion.getField< fields::rateAndState::stateVariable >() ),
    m_stateVariable_n( subRegion.getField< fields::rateAndState::stateVariable_n >() ),
    m_slipRate( subRegion.getField< fields::rateAndState::slipRate >() ),
    m_slipVelocity( subRegion.getField< fields::rateAndState::slipVelocity >() ),
    m_slipVelocity_n( subRegion.getField< fields::rateAndState::slipVelocity_n >() ),
    m_deltaSlip( subRegion.getField< fields::contact::deltaSlip >() ),
    m_deltaSlip_n( subRegion.getField< fields::contact::deltaSlip_n >() ),
    m_dispJump( subRegion.getField< fields::contact::dispJump >() ),
    m_dispJump_n( subRegion.getField< fields::contact::dispJump_n >() ),
    m_error( subRegion.getField< fields::rateAndState::error >() ),
    m_stageRates( subRegion.getField< fields::rateAndState::rungeKuttaStageRates >() ),
    m_frictionLaw( frictionLaw.createKernelUpdates() ),
    m_butcherTable( butcherTable )
  {}

  /**
   * @brief Initialize slip and state buffers
   */
  GEOS_HOST_DEVICE
  void initialize( localIndex const k ) const
  {
    LvArray::tensorOps::copy< 2 >( m_slipVelocity[k], m_slipVelocity_n[k] );
    m_slipRate[k] =  LvArray::tensorOps::l2Norm< 2 >( m_slipVelocity_n[k] );
    m_stateVariable[k] = m_stateVariable_n[k];
  }

  /**
   * @brief Re-uses the last stage rate from the previous time step as the first
   * in the next update. Only valid for FSAL (first-same-as-last) Runge-Kutta methods.
   */
  GEOS_HOST_DEVICE
  void updateStageRatesFSAL( localIndex const k ) const
  {
    LvArray::tensorOps::copy< 3 >( m_stageRates[k][0], m_stageRates[k][m_butcherTable.numStages-1] );
  }

  /**
   * @brief Updates the stage rates rates (the right-hand-side of the ODEs for slip and state)
   */
  GEOS_HOST_DEVICE
  void updateStageRates( localIndex const k, integer const stageIndex ) const
  {
    m_stageRates[k][stageIndex][0] =  m_slipVelocity[k][0];
    m_stageRates[k][stageIndex][1] =  m_slipVelocity[k][1];
    m_stageRates[k][stageIndex][2] =  m_frictionLaw.stateEvolution( k, m_slipRate[k], m_stateVariable[k] );
  }

  /**
   * @brief Update stage values (slip, state and displacement jump) to a Runge-Kutta substage.
   */
  GEOS_HOST_DEVICE
  void updateStageValues( localIndex const k, integer const stageIndex, real64 const dt ) const
  {

    real64 stateVariableIncrement = 0.0;
    real64 deltaSlipIncrement[2] = {0.0, 0.0};

    for( integer i = 0; i < stageIndex; i++ )
    {
      deltaSlipIncrement[0] += m_butcherTable.a[stageIndex-1][i] * m_stageRates[k][i][0];
      deltaSlipIncrement[1] += m_butcherTable.a[stageIndex-1][i] * m_stageRates[k][i][1];
      stateVariableIncrement += m_butcherTable.a[stageIndex-1][i] * m_stageRates[k][i][2];
    }
    m_deltaSlip[k][0] = m_deltaSlip_n[k][0] + dt*deltaSlipIncrement[0];
    m_deltaSlip[k][1] = m_deltaSlip_n[k][1] + dt*deltaSlipIncrement[1];
    m_stateVariable[k] = m_stateVariable_n[k] + dt*stateVariableIncrement;

    m_dispJump[k][1] = m_dispJump_n[k][1] + m_deltaSlip[k][0];
    m_dispJump[k][2] = m_dispJump_n[k][2] + m_deltaSlip[k][1];
  }

  /**
   * @brief Updates slip, state and displacement jump to the next time computes error the local error
   * in the time step
   */
  GEOS_HOST_DEVICE
  void updateSolutionAndLocalError( localIndex const k, real64 const dt, real64 const absTol, real64 const relTol ) const
  {

    real64 deltaSlipIncrement[2] = {0.0, 0.0};
    real64 deltaSlipIncrementLowOrder[2] = {0.0, 0.0};

    real64 stateVariableIncrement = 0.0;
    real64 stateVariableIncrementLowOrder = 0.0;

    for( localIndex i = 0; i < m_butcherTable.numStages; i++ )
    {

      // High order update of solution
      deltaSlipIncrement[0]   += m_butcherTable.b[i] * m_stageRates[k][i][0];
      deltaSlipIncrement[1]   += m_butcherTable.b[i] * m_stageRates[k][i][1];
      stateVariableIncrement  += m_butcherTable.b[i] * m_stageRates[k][i][2];

      // Low order update for error
      deltaSlipIncrementLowOrder[0]   += m_butcherTable.bStar[i] * m_stageRates[k][i][0];
      deltaSlipIncrementLowOrder[1]   += m_butcherTable.bStar[i] * m_stageRates[k][i][1];
      stateVariableIncrementLowOrder  += m_butcherTable.bStar[i] * m_stageRates[k][i][2];
    }

    m_deltaSlip[k][0]  = m_deltaSlip_n[k][0]  + dt * deltaSlipIncrement[0];
    m_deltaSlip[k][1]  = m_deltaSlip_n[k][1]  + dt * deltaSlipIncrement[1];
    m_stateVariable[k] = m_stateVariable_n[k] + dt * stateVariableIncrement;

    real64 const deltaSlipLowOrder[2]  = {m_deltaSlip_n[k][0]  + dt * deltaSlipIncrementLowOrder[0],
                                          m_deltaSlip_n[k][1]  + dt * deltaSlipIncrementLowOrder[1]};
    real64 const stateVariableLowOrder = m_stateVariable_n[k] + dt * stateVariableIncrementLowOrder;

    m_dispJump[k][1] = m_dispJump_n[k][1] + m_deltaSlip[k][0];
    m_dispJump[k][2] = m_dispJump_n[k][2] + m_deltaSlip[k][1];

    // Compute error
    m_error[k][0] = computeError( m_deltaSlip[k][0], deltaSlipLowOrder[0], absTol, relTol );
    m_error[k][1] = computeError( m_deltaSlip[k][1], deltaSlipLowOrder[1], absTol, relTol );
    m_error[k][2] = computeError( m_stateVariable[k], stateVariableLowOrder, absTol, relTol );
  }

  /**
   * @brief Updates slip, state and displacement jump to the next time computes error the local error
   * in the time step. Uses the FSAL (first-same-as-last) property.
   */
  GEOS_HOST_DEVICE
  void updateSolutionAndLocalErrorFSAL( localIndex const k, real64 const dt, real64 const absTol, real64 const relTol ) const
  {

    real64 deltaSlipIncrementLowOrder[2] = {0.0, 0.0};
    real64 stateVariableIncrementLowOrder = 0.0;

    for( localIndex i = 0; i < m_butcherTable.numStages; i++ )
    {
      // In FSAL algorithms the last RK substage update coincides with the
      // high-order update. Only need to compute increments for the the
      // low-order updates for error computation.
      deltaSlipIncrementLowOrder[0]   += m_butcherTable.bStar[i] * m_stageRates[k][i][0];
      deltaSlipIncrementLowOrder[1]   += m_butcherTable.bStar[i] * m_stageRates[k][i][1];
      stateVariableIncrementLowOrder  += m_butcherTable.bStar[i] * m_stageRates[k][i][2];
    }

    real64 const deltaSlipLowOrder[2]  = {m_deltaSlip_n[k][0]  + dt * deltaSlipIncrementLowOrder[0],
                                          m_deltaSlip_n[k][1]  + dt * deltaSlipIncrementLowOrder[1]};
    real64 const stateVariableLowOrder = m_stateVariable_n[k] + dt * stateVariableIncrementLowOrder;

    m_dispJump[k][1] = m_dispJump_n[k][1] + m_deltaSlip[k][0];
    m_dispJump[k][2] = m_dispJump_n[k][2] + m_deltaSlip[k][1];

    // Compute error
    m_error[k][0] = computeError( m_deltaSlip[k][0], deltaSlipLowOrder[0], absTol, relTol );
    m_error[k][1] = computeError( m_deltaSlip[k][1], deltaSlipLowOrder[1], absTol, relTol );
    m_error[k][2] = computeError( m_stateVariable[k], stateVariableLowOrder, absTol, relTol );
  }

  /**
   * @brief Computes the relative error scaled by error tolerances
   */
  GEOS_HOST_DEVICE
  real64 computeError( real64 const highOrderApprox, real64 const lowOrderApprox, real64 const absTol, real64 const relTol ) const
  {
    return (highOrderApprox - lowOrderApprox) /
           ( absTol + relTol * LvArray::math::max( LvArray::math::abs( highOrderApprox ), LvArray::math::abs( lowOrderApprox ) ));
  }

private:

  /// Current state variable
  arrayView1d< real64 > const m_stateVariable;

  /// State variable at t = t_n
  arrayView1d< real64 > const m_stateVariable_n;

  /// Current slip rate (magnitude of slip velocity)
  arrayView1d< real64 > const m_slipRate;

  /// Current slip velocity
  arrayView2d< real64 > const m_slipVelocity;

  /// Slip velocity at time t_n
  arrayView2d< real64 > const m_slipVelocity_n;

  /// Current slip change
  arrayView2d< real64 > const m_deltaSlip;

  /// Slip change at time t_n
  arrayView2d< real64 > const m_deltaSlip_n;

  /// Current displacment jump
  arrayView2d< real64 > const m_dispJump;

  /// Displacment jump at time t_n
  arrayView2d< real64 > const m_dispJump_n;

  /// Local error for each solution component stored as slip1, slip2, state
  arrayView2d< real64 > const m_error;

  /// Stage rates for each solution component stored as slip1, slip2, state
  arrayView3d< real64 > const m_stageRates;

  /// Friction law used for rate-and-state updates
  typename FRICTION_LAW_TYPE::KernelWrapper m_frictionLaw;

  /// Butcher table used for explicit time stepping of slip and state
  TABLE_TYPE m_butcherTable;
};

template< typename FRICTION_TYPE, typename BUTCHER_TABLE_TYPE >
void createAndlaunchODEInitialSubStage( SurfaceElementSubRegion & subRegion,
                                        FRICTION_TYPE & frictionLaw,
                                        BUTCHER_TABLE_TYPE const & butcherTable,
                                        real64 const dt,
                                        bool const successfulStep )
{
  rateAndStateKernels::EmbeddedRungeKuttaKernel< BUTCHER_TABLE_TYPE, FRICTION_TYPE > rkKernel( subRegion, frictionLaw, butcherTable );
  if( butcherTable.FSAL && successfulStep )
  {
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      rkKernel.updateStageRatesFSAL( k );
      rkKernel.updateStageValues( k, 1, dt );
    } );
  }
  else
  {
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      rkKernel.initialize( k );
      rkKernel.updateStageRates( k, 0 );
      rkKernel.updateStageValues( k, 1, dt );
    } );
  }
}

template< typename BUTCHER_TABLE_TYPE, typename FRICTION_TYPE >
void createAndlaunchStepRateStateODESubstage( SurfaceElementSubRegion & subRegion,
                                              FRICTION_TYPE & frictionLaw,
                                              BUTCHER_TABLE_TYPE const & butcherTable,
                                              integer const stageIndex,
                                              real64 const dt )
{

  rateAndStateKernels::EmbeddedRungeKuttaKernel< BUTCHER_TABLE_TYPE, FRICTION_TYPE > rkKernel( subRegion, frictionLaw, butcherTable );
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    rkKernel.updateStageRates( k, stageIndex );
    rkKernel.updateStageValues( k, stageIndex+1, dt );
  } );
}

template< typename BUTCHER_TABLE_TYPE, typename FRICTION_TYPE >
void createAndlaunchStepRateStateODEAndComputeError( SurfaceElementSubRegion & subRegion,
                                                     FRICTION_TYPE & frictionLaw,
                                                     BUTCHER_TABLE_TYPE const & butcherTable,
                                                     real64 const relTolerance,
                                                     real64 const absTolerance,
                                                     real64 const dt )
{

  rateAndStateKernels::EmbeddedRungeKuttaKernel< BUTCHER_TABLE_TYPE, FRICTION_TYPE > rkKernel( subRegion, frictionLaw, butcherTable );
  if( butcherTable.FSAL )
  {
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // Perform last stage rate update
      rkKernel.updateStageRates( k, butcherTable.numStages-1 );
      // Update solution to final time and compute errors
      rkKernel.updateSolutionAndLocalErrorFSAL( k, dt, absTolerance, relTolerance );
    } );
  }
  else
  {
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      // Perform last stage rate update
      rkKernel.updateStageRates( k, butcherTable.numStages-1 );
      // Update solution to final time and compute errors
      rkKernel.updateSolutionAndLocalError( k, dt, absTolerance, relTolerance );
    } );
  }
}

} /* namespace rateAndStateKernels */

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_EMBEDDEDRUNGEKUTTAKERNELS_HPP_ */
