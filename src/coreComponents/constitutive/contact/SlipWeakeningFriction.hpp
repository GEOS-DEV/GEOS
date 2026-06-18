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
 * @file SlipWeakeningFriction.hpp
 */

#ifndef GEOS_CONSTITUTIVE_CONTACT_SLIPWEAKENINGFRICTION_HPP_
#define GEOS_CONSTITUTIVE_CONTACT_SLIPWEAKENINGFRICTION_HPP_

#include "FrictionBase.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "physicsSolvers/solidMechanics/contact/FractureState.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

namespace constitutive
{

/**
 * @class SlipWeakeningFrictionUpdates
 *
 * In-kernel wrapper for slip-weakening friction updates.
 * The friction coefficient decreases linearly from muPeak to muResidual
 * over a characteristic slip distance Dc:
 *   mu(delta) = muResidual + (muPeak - muResidual) * max(0, 1 - delta/Dc)
 */
class SlipWeakeningFrictionUpdates : public FrictionBaseUpdates
{
public:

  SlipWeakeningFrictionUpdates( real64 const & displacementJumpThreshold,
                                real64 const & shearStiffness,
                                arrayView1d< real64 const > const & cohesion,
                                arrayView1d< real64 const > const & initialFrictionCoefficient,
                                arrayView1d< real64 const > const & residualFrictionCoefficient,
                                arrayView1d< real64 const > const & Dc,
                                arrayView1d< real64 > const & cumulativeSlip,
                                arrayView1d< real64 > const & cumulativeSlipSaved,
                                arrayView2d< real64 > const & elasticSlip )
    : FrictionBaseUpdates( displacementJumpThreshold ),
    m_shearStiffness( shearStiffness ),
    m_cohesion( cohesion ),
    m_initialFrictionCoefficient( initialFrictionCoefficient ),
    m_residualFrictionCoefficient( residualFrictionCoefficient ),
    m_Dc( Dc ),
    m_cumulativeSlip( cumulativeSlip ),
    m_cumulativeSlipSaved( cumulativeSlipSaved ),
    m_elasticSlip( elasticSlip )
  {}

  /// Default copy constructor
  SlipWeakeningFrictionUpdates( SlipWeakeningFrictionUpdates const & ) = default;

  /// Default move constructor
  SlipWeakeningFrictionUpdates( SlipWeakeningFrictionUpdates && ) = default;

  /// Deleted default constructor
  SlipWeakeningFrictionUpdates() = delete;

  /// Deleted copy assignment operator
  SlipWeakeningFrictionUpdates & operator=( SlipWeakeningFrictionUpdates const & ) = delete;

  /// Deleted move assignment operator
  SlipWeakeningFrictionUpdates & operator=( SlipWeakeningFrictionUpdates && ) = delete;

  GEOS_HOST_DEVICE
  inline
  virtual real64 computeLimitTangentialTractionNorm( localIndex const k,
                                                     real64 const & normalTraction,
                                                     real64 & dLimitTangentialTractionNorm_dTraction ) const override final;

  GEOS_HOST_DEVICE
  inline
  virtual void computeShearTraction( localIndex const k,
                                     arraySlice1d< real64 const > const & oldDispJump,
                                     arraySlice1d< real64 const > const & dispJump,
                                     integer const & fractureState,
                                     arraySlice1d< real64 > const & tractionVector,
                                     arraySlice2d< real64 > const & dTractionVector_dJump ) const override final;

  GEOS_HOST_DEVICE
  inline
  virtual void updateFractureState( localIndex const k,
                                    arraySlice1d< real64 const > const & dispJump,
                                    arraySlice1d< real64 const > const & tractionVector,
                                    integer & fractureState ) const override final;

  GEOS_HOST_DEVICE
  inline
  virtual void updateElasticSlip( localIndex const k,
                                  arraySlice1d< real64 const > const & dispJump,
                                  arraySlice1d< real64 const > const & oldDispJump,
                                  arraySlice1d< real64 const > const & tractionVector,
                                  integer const & fractureState ) const override final;

  GEOS_HOST_DEVICE
  inline
  virtual void saveState( localIndex const k ) const override final
  {
    m_cumulativeSlipSaved[k] = m_cumulativeSlip[k];
  }

  GEOS_HOST_DEVICE
  inline
  virtual void updateTraction( localIndex const k,
                               arraySlice1d< real64 const > const & oldDispJump,
                               arraySlice1d< real64 const > const & dispJump,
                               arraySlice1d< real64 const > const & penalty,
                               arraySlice1d< real64 const > const & traction,
                               bool const symmetric,
                               bool const fixedLimitTau,
                               real64 const normalTractionTolerance,
                               real64 const tangentialTractionTolerance,
                               real64 ( & dTraction_dDispJump )[3][3],
                               real64 ( & tractionNew )[3],
                               integer & fractureState ) const override final;

  GEOS_HOST_DEVICE
  inline
  virtual void updateTractionOnly( localIndex const k,
                                   arraySlice1d< real64 const > const & dispJump,
                                   arraySlice1d< real64 const > const & deltaDispJump,
                                   arraySlice1d< real64 const > const & penalty,
                                   arraySlice1d< real64 const > const & traction,
                                   arraySlice1d< real64 > const & tractionNew ) const override final;

  GEOS_HOST_DEVICE
  inline
  virtual void constraintCheck( localIndex const k,
                                arraySlice1d< real64 const > const & dispJump,
                                arraySlice1d< real64 const > const & deltaDispJump,
                                arraySlice1d< real64 > const & tractionVector,
                                integer const fractureState,
                                real64 const normalTractionTolerance,
                                real64 const normalDisplacementTolerance,
                                real64 const slidingTolerance,
                                real64 const slidingCheckTolerance,
                                integer & condConv ) const override final;

private:

  /// Linearly interpolated friction coefficient based on cumulative slip
  GEOS_HOST_DEVICE
  inline
  real64 currentFrictionCoefficient( localIndex const k ) const
  {
    real64 const normalizedSlip = m_cumulativeSlip[k] / m_Dc[k];
    real64 const weight = ( normalizedSlip < 1.0 ) ? ( 1.0 - normalizedSlip ) : 0.0;
    return m_residualFrictionCoefficient[k] + ( m_initialFrictionCoefficient[k] - m_residualFrictionCoefficient[k] ) * weight;
  }

  /// The shear stiffness (used for elastic slip tracking and explicit solver path)
  real64 m_shearStiffness;

  arrayView1d< real64 const > m_cohesion;
  arrayView1d< real64 const > m_initialFrictionCoefficient;
  arrayView1d< real64 const > m_residualFrictionCoefficient;
  arrayView1d< real64 const > m_Dc;

  /// Cumulative slip (state variable; semi-implicitly updated during Newton in updateTraction)
  arrayView1d< real64 > m_cumulativeSlip;

  /// Snapshot of cumulativeSlip at the start of the current time step (written by saveState)
  arrayView1d< real64 > m_cumulativeSlipSaved;

  arrayView2d< real64 > m_elasticSlip;
};


/**
 * @class SlipWeakeningFriction
 *
 * Slip-weakening friction constitutive model.
 * mu decreases linearly from muPeak to muResidual over distance Dc.
 */
class SlipWeakeningFriction : public FrictionBase
{
public:

  SlipWeakeningFriction( string const & name, Group * const parent );

  static string catalogName() { return "SlipWeakening"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                        localIndex const numPts ) override final;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = SlipWeakeningFrictionUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates() const;

  struct viewKeyStruct : public FrictionBase::viewKeyStruct
  {
    static constexpr char const * shearStiffnessString() { return "shearStiffness"; }
    static constexpr char const * elasticSlipString() { return "elasticSlip"; }
    static constexpr char const * defaultCohesionString() { return "defaultCohesion"; }
    static constexpr char const * defaultInitialFrictionCoefficientString() { return "defaultInitialFrictionCoefficient"; }
    static constexpr char const * defaultResidualFrictionCoefficientString() { return "defaultResidualFrictionCoefficient"; }
    static constexpr char const * defaultDcString() { return "defaultDc"; }
  };

protected:

  virtual void postInputInitialization() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

private:

  real64 m_shearStiffness;

  array1d< real64 > m_cohesion;
  array1d< real64 > m_initialFrictionCoefficient;
  array1d< real64 > m_residualFrictionCoefficient;
  array1d< real64 > m_Dc;
  array1d< real64 > m_cumulativeSlip;
  array1d< real64 > m_cumulativeSlipSaved;
  array2d< real64 > m_elasticSlip;

  real64 m_defaultCohesion;
  real64 m_defaultInitialFrictionCoefficient;
  real64 m_defaultResidualFrictionCoefficient;
  real64 m_defaultDc;
};


// ============================================================
// Inline implementations
// ============================================================

GEOS_HOST_DEVICE
real64 SlipWeakeningFrictionUpdates::computeLimitTangentialTractionNorm(
  localIndex const k,
  real64 const & normalTraction,
  real64 & dLimitTangentialTractionNorm_dTraction ) const
{
  real64 const mu = currentFrictionCoefficient( k );
  dLimitTangentialTractionNorm_dTraction = -mu;
  return ( m_cohesion[k] - normalTraction * mu );
}


GEOS_HOST_DEVICE
inline void SlipWeakeningFrictionUpdates::computeShearTraction(
  localIndex const k,
  arraySlice1d< real64 const > const & oldDispJump,
  arraySlice1d< real64 const > const & dispJump,
  integer const & fractureState,
  arraySlice1d< real64 > const & tractionVector,
  arraySlice2d< real64 > const & dTractionVector_dJump ) const
{
  using namespace fields::contact;

  real64 const slip[2] = { dispJump[1] - oldDispJump[1],
                           dispJump[2] - oldDispJump[2] };

  switch( fractureState )
  {
    case FractureState::Stick:
    {
      real64 const tau[2] = { m_shearStiffness * ( slip[0] + m_elasticSlip[k][0] ),
                              m_shearStiffness * ( slip[1] + m_elasticSlip[k][1] ) };

      tractionVector[1] = tau[0];
      tractionVector[2] = tau[1];

      dTractionVector_dJump[1][1] = m_shearStiffness;
      dTractionVector_dJump[2][2] = m_shearStiffness;

      break;
    }
    case FractureState::Slip:
    {
      real64 dLimitTau_dNormalTraction;
      real64 const limitTau = computeLimitTangentialTractionNorm( k, tractionVector[0],
                                                                  dLimitTau_dNormalTraction );

      real64 const slipNorm = LvArray::tensorOps::l2Norm< 2 >( slip );

      tractionVector[1] = limitTau * slip[0] / slipNorm;
      tractionVector[2] = limitTau * slip[1] / slipNorm;

      dTractionVector_dJump[1][0] = dTractionVector_dJump[0][0] * dLimitTau_dNormalTraction * slip[0] / slipNorm;
      dTractionVector_dJump[1][1] = limitTau * pow( slip[1], 2 )  / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );
      dTractionVector_dJump[1][2] = -limitTau * slip[0] * slip[1] / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );

      dTractionVector_dJump[2][0] = dTractionVector_dJump[0][0] * dLimitTau_dNormalTraction * slip[1] / slipNorm;
      dTractionVector_dJump[2][1] = -limitTau * slip[0] * slip[1] / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );
      dTractionVector_dJump[2][2] = limitTau * pow( slip[0], 2 )  / pow( LvArray::tensorOps::l2NormSquared< 2 >( slip ), 1.5 );

      break;
    }
  }
}


GEOS_HOST_DEVICE
inline void SlipWeakeningFrictionUpdates::updateFractureState(
  localIndex const k,
  arraySlice1d< real64 const > const & dispJump,
  arraySlice1d< real64 const > const & tractionVector,
  integer & fractureState ) const
{
  using namespace fields::contact;

  if( dispJump[0] > -m_displacementJumpThreshold )
  {
    fractureState = FractureState::Open;
  }
  else
  {
    real64 const tau[2] = { tractionVector[1], tractionVector[2] };
    real64 const tauNorm = LvArray::tensorOps::l2Norm< 2 >( tau );

    real64 dLimitTau_dNormalTraction;
    real64 const limitTau = computeLimitTangentialTractionNorm( k, tractionVector[0],
                                                                dLimitTau_dNormalTraction );

    real64 const yield = tauNorm - limitTau;

    if( yield < 0 )
    {
      fractureState = FractureState::Stick;
    }
    else
    {
      fractureState = FractureState::Slip;
    }
  }
}


GEOS_HOST_DEVICE
inline void SlipWeakeningFrictionUpdates::updateElasticSlip(
  localIndex const k,
  arraySlice1d< real64 const > const & dispJump,
  arraySlice1d< real64 const > const & oldDispJump,
  arraySlice1d< real64 const > const & tractionVector,
  integer const & fractureState ) const
{
  using namespace fields::contact;

  if( fractureState == FractureState::Open )
  {
    m_elasticSlip[k][0] = 0.0;
    m_elasticSlip[k][1] = 0.0;
  }
  else if( fractureState == FractureState::Stick )
  {
    real64 const slip[2] = { dispJump[1] - oldDispJump[1],
                             dispJump[2] - oldDispJump[2] };
    LvArray::tensorOps::add< 2 >( m_elasticSlip[k], slip );
  }
  else if( fractureState == FractureState::Slip )
  {
    // Finalize cumulativeSlip using the converged increment (idempotent with updateTraction's update)
    real64 const slipIncrement[2] = { dispJump[1] - oldDispJump[1],
                                      dispJump[2] - oldDispJump[2] };
    m_cumulativeSlip[k] = m_cumulativeSlipSaved[k] + LvArray::tensorOps::l2Norm< 2 >( slipIncrement );

    m_elasticSlip[k][0] = tractionVector[1] / m_shearStiffness;
    m_elasticSlip[k][1] = tractionVector[2] / m_shearStiffness;
  }
}


GEOS_HOST_DEVICE
inline void SlipWeakeningFrictionUpdates::updateTraction(
  localIndex const k,
  arraySlice1d< real64 const > const & oldDispJump,
  arraySlice1d< real64 const > const & dispJump,
  arraySlice1d< real64 const > const & penalty,
  arraySlice1d< real64 const > const & traction,
  bool const symmetric,
  bool const fixedLimitTau,
  real64 const normalTractionTolerance,
  real64 const tangentialTractionTolerance,
  real64 ( & dTraction_dDispJump )[3][3],
  real64 ( & tractionNew )[3],
  integer & fractureState ) const
{
  using namespace fields::contact;

  real64 dLimitTangentialTractionNorm_dTraction = 0.0;
  real64 limitTau = 0.0;

  // Compute trial traction
  real64 tractionTrial[3];
  tractionTrial[0] = traction[0] + penalty[0] * dispJump[0];
  tractionTrial[1] = traction[1] + penalty[1] * ( dispJump[1] - oldDispJump[1] );
  tractionTrial[2] = traction[2] + penalty[1] * ( dispJump[2] - oldDispJump[2] );

  real64 const tau[2] = { tractionTrial[1], tractionTrial[2] };
  real64 const tractionTrialNorm = LvArray::tensorOps::l2Norm< 2 >( tau );

  fractureState = FractureState::Stick;
  if( tractionTrial[0] > normalTractionTolerance )
  {
    tractionNew[0] = 0.0;
    dTraction_dDispJump[0][0] = 0.0;
    fractureState = FractureState::Open;
  }
  else
  {
    tractionNew[0] = tractionTrial[0];
    dTraction_dDispJump[0][0] = -penalty[0];
  }

  // Semi-implicit update: estimate cumulativeSlip for this Newton iterate using
  // the current tangential increment (dispJump - oldDispJump).  This eliminates
  // the one-step lag that arises when cumulativeSlip is frozen at the previous step.
  // m_cumulativeSlipSaved holds the converged value from the END of the previous step.
  // If the fracture turns out to be Stick or Open, we revert to saved at the end.
  {
    real64 const slipInc[2] = { dispJump[1] - oldDispJump[1], dispJump[2] - oldDispJump[2] };
    m_cumulativeSlip[k] = m_cumulativeSlipSaved[k] + LvArray::tensorOps::l2Norm< 2 >( slipInc );
  }

  // Compute limit tau using the updated mu(delta)
  if( fixedLimitTau )
  {
    limitTau = computeLimitTangentialTractionNorm( k, traction[0],
                                                   dLimitTangentialTractionNorm_dTraction );
  }
  else
  {
    limitTau = computeLimitTangentialTractionNorm( k, tractionNew[0],
                                                   dLimitTangentialTractionNorm_dTraction );
  }

  if( tractionTrialNorm <= tangentialTractionTolerance )
  {
    dTraction_dDispJump[1][1] = -penalty[1];
    dTraction_dDispJump[2][2] = -penalty[1];

    tractionNew[1] = tractionTrial[1];
    tractionNew[2] = tractionTrial[2];

    if( fractureState != FractureState::Open )
    {
      fractureState = FractureState::Stick;
    }
  }
  else if( limitTau <= tangentialTractionTolerance )
  {
    dTraction_dDispJump[1][1] = 0.0;
    dTraction_dDispJump[2][2] = 0.0;

    tractionNew[1] = ( fixedLimitTau ) ? tractionTrial[1] : 0.0;
    tractionNew[2] = ( fixedLimitTau ) ? tractionTrial[2] : 0.0;

    if( fractureState != FractureState::Open )
    {
      fractureState = FractureState::Slip;
    }
  }
  else
  {
    real64 const psi  = ( tractionTrialNorm > limitTau ) ? 1.0 : tractionTrialNorm / limitTau;
    real64 const dpsi = ( tractionTrialNorm > limitTau ) ? 0.0 : 1.0;

    if( fractureState != FractureState::Open )
    {
      fractureState = ( tractionTrialNorm > limitTau ) ? FractureState::Slip : FractureState::Stick;
    }

    real64 dNormTTdgT[3];
    dNormTTdgT[0] = tractionTrial[1] * tractionTrial[1];
    dNormTTdgT[1] = tractionTrial[2] * tractionTrial[2];
    dNormTTdgT[2] = tractionTrial[1] * tractionTrial[2];

    real64 dTdgT[3];
    dTdgT[0] = ( tractionTrialNorm * tractionTrialNorm - dNormTTdgT[0] );
    dTdgT[1] = ( tractionTrialNorm * tractionTrialNorm - dNormTTdgT[1] );
    dTdgT[2] = -dNormTTdgT[2];

    LvArray::tensorOps::scale< 3 >( dNormTTdgT, 1. / std::pow( tractionTrialNorm, 2 ) );
    LvArray::tensorOps::scale< 3 >( dTdgT, 1. / std::pow( tractionTrialNorm, 3 ) );

    dTraction_dDispJump[1][1] = -penalty[1] * ( dpsi * dNormTTdgT[0] + psi * dTdgT[0] * limitTau );
    dTraction_dDispJump[2][2] = -penalty[1] * ( dpsi * dNormTTdgT[1] + psi * dTdgT[1] * limitTau );
    dTraction_dDispJump[1][2] = -penalty[1] * ( dpsi * dNormTTdgT[2] + psi * dTdgT[2] * limitTau );
    dTraction_dDispJump[2][1] = dTraction_dDispJump[1][2];

    if( !symmetric )
    {
      // dLimitTangentialTractionNorm_dTraction = -mu(delta)
      real64 const mu = -dLimitTangentialTractionNorm_dTraction;
      dTraction_dDispJump[1][0] = -dTraction_dDispJump[0][0] * mu *
                                  tractionTrial[1] * ( psi / tractionTrialNorm - dpsi / limitTau );
      dTraction_dDispJump[2][0] = -dTraction_dDispJump[0][0] * mu *
                                  tractionTrial[2] * ( psi / tractionTrialNorm - dpsi / limitTau );
    }

    LvArray::tensorOps::scale< 3 >( tractionTrial, ( psi * limitTau ) / tractionTrialNorm );
    tractionNew[1] = tractionTrial[1];
    tractionNew[2] = tractionTrial[2];
  }

  // Revert cumulativeSlip for non-slip states: weakening only occurs under irreversible slip.
  if( fractureState != FractureState::Slip )
  {
    m_cumulativeSlip[k] = m_cumulativeSlipSaved[k];
  }
}


GEOS_HOST_DEVICE
inline void SlipWeakeningFrictionUpdates::updateTractionOnly(
  localIndex const k,
  arraySlice1d< real64 const > const & dispJump,
  arraySlice1d< real64 const > const & deltaDispJump,
  arraySlice1d< real64 const > const & penalty,
  arraySlice1d< real64 const > const & traction,
  arraySlice1d< real64 > const & tractionNew ) const
{
  real64 const zero = LvArray::NumericLimits< real64 >::epsilon;

  // deltaDispJump = dispJump - oldDispJump (total increment since start of time step,
  // freshly computed by computeDispJump before this call). Update cumulativeSlip here
  // so that limitTau reflects the current converged slip, not the stale value from
  // the last Newton assembly where m_dispJump had not yet been updated.
  {
    real64 const slipInc[2] = { deltaDispJump[1], deltaDispJump[2] };
    m_cumulativeSlip[k] = m_cumulativeSlipSaved[k] + LvArray::tensorOps::l2Norm< 2 >( slipInc );
  }

  tractionNew[0] = traction[0] + penalty[0] * dispJump[0];
  tractionNew[1] = traction[1] + penalty[1] * deltaDispJump[1];
  tractionNew[2] = traction[2] + penalty[1] * deltaDispJump[2];

  real64 const tau[2] = { tractionNew[1], tractionNew[2] };
  real64 const currentTau = LvArray::tensorOps::l2Norm< 2 >( tau );

  real64 dLimitTangentialTractionNorm_dTraction = 0.0;
  real64 const limitTau = computeLimitTangentialTractionNorm( k, tractionNew[0],
                                                              dLimitTangentialTractionNorm_dTraction );

  real64 psi;
  if( limitTau < zero )
  {
    psi = 1.0;
  }
  else
  {
    psi = ( currentTau > limitTau ) ? 1.0 : currentTau / limitTau;
  }

  if( limitTau > zero && currentTau > zero )
  {
    tractionNew[1] *= limitTau * psi / currentTau;
    tractionNew[2] *= limitTau * psi / currentTau;
  }
  else
  {
    tractionNew[1] = 0.0;
    tractionNew[2] = 0.0;
  }
}


GEOS_HOST_DEVICE
inline void SlipWeakeningFrictionUpdates::constraintCheck(
  localIndex const k,
  arraySlice1d< real64 const > const & dispJump,
  arraySlice1d< real64 const > const & deltaDispJump,
  arraySlice1d< real64 > const & tractionVector,
  integer const fractureState,
  real64 const normalTractionTolerance,
  real64 const normalDisplacementTolerance,
  real64 const slidingTolerance,
  real64 const slidingCheckTolerance,
  integer & condConv ) const
{
  using namespace fields::contact;

  // Update cumulativeSlip from the current converged deltaDispJump so that
  // computeLimitTangentialTractionNorm uses the correct mu. This is critical
  // for the simultaneous ALM path where the assembly kernel never calls
  // updateTraction, leaving m_cumulativeSlip frozen at m_cumulativeSlipSaved.
  // Only update for Slip state: in Stick, deltaDispJump is the penalty residual
  // (not real slip) and would cause spurious cumulativeSlip accumulation.
  if( fractureState == FractureState::Slip )
  {
    real64 const slipInc[2] = { deltaDispJump[1], deltaDispJump[2] };
    m_cumulativeSlip[k] = m_cumulativeSlipSaved[k] + LvArray::tensorOps::l2Norm< 2 >( slipInc );
  }

  real64 const deltaDisp[2] = { deltaDispJump[1], deltaDispJump[2] };
  real64 const deltaDispNorm = LvArray::tensorOps::l2Norm< 2 >( deltaDisp );

  real64 const tau[2] = { tractionVector[1], tractionVector[2] };
  real64 const currentTau = LvArray::tensorOps::l2Norm< 2 >( tau );

  real64 dLimitTangentialTractionNorm_dTraction = 0.0;
  real64 const limitTau = computeLimitTangentialTractionNorm( k, tractionVector[0],
                                                              dLimitTangentialTractionNorm_dTraction );

  condConv = 0;

  if( tractionVector[0] >= normalTractionTolerance )
  {
    if( fractureState != FractureState::Open )
    {
      condConv = 1;
    }
    tractionVector[0] = 0.0;
    tractionVector[1] = 0.0;
    tractionVector[2] = 0.0;
  }
  else
  {
    if( ( LvArray::math::abs( dispJump[0] ) > normalDisplacementTolerance ) &&
        ( fractureState != FractureState::Open ) )
    {
      condConv = 2;
    }

    if( fractureState == FractureState::Stick &&
        deltaDispNorm > slidingTolerance )
    {
      condConv = 3;
    }

    if( currentTau > ( LvArray::math::abs( limitTau ) * ( 1.0 + slidingCheckTolerance ) ) )
    {
      condConv = 4;
    }
  }
}

} /* namespace constitutive */

} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_CONTACT_SLIPWEAKENINGFRICTION_HPP_ */
