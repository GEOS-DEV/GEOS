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


/**
 *  @file ExplicitRateAndStateKernels.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_EXPLICITRATEANDSTATEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_EXPLICITRATEANDSTATEKERNELS_HPP_

#include "RateAndStateKernelsBase.hpp"
#include "denseLinearAlgebra/denseLASolvers.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"

namespace geos
{

namespace rateAndStateKernels
{

/**
 * @class ExplicitRateAndStateKernel
 *
 * @brief
 *
 * @details
 */
template< typename FRICTION_LAW_TYPE >
class ExplicitRateAndStateKernel
{
public:

  ExplicitRateAndStateKernel( SurfaceElementSubRegion & subRegion,
                              FRICTION_LAW_TYPE const & frictionLaw,
                              real64 const shearImpedance ):
    m_slipRate( subRegion.getField< fields::rateAndState::slipRate >() ),
    m_stateVariable( subRegion.getField< fields::rateAndState::stateVariable >() ),
    m_normalTraction( subRegion.getField< fields::rateAndState::normalTraction >() ),
    m_shearTraction( subRegion.getField< fields::rateAndState::shearTraction >() ),
    m_slipVelocity( subRegion.getField< fields::rateAndState::slipVelocity >() ),
    m_shearImpedance( shearImpedance ),
    m_frictionLaw( frictionLaw.createKernelUpdates()  )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables
  {
public:

    StackVariables() = default;

    real64 jacobian{};
    real64 rhs{};

  };

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              real64 const dt,
              StackVariables & stack ) const
  {
    GEOS_UNUSED_VAR( dt );
    real64 const normalTraction = m_normalTraction[k];
    real64 const shearTractionMagnitude = LvArray::tensorOps::l2Norm< 2 >( m_shearTraction[k] );

    // Slip rate is bracketed between [0, shear traction magnitude / shear impedance]
    // If slip rate is outside the bracket, re-initialize to the middle value
    real64 const upperBound = shearTractionMagnitude/m_shearImpedance;
    real64 const bracketedSlipRate = m_slipRate[k] > upperBound ? 0.5*upperBound : m_slipRate[k];

    stack.rhs = shearTractionMagnitude - m_shearImpedance *bracketedSlipRate - normalTraction * m_frictionLaw.frictionCoefficient( k, bracketedSlipRate, m_stateVariable[k] );
    stack.jacobian = -m_shearImpedance - normalTraction * m_frictionLaw.dFrictionCoefficient_dSlipRate( k, bracketedSlipRate, m_stateVariable[k] );
  }

  GEOS_HOST_DEVICE
  void solve( localIndex const k,
              StackVariables & stack ) const
  {
    m_slipRate[k] -= stack.rhs/stack.jacobian;

    // Slip rate is bracketed between [0, shear traction magnitude / shear impedance]
    // Check that the update did not end outside of the bracket.
    real64 const shearTractionMagnitude = LvArray::tensorOps::l2Norm< 2 >( m_shearTraction[k] );
    real64 const upperBound = shearTractionMagnitude/m_shearImpedance;
    if( m_slipRate[k] > upperBound ) m_slipRate[k] = 0.5*upperBound;

  }


  GEOS_HOST_DEVICE
  camp::tuple< int, real64 > checkConvergence( StackVariables const & stack,
                                               real64 const tol ) const
  {
    real64 const residualNorm = LvArray::math::abs( stack.rhs );
    int const converged = residualNorm < tol ? 1 : 0;
    camp::tuple< int, real64 > result { converged, residualNorm };
    return result;
  }

  GEOS_HOST_DEVICE
  void projectSlipRate( localIndex const k ) const
  {
    real64 const frictionCoefficient = m_frictionLaw.frictionCoefficient( k, m_slipRate[k], m_stateVariable[k] );
    projectSlipRateBase( k, frictionCoefficient, m_shearImpedance, m_normalTraction, m_shearTraction, m_slipRate, m_slipVelocity );
  }

  GEOS_HOST_DEVICE
  void udpateVariables( localIndex const k ) const
  {
    projectSlipRate( k );
  }

  GEOS_HOST_DEVICE
  void resetState( localIndex const k ) const
  {
    GEOS_UNUSED_VAR( k );
  }

  /**
   * @brief Performs the kernel launch
   * @tparam KernelType The Rate-and-state kernel to launch
   * @tparam POLICY the policy used in the RAJA kernels
   */
  template< typename POLICY >
  static real64
  solveRateAndStateEquation( SurfaceElementSubRegion & subRegion,
                             ExplicitRateAndStateKernel & kernel,
                             real64 dt,
                             integer const maxNewtonIter,
                             real64 const newtonTol )
  {
    GEOS_MARK_FUNCTION;

    newtonSolve< POLICY >( subRegion, kernel, dt, maxNewtonIter, newtonTol );

    forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      kernel.projectSlipRate( k );
    } );
    return dt;
  }

private:

  arrayView1d< real64 > const m_slipRate;

  arrayView1d< real64 > const m_stateVariable;

  arrayView1d< real64 const > const m_normalTraction;

  arrayView2d< real64 const > const m_shearTraction;

  arrayView2d< real64 > const m_slipVelocity;

  real64 const m_shearImpedance;

  typename FRICTION_LAW_TYPE::KernelWrapper m_frictionLaw;

};



} /* namespace rateAndStateKernels */

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_EXPLICITRATEANDSTATEKERNELS_HPP_ */
