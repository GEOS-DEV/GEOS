/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_PHYSICSSOLVERS_RATEANDSTATEKERNELS_HPP_
#define GEOS_PHYSICSSOLVERS_RATEANDSTATEKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include "physicsSolvers/inducedSeismicity/rateAndStateFields.hpp"
#include "denseLinearAlgebra/denseLASolvers.hpp"

namespace geos
{

namespace rateAndStateKernels
{
/**
 * @class RateAndStateKernel
 *
 * @brief
 *
 * @details
 */
class RateAndStateKernel
{
public:

  RateAndStateKernel( SurfaceElementSubRegion & subRegion,
                      constitutive::RateAndStateFriction const & frictionLaw,
                      real64 const shearImpedance ):
    m_slipRate( subRegion.getField< fields::rateAndState::slipRate >() ),
    m_stateVariable( subRegion.getField< fields::rateAndState::stateVariable >() ),
    m_stateVariable_n( subRegion.getField< fields::rateAndState::stateVariable >() ),
    m_traction( subRegion.getField< fields::contact::traction >() ),
    m_shearImpedance( shearImpedance ),
    m_frictionLaw( frictionLaw.createKernelWrapper()  )
  {}

  /**
   * @struct StackVariables
   * @brief Kernel variables located on the stack
   */
  struct StackVariables
  {
public:

    GEOS_HOST_DEVICE
    StackVariables( )
    {}

    real64 jacobian[2][2]{};

    real64 rhs[2]{};

  };

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              real64 const dt,
              StackVariables & stack ) const
  {
    
    // 
    real64 const normalTraction = m_traction[k][0];
    real64 const shearTraction = LvArray::math::sqrt( m_traction[k][1]*m_traction[k][1] + m_traction[k][2]*m_traction[k][2] );

    // Eq 1: shear stress balance
    real64 const tauFriction     = m_frictionLaw.frictionCoefficient( k, m_slipRate[k], m_stateVariable[k] ) * normalTraction;
    real64 const dTauFriction[2] = { m_frictionLaw.dfrictionCoefficient_dStateVariable( k, m_slipRate[k], m_stateVariable[k] ) * normalTraction,
                                     m_frictionLaw.dfrictionCoefficient_dSlipRate( k, m_slipRate[k], m_stateVariable[k] ) * normalTraction };

    stack.rhs[0] = shearTraction - tauFriction - m_shearImpedance * m_slipRate[k];

    // Eq 2: slip law
    stack.rhs[1] = (m_stateVariable[k] - m_stateVariable_n[k]) / dt - m_frictionLaw.dStateVariabledT( k, m_slipRate[k], m_stateVariable[k] );
    real64 const dStateEvolutionLaw[2] = { 1 / dt - m_frictionLaw.dStateVariabledT_dStateVariable( k, m_slipRate[k], m_stateVariable[k] ),
                                           -m_frictionLaw.dStateVariabledT_dSlipRate( k, m_slipRate[k], m_stateVariable[k] ) };

    // Assemble Jacobian matrix
    // derivative shear stress balance w.r.t. theta
    stack.jacobian[0][0] = -dTauFriction[0];
    // derivative shear stress balance w.r.t. slip_velocity
    stack.jacobian[0][1] = -dTauFriction[1];
    // derivative slip law w.r.t. theta
    stack.jacobian[1][0] = dStateEvolutionLaw[0];
    // derivative slip law w.r.t. slip_velocity
    stack.jacobian[1][1] = dStateEvolutionLaw[1];
  }

  GEOS_HOST_DEVICE
  void solve( localIndex const k,
              StackVariables & stack ) const
  {
    /// Solve 2x2 system
    real64 solution[2] = {0.0, 0.0};

    denseLinearAlgebra::solve< 2 >( stack.jacobian, stack.rhs, solution );

    /// Update variables
    m_stateVariable[k] += solution[0];
    m_slipRate[k]      += solution[1];
  }

private:

  arrayView1d< real64 > const m_slipRate;

  arrayView1d< real64 > const m_stateVariable;

  arrayView1d< real64 const > const m_stateVariable_n;

  arrayView2d< real64 const > const m_traction;

  real64 const m_shearImpedance;

  constitutive::RateAndStateFriction::KernelWrapper m_frictionLaw;

};



/**
 * @brief Performs the kernel launch
 * @tparam POLICY the policy used in the RAJA kernels
 */
template< typename POLICY >
static void
createAndLaunch( SurfaceElementSubRegion & subRegion,
                 string const & frictionLawNameKey,
                 real64 const shearImpedance,
                 integer const maxNewtonIter,
                 real64 const time_n,
                 real64 const dt )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n );

  string const & frictionaLawName = subRegion.getReference< string >( frictionLawNameKey );
  constitutive::RateAndStateFriction const & frictionLaw = subRegion.getConstitutiveModel< constitutive::RateAndStateFriction >( frictionaLawName );
  RateAndStateKernel kernel( subRegion, frictionLaw, shearImpedance );

  // Newton loops outside of the kernel launch
  for( integer iter = 0; iter < maxNewtonIter; iter++ )
  {
    /// Kernel 1: Do a solver for all non converged elements
    forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      RateAndStateKernel::StackVariables stack;
      kernel.setup( k, dt, stack );
      kernel.solve( k, stack );
    } );

    /// Kernel 2: Update set of non-converged elements
    // forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    // {

    // } );
  }
}

} /* namespace rateAndStateKernels */

}/* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_RATEANDSTATEKERNELS_HPP_ */
