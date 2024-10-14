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
    m_stateVariable_n( subRegion.getField< fields::rateAndState::stateVariable_n >() ),
    m_traction( subRegion.getField< fields::contact::traction >() ),
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

    GEOS_HOST_DEVICE
    StackVariables( )
    {}

    real64 jacobian[3][3]{};

    real64 rhs[3]{};

  };

  GEOS_HOST_DEVICE
  void setup( localIndex const k,
              real64 const dt,
              StackVariables & stack ) const
  {
    m_slipRate[k][1] = 0.; //TODO. Remove once solution is correctly initialized.
    real64 const normalTraction = m_traction[k][0];
    real64 const shearTraction[2] = { m_traction[k][1],
                                      m_traction[k][2]};
    real64 const slipRateMagnitude = LvArray::math::sqrt(m_slipRate[k][0] * m_slipRate[k][0] + m_slipRate[k][1] * m_slipRate[k][1] );
    real64 const normalizedSlipRate[2] = { m_slipRate[k][0] / slipRateMagnitude,
                                           m_slipRate[k][1] / slipRateMagnitude };
    // std::cout << "normalTraction: " << normalTraction << std::endl;
    // std::cout << "shearTraction[0]: " << shearTraction[0] << std::endl;
    // std::cout << "shearTraction[1]: " << shearTraction[1] << std::endl;

    // std::cout << "normalizedSlipRate[0]: " << normalizedSlipRate[0] << std::endl;
    // std::cout << "normalizedSlipRate[1]: " << normalizedSlipRate[1] << std::endl;

    // Eq 1: shear stress balance
    // Frictional shear strength
    real64 const frictionForce = normalTraction * m_frictionLaw.frictionCoefficient( k, slipRateMagnitude, m_stateVariable[k] );
    real64 const tauFriction[2]   = { frictionForce * normalizedSlipRate[0],
                                      frictionForce * normalizedSlipRate[1] };
    
    // Jacobian entries
    real64 const dFrictionForce_dStateVariable = normalTraction * m_frictionLaw.dFrictionCoefficient_dStateVariable( k, slipRateMagnitude, m_stateVariable[k] );
    real64 const dFrictionForce_dSlipRateMagnitude = normalTraction * m_frictionLaw.dFrictionCoefficient_dSlipRate( k, slipRateMagnitude, m_stateVariable[k] );
    
    // Tangential component 1
    real64 const dTauFriction1_dStateVariable = dFrictionForce_dStateVariable * normalizedSlipRate[0];
    
    real64 const dTauFriction1_dSlipRate1 = dFrictionForce_dSlipRateMagnitude * normalizedSlipRate[0] * normalizedSlipRate[0]
                                          + frictionForce / slipRateMagnitude * normalizedSlipRate[1] * normalizedSlipRate[1];

    real64 const dTauFriction1_dSlipRate2 = ( dFrictionForce_dSlipRateMagnitude + frictionForce / slipRateMagnitude ) 
                                          * ( normalizedSlipRate[0] * normalizedSlipRate[1]);

    real64 const dTauFriction1[3] = { dTauFriction1_dStateVariable, dTauFriction1_dSlipRate1, dTauFriction1_dSlipRate2 };

    // Tangential component 2
    real64 const dTauFriction2_dStateVariable = dFrictionForce_dStateVariable * normalizedSlipRate[1];
    
    real64 const dTauFriction2_dSlipRate1 = dTauFriction1_dSlipRate2;

    real64 const dTauFriction2_dSlipRate2 = dFrictionForce_dSlipRateMagnitude * normalizedSlipRate[1] * normalizedSlipRate[1]
                                          + frictionForce / slipRateMagnitude * normalizedSlipRate[0] * normalizedSlipRate[0];

    real64 const dTauFriction2[3] = { dTauFriction2_dStateVariable, dTauFriction2_dSlipRate1, dTauFriction2_dSlipRate2 };


    // std::cout << "force balance" << std::endl;
    stack.rhs[0] = shearTraction[0] - tauFriction[0] - m_shearImpedance * m_slipRate[k][0];

    stack.rhs[1] = shearTraction[1] - tauFriction[1] - m_shearImpedance * m_slipRate[k][1];

    // Eq 2: slip law
    stack.rhs[2] = (m_stateVariable[k] - m_stateVariable_n[k]) / dt - m_frictionLaw.stateEvolution( k, slipRateMagnitude, m_stateVariable[k] );
    real64 const dStateEvolutionLaw[3] = { 1 / dt - m_frictionLaw.dStateEvolution_dStateVariable( k, slipRateMagnitude, m_stateVariable[k] ),
                                           -m_frictionLaw.dStateEvolution_dSlipRate( k, slipRateMagnitude, m_stateVariable[k] ) * normalizedSlipRate[0], 
                                           -m_frictionLaw.dStateEvolution_dSlipRate( k, slipRateMagnitude, m_stateVariable[k] ) * normalizedSlipRate[1]} ;
    

    // Assemble Jacobian matrix
    // derivative shear stress balance component 1 w.r.t. theta
    stack.jacobian[0][0] = -dTauFriction1[0];
    // derivative shear stress balance component 1 w.r.t. slip_velocity component 1
    stack.jacobian[0][1] = -dTauFriction1[1]- m_shearImpedance;
    // derivative shear stress balance component 1 w.r.t. slip_velocity component 2
    stack.jacobian[0][2] = -dTauFriction1[2];
    // derivative shear stress balance component 2 w.r.t. theta
    stack.jacobian[1][0] = -dTauFriction2[0];
    // derivative shear stress balance component 2 w.r.t. slip_velocity component 1
    stack.jacobian[1][1] = -dTauFriction2[1];
    // derivative shear stress balance component 2 w.r.t. slip_velocity component 2
    stack.jacobian[1][2] = -dTauFriction2[2] - m_shearImpedance;
    // derivative slip law w.r.t. theta
    stack.jacobian[2][0] = dStateEvolutionLaw[0];
    // derivative slip law w.r.t. slip_velocity component 1
    stack.jacobian[2][1] = dStateEvolutionLaw[1];
    // derivative slip law w.r.t. slip_velocity component 2
    stack.jacobian[2][2] = dStateEvolutionLaw[2];
  }

  GEOS_HOST_DEVICE
  void solve( localIndex const k,
              StackVariables & stack ) const
  {
    /// Solve 3x3 system
    real64 solution[3] = {0.0, 0.0, 0.0};

    denseLinearAlgebra::solve< 3 >( stack.jacobian, stack.rhs, solution );

    /// Update variables
    m_stateVariable[k] -= solution[0];
    m_slipRate[k][0]   -= solution[1];
    m_slipRate[k][1]   -= solution[2];

    // // Matteo: debugging tools.
    // printf("solution[0] = %.10e\n", solution[0]);
    // printf("solution[1] = %.10e\n", solution[1]);
    // printf("solution[2] = %.10e\n", solution[2]);

    // printf("m_stateVariable[%d] = %.10e\n", k, m_stateVariable[k]);
    // printf("m_slipRate[0][%d] = %.10e\n", k, m_slipRate[k][0]);
    // printf("m_slipRate[1][%d] = %.10e\n", k, m_slipRate[k][1]);
  }

  GEOS_HOST_DEVICE
  std::pair< int, real64 > checkConvergence( StackVariables const & stack,
                                             real64 const tol ) const
  {
    real64 const residualNorm = LvArray::tensorOps::l2Norm< 3 >( stack.rhs );
    int const converged = residualNorm < tol ? 1 : 0;
    return std::make_pair( converged, residualNorm );
  }

private:

  arrayView2d< real64 > const m_slipRate;

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

  // Newton loop (outside of the kernel launch)
  bool allConverged = false;
  for( integer iter = 0; iter < maxNewtonIter; iter++ )
  {
    RAJA::ReduceMin< parallelDeviceReduce, int > converged( 1 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > residualNorm( 0.0 );
    forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      RateAndStateKernel::StackVariables stack;
      kernel.setup( k, dt, stack );
      kernel.solve( k, stack );
      auto result = kernel.checkConvergence( stack, 1.0e-6 );
      converged.min( std::get< 0 >( result ) );
      residualNorm.max( std::get< 1 >( result ) );
    } );

    real64 const maxResidualNorm = MpiWrapper::max( residualNorm.get() );
    GEOS_LOG_RANK_0( GEOS_FMT( "-----iter {} : residual = {:.10e} ", iter, maxResidualNorm ) );

    if( converged.get() )
    {
      allConverged = true;
      break;
    }
  }
  if( !allConverged )
  {
    GEOS_ERROR( " Failed to converge" );
  }
}

} /* namespace rateAndStateKernels */

}/* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_RATEANDSTATEKERNELS_HPP_ */
