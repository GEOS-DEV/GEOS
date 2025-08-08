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

#ifndef GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_RATEANDSTATEKERNELSBASE_HPP_
#define GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_RATEANDSTATEKERNELSBASE_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include "physicsSolvers/inducedSeismicity/rateAndStateFields.hpp"
#include "constitutive/ConstitutivePassThru.hpp"

namespace geos
{

namespace rateAndStateKernels
{

// TBD: Pass the kernel and add getters for relevant fields to make this function general purpose and avoid
// wrappers?
GEOS_HOST_DEVICE
static void projectSlipRateBase( localIndex const k,
                                 real64 const frictionCoefficient,
                                 real64 const shearImpedance,
                                 arrayView1d< real64 const > const normalTraction,
                                 arrayView2d< real64 const > const shearTraction,
                                 arrayView1d< real64 const > const slipRate,
                                 arrayView2d< real64 > const slipVelocity )
{
  // Project slip rate onto shear traction to get slip velocity components
  real64 const frictionForce = normalTraction[k] * frictionCoefficient;
  real64 const projectionScaling = 1.0 / ( shearImpedance +  frictionForce / slipRate[k] );
  slipVelocity[k][0] = projectionScaling * shearTraction[k][0];
  slipVelocity[k][1] = projectionScaling * shearTraction[k][1];
}

template< typename POLICY, typename KERNEL_TYPE >
static bool newtonSolve( SurfaceElementSubRegion & subRegion,
                         KERNEL_TYPE & kernel,
                         real64 const dt,
                         integer const maxNewtonIter,
                         real64 const newtonTol )
{
  bool allConverged = false;
  for( integer iter = 0; iter < maxNewtonIter; iter++ )
  {
    RAJA::ReduceMin< parallelDeviceReduce, int > converged( 1 );
    RAJA::ReduceMax< parallelDeviceReduce, real64 > residualNorm( 0.0 );
    forAll< POLICY >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename KERNEL_TYPE::StackVariables stack;
      kernel.setup( k, dt, stack );
      kernel.solve( k, stack );
      auto const [elementConverged, elementResidualNorm] = kernel.checkConvergence( stack, newtonTol );
      converged.min( elementConverged );
      residualNorm.max( elementResidualNorm );
    } );

    real64 const maxResidualNorm = MpiWrapper::max( residualNorm.get() );
    GEOS_LOG_RANK_0( GEOS_FMT( "   Newton iter {} : residual = {:.10e} ", iter, maxResidualNorm ) );

    if( converged.get() )
    {
      allConverged = true;
      break;
    }
  }
  return allConverged;
}
template< typename FRICTION_TYPE >
void enforceRateAndVelocityConsistency( FRICTION_TYPE const & frictionLawKernelWrapper,
                                        SurfaceElementSubRegion & subRegion,
                                        real64 const & shearImpedance )
{
  arrayView2d< real64 > const slipVelocity = subRegion.getField< fields::rateAndState::slipVelocity >();
  arrayView1d< real64 > const slipRate  = subRegion.getField< fields::rateAndState::slipRate >();
  arrayView1d< real64 const > const stateVariable  = subRegion.getField< fields::rateAndState::stateVariable >();

  arrayView2d< real64 > const backgroundShearStress = subRegion.getField< fields::rateAndState::backgroundShearStress >();
  arrayView1d< real64 > const backgroundNormalStress = subRegion.getField< fields::rateAndState::backgroundNormalStress >();

  RAJA::ReduceMax< parallelDeviceReduce, int > negativeSlipRate( 0 );
  RAJA::ReduceMax< parallelDeviceReduce, int > bothNonZero( 0 );

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    if( slipRate[k] < 0.0 )
    {
      negativeSlipRate.max( 1 );
    }
    else if( LvArray::tensorOps::l2Norm< 2 >( slipVelocity[k] ) > 0.0 && slipRate[k] > 0.0 )
    {
      bothNonZero.max( 1 );
    }
    else if( LvArray::tensorOps::l2Norm< 2 >( slipVelocity[k] ) > 0.0 )
    {
      slipRate[k] = LvArray::tensorOps::l2Norm< 2 >( slipVelocity[k] );
    }
    else if( slipRate[k] > 0.0 )
    {
      real64 const frictionCoefficient = frictionLawKernelWrapper.frictionCoefficient( k, slipRate[k], stateVariable[k] );
      projectSlipRateBase( k,
                           frictionCoefficient,
                           shearImpedance,
                           backgroundNormalStress,
                           backgroundShearStress,
                           slipRate,
                           slipVelocity );
    }
  } );

  GEOS_ERROR_IF( negativeSlipRate.get() > 0, "SlipRate cannot be negative." );
  GEOS_ERROR_IF( bothNonZero.get() > 0, "Only one between slipRate and slipVelocity can be specified as i.c." );
}

/**
 * @brief Performs the kernel launch
 * @tparam POLICY the policy used in the RAJA kernels
 */
template< template< typename FRICTION_TYPE > class KERNEL_TYPE,
          typename POLICY,
          typename FRICTION_TYPE >
static void
createAndLaunch( SurfaceElementSubRegion & subRegion,
                 FRICTION_TYPE & frictionLaw,
                 real64 const shearImpedance,
                 integer const maxNewtonIter,
                 real64 const newtonTol,
                 real64 const time_n,
                 real64 const totalDt )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n );

  KERNEL_TYPE kernel( subRegion, frictionLaw, shearImpedance );

  real64 dtRemaining = totalDt;
  real64 dt = totalDt;
  for( integer subStep = 0; subStep < 5 && dtRemaining > 0.0; ++subStep )
  {
    real64 dtAccepted = KERNEL_TYPE< FRICTION_TYPE >::template solveRateAndStateEquation< POLICY >( subRegion, kernel, dt, maxNewtonIter, newtonTol );
    dtRemaining -= dtAccepted;

    if( dtRemaining > 0.0 )
    {
      dt = dtAccepted;
    }
    GEOS_LOG_RANK_0( GEOS_FMT( "  sub-step = {} completed, dt = {}, remaining dt = {}", subStep, dt, dtRemaining ) );
  }
}

} /* namespace rateAndStateKernels */

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_INDUCEDSEISMICITY_KERNELS_RATEANDSTATEKERNELSBASE_HPP_ */
