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
 * @file ReactionUpdateKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_REACTIONUPDATEKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_REACTIONUPDATEKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "constitutive/fluid/reactivefluid/ReactiveFluidLayouts.hpp"

namespace geos
{

namespace singlePhaseReactiveBaseKernels
{

/******************************** EquilibriumReactionUpdateKernel ********************************/

struct EquilibriumReactionUpdateKernel
{

  /**
   * @return whether the equilibrium solve converged in every cell
   */
  template< typename REACTION_WRAPPER_TYPE >
  static bool helper( REACTION_WRAPPER_TYPE const & reactionWrapper,
                      arrayView1d< real64 const > const & pres,
                      arrayView1d< real64 const > const & temp,
                      arrayView2d< real64, compflow::USD_COMP > const logPrimaryConc )
  {
    RAJA::ReduceMax< parallelDeviceReduce, integer > anyNotConverged( 0 );

    forAll< parallelDevicePolicy<> >( reactionWrapper.numElems(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      if( !reactionWrapper.updateEquilibriumReaction( k, pres[k], temp[k], logPrimaryConc[k] ) )
      {
        anyNotConverged.max( 1 );
      }
    } );

    return anyNotConverged.get() == 0;
  }

  /**
   * @return whether the equilibrium solve converged in every cell, for the caller to report
   */
  template< typename REACTIVE_FLUID >
  static bool launch( REACTIVE_FLUID const & fluid,
                      arrayView1d< real64 const > const & pres,
                      arrayView1d< real64 const > const & temp,
                      arrayView2d< real64, compflow::USD_COMP > const logPrimaryConc )
  {
    bool converged = true;

    std::visit( [&]( auto const & reactionWrapper )
    {
      converged = helper( reactionWrapper, pres, temp, logPrimaryConc );
    }, fluid.createReactionKernelWrapper());

    return converged;
  }
};

/******************************** MixedSystemReactionUpdateKernel ********************************/

struct MixedSystemReactionUpdateKernel
{

  template< typename REACTION_WRAPPER_TYPE >
  static void helper( REACTION_WRAPPER_TYPE const & reactionWrapper,
                      arrayView1d< real64 const > const & pres,
                      arrayView1d< real64 const > const & temp,
                      arrayView2d< real64 const, compflow::USD_COMP > const logPrimaryConc,
                      arrayView2d< real64 const, compflow::USD_COMP > const surfaceArea )
  {
    forAll< parallelDevicePolicy<> >( reactionWrapper.numElems(), [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      reactionWrapper.updateMixedReactionSystem( k, pres[k], temp[k], logPrimaryConc[k], surfaceArea[k] );
    } );
  }

  template< typename REACTIVE_FLUID >
  static void launch( REACTIVE_FLUID const & fluid,
                      arrayView1d< real64 const > const & pres,
                      arrayView1d< real64 const > const & temp,
                      arrayView2d< real64 const, compflow::USD_COMP > const logPrimaryConc,
                      arrayView2d< real64 const, compflow::USD_COMP > const surfaceArea )
  {
    std::visit( [&]( auto const & reactionWrapper )
    {
      helper( reactionWrapper, pres, temp, logPrimaryConc, surfaceArea );
    }, fluid.createReactionKernelWrapper());
  }
};

} // namespace singlePhaseReactiveBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_FLUIDUPDATEKERNEL_HPP
