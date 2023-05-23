/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ThermalSinglePhasePoromechanicsEFEM_impl.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/contact/ContactFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsEFEM.hpp"
#include "physicsSolvers/contact/SolidMechanicsEFEMKernelsHelper.hpp"

namespace geos
{

namespace thermalSinglePhasePoromechanicsEmbeddedFracturesKernels
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
ThermalSinglePhasePoromechanicsEFEM( NodeManager const & nodeManager,
                                     EdgeManager const & edgeManager,
                                     FaceManager const & faceManager,
                                     localIndex const targetRegionIndex,
                                     SUBREGION_TYPE const & elementSubRegion,
                                     FE_TYPE const & finiteElementSpace,
                                     CONSTITUTIVE_TYPE & inputConstitutiveType,
                                     EmbeddedSurfaceSubRegion const & embeddedSurfSubRegion,
                                     arrayView1d< globalIndex const > const dispDofNumber,
                                     arrayView1d< globalIndex const > const jumpDofNumber,
                                     string const inputFlowDofKey,
                                     globalIndex const rankOffset,
                                     CRSMatrixView< real64, globalIndex const > const inputMatrix,
                                     arrayView1d< real64 > const inputRhs,
                                     real64 const (&inputGravityVector)[3],
                                     string const fluidModelKey ):
  Base( nodeManager,
        edgeManager,
        faceManager,
        targetRegionIndex,
        elementSubRegion,
        finiteElementSpace,
        inputConstitutiveType,
        embeddedSurfSubRegion,
        dispDofNumber,
        jumpDofNumber,
        inputFlowDofKey,
        rankOffset,
        inputMatrix,
        inputRhs,
        inputGravityVector,
        fluidModelKey )
{}


//START_kernelLauncher
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
template< typename POLICY,
          typename KERNEL_TYPE >
real64
ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
kernelLaunch( localIndex const numElems,
              KERNEL_TYPE const & kernelComponent )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( numElems );

  // Define a RAJA reduction variable to get the maximum residual contribution.
  RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

  forAll< POLICY >( kernelComponent.m_fracturedElems.size(),
                    [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localIndex k = kernelComponent.m_fracturedElems[i];
    typename KERNEL_TYPE::StackVariables stack;

    kernelComponent.setup( k, stack );
    for( integer q=0; q<numQuadraturePointsPerElem; ++q )
    {
      kernelComponent.quadraturePointKernel( k, q, stack );
    }
    maxResidual.max( kernelComponent.complete( k, stack ) );
  } );

  return maxResidual.get();
}
//END_kernelLauncher


template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
setup( localIndex const k,
       StackVariables & stack ) const
{
  Base::setup(k, stack);
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
quadraturePointKernel( localIndex const k,
                       localIndex const q,
                       StackVariables & stack ) const
{

  Base::quadraturePointKernel(k, q, stack);

  // assemble KwTmLocal
}

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
real64 ThermalSinglePhasePoromechanicsEFEM< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::
complete( localIndex const k,
          StackVariables & stack ) const
{
  real64 const maxForce = Base::complete( k, stack );

  // add thermal derivatives

  return maxForce;
}


} // namespace thermalSinglePhasePoromechanicsEmbeddedFracturesKernels

} /* namespace geos */

#endif // GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSEFEM_IMPL_HPP_
