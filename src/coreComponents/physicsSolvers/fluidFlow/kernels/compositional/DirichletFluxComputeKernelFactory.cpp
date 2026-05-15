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
 * @file DirichletFluxComputeKernelFactory.cpp
 */

#include "DirichletFluxComputeKernel.hpp"
#include "ThermalDirichletFluxComputeKernel.hpp"

#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"
#include "finiteVolume/BoundaryStencil.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

template< typename POLICY >
void
DirichletFluxComputeKernelFactory< POLICY >::createAndLaunch( integer const numComps,
                                                              integer const numPhases,
                                                              globalIndex const rankOffset,
                                                              BitFlags< KernelFlags > kernelFlags,
                                                              string const & dofKey,
                                                              string const & solverName,
                                                              FaceManager const & faceManager,
                                                              ElementRegionManager const & elemManager,
                                                              BoundaryStencilWrapper const & stencilWrapper,
                                                              constitutive::MultiFluidBase & fluidBase,
                                                              real64 const dt,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  constitutive::constitutiveComponentUpdatePassThru( fluidBase, numComps, [&]( auto & fluid, auto NC )
  {
    using FluidType = TYPEOFREF( fluid );
    typename FluidType::KernelWrapper const fluidWrapper = fluid.createKernelWrapper();

    integer constexpr NUM_COMP = NC();
    integer constexpr NUM_DOF = NC() + 1;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using kernelType = DirichletFluxComputeKernel< typename FluidType::KernelWrapper, NUM_COMP, NUM_DOF >;
    typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
    typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

    kernelType kernel( numPhases, rankOffset, faceManager, stencilWrapper, fluidWrapper,
                       dofNumberAccessor, compFlowAccessors, multiFluidAccessors, capPressureAccessors, permeabilityAccessors,
                       dt, localMatrix, localRhs, kernelFlags );
    kernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  } );
}

} // namespace isothermalCompositionalMultiphaseFVMKernels

namespace thermalCompositionalMultiphaseFVMKernels
{
template< typename POLICY, typename STENCILWRAPPER >
void
DirichletFluxComputeKernelFactory< POLICY, STENCILWRAPPER >::
createAndLaunch( integer const numComps,
                 integer const numPhases,
                 globalIndex const rankOffset,
                 BitFlags< isothermalCompositionalMultiphaseFVMKernels::KernelFlags > kernelFlags,
                 string const & dofKey,
                 string const & solverName,
                 FaceManager const & faceManager,
                 ElementRegionManager const & elemManager,
                 STENCILWRAPPER const & stencilWrapper,
                 constitutive::MultiFluidBase & fluidBase,
                 real64 const dt,
                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                 arrayView1d< real64 > const & localRhs )
{
  constitutive::constitutiveComponentUpdatePassThru< true >( fluidBase, numComps, [&]( auto & fluid, auto NC )
  {
    using FluidType = TYPEOFREF( fluid );
    typename FluidType::KernelWrapper const fluidWrapper = fluid.createKernelWrapper();

    integer constexpr NUM_COMP = NC();
    integer constexpr NUM_DOF = NC() + 2;

    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > dofNumberAccessor =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( dofKey );
    dofNumberAccessor.setName( solverName + "/accessors/" + dofKey );

    using KernelType = DirichletFluxComputeKernel< typename FluidType::KernelWrapper, NUM_COMP, NUM_DOF >;
    typename KernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
    typename KernelType::ThermalCompFlowAccessors thermalCompFlowAccessors( elemManager, solverName );
    typename KernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
    typename KernelType::ThermalMultiFluidAccessors thermalMultiFluidAccessors( elemManager, solverName );
    typename KernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
    typename KernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );
    typename KernelType::ThermalConductivityAccessors thermalConductivityAccessors( elemManager, solverName );

    KernelType kernel( numPhases, rankOffset, faceManager, stencilWrapper, fluidWrapper,
                       dofNumberAccessor, compFlowAccessors, thermalCompFlowAccessors, multiFluidAccessors, thermalMultiFluidAccessors,
                       capPressureAccessors, permeabilityAccessors, thermalConductivityAccessors,
                       dt, localMatrix, localRhs, kernelFlags );
    KernelType::template launch< POLICY >( stencilWrapper.size(), kernel );
  } );
}

} // namespace thermalCompositionalMultiphaseFVMKernels

template class isothermalCompositionalMultiphaseFVMKernels::DirichletFluxComputeKernelFactory< parallelDevicePolicy<> >;
template class isothermalCompositionalMultiphaseFVMKernels::DirichletFluxComputeKernelFactory< serialPolicy >;
template class thermalCompositionalMultiphaseFVMKernels::DirichletFluxComputeKernelFactory< parallelDevicePolicy<>, BoundaryStencilWrapper >;
template class thermalCompositionalMultiphaseFVMKernels::DirichletFluxComputeKernelFactory< serialPolicy, BoundaryStencilWrapper >;

} // namespace geos
