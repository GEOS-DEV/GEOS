/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DirichletFluxComputeZFormulationKernelFactory.cpp
 */

#include "DirichletFluxComputeZFormulationKernel.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernels
{

template< typename POLICY >
void
DirichletFluxComputeZFormulationKernelFactory< POLICY >::
createAndLaunch( integer const numComps,
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

    using kernelType = DirichletFluxComputeZFormulationKernel< typename FluidType::KernelWrapper, POLICY, NUM_COMP, NUM_DOF >;
    typename kernelType::CompFlowAccessors compFlowAccessors( elemManager, solverName );
    typename kernelType::MultiFluidAccessors multiFluidAccessors( elemManager, solverName );
    typename kernelType::CapPressureAccessors capPressureAccessors( elemManager, solverName );
    typename kernelType::PermeabilityAccessors permeabilityAccessors( elemManager, solverName );

    kernelType kernel( numPhases, rankOffset, faceManager, stencilWrapper, fluidWrapper,
                       dofNumberAccessor, compFlowAccessors, multiFluidAccessors, capPressureAccessors, permeabilityAccessors,
                       dt, localMatrix, localRhs, kernelFlags );
    kernel.launchKernel( stencilWrapper.size() );
  } );
}

template class DirichletFluxComputeZFormulationKernelFactory< parallelDevicePolicy<> >;

} // namespace isothermalCompositionalMultiphaseFVMKernels

} // namespace geos
