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
 * @file CompositionalMultiphaseHybridFVMKernelsInstantiations.cpp
 * @brief Centralized explicit template instantiations for compositional hybrid FVM kernels.
 *
 * Define GEOS_ENABLE_MANUAL_HYBRID_FVM_INST to activate the explicit instantiations below.
 * When enabled, ensure overlapping template-generated sources are not compiled to avoid duplicates.
 */

#include "physicsSolvers/fluidFlow/kernels/compositional/CompositionalMultiphaseHybridFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/CompositionalMultiphaseHybridFVMKernels_impl.hpp"

#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"

namespace geos
{
namespace compositionalMultiphaseHybridFVMKernels
{

#if defined(GEOS_ENABLE_MANUAL_HYBRID_FVM_INST)

// -----------------------------------------------------------------------------
// Range driver
// -----------------------------------------------------------------------------
#define GEOS_NFACES_LIST \
  X_NF( 4 )                \
  X_NF( 5 )                \
  X_NF( 6 )                \
  X_NF( 7 )                \
  X_NF( 8 )                \
  X_NF( 9 )                \
  X_NF( 10 )               \
  X_NF( 11 )               \
  X_NF( 12 )               \
  X_NF( 13 )

// -----------------------------------------------------------------------------
// Common parameter macros (kept in sync with kernel signatures)
// -----------------------------------------------------------------------------
#define GEOS_DIRICHLET_LAUNCH_PARAMS \
  integer const numPhases, \
  localIndex const er, \
  localIndex const esr, \
  CellElementSubRegion const & subRegion, \
  constitutive::PermeabilityBase const & permeabilityModel, \
  SortedArrayView< localIndex const > const & boundaryFaceSet, \
  arrayView1d< real64 const > const & facePres, \
  arrayView1d< real64 const > const & faceTemp, \
  arrayView2d< real64 const, compflow::USD_PHASE > const & facePhaseMob, \
  arrayView2d< real64 const, compflow::USD_PHASE > const & facePhaseMassDens, \
  arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & facePhaseCompFrac, \
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
  ArrayOfArraysView< localIndex const > const & faceToNodes, \
  arrayView2d< localIndex const > const & elemRegionList, \
  arrayView2d< localIndex const > const & elemSubRegionList, \
  arrayView2d< localIndex const > const & elemList, \
  arrayView1d< globalIndex const > const & faceDofNumber, \
  arrayView1d< real64 const > const & faceGravCoef, \
  arrayView1d< real64 const > const & transMultiplier, \
  real64 const lengthTolerance, \
  real64 const dt, \
  globalIndex const rankOffset, \
  integer const useTotalMassEquation, \
  DirichletFluxKernel::CompFlowAccessors const & compFlowAccessors, \
  DirichletFluxKernel::MultiFluidAccessors const & multiFluidAccessors, \
  DirichletFluxKernel::ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
  DefaultGlobalMatrixView const & localMatrix, \
  arrayView1d< real64 > const & localRhs

#define GEOS_FLUX_LAUNCH_PARAMS \
  localIndex er, localIndex esr, \
  CellElementSubRegion const & subRegion, \
  constitutive::PermeabilityBase const & permeabilityModel, \
  SortedArrayView< localIndex const > const & regionFilter, \
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
  arrayView2d< localIndex const > const & elemRegionList, \
  arrayView2d< localIndex const > const & elemSubRegionList, \
  arrayView2d< localIndex const > const & elemList, \
  ArrayOfArraysView< localIndex const > const & faceToNodes, \
  arrayView1d< globalIndex const > const & faceDofNumber, \
  arrayView1d< integer const > const & faceGhostRank, \
  arrayView1d< integer const > const & isBoundaryFace, \
  arrayView1d< real64 const > const & facePres, \
  arrayView1d< real64 const > const & faceGravCoef, \
  arrayView1d< real64 const > const & mimFaceGravCoef, \
  arrayView1d< real64 const > const & transMultiplier, \
  FluxKernel::ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
  FluxKernel::ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob, \
  FluxKernel::ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
  FluxKernel::ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
  FluxKernel::ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens, \
  FluxKernel::ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
  FluxKernel::ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens, \
  FluxKernel::ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
  FluxKernel::ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac, \
  FluxKernel::ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
  globalIndex const rankOffset, \
  real64 const lengthTolerance, \
  real64 const dt, \
  integer const useTotalMassEquation, \
  DefaultGlobalMatrixView const & localMatrix, \
  arrayView1d< real64 > const & localRhs

// -----------------------------------------------------------------------------
// IP list helper
// -----------------------------------------------------------------------------
#define GEOS_FOR_EACH_IP( NF, NC, NP, MACRO ) \
  MACRO( NF, NC, NP, mimeticInnerProduct::TPFAInnerProduct ) \
  MACRO( NF, NC, NP, mimeticInnerProduct::QuasiTPFAInnerProduct ) \
  MACRO( NF, NC, NP, mimeticInnerProduct::BdVLMInnerProduct )

// -----------------------------------------------------------------------------
// DirichletFluxKernel explicit instantiations
// -----------------------------------------------------------------------------
#define INSTANTIATE_DIRICHLET( NF, NC, NP, IP ) \
  template void DirichletFluxKernel::launch< NF, NC, NP, IP >( GEOS_DIRICHLET_LAUNCH_PARAMS );

#define INSTANTIATE_DIRICHLET_FOR_NP( NF, NC, NP ) \
  GEOS_FOR_EACH_IP( NF, NC, NP, INSTANTIATE_DIRICHLET )

#define INSTANTIATE_DIRICHLET_FOR_NC( NF, NC ) \
  INSTANTIATE_DIRICHLET_FOR_NP( NF, NC, 2 ) \
  INSTANTIATE_DIRICHLET_FOR_NP( NF, NC, 3 )

#define INSTANTIATE_DIRICHLET_FOR_NF( NF ) \
  INSTANTIATE_DIRICHLET_FOR_NC( NF, 1 ) \
  INSTANTIATE_DIRICHLET_FOR_NC( NF, 2 ) \
  INSTANTIATE_DIRICHLET_FOR_NC( NF, 3 ) \
  INSTANTIATE_DIRICHLET_FOR_NC( NF, 4 ) \
  INSTANTIATE_DIRICHLET_FOR_NC( NF, 5 )

#define X_NF( NF ) INSTANTIATE_DIRICHLET_FOR_NF( NF )
GEOS_NFACES_LIST
#undef X_NF

// -----------------------------------------------------------------------------
// FluxKernel explicit instantiations
// -----------------------------------------------------------------------------
#define INSTANTIATE_FLUX( NF, NC, NP, IP ) \
  template void FluxKernel::launch< NF, NC, NP, IP >( GEOS_FLUX_LAUNCH_PARAMS );

#define INSTANTIATE_FLUX_FOR_NP( NF, NC, NP ) \
  GEOS_FOR_EACH_IP( NF, NC, NP, INSTANTIATE_FLUX )

#define INSTANTIATE_FLUX_FOR_NC( NF, NC ) \
  INSTANTIATE_FLUX_FOR_NP( NF, NC, 2 ) \
  INSTANTIATE_FLUX_FOR_NP( NF, NC, 3 )

#define INSTANTIATE_FLUX_FOR_NF( NF ) \
  INSTANTIATE_FLUX_FOR_NC( NF, 1 ) \
  INSTANTIATE_FLUX_FOR_NC( NF, 2 ) \
  INSTANTIATE_FLUX_FOR_NC( NF, 3 ) \
  INSTANTIATE_FLUX_FOR_NC( NF, 4 ) \
  INSTANTIATE_FLUX_FOR_NC( NF, 5 )

#define X_NF( NF ) INSTANTIATE_FLUX_FOR_NF( NF )
GEOS_NFACES_LIST
#undef X_NF

// Cleanup local helper macros
#undef GEOS_DIRICHLET_LAUNCH_PARAMS
#undef GEOS_FLUX_LAUNCH_PARAMS
#undef GEOS_FOR_EACH_IP
#undef INSTANTIATE_DIRICHLET
#undef INSTANTIATE_DIRICHLET_FOR_NP
#undef INSTANTIATE_DIRICHLET_FOR_NC
#undef INSTANTIATE_DIRICHLET_FOR_NF
#undef INSTANTIATE_FLUX
#undef INSTANTIATE_FLUX_FOR_NP
#undef INSTANTIATE_FLUX_FOR_NC
#undef INSTANTIATE_FLUX_FOR_NF

// -----------------------------------------------------------------------------
// evaluateBCFaceProperties explicit instantiations
// -----------------------------------------------------------------------------
#define INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES( NC, NP ) \
  template void evaluateBCFaceProperties< NC, NP >( \
    integer const numPhases, \
    SortedArrayView< localIndex const > const & boundaryFaceSet, \
    arrayView1d< real64 const > const & facePres, \
    arrayView1d< real64 const > const & faceTemp, \
    arrayView2d< real64 const, compflow::USD_COMP > const & faceCompFrac, \
    arrayView2d< localIndex const > const & elemRegionList, \
    arrayView2d< localIndex const > const & elemSubRegionList, \
    arrayView2d< localIndex const > const & elemList, \
    localIndex const er, \
    localIndex const esr, \
    constitutive::MultiFluidBase & fluid, \
    constitutive::RelativePermeabilityBase & relperm, \
    arrayView2d< real64, compflow::USD_PHASE > const & facePhaseMob, \
    arrayView2d< real64, compflow::USD_PHASE > const & facePhaseMassDens, \
    arrayView3d< real64, compflow::USD_PHASE_COMP > const & facePhaseCompFrac );

#define INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NP( NC, NP ) \
  INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES( NC, NP )

#define INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NC( NC ) \
  INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NP( NC, 2 ) \
  INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NP( NC, 3 )

INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NC( 1 )
INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NC( 2 )
INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NC( 3 )
INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NC( 4 )
INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NC( 5 )

// Cleanup local helper macros
#undef INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES
#undef INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NP
#undef INSTANTIATE_EVALUATE_BC_FACE_PROPERTIES_FOR_NC

// -----------------------------------------------------------------------------
// NOTE: Explicit instantiations for AssemblerKernelHelper and AssemblerKernel
// functions are NOT included here because they are declared as 'inline' in
// CompositionalMultiphaseHybridFVMKernels_impl.hpp. Inline template functions
// are automatically instantiated at the point of use and cannot be explicitly
// instantiated.
//
// The following functions are inline and therefore auto-instantiated:
// - AssemblerKernelHelper::applyGradient
// - AssemblerKernelHelper::assembleFluxDivergence
// - AssemblerKernelHelper::assembleViscousFlux
// - AssemblerKernelHelper::assembleBuoyancyFlux
// - AssemblerKernelHelper::assembleFaceConstraints
// - AssemblerKernel::compute
// -----------------------------------------------------------------------------

#endif // GEOS_ENABLE_MANUAL_HYBRID_FVM_INST

} // namespace compositionalMultiphaseHybridFVMKernels
} // namespace geos
