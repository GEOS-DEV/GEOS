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
 * @file SinglePhaseHybridFVMKernelsInstantiations.cpp
 * @brief Centralized explicit template instantiations for single-phase hybrid FVM kernels.
 *
 * Define GEOS_ENABLE_MANUAL_HYBRID_FVM_INST to activate the explicit instantiations below.
 * When enabled, ensure overlapping template-generated sources are not compiled to avoid duplicates.
 */

#include "physicsSolvers/fluidFlow/kernels/singlePhase/SinglePhaseHybridFVMKernels.hpp"

#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiRTInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"

namespace geos
{
namespace singlePhaseHybridFVMKernels
{

#if defined(GEOS_ENABLE_MANUAL_HYBRID_FVM_INST)

// -----------------------------------------------------------------------------
// Faces-per-element list (kept in sync with kernelLaunchSelectorFaceSwitch)
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
// Inner product list helper
// -----------------------------------------------------------------------------
#define GEOS_FOR_EACH_IP( NF, MACRO ) \
  MACRO( NF, mimeticInnerProduct::TPFAInnerProduct ) \
  MACRO( NF, mimeticInnerProduct::QuasiTPFAInnerProduct ) \
  MACRO( NF, mimeticInnerProduct::QuasiRTInnerProduct ) \
  MACRO( NF, mimeticInnerProduct::SimpleInnerProduct ) \
  MACRO( NF, mimeticInnerProduct::BdVLMInnerProduct )

// -----------------------------------------------------------------------------
// AveragePressureGradientKernel explicit class instantiations
// -----------------------------------------------------------------------------
#define INSTANTIATE_AVG_PRES_GRAD( NF ) \
  template class AveragePressureGradientKernel< NF >;

#define X_NF( NF ) INSTANTIATE_AVG_PRES_GRAD( NF )
GEOS_NFACES_LIST
#undef X_NF
#undef INSTANTIATE_AVG_PRES_GRAD

// -----------------------------------------------------------------------------
// ElementBasedAssemblyKernel explicit class instantiations
// -----------------------------------------------------------------------------
#define INSTANTIATE_ELEM_KERNEL( NF, IP ) \
  template class ElementBasedAssemblyKernel< NF, IP >;

#define INSTANTIATE_ELEM_KERNEL_FOR_IPS( NF ) \
  GEOS_FOR_EACH_IP( NF, INSTANTIATE_ELEM_KERNEL )

#define X_NF( NF ) INSTANTIATE_ELEM_KERNEL_FOR_IPS( NF )
GEOS_NFACES_LIST
#undef X_NF
#undef INSTANTIATE_ELEM_KERNEL_FOR_IPS
#undef INSTANTIATE_ELEM_KERNEL

// Cleanup local helper macros
#undef GEOS_FOR_EACH_IP
#undef GEOS_NFACES_LIST

#endif // GEOS_ENABLE_MANUAL_HYBRID_FVM_INST

} // namespace singlePhaseHybridFVMKernels
} // namespace geos
