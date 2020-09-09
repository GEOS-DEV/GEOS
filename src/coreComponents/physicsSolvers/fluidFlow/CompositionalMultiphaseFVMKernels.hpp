/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFVMKernels
{

/******************************** PhaseMobilityKernel ********************************/

/**
 * @brief Functions to compute phase mobilities and derivatives from density, viscosity and relperm
 */
struct PhaseMobilityKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( arraySlice2d< real64 const > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const > const & phaseDens,
           arraySlice1d< real64 const > const & dPhaseDens_dPres,
           arraySlice2d< real64 const > const & dPhaseDens_dComp,
           arraySlice1d< real64 const > const & phaseVisc,
           arraySlice1d< real64 const > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const > const & phaseRelPerm,
           arraySlice2d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64 > const & phaseMob,
           arraySlice1d< real64 > const & dPhaseMob_dPres,
           arraySlice2d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  Launch( localIndex const size,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  Launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp );
};


/******************************** FluxKernel ********************************/

/**
 * @brief Functions to assemble flux term contributions to residual and Jacobian
 */
struct FluxKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( localIndex const stencilSize,
           localIndex const numPhases,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           arraySlice1d< real64 const > const & stencilWeights,
           ElementView< arrayView1d< real64 const > > const & pres,
           ElementView< arrayView1d< real64 const > > const & dPres,
           ElementView< arrayView1d< real64 const > > const & gravCoef,
           ElementView< arrayView2d< real64 const > > const & phaseMob,
           ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,
           ElementView< arrayView3d< real64 const > > const & dPhaseMob_dComp,
           ElementView< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
           ElementView< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
           ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
           ElementView< arrayView3d< real64 const > > const & phaseDens,
           ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,
           ElementView< arrayView4d< real64 const > > const & dPhaseDens_dComp,
           ElementView< arrayView4d< real64 const > > const & phaseCompFrac,
           ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
           ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
           ElementView< arrayView3d< real64 const > > const & phaseCapPressure,
           ElementView< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
           integer const capPressureFlag,
           real64 const dt,
           arraySlice1d< real64 > const & localFlux,
           arraySlice2d< real64 > const & localFluxJacobian );

  template< localIndex NC, typename STENCIL_TYPE >
  static void
  Launch( localIndex const numPhases,
          STENCIL_TYPE const & stencil,
          globalIndex const rankOffset,
          ElementView< arrayView1d< globalIndex const > > const & dofNumber,
          ElementView< arrayView1d< integer const > > const & ghostRank,
          ElementView< arrayView1d< real64 const > > const & pres,
          ElementView< arrayView1d< real64 const > > const & dPres,
          ElementView< arrayView1d< real64 const > > const & gravCoef,
          ElementView< arrayView2d< real64 const > > const & phaseMob,
          ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,
          ElementView< arrayView3d< real64 const > > const & dPhaseMob_dComp,
          ElementView< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
          ElementView< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
          ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
          ElementView< arrayView3d< real64 const > > const & phaseDens,
          ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,
          ElementView< arrayView4d< real64 const > > const & dPhaseDens_dComp,
          ElementView< arrayView4d< real64 const > > const & phaseCompFrac,
          ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
          ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
          ElementView< arrayView3d< real64 const > > const & phaseCapPressure,
          ElementView< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
          integer const capPressureFlag,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
};

} // namespace CompositionalMultiphaseFVMKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP
