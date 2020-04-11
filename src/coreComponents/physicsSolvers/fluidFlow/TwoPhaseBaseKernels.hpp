/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TwoPhaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/TwoPhaseBase.hpp"

namespace geosx
{

namespace TwoPhaseBaseKernels
{

/******************************** PhaseMobilityKernel ********************************/

struct PhaseMobilityKernel
{

  /**
   * @brief In an element, compute the phase mobilities
   * @param[in] phaseDens phase densities in the element
   * @param[in] dPhaseDens_dPres derivatives of phase densities wrt pressure in the element
   * @param[in] phaseVisc phase viscosities in the element
   * @param[in] dPhaseVisc_dPres derivatives of phase viscosities wrt pressure in the element
   * @param[in] phaseRelPerm phase relative permeabilities in the element
   * @param[in] dPhaseRelPerm_dSat derivatives of the phase relative permeabilities wrt primary saturation in the
   * element
   * @param[inout] phaseMob phase mobilities in the element
   * @param[inout] dPhaseMob_dPres derivatives of phase mobilities wrt pressure in the element
   * @param[inout] dPhaseMob_dSat derivatives of phase mobilities wrt primary saturation in the element
   */
  static void
  Compute( arraySlice1d< real64 const > phaseDens,
           arraySlice1d< real64 const > dPhaseDens_dPres,
           arraySlice1d< real64 const > phaseVisc,
           arraySlice1d< real64 const > dPhaseVisc_dPres,
           arraySlice1d< real64 const > phaseRelPerm,
           arraySlice2d< real64 const > dPhaseRelPerm_dSat,
           arraySlice1d< real64 > phaseMob,
           arraySlice1d< real64 > dPhaseMob_dPres,
           arraySlice1d< real64 > dPhaseMob_dSat );

  /**
   * @brief In a sub region, compute the phase mobilities
   * @param[in] number of elements in the sub region
   * @param[in] phaseDens phase densities in the sub region
   * @param[in] dPhaseDens_dPres derivatives of phase densities wrt pressure in the sub region
   * @param[in] phaseVisc phase viscosities in the sub region
   * @param[in] dPhaseVisc_dPres derivatives of phase viscosities wrt pressure in the sub region
   * @param[in] phaseRelPerm phase relative permeabilities in the sub region
   * @param[in] dPhaseRelPerm_dSat derivatives of the phase relative permeabilities wrt primary saturation in the sub
   * region
   * @param[inout] phaseMob phase mobilities in the sub region
   * @param[inout] dPhaseMob_dPres derivatives of phase mobilities wrt pressure in the sub region
   * @param[inout] dPhaseMob_dSat derivatives of phase mobilities wrt primary saturation in the sub region
   */
  static void
  Launch( localIndex const size,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dSat,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView2d< real64 > const & dPhaseMob_dSat );

};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  /**
   * @brief In an element, compute the accumulation term and its derivatives
   * @param[in] volume volume of the element
   * @param[in] porosityOld old porosity in the element
   * @param[in] porosityRef reference porosity in the element
   * @param[in] pvMult pore volume multiplier in the element
   * @param[in] dPvMult_dPres derivative of pore volume multiplier wrt pressure in the element
   * @param[in] phaseDensOld old phase densities in the element
   * @param[in] phaseDens phase densities in the element
   * @param[in] dPhaseDens_dPres derivatives of phase densities wrt pressure in the element
   * @param[inout] localAccum accumulation in the element
   * @param[inout] localAccumJacobian derivatives of accumulation term in the element
   */
  static void
  Compute( real64 const & volume,
           real64 const & porosityOld,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d< real64 const > const phaseSat,
           arraySlice1d< real64 const > const dPhaseSat,
           arraySlice1d< real64 const > const phaseDensOld,
           arraySlice1d< real64 const > const phaseDens,
           arraySlice1d< real64 const > const dPhaseDens_dPres,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian );

};


} // namespace TwoPhaseBaseKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASEKERNELS_HPP
