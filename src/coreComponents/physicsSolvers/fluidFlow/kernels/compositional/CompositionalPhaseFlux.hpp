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



#pragma once

#include "codingUtilities/Utilities.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/FluxComputeKernelBase.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/PPUPhaseFlux.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/C1PPUPhaseFlux.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/IHUPhaseFlux.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/HU2PhaseFlux.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseFVMKernelUtilities
{

using KernelFlags = isothermalCompositionalMultiphaseFVMKernels::KernelFlags;

/**
 * @brief A stateless, static dispatcher for selecting a compositional phase flux discretization.
 */
struct CompositionalPhaseFlux
{
  /**
   * @brief Selects and calls the appropriate static compute function based on kernelFlags.
   *
   * @tparam NumComp The number of components.
   * @tparam NumSupportPoints The number of support points for the flux.
   * @param kernelFlags The runtime flags used to select the compute kernel.
   * @param ...args The full list of arguments to be forwarded to the compute kernels.
   */
  template<int NumComp, int NumSupportPoints>
  static GEOS_HOST_DEVICE inline void dispatch(
                                               
                                               BitFlags< KernelFlags > const & kernelFlags,
                                               // --- Pass-through arguments for the compute kernels ---
                                               const integer numPhases,
                                               const integer ip,
                                               localIndex const ( &seri )[NumSupportPoints],
                                               localIndex const ( &sesri )[NumSupportPoints],
                                               localIndex const ( &sei )[NumSupportPoints],
                                               real64 const ( &trans )[2],
                                               real64 const ( &dTrans_dPres )[2],
                                               ElementViewConst< arrayView1d< real64 const > > const & pres,
                                               ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                               ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                               ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob,
                                               ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseVolFrac,
                                               ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseVolFrac,
                                               ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                                               ElementViewConst< arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > > const & phaseMassDens,
                                               ElementViewConst< arrayView4d< real64 const, constitutive::multifluid::USD_PHASE_DC > > const & dPhaseMassDens,
                                               ElementViewConst< arrayView3d< real64 const, constitutive::cappres::USD_CAPPRES > > const & phaseCapPressure,
                                               ElementViewConst< arrayView4d< real64 const, constitutive::cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
                                               // --- Output Arguments ---
                                               real64 & potGrad,
                                               real64 & phaseFlux,
                                               real64 ( & dPhaseFlux_dP )[NumSupportPoints],
                                               real64 ( & dPhaseFlux_dC )[NumSupportPoints][NumComp],
                                               real64 & dPhaseFlux_dTrans
                                               )
  {
    // The dispatcher reads the flags and passes the correct boolean arguments down.
    const bool hasCapPressure = kernelFlags.isSet(KernelFlags::CapPressure);
    const bool checkPhasePresenceInGravity = kernelFlags.isSet(KernelFlags::CheckPhasePresenceInGravity);
    
    if (kernelFlags.isSet(KernelFlags::C1PPU))
    {
      C1PPUPhaseFlux::compute<NumComp, NumSupportPoints>(
                                                         numPhases, ip, hasCapPressure, checkPhasePresenceInGravity, seri, sesri, sei, trans, dTrans_dPres,
                                                         pres, gravCoef, phaseMob, dPhaseMob, phaseVolFrac, dPhaseVolFrac, dCompFrac_dCompDens,
                                                         phaseMassDens, dPhaseMassDens, phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                                                         potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC);
    }
    else if (kernelFlags.isSet(KernelFlags::IHU))
    {
      IHUPhaseFlux::compute<NumComp, NumSupportPoints>(
                                                       numPhases, ip, hasCapPressure, checkPhasePresenceInGravity, seri, sesri, sei, trans, dTrans_dPres,
                                                       pres, gravCoef, phaseMob, dPhaseMob, phaseVolFrac, dPhaseVolFrac, dCompFrac_dCompDens,
                                                       phaseMassDens, dPhaseMassDens, phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                                                       potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC);
    }
    else if (kernelFlags.isSet(KernelFlags::HU2PH))
    {
      HU2PhaseFlux::compute<NumComp, NumSupportPoints>(
                                                       numPhases, ip, hasCapPressure, checkPhasePresenceInGravity, seri, sesri, sei, trans, dTrans_dPres,
                                                       pres, gravCoef, phaseMob, dPhaseMob, phaseVolFrac, dPhaseVolFrac, dCompFrac_dCompDens,
                                                       phaseMassDens, dPhaseMassDens, phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                                                       potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC);
    }
    else
    {
      PPUPhaseFlux::compute<NumComp, NumSupportPoints>(
                                                       numPhases, ip, hasCapPressure, checkPhasePresenceInGravity, seri, sesri, sei, trans, dTrans_dPres,
                                                       pres, gravCoef, phaseMob, dPhaseMob, phaseVolFrac, dPhaseVolFrac, dCompFrac_dCompDens,
                                                       phaseMassDens, dPhaseMassDens, phaseCapPressure, dPhaseCapPressure_dPhaseVolFrac,
                                                       potGrad, phaseFlux, dPhaseFlux_dP, dPhaseFlux_dC, dPhaseFlux_dTrans);
    }
  }
};

} // namespace isothermalCompositionalMultiPhaseFVMKernelUtilities

} // namespace geos
