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
 * @file RelativePermeabilityInterpolators.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYINTERPOLATORS_HPP
#define GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYINTERPOLATORS_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"

namespace geosx
{

namespace constitutive
{

namespace relpermInterpolators
{

struct Baker
{

  /**
   * @brief Interpolate the two-phase relperms to compute the three-phase relperm
   * @param[in] shiftedWaterVolFrac
   * @param[in] gasVolFrac
   * @param[out] threePhaseRelPerm
   * @param[out] dThreePhaseRelPerm_dVolFrac
   * @param[in] relPerm_wo
   * @param[in] dRelPerm_wo_dOilVolFrac
   * @param[in] relPerm_go
   * @param[in] dRelPerm_go_dOilVolFrac
   *
   * This function interpolates the two-phase relperms to compute the three-phase relperm
   * The interpolation is based on the modified Baker method, also used as default in Eclipse
   * Reference: Eclipse technical description and PetroWiki
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( real64 const & shiftedWaterVolFrac,
           real64 const & gasVolFrac,
           arraySlice1d< integer const > const & phaseOrder,
           real64 const & relPerm_wo,
           real64 const & dRelPerm_wo_dOilVolFrac,
           real64 const & relPerm_go,
           real64 const & dRelPerm_go_dOilVolFrac,
           real64 & threePhaseRelPerm,
           arraySlice1d< real64 > const & dThreePhaseRelPerm_dVolFrac )
  {
    integer const ip_water = phaseOrder[RelativePermeabilityBase::PhaseType::WATER];
    integer const ip_oil   = phaseOrder[RelativePermeabilityBase::PhaseType::OIL];
    integer const ip_gas   = phaseOrder[RelativePermeabilityBase::PhaseType::GAS];

    // if water phase is immobile, then use the two-phase gas-oil data only
    if( shiftedWaterVolFrac < LvArray::NumericLimits< real64 >::epsilon )
    {
      threePhaseRelPerm = relPerm_go;
      dThreePhaseRelPerm_dVolFrac[ip_oil] = dRelPerm_go_dOilVolFrac;
    }
    // if gas phase is immobile, then use the two-phase water-oil data only
    else if( gasVolFrac < LvArray::NumericLimits< real64 >::epsilon )
    {
      threePhaseRelPerm = relPerm_wo;
      dThreePhaseRelPerm_dVolFrac[ip_oil] = dRelPerm_wo_dOilVolFrac;
    }
    // if both the water phase and the gas phase are mobile,
    // then use a saturation-weighted interpolation of the two-phase oil rel perms
    else
    {
      real64 const sumRelPerm = (shiftedWaterVolFrac * relPerm_wo
                                 + gasVolFrac   * relPerm_go);
      real64 const dSumRelPerm_dWaterVolFrac = relPerm_wo;
      real64 const dSumRelPerm_dOilVolFrac   = shiftedWaterVolFrac * dRelPerm_wo_dOilVolFrac
                                               + gasVolFrac   * dRelPerm_go_dOilVolFrac;
      real64 const dSumRelPerm_dGasVolFrac   = relPerm_go;


      real64 const sumVolFrac    = shiftedWaterVolFrac + gasVolFrac;
      real64 const sumVolFracInv = 1 / sumVolFrac; // div by 0 handled by the if statement above
      real64 const dSumVolFracInv_dWaterVolFrac = -sumVolFracInv * sumVolFracInv;
      real64 const dSumVolFracInv_dGasVolFrac   = dSumVolFracInv_dWaterVolFrac;

      threePhaseRelPerm = sumRelPerm * sumVolFracInv; // three-phase oil rel perm
      dThreePhaseRelPerm_dVolFrac[ip_water] = dSumRelPerm_dWaterVolFrac * sumVolFracInv // derivative w.r.t. Sw
                                              + sumRelPerm                * dSumVolFracInv_dWaterVolFrac;
      dThreePhaseRelPerm_dVolFrac[ip_oil]   = dSumRelPerm_dOilVolFrac   * sumVolFracInv;// derivative w.r.t. So
      dThreePhaseRelPerm_dVolFrac[ip_gas]   = dSumRelPerm_dGasVolFrac   * sumVolFracInv// derivative w.r.t. Sg
                                              + sumRelPerm                * dSumVolFracInv_dGasVolFrac;
    }
  }

};

} // namespace relpermInterpolators

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITYINTERPOLATORS_HPP
