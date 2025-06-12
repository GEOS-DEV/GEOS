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
 * @file CompositionalMultiphaseFluidUtilities.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALMULTIPHASEFLUIDUTILITIES_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALMULTIPHASEFLUIDUTILITIES_HPP_

#include "constitutive/fluid/multifluid/compositional/models/NullModel.hpp"

namespace geos
{
namespace constitutive
{
namespace compositional
{

template< typename FLASH, typename PHASE1, typename PHASE2, typename PHASE3 >
struct CompositionalMultiphaseFluidTraits
{
  static constexpr integer getNumberOfPhases()
  {
    return FLASH::KernelWrapper::getNumberOfPhases();
  }

  static constexpr bool isThermalType()
  {
    if constexpr (getNumberOfPhases() == 3)
    {
      return !( std::is_same_v< typename PHASE1::Enthalpy, compositional::NullModel > ||
                std::is_same_v< typename PHASE2::Enthalpy, compositional::NullModel > ||
                std::is_same_v< typename PHASE3::Enthalpy, compositional::NullModel > );
    }
    else
    {
      return !( std::is_same_v< typename PHASE1::Enthalpy, compositional::NullModel > ||
                std::is_same_v< typename PHASE2::Enthalpy, compositional::NullModel > );
    }
  }
};

} //namespace compositional
} // namespace constitutive
} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_COMPOSITIONALMULTIPHASEFLUIDUTILITIES_HPP_
