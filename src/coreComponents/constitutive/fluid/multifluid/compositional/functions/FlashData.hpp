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
 * @file FlashData.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FLASHDATA_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FLASHDATA_HPP_

#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

/**
 * @brief A POD class to hold data for the negative two phase flash
 * @details This ties all the data that is required for all types of flash
 */
struct FlashData
{
  /// The equation of state for the liquid phase
  EquationOfStateType liquidEos;

  /// The equation of state for the vapour phase
  EquationOfStateType vapourEos;

  /// Salinity for the aqueous phase
  real64 salinity{0.0};
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_FLASHDATA_HPP_
